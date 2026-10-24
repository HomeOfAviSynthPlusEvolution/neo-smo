#include "algorithms/f3kdb_views.hpp"
#include "kernels/f3kdb.hpp"
#include <memory>

namespace neo_smo::f3kdb {
using namespace detail;
namespace {
// Select the sampling layout once per row, not once per pixel/reference.
// Modes 1/3 only use two neighbours; their other references are the center.
template <int Mode>
void collect_row(const double* center, const double* gx, const double* gy, const Distance* distances, std::size_t w,
                 std::size_t h, std::size_t y, int sx, int sy, const std::array<std::int64_t, 256>& vertical_offsets,
                 const std::array<double*, 4>& samples, const std::array<double*, 4>& grad_x,
                 const std::array<double*, 4>& grad_y, std::size_t start = 0) {
  const auto vertical_limit = std::min(y, h - 1 - y);
  const auto step = std::size_t{1} << sx;
  center += start;
  if constexpr (Mode == 7) {
    gx += start;
    gy += start;
  }
  distances += start * step;
  for (std::size_t x = start; x < w; ++x, distances += step, ++center) {
    const std::size_t h1 = distances->first >> sx, v1 = distances->first >> sy;
    const auto horizontal_limit = std::min(x, w - 1 - x);
    if ((Mode != 3 && v1 > vertical_limit) || (Mode != 1 && h1 > horizontal_limit))
      throw std::logic_error("reference outside plane");
    if constexpr (Mode == 1 || Mode == 3) {
      const auto offset = Mode == 1 ? vertical_offsets[v1] : static_cast<std::ptrdiff_t>(h1);
      samples[0][x] = center[offset];
      samples[1][x] = center[-offset];
    } else {
      std::array<std::ptrdiff_t, 4> offsets;
      const auto vertical = static_cast<std::ptrdiff_t>(vertical_offsets[v1]);
      const auto horizontal = static_cast<std::ptrdiff_t>(h1);
      if constexpr (Mode == 2) {
        const std::size_t h2 = distances->second >> sx, v2 = distances->second >> sy;
        if (h2 > horizontal_limit || v2 > vertical_limit)
          throw std::logic_error("reference outside plane");
        const auto a = static_cast<std::ptrdiff_t>(vertical_offsets[v2]) + horizontal;
        const auto b = static_cast<std::ptrdiff_t>(h2) - vertical;
        offsets = {a, b, -a, -b};
      } else {
        offsets = {vertical, -vertical, horizontal, -horizontal};
      }
      for (int j = 0; j < 4; ++j) {
        samples[j][x] = center[offsets[j]];
        if constexpr (Mode == 7) {
          grad_x[j][x] = gx[offsets[j]];
          grad_y[j][x] = gy[offsets[j]];
        }
      }
      if constexpr (Mode == 7) {
        ++gx;
        ++gy;
      }
    }
  }
}
// Modes 1-6 keep source samples in exact 32-bit integer lanes.
// Weighted mode 6 promotes its local references for double-precision weights.
Frame process_integer_frame(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input) {
  const auto& params = filter.parameters();
  const auto& f = params.input;
  const auto& out = params.output;
  const auto& o = params.options;
  const std::size_t planes = f.color == Color::gray ? 1 : 3;
  Frame result{out, {}};
  for (std::size_t plane = 0; plane < planes; ++plane) {
    const int sx = plane ? f.sub_x : 0, sy = plane ? f.sub_y : 0;
    const auto w = f.width >> sx, h = f.height >> sy, n = mul(w, h);
    const auto row_bytes = mul(w, bytes_per_sample(f));
    const auto& view = input[plane];
    if (!view.data || view.stride < row_bytes || view.bytes < add(mul(h - 1, view.stride), row_bytes))
      invalid("plane storage/stride");
    capacity<std::int32_t>(n);
    // Every element is written and validated before sampling; zeroing this
    // full-plane temporary would repeat the unpacker's writes.
    std::unique_ptr<std::int32_t[]> source(new std::int32_t[n]);
    if (!unpack_integer_plane({static_cast<const std::uint8_t*>(view.data), source.get(), w, h, view.stride, f.depth}))
      invalid("integer pixel exceeds declared depth");
    Plane output{w, h, mul(w, bytes_per_sample(out)), {}};
    output.data.resize(mul(output.stride, h));
    const bool use_fs = out.sample == Sample::integer && o.dither_algo == 3;
    const auto stripe_height = use_fs ? std::min<std::size_t>(32, h) : 1;
    std::vector<double> processed(mul(w, stripe_height));
    std::vector<float> boundary(use_fs ? w : 0, 0);
    // Distinct integer codes differ by at least one input code step. If the
    // maximum-neighbour threshold cannot admit that step, all accepted
    // neighbours equal the center and debanding is exactly the identity.
    const bool identity =
        (o.sample_mode == 5 || o.sample_mode == 6) && params.thresholds[plane][1] <= std::ldexp(1.0, 16 - f.depth);
    const auto pitch = mul(add(w, 15) / 16, 16);
    const auto& grain = plane ? filter.grain_c_integer() : filter.grain_y_integer();
    auto center = plane_view(static_cast<const std::int32_t*>(source.get()), w, h).cursor();
    auto distances = plane_view(filter.distances().data(), f.width, f.height).cursor();
    auto noise = plane_view(grain.data() + filter.offsets()[frame], w, h, mul(pitch, sizeof(std::int16_t))).cursor();
    auto target = plane_view(output.data.data(), output.stride, h, output.stride).cursor();
    std::size_t stripe_rows = 0;
    auto* stripe_output = output.data.data();
    auto* processed_row = processed.data();
    for (std::size_t y = 0; y < h; ++y) {
      if (y) {
        ++center;
        ++noise;
        ++target;
        distances += std::size_t{1} << sy;
      }
      if (stripe_rows == 0) {
        stripe_output = target.ptr();
        processed_row = processed.data();
      }
      const auto vlimit = std::min(y, h - 1 - y);
      IntegerRow row{center.ptr(), distances.ptr(), noise.ptr(), processed_row, w, vlimit, y * w, sx, sy};
      const bool direct = out.sample == Sample::integer && o.dither_algo == 1;
      if (direct) {
        row.packed_output = target.ptr();
        row.output_depth = out.depth;
        row.chroma = plane != 0;
        row.limited = o.keep_tv_range;
      }
      const auto count =
          process_integer_row(row, identity ? 0 : o.sample_mode, params.thresholds[plane].data(), o.blur_first);
      for (std::size_t x = count; x < w; ++x) {
        const auto d = distances.ptr()[x << sx];
        const auto h1 = std::size_t(d.first >> sx), v1 = std::size_t(d.first >> sy);
        const auto hlimit = std::min(x, w - 1 - x);
        if ((o.sample_mode != 3 && v1 > vlimit) || (o.sample_mode != 1 && h1 > hlimit))
          throw std::logic_error("reference outside plane");
        const auto* pixel = center.ptr() + x;
        std::array<double, 4> refs{};
        if (o.sample_mode == 1 || o.sample_mode == 3) {
          const auto offset = static_cast<std::ptrdiff_t>(o.sample_mode == 1 ? v1 * w : h1);
          refs[0] = pixel[offset];
          refs[1] = pixel[-offset];
        } else if (o.sample_mode == 2) {
          const auto h2 = std::size_t(d.second >> sx), v2 = std::size_t(d.second >> sy);
          if (h2 > hlimit || v2 > vlimit)
            throw std::logic_error("reference outside plane");
          const auto a = static_cast<std::ptrdiff_t>(v2 * w + h1);
          const auto b = static_cast<std::ptrdiff_t>(h2) - static_cast<std::ptrdiff_t>(v1 * w);
          refs = {double(pixel[a]), double(pixel[b]), double(pixel[-a]), double(pixel[-b])};
        }
        if (o.sample_mode >= 4) {
          const auto vertical = static_cast<std::ptrdiff_t>(v1 * w);
          const auto horizontal = static_cast<std::ptrdiff_t>(h1);
          refs = {double(pixel[vertical]), double(pixel[-vertical]), double(pixel[horizontal]),
                  double(pixel[-horizontal])};
        }
        processed_row[x] =
            deband(*pixel, refs, o.sample_mode, params.thresholds[plane], true, o.blur_first) + noise.ptr()[x];
      }
      ++stripe_rows;
      if (use_fs) {
        if (stripe_rows == stripe_height || y + 1 == h) {
          diffuse_stripe({processed.data(), stripe_output, boundary.data(), w, stripe_rows, output.stride, out.depth,
                          true, plane != 0, o.keep_tv_range});
          stripe_rows = 0;
        } else
          processed_row += w;
      } else {
        const auto start = direct ? count : 0;
        if (start < w)
          pack_row({processed_row + start, target.ptr() + start * bytes_per_sample(out), w - start, y, out.depth,
                    o.dither_algo, true, out.sample == Sample::integer, plane != 0, o.keep_tv_range});
        stripe_rows = 0;
      }
    }
    result.planes.push_back(std::move(output));
  }
  return result;
}
} // namespace

Frame process_highway(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input) {
  const auto& f = filter.parameters().input;
  const auto& out = filter.parameters().output;
  const auto& o = filter.parameters().options;
  const bool integer = f.sample == Sample::integer;
  const std::size_t planes = f.color == Color::gray ? 1 : 3;
  if (frame >= f.frames)
    invalid("frame index out of range");
  if (input.size() != planes)
    invalid("plane count");
  // Signed 32-bit gather indices must cover every source-plane displacement.
  if (integer && o.sample_mode <= 6 && f.sub_x <= 1 && f.sub_y <= 1 &&
      mul(f.width, f.height) <= std::size_t(std::numeric_limits<std::int32_t>::max()))
    return process_integer_frame(filter, frame, input);
  Frame result{out, {}};
  for (std::size_t plane = 0; plane < planes; ++plane) {
    const int sx = plane ? f.sub_x : 0, sy = plane ? f.sub_y : 0;
    const auto w = f.width >> sx, h = f.height >> sy, n = mul(w, h);
    const auto row = mul(w, bytes_per_sample(f));
    const auto& view = input[plane];
    if (!view.data || view.stride < row || view.bytes < add(mul(h - 1, view.stride), row))
      invalid("plane storage/stride");
    capacity<double>(n);
    std::vector<double> source(n);
    const auto* data = static_cast<const std::uint8_t*>(view.data);
    if (!unpack_plane({data, source.data(), w, h, view.stride, f.depth, integer}))
      invalid(integer ? "integer pixel exceeds declared depth" : "nonfinite input pixel");
    Plane output{w, h, mul(w, bytes_per_sample(out)), {}};
    output.data.resize(mul(output.stride, h));
    const bool use_fs = out.sample == Sample::integer && o.dither_algo == 3;
    // Deband fills a bounded stripe before wavefront diffusion. Carry raw signed
    // residuals across stripes so the result is independent of stripe height.
    const auto stripe_height = use_fs ? std::min<std::size_t>(32, h) : 1;
    std::vector<float> fs_boundary(use_fs ? w : 0, 0);
    const auto pitch = mul(add(w, 15) / 16, 16);
    const auto& grain = plane ? filter.grain_c() : filter.grain_y();
    const auto& grain_integer = plane ? filter.grain_c_integer() : filter.grain_y_integer();
    std::vector<double> gx, gy;
    if (o.sample_mode == 7) {
      gx.resize(n);
      gy.resize(n);
      sobel_plane(source.data(), gx.data(), gy.data(), w, h);
    }
    std::array<std::vector<double>, 4> samples, grad_x, grad_y;
    const int references = o.sample_mode == 1 || o.sample_mode == 3 ? 2 : 4;
    for (int j = 0; j < references; ++j) {
      samples[j].resize(w);
      if (o.sample_mode == 7) {
        grad_x[j].resize(w);
        grad_y[j].resize(w);
      }
    }
    std::vector<double> processed(mul(w, stripe_height));
    const auto source_view = plane_view(static_cast<const double*>(source.data()), w, h);
    const auto gx_view =
        o.sample_mode == 7 ? plane_view(static_cast<const double*>(gx.data()), w, h) : span2d::Plane<const double>{};
    const auto gy_view =
        o.sample_mode == 7 ? plane_view(static_cast<const double*>(gy.data()), w, h) : span2d::Plane<const double>{};
    const auto distance_view = plane_view(filter.distances().data(), f.width, f.height);
    const auto processed_view = plane_view(processed.data(), w, stripe_height);
    const auto output_view = plane_view(output.data.data(), output.stride, h, output.stride);
    // Distance stores bytes. Build the reachable vertical displacements once
    // by addition; every sampling mode and both gradient planes reuse them.
    std::array<std::int64_t, 256> vertical_offsets{};
    for (std::size_t i = 1; i <= std::min<std::size_t>(255, h - 1); ++i)
      vertical_offsets[i] = vertical_offsets[i - 1] + source_view.stride();
    auto center = source_view.cursor();
    auto distance_row = distance_view.cursor();
    auto gradient_x = o.sample_mode == 7 ? gx_view.cursor() : span2d::RowCursor<const double>{};
    auto gradient_y = o.sample_mode == 7 ? gy_view.cursor() : span2d::RowCursor<const double>{};
    auto processed_row = processed_view.cursor();
    std::size_t stripe_rows = 0;
    auto* stripe_output = output.data.data();
    auto noise = integer ? span2d::RowCursor<const float>{}
                         : plane_view(grain.data() + filter.offsets()[frame], w, h, mul(pitch, sizeof(float))).cursor();
    auto noise_integer =
        integer ? plane_view(grain_integer.data() + filter.offsets()[frame], w, h, mul(pitch, sizeof(std::int16_t)))
                      .cursor()
                : span2d::RowCursor<const std::int16_t>{};
    auto target = output_view.cursor();
    for (std::size_t y = 0; y < h; ++y) {
      if (y) {
        ++center;
        if (integer)
          ++noise_integer;
        else
          ++noise;
        ++target;
        distance_row += std::size_t{1} << sy;
        if (o.sample_mode == 7) {
          ++gradient_x;
          ++gradient_y;
        }
        if (stripe_rows == stripe_height) {
          stripe_rows = 0;
          processed_row = processed_view.cursor();
          stripe_output = target.ptr();
        } else {
          ++processed_row;
        }
      }
      const std::array<double*, 4> values = {samples[0].data(), samples[1].data(), samples[2].data(),
                                             samples[3].data()};
      const std::array<double*, 4> gradients_x = {grad_x[0].data(), grad_x[1].data(), grad_x[2].data(),
                                                  grad_x[3].data()};
      const std::array<double*, 4> gradients_y = {grad_y[0].data(), grad_y[1].data(), grad_y[2].data(),
                                                  grad_y[3].data()};
      const auto* distances = distance_row.ptr();
      Row batch{};
      batch.center = center.ptr();
      batch.count = w;
      batch.output = processed_row.ptr();
      batch.grain = noise.ptr();
      batch.integer_grain = noise_integer.ptr();
      std::size_t collected = 0;
      std::size_t processed_prefix = 0;
      if (o.sample_mode <= 3) {
        const SamplingRow sampling{center.ptr(),
                                   distances,
                                   vertical_offsets.data(),
                                   {values[0], values[1], values[2], values[3]},
                                   w,
                                   std::min(y, h - 1 - y),
                                   sx,
                                   sy};
        const auto prefix = process_sampled_row(sampling, batch, o.sample_mode,
                                                filter.parameters().thresholds[plane][0], integer, o.blur_first);
        collected = prefix.collected;
        processed_prefix = prefix.processed;
      }
      if (o.sample_mode == 7) {
        const CrossSamplingRow cross{{center.ptr(),
                                      distances,
                                      vertical_offsets.data(),
                                      {values[0], values[1], values[2], values[3]},
                                      w,
                                      std::min(y, h - 1 - y),
                                      sx,
                                      sy},
                                     gradient_x.ptr(),
                                     gradient_y.ptr(),
                                     {gradients_x[0], gradients_x[1], gradients_x[2], gradients_x[3]},
                                     {gradients_y[0], gradients_y[1], gradients_y[2], gradients_y[3]},
                                     y * w,
                                     integer};
        const auto prefix = process_cross_row(cross, batch, filter.parameters().thresholds[plane].data(), integer,
                                              o.max_angle, o.angle_boost);
        collected = prefix.collected;
        processed_prefix = prefix.processed;
      }
      // Keep the full scalar loop starting at a compile-time zero. A variable
      // start here inhibits optimization on targets that do not use gather.
      if (collected == 0) {
        switch (o.sample_mode) {
          case 1:
            collect_row<1>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y);
            break;
          case 2:
            collect_row<2>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y);
            break;
          case 3:
            collect_row<3>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y);
            break;
          case 7:
            collect_row<7>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y);
            break;
          default:
            collect_row<4>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y);
            break;
        }
      } else if (collected < w) {
        switch (o.sample_mode) {
          case 7:
            collect_row<7>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y, collected);
            break;
          case 1:
            collect_row<1>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y, collected);
            break;
          case 2:
            collect_row<2>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y, collected);
            break;
          case 3:
            collect_row<3>(center.ptr(), gradient_x.ptr(), gradient_y.ptr(), distances, w, h, y, sx, sy,
                           vertical_offsets, values, gradients_x, gradients_y, collected);
            break;
        }
      }

      if (o.sample_mode == 7) {
        batch.gx[0] = gradient_x.ptr();
        batch.gy[0] = gradient_y.ptr();
      }
      for (int j = 0; j < 4; ++j) {
        batch.refs[j] = j < references ? samples[j].data() : batch.center;
        batch.gx[j + 1] = grad_x[j].data();
        batch.gy[j + 1] = grad_y[j].data();
      }
      if (processed_prefix) {
        if (o.sample_mode == 7) {
          for (auto& grad : batch.gx)
            grad += processed_prefix;
          for (auto& grad : batch.gy)
            grad += processed_prefix;
        }
        batch.center += processed_prefix;
        batch.output += processed_prefix;
        if (integer)
          batch.integer_grain += processed_prefix;
        else
          batch.grain += processed_prefix;
        batch.count -= processed_prefix;
        for (auto& ref : batch.refs)
          ref += processed_prefix;
      }
      if (batch.count)
        process_row(batch, o.sample_mode, filter.parameters().thresholds[plane].data(), integer, o.blur_first,
                    o.max_angle, o.angle_boost);
      ++stripe_rows;
      if (use_fs) {
        const auto rows = stripe_rows;
        if (rows == stripe_height || y + 1 == h) {
          const FsStripe stripe{processed.data(), stripe_output, fs_boundary.data(), w, rows, output.stride, out.depth,
                                integer,          plane != 0,    o.keep_tv_range};
          diffuse_stripe(stripe);
        }
        continue;
      }
      pack_row({processed.data(), target.ptr(), w, y, out.depth, o.dither_algo, integer, out.sample == Sample::integer,
                plane != 0, o.keep_tv_range});
    }
    result.planes.push_back(std::move(output));
  }
  return result;
}
} // namespace neo_smo::f3kdb
