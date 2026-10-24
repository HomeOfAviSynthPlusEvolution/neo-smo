#include "algorithms/f3kdb_views.hpp"
#include "kernels/f3kdb.hpp"
#include <memory>
#include <cstring>
#include <type_traits>

namespace neo_smo::f3kdb {
using namespace detail;
namespace {
// Select the sampling layout once per row, not once per pixel/reference.
// Modes 1/3 only use two neighbours; their other references are the center.
template <int Mode, class T>
void collect_row(const T* center, const T* gx, const T* gy, const Distance* distances, std::size_t w,
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
template <int Mode, int SX, int SY>
void collect_native8_fixed(const std::uint8_t* source, const Distance* distances, std::size_t w,
                           std::size_t source_stride, const std::array<std::uint8_t*, 4>& refs) {
  for (std::size_t x = 0; x < w; ++x) {
    const auto d = distances[x << SX];
    const std::size_t h1 = d.first >> SX, v1 = d.first >> SY;
    // Filter construction bounds immutable distances for every plane.
    const auto* center = source + x;
    std::ptrdiff_t offsets[4];
    if constexpr (Mode == 1 || Mode == 3) {
      const auto a = std::ptrdiff_t(Mode == 1 ? v1 * source_stride : h1);
      offsets[0] = a;
      offsets[1] = -a;
    } else {
      offsets[0] = std::ptrdiff_t(v1 * source_stride);
      offsets[1] = -offsets[0];
      offsets[2] = std::ptrdiff_t(h1);
      offsets[3] = -offsets[2];
    }
    for (int j = 0; j < (Mode == 1 || Mode == 3 ? 2 : 4); ++j)
      refs[j][x] = center[offsets[j]];
  }
}
template <int Mode>
void collect_byte_row(const ByteRow& row, const std::array<std::uint8_t*, 4>& refs) {
  if (row.sub_x == 0 && row.sub_y == 0)
    collect_native8_fixed<Mode, 0, 0>(row.source, row.distances, row.width, row.stride, refs);
  else if (row.sub_x == 1 && row.sub_y == 0)
    collect_native8_fixed<Mode, 1, 0>(row.source, row.distances, row.width, row.stride, refs);
  else if (row.sub_x == 0 && row.sub_y == 1)
    collect_native8_fixed<Mode, 0, 1>(row.source, row.distances, row.width, row.stride, refs);
  else
    collect_native8_fixed<Mode, 1, 1>(row.source, row.distances, row.width, row.stride, refs);
}

// Eight-bit round-to-even output can operate on borrowed source bytes. Wide
// formats, diffusion, and weighted/directional modes retain their existing path.
Frame process_byte_frame(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input,
                         const std::vector<MutablePlaneView>* external) {
  const auto& params = filter.parameters();
  const auto& f = params.input;
  const auto& out = params.output;
  const auto& o = params.options;
  const std::size_t planes = f.color == Color::gray ? 1 : 3;
  // Validate every borrowed plane before any caller-owned output is written.
  for (std::size_t plane = 0; plane < planes; ++plane) {
    const auto w = f.width >> (plane ? f.sub_x : 0), h = f.height >> (plane ? f.sub_y : 0);
    const auto& view = input[plane];
    const auto extent = add(mul(h - 1, view.stride), w);
    if (!view.data || view.stride < w || view.bytes < extent)
      invalid("plane storage/stride");
    if (extent > std::size_t(std::numeric_limits<std::ptrdiff_t>::max()))
      invalid("plane extent overflow");
  }
  Frame result{out, {}};
  for (std::size_t plane = 0; plane < planes; ++plane) {
    const int sx = plane ? f.sub_x : 0, sy = plane ? f.sub_y : 0;
    const auto w = f.width >> sx, h = f.height >> sy;
    const auto& view = input[plane];
    Plane output{w, h, w, {}};
    std::uint8_t* output_data;
    if (external) {
      output.stride = (*external)[plane].stride;
      output_data = static_cast<std::uint8_t*>((*external)[plane].data);
    } else {
      output.data.resize(mul(w, h));
      output_data = output.data.data();
    }
    const auto pitch = mul(add(w, 15) / 16, 16);
    const auto& grain = plane ? filter.grain_c_integer() : filter.grain_y_integer();
    std::vector<double> tail(w);
    std::array<std::vector<std::uint8_t>, 4> collected;
    const int refs_count = o.sample_mode == 2 ? 0 : (o.sample_mode == 4 ? 4 : 2);
    for (int j = 0; j < refs_count; ++j)
      collected[j].resize(w);
    const std::array<std::uint8_t*, 4> refs{collected[0].data(), collected[1].data(), collected[2].data(),
                                            collected[3].data()};
    for (std::size_t y = 0; y < h; ++y) {
      ByteRow row{static_cast<const std::uint8_t*>(view.data) + y * view.stride,
                  filter.distances().data() + (y << sy) * f.width,
                  grain.data() + filter.offsets()[frame] + y * pitch,
                  output_data + y * output.stride,
                  w,
                  view.stride,
                  sx,
                  sy,
                  plane != 0,
                  o.keep_tv_range,
                  {}};
      row.base = static_cast<const std::uint8_t*>(view.data);
      row.row_offset = y * view.stride;
      row.extent = (h - 1) * view.stride + w;
      auto count = process_byte_row(row, o.sample_mode, params.thresholds[plane].data(), o.blur_first);
      if (count == 0) {
        switch (o.sample_mode) {
          case 1:
            collect_byte_row<1>(row, refs);
            break;
          case 3:
            collect_byte_row<3>(row, refs);
            break;
          case 4:
            collect_byte_row<4>(row, refs);
            break;
        }
        for (int j = 0; j < refs_count; ++j)
          row.collected[j] = refs[j];
        count = process_byte_row(row, o.sample_mode, params.thresholds[plane].data(), o.blur_first);
      }
      for (std::size_t x = count; x < w; ++x) {
        const auto d = row.distances[x << sx];
        const auto h1 = std::ptrdiff_t(d.first >> sx), v1 = std::ptrdiff_t(d.first >> sy);
        const auto* center = row.source + x;
        auto pixel = [&](std::ptrdiff_t offset) {
          return double(center[offset]) * 256;
        };
        const auto stride = std::ptrdiff_t(view.stride);
        std::array<double, 4> samples{};
        if (o.sample_mode == 1 || o.sample_mode == 3) {
          const auto a = o.sample_mode == 1 ? v1 * stride : h1;
          samples[0] = pixel(a);
          samples[1] = pixel(-a);
        } else if (o.sample_mode == 2) {
          const auto h2 = std::ptrdiff_t(d.second >> sx), v2 = std::ptrdiff_t(d.second >> sy);
          const auto a = v2 * stride + h1, b = h2 - v1 * stride;
          samples = {pixel(a), pixel(b), pixel(-a), pixel(-b)};
        } else
          samples = {pixel(v1 * stride), pixel(-v1 * stride), pixel(h1), pixel(-h1)};
        tail[x] = deband(pixel(0), samples, o.sample_mode, params.thresholds[plane], true, o.blur_first) + row.grain[x];
      }
      if (count < w)
        pack_row(
            {tail.data() + count, row.output + count, w - count, y, 8, 1, true, true, plane != 0, o.keep_tv_range});
    }
    if (!external)
      result.planes.push_back(std::move(output));
  }
  return result;
}

// Modes 1-6 widen borrowed source samples into exact 32-bit integer lanes.
// Weighted mode 6 promotes its local references for double-precision weights.
Frame process_integer_frame(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input,
                            const std::vector<MutablePlaneView>* external) {
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
    const auto extent = add(mul(h - 1, view.stride), row_bytes);
    const bool native = row_bytes >= 4 && view.stride <= std::size_t(std::numeric_limits<std::int32_t>::max()) &&
                        extent <= std::size_t(std::numeric_limits<std::int32_t>::max());
    const auto stored = native ? w : n;
    capacity<std::int32_t>(stored);
    std::unique_ptr<std::int32_t[]> source(new std::int32_t[stored]);
    if (native) {
      // Eight/sixteen-bit storage cannot contain out-of-depth codes. Other
      // depths are validated with bounded row scratch, before writing output.
      if (f.depth != 8 && f.depth != 16)
        for (std::size_t y = 0; y < h; ++y)
          if (!unpack_integer_plane({static_cast<const std::uint8_t*>(view.data) + y * view.stride,
                                     source.get(), w, 1, view.stride, f.depth}))
            invalid("integer pixel exceeds declared depth");
    } else if (!unpack_integer_plane({static_cast<const std::uint8_t*>(view.data), source.get(), w, h,
                                      view.stride, f.depth}))
      invalid("integer pixel exceeds declared depth");
    Plane output{w, h, mul(w, bytes_per_sample(out)), {}};
    std::uint8_t* output_data;
    if (external) {
      output.stride = (*external)[plane].stride;
      output_data = static_cast<std::uint8_t*>((*external)[plane].data);
    } else {
      output.data.resize(mul(output.stride, h));
      output_data = output.data.data();
    }
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
    auto center = plane_view(static_cast<const std::int32_t*>(source.get()), w, native ? 1 : h).cursor();
    auto distances = plane_view(filter.distances().data(), f.width, f.height).cursor();
    auto noise = plane_view(grain.data() + filter.offsets()[frame], w, h, mul(pitch, sizeof(std::int16_t))).cursor();
    auto target = plane_view(output_data, output.stride, h, output.stride).cursor();
    std::size_t stripe_rows = 0;
    auto* stripe_output = output_data;
    auto* processed_row = processed.data();
    for (std::size_t y = 0; y < h; ++y) {
      if (y) {
        if (!native)
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
      if (native) {
        row.native = static_cast<const std::uint8_t*>(view.data);
        row.native_stride = view.stride;
        row.native_extent = extent;
        row.input_depth = f.depth;
        row.row_offset = y * view.stride;
      }
      const bool direct = out.sample == Sample::integer && o.dither_algo <= 2;
      if (direct) {
        row.packed_output = target.ptr();
        row.dither = o.dither_algo;
        row.y = y;
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
        const auto pixel = [&](std::ptrdiff_t offset) -> double {
          if (!native)
            return center.ptr()[std::ptrdiff_t(x) + offset];
          const auto* p = static_cast<const std::uint8_t*>(view.data) + y * view.stride + x * bytes_per_sample(f) + offset;
          std::uint16_t code = *p;
          if (f.depth != 8)
            std::memcpy(&code, p, sizeof(code));
          return std::int32_t(code) << (16 - f.depth);
        };
        const auto stride = native ? view.stride : w;
        const auto step = native ? bytes_per_sample(f) : 1;
        std::array<double, 4> refs{};
        if (o.sample_mode == 1 || o.sample_mode == 3) {
          const auto offset = static_cast<std::ptrdiff_t>(o.sample_mode == 1 ? v1 * stride : h1 * step);
          refs[0] = pixel(offset);
          refs[1] = pixel(-offset);
        } else if (o.sample_mode == 2) {
          const auto h2 = std::size_t(d.second >> sx), v2 = std::size_t(d.second >> sy);
          if (h2 > hlimit || v2 > vlimit)
            throw std::logic_error("reference outside plane");
          const auto a = static_cast<std::ptrdiff_t>(v2 * stride + h1 * step);
          const auto b = static_cast<std::ptrdiff_t>(h2 * step) - static_cast<std::ptrdiff_t>(v1 * stride);
          refs = {double(pixel(a)), double(pixel(b)), double(pixel(-a)), double(pixel(-b))};
        }
        if (o.sample_mode >= 4) {
          const auto vertical = static_cast<std::ptrdiff_t>(v1 * stride);
          const auto horizontal = static_cast<std::ptrdiff_t>(h1 * step);
          refs = {double(pixel(vertical)), double(pixel(-vertical)), double(pixel(horizontal)),
                  double(pixel(-horizontal))};
        }
        processed_row[x] =
            deband(pixel(0), refs, o.sample_mode, params.thresholds[plane], true, o.blur_first) + noise.ptr()[x];
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
                    o.dither_algo, true, out.sample == Sample::integer, plane != 0, o.keep_tv_range, start});
        stripe_rows = 0;
      }
    }
    if (!external)
      result.planes.push_back(std::move(output));
  }
  return result;
}
} // namespace

template <class T>
Frame process_highway_typed(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input,
                           const std::vector<MutablePlaneView>* external) {
  const auto& f = filter.parameters().input;
  const auto& out = filter.parameters().output;
  const auto& o = filter.parameters().options;
  const bool integer = f.sample == Sample::integer;
  const std::size_t planes = f.color == Color::gray ? 1 : 3;
  if (frame >= f.frames)
    invalid("frame index out of range");
  if (input.size() != planes)
    invalid("plane count");
  if (integer && f.depth == 8 && out.sample == Sample::integer && out.depth == 8 && o.dither_algo == 1 &&
      o.sample_mode <= 4 && f.sub_x <= 1 && f.sub_y <= 1)
    return process_byte_frame(filter, frame, input, external);
  // Signed 32-bit gather indices must cover every source-plane displacement.
  if (integer && o.sample_mode <= 6 && f.sub_x <= 1 && f.sub_y <= 1 &&
      mul(f.width, f.height) <= std::size_t(std::numeric_limits<std::int32_t>::max()))
    return process_integer_frame(filter, frame, input, external);
  Frame result{out, {}};
  for (std::size_t plane = 0; plane < planes; ++plane) {
    const int sx = plane ? f.sub_x : 0, sy = plane ? f.sub_y : 0;
    const auto w = f.width >> sx, h = f.height >> sy, n = mul(w, h);
    const auto row = mul(w, bytes_per_sample(f));
    const auto& view = input[plane];
    if (!view.data || view.stride < row || view.bytes < add(mul(h - 1, view.stride), row))
      invalid("plane storage/stride");
    capacity<T>(n);
    std::vector<T> source(n);
    const auto* data = static_cast<const std::uint8_t*>(view.data);
    if constexpr (std::is_same_v<T, std::int32_t>) {
      if (!unpack_integer_plane({data, source.data(), w, h, view.stride, f.depth}))
        invalid("integer pixel exceeds declared depth");
    } else if (!unpack_plane({data, source.data(), w, h, view.stride, f.depth, integer}))
      invalid(integer ? "integer pixel exceeds declared depth" : "nonfinite input pixel");
    Plane output{w, h, mul(w, bytes_per_sample(out)), {}};
    std::uint8_t* output_data;
    if (external) {
      output.stride = (*external)[plane].stride;
      output_data = static_cast<std::uint8_t*>((*external)[plane].data);
    } else {
      output.data.resize(mul(output.stride, h));
      output_data = output.data.data();
    }
    const bool use_fs = out.sample == Sample::integer && o.dither_algo == 3;
    // Deband fills a bounded stripe before wavefront diffusion. Carry raw signed
    // residuals across stripes so the result is independent of stripe height.
    const auto stripe_height = use_fs ? std::min<std::size_t>(32, h) : 1;
    std::vector<float> fs_boundary(use_fs ? w : 0, 0);
    const auto pitch = mul(add(w, 15) / 16, 16);
    const auto& grain = plane ? filter.grain_c() : filter.grain_y();
    const auto& grain_integer = plane ? filter.grain_c_integer() : filter.grain_y_integer();
    std::vector<T> gx, gy;
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
    std::vector<double> wide_center, wide_gx, wide_gy;
    if constexpr (std::is_same_v<T, std::int32_t>) {
      wide_center.resize(w);
      wide_gx.resize(w);
      wide_gy.resize(w);
    }
    const auto source_view = plane_view(static_cast<const T*>(source.data()), w, h);
    const auto gx_view =
        o.sample_mode == 7 ? plane_view(static_cast<const T*>(gx.data()), w, h) : span2d::Plane<const T>{};
    const auto gy_view =
        o.sample_mode == 7 ? plane_view(static_cast<const T*>(gy.data()), w, h) : span2d::Plane<const T>{};
    const auto distance_view = plane_view(filter.distances().data(), f.width, f.height);
    const auto processed_view = plane_view(processed.data(), w, stripe_height);
    const auto output_view = plane_view(output_data, output.stride, h, output.stride);
    // Distance stores bytes. Build the reachable vertical displacements once
    // by addition; every sampling mode and both gradient planes reuse them.
    std::array<std::int64_t, 256> vertical_offsets{};
    for (std::size_t i = 1; i <= std::min<std::size_t>(255, h - 1); ++i)
      vertical_offsets[i] = vertical_offsets[i - 1] + source_view.stride();
    auto center = source_view.cursor();
    auto distance_row = distance_view.cursor();
    auto gradient_x = o.sample_mode == 7 ? gx_view.cursor() : span2d::RowCursor<const T>{};
    auto gradient_y = o.sample_mode == 7 ? gy_view.cursor() : span2d::RowCursor<const T>{};
    auto processed_row = processed_view.cursor();
    std::size_t stripe_rows = 0;
    auto* stripe_output = output_data;
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
      const double* center_values;
      const double* gradient_values_x;
      const double* gradient_values_y;
      if constexpr (std::is_same_v<T, std::int32_t>) {
        std::copy_n(center.ptr(), w, wide_center.data());
        std::copy_n(gradient_x.ptr(), w, wide_gx.data());
        std::copy_n(gradient_y.ptr(), w, wide_gy.data());
        center_values = wide_center.data();
        gradient_values_x = wide_gx.data();
        gradient_values_y = wide_gy.data();
      } else {
        center_values = center.ptr();
        gradient_values_x = gradient_x.ptr();
        gradient_values_y = gradient_y.ptr();
      }
      Row batch{};
      batch.center = center_values;
      batch.count = w;
      batch.output = processed_row.ptr();
      batch.grain = noise.ptr();
      batch.integer_grain = noise_integer.ptr();
      std::size_t collected = 0;
      std::size_t processed_prefix = 0;
      if (o.sample_mode <= 3) {
        const SamplingRow sampling{center_values,
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
        CrossSamplingRow cross{{center_values,
                                      distances,
                                      vertical_offsets.data(),
                                      {values[0], values[1], values[2], values[3]},
                                      w,
                                      std::min(y, h - 1 - y),
                                      sx,
                                      sy},
                                     gradient_values_x,
                                     gradient_values_y,
                                     {gradients_x[0], gradients_x[1], gradients_x[2], gradients_x[3]},
                                     {gradients_y[0], gradients_y[1], gradients_y[2], gradients_y[3]},
                                     y * w,
                                     integer};
        if constexpr (std::is_same_v<T, std::int32_t>) {
          cross.compact_source = source.data();
          cross.compact_gx = gx.data();
          cross.compact_gy = gy.data();
        }
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
        batch.gx[0] = gradient_values_x;
        batch.gy[0] = gradient_values_y;
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
    if (!external)
      result.planes.push_back(std::move(output));
  }
  return result;
}
Frame process_highway_impl(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input,
                           const std::vector<MutablePlaneView>* external) {
  const auto& p = filter.parameters();
  // Integer Sobel components are bounded by 4*65535, so their full planes
  // need only signed 32-bit storage. Direction products still use double.
  if (p.input.sample == Sample::integer && p.options.sample_mode == 7 &&
      p.input.sub_x <= 1 && p.input.sub_y <= 1 &&
      mul(p.input.width, p.input.height) <= std::size_t(std::numeric_limits<std::int32_t>::max()))
    return process_highway_typed<std::int32_t>(filter, frame, input, external);
  return process_highway_typed<double>(filter, frame, input, external);
}
Frame process_highway(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input) {
  return process_highway_impl(filter, frame, input, nullptr);
}
void process_highway_into(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input,
                          const std::vector<MutablePlaneView>& output) {
  process_highway_impl(filter, frame, input, &output);
}
} // namespace neo_smo::f3kdb
