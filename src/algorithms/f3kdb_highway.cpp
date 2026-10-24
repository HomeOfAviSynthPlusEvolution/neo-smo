#include "algorithms/f3kdb_views.hpp"
#include "kernels/f3kdb.hpp"

namespace neo_smo::f3kdb {
using namespace detail;
namespace {
// Select the sampling layout once per row, not once per pixel/reference.
// Modes 1/3 only use two neighbours; their other references are the center.
template <int Mode>
void collect_row(span2d::RowCursor<const double> center, span2d::Plane<const double> source,
                 span2d::Plane<const double> gx, span2d::Plane<const double> gy, const Distance* distances,
                 std::size_t y, int sx, int sy, const std::array<double*, 4>& samples,
                 const std::array<double*, 4>& grad_x, const std::array<double*, 4>& grad_y) {
  const auto w = static_cast<std::size_t>(source.width()), h = static_cast<std::size_t>(source.height());
  // The two references in modes 1/3 are symmetric around the center.
  // Validate their shared distance once, then form one linear displacement.
  if constexpr (Mode == 1 || Mode == 3) {
    const auto vertical_limit = std::min(y, h - 1 - y);
    const auto step = std::size_t{1} << sx;
    auto* first = samples[0];
    auto* second = samples[1];
    for (std::size_t x = 0; x < w; ++x, distances += step) {
      const std::size_t distance = distances->first >> (Mode == 1 ? sy : sx);
      const auto limit = Mode == 1 ? vertical_limit : std::min(x, w - 1 - x);
      if (distance > limit)
        throw std::logic_error("reference outside plane");
      const auto offset = Mode == 1 ? distance * w : distance;
      const auto* pixel = center.ptr() + x;
      first[x] = pixel[offset];
      second[x] = *(pixel - offset);
    }
    return;
  }
  constexpr int count = Mode == 1 || Mode == 3 ? 2 : 4;
  for (std::size_t x = 0; x < w; ++x) {
    const auto distance = distances[x << sx];
    const std::size_t h1 = distance.first >> sx, v1 = distance.first >> sy;
    std::array<std::array<std::size_t, 2>, count> positions;
    if constexpr (Mode == 1)
      positions = {{{x, y + v1}, {x, y - v1}}};
    else if constexpr (Mode == 3)
      positions = {{{x + h1, y}, {x - h1, y}}};
    else if constexpr (Mode == 2) {
      const std::size_t h2 = distance.second >> sx, v2 = distance.second >> sy;
      positions = {{{x + h1, y + v2}, {x + h2, y - v1}, {x - h1, y - v2}, {x - h2, y + v1}}};
    } else
      positions = {{{x, y + v1}, {x, y - v1}, {x + h1, y}, {x - h1, y}}};
    for (int j = 0; j < count; ++j) {
      const auto px = positions[j][0], py = positions[j][1];
      if (px >= w || py >= h)
        throw std::logic_error("reference outside plane");
      samples[j][x] = source(py, px);
      if constexpr (Mode == 7) {
        grad_x[j][x] = gx(py, px);
        grad_y[j][x] = gy(py, px);
      }
    }
  }
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
    auto center = source_view.cursor();
    auto noise = plane_view(grain.data() + filter.offsets()[frame], w, h, mul(pitch, sizeof(double))).cursor();
    auto target = output_view.cursor();
    for (std::size_t y = 0; y < h; ++y) {
      if (y) {
        ++center;
        ++noise;
        ++target;
      }
      const std::array<double*, 4> values = {samples[0].data(), samples[1].data(), samples[2].data(),
                                             samples[3].data()};
      const std::array<double*, 4> gradients_x = {grad_x[0].data(), grad_x[1].data(), grad_x[2].data(),
                                                  grad_x[3].data()};
      const std::array<double*, 4> gradients_y = {grad_y[0].data(), grad_y[1].data(), grad_y[2].data(),
                                                  grad_y[3].data()};
      const auto* distances = distance_view.row_ptr(y << sy);
      switch (o.sample_mode) {
        case 1:
          collect_row<1>(center, source_view, gx_view, gy_view, distances, y, sx, sy, values, gradients_x, gradients_y);
          break;
        case 2:
          collect_row<2>(center, source_view, gx_view, gy_view, distances, y, sx, sy, values, gradients_x, gradients_y);
          break;
        case 3:
          collect_row<3>(center, source_view, gx_view, gy_view, distances, y, sx, sy, values, gradients_x, gradients_y);
          break;
        case 7:
          collect_row<7>(center, source_view, gx_view, gy_view, distances, y, sx, sy, values, gradients_x, gradients_y);
          break;
        default:
          collect_row<4>(center, source_view, gx_view, gy_view, distances, y, sx, sy, values, gradients_x, gradients_y);
          break;
      }
      Row batch{};
      batch.center = center.ptr();
      batch.count = w;
      batch.output = processed_view.row_ptr(y % stripe_height);
      batch.grain = noise.ptr();
      if (o.sample_mode == 7) {
        batch.gx[0] = gx_view.row_ptr(y);
        batch.gy[0] = gy_view.row_ptr(y);
      }
      for (int j = 0; j < 4; ++j) {
        batch.refs[j] = j < references ? samples[j].data() : batch.center;
        batch.gx[j + 1] = grad_x[j].data();
        batch.gy[j + 1] = grad_y[j].data();
      }
      process_row(batch, o.sample_mode, filter.parameters().thresholds[plane].data(), integer, o.blur_first,
                  o.max_angle, o.angle_boost);
      if (use_fs) {
        const auto rows = y % stripe_height + 1;
        if (rows == stripe_height || y + 1 == h) {
          const FsStripe stripe{processed.data(),
                                output_view.row_ptr(y + 1 - rows),
                                fs_boundary.data(),
                                w,
                                rows,
                                output.stride,
                                out.depth,
                                integer,
                                plane != 0,
                                o.keep_tv_range};
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
