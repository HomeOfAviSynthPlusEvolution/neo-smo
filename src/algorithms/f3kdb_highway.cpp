#include "algorithms/f3kdb_internal.hpp"
#include "kernels/f3kdb.hpp"

namespace neo_smo::f3kdb {
using namespace detail;
namespace {
// Select the sampling layout once per row, not once per pixel/reference.
// Modes 1/3 only use two neighbours; their other references are the center.
template <int Mode>
void collect_row(const double* source, const double* gx, const double* gy, const Distance* distances, std::size_t w,
                 std::size_t h, std::size_t y, int sx, int sy, const std::array<double*, 4>& samples,
                 const std::array<double*, 4>& grad_x, const std::array<double*, 4>& grad_y) {
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
      const auto index = py * w + px;
      samples[j][x] = source[index];
      if constexpr (Mode == 7) {
        grad_x[j][x] = gx[index];
        grad_y[j][x] = gy[index];
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
    for (std::size_t y = 0; y < h; ++y) {
      const std::array<double*, 4> values = {samples[0].data(), samples[1].data(), samples[2].data(),
                                             samples[3].data()};
      const std::array<double*, 4> gradients_x = {grad_x[0].data(), grad_x[1].data(), grad_x[2].data(),
                                                  grad_x[3].data()};
      const std::array<double*, 4> gradients_y = {grad_y[0].data(), grad_y[1].data(), grad_y[2].data(),
                                                  grad_y[3].data()};
      const auto* distances = filter.distances().data() + (y << sy) * f.width;
      switch (o.sample_mode) {
        case 1:
          collect_row<1>(source.data(), gx.data(), gy.data(), distances, w, h, y, sx, sy, values, gradients_x,
                         gradients_y);
          break;
        case 2:
          collect_row<2>(source.data(), gx.data(), gy.data(), distances, w, h, y, sx, sy, values, gradients_x,
                         gradients_y);
          break;
        case 3:
          collect_row<3>(source.data(), gx.data(), gy.data(), distances, w, h, y, sx, sy, values, gradients_x,
                         gradients_y);
          break;
        case 7:
          collect_row<7>(source.data(), gx.data(), gy.data(), distances, w, h, y, sx, sy, values, gradients_x,
                         gradients_y);
          break;
        default:
          collect_row<4>(source.data(), gx.data(), gy.data(), distances, w, h, y, sx, sy, values, gradients_x,
                         gradients_y);
          break;
      }
      Row batch{};
      batch.center = source.data() + y * w;
      batch.count = w;
      batch.output = processed.data() + (y % stripe_height) * w;
      batch.grain = grain.data() + filter.offsets()[frame] + y * pitch;
      if (o.sample_mode == 7) {
        batch.gx[0] = gx.data() + y * w;
        batch.gy[0] = gy.data() + y * w;
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
                                output.data.data() + (y + 1 - rows) * output.stride,
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
      pack_row({processed.data(), output.data.data() + y * output.stride, w, y, out.depth, o.dither_algo, integer,
                out.sample == Sample::integer, plane != 0, o.keep_tv_range});
    }
    result.planes.push_back(std::move(output));
  }
  return result;
}
} // namespace neo_smo::f3kdb
