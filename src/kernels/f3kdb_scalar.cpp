#include "kernels/f3kdb_scalar.hpp"
#include "algorithms/f3kdb_internal.hpp"
#include <cstring>

namespace neo_smo::f3kdb {
using namespace detail;
namespace {
double plus(double a, double b, bool integer) {
  return integer ? a + b : fp(a + b);
}
double difference(double a, double b, bool integer) {
  return integer ? a - b : fp(a - b);
}
double mean(const double* a, int n, bool integer) {
  if (integer) {
    double sum = 0;
    for (int i = 0; i < n; ++i)
      sum += a[i];
    return std::floor((sum + n / 2) / n);
  }
  double sum = 0;
  // A widened weighted sum prevents overflow and intermediate underflow (e.g.
  // four identical minimum subnormals must retain their constant value).
  for (int i = 0; i < n; ++i)
    sum += a[i] / n;
  return fp(sum);
}
template <class T>
T load(const std::uint8_t* p) {
  T v;
  std::memcpy(&v, p, sizeof(v));
  return v;
}
template <class T>
void store(std::uint8_t* p, T v) {
  std::memcpy(p, &v, sizeof(v));
}
} // namespace
double round_even(double x) {
  if (!std::isfinite(x))
    invalid("nonfinite rounding input");
  if (std::abs(x) >= 0x1p52)
    return x;
  const double lo = std::floor(x), fraction = x - lo;
  return fraction < 0.5 ? lo : fraction > 0.5 ? lo + 1 : std::fmod(lo, 2) == 0 ? lo : lo + 1;
}
double response(double delta, double threshold) {
  if (threshold == 0)
    return delta == 0 ? 1 : 0;
  return std::clamp(3 * (1 - delta / threshold), 0.0, 1.0);
}
double weight(const std::array<double, 4>& c) {
  double sum = 0;
  for (double v : c) {
    bounded(v, 0, 1, "response");
    if (v == 0)
      return 0;
    sum += std::log(v);
  }
  return std::clamp(std::exp(sum / 10), 0.0, 1.0);
}
unsigned bayer16(std::size_t y, std::size_t x) {
  constexpr unsigned b[2][2] = {{0, 2}, {3, 1}};
  unsigned value = 0;
  for (unsigned j = 0; j < 4; ++j)
    value += b[(y >> j) & 1][(x >> j) & 1] << (2 * (3 - j));
  return value;
}
float decode_half(std::uint16_t h) {
  const int exp = (h >> 10) & 31, fraction = h & 1023;
  double v = exp == 0 ? std::ldexp(static_cast<double>(fraction), -24)
             : exp == 31
                 ? (fraction ? std::numeric_limits<double>::quiet_NaN() : std::numeric_limits<double>::infinity())
                 : std::ldexp(static_cast<double>(1024 + fraction), exp - 25);
  return static_cast<float>((h & 0x8000) ? -v : v);
}
std::uint16_t encode_half(double x) {
  if (!std::isfinite(x))
    invalid("nonfinite half input");
  const auto sign = static_cast<std::uint16_t>(std::signbit(x) ? 0x8000 : 0);
  x = std::min(std::abs(x), 65504.0);
  if (x < 0x1p-14)
    return static_cast<std::uint16_t>(sign | static_cast<unsigned>(round_even(std::ldexp(x, 24))));
  int exponent = 0;
  std::frexp(x, &exponent);
  auto mantissa = static_cast<unsigned>(round_even(std::ldexp(x, 11 - exponent)));
  if (mantissa == 2048) {
    mantissa = 1024;
    ++exponent;
  }
  return static_cast<std::uint16_t>(sign | ((exponent + 14) << 10) | (mantissa - 1024));
}
bool aligned(const std::array<double, 2>& g, const std::array<double, 2>& h, double angle) {
  const double gm = std::max(std::abs(g[0]), std::abs(g[1])), hm = std::max(std::abs(h[0]), std::abs(h[1]));
  if (gm == 0 || hm == 0)
    return gm == hm;
  if (angle >= 0.5)
    return true;
  // Integer gradients fit exactly in binary64 products. Normalized F64 vectors
  // also protect tiny F32 gradients and overflow-fallback gradients.
  const bool exact_integer = gm <= 262140 && hm <= 262140 && std::trunc(g[0]) == g[0] && std::trunc(g[1]) == g[1] &&
                             std::trunc(h[0]) == h[0] && std::trunc(h[1]) == h[1];
  const double gx = exact_integer ? g[0] : g[0] / gm, gy = exact_integer ? g[1] : g[1] / gm;
  const double hx = exact_integer ? h[0] : h[0] / hm, hy = exact_integer ? h[1] : h[1] / hm;
  const double c = std::abs(gx * hy - gy * hx), j = std::abs(gx * hx + gy * hy);
  if (angle == 0)
    return c == 0;
  if (angle == 0.25)
    return c <= j;
  constexpr double pi = 3.1415926535897932384626433832795;
  return c * std::cos(pi * angle) <= j * std::sin(pi * angle);
}
double deband(double s, const std::array<double, 4>& r, int mode, const std::array<double, 3>& thresholds, bool integer,
              bool blur, bool directions, double boost) {
  if (mode < 1 || mode > 7)
    invalid("sample mode");
  auto pass = [&](const double* a, int n) {
    const auto average = mean(a, n, integer);
    bool accepted = std::abs(difference(average, s, integer)) < thresholds[0];
    if (!blur) {
      accepted = true;
      for (int i = 0; i < n; ++i)
        accepted = accepted && std::abs(difference(a[i], s, integer)) < thresholds[0];
    }
    return accepted ? average : s;
  };
  if (mode <= 3)
    return pass(r.data(), mode == 2 ? 4 : 2);
  if (mode == 4) {
    const double results[] = {pass(r.data(), 2), pass(r.data() + 2, 2)};
    return mean(results, 2, integer);
  }
  const double a = mean(r.data(), 4, mode == 5 && integer);
  const double da = std::abs(difference(a, s, integer)),
               m = std::max({std::abs(difference(r[0], s, integer)), std::abs(difference(r[1], s, integer)),
                             std::abs(difference(r[2], s, integer)), std::abs(difference(r[3], s, integer))});
  const double qv = std::abs(difference(plus(r[0], r[1], integer), integer ? 2 * s : fp(2 * s), integer));
  const double qh = std::abs(difference(plus(r[2], r[3], integer), integer ? 2 * s : fp(2 * s), integer));
  if (mode == 5)
    return da < thresholds[0] && m < thresholds[1] && qv < thresholds[2] && qh < thresholds[2] ? a : s;
  const double k = mode == 7 && directions ? boost : 1;
  const double f = weight({response(da, thresholds[0] * k), response(m, thresholds[1] * k),
                           response(qv, thresholds[2] * k), response(qh, thresholds[2] * k)});
  if (f == 0)
    return s;
  const double result = integer ? s + (a - s) * f : fp(s + fp(difference(a, s, false) * f));
  return integer ? std::floor(result + 0.5) : result;
}

Frame process_scalar(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input) {
  const auto& params_ = filter.parameters();
  const auto& distances_ = filter.distances();
  const auto& offsets_ = filter.offsets();
  const auto& grain_y_ = filter.grain_y();
  const auto& grain_c_ = filter.grain_c();
  const auto& f = params_.input;
  const auto& out = params_.output;
  const auto& o = params_.options;
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
    for (std::size_t y = 0; y < h; ++y)
      for (std::size_t x = 0; x < w; ++x) {
        const auto* ptr = data + y * view.stride + x * bytes_per_sample(f);
        double v = integer         ? (f.depth == 8 ? *ptr : load<std::uint16_t>(ptr))
                   : f.depth == 16 ? decode_half(load<std::uint16_t>(ptr))
                                   : load<float>(ptr);
        if (!std::isfinite(v))
          invalid("nonfinite input pixel");
        if (integer) {
          if (v >= std::ldexp(1.0, f.depth))
            invalid("integer pixel exceeds declared depth");
          v = std::ldexp(v, 16 - f.depth);
        }
        source[y * w + x] = v;
      }
    Plane output{w, h, mul(w, bytes_per_sample(out)), {}};
    output.data.resize(mul(output.stride, h));
    // Store outgoing residuals; assemble incoming contributions in the contract's
    // fixed order (left, left-up, up, right-up), independent of scheduling.
    std::vector<float> previous(w, 0), current(w, 0);
    auto coordinate = [](std::size_t c, std::size_t extent, int offset) {
      if (offset < 0)
        return c < static_cast<std::size_t>(-offset) ? std::size_t{0} : c - static_cast<std::size_t>(-offset);
      return extent - 1 - c < static_cast<std::size_t>(offset) ? extent - 1 : c + static_cast<std::size_t>(offset);
    };
    auto gradient = [&](std::size_t x, std::size_t y) {
      auto at = [&](int dx, int dy) {
        return source[coordinate(y, h, dy) * w + coordinate(x, w, dx)];
      };
      // Compute exact F64 fallback from original samples if F32 arithmetic
      // overflows. Integer Sobel gradients remain exactly representable.
      const double gx = at(20, -20) + 2 * at(20, 0) + at(20, 20) - at(-20, -20) - 2 * at(-20, 0) - at(-20, 20);
      const double gy = at(-20, -20) + 2 * at(0, -20) + at(20, -20) - at(-20, 20) - 2 * at(0, 20) - at(20, 20);
      return std::array<double, 2>{gx, gy};
    };
    const auto pitch = mul(add(w, 15) / 16, 16);
    const auto& grain = plane ? grain_c_ : grain_y_;
    for (std::size_t y = 0; y < h; ++y) {
      std::fill(current.begin(), current.end(), 0);
      for (std::size_t x = 0; x < w; ++x) {
        const auto distance = distances_[(y << sy) * f.width + (x << sx)];
        const std::size_t h1 = distance.first >> sx, v1 = distance.first >> sy, h2 = distance.second >> sx,
                          v2 = distance.second >> sy;
        std::array<std::array<std::size_t, 2>, 4> positions;
        if (o.sample_mode == 1)
          positions = {{{x, y + v1}, {x, y - v1}, {x, y}, {x, y}}};
        else if (o.sample_mode == 3)
          positions = {{{x + h1, y}, {x - h1, y}, {x, y}, {x, y}}};
        else if (o.sample_mode == 2)
          positions = {{{x + h1, y + v2}, {x + h2, y - v1}, {x - h1, y - v2}, {x - h2, y + v1}}};
        else
          positions = {{{x, y + v1}, {x, y - v1}, {x + h1, y}, {x - h1, y}}};
        std::array<double, 4> refs{};
        bool directions = true;
        const auto g = o.sample_mode == 7 ? gradient(x, y) : std::array<double, 2>{};
        for (std::size_t j = 0; j < 4; ++j) {
          const auto px = positions[j][0], py = positions[j][1];
          if (px >= w || py >= h)
            throw std::logic_error("reference outside plane");
          refs[j] = source[py * w + px];
          if (o.sample_mode == 7)
            directions = directions && aligned(g, gradient(px, py), o.max_angle);
        }
        const double d = deband(source[y * w + x], refs, o.sample_mode, params_.thresholds[plane], integer,
                                o.blur_first, directions, o.angle_boost);
        const double z = plus(d, grain[offsets_[frame] + y * pitch + x], integer);
        auto* target = output.data.data() + y * output.stride + x * bytes_per_sample(out);
        if (out.sample == Sample::floating) {
          double v = integer ? z / 65536 - (plane ? 0.5 : 0) : z;
          if (o.keep_tv_range)
            v = std::clamp(v, 16.0 / 256 - (plane ? 0.5 : 0), (plane ? 240.0 : 235.0) / 256 - (plane ? 0.5 : 0));
          if (out.depth == 16)
            store(target, encode_half(v));
          else
            store(target, static_cast<float>(std::clamp(v, -static_cast<double>(std::numeric_limits<float>::max()),
                                                        static_cast<double>(std::numeric_limits<float>::max()))));
        } else {
          const double scale = std::ldexp(1.0, out.depth), tvscale = std::ldexp(1.0, out.depth - 8);
          const double lo = o.keep_tv_range ? 16 * tvscale : 0,
                       hi = o.keep_tv_range ? (plane ? 240 : 235) * tvscale : scale - 1;
          const double value = integer ? std::ldexp(z, out.depth - 16) : (z + (plane ? 0.5 : 0)) * scale;
          double quantized;
          if (o.dither_algo == 1)
            quantized = round_even(std::clamp(value, lo - 1, hi + 1));
          else if (o.dither_algo == 2)
            quantized = std::floor(std::clamp(value, lo - 1, hi + 1) + (bayer16(y, x) + 0.5) / 256);
          else {
            float error = x ? current[x - 1] * (7.0f / 16) : 0;
            if (y) {
              if (x)
                error += previous[x - 1] * (1.0f / 16);
              error += previous[x] * (5.0f / 16);
              if (x + 1 < w)
                error += previous[x + 1] * (3.0f / 16);
            }
            const double v = std::clamp(fp(std::clamp(value, lo, hi) + error), lo, hi);
            quantized = round_even(v);
            current[x] = static_cast<float>(v - quantized);
          }
          quantized = std::clamp(quantized, lo, hi);
          if (out.depth == 8)
            *target = static_cast<std::uint8_t>(quantized);
          else
            store(target, static_cast<std::uint16_t>(quantized));
        }
      }
      previous.swap(current);
    }
    result.planes.push_back(std::move(output));
  }
  return result;
}
} // namespace neo_smo::f3kdb
