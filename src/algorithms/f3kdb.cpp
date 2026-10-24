#include "algorithms/f3kdb_internal.hpp"
#include "kernels/f3kdb.hpp"

namespace neo_smo::f3kdb {
using namespace detail;
Parameters normalize(const Format& input, const Options& o) {
  valid_format(input);
  Parameters p{input, input, o, {}, {256, 256}};
  if (o.output_sample_type == OutputSample::integer)
    p.output.sample = Sample::integer;
  else if (o.output_sample_type == OutputSample::floating)
    p.output.sample = Sample::floating;
  else if (o.output_sample_type != OutputSample::same)
    invalid("output_sample_type");
  if (o.output_depth < -1)
    invalid("output_depth");
  p.output.depth = o.output_depth != -1                 ? o.output_depth
                   : p.output.sample == input.sample    ? input.depth
                   : p.output.sample == Sample::integer ? 16
                                                        : 32;
  valid_format(p.output);
  if (o.range < 0 || o.range > 255)
    invalid("range");
  if (o.sample_mode < 1 || o.sample_mode > 7)
    invalid("sample_mode");
  if (o.dither_algo < 1 || o.dither_algo > 3)
    invalid("dither_algo");
  if (o.random_algo_ref < 0 || o.random_algo_ref > 2 || o.random_algo_grain < 0 || o.random_algo_grain > 2)
    invalid("random algorithm");
  bounded(o.random_param_ref, 0, std::numeric_limits<double>::max(), "random_param_ref");
  bounded(o.random_param_grain, 0, std::numeric_limits<double>::max(), "random_param_grain");
  bounded(o.angle_boost, 0, std::numeric_limits<double>::max(), "angle_boost");
  bounded(o.max_angle, 0, 1, "max_angle");
  for (auto& t : p.thresholds)
    t = {64, -1, -1};
  p.grain = {64, 64};
  for (const auto& list : o.presets) {
    std::size_t start = 0;
    do {
      const auto end = list.find('/', start);
      const auto token = list.substr(start, end == std::string::npos ? end : end - start);
      double strength = -1;
      if (token == "low")
        strength = 128;
      else if (token == "medium")
        strength = 192;
      else if (token == "high")
        strength = 256;
      else if (token == "veryhigh")
        strength = 320;
      else if (token == "depth")
        strength = 0;
      if (strength >= 0) {
        for (auto& t : p.thresholds)
          t.fill(strength / (o.scale ? 1 : 4));
        p.grain.fill(strength / 4);
      } else if (token == "nograin")
        p.grain.fill(0);
      else if (token == "luma") {
        p.thresholds[1].fill(0);
        p.thresholds[2].fill(0);
        p.grain[1] = 0;
      } else if (token == "chroma") {
        p.thresholds[0].fill(0);
        p.grain[0] = 0;
      } else if (!token.empty())
        invalid("unknown preset");
      if (end == std::string::npos)
        break;
      start = end + 1;
    } while (true);
  }
  for (int plane = 0; plane < 3; ++plane) {
    auto& t = p.thresholds[plane];
    if (o.threshold[plane])
      t[0] = *o.threshold[plane];
    if (o.threshold1[plane])
      t[1] = *o.threshold1[plane];
    if (o.threshold2[plane])
      t[2] = *o.threshold2[plane];
    for (int j = 1; j < 3; ++j)
      if (t[j] == -1)
        t[j] = t[0];
    for (int j = 0; j < 3; ++j) {
      bounded(t[j], 0, o.scale ? 65535 : 511, "threshold");
      t[j] *= o.scale ? 1 : 4;
      if (input.sample == Sample::floating)
        t[j] /= 65536;
      if (!std::isfinite(t[j] * o.angle_boost))
        invalid("boosted threshold overflow");
    }
  }
  if (o.grain_y)
    p.grain[0] = *o.grain_y;
  if (o.grain_c)
    p.grain[1] = *o.grain_c;
  for (auto& g : p.grain) {
    bounded(g, 0, 4096, "grain");
    g *= 4;
    if (input.sample == Sample::floating)
      g /= 65536;
  }
  return p;
}

double Random::next(int algorithm, double sigma) {
  if (algorithm < 0 || algorithm > 2)
    invalid("random algorithm");
  bounded(sigma, 0, std::numeric_limits<double>::max(), "random sigma");
  if (algorithm == 0) {
    state_ ^= state_ << 13;
    state_ ^= state_ >> 17;
    state_ ^= state_ << 5;
    return uniform_value(state_);
  }
  if (algorithm == 1) {
    state_ = 1664525u * state_ + 1013904223u;
    return uniform_value(state_);
  }
  if (sigma == 0) {
    next(1, 0);
    next(1, 0);
    return 0;
  }
  for (int attempt = 0; attempt < 4096; ++attempt) {
    const double x = next(1, 0), y = next(1, 0);
    if (sigma >= 1) {
      if (x > -1 && x < 1 && (y + 1) / 2 < std::exp(-0.5 * (x / sigma) * (x / sigma)))
        return x;
    } else {
      const double r2 = x * x + y * y;
      if (r2 == 0 || r2 > 1)
        continue;
      const double z = sigma * y * std::sqrt(-2 * std::log(r2) / r2);
      if (z > -1 && z < 1)
        return z;
    }
  }
  throw std::runtime_error("Gaussian initialization exceeded 4096 attempts");
}

Filter::Filter(const Format& format, const Options& options) : params_(normalize(format, options)) {
  const auto& f = params_.input;
  const auto& o = params_.options;
  const auto pixels = mul(f.width, f.height);
  const auto length = mul(mul(128, add(f.width, 255) / 128), f.height);
  const auto count = mul(length, o.dynamic_grain ? 3 : 1);
  capacity<Distance>(pixels);
  capacity<double>(count);
  capacity<std::size_t>(f.frames);
  // Prove every plane's final window fits before allocating or advancing RNG.
  for (int plane = 0; plane < (f.color == Color::gray ? 1 : 3); ++plane) {
    const auto w = f.width >> (plane ? f.sub_x : 0), h = f.height >> (plane ? f.sub_y : 0);
    const auto pitch = mul(add(w, 15) / 16, 16);
    if (add(mul(h - 1, pitch), w) > length)
      throw std::overflow_error("grain window exceeds sequence");
  }
  distances_.resize(pixels);
  grain_y_.resize(count);
  if (f.color == Color::yuv)
    grain_c_.resize(count);
  offsets_.resize(f.frames);
  const auto w32 = static_cast<std::uint32_t>(f.width), h32 = static_cast<std::uint32_t>(f.height),
             n32 = static_cast<std::uint32_t>(f.frames);
  Random rng((0x92D68CA2u - static_cast<std::uint32_t>(o.seed)) ^ (w32 << 16) ^ h32 ^ (n32 << 16) ^ n32);
  for (std::size_t y = 0; y < f.height; ++y)
    for (std::size_t x = 0; x < f.width; ++x) {
      rng.next(o.random_algo_grain, o.random_param_grain);
      const auto rx = std::min({static_cast<std::size_t>(o.range), x, f.width - 1 - x}),
                 ry = std::min({static_cast<std::size_t>(o.range), y, f.height - 1 - y});
      const auto radius = o.sample_mode == 1 ? ry : o.sample_mode == 3 ? rx : std::min(rx, ry);
      auto distance = [&] {
        return static_cast<std::uint8_t>(std::abs(static_cast<std::int64_t>(
            away(rng.next(o.random_algo_ref, o.random_param_ref) * static_cast<double>(radius)))));
      };
      auto& d = distances_[y * f.width + x];
      if (radius) {
        d.first = distance();
        if (o.sample_mode == 2)
          d.second = distance();
      }
      if (f.color == Color::yuv && !(x % (std::size_t{1} << f.sub_x)) && !(y % (std::size_t{1} << f.sub_y))) {
        rng.next(o.random_algo_grain, o.random_param_grain);
        rng.next(o.random_algo_grain, o.random_param_grain);
      }
    }
  auto fill = [&](std::vector<double>& table, double amplitude) {
    for (auto& value : table) {
      value = rng.next(o.random_algo_grain, o.random_param_grain) * amplitude;
      value = f.sample == Sample::integer ? away(value) : fp(value);
    }
  };
  fill(grain_y_, params_.grain[0]);
  fill(grain_c_, params_.grain[1]);
  if (o.dynamic_grain)
    for (auto& offset : offsets_) {
      // long double avoids narrowing the potentially wide offset to an int.
      const long double z = static_cast<long double>(rng.next(1, 0)) * length;
      const long double rounded = z > 0 ? std::floor(z + 0.5L) : std::ceil(z - 0.5L);
      const long double aligned_offset = 16 * std::floor((static_cast<long double>(length) + rounded) / 16);
      if (aligned_offset < 0 || aligned_offset > static_cast<long double>(count - length))
        throw std::overflow_error("dynamic grain offset");
      offset = static_cast<std::size_t>(aligned_offset);
      if (offset > count - length)
        throw std::overflow_error("dynamic grain window");
    }
}

Frame Filter::process(std::size_t frame, const std::vector<PlaneView>& input) const {
  return process_frame(*this, frame, input);
}
} // namespace neo_smo::f3kdb
