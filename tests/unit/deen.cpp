#include "hwy/targets.h"
#include "algorithms/deen.hpp"
#include "base/fp16.hpp"
#include "guarded.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using namespace neo_smo;
void check(bool v, const char* why) {
  if (!v)
    throw std::runtime_error(why);
}
std::size_t bytes(DataType t) {
  return t == DataType::U8 ? 1 : t == DataType::F32 ? 4 : 2;
}
double get(const std::uint8_t* p, DataType t) {
  if (t == DataType::U8)
    return *p;
  if (t == DataType::F32) {
    float v;
    std::memcpy(&v, p, 4);
    return v;
  }
  std::uint16_t v;
  std::memcpy(&v, p, 2);
  return t == DataType::F16 ? fp16_to_fp32(v) : v;
}
void put(std::uint8_t* p, DataType t, double x) {
  if (t == DataType::U8) {
    *p = static_cast<std::uint8_t>(x);
    return;
  }
  if (t == DataType::F32) {
    float v = static_cast<float>(x);
    std::memcpy(p, &v, 4);
    return;
  }
  auto v = t == DataType::F16 ? fp32_to_fp16(static_cast<float>(x)) : static_cast<std::uint16_t>(x);
  std::memcpy(p, &v, 2);
}
// Independent direct border mapping and equations, without production tables or padding helpers.
double reference(const DeenOptions& o, const std::array<DeenPlane, 3>& f, int x, int y, bool temporal, bool chroma,
                 double& magnitude, double& low, double& high) {
  const auto& p = f[0];
  const bool integer = p.type == DataType::U8 || p.type == DataType::U16;
  const double peak = integer ? (1u << p.bits) - 1 : 1;
  const double center = get(p.data + y * p.stride + x * bytes(p.type), p.type);
  double numerator = 0, denominator = 0;
  magnitude = 0;
  low = center;
  high = center;
  for (int k = 0; k < (temporal ? 3 : 1); ++k)
    for (int dy = -o.radius; dy <= o.radius; ++dy)
      for (int dx = -o.radius; dx <= o.radius; ++dx) {
        const int xx = std::max(0, std::min(p.width - 1, x + dx)), yy = std::max(0, std::min(p.height - 1, y + dy));
        const double v = get(f[k].data + yy * f[k].stride + xx * bytes(p.type), p.type);
        double g = 1 - (1 - o.minimum) * std::sqrt(double(dx * dx + dy * dy) / (2 * o.radius * o.radius));
        if (std::abs(dx) == o.radius && std::abs(dy) == o.radius)
          g = o.minimum;
        const double t =
            (k == 0 ? (chroma ? o.spatial_uv : o.spatial_y) : (chroma ? o.temporal_uv : o.temporal_y)) * peak / 255;
        // Floating w compares in F32, including differences rounded to the
        // threshold. Keep the independent weight/mean equation in binary64
        // and assess its error for the same selected sample set.
        const bool accepted =
            !integer && o.mode[0] == 'w'
                ? std::abs(static_cast<float>(v) - static_cast<float>(center)) <= static_cast<float>(t)
                : std::abs(v - center) <= t * (o.mode[0] == 'a' ? g : 1);
        if (o.mode[0] == 'a' && !accepted)
          continue;
        const double weight = o.mode[0] == 'w' ? g * (temporal && k == 0 ? 2 : 1) : 1;
        if (weight > 0) {
          const double v0 = accepted ? v : center;
          magnitude = std::max(magnitude, std::abs(v0));
          low = std::min(low, v0);
          high = std::max(high, v0);
        }
        numerator += weight * (accepted ? v : center);
        denominator += weight;
      }
  return numerator / denominator;
}
void run_case(const std::string& mode, DataType type, int bits, int width, int height, int radius, double minimum,
              bool constant, bool chroma, bool temporal, bool at_end, double threshold = -1) {
  const auto b = bytes(type);
  const std::size_t stride = width * b + 3;
  Guarded a(stride * height, at_end), prev((stride + 5) * height, at_end), next((stride + 7) * height, at_end),
      out((stride + 9) * height, at_end);
  std::array<DeenPlane, 3> f{{{a.data, static_cast<std::ptrdiff_t>(stride), width, height, type, bits},
                              {prev.data, static_cast<std::ptrdiff_t>(stride + 5), width, height, type, bits},
                              {next.data, static_cast<std::ptrdiff_t>(stride + 7), width, height, type, bits}}};
  const bool integer = type == DataType::U8 || type == DataType::U16;
  const double peak = integer ? (1u << bits) - 1 : 1;
  for (int k = 0; k < 3; ++k) {
    auto* data = const_cast<std::uint8_t*>(f[k].data);
    std::memset(data, 0xa5, f[k].stride * height);
    for (int y = 0; y < height; ++y)
      for (int x = 0; x < width; ++x) {
        double v = constant ? peak : double((x * 7919 + y * 2797 + k * 431 + 83) % 65536) / 65535 * peak;
        if (integer)
          v = std::floor(v);
        else if (!constant)
          v = v * 4 - 1.5;
        put(data + y * f[k].stride + x * b, type, v);
      }
  }
  DeenOptions o;
  o.mode = mode;
  o.radius = radius;
  o.minimum = minimum;
  o.spatial_y = 91.25;
  o.spatial_uv = 61.75;
  o.temporal_y = 37;
  o.temporal_uv = 127;
  if (threshold >= 0)
    o.spatial_y = o.spatial_uv = o.temporal_y = o.temporal_uv = threshold;
  std::memset(out.data, 0xcc, (stride + 9) * height);
  deen_process(Deen(o), chroma, temporal, f, out.data, stride + 9);
  if (mode == "w2d" && minimum == 1 && integer) {
    std::vector<std::uint8_t> plain((stride + 9) * height);
    auto c = o;
    c.mode = "c2d";
    deen_process(Deen(c), chroma, false, f, plain.data(), stride + 9);
    for (int y = 0; y < height; ++y)
      check(std::memcmp(plain.data() + y * (stride + 9), out.data + y * (stride + 9), width * b) == 0,
            "min=1 w2d differs from c2d");
  }
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      double magnitude = 0, low = 0, high = 0;
      double ref = reference(o, f, x, y, temporal, chroma, magnitude, low, high);
      const double unrounded = ref;
      if (integer)
        ref = std::floor(ref + 0.5);
      const double got = get(out.data + y * (stride + 9) + x * b, type);
      check(std::isfinite(got), "non-finite output");
      check(got >= low && got <= high, "new extrema introduced");
      if (integer && mode[0] == 'w' && got != ref) {
        const double boundary = std::floor(unrounded) + 0.5;
        // Largest-remainder Q15 coefficients have L1 error <=N/65536.
        // Normalization halves its effect over the participating range.
        // Centered Q8 product error is <=(N+1)/512 samples, including the
        // cached bias rounding. U16 uses Q24 coefficients and Q8 products.
        const int side = 2 * radius + 1, taps = side * side * (temporal ? 3 : 1);
        const double budget = type == DataType::U8 ? (taps + 1) / 512.0 + (high - low) * taps / 131072.0
                                                   : (taps + 1) / 512.0 + (high - low) * taps / 33554432.0;
        check(std::abs(unrounded - boundary) <= budget, "integer w difference away from rounding boundary");
      }
      double tolerance = integer ? (mode[0] == 'w' ? 1 : 0) : 2e-6 * std::max(1.0, magnitude);
      if (type == DataType::F16) {
        const double q = std::abs(fp16_to_fp32(fp32_to_fp16(static_cast<float>(ref))));
        tolerance += q < std::ldexp(1.0, -14) ? std::ldexp(1.0, -24) : std::ldexp(1.0, std::ilogb(q) - 10);
      }
      if (std::abs(got - ref) > tolerance) {
        std::cerr << mode << " at " << x << "," << y << " got " << got << " expected " << ref << "\n";
        throw std::runtime_error("reference mismatch");
      }
      if (constant)
        check(got == peak, "constant changed");
    }
    for (std::size_t x = width * b; x < stride + 9; ++x)
      check(out.data[y * (stride + 9) + x] == 0xcc, "padding write");
  }
}
void fixtures(const std::string& mode) {
  std::uint8_t a[9] = {130, 130, 130, 130, 100, 110, 130, 130, 130}, out[9]{};
  DeenOptions o;
  o.mode = mode;
  o.minimum = 1;
  o.spatial_y = 10;
  std::array<DeenPlane, 3> f{{{a, 3, 3, 3, DataType::U8, 8}, {}, {}}};
  deen_process(Deen(o), false, false, f, out, 3);
  check(out[4] == (mode[0] == 'a' ? 105 : 101), "inclusive threshold fixture");
  o.spatial_y = 9.999;
  deen_process(Deen(o), false, false, f, out, 3);
  check(out[4] == 100, "threshold below equality");
  o.spatial_y = 0;
  deen_process(Deen(o), false, false, f, out, 3);
  check(std::memcmp(a, out, 9) == 0, "zero threshold identity");
}
// Check the rounding boundaries of fixed integer means, including large sums.
template <class T = std::uint8_t>
void constant_integer_rounding(const std::string& mode) {
  constexpr int peak = sizeof(T) == 1 ? 255 : 65535;
  const bool temporal = mode[1] == '3';
  const int count = temporal ? 3 : 1;
  for (int radius = 1; radius <= (temporal ? 4 : 7); ++radius) {
    const int side = 2 * radius + 1, divisor = side * side * count;
    std::vector<int> sums;
    for (int q : {0, 1, peak / 2, peak - 1}) {
      sums.push_back(q * divisor + divisor / 2);
      sums.push_back(q * divisor + divisor / 2 + 1);
    }
    sums.push_back(peak * divisor);
    const int width = side * static_cast<int>(sums.size());
    std::array<std::vector<T>, 3> input;
    std::array<DeenPlane, 3> frames{};
    for (int f = 0; f < count; ++f) {
      input[f].resize(width * side);
      frames[f] = {reinterpret_cast<const std::uint8_t*>(input[f].data()),
                   width * static_cast<std::ptrdiff_t>(sizeof(T)),
                   width,
                   side,
                   sizeof(T) == 1 ? DataType::U8 : DataType::U16,
                   int(sizeof(T) * 8)};
    }
    for (std::size_t tile = 0; tile < sums.size(); ++tile) {
      int remaining = sums[tile];
      for (int f = 0; f < count; ++f)
        for (int y = 0; y < side; ++y)
          for (int x = 0; x < side; ++x) {
            const int value = std::min(remaining, peak);
            input[f][y * width + tile * side + x] = static_cast<T>(value);
            remaining -= value;
          }
    }
    DeenOptions o;
    o.mode = mode;
    o.radius = radius;
    o.spatial_y = o.temporal_y = 255;
    std::vector<T> output(width * side);
    deen_process(Deen(o), false, temporal, frames, reinterpret_cast<std::uint8_t*>(output.data()), width * sizeof(T));
    for (std::size_t tile = 0; tile < sums.size(); ++tile)
      check(output[radius * width + tile * side + radius] == (sums[tile] + divisor / 2) / divisor,
            "fixed integer mean rounding boundary");
  }
}
void float_constants(const std::string& mode) {
  for (auto type : {DataType::F16, DataType::F32})
    for (double value : {0.0, -0.125, 0.3, std::ldexp(1.0, -24), -std::ldexp(1.0, -24), 65504.0}) {
      std::uint8_t data[4]{}, out[4]{};
      put(data, type, value);
      DeenOptions o;
      o.mode = mode;
      o.radius = mode[1] == '3' ? 4 : 7;
      o.minimum = 0.25;
      DeenPlane p{data, static_cast<std::ptrdiff_t>(bytes(type)), 1, 1, type, type == DataType::F16 ? 16 : 32};
      deen_process(Deen(o), false, mode[1] == '3', {p, p, p}, out, bytes(type));
      check(get(out, type) == get(data, type), "arbitrary floating constant changed");
    }
}
void weighted_fixture(const std::string& mode) {
  std::uint8_t data[9] = {255, 110, 0, 110, 100, 110, 0, 110, 255}, out[9]{};
  DeenOptions o;
  o.mode = mode;
  o.minimum = 0;
  o.spatial_y = 255;
  DeenPlane p{data, 3, 3, 3, DataType::U8, 8};
  deen_process(Deen(o), false, false, {p, {}, {}}, out, 3);
  check(out[4] == 105, "distance weighted fixture");
}
void weighted_byte_constants(const std::string& mode) {
  constexpr int width = 129;
  const bool temporal = mode[1] == '3';
  for (int radius : {1, temporal ? 4 : 7})
    for (double minimum : {0.0, 0.37, 1.0})
      for (int value : {0, 1, 63, 127, 128, 129, 254, 255}) {
        std::array<std::uint8_t, width> input, output;
        input.fill(static_cast<std::uint8_t>(value));
        DeenPlane plane{input.data(), width, width, 1, DataType::U8, 8};
        DeenOptions o;
        o.mode = mode;
        o.radius = radius;
        o.minimum = minimum;
        o.spatial_y = o.temporal_y = 255;
        deen_process(Deen(o), false, temporal, {plane, plane, plane}, output.data(), width);
        check(input == output, "weighted byte constant changed");
      }
}
void weighted_word_constants(const std::string& mode) {
  constexpr int width = 67;
  for (int bits : {10, 12, 16})
    for (int radius : {1, mode[1] == '3' ? 4 : 7})
      for (double minimum : {0.0, 0.000001, 0.37, 1.0})
        for (unsigned value :
             {0u, 1u, 127u, 255u, 256u, (1u << (bits - 1)) - 1, 1u << (bits - 1), (1u << bits) - 2, (1u << bits) - 1}) {
          std::vector<std::uint16_t> input(width, static_cast<std::uint16_t>(value)), output(width);
          DeenPlane plane{
              reinterpret_cast<const std::uint8_t*>(input.data()), width * 2, width, 1, DataType::U16, bits};
          DeenOptions o;
          o.mode = mode;
          o.radius = radius;
          o.minimum = minimum;
          deen_process(Deen(o), false, mode[1] == '3', {plane, plane, plane},
                       reinterpret_cast<std::uint8_t*>(output.data()), width * 2);
          check(input == output, "weighted word constant changed");
        }
}
void weighted_float_extremes(const std::string& mode) {
  const bool temporal = mode[1] == '3';
  constexpr int width = 65, height = 2;
  for (auto type : {DataType::F16, DataType::F32}) {
    const auto b = bytes(type);
    const int stride = width * b + 3;
    std::array<std::vector<std::uint8_t>, 3> input;
    std::array<DeenPlane, 3> frames{};
    for (int f = 0; f < 3; ++f) {
      input[f].resize((stride + f) * height);
      frames[f] = {input[f].data(), stride + f, width, height, type, type == DataType::F16 ? 16 : 32};
    }
    std::vector<std::uint8_t> output(stride * height);
    DeenOptions o;
    o.mode = mode;
    o.radius = temporal ? 4 : 7;
    o.spatial_y = o.temporal_y = 255;
    o.minimum = 0;
    const double maximum = type == DataType::F16 ? 65504.0 : std::numeric_limits<float>::max();
    const double tiny = type == DataType::F16 ? std::ldexp(1.0, -24) : std::numeric_limits<float>::denorm_min();
    const std::array<double, 7> values{maximum, -maximum, 0, 1, -1, tiny, -tiny};
    for (int f = 0; f < 3; ++f)
      for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
          put(input[f].data() + y * frames[f].stride + x * b, type, values[(x + y + f) % values.size()]);
    deen_process(Deen(o), false, temporal, frames, output.data(), stride);
    for (int y = 0; y < height; ++y)
      for (int x = 0; x < width; ++x) {
        double magnitude = 0, low = 0, high = 0;
        const double ref = reference(o, frames, x, y, temporal, false, magnitude, low, high);
        const double got = get(output.data() + y * stride + x * b, type);
        const double budget = 2e-6 * std::max(1.0, magnitude) + (type == DataType::F16 ? 0x1p-10 : 0);
        check(std::isfinite(got) && got >= low && got <= high, "weighted float extreme range");
        if (std::abs(got - ref) > budget) {
          std::cerr << "weighted float extreme " << static_cast<int>(type) << " at " << x << "," << y << " got " << got
                    << " expected " << ref << " budget " << budget << "\n";
          throw std::runtime_error("weighted float extreme accuracy");
        }
      }
    for (double value : values) {
      for (int f = 0; f < 3; ++f)
        for (int y = 0; y < height; ++y)
          for (int x = 0; x < width; ++x)
            put(input[f].data() + y * frames[f].stride + x * b, type, value);
      deen_process(Deen(o), false, temporal, frames, output.data(), stride);
      for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
          check(get(output.data() + y * stride + x * b, type) == value, "weighted float extreme constant");
    }
    // Test invalid values in every participating frame and in the vector tail.
    for (int f = 0; f < (temporal ? 3 : 1); ++f)
      for (std::uint32_t bad : {0x7f800000u, 0xff800000u, 0x7fc00001u, 0x7f800001u}) {
        auto* ptr = input[f].data() + frames[f].stride + (width - 1) * b;
        const auto saved = get(ptr, type);
        if (type == DataType::F32)
          std::memcpy(ptr, &bad, 4);
        else {
          const std::uint16_t half = bad == 0x7f800000u ? 0x7c00 : bad == 0xff800000u ? 0xfc00 : 0x7c01;
          std::memcpy(ptr, &half, 2);
        }
        bool rejected = false;
        try {
          deen_process(Deen(o), false, temporal, frames, output.data(), stride);
        } catch (const std::invalid_argument&) {
          rejected = true;
        }
        check(rejected, "weighted non-finite input accepted");
        put(ptr, type, saved);
      }
  }
}
void adaptive_fixture(const std::string& mode) {
  std::uint8_t data[9] = {112, 112, 112, 112, 100, 112, 112, 112, 112}, out[9]{};
  DeenOptions o;
  o.mode = mode;
  o.minimum = 0.5;
  o.spatial_y = 20;
  DeenPlane p{data, 3, 3, 3, DataType::U8, 8};
  deen_process(Deen(o), false, false, {p, {}, {}}, out, 3);
  check(out[4] == 110, "adaptive corner vs axial threshold");
  data[0] = 110;
  deen_process(Deen(o), false, false, {p, {}, {}}, out, 3);
  check(out[4] == 110, "adaptive corner equality");
  std::fill_n(data, 9, 130);
  data[4] = 100;
  data[5] = 106;
  o.minimum = 0.010050506338833419;
  deen_process(Deen(o), false, false, {p, {}, {}}, out, 3);
  // Distance evaluation may round this boundary to either side of six.
  check(out[4] == 100 || out[4] == 103, "adaptive near-boundary selection");
  std::fill_n(data, 9, 112);
  data[4] = 100;
  // Four equal corners count even with zero corner threshold; axial samples remain eligible.
  for (int i : {0, 2, 6, 8})
    data[i] = 100;
  o.minimum = 0;
  o.spatial_y = 255;
  deen_process(Deen(o), false, false, {p, {}, {}}, out, 3);
  check(out[4] == 105, "adaptive zero corner threshold still counts equal samples");
}
void rational_adaptive_boundary(const std::string& mode) {
  float data[9] = {1, 1, 1, 1, 0, 1, 1, 1, 0.0625f}, out[9]{};
  DeenOptions o;
  o.mode = mode;
  o.spatial_y = 17;
  o.minimum = 0.9375;
  DeenPlane p{reinterpret_cast<const std::uint8_t*>(data), 12, 3, 3, DataType::F32, 32};
  deen_process(Deen(o), false, false, {p, {}, {}}, reinterpret_cast<std::uint8_t*>(out), 12);
  check(out[4] == 0.03125f, "rational corner threshold equality");
  o.spatial_y = 16.999;
  deen_process(Deen(o), false, false, {p, {}, {}}, reinterpret_cast<std::uint8_t*>(out), 12);
  check(out[4] == 0, "rational corner threshold just below equality");
}
void scaled_boundary(const std::string& mode) {
  std::uint16_t data[9] = {200, 200, 200, 200, 100, 117, 200, 200, 200}, out[9]{};
  DeenOptions o;
  o.mode = mode;
  o.minimum = 1;
  o.spatial_y = 4.237536656891495;
  DeenPlane p{reinterpret_cast<const std::uint8_t*>(data), 6, 3, 3, DataType::U16, 10};
  deen_process(Deen(o), false, false, {p, {}, {}}, reinterpret_cast<std::uint8_t*>(out), 6);
  check(out[4] == (mode[0] == 'a' ? 109 : 102), "scaled threshold accepts rounded equality");
  if (mode[1] == '3') {
    std::uint16_t center = 100, neighbor = 117, result = 0;
    DeenPlane cp{reinterpret_cast<const std::uint8_t*>(&center), 2, 1, 1, DataType::U16, 10};
    DeenPlane np{reinterpret_cast<const std::uint8_t*>(&neighbor), 2, 1, 1, DataType::U16, 10};
    auto temporal = o;
    temporal.temporal_y = o.spatial_y;
    deen_process(Deen(temporal), false, true, {cp, np, np}, reinterpret_cast<std::uint8_t*>(&result), 2);
    check(result == (mode[0] == 'w' ? 109 : 111), "temporal scaled threshold accepts rounded equality");
  }
  o.spatial_y = 4.237;
  deen_process(Deen(o), false, false, {p, {}, {}}, reinterpret_cast<std::uint8_t*>(out), 6);
  check(out[4] == 100, "scaled threshold below difference");
}
void temporal_fixtures(const std::string& mode) {
  std::uint8_t a = 100, b = 80, c = 140, out = 0;
  std::array<DeenPlane, 3> f{
      {{&a, 1, 1, 1, DataType::U8, 8}, {&b, 1, 1, 1, DataType::U8, 8}, {&c, 1, 1, 1, DataType::U8, 8}}};
  DeenOptions o;
  o.mode = mode;
  o.minimum = 1;
  o.spatial_y = o.temporal_y = 255;
  deen_process(Deen(o), false, true, f, &out, 1);
  check(out == (mode[0] == 'w' ? 105 : 107), "three-frame mean fixture");
  o.temporal_y = 20;
  deen_process(Deen(o), false, true, f, &out, 1);
  check(out == (mode[0] == 'a' ? 90 : mode[0] == 'w' ? 95 : 93), "temporal rejection fixture");
  o.scene_threshold = 20;
  check(!deen_scene_cut(Deen(o), {f[0]}, {f[1]}), "scene equality");
  o.scene_threshold = 19.999;
  check(deen_scene_cut(Deen(o), {f[0]}, {f[1]}), "scene above threshold");
  std::uint8_t zeros[10]{}, ones[10]{1}, twos[10]{2};
  DeenPlane zp{zeros, 10, 10, 1, DataType::U8, 8}, op{ones, 10, 10, 1, DataType::U8, 8},
      tp{twos, 10, 10, 1, DataType::U8, 8};
  o.scene_threshold = 0.101;
  check(!deen_scene_cut(Deen(o), {zp, zp, zp}, {op, tp, zp}), "scene mean below threshold");
  o.scene_threshold = 0.099;
  check(deen_scene_cut(Deen(o), {zp, zp, zp}, {op, tp, zp}), "scene mean above threshold");
  std::uint8_t zero5[5]{}, two5[5]{2}, one5[5]{1}, twelve5[5]{12};
  DeenPlane z5{zero5, 5, 5, 1, DataType::U8, 8}, t5{two5, 5, 5, 1, DataType::U8, 8}, o5{one5, 5, 5, 1, DataType::U8, 8},
      d5{twelve5, 5, 5, 1, DataType::U8, 8};
  o.scene_threshold = 1;
  check(!deen_scene_cut(Deen(o), {z5, z5, z5}, {t5, o5, d5}), "scene equality after plane means");
  o.scene_threshold = 0.999;
  check(deen_scene_cut(Deen(o), {z5, z5, z5}, {t5, o5, d5}), "scene mean above threshold");
  o.scenechange = false;
  check(!deen_scene_cut(Deen(o), {}, {}), "disabled scene reads");
  o.scenechange = true;
  o.scene_threshold = 4;
  std::uint8_t zero = 0, twelve = 12;
  DeenPlane z{&zero, 1, 1, 1, DataType::U8, 8}, t{&twelve, 1, 1, 1, DataType::U8, 8};
  check(!deen_scene_cut(Deen(o), {z, z, z}, {z, t, z}), "plane mean equality");
  o.scene_threshold = 3.999;
  check(deen_scene_cut(Deen(o), {z, z, z}, {z, t, z}), "equal plane weighting");
}
void scene_integer_vectors() {
  for (int bits : {8, 10, 16}) {
    const int bytes = bits == 8 ? 1 : 2;
    for (int width : {1, 31, 32, 33, 65, 4101}) {
      const int stride_a = width * bytes + 3, stride_b = width * bytes + 7;
      std::vector<std::uint8_t> a(stride_a * 2 + 1), b(stride_b * 2 + 1);
      const auto type = bits == 8 ? DataType::U8 : DataType::U16;
      for (int y = 0; y < 2; ++y)
        for (int x = 0; x < width; ++x)
          put(b.data() + 1 + y * stride_b + x * bytes, type, (1u << bits) - 1);
      DeenPlane av{a.data() + 1, stride_a, width, 2, type, bits};
      DeenPlane bv{b.data() + 1, stride_b, width, 2, type, bits};
      DeenOptions o;
      o.scene_threshold = 255;
      check(!deen_scene_cut(Deen(o), {av}, {bv}), "maximum scene equality");
      o.scene_threshold = 254;
      check(deen_scene_cut(Deen(o), {av}, {bv}), "maximum scene difference");
      if (bits == 10) {
        put(b.data() + 1, type, 1024);
        bool rejected = false;
        try {
          (void)deen_scene_cut(Deen(o), {av}, {bv});
        } catch (const std::invalid_argument&) {
          rejected = true;
        }
        check(rejected, "scene accepted sample beyond bit depth");
      }
    }
  }
}
void invalid() {
  DeenOptions o;
  for (double v : {-1.0, 256.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    o.spatial_y = v;
    bool caught = false;
    try {
      Deen d(o);
    } catch (const std::exception&) {
      caught = true;
    }
    check(caught, "invalid threshold accepted");
  }
  float nan = std::numeric_limits<float>::quiet_NaN(), output = 0;
  o = DeenOptions{};
  o.mode = "c2d";
  std::array<DeenPlane, 3> f{{{reinterpret_cast<std::uint8_t*>(&nan), 4, 1, 1, DataType::F32, 32}, {}, {}}};
  bool caught = false;
  try {
    deen_process(Deen(o), false, false, f, reinterpret_cast<std::uint8_t*>(&output), 4);
  } catch (const std::exception&) {
    caught = true;
  }
  check(caught, "non-finite sample accepted");
}
} // namespace
int main(int argc, char** argv) {
  try {
    for (const auto target : hwy::SupportedAndGeneratedTargets()) {
      hwy::SetSupportedTargetsForTest(target);
      std::printf("target %s\n", hwy::TargetName(target));

      const std::string mode = argc > 1 ? argv[1] : "c2d";
      fixtures(mode);
      if (mode[0] == 'c') {
        constant_integer_rounding(mode);
        constant_integer_rounding<std::uint16_t>(mode);
      }
      // Exercise both halves of the byte batches and their boundary tails.
      for (int width : {63, 64, 127, 128, 129})
        for (double minimum : {0.000001, 0.37})
          for (bool at_end : {false, true})
            run_case(mode, DataType::U8, 8, width, 3, mode[1] == '3' ? 4 : 7, minimum, false, false, mode[1] == '3',
                     at_end);
      if (mode[0] == 'a' || mode[0] == 'c')
        for (auto type : {DataType::U8, DataType::U16})
          for (int width : {15, 16, 17, 63, 64, 65, 127, 128, 129})
            for (double threshold : {0.0, 7.0, 255.0})
              for (bool at_end : {false, true})
                run_case(mode, type, type == DataType::U8 ? 8 : 16, width, 3, mode[1] == '3' ? 4 : 7, 0.37, false,
                         false, mode[1] == '3', at_end, threshold);
      invalid();
      scene_integer_vectors();
      float_constants(mode);
      scaled_boundary(mode);
      if (mode[0] == 'w') {
        weighted_fixture(mode);
        weighted_byte_constants(mode);
        weighted_word_constants(mode);
        for (int width : {63, 64, 65, 127, 128, 129})
          for (double minimum : {0.000001, 0.37})
            for (double threshold : {0.0, 7.0, 255.0})
              run_case(mode, DataType::U16, 16, width, 3, mode[1] == '3' ? 4 : 7, minimum, false, false, mode[1] == '3',
                       true, threshold);
        weighted_float_extremes(mode);
      }
      if (mode[0] == 'a') {
        adaptive_fixture(mode);
        rational_adaptive_boundary(mode);
      }
      if (mode[1] == '3')
        temporal_fixtures(mode);
      for (auto type : {DataType::U8, DataType::U16, DataType::F16, DataType::F32})
        for (int bits : {8, 10, 12, 16, 32}) {
          if ((type == DataType::U8 && bits != 8) || (type == DataType::U16 && (bits == 8 || bits == 32)) ||
              (type == DataType::F16 && bits != 16) || (type == DataType::F32 && bits != 32))
            continue;
          for (int r = 1; r <= (mode[1] == '3' ? 4 : 7); ++r)
            for (int w : {1, 7, 8, 9, 17, 31, 32, 33, 65})
              for (double m : {0.0, 0.5, 1.0}) {
                run_case(mode, type, bits, w, 3, r, m, false, true, mode[1] == '3', true);
                run_case(mode, type, bits, w, 1, r, m, true, false, false, false);
              }
        }
      std::cout << mode << " passed\n";
    }
    hwy::SetSupportedTargetsForTest(0);
  } catch (const std::exception& e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
