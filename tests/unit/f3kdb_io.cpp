#include "algorithms/f3kdb_views.hpp"
#include "kernels/f3kdb.hpp"
#include "algorithms/f3kdb.hpp"
#ifndef HWY_COMPILE_ALL_ATTAINABLE
#define HWY_COMPILE_ALL_ATTAINABLE
#endif
#include "hwy/targets.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

namespace ref = neo_smo::f3kdb;
namespace simd = neo_smo::f3kdb;
namespace {
std::size_t checks = 0, converted = 0, gradients = 0, gradient_differences = 0;
double gradient_max_error = 0;
void check(bool ok, const std::string& message) {
  ++checks;
  if (!ok)
    throw std::runtime_error(message);
}
template <class T>
void store(void* p, T value) {
  std::memcpy(p, &value, sizeof(value));
}
template <class T>
T load(const void* p) {
  T value;
  std::memcpy(&value, p, sizeof(value));
  return value;
}
void check_plane_views() {
  // Minimal final row: there is no padding after its last active byte.
  std::array<std::uint8_t, 13> storage{};
  auto rows = simd::detail::plane_view(storage.data() + 1, 5, 2, 7).cursor();
  check(rows.ptr() == storage.data() + 1, "byte view origin");
  ++rows;
  rows.row()[4] = 42;
  check(storage[12] == 42, "odd byte stride and minimal final row");
  std::array<double, 8> values{};
  const auto plane = simd::detail::plane_view(values.data(), 3, 2, 4 * sizeof(double));
  auto current = plane.cursor();
  ++current;
  current.row()[2] = 17;
  check(values[6] == 17 && current[-1][2] == 0, "typed element stride and backward row");
  const auto fails = [&](std::size_t w, std::size_t h, std::size_t stride) {
    bool rejected = false;
    try {
      (void)simd::detail::plane_view(values.data(), w, h, stride);
    } catch (const std::exception&) {
      rejected = true;
    }
    check(rejected, "invalid view dimensions/stride rejected before narrowing");
  };
  const auto too_large = static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) + 1;
  fails(too_large, 1, 8);
  fails(1, too_large, 8);
  fails(1, 1, too_large * sizeof(double));
  fails(1, 1, 9);
  fails(2, 1, 8);
}
bool identical(double a, double b) {
  return a == b && (a != 0 || std::signbit(a) == std::signbit(b));
}

void check_unpack(const std::vector<std::uint32_t>& bits, int depth, bool integer, std::size_t width) {
  const std::size_t bytes = depth == 8 ? 1 : depth == 16 || integer ? 2 : 4;
  const auto height = (bits.size() + width - 1) / width, stride = width * bytes + 3;
  // Exact end of the last row: sanitizers can detect a tail over-read.
  std::vector<std::uint8_t> input(1 + (height - 1) * stride + width * bytes, 0xa5);
  std::vector<double> expected(width * height), actual(width * height + 2, 12345.5);
  for (std::size_t i = 0; i < expected.size(); ++i) {
    const auto raw = i < bits.size() ? bits[i] : 0;
    auto* p = input.data() + 1 + (i / width) * stride + (i % width) * bytes;
    if (bytes == 1)
      *p = static_cast<std::uint8_t>(raw);
    else if (bytes == 2)
      store(p, static_cast<std::uint16_t>(raw));
    else
      store(p, raw);
    expected[i] = integer       ? std::ldexp(double(raw), 16 - depth)
                  : depth == 16 ? double(ref::decode_half(static_cast<std::uint16_t>(raw)))
                                : double(load<float>(&raw));
  }
  check(simd::unpack_plane({input.data() + 1, actual.data() + 1, width, height, stride, depth, integer}),
        "valid input rejected");
  for (std::size_t i = 0; i < expected.size(); ++i) {
    ++converted;
    check(identical(actual[i + 1], expected[i]), "unpack depth=" + std::to_string(depth) + " i=" + std::to_string(i));
  }
  check(actual.front() == 12345.5 && actual.back() == 12345.5, "unpack output guards");
  if (integer) {
    std::vector<std::int32_t> packed(width * height + 2, -123456);
    check(simd::unpack_integer_plane({input.data() + 1, packed.data() + 1, width, height, stride, depth}), "valid packed input rejected");
    for (std::size_t i = 0; i < expected.size(); ++i) check(packed[i + 1] == expected[i], "packed unpack value");
    check(packed.front() == -123456 && packed.back() == -123456, "packed unpack guards");
  }
}

void unpack_tests() {
  std::vector<std::uint32_t> bits;
  for (int depth = 8; depth <= 16; ++depth) {
    const auto maximum = (std::uint32_t{1} << depth) - 1;
    bits = {0, 1, 2, maximum / 2, maximum / 2 + 1, maximum - 1, maximum};
    for (unsigned code : {16u, 128u, 235u, 240u}) {
      const auto edge = code << (depth - 8);
      bits.insert(bits.end(), {edge - 1, edge, edge + 1});
    }
    check_unpack(bits, depth, true, 17);
  }
  bits.clear();
  for (unsigned sign : {0u, 0x8000u})
    for (unsigned exponent : {0u, 1u, 2u, 14u, 15u, 16u, 29u, 30u})
      for (unsigned mantissa : {0u, 1u, 2u, 511u, 512u, 1022u, 1023u})
        bits.push_back(sign | (exponent << 10) | mantissa);
  check_unpack(bits, 16, false, 17);
  bits.clear();
  std::mt19937 rng(0x494f3332);
  for (unsigned sign : {0u, 0x80000000u})
    for (unsigned exponent : {0u, 1u, 2u, 126u, 127u, 128u, 253u, 254u})
      for (unsigned mantissa : {0u, 1u, 0x7fffffu, static_cast<unsigned>(rng() & 0x7fffff)})
        bits.push_back(sign | (exponent << 23) | mantissa);
  check_unpack(bits, 32, false, 17);
  for (std::size_t width : {1, 2, 3, 7, 15, 17, 31, 65}) {
    check_unpack({0, 0x8000, 1, 0x8001, 0x7bff, 0xfbff, 0x3c00}, 16, false, width);
    check_unpack({0, 0x80000000, 1, 0x80000001, 0x7f7fffff, 0xff7fffff, 0x3f800000}, 32, false, width);
    for (int depth : {8, 10, 16}) check_unpack({0, 1, 17, (1u << depth) - 1}, depth, true, width);
    for (int depth = 9; depth < 16; ++depth) {
      std::vector<std::uint8_t> raw(width * 2 + 1, 0);
      std::vector<double> out(width);
      store(raw.data() + 1 + (width - 1) * 2, static_cast<std::uint16_t>(1u << depth));
      check(!simd::unpack_plane({raw.data() + 1, out.data(), width, 1, width * 2, depth, true}),
            "out-of-depth integer accepted");
      std::vector<std::int32_t> packed(width);
      check(!simd::unpack_integer_plane({raw.data() + 1, packed.data(), width, 1, width * 2, depth}), "out-of-depth packed integer accepted");
    }
    for (std::uint32_t bad : {0x7f800000u, 0xff800000u, 0x7fc00000u, 0x7f800001u, 0xffc12345u}) {
      std::vector<std::uint8_t> raw(width * 4 + 1, 0);
      std::vector<double> out(width);
      for (std::size_t at = 0; at < width; ++at) {
        store(raw.data() + 1 + 4 * at, bad);
        check(!simd::unpack_plane({raw.data() + 1, out.data(), width, 1, width * 4, 32, false}),
              "nonfinite F32 accepted");
        store(raw.data() + 1 + 4 * at, std::uint32_t{0});
      }
    }
  }
  for (unsigned raw : {0x7c00u, 0xfc00u, 0x7c01u, 0x7e00u, 0x7fffu, 0xfc01u, 0xfe00u, 0xffffu}) {
    std::uint8_t bytes[3];
    store(bytes + 1, static_cast<std::uint16_t>(raw));
    double out;
    check(!simd::unpack_plane({bytes + 1, &out, 1, 1, 2, 16, false}), "nonfinite half accepted");
  }
}

void check_pack(const std::vector<double>& values, int depth, bool integer_input, bool integer_output, bool chroma,
                bool limited, int dither, std::size_t y, std::size_t origin = 0) {
  const std::size_t bytes = depth == 8 ? 1 : depth <= 16 ? 2 : 4;
  std::vector<std::uint8_t> expected(values.size() * bytes + 9, 0xa5), actual = expected;
  const double offset = chroma ? 0.5 : 0, scale = std::ldexp(1.0, depth), tvscale = std::ldexp(1.0, depth - 8);
  const double lo = limited ? 16 * tvscale : 0, hi = limited ? (chroma ? 240 : 235) * tvscale : scale - 1;
  for (std::size_t x = 0; x < values.size(); ++x) {
    auto* p = expected.data() + 1 + x * bytes;
    double v = values[x];
    if (!integer_output) {
      if (integer_input)
        v = v / 65536 - offset;
      if (limited)
        v = std::clamp(v, 16.0 / 256 - offset, (chroma ? 240.0 : 235.0) / 256 - offset);
      if (depth == 16)
        store(p, ref::encode_half(v));
      else
        store(p, static_cast<float>(std::clamp(v, -double(std::numeric_limits<float>::max()),
                                               double(std::numeric_limits<float>::max()))));
    } else {
      v = integer_input ? std::ldexp(v, depth - 16) : (v + offset) * scale;
      v = std::clamp(v, lo - 1, hi + 1);
      v = dither == 1 ? ref::round_even(v) : std::floor(v + (ref::bayer16(y, x + origin) + 0.5) / 256);
      v = std::clamp(v, lo, hi);
      if (bytes == 1)
        *p = static_cast<std::uint8_t>(v);
      else
        store(p, static_cast<std::uint16_t>(v));
    }
  }
  simd::pack_row({values.data(), actual.data() + 1, values.size(), y, depth, dither, integer_input, integer_output,
                  chroma, limited, origin});
  converted += values.size();
  if (actual != expected) {
    std::size_t i = 0;
    while (i < actual.size() && actual[i] == expected[i])
      ++i;
    throw std::runtime_error("pack depth=" + std::to_string(depth) + " integer_in=" + std::to_string(integer_input) +
                             " integer_out=" + std::to_string(integer_output) + " chroma=" + std::to_string(chroma) +
                             " limited=" + std::to_string(limited) + " dither=" + std::to_string(dither) +
                             " y=" + std::to_string(y) + " byte=" + std::to_string(i));
  }
  check(true, "pack bytes including guards");
}

void pack_tests() {
  for (std::size_t origin = 1; origin < 16; ++origin)
    check_pack({-17.0, 0.0, 127.0, 128.0, 255.0, 256.0, 32767.0, 65535.0}, 8, true, true, false, false, 2, 7, origin);
  // Representative exponent/mantissa transitions, even/odd ties and both
  // neighbouring F64 values catch F64->F32->F16 double rounding without enumeration.
  std::vector<double> halves = {0.0, -0.0, 65504, -65504, 1e300, -1e300};
  for (unsigned exponent : {0u, 1u, 2u, 14u, 15u, 16u, 29u, 30u})
    for (unsigned mantissa : {0u, 1u, 511u, 512u, 1022u, 1023u}) {
      const unsigned h = (exponent << 10) | mantissa;
      if (h == 0x7bff)
        continue; // Saturation is covered separately.
      const double midpoint = (double(ref::decode_half(static_cast<std::uint16_t>(h))) +
                               double(ref::decode_half(static_cast<std::uint16_t>(h + 1)))) *
                              0.5;
      for (double v : {std::nextafter(midpoint, 0.0), midpoint,
                       std::nextafter(midpoint, std::numeric_limits<double>::infinity())}) {
        halves.push_back(v);
        halves.push_back(-v);
      }
    }
  check_pack(halves, 16, false, false, false, false, 1, 0);
  std::vector<double> singles = {
      0.0, -0.0, 1e300, -1e300, double(std::numeric_limits<float>::max()), -double(std::numeric_limits<float>::max())};
  for (unsigned exponent : {0u, 1u, 2u, 126u, 127u, 128u, 253u})
    for (unsigned mantissa : {0u, 1u, 0x7ffffeu}) {
      const std::uint32_t bits = (exponent << 23) | mantissa, next = bits + 1;
      const double midpoint = (double(load<float>(&bits)) + double(load<float>(&next))) * 0.5;
      for (double v : {std::nextafter(midpoint, 0.0), midpoint,
                       std::nextafter(midpoint, std::numeric_limits<double>::infinity())}) {
        singles.push_back(v);
        singles.push_back(-v);
      }
    }
  check_pack(singles, 32, false, false, false, false, 1, 0);
  std::mt19937 rng(0x42415945);
  for (bool integer_output : {false, true})
    for (int depth = 8; depth <= 32; ++depth) {
      if (integer_output ? depth > 16 : (depth != 16 && depth != 32))
        continue;
      for (bool integer_input : {false, true})
        for (bool chroma : {false, true})
          for (bool limited : {false, true})
            for (int dither : {1, 2}) { // Full Bayer phase coverage lives in the focused boundary test below.
              const std::size_t y = (depth + integer_input + 2 * chroma + 4 * limited) % 16;
              std::vector<double> values(17);
              for (std::size_t x = 0; x < values.size(); ++x) {
                double code = double(rng() % 65536) + ((x % 2) ? 0.5 : double(rng() % 256) / 256);
                if (x % 11 == 0)
                  code = -1e30;
                if (x % 11 == 1)
                  code = 1e30;
                if (x % 11 == 2)
                  code = 0;
                if (x % 11 == 3)
                  code = 65535;
                values[x] = integer_input ? code : code / 65536 - (chroma ? 0.5 : 0);
              }
              check_pack(values, depth, integer_input, integer_output, chroma, limited, dither, y);
            }
    }
  for (std::size_t width : {1, 2, 3, 7, 15, 17, 31, 65}) {
    std::vector<double> values(width);
    for (std::size_t i = 0; i < width; ++i)
      values[i] = static_cast<int>(i % 7) - 3.5;
    for (int depth : {8, 10, 16, 32})
      for (int dither : {1, 2})
        check_pack(values, depth, false, depth != 32, false, false, dither, std::numeric_limits<std::size_t>::max());
    check_pack(values, 16, false, false, false, false, 1, 0);
  }
  // Exact quantization boundaries and both neighbouring F64 values, including
  // every Bayer phase at depth 8, with representative phases at other depths.
  // Compare after the same domain/chroma conversion as N07.
  for (int depth = 8; depth <= 16; ++depth)
    for (bool integer_input : {false, true})
      for (bool chroma : {false, true})
        for (int dither : {1, 2})
          for (std::size_t y = 0; y < (dither == 2 && depth == 8 ? 16u : 1u); ++y)
            for (int side : {-1, 0, 1}) {
              std::vector<double> values(17);
              for (std::size_t x = 0; x < values.size(); ++x) {
                double code = 100 + double(x % 3) + (dither == 1 ? 0.5 : -(ref::bayer16(y, x) + 0.5) / 256);
                if (side)
                  code = std::nextafter(code, side < 0 ? 0.0 : 1000.0);
                values[x] =
                    integer_input ? std::ldexp(code, 16 - depth) : std::ldexp(code, -depth) - (chroma ? 0.5 : 0);
              }
              check_pack(values, depth, integer_input, true, chroma, false, dither, y);
            }
}

void sobel_tests() {
  std::mt19937 rng(0x534f4245);
  // Pair dimensions: narrow axes, both sides of the clamp distance, and one
  // image with an unclamped interior in both axes. Avoid the Cartesian product.
  for (const auto shape : {std::pair<std::size_t, std::size_t>{1, 1},
                           {2, 3},
                           {3, 2},
                           {19, 21},
                           {20, 20},
                           {21, 19},
                           {39, 2},
                           {2, 40},
                           {41, 41},
                           {65, 1}})
    for (int pattern = 0; pattern < 6; ++pattern) {
      const auto width = shape.first, height = shape.second;
      const auto count = width * height;
      std::vector<double> source(count + 1), gx(count + 2, 12345.5), gy(count + 2, 12345.5);
      auto* data = source.data() + 1;
      for (std::size_t i = 0; i < count; ++i) {
        if (pattern == 0)
          data[i] = rng() % 65536;
        else if (pattern == 1) {
          const std::uint32_t bits = rng() & 0xff7fffffu;
          data[i] = load<float>(&bits);
        } else if (pattern == 2)
          data[i] = (i % 2 ? -1.0 : 1.0) * std::numeric_limits<float>::max();
        else if (pattern == 3)
          data[i] = (i % 3 - 1.0) * std::numeric_limits<float>::denorm_min();
        else
          data[i] = (pattern == 4 ? i % width >= width / 2 : i / width >= height / 2) ? 65535 : 0;
      }
      simd::sobel_plane(data, gx.data() + 1, gy.data() + 1, width, height);
      if (pattern == 0 || pattern >= 4) {
        std::vector<std::int32_t> compact(data, data + count), cx(count + 2, -999999), cy(count + 2, -999999);
        simd::sobel_plane(compact.data(), cx.data() + 1, cy.data() + 1, width, height);
        for (std::size_t i = 0; i < count; ++i)
          check(cx[i + 1] == gx[i + 1] && cy[i + 1] == gy[i + 1], "compact integer Sobel exact");
        check(cx.front() == -999999 && cx.back() == -999999 && cy.front() == -999999 && cy.back() == -999999,
              "compact integer Sobel guards");
      }
      for (std::size_t y = 0; y < height; ++y)
        for (std::size_t x = 0; x < width; ++x) {
          auto at = [&](int dx, int dy) {
            const auto xx = std::clamp<std::int64_t>(static_cast<std::int64_t>(x) + dx, 0, width - 1);
            const auto yy = std::clamp<std::int64_t>(static_cast<std::int64_t>(y) + dy, 0, height - 1);
            return data[std::size_t(yy) * width + std::size_t(xx)];
          };
          const double expected[2] = {
              at(20, -20) + 2 * at(20, 0) + at(20, 20) - at(-20, -20) - 2 * at(-20, 0) - at(-20, 20),
              at(-20, -20) + 2 * at(0, -20) + at(20, -20) - at(-20, 20) - 2 * at(0, 20) - at(20, 20)};
          const double magnitude[2] = {std::abs(at(20, -20)) + 2 * std::abs(at(20, 0)) + std::abs(at(20, 20)) +
                                           std::abs(at(-20, -20)) + 2 * std::abs(at(-20, 0)) + std::abs(at(-20, 20)),
                                       std::abs(at(-20, -20)) + 2 * std::abs(at(0, -20)) + std::abs(at(20, -20)) +
                                           std::abs(at(-20, 20)) + 2 * std::abs(at(0, 20)) + std::abs(at(20, 20))};
          const double actual[2] = {gx[y * width + x + 1], gy[y * width + x + 1]};
          for (int axis = 0; axis < 2; ++axis) {
            const double error = std::abs(actual[axis] - expected[axis]);
            ++gradients;
            gradient_differences += error != 0;
            gradient_max_error = std::max(gradient_max_error, error);
            // Five rounded additions/subtractions; scaling F32 inputs by two is
            // exact in F64. Allow a conservative bound for permitted reassociation.
            const double limit = (pattern == 0 || pattern >= 4) ? 0 : 8 * std::numeric_limits<double>::epsilon() * magnitude[axis];
            check(std::isfinite(actual[axis]) && error <= limit,
                  "Sobel w=" + std::to_string(width) + " h=" + std::to_string(height));
          }
        }
      check(gx.front() == 12345.5 && gx.back() == 12345.5 && gy.front() == 12345.5 && gy.back() == 12345.5,
            "Sobel output guards");
    }
}
} // namespace

int main() {
  try {
    check_plane_views();
    const auto targets = hwy::SupportedAndGeneratedTargets();
    for (auto requested : {HWY_SSE4, HWY_AVX2, HWY_AVX3})
      if (std::find(targets.begin(), targets.end(), requested) == targets.end())
        std::cout << "Skipping unavailable target " << hwy::TargetName(requested) << '\n';
    for (auto target : targets) {
      // Normal x86 builds retain the reduced SIMD selection. Dedicated
      // builds exercise SCALAR and EMU128 without duplicating them here.
#if HWY_ARCH_X86 && !defined(HWY_COMPILE_ONLY_SCALAR) && !defined(HWY_COMPILE_ONLY_EMU128)
      if (target != HWY_SSE4 && target != HWY_AVX2 && target != HWY_AVX3)
        continue;
#endif
      hwy::SetSupportedTargetsForTest(target);
      std::cout << "Testing " << hwy::TargetName(target) << std::endl;
      unpack_tests();
      pack_tests();
      sobel_tests();
      std::cout << hwy::TargetName(target) << ": checks=" << checks << " converted=" << converted
                << " gradients=" << gradients << " gradient_differences=" << gradient_differences
                << " gradient_max_error=" << gradient_max_error << std::endl;
    }
    if (!checks)
      throw std::runtime_error("No supported F3kdb test target executed");
    hwy::SetSupportedTargetsForTest(0);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
