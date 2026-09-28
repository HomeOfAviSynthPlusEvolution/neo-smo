#include "algorithms/dct_filter.hpp"
#include "base/fp16.hpp"
#include "kernels/dct_pixels.hpp"
#include "guarded.hpp"
#include "hwy/targets.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace {
using namespace neo_smo;
void check(bool ok, const char* message) {
  if (!ok)
    throw std::runtime_error(message);
}

std::size_t extended(std::size_t index, std::size_t length) {
  const auto i = static_cast<std::int64_t>(index);
  const auto n = static_cast<std::int64_t>(length);
  return static_cast<std::size_t>(i < n ? i : std::max(std::int64_t{0}, 2 * n - 1 - i));
}

// Independent orthonormal DCT definition, computed in double precision.
void oracle(const double* input, double* output, const std::vector<double>& factors) {
  double matrix[8][8];
  for (int k = 0; k < 8; ++k)
    for (int x = 0; x < 8; ++x)
      matrix[k][x] = (k == 0 ? std::sqrt(0.125) : 0.5) * std::cos((2 * x + 1) * k * std::acos(-1.0) / 16);
  double a[64]{}, b[64]{}, c[64]{};
  for (int y = 0; y < 8; ++y)
    for (int u = 0; u < 8; ++u)
      for (int x = 0; x < 8; ++x)
        a[y * 8 + u] += input[y * 8 + x] * matrix[u][x];
  for (int v = 0; v < 8; ++v)
    for (int u = 0; u < 8; ++u) {
      for (int y = 0; y < 8; ++y)
        b[v * 8 + u] += a[y * 8 + u] * matrix[v][y];
      b[v * 8 + u] *= static_cast<float>(factors[v] * factors[u]);
    }
  for (int y = 0; y < 8; ++y)
    for (int u = 0; u < 8; ++u)
      for (int v = 0; v < 8; ++v)
        c[y * 8 + u] += b[v * 8 + u] * matrix[v][y];
  std::fill(output, output + 64, 0);
  for (int y = 0; y < 8; ++y)
    for (int x = 0; x < 8; ++x)
      for (int u = 0; u < 8; ++u)
        output[y * 8 + x] += c[y * 8 + u] * matrix[u][x];
}

double read(const unsigned char* p, DataType type) {
  if (type == DataType::U8)
    return *p;
  if (type == DataType::F32) {
    float v;
    std::memcpy(&v, p, 4);
    return v;
  }
  std::uint16_t v;
  std::memcpy(&v, p, 2);
  return type == DataType::F16 ? fp16_to_fp32(v) : v;
}

void test_plane(DataType type, int bits, std::size_t width, std::size_t height, const std::vector<double>& factors,
                bool at_end) {
  const bool integer = type == DataType::U8 || type == DataType::U16;
  const std::size_t bytes = type == DataType::U8 ? 1 : (type == DataType::F32 ? 4 : 2);
  const auto ss = width * bytes, ds = (width + 3) * bytes;
  Guarded src(ss * height, at_end), dst(ds * height, at_end);
  std::memset(dst.data, 0xa5, ds * height);
  for (std::size_t i = 0; i < width * height; ++i) {
    const auto n = (i * 7919 + i * i * 13 + 83) % 65536;
    auto* p = src.data + i * bytes;
    if (integer) {
      const auto v = static_cast<std::uint16_t>(n % (1u << bits));
      if (bytes == 1)
        *p = static_cast<std::uint8_t>(v);
      else
        std::memcpy(p, &v, 2);
    } else {
      const auto v = static_cast<float>(n) / 16384.0f - 1.5f;
      if (type == DataType::F32)
        std::memcpy(p, &v, 4);
      else {
        const auto h = fp32_to_fp16(v);
        std::memcpy(p, &h, 2);
      }
    }
  }
  DctFilter filter(factors);
  DctScratch scratch(filter, width);
  scratch.process(filter, type, bits, src.data, dst.data, width, height, ss, ds);
  for (std::size_t by = 0; by < height; by += 8)
    for (std::size_t bx = 0; bx < width; bx += 8) {
      double input[64], expected[64];
      for (std::size_t y = 0; y < 8; ++y)
        for (std::size_t x = 0; x < 8; ++x)
          input[y * 8 + x] = read(src.data + extended(by + y, height) * ss + extended(bx + x, width) * bytes, type);
      oracle(input, expected, factors);
      for (std::size_t y = 0; y < std::min(std::size_t{8}, height - by); ++y)
        for (std::size_t x = 0; x < std::min(std::size_t{8}, width - bx); ++x) {
          auto want = expected[y * 8 + x];
          if (integer)
            want = std::round(std::clamp(want, 0.0, double((1u << bits) - 1)));
          const auto got = read(dst.data + (by + y) * ds + (bx + x) * bytes, type);
          const double tolerance = integer ? 1.0 : (type == DataType::F16 ? 0.002 : 4e-6);
          check(std::isfinite(got) && std::abs(got - want) <= tolerance, "DCT oracle mismatch");
        }
    }
  for (std::size_t y = 0; y < height; ++y)
    for (std::size_t x = ss; x < ds; ++x)
      check(dst.data[y * ds + x] == 0xa5, "DCT wrote row padding");
  const std::vector<unsigned char> first(dst.data, dst.data + ds * height);
  std::memset(dst.data, 0xa5, ds * height);
  scratch.process(filter, type, bits, src.data, dst.data, width, height, ss, ds);
  check(std::memcmp(first.data(), dst.data, ds * height) == 0, "DCT workspace reuse mismatch");
}
} // namespace

int main() {
  try {
    for (const auto bad :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(), -0.1, 1.1}) {
      std::vector<double> factors(8, 1);
      factors[3] = bad;
      bool rejected = false;
      try {
        DctFilter filter(factors);
      } catch (const std::invalid_argument&) {
        rejected = true;
      }
      check(rejected, "DCT accepted invalid factors");
    }
    for (const auto count : {0, 7, 9}) {
      bool rejected = false;
      try {
        DctFilter filter(std::vector<double>(count, 1));
      } catch (const std::invalid_argument&) {
        rejected = true;
      }
      check(rejected, "DCT accepted invalid factor count");
    }
    for (auto target : hwy::SupportedAndGeneratedTargets()) {
      hwy::SetSupportedTargetsForTest(target);
      // Guard both ends; cover intra-block tails and cross-block batches of 8.
      for (bool at_end : {false, true})
        for (const auto width : {1u, 7u, 8u, 9u, 17u, 63u, 64u, 65u, 129u})
          for (auto type : {DataType::U8, DataType::U16, DataType::F16, DataType::F32}) {
            const int bits = type == DataType::U8 ? 8 : (type == DataType::F32 ? 32 : 16);
            test_plane(type, bits, width, width == 1 ? 1 : 11, {1, .9, .8, .7, .5, .3, .1, 0}, at_end);
          }
      for (const auto& factors :
           {std::vector<double>(8, 1), std::vector<double>(8, 0), std::vector<double>{1, 0, 0, 0, 0, 0, 0, 0}})
        for (auto type : {DataType::U8, DataType::U16, DataType::F16, DataType::F32})
          test_plane(type, type == DataType::U8 ? 8 : (type == DataType::F32 ? 32 : 16), 72, 9, factors, true);
      test_plane(DataType::U16, 10, 33, 17, {1, .9, .8, .7, .5, .3, .1, 0}, true);
      const float ties[]{-1, std::nextafter(0.5f, 0.0f), 0.5f, 1.5f, 254.5f, 65535};
      unsigned char out[6]{};
      dct_store_row(DataType::U8, ties, out, 6, 8);
      const unsigned char want[]{0, 0, 1, 2, 255, 255};
      check(std::memcmp(out, want, 6) == 0, "DCT integer rounding");
      std::printf("%s: DCT oracle, guards, formats and rounding passed\n", hwy::TargetName(target));
    }
    hwy::SetSupportedTargetsForTest(0);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
