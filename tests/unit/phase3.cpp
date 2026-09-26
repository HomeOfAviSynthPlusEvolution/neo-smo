#include "base/fp16.hpp"
#include "kernels/dispatch.hpp"
#include "hwy/targets.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
int reflect(int x, int n) {
  if (x < 0) return std::min(-x, n - 1);
  if (x >= n) return std::max(0, 2 * n - 2 - x);
  return x;
}

template <class T>
void check_integer(neo_smo::DataType type, int width) {
  constexpr int height = 7;
  const size_t stride = width * sizeof(T);
  std::vector<T> src(width * height), dst(src.size());
  uint32_t state = 193;
  for (size_t i = 0; i < src.size(); ++i) {
    state = state * 1664525u + 1013904223u;
    src[i] = i % 5 == 0 ? std::numeric_limits<T>::max() : i % 5 == 1 ? 0 : static_cast<T>(state >> 8);
  }
  const auto* s = reinterpret_cast<const uint8_t*>(src.data());
  auto* d = reinterpret_cast<uint8_t*>(dst.data());
  for (int radius = 1; radius <= 3; ++radius) {
    neo_smo::process_inter_quartile_mean_plane(type, radius, s, d, width, height, stride, stride);
    for (int y = 0; y < height; ++y)
      for (int x = 0; x < width; ++x) {
        std::vector<uint32_t> values;
        for (int dy = -radius; dy <= radius; ++dy)
          for (int dx = -radius; dx <= radius; ++dx)
            values.push_back(src[reflect(y + dy, height) * width + reflect(x + dx, width)]);
        std::sort(values.begin(), values.end());
        const size_t n = values.size(), first = n / 4 + 1, last = 3 * n / 4;
        uint64_t sum = 0;
        for (size_t i = first; i < last; ++i) sum += values[i];
        sum += ((values[first - 1] + values[last]) * 3 + 2) / 4;
        const auto expected = static_cast<T>((sum * 2 + n / 2) / n);
        if (dst[y * width + x] != expected) throw std::runtime_error("IQM integer oracle mismatch");
      }
    for (float threshold : {17.75f, static_cast<float>(std::numeric_limits<T>::max())}) {
      neo_smo::process_smart_median_plane(type, radius, threshold, s, d, width, height, stride, stride);
      for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
          std::vector<int64_t> values;
          for (int dy = -radius; dy <= radius; ++dy)
            for (int dx = -radius; dx <= radius; ++dx)
              if (dx || dy) values.push_back(src[reflect(y + dy, height) * width + reflect(x + dx, width)]);
          std::sort(values.begin(), values.end());
          int64_t sum = 0;
          for (auto v : values) sum += v;
          const int64_t average = (sum + values.size() / 2) / values.size();
          const auto left = values[values.size() / 2 - 1], right = values[values.size() / 2];
          const auto dl = left - average, dr = right - average;
          const float variance = std::sqrt(static_cast<float>(dl * dl + dr * dr)) * 13.0f;
          const auto center = src[y * width + x];
          const auto expected = std::round(variance) <= static_cast<unsigned>(threshold)
              ? static_cast<T>(std::clamp<int64_t>(center, left, right)) : center;
          if (dst[y * width + x] != expected) throw std::runtime_error("SmartMedian integer oracle mismatch");
        }
    }
  }
}

template <class T>
void check_threshold(neo_smo::DataType type, float threshold, const std::array<T, 9>& patch, T expected) {
  constexpr int width = 65, height = 9;
  std::vector<T> src(width * height), dst(src.size());
  for (int y = 0; y < height; ++y)
    for (int x = 0; x < width; ++x) src[y * width + x] = patch[(y % 3) * 3 + x % 3];
  neo_smo::process_smart_median_plane(type, 1, threshold, reinterpret_cast<const uint8_t*>(src.data()),
      reinterpret_cast<uint8_t*>(dst.data()), width, height, width * sizeof(T), width * sizeof(T));
  for (int y = 1; y < height - 1; y += 3)
    for (int x = 1; x < width - 1; x += 3)
      if (dst[y * width + x] != expected) throw std::runtime_error("SmartMedian threshold decision mismatch");
}
} // namespace

int main() {
  using neo_smo::DataType;
  for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
    hwy::SetSupportedTargetsForTest(target);
    for (int width : {1, 2, 3, 7, 17, 65}) {
      check_integer<uint8_t>(DataType::U8, width);
      check_integer<uint16_t>(DataType::U16, width);
    }
    check_threshold<uint8_t>(DataType::U8, 17.75f, {0, 0, 0, 0, 255, 2, 2, 2, 2}, 255);
    check_threshold<uint8_t>(DataType::U8, 18.0f, {0, 0, 0, 0, 255, 2, 2, 2, 2}, 2);
    check_threshold<uint16_t>(DataType::U16, 17.75f, {0, 0, 0, 0, 65535, 2, 2, 2, 2}, 65535);
    check_threshold<uint16_t>(DataType::U16, 2704.0f, {724, 724, 724, 996, 65535, 1208, 1208, 1208, 1208}, 65535);
    check_threshold<uint16_t>(DataType::U16, 2705.0f, {724, 724, 724, 996, 65535, 1208, 1208, 1208, 1208}, 1208);
    const uint16_t tenth = neo_smo::fp32_to_fp16(0.1f), one = neo_smo::fp32_to_fp16(1.0f);
    check_threshold<uint16_t>(DataType::F16, 0.9188f, {0, 0, 0, 0, one, tenth, tenth, tenth, tenth}, tenth);
    check_threshold<uint16_t>(DataType::F16, 0.9186f, {0, 0, 0, 0, one, tenth, tenth, tenth, tenth}, one);
    std::printf("%s: phase3 integer oracles and threshold decisions passed\n", hwy::TargetName(target));
  }
  hwy::SetSupportedTargetsForTest(0);
}
