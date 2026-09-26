#include "base/fp16.hpp"
#include "kernels/dispatch.hpp"
#include "hwy/targets.h"
#include "hwy/per_target.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {
float half(float f) { return neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(f)); }
template<class T>
void temporal(neo_smo::DataType type, int width) {
  constexpr int height = 3;
  const size_t size = width * height, pitch = width * sizeof(T);
  const bool f16 = type == neo_smo::DataType::F16;
  const bool integral = type == neo_smo::DataType::U8 || type == neo_smo::DataType::U16;
  auto decode = [&](T v) { return f16 ? neo_smo::fp16_to_fp32(static_cast<uint16_t>(v)) : static_cast<float>(v); };
  auto encode = [&](float v) { return f16 ? static_cast<T>(neo_smo::fp32_to_fp16(v)) : static_cast<T>(v); };
  std::array<std::vector<T>, 21> frames;
  std::array<const uint8_t*, 21> pointers{};
  std::vector<T> dst(size);
  uint32_t state = 91;
  for (size_t j = 0; j < frames.size(); ++j) {
    frames[j].resize(size);
    for (size_t i = 0; i < size; ++i) {
      state = 1664525 * state + 1013904223;
      float value = integral ? static_cast<float>((state >> 8) % (sizeof(T) == 1 ? 256 : 65536)) :
                              static_cast<float>(static_cast<int>((state >> 8) % 2049) - 1024) / 2048;
      frames[j][i] = encode(value);
    }
    pointers[j] = reinterpret_cast<const uint8_t*>(frames[j].data());
  }
  for (int count = 1; count <= 21; ++count) {
    neo_smo::process_temporal_median_plane(type, count, pointers.data(), reinterpret_cast<uint8_t*>(dst.data()), width, height, pitch, pitch);
    for (size_t i = 0; i < size; ++i) {
      std::vector<float> values;
      for (int j = 0; j < count; ++j) values.push_back(decode(frames[j][i]));
      std::sort(values.begin(), values.end());
      float expected = values[count / 2];
      if (count % 2 == 0) {
        float sum = values[count / 2 - 1] + expected;
        if (f16 && hwy::HaveFloat16()) sum = half(sum);
        expected = sum * 0.5f;
        if (integral) expected = std::floor(expected);
      }
      if (dst[i] != encode(expected)) {
        std::fprintf(stderr, "median type=%d count=%d width=%d got=%g expected=%g\n", int(type), count, width, decode(dst[i]), expected);
        throw std::runtime_error("TemporalMedian sorted oracle mismatch");
      }
    }
    const float parameter = integral ? 73.75f : 0.18337f;
    const float threshold = integral ? std::floor(parameter) : f16 && hwy::HaveFloat16() ? half(parameter) : parameter;
    neo_smo::process_temporal_soften_plane(type, count, parameter, pointers.data(), reinterpret_cast<uint8_t*>(dst.data()), width, height, pitch, pitch);
    for (size_t i = 0; i < size; ++i) {
      const float current = decode(frames[0][i]);
      float sum = current;
      for (int j = 1; j < count; ++j) {
        const float v = decode(frames[j][i]);
        sum += std::abs(current - v) <= threshold ? v : current;
        if (f16 && hwy::HaveFloat16()) sum = half(sum);
      }
      float expected = sum / count;
      if (integral) {
        const uint64_t scale = uint64_t{1} << (sizeof(T) == 1 ? 16 : 32);
        expected = static_cast<float>((static_cast<uint64_t>(sum + count / 2) * (scale / count)) / scale);
      }
      if (std::abs(decode(dst[i]) - decode(encode(expected))) > (integral || f16 ? 0 : 1e-6f))
        throw std::runtime_error("TemporalSoften scalar oracle mismatch");
    }
  }
}
}
int main() {
  for (const auto target : hwy::SupportedAndGeneratedTargets()) {
    hwy::SetSupportedTargetsForTest(target);
    for (const int width : {1, 2, 7, 16, 17, 31, 32, 33, 63, 64, 65}) {
      temporal<uint8_t>(neo_smo::DataType::U8, width);
      temporal<uint16_t>(neo_smo::DataType::U16, width);
      temporal<uint16_t>(neo_smo::DataType::F16, width);
      temporal<float>(neo_smo::DataType::F32, width);
    }
    std::printf("%s: Phase4 temporal oracles passed\n", hwy::TargetName(target));
  }
  hwy::SetSupportedTargetsForTest(0);
}
