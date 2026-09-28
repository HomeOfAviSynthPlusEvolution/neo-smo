#include "base/fp16.hpp"
#include "kernels/dispatch.hpp"
#include "hwy/targets.h"
#include "hwy/per_target.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {
constexpr size_t width = 65, height = 17, count = width * height;

template <class T>
std::vector<float> apply(const std::vector<float>& input, bool half, bool chroma, int filter, int mode) {
  std::vector<T> src(count), dst(count);
  for (size_t i = 0; i < count; ++i) {
    if constexpr (sizeof(T) == 2)
      src[i] = neo_smo::fp32_to_fp16(input[i]);
    else
      src[i] = input[i];
  }
  const auto type = half ? neo_smo::DataType::F16 : neo_smo::DataType::F32;
  const auto* s = reinterpret_cast<const uint8_t*>(src.data());
  auto* d = reinterpret_cast<uint8_t*>(dst.data());
  const size_t stride = width * sizeof(T);
  if (filter == 0)
    neo_smo::process_median_plane(type, mode, s, d, width, height, stride, stride);
  if (filter == 1)
    neo_smo::process_remove_grain_plane(type, mode, chroma, s, d, width, height, stride, stride);
  if (filter == 2)
    neo_smo::process_vertical_cleaner_plane(type, mode, chroma, 16, s, d, width, height, stride, stride);
  if (filter >= 3 && filter <= 5) {
    // Independent spatial permutations retain the input's subnormals and ties.
    auto ref1 = src, ref2 = src;
    std::rotate(ref1.begin(), ref1.begin() + 13, ref1.end());
    std::rotate(ref2.begin(), ref2.begin() + 37, ref2.end());
    const auto* r1 = reinterpret_cast<const uint8_t*>(ref1.data());
    const auto* r2 = reinterpret_cast<const uint8_t*>(ref2.data());
    if (filter == 3)
      neo_smo::process_repair_plane(type, mode, chroma, s, r1, d, width, height, stride, stride, stride);
    else if (filter == 4)
      neo_smo::process_clense_plane(type, s, r1, r2, d, width, height, stride, stride, stride, stride);
    else
      neo_smo::process_clense_forward_backward_plane(type, s, r1, r2, d, width, height, stride, stride, stride, stride);
  }
  if (filter == 6)
    neo_smo::process_inter_quartile_mean_plane(type, mode, s, d, width, height, stride, stride);
  if (filter == 7)
    neo_smo::process_smart_median_plane(type, mode, 0.2f, s, d, width, height, stride, stride);
  std::vector<float> out(count);
  for (size_t i = 0; i < count; ++i) {
    if constexpr (sizeof(T) == 2)
      out[i] = neo_smo::fp16_to_fp32(dst[i]);
    else
      out[i] = dst[i];
  }
  return out;
}

// Mode 8 can select a different direction when half-rounded costs tie.
// Check each arithmetic policy against its own scalar calculation.
bool check_rg8_precision(const std::vector<int64_t>& targets) {
  std::vector<float> input(count);
  unsigned state = 41;
  for (auto& value : input) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    value = (1000 + state % 48) * 0x1p-24f;
    if (state & 1024)
      value = -value;
    value = neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(value));
  }
  const int offsets[] = {-int(width) - 1, -int(width), -int(width) + 1, -1};
  for (auto target : targets) {
    hwy::SetSupportedTargetsForTest(target);
    const bool native = hwy::HaveFloat16();
    const auto round = [=](float x) {
      return native ? neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(x)) : x;
    };
    const auto output = apply<uint16_t>(input, true, true, 1, 8);
    for (size_t y = 1; y + 1 < height; ++y)
      for (size_t x = 1; x + 1 < width; ++x) {
        const int i = static_cast<int>(y * width + x);
        float clamped[4], costs[4];
        for (int k = 0; k < 4; ++k) {
          const auto lo = std::min(input[i + offsets[k]], input[i - offsets[k]]);
          const auto hi = std::max(input[i + offsets[k]], input[i - offsets[k]]);
          clamped[k] = std::clamp(input[i], lo, hi);
          costs[k] = std::clamp(round(std::abs(round(input[i] - clamped[k])) + round(round(hi - lo) * 2)), -0.5f, 0.5f);
        }
        const auto best = *std::min_element(costs, costs + 4);
        float expected = clamped[0];
        for (int k : {2, 1, 3})
          if (costs[k] == best)
            expected = clamped[k];
        if (output[i] != expected) {
          std::fprintf(stderr, "%s RG8 native_half=%d pixel=%d ref=%g out=%g\n", hwy::TargetName(target), native, i,
                       expected, output[i]);
          return false;
        }
      }
  }
  return true;
}

} // namespace

int main() {
  const auto targets = hwy::SupportedAndGeneratedTargets();
  if (!check_rg8_precision(targets))
    return 1;
  int64_t references[2]{};
  for (auto target : targets) {
    hwy::SetSupportedTargetsForTest(target);
    const int group = hwy::HaveFloat16() ? 1 : 0;
    if (!references[group]) references[group] = target;
  }
  float maxima[2]{};
  // Four deterministic patterns, including signed subnormals and dense ties.
  // Compare three successive applications, not just isolated helper rounding.
  for (bool half : {false, true})
    for (unsigned seed : {41u, 81u, 193u})
      for (int pattern = 0; pattern < 4; ++pattern)
        for (bool chroma : {false, true}) {
          std::vector<float> input(count);
          unsigned state = seed;
          for (size_t i = 0; i < count; ++i) {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            const unsigned v = state;
            if (pattern == 0)
              input[i] = (v % 1024) * 0x1p-24f;
            if (pattern == 1)
              input[i] = (1000 + v % 48) * 0x1p-24f;
            if (pattern == 2)
              input[i] = (v % 5) * 0.25f;
            if (pattern == 3)
              input[i] = (v % 65536) / 65535.0f;
            if (chroma)
              input[i] = pattern < 2 ? ((v & 1024) ? -input[i] : input[i]) : input[i] - 0.5f;
          }
          for (int filter = 0; filter < 8; ++filter)
            for (int mode = 1; mode <= (filter == 0 || filter >= 6 ? 3 : (filter == 1 || filter == 3) ? 24 : filter == 2 ? 2 : 1); ++mode) {
              std::vector<std::vector<float>> expected[2];
              for (int group = 0; group < (half ? 2 : 1); ++group) {
                const auto reference = half ? references[group] : targets.front();
                if (!reference) continue;
                auto current = input;
                hwy::SetSupportedTargetsForTest(reference);
                for (int step = 0; step < 3; ++step) {
                  current = half ? apply<uint16_t>(current, half, chroma, filter, mode)
                                 : apply<float>(current, half, chroma, filter, mode);
                  expected[group].push_back(current);
                }
              }
              for (auto target : targets) {
                hwy::SetSupportedTargetsForTest(target);
                // F16 arithmetic and F32 arithmetic with F16 storage are separate
                // production policies, not numerically identical implementations.
                const int group = half && hwy::HaveFloat16() ? 1 : 0;
                auto current = input;
                for (int step = 0; step < 3; ++step) {
                  current = half ? apply<uint16_t>(current, half, chroma, filter, mode)
                                 : apply<float>(current, half, chroma, filter, mode);
                  for (size_t i = 0; i < count; ++i) {
                    const float ref = expected[group][step][i], value = current[i];
                    const float error = std::abs(value - ref);
                    // Normalized-video budget: about one half ULP for F16;
                    // for F32, low single-digit millionths of full scale.
                    const float limit = half ? 0x1p-20f + 0x1p-10f * std::abs(ref) : 1e-6f + 1e-6f * std::abs(ref);
                    if (!std::isfinite(ref) || !std::isfinite(value) || error > limit) {
                      std::fprintf(stderr,
                                   "%s half=%d seed=%u pattern=%d chroma=%d filter=%d mode=%d step=%d pixel=%zu ref=%g "
                                   "out=%g error=%g limit=%g\n",
                                   hwy::TargetName(target), half, seed, pattern, chroma, filter, mode, step, i, ref,
                                   value, error, limit);
                      return 1;
                    }
                    maxima[half] = std::max(maxima[half], error);
                  }
                }
              }
            }
        }
  hwy::SetSupportedTargetsForTest(0);
  std::printf("%zu targets; 3 successive passes within each arithmetic policy; max_abs F32=%g F16=%g\n",
              targets.size(), maxima[0], maxima[1]);
}
