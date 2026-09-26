#include "base/fp16.hpp"
#include "hwy/targets.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace test_fp16 {
void run(const float* in, float* out, size_t count);
}

int main() {
  std::vector<float> inputs, expected;
  const auto add = [&](float v, float e) {
    inputs.push_back(v);
    expected.push_back(e);
    inputs.push_back(-v);
    expected.push_back(-e);
  };
  // Every positive finite half value and each rounding boundary, including
  // the subnormal/normal transition. Expected endpoints are known half codes.
  for (unsigned h = 0; h < 0x7bff; ++h) {
    const float a = neo_smo::fp16_to_fp32(static_cast<uint16_t>(h));
    const float b = neo_smo::fp16_to_fp32(static_cast<uint16_t>(h + 1));
    const float mid = (a + b) * 0.5f;
    add(a, a);
    add(std::nextafter(mid, a), a);
    add(mid, (h & 1) ? b : a);
    add(std::nextafter(mid, b), b);
  }
  add(65504.0f, 65504.0f);
  add(65519.0f, 65504.0f);
  add(65520.0f, INFINITY);
  add(INFINITY, INFINITY);
  add(NAN, NAN);
  std::vector<float> out(inputs.size());
  for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
    hwy::SetSupportedTargetsForTest(target);
    test_fp16::run(inputs.data(), out.data(), inputs.size());
    size_t inexact = 0;
    float max_error = 0.0f;
    for (size_t i = 0; i < inputs.size(); ++i) {
      bool ok;
      if (std::isnan(expected[i])) {
        ok = std::isnan(out[i]);
      } else if (std::isinf(expected[i])) {
        ok = out[i] == expected[i];
      } else {
        const float error = std::abs(out[i] - expected[i]);
        ok = std::isfinite(out[i]) && error == 0.0f && std::signbit(out[i]) == std::signbit(expected[i]);
        if (out[i] != expected[i])
          ++inexact;
        max_error = std::max(max_error, error);
      }
      if (!ok) {
        std::fprintf(stderr, "%s input=%a output=%a expected=%a\n", hwy::TargetName(target), inputs[i], out[i],
                     expected[i]);
        return 1;
      }
    }
    std::printf("%s: %zu half rounding cases passed; inexact=%zu max_abs=%g\n", hwy::TargetName(target), inputs.size(),
                inexact, max_error);
  }
  hwy::SetSupportedTargetsForTest(0);
}
