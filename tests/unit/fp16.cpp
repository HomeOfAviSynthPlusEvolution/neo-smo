#include "base/fp16.hpp"
#include "guarded.hpp"
#include "hwy/targets.h"
#include "hwy/per_target.h"
#include <cstring>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace test_fp16 {
void run(const float* in, float* out, size_t count);
void decode(const std::uint16_t* in, float* out, size_t count);
void encode(const float* in, std::uint16_t* out, size_t count);
}

int main(int argc, char** argv) {
  if (argc == 2 && std::strcmp(argv[1], "--capability") == 0) {
    std::puts(hwy::HaveFloat16() ? "native" : "fp32");
    return 0;
  }
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
    std::vector<std::uint16_t> codes(65536), encoded(inputs.size());
    std::vector<float> decoded(65536);
    for (size_t i = 0; i < codes.size(); ++i) codes[i] = static_cast<std::uint16_t>(i);
    test_fp16::decode(codes.data(), decoded.data(), codes.size());
    for (size_t i = 0; i < codes.size(); ++i) {
      const auto e = neo_smo::fp16_to_fp32(codes[i]);
      if (!(std::isnan(e) ? std::isnan(decoded[i]) :
            decoded[i] == e && std::signbit(decoded[i]) == std::signbit(e))) {
        std::fprintf(stderr, "%s decode half=%04x failed\n", hwy::TargetName(target), unsigned(i));
        return 1;
      }
    }
    test_fp16::encode(inputs.data(), encoded.data(), inputs.size());
    for (size_t i = 0; i < inputs.size(); ++i) {
      const auto h = encoded[i];
      const bool ok = std::isnan(inputs[i]) ? ((h & 0x7c00) == 0x7c00 && (h & 0x3ff))
                                           : h == neo_smo::fp32_to_fp16(inputs[i]);
      if (!ok) {
        std::fprintf(stderr, "%s encode input=%a half=%04x failed\n", hwy::TargetName(target), inputs[i], unsigned(h));
        return 1;
      }
    }
    for (bool end : {false, true}) for (size_t n : {1, 3, 7, 15, 16, 17, 31, 33, 65}) {
      Guarded half(n * sizeof(std::uint16_t), end), floats(n * sizeof(float), end);
      auto* hp = reinterpret_cast<std::uint16_t*>(half.data);
      auto* fp = reinterpret_cast<float*>(floats.data);
      for (size_t i = 0; i < n; ++i) hp[i] = static_cast<std::uint16_t>(i + 1);
      test_fp16::decode(hp, fp, n);
      test_fp16::encode(fp, hp, n);
      for (size_t i = 0; i < n; ++i) if (hp[i] != i + 1) return 1;
    }
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
