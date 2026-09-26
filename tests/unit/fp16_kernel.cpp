#include <cstddef>
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "unit/fp16_kernel.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace test_fp16 {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;
#include "common/fp16_simd.hpp"
void round_values(const float* in, float* out, size_t count) {
  const hn::ScalableTag<float> d;
  const size_t lanes = hn::Lanes(d);
  for (size_t i = 0; i < count; i += lanes) {
    const size_t n = std::min(lanes, count - i);
    hn::StoreN(round_f16(d, hn::LoadN(d, in + i, n)), d, out + i, n);
  }
}
} // namespace HWY_NAMESPACE
} // namespace test_fp16
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace test_fp16 {
HWY_EXPORT(round_values);
void run(const float* in, float* out, size_t count) {
  HWY_DYNAMIC_DISPATCH(round_values)(in, out, count);
}
} // namespace test_fp16
#endif
