#include <cstddef>
#include "common/padded_row.hpp"
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
void decode_values(const std::uint16_t* in, float* out, size_t count) {
  const hn::ScalableTag<float> d;
  for (size_t i = 0; i < count; i += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), count - i);
    hn::StoreN(load_f16(d, in + i, n), d, out + i, n);
  }
}
void encode_values(const float* in, std::uint16_t* out, size_t count) {
  const hn::ScalableTag<float> d;
  for (size_t i = 0; i < count; i += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), count - i);
    store_f16(d, hn::LoadN(d, in + i, n), out + i, n);
  }
}
} // namespace HWY_NAMESPACE
} // namespace test_fp16
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace test_fp16 {
HWY_EXPORT(round_values);
HWY_EXPORT(decode_values);
HWY_EXPORT(encode_values);
void decode(const std::uint16_t* in, float* out, size_t count) {
  HWY_DYNAMIC_DISPATCH(decode_values)(in, out, count);
}
void encode(const float* in, std::uint16_t* out, size_t count) {
  HWY_DYNAMIC_DISPATCH(encode_values)(in, out, count);
}
void run(const float* in, float* out, size_t count) {
  HWY_DYNAMIC_DISPATCH(round_values)(in, out, count);
}
} // namespace test_fp16
#endif
