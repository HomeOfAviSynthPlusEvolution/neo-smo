#include "kernels/dct_pixels.hpp"
#include "base/fp16.hpp"
#include <algorithm>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/dct_pixels.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;
#include "common/fp16_simd.hpp"

template <typename T>
void load_integer(const std::uint8_t* src, float* dst, std::size_t width) {
  const hn::ScalableTag<float> d;
  const hn::Rebind<std::int32_t, decltype(d)> di;
  const hn::Rebind<T, decltype(d)> ds;
  for (std::size_t x = 0; x < width; x += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), width - x);
    hn::StoreN(hn::ConvertTo(d, hn::PromoteTo(di, hn::LoadN(ds, reinterpret_cast<const T*>(src) + x, n))), d, dst + x,
               n);
  }
}

template <typename T>
void store_integer(const float* src, std::uint8_t* dst, std::size_t width, int bits) {
  const hn::ScalableTag<float> d;
  const hn::Rebind<std::int32_t, decltype(d)> di;
  const hn::Rebind<T, decltype(d)> ds;
  for (std::size_t x = 0; x < width; x += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), width - x);
    auto v = hn::LoadN(d, src + x, n);
    v = hn::Min(hn::Max(v, hn::Zero(d)), hn::Set(d, static_cast<float>((1u << bits) - 1)));
    // Zig's @round rounds positive ties upward, not to even.
    // Compare the fractional part to avoid rounding a value just below 0.5
    // up to 1.0 while adding 0.5 in single precision.
    const auto whole = hn::ConvertTo(di, v);
    const auto up = hn::RebindMask(di, hn::Ge(hn::Sub(v, hn::ConvertTo(d, whole)), hn::Set(d, 0.5f)));
    const auto rounded = hn::Add(whole, hn::IfThenElse(up, hn::Set(di, 1), hn::Zero(di)));
    hn::StoreN(hn::DemoteTo(ds, rounded), ds, reinterpret_cast<T*>(dst) + x, n);
  }
}

void LoadDctRow(DataType type, const std::uint8_t* src, float* dst, std::size_t width) {
  if (type == DataType::U8)
    return load_integer<std::uint8_t>(src, dst, width);
  if (type == DataType::U16)
    return load_integer<std::uint16_t>(src, dst, width);
  const hn::ScalableTag<float> d;
  for (std::size_t x = 0; x < width; x += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), width - x);
    const auto v = type == DataType::F16 ? load_f16(d, reinterpret_cast<const std::uint16_t*>(src) + x, n)
                                         : hn::LoadN(d, reinterpret_cast<const float*>(src) + x, n);
    hn::StoreN(v, d, dst + x, n);
  }
}

void StoreDctRow(DataType type, const float* src, std::uint8_t* dst, std::size_t width, int bits) {
  if (type == DataType::U8)
    return store_integer<std::uint8_t>(src, dst, width, bits);
  if (type == DataType::U16)
    return store_integer<std::uint16_t>(src, dst, width, bits);
  const hn::ScalableTag<float> d;
  for (std::size_t x = 0; x < width; x += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), width - x);
    const auto v = hn::LoadN(d, src + x, n);
    if (type == DataType::F16)
      store_f16(d, v, reinterpret_cast<std::uint16_t*>(dst) + x, n);
    else
      hn::StoreN(v, d, reinterpret_cast<float*>(dst) + x, n);
  }
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(LoadDctRow);
HWY_EXPORT(StoreDctRow);
void dct_load_row(DataType type, const std::uint8_t* src, float* dst, std::size_t width) {
  HWY_DYNAMIC_DISPATCH(LoadDctRow)(type, src, dst, width);
}
void dct_store_row(DataType type, const float* src, std::uint8_t* dst, std::size_t width, int bits) {
  HWY_DYNAMIC_DISPATCH(StoreDctRow)(type, src, dst, width, bits);
}
} // namespace neo_smo
#endif
