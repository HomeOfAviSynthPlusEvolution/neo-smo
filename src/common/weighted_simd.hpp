// Included inside each Highway target namespace. No include guard.
#include "common/float_arithmetic.hpp"

// Loads never cross a row boundary; only FP16 conversion and mirrored borders
// need staging. The arithmetic remains vectorized for full blocks and tails.
template <bool Half, class D, class T>
HWY_INLINE hn::Vec<D> weighted_load(D d, const T* ptr, std::size_t count) {
  const hn::Rebind<float, D> df;
  if constexpr (Half) {
    HWY_ALIGN float values[hn::MaxLanes(df)]{};
    for (std::size_t i = 0; i < count; ++i) values[i] = fp16_to_fp32(ptr[i]);
    return hn::LoadU(d, values);
  } else if constexpr (std::is_same_v<T, float>) {
    return hn::LoadN(d, ptr, count);
  } else {
    const hn::Rebind<T, D> ds;
    const hn::Rebind<std::uint32_t, D> du;
    const auto value = hn::ConvertTo(df, hn::PromoteTo(du, hn::LoadN(ds, ptr, count)));
    if constexpr (std::is_same_v<hn::TFromD<D>, double>) return hn::PromoteTo(d, value);
    else return value;
  }
}

template <class D>
HWY_INLINE hn::Vec<D> weighted_round(D d, hn::Vec<D> value) {
  const auto floor = hn::Floor(value);
  return hn::Add(floor, hn::IfThenElse(hn::Ge(hn::Sub(value, floor), hn::Set(d, 0.5)),
                                     hn::Set(d, 1), hn::Zero(d)));
}

template <bool Half, class D, class T>
HWY_INLINE void weighted_store(D d, hn::Vec<D> value, T* ptr, std::size_t count) {
  HWY_ALIGN hn::TFromD<D> values[hn::MaxLanes(d)];
  hn::StoreU(value, d, values);
  for (std::size_t i = 0; i < count; ++i) {
    if constexpr (Half) ptr[i] = fp32_to_fp16(values[i]);
    else ptr[i] = static_cast<T>(values[i]);
  }
}
