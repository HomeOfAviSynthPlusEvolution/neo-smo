// Included inside the per-target Highway namespace. No include guard.
namespace hn = hwy::HWY_NAMESPACE;
#include "common/fp16_simd.hpp"

// Convert explicitly to the lane type: float16_t need not be a built-in type.
template <class D, class T>
HWY_INLINE hn::Vec<D> float_set(D d, T value) {
  return hn::Set(d, hwy::ConvertScalarTo<hn::TFromD<D>>(value));
}

// Keep SIMD operations in named functions so each Highway target supplies
// the matching function attributes. Native half lanes round in hardware;
// other targets retain FP32 intermediates until the final store.
template <bool IsF16, class D>
inline hn::Vec<D> float_add([[maybe_unused]] D d, hn::Vec<D> a, hn::Vec<D> b) {
  const auto r = hn::Add(a, b);
  return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_sub([[maybe_unused]] D d, hn::Vec<D> a, hn::Vec<D> b) {
  const auto r = hn::Sub(a, b);
  return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_mul([[maybe_unused]] D d, hn::Vec<D> a, hn::Vec<D> b) {
  const auto r = hn::Mul(a, b);
  return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_div9([[maybe_unused]] D d, hn::Vec<D> a) {
  const auto r = hn::Div(a, float_set(d, 9.0f));
  return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_div([[maybe_unused]] D d, hn::Vec<D> a, hn::Vec<D> b) {
  const auto r = hn::Div(a, b);
  return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_sqrt([[maybe_unused]] D d, hn::Vec<D> a) {
  const auto r = hn::Sqrt(a);
  return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_abs_diff([[maybe_unused]] D d, hn::Vec<D> a, hn::Vec<D> b) {
  return hn::Abs(float_sub<IsF16>(d, a, b));
}
