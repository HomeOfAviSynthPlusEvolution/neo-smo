// Included inside the per-target Highway namespace. No include guard.
#include "common/fp16_simd.hpp"

// Keep SIMD operations in named functions so each Highway target supplies
// the matching function attributes. Half arithmetic rounds after each operation.
template <bool IsF16, class D>
inline hn::Vec<D> float_add(D d, hn::Vec<D> a, hn::Vec<D> b) {
  const auto r = hn::Add(a, b);
  if constexpr (IsF16) return round_f16(d, r);
  else return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_sub(D d, hn::Vec<D> a, hn::Vec<D> b) {
  const auto r = hn::Sub(a, b);
  if constexpr (IsF16) return round_f16(d, r);
  else return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_mul(D d, hn::Vec<D> a, hn::Vec<D> b) {
  const auto r = hn::Mul(a, b);
  if constexpr (IsF16) return round_f16(d, r);
  else return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_div9(D d, hn::Vec<D> a) {
  const auto r = hn::Div(a, hn::Set(d, 9.0f));
  if constexpr (IsF16) return round_f16(d, r);
  else return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_div(D d, hn::Vec<D> a, hn::Vec<D> b) {
  const auto r = hn::Div(a, b);
  if constexpr (IsF16) return round_f16(d, r);
  else return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_sqrt(D d, hn::Vec<D> a) {
  const auto r = hn::Sqrt(a);
  if constexpr (IsF16) return round_f16(d, r);
  else return r;
}

template <bool IsF16, class D>
inline hn::Vec<D> float_abs_diff(D d, hn::Vec<D> a, hn::Vec<D> b) {
  return hn::Abs(float_sub<IsF16>(d, a, b));
}
