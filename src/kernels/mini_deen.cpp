// Algorithm: dubhater/vapoursynth-minideen (ISC), src/minideen.cpp,
// revision 8b915c5134352adc19539560ab5debbbcdbe0ea1.
// Retains the scalar formula, strict threshold, clipped window and extra center weight.
#include "kernels/mini_deen.hpp"
#include <algorithm>
#include <cstring>
#include <type_traits>
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/mini_deen.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"
HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;
template <class T>
T load(const std::uint8_t* row, int x) {
  T value;
  std::memcpy(&value, row + static_cast<std::size_t>(x) * sizeof(T), sizeof(T));
  return value;
}
template <class D>
auto rounded_mean(D d, hn::VFromD<D> sum, hn::VFromD<D> count) {
  const hn::Rebind<float, D> df;
  const hn::Rebind<std::int32_t, D> di;
  auto q = hn::BitCast(d, hn::ConvertTo(di, hn::Div(hn::ConvertTo(df, sum), hn::ConvertTo(df, count))));
  q = hn::Sub(q, hn::IfThenElseZero(hn::Gt(hn::Mul(q, count), sum), hn::Set(d, 1)));
  const auto remainder = hn::Sub(sum, hn::Mul(q, count));
  return hn::Add(q, hn::IfThenElseZero(hn::Ge(hn::Add(remainder, remainder), count), hn::Set(d, 1)));
}
template <class D>
auto load_samples(D d, const hn::TFromD<D>* ptr, std::size_t active) {
  return active == hn::Lanes(d) ? hn::LoadU(d, ptr) : hn::LoadN(d, ptr, active);
}
template <class T>
void process(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
             int width, int height, int radius, unsigned threshold) {
  // 227*255=57885 fits uint16, doubling the lane count for byte input.
  using Acc = std::conditional_t<sizeof(T) == 1 && HWY_TARGET != HWY_SCALAR, std::uint16_t, std::uint32_t>;
  const hn::ScalableTag<Acc> d;
  const hn::Rebind<T, decltype(d)> narrow;
  const int lanes = static_cast<int>(hn::Lanes(d));
  for (int y = 0; y < height; ++y) {
    const auto* row = src + y * src_stride;
    auto* out = dst + y * dst_stride;
    // Test AFTER scaling: threshold=1 can filter high-bit-depth input.
    if (threshold <= 1) {
      std::memcpy(out, row, static_cast<std::size_t>(width) * sizeof(T));
      continue;
    }
    for (int x = 0; x < width;) {
      if (x >= radius && width - x - radius > 0) {
        const int active = std::min(lanes, width - x - radius);
        const auto center = hn::PromoteTo(
            d, load_samples(narrow, reinterpret_cast<const T*>(row + static_cast<std::size_t>(x) * sizeof(T)), active));
        auto sum = hn::Add(center, center), count = hn::Set(d, 2);
        for (int dy = -std::min(y, radius); dy <= std::min(radius, height - 1 - y); ++dy) {
          const auto* neighbor_row = src + (y + dy) * src_stride;
          for (int dx = -radius; dx <= radius; ++dx) {
            const auto sample = hn::PromoteTo(
                d, load_samples(narrow,
                                reinterpret_cast<const T*>(neighbor_row + static_cast<std::size_t>(x + dx) * sizeof(T)),
                                active));
            const auto difference = hn::Sub(hn::Max(center, sample), hn::Min(center, sample));
            const auto pass = hn::Lt(difference, hn::Set(d, threshold));
            sum = hn::Add(sum, hn::IfThenElseZero(pass, sample));
            count = hn::Add(count, hn::IfThenElseZero(pass, hn::Set(d, 1)));
          }
        }
        // Float division gives a quotient estimate; integer correction makes half-up exact.
        hn::VFromD<decltype(d)> q;
        if constexpr (sizeof(Acc) == 2) {
          const hn::Half<decltype(d)> half;
          const hn::Repartition<std::uint32_t, decltype(d)> wide;
          const auto lo = rounded_mean(wide, hn::PromoteTo(wide, hn::LowerHalf(half, sum)),
                                       hn::PromoteTo(wide, hn::LowerHalf(half, count)));
          const auto hi = rounded_mean(wide, hn::PromoteTo(wide, hn::UpperHalf(half, sum)),
                                       hn::PromoteTo(wide, hn::UpperHalf(half, count)));
          q = hn::Combine(d, hn::DemoteTo(half, hi), hn::DemoteTo(half, lo));
        } else {
          q = rounded_mean(d, sum, count);
        }
        hn::StoreN(hn::DemoteTo(narrow, q), narrow, reinterpret_cast<T*>(out + static_cast<std::size_t>(x) * sizeof(T)),
                   active);
        x += active;
        continue;
      }
      const unsigned center = load<T>(row, x);
      unsigned sum = 2 * center, count = 2;
      // Clip offsets, not coordinates: avoid overflow even for large dimensions.
      for (int dy = -std::min(y, radius); dy <= std::min(radius, height - 1 - y); ++dy) {
        const auto* neighbor_row = src + (y + dy) * src_stride;
        for (int dx = -std::min(x, radius); dx <= std::min(radius, width - 1 - x); ++dx) {
          const unsigned sample = load<T>(neighbor_row, x + dx);
          const unsigned difference = center > sample ? center - sample : sample - center;
          if (difference < threshold) {
            sum += sample;
            ++count;
          }
        }
      }
      // At most 227 * 65535; both the sum and rounding numerator fit uint32.
      const T result = static_cast<T>((2 * sum + count) / (2 * count));
      std::memcpy(out + static_cast<std::size_t>(x) * sizeof(T), &result, sizeof(T));
      ++x;
    }
  }
}
void mini_deen_target(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
                      int width, int height, bool byte_samples, int radius, unsigned threshold) {
  if (byte_samples)
    process<std::uint8_t>(src, src_stride, dst, dst_stride, width, height, radius, threshold);
  else
    process<std::uint16_t>(src, src_stride, dst, dst_stride, width, height, radius, threshold);
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(mini_deen_target);
void mini_deen_kernel(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
                      int width, int height, bool byte_samples, int radius, unsigned threshold) {
  HWY_DYNAMIC_DISPATCH(mini_deen_target)(src, src_stride, dst, dst_stride, width, height, byte_samples, radius,
                                         threshold);
}
} // namespace neo_smo
#endif
