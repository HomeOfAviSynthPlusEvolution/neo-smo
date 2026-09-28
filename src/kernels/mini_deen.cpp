// Algorithm: dubhater/vapoursynth-minideen (ISC), src/minideen.cpp,
// revision 8b915c5134352adc19539560ab5debbbcdbe0ea1.
// Retains the scalar formula, strict threshold, clipped window and extra center weight.
#include "kernels/mini_deen.hpp"
#include <algorithm>
#include <cstring>
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
template <class T>
void process(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
             int width, int height, int radius, unsigned threshold) {
  const hn::ScalableTag<std::uint32_t> d;
  const hn::Rebind<T, decltype(d)> narrow;
  const hn::Rebind<float, decltype(d)> df;
  const hn::Rebind<std::int32_t, decltype(d)> di;
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
      if (x >= radius && width - x - radius >= lanes) {
        const auto center = hn::PromoteTo(
            d, hn::LoadU(narrow, reinterpret_cast<const T*>(row + static_cast<std::size_t>(x) * sizeof(T))));
        auto sum = hn::Add(center, center), count = hn::Set(d, 2);
        for (int dy = -std::min(y, radius); dy <= std::min(radius, height - 1 - y); ++dy) {
          const auto* neighbor_row = src + (y + dy) * src_stride;
          for (int dx = -radius; dx <= radius; ++dx) {
            const auto sample = hn::PromoteTo(
                d, hn::LoadU(narrow,
                             reinterpret_cast<const T*>(neighbor_row + static_cast<std::size_t>(x + dx) * sizeof(T))));
            const auto difference = hn::Sub(hn::Max(center, sample), hn::Min(center, sample));
            const auto pass = hn::Lt(difference, hn::Set(d, threshold));
            sum = hn::Add(sum, hn::IfThenElseZero(pass, sample));
            count = hn::Add(count, hn::IfThenElseZero(pass, hn::Set(d, 1)));
          }
        }
        // Float division gives a quotient estimate; integer correction makes half-up exact.
        auto q = hn::BitCast(d, hn::ConvertTo(di, hn::Div(hn::ConvertTo(df, sum), hn::ConvertTo(df, count))));
        q = hn::Sub(q, hn::IfThenElseZero(hn::Gt(hn::Mul(q, count), sum), hn::Set(d, 1)));
        const auto remainder = hn::Sub(sum, hn::Mul(q, count));
        q = hn::Add(q, hn::IfThenElseZero(hn::Ge(hn::Add(remainder, remainder), count), hn::Set(d, 1)));
        hn::StoreU(hn::DemoteTo(narrow, q), narrow,
                   reinterpret_cast<T*>(out + static_cast<std::size_t>(x) * sizeof(T)));
        x += lanes;
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
