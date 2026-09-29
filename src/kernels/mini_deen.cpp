#include "kernels/mini_deen_scalar.hpp"
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
#if HWY_TARGET == HWY_SCALAR || HWY_TARGET == HWY_EMU128
void mini_deen_target(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
                      int width, int height, bool byte_samples, int radius, unsigned threshold) {
  mini_deen_scalar(src, src_stride, dst, dst_stride, width, height, byte_samples, radius, threshold);
}
void mini_deen_float_target(const std::uint8_t* src, std::ptrdiff_t ss, std::uint8_t* dst, std::ptrdiff_t ds, int width,
                            int height, bool half, int radius, float threshold) {
  mini_deen_float_scalar(src, ss, dst, ds, width, height, half, radius, threshold);
}
#else
namespace hn = hwy::HWY_NAMESPACE;
template <class D>
auto rounded_mean(D d, hn::VFromD<D> sum, hn::VFromD<D> count) {
  const hn::Rebind<float, D> df;
  const hn::Rebind<std::int32_t, D> di;
  auto q = hn::BitCast(d, hn::ConvertTo(di, hn::Div(hn::ConvertTo(df, sum), hn::ConvertTo(df, count))));
  q = hn::Sub(q, hn::IfThenElseZero(hn::Gt(hn::Mul(q, count), sum), hn::Set(d, 1)));
  const auto remainder = hn::Sub(sum, hn::Mul(q, count));
  return hn::Add(q, hn::IfThenElseZero(hn::Ge(hn::Add(remainder, remainder), count), hn::Set(d, 1)));
}
#include "kernels/deen_io-inl.hpp"

#include "common/fp16_simd.hpp"
#include "kernels/deen_float_io-inl.hpp"

template <bool Half>
void float_process(const std::uint8_t* src, std::ptrdiff_t ss, std::uint8_t* dst, std::ptrdiff_t ds, int width,
                   int height, int radius, float threshold) {
  const hn::ScalableTag<float> d;
  const hn::Rebind<std::conditional_t<Half, std::uint16_t, std::uint32_t>, decltype(d)> du;
  constexpr std::size_t bytes = Half ? 2 : 4;
  const int lanes = static_cast<int>(hn::Lanes(d));
  const auto exponent = hn::Set(du, Half ? 0x7c00u : 0x7f800000u);
  for (int y = 0; y < height; ++y)
    for (int x = 0; x < width;) {
      const int active = std::min(lanes, width - x);
      const auto bits = deen_raw_load(du, src + y * ss + std::size_t(x) * bytes, active);
      if (!hn::AllFalse(du, hn::Eq(hn::And(bits, exponent), exponent)))
        throw std::invalid_argument("MiniDeen: non-finite sample");
      x += active;
    }
  for (int y = 0; y < height; ++y) {
    int x = 0;
    while (x < width) {
      if (x >= radius && width - x - radius >= lanes) {
        const auto center = deen_float_load<Half, true>(d, src + y * ss, width, x, lanes);
        auto sum = hn::Zero(d), count = hn::Set(d, 2.0f), low = center, high = center;
        for (int dy = -std::min(y, radius); dy <= std::min(radius, height - 1 - y); ++dy)
          for (int dx = -radius; dx <= radius; ++dx) {
            const auto value = deen_float_load<Half, true>(d, src + (y + dy) * ss, width, x + dx, lanes);
            const auto diff = hn::Sub(value, center);
            const auto pass = hn::Lt(hn::Abs(diff), hn::Set(d, threshold));
            sum = hn::Add(sum, hn::IfThenElseZero(pass, diff));
            count = hn::Add(count, hn::IfThenElseZero(pass, hn::Set(d, 1.0f)));
            const auto accepted = hn::IfThenElse(pass, value, center);
            low = hn::Min(low, accepted);
            high = hn::Max(high, accepted);
          }
        const auto result = hn::Min(high, hn::Max(low, hn::Add(center, hn::Div(sum, count))));
        deen_float_store<Half>(d, result, dst + y * ds + std::size_t(x) * bytes, lanes);
        x += lanes;
      } else {
        mini_deen_detail::float_store(
            dst + y * ds, x, Half,
            mini_deen_detail::float_pixel(src, ss, width, height, Half, radius, threshold, x, y));
        ++x;
      }
    }
  }
}
void mini_deen_float_target(const std::uint8_t* src, std::ptrdiff_t ss, std::uint8_t* dst, std::ptrdiff_t ds, int width,
                            int height, bool half, int radius, float threshold) {
  if (half)
    float_process<true>(src, ss, dst, ds, width, height, radius, threshold);
  else
    float_process<false>(src, ss, dst, ds, width, height, radius, threshold);
}

template <class T>
void process(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
             int width, int height, int radius, unsigned threshold) {
  // Select full input vectors, then widen into two independent sums.
  // At most 227 contributions: U8 sums fit U16; U16 sums fit U32.
  using Acc = std::conditional_t<sizeof(T) == 1, std::uint16_t, std::uint32_t>;
  const hn::ScalableTag<T> narrow;
  const hn::Repartition<Acc, decltype(narrow)> d;
#if !HWY_ARCH_X86
  const hn::Half<decltype(narrow)> half;
#endif
  const int lanes = static_cast<int>(hn::Lanes(narrow));
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
        const auto center = deen_raw_load(narrow, row + static_cast<std::size_t>(x) * sizeof(T), active);
        const auto widen_lo = [&](auto value) HWY_ATTR {
#if HWY_ARCH_X86
          return hn::BitCast(d, hn::InterleaveLower(narrow, value, hn::Zero(narrow)));
#else
          return hn::PromoteTo(d, hn::LowerHalf(half, value));
#endif
        };
        const auto widen_hi = [&](auto value) HWY_ATTR {
#if HWY_ARCH_X86
          return hn::BitCast(d, hn::InterleaveUpper(narrow, value, hn::Zero(narrow)));
#else
          return hn::PromoteTo(d, hn::UpperHalf(half, value));
#endif
        };
        auto lo = hn::ShiftLeft<1>(widen_lo(center)), hi = hn::ShiftLeft<1>(widen_hi(center));
        auto count = hn::Set(narrow, 2);
        // Inclusive interval threshold-1 implements the strict integer test.
        const auto limit = hn::Set(narrow, threshold - 1);
        const auto lower = hn::SaturatedSub(center, limit), upper = hn::SaturatedAdd(center, limit);
        for (int dy = -std::min(y, radius); dy <= std::min(radius, height - 1 - y); ++dy) {
          const auto* neighbor_row = src + (y + dy) * src_stride;
          for (int dx = -radius; dx <= radius; ++dx) {
            const auto sample =
                deen_raw_load(narrow, neighbor_row + static_cast<std::size_t>(x + dx) * sizeof(T), active);
            const auto outside = hn::Or(hn::SaturatedSub(sample, upper), hn::SaturatedSub(lower, sample));
            const auto pass = hn::Eq(outside, hn::Zero(narrow));
            const auto value = hn::IfThenElseZero(pass, sample);
            lo = hn::Add(lo, widen_lo(value));
            hi = hn::Add(hi, widen_hi(value));
            count = hn::Add(count, hn::IfThenElseZero(pass, hn::Set(narrow, 1)));
          }
        }
        const auto mean = [&](auto sum, auto counts) HWY_ATTR {
          if constexpr (sizeof(Acc) == 2) {
            const hn::Half<decltype(d)> ah;
            const hn::Repartition<std::uint32_t, decltype(d)> wide;
            const auto qlo = rounded_mean(wide, hn::PromoteTo(wide, hn::LowerHalf(ah, sum)),
                                          hn::PromoteTo(wide, hn::LowerHalf(ah, counts)));
            const auto qhi = rounded_mean(wide, hn::PromoteTo(wide, hn::UpperHalf(ah, sum)),
                                          hn::PromoteTo(wide, hn::UpperHalf(ah, counts)));
            return hn::Combine(d, hn::DemoteTo(ah, qhi), hn::DemoteTo(ah, qlo));
          } else {
            return rounded_mean(d, sum, counts);
          }
        };
        lo = mean(lo, widen_lo(count));
        hi = mean(hi, widen_hi(count));
#if HWY_ARCH_X86
        const auto result = hn::ReorderDemote2To(narrow, lo, hi);
#else
        const auto result = hn::Combine(narrow, hn::DemoteTo(half, hi), hn::DemoteTo(half, lo));
#endif
        deen_store(narrow, result, out + static_cast<std::size_t>(x) * sizeof(T), active);
        x += active;
        continue;
      }
      const T result = mini_deen_detail::pixel<T>(src, src_stride, width, height, radius, threshold, x, y);
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
#endif
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(mini_deen_target);
HWY_EXPORT(mini_deen_float_target);
void mini_deen_float_kernel(const std::uint8_t* src, std::ptrdiff_t ss, std::uint8_t* dst, std::ptrdiff_t ds, int width,
                            int height, bool half, int radius, float threshold) {
  HWY_DYNAMIC_DISPATCH(mini_deen_float_target)(src, ss, dst, ds, width, height, half, radius, threshold);
}
void mini_deen_kernel(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
                      int width, int height, bool byte_samples, int radius, unsigned threshold) {
  HWY_DYNAMIC_DISPATCH(mini_deen_target)(src, src_stride, dst, dst_stride, width, height, byte_samples, radius,
                                         threshold);
}
} // namespace neo_smo
#endif
