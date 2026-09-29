#include "kernels/deen_dispatch.hpp"
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/deen_w.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"
HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
#include "kernels/deen_simd-inl.hpp"
#include "kernels/deen_io-inl.hpp"
#include "common/fp16_simd.hpp"

template <bool Half, bool Interior, class D>
auto weighted_load(D d, const std::uint8_t* row, int width, std::int64_t x, std::size_t active) {
  if constexpr (!Half)
    return deen_row_load<Interior>(d, row, width, x, active);
  else {
    const hn::Rebind<std::uint16_t, D> du;
    const auto bits = deen_row_load<Interior>(du, row, width, x, active);
#if HWY_HAVE_FLOAT16 || (HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX2 && !defined(HWY_DISABLE_F16C))
    const hn::Rebind<hwy::float16_t, D> dh;
    return hn::PromoteTo(d, hn::BitCast(dh, bits));
#else
    HWY_ALIGN std::uint16_t values[hn::MaxLanes(d)]{};
    hn::StoreU(bits, du, values);
    return load_f16(d, values, active);
#endif
  }
}

template <bool Half, class D>
void weighted_store(D d, hn::Vec<D> value, std::uint8_t* dst, std::size_t active) {
  if constexpr (!Half)
    deen_store(d, value, dst, active);
  else {
#if HWY_HAVE_FLOAT16 || (HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX2 && !defined(HWY_DISABLE_F16C))
    const hn::Rebind<hwy::float16_t, D> dh;
    deen_store(dh, hn::DemoteTo(dh, value), dst, active);
#else
    HWY_ALIGN std::uint16_t values[hn::MaxLanes(d)]{};
    store_f16(d, value, values, active);
    std::memcpy(dst, values, active * sizeof(std::uint16_t));
#endif
  }
}

template <bool Half>
void weighted_float_process(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                            double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  const hn::ScalableTag<float> d;
  const hn::Rebind<std::conditional_t<Half, std::uint16_t, std::uint32_t>, decltype(d)> du;
  constexpr std::size_t bytes = Half ? 2 : 4;
  const auto exponent = hn::Set(du, Half ? 0x7c00u : 0x7f800000u);
  const int lanes = static_cast<int>(hn::Lanes(d));
  const auto& p = frames[0];
  // Check all visible inputs, even taps with zero weight. Integer exponent
  // checks remain valid under fast floating-point compilation.
  for (int f = 0; f < count; ++f)
    for (int y = 0; y < p.height; ++y)
      for (int x = 0; x < p.width;) {
        const int active = std::min(lanes, p.width - x);
        const auto bits = deen_raw_load(du, frames[f].data + y * frames[f].stride + x * bytes, active);
        if (!hn::AllFalse(du, hn::Eq(hn::And(bits, exponent), exponent)))
          throw std::invalid_argument("Deen: non-finite sample.");
        x += active;
      }
  const int side = 2 * radius + 1;
  double denominator = 0;
  for (int i = 0; i < side * side; ++i)
    denominator += weights[i];
  denominator *= count == 3 ? 4 : 1;
  const auto reciprocal = hn::Set(d, static_cast<float>(1 / denominator));
  std::array<std::array<float, 225>, 2> coefficients{};
  for (int i = 0; i < side * side; ++i) {
    coefficients[0][i] = static_cast<float>(weights[i] * (count == 3 ? 2 : 1));
    coefficients[1][i] = static_cast<float>(weights[i]);
  }
  std::array<std::array<const std::uint8_t*, 15>, 3> rows{};
  for (int y = 0; y < p.height; ++y) {
    for (int f = 0; f < count; ++f)
      for (int dy = 0; dy < side; ++dy)
        rows[f][dy] = frames[f].data + std::clamp(static_cast<std::int64_t>(y) + dy - radius, std::int64_t{0},
                                                  static_cast<std::int64_t>(p.height - 1)) *
                                           frames[f].stride;
    const auto process_range = [&](auto interior, int begin, int end) HWY_ATTR {
      for (int x = begin; x < end;) {
        const int active = decltype(interior)::value ? lanes : std::min(lanes, end - x);
        const auto center = weighted_load<Half, false>(d, p.data + y * p.stride, p.width, x, active);
        auto sum0 = hn::Zero(d), sum1 = hn::Zero(d);
        auto low0 = center, low1 = center, high0 = center, high1 = center;
        for (int f = 0; f < count; ++f) {
          const auto threshold = hn::Set(d, static_cast<float>(f ? temporal : spatial));
          for (int dy = 0; dy < side; ++dy) {
            const auto* row = rows[f][dy];
            const auto* coeff = coefficients[f != 0].data() + dy * side;
            const auto accumulate = [&](int dx, auto& sum, auto& low, auto& high) HWY_ATTR {
              if (coeff[dx] == 0)
                return;
              const auto sample = weighted_load<Half, decltype(interior)::value>(
                  d, row, p.width, static_cast<std::int64_t>(x) + dx - radius, active);
              const auto difference = hn::Sub(sample, center);
              const auto pass = hn::Le(hn::Abs(difference), threshold);
              // Accepted differences are <=1, even for large finite F32
              // samples. Mask rejected (possibly overflowing) differences
              // before multiplication; the accumulated residual stays small.
              sum = hn::MulAdd(hn::Set(d, coeff[dx]), hn::IfThenElseZero(pass, difference), sum);
              const auto value = hn::IfThenElse(pass, sample, center);
              low = hn::Min(low, value);
              high = hn::Max(high, value);
            };
            int dx = 0;
            for (; dx + 1 < side; dx += 2) {
              accumulate(dx, sum0, low0, high0);
              accumulate(dx + 1, sum1, low1, high1);
            }
            accumulate(dx, sum0, low0, high0);
          }
        }
        const auto result =
            hn::Clamp(hn::MulAdd(hn::Add(sum0, sum1), reciprocal, center), hn::Min(low0, low1), hn::Max(high0, high1));
        weighted_store<Half>(d, result, dst + y * stride + x * bytes, active);
        x += active;
      }
    };
    const int left = std::min(radius, p.width);
    const int vector_end = left + std::max(0, (p.width - radius - left) / lanes) * lanes;
    process_range(std::false_type{}, 0, left);
    process_range(std::true_type{}, left, vector_end);
    process_range(std::false_type{}, vector_end, p.width);
  }
}

void deen_weighted_float_target(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                                double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  if (frames[0].type == DataType::F16)
    weighted_float_process<true>(frames, count, radius, spatial, temporal, weights, dst, stride);
  else
    weighted_float_process<false>(frames, count, radius, spatial, temporal, weights, dst, stride);
}

void deen_w_target(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, double spatial, double temporal, const double* weights, double* dst) {
  const hn::ScalableTag<double> d;
  const auto lanes = hn::Lanes(d);
  const int side = 2 * radius + 1;
  double denominator = 0;
  for (int i = 0; i < side * side; ++i)
    denominator += weights[i];
  denominator *= count == 3 ? 4 : 1;
  for (int y = 0; y < height; ++y)
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto active = std::min(lanes, static_cast<std::size_t>(width) - x);
      const auto center_index = (static_cast<std::size_t>(y) + radius) * pitch + x + radius;
      const auto c = hn::LoadN(d, src[0] + center_index, active);
      auto sum = hn::Zero(d);
      auto low = c, high = c;
      for (int f = 0; f < count; ++f) {
        const auto threshold = f == 0 ? spatial : temporal;
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const double weight = weights[dy * side + dx] * (count == 3 && f == 0 ? 2 : 1);
            if (weight == 0)
              continue;
            const auto s = hn::LoadN(d, src[f] + (static_cast<std::size_t>(y) + dy) * pitch + x + dx, active);
            const auto pass = deen_selection(d, s, c, threshold);
            sum = hn::MulAdd(hn::Set(d, weight), hn::IfThenElseZero(pass, hn::Sub(s, c)), sum);
            const auto value = hn::IfThenElse(pass, s, c);
            low = hn::Min(low, value);
            high = hn::Max(high, value);
          }
      }
      const auto result = hn::Clamp(hn::Add(c, hn::Div(sum, hn::Set(d, denominator))), low, high);
      hn::StoreN(result, d, dst + static_cast<std::size_t>(y) * width + x, active);
    }
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(deen_weighted_float_target);
void deen_weighted_float_kernel(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                                double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  HWY_DYNAMIC_DISPATCH(deen_weighted_float_target)(frames, count, radius, spatial, temporal, weights, dst, stride);
}
HWY_EXPORT(deen_w_target);
void deen_w_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, double spatial, double temporal, const double* weights, double* dst) {
  HWY_DYNAMIC_DISPATCH(deen_w_target)(src, count, pitch, width, height, radius, spatial, temporal, weights, dst);
}
} // namespace neo_smo
#endif
