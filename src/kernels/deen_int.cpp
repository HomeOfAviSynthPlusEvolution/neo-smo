#include "kernels/deen.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/deen_int.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"
HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;
// Access byte-backed host planes without assuming natural sample alignment.
template <class D>
auto integer_raw_load(D d, const std::uint8_t* ptr, std::size_t active) {
#if HWY_TARGET == HWY_SCALAR
  hn::TFromD<D> value;
  std::memcpy(&value, ptr, sizeof(value));
  return hn::Set(d, value);
#else
  const hn::Repartition<std::uint8_t, D> db;
  return hn::BitCast(d,
                     active == hn::Lanes(d) ? hn::LoadU(db, ptr) : hn::LoadN(db, ptr, active * sizeof(hn::TFromD<D>)));
#endif
}
template <class D>
void integer_store(D d, hn::VFromD<D> value, std::uint8_t* ptr, std::size_t active) {
#if HWY_TARGET == HWY_SCALAR
  const auto lane = hn::GetLane(value);
  std::memcpy(ptr, &lane, sizeof(lane));
#else
  const hn::Repartition<std::uint8_t, D> db;
  hn::StoreN(hn::BitCast(db, value), db, ptr, active * sizeof(hn::TFromD<D>));
#endif
}
// Clamp coordinates only at the actual image edges. Interior taps load original pixels.
template <bool Interior = false, class D>
auto integer_load(D d, const DeenPlane& p, std::int64_t y, std::int64_t x, std::size_t active) {
  using T = hn::TFromD<D>;
  const auto* row =
      p.data +
      static_cast<std::ptrdiff_t>(std::clamp(y, std::int64_t{0}, static_cast<std::int64_t>(p.height - 1))) * p.stride;
  if constexpr (Interior) {
    return integer_raw_load(d, row + static_cast<std::size_t>(x) * sizeof(T), hn::Lanes(d));
  }
  if (x >= 0 && x <= p.width - static_cast<int>(active)) {
    return integer_raw_load(d, row + static_cast<std::size_t>(x) * sizeof(T), active);
  }
  HWY_ALIGN T values[hn::MaxLanes(d)]{};
  for (std::size_t i = 0; i < active; ++i)
    std::memcpy(values + i,
                row + static_cast<std::size_t>(std::clamp(static_cast<std::int64_t>(x) + static_cast<std::int64_t>(i),
                                                          std::int64_t{0}, static_cast<std::int64_t>(p.width - 1))) *
                          sizeof(T),
                sizeof(T));
  return hn::Load(d, values);
}
template <bool Constant, class D>
auto integer_mean(D d, hn::VFromD<D> sum, hn::VFromD<D> count, float reciprocal) {
  const hn::Rebind<float, D> df;
  const hn::Rebind<std::int32_t, D> di;
  // Sums <= 243*65535 are exactly representable in F32. Correct the quotient
  // estimate with integer remainders to preserve exact half-up rounding.
  const auto estimate = Constant ? hn::Mul(hn::ConvertTo(df, sum), hn::Set(df, reciprocal))
                                 : hn::Div(hn::ConvertTo(df, sum), hn::ConvertTo(df, count));
  auto q = hn::BitCast(d, hn::ConvertTo(di, estimate));
  q = hn::Sub(q, hn::IfThenElseZero(hn::Gt(hn::Mul(q, count), sum), hn::Set(d, 1)));
  const auto remainder = hn::Sub(sum, hn::Mul(q, count));
  return hn::Add(q, hn::IfThenElseZero(hn::Ge(hn::Add(remainder, remainder), count), hn::Set(d, 1)));
}
template <class T, DeenFamily Family>
void integer_process(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial, double temporal,
                     const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  constexpr bool weighted = Family == DeenFamily::Weighted;
  using Acc = std::conditional_t<sizeof(T) == 1 && !weighted && HWY_TARGET != HWY_SCALAR, std::uint16_t, std::uint32_t>;
  const hn::ScalableTag<Acc> d;
  const hn::Rebind<T, decltype(d)> dn;
  const hn::Rebind<float, decltype(d)> df;
  const hn::Rebind<std::int32_t, decltype(d)> di;
  const int side = radius * 2 + 1, taps = side * side;
  const auto& p = frames[0];
  const int lanes = static_cast<int>(hn::Lanes(d));
  const float reciprocal = 1.0f / (taps * count);
  std::array<std::array<Acc, 225>, 2> limits{};
  double denominator = 0;
  for (int i = 0; i < taps; ++i) {
    denominator += weights[i];
    for (int f = 0; f < 2; ++f)
      limits[f][i] =
          static_cast<Acc>(std::floor((f ? temporal : spatial) * (Family == DeenFamily::Adaptive ? weights[i] : 1)));
  }
  denominator *= count == 3 ? 4 : 1;
  if constexpr (sizeof(T) == 2) {
    if (p.bits < 16) {
      const auto peak = hn::Set(dn, (1u << p.bits) - 1);
      for (int f = 0; f < count; ++f)
        for (int y = 0; y < p.height; ++y)
          for (int x = 0; x < p.width;) {
            const int active = std::min(lanes, p.width - x);
            if (!hn::AllFalse(dn, hn::Gt(integer_load(dn, frames[f], y, x, active), peak)))
              throw std::invalid_argument("Deen: sample exceeds bit depth.");
            x += active;
          }
    }
  }
  for (int y = 0; y < p.height; ++y)
    for (int x = 0; x < p.width;) {
      const int remaining = x < radius             ? std::min(radius - x, p.width - x)
                            : x < p.width - radius ? p.width - radius - x
                                                   : p.width - x;
      const int active = std::min(lanes, remaining);
      const auto batch = [&](auto interior) HWY_ATTR {
        const auto center = hn::PromoteTo(d, integer_load(dn, p, y, x, active));
        auto sum = hn::Zero(d), accepted = hn::Zero(d);
        // Float operations are instantiated only for weighted (u32 accumulator) paths.
        if constexpr (weighted) {
          auto residual = hn::Zero(df);
          const auto cf = hn::ConvertTo(df, center);
          auto low = center, high = center;
          for (int f = 0; f < count; ++f)
            for (int dy = 0; dy < side; ++dy)
              for (int dx = 0; dx < side; ++dx) {
                const float weight = static_cast<float>(weights[dy * side + dx] * (count == 3 && f == 0 ? 2 : 1));
                if (weight == 0)
                  continue;
                const auto sample = hn::PromoteTo(d, integer_load<decltype(interior)::value>(
                                                         dn, frames[f], static_cast<std::int64_t>(y) + dy - radius,
                                                         static_cast<std::int64_t>(x) + dx - radius, active));
                const auto difference = hn::Sub(hn::Max(sample, center), hn::Min(sample, center));
                const auto pass = hn::Le(difference, hn::Set(d, limits[f == 0 ? 0 : 1][0]));
                const auto value = hn::IfThenElse(pass, sample, center);
                residual = hn::MulAdd(hn::Set(df, weight), hn::Sub(hn::ConvertTo(df, value), cf), residual);
                low = hn::Min(low, value);
                high = hn::Max(high, value);
              }
          auto result = hn::MulAdd(residual, hn::Set(df, static_cast<float>(1 / denominator)), cf);
          result = hn::Clamp(result, hn::ConvertTo(df, low), hn::ConvertTo(df, high));
          const auto q = hn::ConvertTo(di, hn::Add(result, hn::Set(df, 0.5f)));
          integer_store(dn, hn::DemoteTo(dn, q), dst + y * stride + static_cast<std::size_t>(x) * sizeof(T), active);
        } else {
          for (int f = 0; f < count; ++f)
            for (int dy = 0; dy < side; ++dy)
              for (int dx = 0; dx < side; ++dx) {
                const auto sample = hn::PromoteTo(d, integer_load<decltype(interior)::value>(
                                                         dn, frames[f], static_cast<std::int64_t>(y) + dy - radius,
                                                         static_cast<std::int64_t>(x) + dx - radius, active));
                const auto difference = hn::Sub(hn::Max(sample, center), hn::Min(sample, center));
                const auto pass = hn::Le(difference, hn::Set(d, limits[f == 0 ? 0 : 1][dy * side + dx]));
                if constexpr (Family == DeenFamily::Adaptive) {
                  sum = hn::Add(sum, hn::IfThenElseZero(pass, sample));
                  accepted = hn::Add(accepted, hn::IfThenElseZero(pass, hn::Set(d, 1)));
                } else
                  sum = hn::Add(sum, hn::IfThenElse(pass, sample, center));
              }
          if constexpr (Family == DeenFamily::Constant)
            accepted = hn::Set(d, taps * count);
          hn::VFromD<decltype(d)> q;
          if constexpr (sizeof(Acc) == 2) {
            const hn::Half<decltype(d)> half;
            const hn::Repartition<std::uint32_t, decltype(d)> wide;
            const auto lo = integer_mean<Family == DeenFamily::Constant>(
                wide, hn::PromoteTo(wide, hn::LowerHalf(half, sum)), hn::PromoteTo(wide, hn::LowerHalf(half, accepted)),
                reciprocal);
            const auto hi = integer_mean<Family == DeenFamily::Constant>(
                wide, hn::PromoteTo(wide, hn::UpperHalf(half, sum)), hn::PromoteTo(wide, hn::UpperHalf(half, accepted)),
                reciprocal);
            q = hn::Combine(d, hn::DemoteTo(half, hi), hn::DemoteTo(half, lo));
          } else
            q = integer_mean<Family == DeenFamily::Constant>(d, sum, accepted, reciprocal);
          integer_store(dn, hn::DemoteTo(dn, q), dst + y * stride + static_cast<std::size_t>(x) * sizeof(T), active);
        }
      };
      if (x >= radius && active == lanes && x <= p.width - radius - lanes)
        batch(std::true_type{});
      else
        batch(std::false_type{});
      x += active;
    }
}
template <class T>
void integer_dispatch(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                      double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  switch (family) {
    case DeenFamily::Constant:
      return integer_process<T, DeenFamily::Constant>(frames, count, radius, spatial, temporal, weights, dst, stride);
    case DeenFamily::Adaptive:
      return integer_process<T, DeenFamily::Adaptive>(frames, count, radius, spatial, temporal, weights, dst, stride);
    case DeenFamily::Weighted:
      return integer_process<T, DeenFamily::Weighted>(frames, count, radius, spatial, temporal, weights, dst, stride);
  }
}
void deen_integer_target(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius,
                         double spatial, double temporal, const double* weights, std::uint8_t* dst,
                         std::ptrdiff_t stride) {
  if (frames[0].type == DataType::U8)
    integer_dispatch<std::uint8_t>(family, frames, count, radius, spatial, temporal, weights, dst, stride);
  else
    integer_dispatch<std::uint16_t>(family, frames, count, radius, spatial, temporal, weights, dst, stride);
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(deen_integer_target);
void deen_integer_kernel(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius,
                         double spatial, double temporal, const double* weights, std::uint8_t* dst,
                         std::ptrdiff_t stride) {
  HWY_DYNAMIC_DISPATCH(deen_integer_target)(family, frames, count, radius, spatial, temporal, weights, dst, stride);
}
} // namespace neo_smo
#endif
