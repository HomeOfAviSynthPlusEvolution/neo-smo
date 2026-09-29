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
auto integer_row_load(D d, const std::uint8_t* row, int width, std::int64_t x, std::size_t active) {
  using T = hn::TFromD<D>;
  if constexpr (Interior) {
    return integer_raw_load(d, row + static_cast<std::size_t>(x) * sizeof(T), hn::Lanes(d));
  }
  if (x >= 0 && x <= width - static_cast<int>(active)) {
    return integer_raw_load(d, row + static_cast<std::size_t>(x) * sizeof(T), active);
  }
  HWY_ALIGN T values[hn::MaxLanes(d)]{};
  for (std::size_t i = 0; i < active; ++i)
    std::memcpy(values + i,
                row + static_cast<std::size_t>(std::clamp(static_cast<std::int64_t>(x) + static_cast<std::int64_t>(i),
                                                          std::int64_t{0}, static_cast<std::int64_t>(width - 1))) *
                          sizeof(T),
                sizeof(T));
  return hn::Load(d, values);
}
template <bool Interior = false, class D>
auto integer_load(D d, const DeenPlane& p, std::int64_t y, std::int64_t x, std::size_t active) {
  const auto* row =
      p.data +
      static_cast<std::ptrdiff_t>(std::clamp(y, std::int64_t{0}, static_cast<std::int64_t>(p.height - 1))) * p.stride;
  return integer_row_load<Interior>(d, row, p.width, x, active);
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
#if HWY_TARGET != HWY_SCALAR
template <DeenFamily Family>
void integer_byte_process(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                          double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  // Compare a full byte vector before widening into two independent sums.
  // At most 243 taps: counts fit u8 and sums (including rounding) fit u16.
  const hn::ScalableTag<std::uint8_t> db;
#if !HWY_ARCH_X86
  const hn::Half<decltype(db)> half;
#endif
  const hn::Repartition<std::uint16_t, decltype(db)> dw;
  const hn::Half<decltype(dw)> wh;
  const hn::Repartition<std::uint32_t, decltype(dw)> dd;
  const int lanes = static_cast<int>(hn::Lanes(db));
  const int side = 2 * radius + 1, taps = side * side, divisor = taps * count;
  const auto reciprocal = hn::Set(dw, 65535 / divisor + 1);
  const auto& p = frames[0];
  struct Tap {
    const std::uint8_t* row;
    int dx;
    std::uint8_t limit;
  };
  std::array<Tap, 243> neighbors;
  std::array<std::array<std::uint8_t, 225>, 2> limits{};
  for (int f = 0; f < 2; ++f)
    for (int i = 0; i < taps; ++i)
      limits[f][i] = static_cast<std::uint8_t>(
          std::floor((f ? temporal : spatial) * (Family == DeenFamily::Adaptive ? weights[i] : 1)));
  for (int y = 0; y < p.height; ++y) {
    int n = 0;
    for (int f = 0; f < count; ++f)
      for (int dy = 0; dy < side; ++dy) {
        const auto* row = frames[f].data + std::clamp(static_cast<std::int64_t>(y) + dy - radius, std::int64_t{0},
                                                      static_cast<std::int64_t>(p.height - 1)) *
                                               frames[f].stride;
        for (int dx = 0; dx < side; ++dx)
          neighbors[n++] = {row, dx - radius, limits[f != 0][dy * side + dx]};
      }
    for (int x = 0; x < p.width;) {
      const int remaining = x < radius             ? std::min(radius - x, p.width - x)
                            : x < p.width - radius ? p.width - radius - x
                                                   : p.width - x;
      const int active = std::min(lanes, remaining);
      const auto batch = [&](auto interior) HWY_ATTR {
        const auto center = integer_load(db, p, y, x, active);
        auto lo = hn::Zero(dw), hi = hn::Zero(dw);
        auto accepted = hn::Zero(db);
        for (int i = 0; i < n; ++i) {
          const auto& tap = neighbors[i];
          const auto sample = integer_row_load<decltype(interior)::value>(
              db, tap.row, p.width, static_cast<std::int64_t>(x) + tap.dx, active);
          const auto difference = hn::Or(hn::SaturatedSub(sample, center), hn::SaturatedSub(center, sample));
          const auto pass = hn::Eq(hn::SaturatedSub(difference, hn::Set(db, tap.limit)), hn::Zero(db));
          const auto value =
              Family == DeenFamily::Adaptive ? hn::IfThenElseZero(pass, sample) : hn::IfThenElse(pass, sample, center);
#if HWY_ARCH_X86
          // Keep each 128-bit block in unpack order until the final pack.
          lo = hn::Add(lo, hn::BitCast(dw, hn::InterleaveLower(db, value, hn::Zero(db))));
          hi = hn::Add(hi, hn::BitCast(dw, hn::InterleaveUpper(db, value, hn::Zero(db))));
#else
          lo = hn::Add(lo, hn::PromoteTo(dw, hn::LowerHalf(half, value)));
          hi = hn::Add(hi, hn::PromoteTo(dw, hn::UpperHalf(half, value)));
#endif
          if constexpr (Family == DeenFamily::Adaptive)
            accepted = hn::Add(accepted, hn::IfThenElseZero(pass, hn::Set(db, 1)));
        }
        const auto mean = [&](auto sum, auto counts) HWY_ATTR {
          if constexpr (Family == DeenFamily::Constant) {
            const auto numerator = hn::Add(sum, hn::Set(dw, divisor / 2));
            auto q = hn::MulHigh(numerator, reciprocal);
            return hn::Sub(q, hn::IfThenElseZero(hn::Gt(hn::Mul(q, hn::Set(dw, divisor)), numerator), hn::Set(dw, 1)));
          } else {
            const auto qlo = integer_mean<false>(dd, hn::PromoteTo(dd, hn::LowerHalf(wh, sum)),
                                                 hn::PromoteTo(dd, hn::LowerHalf(wh, counts)), 0);
            const auto qhi = integer_mean<false>(dd, hn::PromoteTo(dd, hn::UpperHalf(wh, sum)),
                                                 hn::PromoteTo(dd, hn::UpperHalf(wh, counts)), 0);
            return hn::Combine(dw, hn::DemoteTo(wh, qhi), hn::DemoteTo(wh, qlo));
          }
        };
#if HWY_ARCH_X86
        lo = mean(lo, hn::BitCast(dw, hn::InterleaveLower(db, accepted, hn::Zero(db))));
        hi = mean(hi, hn::BitCast(dw, hn::InterleaveUpper(db, accepted, hn::Zero(db))));
        integer_store(db, hn::ReorderDemote2To(db, lo, hi), dst + y * stride + x, active);
#else
        lo = mean(lo, hn::PromoteTo(dw, hn::LowerHalf(half, accepted)));
        hi = mean(hi, hn::PromoteTo(dw, hn::UpperHalf(half, accepted)));
        integer_store(db, hn::Combine(db, hn::DemoteTo(half, hi), hn::DemoteTo(half, lo)), dst + y * stride + x,
                      active);
#endif
      };
      if (x >= radius && active == lanes && x <= p.width - radius - lanes)
        batch(std::true_type{});
      else
        batch(std::false_type{});
      x += active;
    }
  }
}
// Two float accumulators share one narrow comparison batch. With <=243
// nonnegative weighted U8 taps, F32 error is far below half a sample, so
// rounding cannot leave the accepted integer range; per-tap extrema are unnecessary.
void integer_weighted_byte_process(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                                   double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  const hn::ScalableTag<std::uint32_t> d;
  const hn::Rebind<float, decltype(d)> df;
  const hn::Rebind<std::int32_t, decltype(d)> di;
  const hn::Rebind<std::uint8_t, decltype(d)> dh;
  const hn::Twice<decltype(dh)> db;
  const int lanes = static_cast<int>(hn::Lanes(db));
  const int side = 2 * radius + 1, taps = side * side;
  double denominator = 0;
  for (int i = 0; i < taps; ++i)
    denominator += weights[i];
  const auto reciprocal = hn::Set(df, static_cast<float>(1 / (denominator * (count == 3 ? 4 : 1))));
  const auto& p = frames[0];
  struct Tap {
    const std::uint8_t* row;
    int dx;
    std::uint8_t limit;
    float weight;
  };
  std::array<Tap, 243> neighbors;
  for (int y = 0; y < p.height; ++y) {
    int n = 0;
    for (int f = 0; f < count; ++f)
      for (int dy = 0; dy < side; ++dy) {
        const auto* row = frames[f].data + std::clamp(static_cast<std::int64_t>(y) + dy - radius, std::int64_t{0},
                                                      static_cast<std::int64_t>(p.height - 1)) *
                                               frames[f].stride;
        for (int dx = 0; dx < side; ++dx) {
          const auto weight = static_cast<float>(weights[dy * side + dx] * (count == 3 && f == 0 ? 2 : 1));
          if (weight != 0)
            neighbors[n++] = {row, dx - radius, static_cast<std::uint8_t>(std::floor(f ? temporal : spatial)), weight};
        }
      }
    for (int x = 0; x < p.width;) {
      const int remaining = x < radius             ? std::min(radius - x, p.width - x)
                            : x < p.width - radius ? p.width - radius - x
                                                   : p.width - x;
      const int active = std::min(lanes, remaining);
      const auto batch = [&](auto interior) HWY_ATTR {
        const auto center = integer_load(db, p, y, x, active);
        const auto cl = hn::ConvertTo(df, hn::PromoteTo(d, hn::LowerHalf(dh, center)));
        const auto ch = hn::ConvertTo(df, hn::PromoteTo(d, hn::UpperHalf(dh, center)));
        auto lo = hn::Zero(df), hi = hn::Zero(df);
        for (int i = 0; i < n; ++i) {
          const auto& tap = neighbors[i];
          const auto sample = integer_row_load<decltype(interior)::value>(
              db, tap.row, p.width, static_cast<std::int64_t>(x) + tap.dx, active);
          const auto diff = hn::Or(hn::SaturatedSub(sample, center), hn::SaturatedSub(center, sample));
          const auto pass = hn::Eq(hn::SaturatedSub(diff, hn::Set(db, tap.limit)), hn::Zero(db));
          const auto value = hn::IfThenElse(pass, sample, center);
          const auto weight = hn::Set(df, tap.weight);
          lo = hn::MulAdd(weight, hn::Sub(hn::ConvertTo(df, hn::PromoteTo(d, hn::LowerHalf(dh, value))), cl), lo);
          hi = hn::MulAdd(weight, hn::Sub(hn::ConvertTo(df, hn::PromoteTo(d, hn::UpperHalf(dh, value))), ch), hi);
        }
        const auto ql = hn::ConvertTo(di, hn::Add(hn::MulAdd(lo, reciprocal, cl), hn::Set(df, 0.5f)));
        const auto qh = hn::ConvertTo(di, hn::Add(hn::MulAdd(hi, reciprocal, ch), hn::Set(df, 0.5f)));
        integer_store(db, hn::Combine(db, hn::DemoteTo(dh, qh), hn::DemoteTo(dh, ql)), dst + y * stride + x, active);
      };
      if (x >= radius && active == lanes && x <= p.width - radius - lanes)
        batch(std::true_type{});
      else
        batch(std::false_type{});
      x += active;
    }
  }
}

#endif
template <class T, DeenFamily Family>
void integer_process(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial, double temporal,
                     const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
#if HWY_TARGET != HWY_SCALAR
  if constexpr (sizeof(T) == 1 && Family != DeenFamily::Weighted)
    return integer_byte_process<Family>(frames, count, radius, spatial, temporal, weights, dst, stride);
#endif
#if HWY_TARGET != HWY_SCALAR
  if constexpr (sizeof(T) == 1 && Family == DeenFamily::Weighted)
    return integer_weighted_byte_process(frames, count, radius, spatial, temporal, weights, dst, stride);
#endif
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
  struct Tap {
    const std::uint8_t* row;
    int dx;
    Acc limit;
    float weight;
  };
  std::array<Tap, 243> neighbors;
  for (int y = 0; y < p.height; ++y) {
    int num_neighbors = 0;
    if constexpr (weighted) {
      for (int f = 0; f < count; ++f)
        for (int dy = 0; dy < side; ++dy) {
          const auto row_y = std::clamp(static_cast<std::int64_t>(y) + dy - radius, std::int64_t{0},
                                        static_cast<std::int64_t>(p.height - 1));
          const auto* row = frames[f].data + row_y * frames[f].stride;
          for (int dx = 0; dx < side; ++dx) {
            const float weight = static_cast<float>(weights[dy * side + dx] * (count == 3 && f == 0 ? 2 : 1));
            if constexpr (weighted) {
              if (weight == 0)
                continue;
            }
            neighbors[num_neighbors++] = {row, dx - radius, limits[f == 0 ? 0 : 1][dy * side + dx], weight};
          }
        }
    }
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
          for (int i = 0; i < num_neighbors; ++i) {
            const auto& tap = neighbors[i];
            const auto sample =
                hn::PromoteTo(d, integer_row_load<decltype(interior)::value>(
                                     dn, tap.row, p.width, static_cast<std::int64_t>(x) + tap.dx, active));
            const auto difference = hn::Sub(hn::Max(sample, center), hn::Min(sample, center));
            const auto pass = hn::Le(difference, hn::Set(d, tap.limit));
            const auto value = hn::IfThenElse(pass, sample, center);
            residual = hn::MulAdd(hn::Set(df, tap.weight), hn::Sub(hn::ConvertTo(df, value), cf), residual);
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
          if constexpr (sizeof(Acc) == 2 && Family == DeenFamily::Constant) {
            // The largest byte sum plus half the divisor fits u16. A ceiling
            // reciprocal overestimates division by at most one; correct exactly.
            const auto divisor = hn::Set(d, taps * count);
            const auto numerator = hn::Add(sum, hn::Set(d, (taps * count) / 2));
            q = hn::MulHigh(numerator, hn::Set(d, 65535 / (taps * count) + 1));
            q = hn::Sub(q, hn::IfThenElseZero(hn::Gt(hn::Mul(q, divisor), numerator), hn::Set(d, 1)));
          } else if constexpr (sizeof(Acc) == 2) {
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
