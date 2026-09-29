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
#include "kernels/deen_io-inl.hpp"
template <bool Interior = false, class D>
auto integer_load(D d, const DeenPlane& p, std::int64_t y, std::int64_t x, std::size_t active) {
  const auto* row =
      p.data +
      static_cast<std::ptrdiff_t>(std::clamp(y, std::int64_t{0}, static_cast<std::int64_t>(p.height - 1))) * p.stride;
  return deen_row_load<Interior>(d, row, p.width, x, active);
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
  // Exact division over the bounded U8 numerator range. Pick a rounded-down
  // reciprocal plus one in the numerator, or a rounded-up reciprocal. The
  // error bound at quotient 255 proves all smaller quotients exact as well.
  int shift = 0;
  while ((1 << (shift + 1)) <= divisor)
    ++shift;
  const unsigned scale = 1u << (16 + shift);
  unsigned multiplier = scale / divisor;
  unsigned rounding = divisor / 2;
  if (255u * (scale - multiplier * divisor) <= multiplier)
    ++rounding;
  else
    ++multiplier;
  const auto reciprocal = hn::Set(dw, multiplier);
  const auto rounding_bias = hn::Set(dw, rounding);
  const auto& p = frames[0];
  struct Row {
    const std::uint8_t* data;
    const std::uint8_t* limits;
  };
  std::array<Row, 27> rows;
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
        rows[n++] = {row, limits[f != 0].data() + dy * side};
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
          const auto& row = rows[i];
          const auto row_limit = hn::Set(db, row.limits[0]);
          const auto lower = hn::SaturatedSub(center, row_limit), upper = hn::SaturatedAdd(center, row_limit);
          for (int dx = 0; dx < side; ++dx) {
            const auto sample = deen_row_load<decltype(interior)::value>(
                db, row.data, p.width, static_cast<std::int64_t>(x) + dx - radius, active);
            const auto pass = [&]() HWY_ATTR {
              if constexpr (Family == DeenFamily::Constant)
                return hn::Eq(hn::Or(hn::SaturatedSub(sample, upper), hn::SaturatedSub(lower, sample)), hn::Zero(db));
              else {
                const auto difference = hn::Or(hn::SaturatedSub(sample, center), hn::SaturatedSub(center, sample));
                return hn::Eq(hn::SaturatedSub(difference, hn::Set(db, row.limits[dx])), hn::Zero(db));
              }
            }();
            const auto value = Family == DeenFamily::Adaptive ? hn::IfThenElseZero(pass, sample)
                                                              : hn::IfThenElse(pass, sample, center);
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
        }
        const auto mean = [&](auto sum, auto counts) HWY_ATTR {
          if constexpr (Family == DeenFamily::Constant) {
            return hn::ShiftRightSame(hn::MulHigh(hn::Add(sum, rounding_bias), reciprocal), shift);
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
        deen_store(db, hn::ReorderDemote2To(db, lo, hi), dst + y * stride + x, active);
#else
        lo = mean(lo, hn::PromoteTo(dw, hn::LowerHalf(half, accepted)));
        hi = mean(hi, hn::PromoteTo(dw, hn::UpperHalf(half, accepted)));
        deen_store(db, hn::Combine(db, hn::DemoteTo(half, hi), hn::DemoteTo(half, lo)), dst + y * stride + x, active);
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
#if HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX3
// Equal-weight taps can be summed exactly in u16 before conversion to F32.
// A spatial window has at most 225 U8 samples (sum <=57375). Each frame is
// accumulated separately; the current frame's temporal factor is applied in F32.
void integer_weighted_byte_grouped(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                                   double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  const hn::ScalableTag<std::uint8_t> db;
  const hn::Repartition<std::uint16_t, decltype(db)> dw;
  const hn::Repartition<std::uint32_t, decltype(db)> d;
  const hn::Rebind<float, decltype(d)> df;
  const hn::Rebind<std::int32_t, decltype(d)> di;
  const int lanes = static_cast<int>(hn::Lanes(db));
  const int side = 2 * radius + 1, taps = side * side;
  std::array<int, 225> order;
  double denominator = 0;
  for (int i = 0; i < taps; ++i) {
    order[i] = i;
    denominator += weights[i];
  }
  std::sort(order.begin(), order.begin() + taps, [&](int a, int b) { return weights[a] < weights[b]; });
  struct Group {
    int begin, end;
    std::array<float, 2> weight;
  };
  std::array<Group, 225> groups;
  int num_groups = 0;
  for (int begin = 0; begin < taps;) {
    int end = begin + 1;
    while (end < taps && weights[order[end]] == weights[order[begin]])
      ++end;
    if (weights[order[begin]] != 0)
      groups[num_groups++] = {begin,
                              end,
                              {static_cast<float>(weights[order[begin]] * (count == 3 ? 2 : 1)),
                               static_cast<float>(weights[order[begin]])}};
    begin = end;
  }
  const auto reciprocal = hn::Set(df, static_cast<float>(1 / (denominator * (count == 3 ? 4 : 1))));
  const std::array<std::uint8_t, 2> limits{static_cast<std::uint8_t>(std::floor(spatial)),
                                           static_cast<std::uint8_t>(std::floor(temporal))};
  const auto& p = frames[0];
  std::array<std::array<std::ptrdiff_t, 225>, 3> offsets;
  std::array<int, 225> dxs;
  for (int i = 0; i < taps; ++i)
    dxs[i] = order[i] % side - radius;
  for (int y = 0; y < p.height; ++y) {
    for (int f = 0; f < count; ++f)
      for (int i = 0; i < taps; ++i) {
        const int dy = order[i] / side - radius;
        offsets[f][i] =
            std::clamp(static_cast<std::int64_t>(y) + dy, std::int64_t{0}, static_cast<std::int64_t>(p.height - 1)) *
                frames[f].stride +
            dxs[i];
      }
    for (int x = 0; x < p.width;) {
      const int remaining = x < radius             ? std::min(radius - x, p.width - x)
                            : x < p.width - radius ? p.width - radius - x
                                                   : p.width - x;
      const int active = std::min(lanes, remaining);
      const auto batch = [&](auto interior) HWY_ATTR {
        const auto center = integer_load(db, p, y, x, active);
        auto s0 = hn::Zero(df), s1 = hn::Zero(df), s2 = hn::Zero(df), s3 = hn::Zero(df);
        for (int f = 0; f < count; ++f) {
          const auto limit = hn::Set(db, limits[f != 0]);
          const auto lower = hn::SaturatedSub(center, limit), upper = hn::SaturatedAdd(center, limit);
          for (int g = 0; g < num_groups; ++g) {
            const auto& group = groups[g];
            auto lo = hn::Zero(dw), hi = hn::Zero(dw);
            const auto select = [&](int i) HWY_ATTR {
              const auto sample = [&]() HWY_ATTR {
                if constexpr (decltype(interior)::value)
                  return deen_raw_load(db, frames[f].data + (offsets[f][i] + x), lanes);
                else
                  return deen_row_load(db, frames[f].data + (offsets[f][i] - dxs[i]), p.width,
                                       static_cast<std::int64_t>(x) + dxs[i], active);
              }();
              const auto outside = hn::Or(hn::SaturatedSub(sample, upper), hn::SaturatedSub(lower, sample));
              const auto pass = hn::Eq(outside, hn::Zero(db));
              return hn::IfThenElse(pass, sample, center);
            };
            int i = group.begin;
            const hn::Rebind<std::int16_t, decltype(dw)> ds;
            const hn::Rebind<std::int8_t, decltype(db)> dsb;
            const auto ones = hn::Set(dsb, 1);
            for (; i + 1 < group.end; i += 2) {
              const auto a = select(i), b = select(i + 1);
              lo = hn::Add(lo, hn::BitCast(dw, hn::SatWidenMulPairwiseAdd(ds, hn::InterleaveLower(db, a, b), ones)));
              hi = hn::Add(hi, hn::BitCast(dw, hn::SatWidenMulPairwiseAdd(ds, hn::InterleaveUpper(db, a, b), ones)));
            }
            for (; i < group.end; ++i) {
              const auto value = select(i);
              lo = hn::Add(lo, hn::BitCast(dw, hn::InterleaveLower(db, value, hn::Zero(db))));
              hi = hn::Add(hi, hn::BitCast(dw, hn::InterleaveUpper(db, value, hn::Zero(db))));
            }
            const auto weight = hn::Set(df, group.weight[f != 0]);
            s0 = hn::MulAdd(weight, hn::ConvertTo(df, hn::BitCast(d, hn::InterleaveLower(dw, lo, hn::Zero(dw)))), s0);
            s1 = hn::MulAdd(weight, hn::ConvertTo(df, hn::BitCast(d, hn::InterleaveUpper(dw, lo, hn::Zero(dw)))), s1);
            s2 = hn::MulAdd(weight, hn::ConvertTo(df, hn::BitCast(d, hn::InterleaveLower(dw, hi, hn::Zero(dw)))), s2);
            s3 = hn::MulAdd(weight, hn::ConvertTo(df, hn::BitCast(d, hn::InterleaveUpper(dw, hi, hn::Zero(dw)))), s3);
          }
        }
        const auto finish = [&](auto sum) HWY_ATTR {
          return hn::ConvertTo(di, hn::MulAdd(sum, reciprocal, hn::Set(df, 0.5f)));
        };
        const auto lo = hn::ReorderDemote2To(dw, finish(s0), finish(s1));
        const auto hi = hn::ReorderDemote2To(dw, finish(s2), finish(s3));
        deen_store(db, hn::ReorderDemote2To(db, lo, hi), dst + y * stride + x, active);
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
// U8 uses normalized Q15 coefficients stored doubled for unsigned MulHigh.
// Unpacking bytes into the high halves folds in Q8 scaling. The cached bias
// centers the products' truncation interval before the final output rounding.
void integer_weighted_byte_process(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                                   double temporal, const double* weights, const DeenByteWeights& byte_weights,
                                   std::uint8_t* dst, std::ptrdiff_t stride) {
#if HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX3
  // Wide targets favor fewer multiplies for larger distance windows.
  if (radius >= 3)
    return integer_weighted_byte_grouped(frames, count, radius, spatial, temporal, weights, dst, stride);
#endif
  const hn::ScalableTag<std::uint8_t> db;
#if !HWY_ARCH_X86
  const hn::Half<decltype(db)> half;
#endif
  const hn::Repartition<std::uint16_t, decltype(db)> dw;
  const int lanes = static_cast<int>(hn::Lanes(db));
  const int side = 2 * radius + 1;
  // A uniform spatial window is already an exact constant mean.
  if (count == 1 && weights[0] == 1)
    return integer_byte_process<DeenFamily::Constant>(frames, count, radius, spatial, temporal, weights, dst, stride);
  const std::array<std::uint8_t, 2> limits{static_cast<std::uint8_t>(std::floor(spatial)),
                                           static_cast<std::uint8_t>(std::floor(temporal))};
  const auto& p = frames[0];
  std::array<std::array<const std::uint8_t*, 15>, 3> rows;
  for (int y = 0; y < p.height; ++y) {
    for (int f = 0; f < count; ++f)
      for (int dy = 0; dy < side; ++dy)
        rows[f][dy] = frames[f].data + std::clamp(static_cast<std::int64_t>(y) + dy - radius, std::int64_t{0},
                                                  static_cast<std::int64_t>(p.height - 1)) *
                                           frames[f].stride;
    const auto process_range = [&](auto interior, int begin, int end) HWY_ATTR {
      for (int x = begin; x < end;) {
        const int active = decltype(interior)::value ? lanes : std::min(lanes, end - x);
        const auto center = deen_raw_load(db, p.data + y * p.stride + x, active);
        auto lo = hn::Set(dw, byte_weights.rounding), hi = hn::Set(dw, byte_weights.rounding);
        for (int f = 0; f < count; ++f) {
          const auto limit = hn::Set(db, limits[f != 0]);
          const auto lower = hn::SaturatedSub(center, limit), upper = hn::SaturatedAdd(center, limit);
          for (int dy = 0; dy < side; ++dy) {
            const auto* row = rows[f][dy];
            const auto* coeff = byte_weights.coefficients.data() + (f != 0 ? 225 : 0) + dy * side;
            for (int dx = 0; dx < side; ++dx) {
              const auto sample = deen_row_load<decltype(interior)::value>(
                  db, row, p.width, static_cast<std::int64_t>(x) + dx - radius, active);
              const auto outside = hn::Or(hn::SaturatedSub(sample, upper), hn::SaturatedSub(lower, sample));
              const auto value = hn::IfThenElse(hn::Eq(outside, hn::Zero(db)), sample, center);
#if HWY_ARCH_X86
              const auto vl = hn::BitCast(dw, hn::InterleaveLower(db, hn::Zero(db), value));
              const auto vh = hn::BitCast(dw, hn::InterleaveUpper(db, hn::Zero(db), value));
#else
              const auto vl = hn::ShiftLeft<8>(hn::PromoteTo(dw, hn::LowerHalf(half, value)));
              const auto vh = hn::ShiftLeft<8>(hn::PromoteTo(dw, hn::UpperHalf(half, value)));
#endif
              const auto weight = hn::Set(dw, coeff[dx]);
              lo = hn::Add(lo, hn::MulHigh(vl, weight));
              hi = hn::Add(hi, hn::MulHigh(vh, weight));
            }
          }
        }
        // The coefficient sum is 65536 and the centered product error is
        // below 122 Q8 units. Including output rounding, the total stays
        // within U16 and less than half a sample outside the weighted mean.
        lo = hn::ShiftRight<8>(lo);
        hi = hn::ShiftRight<8>(hi);
#if HWY_ARCH_X86
        deen_store(db, hn::ReorderDemote2To(db, lo, hi), dst + y * stride + x, active);
#else
        deen_store(db, hn::Combine(db, hn::DemoteTo(half, hi), hn::DemoteTo(half, lo)), dst + y * stride + x, active);
#endif
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

#endif
template <class T, DeenFamily Family>
void integer_process(const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial, double temporal,
                     const double* weights, const DeenByteWeights& byte_weights, std::uint8_t* dst,
                     std::ptrdiff_t stride) {
#if HWY_TARGET != HWY_SCALAR
  if constexpr (sizeof(T) == 1 && Family != DeenFamily::Weighted)
    return integer_byte_process<Family>(frames, count, radius, spatial, temporal, weights, dst, stride);
#endif
#if HWY_TARGET != HWY_SCALAR
  if constexpr (sizeof(T) == 1 && Family == DeenFamily::Weighted)
    return integer_weighted_byte_process(frames, count, radius, spatial, temporal, weights, byte_weights, dst, stride);
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
                hn::PromoteTo(d, deen_row_load<decltype(interior)::value>(
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
          deen_store(dn, hn::DemoteTo(dn, q), dst + y * stride + static_cast<std::size_t>(x) * sizeof(T), active);
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
          deen_store(dn, hn::DemoteTo(dn, q), dst + y * stride + static_cast<std::size_t>(x) * sizeof(T), active);
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
                      double temporal, const double* weights, const DeenByteWeights& byte_weights, std::uint8_t* dst,
                      std::ptrdiff_t stride) {
  switch (family) {
    case DeenFamily::Constant:
      return integer_process<T, DeenFamily::Constant>(frames, count, radius, spatial, temporal, weights, byte_weights,
                                                      dst, stride);
    case DeenFamily::Adaptive:
      return integer_process<T, DeenFamily::Adaptive>(frames, count, radius, spatial, temporal, weights, byte_weights,
                                                      dst, stride);
    case DeenFamily::Weighted:
      return integer_process<T, DeenFamily::Weighted>(frames, count, radius, spatial, temporal, weights, byte_weights,
                                                      dst, stride);
  }
}
void deen_integer_target(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius,
                         double spatial, double temporal, const double* weights, const DeenByteWeights& byte_weights,
                         std::uint8_t* dst, std::ptrdiff_t stride) {
  if (frames[0].type == DataType::U8)
    integer_dispatch<std::uint8_t>(family, frames, count, radius, spatial, temporal, weights, byte_weights, dst,
                                   stride);
  else
    integer_dispatch<std::uint16_t>(family, frames, count, radius, spatial, temporal, weights, byte_weights, dst,
                                    stride);
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(deen_integer_target);
void deen_integer_kernel(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius,
                         double spatial, double temporal, const double* weights, const DeenByteWeights& byte_weights,
                         std::uint8_t* dst, std::ptrdiff_t stride) {
  HWY_DYNAMIC_DISPATCH(deen_integer_target)(family, frames, count, radius, spatial, temporal, weights, byte_weights,
                                            dst, stride);
}
} // namespace neo_smo
#endif
