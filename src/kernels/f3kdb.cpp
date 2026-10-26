#include "kernels/f3kdb.hpp"
#include "kernels/f3kdb_scalar.hpp"
#include "kernels/f3kdb_weight_table.hpp"
#include <algorithm>
#include <cstring>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#ifndef HWY_COMPILE_ALL_ATTAINABLE
#define HWY_COMPILE_ALL_ATTAINABLE
#endif
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/f3kdb.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"
HWY_BEFORE_NAMESPACE();
namespace neo_smo::f3kdb {
namespace HWY_NAMESPACE {
#if HWY_TARGET == HWY_SCALAR || HWY_TARGET == HWY_EMU128
std::size_t collect_cross_target(const CrossSamplingRow&) {
  return 0;
}
SamplingResult cross_sampled_target(const CrossSamplingRow&, const Row&, const double*, bool, double, double) {
  return {0, 0};
}
std::size_t collect_simple_target(const SamplingRow&, int) {
  return 0;
}
SamplingResult sampled_target(const SamplingRow&, const Row&, int, double, bool, bool) {
  return {0, 0};
}
std::size_t byte_target(const ByteRow&, int, const double*, bool) {
  return 0;
}
std::size_t integer_target(const IntegerRow&, int, const double*, bool) {
  return 0;
}
void row_target(const Row& row, int mode, const double* t, bool integer, bool blur, double angle, double boost) {
  process_row_scalar(row, mode, t, integer, blur, angle, boost);
}
void weights_target(const double* c[4], double* output, std::size_t count) {
  weights_scalar(c, output, count);
}
#else
namespace hn = hwy::HWY_NAMESPACE;
using D = hn::ScalableTag<double>;
using V = hn::VFromD<D>;
using M = hn::MFromD<D>;
HWY_INLINE V load_integer_grain(D d, const std::int16_t* p, std::size_t n) {
  const hn::Rebind<std::int16_t, D> di16;
  const hn::Rebind<std::int32_t, D> di32;
  return hn::PromoteTo(d, hn::PromoteTo(di32, n == hn::Lanes(d) ? hn::LoadU(di16, p) : hn::LoadN(di16, p, n)));
}
HWY_INLINE V load_grain(D d, const float* p, std::size_t n) {
  const hn::Rebind<float, D> df;
  return hn::PromoteTo(d, n == hn::Lanes(d) ? hn::LoadU(df, p) : hn::LoadN(df, p, n));
}
template <int N, bool Blur, bool Integer = true>
HWY_INLINE V sampled_value(D d, V center, const V* refs, double threshold);
HWY_INLINE V fp(D d, V v);
// Read the interleaved distance bytes directly. Chroma subsampling by two
// skips every other pair; wider subsampling uses the existing scalar collector.
template <int Mode, bool Fused = false, bool Blur = true, bool Integer = true>
std::size_t collect_simple(const SamplingRow& row, const Row* output = nullptr, double threshold = 0) {
  const D d;
  const hn::Rebind<std::int64_t, D> di;
  const hn::Rebind<std::uint8_t, D> db;
  static_assert(sizeof(Distance) == 2);
  const auto lanes = hn::Lanes(d);
  const auto full = row.width / lanes * lanes;
  const auto step = std::size_t{1} << row.sub_x;
  const auto* distances = reinterpret_cast<const std::uint8_t*>(row.distances);
  const auto index = hn::Iota(di, 0);
  const auto vertical_limit = hn::Set(di, static_cast<std::int64_t>(row.vertical_limit));
  const auto last = hn::Set(di, static_cast<std::int64_t>(row.width - 1));
  for (std::size_t x = 0; x < full; x += lanes, distances += 2 * step * lanes) {
    hn::VFromD<decltype(db)> first, second;
    if (row.sub_x == 0) {
      hn::LoadInterleaved2(db, distances, first, second);
    } else {
      hn::VFromD<decltype(db)> skip_first, skip_second;
      hn::LoadInterleaved4(db, distances, first, second, skip_first, skip_second);
    }
    const auto distance1 = hn::PromoteTo(di, first);
    const auto h1 = hn::ShiftRightSame(distance1, row.sub_x);
    const auto v1 = hn::ShiftRightSame(distance1, row.sub_y);
    const auto position = hn::Add(index, hn::Set(di, static_cast<std::int64_t>(x)));
    const auto horizontal_limit = hn::Min(position, hn::Sub(last, position));
    auto valid = hn::Eq(index, index);
    if constexpr (Mode != 3)
      valid = hn::And(valid, hn::Le(v1, vertical_limit));
    if constexpr (Mode != 1)
      valid = hn::And(valid, hn::Le(h1, horizontal_limit));
    if (!hn::AllTrue(di, valid))
      throw std::logic_error("reference outside plane");
    const auto* center = row.center + x;
    V refs[Mode == 2 ? 4 : 2];
    if constexpr (Mode == 1 || Mode == 3) {
      const auto offset = Mode == 1 ? hn::GatherIndex(di, row.vertical_offsets, v1) : h1;
      refs[0] = hn::GatherIndex(d, center, hn::Add(index, offset));
      if constexpr (!Fused)
        hn::StoreU(refs[0], d, row.samples[0] + x);
      refs[1] = hn::GatherIndex(d, center, hn::Sub(index, offset));
      if constexpr (!Fused)
        hn::StoreU(refs[1], d, row.samples[1] + x);
    } else {
      const auto distance2 = hn::PromoteTo(di, second);
      const auto h2 = hn::ShiftRightSame(distance2, row.sub_x);
      const auto v2 = hn::ShiftRightSame(distance2, row.sub_y);
      valid = hn::And(hn::Le(v2, vertical_limit), hn::Le(h2, horizontal_limit));
      if (!hn::AllTrue(di, valid))
        throw std::logic_error("reference outside plane");
      const auto a = hn::Add(hn::GatherIndex(di, row.vertical_offsets, v2), h1);
      const auto b = hn::Sub(h2, hn::GatherIndex(di, row.vertical_offsets, v1));
      refs[0] = hn::GatherIndex(d, center, hn::Add(index, a));
      if constexpr (!Fused)
        hn::StoreU(refs[0], d, row.samples[0] + x);
      refs[1] = hn::GatherIndex(d, center, hn::Add(index, b));
      if constexpr (!Fused)
        hn::StoreU(refs[1], d, row.samples[1] + x);
      refs[2] = hn::GatherIndex(d, center, hn::Sub(index, a));
      if constexpr (!Fused)
        hn::StoreU(refs[2], d, row.samples[2] + x);
      refs[3] = hn::GatherIndex(d, center, hn::Sub(index, b));
      if constexpr (!Fused)
        hn::StoreU(refs[3], d, row.samples[3] + x);
    }
    if constexpr (Fused) {
      auto result = sampled_value<Mode == 2 ? 4 : 2, Blur, Integer>(d, hn::LoadU(d, center), refs, threshold);
      if constexpr (Integer)
        result = hn::Add(result, load_integer_grain(d, output->integer_grain + x, hn::Lanes(d)));
      else
        result = fp(d, hn::Add(result, load_grain(d, output->grain + x, hn::Lanes(d))));
      hn::StoreU(result, d, output->output + x);
    }
  }
  return full;
}
std::size_t collect_simple_target(const SamplingRow& row, int mode) {
  // Two-lane gathers cost more than the scalar collector on narrow targets.
  if (hn::Lanes(D{}) < 4 || row.sub_x > 1 || row.width == 0)
    return 0;
  switch (mode) {
    case 1:
      return collect_simple<1>(row);
    case 2:
      return collect_simple<2>(row);
    case 3:
      return collect_simple<3>(row);
    default:
      return 0;
  }
}
// Full vectors use ordinary unaligned accesses. Only the last partial vector
// needs bounded loads/stores (notably important on pre-AVX-512 targets).
HWY_INLINE V load_row(D d, const double* p, std::size_t n) {
  return n == hn::Lanes(d) ? hn::LoadU(d, p) : hn::LoadN(d, p, n);
}
HWY_INLINE void store_row(D d, V v, double* p, std::size_t n) {
  if (n == hn::Lanes(d))
    hn::StoreU(v, d, p);
  else
    hn::StoreN(v, d, p, n);
}
// F32 operations with a widened overflow fallback. Widening before arithmetic
// also preserves subnormals in weighted means and tiny gradient products.
HWY_INLINE V fp(D d, V v) {
  const hn::Rebind<float, D> df;
  const auto finite = hn::Le(hn::Abs(v), hn::Set(d, double(std::numeric_limits<float>::max())));
#if HWY_TARGET == HWY_AVX2
  // AVX2 benefits from avoiding per-lane overflow selection for full blocks.
  if (hn::AllTrue(d, finite))
    return hn::PromoteTo(d, hn::DemoteTo(df, v));
#endif
  const auto safe = hn::IfThenElseZero(finite, v);
  return hn::IfThenElse(finite, hn::PromoteTo(d, hn::DemoteTo(df, safe)), v);
}
HWY_INLINE V diff(D d, V a, V b, bool integer) {
  const auto v = hn::Sub(a, b);
  return integer ? v : fp(d, v);
}
HWY_INLINE V mean(D d, const V* r, int n, bool integer) {
  auto sum = hn::Zero(d);
  const auto inverse = hn::Set(d, 1.0 / n);
  for (int j = 0; j < n; ++j)
    sum = hn::Add(sum, integer ? r[j] : hn::Mul(r[j], inverse));
  return integer ? hn::Floor(hn::Mul(hn::Add(sum, hn::Set(d, n / 2.0)), inverse)) : fp(d, sum);
}
template <int N, bool Blur, bool Integer>
HWY_INLINE V sampled_value(D d, V center, const V* refs, double threshold) {
  const auto average = mean(d, refs, N, Integer);
  const auto limit = hn::Set(d, threshold);
  auto accept = hn::Lt(hn::Abs(diff(d, average, center, Integer)), limit);
  if constexpr (!Blur) {
    accept = hn::Eq(center, center);
    for (int j = 0; j < N; ++j)
      accept = hn::And(accept, hn::Lt(hn::Abs(diff(d, refs[j], center, Integer)), limit));
  }
  return hn::IfThenElse(accept, average, center);
}
SamplingResult sampled_target(const SamplingRow& sampling, const Row& row, int mode, double threshold, bool integer,
                              bool blur) {
  // Four-reference fusion is slower than collection on four-double-lane targets.
  if (hn::Lanes(D{}) < 4 || sampling.sub_x > 1 || sampling.width == 0 || (hn::Lanes(D{}) == 4 && mode == 2))
    return {collect_simple_target(sampling, mode), 0};
  if (!integer) {
    constexpr bool fused_float_sampling = HWY_TARGET == HWY_AVX3_SPR;
    if (!fused_float_sampling)
      return {collect_simple_target(sampling, mode), 0};
    std::size_t n = 0;
    if (mode == 1)
      n = blur ? collect_simple<1, true, true, false>(sampling, &row, threshold)
               : collect_simple<1, true, false, false>(sampling, &row, threshold);
    else if (mode == 3)
      n = blur ? collect_simple<3, true, true, false>(sampling, &row, threshold)
               : collect_simple<3, true, false, false>(sampling, &row, threshold);
    else if (mode == 2)
      n = blur ? collect_simple<2, true, true, false>(sampling, &row, threshold)
               : collect_simple<2, true, false, false>(sampling, &row, threshold);
    else
      return {collect_simple_target(sampling, mode), 0};
    return {n, n};
  }
  std::size_t count = 0;
  switch (mode) {
    case 1:
      count = blur ? collect_simple<1, true, true>(sampling, &row, threshold)
                   : collect_simple<1, true, false>(sampling, &row, threshold);
      break;
    case 2:
      count = blur ? collect_simple<2, true, true>(sampling, &row, threshold)
                   : collect_simple<2, true, false>(sampling, &row, threshold);
      break;
    case 3:
      count = blur ? collect_simple<3, true, true>(sampling, &row, threshold)
                   : collect_simple<3, true, false>(sampling, &row, threshold);
      break;
    default:
      break;
  }
  return {count, count};
}
HWY_INLINE V weight(D d, const V* c);
HWY_INLINE V response(D d, V delta, V threshold);
HWY_INLINE V response_with_inverse(D d, V delta, V threshold, V inverse);
HWY_INLINE M aligned(D d, V gx, V gy, V hx, V hy, double angle, V cosine, V sine);
struct DirectionCenter {
  V x, y, magnitude, normalized_x, normalized_y;
  M integral;
};
HWY_INLINE DirectionCenter prepare_direction(D d, V x, V y) {
  const auto zero = hn::Zero(d), one = hn::Set(d, 1.0);
  const auto magnitude = hn::Max(hn::Abs(x), hn::Abs(y));
  const auto integral = hn::And(hn::Le(magnitude, hn::Set(d, 262140.0)),
                                hn::And(hn::Eq(hn::Trunc(x), x), hn::Eq(hn::Trunc(y), y)));
  const auto divisor = hn::IfThenElse(hn::Eq(magnitude, zero), one, magnitude);
  return {x, y, magnitude, hn::Div(x, divisor), hn::Div(y, divisor), integral};
}
HWY_INLINE M aligned_prepared(D d, const DirectionCenter& g, V hx, V hy,
                              double angle, V cosine, V sine);
// Bounded integer Sobel components need no per-lane integrality/range checks.
HWY_INLINE M aligned_integer(D d, V gx, V gy, V hx, V hy, double angle, V cosine, V sine) {
  const auto zero = hn::Zero(d);
  const auto gz = hn::And(hn::Eq(gx, zero), hn::Eq(gy, zero)), hz = hn::And(hn::Eq(hx, zero), hn::Eq(hy, zero));
  const auto bothzero = hn::And(gz, hz), nonzero = hn::Not(hn::Or(gz, hz));
  if (angle >= 0.5)
    return hn::Or(bothzero, nonzero);
  const auto cross = hn::Abs(hn::Sub(hn::Mul(gx, hy), hn::Mul(gy, hx)));
  const auto dot = hn::Abs(hn::Add(hn::Mul(gx, hx), hn::Mul(gy, hy)));
  auto accept = hn::Eq(cross, zero);
  if (angle == 0.25)
    accept = hn::Le(cross, dot);
  else if (angle != 0)
    accept = hn::Le(hn::Mul(cross, cosine), hn::Mul(dot, sine));
  return hn::Or(bothzero, hn::And(nonzero, accept));
}
template <bool Compact>
HWY_INLINE V cross_center(D d, const double* wide, const std::int32_t* compact,
                          std::size_t row_offset, std::size_t x) {
  if constexpr (Compact) {
    const hn::Rebind<std::int32_t, D> d32;
    return hn::PromoteTo(d, hn::LoadU(d32, compact + row_offset + x));
  } else {
    return hn::LoadU(d, wide + x);
  }
}
template <bool Compact, class Index>
HWY_INLINE V cross_gather(D d, const double* wide, const std::int32_t* compact, Index index) {
  const hn::Rebind<std::int32_t, D> d32;
  if constexpr (Compact)
    return hn::PromoteTo(d, hn::GatherIndex(d32, compact, index));
  else {
    const hn::Rebind<std::int64_t, D> di64;
    if (compact) {
      if constexpr (sizeof(hn::TFromV<Index>) == 4)
        return hn::PromoteTo(d, hn::GatherIndex(d32, compact, index));
      else
        return hn::PromoteTo(d, hn::GatherIndex(d32, compact, hn::DemoteTo(d32, index)));
    }
    if constexpr (sizeof(hn::TFromV<Index>) == 4)
      return hn::GatherIndex(d, wide, hn::PromoteTo(di64, index));
    else
      return hn::GatherIndex(d, wide, index);
  }
}
template <int SX, int SY, bool Fused = false, bool Integer = false, bool Compact = false, bool Narrow = Compact>
std::size_t collect_cross(const CrossSamplingRow& cross, const Row* output = nullptr,
                          const double* thresholds = nullptr, double angle = 0, double boost = 1) {
  // Output doubles may alias caller-owned thresholds. Snapshot the immutable
  // parameters so their conversions and broadcasts can stay outside the loop.
  double row_thresholds[3]{};
  if constexpr (Fused) {
    for (int j = 0; j < 3; ++j)
      row_thresholds[j] = thresholds[j];
    thresholds = row_thresholds;
  }
  double inverse_thresholds[3]{}, inverse_boosted[3]{};
  if constexpr (Fused && !Integer) {
    for (int j = 0; j < 3; ++j) {
      inverse_thresholds[j] = thresholds[j] == 0 ? 1 : 1 / thresholds[j];
      const double boosted = thresholds[j] * boost;
      inverse_boosted[j] = boosted == 0 ? 1 : 1 / boosted;
    }
  }
  const auto& row = cross.sampling;
  const D d;
  const hn::Rebind<std::conditional_t<Narrow, std::int32_t, std::int64_t>, D> di;
  const hn::Rebind<std::int32_t, D> d32;
  const hn::Rebind<std::uint8_t, D> db;
  const auto lanes = hn::Lanes(d), full = row.width / lanes * lanes;
  const auto* distances = reinterpret_cast<const std::uint8_t*>(row.distances);
  const auto index = hn::Iota(di, 0), last = hn::Set(di, static_cast<hn::TFromD<decltype(di)>>(row.width - 1));
  const auto vlimit = hn::Set(di, static_cast<hn::TFromD<decltype(di)>>(row.vertical_limit)),
             pitch = hn::Set(di, static_cast<hn::TFromD<decltype(di)>>(row.width));
  const auto* source = cross.compact_source ? nullptr : row.center - cross.row_offset;
  const auto* gx = cross.compact_gx ? nullptr : cross.gx - cross.row_offset;
  const auto* gy = cross.compact_gy ? nullptr : cross.gy - cross.row_offset;
  constexpr double pi = 3.1415926535897932384626433832795;
  const bool trig = Fused && angle != 0 && angle != 0.25 && angle < 0.5;
  const auto cosine = hn::Set(d, trig ? std::cos(pi * angle) : 0), sine = hn::Set(d, trig ? std::sin(pi * angle) : 0);

  for (std::size_t x = 0; x < full; x += lanes, distances += 2 * (std::size_t{1} << SX) * lanes) {
    hn::VFromD<decltype(db)> first, second;
    if constexpr (SX == 0)
      hn::LoadInterleaved2(db, distances, first, second);
    else {
      hn::VFromD<decltype(db)> skip1, skip2;
      hn::LoadInterleaved4(db, distances, first, second, skip1, skip2);
    }
    const auto small_distance = hn::PromoteTo(d32, first);
    auto distance = hn::Zero(di);
    if constexpr (Narrow)
      distance = small_distance;
    else
      distance = hn::PromoteTo(di, small_distance);
    const auto h = hn::ShiftRightSame(distance, SX), v = hn::ShiftRightSame(distance, SY);
    const auto position = hn::Add(index, hn::Set(di, static_cast<hn::TFromD<decltype(di)>>(x))),
               hlimit = hn::Min(position, hn::Sub(last, position));
    if (!hn::AllTrue(di, hn::And(hn::Le(h, hlimit), hn::Le(v, vlimit))))
      throw std::logic_error("reference outside plane");
    const auto absolute = hn::Add(position, hn::Set(di, static_cast<hn::TFromD<decltype(di)>>(cross.row_offset)));
    const auto vertical = hn::Mul(v, pitch);
    const hn::VFromD<decltype(di)> offsets[4] = {hn::Add(absolute, vertical), hn::Sub(absolute, vertical),
                                                 hn::Add(absolute, h), hn::Sub(absolute, h)};
    V r[4];
    for (int j = 0; j < 4; ++j) {
      r[j] = cross_gather<Compact>(d, source, cross.compact_source, offsets[j]);
      if constexpr (!Fused) {
        hn::StoreU(r[j], d, row.samples[j] + x);
        hn::StoreU(cross_gather<Compact>(d, gx, cross.compact_gx, offsets[j]), d, cross.grad_x[j] + x);
        hn::StoreU(cross_gather<Compact>(d, gy, cross.compact_gy, offsets[j]), d, cross.grad_y[j] + x);
      }
    }
    if constexpr (Fused) {
      const auto s = cross_center<Compact>(d, row.center, cross.compact_source, cross.row_offset, x),
                 zero = hn::Zero(d), one = hn::Set(d, 1.0);
      const auto average =
          Integer ? hn::Mul(hn::Add(hn::Add(r[0], r[1]), hn::Add(r[2], r[3])), hn::Set(d, 0.25)) : mean(d, r, 4, false);
      const auto da = hn::Abs(diff(d, average, s, Integer));
      auto maximum = zero;
      for (int j = 0; j < 4; ++j)
        maximum = hn::Max(maximum, hn::Abs(diff(d, r[j], s, Integer)));
      // A zero maximum-difference response rejects both possible boosts.
      // An integer mean equal to the center also makes every blend identical.
      // Neither case needs direction tests or gradient memory accesses.
      auto unchanged =
          hn::And(hn::Gt(maximum, zero), hn::Ge(maximum, hn::Set(d, std::max(thresholds[1], thresholds[1] * boost))));
      if constexpr (Integer)
        unchanged = hn::Or(unchanged, hn::Eq(average, s));
      if (hn::AllTrue(d, unchanged)) {
        const auto result = Integer ? hn::Add(s, load_integer_grain(d, output->integer_grain + x, lanes))
                                    : fp(d, hn::Add(s, load_grain(d, output->grain + x, lanes)));
        hn::StoreU(result, d, output->output + x);
        continue;
      }
      auto directions = hn::Eq(zero, zero);
      const auto cgx = cross_center<Compact>(d, cross.gx, cross.compact_gx, cross.row_offset, x),
                 cgy = cross_center<Compact>(d, cross.gy, cross.compact_gy, cross.row_offset, x);
      DirectionCenter prepared{};
      if constexpr (!Integer)
        prepared = prepare_direction(d, cgx, cgy);
      for (int j = 0; j < 4; ++j) {
        const auto hx = cross_gather<Compact>(d, gx, cross.compact_gx, offsets[j]),
                   hy = cross_gather<Compact>(d, gy, cross.compact_gy, offsets[j]);
        M direction;
        if constexpr (Compact)
          direction = aligned_integer(d, cgx, cgy, hx, hy, angle, cosine, sine);
        else if constexpr (!Integer)
          direction = cross.bounded_integer_gradients
                          ? aligned_integer(d, cgx, cgy, hx, hy, angle, cosine, sine)
                          : aligned_prepared(d, prepared, hx, hy, angle, cosine, sine);
        else
          direction = cross.bounded_integer_gradients
                          ? aligned_integer(d, cgx, cgy, hx, hy, angle, cosine, sine)
                          : aligned(d, cgx, cgy, hx, hy, angle, cosine, sine);
        directions = hn::And(directions, direction);
      }
      auto vertical = hn::Add(r[0], r[1]), horizontal = hn::Add(r[2], r[3]), twice = hn::Add(s, s);
      if constexpr (!Integer) {
        vertical = fp(d, vertical);
        horizontal = fp(d, horizontal);
        twice = fp(d, twice);
      }
      const auto qv = hn::Abs(diff(d, vertical, twice, Integer)), qh = hn::Abs(diff(d, horizontal, twice, Integer));
      const auto k = hn::IfThenElse(directions, hn::Set(d, boost), one);
      V c[4];
      if constexpr (!Integer) {
        V inverse[3];
        for (int j = 0; j < 3; ++j)
          inverse[j] = hn::IfThenElse(directions, hn::Set(d, inverse_boosted[j]), hn::Set(d, inverse_thresholds[j]));
        c[0] = response_with_inverse(d, da, hn::Mul(k, hn::Set(d, thresholds[0])), inverse[0]);
        c[1] = response_with_inverse(d, maximum, hn::Mul(k, hn::Set(d, thresholds[1])), inverse[1]);
        c[2] = response_with_inverse(d, qv, hn::Mul(k, hn::Set(d, thresholds[2])), inverse[2]);
        c[3] = response_with_inverse(d, qh, hn::Mul(k, hn::Set(d, thresholds[2])), inverse[2]);
      } else {
        c[0] = response(d, da, hn::Mul(k, hn::Set(d, thresholds[0])));
        c[1] = response(d, maximum, hn::Mul(k, hn::Set(d, thresholds[1])));
        c[2] = response(d, qv, hn::Mul(k, hn::Set(d, thresholds[2])));
        c[3] = response(d, qh, hn::Mul(k, hn::Set(d, thresholds[2])));
      }
      const auto f = weight(d, c);
      const auto blend = Integer ? hn::Floor(hn::Add(hn::MulAdd(hn::Sub(average, s), f, s), hn::Set(d, 0.5)))
                                 : fp(d, hn::Add(s, fp(d, hn::Mul(diff(d, average, s, false), f))));
      auto result = hn::IfThenElse(hn::Eq(f, zero), s, blend);
      if constexpr (Integer)
        result = hn::Add(result, load_integer_grain(d, output->integer_grain + x, lanes));
      else
        result = fp(d, hn::Add(result, load_grain(d, output->grain + x, lanes)));
      hn::StoreU(result, d, output->output + x);
    }
  }
  return full;
}
std::size_t collect_cross_target(const CrossSamplingRow& row) {
  if (row.sampling.sub_x > 1 || row.sampling.sub_y > 1 || row.sampling.width == 0)
    return 0;
  if (row.sampling.sub_x == 0)
    return row.sampling.sub_y == 0 ? collect_cross<0, 0>(row) : collect_cross<0, 1>(row);
  return row.sampling.sub_y == 0 ? collect_cross<1, 0>(row) : collect_cross<1, 1>(row);
}
template <bool Integer, bool Compact = false, bool Narrow = Compact>
std::size_t cross_dispatch(const CrossSamplingRow& row, const Row& output, const double* t, double angle,
                           double boost) {
  if (row.sampling.sub_x > 1 || row.sampling.sub_y > 1 || row.sampling.width == 0)
    return 0;
  if (row.sampling.sub_x == 0)
    return row.sampling.sub_y == 0 ? collect_cross<0, 0, true, Integer, Compact, Narrow>(row, &output, t, angle, boost)
                                   : collect_cross<0, 1, true, Integer, Compact, Narrow>(row, &output, t, angle, boost);
  return row.sampling.sub_y == 0 ? collect_cross<1, 0, true, Integer, Compact, Narrow>(row, &output, t, angle, boost)
                                 : collect_cross<1, 1, true, Integer, Compact, Narrow>(row, &output, t, angle, boost);
}
SamplingResult cross_sampled_target(const CrossSamplingRow& row, const Row& output, const double* t, bool integer,
                                    double angle, double boost) {
  // Compact planes fit signed 32-bit indices. Prove the whole possible
  // reference interval before selecting narrow coordinate arithmetic.
  const auto width = row.sampling.width;
  constexpr auto limit = std::size_t(std::numeric_limits<std::int32_t>::max());
  if (integer && row.bounded_integer_gradients && row.compact_source && row.compact_gx && row.compact_gy &&
      width && width <= limit &&
      row.row_offset <= limit - (width - 1) &&
      row.sampling.vertical_limit <= row.row_offset / width &&
      row.sampling.vertical_limit <= (limit - (width - 1) - row.row_offset) / width) {
    const auto n = cross_dispatch<true, true>(row, output, t, angle, boost);
    return {n, n};
  }
#if HWY_TARGET == HWY_AVX3_SPR
  if (!integer && width && width <= limit && row.row_offset <= limit - (width - 1) &&
      row.sampling.vertical_limit <= row.row_offset / width &&
      row.sampling.vertical_limit <= (limit - (width - 1) - row.row_offset) / width) {
    const auto n = cross_dispatch<false, false, true>(row, output, t, angle, boost);
    return {n, n};
  }
#endif
  const auto n = integer ? cross_dispatch<true>(row, output, t, angle, boost)
                         : cross_dispatch<false>(row, output, t, angle, boost);
  return {n, n};
}
HWY_INLINE V weight(D d, const V* c);
HWY_INLINE V response(D d, V delta, V threshold);
HWY_INLINE V weighted_half(D d, V center, const V* refs, const double* thresholds) {
  auto sum = hn::Zero(d), maximum = hn::Zero(d);
  for (int j = 0; j < 4; ++j) {
    sum = hn::Add(sum, refs[j]);
    maximum = hn::Max(maximum, hn::Abs(hn::Sub(refs[j], center)));
  }
  const auto average = hn::Mul(sum, hn::Set(d, 0.25));
  const auto da = hn::Abs(hn::Sub(average, center));
  const auto twice = hn::Add(center, center);
  const auto qv = hn::Abs(hn::Sub(hn::Add(refs[0], refs[1]), twice));
  const auto qh = hn::Abs(hn::Sub(hn::Add(refs[2], refs[3]), twice));
  const V c[4] = {response(d, da, hn::Set(d, thresholds[0])), response(d, maximum, hn::Set(d, thresholds[1])),
                  response(d, qv, hn::Set(d, thresholds[2])), response(d, qh, hn::Set(d, thresholds[2]))};
  const auto f = weight(d, c);
  return hn::IfThenElse(hn::Eq(f, hn::Zero(d)), center,
                        hn::Floor(hn::Add(hn::MulAdd(hn::Sub(average, center), f, center), hn::Set(d, 0.5))));
}
template <class DI>
HWY_INLINE auto weighted_integer_fast(DI di, hn::VFromD<DI> center,
                                     const hn::VFromD<DI>* refs, const float* thresholds,
                                     const float* slopes, const std::uint64_t* roots, float input_scale) {
  const hn::Rebind<float, DI> d;
  const auto zero = hn::Zero(d), one = hn::Set(d, 1.0f);
  const auto s = hn::ConvertTo(d, center);
  // Keep exact 16-bit sums and differences in int32 until the responses
  // need floats. The largest magnitude is 4 * 65535, well within int32.
  hn::VFromD<DI> signed_delta[4];
  for (int j = 0; j < 4; ++j) signed_delta[j] = hn::Sub(refs[j], center);
  const auto vertical_delta = hn::Add(signed_delta[0], signed_delta[1]);
  const auto horizontal_delta = hn::Add(signed_delta[2], signed_delta[3]);
  const auto maximum = hn::Max(hn::Max(hn::Abs(signed_delta[0]), hn::Abs(signed_delta[1])),
                               hn::Max(hn::Abs(signed_delta[2]), hn::Abs(signed_delta[3])));
  const auto blend_delta = hn::Mul(hn::ConvertTo(d, hn::Add(vertical_delta, horizontal_delta)), hn::Set(d, 0.25f));
  const hn::VFromD<decltype(d)> delta[4] = {
      hn::Abs(blend_delta), hn::ConvertTo(d, maximum),
      hn::ConvertTo(d, hn::Abs(vertical_delta)),
      hn::ConvertTo(d, hn::Abs(horizontal_delta))};
  hn::VFromD<decltype(d)> c[4];
  for (int j = 0; j < 4; ++j) {
    const int k = j < 2 ? j : 2;
    // Subtract before scaling: near the upper endpoint this difference is
    // exact, whereas 1 - delta / threshold loses relative accuracy.
    c[j] = hn::Min(one, hn::Max(zero, hn::Mul(
        hn::Sub(hn::Set(d, thresholds[k]), delta[j]), hn::Set(d, slopes[k]))));
  }
  const auto product = hn::Mul(hn::Mul(c[0], c[2]), hn::Mul(c[1], c[3]));
  // Integral thresholds in [1,65535] and quarter-integral differences make
  // every positive response >= 0.75/65535. Their product is normal binary32.
  const auto positive = hn::Gt(product, zero);
  const auto bits = hn::BitCast(di, hn::IfThenElse(positive, product, one));
  const auto index = hn::ShiftRight<15>(hn::Sub(bits, hn::Set(di, integer_weight_first_bits)));
  const auto fraction = hn::Mul(hn::ConvertTo(d, hn::And(bits, hn::Set(di, 0x7fff))), hn::Set(d, 0x1p-15f));
  // Read each base/slope pair once. Two half-width qword gathers service N
  // pixels; the exact endpoint difference was computed when building the table.
  const hn::Half<DI> half;
  const hn::Repartition<std::uint64_t, DI> dq;
  const hn::Rebind<std::int64_t, decltype(dq)> diq;
  const auto low = hn::BitCast(d, hn::GatherIndex(dq, roots, hn::PromoteTo(diq, hn::LowerHalf(half, index))));
  const auto high = hn::BitCast(d, hn::GatherIndex(dq, roots, hn::PromoteTo(diq, hn::UpperHalf(half, index))));
  const auto lower = hn::ConcatEven(d, high, low);
  const auto slope = hn::ConcatOdd(d, high, low);
  // The table shares endpoints, and positive-slope interpolation is monotone
  // even after rounding. The binary32 root error is below 2.4e-7; response
  // and product rounding keep the combined weight error below 2^-20.
  const auto f = hn::IfThenElseZero(positive, hn::MulAdd(fraction, slope, lower));
  // The rounded 16-bit blend fits int32; saturation is redundant.
  const auto blended = hn::Mul(hn::MulAdd(blend_delta, f, s), hn::Set(d, input_scale));
  return hn::ConvertInRangeTo(di, hn::Add(blended, hn::Set(d, 0.5f)));
}
template <class DI>
HWY_INLINE auto weighted_integer(DI d, hn::VFromD<DI> center, const hn::VFromD<DI>* refs, const double* thresholds) {
  const D dd;
  const hn::Half<DI> half;
  V low[4], high[4];
  for (int j = 0; j < 4; ++j) {
    low[j] = hn::PromoteTo(dd, hn::LowerHalf(half, refs[j]));
    high[j] = hn::PromoteTo(dd, hn::UpperHalf(half, refs[j]));
  }
  const auto lo = weighted_half(dd, hn::PromoteTo(dd, hn::LowerHalf(half, center)), low, thresholds);
  const auto hi = weighted_half(dd, hn::PromoteTo(dd, hn::UpperHalf(half, center)), high, thresholds);
  return hn::Combine(d, hn::DemoteTo(half, hi), hn::DemoteTo(half, lo));
}
template <int Mode, bool Blur, int SX = 0, int SY = 0, bool WordGather = false>
std::size_t byte_reference_samples(const ByteRow& input_row, const double* thresholds) {
  const ByteRow row = input_row;
  using D = std::conditional_t<WordGather, hn::Rebind<std::uint16_t, hn::ScalableTag<std::int32_t>>,
                               std::conditional_t<Mode == 2 && HWY_TARGET == HWY_AVX3_ZEN4,
                                                  hn::CappedTag<std::uint16_t, 16>, hn::ScalableTag<std::uint16_t>>>;
  const D d;
  const hn::Rebind<std::uint8_t, decltype(d)> db;
  const hn::Rebind<std::int16_t, decltype(d)> ds;
  const auto lanes = hn::Lanes(d), full = row.width / lanes * lanes;
  const auto* source = row.source;
  if constexpr (Mode == 3) {
    if (lanes < 32 && !row.collected[0])
      return 0;
  }

  constexpr int N = Mode == 1 || Mode == 3 ? 2 : 4;
  const auto main_limit = hn::Set(d, std::uint16_t(std::ceil(thresholds[0] / (N == 2 ? 128.0 : 64.0))));
  const auto pixel_limit = hn::Set(d, std::uint16_t(std::ceil(thresholds[0] / 256.0)));
  const auto pair_limit = hn::Set(d, std::uint16_t(std::ceil(thresholds[0] / 128.0)));
  for (std::size_t x = 0; x < full; x += lanes) {
    const auto center = hn::PromoteTo(d, hn::LoadU(db, source + x));
    hn::VFromD<decltype(d)> r[N];
    auto sum = hn::Zero(d);
    if constexpr (Mode == 2 || WordGather) {
      if constexpr (WordGather) {
        const hn::Rebind<std::int32_t, D> di;
        const hn::Rebind<std::uint32_t, D> du;
        hn::VFromD<decltype(db)> first, second;
        const auto* distance = reinterpret_cast<const std::uint8_t*>(row.distances) + 2 * (x << SX);
        if constexpr (SX == 0)
          hn::LoadInterleaved2(db, distance, first, second);
        else {
          hn::VFromD<decltype(db)> skip1, skip2;
          hn::LoadInterleaved4(db, distance, first, second, skip1, skip2);
        }
        const auto one = hn::PromoteTo(di, first), two = hn::PromoteTo(di, second);
        const auto stride = hn::Set(di, static_cast<std::int32_t>(row.stride));
        const auto h1 = hn::ShiftRight<SX>(one), h2 = hn::ShiftRight<SX>(two);
        const auto v1 = hn::Mul(hn::ShiftRight<SY>(one), stride), v2 = hn::Mul(hn::ShiftRight<SY>(two), stride);
        const auto column = hn::Iota(di, static_cast<std::int32_t>(x));
        const auto row_start = hn::Set(di, static_cast<std::int32_t>(row.row_offset));
        hn::VFromD<decltype(di)> columns[4], rows[4];
        if constexpr (Mode == 2) {
          columns[0] = hn::Add(column, h1);
          columns[1] = hn::Add(column, h2);
          columns[2] = hn::Sub(column, h1);
          columns[3] = hn::Sub(column, h2);
          rows[0] = hn::Add(row_start, v2);
          rows[1] = hn::Sub(row_start, v1);
          rows[2] = hn::Sub(row_start, v2);
          rows[3] = hn::Add(row_start, v1);
        } else {
          columns[0] = column;
          columns[1] = column;
          columns[2] = hn::Add(column, h1);
          columns[3] = hn::Sub(column, h1);
          rows[0] = hn::Add(row_start, v1);
          rows[1] = hn::Sub(row_start, v1);
          rows[2] = row_start;
          rows[3] = row_start;
        }
        for (int j = 0; j < N; ++j) {
          // Slide the final word inward within its visible row. Every load
          // covers four actual pixels, without requiring input padding.
          const auto first = hn::Min(columns[j], hn::Set(di, static_cast<std::int32_t>(row.width - 4)));
          const auto start = hn::Add(rows[j], first);
          const auto shift = hn::BitCast(du, hn::ShiftLeft<3>(hn::Sub(columns[j], first)));
          const auto bits = hn::GatherOffset(du, reinterpret_cast<const std::uint32_t*>(row.base), start);
          const auto byte = HWY_IS_LITTLE_ENDIAN
                                ? hn::And(hn::Shr(bits, shift), hn::Set(du, 255))
                                : hn::And(hn::Shr(bits, hn::Sub(hn::Set(du, 24), shift)), hn::Set(du, 255));
          r[j] = hn::DemoteTo(d, byte);
          sum = hn::Add(sum, r[j]);
        }
      } else {
        HWY_ALIGN std::uint8_t gathered[4][hn::MaxLanes(d)];
        for (std::size_t i = 0; i < lanes; ++i) {
          const auto distance = row.distances[(x + i) << SX];
          const auto h1 = std::ptrdiff_t(distance.first >> SX), v1 = std::ptrdiff_t(distance.first >> SY);
          const auto h2 = std::ptrdiff_t(distance.second >> SX), v2 = std::ptrdiff_t(distance.second >> SY);
          const auto a = v2 * std::ptrdiff_t(row.stride) + h1, b = h2 - v1 * std::ptrdiff_t(row.stride);
          const auto* pixel = source + x + i;
          gathered[0][i] = pixel[a];
          gathered[1][i] = pixel[b];
          gathered[2][i] = pixel[-a];
          gathered[3][i] = pixel[-b];
        }
        for (int j = 0; j < 4; ++j) {
          r[j] = hn::PromoteTo(d, hn::LoadU(db, gathered[j]));
          sum = hn::Add(sum, r[j]);
        }
      }
    } else if constexpr (Mode == 3) {
      if (lanes >= 32) {
        hn::VFromD<decltype(db)> first, second;
        const auto* distances = reinterpret_cast<const std::uint8_t*>(row.distances) + (x << SX) * 2;
        if constexpr (SX == 0)
          hn::LoadInterleaved2(db, distances, first, second);
        else {
          hn::VFromD<decltype(db)> skip1, skip2;
          hn::LoadInterleaved4(db, distances, first, second, skip1, skip2);
        }
        const auto h = hn::ShiftRight<SX>(hn::PromoteTo(d, first));
        // Both complete neighbouring vectors stay inside the source row.
        if (x >= lanes && row.width - x >= 2 * lanes && hn::AllTrue(d, hn::Le(h, hn::Set(d, std::uint16_t(lanes))))) {
          const auto left = hn::PromoteTo(d, hn::LoadU(db, source + x - lanes));
          const auto right = hn::PromoteTo(d, hn::LoadU(db, source + x + lanes));
          const auto index = hn::Iota(d, 0);
          r[0] = hn::TwoTablesLookupLanes(center, right, hn::IndicesFromVec(d, hn::Add(index, h)));
          r[1] = hn::TwoTablesLookupLanes(
              left, center, hn::IndicesFromVec(d, hn::Sub(hn::Add(index, hn::Set(d, std::uint16_t(lanes))), h)));
        } else {
          HWY_ALIGN std::uint8_t plus[hn::MaxLanes(d)], minus[hn::MaxLanes(d)];
          for (std::size_t i = 0; i < lanes; ++i) {
            const auto at = x + i, offset = std::size_t(row.distances[at << SX].first >> SX);
            plus[i] = source[at + offset];
            minus[i] = source[at - offset];
          }
          r[0] = hn::PromoteTo(d, hn::LoadU(db, plus));
          r[1] = hn::PromoteTo(d, hn::LoadU(db, minus));
        }
      } else {
        r[0] = hn::PromoteTo(d, hn::LoadU(db, row.collected[0] + x));
        r[1] = hn::PromoteTo(d, hn::LoadU(db, row.collected[1] + x));
      }
      sum = hn::Add(r[0], r[1]);
    } else {
      for (int j = 0; j < N; ++j) {
        r[j] = hn::PromoteTo(d, hn::LoadU(db, row.collected[j] + x));
        sum = hn::Add(sum, r[j]);
      }
    }
    const auto scaled_center = hn::ShiftLeft<(N == 2 ? 1 : 2)>(center);
    const auto delta = hn::Sub(hn::Max(sum, scaled_center), hn::Min(sum, scaled_center));
    auto accepted = hn::Lt(delta, main_limit);
    if constexpr (!Blur && Mode <= 3) {
      accepted = hn::Eq(center, center);
      for (int j = 0; j < N; ++j)
        accepted = hn::And(accepted, hn::Lt(hn::Sub(hn::Max(r[j], center), hn::Min(r[j], center)), pixel_limit));
    }
    if constexpr (Mode == 4) {
      hn::VFromD<decltype(d)> pairs[2];
      const auto twice = hn::ShiftLeft<1>(center);
      for (int j : {0, 2}) {
        const auto pair = hn::Add(r[j], r[j + 1]);
        auto accept = hn::Lt(hn::Sub(hn::Max(pair, twice), hn::Min(pair, twice)), pair_limit);
        if constexpr (!Blur) {
          accept = hn::And(hn::Lt(hn::Sub(hn::Max(r[j], center), hn::Min(r[j], center)), pixel_limit),
                           hn::Lt(hn::Sub(hn::Max(r[j + 1], center), hn::Min(r[j + 1], center)), pixel_limit));
        }
        pairs[j / 2] = hn::IfThenElse(accept, pair, twice);
      }
      sum = hn::Add(pairs[0], pairs[1]);
      accepted = hn::Eq(center, center);
    }
    auto value = hn::ShiftLeft<(N == 2 ? 7 : 6)>(hn::IfThenElse(accepted, sum, scaled_center));
    const auto grain = hn::LoadU(ds, row.grain + x);
    value = hn::SaturatedAdd(value, hn::BitCast(d, hn::Max(hn::Zero(ds), grain)));
    value = hn::SaturatedSub(value, hn::BitCast(d, hn::Neg(hn::Min(hn::Zero(ds), grain))));
    auto code = hn::ShiftRight<8>(value);
    const auto remainder = hn::And(value, hn::Set(d, 255));
    const auto up =
        hn::Or(hn::Gt(remainder, hn::Set(d, 128)),
               hn::And(hn::Eq(remainder, hn::Set(d, 128)), hn::Ne(hn::And(code, hn::Set(d, 1)), hn::Zero(d))));
    code = hn::Add(code, hn::IfThenElseZero(up, hn::Set(d, 1)));
    if (row.limited)
      code = hn::Min(hn::Set(d, row.chroma ? 240 : 235), hn::Max(hn::Set(d, 16), code));
    hn::StoreU(hn::DemoteTo(db, code), db, row.output + x);
  }
  return full;
}
template <int Mode, bool Blur, int SX, int SY>
std::size_t byte_cross(const ByteRow& row, const double* thresholds) {
  // Wide gathers avoid serial byte insertion. Zen 4 benefits from compact
  // byte loads in 16-pixel blocks; narrow targets also retain byte collection.
  if constexpr (HWY_TARGET != HWY_AVX3_ZEN4 && hn::MaxLanes(hn::ScalableTag<std::int32_t>()) >= 16) {
    if (row.width >= 4 && row.extent <= std::size_t(std::numeric_limits<std::int32_t>::max()) &&
        row.stride <= std::size_t(std::numeric_limits<std::int32_t>::max()))
      return byte_reference_samples<Mode, Blur, SX, SY, true>(row, thresholds);
  }
  if constexpr (Mode != 2) {
    if (!row.collected[0])
      return 0;
  }
  return byte_reference_samples<Mode, Blur, SX, SY>(row, thresholds);
}
template <int Mode, bool Blur>
std::size_t byte_cross_layout(const ByteRow& row, const double* t) {
  if (row.sub_x == 0 && row.sub_y == 0)
    return byte_cross<Mode, Blur, 0, 0>(row, t);
  if (row.sub_x == 1 && row.sub_y == 0)
    return byte_cross<Mode, Blur, 1, 0>(row, t);
  if (row.sub_x == 0 && row.sub_y == 1)
    return byte_cross<Mode, Blur, 0, 1>(row, t);
  return byte_cross<Mode, Blur, 1, 1>(row, t);
}
std::size_t byte_target(const ByteRow& row, int mode, const double* t, bool blur) {
  switch (mode) {
    case 1:
      return blur ? byte_cross_layout<1, true>(row, t) : byte_cross_layout<1, false>(row, t);
    case 2:
      return blur ? byte_cross_layout<2, true>(row, t) : byte_cross_layout<2, false>(row, t);
    case 3:
      if (row.sub_x == 0)
        return blur ? byte_reference_samples<3, true, 0>(row, t) : byte_reference_samples<3, false, 0>(row, t);
      return blur ? byte_reference_samples<3, true, 1>(row, t) : byte_reference_samples<3, false, 1>(row, t);
    case 4:
      return blur ? byte_cross_layout<4, true>(row, t) : byte_cross_layout<4, false>(row, t);
    default:
      return 0;
  }
}
#if HWY_TARGET == HWY_AVX3 || HWY_TARGET == HWY_AVX3_DL || HWY_TARGET == HWY_AVX3_ZEN4 || HWY_TARGET == HWY_AVX3_SPR || \
    HWY_TARGET == HWY_AVX2 || HWY_TARGET == HWY_SSE4 || HWY_TARGET == HWY_SSSE3
template <class V>
HWY_INLINE V word_absdiff(V a, V b) {
  return hn::Sub(hn::Max(a, b), hn::Min(a, b));
}
template <class D>
HWY_INLINE hn::VFromD<D> native_word_load(D d, const std::uint8_t* source, int depth) {
  const hn::Rebind<std::uint8_t, D> db;
  const hn::Repartition<std::uint8_t, D> raw;
  return hn::ShiftLeftSame(depth == 8 ? hn::PromoteTo(d, hn::LoadU(db, source))
                                     : hn::BitCast(d, hn::LoadU(raw, source)), 16 - depth);
}
template <int Mode, bool Blur, int SX, int SY>
std::size_t native_word_pairs(const IntegerRow& input_row, const double* thresholds) {
  static_assert(Mode == 1 || Mode == 3);
  const IntegerRow row = input_row;
  const int bytes = row.input_depth == 8 ? 1 : 2;
  const hn::ScalableTag<std::uint16_t> d;
  const hn::Repartition<std::int32_t, decltype(d)> di;
  const hn::Rebind<std::uint32_t, decltype(di)> du;
  const hn::Rebind<std::uint16_t, decltype(di)> dh;
  const hn::Rebind<std::uint8_t, decltype(d)> db;
  const hn::Repartition<std::uint8_t, decltype(d)> raw;
  const hn::Rebind<std::int16_t, decltype(d)> dg;
  const auto lanes = hn::Lanes(d), half = hn::Lanes(di), full = row.width / lanes * lanes;
  const auto threshold = hn::Set(d, std::uint16_t(std::ceil(thresholds[0])));
  const auto pitch = hn::Set(di, std::int32_t(row.native_stride)),
             vlimit = hn::Set(di, std::int32_t(row.vertical_limit));
  const auto last = hn::Set(di, std::int32_t(row.width - 1));
  constexpr int N = 2;
  for (std::size_t x = 0; x < full; x += lanes) {
    hn::VFromD<decltype(db)> first, second;
    const auto* distances = reinterpret_cast<const std::uint8_t*>(row.distances) + (x << SX) * 2;
    if constexpr (SX == 0)
      hn::LoadInterleaved2(db, distances, first, second);
    else {
      hn::VFromD<decltype(db)> a, b;
      hn::LoadInterleaved4(db, distances, first, second, a, b);
    }
    const auto d1 = hn::PromoteTo(d, first);
    hn::VFromD<decltype(d)> refs[N];
    bool tabled = false;
    if constexpr (Mode == 3) {
      const auto h = hn::ShiftRight<SX>(d1);
      // The three loads cover the row interior; the distance check proves
      // both selected neighbours lie in those tables. Edges use gathers.
      if (x >= lanes && row.width - x >= 2 * lanes && hn::AllTrue(d, hn::Le(h, hn::Set(d, std::uint16_t(lanes))))) {
        const auto* source = row.native + row.row_offset + x * bytes;
        const auto left = native_word_load(d, source - lanes * bytes, row.input_depth),
                   center = native_word_load(d, source, row.input_depth),
                   right = native_word_load(d, source + lanes * bytes, row.input_depth);
        const auto index = hn::Iota(d, 0);
        refs[0] = hn::TwoTablesLookupLanes(center, right, hn::IndicesFromVec(d, hn::Add(index, h)));
        refs[1] = hn::TwoTablesLookupLanes(
            left, center, hn::IndicesFromVec(d, hn::Sub(hn::Add(index, hn::Set(d, std::uint16_t(lanes))), h)));
        tabled = true;
      }
    }
    if (!tabled) {
      hn::VFromD<decltype(dh)> gathered[2][N];
      for (int part = 0; part < 2; ++part) {
        const auto dist1 = hn::PromoteTo(di, part == 0 ? hn::LowerHalf(dh, d1) : hn::UpperHalf(dh, d1));
        const auto h1 = hn::ShiftRight<SX>(dist1), v1 = hn::ShiftRight<SY>(dist1);
        const auto pos = hn::Iota(di, std::int32_t(x + part * half)), hlimit = hn::Min(pos, hn::Sub(last, pos));
        auto valid = hn::Eq(pos, pos);
        if constexpr (Mode != 3)
          valid = hn::And(valid, hn::Le(v1, vlimit));
        if constexpr (Mode != 1)
          valid = hn::And(valid, hn::Le(h1, hlimit));
        if (!hn::AllTrue(di, valid))
          throw std::logic_error("reference outside plane");
        const auto absolute = hn::Add(hn::Mul(pos, hn::Set(di, bytes)), hn::Set(di, std::int32_t(row.row_offset)));
        hn::VFromD<decltype(di)> offsets[N];
        const auto off = Mode == 1 ? hn::Mul(v1, pitch) : hn::Mul(h1, hn::Set(di, bytes));
        offsets[0] = hn::Add(absolute, off);
        offsets[1] = hn::Sub(absolute, off);
        for (int j = 0; j < N; ++j) {
          // Clamp the four-byte load, then extract the requested sample. This
          // includes the last visible sample without requiring source padding.
          const auto start = hn::Min(offsets[j], hn::Set(di, std::int32_t(row.native_extent - 4)));
          const auto shift = hn::BitCast(du, hn::ShiftLeft<3>(hn::Sub(offsets[j], start)));
          const auto words = hn::GatherOffset(du, reinterpret_cast<const std::uint32_t*>(row.native), start);
          const auto native_shift = HWY_IS_LITTLE_ENDIAN ? shift : hn::Sub(hn::Set(du, 32 - 8 * bytes), shift);
          const auto value = hn::And(hn::Shr(words, native_shift), hn::Set(du, bytes == 1 ? 255u : 65535u));
          gathered[part][j] = hn::DemoteTo(dh, hn::ShiftLeftSame(hn::BitCast(di, value), 16 - row.input_depth));
        }
      }
      for (int j = 0; j < N; ++j)
        refs[j] = hn::Combine(d, gathered[1][j], gathered[0][j]);
    }
    const auto* source = row.native + row.row_offset + x * bytes;
    const auto center = native_word_load(d, source, row.input_depth);
    auto mean = hn::AverageRound(refs[0], refs[1]);
    auto accepted = hn::Lt(word_absdiff(mean, center), threshold);
    if constexpr (!Blur && Mode <= 3) {
      accepted = hn::Eq(center, center);
      for (int j = 0; j < N; ++j)
        accepted = hn::And(accepted, hn::Lt(word_absdiff(refs[j], center), threshold));
    }
    auto value = hn::IfThenElse(accepted, mean, center);
    // Saturation to [0,65535] commutes with the final output clamp.
    const auto grain = hn::LoadU(dg, row.grain + x);
    value = hn::SaturatedAdd(value, hn::BitCast(d, hn::Max(hn::Zero(dg), grain)));
    value = hn::SaturatedSub(value, hn::BitCast(d, hn::Neg(hn::Min(hn::Zero(dg), grain))));
    const int shift = 16 - row.output_depth;
    auto code = value;
    if (shift) {
      const auto q = hn::ShiftRightSame(value, shift), rem = hn::Sub(value, hn::ShiftLeftSame(q, shift));
      const auto midpoint = hn::Set(d, std::uint16_t(1u << (shift - 1)));
      const auto up =
          hn::Or(hn::Gt(rem, midpoint), hn::And(hn::Eq(rem, midpoint), hn::Ne(hn::And(q, hn::Set(d, 1)), hn::Zero(d))));
      code = hn::Add(q, hn::IfThenElseZero(up, hn::Set(d, 1)));
    }
    const auto low = std::uint16_t(row.limited ? 16 << (row.output_depth - 8) : 0),
               high = std::uint16_t(row.limited ? (row.chroma ? 240 : 235) << (row.output_depth - 8)
                                                : (1u << row.output_depth) - 1);
    code = hn::Min(hn::Set(d, high), hn::Max(hn::Set(d, low), code));
    if (row.output_depth == 8)
      hn::StoreU(hn::DemoteTo(db, code), db, row.packed_output + x);
    else
      hn::StoreU(hn::BitCast(raw, code), raw, row.packed_output + 2 * x);
  }
  return full;
}
#endif

// Keep SIMD operations in the per-target namespace so their target attributes
// apply to every instantiation, including clang-cl builds.
template <bool DirectNative, bool FastWeights, bool RelativeGather, class D32, class Index>
HWY_INLINE hn::VFromD<D32> integer_gather(D32 d, const IntegerRow& row, Index position, int input_depth, int bytes) {
  const hn::Rebind<std::int32_t, D32> di;
  const hn::Rebind<std::uint32_t, D32> du;
  if constexpr (DirectNative) {
    // The row dispatcher proved that these four-byte reads fit the plane.
    // Coordinate validation below still precedes every reference read.
    auto bits = hn::GatherOffset(
        du, reinterpret_cast<const std::uint32_t*>(row.native + (RelativeGather ? row.row_offset : 0)), position);
    if constexpr (!HWY_IS_LITTLE_ENDIAN)
      bits = hn::ShiftRightSame(bits, 32 - 8 * bytes);
    const auto code = hn::And(bits, hn::Set(du, bytes == 1 ? 255u : 65535u));
    if constexpr (FastWeights)
      return hn::BitCast(d, code);
    return hn::ShiftLeftSame(hn::BitCast(d, code), 16 - input_depth);
  }
  if (!row.native)
    return hn::GatherIndex(d, row.center - row.row_offset, position);
  const auto start = hn::Min(position, hn::Set(di, std::int32_t(row.native_extent - 4)));
  const auto shift = hn::BitCast(du, hn::ShiftLeft<3>(hn::Sub(position, start)));
  const auto bits = hn::GatherOffset(du, reinterpret_cast<const std::uint32_t*>(row.native), start);
  const auto native_shift = HWY_IS_LITTLE_ENDIAN ? shift : hn::Sub(hn::Set(du, 32 - 8 * bytes), shift);
  const auto code = hn::And(hn::Shr(bits, native_shift), hn::Set(du, bytes == 1 ? 255u : 65535u));
  return hn::ShiftLeftSame(hn::BitCast(d, code), 16 - input_depth);
}
template <bool DirectNative, bool FastWeights, class D32>
HWY_INLINE hn::VFromD<D32> integer_center(D32 d, const IntegerRow& row, std::size_t x, int input_depth, int bytes) {
  const hn::Rebind<std::uint8_t, D32> db;
  if constexpr (!DirectNative) {
    if (!row.native)
      return hn::LoadU(d, row.center + x);
  }
  const auto* source = row.native + row.row_offset + x * bytes;
  const hn::Rebind<std::uint16_t, D32> dw;
  const hn::Repartition<std::uint8_t, decltype(dw)> dwb;
  const auto code =
      bytes == 1 ? hn::PromoteTo(d, hn::LoadU(db, source)) : hn::PromoteTo(d, hn::BitCast(dw, hn::LoadU(dwb, source)));
  if constexpr (FastWeights && DirectNative)
    return code;
  return hn::ShiftLeftSame(code, 16 - input_depth);
}

template <int Mode, bool Blur, int SX, int SY, bool DirectNative = false, bool FastWeights = false,
          int PackedDither = -1, bool ValidatedDistances = false, bool CachedOffsets = false, int NativeDepth = 0,
          int OutputDepth = 0>
std::size_t integer_samples(const IntegerRow& input_row, const double* thresholds) {
  // Output byte stores may alias caller-owned metadata. Keep the immutable
  // descriptor local so the sampling loop can retain its invariants.
  const IntegerRow row = input_row;
  const int input_depth = NativeDepth ? NativeDepth : row.input_depth;
  const int output_depth = OutputDepth ? OutputDepth : row.output_depth;
  const hn::ScalableTag<std::int32_t> d;
  const hn::Rebind<std::int32_t, decltype(d)> di;
  const hn::Rebind<std::uint8_t, decltype(d)> db;
  const auto lanes = hn::Lanes(d), full = row.width / lanes * lanes;
  const auto* distances = reinterpret_cast<const std::uint8_t*>(row.distances);
  const auto step = std::size_t{1} << SX;
  const auto index = hn::Iota(di, 0), last = hn::Set(di, std::int32_t(row.width - 1)),
             vlimit = hn::Set(di, std::int32_t(row.vertical_limit));
  const auto bytes = DirectNative || row.native ? (input_depth == 8 ? 1 : 2) : 1;
  const auto pitch = hn::Set(di, std::int32_t(row.native ? row.native_stride : row.width));
  // These 512-bit targets use native gathers with signed byte offsets.
  // Generic Highway gathers require nonnegative offsets, so all other
  // targets retain the plane origin even when reference offsets are cached.
  constexpr bool relative_gather = CachedOffsets && (HWY_TARGET == HWY_AVX3 || HWY_TARGET == HWY_AVX3_DL ||
                                                     HWY_TARGET == HWY_AVX3_ZEN4 || HWY_TARGET == HWY_AVX3_SPR);
  // Integer differences satisfy delta < t exactly when delta < ceil(t).
  const auto limit = hn::Set(d, std::int32_t(std::ceil(thresholds[0])));
  auto maximum_limit = hn::Zero(d), pair_limit = hn::Zero(d);
  if constexpr (Mode == 5) {
    maximum_limit = hn::Set(d, std::int32_t(std::ceil(thresholds[1])));
    pair_limit = hn::Set(d, std::int32_t(std::ceil(thresholds[2])));
  }
  // Power-of-two unit changes are exact. Borrow native codes directly and
  // scale the blend once, instead of shifting five input vectors per block.
  const float input_scale = FastWeights && DirectNative ? float(1u << (16 - input_depth)) : 1.0f;
  float weight_thresholds[3]{}, weight_slopes[3]{};
  if constexpr (FastWeights) {
    for (int j = 0; j < 3; ++j) {
      weight_thresholds[j] = static_cast<float>(thresholds[j]) / input_scale;
      weight_slopes[j] = static_cast<float>(3.0 / thresholds[j]) * input_scale;
    }
  }
  const std::uint64_t* weight_roots = FastWeights ? integer_weight_table().data() : nullptr;
  // The Bayer matrix repeats every 16 columns. Compute the exact integer
  // bias once per row: floor((2*bayer+1)*step/512).
  std::int32_t bayer_bias[16]{};
  if ((PackedDither >= 0 ? PackedDither : row.dither) == 2 && output_depth < 16) {
    for (int column = 0; column < 16; ++column) {
      int value = 0;
      for (int bit = 0; bit < 4; ++bit)
        value |= ((2 * ((column >> bit) & 1)) ^ (3 * ((row.y >> bit) & 1))) << (2 * (3 - bit));
      bayer_bias[column] = ((2 * value + 1) << (16 - output_depth)) >> 9;
    }
  }
  for (std::size_t x = 0; x < full; x += lanes, distances += 2 * step * lanes) {
    const auto value = integer_center<DirectNative, FastWeights>(d, row, x, input_depth, bytes);
    auto debanded = value;
    if constexpr (Mode != 0) {
      constexpr int N = Mode == 1 || Mode == 3 ? 2 : 4;
      hn::VFromD<decltype(d)> refs[N];
      if constexpr (CachedOffsets) {
        static_assert(Mode == 6 && DirectNative && ValidatedDistances);
        const auto pos = hn::Add(index, hn::Set(di, std::int32_t(x)));
        const auto absolute =
            hn::Add(hn::Mul(pos, hn::Set(di, bytes)), hn::Set(di, relative_gather ? 0 : std::int32_t(row.row_offset)));
        const auto vertical = hn::LoadU(d, row.vertical_offsets + x);
        const auto horizontal = hn::LoadU(d, row.horizontal_offsets + x);
        refs[0] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Add(absolute, vertical),
                                                                             input_depth, bytes);
        refs[1] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Sub(absolute, vertical),
                                                                             input_depth, bytes);
        refs[2] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Add(absolute, horizontal),
                                                                             input_depth, bytes);
        refs[3] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Sub(absolute, horizontal),
                                                                             input_depth, bytes);
      } else {
        hn::VFromD<decltype(db)> first, second;
        if (SX == 0)
          hn::LoadInterleaved2(db, distances, first, second);
        else {
          hn::VFromD<decltype(db)> skip1, skip2;
          hn::LoadInterleaved4(db, distances, first, second, skip1, skip2);
        }
        const auto distance1 = hn::PromoteTo(di, first), h1 = hn::ShiftRightSame(distance1, SX),
                   v1 = hn::ShiftRightSame(distance1, SY);
        const auto pos = hn::Add(index, hn::Set(di, std::int32_t(x))), hlimit = hn::Min(pos, hn::Sub(last, pos));
        auto valid = hn::Eq(index, index);
        if constexpr (Mode != 3)
          valid = hn::And(valid, hn::Le(v1, vlimit));
        if constexpr (Mode != 1)
          valid = hn::And(valid, hn::Le(h1, hlimit));
        if constexpr (!ValidatedDistances) {
          if (!hn::AllTrue(di, valid))
            throw std::logic_error("reference outside plane");
        }
        const auto absolute = hn::Add(hn::Mul(pos, hn::Set(di, bytes)), hn::Set(di, std::int32_t(row.row_offset)));
        if constexpr (Mode == 1 || Mode == 3) {
          bool tabled = false;
          if constexpr (Mode == 3) {
            // Two adjacent vectors cover each side when every distance fits one
            // vector. Loads stay inside this row, including at image borders.
            if (lanes >= 16 && x >= lanes && row.width - x >= 2 * lanes &&
                hn::AllTrue(di, hn::Le(h1, hn::Set(di, std::int32_t(lanes))))) {
              const auto left = integer_center<DirectNative, FastWeights>(d, row, x - lanes, input_depth, bytes),
                         right = integer_center<DirectNative, FastWeights>(d, row, x + lanes, input_depth, bytes);
              refs[0] = hn::TwoTablesLookupLanes(value, right, hn::IndicesFromVec(d, hn::Add(index, h1)));
              refs[1] = hn::TwoTablesLookupLanes(
                  left, value, hn::IndicesFromVec(d, hn::Sub(hn::Add(index, hn::Set(di, std::int32_t(lanes))), h1)));
              tabled = true;
            }
          }
          if (!tabled) {
            const auto off = Mode == 1 ? hn::Mul(v1, pitch) : hn::Mul(h1, hn::Set(di, bytes));
            refs[0] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Add(absolute, off),
                                                                                 input_depth, bytes);
            refs[1] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Sub(absolute, off),
                                                                                 input_depth, bytes);
          }
        } else if constexpr (Mode == 2) {
          const auto distance2 = hn::PromoteTo(di, second), h2 = hn::ShiftRightSame(distance2, SX),
                     v2 = hn::ShiftRightSame(distance2, SY);
          if (!hn::AllTrue(di, hn::And(hn::Le(h2, hlimit), hn::Le(v2, vlimit))))
            throw std::logic_error("reference outside plane");
          const auto a = hn::Add(hn::Mul(v2, pitch), hn::Mul(h1, hn::Set(di, bytes))),
                     b = hn::Sub(hn::Mul(h2, hn::Set(di, bytes)), hn::Mul(v1, pitch));
          refs[0] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Add(absolute, a),
                                                                               input_depth, bytes);
          refs[1] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Add(absolute, b),
                                                                               input_depth, bytes);
          refs[2] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Sub(absolute, a),
                                                                               input_depth, bytes);
          refs[3] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Sub(absolute, b),
                                                                               input_depth, bytes);
        }
        if constexpr (Mode >= 4) {
          const auto vertical = hn::Mul(v1, pitch);
          refs[0] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Add(absolute, vertical),
                                                                               input_depth, bytes);
          refs[1] = integer_gather<DirectNative, FastWeights, relative_gather>(d, row, hn::Sub(absolute, vertical),
                                                                               input_depth, bytes);
          refs[2] = integer_gather<DirectNative, FastWeights, relative_gather>(
              d, row, hn::Add(absolute, hn::Mul(h1, hn::Set(di, bytes))), input_depth, bytes);
          refs[3] = integer_gather<DirectNative, FastWeights, relative_gather>(
              d, row, hn::Sub(absolute, hn::Mul(h1, hn::Set(di, bytes))), input_depth, bytes);
        }
      }
      auto sum = hn::Zero(d);
      for (int j = 0; j < N; ++j)
        sum = hn::Add(sum, refs[j]);
      auto mean = hn::ShiftRight<(N == 4 ? 2 : 1)>(hn::Add(sum, hn::Set(d, N / 2)));
      auto accept = hn::Lt(hn::Abs(hn::Sub(mean, value)), limit);
      if constexpr (!Blur && Mode <= 3) {
        accept = hn::Eq(value, value);
        for (int j = 0; j < N; ++j)
          accept = hn::And(accept, hn::Lt(hn::Abs(hn::Sub(refs[j], value)), limit));
      }
      if constexpr (Mode == 4) {
        hn::VFromD<decltype(d)> results[2];
        for (int j : {0, 2}) {
          const auto average = hn::ShiftRight<1>(hn::Add(hn::Add(refs[j], refs[j + 1]), hn::Set(d, 1)));
          auto accepted = hn::Lt(hn::Abs(hn::Sub(average, value)), limit);
          if constexpr (!Blur)
            accepted = hn::And(hn::Lt(hn::Abs(hn::Sub(refs[j], value)), limit),
                               hn::Lt(hn::Abs(hn::Sub(refs[j + 1], value)), limit));
          results[j / 2] = hn::IfThenElse(accepted, average, value);
        }
        mean = hn::ShiftRight<1>(hn::Add(hn::Add(results[0], results[1]), hn::Set(d, 1)));
        accept = hn::Eq(value, value);
      } else if constexpr (Mode == 5) {
        for (int j = 0; j < 4; ++j)
          accept = hn::And(accept, hn::Lt(hn::Abs(hn::Sub(refs[j], value)), maximum_limit));
        const auto twice = hn::Add(value, value);
        for (int j : {0, 2})
          accept = hn::And(accept, hn::Lt(hn::Abs(hn::Sub(hn::Add(refs[j], refs[j + 1]), twice)), pair_limit));
      }
      if constexpr (Mode == 6) {
        if constexpr (FastWeights)
          mean = weighted_integer_fast(d, value, refs, weight_thresholds, weight_slopes, weight_roots, input_scale);
        else
          mean = weighted_integer(d, value, refs, thresholds);
        accept = hn::Eq(value, value);
      }
      debanded = hn::IfThenElse(accept, mean, value);
    }
    const hn::Rebind<std::int16_t, decltype(d)> dg;
    const auto result = hn::Add(debanded, hn::PromoteTo(d, hn::LoadU(dg, row.grain + x)));
    if (PackedDither >= 0 ? PackedDither != 0 : row.packed_output != nullptr) {
      const int shift = 16 - output_depth;
      auto code = result;
      if (shift && (PackedDither >= 0 ? PackedDither : row.dither) == 2) {
        hn::VFromD<decltype(d)> bias;
        // Fixed-width blocks divide the 16-column period and never wrap
        // within a block. Wider/scalable blocks retain the indexed lookup.
        if (lanes <= 16 && 16 % lanes == 0)
          bias = hn::LoadU(d, bayer_bias + (x & 15));
        else {
          const auto positions = hn::And(hn::Iota(d, static_cast<std::int32_t>(x)), hn::Set(d, 15));
          bias = hn::GatherIndex(d, bayer_bias, positions);
        }
        code = hn::ShiftRightSame(hn::Add(result, bias), shift);
      } else if (shift) {
        const auto quotient = hn::ShiftRightSame(result, shift);
        // Round ties to even using an integer bias, including negative grain.
        // For result=q*2^shift+r, add half-1+(q&1) before the final shift.
        // The processing range plus bias is far inside int32.
        const auto bias = hn::Add(hn::Set(d, (std::int32_t{1} << (shift - 1)) - 1),
                                  hn::And(quotient, hn::Set(d, 1)));
        code = hn::ShiftRightSame(hn::Add(result, bias), shift);
      }
      // DemoteTo saturates to the unsigned storage range. Only limited
      // range and sub-16-bit word formats require an explicit bound.
      if (row.limited) {
        const int low = 16 << (output_depth - 8);
        const int high = (row.chroma ? 240 : 235) << (output_depth - 8);
        code = hn::Min(hn::Set(d, high), hn::Max(hn::Set(d, low), code));
      } else if (output_depth != 8 && output_depth != 16) {
        code = hn::Min(hn::Set(d, (1 << output_depth) - 1), code);
      }
      const hn::Rebind<std::uint16_t, decltype(d)> du16;
      if (output_depth == 8) {
        const hn::Rebind<std::uint8_t, decltype(d)> du8;
        hn::StoreU(hn::DemoteTo(du8, code), du8, row.packed_output + x);
      } else {
        const auto words = hn::DemoteTo(du16, code);
        const hn::Repartition<std::uint8_t, decltype(du16)> bytes;
        hn::StoreU(hn::BitCast(bytes, words), bytes, row.packed_output + 2 * x);
      }
      continue;
    }
    const hn::ScalableTag<double> dd;
    const hn::Half<decltype(d)> half;
    const auto low = hn::PromoteTo(dd, hn::LowerHalf(half, result));
    const auto high = hn::PromoteTo(dd, hn::UpperHalf(half, result));
    hn::StoreU(low, dd, row.output + x);
    hn::StoreU(high, dd, row.output + x + hn::Lanes(dd));
  }
  return full;
}
template <int SX, int SY, bool ValidatedDistances = false, bool CachedOffsets = false>
std::size_t integer_native_weighted(const IntegerRow& row, const double* thresholds) {
  // Keep output selection out of the general weighted loop. Zero denotes
  // the intermediate output used by floating conversion and error diffusion.
  if (!row.packed_output)
    return integer_samples<6, true, SX, SY, true, true, 0, ValidatedDistances, CachedOffsets>(row, thresholds);
  if (row.dither == 2)
    return integer_samples<6, true, SX, SY, true, true, 2, ValidatedDistances, CachedOffsets>(row, thresholds);
  return integer_samples<6, true, SX, SY, true, true, 1, ValidatedDistances, CachedOffsets>(row, thresholds);
}
std::size_t integer_cached_weighted(const IntegerRow& row, const double* thresholds) {
  // Cached offsets already encode subsampling. Share one layout-independent
  // kernel and specialize common storage pairs without duplicating layouts.
  if (row.packed_output && row.dither == 1 && row.input_depth == row.output_depth) {
    switch (row.input_depth) {
      case 8: return integer_samples<6, true, 0, 0, true, true, 1, true, true, 8, 8>(row, thresholds);
      case 10: return integer_samples<6, true, 0, 0, true, true, 1, true, true, 10, 10>(row, thresholds);
      case 16: return integer_samples<6, true, 0, 0, true, true, 1, true, true, 16, 16>(row, thresholds);
    }
  }
  if (row.packed_output && row.dither == 2 && row.input_depth == 16 && row.output_depth == 8)
    return integer_samples<6, true, 0, 0, true, true, 2, true, true, 16, 8>(row, thresholds);
  return integer_native_weighted<0, 0, true, true>(row, thresholds);
}
template <int SX, int SY>
std::size_t integer_dispatch(const IntegerRow& row, int mode, const double* thresholds, bool blur) {
  if (row.sub_x > 1 || row.width == 0)
    return 0;
#if HWY_TARGET == HWY_AVX3 || HWY_TARGET == HWY_AVX3_DL || HWY_TARGET == HWY_AVX3_ZEN4 || HWY_TARGET == HWY_AVX3_SPR || \
    HWY_TARGET == HWY_AVX2 || HWY_TARGET == HWY_SSE4 || HWY_TARGET == HWY_SSSE3
  // Extend compact-word arithmetic to narrow targets with measured gains.
  // Other output contracts and threshold ranges keep the general kernel.
  if (row.native && row.native_extent >= 4 && row.packed_output && row.dither == 1 &&
      thresholds[0] > 0 && thresholds[0] <= 65535) {
    // SSSE3 mode 1 keeps I32 after mixed/negative full-frame results.
    if constexpr (HWY_TARGET != HWY_SSSE3) {
      if (mode == 1)
        return blur ? native_word_pairs<1, true, SX, SY>(row, thresholds)
                    : native_word_pairs<1, false, SX, SY>(row, thresholds);
    }
    if (mode == 3)
      return blur ? native_word_pairs<3, true, SX, SY>(row, thresholds)
                  : native_word_pairs<3, false, SX, SY>(row, thresholds);
  }
#endif
  switch (mode) {
    case 0:
      return integer_samples<0, true, SX, SY>(row, thresholds);
    case 1:
      return blur ? integer_samples<1, true, SX, SY>(row, thresholds)
                  : integer_samples<1, false, SX, SY>(row, thresholds);
    case 2:
      return blur ? integer_samples<2, true, SX, SY>(row, thresholds)
                  : integer_samples<2, false, SX, SY>(row, thresholds);
    case 3:
      return blur ? integer_samples<3, true, SX, SY>(row, thresholds)
                  : integer_samples<3, false, SX, SY>(row, thresholds);
    case 4:
      return blur ? integer_samples<4, true, SX, SY>(row, thresholds)
                  : integer_samples<4, false, SX, SY>(row, thresholds);
    case 5:
      return integer_samples<5, true, SX, SY>(row, thresholds);
    case 6: {
      // Select the numerical contract once per row. Keeping the F64 fallback
      // out of the F32 loop reduces live registers and loop-invariant spills.
      bool fast_weights = true;
      for (int j = 0; j < 3; ++j)
        fast_weights &= thresholds[j] >= 1 && thresholds[j] <= 65535 &&
                        thresholds[j] == std::floor(thresholds[j]);
      if (row.native) {
        const std::size_t bytes = row.input_depth == 8 ? 1 : 2;
        const auto row_bytes = row.width * bytes;
        // Only the final (4/bytes - 1) pixels can cross the plane end.
        // Their horizontal limit also bounds the shared vertical distance:
        // d < (4/bytes - 1)*2^SX, hence dy <= (d >> SY). Outside these
        // final rows, all four-byte gathers fit even without trailing padding.
        const std::size_t last_rows = (((4 / bytes - 1) << SX) - 1) >> SY;
        if (row.native_extent >= row_bytes && row.row_offset <= row.native_extent - row_bytes &&
            row.native_extent - row_bytes - row.row_offset > last_rows * row.native_stride) {
          if (fast_weights && row.validated_distances && row.vertical_offsets && row.horizontal_offsets)
            return integer_cached_weighted(row, thresholds);
          return fast_weights ? (row.validated_distances ? integer_native_weighted<SX, SY, true>(row, thresholds)
                                                        : integer_native_weighted<SX, SY>(row, thresholds))
                              : integer_samples<6, true, SX, SY, true, false>(row, thresholds);
        }
      }
      return fast_weights ? integer_samples<6, true, SX, SY, false, true>(row, thresholds)
                          : integer_samples<6, true, SX, SY, false, false>(row, thresholds);
    }
    default:
      return 0;
  }
}
std::size_t integer_target(const IntegerRow& row, int mode, const double* thresholds, bool blur) {
  if (row.sub_x > 1 || row.sub_y > 1)
    return 0;
  if (row.sub_x == 0)
    return row.sub_y == 0 ? integer_dispatch<0, 0>(row, mode, thresholds, blur)
                          : integer_dispatch<0, 1>(row, mode, thresholds, blur);
  return row.sub_y == 0 ? integer_dispatch<1, 0>(row, mode, thresholds, blur)
                        : integer_dispatch<1, 1>(row, mode, thresholds, blur);
}
// x is a positive normal binary64 value in (0, 1]. Split x = m * 2^(10q+r),
// where m is in [1,2) and r in [0,9]. Only m needs a polynomial; powers of two
// are exact and the eleven-entry table joins adjacent exponent intervals.
HWY_INLINE V tenth_root_normal(D d, V x) {
  const hn::Rebind<std::int64_t, D> di;
  const auto bits = hn::BitCast(di, x);
  const auto exponent = hn::Sub(hn::ShiftRight<52>(bits), hn::Set(di, 1023));
#if HWY_TARGET == HWY_AVX2
  // Exact floor(e/10) for every normal exponent e in [-1022, 0].
  const hn::Rebind<std::int32_t, D> di32;
  const auto e32 = hn::DemoteTo(di32, exponent);
  const auto q32 = hn::Neg(hn::ShiftRight<15>(hn::Mul(hn::Sub(hn::Set(di32, 9), e32), hn::Set(di32, 3277))));
  const auto q = hn::PromoteTo(di, q32);
  const auto r = hn::PromoteTo(di, hn::Sub(e32, hn::Mul(q32, hn::Set(di32, 10))));
#else
  // There are only 1023 possible exponents. Multiplication by binary64 0.1
  // followed by floor gives floor(e/10) exactly throughout [-1022, 0].
  const auto q = hn::ConvertTo(di, hn::Floor(hn::Mul(hn::ConvertTo(d, exponent), hn::Set(d, 0.1))));
  const auto r = hn::Sub(exponent, hn::Mul(q, hn::Set(di, 10)));
#endif
#if HWY_TARGET == HWY_AVX3_SPR
  alignas(64) static constexpr double powers[16] = {
#else
  static constexpr double powers[11] = {
#endif
      0x1.0000000000000p+0, 0x1.125fbee250664p+0, 0x1.2611186bae675p+0,
      0x1.3b2c47bff8329p+0, 0x1.51cb453b9536cp+0, 0x1.6a09e667f3bcdp+0,
      0x1.8406003b2ae5cp+0, 0x1.9fdf8bcce533dp+0, 0x1.bdb8cdadbe120p+0,
      0x1.ddb680117ab12p+0, 0x1.0000000000000p+1};
#if HWY_TARGET == HWY_AVX3_SPR
  const auto table_lo = hn::Load(d, powers), table_hi = hn::Load(d, powers + 8);
  const auto lower = hn::TwoTablesLookupLanes(table_lo, table_hi, hn::IndicesFromVec(d, r));
  const auto upper = hn::TwoTablesLookupLanes(table_lo, table_hi,
      hn::IndicesFromVec(d, hn::Add(r, hn::Set(di, 1))));
#else
  const auto lower = hn::GatherIndex(d, powers, r);
  const auto upper = hn::GatherIndex(d, powers, hn::Add(r, hn::Set(di, 1)));
#endif
  const auto scale = hn::BitCast(d, hn::ShiftLeft<52>(hn::Add(q, hn::Set(di, 1023))));
  const auto mantissa = hn::BitCast(
      d, hn::Or(hn::And(bits, hn::Set(di, 0x000fffffffffffffLL)), hn::Set(di, 0x3ff0000000000000LL)));
  const auto t = hn::Sub(mantissa, hn::Set(d, 1.0));
  const auto z = hn::Sub(t, hn::Set(d, 0.5));
  // Endpoint-factored Chebyshev approximation on t in [0,1]:
  // (1+t)^0.1 = 1 + t * (2^0.1-1 + (1-t) * P(t-0.5)).
#if HWY_TARGET == HWY_AVX3_SPR
  // Degree-five endpoint-factored fit; max measured error < 7e-8.
  // A coefficient interval bound keeps the derivative positive.
  constexpr double coefficients[6] = {
      0x1.67fe7bd9c4501p-6, -0x1.33cb006754197p-7, 0x1.2ae1af7536712p-8,
      -0x1.37cbab919fdb6p-9, 0x1.85f25f1bcbafcp-10, -0x1.c3e926c09fdb9p-11};
  auto polynomial = hn::Set(d, coefficients[5]);
  for (int j = 4; j >= 0; --j)
    polynomial = hn::MulAdd(polynomial, z, hn::Set(d, coefficients[j]));
#else
  // Max measured error < 2.2e-10; a coefficient interval bound gives a
  // derivative > 0.039. Endpoint factoring keeps exponent joins continuous.
  constexpr double coefficients[9] = {
      0x1.67fd762703c62p-6, -0x1.33c6958706f00p-7, 0x1.2c02b1c129546p-8,
      -0x1.3a011bb834600p-9, 0x1.57944156155c1p-10, -0x1.7d3b21ee12484p-11,
      0x1.b4f8b68efa72ep-12, -0x1.431b83c280ea9p-12, 0x1.8b15e6f1fa946p-13};
  auto polynomial = hn::Set(d, coefficients[8]);
  for (int j = 7; j >= 0; --j)
    polynomial = hn::MulAdd(polynomial, z, hn::Set(d, coefficients[j]));
#endif
  const auto correction = hn::Mul(hn::Sub(hn::Set(d, 1.0), t), polynomial);
  const auto root = hn::MulAdd(t, hn::MulAdd(lower, correction, hn::Sub(upper, lower)), lower);
  return hn::Mul(scale, hn::Min(upper, hn::Max(lower, root)));
}

HWY_INLINE V weight(D d, const V* c) {
  auto all_positive = hn::Gt(c[0], hn::Zero(d));
  auto all_one = hn::Eq(c[0], hn::Set(d, 1.0));
  for (int j = 1; j < 4; ++j) {
    all_positive = hn::And(all_positive, hn::Gt(c[j], hn::Zero(d)));
    all_one = hn::And(all_one, hn::Eq(c[j], hn::Set(d, 1.0)));
  }
  // Each lane can independently be an exact zero or one endpoint.
  if (hn::AllTrue(d, hn::Or(hn::Not(all_positive), all_one)))
    return hn::IfThenElseZero(all_positive, hn::Set(d, 1.0));
  // Ordinary products use a bounded tenth-root approximation. A normal final
  // product guarantees normal intermediate products because every response is
  // in [0, 1]. Preserve the logarithmic sum below for underflow/subnormal cases.
  const auto product = hn::Mul(hn::Mul(c[0], c[1]), hn::Mul(c[2], c[3]));
  if (hn::AllTrue(d, hn::Or(hn::Not(all_positive),
                           hn::Ge(product, hn::Set(d, std::numeric_limits<double>::min()))))) {
    const auto argument = hn::IfThenElse(all_positive, product, hn::Set(d, 1.0));
    const auto result = hn::Min(hn::Set(d, 1.0), tenth_root_normal(d, argument));
    return hn::IfThenElseZero(all_positive, result);
  }
  auto sum = hn::Zero(d);
  auto positive = hn::Gt(c[0], hn::Zero(d));
  for (int j = 0; j < 4; ++j) {
    positive = hn::And(positive, hn::Gt(c[j], hn::Zero(d)));
    // Avoid log(0), including inactive tail lanes. Exact endpoints are explicit.
    const auto value = hn::IfThenElse(hn::Gt(c[j], hn::Zero(d)), c[j], hn::Set(d, 1.0));
    // Normalize subnormal responses before Highway Log, whose polynomial
    // expects a normal argument. The logarithm remains finite even when the
    // product of four responses would underflow.
    const auto tiny = hn::Lt(value, hn::Set(d, std::numeric_limits<double>::min()));
    const auto argument = hn::IfThenElse(tiny, hn::Mul(value, hn::Set(d, 0x1p52)), value);
    const auto logarithm =
        hn::Sub(hn::Log(d, argument), hn::IfThenElseZero(tiny, hn::Set(d, 52 * 0.693147180559945309417232121458)));
    sum = hn::Add(sum, logarithm);
  }
  return hn::IfThenElseZero(positive, hn::Min(hn::Set(d, 1.0), hn::Exp(d, hn::Mul(sum, hn::Set(d, 0.1)))));
}
HWY_INLINE V response(D d, V delta, V threshold) {
  const auto zero = hn::Zero(d), one = hn::Set(d, 1.0);
  const auto iszero = hn::Eq(threshold, zero);
  const auto ratio = hn::Div(delta, hn::IfThenElse(iszero, one, threshold));
  return hn::IfThenElse(iszero, hn::IfThenElse(hn::Eq(delta, zero), one, zero),
                        hn::Max(zero, hn::Min(one, hn::Mul(hn::Set(d, 3.0), hn::Sub(one, ratio)))));
}
HWY_INLINE V response_with_inverse(D d, V delta, V threshold, V inverse) {
  const auto zero = hn::Zero(d), one = hn::Set(d, 1.0);
  if (!hn::AllTrue(d, hn::Le(hn::Abs(inverse), hn::Set(d, std::numeric_limits<double>::max()))))
    return response(d, delta, threshold);
  const auto iszero = hn::Eq(threshold, zero);
  auto ratio = hn::Mul(delta, inverse);
  // Preserve division at the zero-response boundary, where a tiny absolute
  // ratio error can become a visible difference after the tenth root.
  const auto boundary = hn::Lt(hn::Abs(hn::Sub(one, ratio)), hn::Set(d, 0x1p-40));
  if (!hn::AllFalse(d, boundary))
    ratio = hn::IfThenElse(boundary, hn::Div(delta, hn::IfThenElse(iszero, one, threshold)), ratio);
  return hn::IfThenElse(iszero, hn::IfThenElse(hn::Eq(delta, zero), one, zero),
                        hn::Max(zero, hn::Min(one, hn::Mul(hn::Set(d, 3.0), hn::Sub(one, ratio)))));
}
HWY_INLINE M aligned(D d, V gx, V gy, V hx, V hy, double angle, V cosine, V sine) {
  const auto zero = hn::Zero(d), one = hn::Set(d, 1.0);
  const auto gm = hn::Max(hn::Abs(gx), hn::Abs(gy)), hm = hn::Max(hn::Abs(hx), hn::Abs(hy));
  const auto bothzero = hn::And(hn::Eq(gm, zero), hn::Eq(hm, zero));
  const auto nonzero = hn::And(hn::Gt(gm, zero), hn::Gt(hm, zero));
  if (angle >= 0.5)
    return hn::Or(bothzero, nonzero);
  auto integral = hn::And(hn::Le(gm, hn::Set(d, 262140.0)), hn::Le(hm, hn::Set(d, 262140.0)));
  integral = hn::And(integral, hn::And(hn::Eq(hn::Trunc(gx), gx), hn::Eq(hn::Trunc(gy), gy)));
  integral = hn::And(integral, hn::And(hn::Eq(hn::Trunc(hx), hx), hn::Eq(hn::Trunc(hy), hy)));
  // These bounded integer products and sums are exact in double. Their
  // normalization factors are one, so the four divisions have no effect.
  if (!hn::AllTrue(d, integral)) {
    const auto gs = hn::IfThenElse(hn::Or(integral, hn::Eq(gm, zero)), one, gm);
    const auto hs = hn::IfThenElse(hn::Or(integral, hn::Eq(hm, zero)), one, hm);
    gx = hn::Div(gx, gs);
    gy = hn::Div(gy, gs);
    hx = hn::Div(hx, hs);
    hy = hn::Div(hy, hs);
  }
  const auto cross = hn::Abs(hn::Sub(hn::Mul(gx, hy), hn::Mul(gy, hx)));
  const auto dot = hn::Abs(hn::Add(hn::Mul(gx, hx), hn::Mul(gy, hy)));
  auto accept = hn::Eq(cross, zero);
  if (angle == 0.25)
    accept = hn::Le(cross, dot);
  else if (angle != 0) {
    accept = hn::Le(hn::Mul(cross, cosine), hn::Mul(dot, sine));
  }
  return hn::Or(bothzero, hn::And(nonzero, accept));
}
HWY_INLINE M aligned_prepared(D d, const DirectionCenter& g, V hx, V hy, double angle, V cosine, V sine) {
  const auto zero = hn::Zero(d), one = hn::Set(d, 1.0);
  const auto gm = g.magnitude, hm = hn::Max(hn::Abs(hx), hn::Abs(hy));
  auto gx = g.x, gy = g.y;
  const auto bothzero = hn::And(hn::Eq(gm, zero), hn::Eq(hm, zero));
  const auto nonzero = hn::And(hn::Gt(gm, zero), hn::Gt(hm, zero));
  if (angle >= 0.5)
    return hn::Or(bothzero, nonzero);
  auto integral = hn::And(g.integral, hn::Le(hm, hn::Set(d, 262140.0)));
  integral = hn::And(integral, hn::And(hn::Eq(hn::Trunc(hx), hx), hn::Eq(hn::Trunc(hy), hy)));
  // These bounded integer products and sums are exact in double. Their
  // normalization factors are one, so the four divisions have no effect.
  if (!hn::AllTrue(d, integral)) {
    const auto hs = hn::IfThenElse(hn::Or(integral, hn::Eq(hm, zero)), one, hm);
    gx = hn::IfThenElse(integral, gx, g.normalized_x);
    gy = hn::IfThenElse(integral, gy, g.normalized_y);
    hx = hn::Div(hx, hs);
    hy = hn::Div(hy, hs);
  }
  const auto cross = hn::Abs(hn::Sub(hn::Mul(gx, hy), hn::Mul(gy, hx)));
  const auto dot = hn::Abs(hn::Add(hn::Mul(gx, hx), hn::Mul(gy, hy)));
  auto accept = hn::Eq(cross, zero);
  if (angle == 0.25)
    accept = hn::Le(cross, dot);
  else if (angle != 0) {
    accept = hn::Le(hn::Mul(cross, cosine), hn::Mul(dot, sine));
  }
  return hn::Or(bothzero, hn::And(nonzero, accept));
}
// Modes 1/3 use two references and mode 2 uses four. Resolve these uniform
// choices before the pixel loop so the common path has fixed arithmetic and
// only loads references that participate in the mean.
template <int N, bool Integer, bool Blur>
void simple_row(const Row& row, double threshold) {
  const D d;
  const auto limit = hn::Set(d, threshold);
  const auto zero = hn::Zero(d);
  for (std::size_t x = 0; x < row.count; x += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), row.count - x);
    const auto center = load_row(d, row.center + x, n);
    V refs[N];
    for (int j = 0; j < N; ++j)
      refs[j] = load_row(d, row.refs[j] + x, n);
    const auto average = mean(d, refs, N, Integer);
    auto accepted = hn::Lt(hn::Abs(diff(d, average, center, Integer)), limit);
    if constexpr (!Blur) {
      accepted = hn::Eq(zero, zero);
      for (int j = 0; j < N; ++j)
        accepted = hn::And(accepted, hn::Lt(hn::Abs(diff(d, refs[j], center, Integer)), limit));
    }
    auto result =
        hn::Add(hn::IfThenElse(accepted, average, center),
                (Integer ? load_integer_grain(d, row.integer_grain + x, n) : load_grain(d, row.grain + x, n)));
    if constexpr (!Integer)
      result = fp(d, result);
    store_row(d, result, row.output + x, n);
  }
}
template <int N>
void simple_row_dispatch(const Row& row, double threshold, bool integer, bool blur) {
  if (integer) {
    if (blur)
      simple_row<N, true, true>(row, threshold);
    else
      simple_row<N, true, false>(row, threshold);
  } else {
    if (blur)
      simple_row<N, false, true>(row, threshold);
    else
      simple_row<N, false, false>(row, threshold);
  }
}
#if HWY_TARGET == HWY_AVX2
HWY_INLINE V row_pass(D d, V s, const V* values, int n, double threshold, bool integer, bool blur) {
  const auto a = mean(d, values, n, integer), zero = hn::Zero(d);
  auto accepted = hn::Lt(hn::Abs(diff(d, a, s, integer)), hn::Set(d, threshold));
  if (!blur) {
    accepted = hn::Eq(zero, zero);
    for (int j = 0; j < n; ++j)
      accepted = hn::And(accepted, hn::Lt(hn::Abs(diff(d, values[j], s, integer)), hn::Set(d, threshold)));
  }
  return hn::IfThenElse(accepted, a, s);
}
HWY_INLINE V row_pair(D d, V s, V v, V w, bool integer) {
  auto sum = hn::Add(v, w), twice = hn::Add(s, s);
  if (!integer) {
    sum = fp(d, sum);
    twice = fp(d, twice);
  }
  return hn::Abs(diff(d, sum, twice, integer));
}
template <int Mode, bool Integer>
void row_specialized(const Row& row, const double* thresholds, bool blur, double angle,
                     double boost) {
  constexpr bool integer = Integer;
  constexpr int mode = Mode;
  if (mode <= 3) {
    if (mode == 2)
      simple_row_dispatch<4>(row, thresholds[0], integer, blur);
    else
      simple_row_dispatch<2>(row, thresholds[0], integer, blur);
    return;
  }
  double inverse_thresholds[3]{};
  if (!integer && mode == 6)
    for (int j = 0; j < 3; ++j)
      inverse_thresholds[j] = thresholds[j] == 0 ? 1 : 1 / thresholds[j];
  const D d;
  const auto zero = hn::Zero(d), half = hn::Set(d, 0.5), one = hn::Set(d, 1.0);
  // The direction limit is uniform across the row and all four neighbours.
  constexpr double pi = 3.1415926535897932384626433832795;
  const bool trig = mode == 7 && angle != 0 && angle != 0.25 && angle < 0.5;
  const auto cosine = hn::Set(d, trig ? std::cos(pi * angle) : 0);
  const auto sine = hn::Set(d, trig ? std::sin(pi * angle) : 0);
  for (std::size_t x = 0; x < row.count; x += hn::Lanes(d)) {
    const auto count = std::min(hn::Lanes(d), row.count - x);
    const auto s = load_row(d, row.center + x, count);
    V r[4];
    for (int j = 0; j < 4; ++j)
      r[j] = load_row(d, row.refs[j] + x, count);
    auto result = s;
    if (mode == 4) {
      const V results[2] = {row_pass(d, s, r, 2, thresholds[0], integer, blur), row_pass(d, s, r + 2, 2, thresholds[0], integer, blur)};
      result = mean(d, results, 2, integer);
    } else {
      const auto a = mean(d, r, 4, mode == 5 && integer);
      const auto da = hn::Abs(diff(d, a, s, integer));
      auto maximum = zero;
      for (int j = 0; j < 4; ++j)
        maximum = hn::Max(maximum, hn::Abs(diff(d, r[j], s, integer)));
      const auto qv = row_pair(d, s, r[0], r[1], integer), qh = row_pair(d, s, r[2], r[3], integer);
      if (mode == 5) {
        auto accepted = hn::And(hn::Lt(da, hn::Set(d, thresholds[0])), hn::Lt(maximum, hn::Set(d, thresholds[1])));
        accepted =
            hn::And(accepted, hn::And(hn::Lt(qv, hn::Set(d, thresholds[2])), hn::Lt(qh, hn::Set(d, thresholds[2]))));
        result = hn::IfThenElse(accepted, a, s);
      } else {
        auto k = one;
        if (mode == 7) {
          auto directions = hn::Eq(zero, zero);
          const auto gx = load_row(d, row.gx[0] + x, count), gy = load_row(d, row.gy[0] + x, count);
          for (int j = 1; j < 5; ++j)
            directions = hn::And(directions, aligned(d, gx, gy, load_row(d, row.gx[j] + x, count),
                                                     load_row(d, row.gy[j] + x, count), angle, cosine, sine));
          k = hn::IfThenElse(directions, hn::Set(d, boost), one);
        }
        V c[4];
        if (!integer && mode == 6) {
          c[0] = response_with_inverse(d, da, hn::Set(d, thresholds[0]), hn::Set(d, inverse_thresholds[0]));
          c[1] = response_with_inverse(d, maximum, hn::Set(d, thresholds[1]), hn::Set(d, inverse_thresholds[1]));
          c[2] = response_with_inverse(d, qv, hn::Set(d, thresholds[2]), hn::Set(d, inverse_thresholds[2]));
          c[3] = response_with_inverse(d, qh, hn::Set(d, thresholds[2]), hn::Set(d, inverse_thresholds[2]));
        } else {
          c[0] = response(d, da, hn::Mul(k, hn::Set(d, thresholds[0])));
          c[1] = response(d, maximum, hn::Mul(k, hn::Set(d, thresholds[1])));
          c[2] = response(d, qv, hn::Mul(k, hn::Set(d, thresholds[2])));
          c[3] = response(d, qh, hn::Mul(k, hn::Set(d, thresholds[2])));
        }
        const auto f = weight(d, c);
        // FMA in the integer blend is permitted; final integer error is tested.
        const auto blend = integer ? hn::Floor(hn::Add(hn::MulAdd(hn::Sub(a, s), f, s), half))
                                   : fp(d, hn::Add(s, fp(d, hn::Mul(diff(d, a, s, false), f))));
        result = hn::IfThenElse(hn::Eq(f, zero), s, blend);
      }
    }
    result = hn::Add(
        result, (integer ? load_integer_grain(d, row.integer_grain + x, count) : load_grain(d, row.grain + x, count)));
    if (!integer)
      result = fp(d, result);
    store_row(d, result, row.output + x, count);
  }
}
#endif
void row_target(const Row& row, int mode, const double* thresholds, bool integer, bool blur, double angle,
                double boost) {
#if HWY_TARGET == HWY_AVX2
  if (!integer && (mode == 5 || mode == 6)) {
    if (mode == 5) row_specialized<5, false>(row, thresholds, blur, angle, boost);
    else row_specialized<6, false>(row, thresholds, blur, angle, boost);
    return;
  }
#endif

  if (mode <= 3) {
    if (mode == 2)
      simple_row_dispatch<4>(row, thresholds[0], integer, blur);
    else
      simple_row_dispatch<2>(row, thresholds[0], integer, blur);
    return;
  }
  double inverse_thresholds[3]{};
  if (!integer && mode == 6)
    for (int j = 0; j < 3; ++j)
      inverse_thresholds[j] = thresholds[j] == 0 ? 1 : 1 / thresholds[j];
  const D d;
  const auto zero = hn::Zero(d), half = hn::Set(d, 0.5), one = hn::Set(d, 1.0);
  // The direction limit is uniform across the row and all four neighbours.
  constexpr double pi = 3.1415926535897932384626433832795;
  const bool trig = mode == 7 && angle != 0 && angle != 0.25 && angle < 0.5;
  const auto cosine = hn::Set(d, trig ? std::cos(pi * angle) : 0);
  const auto sine = hn::Set(d, trig ? std::sin(pi * angle) : 0);
  for (std::size_t x = 0; x < row.count; x += hn::Lanes(d)) {
    const auto count = std::min(hn::Lanes(d), row.count - x);
    const auto s = load_row(d, row.center + x, count);
    V r[4];
    for (int j = 0; j < 4; ++j)
      r[j] = load_row(d, row.refs[j] + x, count);
    auto pass = [&](const V* values, int n) {
      const auto a = mean(d, values, n, integer);
      auto accepted = hn::Lt(hn::Abs(diff(d, a, s, integer)), hn::Set(d, thresholds[0]));
      if (!blur) {
        accepted = hn::Eq(zero, zero);
        for (int j = 0; j < n; ++j)
          accepted = hn::And(accepted, hn::Lt(hn::Abs(diff(d, values[j], s, integer)), hn::Set(d, thresholds[0])));
      }
      return hn::IfThenElse(accepted, a, s);
    };
    auto result = s;
    if (mode == 4) {
      const V results[2] = {pass(r, 2), pass(r + 2, 2)};
      result = mean(d, results, 2, integer);
    } else {
      const auto a = mean(d, r, 4, mode == 5 && integer);
      const auto da = hn::Abs(diff(d, a, s, integer));
      auto maximum = zero;
      for (int j = 0; j < 4; ++j)
        maximum = hn::Max(maximum, hn::Abs(diff(d, r[j], s, integer)));
      auto pair = [&](V v, V w) {
        auto sum = hn::Add(v, w), twice = hn::Add(s, s);
        if (!integer) {
          sum = fp(d, sum);
          twice = fp(d, twice);
        }
        return hn::Abs(diff(d, sum, twice, integer));
      };
      const auto qv = pair(r[0], r[1]), qh = pair(r[2], r[3]);
      if (mode == 5) {
        auto accepted = hn::And(hn::Lt(da, hn::Set(d, thresholds[0])), hn::Lt(maximum, hn::Set(d, thresholds[1])));
        accepted =
            hn::And(accepted, hn::And(hn::Lt(qv, hn::Set(d, thresholds[2])), hn::Lt(qh, hn::Set(d, thresholds[2]))));
        result = hn::IfThenElse(accepted, a, s);
      } else {
        auto k = one;
        if (mode == 7) {
          auto directions = hn::Eq(zero, zero);
          const auto gx = load_row(d, row.gx[0] + x, count), gy = load_row(d, row.gy[0] + x, count);
          for (int j = 1; j < 5; ++j)
            directions = hn::And(directions, aligned(d, gx, gy, load_row(d, row.gx[j] + x, count),
                                                     load_row(d, row.gy[j] + x, count), angle, cosine, sine));
          k = hn::IfThenElse(directions, hn::Set(d, boost), one);
        }
        V c[4];
        if (!integer && mode == 6) {
          c[0] = response_with_inverse(d, da, hn::Set(d, thresholds[0]), hn::Set(d, inverse_thresholds[0]));
          c[1] = response_with_inverse(d, maximum, hn::Set(d, thresholds[1]), hn::Set(d, inverse_thresholds[1]));
          c[2] = response_with_inverse(d, qv, hn::Set(d, thresholds[2]), hn::Set(d, inverse_thresholds[2]));
          c[3] = response_with_inverse(d, qh, hn::Set(d, thresholds[2]), hn::Set(d, inverse_thresholds[2]));
        } else {
          c[0] = response(d, da, hn::Mul(k, hn::Set(d, thresholds[0])));
          c[1] = response(d, maximum, hn::Mul(k, hn::Set(d, thresholds[1])));
          c[2] = response(d, qv, hn::Mul(k, hn::Set(d, thresholds[2])));
          c[3] = response(d, qh, hn::Mul(k, hn::Set(d, thresholds[2])));
        }
        const auto f = weight(d, c);
        // FMA in the integer blend is permitted; final integer error is tested.
        const auto blend = integer ? hn::Floor(hn::Add(hn::MulAdd(hn::Sub(a, s), f, s), half))
                                   : fp(d, hn::Add(s, fp(d, hn::Mul(diff(d, a, s, false), f))));
        result = hn::IfThenElse(hn::Eq(f, zero), s, blend);
      }
    }
    result = hn::Add(
        result, (integer ? load_integer_grain(d, row.integer_grain + x, count) : load_grain(d, row.grain + x, count)));
    if (!integer)
      result = fp(d, result);
    store_row(d, result, row.output + x, count);
  }
}
void weights_target(const double* c[4], double* output, std::size_t count) {
  const D d;
  for (std::size_t x = 0; x < count; x += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), count - x);
    const V values[4] = {load_row(d, c[0] + x, n), load_row(d, c[1] + x, n), load_row(d, c[2] + x, n),
                         load_row(d, c[3] + x, n)};
    store_row(d, weight(d, values), output + x, n);
  }
}
#endif
void frame_into_target(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input,
                       const std::vector<MutablePlaneView>& output) {
#if HWY_TARGET == HWY_SCALAR || HWY_TARGET == HWY_EMU128
  const auto result = process_scalar(filter, frame, input);
  for (std::size_t plane = 0; plane < result.planes.size(); ++plane) {
    const auto& p = result.planes[plane];
    for (std::size_t y = 0; y < p.height; ++y)
      std::memcpy(static_cast<std::uint8_t*>(output[plane].data) + y * output[plane].stride,
                  p.data.data() + y * p.stride, p.stride);
  }
#else
  process_highway_into(filter, frame, input, output);
#endif
}
Frame frame_target(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input) {
#if HWY_TARGET == HWY_SCALAR || HWY_TARGET == HWY_EMU128
  return process_scalar(filter, frame, input);
#else
  return process_highway(filter, frame, input);
#endif
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo::f3kdb
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo::f3kdb {
HWY_EXPORT(frame_into_target);
void process_frame_into(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input,
                        const std::vector<MutablePlaneView>& output) {
  HWY_DYNAMIC_DISPATCH(frame_into_target)(filter, frame, input, output);
}
HWY_EXPORT(frame_target);
Frame process_frame(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input) {
  return HWY_DYNAMIC_DISPATCH(frame_target)(filter, frame, input);
}
HWY_EXPORT(cross_sampled_target);
SamplingResult process_cross_row(const CrossSamplingRow& row, const Row& output, const double* t, bool integer,
                                 double angle, double boost) {
  return HWY_DYNAMIC_DISPATCH(cross_sampled_target)(row, output, t, integer, angle, boost);
}
HWY_EXPORT(collect_cross_target);
std::size_t collect_cross_row(const CrossSamplingRow& row) {
  return HWY_DYNAMIC_DISPATCH(collect_cross_target)(row);
}
HWY_EXPORT(byte_target);
std::size_t process_byte_row(const ByteRow& row, int mode, const double* thresholds, bool blur) {
  return HWY_DYNAMIC_DISPATCH(byte_target)(row, mode, thresholds, blur);
}
HWY_EXPORT(integer_target);
std::size_t process_integer_row(const IntegerRow& row, int mode, const double* thresholds, bool blur) {
  return HWY_DYNAMIC_DISPATCH(integer_target)(row, mode, thresholds, blur);
}
HWY_EXPORT(sampled_target);
SamplingResult process_sampled_row(const SamplingRow& sampling, const Row& row, int mode, double threshold,
                                   bool integer, bool blur) {
  return HWY_DYNAMIC_DISPATCH(sampled_target)(sampling, row, mode, threshold, integer, blur);
}
HWY_EXPORT(collect_simple_target);
std::size_t collect_simple_row(const SamplingRow& row, int mode) {
  return HWY_DYNAMIC_DISPATCH(collect_simple_target)(row, mode);
}
HWY_EXPORT(row_target);
HWY_EXPORT(weights_target);
void process_row(const Row& row, int mode, const double* t, bool integer, bool blur, double angle, double boost) {
  HWY_DYNAMIC_DISPATCH(row_target)(row, mode, t, integer, blur, angle, boost);
}
void weights(const double* c[4], double* output, std::size_t n) {
  HWY_DYNAMIC_DISPATCH(weights_target)(c, output, n);
}
} // namespace neo_smo::f3kdb
#endif
