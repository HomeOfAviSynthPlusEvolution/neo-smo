#include "kernels/f3kdb.hpp"
#include "kernels/f3kdb_scalar.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
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
std::size_t collect_simple_target(const SamplingRow&, int) {
  return 0;
}
SamplingResult sampled_target(const SamplingRow&, const Row&, int, double, bool, bool) {
  return {0, 0};
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
template <int N, bool Blur>
HWY_INLINE V sampled_value(D d, V center, const V* refs, double threshold);
// Read the interleaved distance bytes directly. Chroma subsampling by two
// skips every other pair; wider subsampling uses the existing scalar collector.
template <int Mode, bool Fused = false, bool Blur = true>
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
      auto result = sampled_value<Mode == 2 ? 4 : 2, Blur>(d, hn::LoadU(d, center), refs, threshold);
      result = hn::Add(result, load_integer_grain(d, output->integer_grain + x, hn::Lanes(d)));
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
template <int N, bool Blur>
HWY_INLINE V sampled_value(D d, V center, const V* refs, double threshold) {
  const auto average = mean(d, refs, N, true);
  const auto limit = hn::Set(d, threshold);
  auto accept = hn::Lt(hn::Abs(diff(d, average, center, true)), limit);
  if constexpr (!Blur) {
    accept = hn::Eq(center, center);
    for (int j = 0; j < N; ++j)
      accept = hn::And(accept, hn::Lt(hn::Abs(diff(d, refs[j], center, true)), limit));
  }
  return hn::IfThenElse(accept, average, center);
}
SamplingResult sampled_target(const SamplingRow& sampling, const Row& row, int mode, double threshold, bool integer,
                              bool blur) {
  // Four-reference fusion is slower than collection on four-double-lane targets.
  if (!integer || hn::Lanes(D{}) < 4 || sampling.sub_x > 1 || sampling.width == 0 || (hn::Lanes(D{}) == 4 && mode == 2))
    return {collect_simple_target(sampling, mode), 0};
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
template <int Mode, bool Blur, int SX, int SY>
std::size_t integer_samples(const IntegerRow& row, const double* thresholds) {
  const hn::ScalableTag<std::int32_t> d;
  const hn::Rebind<std::int32_t, decltype(d)> di;
  const hn::Rebind<std::uint8_t, decltype(d)> db;
  const auto lanes = hn::Lanes(d), full = row.width / lanes * lanes;
  const auto* distances = reinterpret_cast<const std::uint8_t*>(row.distances);
  const auto step = std::size_t{1} << SX;
  const auto index = hn::Iota(di, 0), last = hn::Set(di, std::int32_t(row.width - 1)),
             vlimit = hn::Set(di, std::int32_t(row.vertical_limit));
  const auto pitch = hn::Set(di, std::int32_t(row.width));
  // Integer differences satisfy delta < t exactly when delta < ceil(t).
  const auto limit = hn::Set(d, std::int32_t(std::ceil(thresholds[0])));
  for (std::size_t x = 0; x < full; x += lanes, distances += 2 * step * lanes) {
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
    if (!hn::AllTrue(di, valid))
      throw std::logic_error("reference outside plane");
    const auto* center = row.center + x;
    const auto* source = row.center - row.row_offset;
    const auto absolute = hn::Add(pos, hn::Set(di, std::int32_t(row.row_offset)));
    constexpr int N = Mode == 1 || Mode == 3 ? 2 : 4;
    hn::VFromD<decltype(d)> refs[N];
    if constexpr (Mode == 1 || Mode == 3) {
      const auto off = Mode == 1 ? hn::Mul(v1, pitch) : h1;
      refs[0] = hn::GatherIndex(d, source, hn::Add(absolute, off));
      refs[1] = hn::GatherIndex(d, source, hn::Sub(absolute, off));
    } else if constexpr (Mode == 2) {
      const auto distance2 = hn::PromoteTo(di, second), h2 = hn::ShiftRightSame(distance2, SX),
                 v2 = hn::ShiftRightSame(distance2, SY);
      if (!hn::AllTrue(di, hn::And(hn::Le(h2, hlimit), hn::Le(v2, vlimit))))
        throw std::logic_error("reference outside plane");
      const auto a = hn::Add(hn::Mul(v2, pitch), h1), b = hn::Sub(h2, hn::Mul(v1, pitch));
      refs[0] = hn::GatherIndex(d, source, hn::Add(absolute, a));
      refs[1] = hn::GatherIndex(d, source, hn::Add(absolute, b));
      refs[2] = hn::GatherIndex(d, source, hn::Sub(absolute, a));
      refs[3] = hn::GatherIndex(d, source, hn::Sub(absolute, b));
    }
    if constexpr (Mode >= 4) {
      const auto vertical = hn::Mul(v1, pitch);
      refs[0] = hn::GatherIndex(d, source, hn::Add(absolute, vertical));
      refs[1] = hn::GatherIndex(d, source, hn::Sub(absolute, vertical));
      refs[2] = hn::GatherIndex(d, source, hn::Add(absolute, h1));
      refs[3] = hn::GatherIndex(d, source, hn::Sub(absolute, h1));
    }
    auto sum = hn::Zero(d);
    for (int j = 0; j < N; ++j)
      sum = hn::Add(sum, refs[j]);
    auto mean = hn::ShiftRight<(N == 4 ? 2 : 1)>(hn::Add(sum, hn::Set(d, N / 2)));
    const auto value = hn::LoadU(d, center);
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
      const auto maximum_limit = hn::Set(d, std::int32_t(std::ceil(thresholds[1])));
      for (int j = 0; j < 4; ++j)
        accept = hn::And(accept, hn::Lt(hn::Abs(hn::Sub(refs[j], value)), maximum_limit));
      const auto pair_limit = hn::Set(d, std::int32_t(std::ceil(thresholds[2])));
      const auto twice = hn::Add(value, value);
      for (int j : {0, 2})
        accept = hn::And(accept, hn::Lt(hn::Abs(hn::Sub(hn::Add(refs[j], refs[j + 1]), twice)), pair_limit));
    }
    const hn::Rebind<std::int16_t, decltype(d)> dg;
    const auto result = hn::Add(hn::IfThenElse(accept, mean, value), hn::PromoteTo(d, hn::LoadU(dg, row.grain + x)));
    const hn::ScalableTag<double> dd;
    const hn::Half<decltype(d)> half;
    const auto low = hn::PromoteTo(dd, hn::LowerHalf(half, result));
    const auto high = hn::PromoteTo(dd, hn::UpperHalf(half, result));
    hn::StoreU(low, dd, row.output + x);
    hn::StoreU(high, dd, row.output + x + hn::Lanes(dd));
  }
  return full;
}
template <int SX, int SY>
std::size_t integer_dispatch(const IntegerRow& row, int mode, const double* thresholds, bool blur) {
  if (row.sub_x > 1 || row.width == 0)
    return 0;
  switch (mode) {
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
  const auto gs = hn::IfThenElse(hn::Or(integral, hn::Eq(gm, zero)), one, gm);
  const auto hs = hn::IfThenElse(hn::Or(integral, hn::Eq(hm, zero)), one, hm);
  gx = hn::Div(gx, gs);
  gy = hn::Div(gy, gs);
  hx = hn::Div(hx, hs);
  hy = hn::Div(hy, hs);
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
void row_target(const Row& row, int mode, const double* thresholds, bool integer, bool blur, double angle,
                double boost) {
  if (mode <= 3) {
    if (mode == 2)
      simple_row_dispatch<4>(row, thresholds[0], integer, blur);
    else
      simple_row_dispatch<2>(row, thresholds[0], integer, blur);
    return;
  }
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
        const V c[4] = {response(d, da, hn::Mul(k, hn::Set(d, thresholds[0]))),
                        response(d, maximum, hn::Mul(k, hn::Set(d, thresholds[1]))),
                        response(d, qv, hn::Mul(k, hn::Set(d, thresholds[2]))),
                        response(d, qh, hn::Mul(k, hn::Set(d, thresholds[2])))};
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
HWY_EXPORT(frame_target);
Frame process_frame(const Filter& filter, std::size_t frame, const std::vector<PlaneView>& input) {
  return HWY_DYNAMIC_DISPATCH(frame_target)(filter, frame, input);
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
