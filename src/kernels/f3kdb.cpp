#include "kernels/f3kdb.hpp"
#include "kernels/f3kdb_scalar.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
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
HWY_INLINE V weight(D d, const V* c) {
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
HWY_INLINE M aligned(D d, V gx, V gy, V hx, V hy, double angle) {
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
    constexpr double pi = 3.1415926535897932384626433832795;
    accept = hn::Le(hn::Mul(cross, hn::Set(d, std::cos(pi * angle))), hn::Mul(dot, hn::Set(d, std::sin(pi * angle))));
  }
  return hn::Or(bothzero, hn::And(nonzero, accept));
}
void row_target(const Row& row, int mode, const double* thresholds, bool integer, bool blur, double angle,
                double boost) {
  const D d;
  const auto zero = hn::Zero(d), half = hn::Set(d, 0.5), one = hn::Set(d, 1.0);
  for (std::size_t x = 0; x < row.count; x += hn::Lanes(d)) {
    const auto count = std::min(hn::Lanes(d), row.count - x);
    const auto s = hn::LoadN(d, row.center + x, count);
    V r[4];
    for (int j = 0; j < 4; ++j)
      r[j] = hn::LoadN(d, row.refs[j] + x, count);
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
    if (mode <= 3)
      result = pass(r, mode == 2 ? 4 : 2);
    else if (mode == 4) {
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
          const auto gx = hn::LoadN(d, row.gx[0] + x, count), gy = hn::LoadN(d, row.gy[0] + x, count);
          for (int j = 1; j < 5; ++j)
            directions = hn::And(directions, aligned(d, gx, gy, hn::LoadN(d, row.gx[j] + x, count),
                                                     hn::LoadN(d, row.gy[j] + x, count), angle));
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
    result = hn::Add(result, hn::LoadN(d, row.grain + x, count));
    if (!integer)
      result = fp(d, result);
    hn::StoreN(result, d, row.output + x, count);
  }
}
void weights_target(const double* c[4], double* output, std::size_t count) {
  const D d;
  for (std::size_t x = 0; x < count; x += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), count - x);
    const V values[4] = {hn::LoadN(d, c[0] + x, n), hn::LoadN(d, c[1] + x, n), hn::LoadN(d, c[2] + x, n),
                         hn::LoadN(d, c[3] + x, n)};
    hn::StoreN(weight(d, values), d, output + x, n);
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
