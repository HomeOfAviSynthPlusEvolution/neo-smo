#include "algorithms/f3kdb_views.hpp"
#include "kernels/f3kdb.hpp"
#include "kernels/f3kdb_scalar.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <type_traits>
#include <vector>

#ifndef HWY_COMPILE_ALL_ATTAINABLE
#define HWY_COMPILE_ALL_ATTAINABLE
#endif
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/f3kdb_fs.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo::f3kdb {
namespace HWY_NAMESPACE {
#if HWY_TARGET == HWY_SCALAR || HWY_TARGET == HWY_EMU128
void diffuse_stripe_target(const FsStripe& stripe) {
  diffuse_stripe_scalar(stripe);
}
#else
namespace hn = hwy::HWY_NAMESPACE;

// Adapted from the author's neo-f3kdb process_stripe_fs_wavefront_hwy schedule.
// Each lane owns one row: x = t - 2 * lane. Left/right-up, up and left-up
// dependencies then belong to t-1, t-2 and t-3 respectively. Unlike the old
// unsigned-remainder algorithm, N07 uses signed F32 residuals and ties-to-even.
void diffuse_stripe_target(const FsStripe& stripe) {
  if (!stripe.width || !stripe.height)
    return;
  const hn::ScalableTag<double> d;
  const hn::Rebind<float, decltype(d)> df;
  const std::size_t lanes = hn::Lanes(d);
  const auto zero = hn::Zero(d);
  const double scale = std::ldexp(1.0, stripe.output_depth);
  const double tvscale = std::ldexp(1.0, stripe.output_depth - 8);
  const auto lower = hn::Set(d, stripe.limited ? 16 * tvscale : 0);
  const auto upper = hn::Set(d, stripe.limited ? (stripe.chroma ? 240 : 235) * tvscale : scale - 1);
  const auto factor = hn::Set(d, stripe.integer_input ? std::ldexp(1.0, stripe.output_depth - 16) : scale);
  const auto offset = hn::Set(d, !stripe.integer_input && stripe.chroma ? 0.5 : 0.0);
  const auto w_left = hn::Set(d, 7.0 / 16);
  const auto w_left_up = hn::Set(d, 1.0 / 16);
  const auto w_up = hn::Set(d, 5.0 / 16);
  const auto w_right_up = hn::Set(d, 3.0 / 16);
  // Explicit F32 roundings prevent FMA contraction from changing the specified
  // incoming-error accumulation. All arguments here are finite and bounded.
  const auto f32 = [&](auto v) {
    return hn::PromoteTo(d, hn::DemoteTo(df, v));
  };
  const auto clip = [&](auto v) {
    return hn::Min(upper, hn::Max(lower, v));
  };

  HWY_ALIGN double source[hn::MaxLanes(d)];
  HWY_ALIGN double valid[hn::MaxLanes(d)];
  HWY_ALIGN double output[hn::MaxLanes(d)];
  HWY_ALIGN double residual[hn::MaxLanes(d)];
  std::vector<float> next_boundary(stripe.width);

  const auto values = detail::plane_view(stripe.values, stripe.width, stripe.height);
  const auto target = detail::plane_view(stripe.output, detail::mul(stripe.width, stripe.output_depth == 8 ? 1 : 2),
                                         stripe.height, stripe.output_stride);
  auto value_row = values.cursor();
  auto target_row = target.cursor();
  const double* inputs[hn::MaxLanes(d)];
  std::uint8_t* outputs[hn::MaxLanes(d)];
  for (std::size_t y0 = 0; y0 < stripe.height; y0 += lanes) {
    const auto rows = std::min(lanes, stripe.height - y0);
    // Each lane stays on this row throughout the wavefront's horizontal scan.
    for (std::size_t k = 0; k < rows; ++k) {
      if (y0 || k) {
        ++value_row;
        ++target_row;
      }
      inputs[k] = value_row.ptr();
      outputs[k] = target_row.ptr();
    }
    auto prev1 = zero, prev2 = zero, prev3 = zero;
    const auto end = stripe.width + 2 * (rows - 1);
    const hn::Rebind<std::int64_t, decltype(d)> di;
    const auto offsets = hn::Mul(hn::Iota(di, 0), hn::Set(di, std::int64_t(stripe.width) - 2));
    const auto step = [&](std::size_t t, auto interior_tag) {
      constexpr bool interior = decltype(interior_tag)::value;
      auto left_up = hn::Slide1Up(d, prev3);
      auto up = hn::Slide1Up(d, prev2);
      auto right_up = hn::Slide1Up(d, prev1);
      // Keep bridge contributions separate: pre-summing them would change the
      // mandated order left, left-up, up, right-up at every stripe boundary.
      left_up = hn::InsertLane(left_up, 0, interior || (t > 0 && t <= stripe.width) ? stripe.boundary[t - 1] : 0.0);
      up = hn::InsertLane(up, 0, interior || t < stripe.width ? stripe.boundary[t] : 0.0);
      right_up = hn::InsertLane(right_up, 0, interior || t + 1 < stripe.width ? stripe.boundary[t + 1] : 0.0);
      auto error = f32(hn::Mul(prev1, w_left));
      error = f32(hn::Add(error, f32(hn::Mul(left_up, w_left_up))));
      error = f32(hn::Add(error, f32(hn::Mul(up, w_up))));
      error = f32(hn::Add(error, f32(hn::Mul(right_up, w_right_up))));

      auto samples = zero;
      auto active = hn::Eq(zero, zero);
      if constexpr (interior) {
        samples = hn::GatherIndex(d, inputs[0] + t, offsets);
      } else {
        for (std::size_t k = 0; k < lanes; ++k) {
          const bool enabled = k < rows && t >= 2 * k && t - 2 * k < stripe.width;
          source[k] = enabled ? inputs[k][t - 2 * k] : 0;
          valid[k] = enabled ? 1 : 0;
        }
        active = hn::Gt(hn::LoadU(d, valid), zero);
        samples = hn::LoadU(d, source);
      }
      const auto value = hn::Mul(hn::Add(samples, offset), factor);
      const auto v = clip(f32(hn::Add(clip(value), error)));
      const auto quantized = hn::Round(v);
      const auto next = hn::IfThenElseZero(active, f32(hn::Sub(v, quantized)));
      hn::StoreU(quantized, d, output);
      hn::StoreU(next, d, residual);
      for (std::size_t k = 0; k < rows; ++k) {
        if constexpr (!interior) {
          if (t < 2 * k || t - 2 * k >= stripe.width)
            continue;
        }
        const auto x = t - 2 * k;
        auto* dst = outputs[k];
        if (stripe.output_depth == 8) {
          dst[x] = static_cast<std::uint8_t>(output[k]);
        } else {
          const auto code = static_cast<std::uint16_t>(output[k]);
          std::memcpy(dst + 2 * x, &code, sizeof(code));
        }
        if (k + 1 == rows)
          next_boundary[x] = static_cast<float>(residual[k]);
      }
      prev3 = prev2;
      prev2 = prev1;
      prev1 = next;
    };
    // AVX2 benefits from gathering full interior stripes; wider targets retain
    // the scalar collection schedule to avoid a gather throughput regression.
    const auto first = std::max<std::size_t>(1, 2 * (rows - 1));
    const auto last = stripe.width - 1;
    if (HWY_TARGET == HWY_AVX2 && rows == lanes && first < last) {
      for (std::size_t t = 0; t < first; ++t)
        step(t, std::false_type{});
      for (std::size_t t = first; t < last; ++t)
        step(t, std::true_type{});
      for (std::size_t t = last; t < end; ++t)
        step(t, std::false_type{});
    } else {
      for (std::size_t t = 0; t < end; ++t)
        step(t, std::false_type{});
    }
    std::copy(next_boundary.begin(), next_boundary.end(), stripe.boundary);
  }
}
#endif
} // namespace HWY_NAMESPACE
} // namespace neo_smo::f3kdb
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo::f3kdb {
HWY_EXPORT(diffuse_stripe_target);
void diffuse_stripe(const FsStripe& stripe) {
  HWY_DYNAMIC_DISPATCH(diffuse_stripe_target)(stripe);
}
} // namespace neo_smo::f3kdb
#endif
