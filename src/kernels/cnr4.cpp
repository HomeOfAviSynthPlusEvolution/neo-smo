#include "kernels/dispatch.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <array>
#include "base/fp16.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/cnr4.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/weighted_simd.hpp"
namespace {

inline std::vector<float> calculate_temporal_weights(int wmode, int radius) {
  std::vector<float> weights(static_cast<std::size_t>(radius * 2), 1.0f);
  if (wmode == 0) { // equal
    return weights;
  }
  if (wmode == 1) { // sqrt
    for (int i = 0; i < radius; ++i) {
      const float val = std::sqrt(static_cast<float>(radius - i) / static_cast<float>(radius * 2));
      weights[static_cast<std::size_t>(radius - 1 - i)] = val;
      weights[static_cast<std::size_t>(radius + i)] = val;
    }
    return weights;
  }
  if (wmode == 2) { // sin
    for (int i = 0; i < radius; ++i) {
      const float val = std::sin(static_cast<float>(radius + 1 - i) / static_cast<float>(radius * 2));
      weights[static_cast<std::size_t>(radius - 1 - i)] = val;
      weights[static_cast<std::size_t>(radius + i)] = val;
    }
    return weights;
  }
  if (wmode == 3) { // linear_decrease
    for (int i = 0; i < radius * 2; ++i) {
      weights[static_cast<std::size_t>(i)] = 1.0f / static_cast<float>(i < radius ? (radius - i + 1) : (i - radius + 2));
    }
    return weights;
  }
  return weights;
}

template <typename T>
HWY_NOINLINE void cnr4_process_frame_impl(
    int radius, int depth, int width, int height, std::size_t stride_y, std::size_t stride_uv,
    const T* const curr[3], const T* const curr_ref[3],
    const std::array<const T*, 3>* src, const std::array<const T*, 3>* ref,
    const float* temporal_weights, const std::uint8_t* table_y,
    const std::uint8_t* table_u, const std::uint8_t* table_v, T* dst_u, T* dst_v) {
  // Every product and sum is an integer below 2^53, including 16-bit/radius-10
  // accumulation. Binary64 SIMD therefore keeps exact integer arithmetic while
  // supporting division on all of the generated x86 targets.
  const hn::ScalableTag<double> d;
  const hn::Rebind<float, decltype(d)> df;
  const hn::Rebind<std::int32_t, decltype(d)> di;
  const std::size_t lanes = hn::Lanes(d);
  const double max_value = static_cast<double>(std::uint64_t{1} << (depth * 2));
  const auto max = hn::Set(d, max_value), half = hn::Set(d, max_value / 2);
  const auto divisor = hn::Set(d, max_value * radius * 2);
  const auto weight_scale = hn::Set(d, static_cast<double>(1u << ((depth - 8) * 2)));
  const auto index_scale = hn::Set(d, 1.0 / (1u << (depth - 8)));
  // Widen the tiny LUTs once per pass so Highway can gather directly.
  std::array<float, 256> tables[3];
  const std::uint8_t* input_tables[3] = {table_y, table_u, table_v};
  for (int p = 0; p < 3; ++p)
    for (int i = 0; i < 256; ++i) tables[p][i] = input_tables[p][i];
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; x += static_cast<int>(lanes)) {
      const auto count = std::min(lanes, static_cast<std::size_t>(width - x));
      auto load = [&](const T* const* planes, int p) HWY_ATTR {
        return weighted_load<false>(d, planes[p] + y * (p == 0 ? stride_y : stride_uv) + x, count);
      };
      const auto current_u = load(curr, 1), current_v = load(curr, 2);
      const auto ry = load(curr_ref, 0), ru = load(curr_ref, 1), rv = load(curr_ref, 2);
      auto total_u = hn::Zero(d), total_v = hn::Zero(d);
      for (int i = 0; i < radius * 2; ++i) {
        const auto dy = hn::Abs(hn::Sub(ry, load(ref[i].data(), 0)));
        const auto du = hn::Abs(hn::Sub(ru, load(ref[i].data(), 1)));
        const auto dv = hn::Abs(hn::Sub(rv, load(ref[i].data(), 2)));
        auto lookup = [&](auto diff, int p) HWY_ATTR {
          const auto idx = hn::ConvertTo(di, hn::DemoteTo(df, hn::Min(hn::Mul(diff, index_scale), hn::Set(d, 255))));
          return hn::PromoteTo(d, hn::GatherIndex(df, tables[p].data(), idx));
        };
        const auto wy = lookup(dy, 0);
        auto blend = [&](auto diff, auto center, int p) HWY_ATTR {
          auto weight = hn::Mul(hn::Mul(wy, lookup(diff, p)), weight_scale);
          const auto weighted = hn::Mul(hn::DemoteTo(df, weight), hn::Set(df, temporal_weights[i]));
          weight = hn::PromoteTo(d, weighted_round(df, weighted));
          const auto numerator = hn::Add(hn::Add(hn::Mul(weight, load(src[i].data(), p)),
                                                hn::Mul(hn::Sub(max, weight), center)), half);
          return hn::Mul(hn::Sub(max, hn::Add(dy, diff)), hn::Floor(hn::Div(numerator, max)));
        };
        total_u = hn::Add(total_u, blend(du, current_u, 1));
        total_v = hn::Add(total_v, blend(dv, current_v, 2));
      }
      const auto rounding = hn::Set(d, max_value * radius);
      weighted_store<false>(d, hn::Floor(hn::Div(hn::Add(total_u, rounding), divisor)), dst_u + y * stride_uv + x, count);
      weighted_store<false>(d, hn::Floor(hn::Div(hn::Add(total_v, rounding), divisor)), dst_v + y * stride_uv + x, count);
    }
  }
}

template <typename T>
void cnr4_driver(
    int width,
    int height,
    std::size_t stride_bytes,
    int depth,
    int radius,
    int tmode,
    int wmode,
    const std::uint8_t* const curr_bytes[3],
    const std::uint8_t* const curr_ref_bytes[3],
    const std::array<const std::uint8_t*, 3>* src_bytes,
    const std::array<const std::uint8_t*, 3>* ref_bytes,
    std::size_t num_frames,
    std::uint8_t* dst_u_bytes,
    std::uint8_t* dst_v_bytes,
    const std::uint8_t* table_y,
    const std::uint8_t* table_u,
    const std::uint8_t* table_v
) {
  const std::size_t stride = stride_bytes / sizeof(T);
  const T* const curr[3] = {
      reinterpret_cast<const T*>(curr_bytes[0]),
      reinterpret_cast<const T*>(curr_bytes[1]),
      reinterpret_cast<const T*>(curr_bytes[2]),
  };
  const T* const curr_ref[3] = {
      reinterpret_cast<const T*>(curr_ref_bytes[0]),
      reinterpret_cast<const T*>(curr_ref_bytes[1]),
      reinterpret_cast<const T*>(curr_ref_bytes[2]),
  };

  auto* dst_u = reinterpret_cast<T*>(dst_u_bytes);
  auto* dst_v = reinterpret_cast<T*>(dst_v_bytes);

  if (tmode == 0) { // inv_diff
    const auto t_weights = calculate_temporal_weights(wmode, radius);
    std::vector<std::array<const T*, 3>> src_ptrs(static_cast<std::size_t>(radius * 2));
    std::vector<std::array<const T*, 3>> ref_ptrs(static_cast<std::size_t>(radius * 2));
    for (int i = 0; i < radius * 2; ++i) {
      src_ptrs[static_cast<std::size_t>(i)][0] = reinterpret_cast<const T*>(src_bytes[i][0]);
      src_ptrs[static_cast<std::size_t>(i)][1] = reinterpret_cast<const T*>(src_bytes[i][1]);
      src_ptrs[static_cast<std::size_t>(i)][2] = reinterpret_cast<const T*>(src_bytes[i][2]);
      ref_ptrs[static_cast<std::size_t>(i)][0] = reinterpret_cast<const T*>(ref_bytes[i][0]);
      ref_ptrs[static_cast<std::size_t>(i)][1] = reinterpret_cast<const T*>(ref_bytes[i][1]);
      ref_ptrs[static_cast<std::size_t>(i)][2] = reinterpret_cast<const T*>(ref_bytes[i][2]);
    }

    cnr4_process_frame_impl(
        radius, depth, width, height, stride, stride,
        curr, curr_ref,
        src_ptrs.data(),
        ref_ptrs.data(),
        t_weights.data(), table_y, table_u, table_v, dst_u, dst_v);
  } else {
    // CNR2 multi-pass modes (1..4)
    const bool should_update_src = (tmode == 2 || tmode == 4); // cnr2, cnr2_expanding
    const bool use_expanding_radius = (tmode == 3 || tmode == 4); // expanding modes

    std::vector<std::array<const T*, 3>> srcs(num_frames);
    std::vector<std::array<const T*, 3>> refs(num_frames);
    for (std::size_t i = 0; i < num_frames; ++i) {
      srcs[i][0] = reinterpret_cast<const T*>(src_bytes[i][0]);
      srcs[i][1] = reinterpret_cast<const T*>(src_bytes[i][1]);
      srcs[i][2] = reinterpret_cast<const T*>(src_bytes[i][2]);
      refs[i][0] = reinterpret_cast<const T*>(ref_bytes[i][0]);
      refs[i][1] = reinterpret_cast<const T*>(ref_bytes[i][1]);
      refs[i][2] = reinterpret_cast<const T*>(ref_bytes[i][2]);
    }

    const std::size_t plane_size = stride * static_cast<std::size_t>(height);
    std::vector<T> left_u(plane_size), left_v(plane_size);
    std::vector<T> right_u(plane_size), right_v(plane_size);

    for (int i = 1; i < radius; ++i) {
      const int l_idx = i;
      const int r_idx = static_cast<int>(num_frames) - 1 - i;
      const int cur_radius = use_expanding_radius ? std::min(radius, l_idx) : 1;
      const auto t_weights = calculate_temporal_weights(wmode, cur_radius);

      // Left
      std::vector<std::array<const T*, 3>> tmp_srcs(static_cast<std::size_t>(cur_radius * 2));
      std::vector<std::array<const T*, 3>> tmp_refs(static_cast<std::size_t>(cur_radius * 2));
      if (!use_expanding_radius) {
        tmp_srcs[0][0] = srcs[static_cast<std::size_t>(l_idx - 1)][0];
        tmp_srcs[0][1] = srcs[static_cast<std::size_t>(l_idx - 1)][1];
        tmp_srcs[0][2] = srcs[static_cast<std::size_t>(l_idx - 1)][2];
        tmp_srcs[1][0] = srcs[static_cast<std::size_t>(l_idx + 1)][0];
        tmp_srcs[1][1] = srcs[static_cast<std::size_t>(l_idx + 1)][1];
        tmp_srcs[1][2] = srcs[static_cast<std::size_t>(l_idx + 1)][2];

        tmp_refs[0][0] = refs[static_cast<std::size_t>(l_idx - 1)][0];
        tmp_refs[0][1] = refs[static_cast<std::size_t>(l_idx - 1)][1];
        tmp_refs[0][2] = refs[static_cast<std::size_t>(l_idx - 1)][2];
        tmp_refs[1][0] = refs[static_cast<std::size_t>(l_idx + 1)][0];
        tmp_refs[1][1] = refs[static_cast<std::size_t>(l_idx + 1)][1];
        tmp_refs[1][2] = refs[static_cast<std::size_t>(l_idx + 1)][2];
      } else {
        for (int j = 0; j < cur_radius; ++j) {
          tmp_srcs[static_cast<std::size_t>(j)][0] = srcs[static_cast<std::size_t>(l_idx - cur_radius + j)][0];
          tmp_srcs[static_cast<std::size_t>(j)][1] = srcs[static_cast<std::size_t>(l_idx - cur_radius + j)][1];
          tmp_srcs[static_cast<std::size_t>(j)][2] = srcs[static_cast<std::size_t>(l_idx - cur_radius + j)][2];

          tmp_refs[static_cast<std::size_t>(j)][0] = refs[static_cast<std::size_t>(l_idx - cur_radius + j)][0];
          tmp_refs[static_cast<std::size_t>(j)][1] = refs[static_cast<std::size_t>(l_idx - cur_radius + j)][1];
          tmp_refs[static_cast<std::size_t>(j)][2] = refs[static_cast<std::size_t>(l_idx - cur_radius + j)][2];

          tmp_srcs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][0] = srcs[static_cast<std::size_t>(l_idx + cur_radius - j)][0];
          tmp_srcs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][1] = srcs[static_cast<std::size_t>(l_idx + cur_radius - j)][1];
          tmp_srcs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][2] = srcs[static_cast<std::size_t>(l_idx + cur_radius - j)][2];

          tmp_refs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][0] = refs[static_cast<std::size_t>(l_idx + cur_radius - j)][0];
          tmp_refs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][1] = refs[static_cast<std::size_t>(l_idx + cur_radius - j)][1];
          tmp_refs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][2] = refs[static_cast<std::size_t>(l_idx + cur_radius - j)][2];
        }
      }

      cnr4_process_frame_impl(
          cur_radius, depth, width, height, stride, stride,
          srcs[static_cast<std::size_t>(l_idx)].data(), refs[static_cast<std::size_t>(l_idx)].data(),
          tmp_srcs.data(),
          tmp_refs.data(),
          t_weights.data(), table_y, table_u, table_v, left_u.data(), left_v.data());

      srcs[static_cast<std::size_t>(l_idx)][1] = should_update_src ? left_u.data() : srcs[static_cast<std::size_t>(l_idx)][1];
      srcs[static_cast<std::size_t>(l_idx)][2] = should_update_src ? left_v.data() : srcs[static_cast<std::size_t>(l_idx)][2];
      refs[static_cast<std::size_t>(l_idx)][1] = left_u.data();
      refs[static_cast<std::size_t>(l_idx)][2] = left_v.data();

      // Right
      if (!use_expanding_radius) {
        tmp_srcs[0][0] = srcs[static_cast<std::size_t>(r_idx - 1)][0];
        tmp_srcs[0][1] = srcs[static_cast<std::size_t>(r_idx - 1)][1];
        tmp_srcs[0][2] = srcs[static_cast<std::size_t>(r_idx - 1)][2];
        tmp_srcs[1][0] = srcs[static_cast<std::size_t>(r_idx + 1)][0];
        tmp_srcs[1][1] = srcs[static_cast<std::size_t>(r_idx + 1)][1];
        tmp_srcs[1][2] = srcs[static_cast<std::size_t>(r_idx + 1)][2];

        tmp_refs[0][0] = refs[static_cast<std::size_t>(r_idx - 1)][0];
        tmp_refs[0][1] = refs[static_cast<std::size_t>(r_idx - 1)][1];
        tmp_refs[0][2] = refs[static_cast<std::size_t>(r_idx - 1)][2];
        tmp_refs[1][0] = refs[static_cast<std::size_t>(r_idx + 1)][0];
        tmp_refs[1][1] = refs[static_cast<std::size_t>(r_idx + 1)][1];
        tmp_refs[1][2] = refs[static_cast<std::size_t>(r_idx + 1)][2];
      } else {
        for (int j = 0; j < cur_radius; ++j) {
          tmp_srcs[static_cast<std::size_t>(j)][0] = srcs[static_cast<std::size_t>(r_idx - cur_radius + j)][0];
          tmp_srcs[static_cast<std::size_t>(j)][1] = srcs[static_cast<std::size_t>(r_idx - cur_radius + j)][1];
          tmp_srcs[static_cast<std::size_t>(j)][2] = srcs[static_cast<std::size_t>(r_idx - cur_radius + j)][2];

          tmp_refs[static_cast<std::size_t>(j)][0] = refs[static_cast<std::size_t>(r_idx - cur_radius + j)][0];
          tmp_refs[static_cast<std::size_t>(j)][1] = refs[static_cast<std::size_t>(r_idx - cur_radius + j)][1];
          tmp_refs[static_cast<std::size_t>(j)][2] = refs[static_cast<std::size_t>(r_idx - cur_radius + j)][2];

          tmp_srcs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][0] = srcs[static_cast<std::size_t>(r_idx + cur_radius - j)][0];
          tmp_srcs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][1] = srcs[static_cast<std::size_t>(r_idx + cur_radius - j)][1];
          tmp_srcs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][2] = srcs[static_cast<std::size_t>(r_idx + cur_radius - j)][2];

          tmp_refs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][0] = refs[static_cast<std::size_t>(r_idx + cur_radius - j)][0];
          tmp_refs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][1] = refs[static_cast<std::size_t>(r_idx + cur_radius - j)][1];
          tmp_refs[static_cast<std::size_t>(cur_radius * 2 - 1 - j)][2] = refs[static_cast<std::size_t>(r_idx + cur_radius - j)][2];
        }
      }

      cnr4_process_frame_impl(
          cur_radius, depth, width, height, stride, stride,
          srcs[static_cast<std::size_t>(r_idx)].data(), refs[static_cast<std::size_t>(r_idx)].data(),
          tmp_srcs.data(),
          tmp_refs.data(),
          t_weights.data(), table_y, table_u, table_v, right_u.data(), right_v.data());

      srcs[static_cast<std::size_t>(r_idx)][1] = should_update_src ? right_u.data() : srcs[static_cast<std::size_t>(r_idx)][1];
      srcs[static_cast<std::size_t>(r_idx)][2] = should_update_src ? right_v.data() : srcs[static_cast<std::size_t>(r_idx)][2];
      refs[static_cast<std::size_t>(r_idx)][1] = right_u.data();
      refs[static_cast<std::size_t>(r_idx)][2] = right_v.data();
    }

    // Final combination for center frame
    const auto t_weights = calculate_temporal_weights(wmode, radius);
    std::vector<std::array<const T*, 3>> tmp_srcs(static_cast<std::size_t>(radius * 2));
    std::vector<std::array<const T*, 3>> tmp_refs(static_cast<std::size_t>(radius * 2));
    for (int j = 0; j < radius; ++j) {
      tmp_srcs[static_cast<std::size_t>(j)][0] = srcs[static_cast<std::size_t>(j)][0];
      tmp_srcs[static_cast<std::size_t>(j)][1] = srcs[static_cast<std::size_t>(j)][1];
      tmp_srcs[static_cast<std::size_t>(j)][2] = srcs[static_cast<std::size_t>(j)][2];
      tmp_refs[static_cast<std::size_t>(j)][0] = refs[static_cast<std::size_t>(j)][0];
      tmp_refs[static_cast<std::size_t>(j)][1] = refs[static_cast<std::size_t>(j)][1];
      tmp_refs[static_cast<std::size_t>(j)][2] = refs[static_cast<std::size_t>(j)][2];

      tmp_srcs[static_cast<std::size_t>(radius * 2 - 1 - j)][0] = srcs[num_frames - 1 - static_cast<std::size_t>(j)][0];
      tmp_srcs[static_cast<std::size_t>(radius * 2 - 1 - j)][1] = srcs[num_frames - 1 - static_cast<std::size_t>(j)][1];
      tmp_srcs[static_cast<std::size_t>(radius * 2 - 1 - j)][2] = srcs[num_frames - 1 - static_cast<std::size_t>(j)][2];
      tmp_refs[static_cast<std::size_t>(radius * 2 - 1 - j)][0] = refs[num_frames - 1 - static_cast<std::size_t>(j)][0];
      tmp_refs[static_cast<std::size_t>(radius * 2 - 1 - j)][1] = refs[num_frames - 1 - static_cast<std::size_t>(j)][1];
      tmp_refs[static_cast<std::size_t>(radius * 2 - 1 - j)][2] = refs[num_frames - 1 - static_cast<std::size_t>(j)][2];
    }

    cnr4_process_frame_impl(
        radius, depth, width, height, stride, stride,
        curr, curr_ref,
        tmp_srcs.data(),
        tmp_refs.data(),
        t_weights.data(), table_y, table_u, table_v, dst_u, dst_v);
  }
}

} // namespace

void dispatch_cnr4_target(
    DataType dtype,
    int width,
    int height,
    std::size_t stride_bytes,
    int depth,
    int radius,
    int tmode,
    int wmode,
    const std::uint8_t* const curr[3],
    const std::uint8_t* const curr_ref[3],
    const std::array<const std::uint8_t*, 3>* src,
    const std::array<const std::uint8_t*, 3>* ref,
    std::size_t num_frames,
    std::uint8_t* dst_u,
    std::uint8_t* dst_v,
    const std::uint8_t* table_y,
    const std::uint8_t* table_u,
    const std::uint8_t* table_v
) {
  if (dtype == DataType::U8) {
    cnr4_driver<std::uint8_t>(
        width, height, stride_bytes, depth, radius, tmode, wmode,
        curr, curr_ref, src, ref, num_frames, dst_u, dst_v,
        table_y, table_u, table_v);
  } else if (dtype == DataType::U16) {
    cnr4_driver<std::uint16_t>(
        width, height, stride_bytes, depth, radius, tmode, wmode,
        curr, curr_ref, src, ref, num_frames, dst_u, dst_v,
        table_y, table_u, table_v);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {

HWY_EXPORT(dispatch_cnr4_target);

void process_cnr4_frame(
    DataType dtype,
    int width,
    int height,
    std::size_t stride_bytes,
    int depth,
    int radius,
    int tmode,
    int wmode,
    const std::uint8_t* const curr[3],
    const std::uint8_t* const curr_ref[3],
    const std::array<const std::uint8_t*, 3>* src,
    const std::array<const std::uint8_t*, 3>* ref,
    std::size_t num_frames,
    std::uint8_t* dst_u,
    std::uint8_t* dst_v,
    const std::uint8_t* table_y,
    const std::uint8_t* table_u,
    const std::uint8_t* table_v
) {
  HWY_DYNAMIC_DISPATCH(dispatch_cnr4_target)(
      dtype, width, height, stride_bytes, depth, radius, tmode, wmode,
      curr, curr_ref, src, ref, num_frames, dst_u, dst_v,
      table_y, table_u, table_v);
}

} // namespace neo_smo
#endif
