#include "kernels/dispatch.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <array>
#include "base/fp16.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/ttempsmooth.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
#include "common/weighted_simd.hpp"

namespace {

template <typename T>
HWY_NOINLINE void ttempsmooth_row_impl(
    int width,
    std::size_t src_stride,
    std::size_t ref_stride,
    std::size_t dst_stride,
    int y,
    const T* curr_plane,
    const T* curr_ref_plane,
    const T* const* prev_planes,
    const T* const* prev_ref_planes,
    const T* const* next_planes,
    const T* const* next_ref_planes,
    T* dst_plane,
    int maxr,
    int num_prev,
    int num_next,
    float threshold,
    bool fp,
    int bits_per_sample,
    int weight_mode,
    float center_weight,
    const float* temporal_weights,
    const float* temporal_difference_weights
) {
  const hn::ScalableTag<float> df;
  const std::size_t lanes = hn::Lanes(df);
  constexpr bool is_float = std::is_floating_point_v<T>;
  const int shift = is_float ? 0 : (bits_per_sample - 8);
  const float format_max = is_float ? 1.0f : static_cast<float>((1u << bits_per_sample) - 1u);

  const T* curr_row = curr_plane + static_cast<std::size_t>(y) * src_stride;
  const T* curr_ref_row = curr_ref_plane + static_cast<std::size_t>(y) * ref_stride;
  T* dst_row = dst_plane + static_cast<std::size_t>(y) * dst_stride;

  const auto thresh_vec = hn::Set(df, threshold);
  const auto one_vec = hn::Set(df, 1.0f);
  const auto center_weight_vec = hn::Set(df, center_weight);

  const hn::RebindToSigned<decltype(df)> di;
  for (int x = 0; x < width; x += static_cast<int>(lanes)) {
    const std::size_t count = std::min(lanes, static_cast<std::size_t>(width - x));
    const auto current_pixel = weighted_load<false>(df, curr_ref_row + x, count);
    const auto currv = weighted_load<false>(df, curr_row + x, count);

    auto weight_sum = center_weight_vec;
    auto sum = hn::Mul(currv, center_weight_vec);

    // Two directions: 0 = past (prev), 1 = future (next)
    for (int dir = 0; dir < 2; ++dir) {
      const int num_frames = (dir == 0) ? num_prev : num_next;
      if (num_frames <= 0) {
        break;
      }

      const T* const* dir_src = (dir == 0) ? prev_planes : next_planes;
      const T* const* dir_ref = (dir == 0) ? prev_ref_planes : next_ref_planes;

      auto active_mask = hn::Eq(one_vec, one_vec);
      auto prev_temporal_pixel = current_pixel; // for step 0, prev is itself

      for (int i = 0; i < num_frames; ++i) {
        const T* ref_row = dir_ref[i] + static_cast<std::size_t>(y) * ref_stride;
        const T* src_row = dir_src[i] + static_cast<std::size_t>(y) * src_stride;

        const auto temporal_pixel1 = weighted_load<false>(df, ref_row + x, count);
        const auto srcv = weighted_load<false>(df, src_row + x, count);

        auto diff = hn::Abs(hn::Sub(current_pixel, temporal_pixel1));
        if constexpr (is_float) {
          diff = hn::Min(diff, one_vec);
        }

        auto lt_diff = hn::Lt(diff, thresh_vec);
        if (i == 0) {
          active_mask = hn::And(active_mask, lt_diff);
        } else {
          auto tdiff = hn::Abs(hn::Sub(temporal_pixel1, prev_temporal_pixel));
          if constexpr (is_float) {
            tdiff = hn::Min(tdiff, one_vec);
          }
          auto lt_tdiff = hn::Lt(tdiff, thresh_vec);
          active_mask = hn::And(active_mask, hn::And(lt_diff, lt_tdiff));
        }

        if (hn::AllFalse(df, active_mask)) {
          break;
        }

        auto weight = one_vec;
        if (weight_mode == 1) { // Temporal
          weight = hn::Set(df, temporal_weights[1 + i]);
        } else { // InverseDifference
          const float* table = temporal_difference_weights + static_cast<std::size_t>(i) * 256;
          auto index_value = diff;
          if constexpr (is_float) index_value = hn::Mul(diff, hn::Set(df, 255.0f));
          else if constexpr (sizeof(T) == 2) index_value = hn::Mul(diff, hn::Set(df, 1.0f / (1u << shift)));
          // Clamp before conversion/gather, including inactive padded lanes.
          index_value = hn::Min(hn::Max(index_value, hn::Zero(df)), hn::Set(df, 255.0f));
          weight = hn::GatherIndex(df, table, hn::ConvertTo(di, index_value));
        }

        weight_sum = hn::IfThenElse(active_mask, hn::Add(weight_sum, weight), weight_sum);
        sum = hn::IfThenElse(active_mask, hn::Add(sum, hn::Mul(srcv, weight)), sum);

        prev_temporal_pixel = temporal_pixel1;
      }
    }

    auto res = one_vec;
    if (fp) {
      res = hn::Add(hn::Mul(currv, hn::Sub(one_vec, weight_sum)), sum);
    } else {
      res = hn::Div(sum, weight_sum);
    }

    if constexpr (!is_float)
      res = hn::Min(hn::Max(weighted_round(df, res), hn::Zero(df)), hn::Set(df, format_max));
    weighted_store<false>(df, res, dst_row + x, count);
  }
}

template <typename T>
void ttempsmooth_plane_impl(
    int width,
    int height,
    std::size_t src_stride,
    std::size_t ref_stride,
    std::size_t dst_stride,
    const T* curr_plane,
    const T* curr_ref_plane,
    const T* const* prev_planes,
    const T* const* prev_ref_planes,
    const T* const* next_planes,
    const T* const* next_ref_planes,
    T* dst_plane,
    int maxr,
    int num_prev,
    int num_next,
    float threshold,
    bool fp,
    int bits_per_sample,
    int weight_mode,
    float center_weight,
    const float* temporal_weights,
    const float* temporal_difference_weights
) {
  for (int y = 0; y < height; ++y) {
    ttempsmooth_row_impl(
        width, src_stride, ref_stride, dst_stride, y,
        curr_plane, curr_ref_plane,
        prev_planes, prev_ref_planes,
        next_planes, next_ref_planes,
        dst_plane,
        maxr, num_prev, num_next,
        threshold, fp, bits_per_sample, weight_mode,
        center_weight, temporal_weights, temporal_difference_weights);
  }
}

} // namespace

void dispatch_ttempsmooth_target(
    DataType dtype,
    int width,
    int height,
    std::size_t src_stride_bytes,
    std::size_t ref_stride_bytes,
    std::size_t dst_stride_bytes,
    const std::uint8_t* curr,
    const std::uint8_t* curr_ref,
    const std::uint8_t* const* prev,
    const std::uint8_t* const* prev_ref,
    const std::uint8_t* const* next,
    const std::uint8_t* const* next_ref,
    std::uint8_t* dstp,
    int maxr,
    int num_prev,
    int num_next,
    float threshold,
    bool fp,
    int bits_per_sample,
    int weight_mode,
    float center_weight,
    const float* temporal_weights,
    const float* temporal_difference_weights
) {
  std::array<const std::uint16_t*, 7> prev16{}, prev_ref16{}, next16{}, next_ref16{};
  std::array<const float*, 7> prev32{}, prev_ref32{}, next32{}, next_ref32{};
  for (int i = 0; i < num_prev; ++i) {
    prev16[i] = reinterpret_cast<const std::uint16_t*>(prev[i]);
    prev_ref16[i] = reinterpret_cast<const std::uint16_t*>(prev_ref[i]);
    prev32[i] = reinterpret_cast<const float*>(prev[i]);
    prev_ref32[i] = reinterpret_cast<const float*>(prev_ref[i]);
  }
  for (int i = 0; i < num_next; ++i) {
    next16[i] = reinterpret_cast<const std::uint16_t*>(next[i]);
    next_ref16[i] = reinterpret_cast<const std::uint16_t*>(next_ref[i]);
    next32[i] = reinterpret_cast<const float*>(next[i]);
    next_ref32[i] = reinterpret_cast<const float*>(next_ref[i]);
  }
  if (dtype == DataType::U8) {
    ttempsmooth_plane_impl(
        width, height,
        src_stride_bytes, ref_stride_bytes, dst_stride_bytes,
        curr, curr_ref,
        reinterpret_cast<const std::uint8_t* const*>(prev),
        reinterpret_cast<const std::uint8_t* const*>(prev_ref),
        reinterpret_cast<const std::uint8_t* const*>(next),
        reinterpret_cast<const std::uint8_t* const*>(next_ref),
        dstp,
        maxr, num_prev, num_next,
        threshold, fp, bits_per_sample, weight_mode,
        center_weight, temporal_weights, temporal_difference_weights);
  } else if (dtype == DataType::U16) {
    ttempsmooth_plane_impl(
        width, height,
        src_stride_bytes / sizeof(std::uint16_t),
        ref_stride_bytes / sizeof(std::uint16_t),
        dst_stride_bytes / sizeof(std::uint16_t),
        reinterpret_cast<const std::uint16_t*>(curr),
        reinterpret_cast<const std::uint16_t*>(curr_ref),
        prev16.data(),
        prev_ref16.data(),
        next16.data(),
        next_ref16.data(),
        reinterpret_cast<std::uint16_t*>(dstp),
        maxr, num_prev, num_next,
        threshold, fp, bits_per_sample, weight_mode,
        center_weight, temporal_weights, temporal_difference_weights);
  } else if (dtype == DataType::F32) {
    ttempsmooth_plane_impl(
        width, height,
        src_stride_bytes / sizeof(float),
        ref_stride_bytes / sizeof(float),
        dst_stride_bytes / sizeof(float),
        reinterpret_cast<const float*>(curr),
        reinterpret_cast<const float*>(curr_ref),
        prev32.data(),
        prev_ref32.data(),
        next32.data(),
        next_ref32.data(),
        reinterpret_cast<float*>(dstp),
        maxr, num_prev, num_next,
        threshold, fp, bits_per_sample, weight_mode,
        center_weight, temporal_weights, temporal_difference_weights);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {

HWY_EXPORT(dispatch_ttempsmooth_target);

void process_ttempsmooth_plane(
    DataType dtype,
    int width,
    int height,
    std::size_t src_stride_bytes,
    std::size_t ref_stride_bytes,
    std::size_t dst_stride_bytes,
    const std::uint8_t* curr,
    const std::uint8_t* curr_ref,
    const std::uint8_t* const* prev,
    const std::uint8_t* const* prev_ref,
    const std::uint8_t* const* next,
    const std::uint8_t* const* next_ref,
    std::uint8_t* dstp,
    int maxr,
    int num_prev,
    int num_next,
    float threshold,
    bool fp,
    int bits_per_sample,
    int weight_mode,
    float center_weight,
    const float* temporal_weights,
    const float* temporal_difference_weights
) {
  HWY_DYNAMIC_DISPATCH(dispatch_ttempsmooth_target)(
      dtype, width, height,
      src_stride_bytes, ref_stride_bytes, dst_stride_bytes,
      curr, curr_ref, prev, prev_ref, next, next_ref, dstp,
      maxr, num_prev, num_next,
      threshold, fp, bits_per_sample, weight_mode,
      center_weight, temporal_weights, temporal_difference_weights);
}

} // namespace neo_smo
#endif
