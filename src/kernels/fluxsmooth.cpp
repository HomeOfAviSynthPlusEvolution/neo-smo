#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/fluxsmooth.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/grid.hpp"
#include "common/float_arithmetic.hpp"

// ---------------------------------------------------------------------------
// FluxSmoothT
// ---------------------------------------------------------------------------
template <typename T>
void fluxsmooth_t_int_impl(T temporal_threshold, const T* prevp, const T* currp, const T* nextp,
                           T* dstp, int width, int height, std::size_t prev_stride,
                           std::size_t curr_stride, std::size_t next_stride, std::size_t dst_stride) {
  hn::ScalableTag<std::int32_t> d;
  const std::size_t lanes = hn::Lanes(d);
  const auto thresh = hn::Set(d, static_cast<std::int32_t>(temporal_threshold));
  const auto one = hn::Set(d, 1);
  const auto zero = hn::Zero(d);

  std::vector<std::int32_t> b_prev(static_cast<std::size_t>(width) + lanes, 0);
  std::vector<std::int32_t> b_curr(static_cast<std::size_t>(width) + lanes, 0);
  std::vector<std::int32_t> b_next(static_cast<std::size_t>(width) + lanes, 0);
  std::vector<std::int32_t> out_i32(static_cast<std::size_t>(width) + lanes, 0);

  for (int y = 0; y < height; ++y) {
    const T* p_row = prevp + static_cast<std::size_t>(y) * prev_stride;
    const T* c_row = currp + static_cast<std::size_t>(y) * curr_stride;
    const T* n_row = nextp + static_cast<std::size_t>(y) * next_stride;
    for (int x = 0; x < width; ++x) {
      b_prev[static_cast<std::size_t>(x)] = static_cast<std::int32_t>(p_row[x]);
      b_curr[static_cast<std::size_t>(x)] = static_cast<std::int32_t>(c_row[x]);
      b_next[static_cast<std::size_t>(x)] = static_cast<std::int32_t>(n_row[x]);
    }

    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto p = hn::LoadU(d, b_prev.data() + x);
      const auto c = hn::LoadU(d, b_curr.data() + x);
      const auto n = hn::LoadU(d, b_next.data() + x);

      const auto prevnextless = hn::And(hn::Lt(p, c), hn::Lt(n, c));
      const auto prevnextmore = hn::And(hn::Gt(p, c), hn::Gt(n, c));
      const auto mask_either = hn::Or(prevnextless, prevnextmore);

      const auto p_diff = hn::AbsDiff(p, c);
      const auto n_diff = hn::AbsDiff(n, c);

      const auto p_match = hn::Le(p_diff, thresh);
      const auto n_match = hn::Le(n_diff, thresh);

      auto sum = c;
      auto count = one;

      sum = hn::Add(sum, hn::IfThenElse(p_match, p, zero));
      count = hn::Add(count, hn::IfThenElse(p_match, one, zero));

      sum = hn::Add(sum, hn::IfThenElse(n_match, n, zero));
      count = hn::Add(count, hn::IfThenElse(n_match, one, zero));

      // (sum * 2 + count) / (count * 2)
      const auto numerator = hn::Add(hn::ShiftLeft<1>(sum), count);
      const auto denominator = hn::ShiftLeft<1>(count);
      const auto filtered = hn::Div(numerator, denominator);

      const auto res = hn::IfThenElse(mask_either, filtered, c);
      hn::StoreU(res, d, out_i32.data() + x);
    }

    for (int x = 0; x < width; ++x) {
      dst_row[x] = static_cast<T>(out_i32[static_cast<std::size_t>(x)]);
    }
  }
}

template <bool IsF16, typename StorageT>
void fluxsmooth_t_float_impl(float temporal_threshold, const StorageT* prevp, const StorageT* currp,
                             const StorageT* nextp, StorageT* dstp, int width, int height,
                             std::size_t prev_stride, std::size_t curr_stride,
                             std::size_t next_stride, std::size_t dst_stride) {
  const hn::ScalableTag<float> d;
  const std::size_t lanes = hn::Lanes(d);
  const auto thresh = hn::Set(d, temporal_threshold);
  const auto one = hn::Set(d, 1.0f);
  const auto zero = hn::Zero(d);

  if constexpr (IsF16) {
    std::vector<float> b_prev(static_cast<std::size_t>(width) + lanes);
    std::vector<float> b_curr(static_cast<std::size_t>(width) + lanes);
    std::vector<float> b_next(static_cast<std::size_t>(width) + lanes);
    std::vector<float> b_out(static_cast<std::size_t>(width) + lanes);

    for (int y = 0; y < height; ++y) {
      fill_mirrored_row_fp16_to_fp32(b_prev.data(), prevp + static_cast<std::size_t>(y) * prev_stride, width, 0);
      fill_mirrored_row_fp16_to_fp32(b_curr.data(), currp + static_cast<std::size_t>(y) * curr_stride, width, 0);
      fill_mirrored_row_fp16_to_fp32(b_next.data(), nextp + static_cast<std::size_t>(y) * next_stride, width, 0);

      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto p = hn::LoadU(d, b_prev.data() + x);
        const auto c = hn::LoadU(d, b_curr.data() + x);
        const auto n = hn::LoadU(d, b_next.data() + x);

        const auto prevnextless = hn::And(hn::Lt(p, c), hn::Lt(n, c));
        const auto prevnextmore = hn::And(hn::Gt(p, c), hn::Gt(n, c));
        const auto mask_either = hn::Or(prevnextless, prevnextmore);

        const auto p_diff = float_abs_diff<true>(d, p, c);
        const auto n_diff = float_abs_diff<true>(d, n, c);

        const auto p_match = hn::Le(p_diff, thresh);
        const auto n_match = hn::Le(n_diff, thresh);

        auto sum = c;
        auto count = one;

        sum = float_add<true>(d, sum, hn::IfThenElse(p_match, p, zero));
        count = float_add<true>(d, count, hn::IfThenElse(p_match, one, zero));

        sum = float_add<true>(d, sum, hn::IfThenElse(n_match, n, zero));
        count = float_add<true>(d, count, hn::IfThenElse(n_match, one, zero));

        const auto filtered = float_div<true>(d, sum, count);
        const auto res = hn::IfThenElse(mask_either, filtered, c);

        hn::StoreU(res, d, b_out.data() + x);
      }
      convert_row_fp32_to_fp16(dst_row, b_out.data(), width);
    }
  } else {
    for (int y = 0; y < height; ++y) {
      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const std::size_t rem = static_cast<std::size_t>(width - x);
        const std::size_t count_lanes = std::min(lanes, rem);

        const auto p = hn::LoadN(d, prevp + static_cast<std::size_t>(y) * prev_stride + x, count_lanes);
        const auto c = hn::LoadN(d, currp + static_cast<std::size_t>(y) * curr_stride + x, count_lanes);
        const auto n = hn::LoadN(d, nextp + static_cast<std::size_t>(y) * next_stride + x, count_lanes);

        const auto prevnextless = hn::And(hn::Lt(p, c), hn::Lt(n, c));
        const auto prevnextmore = hn::And(hn::Gt(p, c), hn::Gt(n, c));
        const auto mask_either = hn::Or(prevnextless, prevnextmore);

        const auto p_diff = hn::Abs(hn::Sub(p, c));
        const auto n_diff = hn::Abs(hn::Sub(n, c));

        const auto p_match = hn::Le(p_diff, thresh);
        const auto n_match = hn::Le(n_diff, thresh);

        auto sum = c;
        auto count = one;

        sum = hn::Add(sum, hn::IfThenElse(p_match, p, zero));
        count = hn::Add(count, hn::IfThenElse(p_match, one, zero));

        sum = hn::Add(sum, hn::IfThenElse(n_match, n, zero));
        count = hn::Add(count, hn::IfThenElse(n_match, one, zero));

        const auto filtered = hn::Div(sum, count);
        const auto res = hn::IfThenElse(mask_either, filtered, c);

        if (rem >= lanes) hn::StoreU(res, d, dst_row + x);
        else hn::StoreN(res, d, dst_row + x, rem);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// FluxSmoothST
// ---------------------------------------------------------------------------
template <typename T>
void fluxsmooth_st_int_impl(std::int32_t temporal_threshold, std::int32_t spatial_threshold, const T* prevp, const T* currp,
                            const T* nextp, T* dstp, int width, int height, std::size_t prev_stride,
                            std::size_t curr_stride, std::size_t next_stride, std::size_t dst_stride) {
  constexpr int kRadius = 1;
  hn::ScalableTag<std::int32_t> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + lanes;
  const auto t_thresh = hn::Set(d, static_cast<std::int32_t>(temporal_threshold));
  const auto s_thresh = hn::Set(d, static_cast<std::int32_t>(spatial_threshold));
  const auto one = hn::Set(d, 1);
  const auto zero = hn::Zero(d);

  copy_first_n_lines(dstp, currp, static_cast<std::size_t>(width), dst_stride, curr_stride, 1);

  std::vector<std::int32_t> b_prev(static_cast<std::size_t>(width) + lanes, 0);
  std::vector<std::int32_t> b_next(static_cast<std::size_t>(width) + lanes, 0);
  std::vector<std::int32_t> b_curr(3 * padded_len, 0);
  std::array<std::int32_t*, 3> r_curr{b_curr.data() + 0 * padded_len + kRadius,
                                      b_curr.data() + 1 * padded_len + kRadius,
                                      b_curr.data() + 2 * padded_len + kRadius};

  auto fill_i32 = [&](std::int32_t* dst, const T* srow) {
    for (std::int64_t x = -kRadius; x < static_cast<std::int64_t>(width) + kRadius; ++x) {
      const std::size_t idx = (x >= 0 && x < width) ? static_cast<std::size_t>(x) : mirror_index(x, width);
      dst[x] = static_cast<std::int32_t>(srow[idx]);
    }
  };

  std::vector<std::int32_t> out_i32(static_cast<std::size_t>(width) + lanes);

  for (int y = 1; y < height - 1; ++y) {
    fill_i32(r_curr[0], currp + static_cast<std::size_t>(y - 1) * curr_stride);
    fill_i32(r_curr[1], currp + static_cast<std::size_t>(y) * curr_stride);
    fill_i32(r_curr[2], currp + static_cast<std::size_t>(y + 1) * curr_stride);

    for (int x = 0; x < width; ++x) {
      b_prev[static_cast<std::size_t>(x)] = static_cast<std::int32_t>(prevp[static_cast<std::size_t>(y) * prev_stride + x]);
      b_next[static_cast<std::size_t>(x)] = static_cast<std::int32_t>(nextp[static_cast<std::size_t>(y) * next_stride + x]);
    }

    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto p = hn::LoadU(d, b_prev.data() + x);
      const auto n = hn::LoadU(d, b_next.data() + x);
      const auto gc = Grid3x3<decltype(d)>::load(d, r_curr[0], r_curr[1], r_curr[2], x);
      const auto c = gc.center_center;

      const auto prevnextless = hn::And(hn::Lt(p, c), hn::Lt(n, c));
      const auto prevnextmore = hn::And(hn::Gt(p, c), hn::Gt(n, c));
      const auto mask_either = hn::Or(prevnextless, prevnextmore);

      auto sum = c;
      auto count = one;

      // Temporal neighbors
      const auto p_match = hn::Le(hn::AbsDiff(p, c), t_thresh);
      const auto n_match = hn::Le(hn::AbsDiff(n, c), t_thresh);
      sum = hn::Add(sum, hn::IfThenElse(p_match, p, zero));
      count = hn::Add(count, hn::IfThenElse(p_match, one, zero));
      sum = hn::Add(sum, hn::IfThenElse(n_match, n, zero));
      count = hn::Add(count, hn::IfThenElse(n_match, one, zero));

      // 8 spatial neighbors in current frame
      const hn::Vec<decltype(d)> spatial_neighbors[8] = {
          gc.top_left, gc.top_center, gc.top_right,
          gc.center_left, gc.center_right,
          gc.bottom_left, gc.bottom_center, gc.bottom_right,
      };

      for (int i = 0; i < 8; ++i) {
        const auto sn = spatial_neighbors[i];
        const auto s_match = hn::Le(hn::AbsDiff(sn, c), s_thresh);
        sum = hn::Add(sum, hn::IfThenElse(s_match, sn, zero));
        count = hn::Add(count, hn::IfThenElse(s_match, one, zero));
      }

      const auto numerator = hn::Add(hn::ShiftLeft<1>(sum), count);
      const auto denominator = hn::ShiftLeft<1>(count);
      const auto filtered = hn::Div(numerator, denominator);

      const auto res = hn::IfThenElse(mask_either, filtered, c);
      hn::StoreU(res, d, out_i32.data() + x);
    }

    dst_row[0] = currp[static_cast<std::size_t>(y) * curr_stride];
    for (int x = 1; x < width - 1; ++x) {
      dst_row[x] = static_cast<T>(out_i32[static_cast<std::size_t>(x)]);
    }
    dst_row[width - 1] = currp[static_cast<std::size_t>(y) * curr_stride + width - 1];
  }

  copy_last_n_lines(dstp, currp, static_cast<std::size_t>(width), static_cast<std::size_t>(height), dst_stride, curr_stride, 1);
}

template <bool IsF16, typename StorageT>
void fluxsmooth_st_float_impl(float temporal_threshold, float spatial_threshold,
                             const StorageT* prevp, const StorageT* currp, const StorageT* nextp,
                             StorageT* dstp, int width, int height, std::size_t prev_stride,
                             std::size_t curr_stride, std::size_t next_stride, std::size_t dst_stride) {
  constexpr int kRadius = 1;
  const hn::ScalableTag<float> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + lanes;
  const auto t_thresh = hn::Set(d, temporal_threshold);
  const auto s_thresh = hn::Set(d, spatial_threshold);
  const auto one = hn::Set(d, 1.0f);
  const auto zero = hn::Zero(d);

  copy_first_n_lines(dstp, currp, static_cast<std::size_t>(width), dst_stride, curr_stride, 1);

  std::vector<float> b_curr(3 * padded_len);
  std::array<float*, 3> r_curr{b_curr.data() + 0 * padded_len + kRadius,
                              b_curr.data() + 1 * padded_len + kRadius,
                              b_curr.data() + 2 * padded_len + kRadius};

  std::vector<float> b_prev(static_cast<std::size_t>(width) + lanes);
  std::vector<float> b_next(static_cast<std::size_t>(width) + lanes);
  std::vector<float> b_out(static_cast<std::size_t>(width) + lanes);

  for (int y = 1; y < height - 1; ++y) {
    if constexpr (IsF16) {
      fill_mirrored_row_fp16_to_fp32(r_curr[0] - kRadius, currp + static_cast<std::size_t>(y - 1) * curr_stride, width, kRadius);
      fill_mirrored_row_fp16_to_fp32(r_curr[1] - kRadius, currp + static_cast<std::size_t>(y) * curr_stride, width, kRadius);
      fill_mirrored_row_fp16_to_fp32(r_curr[2] - kRadius, currp + static_cast<std::size_t>(y + 1) * curr_stride, width, kRadius);

      fill_mirrored_row_fp16_to_fp32(b_prev.data(), prevp + static_cast<std::size_t>(y) * prev_stride, width, 0);
      fill_mirrored_row_fp16_to_fp32(b_next.data(), nextp + static_cast<std::size_t>(y) * next_stride, width, 0);
    } else {
      fill_mirrored_row(r_curr[0] - kRadius, currp + static_cast<std::size_t>(y - 1) * curr_stride, width, kRadius);
      fill_mirrored_row(r_curr[1] - kRadius, currp + static_cast<std::size_t>(y) * curr_stride, width, kRadius);
      fill_mirrored_row(r_curr[2] - kRadius, currp + static_cast<std::size_t>(y + 1) * curr_stride, width, kRadius);
    }

    StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      hn::Vec<decltype(d)> p, n;
      if constexpr (IsF16) {
        p = hn::LoadU(d, b_prev.data() + x);
        n = hn::LoadU(d, b_next.data() + x);
      } else {
        const std::size_t rem = static_cast<std::size_t>(width - x);
        const std::size_t count_lanes = std::min(lanes, rem);
        p = hn::LoadN(d, prevp + static_cast<std::size_t>(y) * prev_stride + x, count_lanes);
        n = hn::LoadN(d, nextp + static_cast<std::size_t>(y) * next_stride + x, count_lanes);
      }
      const auto gc = Grid3x3<decltype(d)>::load(d, r_curr[0], r_curr[1], r_curr[2], x);
      const auto c = gc.center_center;

      const auto prevnextless = hn::And(hn::Lt(p, c), hn::Lt(n, c));
      const auto prevnextmore = hn::And(hn::Gt(p, c), hn::Gt(n, c));
      const auto mask_either = hn::Or(prevnextless, prevnextmore);

      auto sum = c;
      auto count = one;

      // Temporal neighbors
      const auto p_diff = float_abs_diff<IsF16>(d, p, c);
      const auto n_diff = float_abs_diff<IsF16>(d, n, c);
      const auto p_match = hn::Le(p_diff, t_thresh);
      const auto n_match = hn::Le(n_diff, t_thresh);

      sum = float_add<IsF16>(d, sum, hn::IfThenElse(p_match, p, zero));
      count = float_add<IsF16>(d, count, hn::IfThenElse(p_match, one, zero));
      sum = float_add<IsF16>(d, sum, hn::IfThenElse(n_match, n, zero));
      count = float_add<IsF16>(d, count, hn::IfThenElse(n_match, one, zero));

      // 8 spatial neighbors
      const hn::Vec<decltype(d)> spatial_neighbors[8] = {
          gc.top_left, gc.top_center, gc.top_right,
          gc.center_left, gc.center_right,
          gc.bottom_left, gc.bottom_center, gc.bottom_right,
      };

      for (int i = 0; i < 8; ++i) {
        const auto sn = spatial_neighbors[i];
        const auto s_diff = float_abs_diff<IsF16>(d, sn, c);
        const auto s_match = hn::Le(s_diff, s_thresh);
        sum = float_add<IsF16>(d, sum, hn::IfThenElse(s_match, sn, zero));
        count = float_add<IsF16>(d, count, hn::IfThenElse(s_match, one, zero));
      }

      const auto filtered = float_div<IsF16>(d, sum, count);
      const auto res = hn::IfThenElse(mask_either, filtered, c);

      if constexpr (IsF16) {
        hn::StoreU(res, d, b_out.data() + x);
      } else {
        const std::size_t rem = static_cast<std::size_t>(width - x);
        if (rem >= lanes) hn::StoreU(res, d, dst_row + x);
        else hn::StoreN(res, d, dst_row + x, rem);
      }
    }

    if constexpr (IsF16) {
      convert_row_fp32_to_fp16(dst_row, b_out.data(), width);
    }

    dst_row[0] = currp[static_cast<std::size_t>(y) * curr_stride];
    dst_row[width - 1] = currp[static_cast<std::size_t>(y) * curr_stride + width - 1];
  }

  copy_last_n_lines(dstp, currp, static_cast<std::size_t>(width), static_cast<std::size_t>(height), dst_stride, curr_stride, 1);
}

void dispatch_fluxsmooth_t_target(DataType dtype, float temporal_threshold, const std::uint8_t* prevp,
                                  const std::uint8_t* currp, const std::uint8_t* nextp,
                                  std::uint8_t* dstp, std::size_t width, std::size_t height,
                                  std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                  std::size_t next_stride_bytes, std::size_t dst_stride_bytes) {
  if (dtype == DataType::F16) {
    temporal_threshold = fp16_to_fp32(fp32_to_fp16(temporal_threshold));
  }
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);

  if (dtype == DataType::U8) {
    const auto th = static_cast<std::uint8_t>(std::clamp(temporal_threshold, 0.0f, 255.0f));
    fluxsmooth_t_int_impl<std::uint8_t>(th, prevp, currp, nextp, dstp, w, h,
                                       prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    const auto th = static_cast<std::uint16_t>(std::clamp(temporal_threshold, 0.0f, 65535.0f));
    fluxsmooth_t_int_impl<std::uint16_t>(
        th, reinterpret_cast<const std::uint16_t*>(prevp),
        reinterpret_cast<const std::uint16_t*>(currp), reinterpret_cast<const std::uint16_t*>(nextp),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, prev_stride_bytes / 2, curr_stride_bytes / 2,
        next_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    fluxsmooth_t_float_impl<true, std::uint16_t>(
        temporal_threshold, reinterpret_cast<const std::uint16_t*>(prevp),
        reinterpret_cast<const std::uint16_t*>(currp), reinterpret_cast<const std::uint16_t*>(nextp),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, prev_stride_bytes / 2, curr_stride_bytes / 2,
        next_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    fluxsmooth_t_float_impl<false, float>(
        temporal_threshold, reinterpret_cast<const float*>(prevp),
        reinterpret_cast<const float*>(currp), reinterpret_cast<const float*>(nextp),
        reinterpret_cast<float*>(dstp), w, h, prev_stride_bytes / 4, curr_stride_bytes / 4,
        next_stride_bytes / 4, dst_stride_bytes / 4);
  }
}

void dispatch_fluxsmooth_st_target(DataType dtype, float temporal_threshold, float spatial_threshold,
                                   const std::uint8_t* prevp, const std::uint8_t* currp,
                                   const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                                   std::size_t height, std::size_t prev_stride_bytes,
                                   std::size_t curr_stride_bytes, std::size_t next_stride_bytes,
                                   std::size_t dst_stride_bytes) {
  if (dtype == DataType::F16) spatial_threshold = fp16_to_fp32(fp32_to_fp16(spatial_threshold));
  if (dtype == DataType::F16) {
    temporal_threshold = fp16_to_fp32(fp32_to_fp16(temporal_threshold));
  }
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);

  if (dtype == DataType::U8) {
    const auto t_th = static_cast<std::int32_t>(std::clamp(temporal_threshold, -2147483648.0f, 255.0f));
    const auto s_th = static_cast<std::int32_t>(std::clamp(spatial_threshold, -2147483648.0f, 255.0f));
    fluxsmooth_st_int_impl<std::uint8_t>(t_th, s_th, prevp, currp, nextp, dstp, w, h,
                                        prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    const auto t_th = static_cast<std::int32_t>(std::clamp(temporal_threshold, -2147483648.0f, 65535.0f));
    const auto s_th = static_cast<std::int32_t>(std::clamp(spatial_threshold, -2147483648.0f, 65535.0f));
    fluxsmooth_st_int_impl<std::uint16_t>(
        t_th, s_th, reinterpret_cast<const std::uint16_t*>(prevp),
        reinterpret_cast<const std::uint16_t*>(currp), reinterpret_cast<const std::uint16_t*>(nextp),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, prev_stride_bytes / 2, curr_stride_bytes / 2,
        next_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    fluxsmooth_st_float_impl<true, std::uint16_t>(
        temporal_threshold, spatial_threshold, reinterpret_cast<const std::uint16_t*>(prevp),
        reinterpret_cast<const std::uint16_t*>(currp), reinterpret_cast<const std::uint16_t*>(nextp),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, prev_stride_bytes / 2, curr_stride_bytes / 2,
        next_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    fluxsmooth_st_float_impl<false, float>(
        temporal_threshold, spatial_threshold, reinterpret_cast<const float*>(prevp),
        reinterpret_cast<const float*>(currp), reinterpret_cast<const float*>(nextp),
        reinterpret_cast<float*>(dstp), w, h, prev_stride_bytes / 4, curr_stride_bytes / 4,
        next_stride_bytes / 4, dst_stride_bytes / 4);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_fluxsmooth_t_target);
HWY_EXPORT(dispatch_fluxsmooth_st_target);

void process_fluxsmooth_t_plane(DataType dtype, float temporal_threshold, const std::uint8_t* prevp,
                                const std::uint8_t* currp, const std::uint8_t* nextp,
                                std::uint8_t* dstp, std::size_t width, std::size_t height,
                                std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                std::size_t next_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_fluxsmooth_t_target)(
      dtype, temporal_threshold, prevp, currp, nextp, dstp, width, height,
      prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
}

void process_fluxsmooth_st_plane(DataType dtype, float temporal_threshold, float spatial_threshold,
                                 const std::uint8_t* prevp, const std::uint8_t* currp,
                                 const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                                 std::size_t height, std::size_t prev_stride_bytes,
                                 std::size_t curr_stride_bytes, std::size_t next_stride_bytes,
                                 std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_fluxsmooth_st_target)(
      dtype, temporal_threshold, spatial_threshold, prevp, currp, nextp, dstp, width, height,
      prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
