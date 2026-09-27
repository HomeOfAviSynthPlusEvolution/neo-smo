#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/degrain_median.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/grid.hpp"
#include "common/float_arithmetic.hpp"
#include "common/fp16_rows.hpp"

// Avoid expanding the repeated scalar fallback neighbor checks in native MSVC.
// Hardware SIMD and other compilers retain their existing inlining.
#if defined(NEO_SMO_MSVC_OUTLINE_DGM_EMU) && HWY_TARGET == HWY_EMU128
#define NEO_SMO_DGM_NEIGHBOR_INLINE HWY_NOINLINE
#else
#define NEO_SMO_DGM_NEIGHBOR_INLINE HWY_INLINE
#endif

template <class D, class V = hn::Vec<D>>
NEO_SMO_DGM_NEIGHBOR_INLINE void check_better_neighbors_int(D d, V a, V b, V& diff, V& min_v, V& max_v) {
  const auto newdiff = hn::AbsDiff(a, b);
  const auto le = hn::Le(newdiff, diff);
  diff = hn::Min(newdiff, diff);
  min_v = hn::IfThenElse(le, hn::Min(a, b), min_v);
  max_v = hn::IfThenElse(le, hn::Max(a, b), max_v);
}

template <bool IsF16, class D, class V = hn::Vec<D>>
NEO_SMO_DGM_NEIGHBOR_INLINE void check_better_neighbors_float(D d, V a, V b, V& diff, V& min_v, V& max_v) {
  const auto newdiff = float_abs_diff<IsF16>(d, a, b);
  const auto le = hn::Le(newdiff, diff);
  diff = hn::Min(newdiff, diff);
  min_v = hn::IfThenElse(le, hn::Min(a, b), min_v);
  max_v = hn::IfThenElse(le, hn::Max(a, b), max_v);
}

template <int Mode, class D, class V = hn::Vec<D>>
NEO_SMO_DGM_NEIGHBOR_INLINE void diag_weight_int(D d, V old_pixel, V a, V b, V& old_result, V& old_weight, V pixel_max) {
  const auto np = hn::Clamp(old_pixel, hn::Min(a, b), hn::Max(a, b));
  auto w = hn::Max(hn::SaturatedSub(old_pixel, hn::Max(a, b)), hn::SaturatedSub(hn::Min(a, b), old_pixel));

  if constexpr (Mode != 5) {
    auto neighbor_diff = hn::AbsDiff(a, b);
    if constexpr (Mode == 4) {
      w = hn::SaturatedAdd(w, w);
    } else if constexpr (Mode == 2) {
      neighbor_diff = hn::SaturatedAdd(neighbor_diff, neighbor_diff);
    } else if constexpr (Mode == 1) {
      neighbor_diff = hn::SaturatedAdd(hn::SaturatedAdd(neighbor_diff, neighbor_diff),
                                      hn::SaturatedAdd(neighbor_diff, neighbor_diff));
    }
    w = hn::Min(hn::SaturatedAdd(w, neighbor_diff), pixel_max);
  }

  const auto le = hn::Le(w, old_weight);
  old_weight = hn::Min(w, old_weight);
  old_result = hn::IfThenElse(le, np, old_result);
}

template <int Mode, bool IsF16, class D, class V = hn::Vec<D>>
NEO_SMO_DGM_NEIGHBOR_INLINE void diag_weight_float(D d, V old_pixel, V a, V b, V& old_result, V& old_weight, V pixel_min, V pixel_max) {
  const auto min_ab = hn::Min(a, b);
  const auto max_ab = hn::Max(a, b);
  const auto new_pixel = hn::Clamp(old_pixel, min_ab, max_ab);

  const auto pixel_clamped_diff = hn::Max(float_sub<IsF16>(d, old_pixel, max_ab), pixel_min);
  const auto weight_diff = hn::Max(float_sub<IsF16>(d, min_ab, old_pixel), pixel_min);
  auto w = hn::Max(weight_diff, pixel_clamped_diff);

  if constexpr (Mode != 5) {
    auto neighbor_diff = float_abs_diff<IsF16>(d, a, b);
    if constexpr (Mode == 4) {
      w = float_add<IsF16>(d, w, w);
    } else if constexpr (Mode == 2) {
      neighbor_diff = float_add<IsF16>(d, neighbor_diff, neighbor_diff);
    } else if constexpr (Mode == 1) {
      const auto two_d = float_add<IsF16>(d, neighbor_diff, neighbor_diff);
      neighbor_diff = float_add<IsF16>(d, two_d, two_d);
    }
    w = hn::Min(float_add<IsF16>(d, w, neighbor_diff), pixel_max);
  }

  const auto le = hn::Le(w, old_weight);
  old_weight = hn::Min(w, old_weight);
  old_result = hn::IfThenElse(le, new_pixel, old_result);
}

#undef NEO_SMO_DGM_NEIGHBOR_INLINE

template <class D, class V = hn::Vec<D>>
HWY_INLINE V limit_pixel_correction_int(D d, V old_pixel, V new_pixel, V limit, V pixel_max) {
  (void)d;
  const auto lower = hn::SaturatedSub(old_pixel, limit);
  const auto upper = hn::Min(hn::SaturatedAdd(old_pixel, limit), pixel_max);
  return hn::Clamp(new_pixel, lower, upper);
}

template <bool IsF16, class D, class V = hn::Vec<D>>
HWY_INLINE V limit_pixel_correction_float(D d, V old_pixel, V new_pixel, V limit, V pixel_min, V pixel_max) {
  const auto lower = hn::Max(pixel_min, float_sub<IsF16>(d, old_pixel, limit));
  const auto upper = hn::Min(pixel_max, float_add<IsF16>(d, old_pixel, limit));
  return hn::Clamp(new_pixel, lower, upper);
}

template <int Mode, bool NoRow, class D, class Grid, class V = hn::Vec<D>>
// Keep native MSVC's compile-cost boundary; other compilers can keep the
// neighborhood in registers instead of passing three grids through memory.
#if defined(_MSC_VER) && !defined(__clang__)
HWY_NOINLINE
#else
HWY_INLINE
#endif
V eval_dgm_int(D d, const Grid& prev, const Grid& curr, const Grid& next,
                            V limit, V pixel_max) {
  if constexpr (Mode == 0) {
    V diff = pixel_max;
    V max_v = pixel_max;
    V min_v = hn::Zero(d);

    check_better_neighbors_int(d, next.top_left, prev.bottom_right, diff, min_v, max_v);
    check_better_neighbors_int(d, next.top_right, prev.bottom_left, diff, min_v, max_v);
    check_better_neighbors_int(d, next.bottom_left, prev.top_right, diff, min_v, max_v);
    check_better_neighbors_int(d, next.bottom_right, prev.top_left, diff, min_v, max_v);

    check_better_neighbors_int(d, next.bottom_center, prev.top_center, diff, min_v, max_v);
    check_better_neighbors_int(d, next.top_center, prev.bottom_center, diff, min_v, max_v);

    check_better_neighbors_int(d, next.center_left, prev.center_right, diff, min_v, max_v);
    check_better_neighbors_int(d, next.center_right, prev.center_left, diff, min_v, max_v);

    check_better_neighbors_int(d, next.center_center, prev.center_center, diff, min_v, max_v);

    check_better_neighbors_int(d, curr.top_left, curr.bottom_right, diff, min_v, max_v);
    check_better_neighbors_int(d, curr.top_right, curr.bottom_left, diff, min_v, max_v);

    check_better_neighbors_int(d, curr.top_center, curr.bottom_center, diff, min_v, max_v);

    if constexpr (!NoRow) {
      check_better_neighbors_int(d, curr.center_left, curr.center_right, diff, min_v, max_v);
    }

    const auto res = hn::Clamp(curr.center_center, min_v, max_v);
    return limit_pixel_correction_int(d, curr.center_center, res, limit, pixel_max);
  } else {
    V result = hn::Zero(d);
    V weight = pixel_max;
    const V old_pixel = curr.center_center;

    diag_weight_int<Mode>(d, old_pixel, curr.top_left, curr.bottom_right, result, weight, pixel_max);
    diag_weight_int<Mode>(d, old_pixel, curr.bottom_left, curr.top_right, result, weight, pixel_max);
    diag_weight_int<Mode>(d, old_pixel, curr.bottom_center, curr.top_center, result, weight, pixel_max);

    if constexpr (!NoRow) {
      diag_weight_int<Mode>(d, old_pixel, curr.center_left, curr.center_right, result, weight, pixel_max);
    }

    diag_weight_int<Mode>(d, old_pixel, next.top_left, prev.bottom_right, result, weight, pixel_max);
    diag_weight_int<Mode>(d, old_pixel, next.top_right, prev.bottom_left, result, weight, pixel_max);
    diag_weight_int<Mode>(d, old_pixel, next.bottom_left, prev.top_right, result, weight, pixel_max);
    diag_weight_int<Mode>(d, old_pixel, next.bottom_right, prev.top_left, result, weight, pixel_max);

    diag_weight_int<Mode>(d, old_pixel, next.bottom_center, prev.top_center, result, weight, pixel_max);
    diag_weight_int<Mode>(d, old_pixel, next.top_center, prev.bottom_center, result, weight, pixel_max);

    diag_weight_int<Mode>(d, old_pixel, next.center_left, prev.center_right, result, weight, pixel_max);
    diag_weight_int<Mode>(d, old_pixel, next.center_right, prev.center_left, result, weight, pixel_max);
    diag_weight_int<Mode>(d, old_pixel, next.center_center, prev.center_center, result, weight, pixel_max);

    return limit_pixel_correction_int(d, old_pixel, result, limit, pixel_max);
  }
}

template <int Mode, bool NoRow, bool IsF16, class D, class Grid, class V = hn::Vec<D>>
HWY_NOINLINE V eval_dgm_float(D d, const Grid& prev, const Grid& curr, const Grid& next,
                              V limit, V pixel_min, V pixel_max) {
  if constexpr (Mode == 0) {
    V diff = pixel_max;
    V max_v = pixel_max;
    V min_v = hn::Zero(d);

    check_better_neighbors_float<IsF16>(d, next.top_left, prev.bottom_right, diff, min_v, max_v);
    check_better_neighbors_float<IsF16>(d, next.top_right, prev.bottom_left, diff, min_v, max_v);
    check_better_neighbors_float<IsF16>(d, next.bottom_left, prev.top_right, diff, min_v, max_v);
    check_better_neighbors_float<IsF16>(d, next.bottom_right, prev.top_left, diff, min_v, max_v);

    check_better_neighbors_float<IsF16>(d, next.bottom_center, prev.top_center, diff, min_v, max_v);
    check_better_neighbors_float<IsF16>(d, next.top_center, prev.bottom_center, diff, min_v, max_v);

    check_better_neighbors_float<IsF16>(d, next.center_left, prev.center_right, diff, min_v, max_v);
    check_better_neighbors_float<IsF16>(d, next.center_right, prev.center_left, diff, min_v, max_v);

    check_better_neighbors_float<IsF16>(d, next.center_center, prev.center_center, diff, min_v, max_v);

    check_better_neighbors_float<IsF16>(d, curr.top_left, curr.bottom_right, diff, min_v, max_v);
    check_better_neighbors_float<IsF16>(d, curr.top_right, curr.bottom_left, diff, min_v, max_v);

    check_better_neighbors_float<IsF16>(d, curr.top_center, curr.bottom_center, diff, min_v, max_v);

    if constexpr (!NoRow) {
      check_better_neighbors_float<IsF16>(d, curr.center_left, curr.center_right, diff, min_v, max_v);
    }

    const auto res = hn::Clamp(curr.center_center, min_v, max_v);
    return limit_pixel_correction_float<IsF16>(d, curr.center_center, res, limit, pixel_min, pixel_max);
  } else {
    V result = hn::Zero(d);
    V weight = pixel_max;
    const V old_pixel = curr.center_center;

    diag_weight_float<Mode, IsF16>(d, old_pixel, curr.top_left, curr.bottom_right, result, weight, pixel_min, pixel_max);
    diag_weight_float<Mode, IsF16>(d, old_pixel, curr.bottom_left, curr.top_right, result, weight, pixel_min, pixel_max);
    diag_weight_float<Mode, IsF16>(d, old_pixel, curr.bottom_center, curr.top_center, result, weight, pixel_min, pixel_max);

    if constexpr (!NoRow) {
      diag_weight_float<Mode, IsF16>(d, old_pixel, curr.center_left, curr.center_right, result, weight, pixel_min, pixel_max);
    }

    diag_weight_float<Mode, IsF16>(d, old_pixel, next.top_left, prev.bottom_right, result, weight, pixel_min, pixel_max);
    diag_weight_float<Mode, IsF16>(d, old_pixel, next.top_right, prev.bottom_left, result, weight, pixel_min, pixel_max);
    diag_weight_float<Mode, IsF16>(d, old_pixel, next.bottom_left, prev.top_right, result, weight, pixel_min, pixel_max);
    diag_weight_float<Mode, IsF16>(d, old_pixel, next.bottom_right, prev.top_left, result, weight, pixel_min, pixel_max);

    diag_weight_float<Mode, IsF16>(d, old_pixel, next.bottom_center, prev.top_center, result, weight, pixel_min, pixel_max);
    diag_weight_float<Mode, IsF16>(d, old_pixel, next.top_center, prev.bottom_center, result, weight, pixel_min, pixel_max);

    diag_weight_float<Mode, IsF16>(d, old_pixel, next.center_left, prev.center_right, result, weight, pixel_min, pixel_max);
    diag_weight_float<Mode, IsF16>(d, old_pixel, next.center_right, prev.center_left, result, weight, pixel_min, pixel_max);
    diag_weight_float<Mode, IsF16>(d, old_pixel, next.center_center, prev.center_center, result, weight, pixel_min, pixel_max);

    return limit_pixel_correction_float<IsF16>(d, old_pixel, result, limit, pixel_min, pixel_max);
  }
}

template <int Mode, bool NoRow, typename T>
void degrain_median_int_impl(int bits_per_sample, T limit, bool interlaced, const T* prevp,
                             const T* currp, const T* nextp, T* dstp, int width, int height,
                             std::size_t prev_stride, std::size_t curr_stride,
                             std::size_t next_stride, std::size_t dst_stride) {
  constexpr int kRadius = 1;
  const hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + lanes;
  const auto pixel_max = hn::Set(d, get_format_maximum<T>(bits_per_sample, false));
  const auto limit_vec = hn::Set(d, limit);

  const int skip_rows = interlaced ? 2 : 1;
  if (height <= 2 * skip_rows || width < 3) {
    copy_plane(dstp, currp, width, height, dst_stride, curr_stride);
    return;
  }

  copy_first_n_lines(dstp, currp, static_cast<std::size_t>(width), dst_stride, curr_stride, skip_rows);

  std::vector<T> b_prev(3 * padded_len, 0), b_curr(3 * padded_len, 0), b_next(3 * padded_len, 0);
  std::array<T*, 3> r_prev{b_prev.data() + 0 * padded_len + kRadius, b_prev.data() + 1 * padded_len + kRadius, b_prev.data() + 2 * padded_len + kRadius};
  std::array<T*, 3> r_curr{b_curr.data() + 0 * padded_len + kRadius, b_curr.data() + 1 * padded_len + kRadius, b_curr.data() + 2 * padded_len + kRadius};
  std::array<T*, 3> r_next{b_next.data() + 0 * padded_len + kRadius, b_next.data() + 1 * padded_len + kRadius, b_next.data() + 2 * padded_len + kRadius};

  for (int y = skip_rows; y < height - skip_rows; ++y) {
    fill_mirrored_row(r_prev[0] - kRadius, prevp + static_cast<std::size_t>(y - skip_rows) * prev_stride, width, kRadius);
    fill_mirrored_row(r_prev[1] - kRadius, prevp + static_cast<std::size_t>(y) * prev_stride, width, kRadius);
    fill_mirrored_row(r_prev[2] - kRadius, prevp + static_cast<std::size_t>(y + skip_rows) * prev_stride, width, kRadius);

    fill_mirrored_row(r_curr[0] - kRadius, currp + static_cast<std::size_t>(y - skip_rows) * curr_stride, width, kRadius);
    fill_mirrored_row(r_curr[1] - kRadius, currp + static_cast<std::size_t>(y) * curr_stride, width, kRadius);
    fill_mirrored_row(r_curr[2] - kRadius, currp + static_cast<std::size_t>(y + skip_rows) * curr_stride, width, kRadius);

    fill_mirrored_row(r_next[0] - kRadius, nextp + static_cast<std::size_t>(y - skip_rows) * next_stride, width, kRadius);
    fill_mirrored_row(r_next[1] - kRadius, nextp + static_cast<std::size_t>(y) * next_stride, width, kRadius);
    fill_mirrored_row(r_next[2] - kRadius, nextp + static_cast<std::size_t>(y + skip_rows) * next_stride, width, kRadius);

    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto gp = Grid3x3<decltype(d)>::load(d, r_prev[0], r_prev[1], r_prev[2], x);
      const auto gc = Grid3x3<decltype(d)>::load(d, r_curr[0], r_curr[1], r_curr[2], x);
      const auto gn = Grid3x3<decltype(d)>::load(d, r_next[0], r_next[1], r_next[2], x);

      const auto res = eval_dgm_int<Mode, NoRow>(d, gp, gc, gn, limit_vec, pixel_max);

      const std::size_t rem = static_cast<std::size_t>(width - x);
      if (rem >= lanes) hn::StoreU(res, d, dst_row + x);
      else hn::StoreN(res, d, dst_row + x, rem);
    }

    dst_row[0] = currp[static_cast<std::size_t>(y) * curr_stride];
    dst_row[width - 1] = currp[static_cast<std::size_t>(y) * curr_stride + width - 1];
  }

  copy_last_n_lines(dstp, currp, static_cast<std::size_t>(width), static_cast<std::size_t>(height), dst_stride, curr_stride, skip_rows);
}

template <int Mode, bool NoRow, bool IsF16, typename StorageT>
void degrain_median_float_impl(float limit, bool interlaced, bool chroma, const StorageT* prevp,
                               const StorageT* currp, const StorageT* nextp, StorageT* dstp,
                               int width, int height, std::size_t prev_stride,
                               std::size_t curr_stride, std::size_t next_stride, std::size_t dst_stride) {
  using ComputeT = FloatLane<IsF16>;
  constexpr int kRadius = 1;
  const hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + lanes;
  const auto pixel_min = hn::Set(d, chroma ? -0.5f : 0.0f);
  const auto pixel_max = hn::Set(d, chroma ? 0.5f : 1.0f);
  const auto limit_vec = hn::Set(d, limit);

  const int skip_rows = interlaced ? 2 : 1;
  if (height <= 2 * skip_rows || width < 3) {
    copy_plane(dstp, currp, width, height, dst_stride, curr_stride);
    return;
  }

  copy_first_n_lines(dstp, currp, static_cast<std::size_t>(width), dst_stride, curr_stride, skip_rows);

  std::vector<ComputeT> b_prev(3 * padded_len), b_curr(3 * padded_len), b_next(3 * padded_len);
  std::vector<ComputeT> b_out(static_cast<std::size_t>(width) + lanes);
  std::array<ComputeT*, 3> r_prev{b_prev.data() + 0 * padded_len + kRadius, b_prev.data() + 1 * padded_len + kRadius, b_prev.data() + 2 * padded_len + kRadius};
  std::array<ComputeT*, 3> r_curr{b_curr.data() + 0 * padded_len + kRadius, b_curr.data() + 1 * padded_len + kRadius, b_curr.data() + 2 * padded_len + kRadius};
  std::array<ComputeT*, 3> r_next{b_next.data() + 0 * padded_len + kRadius, b_next.data() + 1 * padded_len + kRadius, b_next.data() + 2 * padded_len + kRadius};

  for (int y = skip_rows; y < height - skip_rows; ++y) {
    if constexpr (IsF16) {
      fill_mirrored_row_f16(r_prev[0] - kRadius, prevp + static_cast<std::size_t>(y - skip_rows) * prev_stride, width, kRadius);
      fill_mirrored_row_f16(r_prev[1] - kRadius, prevp + static_cast<std::size_t>(y) * prev_stride, width, kRadius);
      fill_mirrored_row_f16(r_prev[2] - kRadius, prevp + static_cast<std::size_t>(y + skip_rows) * prev_stride, width, kRadius);

      fill_mirrored_row_f16(r_curr[0] - kRadius, currp + static_cast<std::size_t>(y - skip_rows) * curr_stride, width, kRadius);
      fill_mirrored_row_f16(r_curr[1] - kRadius, currp + static_cast<std::size_t>(y) * curr_stride, width, kRadius);
      fill_mirrored_row_f16(r_curr[2] - kRadius, currp + static_cast<std::size_t>(y + skip_rows) * curr_stride, width, kRadius);

      fill_mirrored_row_f16(r_next[0] - kRadius, nextp + static_cast<std::size_t>(y - skip_rows) * next_stride, width, kRadius);
      fill_mirrored_row_f16(r_next[1] - kRadius, nextp + static_cast<std::size_t>(y) * next_stride, width, kRadius);
      fill_mirrored_row_f16(r_next[2] - kRadius, nextp + static_cast<std::size_t>(y + skip_rows) * next_stride, width, kRadius);
    } else {
      fill_mirrored_row(r_prev[0] - kRadius, prevp + static_cast<std::size_t>(y - skip_rows) * prev_stride, width, kRadius);
      fill_mirrored_row(r_prev[1] - kRadius, prevp + static_cast<std::size_t>(y) * prev_stride, width, kRadius);
      fill_mirrored_row(r_prev[2] - kRadius, prevp + static_cast<std::size_t>(y + skip_rows) * prev_stride, width, kRadius);

      fill_mirrored_row(r_curr[0] - kRadius, currp + static_cast<std::size_t>(y - skip_rows) * curr_stride, width, kRadius);
      fill_mirrored_row(r_curr[1] - kRadius, currp + static_cast<std::size_t>(y) * curr_stride, width, kRadius);
      fill_mirrored_row(r_curr[2] - kRadius, currp + static_cast<std::size_t>(y + skip_rows) * curr_stride, width, kRadius);

      fill_mirrored_row(r_next[0] - kRadius, nextp + static_cast<std::size_t>(y - skip_rows) * next_stride, width, kRadius);
      fill_mirrored_row(r_next[1] - kRadius, nextp + static_cast<std::size_t>(y) * next_stride, width, kRadius);
      fill_mirrored_row(r_next[2] - kRadius, nextp + static_cast<std::size_t>(y + skip_rows) * next_stride, width, kRadius);
    }

    StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto gp = Grid3x3<decltype(d)>::load(d, r_prev[0], r_prev[1], r_prev[2], x);
      const auto gc = Grid3x3<decltype(d)>::load(d, r_curr[0], r_curr[1], r_curr[2], x);
      const auto gn = Grid3x3<decltype(d)>::load(d, r_next[0], r_next[1], r_next[2], x);

      const auto res = eval_dgm_float<Mode, NoRow, IsF16>(d, gp, gc, gn, limit_vec, pixel_min, pixel_max);

      if constexpr (IsF16) {
        hn::StoreU(res, d, b_out.data() + x);
      } else {
        const std::size_t rem = static_cast<std::size_t>(width - x);
        if (rem >= lanes) hn::StoreU(res, d, dst_row + x);
        else hn::StoreN(res, d, dst_row + x, rem);
      }
    }

    if constexpr (IsF16) {
      store_row_f16(dst_row, b_out.data(), width);
    }
    dst_row[0] = currp[static_cast<std::size_t>(y) * curr_stride];
    dst_row[width - 1] = currp[static_cast<std::size_t>(y) * curr_stride + width - 1];
  }

  copy_last_n_lines(dstp, currp, static_cast<std::size_t>(width), static_cast<std::size_t>(height), dst_stride, curr_stride, skip_rows);
}

void dispatch_degrain_median_target(DataType dtype, int mode, float limit, bool interlaced, bool norow,
                                    bool chroma, int bits_per_sample, const std::uint8_t* prevp,
                                    const std::uint8_t* currp, const std::uint8_t* nextp,
                                    std::uint8_t* dstp, std::size_t width, std::size_t height,
                                    std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                    std::size_t next_stride_bytes, std::size_t dst_stride_bytes) {
  if (dtype == DataType::F16) {
    if constexpr (HWY_HAVE_FLOAT16) limit = fp16_to_fp32(fp32_to_fp16(limit));
  }
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);

  auto dispatch_all = [&](auto mode_tag, auto norow_tag) {
    constexpr int M = decltype(mode_tag)::value;
    constexpr bool NR = decltype(norow_tag)::value;

    if (dtype == DataType::U8) {
      const auto lim = static_cast<std::uint8_t>(std::clamp(limit, 0.0f, 255.0f));
      degrain_median_int_impl<M, NR, std::uint8_t>(bits_per_sample, lim, interlaced, prevp, currp, nextp, dstp, w, h,
                                                   prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
    } else if (dtype == DataType::U16) {
      const auto max_val = static_cast<float>((1 << bits_per_sample) - 1);
      const auto lim = static_cast<std::uint16_t>(std::clamp(limit, 0.0f, max_val));
      degrain_median_int_impl<M, NR, std::uint16_t>(
          bits_per_sample, lim, interlaced, reinterpret_cast<const std::uint16_t*>(prevp),
          reinterpret_cast<const std::uint16_t*>(currp), reinterpret_cast<const std::uint16_t*>(nextp),
          reinterpret_cast<std::uint16_t*>(dstp), w, h, prev_stride_bytes / 2, curr_stride_bytes / 2,
          next_stride_bytes / 2, dst_stride_bytes / 2);
    } else if (dtype == DataType::F16) {
      degrain_median_float_impl<M, NR, true, std::uint16_t>(
          limit, interlaced, chroma, reinterpret_cast<const std::uint16_t*>(prevp),
          reinterpret_cast<const std::uint16_t*>(currp), reinterpret_cast<const std::uint16_t*>(nextp),
          reinterpret_cast<std::uint16_t*>(dstp), w, h, prev_stride_bytes / 2, curr_stride_bytes / 2,
          next_stride_bytes / 2, dst_stride_bytes / 2);
    } else if (dtype == DataType::F32) {
      degrain_median_float_impl<M, NR, false, float>(
          limit, interlaced, chroma, reinterpret_cast<const float*>(prevp),
          reinterpret_cast<const float*>(currp), reinterpret_cast<const float*>(nextp),
          reinterpret_cast<float*>(dstp), w, h, prev_stride_bytes / 4, curr_stride_bytes / 4,
          next_stride_bytes / 4, dst_stride_bytes / 4);
    }
  };

  #define DISPATCH_MODE_NOROW(m, nr) \
    case m * 2 + (nr ? 1 : 0): \
      dispatch_all(std::integral_constant<int, m>{}, std::bool_constant<nr>{}); \
      break;

  switch (mode * 2 + (norow ? 1 : 0)) {
    DISPATCH_MODE_NOROW(0, false)
    DISPATCH_MODE_NOROW(0, true)
    DISPATCH_MODE_NOROW(1, false)
    DISPATCH_MODE_NOROW(1, true)
    DISPATCH_MODE_NOROW(2, false)
    DISPATCH_MODE_NOROW(2, true)
    DISPATCH_MODE_NOROW(3, false)
    DISPATCH_MODE_NOROW(3, true)
    DISPATCH_MODE_NOROW(4, false)
    DISPATCH_MODE_NOROW(4, true)
    DISPATCH_MODE_NOROW(5, false)
    DISPATCH_MODE_NOROW(5, true)
    default: break;
  }
  #undef DISPATCH_MODE_NOROW
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_degrain_median_target);

void process_degrain_median_plane(DataType dtype, int mode, float limit, bool interlaced, bool norow,
                                  bool chroma, int bits_per_sample, const std::uint8_t* prevp,
                                  const std::uint8_t* currp, const std::uint8_t* nextp,
                                  std::uint8_t* dstp, std::size_t width, std::size_t height,
                                  std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                  std::size_t next_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_degrain_median_target)(
      dtype, mode, limit, interlaced, norow, chroma, bits_per_sample, prevp, currp, nextp, dstp,
      width, height, prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
