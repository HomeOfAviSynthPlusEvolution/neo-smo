#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/temporal_repair.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/float_arithmetic.hpp"
#include "common/fp16_rows.hpp"

// Saturated subtraction and addition helpers for vectors
template <class D, class V = hn::Vec<D>>
HWY_INLINE V sub_sat_int(D d, V a, V b) {
  (void)d;
  return hn::SaturatedSub(a, b);
}

template <class D, class V = hn::Vec<D>>
HWY_INLINE V add_sat_int(D d, V a, V b, V max_val) {
  return hn::Min(hn::SaturatedAdd(a, b), max_val);
}

template <bool IsF16, class D, class V = hn::Vec<D>>
HWY_INLINE V sub_sat_float(D d, V a, V b, V min_val) {
  return hn::Max(min_val, float_sub<IsF16>(d, a, b));
}

template <bool IsF16, class D, class V = hn::Vec<D>>
HWY_INLINE V add_sat_float(D d, V a, V b, V max_val) {
  return hn::Min(max_val, float_add<IsF16>(d, a, b));
}

// ---------------------------------------------------------------------------
// Pure temporal modes: 0 and 4
// ---------------------------------------------------------------------------
template <typename T>
void temporal_repair_pt_int_impl(int mode, int bits_per_sample, const T* srcp, const T* prevp,
                                 const T* currp, const T* nextp, T* dstp, int width, int height,
                                 std::size_t src_stride, std::size_t prev_stride,
                                 std::size_t curr_stride, std::size_t next_stride, std::size_t dst_stride) {
  const hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);
  const auto format_max = hn::Set(d, get_format_maximum<T>(bits_per_sample, false));

  for (int y = 0; y < height; ++y) {
    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const std::size_t rem = static_cast<std::size_t>(width - x);
      const std::size_t count = std::min(lanes, rem);

      const auto s = hn::LoadN(d, srcp + static_cast<std::size_t>(y) * src_stride + x, count);
      const auto p = hn::LoadN(d, prevp + static_cast<std::size_t>(y) * prev_stride + x, count);
      const auto c = hn::LoadN(d, currp + static_cast<std::size_t>(y) * curr_stride + x, count);
      const auto n = hn::LoadN(d, nextp + static_cast<std::size_t>(y) * next_stride + x, count);

      auto res = s;
      if (mode == 0) {
        const auto mn = hn::Min(hn::Min(p, c), n);
        const auto mx = hn::Max(hn::Max(p, c), n);
        res = hn::Clamp(s, mn, mx);
      } else { // mode == 4
        const auto brightest_neighbor = hn::Max(p, n);
        const auto darkest_neighbor = hn::Min(p, n);

        const auto diff_curr_darkest = sub_sat_int(d, c, darkest_neighbor);
        const auto darkest_plus_weighted_diff = add_sat_int(
            d, add_sat_int(d, diff_curr_darkest, diff_curr_darkest, format_max), darkest_neighbor, format_max);

        const auto diff_curr_brightest = sub_sat_int(d, brightest_neighbor, c);
        const auto brightest_minus_weighted_diff =
            sub_sat_int(d, brightest_neighbor, add_sat_int(d, diff_curr_brightest, diff_curr_brightest, format_max));

        const auto upper = hn::Min(darkest_plus_weighted_diff, brightest_neighbor);
        const auto lower = hn::Max(brightest_minus_weighted_diff, darkest_neighbor);

        const auto skip_cond = hn::Or(hn::Eq(darkest_neighbor, upper), hn::Eq(brightest_neighbor, lower));
        res = hn::IfThenElse(skip_cond, c, hn::Clamp(s, lower, upper));
      }

      if (rem >= lanes) hn::StoreU(res, d, dst_row + x);
      else hn::StoreN(res, d, dst_row + x, rem);
    }
  }
}

template <bool IsF16, typename StorageT>
void temporal_repair_pt_float_impl(int mode, bool chroma, const StorageT* srcp, const StorageT* prevp,
                                   const StorageT* currp, const StorageT* nextp, StorageT* dstp,
                                   int width, int height, std::size_t src_stride,
                                   std::size_t prev_stride, std::size_t curr_stride,
                                   std::size_t next_stride, std::size_t dst_stride) {
  using ComputeT = FloatLane<IsF16>;
  const hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);
  const auto format_min = hn::Set(d, chroma ? -0.5f : 0.0f);
  const auto format_max = hn::Set(d, chroma ? 0.5f : 1.0f);

  if constexpr (IsF16) {
    for (int y = 0; y < height; ++y) {
      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto count_lanes = std::min(lanes, static_cast<std::size_t>(width) - x);
        const auto s = load_f16(d, srcp + static_cast<std::size_t>(y) * src_stride + x, count_lanes);
        const auto p = load_f16(d, prevp + static_cast<std::size_t>(y) * prev_stride + x, count_lanes);
        const auto c = load_f16(d, currp + static_cast<std::size_t>(y) * curr_stride + x, count_lanes);
        const auto n = load_f16(d, nextp + static_cast<std::size_t>(y) * next_stride + x, count_lanes);

        auto res = s;
        if (mode == 0) {
          const auto mn = hn::Min(hn::Min(p, c), n);
          const auto mx = hn::Max(hn::Max(p, c), n);
          res = hn::Clamp(s, mn, mx);
        } else {
          const auto brightest_neighbor = hn::Max(p, n);
          const auto darkest_neighbor = hn::Min(p, n);

          const auto diff_curr_darkest = sub_sat_float<true>(d, c, darkest_neighbor, hn::Zero(d));
          const auto darkest_plus_weighted_diff = add_sat_float<true>(
              d, add_sat_float<true>(d, diff_curr_darkest, diff_curr_darkest, format_max), darkest_neighbor, format_max);

          const auto diff_curr_brightest = sub_sat_float<true>(d, brightest_neighbor, c, hn::Zero(d));
          const auto brightest_minus_weighted_diff = sub_sat_float<true>(
              d, brightest_neighbor, add_sat_float<true>(d, diff_curr_brightest, diff_curr_brightest, format_max), format_min);

          auto upper = hn::Min(darkest_plus_weighted_diff, brightest_neighbor);
          const auto lower = hn::Max(brightest_minus_weighted_diff, darkest_neighbor);
          upper = hn::Max(upper, lower);

          const auto skip_cond = hn::Or(hn::Eq(darkest_neighbor, upper), hn::Eq(brightest_neighbor, lower));
          res = hn::IfThenElse(skip_cond, c, hn::Clamp(s, lower, upper));
        }
        store_f16(d, res, dst_row + x, count_lanes);
      }
    }
  } else {
    for (int y = 0; y < height; ++y) {
      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const std::size_t rem = static_cast<std::size_t>(width - x);
        const std::size_t count = std::min(lanes, rem);

        const auto s = hn::LoadN(d, srcp + static_cast<std::size_t>(y) * src_stride + x, count);
        const auto p = hn::LoadN(d, prevp + static_cast<std::size_t>(y) * prev_stride + x, count);
        const auto c = hn::LoadN(d, currp + static_cast<std::size_t>(y) * curr_stride + x, count);
        const auto n = hn::LoadN(d, nextp + static_cast<std::size_t>(y) * next_stride + x, count);

        auto res = s;
        if (mode == 0) {
          const auto mn = hn::Min(hn::Min(p, c), n);
          const auto mx = hn::Max(hn::Max(p, c), n);
          res = hn::Clamp(s, mn, mx);
        } else {
          const auto brightest_neighbor = hn::Max(p, n);
          const auto darkest_neighbor = hn::Min(p, n);

          const auto diff_curr_darkest = sub_sat_float<false>(d, c, darkest_neighbor, hn::Zero(d));
          const auto darkest_plus_weighted_diff = add_sat_float<false>(
              d, add_sat_float<false>(d, diff_curr_darkest, diff_curr_darkest, format_max), darkest_neighbor, format_max);

          const auto diff_curr_brightest = sub_sat_float<false>(d, brightest_neighbor, c, hn::Zero(d));
          const auto brightest_minus_weighted_diff = sub_sat_float<false>(
              d, brightest_neighbor, add_sat_float<false>(d, diff_curr_brightest, diff_curr_brightest, format_max), format_min);

          auto upper = hn::Min(darkest_plus_weighted_diff, brightest_neighbor);
          const auto lower = hn::Max(brightest_minus_weighted_diff, darkest_neighbor);
          upper = hn::Max(upper, lower);

          const auto skip_cond = hn::Or(hn::Eq(darkest_neighbor, upper), hn::Eq(brightest_neighbor, lower));
          res = hn::IfThenElse(skip_cond, c, hn::Clamp(s, lower, upper));
        }

        if (rem >= lanes) hn::StoreU(res, d, dst_row + x);
        else hn::StoreN(res, d, dst_row + x, rem);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Spatio-temporal modes: 1, 2, 3
// ---------------------------------------------------------------------------
template <typename T>
void temporal_repair_st_int_impl(int mode, int bits_per_sample, const T* srcp, const T* prevp,
                                 const T* currp, const T* nextp, T* dstp, int width, int height,
                                 std::size_t src_stride, std::size_t prev_stride,
                                 std::size_t curr_stride, std::size_t next_stride, std::size_t dst_stride) {
  constexpr int kRadius = 1;
  hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + lanes;
  const auto format_max = hn::Set(d, get_format_maximum<T>(bits_per_sample, false));

  std::vector<T> b_prev(3 * padded_len, 0), b_curr(3 * padded_len, 0), b_next(3 * padded_len, 0);
  std::array<T*, 3> r_prev{b_prev.data() + 0 * padded_len + kRadius, b_prev.data() + 1 * padded_len + kRadius, b_prev.data() + 2 * padded_len + kRadius};
  std::array<T*, 3> r_curr{b_curr.data() + 0 * padded_len + kRadius, b_curr.data() + 1 * padded_len + kRadius, b_curr.data() + 2 * padded_len + kRadius};
  std::array<T*, 3> r_next{b_next.data() + 0 * padded_len + kRadius, b_next.data() + 1 * padded_len + kRadius, b_next.data() + 2 * padded_len + kRadius};

  std::array<int, 3> cached_y{-1, -1, -1};
  for (int y = 0; y < height; ++y) {
    for (int dy = -1; dy <= 1; ++dy) {
      const auto my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      const auto slot = my % 3;
      auto* p = b_prev.data() + slot * padded_len + kRadius;
      auto* c = b_curr.data() + slot * padded_len + kRadius;
      auto* n = b_next.data() + slot * padded_len + kRadius;
      if (cached_y[slot] != static_cast<int>(my)) {
        fill_mirrored_row(p - kRadius, prevp + my * prev_stride, width, kRadius);
        fill_mirrored_row(c - kRadius, currp + my * curr_stride, width, kRadius);
        fill_mirrored_row(n - kRadius, nextp + my * next_stride, width, kRadius);
        cached_y[slot] = static_cast<int>(my);
      }
      r_prev[dy + 1] = p; r_curr[dy + 1] = c; r_next[dy + 1] = n;
    }

    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const std::size_t rem = static_cast<std::size_t>(width - x);
      const std::size_t count = std::min(lanes, rem);

      const auto s = hn::LoadN(d, srcp + static_cast<std::size_t>(y) * src_stride + x, count);

      hn::Vec<decltype(d)> gp[9], gc[9], gn[9];
      int idx = 0;
      for (int ky = 0; ky < 3; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          gp[idx] = hn::LoadU(d, r_prev[static_cast<std::size_t>(ky)] + x + kx);
          gc[idx] = hn::LoadU(d, r_curr[static_cast<std::size_t>(ky)] + x + kx);
          gn[idx] = hn::LoadU(d, r_next[static_cast<std::size_t>(ky)] + x + kx);
          ++idx;
        }
      }

      auto res = s;
      if (mode == 1) {
        auto brightest_diff_max = hn::Zero(d);
        auto darkest_diff_max = hn::Zero(d);
        for (int i = 0; i < 9; ++i) {
          if (i == 4) continue; // skip center
          const auto b_diff = sub_sat_int(d, hn::Max(gp[i], gn[i]), gc[i]);
          const auto d_diff = sub_sat_int(d, gc[i], hn::Min(gp[i], gn[i]));
          brightest_diff_max = hn::Max(brightest_diff_max, b_diff);
          darkest_diff_max = hn::Max(darkest_diff_max, d_diff);
        }
        const auto brightest_curr_diff = add_sat_int(d, gc[4], brightest_diff_max, format_max);
        const auto darkest_curr_diff = sub_sat_int(d, gc[4], darkest_diff_max);

        const auto mx = hn::Max(hn::Max(brightest_curr_diff, gp[4]), gn[4]);
        const auto mn = hn::Min(hn::Min(darkest_curr_diff, gp[4]), gn[4]);
        res = hn::Clamp(s, mn, mx);
      } else if (mode == 2) {
        auto brightest_diff_max = hn::Zero(d);
        auto darkest_diff_max = hn::Zero(d);
        for (int i = 0; i < 9; ++i) {
          const auto b_diff = sub_sat_int(d, hn::Max(gp[i], gn[i]), gc[i]);
          const auto d_diff = sub_sat_int(d, gc[i], hn::Min(gp[i], gn[i]));
          brightest_diff_max = hn::Max(brightest_diff_max, b_diff);
          darkest_diff_max = hn::Max(darkest_diff_max, d_diff);
        }
        const auto diff_max = hn::Max(brightest_diff_max, darkest_diff_max);
        const auto curr_diff_upper = add_sat_int(d, gc[4], diff_max, format_max);
        const auto curr_diff_lower = sub_sat_int(d, gc[4], diff_max);
        res = hn::Clamp(s, curr_diff_lower, curr_diff_upper);
      } else if (mode == 3) {
        auto prev_diff_max = hn::Zero(d);
        auto next_diff_max = hn::Zero(d);
        for (int i = 0; i < 9; ++i) {
          const auto pdiff = hn::AbsDiff(gc[i], gp[i]);
          const auto ndiff = hn::AbsDiff(gc[i], gn[i]);
          prev_diff_max = hn::Max(prev_diff_max, pdiff);
          next_diff_max = hn::Max(next_diff_max, ndiff);
        }
        const auto diff_min = hn::Min(prev_diff_max, next_diff_max);
        const auto curr_diff_upper = add_sat_int(d, gc[4], diff_min, format_max);
        const auto curr_diff_lower = sub_sat_int(d, gc[4], diff_min);
        res = hn::Clamp(s, curr_diff_lower, curr_diff_upper);
      }

      if (rem >= lanes) hn::StoreU(res, d, dst_row + x);
      else hn::StoreN(res, d, dst_row + x, rem);
    }
  }
}

template <bool IsF16, typename StorageT>
void temporal_repair_st_float_impl(int mode, bool chroma, const StorageT* srcp, const StorageT* prevp,
                                   const StorageT* currp, const StorageT* nextp, StorageT* dstp,
                                   int width, int height, std::size_t src_stride,
                                   std::size_t prev_stride, std::size_t curr_stride,
                                   std::size_t next_stride, std::size_t dst_stride) {
  using ComputeT = FloatLane<IsF16>;
  constexpr int kRadius = 1;
  const hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + lanes;
  const auto format_min = hn::Set(d, chroma ? -0.5f : 0.0f);
  const auto format_max = hn::Set(d, chroma ? 0.5f : 1.0f);

  std::vector<ComputeT> b_prev(3 * padded_len), b_curr(3 * padded_len), b_next(3 * padded_len);
  std::vector<ComputeT> b_src(IsF16 ? static_cast<std::size_t>(width) + lanes : 0),
      b_out(IsF16 ? static_cast<std::size_t>(width) + lanes : 0);
  std::array<ComputeT*, 3> r_prev{b_prev.data() + 0 * padded_len + kRadius, b_prev.data() + 1 * padded_len + kRadius, b_prev.data() + 2 * padded_len + kRadius};
  std::array<ComputeT*, 3> r_curr{b_curr.data() + 0 * padded_len + kRadius, b_curr.data() + 1 * padded_len + kRadius, b_curr.data() + 2 * padded_len + kRadius};
  std::array<ComputeT*, 3> r_next{b_next.data() + 0 * padded_len + kRadius, b_next.data() + 1 * padded_len + kRadius, b_next.data() + 2 * padded_len + kRadius};

  std::array<int, 3> cached_y{-1, -1, -1};
  for (int y = 0; y < height; ++y) {
    for (int dy = -1; dy <= 1; ++dy) {
      const auto my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      const auto slot = my % 3;
      auto* p = b_prev.data() + slot * padded_len + kRadius;
      auto* c = b_curr.data() + slot * padded_len + kRadius;
      auto* n = b_next.data() + slot * padded_len + kRadius;
      if (cached_y[slot] != static_cast<int>(my)) {
        if constexpr (IsF16) {
          fill_mirrored_row_f16(p - kRadius, prevp + my * prev_stride, width, kRadius);
          fill_mirrored_row_f16(c - kRadius, currp + my * curr_stride, width, kRadius);
          fill_mirrored_row_f16(n - kRadius, nextp + my * next_stride, width, kRadius);
        } else {
          fill_mirrored_row(p - kRadius, prevp + my * prev_stride, width, kRadius);
          fill_mirrored_row(c - kRadius, currp + my * curr_stride, width, kRadius);
          fill_mirrored_row(n - kRadius, nextp + my * next_stride, width, kRadius);
        }
        cached_y[slot] = static_cast<int>(my);
      }
      r_prev[dy + 1] = p; r_curr[dy + 1] = c; r_next[dy + 1] = n;
    }

    if constexpr (IsF16) {
      fill_mirrored_row_f16(b_src.data(), srcp + static_cast<std::size_t>(y) * src_stride, width, 0);
    }

    StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      hn::Vec<decltype(d)> s;
      if constexpr (IsF16) {
        s = hn::LoadU(d, b_src.data() + x);
      } else {
        const std::size_t rem = static_cast<std::size_t>(width - x);
        s = hn::LoadN(d, srcp + static_cast<std::size_t>(y) * src_stride + x, std::min(lanes, rem));
      }

      hn::Vec<decltype(d)> gp[9], gc[9], gn[9];
      int idx = 0;
      for (int ky = 0; ky < 3; ++ky) {
        for (int kx = -1; kx <= 1; ++kx) {
          gp[idx] = hn::LoadU(d, r_prev[static_cast<std::size_t>(ky)] + x + kx);
          gc[idx] = hn::LoadU(d, r_curr[static_cast<std::size_t>(ky)] + x + kx);
          gn[idx] = hn::LoadU(d, r_next[static_cast<std::size_t>(ky)] + x + kx);
          ++idx;
        }
      }

      auto res = s;
      if (mode == 1) {
        auto brightest_diff_max = hn::Zero(d);
        auto darkest_diff_max = hn::Zero(d);
        for (int i = 0; i < 9; ++i) {
          if (i == 4) continue;
          const auto b_diff = sub_sat_float<IsF16>(d, hn::Max(gp[i], gn[i]), gc[i], hn::Zero(d));
          const auto d_diff = sub_sat_float<IsF16>(d, gc[i], hn::Min(gp[i], gn[i]), hn::Zero(d));
          brightest_diff_max = hn::Max(brightest_diff_max, b_diff);
          darkest_diff_max = hn::Max(darkest_diff_max, d_diff);
        }
        const auto brightest_curr_diff = add_sat_float<IsF16>(d, gc[4], brightest_diff_max, format_max);
        const auto darkest_curr_diff = sub_sat_float<IsF16>(d, gc[4], darkest_diff_max, format_min);

        const auto mx = hn::Max(hn::Max(brightest_curr_diff, gp[4]), gn[4]);
        const auto mn = hn::Min(hn::Min(darkest_curr_diff, gp[4]), gn[4]);
        res = hn::Clamp(s, mn, mx);
      } else if (mode == 2) {
        auto brightest_diff_max = hn::Zero(d);
        auto darkest_diff_max = hn::Zero(d);
        for (int i = 0; i < 9; ++i) {
          const auto b_diff = sub_sat_float<IsF16>(d, hn::Max(gp[i], gn[i]), gc[i], hn::Zero(d));
          const auto d_diff = sub_sat_float<IsF16>(d, gc[i], hn::Min(gp[i], gn[i]), hn::Zero(d));
          brightest_diff_max = hn::Max(brightest_diff_max, b_diff);
          darkest_diff_max = hn::Max(darkest_diff_max, d_diff);
        }
        const auto diff_max = hn::Max(brightest_diff_max, darkest_diff_max);
        auto curr_diff_upper = add_sat_float<IsF16>(d, gc[4], diff_max, format_max);
        const auto curr_diff_lower = sub_sat_float<IsF16>(d, gc[4], diff_max, format_min);
        curr_diff_upper = hn::Max(curr_diff_upper, curr_diff_lower);
        res = hn::Clamp(s, curr_diff_lower, curr_diff_upper);
      } else if (mode == 3) {
        auto prev_diff_max = hn::Zero(d);
        auto next_diff_max = hn::Zero(d);
        for (int i = 0; i < 9; ++i) {
          const auto pdiff = float_abs_diff<IsF16>(d, gc[i], gp[i]);
          const auto ndiff = float_abs_diff<IsF16>(d, gc[i], gn[i]);
          prev_diff_max = hn::Max(prev_diff_max, pdiff);
          next_diff_max = hn::Max(next_diff_max, ndiff);
        }
        const auto diff_min = hn::Min(prev_diff_max, next_diff_max);
        auto curr_diff_upper = add_sat_float<IsF16>(d, gc[4], diff_min, format_max);
        const auto curr_diff_lower = sub_sat_float<IsF16>(d, gc[4], diff_min, format_min);
        curr_diff_upper = hn::Max(curr_diff_upper, curr_diff_lower);
        res = hn::Clamp(s, curr_diff_lower, curr_diff_upper);
      }

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
  }
}

void dispatch_temporal_repair_target(DataType dtype, int mode, bool chroma, int bits_per_sample,
                                     const std::uint8_t* srcp, const std::uint8_t* prevp,
                                     const std::uint8_t* currp, const std::uint8_t* nextp,
                                     std::uint8_t* dstp, std::size_t width, std::size_t height,
                                     std::size_t src_stride_bytes, std::size_t prev_stride_bytes,
                                     std::size_t curr_stride_bytes, std::size_t next_stride_bytes,
                                     std::size_t dst_stride_bytes) {
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);

  if (mode == 0 || mode == 4) {
    if (dtype == DataType::U8) {
      temporal_repair_pt_int_impl<std::uint8_t>(
          mode, bits_per_sample, srcp, prevp, currp, nextp, dstp, w, h,
          src_stride_bytes, prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
    } else if (dtype == DataType::U16) {
      temporal_repair_pt_int_impl<std::uint16_t>(
          mode, bits_per_sample, reinterpret_cast<const std::uint16_t*>(srcp),
          reinterpret_cast<const std::uint16_t*>(prevp), reinterpret_cast<const std::uint16_t*>(currp),
          reinterpret_cast<const std::uint16_t*>(nextp), reinterpret_cast<std::uint16_t*>(dstp),
          w, h, src_stride_bytes / 2, prev_stride_bytes / 2, curr_stride_bytes / 2,
          next_stride_bytes / 2, dst_stride_bytes / 2);
    } else if (dtype == DataType::F16) {
      temporal_repair_pt_float_impl<true, std::uint16_t>(
          mode, chroma, reinterpret_cast<const std::uint16_t*>(srcp),
          reinterpret_cast<const std::uint16_t*>(prevp), reinterpret_cast<const std::uint16_t*>(currp),
          reinterpret_cast<const std::uint16_t*>(nextp), reinterpret_cast<std::uint16_t*>(dstp),
          w, h, src_stride_bytes / 2, prev_stride_bytes / 2, curr_stride_bytes / 2,
          next_stride_bytes / 2, dst_stride_bytes / 2);
    } else if (dtype == DataType::F32) {
      temporal_repair_pt_float_impl<false, float>(
          mode, chroma, reinterpret_cast<const float*>(srcp),
          reinterpret_cast<const float*>(prevp), reinterpret_cast<const float*>(currp),
          reinterpret_cast<const float*>(nextp), reinterpret_cast<float*>(dstp),
          w, h, src_stride_bytes / 4, prev_stride_bytes / 4, curr_stride_bytes / 4,
          next_stride_bytes / 4, dst_stride_bytes / 4);
    }
  } else { // mode 1, 2, 3
    if (dtype == DataType::U8) {
      temporal_repair_st_int_impl<std::uint8_t>(
          mode, bits_per_sample, srcp, prevp, currp, nextp, dstp, w, h,
          src_stride_bytes, prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
    } else if (dtype == DataType::U16) {
      temporal_repair_st_int_impl<std::uint16_t>(
          mode, bits_per_sample, reinterpret_cast<const std::uint16_t*>(srcp),
          reinterpret_cast<const std::uint16_t*>(prevp), reinterpret_cast<const std::uint16_t*>(currp),
          reinterpret_cast<const std::uint16_t*>(nextp), reinterpret_cast<std::uint16_t*>(dstp),
          w, h, src_stride_bytes / 2, prev_stride_bytes / 2, curr_stride_bytes / 2,
          next_stride_bytes / 2, dst_stride_bytes / 2);
    } else if (dtype == DataType::F16) {
      temporal_repair_st_float_impl<true, std::uint16_t>(
          mode, chroma, reinterpret_cast<const std::uint16_t*>(srcp),
          reinterpret_cast<const std::uint16_t*>(prevp), reinterpret_cast<const std::uint16_t*>(currp),
          reinterpret_cast<const std::uint16_t*>(nextp), reinterpret_cast<std::uint16_t*>(dstp),
          w, h, src_stride_bytes / 2, prev_stride_bytes / 2, curr_stride_bytes / 2,
          next_stride_bytes / 2, dst_stride_bytes / 2);
    } else if (dtype == DataType::F32) {
      temporal_repair_st_float_impl<false, float>(
          mode, chroma, reinterpret_cast<const float*>(srcp),
          reinterpret_cast<const float*>(prevp), reinterpret_cast<const float*>(currp),
          reinterpret_cast<const float*>(nextp), reinterpret_cast<float*>(dstp),
          w, h, src_stride_bytes / 4, prev_stride_bytes / 4, curr_stride_bytes / 4,
          next_stride_bytes / 4, dst_stride_bytes / 4);
    }
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_temporal_repair_target);

void process_temporal_repair_plane(DataType dtype, int mode, bool chroma, int bits_per_sample,
                                   const std::uint8_t* srcp, const std::uint8_t* prevp,
                                   const std::uint8_t* currp, const std::uint8_t* nextp,
                                   std::uint8_t* dstp, std::size_t width, std::size_t height,
                                   std::size_t src_stride_bytes, std::size_t prev_stride_bytes,
                                   std::size_t curr_stride_bytes, std::size_t next_stride_bytes,
                                   std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_temporal_repair_target)(
      dtype, mode, chroma, bits_per_sample, srcp, prevp, currp, nextp, dstp, width, height,
      src_stride_bytes, prev_stride_bytes, curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
