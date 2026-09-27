#include "kernels/repair_dispatch.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/repair_int.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/grid.hpp"
#include "common/rank_clamp.hpp"

template <class V>
HWY_INLINE V repair_abs_diff32(V a, V b) {
#if defined(NEO_SMO_MSVC_REPAIR_ABS_DIFF)
  // Reduce MSVC optimization cost. Both operands are U8/U16 samples or clamps
  // within that range, so their difference fits in [-65535, 65535], tails included.
  return hn::Abs(hn::Sub(a, b));
#else
  return hn::AbsDiff(a, b);
#endif
}

// U8 intermediates fit signed 16-bit lanes; U16 requires signed 32-bit lanes.
template <class D, class V = hn::Vec<D>>
#if defined(_MSC_VER) && !defined(__clang__)
HWY_NOINLINE
#else
HWY_INLINE
#endif
V eval_repair_int(D d, int mode, V src, const Grid3x3<D>& g, std::int32_t type_max) {
  const V c = g.center_center;
  const V zero = hn::Zero(d);
  const V vmax = hn::Set(d, type_max);

  switch (mode) {
    case 1: {
      const V mn = g.min_with_center(d);
      const V mx = g.max_with_center(d);
      return hn::Clamp(src, mn, mx);
    }
    case 2:
    case 3:
    case 4: {
      V a[9];
      g.sort_with_center(d, a);
      const int lo = mode - 1;
      const int hi = 9 - mode;
      return hn::Clamp(src, a[lo], a[hi]);
    }
    case 5:
    case 6:
    case 7:
    case 8:
    case 9: {
      const auto s = g.min_max_opposites_with_center(d);
      const V clamp1 = hn::Clamp(src, s.min1, s.max1);
      const V clamp2 = hn::Clamp(src, s.min2, s.max2);
      const V clamp3 = hn::Clamp(src, s.min3, s.max3);
      const V clamp4 = hn::Clamp(src, s.min4, s.max4);

      V c1, c2, c3, c4;
      if (mode == 5) {
        c1 = repair_abs_diff32(src, clamp1);
        c2 = repair_abs_diff32(src, clamp2);
        c3 = repair_abs_diff32(src, clamp3);
        c4 = repair_abs_diff32(src, clamp4);
      } else {
        const V d1 = hn::Sub(s.max1, s.min1);
        const V d2 = hn::Sub(s.max2, s.min2);
        const V d3 = hn::Sub(s.max3, s.min3);
        const V d4 = hn::Sub(s.max4, s.min4);
        if (mode == 6) {
          c1 = hn::Min(hn::Add(hn::ShiftLeft<1>(repair_abs_diff32(src, clamp1)), d1), vmax);
          c2 = hn::Min(hn::Add(hn::ShiftLeft<1>(repair_abs_diff32(src, clamp2)), d2), vmax);
          c3 = hn::Min(hn::Add(hn::ShiftLeft<1>(repair_abs_diff32(src, clamp3)), d3), vmax);
          c4 = hn::Min(hn::Add(hn::ShiftLeft<1>(repair_abs_diff32(src, clamp4)), d4), vmax);
        } else if (mode == 7) {
          c1 = hn::Add(repair_abs_diff32(src, clamp1), d1);
          c2 = hn::Add(repair_abs_diff32(src, clamp2), d2);
          c3 = hn::Add(repair_abs_diff32(src, clamp3), d3);
          c4 = hn::Add(repair_abs_diff32(src, clamp4), d4);
        } else if (mode == 8) {
          c1 = hn::Clamp(hn::Add(repair_abs_diff32(src, clamp1), hn::ShiftLeft<1>(d1)), zero, vmax);
          c2 = hn::Clamp(hn::Add(repair_abs_diff32(src, clamp2), hn::ShiftLeft<1>(d2)), zero, vmax);
          c3 = hn::Clamp(hn::Add(repair_abs_diff32(src, clamp3), hn::ShiftLeft<1>(d3)), zero, vmax);
          c4 = hn::Clamp(hn::Add(repair_abs_diff32(src, clamp4), hn::ShiftLeft<1>(d4)), zero, vmax);
        } else { // mode == 9
          c1 = d1;
          c2 = d2;
          c3 = d3;
          c4 = d4;
        }
      }

      const V mindiff = hn::Min(hn::Min(c1, c2), hn::Min(c3, c4));
      V res = clamp1;
      res = hn::IfThenElse(hn::Eq(mindiff, c3), clamp3, res);
      res = hn::IfThenElse(hn::Eq(mindiff, c2), clamp2, res);
      res = hn::IfThenElse(hn::Eq(mindiff, c4), clamp4, res);
      return res;
    }
    case 10: {
      const V d1 = repair_abs_diff32(src, g.top_left);
      const V d2 = repair_abs_diff32(src, g.top_center);
      const V d3 = repair_abs_diff32(src, g.top_right);
      const V d4 = repair_abs_diff32(src, g.center_left);
      const V d5 = repair_abs_diff32(src, g.center_right);
      const V d6 = repair_abs_diff32(src, g.bottom_left);
      const V d7 = repair_abs_diff32(src, g.bottom_center);
      const V d8 = repair_abs_diff32(src, g.bottom_right);
      const V dc = repair_abs_diff32(src, c);
      const V mindiff = hn::Min(hn::Min(hn::Min(d1, d2), hn::Min(d3, d4)),
                                hn::Min(hn::Min(d5, d6), hn::Min(hn::Min(d7, d8), dc)));
      V res = g.center_left;
      res = hn::IfThenElse(hn::Eq(mindiff, dc), c, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d5), g.center_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d1), g.top_left, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d3), g.top_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), g.top_center, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d6), g.bottom_left, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d8), g.bottom_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d7), g.bottom_center, res);
      return res;
    }
    case 11: {
      const V mn = g.min_with_center(d);
      const V mx = g.max_with_center(d);
      return hn::Clamp(src, mn, mx);
    }
    case 12:
    case 13:
    case 14: {
      V a[8];
      g.sort_without_center(d, a);
      const int lo = mode - 11;
      const int hi = 18 - mode;
      const V mn = hn::Min(a[lo], c);
      const V mx = hn::Max(a[hi], c);
      return hn::Clamp(src, mn, mx);
    }
    case 15:
    case 16: {
      const auto s = g.min_max_opposites_without_center(d);
      V c1, c2, c3, c4;
      if (mode == 15) {
        c1 = repair_abs_diff32(c, hn::Clamp(c, s.min1, s.max1));
        c2 = repair_abs_diff32(c, hn::Clamp(c, s.min2, s.max2));
        c3 = repair_abs_diff32(c, hn::Clamp(c, s.min3, s.max3));
        c4 = repair_abs_diff32(c, hn::Clamp(c, s.min4, s.max4));
      } else {
        const V d1 = hn::Sub(s.max1, s.min1);
        const V d2 = hn::Sub(s.max2, s.min2);
        const V d3 = hn::Sub(s.max3, s.min3);
        const V d4 = hn::Sub(s.max4, s.min4);
        c1 = hn::Min(hn::Add(hn::ShiftLeft<1>(repair_abs_diff32(c, hn::Clamp(c, s.min1, s.max1))), d1), vmax);
        c2 = hn::Min(hn::Add(hn::ShiftLeft<1>(repair_abs_diff32(c, hn::Clamp(c, s.min2, s.max2))), d2), vmax);
        c3 = hn::Min(hn::Add(hn::ShiftLeft<1>(repair_abs_diff32(c, hn::Clamp(c, s.min3, s.max3))), d3), vmax);
        c4 = hn::Min(hn::Add(hn::ShiftLeft<1>(repair_abs_diff32(c, hn::Clamp(c, s.min4, s.max4))), d4), vmax);
      }
      const V mindiff = hn::Min(hn::Min(c1, c2), hn::Min(c3, c4));
      V pmin = s.min1;
      V pmax = s.max1;
      pmin = hn::IfThenElse(hn::Eq(mindiff, c3), s.min3, pmin);
      pmax = hn::IfThenElse(hn::Eq(mindiff, c3), s.max3, pmax);
      pmin = hn::IfThenElse(hn::Eq(mindiff, c2), s.min2, pmin);
      pmax = hn::IfThenElse(hn::Eq(mindiff, c2), s.max2, pmax);
      pmin = hn::IfThenElse(hn::Eq(mindiff, c4), s.min4, pmin);
      pmax = hn::IfThenElse(hn::Eq(mindiff, c4), s.max4, pmax);
      const V mn = hn::Min(pmin, c);
      const V mx = hn::Max(pmax, c);
      return hn::Clamp(src, mn, mx);
    }
    case 17: {
      const auto s = g.min_max_opposites_without_center(d);
      const V l = hn::Max(hn::Max(s.min1, s.min2), hn::Max(s.min3, s.min4));
      const V u = hn::Min(hn::Min(s.max1, s.max2), hn::Min(s.max3, s.max4));
      const V mn = hn::Min(hn::Min(l, u), c);
      const V mx = hn::Max(hn::Max(l, u), c);
      return hn::Clamp(src, mn, mx);
    }
    case 18: {
      const V d1 = hn::Max(repair_abs_diff32(c, g.top_left), repair_abs_diff32(c, g.bottom_right));
      const V d2 = hn::Max(repair_abs_diff32(c, g.top_center), repair_abs_diff32(c, g.bottom_center));
      const V d3 = hn::Max(repair_abs_diff32(c, g.top_right), repair_abs_diff32(c, g.bottom_left));
      const V d4 = hn::Max(repair_abs_diff32(c, g.center_left), repair_abs_diff32(c, g.center_right));
      const V mindiff = hn::Min(hn::Min(d1, d2), hn::Min(d3, d4));
      V pmin = hn::Min(g.top_left, g.bottom_right);
      V pmax = hn::Max(g.top_left, g.bottom_right);
      pmin = hn::IfThenElse(hn::Eq(mindiff, d3), hn::Min(g.top_right, g.bottom_left), pmin);
      pmax = hn::IfThenElse(hn::Eq(mindiff, d3), hn::Max(g.top_right, g.bottom_left), pmax);
      pmin = hn::IfThenElse(hn::Eq(mindiff, d2), hn::Min(g.top_center, g.bottom_center), pmin);
      pmax = hn::IfThenElse(hn::Eq(mindiff, d2), hn::Max(g.top_center, g.bottom_center), pmax);
      pmin = hn::IfThenElse(hn::Eq(mindiff, d4), hn::Min(g.center_left, g.center_right), pmin);
      pmax = hn::IfThenElse(hn::Eq(mindiff, d4), hn::Max(g.center_left, g.center_right), pmax);
      const V mn = hn::Min(pmin, c);
      const V mx = hn::Max(pmax, c);
      return hn::Clamp(src, mn, mx);
    }
    case 19: {
      const V d1 = repair_abs_diff32(c, g.top_left);
      const V d2 = repair_abs_diff32(c, g.top_center);
      const V d3 = repair_abs_diff32(c, g.top_right);
      const V d4 = repair_abs_diff32(c, g.center_left);
      const V d5 = repair_abs_diff32(c, g.center_right);
      const V d6 = repair_abs_diff32(c, g.bottom_left);
      const V d7 = repair_abs_diff32(c, g.bottom_center);
      const V d8 = repair_abs_diff32(c, g.bottom_right);
      const V mindiff = hn::Min(hn::Min(hn::Min(d1, d2), hn::Min(d3, d4)),
                                hn::Min(hn::Min(d5, d6), hn::Min(d7, d8)));
      const V lo = hn::Clamp(hn::Sub(c, mindiff), zero, vmax);
      const V hi = hn::Clamp(hn::Add(c, mindiff), zero, vmax);
      return hn::Clamp(src, lo, hi);
    }
    case 20: {
      const V d1 = repair_abs_diff32(c, g.top_left);
      const V d2 = repair_abs_diff32(c, g.top_center);
      const V d3 = repair_abs_diff32(c, g.top_right);
      const V d4 = repair_abs_diff32(c, g.center_left);
      const V d5 = repair_abs_diff32(c, g.center_right);
      const V d6 = repair_abs_diff32(c, g.bottom_left);
      const V d7 = repair_abs_diff32(c, g.bottom_center);
      const V d8 = repair_abs_diff32(c, g.bottom_right);

      V maxdiff = hn::Max(d1, d2);
      V mindiff = hn::Min(d1, d2);

      #define NEO_SMO_STEP_DIFF(di) \
        maxdiff = hn::Clamp(maxdiff, hn::Min(mindiff, di), hn::Max(mindiff, di)); \
        mindiff = hn::Min(mindiff, di)

      NEO_SMO_STEP_DIFF(d3);
      NEO_SMO_STEP_DIFF(d4);
      NEO_SMO_STEP_DIFF(d5);
      NEO_SMO_STEP_DIFF(d6);
      NEO_SMO_STEP_DIFF(d7);
      NEO_SMO_STEP_DIFF(d8);
      #undef NEO_SMO_STEP_DIFF

      const V lo = hn::Clamp(hn::Sub(c, maxdiff), zero, vmax);
      const V hi = hn::Clamp(hn::Add(c, maxdiff), zero, vmax);
      return hn::Clamp(src, lo, hi);
    }
    case 21: {
      const auto s = g.min_max_opposites_without_center(d);
      const V d1 = hn::Clamp(hn::Sub(s.max1, c), zero, vmax);
      const V d2 = hn::Clamp(hn::Sub(s.max2, c), zero, vmax);
      const V d3 = hn::Clamp(hn::Sub(s.max3, c), zero, vmax);
      const V d4 = hn::Clamp(hn::Sub(s.max4, c), zero, vmax);
      const V rd1 = hn::Clamp(hn::Sub(c, s.min1), zero, vmax);
      const V rd2 = hn::Clamp(hn::Sub(c, s.min2), zero, vmax);
      const V rd3 = hn::Clamp(hn::Sub(c, s.min3), zero, vmax);
      const V rd4 = hn::Clamp(hn::Sub(c, s.min4), zero, vmax);
      const V u = hn::Min(hn::Min(hn::Max(d1, rd1), hn::Max(d2, rd2)),
                          hn::Min(hn::Max(d3, rd3), hn::Max(d4, rd4)));
      const V lo = hn::Clamp(hn::Sub(c, u), zero, vmax);
      const V hi = hn::Clamp(hn::Add(c, u), zero, vmax);
      return hn::Clamp(src, lo, hi);
    }
    case 22: {
      const V d1 = repair_abs_diff32(src, g.top_left);
      const V d2 = repair_abs_diff32(src, g.top_center);
      const V d3 = repair_abs_diff32(src, g.top_right);
      const V d4 = repair_abs_diff32(src, g.center_left);
      const V d5 = repair_abs_diff32(src, g.center_right);
      const V d6 = repair_abs_diff32(src, g.bottom_left);
      const V d7 = repair_abs_diff32(src, g.bottom_center);
      const V d8 = repair_abs_diff32(src, g.bottom_right);
      const V mindiff = hn::Min(hn::Min(hn::Min(d1, d2), hn::Min(d3, d4)),
                                hn::Min(hn::Min(d5, d6), hn::Min(d7, d8)));
      const V lo = hn::Clamp(hn::Sub(src, mindiff), zero, vmax);
      const V hi = hn::Clamp(hn::Add(src, mindiff), zero, vmax);
      return hn::Clamp(c, lo, hi);
    }
    case 23: {
      const V d1 = repair_abs_diff32(src, g.top_left);
      const V d2 = repair_abs_diff32(src, g.top_center);
      const V d3 = repair_abs_diff32(src, g.top_right);
      const V d4 = repair_abs_diff32(src, g.center_left);
      const V d5 = repair_abs_diff32(src, g.center_right);
      const V d6 = repair_abs_diff32(src, g.bottom_left);
      const V d7 = repair_abs_diff32(src, g.bottom_center);
      const V d8 = repair_abs_diff32(src, g.bottom_right);

      V maxdiff = hn::Max(d1, d2);
      V mindiff = hn::Min(d1, d2);

      #define NEO_SMO_STEP_DIFF(di) \
        maxdiff = hn::Clamp(maxdiff, hn::Min(mindiff, di), hn::Max(mindiff, di)); \
        mindiff = hn::Min(mindiff, di)

      NEO_SMO_STEP_DIFF(d3);
      NEO_SMO_STEP_DIFF(d4);
      NEO_SMO_STEP_DIFF(d5);
      NEO_SMO_STEP_DIFF(d6);
      NEO_SMO_STEP_DIFF(d7);
      NEO_SMO_STEP_DIFF(d8);
      #undef NEO_SMO_STEP_DIFF

      const V lo = hn::Clamp(hn::Sub(src, maxdiff), zero, vmax);
      const V hi = hn::Clamp(hn::Add(src, maxdiff), zero, vmax);
      return hn::Clamp(c, lo, hi);
    }
    case 24: {
      const auto s = g.min_max_opposites_without_center(d);
      const V d1 = hn::Clamp(hn::Sub(s.max1, src), zero, vmax);
      const V d2 = hn::Clamp(hn::Sub(s.max2, src), zero, vmax);
      const V d3 = hn::Clamp(hn::Sub(s.max3, src), zero, vmax);
      const V d4 = hn::Clamp(hn::Sub(s.max4, src), zero, vmax);
      const V rd1 = hn::Clamp(hn::Sub(src, s.min1), zero, vmax);
      const V rd2 = hn::Clamp(hn::Sub(src, s.min2), zero, vmax);
      const V rd3 = hn::Clamp(hn::Sub(src, s.min3), zero, vmax);
      const V rd4 = hn::Clamp(hn::Sub(src, s.min4), zero, vmax);
      const V u = hn::Min(hn::Min(hn::Max(d1, rd1), hn::Max(d2, rd2)),
                          hn::Min(hn::Max(d3, rd3), hn::Max(d4, rd4)));
      const V lo = hn::Clamp(hn::Sub(src, u), zero, vmax);
      const V hi = hn::Clamp(hn::Add(src, u), zero, vmax);
      return hn::Clamp(c, lo, hi);
    }
    default:
      return src;
  }
}

template <typename T>
void repair_int_impl(int mode, const T* srcp, const T* repairp, T* dstp, int width, int height,
                     std::size_t src_stride, std::size_t repair_stride, std::size_t dst_stride) {
  // These modes only select ranks and clamp; they need no widened arithmetic.
#define NEO_SMO_REPAIR_RANK(M) case M: rank_clamp_plane<true, M>(srcp, repairp, dstp, width, height, src_stride, repair_stride, dst_stride); return
  switch (mode) {
    case 11:
    NEO_SMO_REPAIR_RANK(1);
    NEO_SMO_REPAIR_RANK(2);
    NEO_SMO_REPAIR_RANK(3);
    NEO_SMO_REPAIR_RANK(4);
    NEO_SMO_REPAIR_RANK(5);
    NEO_SMO_REPAIR_RANK(6);
    NEO_SMO_REPAIR_RANK(8);
    NEO_SMO_REPAIR_RANK(15);
    NEO_SMO_REPAIR_RANK(16);
    NEO_SMO_REPAIR_RANK(18);
    NEO_SMO_REPAIR_RANK(9);
    NEO_SMO_REPAIR_RANK(17);
    NEO_SMO_REPAIR_RANK(19);
    NEO_SMO_REPAIR_RANK(20);
    NEO_SMO_REPAIR_RANK(21);
    NEO_SMO_REPAIR_RANK(22);
    NEO_SMO_REPAIR_RANK(23);
    NEO_SMO_REPAIR_RANK(24);
    NEO_SMO_REPAIR_RANK(12);
    NEO_SMO_REPAIR_RANK(13);
    NEO_SMO_REPAIR_RANK(14);
  }
#undef NEO_SMO_REPAIR_RANK
  constexpr int kRadius = 1;
  using ComputeT = std::conditional_t<sizeof(T) == 1, std::int16_t, std::int32_t>;
  hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + kSimdPad;
  std::vector<ComputeT> row_buffers(checked_product(3, padded_len));
  std::array<ComputeT*, 3> rows{
      row_buffers.data() + 0 * padded_len + kRadius,
      row_buffers.data() + 1 * padded_len + kRadius,
      row_buffers.data() + 2 * padded_len + kRadius,
  };
  const hn::Rebind<T, decltype(d)> ds;
  std::array<int, 3> cached_y{-1, -1, -1};

  const std::int32_t type_max = static_cast<std::int32_t>(std::numeric_limits<T>::max());
  auto fill_i32_row = [&](ComputeT* dst, const T* srow) {
    dst[-1] = srow[mirror_index(-1, width)];
    for (int x = 0; x < width; ++x) dst[x] = srow[x];
    dst[width] = srow[mirror_index(width, width)];
  };

  for (int y = 0; y < height; ++y) {
    const T* src_row = srcp + static_cast<std::size_t>(y) * src_stride;
    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;

    for (int dy = -1; dy <= 1; ++dy) {
      const std::size_t my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      const std::size_t slot = my % 3;
      auto* row = row_buffers.data() + slot * padded_len + kRadius;
      if (cached_y[slot] != static_cast<int>(my)) {
        fill_i32_row(row, repairp + my * repair_stride);
        cached_y[slot] = static_cast<int>(my);
      }
      rows[static_cast<std::size_t>(dy + 1)] = row;
    }

    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto g = Grid3x3<decltype(d)>::load(d, rows[0], rows[1], rows[2], x);
      const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
      const auto s = hn::PromoteTo(d, hn::LoadN(ds, src_row + x, count));
      const auto res = eval_repair_int(d, mode, s, g, type_max);
      hn::StoreN(hn::DemoteTo(ds, res), ds, dst_row + x, count);
    }

  }
}

void dispatch_repair_int_target(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp,
    const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width,
    std::size_t height, std::size_t src_stride_bytes,
    std::size_t repair_stride_bytes, std::size_t dst_stride_bytes) {
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);
  (void)chroma;
  if (dtype == DataType::U8) {
    repair_int_impl<std::uint8_t>(mode, srcp, repairp, dstp, w, h, src_stride_bytes,
                                  repair_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    repair_int_impl<std::uint16_t>(mode, reinterpret_cast<const std::uint16_t*>(srcp),
                                   reinterpret_cast<const std::uint16_t*>(repairp),
                                   reinterpret_cast<std::uint16_t*>(dstp), w, h,
                                   src_stride_bytes / 2, repair_stride_bytes / 2,
                                   dst_stride_bytes / 2);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_repair_int_target);
void process_repair_int_plane(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp,
    const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width,
    std::size_t height, std::size_t src_stride_bytes,
    std::size_t repair_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_repair_int_target)(dtype, mode, chroma, srcp, repairp, dstp, width, height, src_stride_bytes, repair_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
