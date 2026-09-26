#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/remove_grain.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/grid.hpp"

#include "common/float_arithmetic.hpp"
#include "common/fp16_rows.hpp"

inline bool rg_should_skip_line(int mode, int line) noexcept {
  if (mode == 13 || mode == 15) {
    return (line & 1) != 0;
  }
  if (mode == 14 || mode == 16) {
    return (line & 1) == 0;
  }
  return false;
}

// Evaluate RemoveGrain mode on 32-bit signed integer lanes (covers both u8 and u16 without overflow)
template <class D, class V = hn::Vec<D>>
HWY_INLINE V eval_rg_int32(D d, int mode, const Grid3x3<D>& g, std::int32_t type_max) {
  const V c = g.center_center;
  const V zero = hn::Zero(d);
  const V one = hn::Set(d, 1);
  const V four = hn::Set(d, 4);
  const V eight = hn::Set(d, 8);
  const V vmax = hn::Set(d, type_max);

  switch (mode) {
    case 1: {
      const V mn = g.min_without_center(d);
      const V mx = g.max_without_center(d);
      return hn::Clamp(c, mn, mx);
    }
    case 2:
    case 3:
    case 4: {
      V a[8];
      g.sort_without_center(d, a);
      const int lo = mode - 1;
      const int hi = 8 - mode;
      return hn::Clamp(c, a[lo], a[hi]);
    }
    case 5:
    case 6:
    case 7:
    case 8:
    case 9: {
      const auto s = g.min_max_opposites_without_center(d);
      const V clamp1 = hn::Clamp(c, s.min1, s.max1);
      const V clamp2 = hn::Clamp(c, s.min2, s.max2);
      const V clamp3 = hn::Clamp(c, s.min3, s.max3);
      const V clamp4 = hn::Clamp(c, s.min4, s.max4);

      V c1, c2, c3, c4;
      if (mode == 5) {
        c1 = hn::AbsDiff(c, clamp1);
        c2 = hn::AbsDiff(c, clamp2);
        c3 = hn::AbsDiff(c, clamp3);
        c4 = hn::AbsDiff(c, clamp4);
      } else {
        const V d1 = hn::Sub(s.max1, s.min1);
        const V d2 = hn::Sub(s.max2, s.min2);
        const V d3 = hn::Sub(s.max3, s.min3);
        const V d4 = hn::Sub(s.max4, s.min4);
        if (mode == 6) {
          c1 = hn::Min(hn::Add(hn::ShiftLeft<1>(hn::AbsDiff(c, clamp1)), d1), vmax);
          c2 = hn::Min(hn::Add(hn::ShiftLeft<1>(hn::AbsDiff(c, clamp2)), d2), vmax);
          c3 = hn::Min(hn::Add(hn::ShiftLeft<1>(hn::AbsDiff(c, clamp3)), d3), vmax);
          c4 = hn::Min(hn::Add(hn::ShiftLeft<1>(hn::AbsDiff(c, clamp4)), d4), vmax);
        } else if (mode == 7) {
          c1 = hn::Add(hn::AbsDiff(c, clamp1), d1);
          c2 = hn::Add(hn::AbsDiff(c, clamp2), d2);
          c3 = hn::Add(hn::AbsDiff(c, clamp3), d3);
          c4 = hn::Add(hn::AbsDiff(c, clamp4), d4);
        } else if (mode == 8) {
          c1 = hn::Clamp(hn::Add(hn::AbsDiff(c, clamp1), hn::ShiftLeft<1>(d1)), zero, vmax);
          c2 = hn::Clamp(hn::Add(hn::AbsDiff(c, clamp2), hn::ShiftLeft<1>(d2)), zero, vmax);
          c3 = hn::Clamp(hn::Add(hn::AbsDiff(c, clamp3), hn::ShiftLeft<1>(d3)), zero, vmax);
          c4 = hn::Clamp(hn::Add(hn::AbsDiff(c, clamp4), hn::ShiftLeft<1>(d4)), zero, vmax);
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
      const V d1 = hn::AbsDiff(c, g.top_left);
      const V d2 = hn::AbsDiff(c, g.top_center);
      const V d3 = hn::AbsDiff(c, g.top_right);
      const V d4 = hn::AbsDiff(c, g.center_left);
      const V d5 = hn::AbsDiff(c, g.center_right);
      const V d6 = hn::AbsDiff(c, g.bottom_left);
      const V d7 = hn::AbsDiff(c, g.bottom_center);
      const V d8 = hn::AbsDiff(c, g.bottom_right);
      const V mindiff = hn::Min(hn::Min(hn::Min(d1, d2), hn::Min(d3, d4)), hn::Min(hn::Min(d5, d6), hn::Min(d7, d8)));
      V res = g.center_left;
      res = hn::IfThenElse(hn::Eq(mindiff, d5), g.center_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d1), g.top_left, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d3), g.top_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), g.top_center, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d6), g.bottom_left, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d8), g.bottom_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d7), g.bottom_center, res);
      return res;
    }
    case 11:
    case 12: {
      const V cross = hn::Add(hn::Add(g.top_center, g.center_left), hn::Add(g.center_right, g.bottom_center));
      const V diag = hn::Add(hn::Add(g.top_left, g.top_right), hn::Add(g.bottom_left, g.bottom_right));
      const V sum = hn::Add(hn::Add(hn::ShiftLeft<2>(c), hn::ShiftLeft<1>(cross)), diag);
      return hn::ShiftRight<4>(hn::Add(sum, eight));
    }
    case 13:
    case 14: {
      const V d1 = hn::AbsDiff(g.top_left, g.bottom_right);
      const V d2 = hn::AbsDiff(g.top_center, g.bottom_center);
      const V d3 = hn::AbsDiff(g.top_right, g.bottom_left);
      const V mindiff = hn::Min(hn::Min(d1, d2), d3);
      const V a1 = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_left, g.bottom_right), one));
      const V a2 = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_center, g.bottom_center), one));
      const V a3 = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_right, g.bottom_left), one));
      V res = a1;
      res = hn::IfThenElse(hn::Eq(mindiff, d3), a3, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), a2, res);
      return res;
    }
    case 15:
    case 16: {
      const V d1 = hn::AbsDiff(g.top_left, g.bottom_right);
      const V d2 = hn::AbsDiff(g.top_center, g.bottom_center);
      const V d3 = hn::AbsDiff(g.top_right, g.bottom_left);
      const V mindiff = hn::Min(hn::Min(d1, d2), d3);
      const V sum = hn::Add(hn::ShiftLeft<1>(hn::Add(g.top_center, g.bottom_center)),
                            hn::Add(hn::Add(g.top_left, g.top_right), hn::Add(g.bottom_left, g.bottom_right)));
      const V avg = hn::ShiftRight<3>(hn::Add(sum, four));
      const V c1 = hn::Clamp(avg, hn::Min(g.top_left, g.bottom_right), hn::Max(g.top_left, g.bottom_right));
      const V c2 = hn::Clamp(avg, hn::Min(g.top_center, g.bottom_center), hn::Max(g.top_center, g.bottom_center));
      const V c3 = hn::Clamp(avg, hn::Min(g.top_right, g.bottom_left), hn::Max(g.top_right, g.bottom_left));
      V res = c1;
      res = hn::IfThenElse(hn::Eq(mindiff, d3), c3, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), c2, res);
      return res;
    }
    case 17: {
      const auto s = g.min_max_opposites_without_center(d);
      const V l = hn::Max(hn::Max(s.min1, s.min2), hn::Max(s.min3, s.min4));
      const V u = hn::Min(hn::Min(s.max1, s.max2), hn::Min(s.max3, s.max4));
      return hn::Clamp(c, hn::Min(l, u), hn::Max(l, u));
    }
    case 18: {
      const V d1 = hn::Max(hn::AbsDiff(c, g.top_left), hn::AbsDiff(c, g.bottom_right));
      const V d2 = hn::Max(hn::AbsDiff(c, g.top_center), hn::AbsDiff(c, g.bottom_center));
      const V d3 = hn::Max(hn::AbsDiff(c, g.top_right), hn::AbsDiff(c, g.bottom_left));
      const V d4 = hn::Max(hn::AbsDiff(c, g.center_left), hn::AbsDiff(c, g.center_right));
      const V mindiff = hn::Min(hn::Min(d1, d2), hn::Min(d3, d4));
      const V c1 = hn::Clamp(c, hn::Min(g.top_left, g.bottom_right), hn::Max(g.top_left, g.bottom_right));
      const V c2 = hn::Clamp(c, hn::Min(g.top_center, g.bottom_center), hn::Max(g.top_center, g.bottom_center));
      const V c3 = hn::Clamp(c, hn::Min(g.top_right, g.bottom_left), hn::Max(g.top_right, g.bottom_left));
      const V c4 = hn::Clamp(c, hn::Min(g.center_left, g.center_right), hn::Max(g.center_left, g.center_right));
      V res = c1;
      res = hn::IfThenElse(hn::Eq(mindiff, d3), c3, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), c2, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d4), c4, res);
      return res;
    }
    case 19: {
      const V sum = hn::Add(hn::Add(hn::Add(g.top_left, g.top_center), hn::Add(g.top_right, g.center_left)),
                            hn::Add(hn::Add(g.center_right, g.bottom_left), hn::Add(g.bottom_center, g.bottom_right)));
      return hn::ShiftRight<3>(hn::Add(sum, four));
    }
    case 20: {
      const V sum8 = hn::Add(hn::Add(hn::Add(g.top_left, g.top_center), hn::Add(g.top_right, g.center_left)),
                             hn::Add(hn::Add(g.center_right, g.bottom_left), hn::Add(g.bottom_center, g.bottom_right)));
      const V sum = hn::Add(hn::Add(sum8, c), four);
      // sum is in [4, 9 * 65535 + 4]. Multiplication by ceil(2^33 / 9)
      // followed by a 33-bit shift gives the exact nonnegative quotient.
      return hn::ShiftRight<1>(hn::MulHigh(sum, hn::Set(d, 0x38E38E39)));
    }
    case 21: {
      const V l1l = hn::ShiftRight<1>(hn::Add(g.top_left, g.bottom_right));
      const V l2l = hn::ShiftRight<1>(hn::Add(g.top_center, g.bottom_center));
      const V l3l = hn::ShiftRight<1>(hn::Add(g.top_right, g.bottom_left));
      const V l4l = hn::ShiftRight<1>(hn::Add(g.center_left, g.center_right));
      const V l1h = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_left, g.bottom_right), one));
      const V l2h = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_center, g.bottom_center), one));
      const V l3h = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_right, g.bottom_left), one));
      const V l4h = hn::ShiftRight<1>(hn::Add(hn::Add(g.center_left, g.center_right), one));
      const V mn = hn::Min(hn::Min(l1l, l2l), hn::Min(l3l, l4l));
      const V mx = hn::Max(hn::Max(l1h, l2h), hn::Max(l3h, l4h));
      return hn::Clamp(c, mn, mx);
    }
    case 22: {
      const V l1 = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_left, g.bottom_right), one));
      const V l2 = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_center, g.bottom_center), one));
      const V l3 = hn::ShiftRight<1>(hn::Add(hn::Add(g.top_right, g.bottom_left), one));
      const V l4 = hn::ShiftRight<1>(hn::Add(hn::Add(g.center_left, g.center_right), one));
      const V mn = hn::Min(hn::Min(l1, l2), hn::Min(l3, l4));
      const V mx = hn::Max(hn::Max(l1, l2), hn::Max(l3, l4));
      return hn::Clamp(c, mn, mx);
    }
    case 23: {
      const auto s = g.min_max_opposites_without_center(d);
      const V ld1 = hn::Sub(s.max1, s.min1);
      const V ld2 = hn::Sub(s.max2, s.min2);
      const V ld3 = hn::Sub(s.max3, s.min3);
      const V ld4 = hn::Sub(s.max4, s.min4);
      const V h1 = hn::Min(hn::Sub(c, s.max1), ld1);
      const V h2 = hn::Min(hn::Sub(c, s.max2), ld2);
      const V h3 = hn::Min(hn::Sub(c, s.max3), ld3);
      const V h4 = hn::Min(hn::Sub(c, s.max4), ld4);
      const V h = hn::Max(zero, hn::Max(hn::Max(h1, h2), hn::Max(h3, h4)));
      const V l1 = hn::Min(hn::Sub(s.min1, c), ld1);
      const V l2 = hn::Min(hn::Sub(s.min2, c), ld2);
      const V l3 = hn::Min(hn::Sub(s.min3, c), ld3);
      const V l4 = hn::Min(hn::Sub(s.min4, c), ld4);
      const V l = hn::Max(zero, hn::Max(hn::Max(l1, l2), hn::Max(l3, l4)));
      return hn::Add(hn::Sub(c, h), l);
    }
    case 24: {
      const auto s = g.min_max_opposites_without_center(d);
      const V ld1 = hn::Sub(s.max1, s.min1);
      const V ld2 = hn::Sub(s.max2, s.min2);
      const V ld3 = hn::Sub(s.max3, s.min3);
      const V ld4 = hn::Sub(s.max4, s.min4);
      const V th1 = hn::Sub(c, s.max1);
      const V th2 = hn::Sub(c, s.max2);
      const V th3 = hn::Sub(c, s.max3);
      const V th4 = hn::Sub(c, s.max4);
      const V h1 = hn::Min(th1, hn::Sub(ld1, th1));
      const V h2 = hn::Min(th2, hn::Sub(ld2, th2));
      const V h3 = hn::Min(th3, hn::Sub(ld3, th3));
      const V h4 = hn::Min(th4, hn::Sub(ld4, th4));
      const V h = hn::Max(zero, hn::Max(hn::Max(h1, h2), hn::Max(h3, h4)));
      const V tl1 = hn::Sub(s.min1, c);
      const V tl2 = hn::Sub(s.min2, c);
      const V tl3 = hn::Sub(s.min3, c);
      const V tl4 = hn::Sub(s.min4, c);
      const V l1 = hn::Min(tl1, hn::Sub(ld1, tl1));
      const V l2 = hn::Min(tl2, hn::Sub(ld2, tl2));
      const V l3 = hn::Min(tl3, hn::Sub(ld3, tl3));
      const V l4 = hn::Min(tl4, hn::Sub(ld4, tl4));
      const V l = hn::Max(zero, hn::Max(hn::Max(l1, l2), hn::Max(l3, l4)));
      return hn::Add(hn::Sub(c, h), l);
    }
    default:
      return c;
  }
}

template <bool IsF16, class D, class V = hn::Vec<D>>
HWY_INLINE V eval_rg_float(D d, int mode, const Grid3x3<D>& g, bool chroma) {
  const V c = g.center_center;
  const V zero = hn::Zero(d);
  const V two = hn::Set(d, 2.0f);
  const V four = hn::Set(d, 4.0f);
  const V inv2 = hn::Set(d, 0.5f);
  const V inv8 = hn::Set(d, 0.125f);
  const V inv16 = hn::Set(d, 0.0625f);
  const V vmin = hn::Set(d, chroma ? -0.5f : 0.0f);
  const V vmax = hn::Set(d, chroma ? 0.5f : 1.0f);

  switch (mode) {
    case 1: {
      const V mn = g.min_without_center(d);
      const V mx = g.max_without_center(d);
      return hn::Clamp(c, mn, mx);
    }
    case 2:
    case 3:
    case 4: {
      V a[8];
      g.sort_without_center(d, a);
      const int lo = mode - 1;
      const int hi = 8 - mode;
      return hn::Clamp(c, a[lo], a[hi]);
    }
    case 5:
    case 6:
    case 7:
    case 8:
    case 9: {
      const auto s = g.min_max_opposites_without_center(d);
      const V clamp1 = hn::Clamp(c, s.min1, s.max1);
      const V clamp2 = hn::Clamp(c, s.min2, s.max2);
      const V clamp3 = hn::Clamp(c, s.min3, s.max3);
      const V clamp4 = hn::Clamp(c, s.min4, s.max4);

      V c1, c2, c3, c4;
      if (mode == 5) {
        c1 = float_abs_diff<IsF16>(d, c, clamp1);
        c2 = float_abs_diff<IsF16>(d, c, clamp2);
        c3 = float_abs_diff<IsF16>(d, c, clamp3);
        c4 = float_abs_diff<IsF16>(d, c, clamp4);
      } else {
        const V d1 = float_sub<IsF16>(d, s.max1, s.min1);
        const V d2 = float_sub<IsF16>(d, s.max2, s.min2);
        const V d3 = float_sub<IsF16>(d, s.max3, s.min3);
        const V d4 = float_sub<IsF16>(d, s.max4, s.min4);
        if (mode == 6) {
          c1 = hn::Min(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp1), two), d1), vmax);
          c2 = hn::Min(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp2), two), d2), vmax);
          c3 = hn::Min(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp3), two), d3), vmax);
          c4 = hn::Min(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp4), two), d4), vmax);
        } else if (mode == 7) {
          c1 = float_add<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp1), d1);
          c2 = float_add<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp2), d2);
          c3 = float_add<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp3), d3);
          c4 = float_add<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp4), d4);
        } else if (mode == 8) {
          c1 = hn::Clamp(float_add<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp1), float_mul<IsF16>(d, d1, two)), vmin, vmax);
          c2 = hn::Clamp(float_add<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp2), float_mul<IsF16>(d, d2, two)), vmin, vmax);
          c3 = hn::Clamp(float_add<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp3), float_mul<IsF16>(d, d3, two)), vmin, vmax);
          c4 = hn::Clamp(float_add<IsF16>(d, float_abs_diff<IsF16>(d, c, clamp4), float_mul<IsF16>(d, d4, two)), vmin, vmax);
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
      const V d1 = float_abs_diff<IsF16>(d, c, g.top_left);
      const V d2 = float_abs_diff<IsF16>(d, c, g.top_center);
      const V d3 = float_abs_diff<IsF16>(d, c, g.top_right);
      const V d4 = float_abs_diff<IsF16>(d, c, g.center_left);
      const V d5 = float_abs_diff<IsF16>(d, c, g.center_right);
      const V d6 = float_abs_diff<IsF16>(d, c, g.bottom_left);
      const V d7 = float_abs_diff<IsF16>(d, c, g.bottom_center);
      const V d8 = float_abs_diff<IsF16>(d, c, g.bottom_right);
      const V mindiff = hn::Min(hn::Min(hn::Min(d1, d2), hn::Min(d3, d4)), hn::Min(hn::Min(d5, d6), hn::Min(d7, d8)));
      V res = g.center_left;
      res = hn::IfThenElse(hn::Eq(mindiff, d5), g.center_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d1), g.top_left, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d3), g.top_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), g.top_center, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d6), g.bottom_left, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d8), g.bottom_right, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d7), g.bottom_center, res);
      return res;
    }
    case 11:
    case 12: {
      // Exact left-to-right evaluation matching Zig:
      // 4 * center_center + 2 * (top_center + center_left + center_right + bottom_center) + top_left + top_right + bottom_left + bottom_right
      const V cross = float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, g.top_center, g.center_left), g.center_right), g.bottom_center);
      const V term1 = float_add<IsF16>(d, float_mul<IsF16>(d, four, c), float_mul<IsF16>(d, two, cross));
      const V sum = float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, term1, g.top_left), g.top_right), g.bottom_left), g.bottom_right);
      return float_mul<IsF16>(d, sum, inv16);
    }
    case 13:
    case 14: {
      const V d1 = float_abs_diff<IsF16>(d, g.top_left, g.bottom_right);
      const V d2 = float_abs_diff<IsF16>(d, g.top_center, g.bottom_center);
      const V d3 = float_abs_diff<IsF16>(d, g.top_right, g.bottom_left);
      const V mindiff = hn::Min(hn::Min(d1, d2), d3);
      const V a1 = float_mul<IsF16>(d, float_add<IsF16>(d, g.top_left, g.bottom_right), inv2);
      const V a2 = float_mul<IsF16>(d, float_add<IsF16>(d, g.top_center, g.bottom_center), inv2);
      const V a3 = float_mul<IsF16>(d, float_add<IsF16>(d, g.top_right, g.bottom_left), inv2);
      V res = a1;
      res = hn::IfThenElse(hn::Eq(mindiff, d3), a3, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), a2, res);
      return res;
    }
    case 15:
    case 16: {
      const V d1 = float_abs_diff<IsF16>(d, g.top_left, g.bottom_right);
      const V d2 = float_abs_diff<IsF16>(d, g.top_center, g.bottom_center);
      const V d3 = float_abs_diff<IsF16>(d, g.top_right, g.bottom_left);
      const V mindiff = hn::Min(hn::Min(d1, d2), d3);
      // (2 * (top_center + bottom_center) + top_left + top_right + bottom_left + bottom_right) / 8
      const V term1 = float_mul<IsF16>(d, two, float_add<IsF16>(d, g.top_center, g.bottom_center));
      const V sum = float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, term1, g.top_left), g.top_right), g.bottom_left), g.bottom_right);
      const V avg = float_mul<IsF16>(d, sum, inv8);
      const V c1 = hn::Clamp(avg, hn::Min(g.top_left, g.bottom_right), hn::Max(g.top_left, g.bottom_right));
      const V c2 = hn::Clamp(avg, hn::Min(g.top_center, g.bottom_center), hn::Max(g.top_center, g.bottom_center));
      const V c3 = hn::Clamp(avg, hn::Min(g.top_right, g.bottom_left), hn::Max(g.top_right, g.bottom_left));
      V res = c1;
      res = hn::IfThenElse(hn::Eq(mindiff, d3), c3, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), c2, res);
      return res;
    }
    case 17: {
      const auto s = g.min_max_opposites_without_center(d);
      const V l = hn::Max(hn::Max(s.min1, s.min2), hn::Max(s.min3, s.min4));
      const V u = hn::Min(hn::Min(s.max1, s.max2), hn::Min(s.max3, s.max4));
      return hn::Clamp(c, hn::Min(l, u), hn::Max(l, u));
    }
    case 18: {
      const V d1 = hn::Max(float_abs_diff<IsF16>(d, c, g.top_left), float_abs_diff<IsF16>(d, c, g.bottom_right));
      const V d2 = hn::Max(float_abs_diff<IsF16>(d, c, g.top_center), float_abs_diff<IsF16>(d, c, g.bottom_center));
      const V d3 = hn::Max(float_abs_diff<IsF16>(d, c, g.top_right), float_abs_diff<IsF16>(d, c, g.bottom_left));
      const V d4 = hn::Max(float_abs_diff<IsF16>(d, c, g.center_left), float_abs_diff<IsF16>(d, c, g.center_right));
      const V mindiff = hn::Min(hn::Min(d1, d2), hn::Min(d3, d4));
      const V c1 = hn::Clamp(c, hn::Min(g.top_left, g.bottom_right), hn::Max(g.top_left, g.bottom_right));
      const V c2 = hn::Clamp(c, hn::Min(g.top_center, g.bottom_center), hn::Max(g.top_center, g.bottom_center));
      const V c3 = hn::Clamp(c, hn::Min(g.top_right, g.bottom_left), hn::Max(g.top_right, g.bottom_left));
      const V c4 = hn::Clamp(c, hn::Min(g.center_left, g.center_right), hn::Max(g.center_left, g.center_right));
      V res = c1;
      res = hn::IfThenElse(hn::Eq(mindiff, d3), c3, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d2), c2, res);
      res = hn::IfThenElse(hn::Eq(mindiff, d4), c4, res);
      return res;
    }
    case 19: {
      // top_left + top_center + top_right + center_left + center_right + bottom_left + bottom_center + bottom_right
      const V s1 = float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, g.top_left, g.top_center), g.top_right), g.center_left);
      const V sum = float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, s1, g.center_right), g.bottom_left), g.bottom_center), g.bottom_right);
      return float_mul<IsF16>(d, sum, inv8);
    }
    case 20: {
      // top_left + top_center + top_right + center_center + center_left + center_right + bottom_left + bottom_center + bottom_right
      const V s1 = float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, g.top_left, g.top_center), g.top_right), c), g.center_left);
      const V sum = float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, float_add<IsF16>(d, s1, g.center_right), g.bottom_left), g.bottom_center), g.bottom_right);
      return float_div9<IsF16>(d, sum);
    }
    case 21:
    case 22: {
      const V l1 = float_mul<IsF16>(d, float_add<IsF16>(d, g.top_left, g.bottom_right), inv2);
      const V l2 = float_mul<IsF16>(d, float_add<IsF16>(d, g.top_center, g.bottom_center), inv2);
      const V l3 = float_mul<IsF16>(d, float_add<IsF16>(d, g.top_right, g.bottom_left), inv2);
      const V l4 = float_mul<IsF16>(d, float_add<IsF16>(d, g.center_left, g.center_right), inv2);
      const V mn = hn::Min(hn::Min(l1, l2), hn::Min(l3, l4));
      const V mx = hn::Max(hn::Max(l1, l2), hn::Max(l3, l4));
      return hn::Clamp(c, mn, mx);
    }
    case 23: {
      const auto s = g.min_max_opposites_without_center(d);
      const V ld1 = float_sub<IsF16>(d, s.max1, s.min1);
      const V ld2 = float_sub<IsF16>(d, s.max2, s.min2);
      const V ld3 = float_sub<IsF16>(d, s.max3, s.min3);
      const V ld4 = float_sub<IsF16>(d, s.max4, s.min4);
      const V h1 = hn::Min(float_sub<IsF16>(d, c, s.max1), ld1);
      const V h2 = hn::Min(float_sub<IsF16>(d, c, s.max2), ld2);
      const V h3 = hn::Min(float_sub<IsF16>(d, c, s.max3), ld3);
      const V h4 = hn::Min(float_sub<IsF16>(d, c, s.max4), ld4);
      const V h = hn::Max(zero, hn::Max(hn::Max(h1, h2), hn::Max(h3, h4)));
      const V l1 = hn::Min(float_sub<IsF16>(d, s.min1, c), ld1);
      const V l2 = hn::Min(float_sub<IsF16>(d, s.min2, c), ld2);
      const V l3 = hn::Min(float_sub<IsF16>(d, s.min3, c), ld3);
      const V l4 = hn::Min(float_sub<IsF16>(d, s.min4, c), ld4);
      const V l = hn::Max(zero, hn::Max(hn::Max(l1, l2), hn::Max(l3, l4)));
      return float_add<IsF16>(d, float_sub<IsF16>(d, c, h), l);
    }
    case 24: {
      const auto s = g.min_max_opposites_without_center(d);
      const V ld1 = float_sub<IsF16>(d, s.max1, s.min1);
      const V ld2 = float_sub<IsF16>(d, s.max2, s.min2);
      const V ld3 = float_sub<IsF16>(d, s.max3, s.min3);
      const V ld4 = float_sub<IsF16>(d, s.max4, s.min4);
      const V th1 = float_sub<IsF16>(d, c, s.max1);
      const V th2 = float_sub<IsF16>(d, c, s.max2);
      const V th3 = float_sub<IsF16>(d, c, s.max3);
      const V th4 = float_sub<IsF16>(d, c, s.max4);
      const V h1 = hn::Min(th1, float_sub<IsF16>(d, ld1, th1));
      const V h2 = hn::Min(th2, float_sub<IsF16>(d, ld2, th2));
      const V h3 = hn::Min(th3, float_sub<IsF16>(d, ld3, th3));
      const V h4 = hn::Min(th4, float_sub<IsF16>(d, ld4, th4));
      const V h = hn::Max(zero, hn::Max(hn::Max(h1, h2), hn::Max(h3, h4)));
      const V tl1 = float_sub<IsF16>(d, s.min1, c);
      const V tl2 = float_sub<IsF16>(d, s.min2, c);
      const V tl3 = float_sub<IsF16>(d, s.min3, c);
      const V tl4 = float_sub<IsF16>(d, s.min4, c);
      const V l1 = hn::Min(tl1, float_sub<IsF16>(d, ld1, tl1));
      const V l2 = hn::Min(tl2, float_sub<IsF16>(d, ld2, tl2));
      const V l3 = hn::Min(tl3, float_sub<IsF16>(d, ld3, tl3));
      const V l4 = hn::Min(tl4, float_sub<IsF16>(d, ld4, tl4));
      const V l = hn::Max(zero, hn::Max(hn::Max(l1, l2), hn::Max(l3, l4)));
      return float_add<IsF16>(d, float_sub<IsF16>(d, c, h), l);
    }
    default:
      return c;
  }
}

template <typename T>
void remove_grain_int_impl(int mode, const T* srcp, T* dstp, int width, int height, std::size_t src_stride,
                           std::size_t dst_stride) {
  constexpr int kRadius = 1;
  hn::ScalableTag<std::int32_t> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + kSimdPad;
  std::vector<std::int32_t> row_buffers(checked_product(3, padded_len));
  std::array<std::int32_t*, 3> rows{
      row_buffers.data() + 0 * padded_len + kRadius,
      row_buffers.data() + 1 * padded_len + kRadius,
      row_buffers.data() + 2 * padded_len + kRadius,
  };
  std::array<int, 3> cached_y{-1, -1, -1};
  const hn::Rebind<T, decltype(d)> ds;

  const std::int32_t type_max = static_cast<std::int32_t>(std::numeric_limits<T>::max());
  auto fill_i32_row = [&](std::int32_t* dst, const T* srow) {
    dst[-1] = static_cast<std::int32_t>(srow[mirror_index(-1, width)]);
    // Keep mirror indexing out of the interior so widening is a contiguous load.
    for (int x = 0; x < width; ++x) {
      dst[x] = static_cast<std::int32_t>(srow[x]);
    }
    dst[width] = static_cast<std::int32_t>(srow[mirror_index(width, width)]);
  };

  for (int y = 0; y < height; ++y) {
    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    if (y > 0 && y < height - 1 && rg_should_skip_line(mode, y)) {
      std::memcpy(dst_row, srcp + static_cast<std::size_t>(y) * src_stride,
                  static_cast<std::size_t>(width) * sizeof(T));
      continue;
    }

    for (int dy = -1; dy <= 1; ++dy) {
      const std::size_t my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      const std::size_t slot = my % 3;
      auto* row = row_buffers.data() + slot * padded_len + kRadius;
      if (cached_y[slot] != static_cast<int>(my)) {
        fill_i32_row(row, srcp + my * src_stride);
        cached_y[slot] = static_cast<int>(my);
      }
      rows[static_cast<std::size_t>(dy + 1)] = row;
    }

    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto g = Grid3x3<decltype(d)>::load(d, rows[0], rows[1], rows[2], static_cast<std::size_t>(x));
      const auto res = eval_rg_int32(d, mode, g, type_max);
      hn::StoreN(hn::DemoteTo(ds, res), ds, dst_row + x,
                 std::min(lanes, static_cast<std::size_t>(width) - x));
    }

  }
}

template <bool IsF16, typename StorageT>
void remove_grain_float_impl(int mode, bool chroma, const StorageT* srcp, StorageT* dstp, int width, int height,
                             std::size_t src_stride, std::size_t dst_stride) {
  using ComputeT = FloatLane<IsF16>;
  constexpr int kRadius = 1;
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
  std::vector<ComputeT> out_f32(static_cast<std::size_t>(width) + kSimdPad);

  for (int y = 0; y < height; ++y) {
    StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    if (y > 0 && y < height - 1 && rg_should_skip_line(mode, y)) {
      std::memcpy(dst_row, srcp + static_cast<std::size_t>(y) * src_stride,
                  static_cast<std::size_t>(width) * sizeof(StorageT));
      continue;
    }

    for (int dy = -1; dy <= 1; ++dy) {
      const std::size_t my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      if constexpr (IsF16) {
        fill_mirrored_row_f16(rows[static_cast<std::size_t>(dy + 1)] - kRadius, srcp + my * src_stride, width,
                                       kRadius);
      } else {
        fill_mirrored_row(rows[static_cast<std::size_t>(dy + 1)] - kRadius, srcp + my * src_stride, width, kRadius);
      }
    }

    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto g = Grid3x3<decltype(d)>::load(d, rows[0], rows[1], rows[2], static_cast<std::size_t>(x));
      const auto res = eval_rg_float<IsF16>(d, mode, g, chroma);
      hn::StoreU(res, d, out_f32.data() + x);
    }

    if constexpr (IsF16) {
      store_row_f16(dst_row, out_f32.data(), width);
    } else {
      std::memcpy(dst_row, out_f32.data(), static_cast<std::size_t>(width) * sizeof(float));
    }
  }
}

void dispatch_remove_grain_target(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp, std::uint8_t* dstp,
                                  std::size_t width, std::size_t height, std::size_t src_stride_bytes,
                                  std::size_t dst_stride_bytes) {
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);
  if (dtype == DataType::U8) {
    remove_grain_int_impl<std::uint8_t>(mode, srcp, dstp, w, h, src_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    remove_grain_int_impl<std::uint16_t>(mode, reinterpret_cast<const std::uint16_t*>(srcp),
                                         reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2,
                                         dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    remove_grain_float_impl<true, std::uint16_t>(mode, chroma, reinterpret_cast<const std::uint16_t*>(srcp),
                                                 reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2,
                                                 dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    remove_grain_float_impl<false, float>(mode, chroma, reinterpret_cast<const float*>(srcp),
                                          reinterpret_cast<float*>(dstp), w, h, src_stride_bytes / 4,
                                          dst_stride_bytes / 4);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_remove_grain_target);
void process_remove_grain_plane(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp, std::uint8_t* dstp,
                                std::size_t width, std::size_t height, std::size_t src_stride_bytes,
                                std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_remove_grain_target)(dtype, mode, chroma, srcp, dstp, width, height, src_stride_bytes,
                                                     dst_stride_bytes);
}
} // namespace neo_smo
#endif
