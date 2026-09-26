// Included inside the per-target Highway namespace. No include guard.
#include "common/grid.hpp"
#include "common/float_arithmetic.hpp"

// Evaluate Repair mode on floating-point lanes
template <bool IsF16, class D, class V = hn::Vec<D>>
HWY_NOINLINE V eval_repair_float(D d, int mode, V src, const Grid3x3<D>& g, V vmin, V vmax) {
  const V c = g.center_center;
  const V two = hn::Set(d, 2.0f);

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
        c1 = float_abs_diff<IsF16>(d, src, clamp1);
        c2 = float_abs_diff<IsF16>(d, src, clamp2);
        c3 = float_abs_diff<IsF16>(d, src, clamp3);
        c4 = float_abs_diff<IsF16>(d, src, clamp4);
      } else {
        const V d1 = float_sub<IsF16>(d, s.max1, s.min1);
        const V d2 = float_sub<IsF16>(d, s.max2, s.min2);
        const V d3 = float_sub<IsF16>(d, s.max3, s.min3);
        const V d4 = float_sub<IsF16>(d, s.max4, s.min4);
        if (mode == 6) {
          c1 = hn::Min(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp1), two), d1), vmax);
          c2 = hn::Min(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp2), two), d2), vmax);
          c3 = hn::Min(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp3), two), d3), vmax);
          c4 = hn::Min(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp4), two), d4), vmax);
        } else if (mode == 7) {
          c1 = float_add<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp1), d1);
          c2 = float_add<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp2), d2);
          c3 = float_add<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp3), d3);
          c4 = float_add<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp4), d4);
        } else if (mode == 8) {
          c1 = hn::Clamp(float_add<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp1), float_mul<IsF16>(d, d1, two)), vmin, vmax);
          c2 = hn::Clamp(float_add<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp2), float_mul<IsF16>(d, d2, two)), vmin, vmax);
          c3 = hn::Clamp(float_add<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp3), float_mul<IsF16>(d, d3, two)), vmin, vmax);
          c4 = hn::Clamp(float_add<IsF16>(d, float_abs_diff<IsF16>(d, src, clamp4), float_mul<IsF16>(d, d4, two)), vmin, vmax);
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
      const V d1 = float_abs_diff<IsF16>(d, src, g.top_left);
      const V d2 = float_abs_diff<IsF16>(d, src, g.top_center);
      const V d3 = float_abs_diff<IsF16>(d, src, g.top_right);
      const V d4 = float_abs_diff<IsF16>(d, src, g.center_left);
      const V d5 = float_abs_diff<IsF16>(d, src, g.center_right);
      const V d6 = float_abs_diff<IsF16>(d, src, g.bottom_left);
      const V d7 = float_abs_diff<IsF16>(d, src, g.bottom_center);
      const V d8 = float_abs_diff<IsF16>(d, src, g.bottom_right);
      const V dc = float_abs_diff<IsF16>(d, src, c);
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
        c1 = float_abs_diff<IsF16>(d, c, hn::Clamp(c, s.min1, s.max1));
        c2 = float_abs_diff<IsF16>(d, c, hn::Clamp(c, s.min2, s.max2));
        c3 = float_abs_diff<IsF16>(d, c, hn::Clamp(c, s.min3, s.max3));
        c4 = float_abs_diff<IsF16>(d, c, hn::Clamp(c, s.min4, s.max4));
      } else {
        const V d1 = float_sub<IsF16>(d, s.max1, s.min1);
        const V d2 = float_sub<IsF16>(d, s.max2, s.min2);
        const V d3 = float_sub<IsF16>(d, s.max3, s.min3);
        const V d4 = float_sub<IsF16>(d, s.max4, s.min4);
        c1 = hn::Clamp(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, c, hn::Clamp(c, s.min1, s.max1)), two), d1), vmin, vmax);
        c2 = hn::Clamp(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, c, hn::Clamp(c, s.min2, s.max2)), two), d2), vmin, vmax);
        c3 = hn::Clamp(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, c, hn::Clamp(c, s.min3, s.max3)), two), d3), vmin, vmax);
        c4 = hn::Clamp(float_add<IsF16>(d, float_mul<IsF16>(d, float_abs_diff<IsF16>(d, c, hn::Clamp(c, s.min4, s.max4)), two), d4), vmin, vmax);
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
      const V d1 = hn::Max(float_abs_diff<IsF16>(d, c, g.top_left), float_abs_diff<IsF16>(d, c, g.bottom_right));
      const V d2 = hn::Max(float_abs_diff<IsF16>(d, c, g.top_center), float_abs_diff<IsF16>(d, c, g.bottom_center));
      const V d3 = hn::Max(float_abs_diff<IsF16>(d, c, g.top_right), float_abs_diff<IsF16>(d, c, g.bottom_left));
      const V d4 = hn::Max(float_abs_diff<IsF16>(d, c, g.center_left), float_abs_diff<IsF16>(d, c, g.center_right));
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
      const V d1 = float_abs_diff<IsF16>(d, c, g.top_left);
      const V d2 = float_abs_diff<IsF16>(d, c, g.top_center);
      const V d3 = float_abs_diff<IsF16>(d, c, g.top_right);
      const V d4 = float_abs_diff<IsF16>(d, c, g.center_left);
      const V d5 = float_abs_diff<IsF16>(d, c, g.center_right);
      const V d6 = float_abs_diff<IsF16>(d, c, g.bottom_left);
      const V d7 = float_abs_diff<IsF16>(d, c, g.bottom_center);
      const V d8 = float_abs_diff<IsF16>(d, c, g.bottom_right);
      const V mindiff = hn::Min(hn::Min(hn::Min(d1, d2), hn::Min(d3, d4)),
                                hn::Min(hn::Min(d5, d6), hn::Min(d7, d8)));
      const V lo = hn::Clamp(float_sub<IsF16>(d, c, mindiff), vmin, vmax);
      const V hi = hn::Clamp(float_add<IsF16>(d, c, mindiff), vmin, vmax);
      return hn::Clamp(src, lo, hi);
    }
    case 20: {
      const V d1 = float_abs_diff<IsF16>(d, c, g.top_left);
      const V d2 = float_abs_diff<IsF16>(d, c, g.top_center);
      const V d3 = float_abs_diff<IsF16>(d, c, g.top_right);
      const V d4 = float_abs_diff<IsF16>(d, c, g.center_left);
      const V d5 = float_abs_diff<IsF16>(d, c, g.center_right);
      const V d6 = float_abs_diff<IsF16>(d, c, g.bottom_left);
      const V d7 = float_abs_diff<IsF16>(d, c, g.bottom_center);
      const V d8 = float_abs_diff<IsF16>(d, c, g.bottom_right);

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

      const V lo = hn::Clamp(float_sub<IsF16>(d, c, maxdiff), vmin, vmax);
      const V hi = hn::Clamp(float_add<IsF16>(d, c, maxdiff), vmin, vmax);
      return hn::Clamp(src, lo, hi);
    }
    case 21: {
      const auto s = g.min_max_opposites_without_center(d);
      const V d1 = hn::Clamp(float_sub<IsF16>(d, s.max1, c), vmin, vmax);
      const V d2 = hn::Clamp(float_sub<IsF16>(d, s.max2, c), vmin, vmax);
      const V d3 = hn::Clamp(float_sub<IsF16>(d, s.max3, c), vmin, vmax);
      const V d4 = hn::Clamp(float_sub<IsF16>(d, s.max4, c), vmin, vmax);
      const V rd1 = hn::Clamp(float_sub<IsF16>(d, c, s.min1), vmin, vmax);
      const V rd2 = hn::Clamp(float_sub<IsF16>(d, c, s.min2), vmin, vmax);
      const V rd3 = hn::Clamp(float_sub<IsF16>(d, c, s.min3), vmin, vmax);
      const V rd4 = hn::Clamp(float_sub<IsF16>(d, c, s.min4), vmin, vmax);
      const V u = hn::Min(hn::Min(hn::Max(d1, rd1), hn::Max(d2, rd2)),
                          hn::Min(hn::Max(d3, rd3), hn::Max(d4, rd4)));
      const V lo = hn::Clamp(float_sub<IsF16>(d, c, u), vmin, vmax);
      const V hi = hn::Clamp(float_add<IsF16>(d, c, u), vmin, vmax);
      return hn::Clamp(src, lo, hi);
    }
    case 22: {
      const V d1 = float_abs_diff<IsF16>(d, src, g.top_left);
      const V d2 = float_abs_diff<IsF16>(d, src, g.top_center);
      const V d3 = float_abs_diff<IsF16>(d, src, g.top_right);
      const V d4 = float_abs_diff<IsF16>(d, src, g.center_left);
      const V d5 = float_abs_diff<IsF16>(d, src, g.center_right);
      const V d6 = float_abs_diff<IsF16>(d, src, g.bottom_left);
      const V d7 = float_abs_diff<IsF16>(d, src, g.bottom_center);
      const V d8 = float_abs_diff<IsF16>(d, src, g.bottom_right);
      const V mindiff = hn::Min(hn::Min(hn::Min(d1, d2), hn::Min(d3, d4)),
                                hn::Min(hn::Min(d5, d6), hn::Min(d7, d8)));
      const V lo = hn::Clamp(float_sub<IsF16>(d, src, mindiff), vmin, vmax);
      const V hi = hn::Clamp(float_add<IsF16>(d, src, mindiff), vmin, vmax);
      return hn::Clamp(c, lo, hi);
    }
    case 23: {
      const V d1 = float_abs_diff<IsF16>(d, src, g.top_left);
      const V d2 = float_abs_diff<IsF16>(d, src, g.top_center);
      const V d3 = float_abs_diff<IsF16>(d, src, g.top_right);
      const V d4 = float_abs_diff<IsF16>(d, src, g.center_left);
      const V d5 = float_abs_diff<IsF16>(d, src, g.center_right);
      const V d6 = float_abs_diff<IsF16>(d, src, g.bottom_left);
      const V d7 = float_abs_diff<IsF16>(d, src, g.bottom_center);
      const V d8 = float_abs_diff<IsF16>(d, src, g.bottom_right);

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

      const V lo = hn::Clamp(float_sub<IsF16>(d, src, maxdiff), vmin, vmax);
      const V hi = hn::Clamp(float_add<IsF16>(d, src, maxdiff), vmin, vmax);
      return hn::Clamp(c, lo, hi);
    }
    case 24: {
      const auto s = g.min_max_opposites_without_center(d);
      const V d1 = hn::Clamp(float_sub<IsF16>(d, s.max1, src), vmin, vmax);
      const V d2 = hn::Clamp(float_sub<IsF16>(d, s.max2, src), vmin, vmax);
      const V d3 = hn::Clamp(float_sub<IsF16>(d, s.max3, src), vmin, vmax);
      const V d4 = hn::Clamp(float_sub<IsF16>(d, s.max4, src), vmin, vmax);
      const V rd1 = hn::Clamp(float_sub<IsF16>(d, src, s.min1), vmin, vmax);
      const V rd2 = hn::Clamp(float_sub<IsF16>(d, src, s.min2), vmin, vmax);
      const V rd3 = hn::Clamp(float_sub<IsF16>(d, src, s.min3), vmin, vmax);
      const V rd4 = hn::Clamp(float_sub<IsF16>(d, src, s.min4), vmin, vmax);
      const V u = hn::Min(hn::Min(hn::Max(d1, rd1), hn::Max(d2, rd2)),
                          hn::Min(hn::Max(d3, rd3), hn::Max(d4, rd4)));
      const V lo = hn::Clamp(float_sub<IsF16>(d, src, u), vmin, vmax);
      const V hi = hn::Clamp(float_add<IsF16>(d, src, u), vmin, vmax);
      return hn::Clamp(c, lo, hi);
    }
    default:
      return src;
  }
}

template <bool IsF16, typename StorageT>
void repair_float_impl(int mode, bool chroma, const StorageT* srcp, const StorageT* repairp, StorageT* dstp,
                       int width, int height, std::size_t src_stride, std::size_t repair_stride,
                       std::size_t dst_stride) {
  constexpr int kRadius = 1;
  hn::ScalableTag<float> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * kRadius + kSimdPad;
  std::vector<float> row_buffers(checked_product(3, padded_len));
  std::array<float*, 3> rows{
      row_buffers.data() + 0 * padded_len + kRadius,
      row_buffers.data() + 1 * padded_len + kRadius,
      row_buffers.data() + 2 * padded_len + kRadius,
  };
  std::vector<float> src_f32(static_cast<std::size_t>(width) + kSimdPad, 0.0f);
  std::vector<float> out_f32(static_cast<std::size_t>(width) + kSimdPad);

  const float max_val = chroma ? 0.5f : 1.0f;
  const float min_val = chroma ? -0.5f : 0.0f;
  const auto vmax = hn::Set(d, max_val);
  const auto vmin = hn::Set(d, min_val);

  for (int y = 0; y < height; ++y) {
    const StorageT* src_row = srcp + static_cast<std::size_t>(y) * src_stride;
    StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;

    for (int dy = -1; dy <= 1; ++dy) {
      const std::size_t my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      if constexpr (IsF16) {
        fill_mirrored_row_fp16_to_fp32(rows[static_cast<std::size_t>(dy + 1)] - kRadius,
                                       repairp + my * repair_stride, width, kRadius);
      } else {
        fill_mirrored_row(rows[static_cast<std::size_t>(dy + 1)] - kRadius,
                          repairp + my * repair_stride, width, kRadius);
      }
    }

    if constexpr (IsF16) {
      fill_mirrored_row_fp16_to_fp32(src_f32.data(), src_row, width, 0);
    } else {
      std::memcpy(src_f32.data(), src_row, static_cast<std::size_t>(width) * sizeof(float));
    }

    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto g = Grid3x3<decltype(d)>::load(d, rows[0], rows[1], rows[2], x);
      const auto s = hn::LoadU(d, src_f32.data() + x);
      const auto res = eval_repair_float<IsF16>(d, mode, s, g, vmin, vmax);
      hn::StoreU(res, d, out_f32.data() + x);
    }

    if constexpr (IsF16) {
      convert_row_fp32_to_fp16(dst_row, out_f32.data(), width);
    } else {
      std::memcpy(dst_row, out_f32.data(), static_cast<std::size_t>(width) * sizeof(float));
    }
  }
}
