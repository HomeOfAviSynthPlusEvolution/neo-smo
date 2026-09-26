// Included inside `namespace neo_smo { namespace HWY_NAMESPACE { ... } }` per Highway target.
#include "common/sorting_networks.hpp"

template <class D, class V = hn::Vec<D>>
struct Grid3x3 {
  V top_left;
  V top_center;
  V top_right;
  V center_left;
  V center_center;
  V center_right;
  V bottom_left;
  V bottom_center;
  V bottom_right;

  template <typename T>
  static HWY_INLINE Grid3x3 load(D d, const T* row_top, const T* row_center, const T* row_bottom, std::size_t x) {
    return Grid3x3{
      hn::LoadU(d, row_top + x - 1),
      hn::LoadU(d, row_top + x),
      hn::LoadU(d, row_top + x + 1),
      hn::LoadU(d, row_center + x - 1),
      hn::LoadU(d, row_center + x),
      hn::LoadU(d, row_center + x + 1),
      hn::LoadU(d, row_bottom + x - 1),
      hn::LoadU(d, row_bottom + x),
      hn::LoadU(d, row_bottom + x + 1),
    };
  }

  HWY_INLINE V min_without_center(D d) const {
    (void)d;
    return hn::Min(
      hn::Min(hn::Min(top_left, top_center), hn::Min(top_right, center_left)),
      hn::Min(hn::Min(center_right, bottom_left), hn::Min(bottom_center, bottom_right))
    );
  }

  HWY_INLINE V min_with_center(D d) const {
    return hn::Min(min_without_center(d), center_center);
  }

  HWY_INLINE V max_without_center(D d) const {
    (void)d;
    return hn::Max(
      hn::Max(hn::Max(top_left, top_center), hn::Max(top_right, center_left)),
      hn::Max(hn::Max(center_right, bottom_left), hn::Max(bottom_center, bottom_right))
    );
  }

  HWY_INLINE V max_with_center(D d) const {
    return hn::Max(max_without_center(d), center_center);
  }

  HWY_INLINE void sort_without_center(D d, V (&out)[8]) const {
    out[0] = top_left;
    out[1] = top_center;
    out[2] = top_right;
    out[3] = center_left;
    out[4] = center_right;
    out[5] = bottom_left;
    out[6] = bottom_center;
    out[7] = bottom_right;
    sort8(d, out);
  }

  HWY_INLINE void sort_with_center(D d, V (&out)[9]) const {
    out[0] = top_left;
    out[1] = top_center;
    out[2] = top_right;
    out[3] = center_left;
    out[4] = center_center;
    out[5] = center_right;
    out[6] = bottom_left;
    out[7] = bottom_center;
    out[8] = bottom_right;
    sort9(d, out);
  }

  struct Opposites {
    V max1, min1;
    V max2, min2;
    V max3, min3;
    V max4, min4;
  };

  HWY_INLINE Opposites min_max_opposites_without_center(D d) const {
    (void)d;
    return Opposites{
      hn::Max(top_left, bottom_right), hn::Min(top_left, bottom_right),
      hn::Max(top_center, bottom_center), hn::Min(top_center, bottom_center),
      hn::Max(top_right, bottom_left), hn::Min(top_right, bottom_left),
      hn::Max(center_left, center_right), hn::Min(center_left, center_right),
    };
  }

  HWY_INLINE Opposites min_max_opposites_with_center(D d) const {
    (void)d;
    return Opposites{
      hn::Max(hn::Max(top_left, bottom_right), center_center), hn::Min(hn::Min(top_left, bottom_right), center_center),
      hn::Max(hn::Max(top_center, bottom_center), center_center), hn::Min(hn::Min(top_center, bottom_center), center_center),
      hn::Max(hn::Max(top_right, bottom_left), center_center), hn::Min(hn::Min(top_right, bottom_left), center_center),
      hn::Max(hn::Max(center_left, center_right), center_center), hn::Min(hn::Min(center_left, center_right), center_center),
    };
  }
};
