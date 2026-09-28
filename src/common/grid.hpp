// Included inside `namespace neo_smo { namespace HWY_NAMESPACE { ... } }` per Highway target.
#include "common/sorting_networks.hpp"

template <class D, class V = hn::Vec<D>>
struct Grid3x3 {
  // Direct brace initialization extends loaded-vector temporaries to the grid's
  // lifetime on scalable targets. Do not return these reference-backed grids.
  using Member = std::conditional_t<NEO_SMO_SIZELESS, const V&, V>;
  Member top_left;
  Member top_center;
  Member top_right;
  Member center_left;
  Member center_center;
  Member center_right;
  Member bottom_left;
  Member bottom_center;
  Member bottom_right;

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

  HWY_INLINE void sort_without_center(D d, VectorArrayView<V> out) const {
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

  HWY_INLINE void sort_with_center(D d, VectorArrayView<V> out) const {
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

#if NEO_SMO_SIZELESS
  struct Opposites {
    const Grid3x3* grid;
    bool with_center;
    HWY_INLINE V max1() const {
      const auto value = hn::Max(grid->top_left, grid->bottom_right);
      return with_center ? hn::Max(value, grid->center_center) : value;
    }
    HWY_INLINE V min1() const {
      const auto value = hn::Min(grid->top_left, grid->bottom_right);
      return with_center ? hn::Min(value, grid->center_center) : value;
    }
    HWY_INLINE V max2() const {
      const auto value = hn::Max(grid->top_center, grid->bottom_center);
      return with_center ? hn::Max(value, grid->center_center) : value;
    }
    HWY_INLINE V min2() const {
      const auto value = hn::Min(grid->top_center, grid->bottom_center);
      return with_center ? hn::Min(value, grid->center_center) : value;
    }
    HWY_INLINE V max3() const {
      const auto value = hn::Max(grid->top_right, grid->bottom_left);
      return with_center ? hn::Max(value, grid->center_center) : value;
    }
    HWY_INLINE V min3() const {
      const auto value = hn::Min(grid->top_right, grid->bottom_left);
      return with_center ? hn::Min(value, grid->center_center) : value;
    }
    HWY_INLINE V max4() const {
      const auto value = hn::Max(grid->center_left, grid->center_right);
      return with_center ? hn::Max(value, grid->center_center) : value;
    }
    HWY_INLINE V min4() const {
      const auto value = hn::Min(grid->center_left, grid->center_right);
      return with_center ? hn::Min(value, grid->center_center) : value;
    }
  };
  HWY_INLINE Opposites min_max_opposites_without_center(D) const { return {this, false}; }
  HWY_INLINE Opposites min_max_opposites_with_center(D) const { return {this, true}; }
#else
  struct Opposites {
    V max1_, min1_;
    V max2_, min2_;
    V max3_, min3_;
    V max4_, min4_;
    HWY_INLINE V max1() const { return max1_; }
    HWY_INLINE V min1() const { return min1_; }
    HWY_INLINE V max2() const { return max2_; }
    HWY_INLINE V min2() const { return min2_; }
    HWY_INLINE V max3() const { return max3_; }
    HWY_INLINE V min3() const { return min3_; }
    HWY_INLINE V max4() const { return max4_; }
    HWY_INLINE V min4() const { return min4_; }
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
#endif
};

// Expand at the declaration, rather than returning a reference-backed aggregate
// from a helper: the vector temporaries must outlive all uses of the grid->
#undef NEO_SMO_LOAD_GRID
#define NEO_SMO_LOAD_GRID(d, top, center, bottom, x)                                                                   \
  {hn::LoadU(d, (top) + (x) - 1),    hn::LoadU(d, (top) + (x)),    hn::LoadU(d, (top) + (x) + 1),                      \
   hn::LoadU(d, (center) + (x) - 1), hn::LoadU(d, (center) + (x)), hn::LoadU(d, (center) + (x) + 1),                   \
   hn::LoadU(d, (bottom) + (x) - 1), hn::LoadU(d, (bottom) + (x)), hn::LoadU(d, (bottom) + (x) + 1)}
