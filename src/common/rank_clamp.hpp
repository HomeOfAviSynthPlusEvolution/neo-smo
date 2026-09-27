// Included after grid.hpp inside each Highway target namespace.
// Ranking and clamping need no arithmetic headroom: retain the sample width.
template <bool WithCenter, int Mode, class T>
void rank_clamp_plane(const T* src, const T* reference, T* dst,
    int width, int height, std::size_t src_stride, std::size_t reference_stride,
    std::size_t dst_stride) {
  static_assert((Mode >= 1 && Mode <= 4) || Mode == 5 || Mode == 6 || Mode == 8 || Mode == 9 || Mode == 17 || Mode == 18 || (!WithCenter && Mode == 22) || (WithCenter && ((Mode >= 12 && Mode <= 16) || (Mode >= 19 && Mode <= 24))));
  const hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 + lanes;
  std::vector<T> buffer(checked_product(3, padded_len));
  std::array<int, 3> cached_y{-1, -1, -1};
  std::array<T*, 3> rows{};
  for (int y = 0; y < height; ++y) {
    for (int dy = -1; dy <= 1; ++dy) {
      const auto my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      const auto slot = my % 3;
      auto* row = buffer.data() + slot * padded_len + 1;
      if (cached_y[slot] != static_cast<int>(my)) {
        const auto* input = reference + my * reference_stride;
        std::copy_n(input, width, row);
        row[-1] = input[mirror_index(-1, width)];
        row[width] = input[mirror_index(width, width)];
        cached_y[slot] = static_cast<int>(my);
      }
      rows[dy + 1] = row;
    }
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
      const auto g = Grid3x3<decltype(d)>::load(d, rows[0], rows[1], rows[2], x);
      auto value = g.center_center;
      if constexpr (WithCenter) value = hn::LoadN(d, src + y * src_stride + x, count);
      auto lo = hn::Zero(d), hi = hn::Zero(d);
      if constexpr (Mode == 1) {
        lo = g.min_without_center(d);
        hi = g.max_without_center(d);
        if constexpr (WithCenter) {
          lo = hn::Min(lo, g.center_center);
          hi = hn::Max(hi, g.center_center);
        }
      } else if constexpr (Mode == 5 || Mode == 6 || Mode == 8 || Mode == 9 || Mode == 15 || Mode == 16 || Mode == 18) {
        auto pairs = g.min_max_opposites_without_center(d);
        if constexpr (WithCenter && Mode <= 9) pairs = g.min_max_opposites_with_center(d);
        const auto pivot = Mode >= 15 ? g.center_center : value;
        auto pair_cost = [&](auto low, auto high) HWY_ATTR {
          const auto range = hn::Sub(high, low);
          if constexpr (Mode == 9) return range;
          else if constexpr (Mode == 18) return hn::Max(
              hn::Sub(hn::Max(pivot, low), hn::Min(pivot, low)),
              hn::Sub(hn::Max(pivot, high), hn::Min(pivot, high)));
          else {
            const auto clamped = hn::Clamp(pivot, low, high);
            const auto delta = hn::Sub(hn::Max(pivot, clamped), hn::Min(pivot, clamped));
            if constexpr (Mode == 5 || Mode == 15) return delta;
            else if constexpr (Mode == 8) return hn::SaturatedAdd(delta, hn::SaturatedAdd(range, range));
            else return hn::SaturatedAdd(hn::SaturatedAdd(delta, delta), range);
          }
        };
        const auto c1 = pair_cost(pairs.min1, pairs.max1), c2 = pair_cost(pairs.min2, pairs.max2);
        const auto c3 = pair_cost(pairs.min3, pairs.max3), c4 = pair_cost(pairs.min4, pairs.max4);
        const auto cost = hn::Min(hn::Min(c1, c2), hn::Min(c3, c4));
        // Match upstream tie priority: horizontal, vertical, diagonal 3, diagonal 1.
        lo = pairs.min1; hi = pairs.max1;
        lo = hn::IfThenElse(hn::Eq(cost, c3), pairs.min3, lo);
        hi = hn::IfThenElse(hn::Eq(cost, c3), pairs.max3, hi);
        lo = hn::IfThenElse(hn::Eq(cost, c2), pairs.min2, lo);
        hi = hn::IfThenElse(hn::Eq(cost, c2), pairs.max2, hi);
        lo = hn::IfThenElse(hn::Eq(cost, c4), pairs.min4, lo);
        hi = hn::IfThenElse(hn::Eq(cost, c4), pairs.max4, hi);
        if constexpr (WithCenter && Mode >= 15) {
          lo = hn::Min(lo, g.center_center);
          hi = hn::Max(hi, g.center_center);
        }
      } else if constexpr (!WithCenter && Mode == 22) {
        const auto a = hn::AverageRound(g.top_left, g.bottom_right);
        const auto b = hn::AverageRound(g.top_center, g.bottom_center);
        const auto c = hn::AverageRound(g.top_right, g.bottom_left);
        const auto e = hn::AverageRound(g.center_left, g.center_right);
        lo = hn::Min(hn::Min(a, b), hn::Min(c, e));
        hi = hn::Max(hn::Max(a, b), hn::Max(c, e));
      } else if constexpr (Mode == 17) {
        const auto pairs = g.min_max_opposites_without_center(d);
        const auto lower = hn::Max(hn::Max(pairs.min1, pairs.min2), hn::Max(pairs.min3, pairs.min4));
        const auto upper = hn::Min(hn::Min(pairs.max1, pairs.max2), hn::Min(pairs.max3, pairs.max4));
        lo = hn::Min(lower, upper);
        hi = hn::Max(lower, upper);
        if constexpr (WithCenter) {
          lo = hn::Min(lo, g.center_center);
          hi = hn::Max(hi, g.center_center);
        }
      } else if constexpr (Mode >= 19) {
        const auto center = Mode >= 22 ? value : g.center_center;
        const hn::Vec<decltype(d)> neighbors[8] = {g.top_left, g.top_center, g.top_right,
            g.center_left, g.center_right, g.bottom_left, g.bottom_center, g.bottom_right};
        hn::Vec<decltype(d)> diff[8];
        for (int i = 0; i < 8; ++i) diff[i] = hn::Sub(hn::Max(center, neighbors[i]), hn::Min(center, neighbors[i]));
        auto radius = hn::Zero(d);
        if constexpr (Mode == 21 || Mode == 24) {
          radius = hn::Min(hn::Min(hn::Max(diff[0], diff[7]), hn::Max(diff[1], diff[6])),
                           hn::Min(hn::Max(diff[2], diff[5]), hn::Max(diff[3], diff[4])));
        } else {
          const auto a = hn::Min(diff[0], diff[1]), b = hn::Min(diff[2], diff[3]);
          const auto c = hn::Min(diff[4], diff[5]), e = hn::Min(diff[6], diff[7]);
          if constexpr (Mode == 20 || Mode == 23) {
            // Merge pairs of two smallest distances; ties count as distinct samples.
            const auto second_ab = hn::Min(hn::Max(a, b),
                hn::Min(hn::Max(diff[0], diff[1]), hn::Max(diff[2], diff[3])));
            const auto second_ce = hn::Min(hn::Max(c, e),
                hn::Min(hn::Max(diff[4], diff[5]), hn::Max(diff[6], diff[7])));
            radius = hn::Min(hn::Max(hn::Min(a, b), hn::Min(c, e)), hn::Min(second_ab, second_ce));
          } else {
            radius = hn::Min(hn::Min(a, b), hn::Min(c, e));
          }
        }
        lo = hn::SaturatedSub(center, radius);
        hi = hn::SaturatedAdd(center, radius);
        if constexpr (Mode >= 22) value = g.center_center;
      } else if constexpr (Mode >= 12) {
        hn::Vec<decltype(d)> sorted[8];
        g.sort_without_center(d, sorted);
        lo = hn::Min(sorted[Mode - 11], g.center_center);
        hi = hn::Max(sorted[18 - Mode], g.center_center);
      } else if constexpr (WithCenter) {
        hn::Vec<decltype(d)> sorted[9];
        g.sort_with_center(d, sorted);
        lo = sorted[Mode - 1]; hi = sorted[9 - Mode];
      } else {
        hn::Vec<decltype(d)> sorted[8];
        g.sort_without_center(d, sorted);
        lo = sorted[Mode - 1]; hi = sorted[8 - Mode];
      }
      hn::StoreN(hn::Clamp(value, lo, hi), d, dst + y * dst_stride + x, count);
    }
  }
}
