// Included after grid.hpp inside each Highway target namespace.
// Ranking and clamping need no arithmetic headroom: retain the sample width.
template <bool WithCenter, int Mode, class T>
void rank_clamp_plane(const T* src, const T* reference, T* dst,
    int width, int height, std::size_t src_stride, std::size_t reference_stride,
    std::size_t dst_stride) {
  static_assert((Mode >= 1 && Mode <= 4) || (WithCenter && Mode >= 12 && Mode <= 14));
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
      auto lo = hn::Zero(d), hi = hn::Zero(d);
      if constexpr (Mode == 1) {
        lo = g.min_without_center(d);
        hi = g.max_without_center(d);
        if constexpr (WithCenter) {
          lo = hn::Min(lo, g.center_center);
          hi = hn::Max(hi, g.center_center);
        }
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
      auto value = g.center_center;
      if constexpr (WithCenter) value = hn::LoadN(d, src + y * src_stride + x, count);
      hn::StoreN(hn::Clamp(value, lo, hi), d, dst + y * dst_stride + x, count);
    }
  }
}
