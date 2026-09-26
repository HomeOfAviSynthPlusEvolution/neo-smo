// Included after fp16_simd.hpp inside each Highway target namespace.
// Only the few reflected border pixels need scalar conversion. The interior
// and tail use the same bounded vector loads as temporal/weighted kernels.
template <class T>
HWY_INLINE void fill_mirrored_row_f16(T* dst, const std::uint16_t* src,
                                            int width, int radius) {
  const hn::ScalableTag<T> d;
  const auto lanes = hn::Lanes(d);
  for (int x = -radius; x < 0; ++x)
    dst[radius + x] = static_cast<T>(neo_smo::fp16_to_fp32(src[neo_smo::mirror_index(x, width)]));
  for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
    const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
    hn::StoreN(load_f16(d, src + x, count), d, dst + radius + x, count);
  }
  for (std::int64_t x = width; x < static_cast<std::int64_t>(width) + radius; ++x)
    dst[radius + x] = static_cast<T>(neo_smo::fp16_to_fp32(src[neo_smo::mirror_index(x, width)]));
}

template <class T>
HWY_INLINE void store_row_f16(std::uint16_t* dst, const T* src, int width) {
  const hn::ScalableTag<T> d;
  const auto lanes = hn::Lanes(d);
  for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
    const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
    store_f16(d, hn::LoadN(d, src + x, count), dst + x, count);
  }
}
