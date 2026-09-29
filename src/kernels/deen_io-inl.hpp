// Included only inside hardware SIMD target namespaces; software targets use plain C++.
// Access byte-backed host planes without assuming natural sample alignment.
template <class D>
auto deen_raw_load(D d, const std::uint8_t* ptr, std::size_t active) {
  const hn::Repartition<std::uint8_t, D> db;
  return hn::BitCast(d,
                     active == hn::Lanes(d) ? hn::LoadU(db, ptr) : hn::LoadN(db, ptr, active * sizeof(hn::TFromD<D>)));
}
template <class D>
void deen_store([[maybe_unused]] D d, hn::VFromD<D> value, std::uint8_t* ptr, std::size_t active) {
  const hn::Repartition<std::uint8_t, D> db;
  hn::StoreN(hn::BitCast(db, value), db, ptr, active * sizeof(hn::TFromD<D>));
}
// Clamp coordinates only at the actual image edges. Interior taps load original pixels.
template <bool Interior = false, class D>
auto deen_row_load(D d, const std::uint8_t* row, int width, std::int64_t x, std::size_t active) {
  using T = hn::TFromD<D>;
  if constexpr (Interior) {
    return deen_raw_load(d, row + static_cast<std::size_t>(x) * sizeof(T), hn::Lanes(d));
  }
  if (x >= 0 && x <= width - static_cast<int>(active)) {
    return deen_raw_load(d, row + static_cast<std::size_t>(x) * sizeof(T), active);
  }
  HWY_ALIGN T values[hn::MaxLanes(d)]{};
  for (std::size_t i = 0; i < active; ++i)
    std::memcpy(values + i,
                row + static_cast<std::size_t>(std::clamp(static_cast<std::int64_t>(x) + static_cast<std::int64_t>(i),
                                                          std::int64_t{0}, static_cast<std::int64_t>(width - 1))) *
                          sizeof(T),
                sizeof(T));
  return hn::Load(d, values);
}
