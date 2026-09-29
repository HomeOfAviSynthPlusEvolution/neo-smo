// Included inside a Highway target namespace.
template <bool Half, bool Interior, class D>
auto deen_float_load(D d, const std::uint8_t* row, int width, std::int64_t x, std::size_t active) {
  if constexpr (!Half)
    return deen_row_load<Interior>(d, row, width, x, active);
  else {
    const hn::Rebind<std::uint16_t, D> du;
    const auto bits = deen_row_load<Interior>(du, row, width, x, active);
#if HWY_HAVE_FLOAT16 || (HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX2 && !defined(HWY_DISABLE_F16C))
    const hn::Rebind<hwy::float16_t, D> dh;
    return hn::PromoteTo(d, hn::BitCast(dh, bits));
#else
    HWY_ALIGN std::uint16_t values[hn::MaxLanes(d)]{};
    hn::StoreU(bits, du, values);
    return load_f16(d, values, active);
#endif
  }
}

template <bool Half, class D>
void deen_float_store(D d, hn::Vec<D> value, std::uint8_t* dst, std::size_t active) {
  if constexpr (!Half)
    deen_store(d, value, dst, active);
  else {
#if HWY_HAVE_FLOAT16 || (HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX2 && !defined(HWY_DISABLE_F16C))
    const hn::Rebind<hwy::float16_t, D> dh;
    deen_store(dh, hn::DemoteTo(dh, value), dst, active);
#else
    HWY_ALIGN std::uint16_t values[hn::MaxLanes(d)]{};
    store_f16(d, value, values, active);
    std::memcpy(dst, values, active * sizeof(std::uint16_t));
#endif
  }
}
