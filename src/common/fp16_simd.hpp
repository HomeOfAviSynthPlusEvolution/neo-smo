// FP16 arithmetic is native only when the current Highway target supports it.
template <bool Half>
using FloatLane = std::conditional_t<Half && HWY_HAVE_FLOAT16, hwy::float16_t, float>;

// Included inside each Highway target namespace. No include guard.
// Native conversion is available on F16C/FP16 targets. The portable vector path
// explicitly handles the full binary16 domain, including subnormals and NaNs.
template <class D>
#if defined(NEO_SMO_OUTLINE_EMU_ROUND_F16) && HWY_TARGET == HWY_EMU128
HWY_NOINLINE
#else
HWY_INLINE
#endif
hn::Vec<D> round_f16(D d, hn::Vec<D> v) {
#if HWY_HAVE_FLOAT16 || (HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX2 && !defined(HWY_DISABLE_F16C))
  const hn::Rebind<hwy::float16_t, D> dh;
  return hn::PromoteTo(d, hn::DemoteTo(dh, v));
#else
  const hn::RebindToUnsigned<D> du;
  const auto bits = hn::BitCast(du, v);
  const auto sign = hn::And(bits, hn::Set(du, 0x80000000u));
  const auto magnitude = hn::And(bits, hn::Set(du, 0x7fffffffu));
  const auto odd = hn::And(hn::ShiftRight<13>(magnitude), hn::Set(du, 1u));
  auto rounded = hn::And(hn::Add(magnitude, hn::Add(hn::Set(du, 0x0fffu), odd)), hn::Set(du, 0xffffe000u));

  // Half subnormals have a fixed 2^-24 quantum (not ten fractional bits).
  // Mask first so inactive infinity/NaN/large lanes do not enter arithmetic.
  const auto tiny = hn::Lt(magnitude, hn::Set(du, 0x38800000u));
  const auto small = hn::BitCast(d, hn::IfThenElseZero(tiny, magnitude));
  // An integer conversion keeps this quantization explicit under fast math;
  // a floating-point Round implementation can lose it through reassociation.
  const auto quantum = hn::NearestInt(hn::Mul(small, hn::Set(d, 0x1p24f)));
  const auto subnormal = hn::Mul(hn::ConvertTo(d, quantum), hn::Set(d, 0x1p-24f));
  rounded = hn::IfThenElse(tiny, hn::BitCast(du, subnormal), rounded);
  rounded = hn::IfThenElse(hn::Ge(magnitude, hn::Set(du, 0x477ff000u)), hn::Set(du, 0x7f800000u), rounded);
  // Keep NaNs as NaNs, including payloads that would otherwise round to inf.
  rounded =
      hn::IfThenElse(hn::Gt(magnitude, hn::Set(du, 0x7f800000u)), hn::Or(magnitude, hn::Set(du, 0x00400000u)), rounded);
  return hn::BitCast(d, hn::Or(sign, rounded));
#endif
}

template <class D>
HWY_INLINE hn::Vec<D> load_f16(D d, const std::uint16_t* src, std::size_t count) {
  const hn::Rebind<std::uint16_t, D> du16;
#if HWY_HAVE_FLOAT16
  if constexpr (std::is_same_v<hn::TFromD<D>, hwy::float16_t>) {
    return hn::BitCast(d, hn::LoadN(du16, src, count));
  } else
#endif
  {
  const auto half = hn::LoadN(du16, src, count);
#if HWY_HAVE_FLOAT16 || (HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX2 && !defined(HWY_DISABLE_F16C))
  const hn::Rebind<hwy::float16_t, D> dh;
  return hn::PromoteTo(d, hn::BitCast(dh, half));
#else
  const hn::Rebind<std::uint32_t, D> du;
  const auto bits = hn::PromoteTo(du, half);
  const auto sign = hn::ShiftLeft<16>(hn::And(bits, hn::Set(du, 0x8000)));
  const auto exp = hn::And(bits, hn::Set(du, 0x7c00));
  const auto mant = hn::And(bits, hn::Set(du, 0x03ff));
  auto out = hn::Add(hn::ShiftLeft<13>(hn::And(bits, hn::Set(du, 0x7fff))), hn::Set(du, 0x38000000));
  const auto subnormal = hn::Mul(hn::ConvertTo(d, mant), hn::Set(d, 0x1p-24f));
  out = hn::IfThenElse(hn::Eq(exp, hn::Zero(du)), hn::BitCast(du, subnormal), out);
  out = hn::IfThenElse(hn::Eq(exp, hn::Set(du, 0x7c00)),
                      hn::Or(hn::Set(du, 0x7f800000), hn::ShiftLeft<13>(mant)), out);
  return hn::BitCast(d, hn::Or(sign, out));
#endif
  }
}

template <class D>
HWY_INLINE void store_f16([[maybe_unused]] D d, hn::Vec<D> value, std::uint16_t* dst, std::size_t count) {
  const hn::Rebind<std::uint16_t, D> du16;
#if HWY_HAVE_FLOAT16
  if constexpr (std::is_same_v<hn::TFromD<D>, hwy::float16_t>) {
    hn::StoreN(hn::BitCast(du16, value), du16, dst, count);
  } else
#endif
  {
#if HWY_HAVE_FLOAT16 || (HWY_ARCH_X86 && HWY_TARGET <= HWY_AVX2 && !defined(HWY_DISABLE_F16C))
  const hn::Rebind<hwy::float16_t, D> dh;
  const auto half = hn::BitCast(du16, hn::DemoteTo(dh, value));
#else
  const hn::Rebind<std::uint32_t, D> du;
  const auto bits = hn::BitCast(du, round_f16(d, value));
  const auto sign = hn::And(hn::ShiftRight<16>(bits), hn::Set(du, 0x8000));
  const auto mag = hn::And(bits, hn::Set(du, 0x7fffffff));
  const auto tiny = hn::Lt(mag, hn::Set(du, 0x38800000));
  const auto small = hn::BitCast(d, hn::IfThenElseZero(tiny, mag));
  const auto subnormal = hn::ConvertTo(du, hn::Mul(small, hn::Set(d, 0x1p24f)));
  auto out = hn::Sub(hn::ShiftRight<13>(mag), hn::Set(du, 0x1c000));
  out = hn::IfThenElse(tiny, subnormal, out);
  const auto special = hn::Or(hn::Set(du, 0x7c00),
      hn::And(hn::ShiftRight<13>(mag), hn::Set(du, 0x03ff)));
  out = hn::IfThenElse(hn::Ge(mag, hn::Set(du, 0x7f800000)), special, out);
  const auto half = hn::DemoteTo(du16, hn::Or(sign, out));
#endif
  hn::StoreN(half, du16, dst, count);
  }
}
