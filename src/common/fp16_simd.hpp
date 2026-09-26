// Included inside the per-target Highway namespace.
// Round to binary16, keeping binary32 storage for subsequent arithmetic.
template <class D>
#if defined(NEO_SMO_OUTLINE_EMU_ROUND_F16) && HWY_TARGET == HWY_EMU128
HWY_NOINLINE
#else
HWY_INLINE
#endif
hn::Vec<D> round_f16(D d, hn::Vec<D> v) {
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
}
