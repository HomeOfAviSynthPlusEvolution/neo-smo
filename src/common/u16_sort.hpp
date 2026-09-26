// Included inside each Highway target namespace, after sorting_networks.hpp.
#if defined(NEO_SMO_MSVC_U16_SORT) && (HWY_TARGET == HWY_SSE2 || HWY_TARGET == HWY_SSSE3)
// These networks receive only U8/U16 samples widened to int32. Bias into the
// signed 16-bit range so old SSE can use native word Min/Max. The upper word
// is the sign extension and remains consistent with the sorted lower word.
// Converting once around the network avoids costly repeated int32 Min/Max
// expansion in MSVC; all subsequent filter arithmetic remains int32.
template <int Count, class D, class Sort>
HWY_INLINE void sort_u16_samples(D d, hn::Vec<D>* values, Sort sort) {
  const hn::Repartition<std::int16_t, D> d16;
  hn::Vec<decltype(d16)> words[Count];
  const auto bias = hn::Set(d, 32768);
  for (int i = 0; i < Count; ++i)
    words[i] = hn::BitCast(d16, hn::Sub(values[i], bias));
  sort(d16, words);
  for (int i = 0; i < Count; ++i)
    values[i] = hn::Add(hn::BitCast(d, words[i]), bias);
}
#endif
