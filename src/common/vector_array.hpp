// Included once inside each Highway target namespace, before sorting helpers.
// Scalable vectors cannot be array elements. Keep them as local variables and
// index a non-owning view of their addresses instead. Constant-index sorting
// networks inline through the view, allowing the vectors to stay in registers.
// Fixed-width targets retain their original native vector arrays.
#undef NEO_SMO_SIZELESS
#define NEO_SMO_SIZELESS (HWY_HAVE_SCALABLE || HWY_TARGET == HWY_SVE_256 || HWY_TARGET == HWY_SVE2_128)
#undef NEO_SMO_VECTOR_ARRAY
#if NEO_SMO_SIZELESS
template <class V>
struct VectorArrayView {
  V* const* entries;
  HWY_INLINE V& operator[](std::size_t i) const { return *entries[i]; }
};

template <std::size_t Count, class V, std::size_t Capacity>
HWY_INLINE std::array<V*, Count> vector_pointer_prefix(const std::array<V*, Capacity>& pointers) {
  static_assert(Count <= Capacity, "Vector array exceeds the largest sorting network");
  std::array<V*, Count> result;
  for (std::size_t i = 0; i < Count; ++i)
    result[i] = pointers[i];
  return result;
}

// The backing variables must live in the caller's scope. Never return this view
// or retain it beyond that scope. Only the first Count addresses escape; unused
// variables and the temporary address list are removed by optimization.
#define NEO_SMO_VECTOR_ARRAY(D, name, Count)                                                                           \
  hn::Vec<D> name##_v0;                                                                                                \
  hn::Vec<D> name##_v1;                                                                                                \
  hn::Vec<D> name##_v2;                                                                                                \
  hn::Vec<D> name##_v3;                                                                                                \
  hn::Vec<D> name##_v4;                                                                                                \
  hn::Vec<D> name##_v5;                                                                                                \
  hn::Vec<D> name##_v6;                                                                                                \
  hn::Vec<D> name##_v7;                                                                                                \
  hn::Vec<D> name##_v8;                                                                                                \
  hn::Vec<D> name##_v9;                                                                                                \
  hn::Vec<D> name##_v10;                                                                                               \
  hn::Vec<D> name##_v11;                                                                                               \
  hn::Vec<D> name##_v12;                                                                                               \
  hn::Vec<D> name##_v13;                                                                                               \
  hn::Vec<D> name##_v14;                                                                                               \
  hn::Vec<D> name##_v15;                                                                                               \
  hn::Vec<D> name##_v16;                                                                                               \
  hn::Vec<D> name##_v17;                                                                                               \
  hn::Vec<D> name##_v18;                                                                                               \
  hn::Vec<D> name##_v19;                                                                                               \
  hn::Vec<D> name##_v20;                                                                                               \
  hn::Vec<D> name##_v21;                                                                                               \
  hn::Vec<D> name##_v22;                                                                                               \
  hn::Vec<D> name##_v23;                                                                                               \
  hn::Vec<D> name##_v24;                                                                                               \
  hn::Vec<D> name##_v25;                                                                                               \
  hn::Vec<D> name##_v26;                                                                                               \
  hn::Vec<D> name##_v27;                                                                                               \
  hn::Vec<D> name##_v28;                                                                                               \
  hn::Vec<D> name##_v29;                                                                                               \
  hn::Vec<D> name##_v30;                                                                                               \
  hn::Vec<D> name##_v31;                                                                                               \
  hn::Vec<D> name##_v32;                                                                                               \
  hn::Vec<D> name##_v33;                                                                                               \
  hn::Vec<D> name##_v34;                                                                                               \
  hn::Vec<D> name##_v35;                                                                                               \
  hn::Vec<D> name##_v36;                                                                                               \
  hn::Vec<D> name##_v37;                                                                                               \
  hn::Vec<D> name##_v38;                                                                                               \
  hn::Vec<D> name##_v39;                                                                                               \
  hn::Vec<D> name##_v40;                                                                                               \
  hn::Vec<D> name##_v41;                                                                                               \
  hn::Vec<D> name##_v42;                                                                                               \
  hn::Vec<D> name##_v43;                                                                                               \
  hn::Vec<D> name##_v44;                                                                                               \
  hn::Vec<D> name##_v45;                                                                                               \
  hn::Vec<D> name##_v46;                                                                                               \
  hn::Vec<D> name##_v47;                                                                                               \
  hn::Vec<D> name##_v48;                                                                                               \
  auto name##_pointers = vector_pointer_prefix<Count>(                                                                 \
      std::array{&name##_v0,  &name##_v1,  &name##_v2,  &name##_v3,  &name##_v4,  &name##_v5,  &name##_v6,             \
                 &name##_v7,  &name##_v8,  &name##_v9,  &name##_v10, &name##_v11, &name##_v12, &name##_v13,            \
                 &name##_v14, &name##_v15, &name##_v16, &name##_v17, &name##_v18, &name##_v19, &name##_v20,            \
                 &name##_v21, &name##_v22, &name##_v23, &name##_v24, &name##_v25, &name##_v26, &name##_v27,            \
                 &name##_v28, &name##_v29, &name##_v30, &name##_v31, &name##_v32, &name##_v33, &name##_v34,            \
                 &name##_v35, &name##_v36, &name##_v37, &name##_v38, &name##_v39, &name##_v40, &name##_v41,            \
                 &name##_v42, &name##_v43, &name##_v44, &name##_v45, &name##_v46, &name##_v47, &name##_v48});          \
  VectorArrayView<hn::Vec<D>> name {                                                                                   \
    name##_pointers.data()                                                                                             \
  }
#else
template <class V>
using VectorArrayView = V*;
#define NEO_SMO_VECTOR_ARRAY(D, name, Count) hn::Vec<D> name[Count]
#endif
