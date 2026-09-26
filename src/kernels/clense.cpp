#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/clense.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/sorting_networks.hpp"
#include "common/float_arithmetic.hpp"

template <typename T>
void clense_int_impl(const T* srcp, const T* prevp, const T* nextp, T* dstp,
                     std::size_t width, std::size_t height,
                     std::size_t src_stride, std::size_t prev_stride,
                     std::size_t next_stride, std::size_t dst_stride) {
  const hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);

  for (std::size_t y = 0; y < height; ++y) {
    const T* r_src = srcp + y * src_stride;
    const T* r_prev = prevp + y * prev_stride;
    const T* r_next = nextp + y * next_stride;
    T* r_dst = dstp + y * dst_stride;

    for (std::size_t x = 0; x < width; x += lanes) {
      const std::size_t rem = width - x;
      const std::size_t count = std::min(lanes, rem);
      const auto s = hn::LoadN(d, r_src + x, count);
      const auto p = hn::LoadN(d, r_prev + x, count);
      const auto n = hn::LoadN(d, r_next + x, count);
      const auto res = median3(d, p, s, n);
      if (rem >= lanes)
        hn::StoreU(res, d, r_dst + x);
      else
        hn::StoreN(res, d, r_dst + x, rem);
    }
  }
}

template <bool IsF16, typename StorageT>
void clense_float_impl(const StorageT* srcp, const StorageT* prevp, const StorageT* nextp, StorageT* dstp,
                       std::size_t width, std::size_t height,
                       std::size_t src_stride, std::size_t prev_stride,
                       std::size_t next_stride, std::size_t dst_stride) {
  using ComputeT = FloatLane<IsF16>;
  const hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);

  if constexpr (IsF16) {
    for (std::size_t y = 0; y < height; ++y) {
      const StorageT* r_src = srcp + y * src_stride;
      const StorageT* r_prev = prevp + y * prev_stride;
      const StorageT* r_next = nextp + y * next_stride;
      StorageT* r_dst = dstp + y * dst_stride;

      for (std::size_t x = 0; x < width; x += lanes) {
        const auto count_lanes = std::min(lanes, static_cast<std::size_t>(width) - x);
        const auto s = load_f16(d, r_src + x, count_lanes);
        const auto p = load_f16(d, r_prev + x, count_lanes);
        const auto n = load_f16(d, r_next + x, count_lanes);
        const auto res = median3(d, p, s, n);
        store_f16(d, res, r_dst + x, count_lanes);
      }
    }
  } else {
    for (std::size_t y = 0; y < height; ++y) {
      const float* r_src = srcp + y * src_stride;
      const float* r_prev = prevp + y * prev_stride;
      const float* r_next = nextp + y * next_stride;
      float* r_dst = dstp + y * dst_stride;

      for (std::size_t x = 0; x < width; x += lanes) {
        const std::size_t rem = width - x;
        const std::size_t count = std::min(lanes, rem);
        const auto s = hn::LoadN(d, r_src + x, count);
        const auto p = hn::LoadN(d, r_prev + x, count);
        const auto n = hn::LoadN(d, r_next + x, count);
        const auto res = median3(d, p, s, n);
        if (rem >= lanes)
          hn::StoreU(res, d, r_dst + x);
        else
          hn::StoreN(res, d, r_dst + x, rem);
      }
    }
  }
}

template <typename T>
void clense_forward_backward_int_impl(const T* srcp, const T* ref1p, const T* ref2p, T* dstp,
                                      std::size_t width, std::size_t height,
                                      std::size_t src_stride, std::size_t ref1_stride,
                                      std::size_t ref2_stride, std::size_t dst_stride) {
  const hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);

  for (std::size_t y = 0; y < height; ++y) {
    const T* r_src = srcp + y * src_stride;
    const T* r_ref1 = ref1p + y * ref1_stride;
    const T* r_ref2 = ref2p + y * ref2_stride;
    T* r_dst = dstp + y * dst_stride;

    for (std::size_t x = 0; x < width; x += lanes) {
      const std::size_t rem = width - x;
      const std::size_t count = std::min(lanes, rem);
      const auto s = hn::LoadN(d, r_src + x, count);
      const auto ref1 = hn::LoadN(d, r_ref1 + x, count);
      const auto ref2 = hn::LoadN(d, r_ref2 + x, count);

      const auto minref = hn::Min(ref1, ref2);
      const auto maxref = hn::Max(ref1, ref2);

      const auto lowref = hn::SaturatedSub(minref, hn::SaturatedSub(ref2, minref));
      const auto highref = hn::SaturatedAdd(hn::SaturatedSub(maxref, ref2), maxref);

      const auto res = hn::Clamp(s, lowref, highref);

      if (rem >= lanes)
        hn::StoreU(res, d, r_dst + x);
      else
        hn::StoreN(res, d, r_dst + x, rem);
    }
  }
}

template <bool IsF16, typename StorageT>
void clense_forward_backward_float_impl(const StorageT* srcp, const StorageT* ref1p, const StorageT* ref2p, StorageT* dstp,
                                        std::size_t width, std::size_t height,
                                        std::size_t src_stride, std::size_t ref1_stride,
                                        std::size_t ref2_stride, std::size_t dst_stride) {
  using ComputeT = FloatLane<IsF16>;
  const hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);
  const auto two = hn::Set(d, 2.0f);

  if constexpr (IsF16) {
    for (std::size_t y = 0; y < height; ++y) {
      const StorageT* r_src = srcp + y * src_stride;
      const StorageT* r_ref1 = ref1p + y * ref1_stride;
      const StorageT* r_ref2 = ref2p + y * ref2_stride;
      StorageT* r_dst = dstp + y * dst_stride;

      for (std::size_t x = 0; x < width; x += lanes) {
        const auto count_lanes = std::min(lanes, static_cast<std::size_t>(width) - x);
        const auto s = load_f16(d, r_src + x, count_lanes);
        const auto ref1 = load_f16(d, r_ref1 + x, count_lanes);
        const auto ref2 = load_f16(d, r_ref2 + x, count_lanes);

        const auto minref = hn::Min(ref1, ref2);
        const auto maxref = hn::Max(ref1, ref2);

        const auto lowref = float_sub<IsF16>(d, float_mul<IsF16>(d, minref, two), ref2);
        const auto highref = float_sub<IsF16>(d, float_mul<IsF16>(d, maxref, two), ref2);

        const auto res = hn::Clamp(s, lowref, highref);
        store_f16(d, res, r_dst + x, count_lanes);
      }
    }
  } else {
    for (std::size_t y = 0; y < height; ++y) {
      const float* r_src = srcp + y * src_stride;
      const float* r_ref1 = ref1p + y * ref1_stride;
      const float* r_ref2 = ref2p + y * ref2_stride;
      float* r_dst = dstp + y * dst_stride;

      for (std::size_t x = 0; x < width; x += lanes) {
        const std::size_t rem = width - x;
        const std::size_t count = std::min(lanes, rem);
        const auto s = hn::LoadN(d, r_src + x, count);
        const auto ref1 = hn::LoadN(d, r_ref1 + x, count);
        const auto ref2 = hn::LoadN(d, r_ref2 + x, count);

        const auto minref = hn::Min(ref1, ref2);
        const auto maxref = hn::Max(ref1, ref2);

        const auto lowref = hn::Sub(hn::Mul(minref, two), ref2);
        const auto highref = hn::Sub(hn::Mul(maxref, two), ref2);

        const auto res = hn::Clamp(s, lowref, highref);

        if (rem >= lanes)
          hn::StoreU(res, d, r_dst + x);
        else
          hn::StoreN(res, d, r_dst + x, rem);
      }
    }
  }
}

void dispatch_clense_target(DataType dtype, const std::uint8_t* srcp, const std::uint8_t* prevp,
                            const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                            std::size_t height, std::size_t src_stride_bytes,
                            std::size_t prev_stride_bytes, std::size_t next_stride_bytes,
                            std::size_t dst_stride_bytes) {
  if (dtype == DataType::U8) {
    clense_int_impl<std::uint8_t>(srcp, prevp, nextp, dstp, width, height,
                                  src_stride_bytes, prev_stride_bytes, next_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    clense_int_impl<std::uint16_t>(reinterpret_cast<const std::uint16_t*>(srcp),
                                   reinterpret_cast<const std::uint16_t*>(prevp),
                                   reinterpret_cast<const std::uint16_t*>(nextp),
                                   reinterpret_cast<std::uint16_t*>(dstp),
                                   width, height,
                                   src_stride_bytes / 2, prev_stride_bytes / 2,
                                   next_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    clense_float_impl<true, std::uint16_t>(reinterpret_cast<const std::uint16_t*>(srcp),
                                          reinterpret_cast<const std::uint16_t*>(prevp),
                                          reinterpret_cast<const std::uint16_t*>(nextp),
                                          reinterpret_cast<std::uint16_t*>(dstp),
                                          width, height,
                                          src_stride_bytes / 2, prev_stride_bytes / 2,
                                          next_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    clense_float_impl<false, float>(reinterpret_cast<const float*>(srcp),
                                    reinterpret_cast<const float*>(prevp),
                                    reinterpret_cast<const float*>(nextp),
                                    reinterpret_cast<float*>(dstp),
                                    width, height,
                                    src_stride_bytes / 4, prev_stride_bytes / 4,
                                    next_stride_bytes / 4, dst_stride_bytes / 4);
  }
}

void dispatch_clense_forward_backward_target(DataType dtype, const std::uint8_t* srcp,
                                            const std::uint8_t* ref1p, const std::uint8_t* ref2p,
                                            std::uint8_t* dstp, std::size_t width, std::size_t height,
                                            std::size_t src_stride_bytes, std::size_t ref1_stride_bytes,
                                            std::size_t ref2_stride_bytes, std::size_t dst_stride_bytes) {
  if (dtype == DataType::U8) {
    clense_forward_backward_int_impl<std::uint8_t>(srcp, ref1p, ref2p, dstp, width, height,
                                                   src_stride_bytes, ref1_stride_bytes, ref2_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    clense_forward_backward_int_impl<std::uint16_t>(reinterpret_cast<const std::uint16_t*>(srcp),
                                                    reinterpret_cast<const std::uint16_t*>(ref1p),
                                                    reinterpret_cast<const std::uint16_t*>(ref2p),
                                                    reinterpret_cast<std::uint16_t*>(dstp),
                                                    width, height,
                                                    src_stride_bytes / 2, ref1_stride_bytes / 2,
                                                    ref2_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    clense_forward_backward_float_impl<true, std::uint16_t>(reinterpret_cast<const std::uint16_t*>(srcp),
                                                           reinterpret_cast<const std::uint16_t*>(ref1p),
                                                           reinterpret_cast<const std::uint16_t*>(ref2p),
                                                           reinterpret_cast<std::uint16_t*>(dstp),
                                                           width, height,
                                                           src_stride_bytes / 2, ref1_stride_bytes / 2,
                                                           ref2_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    clense_forward_backward_float_impl<false, float>(reinterpret_cast<const float*>(srcp),
                                                     reinterpret_cast<const float*>(ref1p),
                                                     reinterpret_cast<const float*>(ref2p),
                                                     reinterpret_cast<float*>(dstp),
                                                     width, height,
                                                     src_stride_bytes / 4, ref1_stride_bytes / 4,
                                                     ref2_stride_bytes / 4, dst_stride_bytes / 4);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_clense_target);
HWY_EXPORT(dispatch_clense_forward_backward_target);

void process_clense_plane(DataType dtype, const std::uint8_t* srcp, const std::uint8_t* prevp,
                          const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                          std::size_t height, std::size_t src_stride_bytes,
                          std::size_t prev_stride_bytes, std::size_t next_stride_bytes,
                          std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_clense_target)(dtype, srcp, prevp, nextp, dstp, width, height,
                                               src_stride_bytes, prev_stride_bytes, next_stride_bytes,
                                               dst_stride_bytes);
}

void process_clense_forward_backward_plane(DataType dtype, const std::uint8_t* srcp,
                                          const std::uint8_t* ref1p, const std::uint8_t* ref2p,
                                          std::uint8_t* dstp, std::size_t width, std::size_t height,
                                          std::size_t src_stride_bytes, std::size_t ref1_stride_bytes,
                                          std::size_t ref2_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_clense_forward_backward_target)(dtype, srcp, ref1p, ref2p, dstp, width, height,
                                                                src_stride_bytes, ref1_stride_bytes,
                                                                ref2_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
