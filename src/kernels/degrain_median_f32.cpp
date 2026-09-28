#include "kernels/degrain_median_dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/degrain_median_f32.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
#include "kernels/degrain_median-inl.hpp"

void dispatch_degrain_median_f32_target(int mode, float limit, bool interlaced, bool norow, bool chroma,
                                        int bits_per_sample, const std::uint8_t* prevp, const std::uint8_t* currp,
                                        const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                                        std::size_t height, std::size_t prev_stride_bytes,
                                        std::size_t curr_stride_bytes, std::size_t next_stride_bytes,
                                        std::size_t dst_stride_bytes) {
  dispatch_degrain_median_impl<DataType::F32>(mode, limit, interlaced, norow, chroma, bits_per_sample, prevp, currp,
                                              nextp, dstp, width, height, prev_stride_bytes, curr_stride_bytes,
                                              next_stride_bytes, dst_stride_bytes);
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_degrain_median_f32_target);
void process_degrain_median_f32_plane(int mode, float limit, bool interlaced, bool norow, bool chroma,
                                      int bits_per_sample, const std::uint8_t* prevp, const std::uint8_t* currp,
                                      const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                                      std::size_t height, std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                      std::size_t next_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_degrain_median_f32_target)(mode, limit, interlaced, norow, chroma, bits_per_sample,
                                                           prevp, currp, nextp, dstp, width, height, prev_stride_bytes,
                                                           curr_stride_bytes, next_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
