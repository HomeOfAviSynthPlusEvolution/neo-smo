#include "kernels/repair_dispatch.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/repair_f16.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "kernels/repair_float-inl.hpp"

void dispatch_repair_f16_target(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp,
    const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width,
    std::size_t height, std::size_t src_stride_bytes,
    std::size_t repair_stride_bytes, std::size_t dst_stride_bytes) {
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);
  (void)dtype;
  repair_float_impl<true, std::uint16_t>(mode, chroma,
      reinterpret_cast<const std::uint16_t*>(srcp), reinterpret_cast<const std::uint16_t*>(repairp),
      reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / sizeof(std::uint16_t),
      repair_stride_bytes / sizeof(std::uint16_t), dst_stride_bytes / sizeof(std::uint16_t));
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_repair_f16_target);
void process_repair_f16_plane(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp,
    const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width,
    std::size_t height, std::size_t src_stride_bytes,
    std::size_t repair_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_repair_f16_target)(dtype, mode, chroma, srcp, repairp, dstp, width, height, src_stride_bytes, repair_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
