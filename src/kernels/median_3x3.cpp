#include "kernels/median_dispatch.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/median_3x3.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "kernels/median-inl.hpp"

void dispatch_median_3x3_target(DataType dtype, const std::uint8_t* srcp, std::uint8_t* dstp, std::size_t width,
                            std::size_t height, std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);
  if (dtype == DataType::U8) {
    const auto* s = reinterpret_cast<const std::uint8_t*>(srcp);
    auto* d = reinterpret_cast<std::uint8_t*>(dstp);
    median_plane_impl<std::uint8_t, 1>(s, d, w, h, src_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    const auto* s = reinterpret_cast<const std::uint16_t*>(srcp);
    auto* d = reinterpret_cast<std::uint16_t*>(dstp);
    const std::size_t ss = src_stride_bytes / sizeof(std::uint16_t);
    const std::size_t ds = dst_stride_bytes / sizeof(std::uint16_t);
    median_plane_impl<std::uint16_t, 1>(s, d, w, h, ss, ds);
  } else if (dtype == DataType::F16) {
    const auto* s = reinterpret_cast<const std::uint16_t*>(srcp);
    auto* d = reinterpret_cast<std::uint16_t*>(dstp);
    const std::size_t ss = src_stride_bytes / sizeof(std::uint16_t);
    const std::size_t ds = dst_stride_bytes / sizeof(std::uint16_t);
    median_plane_f16_impl<1>(s, d, w, h, ss, ds);
  } else if (dtype == DataType::F32) {
    const auto* s = reinterpret_cast<const float*>(srcp);
    auto* d = reinterpret_cast<float*>(dstp);
    const std::size_t ss = src_stride_bytes / sizeof(float);
    const std::size_t ds = dst_stride_bytes / sizeof(float);
    median_plane_impl<float, 1>(s, d, w, h, ss, ds);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_median_3x3_target);
void process_median_3x3_plane(DataType dtype, const std::uint8_t* srcp, std::uint8_t* dstp, std::size_t width,
                          std::size_t height, std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_median_3x3_target)(dtype, srcp, dstp, width, height, src_stride_bytes,
                                               dst_stride_bytes);
}
} // namespace neo_smo
#endif
