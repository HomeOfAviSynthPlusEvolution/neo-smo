#include "kernels/repair_dispatch.hpp"

namespace neo_smo {
void process_repair_plane(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp,
    const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width,
    std::size_t height, std::size_t src_stride_bytes,
    std::size_t repair_stride_bytes, std::size_t dst_stride_bytes) {
  switch (dtype) {
    case DataType::U8:
    case DataType::U16:
      process_repair_int_plane(dtype, mode, chroma, srcp, repairp, dstp, width, height, src_stride_bytes, repair_stride_bytes, dst_stride_bytes);
      break;
    case DataType::F16:
      process_repair_f16_plane(dtype, mode, chroma, srcp, repairp, dstp, width, height, src_stride_bytes, repair_stride_bytes, dst_stride_bytes);
      break;
    case DataType::F32:
      process_repair_f32_plane(dtype, mode, chroma, srcp, repairp, dstp, width, height, src_stride_bytes, repair_stride_bytes, dst_stride_bytes);
      break;
  }
}
} // namespace neo_smo
