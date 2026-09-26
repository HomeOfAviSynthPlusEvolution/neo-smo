#include "kernels/median_dispatch.hpp"

namespace neo_smo {
void process_median_plane(DataType dtype, int radius, const std::uint8_t* srcp, std::uint8_t* dstp, std::size_t width,
                          std::size_t height, std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  switch (radius) {
    case 1: process_median_3x3_plane(dtype, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes); break;
    case 2: process_median_5x5_plane(dtype, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes); break;
    case 3: process_median_7x7_plane(dtype, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes); break;
  }
}
} // namespace neo_smo
