#include "kernels/degrain_median_dispatch.hpp"

namespace neo_smo {
void process_degrain_median_plane(DataType dtype, int mode, float limit, bool interlaced, bool norow, bool chroma,
                                  int bits_per_sample, const std::uint8_t* prevp, const std::uint8_t* currp,
                                  const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width, std::size_t height,
                                  std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                  std::size_t next_stride_bytes, std::size_t dst_stride_bytes) {
  switch (dtype) {
    case DataType::U8:
      process_degrain_median_u8_plane(mode, limit, interlaced, norow, chroma, bits_per_sample, prevp, currp, nextp,
                                      dstp, width, height, prev_stride_bytes, curr_stride_bytes, next_stride_bytes,
                                      dst_stride_bytes);
      break;
    case DataType::U16:
      process_degrain_median_u16_plane(mode, limit, interlaced, norow, chroma, bits_per_sample, prevp, currp, nextp,
                                       dstp, width, height, prev_stride_bytes, curr_stride_bytes, next_stride_bytes,
                                       dst_stride_bytes);
      break;
    case DataType::F16:
      process_degrain_median_f16_plane(mode, limit, interlaced, norow, chroma, bits_per_sample, prevp, currp, nextp,
                                       dstp, width, height, prev_stride_bytes, curr_stride_bytes, next_stride_bytes,
                                       dst_stride_bytes);
      break;
    case DataType::F32:
      process_degrain_median_f32_plane(mode, limit, interlaced, norow, chroma, bits_per_sample, prevp, currp, nextp,
                                       dstp, width, height, prev_stride_bytes, curr_stride_bytes, next_stride_bytes,
                                       dst_stride_bytes);
      break;
  }
}
} // namespace neo_smo
