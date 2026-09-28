#pragma once
#include "kernels/dispatch.hpp"

namespace neo_smo {
void process_degrain_median_u8_plane(int mode, float limit, bool interlaced, bool norow, bool chroma,
                                     int bits_per_sample, const std::uint8_t* prevp, const std::uint8_t* currp,
                                     const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                                     std::size_t height, std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                     std::size_t next_stride_bytes, std::size_t dst_stride_bytes);
void process_degrain_median_u16_plane(int mode, float limit, bool interlaced, bool norow, bool chroma,
                                      int bits_per_sample, const std::uint8_t* prevp, const std::uint8_t* currp,
                                      const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                                      std::size_t height, std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                      std::size_t next_stride_bytes, std::size_t dst_stride_bytes);
void process_degrain_median_f16_plane(int mode, float limit, bool interlaced, bool norow, bool chroma,
                                      int bits_per_sample, const std::uint8_t* prevp, const std::uint8_t* currp,
                                      const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                                      std::size_t height, std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                      std::size_t next_stride_bytes, std::size_t dst_stride_bytes);
void process_degrain_median_f32_plane(int mode, float limit, bool interlaced, bool norow, bool chroma,
                                      int bits_per_sample, const std::uint8_t* prevp, const std::uint8_t* currp,
                                      const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width,
                                      std::size_t height, std::size_t prev_stride_bytes, std::size_t curr_stride_bytes,
                                      std::size_t next_stride_bytes, std::size_t dst_stride_bytes);
} // namespace neo_smo
