#pragma once
#include "kernels/dispatch.hpp"
namespace neo_smo {
void process_median_3x3_plane(DataType dtype, const std::uint8_t* srcp, std::uint8_t* dstp, std::size_t width,
                          std::size_t height, std::size_t src_stride_bytes, std::size_t dst_stride_bytes);
void process_median_5x5_plane(DataType dtype, const std::uint8_t* srcp, std::uint8_t* dstp, std::size_t width,
                          std::size_t height, std::size_t src_stride_bytes, std::size_t dst_stride_bytes);
void process_median_7x7_plane(DataType dtype, const std::uint8_t* srcp, std::uint8_t* dstp, std::size_t width,
                          std::size_t height, std::size_t src_stride_bytes, std::size_t dst_stride_bytes);
} // namespace neo_smo
