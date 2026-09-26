#pragma once
#include "kernels/dispatch.hpp"

namespace neo_smo {
void process_repair_int_plane(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp,
    const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width,
    std::size_t height, std::size_t src_stride_bytes,
    std::size_t repair_stride_bytes, std::size_t dst_stride_bytes);
void process_repair_f16_plane(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp,
    const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width,
    std::size_t height, std::size_t src_stride_bytes,
    std::size_t repair_stride_bytes, std::size_t dst_stride_bytes);
void process_repair_f32_plane(DataType dtype, int mode, bool chroma, const std::uint8_t* srcp,
    const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width,
    std::size_t height, std::size_t src_stride_bytes,
    std::size_t repair_stride_bytes, std::size_t dst_stride_bytes);
} // namespace neo_smo
