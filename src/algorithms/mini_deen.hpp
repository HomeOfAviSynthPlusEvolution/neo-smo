#pragma once
#include <cstddef>
#include <cstdint>
namespace neo_smo {
// Integer samples only. Threshold is in 8-bit units; source and destination do not overlap.
void mini_deen_process(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
                       int width, int height, int bits, int radius, int threshold);
} // namespace neo_smo
