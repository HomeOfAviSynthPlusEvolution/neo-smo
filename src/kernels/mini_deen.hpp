#pragma once
#include <cstddef>
#include <cstdint>
namespace neo_smo {
void mini_deen_kernel(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
                      int width, int height, bool byte_samples, int radius, unsigned threshold);
}
