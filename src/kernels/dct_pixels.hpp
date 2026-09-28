#pragma once

#include "base/types.hpp"
#include <cstddef>

namespace neo_smo {
void dct_load_row(DataType type, const std::uint8_t* src, float* dst, std::size_t width);
void dct_store_row(DataType type, const float* src, std::uint8_t* dst, std::size_t width, int bits);
} // namespace neo_smo
