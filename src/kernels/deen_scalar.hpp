#pragma once
#include "kernels/deen.hpp"

namespace neo_smo {
// Plain C++ fallbacks; no Highway vector types or operations.
void deen_integer_scalar(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius,
                         double spatial, double temporal, const double* weights, const DeenByteWeights& table,
                         std::uint8_t* dst, std::ptrdiff_t stride);
void deen_float_scalar(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                       double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride);
std::uint64_t deen_sad_scalar(const DeenPlane& a, const DeenPlane& b);
void deen_padded_scalar(DeenFamily family, const std::array<const double*, 3>& src, int count, std::size_t pitch,
                        int width, int height, int radius, double spatial, double temporal, const double* weights,
                        double* dst);
} // namespace neo_smo
