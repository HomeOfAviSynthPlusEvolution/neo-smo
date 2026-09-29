#pragma once
#include "algorithms/deen.hpp"
namespace neo_smo {
std::uint64_t deen_integer_sad(const DeenPlane& a, const DeenPlane& b);
// Integer planes are read directly; each frame retains its own stride.
void deen_integer_kernel(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius,
                         double spatial, double temporal, const double* weights, const DeenByteWeights& byte_weights,
                         std::uint8_t* dst, std::ptrdiff_t stride);
void deen_float_kernel(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                       double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride);
// Each input is an edge-extended plane. Output contains visible samples only.
void deen_kernel(DeenFamily family, const std::array<const double*, 3>& src, int count, std::size_t pitch, int width,
                 int height, int radius, double spatial, double temporal, const double* weights, double* dst);
} // namespace neo_smo
