#pragma once

#include "kernels/deen.hpp"
#include <array>
#include <cstddef>

namespace neo_smo {
// Highway-dispatched implementations; count is one or three.
// Inputs have radius samples of replicated padding on every side.
void deen_c_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, DeenThreshold spatial, DeenThreshold temporal, double* dst);
void deen_w_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, DeenThreshold spatial, DeenThreshold temporal, const double* weights, double* dst);
void deen_a_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, DeenThreshold spatial, DeenThreshold temporal, const double* weights, double* dst);
} // namespace neo_smo
