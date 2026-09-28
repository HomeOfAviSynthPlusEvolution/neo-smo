#pragma once
#include "algorithms/deen.hpp"
#include "kernels/deen_precision.hpp"
namespace neo_smo {
using DeenThreshold = deen_detail::Threshold;
// Each input is an edge-extended plane. Output contains visible samples only.
void deen_kernel(DeenFamily family, const std::array<const double*, 3>& src, int count, std::size_t pitch, int width,
                 int height, int radius, DeenThreshold spatial, DeenThreshold temporal, const double* weights,
                 double* dst);
} // namespace neo_smo
