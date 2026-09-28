#pragma once
#include "algorithms/deen.hpp"
namespace neo_smo {
// Each input is an edge-extended plane. Output contains visible samples only.
void deen_kernel(DeenFamily family, const std::array<const double*, 3>& src, int count, std::size_t pitch, int width,
                 int height, int radius, double spatial, double temporal, const double* weights, double* dst);
} // namespace neo_smo
