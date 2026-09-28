#include "kernels/deen_dispatch.hpp"
#include <algorithm>
#include <cmath>

namespace neo_smo {
void deen_w_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, DeenThreshold spatial, DeenThreshold temporal, const double* weights, double* dst) {
  const int side = 2 * radius + 1;
  double denominator = 0;
  for (int i = 0; i < side * side; ++i)
    denominator += weights[i];
  denominator *= count == 3 ? 4 : 1;
  for (int y = 0; y < height; ++y)
    for (int x = 0; x < width; ++x) {
      const auto center_index = (static_cast<std::size_t>(y) + radius) * pitch + x + radius;
      const double c = src[0][center_index];
      double residual = 0, low = c, high = c;
      for (int f = 0; f < count; ++f) {
        const auto threshold = f == 0 ? spatial : temporal;
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const double weight = weights[dy * side + dx] * (count == 3 && f == 0 ? 2 : 1);
            if (weight == 0)
              continue;
            const double s = src[f][(static_cast<std::size_t>(y) + dy) * pitch + x + dx];
            if (deen_detail::within_threshold(s, c, threshold)) {
              // Center-relative accumulation preserves constants without a rounding bias.
              residual += weight * (s - c);
              low = std::min(low, s);
              high = std::max(high, s);
            }
          }
      }
      dst[static_cast<std::size_t>(y) * width + x] = std::clamp(c + residual / denominator, low, high);
    }
}
} // namespace neo_smo
