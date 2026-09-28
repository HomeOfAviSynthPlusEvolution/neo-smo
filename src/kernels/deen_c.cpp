#include "kernels/deen_dispatch.hpp"
#include <cmath>

namespace neo_smo {
void deen_c_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, DeenThreshold spatial, DeenThreshold temporal, double* dst) {
  const int side = 2 * radius + 1;
  for (int y = 0; y < height; ++y)
    for (int x = 0; x < width; ++x) {
      const auto center_index = (static_cast<std::size_t>(y) + radius) * pitch + x + radius;
      const double c = src[0][center_index];
      double sum = 0;
      for (int f = 0; f < count; ++f) {
        const auto threshold = f == 0 ? spatial : temporal;
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const double s = src[f][(static_cast<std::size_t>(y) + dy) * pitch + x + dx];
            sum += deen_detail::within_threshold(s, c, threshold) ? s : c;
          }
      }
      dst[static_cast<std::size_t>(y) * width + x] = sum / (count * side * side);
    }
}
} // namespace neo_smo
