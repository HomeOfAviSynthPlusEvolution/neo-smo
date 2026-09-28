#include "kernels/deen_dispatch.hpp"
#include "kernels/deen_precision.hpp"
#include <limits>
#include <algorithm>
#include <cmath>

namespace neo_smo {
void deen_a_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, DeenThreshold spatial, DeenThreshold temporal, const double* weights, double* dst) {
  const int side = 2 * radius + 1;
  for (int y = 0; y < height; ++y)
    for (int x = 0; x < width; ++x) {
      const auto center_index = (static_cast<std::size_t>(y) + radius) * pitch + x + radius;
      const double c = src[0][center_index];
      double sum = 0, low = c, high = c;
      int accepted = 0;
      for (int f = 0; f < count; ++f) {
        const auto threshold = f == 0 ? spatial : temporal;
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const double s = src[f][(static_cast<std::size_t>(y) + dy) * pitch + x + dx];
            const double limit = threshold.hi * weights[dy * side + dx];
            const double difference = std::abs(s - c);
            bool include = difference <= limit;
            const double uncertainty = 16 * std::numeric_limits<double>::epsilon() * std::max(threshold.hi, difference);
            if (std::abs(difference - limit) <= uncertainty)
              include = deen_detail::adaptive_boundary(s, c, threshold, weights[0], dx - radius, dy - radius, radius);
            if (include) {
              sum += s;
              ++accepted;
              low = std::min(low, s);
              high = std::max(high, s);
            }
          }
      }
      // Current-frame offset (0,0) always passes; no extra center seeding.
      dst[static_cast<std::size_t>(y) * width + x] = std::clamp(sum / accepted, low, high);
    }
}
} // namespace neo_smo
