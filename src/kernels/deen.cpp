#include "kernels/deen.hpp"
#include "kernels/deen_dispatch.hpp"
#include "base/checked.hpp"

namespace neo_smo {
void deen_kernel(DeenFamily family, const std::array<const double*, 3>& src, int count, std::size_t pitch, int width,
                 int height, int radius, DeenThreshold spatial, DeenThreshold temporal, const double* weights,
                 double* dst) {
  switch (family) {
    case DeenFamily::Constant:
      deen_c_kernel(src, count, pitch, width, height, radius, spatial, temporal, dst);
      return;
    case DeenFamily::Weighted:
      deen_w_kernel(src, count, pitch, width, height, radius, spatial, temporal, weights, dst);
      return;
    case DeenFamily::Adaptive:
      deen_a_kernel(src, count, pitch, width, height, radius, spatial, temporal, weights, dst);
      return;
  }
  require(false, "Deen: invalid kernel family.");
}
} // namespace neo_smo
