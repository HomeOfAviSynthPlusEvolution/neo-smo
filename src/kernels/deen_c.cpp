#include "kernels/deen_dispatch.hpp"
#include <algorithm>
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/deen_c.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"
HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
#include "kernels/deen_simd-inl.hpp"
void deen_c_target(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, double spatial, double temporal, double* dst) {
  const hn::ScalableTag<double> d;
  const auto lanes = hn::Lanes(d);
  const int side = 2 * radius + 1;
  for (int y = 0; y < height; ++y)
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto active = std::min(lanes, static_cast<std::size_t>(width) - x);
      const auto center_index = (static_cast<std::size_t>(y) + radius) * pitch + x + radius;
      const auto c = hn::LoadN(d, src[0] + center_index, active);
      auto sum = hn::Zero(d);
      for (int f = 0; f < count; ++f) {
        const auto threshold = f == 0 ? spatial : temporal;
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const auto s = hn::LoadN(d, src[f] + (static_cast<std::size_t>(y) + dy) * pitch + x + dx, active);
            const auto pass = deen_selection(d, s, c, threshold);
            sum = hn::Add(sum, hn::IfThenElse(pass, s, c));
          }
      }
      const auto result = hn::Div(sum, hn::Set(d, count * side * side));
      hn::StoreN(result, d, dst + static_cast<std::size_t>(y) * width + x, active);
    }
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(deen_c_target);
void deen_c_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, double spatial, double temporal, double* dst) {
  HWY_DYNAMIC_DISPATCH(deen_c_target)(src, count, pitch, width, height, radius, spatial, temporal, dst);
}
} // namespace neo_smo
#endif
