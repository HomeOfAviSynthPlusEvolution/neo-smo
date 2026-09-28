#include "kernels/deen_dispatch.hpp"
#include <algorithm>
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/deen_w.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"
HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
#include "kernels/deen_simd-inl.hpp"
void deen_w_target(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, double spatial, double temporal, const double* weights, double* dst) {
  const hn::ScalableTag<double> d;
  const auto lanes = hn::Lanes(d);
  const int side = 2 * radius + 1;
  double denominator = 0;
  for (int i = 0; i < side * side; ++i)
    denominator += weights[i];
  denominator *= count == 3 ? 4 : 1;
  for (int y = 0; y < height; ++y)
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto active = std::min(lanes, static_cast<std::size_t>(width) - x);
      const auto center_index = (static_cast<std::size_t>(y) + radius) * pitch + x + radius;
      const auto c = hn::LoadN(d, src[0] + center_index, active);
      auto sum = hn::Zero(d);
      auto low = c, high = c;
      for (int f = 0; f < count; ++f) {
        const auto threshold = f == 0 ? spatial : temporal;
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const double weight = weights[dy * side + dx] * (count == 3 && f == 0 ? 2 : 1);
            if (weight == 0)
              continue;
            const auto s = hn::LoadN(d, src[f] + (static_cast<std::size_t>(y) + dy) * pitch + x + dx, active);
            const auto pass = deen_selection(d, s, c, threshold);
            sum = hn::MulAdd(hn::Set(d, weight), hn::IfThenElseZero(pass, hn::Sub(s, c)), sum);
            const auto value = hn::IfThenElse(pass, s, c);
            low = hn::Min(low, value);
            high = hn::Max(high, value);
          }
      }
      const auto result = hn::Clamp(hn::Add(c, hn::Div(sum, hn::Set(d, denominator))), low, high);
      hn::StoreN(result, d, dst + static_cast<std::size_t>(y) * width + x, active);
    }
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(deen_w_target);
void deen_w_kernel(const std::array<const double*, 3>& src, int count, std::size_t pitch, int width, int height,
                   int radius, double spatial, double temporal, const double* weights, double* dst) {
  HWY_DYNAMIC_DISPATCH(deen_w_target)(src, count, pitch, width, height, radius, spatial, temporal, weights, dst);
}
} // namespace neo_smo
#endif
