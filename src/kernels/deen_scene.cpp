#include "kernels/deen_scalar.hpp"
#include "kernels/deen.hpp"
#include <algorithm>
#include <cstring>
#include <stdexcept>
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/deen_scene.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"
HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
#if HWY_TARGET == HWY_SCALAR || HWY_TARGET == HWY_EMU128
std::uint64_t deen_sad_target(const DeenPlane& a, const DeenPlane& b) {
  return deen_sad_scalar(a, b);
}
#else

namespace hn = hwy::HWY_NAMESPACE;
std::uint64_t deen_sad_target(const DeenPlane& a, const DeenPlane& b) {
  if (a.type == DataType::U8) {
    const hn::ScalableTag<std::uint8_t> d;
    const hn::Repartition<std::uint64_t, decltype(d)> ds;
    auto sum = hn::Zero(ds);
    for (int y = 0; y < a.height; ++y)
      for (int x = 0; x < a.width;) {
        const auto active = std::min(hn::Lanes(d), static_cast<std::size_t>(a.width - x));
        const auto av = hn::LoadN(d, a.data + y * a.stride + x, active);
        const auto bv = hn::LoadN(d, b.data + y * b.stride + x, active);
        sum = hn::Add(sum, hn::SumsOf8(hn::Sub(hn::Max(av, bv), hn::Min(av, bv))));
        x += static_cast<int>(active);
      }
    return hn::ReduceSum(ds, sum);
  }
  const hn::ScalableTag<std::uint32_t> d;
  const hn::Rebind<std::uint16_t, decltype(d)> dn;
  const hn::Repartition<std::uint8_t, decltype(dn)> db;
  const auto peak = hn::Set(d, (1u << a.bits) - 1);
  std::uint64_t total = 0;
  auto maximum = hn::Zero(d);
  for (int y = 0; y < a.height; ++y)
    for (int x = 0; x < a.width;) {
      // Reduce in bounded chunks: even the horizontal u32 sum cannot overflow.
      const int end = x + std::min(4096, a.width - x);
      auto sum = hn::Zero(d);
      while (x < end) {
        const auto active = std::min(hn::Lanes(d), static_cast<std::size_t>(end - x));
        const auto load = [&](const DeenPlane& p) HWY_ATTR {
          const auto* ptr = p.data + y * p.stride + static_cast<std::size_t>(x) * 2;
          return hn::PromoteTo(d, hn::BitCast(dn, hn::LoadN(db, ptr, active * 2)));
        };
        const auto av = load(a), bv = load(b);
        maximum = hn::Max(maximum, hn::Max(av, bv));
        sum = hn::Add(sum, hn::Sub(hn::Max(av, bv), hn::Min(av, bv)));
        x += static_cast<int>(active);
      }
      total += hn::ReduceSum(d, sum);
    }
  if (!hn::AllFalse(d, hn::Gt(maximum, peak)))
    throw std::invalid_argument("Deen: sample exceeds bit depth.");
  return total;
}
#endif
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();
#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(deen_sad_target);
std::uint64_t deen_integer_sad(const DeenPlane& a, const DeenPlane& b) {
  return HWY_DYNAMIC_DISPATCH(deen_sad_target)(a, b);
}
} // namespace neo_smo
#endif
