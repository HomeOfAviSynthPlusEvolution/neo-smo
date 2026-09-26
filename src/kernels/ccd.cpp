#include "kernels/dispatch.hpp"
#include "common/padded_row.hpp"
#include <array>
#include <type_traits>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/ccd.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {
#include "common/weighted_simd.hpp"

template <bool Half, bool Interior, class D, class T>
HWY_INLINE hn::Vec<D> ccd_load(D d, const T* plane, std::size_t stride,
                              int width, int height, int x, int y, std::size_t count) {
  if constexpr (Interior) {
    return weighted_load<Half>(d, plane + static_cast<std::size_t>(y) * stride + x, hn::Lanes(d));
  }
  const auto* row = plane + neo_smo::mirror_index(y, height) * stride;
  if (x >= 0 && static_cast<std::size_t>(x) + count <= static_cast<std::size_t>(width))
    return weighted_load<Half>(d, row + x, count);
  HWY_ALIGN T values[hn::MaxLanes(d)]{};
  for (std::size_t i = 0; i < count; ++i)
    values[i] = row[neo_smo::mirror_index(static_cast<std::int64_t>(x) + i, width)];
  return weighted_load<Half>(d, values, count);
}

template <bool Half, bool RGB, class D>
HWY_INLINE hn::Vec<D> ccd_distance(D d, hn::Vec<D> a, hn::Vec<D> b, hn::Vec<D> c,
                                  hn::Vec<D> ca, hn::Vec<D> cb, hn::Vec<D> cc) {
  const auto da = float_sub<Half>(d, a, ca);
  const auto db = float_sub<Half>(d, b, cb);
  const auto dc = float_sub<Half>(d, c, cc);
  auto aa = float_mul<Half>(d, da, da);
  if constexpr (!RGB) aa = float_mul<Half>(d, aa, hn::Set(d, 4));
  return float_add<Half>(d, float_add<Half>(d, aa, float_mul<Half>(d, db, db)), float_mul<Half>(d, dc, dc));
}

template <class T, bool Half, bool RGB>
HWY_NOINLINE void ccd_impl(int width, int height, std::size_t stride_bytes,
    const std::uint8_t* const* src, const std::uint8_t* const* ref,
    std::uint8_t* dst_r, std::uint8_t* dst_g, std::uint8_t* dst_b,
    float threshold, int radius, const float* weights, const Point* points, int num_points, int bits) {
  constexpr bool integer = !Half && !std::is_same_v<T, float>;
  // U8 squared-distance accumulation is bounded by 6 * 255^2 * 21 < 2^23,
  // so FP32 represents each integer exactly. U16 still requires binary64.
  constexpr bool wide_integer = integer && sizeof(T) > 1;
  const hn::ScalableTag<std::conditional_t<wide_integer, double, FloatLane<Half>>> d;
  const hn::Rebind<float, decltype(d)> df;
  const std::size_t stride = stride_bytes / sizeof(T), lanes = hn::Lanes(d);
  // floor((ssd + radius) / diameter) < floor(threshold) is exactly
  // ssd < floor(threshold) * diameter - radius for integer SSD values.
  const double cutoff = integer ? std::floor(static_cast<double>(threshold)) * (radius * 2 + 1) - radius : threshold;
  const auto threshold_v = hn::Set(d, cutoff);
  const auto div = hn::Set(d, radius * 2 + 1);
  T* destinations[3] = {reinterpret_cast<T*>(dst_r), reinterpret_cast<T*>(dst_g), reinterpret_cast<T*>(dst_b)};
  int min_x = 0, max_x = 0, min_y = 0, max_y = 0;
  for (int p = 0; p < num_points; ++p) {
    min_x = std::min(min_x, points[p].x); max_x = std::max(max_x, points[p].x);
    min_y = std::min(min_y, points[p].y); max_y = std::max(max_y, points[p].y);
  }
  for (int y = 0; y < height; ++y) {
    auto block = [&](int x, auto interior_tag) HWY_ATTR {
      constexpr bool Interior = decltype(interior_tag)::value;
      const auto count = std::min(lanes, static_cast<std::size_t>(width - x));
      auto load_ref = [&](int frame, int plane, int px, int py) HWY_ATTR {
        return ccd_load<Half, Interior>(d, reinterpret_cast<const T*>(ref[frame * 3 + plane]), stride, width, height, px, py, count);
      };
      auto load_src = [&](int plane, int px, int py) HWY_ATTR {
        return ccd_load<Half, Interior>(d, reinterpret_cast<const T*>(src[radius * 3 + plane]), stride, width, height, px, py, count);
      };
      const auto ca = load_ref(radius, 0, x, y), cb = load_ref(radius, 1, x, y), cc = load_ref(radius, 2, x, y);
      auto ta = hn::Zero(d), tb = load_src(1, x, y), tc = load_src(2, x, y);
      if constexpr (RGB) ta = load_src(0, x, y);
      auto accepted = hn::Set(d, 1);
      for (int p = 0; p < num_points; ++p) {
        const int px = x + points[p].x, py = y + points[p].y;
        auto distance = [&](int frame) HWY_ATTR {
          return ccd_distance<Half, RGB>(d, load_ref(frame, 0, px, py), load_ref(frame, 1, px, py), load_ref(frame, 2, px, py), ca, cb, cc);
        };
        auto ssd = distance(radius);
        for (int i = 0; i < radius; ++i) {
          const int prev = radius - 1 - i, next = radius + 1 + i;
          const auto a = distance(prev), b = distance(next);
          if constexpr (integer) {
            if constexpr (wide_integer) {
              const auto term = hn::Add(hn::Mul(hn::DemoteTo(df, a), hn::Set(df, weights[prev])),
                                        hn::Mul(hn::DemoteTo(df, b), hn::Set(df, weights[next])));
              ssd = hn::Add(ssd, hn::PromoteTo(d, weighted_round(df, term)));
            } else {
              const auto term = hn::Add(hn::Mul(a, hn::Set(d, weights[prev])),
                                        hn::Mul(b, hn::Set(d, weights[next])));
              ssd = hn::Add(ssd, weighted_round(d, term));
            }
          } else {
            auto wp = hn::Set(d, weights[prev]), wn = hn::Set(d, weights[next]);
            ssd = float_add<Half>(d, ssd, float_add<Half>(d, float_mul<Half>(d, a, wp), float_mul<Half>(d, b, wn)));
          }
        }
        if constexpr (!integer) {
          if (radius) ssd = float_div<Half>(d, ssd, div);
        }
        const auto mask = hn::Lt(ssd, threshold_v);
        if constexpr (RGB) ta = hn::IfThenElse(mask, float_add<Half>(d, ta, load_src(0, px, py)), ta);
        tb = hn::IfThenElse(mask, float_add<Half>(d, tb, load_src(1, px, py)), tb);
        tc = hn::IfThenElse(mask, float_add<Half>(d, tc, load_src(2, px, py)), tc);
        accepted = hn::Add(accepted, hn::IfThenElse(mask, hn::Set(d, 1), hn::Zero(d)));
      }
      auto store = [&](hn::Vec<decltype(d)> total, int plane) HWY_ATTR {
        if constexpr (!Half && !std::is_same_v<T, float>) {
          // Round the rational integer average without an FP32 reciprocal
          // moving exact half ties down by one code value.
          const auto numerator = hn::Add(hn::Add(total, total), accepted);
          const auto result = hn::Floor(hn::Div(numerator, hn::Add(accepted, accepted)));
          weighted_store<false>(d, hn::Min(result, hn::Set(d, (1u << bits) - 1)), destinations[plane] + y * stride + x, count);
        } else {
          weighted_store<Half>(d, float_div<Half>(d, total, accepted), destinations[plane] + y * stride + x, count);
        }
      };
      if constexpr (RGB) store(ta, 0);
      store(tb, 1); store(tc, 2);
    };
    for (int x = 0; x < width; x += static_cast<int>(lanes)) {
      const bool interior = static_cast<int64_t>(x) + min_x >= 0 &&
          static_cast<int64_t>(x) + lanes + max_x <= static_cast<std::size_t>(width) &&
          static_cast<int64_t>(y) + min_y >= 0 && static_cast<int64_t>(y) + max_y < height;
      if (interior) block(x, std::true_type{});
      else block(x, std::false_type{});
    }
  }
}

void dispatch_ccd_target(DataType dtype, bool is_rgb, int width, int height, std::size_t stride_bytes,
    const std::uint8_t* const* src, const std::uint8_t* const* ref,
    std::uint8_t* dst_r, std::uint8_t* dst_g, std::uint8_t* dst_b,
    float threshold, int temporal_radius, const float* weights, const Point* points,
    int num_points, int, float, int bits_per_sample) {
#define CCD_RUN(T, HALF) \
  if (is_rgb) ccd_impl<T, HALF, true>(width, height, stride_bytes, src, ref, dst_r, dst_g, dst_b, threshold, temporal_radius, weights, points, num_points, bits_per_sample); \
  else ccd_impl<T, HALF, false>(width, height, stride_bytes, src, ref, dst_r, dst_g, dst_b, threshold, temporal_radius, weights, points, num_points, bits_per_sample)
  switch (dtype) {
    case DataType::U8: { CCD_RUN(std::uint8_t, false); break; }
    case DataType::U16: { CCD_RUN(std::uint16_t, false); break; }
    case DataType::F16: { CCD_RUN(std::uint16_t, true); break; }
    case DataType::F32: { CCD_RUN(float, false); break; }
  }
#undef CCD_RUN
}
} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {

HWY_EXPORT(dispatch_ccd_target);

void process_ccd_planes(
    DataType dtype,
    bool is_rgb,
    int width,
    int height,
    std::size_t stride_bytes,
    const std::uint8_t* const* src,
    const std::uint8_t* const* ref,
    std::uint8_t* dst_r,
    std::uint8_t* dst_g,
    std::uint8_t* dst_b,
    float threshold,
    int temporal_radius,
    const float* weights,
    const Point* points,
    int num_points,
    int diameter,
    float scale,
    int bits_per_sample
) {
  HWY_DYNAMIC_DISPATCH(dispatch_ccd_target)(
      dtype, is_rgb, width, height, stride_bytes,
      src, ref, dst_r, dst_g, dst_b,
      threshold, temporal_radius, weights,
      points, num_points, diameter, scale, bits_per_sample);
}

} // namespace neo_smo
#endif
