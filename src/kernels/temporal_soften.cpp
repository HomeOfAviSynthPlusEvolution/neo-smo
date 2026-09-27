#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/temporal_soften.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/float_arithmetic.hpp"

template <typename T>
void temporal_soften_int_impl(int diameter, T threshold, const T* const* srcp_planes, T* dstp,
                              int width, int height, std::size_t src_stride, std::size_t dst_stride) {
  if (diameter == 1) {
    copy_plane(dstp, srcp_planes[0], width, height, dst_stride, src_stride);
    return;
  }
  using ComputeT = std::conditional_t<sizeof(T) == 1, std::int16_t, std::int32_t>;
  const hn::ScalableTag<ComputeT> d;
  const hn::Rebind<T, decltype(d)> ds;
  const std::size_t lanes = hn::Lanes(d);
  const auto thresh = hn::Set(d, threshold);
  const auto rounding = hn::Set(d, diameter / 2);
  constexpr unsigned scale_bits = sizeof(T) == 1 ? 16 : 32;
  const hn::RebindToUnsigned<decltype(d)> du;
  const auto multiplier = hn::Set(du, (std::uint64_t{1} << scale_bits) / diameter);
  for (int y = 0; y < height; ++y) {
    const auto offset = static_cast<std::size_t>(y) * src_stride;
    auto* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
      const auto curr = hn::PromoteTo(d, hn::LoadN(ds, srcp_planes[0] + offset + x, count));
      auto sum = curr;
      for (int i = 1; i < diameter; ++i) {
        const auto value = hn::PromoteTo(d, hn::LoadN(ds, srcp_planes[i] + offset + x, count));
        sum = hn::Add(sum, hn::IfThenElse(hn::Le(hn::AbsDiff(curr, value), thresh), value, curr));
      }
      const auto value = hn::Add(sum, rounding);
      // Preserve the upstream fixed-point reciprocal, including its truncation.
      const auto result = hn::BitCast(d, hn::MulHigh(hn::BitCast(du, value), multiplier));
      hn::StoreN(hn::DemoteTo(ds, result), ds, dst_row + x, count);
    }
  }
}

template <bool IsF16, typename StorageT, int FixedDiameter = 0>
void temporal_soften_float_impl(int runtime_diameter, float threshold, const StorageT* const* srcp_planes,
                                StorageT* dstp, int width, int height, std::size_t src_stride,
                                std::size_t dst_stride) {
  const int diameter = FixedDiameter ? FixedDiameter : runtime_diameter;
  using ComputeT = FloatLane<IsF16>;
  const hn::ScalableTag<float> df;
  const hn::Rebind<ComputeT, decltype(df)> d;
  const std::size_t lanes = hn::Lanes(d);
  const auto thresh_vec = hn::Set(d, threshold);
  const auto reciprocal = hn::Set(df, 1.0f / static_cast<float>(diameter));

  if constexpr (IsF16) {
    for (int y = 0; y < height; ++y) {
      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
        const auto curr = load_f16(d, srcp_planes[0] + static_cast<std::size_t>(y) * src_stride + x, count);
        auto sum = curr;

        for (int i = 1; i < diameter; ++i) {
          const auto f_val = load_f16(
              d, srcp_planes[i] + static_cast<std::size_t>(y) * src_stride + x, count);
          auto chosen = curr;
#if HWY_HAVE_FLOAT16
          if constexpr (std::is_same_v<ComputeT, hwy::float16_t>) {
            const auto diff = hn::Abs(hn::Sub(hn::PromoteTo(df, curr), hn::PromoteTo(df, f_val)));
            // Compare in F32, but select the original half lanes without a conversion back.
            chosen = hn::IfThenElse(hn::DemoteMaskTo(d, df, hn::Le(diff, hn::Set(df, threshold))), f_val, curr);
          } else
#endif
          {
            const auto diff = hn::Abs(hn::Sub(curr, f_val));
            chosen = hn::IfThenElse(hn::Le(diff, thresh_vec), f_val, curr);
          }
          sum = float_add<true>(d, sum, chosen);
        }

#if HWY_HAVE_FLOAT16
        if constexpr (std::is_same_v<ComputeT, hwy::float16_t>) {
          const auto res = hn::Mul(hn::PromoteTo(df, sum), reciprocal);
          store_f16(df, res, dst_row + x, count);
        } else
#endif
        {
          store_f16(d, hn::Mul(sum, reciprocal), dst_row + x, count);
        }
      }
    }
  } else {
    for (int y = 0; y < height; ++y) {
      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const std::size_t rem = static_cast<std::size_t>(width - x);
        const std::size_t count = std::min(lanes, rem);

        const auto curr = hn::LoadN(d, srcp_planes[0] + static_cast<std::size_t>(y) * src_stride + x, count);
        auto sum = curr;

        for (int i = 1; i < diameter; ++i) {
          const auto f_val = hn::LoadN(d, srcp_planes[i] + static_cast<std::size_t>(y) * src_stride + x, count);
          const auto diff = hn::Abs(hn::Sub(curr, f_val));
          const auto chosen = hn::IfThenElse(hn::Le(diff, thresh_vec), f_val, curr);
          sum = hn::Add(sum, chosen);
        }

        const auto res = hn::Mul(sum, reciprocal);
        if (rem >= lanes) {
          hn::StoreU(res, d, dst_row + x);
        } else {
          hn::StoreN(res, d, dst_row + x, rem);
        }
      }
    }
  }
}

template <int FixedDiameter>
void temporal_soften_dispatch_impl(DataType dtype, int diameter, float threshold,
                                     const std::uint8_t* const* srcp_planes, std::uint8_t* dstp,
                                     std::size_t width, std::size_t height,
                                     std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  if (dtype == DataType::F16) {
    if constexpr (HWY_HAVE_FLOAT16) threshold = fp16_to_fp32(fp32_to_fp16(threshold));
  }
  std::array<const std::uint16_t*, 21> u16_planes{};
  std::array<const float*, 21> f32_planes{};
  for (int i = 0; i < diameter; ++i) {
    u16_planes[i] = reinterpret_cast<const std::uint16_t*>(srcp_planes[i]);
    f32_planes[i] = reinterpret_cast<const float*>(srcp_planes[i]);
  }
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);

  if (dtype == DataType::U8) {
    const auto th = static_cast<std::uint8_t>(std::clamp(threshold, 0.0f, 255.0f));
    temporal_soften_int_impl<std::uint8_t>(diameter, th, srcp_planes, dstp, w, h, src_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    const auto th = static_cast<std::uint16_t>(std::clamp(threshold, 0.0f, 65535.0f));
    temporal_soften_int_impl<std::uint16_t>(
        diameter, th, u16_planes.data(),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    temporal_soften_float_impl<true, std::uint16_t, FixedDiameter>(
        diameter, threshold, u16_planes.data(),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    temporal_soften_float_impl<false, float, FixedDiameter>(
        diameter, threshold, f32_planes.data(),
        reinterpret_cast<float*>(dstp), w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
  }
}

void dispatch_temporal_soften_target(DataType dtype, int diameter, float threshold,
    const std::uint8_t* const* srcp_planes, std::uint8_t* dstp,
    std::size_t width, std::size_t height, std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  if (diameter == 21 && (dtype == DataType::F16 || dtype == DataType::F32)) {
    temporal_soften_dispatch_impl<21>(dtype, diameter, threshold, srcp_planes, dstp,
        width, height, src_stride_bytes, dst_stride_bytes);
  } else {
    temporal_soften_dispatch_impl<0>(dtype, diameter, threshold, srcp_planes, dstp,
        width, height, src_stride_bytes, dst_stride_bytes);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_temporal_soften_target);

void process_temporal_soften_plane(DataType dtype, int diameter, float threshold,
                                   const std::uint8_t* const* srcp_planes, std::uint8_t* dstp,
                                   std::size_t width, std::size_t height,
                                   std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_temporal_soften_target)(
      dtype, diameter, threshold, srcp_planes, dstp, width, height, src_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
