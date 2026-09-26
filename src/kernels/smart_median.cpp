#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/smart_median.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/sorting_networks.hpp"
#include "common/u16_sort.hpp"
#include "common/float_arithmetic.hpp"
#include "common/fp16_rows.hpp"

template <int Radius, class D, class V = hn::Vec<D>>
HWY_INLINE V eval_smart_median_int32(D d, V center, V* values, V threshold) {
  constexpr int kSide = 2 * Radius + 1;
  constexpr int kCount = kSide * kSide;
  constexpr int kEvenCount = kCount - 1;
  constexpr int kCenterIdx = kEvenCount / 2;

  if constexpr (Radius == 1) {
    median8_even(d, values);
  } else if constexpr (Radius == 2) {
    median24_even(d, values);
  } else {
#if defined(NEO_SMO_MSVC_U16_SORT) && (HWY_TARGET == HWY_SSE2 || HWY_TARGET == HWY_SSSE3)
    sort_u16_samples<48>(d, values, [](auto ds, auto* v) { median48_even(ds, v); });
#else
    median48_even(d, values);
#endif
  }

  const V median_left = values[kCenterIdx - 1];
  const V median_right = values[kCenterIdx];

  V sum = values[0];
  for (int i = 1; i < kEvenCount; ++i) {
    sum = hn::Add(sum, values[i]);
  }

  // After removing the power of two, MulHigh divides the bounded sum by 3 exactly.
  V average;
  if constexpr (Radius == 1) {
    average = hn::ShiftRight<3>(hn::Add(sum, hn::Set(d, 4)));
  } else if constexpr (Radius == 2) {
    average = hn::MulHigh(hn::ShiftRight<3>(hn::Add(sum, hn::Set(d, 12))), hn::Set(d, 0x55555556));
  } else {
    average = hn::MulHigh(hn::ShiftRight<4>(hn::Add(sum, hn::Set(d, 24))), hn::Set(d, 0x55555556));
  }

  const V diff_l = hn::Sub(median_left, average);
  const V diff_r = hn::Sub(median_right, average);
  const V sq_l = hn::Mul(diff_l, diff_l);
  const V sq_r = hn::Mul(diff_r, diff_r);
  const V sq_sum = hn::Add(sq_l, sq_r);

  const hn::RebindToFloat<D> df;
  const auto f_sum = hn::ConvertTo(df, sq_sum);
  const auto f_var = hn::Mul(hn::Sqrt(f_sum), hn::Set(df, 13.0f));
  // Zig rounds nonnegative half-integers upward, not to the nearest even integer.
  const auto curved_variance = hn::ConvertTo(d, hn::Floor(hn::Add(f_var, hn::Set(df, 0.5f))));

  const auto lte = hn::Le(curved_variance, threshold);
  const auto med_lo = hn::Min(median_left, median_right);
  const auto med_hi = hn::Max(median_left, median_right);
  const auto clamped = hn::Clamp(center, med_lo, med_hi);

  return hn::IfThenElse(lte, clamped, center);
}

template <bool IsF16, int Radius, class D, class V = hn::Vec<D>>
HWY_INLINE V eval_smart_median_float(D d, V center, V* values, V threshold) {
  constexpr int kSide = 2 * Radius + 1;
  constexpr int kCount = kSide * kSide;
  constexpr int kEvenCount = kCount - 1;
  constexpr int kCenterIdx = kEvenCount / 2;

  if constexpr (Radius == 1) {
    median8_even(d, values);
  } else if constexpr (Radius == 2) {
    median24_even(d, values);
  } else {
    median48_even(d, values);
  }

  const V median_left = values[kCenterIdx - 1];
  const V median_right = values[kCenterIdx];

  V sum = values[0];
  for (int i = 1; i < kEvenCount; ++i) {
    sum = float_add<IsF16>(d, sum, values[i]);
  }

  const V count_vec = hn::Set(d, static_cast<float>(kEvenCount));
  const V average = float_div<IsF16>(d, sum, count_vec);

  const V diff_l = float_sub<IsF16>(d, median_left, average);
  const V diff_r = float_sub<IsF16>(d, median_right, average);
  const V sq_l = float_mul<IsF16>(d, diff_l, diff_l);
  const V sq_r = float_mul<IsF16>(d, diff_r, diff_r);
  const V sq_sum = float_add<IsF16>(d, sq_l, sq_r);

  const V curved_variance = float_mul<IsF16>(d, float_sqrt<IsF16>(d, sq_sum), hn::Set(d, 13.0f));

  const auto lte = hn::Le(curved_variance, threshold);
  const auto med_lo = hn::Min(median_left, median_right);
  const auto med_hi = hn::Max(median_left, median_right);
  const auto clamped = hn::Clamp(center, med_lo, med_hi);

  return hn::IfThenElse(lte, clamped, center);
}

template <typename T, int Radius>
void smart_median_int_impl(T threshold, const T* srcp, T* dstp, int width, int height,
                           std::size_t src_stride, std::size_t dst_stride) {
  constexpr int kSide = 2 * Radius + 1;
  constexpr int kCount = kSide * kSide;
  constexpr int kEvenCount = kCount - 1;
  constexpr int kCenterOffset = kCount / 2;

  hn::ScalableTag<std::int32_t> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * Radius + kSimdPad;

  std::vector<std::int32_t> row_buffers(checked_product(kSide, padded_len));
  std::array<std::int32_t*, kSide> rows{};
  for (int i = 0; i < kSide; ++i) {
    rows[static_cast<std::size_t>(i)] = row_buffers.data() + static_cast<std::size_t>(i) * padded_len + Radius;
  }
  const hn::Rebind<T, decltype(d)> ds;
  std::array<int, kSide> cached_y;
  cached_y.fill(-1);
  const auto thresh_vec = hn::Set(d, static_cast<std::int32_t>(threshold));

  auto fill_i32_row = [&](std::int32_t* dst, const T* srow) {
    for (int x = -Radius; x < 0; ++x) dst[x] = srow[mirror_index(x, width)];
    for (int x = 0; x < width; ++x) dst[x] = srow[x];
    for (std::int64_t x = width; x < static_cast<std::int64_t>(width) + Radius; ++x) dst[x] = srow[mirror_index(x, width)];
  };

  for (int y = 0; y < height; ++y) {
    for (int dy = -Radius; dy <= Radius; ++dy) {
      const std::size_t my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      const auto slot = my % kSide;
      auto* row = row_buffers.data() + slot * padded_len + Radius;
      if (cached_y[slot] != static_cast<int>(my)) {
        fill_i32_row(row, srcp + my * src_stride);
        cached_y[slot] = static_cast<int>(my);
      }
      rows[static_cast<std::size_t>(dy + Radius)] = row;
    }

    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      hn::Vec<decltype(d)> values[kEvenCount];
      hn::Vec<decltype(d)> center{};
      int flat_idx = 0;
      int val_idx = 0;

      for (int ky = 0; ky < kSide; ++ky) {
        const std::int32_t* rptr = rows[static_cast<std::size_t>(ky)];
        for (int kx = -Radius; kx <= Radius; ++kx) {
          const auto val = hn::LoadU(d, rptr + x + kx);
          if (flat_idx == kCenterOffset) {
            center = val;
          } else {
            values[val_idx++] = val;
          }
          ++flat_idx;
        }
      }

      const auto res = eval_smart_median_int32<Radius>(d, center, values, thresh_vec);
      hn::StoreN(hn::DemoteTo(ds, res), ds, dst_row + x,
                 std::min(lanes, static_cast<std::size_t>(width) - x));
    }

  }
}

template <bool IsF16, typename StorageT, int Radius>
void smart_median_float_impl(float threshold, const StorageT* srcp, StorageT* dstp, int width, int height,
                             std::size_t src_stride, std::size_t dst_stride) {
  using ComputeT = FloatLane<IsF16>;
  constexpr int kSide = 2 * Radius + 1;
  constexpr int kCount = kSide * kSide;
  constexpr int kEvenCount = kCount - 1;
  constexpr int kCenterOffset = kCount / 2;

  hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * Radius + kSimdPad;

  std::vector<ComputeT> row_buffers(checked_product(kSide, padded_len));
  std::vector<ComputeT> out_f32(static_cast<std::size_t>(width) + kSimdPad);
  std::array<ComputeT*, kSide> rows{};
  for (int i = 0; i < kSide; ++i) {
    rows[static_cast<std::size_t>(i)] = row_buffers.data() + static_cast<std::size_t>(i) * padded_len + Radius;
  }
  const auto thresh_vec = hn::Set(d, threshold);

  for (int y = 0; y < height; ++y) {
    for (int dy = -Radius; dy <= Radius; ++dy) {
      const std::size_t my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      if constexpr (IsF16) {
        fill_mirrored_row_f16(rows[static_cast<std::size_t>(dy + Radius)] - Radius,
                                       srcp + my * src_stride, width, Radius);
      } else {
        fill_mirrored_row(rows[static_cast<std::size_t>(dy + Radius)] - Radius,
                          srcp + my * src_stride, width, Radius);
      }
    }

    StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      hn::Vec<decltype(d)> values[kEvenCount];
      hn::Vec<decltype(d)> center{};
      int flat_idx = 0;
      int val_idx = 0;

      for (int ky = 0; ky < kSide; ++ky) {
        const ComputeT* rptr = rows[static_cast<std::size_t>(ky)];
        for (int kx = -Radius; kx <= Radius; ++kx) {
          const auto val = hn::LoadU(d, rptr + x + kx);
          if (flat_idx == kCenterOffset) {
            center = val;
          } else {
            values[val_idx++] = val;
          }
          ++flat_idx;
        }
      }

      const auto res = eval_smart_median_float<IsF16, Radius>(d, center, values, thresh_vec);
      hn::StoreU(res, d, out_f32.data() + x);
    }

    if constexpr (IsF16) {
      store_row_f16(dst_row, out_f32.data(), width);
    } else {
      std::memcpy(dst_row, out_f32.data(), static_cast<std::size_t>(width) * sizeof(float));
    }
  }
}

void dispatch_smart_median_target(DataType dtype, int radius, float threshold, const std::uint8_t* srcp,
                                  std::uint8_t* dstp, std::size_t width, std::size_t height,
                                  std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);

  if (dtype == DataType::U8) {
    const auto th = static_cast<std::uint8_t>(std::clamp(threshold, 0.0f, 255.0f));
    if (radius == 1) smart_median_int_impl<std::uint8_t, 1>(th, srcp, dstp, w, h, src_stride_bytes, dst_stride_bytes);
    else if (radius == 2) smart_median_int_impl<std::uint8_t, 2>(th, srcp, dstp, w, h, src_stride_bytes, dst_stride_bytes);
    else if (radius == 3) smart_median_int_impl<std::uint8_t, 3>(th, srcp, dstp, w, h, src_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    const auto th = static_cast<std::uint16_t>(std::clamp(threshold, 0.0f, 65535.0f));
    auto* s = reinterpret_cast<const std::uint16_t*>(srcp);
    auto* d = reinterpret_cast<std::uint16_t*>(dstp);
    if (radius == 1) smart_median_int_impl<std::uint16_t, 1>(th, s, d, w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
    else if (radius == 2) smart_median_int_impl<std::uint16_t, 2>(th, s, d, w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
    else if (radius == 3) smart_median_int_impl<std::uint16_t, 3>(th, s, d, w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    // The reference converts the parameter to the plane's sample type first.
    if constexpr (HWY_HAVE_FLOAT16) threshold = fp16_to_fp32(fp32_to_fp16(threshold));
    auto* s = reinterpret_cast<const std::uint16_t*>(srcp);
    auto* d = reinterpret_cast<std::uint16_t*>(dstp);
    if (radius == 1) smart_median_float_impl<true, std::uint16_t, 1>(threshold, s, d, w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
    else if (radius == 2) smart_median_float_impl<true, std::uint16_t, 2>(threshold, s, d, w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
    else if (radius == 3) smart_median_float_impl<true, std::uint16_t, 3>(threshold, s, d, w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    auto* s = reinterpret_cast<const float*>(srcp);
    auto* d = reinterpret_cast<float*>(dstp);
    if (radius == 1) smart_median_float_impl<false, float, 1>(threshold, s, d, w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
    else if (radius == 2) smart_median_float_impl<false, float, 2>(threshold, s, d, w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
    else if (radius == 3) smart_median_float_impl<false, float, 3>(threshold, s, d, w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_smart_median_target);

void process_smart_median_plane(DataType dtype, int radius, float threshold, const std::uint8_t* srcp,
                                std::uint8_t* dstp, std::size_t width, std::size_t height,
                                std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_smart_median_target)(dtype, radius, threshold, srcp, dstp, width,
                                                    height, src_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
