#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/inter_quartile_mean.cpp"
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
HWY_INLINE V eval_iqm_int(D d, V* vals) {
  const V three = hn::Set(d, 3);
  const V two = hn::Set(d, 2);

  if constexpr (Radius == 1) {
    sort9(d, vals);
    const V inner_sum = hn::Add(hn::Add(vals[3], vals[4]), vals[5]);
    const V ends = hn::Add(vals[2], vals[6]);
    const V weighted = hn::ShiftRight<2>(hn::Add(hn::Mul(ends, three), two));
    const V sum = hn::Add(inner_sum, weighted);
    const V four = hn::Set(d, 4);
    return hn::Add(hn::Mul(sum, two), four);
  } else if constexpr (Radius == 2) {
#if defined(NEO_SMO_MSVC_U16_SORT) && (HWY_TARGET == HWY_SSE2 || HWY_TARGET == HWY_SSSE3)
    if constexpr (sizeof(hn::TFromD<D>) == 4)
      sort_u16_samples<25>(d, vals, [](auto ds, auto* v) { sort25(ds, v); });
    else sort25(d, vals);
#else
    sort25(d, vals);
#endif
    V inner_sum = vals[7];
    for (int i = 8; i <= 17; ++i) {
      inner_sum = hn::Add(inner_sum, vals[i]);
    }
    const V ends = hn::Add(vals[6], vals[18]);
    const V twelve = hn::Set(d, 12);
    const V weighted = hn::ShiftRight<2>(hn::Add(hn::Mul(ends, three), two));
    const V sum = hn::Add(inner_sum, weighted);
    return hn::Add(hn::Mul(sum, two), twelve);
  } else {
#if defined(NEO_SMO_MSVC_U16_SORT) && (HWY_TARGET == HWY_SSE2 || HWY_TARGET == HWY_SSSE3)
    if constexpr (sizeof(hn::TFromD<D>) == 4)
      sort_u16_samples<49>(d, vals, [](auto ds, auto* v) { sort49(ds, v); });
    else sort49(d, vals);
#else
    sort49(d, vals);
#endif
    V inner_sum = vals[13];
    for (int i = 14; i <= 35; ++i) {
      inner_sum = hn::Add(inner_sum, vals[i]);
    }
    const V ends = hn::Add(vals[12], vals[36]);
    const V twenty_four = hn::Set(d, 24);
    const V weighted = hn::ShiftRight<2>(hn::Add(hn::Mul(ends, three), two));
    const V sum = hn::Add(inner_sum, weighted);
    return hn::Add(hn::Mul(sum, two), twenty_four);
  }
}

template <bool IsF16, int Radius, class D, class V = hn::Vec<D>>
HWY_INLINE V eval_iqm_float(D d, V* vals) {
  const V three_quarters = float_set(d, 0.75f);
  if constexpr (Radius == 1) {
    sort9(d, vals);
    const V inner_sum = float_add<IsF16>(d, float_add<IsF16>(d, vals[3], vals[4]), vals[5]);
    const V ends = float_add<IsF16>(d, vals[2], vals[6]);
    const V weighted = float_mul<IsF16>(d, ends, three_quarters);
    const V sum = float_add<IsF16>(d, inner_sum, weighted);
    const V len_half = float_set(d, 4.5f);
    return float_div<IsF16>(d, sum, len_half);
  } else if constexpr (Radius == 2) {
    sort25(d, vals);
    V inner_sum = vals[7];
    for (int i = 8; i <= 17; ++i) {
      inner_sum = float_add<IsF16>(d, inner_sum, vals[i]);
    }
    const V ends = float_add<IsF16>(d, vals[6], vals[18]);
    const V weighted = float_mul<IsF16>(d, ends, three_quarters);
    const V sum = float_add<IsF16>(d, inner_sum, weighted);
    const V len_half = float_set(d, 12.5f);
    return float_div<IsF16>(d, sum, len_half);
  } else {
    sort49(d, vals);
    V inner_sum = vals[13];
    for (int i = 14; i <= 35; ++i) {
      inner_sum = float_add<IsF16>(d, inner_sum, vals[i]);
    }
    const V ends = float_add<IsF16>(d, vals[12], vals[36]);
    const V weighted = float_mul<IsF16>(d, ends, three_quarters);
    const V sum = float_add<IsF16>(d, inner_sum, weighted);
    const V len_half = float_set(d, 24.5f);
    return float_div<IsF16>(d, sum, len_half);
  }
}

// Sort in native sample lanes, then widen only the retained quartiles.
template <int Radius, class D, class T>
void iqm_native_int(D d, const T* srcp, T* dstp,
    int width, int height, std::size_t src_stride, std::size_t dst_stride) {
  constexpr int side = 2 * Radius + 1, count = side * side, quartile = count / 4;
  using Wide = std::conditional_t<sizeof(T) == 1, std::int16_t, std::int32_t>;
  const hn::Half<D> dh;
  const hn::Rebind<Wide, decltype(dh)> dw;
  const auto lanes = hn::Lanes(d);
  const auto padded_len = static_cast<std::size_t>(width) + 2 * Radius + lanes;
  std::vector<T> buffer(checked_product(side, padded_len));
  std::array<int, side> cached_y;
  cached_y.fill(-1);
  std::array<const T*, side> rows{};
  for (int y = 0; y < height; ++y) {
    for (int dy = -Radius; dy <= Radius; ++dy) {
      const auto my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      const auto slot = my % side;
      auto* row = buffer.data() + slot * padded_len + Radius;
      if (cached_y[slot] != static_cast<int>(my)) {
        fill_mirrored_row(row - Radius, srcp + my * src_stride, width, Radius);
        cached_y[slot] = static_cast<int>(my);
      }
      rows[dy + Radius] = row;
    }
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      hn::Vec<D> values[count];
      for (int ky = 0; ky < side; ++ky) for (int kx = -Radius; kx <= Radius; ++kx)
        values[ky * side + kx + Radius] = hn::LoadU(d, rows[ky] + x + kx);
      if constexpr (Radius == 1) sort9(d, values);
      else sort25(d, values);
      auto average = [&](auto upper) HWY_ATTR {
        auto widen = [&](int i) HWY_ATTR {
          // Defer lookup via ADL: the scalar target has no UpperHalf/Combine.
          if constexpr (decltype(upper)::value) return hn::PromoteTo(dw, UpperHalf(dh, values[i]));
          else return hn::PromoteTo(dw, hn::LowerHalf(dh, values[i]));
        };
        auto sum = widen(quartile + 1);
        for (int i = quartile + 2; i < count - quartile - 1; ++i) sum = hn::Add(sum, widen(i));
        const auto ends = hn::Add(widen(quartile), widen(count - quartile - 1));
        sum = hn::Add(sum, hn::ShiftRight<2>(hn::Add(hn::Mul(ends, hn::Set(dw, 3)), hn::Set(dw, 2))));
        const auto numerator = hn::Add(hn::Add(sum, sum), hn::Set(dw, count / 2));
        const auto divisor = hn::Set(dw, count);
        auto q = hn::MulHigh(numerator, hn::Set(dw, ((std::uint64_t{1} << (sizeof(Wide) * 8)) + count - 1) / count));
        q = hn::Sub(q, hn::IfThenElse(hn::Gt(hn::Mul(q, divisor), numerator), hn::Set(dw, 1), hn::Zero(dw)));
        return hn::DemoteTo(dh, q);
      };
      const auto lo = average(std::false_type{}), hi = average(std::true_type{});
      hn::StoreN(Combine(d, hi, lo), d, dstp + static_cast<std::size_t>(y) * dst_stride + x,
                 std::min(lanes, static_cast<std::size_t>(width) - x));
    }
  }
}

template <typename T, int Radius>
void iqm_int_impl(const T* srcp, T* dstp, int width, int height, std::size_t src_stride, std::size_t dst_stride) {
  if constexpr (Radius <= 2) {
    const hn::ScalableTag<T> packed;
    if constexpr (hn::MaxLanes(packed) >= 2) {
      iqm_native_int<Radius>(packed, srcp, dstp, width, height, src_stride, dst_stride);
      return;
    }
  }
  constexpr int kSide = 2 * Radius + 1;
  constexpr int kCount = kSide * kSide;
  using ComputeT = std::conditional_t<sizeof(T) == 1, std::int16_t, std::int32_t>;
  hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * Radius + kSimdPad;

  std::vector<ComputeT> row_buffers(checked_product(kSide, padded_len));
  std::array<ComputeT*, kSide> rows{};
  for (int i = 0; i < kSide; ++i) {
    rows[static_cast<std::size_t>(i)] = row_buffers.data() + static_cast<std::size_t>(i) * padded_len + Radius;
  }
  std::array<int, kSide> cached_y;
  cached_y.fill(-1);
  const hn::Rebind<T, decltype(d)> ds;

  auto fill_i32_row = [&](ComputeT* dst, const T* srow) {
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
      hn::Vec<decltype(d)> vals[kCount];
      int idx = 0;
      for (int ky = 0; ky < kSide; ++ky) {
        const ComputeT* rptr = rows[static_cast<std::size_t>(ky)];
        for (int kx = -Radius; kx <= Radius; ++kx) {
          vals[idx++] = hn::LoadU(d, rptr + x + kx);
        }
      }

      const auto res = eval_iqm_int<Radius>(d, vals);
      // Rounded-up reciprocal gives at most one excess; correct it exactly.
      constexpr std::uint64_t scale = std::uint64_t{1} << (sizeof(ComputeT) * 8);
      const auto divisor = hn::Set(d, kCount);
      auto quotient = hn::MulHigh(res, hn::Set(d, (scale + kCount - 1) / kCount));
      quotient = hn::Sub(quotient, hn::IfThenElse(hn::Gt(hn::Mul(quotient, divisor), res),
                                                hn::Set(d, 1), hn::Zero(d)));
      hn::StoreN(hn::DemoteTo(ds, quotient), ds, dst_row + x,
                 std::min(lanes, static_cast<std::size_t>(width) - x));
    }

  }
}

template <bool IsF16, typename StorageT, int Radius>
void iqm_float_impl(const StorageT* srcp, StorageT* dstp, int width, int height, std::size_t src_stride,
                    std::size_t dst_stride) {
  using ComputeT = FloatLane<IsF16>;
  constexpr int kSide = 2 * Radius + 1;
  constexpr int kCount = kSide * kSide;
  hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * Radius + kSimdPad;

  std::vector<ComputeT> row_buffers(checked_product(kSide, padded_len));
  std::vector<ComputeT> out_f32(IsF16 ? static_cast<std::size_t>(width) + kSimdPad : 0);
  std::array<int, kSide> cached_y;
  cached_y.fill(-1);
  std::array<ComputeT*, kSide> rows{};
  for (int i = 0; i < kSide; ++i) {
    rows[static_cast<std::size_t>(i)] = row_buffers.data() + static_cast<std::size_t>(i) * padded_len + Radius;
  }

  for (int y = 0; y < height; ++y) {
    for (int dy = -Radius; dy <= Radius; ++dy) {
      const auto my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      const auto slot = my % kSide;
      auto* row = row_buffers.data() + slot * padded_len + Radius;
      if (cached_y[slot] != static_cast<int>(my)) {
        if constexpr (IsF16) fill_mirrored_row_f16(row - Radius, srcp + my * src_stride, width, Radius);
        else fill_mirrored_row(row - Radius, srcp + my * src_stride, width, Radius);
        cached_y[slot] = static_cast<int>(my);
      }
      rows[dy + Radius] = row;
    }

    StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      hn::Vec<decltype(d)> vals[kCount];
      int idx = 0;
      for (int ky = 0; ky < kSide; ++ky) {
        const ComputeT* rptr = rows[static_cast<std::size_t>(ky)];
        for (int kx = -Radius; kx <= Radius; ++kx) {
          vals[idx++] = hn::LoadU(d, rptr + x + kx);
        }
      }

      const auto res = eval_iqm_float<IsF16, Radius>(d, vals);
      if constexpr (IsF16) hn::StoreU(res, d, out_f32.data() + x);
      else hn::StoreN(res, d, dst_row + x, std::min(lanes, static_cast<std::size_t>(width) - x));
    }

    if constexpr (IsF16) {
      store_row_f16(dst_row, out_f32.data(), width);
    }
  }
}

void dispatch_inter_quartile_mean_target(DataType dtype, int radius, const std::uint8_t* srcp,
                                        std::uint8_t* dstp, std::size_t width, std::size_t height,
                                        std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);

  auto dispatch_radius = [&](auto tag_type) {
    using T = decltype(tag_type);
    if constexpr (std::is_same_v<T, float>) {
      if (radius == 1) iqm_float_impl<false, float, 1>(reinterpret_cast<const float*>(srcp), reinterpret_cast<float*>(dstp), w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
      else if (radius == 2) iqm_float_impl<false, float, 2>(reinterpret_cast<const float*>(srcp), reinterpret_cast<float*>(dstp), w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
      else if (radius == 3) iqm_float_impl<false, float, 3>(reinterpret_cast<const float*>(srcp), reinterpret_cast<float*>(dstp), w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
    } else if constexpr (std::is_same_v<T, double>) { // used as tag for F16
      if (radius == 1) iqm_float_impl<true, std::uint16_t, 1>(reinterpret_cast<const std::uint16_t*>(srcp), reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
      else if (radius == 2) iqm_float_impl<true, std::uint16_t, 2>(reinterpret_cast<const std::uint16_t*>(srcp), reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
      else if (radius == 3) iqm_float_impl<true, std::uint16_t, 3>(reinterpret_cast<const std::uint16_t*>(srcp), reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
    } else {
      if (radius == 1) iqm_int_impl<T, 1>(reinterpret_cast<const T*>(srcp), reinterpret_cast<T*>(dstp), w, h, src_stride_bytes / sizeof(T), dst_stride_bytes / sizeof(T));
      else if (radius == 2) iqm_int_impl<T, 2>(reinterpret_cast<const T*>(srcp), reinterpret_cast<T*>(dstp), w, h, src_stride_bytes / sizeof(T), dst_stride_bytes / sizeof(T));
      else if (radius == 3) iqm_int_impl<T, 3>(reinterpret_cast<const T*>(srcp), reinterpret_cast<T*>(dstp), w, h, src_stride_bytes / sizeof(T), dst_stride_bytes / sizeof(T));
    }
  };

  if (dtype == DataType::U8) {
    dispatch_radius(std::uint8_t{});
  } else if (dtype == DataType::U16) {
    dispatch_radius(std::uint16_t{});
  } else if (dtype == DataType::F16) {
    dispatch_radius(double{});
  } else if (dtype == DataType::F32) {
    dispatch_radius(float{});
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_inter_quartile_mean_target);

void process_inter_quartile_mean_plane(DataType dtype, int radius, const std::uint8_t* srcp,
                                       std::uint8_t* dstp, std::size_t width, std::size_t height,
                                       std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_inter_quartile_mean_target)(dtype, radius, srcp, dstp, width, height,
                                                           src_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
