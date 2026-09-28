#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/vertical_cleaner.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/grid.hpp"

#include "common/float_arithmetic.hpp"
#include "common/fp16_rows.hpp"

template <typename T>
void vertical_cleaner_int_impl(int mode, int bits_per_sample, const T* srcp, T* dstp, int width, int height,
                               std::size_t src_stride, std::size_t dst_stride) {
  const hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);

  if (mode == 1) {
    copy_first_n_lines(dstp, srcp, static_cast<std::size_t>(width), dst_stride, src_stride, 1);
    for (int y = 1; y < height - 1; ++y) {
      const T* r_top = srcp + static_cast<std::size_t>(y - 1) * src_stride;
      const T* r_cur = srcp + static_cast<std::size_t>(y) * src_stride;
      const T* r_bot = srcp + static_cast<std::size_t>(y + 1) * src_stride;
      T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto top = hn::LoadN(d, r_top + x, std::min(lanes, static_cast<std::size_t>(width - x)));
        const auto cur = hn::LoadN(d, r_cur + x, std::min(lanes, static_cast<std::size_t>(width - x)));
        const auto bot = hn::LoadN(d, r_bot + x, std::min(lanes, static_cast<std::size_t>(width - x)));
        const auto res = median3(d, top, cur, bot);
        const std::size_t rem = static_cast<std::size_t>(width - x);
        if (rem >= lanes)
          hn::StoreU(res, d, dst_row + x);
        else
          hn::StoreN(res, d, dst_row + x, rem);
      }
    }
    copy_last_n_lines(dstp, srcp, static_cast<std::size_t>(width), static_cast<std::size_t>(height), dst_stride,
                      src_stride, 1);
  } else if (mode == 2) {
    copy_first_n_lines(dstp, srcp, static_cast<std::size_t>(width), dst_stride, src_stride, 2);
    const T max_s = get_format_maximum<T>(bits_per_sample, false);
    const auto vmax = hn::Set(d, max_s);

    for (int y = 2; y < height - 2; ++y) {
      const T* rp2 = srcp + static_cast<std::size_t>(y - 2) * src_stride;
      const T* rp1 = srcp + static_cast<std::size_t>(y - 1) * src_stride;
      const T* rc = srcp + static_cast<std::size_t>(y) * src_stride;
      const T* rn1 = srcp + static_cast<std::size_t>(y + 1) * src_stride;
      const T* rn2 = srcp + static_cast<std::size_t>(y + 2) * src_stride;
      T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;

      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto p2 = hn::LoadN(d, rp2 + x, std::min(lanes, static_cast<std::size_t>(width - x)));
        const auto p1 = hn::LoadN(d, rp1 + x, std::min(lanes, static_cast<std::size_t>(width - x)));
        const auto c = hn::LoadN(d, rc + x, std::min(lanes, static_cast<std::size_t>(width - x)));
        const auto n1 = hn::LoadN(d, rn1 + x, std::min(lanes, static_cast<std::size_t>(width - x)));
        const auto n2 = hn::LoadN(d, rn2 + x, std::min(lanes, static_cast<std::size_t>(width - x)));

        const auto up_p = hn::Min(hn::SaturatedAdd(hn::Min(hn::SaturatedSub(p1, p2), vmax), p1), vmax);
        const auto up_n = hn::Min(hn::SaturatedAdd(hn::Min(hn::SaturatedSub(n1, n2), vmax), n1), vmax);
        const auto upper = hn::Max(hn::Max(hn::Min(up_p, up_n), p1), n1);

        const auto lo_p = hn::Min(hn::SaturatedSub(p1, hn::Min(hn::SaturatedSub(p2, p1), vmax)), vmax);
        const auto lo_n = hn::Min(hn::SaturatedSub(n1, hn::Min(hn::SaturatedSub(n2, n1), vmax)), vmax);
        const auto lower = hn::Min(hn::Min(p1, n1), hn::Max(lo_p, lo_n));

        const auto res = hn::Clamp(c, lower, upper);
        const std::size_t rem = static_cast<std::size_t>(width - x);
        if (rem >= lanes)
          hn::StoreU(res, d, dst_row + x);
        else
          hn::StoreN(res, d, dst_row + x, rem);
      }
    }
    copy_last_n_lines(dstp, srcp, static_cast<std::size_t>(width), static_cast<std::size_t>(height), dst_stride,
                      src_stride, 2);
  }
}

template <bool IsF16, typename StorageT>
void vertical_cleaner_float_impl(int mode, bool chroma, const StorageT* srcp, StorageT* dstp, int width, int height,
                                 std::size_t src_stride, std::size_t dst_stride) {
  const hn::ScalableTag<FloatLane<IsF16>> d;
  const std::size_t lanes = hn::Lanes(d);
  const int radius = mode == 1 ? 1 : 2;
  if (mode != 1 && mode != 2) return;
  copy_first_n_lines(dstp, srcp, static_cast<std::size_t>(width), dst_stride, src_stride, radius);
  const auto vmin = float_set(d, chroma ? -0.5f : 0.0f);
  const auto vmax = float_set(d, chroma ? 0.5f : 1.0f);
  for (int y = radius; y < height - radius; ++y) {
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
      auto load = [&](int dy) HWY_ATTR {
        const auto* row = srcp + static_cast<std::size_t>(y + dy) * src_stride + x;
        if constexpr (IsF16) return load_f16(d, row, count);
        else return hn::LoadN(d, row, count);
      };
      const auto p1 = load(-1), c = load(0), n1 = load(1);
      auto result = c;
      if (mode == 1) {
        result = median3(d, p1, c, n1);
      } else {
        const auto p2 = load(-2), n2 = load(2);
        const auto up_p = hn::Clamp(float_add<IsF16>(d, hn::Clamp(float_sub<IsF16>(d, p1, p2), vmin, vmax), p1), vmin, vmax);
        const auto up_n = hn::Clamp(float_add<IsF16>(d, hn::Clamp(float_sub<IsF16>(d, n1, n2), vmin, vmax), n1), vmin, vmax);
        const auto upper = hn::Max(hn::Max(hn::Min(up_p, up_n), p1), n1);

        const auto lo_p = hn::Clamp(float_sub<IsF16>(d, p1, hn::Clamp(float_sub<IsF16>(d, p2, p1), vmin, vmax)), vmin, vmax);
        const auto lo_n = hn::Clamp(float_sub<IsF16>(d, n1, hn::Clamp(float_sub<IsF16>(d, n2, n1), vmin, vmax)), vmin, vmax);
        const auto lower = hn::Min(hn::Min(p1, n1), hn::Max(lo_p, lo_n));

        result = hn::Clamp(c, lower, upper);
      }
      auto* dst = dstp + static_cast<std::size_t>(y) * dst_stride + x;
      if constexpr (IsF16) store_f16(d, result, dst, count);
      else hn::StoreN(result, d, dst, count);
    }
  }
  copy_last_n_lines(dstp, srcp, static_cast<std::size_t>(width), static_cast<std::size_t>(height),
                    dst_stride, src_stride, radius);
}

void dispatch_vertical_cleaner_target(DataType dtype, int mode, bool chroma, int bits_per_sample,
                                      const std::uint8_t* srcp, std::uint8_t* dstp, std::size_t width,
                                      std::size_t height, std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);
  if (dtype == DataType::U8) {
    vertical_cleaner_int_impl<std::uint8_t>(mode, bits_per_sample, srcp, dstp, w, h, src_stride_bytes,
                                            dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    vertical_cleaner_int_impl<std::uint16_t>(mode, bits_per_sample, reinterpret_cast<const std::uint16_t*>(srcp),
                                             reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2,
                                             dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    vertical_cleaner_float_impl<true, std::uint16_t>(mode, chroma, reinterpret_cast<const std::uint16_t*>(srcp),
                                                     reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2,
                                                     dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    vertical_cleaner_float_impl<false, float>(mode, chroma, reinterpret_cast<const float*>(srcp),
                                              reinterpret_cast<float*>(dstp), w, h, src_stride_bytes / 4,
                                              dst_stride_bytes / 4);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_vertical_cleaner_target);
void process_vertical_cleaner_plane(DataType dtype, int mode, bool chroma, int bits_per_sample,
                                    const std::uint8_t* srcp, std::uint8_t* dstp, std::size_t width, std::size_t height,
                                    std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_vertical_cleaner_target)(dtype, mode, chroma, bits_per_sample, srcp, dstp, width,
                                                         height, src_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
