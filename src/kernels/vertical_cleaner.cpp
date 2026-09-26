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
  const hn::ScalableTag<float> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;

  auto load_row_f32 = [&](std::vector<float>& buf, const StorageT* row_ptr) {
    if constexpr (IsF16) {
      fill_mirrored_row_fp16_to_fp32(buf.data(), row_ptr, width, 0);
    } else {
      std::memcpy(buf.data(), row_ptr, static_cast<std::size_t>(width) * sizeof(float));
    }
  };

  if (mode == 1) {
    copy_first_n_lines(dstp, srcp, static_cast<std::size_t>(width), dst_stride, src_stride, 1);
    std::vector<float> r_top(static_cast<std::size_t>(width) + kSimdPad);
    std::vector<float> r_cur(static_cast<std::size_t>(width) + kSimdPad);
    std::vector<float> r_bot(static_cast<std::size_t>(width) + kSimdPad);
    std::vector<float> out_f32(static_cast<std::size_t>(width) + kSimdPad);

    for (int y = 1; y < height - 1; ++y) {
      load_row_f32(r_top, srcp + static_cast<std::size_t>(y - 1) * src_stride);
      load_row_f32(r_cur, srcp + static_cast<std::size_t>(y) * src_stride);
      load_row_f32(r_bot, srcp + static_cast<std::size_t>(y + 1) * src_stride);

      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto top = hn::LoadU(d, r_top.data() + x);
        const auto cur = hn::LoadU(d, r_cur.data() + x);
        const auto bot = hn::LoadU(d, r_bot.data() + x);
        const auto res = median3(d, top, cur, bot);
        hn::StoreU(res, d, out_f32.data() + x);
      }

      if constexpr (IsF16) {
        convert_row_fp32_to_fp16(dstp + static_cast<std::size_t>(y) * dst_stride, out_f32.data(), width);
      } else {
        std::memcpy(dstp + static_cast<std::size_t>(y) * dst_stride, out_f32.data(),
                    static_cast<std::size_t>(width) * sizeof(float));
      }
    }
    copy_last_n_lines(dstp, srcp, static_cast<std::size_t>(width), static_cast<std::size_t>(height), dst_stride,
                      src_stride, 1);
  } else if (mode == 2) {
    copy_first_n_lines(dstp, srcp, static_cast<std::size_t>(width), dst_stride, src_stride, 2);
    const float min_s = chroma ? -0.5f : 0.0f;
    const float max_s = chroma ? 0.5f : 1.0f;
    const auto vmin = hn::Set(d, min_s);
    const auto vmax = hn::Set(d, max_s);

    std::vector<float> rp2(static_cast<std::size_t>(width) + kSimdPad);
    std::vector<float> rp1(static_cast<std::size_t>(width) + kSimdPad);
    std::vector<float> rc(static_cast<std::size_t>(width) + kSimdPad);
    std::vector<float> rn1(static_cast<std::size_t>(width) + kSimdPad);
    std::vector<float> rn2(static_cast<std::size_t>(width) + kSimdPad);
    std::vector<float> out_f32(static_cast<std::size_t>(width) + kSimdPad);

    for (int y = 2; y < height - 2; ++y) {
      load_row_f32(rp2, srcp + static_cast<std::size_t>(y - 2) * src_stride);
      load_row_f32(rp1, srcp + static_cast<std::size_t>(y - 1) * src_stride);
      load_row_f32(rc, srcp + static_cast<std::size_t>(y) * src_stride);
      load_row_f32(rn1, srcp + static_cast<std::size_t>(y + 1) * src_stride);
      load_row_f32(rn2, srcp + static_cast<std::size_t>(y + 2) * src_stride);

      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto p2 = hn::LoadU(d, rp2.data() + x);
        const auto p1 = hn::LoadU(d, rp1.data() + x);
        const auto c = hn::LoadU(d, rc.data() + x);
        const auto n1 = hn::LoadU(d, rn1.data() + x);
        const auto n2 = hn::LoadU(d, rn2.data() + x);

        const auto up_p = hn::Clamp(float_add<IsF16>(d, hn::Clamp(float_sub<IsF16>(d, p1, p2), vmin, vmax), p1), vmin, vmax);
        const auto up_n = hn::Clamp(float_add<IsF16>(d, hn::Clamp(float_sub<IsF16>(d, n1, n2), vmin, vmax), n1), vmin, vmax);
        const auto upper = hn::Max(hn::Max(hn::Min(up_p, up_n), p1), n1);

        const auto lo_p = hn::Clamp(float_sub<IsF16>(d, p1, hn::Clamp(float_sub<IsF16>(d, p2, p1), vmin, vmax)), vmin, vmax);
        const auto lo_n = hn::Clamp(float_sub<IsF16>(d, n1, hn::Clamp(float_sub<IsF16>(d, n2, n1), vmin, vmax)), vmin, vmax);
        const auto lower = hn::Min(hn::Min(p1, n1), hn::Max(lo_p, lo_n));

        const auto res = hn::Clamp(c, lower, upper);
        hn::StoreU(res, d, out_f32.data() + x);
      }

      if constexpr (IsF16) {
        convert_row_fp32_to_fp16(dstp + static_cast<std::size_t>(y) * dst_stride, out_f32.data(), width);
      } else {
        std::memcpy(dstp + static_cast<std::size_t>(y) * dst_stride, out_f32.data(),
                    static_cast<std::size_t>(width) * sizeof(float));
      }
    }
    copy_last_n_lines(dstp, srcp, static_cast<std::size_t>(width), static_cast<std::size_t>(height), dst_stride,
                      src_stride, 2);
  }
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
