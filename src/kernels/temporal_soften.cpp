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
  hn::ScalableTag<std::int32_t> d;
  const std::size_t lanes = hn::Lanes(d);
  const auto thresh_vec = hn::Set(d, static_cast<std::int32_t>(threshold));
  const int half_frames = diameter / 2;

  std::vector<std::int32_t> curr_i32(static_cast<std::size_t>(width) + lanes, 0);
  std::vector<std::int32_t> frame_i32(static_cast<std::size_t>(width) + lanes, 0);
  std::vector<std::int32_t> sum_i32(static_cast<std::size_t>(width) + lanes, 0);

  for (int y = 0; y < height; ++y) {
    const T* c_row = srcp_planes[0] + static_cast<std::size_t>(y) * src_stride;
    for (int x = 0; x < width; ++x) {
      curr_i32[static_cast<std::size_t>(x)] = static_cast<std::int32_t>(c_row[x]);
    }

    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto curr = hn::LoadU(d, curr_i32.data() + x);
      hn::StoreU(curr, d, sum_i32.data() + x);
    }

    for (int i = 1; i < diameter; ++i) {
      const T* f_row = srcp_planes[i] + static_cast<std::size_t>(y) * src_stride;
      for (int x = 0; x < width; ++x) {
        frame_i32[static_cast<std::size_t>(x)] = static_cast<std::int32_t>(f_row[x]);
      }

      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto curr = hn::LoadU(d, curr_i32.data() + x);
        const auto f_val = hn::LoadU(d, frame_i32.data() + x);
        const auto diff = hn::AbsDiff(curr, f_val);
        const auto chosen = hn::IfThenElse(hn::Le(diff, thresh_vec), f_val, curr);
        const auto old_sum = hn::LoadU(d, sum_i32.data() + x);
        hn::StoreU(hn::Add(old_sum, chosen), d, sum_i32.data() + x);
      }
    }

    for (int x = 0; x < width; ++x) {
      if constexpr (std::is_same_v<T, std::uint8_t>) {
        const std::uint32_t mul = (1u << 16) / static_cast<std::uint32_t>(diameter);
        const auto val = static_cast<std::uint32_t>(sum_i32[static_cast<std::size_t>(x)] + half_frames);
        dst_row[x] = static_cast<T>((val * mul) >> 16);
      } else {
        const std::uint64_t mul = (1ULL << 32) / static_cast<std::uint64_t>(diameter);
        const auto val = static_cast<std::uint64_t>(sum_i32[static_cast<std::size_t>(x)] + half_frames);
        dst_row[x] = static_cast<T>((val * mul) >> 32);
      }
    }
  }
}

template <bool IsF16, typename StorageT>
void temporal_soften_float_impl(int diameter, float threshold, const StorageT* const* srcp_planes,
                                StorageT* dstp, int width, int height, std::size_t src_stride,
                                std::size_t dst_stride) {
  const hn::ScalableTag<float> d;
  const std::size_t lanes = hn::Lanes(d);
  const auto thresh_vec = hn::Set(d, threshold);
  const auto frames_vec = hn::Set(d, static_cast<float>(diameter));

  if constexpr (IsF16) {
    std::vector<float> r_bufs(static_cast<std::size_t>(diameter) * (static_cast<std::size_t>(width) + lanes));
    std::vector<float> out_f32(static_cast<std::size_t>(width) + lanes);

    for (int y = 0; y < height; ++y) {
      for (int i = 0; i < diameter; ++i) {
        fill_mirrored_row_fp16_to_fp32(
            r_bufs.data() + static_cast<std::size_t>(i) * (static_cast<std::size_t>(width) + lanes),
            srcp_planes[i] + static_cast<std::size_t>(y) * src_stride, width, 0);
      }

      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto curr = hn::LoadU(d, r_bufs.data() + x);
        auto sum = curr;

        for (int i = 1; i < diameter; ++i) {
          const auto f_val = hn::LoadU(
              d, r_bufs.data() + static_cast<std::size_t>(i) * (static_cast<std::size_t>(width) + lanes) + x);
          const auto diff = hn::Abs(hn::Sub(curr, f_val)); // Zig compares in f32, including f16 input.
          const auto chosen = hn::IfThenElse(hn::Le(diff, thresh_vec), f_val, curr);
          sum = float_add<true>(d, sum, chosen);
        }

        const auto res = float_div<true>(d, sum, frames_vec);
        hn::StoreU(res, d, out_f32.data() + x);
      }

      convert_row_fp32_to_fp16(dst_row, out_f32.data(), width);
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

        const auto res = hn::Div(sum, frames_vec);
        if (rem >= lanes) {
          hn::StoreU(res, d, dst_row + x);
        } else {
          hn::StoreN(res, d, dst_row + x, rem);
        }
      }
    }
  }
}

void dispatch_temporal_soften_target(DataType dtype, int diameter, float threshold,
                                     const std::uint8_t* const* srcp_planes, std::uint8_t* dstp,
                                     std::size_t width, std::size_t height,
                                     std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  if (dtype == DataType::F16) {
    threshold = fp16_to_fp32(fp32_to_fp16(threshold));
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
    temporal_soften_float_impl<true, std::uint16_t>(
        diameter, threshold, u16_planes.data(),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    temporal_soften_float_impl<false, float>(
        diameter, threshold, f32_planes.data(),
        reinterpret_cast<float*>(dstp), w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
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
