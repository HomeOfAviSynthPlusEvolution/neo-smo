// Included inside the per-target Highway namespace. No include guard.
#include "common/sorting_networks.hpp"

template <typename T, int Radius>
void median_plane_impl(const T* srcp, T* dstp, int width, int height, std::size_t src_stride, std::size_t dst_stride) {
  constexpr int kSide = 2 * Radius + 1;
  constexpr int kCount = kSide * kSide;
  const hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * Radius + kSimdPad;

  std::vector<T> row_buffers(checked_product(kSide, padded_len));
  std::array<T*, kSide> rows{};
  for (int i = 0; i < kSide; ++i) {
    rows[static_cast<std::size_t>(i)] = row_buffers.data() + static_cast<std::size_t>(i) * padded_len + Radius;
  }

  for (int y = 0; y < height; ++y) {
    for (int dy = -Radius; dy <= Radius; ++dy) {
      const std::size_t my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      fill_mirrored_row(rows[static_cast<std::size_t>(dy + Radius)] - Radius, srcp + my * src_stride, width, Radius);
    }

    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      hn::Vec<decltype(d)> vals[kCount];
      int idx = 0;
      for (int ky = 0; ky < kSide; ++ky) {
        const T* rptr = rows[static_cast<std::size_t>(ky)];
        for (int kx = -Radius; kx <= Radius; ++kx) {
          vals[idx++] = hn::LoadU(d, rptr + x + kx);
        }
      }

      hn::Vec<decltype(d)> med;
      if constexpr (Radius == 1) {
        med = median9(d, vals);
      } else if constexpr (Radius == 2) {
        med = median25(d, vals);
      } else {
        med = median49(d, vals);
      }

      const std::size_t remaining = static_cast<std::size_t>(width - x);
      if (remaining >= lanes) {
        hn::StoreU(med, d, dst_row + x);
      } else {
        hn::StoreN(med, d, dst_row + x, remaining);
      }
    }
  }
}

template <int Radius>
void median_plane_f16_impl(const std::uint16_t* srcp, std::uint16_t* dstp, int width, int height,
                           std::size_t src_stride, std::size_t dst_stride) {
  constexpr int kSide = 2 * Radius + 1;
  constexpr int kCount = kSide * kSide;
  const hn::ScalableTag<float> d;
  const std::size_t lanes = hn::Lanes(d);
  const std::size_t kSimdPad = lanes;
  const std::size_t padded_len = static_cast<std::size_t>(width) + 2 * Radius + kSimdPad;

  std::vector<float> row_buffers(checked_product(kSide, padded_len));
  std::vector<float> dst_f32(static_cast<std::size_t>(width) + kSimdPad);
  std::array<float*, kSide> rows{};
  for (int i = 0; i < kSide; ++i) {
    rows[static_cast<std::size_t>(i)] = row_buffers.data() + static_cast<std::size_t>(i) * padded_len + Radius;
  }

  for (int y = 0; y < height; ++y) {
    for (int dy = -Radius; dy <= Radius; ++dy) {
      const std::size_t my = mirror_index(static_cast<std::int64_t>(y) + dy, height);
      fill_mirrored_row_fp16_to_fp32(rows[static_cast<std::size_t>(dy + Radius)] - Radius, srcp + my * src_stride,
                                     width, Radius);
    }

    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      hn::Vec<decltype(d)> vals[kCount];
      int idx = 0;
      for (int ky = 0; ky < kSide; ++ky) {
        const float* rptr = rows[static_cast<std::size_t>(ky)];
        for (int kx = -Radius; kx <= Radius; ++kx) {
          vals[idx++] = hn::LoadU(d, rptr + x + kx);
        }
      }

      hn::Vec<decltype(d)> med;
      if constexpr (Radius == 1) {
        med = median9(d, vals);
      } else if constexpr (Radius == 2) {
        med = median25(d, vals);
      } else {
        med = median49(d, vals);
      }

      hn::StoreU(med, d, dst_f32.data() + x);
    }

    convert_row_fp32_to_fp16(dstp + static_cast<std::size_t>(y) * dst_stride, dst_f32.data(), width);
  }
}
