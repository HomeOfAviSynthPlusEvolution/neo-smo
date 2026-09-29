#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace neo_smo {
namespace mini_deen_detail {
template <class T>
T load(const std::uint8_t* row, int x) {
  T value;
  std::memcpy(&value, row + static_cast<std::size_t>(x) * sizeof(T), sizeof(T));
  return value;
}
// Shared by the plain C++ fallback and SIMD border pixels.
template <class T>
T pixel(const std::uint8_t* src, std::ptrdiff_t src_stride, int width, int height, int radius, unsigned threshold,
        int x, int y) {
  const auto* row = src + y * src_stride;
  const unsigned center = load<T>(row, x);
  unsigned sum = 2 * center, count = 2;
  // Clip offsets, not coordinates: avoid overflow even for large dimensions.
  for (int dy = -std::min(y, radius); dy <= std::min(radius, height - 1 - y); ++dy) {
    const auto* neighbor_row = src + (y + dy) * src_stride;
    for (int dx = -std::min(x, radius); dx <= std::min(radius, width - 1 - x); ++dx) {
      const unsigned sample = load<T>(neighbor_row, x + dx);
      const unsigned difference = center > sample ? center - sample : sample - center;
      if (difference < threshold) {
        sum += sample;
        ++count;
      }
    }
  }
  // At most 227 * 65535; both the sum and rounding numerator fit uint32.
  const T result = static_cast<T>((2 * sum + count) / (2 * count));
  return result;
}
} // namespace mini_deen_detail
inline void mini_deen_scalar(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst,
                             std::ptrdiff_t dst_stride, int width, int height, bool byte_samples, int radius,
                             unsigned threshold) {
  const unsigned bytes = byte_samples ? 1 : 2;
  for (int y = 0; y < height; ++y) {
    auto* out = dst + y * dst_stride;
    if (threshold <= 1) {
      std::memcpy(out, src + y * src_stride, std::size_t(width) * bytes);
      continue;
    }
    for (int x = 0; x < width; ++x) {
      if (byte_samples)
        out[x] = mini_deen_detail::pixel<std::uint8_t>(src, src_stride, width, height, radius, threshold, x, y);
      else {
        const auto value =
            mini_deen_detail::pixel<std::uint16_t>(src, src_stride, width, height, radius, threshold, x, y);
        std::memcpy(out + std::size_t(x) * 2, &value, sizeof(value));
      }
    }
  }
}
} // namespace neo_smo
