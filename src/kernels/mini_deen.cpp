// Algorithm: dubhater/vapoursynth-minideen (ISC), src/minideen.cpp,
// revision 8b915c5134352adc19539560ab5debbbcdbe0ea1.
// Retains the scalar formula, strict threshold, clipped window and extra center weight.
#include "kernels/mini_deen.hpp"
#include <algorithm>
#include <cstring>
namespace neo_smo {
namespace {
template <class T>
T load(const std::uint8_t* row, int x) {
  T value;
  std::memcpy(&value, row + static_cast<std::size_t>(x) * sizeof(T), sizeof(T));
  return value;
}
template <class T>
void process(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
             int width, int height, int radius, unsigned threshold) {
  for (int y = 0; y < height; ++y) {
    const auto* row = src + y * src_stride;
    auto* out = dst + y * dst_stride;
    // Test AFTER scaling: threshold=1 can filter high-bit-depth input.
    if (threshold <= 1) {
      std::memcpy(out, row, static_cast<std::size_t>(width) * sizeof(T));
      continue;
    }
    for (int x = 0; x < width; ++x) {
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
      std::memcpy(out + static_cast<std::size_t>(x) * sizeof(T), &result, sizeof(T));
    }
  }
}
} // namespace
void mini_deen_kernel(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
                      int width, int height, bool byte_samples, int radius, unsigned threshold) {
  if (byte_samples)
    process<std::uint8_t>(src, src_stride, dst, dst_stride, width, height, radius, threshold);
  else
    process<std::uint16_t>(src, src_stride, dst, dst_stride, width, height, radius, threshold);
}
} // namespace neo_smo
