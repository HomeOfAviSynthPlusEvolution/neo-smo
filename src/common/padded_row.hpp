#pragma once

#include "base/fp16.hpp"
#include "base/checked.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

namespace neo_smo {

// Exact mirror indexing matching zsmooth's Grid.initFromCenterMirrored / ArrayGrid.initFromCenterMirrored:
// min_idx = abs(coord), clamped = min(min_idx, 2 * (len - 1) - coord)
inline std::size_t mirror_index(std::int64_t coord, int length) noexcept {
  if (length <= 1) {
    return 0;
  }
  const std::int64_t min_idx = std::abs(coord);
  const std::int64_t reflected = 2 * (static_cast<std::int64_t>(length) - 1) - coord;
  const std::int64_t clamped = std::min(min_idx, reflected);
  return static_cast<std::size_t>(std::clamp<std::int64_t>(clamped, 0, length - 1));
}

// Fill the useful neighborhood of a zero-initialized padded row buffer.
// Callers allocate width + 2 * radius + Lanes(d) elements; extra lanes remain zero.
// Element `dst[radius + x]` corresponds to `src_row[x]` for x in [0, width),
// with mirrored borders for x in [-radius, 0) and [width, width + radius).
template <typename T>
inline void fill_mirrored_row(T* dst, const T* src_row, int width, int radius) noexcept {
  for (int x = -radius; x < 0; ++x) {
    dst[radius + x] = src_row[mirror_index(x, width)];
  }
  std::memcpy(dst + radius, src_row, static_cast<std::size_t>(width) * sizeof(T));
  for (std::int64_t x = width; x < static_cast<std::int64_t>(width) + radius; ++x) {
    dst[radius + x] = src_row[mirror_index(x, width)];
  }
}

// Fill a float padded row buffer from an IEEE-754 FP16 (`uint16_t`) source row.
inline void fill_mirrored_row_fp16_to_fp32(float* dst, const std::uint16_t* src_row, int width, int radius) noexcept {
  for (std::int64_t x = -radius; x < static_cast<std::int64_t>(width) + radius; ++x) {
    const std::size_t idx = (x >= 0 && x < width) ? static_cast<std::size_t>(x) : mirror_index(x, width);
    dst[radius + x] = fp16_to_fp32(src_row[idx]);
  }
}

inline void convert_row_fp32_to_fp16(std::uint16_t* dst_row, const float* src_row, int width) noexcept {
  for (int x = 0; x < width; ++x) {
    dst_row[x] = fp32_to_fp16(src_row[x]);
  }
}

} // namespace neo_smo
