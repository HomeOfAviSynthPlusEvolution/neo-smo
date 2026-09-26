#pragma once

#include <cstddef>
#include <cstring>

namespace neo_smo {

template <typename T>
inline void copy_plane(T* dstp, const T* srcp, std::size_t width, std::size_t height, std::size_t dst_stride, std::size_t src_stride) noexcept {
  if (width == 0 || height == 0 || width > (static_cast<std::size_t>(1) << 30)) {
    return;
  }
  for (std::size_t r = 0; r < height; ++r) {
    std::memcpy(dstp + r * dst_stride, srcp + r * src_stride, width * sizeof(T));
  }
}

template <typename T>
inline void copy_first_n_lines(T* dstp, const T* srcp, std::size_t width, std::size_t dst_stride, std::size_t src_stride, std::size_t n) noexcept {
  if (width == 0 || n == 0 || width > (static_cast<std::size_t>(1) << 30)) {
    return;
  }
  for (std::size_t r = 0; r < n; ++r) {
    std::memcpy(dstp + r * dst_stride, srcp + r * src_stride, width * sizeof(T));
  }
}

template <typename T>
inline void copy_last_n_lines(T* dstp, const T* srcp, std::size_t width, std::size_t height, std::size_t dst_stride, std::size_t src_stride, std::size_t n) noexcept {
  if (width == 0 || height == 0 || n == 0 || width > (static_cast<std::size_t>(1) << 30)) {
    return;
  }
  const std::size_t start = (height > n) ? (height - n) : 0;
  for (std::size_t r = start; r < height; ++r) {
    std::memcpy(dstp + r * dst_stride, srcp + r * src_stride, width * sizeof(T));
  }
}

} // namespace neo_smo
