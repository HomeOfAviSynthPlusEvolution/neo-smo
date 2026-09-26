#pragma once

#include "base/types.hpp"
#include <cstddef>
#include <cstdint>

namespace neo_smo {

void process_median_plane(
  DataType dtype,
  int radius,
  const std::uint8_t* srcp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_vertical_cleaner_plane(
  DataType dtype,
  int mode,
  bool chroma,
  int bits_per_sample,
  const std::uint8_t* srcp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_remove_grain_plane(
  DataType dtype,
  int mode,
  bool chroma,
  const std::uint8_t* srcp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t dst_stride_bytes
);

} // namespace neo_smo
