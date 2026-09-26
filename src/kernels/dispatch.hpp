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

void process_repair_plane(
  DataType dtype,
  int mode,
  bool chroma,
  const std::uint8_t* srcp,
  const std::uint8_t* repairp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t repair_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_clense_plane(
  DataType dtype,
  const std::uint8_t* srcp,
  const std::uint8_t* prevp,
  const std::uint8_t* nextp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t prev_stride_bytes,
  std::size_t next_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_clense_forward_backward_plane(
  DataType dtype,
  const std::uint8_t* srcp,
  const std::uint8_t* ref1p,
  const std::uint8_t* ref2p,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t ref1_stride_bytes,
  std::size_t ref2_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_inter_quartile_mean_plane(
  DataType dtype,
  int radius,
  const std::uint8_t* srcp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_smart_median_plane(
  DataType dtype,
  int radius,
  float threshold,
  const std::uint8_t* srcp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_temporal_median_plane(
  DataType dtype,
  int diameter,
  const std::uint8_t* const* srcp_planes,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_temporal_soften_plane(
  DataType dtype,
  int diameter,
  float threshold,
  const std::uint8_t* const* srcp_planes,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_temporal_repair_plane(
  DataType dtype,
  int mode,
  bool chroma,
  int bits_per_sample,
  const std::uint8_t* srcp,
  const std::uint8_t* prevp,
  const std::uint8_t* currp,
  const std::uint8_t* nextp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t prev_stride_bytes,
  std::size_t curr_stride_bytes,
  std::size_t next_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_degrain_median_plane(
  DataType dtype,
  int mode,
  float limit,
  bool interlaced,
  bool norow,
  bool chroma,
  int bits_per_sample,
  const std::uint8_t* prevp,
  const std::uint8_t* currp,
  const std::uint8_t* nextp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t prev_stride_bytes,
  std::size_t curr_stride_bytes,
  std::size_t next_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_fluxsmooth_t_plane(
  DataType dtype,
  float temporal_threshold,
  const std::uint8_t* prevp,
  const std::uint8_t* currp,
  const std::uint8_t* nextp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t prev_stride_bytes,
  std::size_t curr_stride_bytes,
  std::size_t next_stride_bytes,
  std::size_t dst_stride_bytes
);

void process_fluxsmooth_st_plane(
  DataType dtype,
  float temporal_threshold,
  float spatial_threshold,
  const std::uint8_t* prevp,
  const std::uint8_t* currp,
  const std::uint8_t* nextp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t prev_stride_bytes,
  std::size_t curr_stride_bytes,
  std::size_t next_stride_bytes,
  std::size_t dst_stride_bytes
);

} // namespace neo_smo
