#pragma once

#include "base/types.hpp"
#include <cstddef>
#include <cstdint>
#include <array>

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

void process_ttempsmooth_plane(
  DataType dtype,
  int width,
  int height,
  std::size_t src_stride_bytes,
  std::size_t ref_stride_bytes,
  std::size_t dst_stride_bytes,
  const std::uint8_t* curr,
  const std::uint8_t* curr_ref,
  const std::uint8_t* const* prev,
  const std::uint8_t* const* prev_ref,
  const std::uint8_t* const* next,
  const std::uint8_t* const* next_ref,
  std::uint8_t* dstp,
  int maxr,
  int num_prev,
  int num_next,
  float threshold,
  bool fp,
  int bits_per_sample,
  int weight_mode,
  float center_weight,
  const float* temporal_weights,
  const float* temporal_difference_weights
);

struct Point {
  int x;
  int y;
};

void process_ccd_planes(
  DataType dtype,
  bool is_rgb,
  int width,
  int height,
  std::size_t stride_bytes,
  const std::uint8_t* const* src,
  const std::uint8_t* const* ref,
  std::uint8_t* dst_r,
  std::uint8_t* dst_g,
  std::uint8_t* dst_b,
  float threshold,
  int temporal_radius,
  const float* weights,
  const Point* points,
  int num_points,
  int diameter,
  float scale,
  int bits_per_sample
);

void process_cnr4_frame(
  DataType dtype,
  int width,
  int height,
  std::size_t stride_bytes,
  int depth,
  int radius,
  int tmode,
  int wmode,
  const std::uint8_t* const curr[3],
  const std::uint8_t* const curr_ref[3],
  const std::array<const std::uint8_t*, 3>* src,
  const std::array<const std::uint8_t*, 3>* ref,
  std::size_t num_frames,
  std::uint8_t* dst_u,
  std::uint8_t* dst_v,
  const std::uint8_t* table_y,
  const std::uint8_t* table_u,
  const std::uint8_t* table_v
);

} // namespace neo_smo
