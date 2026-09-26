#pragma once

#include "algorithms/plan.hpp"
#include "kernels/dispatch.hpp"
#include <dualsynth/param.hpp>

namespace neo_smo::plugin {

inline ds::FilterDescriptor descriptor(Algorithm alg) {
  using P = ds::ParamType;
  ds::FilterDescriptor d;
  d.name = algorithm_name(alg);
  auto add = [&](const char* name, P type, bool required, bool is_array) {
    d.params.push_back({name, type, {}, required, is_array, true, true, ds::AvisynthArrayBinding::Native});
  };
  add("clip", P::Clip, true, false);
  if (alg == Algorithm::Median) {
    add("radius", P::Integer, false, true);
    add("planes", P::Integer, false, true);
  } else if (alg == Algorithm::VerticalCleaner) {
    add("mode", P::Integer, true, true);
  } else if (alg == Algorithm::RemoveGrain) {
    add("mode", P::Integer, true, true);
  }
  return d;
}

inline const char* vs_signature(Algorithm alg) noexcept {
  switch (alg) {
    case Algorithm::Median:
      return "clip:vnode;radius:int[]:opt;planes:int[]:opt;";
    case Algorithm::VerticalCleaner:
      return "clip:vnode;mode:int[];";
    case Algorithm::RemoveGrain:
      return "clip:vnode;mode:int[];";
    case Algorithm::Repair:
      return "clip:vnode;repairclip:vnode;mode:int[];";
    case Algorithm::Clense:
      return "clip:vnode;previous:vnode:opt;next:vnode:opt;planes:int[]:opt;";
    case Algorithm::ForwardClense:
    case Algorithm::BackwardClense:
      return "clip:vnode;planes:int[]:opt;";
    case Algorithm::InterQuartileMean:
      return "clip:vnode;radius:int[]:opt;planes:int[]:opt;";
    case Algorithm::SmartMedian:
      return "clip:vnode;radius:int[]:opt;threshold:float[]:opt;scalep:int:opt;planes:int[]:opt;";
  }
  return "clip:vnode;";
}

inline void execute_plane(
  const FilterPlan& plan,
  int plane,
  const std::uint8_t* srcp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t dst_stride_bytes
) {
  const DataType dtype = get_data_type(plan.format.bytes_per_sample, plan.format.is_float);
  const bool chroma = (plan.format.color_family == 3) && (plane > 0);
  const int param = plan.params[static_cast<std::size_t>(plane)];

  switch (plan.algorithm) {
    case Algorithm::Median:
      process_median_plane(dtype, param, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes);
      break;
    case Algorithm::VerticalCleaner:
      process_vertical_cleaner_plane(dtype, param, chroma, plan.format.bits_per_sample, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes);
      break;
    case Algorithm::RemoveGrain:
      process_remove_grain_plane(dtype, param, chroma, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes);
      break;
    case Algorithm::InterQuartileMean:
      process_inter_quartile_mean_plane(dtype, param, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes);
      break;
    case Algorithm::SmartMedian:
      process_smart_median_plane(dtype, param, plan.thresholds[static_cast<std::size_t>(plane)], srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes);
      break;
    default:
      break;
  }
}

inline void execute_repair_plane(
  const FilterPlan& plan,
  int plane,
  const std::uint8_t* srcp,
  const std::uint8_t* repairp,
  std::uint8_t* dstp,
  std::size_t width,
  std::size_t height,
  std::size_t src_stride_bytes,
  std::size_t repair_stride_bytes,
  std::size_t dst_stride_bytes
) {
  const DataType dtype = get_data_type(plan.format.bytes_per_sample, plan.format.is_float);
  const bool chroma = (plan.format.color_family == 3) && (plane > 0);
  const int mode = plan.params[static_cast<std::size_t>(plane)];
  process_repair_plane(dtype, mode, chroma, srcp, repairp, dstp, width, height, src_stride_bytes, repair_stride_bytes, dst_stride_bytes);
}

inline void execute_clense_plane(
  const FilterPlan& plan,
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
) {
  const DataType dtype = get_data_type(plan.format.bytes_per_sample, plan.format.is_float);
  process_clense_plane(dtype, srcp, prevp, nextp, dstp, width, height, src_stride_bytes, prev_stride_bytes, next_stride_bytes, dst_stride_bytes);
}

inline void execute_clense_forward_backward_plane(
  const FilterPlan& plan,
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
) {
  const DataType dtype = get_data_type(plan.format.bytes_per_sample, plan.format.is_float);
  process_clense_forward_backward_plane(dtype, srcp, ref1p, ref2p, dstp, width, height, src_stride_bytes, ref1_stride_bytes, ref2_stride_bytes, dst_stride_bytes);
}

} // namespace neo_smo::plugin
