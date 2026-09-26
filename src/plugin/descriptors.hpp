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
  }
}

} // namespace neo_smo::plugin
