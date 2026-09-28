#pragma once

#include "algorithms/plan.hpp"
#include "kernels/dispatch.hpp"
#include <dualsynth/param.hpp>

namespace neo_smo::plugin {

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
    case Algorithm::TemporalMedian:
      return "clip:vnode;radius:int:opt;planes:int[]:opt;scenechange:int:opt;";
    case Algorithm::TemporalSoften:
      return "clip:vnode;radius:int:opt;threshold:float[]:opt;scenechange:int:opt;scalep:int:opt;planes:int[]:opt;";
    case Algorithm::TemporalRepair:
      return "clip:vnode;repairclip:vnode;mode:int[]:opt;planes:int[]:opt;";
    case Algorithm::DegrainMedian:
      return "clip:vnode;limit:float[]:opt;mode:int[]:opt;interlaced:int:opt;norow:int:opt;scalep:int:opt;";
    case Algorithm::FluxSmoothT:
      return "clip:vnode;temporal_threshold:float[]:opt;planes:int[]:opt;scalep:int:opt;";
    case Algorithm::FluxSmoothST:
      return "clip:vnode;temporal_threshold:float[]:opt;spatial_threshold:float[]:opt;planes:int[]:opt;scalep:int:opt;";
    case Algorithm::TTempSmooth:
      return "clip:vnode;maxr:int:opt;thresh:int[]:opt;mdiff:int[]:opt;strength:int:opt;scthresh:float:opt;fp:int:opt;"
             "pfclip:vnode:opt;planes:int[]:opt;";
    case Algorithm::CCD:
      return "clip:vnode;threshold:float:opt;temporal_radius:int:opt;points:int[]:opt;scale:float:opt;ref:vnode:opt;";
    case Algorithm::Cnr4:
      return "clip:vnode;mode:data:opt;radius:int:opt;sense:int[]:opt;str:int[]:opt;pow:float[]:opt;tmode:int:opt;"
             "wmode:int:opt;scenechange:int:opt;ref:vnode:opt;";
    case Algorithm::DCTFilter:
      return "clip:vnode;factors:float[];planes:int[]:opt;";
  }
  return "clip:vnode;";
}

// One parameter list for both hosts; AviSynth uses native arrays and bools.
inline ds::FilterDescriptor descriptor(Algorithm alg) {
  using P = ds::ParamType;
  ds::FilterDescriptor d;
  d.name = algorithm_name(alg);
  const std::string signature = vs_signature(alg);
  for (std::size_t from = 0; from < signature.size();) {
    const auto end = signature.find(';', from);
    const auto item = signature.substr(from, end - from);
    const auto colon = item.find(':');
    const auto name = item.substr(0, colon);
    const auto type = item.substr(colon + 1);
    P kind = type.find("vnode") == 0   ? P::Clip
             : type.find("float") == 0 ? P::Float
             : type.find("data") == 0  ? P::String
                                       : P::Integer;
    if (name == "scalep" || name == "interlaced" || name == "norow" || name == "fp" ||
        (name == "scenechange" && alg != Algorithm::TemporalSoften))
      kind = P::Boolean;
    d.params.push_back(
        {name,
         kind,
         {},
         type.find(":opt") == std::string::npos,
         type.find("[]") != std::string::npos,
         true,
         true,
         type.find("[]") != std::string::npos ? ds::AvisynthArrayBinding::Native : ds::AvisynthArrayBinding::Legacy});
    from = end + 1;
  }
  return d;
}

inline void execute_plane(const FilterPlan& plan, int plane, const std::uint8_t* srcp, std::uint8_t* dstp,
                          std::size_t width, std::size_t height, std::size_t src_stride_bytes,
                          std::size_t dst_stride_bytes) {
  const DataType dtype = get_data_type(plan.format.bytes_per_sample, plan.format.is_float);
  const bool chroma = (plan.format.color_family == 3) && (plane > 0);
  const int param = plan.params[static_cast<std::size_t>(plane)];

  switch (plan.algorithm) {
    case Algorithm::Median:
      process_median_plane(dtype, param, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes);
      break;
    case Algorithm::VerticalCleaner:
      process_vertical_cleaner_plane(dtype, param, chroma, plan.format.bits_per_sample, srcp, dstp, width, height,
                                     src_stride_bytes, dst_stride_bytes);
      break;
    case Algorithm::RemoveGrain:
      process_remove_grain_plane(dtype, param, chroma, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes);
      break;
    case Algorithm::InterQuartileMean:
      process_inter_quartile_mean_plane(dtype, param, srcp, dstp, width, height, src_stride_bytes, dst_stride_bytes);
      break;
    case Algorithm::SmartMedian:
      process_smart_median_plane(dtype, param, plan.thresholds[static_cast<std::size_t>(plane)], srcp, dstp, width,
                                 height, src_stride_bytes, dst_stride_bytes);
      break;
    default:
      break;
  }
}

inline void execute_repair_plane(const FilterPlan& plan, int plane, const std::uint8_t* srcp,
                                 const std::uint8_t* repairp, std::uint8_t* dstp, std::size_t width, std::size_t height,
                                 std::size_t src_stride_bytes, std::size_t repair_stride_bytes,
                                 std::size_t dst_stride_bytes) {
  const DataType dtype = get_data_type(plan.format.bytes_per_sample, plan.format.is_float);
  const bool chroma = (plan.format.color_family == 3) && (plane > 0);
  const int mode = plan.params[static_cast<std::size_t>(plane)];
  process_repair_plane(dtype, mode, chroma, srcp, repairp, dstp, width, height, src_stride_bytes, repair_stride_bytes,
                       dst_stride_bytes);
}

inline void execute_clense_plane(const FilterPlan& plan, const std::uint8_t* srcp, const std::uint8_t* prevp,
                                 const std::uint8_t* nextp, std::uint8_t* dstp, std::size_t width, std::size_t height,
                                 std::size_t src_stride_bytes, std::size_t prev_stride_bytes,
                                 std::size_t next_stride_bytes, std::size_t dst_stride_bytes) {
  const DataType dtype = get_data_type(plan.format.bytes_per_sample, plan.format.is_float);
  process_clense_plane(dtype, srcp, prevp, nextp, dstp, width, height, src_stride_bytes, prev_stride_bytes,
                       next_stride_bytes, dst_stride_bytes);
}

inline void execute_clense_forward_backward_plane(const FilterPlan& plan, const std::uint8_t* srcp,
                                                  const std::uint8_t* ref1p, const std::uint8_t* ref2p,
                                                  std::uint8_t* dstp, std::size_t width, std::size_t height,
                                                  std::size_t src_stride_bytes, std::size_t ref1_stride_bytes,
                                                  std::size_t ref2_stride_bytes, std::size_t dst_stride_bytes) {
  const DataType dtype = get_data_type(plan.format.bytes_per_sample, plan.format.is_float);
  process_clense_forward_backward_plane(dtype, srcp, ref1p, ref2p, dstp, width, height, src_stride_bytes,
                                        ref1_stride_bytes, ref2_stride_bytes, dst_stride_bytes);
}

} // namespace neo_smo::plugin
