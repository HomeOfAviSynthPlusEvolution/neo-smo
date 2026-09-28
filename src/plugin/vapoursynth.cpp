#include "plugin/descriptors.hpp"
#include "plugin/dct_filter.hpp"
#include "common/temporal_window.hpp"
#include "neo_smo_version.hpp"
#include <vapoursynth/VapourSynth4.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <limits>

namespace neo_smo::plugin {
namespace {

std::array<bool, 3> normalize_planes(int num_planes, const std::vector<int>& planes_list, bool planes_specified, const std::string& prefix) {
  if (!planes_specified) {
    return {true, true, true};
  }
  std::array<bool, 3> out{false, false, false};
  for (int p : planes_list) {
    require(p >= 0 && p < num_planes, prefix + ": Plane index out of range.");
    require(!out[static_cast<std::size_t>(p)], prefix + ": Plane specified twice.");
    out[static_cast<std::size_t>(p)] = true;
  }
  return out;
}

struct VSFilterInstance {
  VSNode* node = nullptr;
  VSVideoInfo vi{};
  FilterPlan plan{};
};

struct VSRepairInstance {
  VSNode* node = nullptr;
  VSNode* repair_node = nullptr;
  VSVideoInfo vi{};
  FilterPlan plan{};
};

struct VSClenseInstance {
  VSNode* cnode = nullptr;
  VSNode* pnode = nullptr;
  VSNode* nnode = nullptr;
  VSVideoInfo vi{};
  FilterPlan plan{};
};

FormatInfo extract_format_info(const VSVideoInfo* vi) {
  FormatInfo fmt{};
  if (vi->format.colorFamily == cfGray)
    fmt.color_family = 1;
  else if (vi->format.colorFamily == cfRGB)
    fmt.color_family = 2;
  else if (vi->format.colorFamily == cfYUV)
    fmt.color_family = 3;
  fmt.is_float = (vi->format.sampleType == stFloat);
  fmt.bits_per_sample = vi->format.bitsPerSample;
  fmt.bytes_per_sample = vi->format.bytesPerSample;
  fmt.num_planes = vi->format.numPlanes;
  fmt.subsampling_w = vi->format.subSamplingW;
  fmt.subsampling_h = vi->format.subSamplingH;
  fmt.width = vi->width;
  fmt.height = vi->height;
  return fmt;
}

bool is_same_video_info(const VSVideoInfo* a, const VSVideoInfo* b) {
  return a->format.colorFamily == b->format.colorFamily &&
         a->format.sampleType == b->format.sampleType &&
         a->format.bitsPerSample == b->format.bitsPerSample &&
         a->format.bytesPerSample == b->format.bytesPerSample &&
         a->format.subSamplingW == b->format.subSamplingW &&
         a->format.subSamplingH == b->format.subSamplingH &&
         a->format.numPlanes == b->format.numPlanes &&
         a->width == b->width &&
         a->height == b->height;
}

int read_int(const VSMap* in, const char* key, int index, const VSAPI* vsapi) {
  int err = 0;
  const std::int64_t value = vsapi->mapGetInt(in, key, index, &err);
  require(!err, "neo-smo: invalid integer parameter.");
  require(value >= std::numeric_limits<int>::min() && value <= std::numeric_limits<int>::max(),
          "neo-smo: integer parameter out of range.");
  return static_cast<int>(value);
}

// ---------------------------------------------------------------------------
// Single-clip spatial filters (Median, VerticalCleaner, RemoveGrain)
// ---------------------------------------------------------------------------
const VSFrame* VS_CC filter_get_frame(int n, int activation_reason, void* instance_data, void**,
                                      VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSFilterInstance*>(instance_data);

  try {
    if (activation_reason == arInitial) {
      vsapi->requestFrameFilter(n, d->node, frame_ctx);
    } else if (activation_reason == arAllFramesReady) {
      const auto free_frame = [vsapi](const VSFrame* f) {
        vsapi->freeFrame(f);
      };
      std::unique_ptr<const VSFrame, decltype(free_frame)> src_owner(vsapi->getFrameFilter(n, d->node, frame_ctx),
                                                                     free_frame);
      const VSFrame* src = src_owner.get();
      if (!src) {
        return nullptr;
      }

      const VSFrame* plane_src[3] = {
          d->plan.process[0] ? nullptr : src,
          d->plan.process[1] ? nullptr : src,
          d->plan.process[2] ? nullptr : src,
      };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, src, core);

      if (!dst)
        throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      for (int plane = 0; plane < d->vi.format.numPlanes; ++plane) {
        if (!d->plan.process[static_cast<std::size_t>(plane)]) {
          continue;
        }
        const auto width = static_cast<std::size_t>(vsapi->getFrameWidth(dst, plane));
        const auto height = static_cast<std::size_t>(vsapi->getFrameHeight(dst, plane));
        const auto src_stride = static_cast<std::size_t>(vsapi->getStride(src, plane));
        const auto dst_stride = static_cast<std::size_t>(vsapi->getStride(dst, plane));
        const auto* srcp = vsapi->getReadPtr(src, plane);
        auto* dstp = vsapi->getWritePtr(dst, plane);

        execute_plane(d->plan, plane, srcp, dstp, width, height, src_stride, dst_stride);
      }

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("neo-smo: frame processing failed", frame_ctx);
  }
  return nullptr;
}

void VS_CC filter_free(void* instance_data, VSCore*, const VSAPI* vsapi) {
  auto* d = static_cast<VSFilterInstance*>(instance_data);
  if (d) {
    if (d->node) {
      vsapi->freeNode(d->node);
    }
    delete d;
  }
}

template <Algorithm Alg>
void VS_CC filter_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* node = nullptr;
  try {
    int err = 0;
    node = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && node != nullptr, std::string(algorithm_name(Alg)) + ": clip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(node);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            std::string(algorithm_name(Alg)) + ": only constant format and dimensions supported.");

    FormatInfo fmt = extract_format_info(vi);

    const char* param_key = (Alg == Algorithm::Median || Alg == Algorithm::InterQuartileMean || Alg == Algorithm::SmartMedian) ? "radius" : "mode";
    std::vector<int> param_list;
    const int num_params = vsapi->mapNumElements(in, param_key);
    for (int i = 0; i < num_params; ++i) {
      param_list.push_back(read_int(in, param_key, i, vsapi));
    }

    std::vector<float> threshold_list;
    bool scalep = false;
    if constexpr (Alg == Algorithm::SmartMedian) {
      const int num_th = vsapi->mapNumElements(in, "threshold");
      for (int i = 0; i < num_th; ++i) {
        int th_err = 0;
        const double th_val = vsapi->mapGetFloat(in, "threshold", i, &th_err);
        require(!th_err, "SmartMedian: invalid threshold parameter.");
        threshold_list.push_back(static_cast<float>(th_val));
      }
      int sc_err = 0;
      const std::int64_t sc_val = vsapi->mapGetInt(in, "scalep", 0, &sc_err);
      if (!sc_err) {
        scalep = (sc_val != 0);
      }
    }

    std::vector<int> planes_list;
    const int num_planes = vsapi->mapNumElements(in, "planes");
    const bool planes_specified = (num_planes >= 0);
    for (int i = 0; i < num_planes; ++i) {
      planes_list.push_back(read_int(in, "planes", i, vsapi));
    }

    FilterPlan plan = build_plan(Alg, fmt, param_list, threshold_list, scalep, planes_list, planes_specified);

    auto instance = std::make_unique<VSFilterInstance>();
    instance->node = node;
    instance->vi = *vi;
    instance->plan = plan;

    VSFilterDependency dep{node, rpStrictSpatial};
    vsapi->createVideoFilter(out, algorithm_name(Alg), vi, filter_get_frame, filter_free, fmParallel, &dep, 1,
                             instance.release(), core);
  } catch (const std::exception& e) {
    if (node) {
      vsapi->freeNode(node);
    }
    vsapi->mapSetError(out, e.what());
  }
}

// ---------------------------------------------------------------------------
// Repair (clip, repairclip, mode)
// ---------------------------------------------------------------------------
const VSFrame* VS_CC repair_get_frame(int n, int activation_reason, void* instance_data, void**,
                                      VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSRepairInstance*>(instance_data);

  try {
    if (activation_reason == arInitial) {
      vsapi->requestFrameFilter(n, d->node, frame_ctx);
      vsapi->requestFrameFilter(n, d->repair_node, frame_ctx);
    } else if (activation_reason == arAllFramesReady) {
      const auto free_frame = [vsapi](const VSFrame* f) {
        vsapi->freeFrame(f);
      };
      std::unique_ptr<const VSFrame, decltype(free_frame)> src_owner(vsapi->getFrameFilter(n, d->node, frame_ctx),
                                                                     free_frame);
      std::unique_ptr<const VSFrame, decltype(free_frame)> repair_owner(vsapi->getFrameFilter(n, d->repair_node, frame_ctx),
                                                                        free_frame);
      const VSFrame* src = src_owner.get();
      const VSFrame* repair = repair_owner.get();
      if (!src || !repair) {
        return nullptr;
      }

      const VSFrame* plane_src[3] = {
          d->plan.process[0] ? nullptr : src,
          d->plan.process[1] ? nullptr : src,
          d->plan.process[2] ? nullptr : src,
      };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, src, core);

      if (!dst)
        throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      for (int plane = 0; plane < d->vi.format.numPlanes; ++plane) {
        if (!d->plan.process[static_cast<std::size_t>(plane)]) {
          continue;
        }
        const auto width = static_cast<std::size_t>(vsapi->getFrameWidth(dst, plane));
        const auto height = static_cast<std::size_t>(vsapi->getFrameHeight(dst, plane));
        const auto src_stride = static_cast<std::size_t>(vsapi->getStride(src, plane));
        const auto repair_stride = static_cast<std::size_t>(vsapi->getStride(repair, plane));
        const auto dst_stride = static_cast<std::size_t>(vsapi->getStride(dst, plane));
        const auto* srcp = vsapi->getReadPtr(src, plane);
        const auto* repairp = vsapi->getReadPtr(repair, plane);
        auto* dstp = vsapi->getWritePtr(dst, plane);

        execute_repair_plane(d->plan, plane, srcp, repairp, dstp, width, height, src_stride, repair_stride, dst_stride);
      }

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("Repair: frame processing failed", frame_ctx);
  }
  return nullptr;
}

void VS_CC repair_free(void* instance_data, VSCore*, const VSAPI* vsapi) {
  auto* d = static_cast<VSRepairInstance*>(instance_data);
  if (d) {
    if (d->node) {
      vsapi->freeNode(d->node);
    }
    if (d->repair_node) {
      vsapi->freeNode(d->repair_node);
    }
    delete d;
  }
}

void VS_CC repair_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* node = nullptr;
  VSNode* repair_node = nullptr;
  try {
    int err = 0;
    node = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && node != nullptr, "Repair: clip is required.");

    repair_node = vsapi->mapGetNode(in, "repairclip", 0, &err);
    require(!err && repair_node != nullptr, "Repair: repairclip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(node);
    const VSVideoInfo* rvi = vsapi->getVideoInfo(repair_node);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            "Repair: only constant format and dimensions supported.");
    require(is_same_video_info(vi, rvi), "Repair: Input clips must have the same format.");

    FormatInfo fmt = extract_format_info(vi);

    std::vector<int> param_list;
    const int num_params = vsapi->mapNumElements(in, "mode");
    for (int i = 0; i < num_params; ++i) {
      param_list.push_back(read_int(in, "mode", i, vsapi));
    }

    FilterPlan plan = build_plan(Algorithm::Repair, fmt, param_list, {}, false);

    auto instance = std::make_unique<VSRepairInstance>();
    instance->node = node;
    instance->repair_node = repair_node;
    instance->vi = *vi;
    instance->plan = plan;

    VSFilterDependency deps[2] = {
        {node, rpStrictSpatial},
        {repair_node, rvi->numFrames >= vi->numFrames ? rpStrictSpatial : rpGeneral},
    };
    vsapi->createVideoFilter(out, "Repair", vi, repair_get_frame, repair_free, fmParallel, deps, 2,
                             instance.release(), core);
  } catch (const std::exception& e) {
    if (node) {
      vsapi->freeNode(node);
    }
    if (repair_node) {
      vsapi->freeNode(repair_node);
    }
    vsapi->mapSetError(out, e.what());
  }
}

// ---------------------------------------------------------------------------
// Clense family (Clense, ForwardClense, BackwardClense)
// ---------------------------------------------------------------------------
template <Algorithm Alg>
const VSFrame* VS_CC clense_get_frame(int n, int activation_reason, void* instance_data, void** frame_data,
                                      VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSClenseInstance*>(instance_data);

  try {
    if (activation_reason == arInitial) {
      if constexpr (Alg == Algorithm::Clense) {
        if (n >= 1 && n <= d->vi.numFrames - 2) {
          *frame_data = reinterpret_cast<void*>(1);
          vsapi->requestFrameFilter(n - 1, d->pnode, frame_ctx);
          vsapi->requestFrameFilter(n, d->cnode, frame_ctx);
          vsapi->requestFrameFilter(n + 1, d->nnode, frame_ctx);
        } else {
          *frame_data = nullptr;
          vsapi->requestFrameFilter(n, d->cnode, frame_ctx);
        }
      } else if constexpr (Alg == Algorithm::ForwardClense) {
        vsapi->requestFrameFilter(n, d->cnode, frame_ctx);
        if (n <= d->vi.numFrames - 3) {
          *frame_data = reinterpret_cast<void*>(1);
          vsapi->requestFrameFilter(n + 1, d->cnode, frame_ctx);
          vsapi->requestFrameFilter(n + 2, d->cnode, frame_ctx);
        } else {
          *frame_data = nullptr;
        }
      } else if constexpr (Alg == Algorithm::BackwardClense) {
        if (n >= 2) {
          *frame_data = reinterpret_cast<void*>(1);
          vsapi->requestFrameFilter(n - 2, d->cnode, frame_ctx);
          vsapi->requestFrameFilter(n - 1, d->cnode, frame_ctx);
        } else {
          *frame_data = nullptr;
        }
        vsapi->requestFrameFilter(n, d->cnode, frame_ctx);
      }
    } else if (activation_reason == arAllFramesReady) {
      if (*frame_data != reinterpret_cast<void*>(1)) {
        return vsapi->getFrameFilter(n, d->cnode, frame_ctx);
      }

      const auto free_frame = [vsapi](const VSFrame* f) {
        vsapi->freeFrame(f);
      };

      const VSFrame* ref1 = nullptr;
      const VSFrame* src = nullptr;
      const VSFrame* ref2 = nullptr;

      if constexpr (Alg == Algorithm::Clense) {
        ref1 = vsapi->getFrameFilter(n - 1, d->pnode, frame_ctx);
        src = vsapi->getFrameFilter(n, d->cnode, frame_ctx);
        ref2 = vsapi->getFrameFilter(n + 1, d->nnode, frame_ctx);
      } else if constexpr (Alg == Algorithm::ForwardClense) {
        ref1 = vsapi->getFrameFilter(n + 1, d->cnode, frame_ctx);
        src = vsapi->getFrameFilter(n, d->cnode, frame_ctx);
        ref2 = vsapi->getFrameFilter(n + 2, d->cnode, frame_ctx);
      } else if constexpr (Alg == Algorithm::BackwardClense) {
        ref1 = vsapi->getFrameFilter(n - 1, d->cnode, frame_ctx);
        src = vsapi->getFrameFilter(n, d->cnode, frame_ctx);
        ref2 = vsapi->getFrameFilter(n - 2, d->cnode, frame_ctx);
      }

      std::unique_ptr<const VSFrame, decltype(free_frame)> ref1_owner(ref1, free_frame);
      std::unique_ptr<const VSFrame, decltype(free_frame)> src_owner(src, free_frame);
      std::unique_ptr<const VSFrame, decltype(free_frame)> ref2_owner(ref2, free_frame);

      if (!ref1 || !src || !ref2) {
        return nullptr;
      }

      const VSFrame* plane_src[3] = {
          d->plan.process[0] ? nullptr : src,
          d->plan.process[1] ? nullptr : src,
          d->plan.process[2] ? nullptr : src,
      };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, src, core);

      if (!dst)
        throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      for (int plane = 0; plane < d->vi.format.numPlanes; ++plane) {
        if (!d->plan.process[static_cast<std::size_t>(plane)]) {
          continue;
        }
        const auto width = static_cast<std::size_t>(vsapi->getFrameWidth(dst, plane));
        const auto height = static_cast<std::size_t>(vsapi->getFrameHeight(dst, plane));
        const auto src_stride = static_cast<std::size_t>(vsapi->getStride(src, plane));
        const auto ref1_stride = static_cast<std::size_t>(vsapi->getStride(ref1, plane));
        const auto ref2_stride = static_cast<std::size_t>(vsapi->getStride(ref2, plane));
        const auto dst_stride = static_cast<std::size_t>(vsapi->getStride(dst, plane));
        const auto* srcp = vsapi->getReadPtr(src, plane);
        const auto* ref1p = vsapi->getReadPtr(ref1, plane);
        const auto* ref2p = vsapi->getReadPtr(ref2, plane);
        auto* dstp = vsapi->getWritePtr(dst, plane);

        if constexpr (Alg == Algorithm::Clense) {
          execute_clense_plane(d->plan, srcp, ref1p, ref2p, dstp, width, height, src_stride, ref1_stride, ref2_stride, dst_stride);
        } else {
          execute_clense_forward_backward_plane(d->plan, srcp, ref1p, ref2p, dstp, width, height, src_stride, ref1_stride, ref2_stride, dst_stride);
        }
      }

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("Clense: frame processing failed", frame_ctx);
  }
  return nullptr;
}

void VS_CC clense_free(void* instance_data, VSCore*, const VSAPI* vsapi) {
  auto* d = static_cast<VSClenseInstance*>(instance_data);
  if (d) {
    if (d->cnode) {
      vsapi->freeNode(d->cnode);
    }
    if (d->pnode) {
      vsapi->freeNode(d->pnode);
    }
    if (d->nnode) {
      vsapi->freeNode(d->nnode);
    }
    delete d;
  }
}

template <Algorithm Alg>
void VS_CC clense_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* cnode = nullptr;
  VSNode* pnode = nullptr;
  VSNode* nnode = nullptr;
  try {
    int err = 0;
    cnode = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && cnode != nullptr, "Clense: clip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(cnode);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            "Clense: only constant format input supported");

    if constexpr (Alg == Algorithm::Clense) {
      pnode = vsapi->mapGetNode(in, "previous", 0, &err);
      if (err || !pnode) {
        pnode = vsapi->addNodeRef(cnode);
      } else {
        const VSVideoInfo* pvi = vsapi->getVideoInfo(pnode);
        require(is_same_video_info(vi, pvi), "Clense: previous clip does not have the same format as the main clip.");
      }

      nnode = vsapi->mapGetNode(in, "next", 0, &err);
      if (err || !nnode) {
        nnode = vsapi->addNodeRef(cnode);
      } else {
        const VSVideoInfo* nvi = vsapi->getVideoInfo(nnode);
        require(is_same_video_info(vi, nvi), "Clense: next clip does not have the same format as the main clip.");
      }
    }

    FormatInfo fmt = extract_format_info(vi);

    std::vector<int> planes_list;
    const int num_planes = vsapi->mapNumElements(in, "planes");
    const bool planes_specified = (num_planes >= 0);
    for (int i = 0; i < num_planes; ++i) {
      planes_list.push_back(read_int(in, "planes", i, vsapi));
    }

    FilterPlan plan = build_plan(Alg, fmt, {}, planes_list, planes_specified);

    auto instance = std::make_unique<VSClenseInstance>();
    instance->cnode = cnode;
    instance->pnode = pnode;
    instance->nnode = nnode;
    instance->vi = *vi;
    instance->plan = plan;

    if constexpr (Alg == Algorithm::Clense) {
      // Shared temporal inputs and short references reuse frames.
      VSFilterDependency deps[3]{};
      int num_deps = 0;
      for (VSNode* source : {cnode, pnode, nnode}) {
        int i = 0;
        while (i < num_deps && deps[i].source != source)
          ++i;
        if (i < num_deps) {
          deps[i].requestPattern = rpGeneral;
        } else {
          const int pattern = source == cnode ? rpStrictSpatial :
              (vsapi->getVideoInfo(source)->numFrames >= vi->numFrames ? rpNoFrameReuse : rpGeneral);
          deps[num_deps++] = {source, pattern};
        }
      }
      vsapi->createVideoFilter(out, "Clense", vi, clense_get_frame<Alg>, clense_free, fmParallel, deps, num_deps,
                               instance.release(), core);
    } else {
      VSFilterDependency dep{cnode, rpGeneral};
      vsapi->createVideoFilter(out, algorithm_name(Alg), vi, clense_get_frame<Alg>, clense_free, fmParallel, &dep, 1,
                               instance.release(), core);
    }
  } catch (const std::exception& e) {
    if (cnode) vsapi->freeNode(cnode);
    if (pnode) vsapi->freeNode(pnode);
    if (nnode) vsapi->freeNode(nnode);
    vsapi->mapSetError(out, e.what());
  }
}

// ---------------------------------------------------------------------------
// Temporal Filters (TemporalMedian, TemporalSoften)
// ---------------------------------------------------------------------------
template <Algorithm Alg>
const VSFrame* VS_CC temporal_filter_get_frame(int n, int activation_reason, void* instance_data, void**,
                                               VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSFilterInstance*>(instance_data);
  const int radius = d->plan.params[0];
  const bool passthrough = Alg == Algorithm::TemporalMedian && (n < radius || n > d->vi.numFrames - 1 - radius);
  const int first = std::max(0, n - radius);
  const int last = n + std::min(radius, d->vi.numFrames - 1 - n);
  const int center = n - first;
  const int diameter = last - first + 1;

  try {
    if (activation_reason == arInitial) {
      if (passthrough) {
        vsapi->requestFrameFilter(n, d->node, frame_ctx);
      } else {
        for (int i = first; i <= last; ++i) {
          vsapi->requestFrameFilter(i, d->node, frame_ctx);
        }
      }
    } else if (activation_reason == arAllFramesReady) {
      if (passthrough) {
        return vsapi->getFrameFilter(n, d->node, frame_ctx);
      }

      const auto free_frame = [vsapi](const VSFrame* f) {
        vsapi->freeFrame(f);
      };

      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> frame_owners;
      frame_owners.reserve(static_cast<std::size_t>(diameter));
      std::vector<const VSFrame*> raw_frames(static_cast<std::size_t>(diameter));
      for (int i = first; i <= last; ++i) {
        const VSFrame* f = vsapi->getFrameFilter(i, d->node, frame_ctx);
        frame_owners.emplace_back(f, free_frame);
        if (!f) return nullptr;
        raw_frames[static_cast<std::size_t>(i - first)] = f;
      }

      auto window = resolve_temporal_window(raw_frames.data(), static_cast<std::size_t>(center),
                                            static_cast<std::size_t>(diameter), d->plan.scenechange, vsapi,
                                            algorithm_name(Alg), Alg == Algorithm::TemporalMedian);
      const int active_diameter = static_cast<int>(window.to_idx - window.from_idx + 1);

      const VSFrame* center_frame = raw_frames[static_cast<std::size_t>(center)];
      const VSFrame* plane_src[3] = {
          d->plan.process[0] ? nullptr : center_frame,
          d->plan.process[1] ? nullptr : center_frame,
          d->plan.process[2] ? nullptr : center_frame,
      };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, center_frame, core);
      if (!dst) throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      const DataType dtype = get_data_type(d->plan.format.bytes_per_sample, d->plan.format.is_float);

      for (int plane = 0; plane < d->vi.format.numPlanes; ++plane) {
        if (!d->plan.process[static_cast<std::size_t>(plane)]) continue;

        const auto width = static_cast<std::size_t>(vsapi->getFrameWidth(dst, plane));
        const auto height = static_cast<std::size_t>(vsapi->getFrameHeight(dst, plane));
        const auto src_stride = static_cast<std::size_t>(vsapi->getStride(center_frame, plane));
        const auto dst_stride = static_cast<std::size_t>(vsapi->getStride(dst, plane));
        auto* dstp = vsapi->getWritePtr(dst, plane);

        std::vector<const std::uint8_t*> plane_ptrs(static_cast<std::size_t>(active_diameter));
        if constexpr (Alg == Algorithm::TemporalSoften) {
          plane_ptrs[0] = vsapi->getReadPtr(center_frame, plane);
          std::size_t idx = 1;
          // Match upstream accumulation order: current, nearest previous first, then next.
          for (std::size_t k = static_cast<std::size_t>(center); k > window.from_idx; --k)
            plane_ptrs[idx++] = vsapi->getReadPtr(raw_frames[k - 1], plane);
          for (std::size_t k = static_cast<std::size_t>(center) + 1; k <= window.to_idx; ++k)
            plane_ptrs[idx++] = vsapi->getReadPtr(raw_frames[k], plane);
          process_temporal_soften_plane(dtype, active_diameter, d->plan.thresholds[static_cast<std::size_t>(plane)],
                                       plane_ptrs.data(), dstp, width, height, src_stride, dst_stride);
        } else {
          for (int k = 0; k < active_diameter; ++k) {
            plane_ptrs[static_cast<std::size_t>(k)] = vsapi->getReadPtr(raw_frames[window.from_idx + k], plane);
          }
          process_temporal_median_plane(dtype, active_diameter, plane_ptrs.data(), dstp, width, height, src_stride, dst_stride);
        }
      }

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("neo-smo: frame processing failed", frame_ctx);
  }
  return nullptr;
}

template <Algorithm Alg>
void VS_CC temporal_filter_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* node = nullptr;
  try {
    int err = 0;
    node = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && node != nullptr, std::string(algorithm_name(Alg)) + ": clip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(node);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            std::string(algorithm_name(Alg)) + ": only constant format input supported");

    FormatInfo fmt = extract_format_info(vi);

    int radius = (Alg == Algorithm::TemporalSoften) ? 4 : 1;
    const int num_r = vsapi->mapNumElements(in, "radius");
    if (num_r > 0) {
      radius = read_int(in, "radius", 0, vsapi);
    }

    std::vector<float> threshold_list;
    if constexpr (Alg == Algorithm::TemporalSoften) {
      const int num_th = vsapi->mapNumElements(in, "threshold");
      for (int i = 0; i < num_th; ++i) {
        int th_err = 0;
        const double th_val = vsapi->mapGetFloat(in, "threshold", i, &th_err);
        require(!th_err && std::isfinite(static_cast<float>(th_val)), "TemporalSoften: invalid threshold parameter.");
        threshold_list.push_back(static_cast<float>(th_val));
      }
    }

    bool scenechange = false;
    if constexpr (Alg == Algorithm::TemporalSoften) {
      int sc_err = 0;
      const std::int64_t sc_val = vsapi->mapGetInt(in, "scenechange", 0, &sc_err);
      if (!sc_err) {
        require(sc_val >= -1 && sc_val <= 254, "TemporalSoften: scenechange must be between -1 and 254 (inclusive)");
        scenechange = (sc_val != 0);
        if (sc_val > 0) {
          require(vi->format.colorFamily != cfRGB, "TemporalSoften: Scene change support does not work with RGB.");
          VSPlugin* misc = vsapi->getPluginByID("com.vapoursynth.misc", core);
          require(misc != nullptr, "TemporalSoften: Miscellaneous filters (https://github.com/vapoursynth/vs-miscfilters-obsolete) plugin is required in order to use scene change detection.");
          const auto free_map = [vsapi](VSMap* map) { vsapi->freeMap(map); };
          std::unique_ptr<VSMap, decltype(free_map)> args(vsapi->createMap(), free_map);
          if (!args) throw std::bad_alloc();
          vsapi->mapSetNode(args.get(), "clip", node, maReplace);
          vsapi->mapSetFloat(args.get(), "threshold", static_cast<double>(sc_val) / 255.0, maReplace);
          std::unique_ptr<VSMap, decltype(free_map)> ret(vsapi->invoke(misc, "SCDetect", args.get()), free_map);
          if (!ret) throw std::bad_alloc();
          const char* ret_msg = vsapi->mapGetError(ret.get());
          if (ret_msg) throw std::runtime_error(ret_msg);
          int ret_err = 0;
          VSNode* detected = vsapi->mapGetNode(ret.get(), "clip", 0, &ret_err);
          if (ret_err || !detected) {
            if (detected) vsapi->freeNode(detected);
            throw std::runtime_error("TemporalSoften: Unexpected error while invoking SCDetect.");
          }
          vsapi->freeNode(node);
          node = detected;
          vi = vsapi->getVideoInfo(node);
        }
      }
    } else {
      int sc_err = 0;
      const std::int64_t sc_val = vsapi->mapGetInt(in, "scenechange", 0, &sc_err);
      if (!sc_err) {
        scenechange = (sc_val != 0);
      }
    }

    std::vector<int> planes_list;
    const int num_planes = vsapi->mapNumElements(in, "planes");
    const bool planes_specified = (num_planes >= 0);
    for (int i = 0; i < num_planes; ++i) {
      planes_list.push_back(read_int(in, "planes", i, vsapi));
    }

    const bool scalep = vsapi->mapGetInt(in, "scalep", 0, &err) != 0 && !err;
    FilterPlan plan = build_plan(Alg, fmt, {radius}, threshold_list, scalep, planes_list, planes_specified);
    plan.scenechange = scenechange;

    auto instance = std::make_unique<VSFilterInstance>();
    instance->node = node;
    instance->vi = *vi;
    instance->plan = plan;

    VSFilterDependency dep{node, rpGeneral};
    vsapi->createVideoFilter(out, algorithm_name(Alg), vi, temporal_filter_get_frame<Alg>, filter_free, fmParallel, &dep, 1,
                             instance.release(), core);
  } catch (const std::exception& e) {
    if (node) vsapi->freeNode(node);
    vsapi->mapSetError(out, e.what());
  }
}

// ---------------------------------------------------------------------------
// TemporalRepair
// ---------------------------------------------------------------------------
const VSFrame* VS_CC temporal_repair_get_frame(int n, int activation_reason, void* instance_data, void**,
                                               VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSRepairInstance*>(instance_data);

  try {
    if (activation_reason == arInitial) {
      if (n == 0 || n == d->vi.numFrames - 1) {
        vsapi->requestFrameFilter(n, d->node, frame_ctx);
      } else {
        vsapi->requestFrameFilter(n, d->node, frame_ctx);
        vsapi->requestFrameFilter(n - 1, d->repair_node, frame_ctx);
        vsapi->requestFrameFilter(n, d->repair_node, frame_ctx);
        vsapi->requestFrameFilter(n + 1, d->repair_node, frame_ctx);
      }
    } else if (activation_reason == arAllFramesReady) {
      if (n == 0 || n == d->vi.numFrames - 1) {
        return vsapi->getFrameFilter(n, d->node, frame_ctx);
      }

      const auto free_frame = [vsapi](const VSFrame* f) {
        vsapi->freeFrame(f);
      };

      std::unique_ptr<const VSFrame, decltype(free_frame)> src_owner(vsapi->getFrameFilter(n, d->node, frame_ctx), free_frame);
      std::unique_ptr<const VSFrame, decltype(free_frame)> prev_owner(vsapi->getFrameFilter(n - 1, d->repair_node, frame_ctx), free_frame);
      std::unique_ptr<const VSFrame, decltype(free_frame)> curr_owner(vsapi->getFrameFilter(n, d->repair_node, frame_ctx), free_frame);
      std::unique_ptr<const VSFrame, decltype(free_frame)> next_owner(vsapi->getFrameFilter(n + 1, d->repair_node, frame_ctx), free_frame);

      const VSFrame* src = src_owner.get();
      const VSFrame* prev = prev_owner.get();
      const VSFrame* curr = curr_owner.get();
      const VSFrame* next = next_owner.get();
      if (!src || !prev || !curr || !next) return nullptr;

      const VSFrame* plane_src[3] = {
          d->plan.process[0] ? nullptr : src,
          d->plan.process[1] ? nullptr : src,
          d->plan.process[2] ? nullptr : src,
      };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, src, core);
      if (!dst) throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      const DataType dtype = get_data_type(d->plan.format.bytes_per_sample, d->plan.format.is_float);

      for (int plane = 0; plane < d->vi.format.numPlanes; ++plane) {
        if (!d->plan.process[static_cast<std::size_t>(plane)]) continue;

        const auto width = static_cast<std::size_t>(vsapi->getFrameWidth(dst, plane));
        const auto height = static_cast<std::size_t>(vsapi->getFrameHeight(dst, plane));
        const auto src_stride = static_cast<std::size_t>(vsapi->getStride(src, plane));
        const auto prev_stride = static_cast<std::size_t>(vsapi->getStride(prev, plane));
        const auto curr_stride = static_cast<std::size_t>(vsapi->getStride(curr, plane));
        const auto next_stride = static_cast<std::size_t>(vsapi->getStride(next, plane));
        const auto dst_stride = static_cast<std::size_t>(vsapi->getStride(dst, plane));

        const auto* srcp = vsapi->getReadPtr(src, plane);
        const auto* prevp = vsapi->getReadPtr(prev, plane);
        const auto* currp = vsapi->getReadPtr(curr, plane);
        const auto* nextp = vsapi->getReadPtr(next, plane);
        auto* dstp = vsapi->getWritePtr(dst, plane);

        const bool chroma = (d->plan.format.color_family == 3) && (plane > 0);
        const int mode = d->plan.params[static_cast<std::size_t>(plane)];

        process_temporal_repair_plane(dtype, mode, chroma, d->plan.format.bits_per_sample,
                                     srcp, prevp, currp, nextp, dstp, width, height,
                                     src_stride, prev_stride, curr_stride, next_stride, dst_stride);
      }

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("TemporalRepair: frame processing failed", frame_ctx);
  }
  return nullptr;
}

void VS_CC temporal_repair_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* node = nullptr;
  VSNode* repair_node = nullptr;
  try {
    int err = 0;
    node = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && node != nullptr, "TemporalRepair: clip is required.");

    repair_node = vsapi->mapGetNode(in, "repairclip", 0, &err);
    require(!err && repair_node != nullptr, "TemporalRepair: repairclip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(node);
    const VSVideoInfo* rvi = vsapi->getVideoInfo(repair_node);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            "TemporalRepair: only constant format and dimensions supported.");
    require(is_same_video_info(vi, rvi), "TemporalRepair: Input clips must have the same format.");

    FormatInfo fmt = extract_format_info(vi);

    std::vector<int> mode_list;
    const int num_m = vsapi->mapNumElements(in, "mode");
    for (int i = 0; i < num_m; ++i) {
      mode_list.push_back(read_int(in, "mode", i, vsapi));
    }

    std::vector<int> planes_list;
    const int num_planes = vsapi->mapNumElements(in, "planes");
    const bool planes_specified = (num_planes >= 0);
    for (int i = 0; i < num_planes; ++i) {
      planes_list.push_back(read_int(in, "planes", i, vsapi));
    }

    FilterPlan plan = build_plan(Algorithm::TemporalRepair, fmt, mode_list, {}, false, planes_list, planes_specified);

    auto instance = std::make_unique<VSRepairInstance>();
    instance->node = node;
    instance->repair_node = repair_node;
    instance->vi = *vi;
    instance->plan = plan;

    VSFilterDependency deps[2] = {
        {node, rpStrictSpatial},
        {repair_node, rpGeneral},
    };
    const int dependency_count = node == repair_node ? 1 : 2;
    if (dependency_count == 1) deps[0].requestPattern = rpGeneral;
    vsapi->createVideoFilter(out, "TemporalRepair", vi, temporal_repair_get_frame, repair_free, fmParallel, deps, dependency_count,
                             instance.release(), core);
  } catch (const std::exception& e) {
    if (node) vsapi->freeNode(node);
    if (repair_node) vsapi->freeNode(repair_node);
    vsapi->mapSetError(out, e.what());
  }
}

// ---------------------------------------------------------------------------
// Trio filters (DegrainMedian, FluxSmoothT, FluxSmoothST)
// ---------------------------------------------------------------------------
template <Algorithm Alg>
const VSFrame* VS_CC trio_filter_get_frame(int n, int activation_reason, void* instance_data, void**,
                                           VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSFilterInstance*>(instance_data);

  try {
    if (activation_reason == arInitial) {
      if (n == 0 || n == d->vi.numFrames - 1) {
        vsapi->requestFrameFilter(n, d->node, frame_ctx);
      } else {
        vsapi->requestFrameFilter(n - 1, d->node, frame_ctx);
        vsapi->requestFrameFilter(n, d->node, frame_ctx);
        vsapi->requestFrameFilter(n + 1, d->node, frame_ctx);
      }
    } else if (activation_reason == arAllFramesReady) {
      if (n == 0 || n == d->vi.numFrames - 1) {
        return vsapi->getFrameFilter(n, d->node, frame_ctx);
      }

      const auto free_frame = [vsapi](const VSFrame* f) {
        vsapi->freeFrame(f);
      };

      std::unique_ptr<const VSFrame, decltype(free_frame)> prev_owner(vsapi->getFrameFilter(n - 1, d->node, frame_ctx), free_frame);
      std::unique_ptr<const VSFrame, decltype(free_frame)> curr_owner(vsapi->getFrameFilter(n, d->node, frame_ctx), free_frame);
      std::unique_ptr<const VSFrame, decltype(free_frame)> next_owner(vsapi->getFrameFilter(n + 1, d->node, frame_ctx), free_frame);

      const VSFrame* prev = prev_owner.get();
      const VSFrame* curr = curr_owner.get();
      const VSFrame* next = next_owner.get();
      if (!prev || !curr || !next) return nullptr;

      const VSFrame* plane_src[3] = {
          d->plan.process[0] ? nullptr : curr,
          d->plan.process[1] ? nullptr : curr,
          d->plan.process[2] ? nullptr : curr,
      };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, curr, core);
      if (!dst) throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      const DataType dtype = get_data_type(d->plan.format.bytes_per_sample, d->plan.format.is_float);

      for (int plane = 0; plane < d->vi.format.numPlanes; ++plane) {
        if (!d->plan.process[static_cast<std::size_t>(plane)]) continue;

        const auto width = static_cast<std::size_t>(vsapi->getFrameWidth(dst, plane));
        const auto height = static_cast<std::size_t>(vsapi->getFrameHeight(dst, plane));
        const auto prev_stride = static_cast<std::size_t>(vsapi->getStride(prev, plane));
        const auto curr_stride = static_cast<std::size_t>(vsapi->getStride(curr, plane));
        const auto next_stride = static_cast<std::size_t>(vsapi->getStride(next, plane));
        const auto dst_stride = static_cast<std::size_t>(vsapi->getStride(dst, plane));

        const auto* prevp = vsapi->getReadPtr(prev, plane);
        const auto* currp = vsapi->getReadPtr(curr, plane);
        const auto* nextp = vsapi->getReadPtr(next, plane);
        auto* dstp = vsapi->getWritePtr(dst, plane);

        const bool chroma = (d->plan.format.color_family == 3) && (plane > 0);

        if constexpr (Alg == Algorithm::DegrainMedian) {
          const int mode = d->plan.params[static_cast<std::size_t>(plane)];
          const float limit = d->plan.thresholds[static_cast<std::size_t>(plane)];
          process_degrain_median_plane(dtype, mode, limit, d->plan.interlaced, d->plan.norow, chroma,
                                       d->plan.format.bits_per_sample, prevp, currp, nextp, dstp,
                                       width, height, prev_stride, curr_stride, next_stride, dst_stride);
        } else if constexpr (Alg == Algorithm::FluxSmoothT) {
          const float t_th = d->plan.thresholds[static_cast<std::size_t>(plane)];
          process_fluxsmooth_t_plane(dtype, t_th, prevp, currp, nextp, dstp, width, height,
                                    prev_stride, curr_stride, next_stride, dst_stride);
        } else if constexpr (Alg == Algorithm::FluxSmoothST) {
          const float t_th = d->plan.thresholds[static_cast<std::size_t>(plane)];
          const float s_th = d->plan.secondary_thresholds[static_cast<std::size_t>(plane)];
          process_fluxsmooth_st_plane(dtype, t_th, s_th, prevp, currp, nextp, dstp, width, height,
                                     prev_stride, curr_stride, next_stride, dst_stride);
        }
      }

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("neo-smo: frame processing failed", frame_ctx);
  }
  return nullptr;
}

template <Algorithm Alg>
void VS_CC trio_filter_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* node = nullptr;
  try {
    int err = 0;
    node = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && node != nullptr, std::string(algorithm_name(Alg)) + ": clip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(node);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            std::string(algorithm_name(Alg)) + ": only constant format input supported");

    FormatInfo fmt = extract_format_info(vi);

    std::vector<int> mode_list;
    std::vector<float> th_list;
    std::vector<float> sec_th_list;
    bool interlaced = false;
    bool norow = false;

    if constexpr (Alg == Algorithm::DegrainMedian) {
      const int num_m = vsapi->mapNumElements(in, "mode");
      for (int i = 0; i < num_m; ++i) {
        mode_list.push_back(read_int(in, "mode", i, vsapi));
      }
      const int num_l = vsapi->mapNumElements(in, "limit");
      for (int i = 0; i < num_l; ++i) {
        int l_err = 0;
        const double l_val = vsapi->mapGetFloat(in, "limit", i, &l_err);
        require(!l_err && std::isfinite(static_cast<float>(l_val)), "DegrainMedian: invalid limit parameter.");
        th_list.push_back(static_cast<float>(l_val));
      }
      int int_err = 0;
      const std::int64_t int_val = vsapi->mapGetInt(in, "interlaced", 0, &int_err);
      if (!int_err) interlaced = (int_val != 0);

      int nr_err = 0;
      const std::int64_t nr_val = vsapi->mapGetInt(in, "norow", 0, &nr_err);
      if (!nr_err) norow = (nr_val != 0);
    } else if constexpr (Alg == Algorithm::FluxSmoothT) {
      const int num_th = vsapi->mapNumElements(in, "temporal_threshold");
      for (int i = 0; i < num_th; ++i) {
        int th_err = 0;
        const double th_val = vsapi->mapGetFloat(in, "temporal_threshold", i, &th_err);
        require(!th_err && std::isfinite(static_cast<float>(th_val)), "FluxSmoothT: invalid temporal_threshold parameter.");
        th_list.push_back(static_cast<float>(th_val));
      }
    } else if constexpr (Alg == Algorithm::FluxSmoothST) {
      const int num_tth = vsapi->mapNumElements(in, "temporal_threshold");
      for (int i = 0; i < num_tth; ++i) {
        int th_err = 0;
        const double th_val = vsapi->mapGetFloat(in, "temporal_threshold", i, &th_err);
        require(!th_err && std::isfinite(static_cast<float>(th_val)), "FluxSmoothST: invalid temporal_threshold parameter.");
        th_list.push_back(static_cast<float>(th_val));
      }
      const int num_sth = vsapi->mapNumElements(in, "spatial_threshold");
      for (int i = 0; i < num_sth; ++i) {
        int th_err = 0;
        const double th_val = vsapi->mapGetFloat(in, "spatial_threshold", i, &th_err);
        require(!th_err && std::isfinite(static_cast<float>(th_val)), "FluxSmoothST: invalid spatial_threshold parameter.");
        sec_th_list.push_back(static_cast<float>(th_val));
      }
    }

    std::vector<int> planes_list;
    const int num_planes = vsapi->mapNumElements(in, "planes");
    const bool planes_specified = (num_planes >= 0);
    for (int i = 0; i < num_planes; ++i) {
      planes_list.push_back(read_int(in, "planes", i, vsapi));
    }

    const bool scalep = vsapi->mapGetInt(in, "scalep", 0, &err) != 0 && !err;
    FilterPlan plan = build_plan(Alg, fmt, mode_list, th_list, scalep, planes_list, planes_specified);
    plan.interlaced = interlaced;
    plan.norow = norow;
    if constexpr (Alg == Algorithm::FluxSmoothST) {
      plan.secondary_thresholds = flux_thresholds(fmt, sec_th_list, scalep);
      for (std::size_t i = 0; i < 3; ++i)
        plan.process[i] = plan.process[i] && (plan.thresholds[i] >= 0 || plan.secondary_thresholds[i] >= 0);
    }

    auto instance = std::make_unique<VSFilterInstance>();
    instance->node = node;
    instance->vi = *vi;
    instance->plan = plan;

    VSFilterDependency dep{node, rpGeneral};
    vsapi->createVideoFilter(out, algorithm_name(Alg), vi, trio_filter_get_frame<Alg>, filter_free, fmParallel, &dep, 1,
                             instance.release(), core);
  } catch (const std::exception& e) {
    if (node) vsapi->freeNode(node);
    vsapi->mapSetError(out, e.what());
  }
}

// ---------------------------------------------------------------------------
// TTempSmooth
// ---------------------------------------------------------------------------
struct VSTTempSmoothInstance {
  VSNode* node = nullptr;
  VSNode* pf_node = nullptr;
  VSVideoInfo vi{};
  int maxr = 3;
  std::array<int, 3> threshold{4, 5, 5};
  std::array<int, 3> mdiff{2, 3, 3};
  int strength = 2;
  float scthresh = 12.0f;
  bool scenechange = true;
  bool fp = true;
  std::array<bool, 3> process{true, true, true};
  std::array<int, 3> weight_mode{0, 0, 0};
  std::array<float, 3> center_weights{0.0f, 0.0f, 0.0f};
  std::array<std::vector<float>, 3> temporal_weights{};
  std::array<std::vector<float>, 3> temporal_difference_weights{};
};

void VS_CC ttempsmooth_free(void* instance_data, VSCore*, const VSAPI* vsapi) {
  auto* d = static_cast<VSTTempSmoothInstance*>(instance_data);
  if (d) {
    if (d->node) vsapi->freeNode(d->node);
    if (d->pf_node) vsapi->freeNode(d->pf_node);
    delete d;
  }
}

const VSFrame* VS_CC ttempsmooth_get_frame(int n, int activation_reason, void* instance_data, void**,
                                           VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSTTempSmoothInstance*>(instance_data);
  try {
    const int first = std::max(0, n - d->maxr);
    const int last = std::min(d->vi.numFrames - 1, n + d->maxr);
    const bool has_ref = (d->pf_node != nullptr);

    if (activation_reason == arInitial) {
      for (int i = first; i <= last; ++i) {
        vsapi->requestFrameFilter(i, d->node, frame_ctx);
        if (has_ref) vsapi->requestFrameFilter(i, d->pf_node, frame_ctx);
      }
    } else if (activation_reason == arAllFramesReady) {
      const auto free_frame = [vsapi](const VSFrame* f) { vsapi->freeFrame(f); };
      const int diameter = d->maxr * 2 + 1;

      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> src_frames;
      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> ref_frames;
      src_frames.reserve(static_cast<std::size_t>(diameter));
      if (has_ref) ref_frames.reserve(static_cast<std::size_t>(diameter));

      for (int i = -d->maxr; i <= d->maxr; ++i) {
        const int frame_num = std::clamp(n + i, 0, d->vi.numFrames - 1);
        src_frames.emplace_back(vsapi->getFrameFilter(frame_num, d->node, frame_ctx), free_frame);
        if (!src_frames.back()) return nullptr;
        if (has_ref) {
          ref_frames.emplace_back(vsapi->getFrameFilter(frame_num, d->pf_node, frame_ctx), free_frame);
        if (!ref_frames.back()) return nullptr;
        }
      }

      int from_frame_idx = 0;
      int to_frame_idx = diameter - 1;
      if (d->scenechange) {
        const auto& frames = has_ref ? ref_frames : src_frames;
        for (int i = d->maxr; i > 0; --i) {
          int err = 0;
          const auto* props = vsapi->getFramePropertiesRO(frames[static_cast<std::size_t>(i)].get());
          if (vsapi->mapGetInt(props, "_SceneChangePrev", 0, &err) == 1 && !err) {
            from_frame_idx = i;
            break;
          }
        }
        for (int i = d->maxr; i < diameter - 1; ++i) {
          int err = 0;
          const auto* props = vsapi->getFramePropertiesRO(frames[static_cast<std::size_t>(i)].get());
          if (vsapi->mapGetInt(props, "_SceneChangeNext", 0, &err) == 1 && !err) {
            to_frame_idx = i;
            break;
          }
        }
      }

      const VSFrame* center_src = src_frames[static_cast<std::size_t>(d->maxr)].get();
      const VSFrame* plane_src[3] = {
          d->process[0] ? nullptr : center_src,
          d->process[1] ? nullptr : center_src,
          d->process[2] ? nullptr : center_src,
      };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, center_src, core);
      if (!dst) throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      const DataType dtype = get_data_type(d->vi.format.bytesPerSample, d->vi.format.sampleType == stFloat);

      for (int plane = 0; plane < d->vi.format.numPlanes; ++plane) {
        if (!d->process[static_cast<std::size_t>(plane)]) continue;

        const auto width = vsapi->getFrameWidth(dst, plane);
        const auto height = vsapi->getFrameHeight(dst, plane);
        const auto src_stride = static_cast<std::size_t>(vsapi->getStride(center_src, plane));
        const auto ref_stride = has_ref ? static_cast<std::size_t>(vsapi->getStride(ref_frames[static_cast<std::size_t>(d->maxr)].get(), plane)) : src_stride;
        const auto dst_stride = static_cast<std::size_t>(vsapi->getStride(dst, plane));

        const auto* curr_ptr = vsapi->getReadPtr(center_src, plane);
        const auto* curr_ref_ptr = has_ref ? vsapi->getReadPtr(ref_frames[static_cast<std::size_t>(d->maxr)].get(), plane) : curr_ptr;
        auto* dst_ptr = vsapi->getWritePtr(dst, plane);

        std::vector<const std::uint8_t*> prev_ptrs;
        std::vector<const std::uint8_t*> prev_ref_ptrs;
        std::vector<const std::uint8_t*> next_ptrs;
        std::vector<const std::uint8_t*> next_ref_ptrs;

        const int num_prev = d->maxr - from_frame_idx;
        const int num_next = to_frame_idx - d->maxr;

        for (int i = 0; i < num_prev; ++i) {
          prev_ptrs.push_back(vsapi->getReadPtr(src_frames[static_cast<std::size_t>(d->maxr - 1 - i)].get(), plane));
          prev_ref_ptrs.push_back(has_ref ? vsapi->getReadPtr(ref_frames[static_cast<std::size_t>(d->maxr - 1 - i)].get(), plane) : prev_ptrs.back());
        }
        for (int i = 0; i < num_next; ++i) {
          next_ptrs.push_back(vsapi->getReadPtr(src_frames[static_cast<std::size_t>(d->maxr + 1 + i)].get(), plane));
          next_ref_ptrs.push_back(has_ref ? vsapi->getReadPtr(ref_frames[static_cast<std::size_t>(d->maxr + 1 + i)].get(), plane) : next_ptrs.back());
        }

        const float threshold = scale_to_format(extract_format_info(&d->vi), static_cast<float>(d->threshold[static_cast<std::size_t>(plane)]));

        process_ttempsmooth_plane(
            dtype, width, height,
            src_stride, ref_stride, dst_stride,
            curr_ptr, curr_ref_ptr,
            prev_ptrs.data(), prev_ref_ptrs.data(),
            next_ptrs.data(), next_ref_ptrs.data(),
            dst_ptr,
            d->maxr, num_prev, num_next,
            threshold, d->fp, d->vi.format.bitsPerSample,
            d->weight_mode[static_cast<std::size_t>(plane)],
            d->center_weights[static_cast<std::size_t>(plane)],
            d->temporal_weights[static_cast<std::size_t>(plane)].data(),
            d->temporal_difference_weights[static_cast<std::size_t>(plane)].data());
      }

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("TTempSmooth: frame processing failed", frame_ctx);
  }
  return nullptr;
}

void VS_CC ttempsmooth_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* node = nullptr;
  VSNode* pf_node = nullptr;
  try {
    int err = 0;
    node = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && node != nullptr, "TTempSmooth: clip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(node);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            "TTempSmooth: only constant format input supported");
    require(!(vi->format.sampleType == stFloat && vi->format.bytesPerSample == 2),
            "TTempSmooth: unsupported format (f16 not supported)");

    int maxr = 3;
    const std::int64_t maxr_val = vsapi->mapGetInt(in, "maxr", 0, &err);
    if (!err) {
      require(maxr_val >= 1 && maxr_val <= 7, "TTempSmooth: maxr must be between 1 and 7 (inclusive)");
      maxr = static_cast<int>(maxr_val);
    }

    std::array<int, 3> threshold{4, 5, 5};
    const int num_thresh = vsapi->mapNumElements(in, "thresh");
    if (num_thresh > 0) {
      for (int plane = 0; plane < 3; ++plane) {
        if (plane < num_thresh) {
          const std::int64_t val = vsapi->mapGetInt(in, "thresh", plane, &err);
          require(!err && val >= 1 && val <= 256, "TTempSmooth: thresh must be between 1 and 256");
          threshold[static_cast<std::size_t>(plane)] = static_cast<int>(val);
        } else {
          threshold[static_cast<std::size_t>(plane)] = threshold[static_cast<std::size_t>(plane - 1)];
        }
      }
    }

    std::array<int, 3> mdiff{2, 3, 3};
    const int num_mdiff = vsapi->mapNumElements(in, "mdiff");
    if (num_mdiff > 0) {
      for (int plane = 0; plane < 3; ++plane) {
        if (plane < num_mdiff) {
          const std::int64_t val = vsapi->mapGetInt(in, "mdiff", plane, &err);
          require(!err && val >= 0 && val <= 255, "TTempSmooth: mdiff must be between 0 and 255");
          mdiff[static_cast<std::size_t>(plane)] = static_cast<int>(val);
        } else {
          mdiff[static_cast<std::size_t>(plane)] = mdiff[static_cast<std::size_t>(plane - 1)];
        }
      }
    }

    int strength = 2;
    const std::int64_t strength_val = vsapi->mapGetInt(in, "strength", 0, &err);
    if (!err) {
      require(strength_val >= 1 && strength_val <= 8, "TTempSmooth: strength must be between 1 and 8 (inclusive)");
      strength = static_cast<int>(strength_val);
    }

    float scthresh = 12.0f;
    const double sc_val = vsapi->mapGetFloat(in, "scthresh", 0, &err);
    if (!err) {
      require(sc_val >= -1.0 && sc_val <= 100.0, "TTempSmooth: scthresh must be between -1 and 100.0 (inclusive)");
      scthresh = static_cast<float>(sc_val);
    }

    bool fp = true;
    const std::int64_t fp_val = vsapi->mapGetInt(in, "fp", 0, &err);
    if (!err) {
      fp = (fp_val != 0);
    }

    pf_node = vsapi->mapGetNode(in, "pfclip", 0, &err);
    if (!err && pf_node) {
      const VSVideoInfo* refvi = vsapi->getVideoInfo(pf_node);
      require(is_same_video_info(vi, refvi), "pfclip must have same format and dimensions as the main clip");
    }

    std::vector<int> planes_list;
    const int num_planes = vsapi->mapNumElements(in, "planes");
    const bool planes_specified = (num_planes >= 0);
    for (int i = 0; i < num_planes; ++i) {
      planes_list.push_back(read_int(in, "planes", i, vsapi));
    }
    const auto process = normalize_planes(vi->format.numPlanes, planes_list, planes_specified, "TTempSmooth");

    if (scthresh > 0.0f) {
      require(vi->format.colorFamily != cfRGB,
              "TTempSmooth: scthresh > 0 does not work with RGB. "
              "Invoke SCDetect (or similar) yourself with an RGB->YUV converted clip, "
              "copy the properties (CopyFrameProps) to your input clip, "
              "and then invoke TTempSmooth with scthresh=-1 to use those properties.");
      VSPlugin* misc = vsapi->getPluginByID("com.vapoursynth.misc", core);
      require(misc != nullptr, "TTempSmooth: Miscellaneous filters (https://github.com/vapoursynth/vs-miscfilters-obsolete) plugin is required in order to use scene change detection.");

      const auto free_map = [vsapi](VSMap* map) { vsapi->freeMap(map); };
      std::unique_ptr<VSMap, decltype(free_map)> args(vsapi->createMap(), free_map);
      if (!args) throw std::bad_alloc();
      VSNode* target_node = pf_node ? pf_node : node;
      vsapi->mapSetNode(args.get(), "clip", target_node, maReplace);
      vsapi->mapSetFloat(args.get(), "threshold", static_cast<double>(scthresh) / 100.0, maReplace);
      std::unique_ptr<VSMap, decltype(free_map)> ret(vsapi->invoke(misc, "SCDetect", args.get()), free_map);
      if (!ret) throw std::bad_alloc();
      const char* ret_msg = vsapi->mapGetError(ret.get());
      if (ret_msg) throw std::runtime_error(ret_msg);
      int ret_err = 0;
      VSNode* detected = vsapi->mapGetNode(ret.get(), "clip", 0, &ret_err);
      require(!ret_err && detected != nullptr, "TTempSmooth: Unexpected error while invoking SCDetect");
      if (pf_node) {
        vsapi->freeNode(pf_node);
        pf_node = detected;
      } else {
        vsapi->freeNode(node);
        node = detected;
        vi = vsapi->getVideoInfo(node);
      }
    }

    auto instance = std::make_unique<VSTTempSmoothInstance>();
    instance->node = node;
    instance->pf_node = pf_node;
    instance->vi = *vi;
    instance->maxr = maxr;
    instance->threshold = threshold;
    instance->mdiff = mdiff;
    instance->strength = strength;
    instance->scthresh = scthresh;
    instance->scenechange = (scthresh != 0.0f);
    instance->fp = fp;
    instance->process = process;

    auto calc_temp_weights = [](int r, int str, std::vector<float>& weights, float& cw) {
      weights.resize(static_cast<std::size_t>(r + 1));
      for (int i = 0; i <= r; ++i) {
        weights[static_cast<std::size_t>(i)] = (i < str) ? 1.0f : 1.0f / static_cast<float>(i - str + 2);
      }
      float sum = weights[0];
      for (int i = 1; i <= r; ++i) {
        sum += weights[static_cast<std::size_t>(i)] * 2.0f;
      }
      for (int i = 0; i <= r; ++i) {
        weights[static_cast<std::size_t>(i)] /= sum;
      }
      cw = weights[0];
    };

    auto calc_diff_weights = [](int thresh, int md, int r, int str, std::vector<float>& diff_weights, float& cw) {
      diff_weights.assign(static_cast<std::size_t>(r) * 256, 0.0f);
      std::vector<float> tw(static_cast<std::size_t>(r + 1), 0.0f);
      for (int i = 0; i <= r; ++i) {
        tw[static_cast<std::size_t>(i)] = (i < str) ? 1.0f : 1.0f / static_cast<float>(i - str + 2);
      }
      std::vector<float> dw(256, 0.0f);
      const float step = 256.0f / static_cast<float>(thresh - std::min(md, thresh - 1));
      float base = 256.0f;
      for (int diff = 0; diff < thresh; ++diff) {
        if (diff < md) {
          dw[static_cast<std::size_t>(diff)] = 256.0f;
        } else {
          if (base > 0.0f) dw[static_cast<std::size_t>(diff)] = base;
          else break;
          base -= step;
        }
      }
      float t_sum = tw[0];
      for (int rad = 1; rad <= r; ++rad) {
        t_sum += tw[static_cast<std::size_t>(rad)] * 2.0f;
        for (int diff = 0; diff < 256; ++diff) {
          diff_weights[static_cast<std::size_t>(rad - 1) * 256 + static_cast<std::size_t>(diff)] =
              tw[static_cast<std::size_t>(rad)] * dw[static_cast<std::size_t>(diff)] / 256.0f;
        }
      }
      for (int rad = 0; rad < r; ++rad) {
        for (int diff = 0; diff < 256; ++diff) {
          diff_weights[static_cast<std::size_t>(rad) * 256 + static_cast<std::size_t>(diff)] /= t_sum;
        }
      }
      cw = tw[0] / t_sum;
    };

    for (int plane = 0; plane < 3; ++plane) {
      if (!process[static_cast<std::size_t>(plane)]) continue;
      if (threshold[static_cast<std::size_t>(plane)] > mdiff[static_cast<std::size_t>(plane)] + 1) {
        instance->weight_mode[static_cast<std::size_t>(plane)] = 0; // inv_diff
        calc_diff_weights(threshold[static_cast<std::size_t>(plane)], mdiff[static_cast<std::size_t>(plane)],
                          maxr, strength, instance->temporal_difference_weights[static_cast<std::size_t>(plane)],
                          instance->center_weights[static_cast<std::size_t>(plane)]);
      } else {
        instance->weight_mode[static_cast<std::size_t>(plane)] = 1; // temporal
        calc_temp_weights(maxr, strength, instance->temporal_weights[static_cast<std::size_t>(plane)],
                          instance->center_weights[static_cast<std::size_t>(plane)]);
      }
    }

    std::vector<VSFilterDependency> deps;
    deps.push_back({instance->node, rpGeneral});
    if (instance->pf_node && instance->pf_node != instance->node) {
      deps.push_back({instance->pf_node, rpGeneral});
    }

    vsapi->createVideoFilter(out, "TTempSmooth", vi, ttempsmooth_get_frame, ttempsmooth_free, fmParallel,
                             deps.data(), static_cast<int>(deps.size()), instance.release(), core);
  } catch (const std::exception& e) {
    if (node) vsapi->freeNode(node);
    if (pf_node) vsapi->freeNode(pf_node);
    vsapi->mapSetError(out, e.what());
  }
}

// ---------------------------------------------------------------------------
// CCD
// ---------------------------------------------------------------------------
struct VSCCDInstance {
  VSNode* node = nullptr;
  VSNode* ref_node = nullptr;
  VSVideoInfo vi{};
  float threshold = 4.0f;
  int temporal_radius = 0;
  std::vector<float> weights;
  std::vector<Point> points;
  int diameter = 0;
  float scale = 1.0f;
};

void VS_CC ccd_free(void* instance_data, VSCore*, const VSAPI* vsapi) {
  auto* d = static_cast<VSCCDInstance*>(instance_data);
  if (d) {
    if (d->node) vsapi->freeNode(d->node);
    if (d->ref_node) vsapi->freeNode(d->ref_node);
    delete d;
  }
}

const VSFrame* VS_CC ccd_get_frame(int n, int activation_reason, void* instance_data, void**,
                                   VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSCCDInstance*>(instance_data);
  try {
    const int first = std::max(0, n - d->temporal_radius);
    const int last = std::min(d->vi.numFrames - 1, n + d->temporal_radius);

    if (activation_reason == arInitial) {
      for (int i = first; i <= last; ++i) {
        vsapi->requestFrameFilter(i, d->node, frame_ctx);
        if (d->ref_node) vsapi->requestFrameFilter(i, d->ref_node, frame_ctx);
      }
    } else if (activation_reason == arAllFramesReady) {
      if (n < d->temporal_radius || n > d->vi.numFrames - 1 - d->temporal_radius) {
        return vsapi->getFrameFilter(n, d->node, frame_ctx);
      }

      const auto free_frame = [vsapi](const VSFrame* f) { vsapi->freeFrame(f); };
      const int temporal_diameter = d->temporal_radius * 2 + 1;

      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> src_frames;
      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> ref_frames;
      src_frames.reserve(static_cast<std::size_t>(temporal_diameter));
      ref_frames.reserve(static_cast<std::size_t>(temporal_diameter));

      for (int i = 0; i < temporal_diameter; ++i) {
        const int frame_num = n - d->temporal_radius + i;
        src_frames.emplace_back(vsapi->getFrameFilter(frame_num, d->node, frame_ctx), free_frame);
        if (!src_frames.back()) return nullptr;
        ref_frames.emplace_back(vsapi->getFrameFilter(frame_num, d->ref_node, frame_ctx), free_frame);
        if (!ref_frames.back()) return nullptr;
      }

      const bool is_rgb = (d->vi.format.colorFamily == cfRGB);
      const VSFrame* center_src = src_frames[static_cast<std::size_t>(d->temporal_radius)].get();

      const VSFrame* plane_src[3] = {
          is_rgb ? nullptr : center_src,
          nullptr,
          nullptr,
      };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, center_src, core);
      if (!dst) throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      std::vector<const std::uint8_t*> src_slices;
      std::vector<const std::uint8_t*> ref_slices;
      src_slices.reserve(static_cast<std::size_t>(temporal_diameter * 3));
      ref_slices.reserve(static_cast<std::size_t>(temporal_diameter * 3));

      for (int i = 0; i < temporal_diameter; ++i) {
        src_slices.push_back(vsapi->getReadPtr(src_frames[static_cast<std::size_t>(i)].get(), 0));
        src_slices.push_back(vsapi->getReadPtr(src_frames[static_cast<std::size_t>(i)].get(), 1));
        src_slices.push_back(vsapi->getReadPtr(src_frames[static_cast<std::size_t>(i)].get(), 2));

        ref_slices.push_back(vsapi->getReadPtr(ref_frames[static_cast<std::size_t>(i)].get(), 0));
        ref_slices.push_back(vsapi->getReadPtr(ref_frames[static_cast<std::size_t>(i)].get(), 1));
        ref_slices.push_back(vsapi->getReadPtr(ref_frames[static_cast<std::size_t>(i)].get(), 2));
      }

      const DataType dtype = get_data_type(d->vi.format.bytesPerSample, d->vi.format.sampleType == stFloat);
      const int target_plane = is_rgb ? 0 : 1;
      const int width = vsapi->getFrameWidth(dst, target_plane);
      const int height = vsapi->getFrameHeight(dst, target_plane);
      const std::size_t stride = static_cast<std::size_t>(vsapi->getStride(dst, target_plane));

      std::uint8_t* dst_r = is_rgb ? vsapi->getWritePtr(dst, 0) : nullptr;
      std::uint8_t* dst_g = vsapi->getWritePtr(dst, 1);
      std::uint8_t* dst_b = vsapi->getWritePtr(dst, 2);

      process_ccd_planes(
          dtype, is_rgb, width, height, stride,
          src_slices.data(), ref_slices.data(),
          dst_r, dst_g, dst_b,
          d->threshold, d->temporal_radius, d->weights.data(),
          d->points.data(), static_cast<int>(d->points.size()),
          d->diameter, d->scale, d->vi.format.bitsPerSample);

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("CCD: frame processing failed", frame_ctx);
  }
  return nullptr;
}

void VS_CC ccd_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* node = nullptr;
  VSNode* ref_node = nullptr;
  try {
    int err = 0;
    node = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && node != nullptr, "CCD: clip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(node);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            "CCD: only constant format input supported");
    require(vi->format.colorFamily == cfRGB || vi->format.colorFamily == cfYUV,
            "CCD: only RGB and YUV color formats are supported");

    const float format_max = (vi->format.sampleType == stFloat) ? 1.0f : static_cast<float>((1u << vi->format.bitsPerSample) - 1u);
    double threshold_user = 4.0;
    const double th_val = vsapi->mapGetFloat(in, "threshold", 0, &err);
    if (!err) threshold_user = th_val;

    require(std::isfinite(threshold_user) && std::isfinite(static_cast<float>(threshold_user)),
            "CCD: threshold must be finite.");
    float threshold = static_cast<float>(threshold_user) * format_max;
    threshold = (threshold * threshold) / (255.0f * 255.0f * 3.0f);
    require(std::isfinite(threshold), "CCD: threshold is too large.");

    int temporal_radius = 0;
    const std::int64_t tr_val = vsapi->mapGetInt(in, "temporal_radius", 0, &err);
    if (!err) {
      require(tr_val >= 0 && tr_val <= 10, "CCD: temporal radius must be <= 10");
      temporal_radius = static_cast<int>(tr_val);
    }

    const int temporal_diameter = temporal_radius * 2 + 1;
    std::vector<float> weights(static_cast<std::size_t>(temporal_diameter), 1.0f);
    for (int r = 0; r < temporal_radius; ++r) {
      const float tr = static_cast<float>(temporal_radius);
      const float fr = static_cast<float>(r);
      weights[static_cast<std::size_t>(temporal_radius - 1 - r)] = std::sqrt((tr + 1.0f - fr) / ((tr + 1.0f) * 2.0f));
      weights[static_cast<std::size_t>(temporal_radius + 1 + r)] = std::sin((tr + 2.0f - fr) / ((tr + 1.0f) * 2.0f));
    }
    weights[static_cast<std::size_t>(temporal_radius)] = 1.0f;

    float scale = static_cast<float>(vi->height >> vi->format.subSamplingH) / 240.0f;
    const double sc_val = vsapi->mapGetFloat(in, "scale", 0, &err);
    if (!err) {
      require(std::isfinite(sc_val) && sc_val >= 1.0 && sc_val <= std::numeric_limits<float>::max(),
              "CCD: scale must be finite and greater than or equal to 1.0");
      scale = static_cast<float>(sc_val);
    }

    require(scale >= 1.0f, "CCD: scale must be greater than or equal to 1.0");
    // Bound coordinates before rounding/narrowing them. The diameter check below
    // applies the reference's stricter geometry rule after points are scaled.
    require(static_cast<double>(scale) * 12.0 <= std::numeric_limits<int>::max() / 4,
            "CCD: scale is too large.");
    const int num_points_arr = vsapi->mapNumElements(in, "points");
    require(num_points_arr < 0 || num_points_arr == 3, "CCD: The points array must have 3 boolean elements.");
    bool low = true, medium = true, high = false;
    if (num_points_arr == 3) {
      low = (vsapi->mapGetInt(in, "points", 0, &err) != 0);
      medium = (vsapi->mapGetInt(in, "points", 1, &err) != 0);
      high = (vsapi->mapGetInt(in, "points", 2, &err) != 0);
    }
    require(low || medium || high, "CCD: A minimum of one set of points must be used.");

    std::vector<Point> points;
    const Point low_pts[] = {{-4, -4}, {4, -4}, {-4, 4}, {4, 4}};
    const Point med_pts[] = {{-8, -8}, {0, -8}, {8, -8}, {-8, 0}, {8, 0}, {-8, 8}, {0, 8}, {8, 8}};
    const Point hi_pts[] = {{-12, -12}, {-4, -12}, {4, -12}, {12, -12},
                            {-12, -4},                         {12, -4},
                            {-12, 4},                          {12, 4},
                            {-12, 12},  {-4, 12},   {4, 12},   {12, 12}};

    if (low) for (const auto& p : low_pts) points.push_back(p);
    if (medium) for (const auto& p : med_pts) points.push_back(p);
    if (high) for (const auto& p : hi_pts) points.push_back(p);

    for (auto& pt : points) {
      pt.x = static_cast<int>(std::round(static_cast<float>(pt.x >> vi->format.subSamplingW) * scale));
      pt.y = static_cast<int>(std::round(static_cast<float>(pt.y >> vi->format.subSamplingH) * scale));
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
      if (a.y < b.y) return true;
      if (a.y == b.y) return a.x < b.x;
      return false;
    });

    int diameter = 0;
    for (const auto& pt : points) {
      diameter = std::max(diameter, std::abs(pt.x) * 2 + 1);
      diameter = std::max(diameter, std::abs(pt.y) * 2 + 1);
    }

    const double scaled_size = std::round(static_cast<double>(diameter) * scale);
    require(scaled_size <= std::numeric_limits<int>::max(), "CCD: scale is too large.");
    const int scaled_diameter = static_cast<int>(scaled_size);
    if (scaled_diameter > vi->width || scaled_diameter > vi->height) {
      char msg[256];
      std::snprintf(msg, sizeof(msg),
                    "CCD: Scale %.1f produces a scaled filter diameter of %d, which is beyond the width %d or height %d. Reduce the scale amount.",
                    scale, scaled_diameter, vi->width, vi->height);
      throw std::runtime_error(msg);
    }

    ref_node = vsapi->mapGetNode(in, "ref", 0, &err);
    if (!err && ref_node) {
      const VSVideoInfo* refvi = vsapi->getVideoInfo(ref_node);
      require(is_same_video_info(vi, refvi), "CCD: ref and source clip format, width, and height must match.");
    } else {
      ref_node = vsapi->addNodeRef(node);
    }

    const auto free_node = [vsapi](VSNode* n) { vsapi->freeNode(n); };
    // YUV luma downscaling
    const bool has_luma = (vi->format.colorFamily == cfYUV);
    const bool should_resize = has_luma && (vi->format.subSamplingW > 0 || vi->format.subSamplingH > 0);
    if (should_resize) {
      VSPlugin* std_plugin = vsapi->getPluginByID("com.vapoursynth.std", core);
      VSPlugin* resize_plugin = vsapi->getPluginByID("com.vapoursynth.resize", core);
      require(std_plugin && resize_plugin, "CCD: std and resize plugins required for YUV subsampled processing.");

      const auto free_map = [vsapi](VSMap* m) { vsapi->freeMap(m); };

      // 1. Extract luma
      std::unique_ptr<VSMap, decltype(free_map)> a1(vsapi->createMap(), free_map);
      if (!a1) throw std::bad_alloc();
      vsapi->mapSetNode(a1.get(), "clips", ref_node, maAppend);
      vsapi->mapSetInt(a1.get(), "planes", 0, maAppend);
      vsapi->mapSetInt(a1.get(), "colorfamily", static_cast<int>(cfGray), maAppend);
      std::unique_ptr<VSMap, decltype(free_map)> r1(vsapi->invoke(std_plugin, "ShufflePlanes", a1.get()), free_map);
      if (!r1) throw std::bad_alloc();
      const char* err1 = vsapi->mapGetError(r1.get());
      if (err1) throw std::runtime_error(err1);
      int rerr1 = 0;
      std::unique_ptr<VSNode, decltype(free_node)> luma_owner(vsapi->mapGetNode(r1.get(), "clip", 0, &rerr1), free_node);
      VSNode* luma = luma_owner.get();
      require(!rerr1 && luma != nullptr, "CCD: failed to extract luma.");

      // 2. Resize luma
      std::unique_ptr<VSMap, decltype(free_map)> a2(vsapi->createMap(), free_map);
      if (!a2) throw std::bad_alloc();
      vsapi->mapSetNode(a2.get(), "clip", luma, maAppend);
      luma_owner.reset();
      vsapi->mapSetInt(a2.get(), "width", vi->width >> vi->format.subSamplingW, maAppend);
      vsapi->mapSetInt(a2.get(), "height", vi->height >> vi->format.subSamplingH, maAppend);
      std::unique_ptr<VSMap, decltype(free_map)> r2(vsapi->invoke(resize_plugin, "Bilinear", a2.get()), free_map);
      if (!r2) throw std::bad_alloc();
      const char* err2 = vsapi->mapGetError(r2.get());
      if (err2) throw std::runtime_error(err2);
      int rerr2 = 0;
      std::unique_ptr<VSNode, decltype(free_node)> resized_owner(vsapi->mapGetNode(r2.get(), "clip", 0, &rerr2), free_node);
      VSNode* resized_luma = resized_owner.get();
      require(!rerr2 && resized_luma != nullptr, "CCD: failed to resize luma.");

      // 3. Merge resized luma with chroma
      std::unique_ptr<VSMap, decltype(free_map)> a3(vsapi->createMap(), free_map);
      if (!a3) throw std::bad_alloc();
      vsapi->mapSetNode(a3.get(), "clips", resized_luma, maAppend);
      resized_owner.reset();
      vsapi->mapSetNode(a3.get(), "clips", ref_node, maAppend);
      vsapi->freeNode(ref_node);
      ref_node = nullptr;
      vsapi->mapSetInt(a3.get(), "planes", 0, maAppend);
      vsapi->mapSetInt(a3.get(), "planes", 1, maAppend);
      vsapi->mapSetInt(a3.get(), "planes", 2, maAppend);
      vsapi->mapSetInt(a3.get(), "colorfamily", static_cast<int>(cfYUV), maAppend);
      std::unique_ptr<VSMap, decltype(free_map)> r3(vsapi->invoke(std_plugin, "ShufflePlanes", a3.get()), free_map);
      if (!r3) throw std::bad_alloc();
      const char* err3 = vsapi->mapGetError(r3.get());
      if (err3) throw std::runtime_error(err3);
      int rerr3 = 0;
      ref_node = vsapi->mapGetNode(r3.get(), "clip", 0, &rerr3);
      require(!rerr3 && ref_node != nullptr, "CCD: failed to merge resized luma.");
    }

    auto instance = std::make_unique<VSCCDInstance>();
    instance->node = node;
    instance->ref_node = ref_node;
    instance->vi = *vi;
    instance->threshold = threshold;
    instance->temporal_radius = temporal_radius;
    instance->weights = std::move(weights);
    instance->points = std::move(points);
    instance->diameter = diameter;
    instance->scale = scale;

    std::vector<VSFilterDependency> deps;
    const int source_pattern = temporal_radius == 0 ? rpStrictSpatial : rpGeneral;
    deps.push_back({instance->node, source_pattern});
    if (instance->ref_node && instance->ref_node != instance->node) {
      // A shorter reference is clamped by the host, so it is not strictly
      // spatial even when the source requests only the current frame.
      const auto* refvi = vsapi->getVideoInfo(instance->ref_node);
      const int reference_pattern = temporal_radius == 0 && refvi->numFrames >= vi->numFrames
          ? rpStrictSpatial : rpGeneral;
      deps.push_back({instance->ref_node, reference_pattern});
    }

    vsapi->createVideoFilter(out, "CCD", vi, ccd_get_frame, ccd_free, fmParallel,
                             deps.data(), static_cast<int>(deps.size()), instance.release(), core);
  } catch (const std::exception& e) {
    if (node) vsapi->freeNode(node);
    if (ref_node) vsapi->freeNode(ref_node);
    vsapi->mapSetError(out, e.what());
  }
}

// ---------------------------------------------------------------------------
// Cnr4
// ---------------------------------------------------------------------------
struct VSCnr4Instance {
  VSNode* node = nullptr;
  VSNode* ref_node = nullptr;
  VSNode* luma_node = nullptr;
  VSNode* ref_luma_node = nullptr;
  VSVideoInfo vi{};
  int radius = 2;
  int tmode = 0;
  int wmode = 0;
  bool scenechange = true;
  std::array<std::uint8_t, 256> table_y{};
  std::array<std::uint8_t, 256> table_u{};
  std::array<std::uint8_t, 256> table_v{};
};

void VS_CC cnr4_free(void* instance_data, VSCore*, const VSAPI* vsapi) {
  auto* d = static_cast<VSCnr4Instance*>(instance_data);
  if (d) {
    if (d->node) vsapi->freeNode(d->node);
    if (d->ref_node) vsapi->freeNode(d->ref_node);
    if (d->luma_node) vsapi->freeNode(d->luma_node);
    if (d->ref_luma_node) vsapi->freeNode(d->ref_luma_node);
    delete d;
  }
}

const VSFrame* VS_CC cnr4_get_frame(int n, int activation_reason, void* instance_data, void**,
                                    VSFrameContext* frame_ctx, VSCore* core, const VSAPI* vsapi) {
  auto* d = static_cast<VSCnr4Instance*>(instance_data);
  try {
    const bool src_includes_center = (d->tmode != 0); // inv_diff does not include center in surrounding array

    if (activation_reason == arInitial) {
      for (int i = -d->radius; i <= d->radius; ++i) {
        const int fn = std::clamp(n + i, 0, d->vi.numFrames - 1);
        vsapi->requestFrameFilter(fn, d->node, frame_ctx);
        vsapi->requestFrameFilter(fn, d->luma_node, frame_ctx);
        if (d->ref_node) {
          vsapi->requestFrameFilter(fn, d->ref_node, frame_ctx);
          vsapi->requestFrameFilter(fn, d->ref_luma_node, frame_ctx);
        }
      }
    } else if (activation_reason == arAllFramesReady) {
      const auto free_frame = [vsapi](const VSFrame* f) { vsapi->freeFrame(f); };

      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> src_frames;
      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> luma_frames;
      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> ref_frames;
      std::vector<std::unique_ptr<const VSFrame, decltype(free_frame)>> ref_luma_frames;

      const std::size_t capacity = static_cast<std::size_t>(d->radius * 2 + 1);
      src_frames.reserve(capacity);
      luma_frames.reserve(capacity);
      ref_frames.reserve(capacity);
      ref_luma_frames.reserve(capacity);

      for (int i = -d->radius; i <= d->radius; ++i) {
        if (!src_includes_center && i == 0) continue;
        const int fn = std::clamp(n + i, 0, d->vi.numFrames - 1);
        src_frames.emplace_back(vsapi->getFrameFilter(fn, d->node, frame_ctx), free_frame);
        if (!src_frames.back()) return nullptr;
        luma_frames.emplace_back(vsapi->getFrameFilter(fn, d->luma_node, frame_ctx), free_frame);
        if (!luma_frames.back()) return nullptr;
        if (d->ref_node) {
          ref_frames.emplace_back(vsapi->getFrameFilter(fn, d->ref_node, frame_ctx), free_frame);
        if (!ref_frames.back()) return nullptr;
          ref_luma_frames.emplace_back(vsapi->getFrameFilter(fn, d->ref_luma_node, frame_ctx), free_frame);
        if (!ref_luma_frames.back()) return nullptr;
        if (!luma_frames.back()) return nullptr;
        }
      }

      const int frame_count = static_cast<int>(src_frames.size());

      std::unique_ptr<const VSFrame, decltype(free_frame)> curr(vsapi->getFrameFilter(n, d->node, frame_ctx), free_frame);
      if (!curr) return nullptr;
      std::unique_ptr<const VSFrame, decltype(free_frame)> curr_luma(vsapi->getFrameFilter(n, d->luma_node, frame_ctx), free_frame);
      if (!curr_luma) return nullptr;
      std::unique_ptr<const VSFrame, decltype(free_frame)> curr_ref(
          d->ref_node ? vsapi->getFrameFilter(n, d->ref_node, frame_ctx) : vsapi->addFrameRef(curr.get()), free_frame);
      if (!curr_ref) return nullptr;
      std::unique_ptr<const VSFrame, decltype(free_frame)> curr_ref_luma(
          d->ref_node ? vsapi->getFrameFilter(n, d->ref_luma_node, frame_ctx) : vsapi->addFrameRef(curr_luma.get()), free_frame);

      if (!curr_ref_luma) return nullptr;
      int start_idx = 0;
      int end_idx = frame_count - 1;

      if (d->scenechange) {
        const auto* props0 = vsapi->getFramePropertiesRO(src_frames[0].get());
        int sc_err1 = 0, sc_err2 = 0;
        vsapi->mapGetInt(props0, "_SceneChangePrev", 0, &sc_err1);
        vsapi->mapGetInt(props0, "_SceneChangeNext", 0, &sc_err2);
        if (sc_err1 || sc_err2) {
          vsapi->setFilterError("Cnr4: Scene change handling enabled, but input frame is missing scene change properties. "
                                "Either set scenechange=False or run scene change detection on your input clip.", frame_ctx);
          return nullptr;
        }

        // Walk backwards on left
        int mid_left = src_includes_center ? (frame_count / 2) : (frame_count / 2 - 1);
        for (int i = mid_left; i > 0; --i) {
          int sc_err = 0;
          const auto* p = vsapi->getFramePropertiesRO(src_frames[static_cast<std::size_t>(i)].get());
          if (vsapi->mapGetInt(p, "_SceneChangePrev", 0, &sc_err) == 1 && !sc_err) {
            start_idx = i;
            break;
          }
        }

        // Walk forwards on right
        for (int i = frame_count / 2; i < frame_count - 1; ++i) {
          int sc_err = 0;
          const auto* p = vsapi->getFramePropertiesRO(src_frames[static_cast<std::size_t>(i)].get());
          if (vsapi->mapGetInt(p, "_SceneChangeNext", 0, &sc_err) == 1 && !sc_err) {
            end_idx = i;
            break;
          }
        }

        if (!src_includes_center) {
          const auto* curr_props = vsapi->getFramePropertiesRO(curr.get());
          int sc_err = 0;
          if (vsapi->mapGetInt(curr_props, "_SceneChangePrev", 0, &sc_err) == 1 && !sc_err) {
            start_idx = frame_count / 2;
          }
          if (vsapi->mapGetInt(curr_props, "_SceneChangeNext", 0, &sc_err) == 1 && !sc_err) {
            end_idx = frame_count / 2 - 1;
          }
        }
      }

      // Replace unusable frames with center frame
      for (int i = 0; i < start_idx; ++i) {
        src_frames[static_cast<std::size_t>(i)].reset(vsapi->addFrameRef(curr.get()));
        luma_frames[static_cast<std::size_t>(i)].reset(vsapi->addFrameRef(curr_luma.get()));
        if (d->ref_node) {
          ref_frames[static_cast<std::size_t>(i)].reset(vsapi->addFrameRef(curr_ref.get()));
          ref_luma_frames[static_cast<std::size_t>(i)].reset(vsapi->addFrameRef(curr_ref_luma.get()));
        }
      }
      for (int i = end_idx + 1; i < frame_count; ++i) {
        src_frames[static_cast<std::size_t>(i)].reset(vsapi->addFrameRef(curr.get()));
        luma_frames[static_cast<std::size_t>(i)].reset(vsapi->addFrameRef(curr_luma.get()));
        if (d->ref_node) {
          ref_frames[static_cast<std::size_t>(i)].reset(vsapi->addFrameRef(curr_ref.get()));
          ref_luma_frames[static_cast<std::size_t>(i)].reset(vsapi->addFrameRef(curr_ref_luma.get()));
        }
      }

      // Build src/ref arrays: plane 0 = luma, plane 1 = U, plane 2 = V
      std::vector<std::array<const std::uint8_t*, 3>> src_planes(static_cast<std::size_t>(frame_count));
      std::vector<std::array<const std::uint8_t*, 3>> ref_planes(static_cast<std::size_t>(frame_count));

      for (int i = 0; i < frame_count; ++i) {
        src_planes[static_cast<std::size_t>(i)][0] = vsapi->getReadPtr(luma_frames[static_cast<std::size_t>(i)].get(), 0);
        src_planes[static_cast<std::size_t>(i)][1] = vsapi->getReadPtr(src_frames[static_cast<std::size_t>(i)].get(), 1);
        src_planes[static_cast<std::size_t>(i)][2] = vsapi->getReadPtr(src_frames[static_cast<std::size_t>(i)].get(), 2);

        ref_planes[static_cast<std::size_t>(i)][0] = d->ref_node ? vsapi->getReadPtr(ref_luma_frames[static_cast<std::size_t>(i)].get(), 0) : src_planes[static_cast<std::size_t>(i)][0];
        ref_planes[static_cast<std::size_t>(i)][1] = d->ref_node ? vsapi->getReadPtr(ref_frames[static_cast<std::size_t>(i)].get(), 1) : src_planes[static_cast<std::size_t>(i)][1];
        ref_planes[static_cast<std::size_t>(i)][2] = d->ref_node ? vsapi->getReadPtr(ref_frames[static_cast<std::size_t>(i)].get(), 2) : src_planes[static_cast<std::size_t>(i)][2];
      }

      const VSFrame* plane_src[3] = { curr.get(), nullptr, nullptr };
      int planes[3] = {0, 1, 2};
      VSFrame* dst = vsapi->newVideoFrame2(&d->vi.format, d->vi.width, d->vi.height, plane_src, planes, curr.get(), core);
      if (!dst) throw std::bad_alloc();
      std::unique_ptr<VSFrame, decltype(free_frame)> dst_owner(dst, free_frame);

      const DataType dtype = get_data_type(d->vi.format.bytesPerSample, false);
      const int width_uv = vsapi->getFrameWidth(dst, 1);
      const int height_uv = vsapi->getFrameHeight(dst, 1);
      const std::size_t stride_bytes = static_cast<std::size_t>(vsapi->getStride(dst, 1));

      const std::uint8_t* const curr_ptrs[3] = {
          vsapi->getReadPtr(curr_luma.get(), 0),
          vsapi->getReadPtr(curr.get(), 1),
          vsapi->getReadPtr(curr.get(), 2),
      };
      const std::uint8_t* const curr_ref_ptrs[3] = {
          d->ref_node ? vsapi->getReadPtr(curr_ref_luma.get(), 0) : curr_ptrs[0],
          d->ref_node ? vsapi->getReadPtr(curr_ref.get(), 1) : curr_ptrs[1],
          d->ref_node ? vsapi->getReadPtr(curr_ref.get(), 2) : curr_ptrs[2],
      };

      std::uint8_t* dst_u = vsapi->getWritePtr(dst, 1);
      std::uint8_t* dst_v = vsapi->getWritePtr(dst, 2);

      process_cnr4_frame(
          dtype, width_uv, height_uv, stride_bytes,
          d->vi.format.bitsPerSample, d->radius, d->tmode, d->wmode,
          curr_ptrs, curr_ref_ptrs,
          src_planes.data(),
          ref_planes.data(),
          static_cast<std::size_t>(frame_count),
          dst_u, dst_v,
          d->table_y.data(), d->table_u.data(), d->table_v.data());

      return dst_owner.release();
    }
  } catch (const std::exception& e) {
    vsapi->setFilterError(e.what(), frame_ctx);
  } catch (...) {
    vsapi->setFilterError("Cnr4: frame processing failed", frame_ctx);
  }
  return nullptr;
}

void VS_CC cnr4_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi) {
  VSNode* node = nullptr;
  VSNode* ref_node = nullptr;
  VSNode* luma_node = nullptr;
  VSNode* ref_luma_node = nullptr;
  try {
    int err = 0;
    node = vsapi->mapGetNode(in, "clip", 0, &err);
    require(!err && node != nullptr, "Cnr4: clip is required.");

    const VSVideoInfo* vi = vsapi->getVideoInfo(node);
    require(vi->format.colorFamily != cfUndefined && vi->width > 0 && vi->height > 0,
            "Cnr4: only constant format input supported");
    require(vi->format.colorFamily == cfYUV, "Cnr4: only YUV color formats are supported");
    require(vi->format.sampleType == stInteger && vi->format.bitsPerSample >= 8 && vi->format.bitsPerSample <= 16,
            "Cnr4: only 8-16 bit integer formats are supported");

    std::string mode = "oxx";
    const char* mode_str = vsapi->mapGetData(in, "mode", 0, &err);
    if (!err && mode_str) {
      const int mode_len = vsapi->mapGetDataSize(in, "mode", 0, &err);
      require(mode_len == 3, "Cnr4: mode must have 3 characters");
      mode = std::string(mode_str, 3);
    }
    for (char c : mode) {
      require(c == 'o' || c == 'x', "Cnr4: Only 'o' and 'x' are recognized characters in mode");
    }

    int radius = 2;
    const std::int64_t rad_val = vsapi->mapGetInt(in, "radius", 0, &err);
    if (!err) {
      require(rad_val >= 1 && rad_val <= 10, "Cnr4: radius must be between 1 and 10");
      radius = static_cast<int>(rad_val);
    }

    int tmode = 0;
    const std::int64_t tm_val = vsapi->mapGetInt(in, "tmode", 0, &err);
    if (!err) {
      require(tm_val >= 0 && tm_val <= 4, "Cnr4: tmode can only be between 0 and 4");
      tmode = static_cast<int>(tm_val);
    }

    int wmode = 0;
    const std::int64_t wm_val = vsapi->mapGetInt(in, "wmode", 0, &err);
    if (!err) {
      require(wm_val >= 0 && wm_val <= 3, "Cnr4: wmode can only be between 0 and 3");
      wmode = static_cast<int>(wm_val);
    }

    bool scenechange = true;
    const std::int64_t sc_val = vsapi->mapGetInt(in, "scenechange", 0, &err);
    if (!err) {
      scenechange = (sc_val != 0);
    }

    // Sense
    std::array<int, 3> sense{35, 47, 47};
    const int num_sense = vsapi->mapNumElements(in, "sense");
    if (num_sense >= 0) {
      require(num_sense == 3, "Cnr4: sense must have 3 elements");
      for (int i = 0; i < 3; ++i) {
        const std::int64_t s_val = vsapi->mapGetInt(in, "sense", i, &err);
        require(!err && s_val >= -1 && s_val <= 255, "Cnr4: sense must be between -1 and 255");
        if (s_val != -1) sense[static_cast<std::size_t>(i)] = static_cast<int>(s_val);
      }
    }

    // Str
    std::array<int, 3> str{192, 255, 255};
    const int num_str = vsapi->mapNumElements(in, "str");
    if (num_str >= 0) {
      require(num_str == 3, "Cnr4: str must have 3 elements");
      const char* const str_names[3] = {"l_str", "u_str", "v_str"};
      for (int i = 0; i < 3; ++i) {
        const std::int64_t st_val = vsapi->mapGetInt(in, "str", i, &err);
        if (err || st_val < -1 || st_val > 255) {
          char msg[128];
          std::snprintf(msg, sizeof(msg), "Cnr4: %s must be between -1 and 255", str_names[i]);
          throw std::runtime_error(msg);
        }
        if (st_val != -1) str[static_cast<std::size_t>(i)] = static_cast<int>(st_val);
      }
    }

    // Pow
    std::array<float, 3> pow_vals{1.0f, 1.0f, 1.0f};
    const int num_pow = vsapi->mapNumElements(in, "pow");
    if (num_pow >= 0) {
      require(num_pow == 3, "Cnr4: pow must contain 3 elements");
      for (int i = 0; i < 3; ++i) {
        const double p_val = vsapi->mapGetFloat(in, "pow", i, &err);
        require(!err && std::isfinite(p_val) && p_val >= 0.0 && p_val <= std::numeric_limits<float>::max(),
                "Cnr4: pow must be finite and >= 0");
        pow_vals[static_cast<std::size_t>(i)] = static_cast<float>(p_val);
      }
    }

    ref_node = vsapi->mapGetNode(in, "ref", 0, &err);
    if (!err && ref_node) {
      const VSVideoInfo* refvi = vsapi->getVideoInfo(ref_node);
      require(is_same_video_info(vi, refvi), "Cnr4: ref must be the same video format, width, and height as the source clip");
    }

    // Extract luma plane
    VSPlugin* std_plugin = vsapi->getPluginByID("com.vapoursynth.std", core);
    VSPlugin* resize_plugin = vsapi->getPluginByID("com.vapoursynth.resize", core);
    require(std_plugin && resize_plugin, "Cnr4: std and resize plugins required.");

    const auto free_map = [vsapi](VSMap* m) { vsapi->freeMap(m); };

    const auto free_node = [vsapi](VSNode* n) { vsapi->freeNode(n); };
    auto extract_and_resize_luma = [&](VSNode* src_clip, VSNode*& out_luma) {
      std::unique_ptr<VSMap, decltype(free_map)> a1(vsapi->createMap(), free_map);
      if (!a1) throw std::bad_alloc();
      vsapi->mapSetNode(a1.get(), "clips", src_clip, maAppend);
      vsapi->mapSetInt(a1.get(), "planes", 0, maAppend);
      vsapi->mapSetInt(a1.get(), "colorfamily", static_cast<int>(cfGray), maAppend);
      std::unique_ptr<VSMap, decltype(free_map)> r1(vsapi->invoke(std_plugin, "ShufflePlanes", a1.get()), free_map);
      if (!r1) throw std::bad_alloc();
      const char* err1 = vsapi->mapGetError(r1.get());
      if (err1) throw std::runtime_error(err1);
      int rerr1 = 0;
      std::unique_ptr<VSNode, decltype(free_node)> luma_owner(vsapi->mapGetNode(r1.get(), "clip", 0, &rerr1), free_node);
      VSNode* luma = luma_owner.get();
      require(!rerr1 && luma != nullptr, "Cnr4: failed to extract luma.");

      if (vi->format.subSamplingW > 0 || vi->format.subSamplingH > 0) {
        std::unique_ptr<VSMap, decltype(free_map)> a2(vsapi->createMap(), free_map);
      if (!a2) throw std::bad_alloc();
        vsapi->mapSetNode(a2.get(), "clip", luma, maAppend);
        luma_owner.reset();
        vsapi->mapSetInt(a2.get(), "width", vi->width >> vi->format.subSamplingW, maAppend);
        vsapi->mapSetInt(a2.get(), "height", vi->height >> vi->format.subSamplingH, maAppend);
        std::unique_ptr<VSMap, decltype(free_map)> r2(vsapi->invoke(resize_plugin, "Bilinear", a2.get()), free_map);
      if (!r2) throw std::bad_alloc();
        const char* err2 = vsapi->mapGetError(r2.get());
        if (err2) throw std::runtime_error(err2);
        int rerr2 = 0;
        out_luma = vsapi->mapGetNode(r2.get(), "clip", 0, &rerr2);
        require(!rerr2 && out_luma != nullptr, "Cnr4: failed to resize luma.");
      } else {
        out_luma = luma_owner.release();
      }
    };

    extract_and_resize_luma(node, luma_node);
    if (ref_node) extract_and_resize_luma(ref_node, ref_luma_node);

    auto instance = std::make_unique<VSCnr4Instance>();
    instance->node = node;
    instance->ref_node = ref_node;
    instance->luma_node = luma_node;
    instance->ref_luma_node = ref_luma_node;
    instance->vi = *vi;
    instance->radius = radius;
    instance->tmode = tmode;
    instance->wmode = wmode;
    instance->scenechange = scenechange;

    // Calculate LUTs
    instance->table_y.fill(0);
    instance->table_u.fill(0);
    instance->table_v.fill(0);

    const float pi = 3.14159265358979323846f;
    auto build_table = [&](char m, int s_sense, int s_str, float p_val, std::array<std::uint8_t, 256>& tbl) {
      const float str_f = static_cast<float>(s_str);
      const float sense_f = static_cast<float>(s_sense);
      // Zero sensitivity admits identical samples only. Avoid 0/0 and NaN to
      // integer conversion; zero power is the limiting infinitely steep curve.
      if (s_sense == 0) { tbl[0] = static_cast<std::uint8_t>(s_str); return; }
      const float inv_p = p_val > 0.0f ? (1.0f / p_val) : std::numeric_limits<float>::infinity();
      for (int l = 0; l <= s_str; ++l) {
        const float lf = static_cast<float>(l);
        float base_val = 0.0f;
        if (m == 'o') {
          base_val = (1.0f + std::cos(lf * lf * pi / (sense_f * sense_f))) / 2.0f;
        } else if (m == 'x') {
          base_val = (1.0f + std::cos(lf * pi / sense_f)) / 2.0f;
        }
        const float shaped = std::isinf(inv_p) ? (base_val >= 1.0f ? 1.0f : 0.0f) : std::pow(base_val, inv_p);
        tbl[static_cast<std::size_t>(l)] = static_cast<std::uint8_t>(std::clamp(str_f * shaped, 0.0f, 255.0f));
      }
    };

    build_table(mode[0], sense[0], str[0], pow_vals[0], instance->table_y);
    build_table(mode[1], sense[1], str[1], pow_vals[1], instance->table_u);
    build_table(mode[2], sense[2], str[2], pow_vals[2], instance->table_v);

    std::vector<VSFilterDependency> deps;
    deps.push_back({instance->node, rpGeneral});
    deps.push_back({instance->luma_node, rpGeneral});
    if (instance->ref_node) {
      if (instance->ref_node != instance->node) deps.push_back({instance->ref_node, rpGeneral});
      if (instance->ref_luma_node != instance->luma_node) deps.push_back({instance->ref_luma_node, rpGeneral});
    }

    vsapi->createVideoFilter(out, "Cnr4", vi, cnr4_get_frame, cnr4_free, fmParallel,
                             deps.data(), static_cast<int>(deps.size()), instance.release(), core);
  } catch (const std::exception& e) {
    if (node) vsapi->freeNode(node);
    if (ref_node) vsapi->freeNode(ref_node);
    if (luma_node) vsapi->freeNode(luma_node);
    if (ref_luma_node) vsapi->freeNode(ref_luma_node);
    vsapi->mapSetError(out, e.what());
  }
}

} // namespace
} // namespace neo_smo::plugin

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin* plugin, const VSPLUGINAPI* vspapi) {
  using namespace neo_smo;
  using namespace neo_smo::plugin;

  vspapi->configPlugin("org.neofilters.neo_smo", "neo_smo", "neo-smo " NEO_SMO_VERSION_STRING " smoothing filters",
                       VS_MAKE_VERSION(NEO_SMO_VERSION_MAJOR, NEO_SMO_VERSION_MINOR), VAPOURSYNTH_API_VERSION, 0,
                       plugin);

  vspapi->registerFunction("Median", vs_signature(Algorithm::Median), "clip:vnode;", filter_create<Algorithm::Median>,
                           nullptr, plugin);
  vspapi->registerFunction("VerticalCleaner", vs_signature(Algorithm::VerticalCleaner), "clip:vnode;",
                           filter_create<Algorithm::VerticalCleaner>, nullptr, plugin);
  vspapi->registerFunction("RemoveGrain", vs_signature(Algorithm::RemoveGrain), "clip:vnode;",
                           filter_create<Algorithm::RemoveGrain>, nullptr, plugin);
  vspapi->registerFunction("Repair", vs_signature(Algorithm::Repair), "clip:vnode;",
                           repair_create, nullptr, plugin);
  vspapi->registerFunction("Clense", vs_signature(Algorithm::Clense), "clip:vnode;",
                           clense_create<Algorithm::Clense>, nullptr, plugin);
  vspapi->registerFunction("ForwardClense", vs_signature(Algorithm::ForwardClense), "clip:vnode;",
                           clense_create<Algorithm::ForwardClense>, nullptr, plugin);
  vspapi->registerFunction("BackwardClense", vs_signature(Algorithm::BackwardClense), "clip:vnode;",
                           clense_create<Algorithm::BackwardClense>, nullptr, plugin);
  vspapi->registerFunction("InterQuartileMean", vs_signature(Algorithm::InterQuartileMean), "clip:vnode;",
                           filter_create<Algorithm::InterQuartileMean>, nullptr, plugin);
  vspapi->registerFunction("SmartMedian", vs_signature(Algorithm::SmartMedian), "clip:vnode;",
                           filter_create<Algorithm::SmartMedian>, nullptr, plugin);
  vspapi->registerFunction("TemporalMedian", vs_signature(Algorithm::TemporalMedian), "clip:vnode;",
                           temporal_filter_create<Algorithm::TemporalMedian>, nullptr, plugin);
  vspapi->registerFunction("TemporalSoften", vs_signature(Algorithm::TemporalSoften), "clip:vnode;",
                           temporal_filter_create<Algorithm::TemporalSoften>, nullptr, plugin);
  vspapi->registerFunction("TemporalRepair", vs_signature(Algorithm::TemporalRepair), "clip:vnode;",
                           temporal_repair_create, nullptr, plugin);
  vspapi->registerFunction("DegrainMedian", vs_signature(Algorithm::DegrainMedian), "clip:vnode;",
                           trio_filter_create<Algorithm::DegrainMedian>, nullptr, plugin);
  vspapi->registerFunction("FluxSmoothT", vs_signature(Algorithm::FluxSmoothT), "clip:vnode;",
                           trio_filter_create<Algorithm::FluxSmoothT>, nullptr, plugin);
  vspapi->registerFunction("FluxSmoothST", vs_signature(Algorithm::FluxSmoothST), "clip:vnode;",
                           trio_filter_create<Algorithm::FluxSmoothST>, nullptr, plugin);
  vspapi->registerFunction("TTempSmooth", vs_signature(Algorithm::TTempSmooth), "clip:vnode;",
                           ttempsmooth_create, nullptr, plugin);
  vspapi->registerFunction("CCD", vs_signature(Algorithm::CCD), "clip:vnode;",
                           ccd_create, nullptr, plugin);
  vspapi->registerFunction("Cnr4", vs_signature(Algorithm::Cnr4), "clip:vnode;",
                           cnr4_create, nullptr, plugin);
  vspapi->registerFunction("DCTFilter", "clip:vnode;factors:float[];planes:int[]:opt;", "clip:vnode;",
                           dct_filter_create, nullptr, plugin);
}
