#include "plugin/descriptors.hpp"
#include "common/temporal_window.hpp"
#include "neo_smo_version.hpp"
#include <vapoursynth/VapourSynth4.h>
#include <cmath>
#include <memory>
#include <limits>

namespace neo_smo::plugin {
namespace {

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
}
