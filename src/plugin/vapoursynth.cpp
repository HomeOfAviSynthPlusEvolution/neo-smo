#include "plugin/descriptors.hpp"
#include "neo_smo_version.hpp"
#include <vapoursynth/VapourSynth4.h>
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

    const char* param_key = (Alg == Algorithm::Median) ? "radius" : "mode";
    std::vector<int> param_list;
    const int num_params = vsapi->mapNumElements(in, param_key);
    for (int i = 0; i < num_params; ++i) {
      param_list.push_back(read_int(in, param_key, i, vsapi));
    }

    std::vector<int> planes_list;
    const int num_planes = vsapi->mapNumElements(in, "planes");
    const bool planes_specified = (num_planes >= 0);
    for (int i = 0; i < num_planes; ++i) {
      planes_list.push_back(read_int(in, "planes", i, vsapi));
    }

    FilterPlan plan = build_plan(Alg, fmt, param_list, planes_list, planes_specified);

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
}
