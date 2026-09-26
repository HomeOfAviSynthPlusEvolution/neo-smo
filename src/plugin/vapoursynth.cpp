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

int read_int(const VSMap* in, const char* key, int index, const VSAPI* vsapi) {
  int err = 0;
  const std::int64_t value = vsapi->mapGetInt(in, key, index, &err);
  require(!err, "neo-smo: invalid integer parameter.");
  require(value >= std::numeric_limits<int>::min() && value <= std::numeric_limits<int>::max(),
          "neo-smo: integer parameter out of range.");
  return static_cast<int>(value);
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
}
