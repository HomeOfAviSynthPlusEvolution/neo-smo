#include "plugin/deen.hpp"
#include "algorithms/deen.hpp"
#include "algorithms/scene_detection.hpp"
#include "base/checked.hpp"
#include "plugin/deen_params_vs.hpp"
#include <memory>

namespace neo_smo::plugin {
namespace {
struct Instance {
  explicit Instance(DeenConfig c) : config(std::move(c)) {
    for (int p = 0; p < 3; ++p)
      filters[p] = std::make_unique<Deen>(deen_plane_filter(config, p));
  }
  DeenConfig config;
  bool internal_scene = false;
  std::array<std::unique_ptr<Deen>, 3> filters;
  VSNode* node = nullptr;
  VSVideoInfo vi{};
};
struct FrameDeleter {
  const VSAPI* api;
  void operator()(const VSFrame* p) const {
    if (p)
      api->freeFrame(p);
  }
};
std::vector<DeenPlane> views(const VSFrame* f, const VSVideoInfo& vi, const VSAPI* api) {
  const auto& a = *api->getVideoFrameFormat(f);
  const auto& b = vi.format;
  require(a.colorFamily == b.colorFamily && a.sampleType == b.sampleType && a.bitsPerSample == b.bitsPerSample &&
              a.bytesPerSample == b.bytesPerSample && a.numPlanes == b.numPlanes && a.subSamplingW == b.subSamplingW &&
              a.subSamplingH == b.subSamplingH,
          "Deen: input format changed.");
  std::vector<DeenPlane> out;
  for (int p = 0; p < b.numPlanes; ++p) {
    const int w = api->getFrameWidth(f, p), h = api->getFrameHeight(f, p);
    const bool uv = b.colorFamily == cfYUV && p > 0;
    require(w == (vi.width >> (uv ? b.subSamplingW : 0)) && h == (vi.height >> (uv ? b.subSamplingH : 0)),
            "Deen: input dimensions changed.");
    out.push_back({api->getReadPtr(f, p), api->getStride(f, p), w, h,
                   get_data_type(b.bytesPerSample, b.sampleType == stFloat), b.bitsPerSample});
  }
  return out;
}
const VSFrame* VS_CC get_frame(int n, int activation, void* data, void**, VSFrameContext* ctx, VSCore* core,
                               const VSAPI* api) {
  const auto& d = *static_cast<Instance*>(data);
  const bool temporal = d.config.temporal && n > 0 && n < d.vi.numFrames - 1;
  const int count = temporal || d.internal_scene ? 3 : 1;
  const int indices[3] = {n, std::max(n - 1, 0), std::min(n + 1, d.vi.numFrames - 1)};
  try {
    if (activation == arInitial) {
      for (int i = 0; i < count; ++i)
        api->requestFrameFilter(indices[i], d.node, ctx);
    } else if (activation == arAllFramesReady) {
      using Frame = std::unique_ptr<const VSFrame, FrameDeleter>;
      std::array<Frame, 3> frames{{Frame(nullptr, {api}), Frame(nullptr, {api}), Frame(nullptr, {api})}};
      std::array<std::vector<DeenPlane>, 3> input;
      for (int i = 0; i < count; ++i) {
        frames[i].reset(api->getFrameFilter(indices[i], d.node, ctx));
        if (!frames[i])
          return nullptr;
        input[i] = views(frames[i].get(), d.vi, api);
      }
      bool use_temporal = temporal;
      bool left = false, right = false;
      if (d.internal_scene) {
        const double threshold = d.config.scenechange / 255.0;
        right = n < d.vi.numFrames - 1 && scene_difference(input[0][0], input[2][0]) > threshold;
        left = n == 0 ? right : scene_difference(input[1][0], input[0][0]) > threshold;
        use_temporal = temporal && !left && !right;
      } else if (temporal && d.config.scenechange != 0) {
        const auto* props = api->getFramePropertiesRO(frames[0].get());
        int err = 0;
        const bool left = api->mapGetInt(props, "_SceneChangePrev", 0, &err) != 0 && !err;
        const bool right = api->mapGetInt(props, "_SceneChangeNext", 0, &err) != 0 && !err;
        use_temporal = !left && !right;
      }
      const VSFrame* plane_src[3] = {d.config.process[0] ? nullptr : frames[0].get(),
                                     d.config.process[1] ? nullptr : frames[0].get(),
                                     d.config.process[2] ? nullptr : frames[0].get()};
      const int plane_indices[3] = {0, 1, 2};
      std::unique_ptr<VSFrame, FrameDeleter> dst(
          api->newVideoFrame2(&d.vi.format, d.vi.width, d.vi.height, plane_src, plane_indices, frames[0].get(), core),
          {api});
      if (!dst)
        throw std::bad_alloc();
      if (d.internal_scene) {
        auto* props = api->getFramePropertiesRW(dst.get());
        api->mapSetInt(props, "_SceneChangePrev", left, maReplace);
        api->mapSetInt(props, "_SceneChangeNext", right, maReplace);
      }
      for (int p = 0; p < d.vi.format.numPlanes; ++p) {
        if (!d.config.process[p])
          continue;
        std::array<DeenPlane, 3> planes{};
        for (int i = 0; i < (use_temporal ? 3 : 1); ++i)
          planes[i] = input[i][p];
        deen_process_native(*d.filters[p], d.config.threshold[p], d.config.temporal_threshold[p], use_temporal, planes,
                            api->getWritePtr(dst.get(), p), api->getStride(dst.get(), p));
      }
      return dst.release();
    }
  } catch (const std::exception& e) {
    api->setFilterError(e.what(), ctx);
  } catch (...) {
    api->setFilterError("Deen: frame processing failed.", ctx);
  }
  return nullptr;
}
void VS_CC free_filter(void* data, VSCore*, const VSAPI* api) {
  auto* d = static_cast<Instance*>(data);
  api->freeNode(d->node);
  delete d;
}
} // namespace
void VS_CC deen_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* api) {
  const auto release = [api](VSNode* p) {
    if (p)
      api->freeNode(p);
  };
  std::unique_ptr<VSNode, decltype(release)> node(nullptr, release);
  try {
    int err = 0;
    node.reset(api->mapGetNode(in, "clip", 0, &err));
    require(!err && node, "Deen: clip is required.");
    const auto vi = *api->getVideoInfo(node.get());
    const auto& f = vi.format;
    require((f.colorFamily == cfGray && f.numPlanes == 1) ||
                ((f.colorFamily == cfYUV || f.colorFamily == cfRGB) && f.numPlanes == 3),
            "Deen: unsupported color format.");
    require(vi.width > 0 && vi.height > 0 && vi.numFrames > 0,
            "Deen: constant dimensions and positive length required.");
    require((f.sampleType == stInteger && ((f.bitsPerSample == 8 && f.bytesPerSample == 1) ||
                                           (f.bitsPerSample >= 9 && f.bitsPerSample <= 16 && f.bytesPerSample == 2))) ||
                (f.sampleType == stFloat && ((f.bitsPerSample == 16 && f.bytesPerSample == 2) ||
                                             (f.bitsPerSample == 32 && f.bytesPerSample == 4))),
            "Deen: unsupported sample format.");
    const FormatInfo format{f.colorFamily == cfGray  ? 1
                            : f.colorFamily == cfRGB ? 2
                                                     : 3,
                            f.sampleType == stFloat,
                            f.bitsPerSample,
                            f.bytesPerSample,
                            f.numPlanes,
                            f.subSamplingW,
                            f.subSamplingH};
    auto d = std::make_unique<Instance>(deen_config(deen_parameters(in, api, false), format, false));
    if (d->config.temporal && d->config.scenechange > 0) {
      VSPlugin* misc = api->getPluginByID("com.vapoursynth.misc", core);
      require(vi.numFrames > 1, "Deen: automatic scene detection requires more than one frame");
      if (!misc || !api->getPluginFunctionByName("SCDetect", misc)) {
        d->internal_scene = true;
      } else {
        const auto free_map = [api](VSMap* map) {
          if (map)
            api->freeMap(map);
        };
        std::unique_ptr<VSMap, decltype(free_map)> args(api->createMap(), free_map);
        if (!args)
          throw std::bad_alloc();
        api->mapSetNode(args.get(), "clip", node.get(), maReplace);
        api->mapSetFloat(args.get(), "threshold", d->config.scenechange / 255.0, maReplace);
        std::unique_ptr<VSMap, decltype(free_map)> ret(api->invoke(misc, "SCDetect", args.get()), free_map);
        if (!ret)
          throw std::bad_alloc();
        if (const char* message = api->mapGetError(ret.get()))
          throw std::runtime_error(message);
        VSNode* detected = api->mapGetNode(ret.get(), "clip", 0, &err);
        if (err || !detected) {
          if (detected)
            api->freeNode(detected);
          throw std::runtime_error("Deen: SCDetect returned no clip");
        }
        node.reset(detected);
      }
    }
    d->node = node.get();
    d->vi = vi;
    const VSFilterDependency dep{node.get(), d->config.temporal ? rpGeneral : rpStrictSpatial};
    api->createVideoFilter(out, "Deen", &vi, get_frame, free_filter, fmParallel, &dep, 1, d.get(), core);
    d.release();
    node.release();
  } catch (const std::exception& e) {
    api->mapSetError(out, e.what());
  } catch (...) {
    api->mapSetError(out, "Deen: creation failed.");
  }
}
} // namespace neo_smo::plugin
