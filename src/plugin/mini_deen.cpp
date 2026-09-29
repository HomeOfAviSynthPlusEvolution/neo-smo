#include "plugin/mini_deen.hpp"
#include "algorithms/mini_deen.hpp"
#include "base/checked.hpp"
#include "plugin/deen_params_vs.hpp"
#include <array>
#include <memory>
namespace neo_smo::plugin {
namespace {
struct Instance {
  VSNode* node = nullptr;
  VSVideoInfo vi{};
  DeenConfig config;
};
struct FrameDeleter {
  const VSAPI* api;
  void operator()(const VSFrame* frame) const {
    if (frame)
      api->freeFrame(frame);
  }
};
const VSFrame* VS_CC get_frame(int n, int activation, void* data, void**, VSFrameContext* ctx, VSCore* core,
                               const VSAPI* api) {
  const auto& d = *static_cast<const Instance*>(data);
  try {
    if (activation == arInitial) {
      api->requestFrameFilter(n, d.node, ctx);
    } else if (activation == arAllFramesReady) {
      std::unique_ptr<const VSFrame, FrameDeleter> src(api->getFrameFilter(n, d.node, ctx), {api});
      if (!src)
        return nullptr;
      const auto& actual = *api->getVideoFrameFormat(src.get());
      const auto& f = d.vi.format;
      require(actual.colorFamily == f.colorFamily && actual.sampleType == f.sampleType &&
                  actual.bitsPerSample == f.bitsPerSample && actual.bytesPerSample == f.bytesPerSample &&
                  actual.numPlanes == f.numPlanes && actual.subSamplingW == f.subSamplingW &&
                  actual.subSamplingH == f.subSamplingH,
              "MiniDeen: input format changed.");
      for (int p = 0; p < f.numPlanes; ++p) {
        const bool uv = f.colorFamily == cfYUV && p > 0;
        require(api->getFrameWidth(src.get(), p) == (d.vi.width >> (uv ? f.subSamplingW : 0)) &&
                    api->getFrameHeight(src.get(), p) == (d.vi.height >> (uv ? f.subSamplingH : 0)),
                "MiniDeen: input dimensions changed.");
      }
      const VSFrame* planes[3] = {d.config.process[0] ? nullptr : src.get(), d.config.process[1] ? nullptr : src.get(),
                                  d.config.process[2] ? nullptr : src.get()};
      const int indices[3] = {0, 1, 2};
      std::unique_ptr<VSFrame, FrameDeleter> dst(
          api->newVideoFrame2(&f, d.vi.width, d.vi.height, planes, indices, src.get(), core), {api});
      if (!dst)
        throw std::bad_alloc();
      for (int p = 0; p < f.numPlanes; ++p)
        if (d.config.process[p])
          mini_deen_process_native(api->getReadPtr(src.get(), p), api->getStride(src.get(), p),
                                   api->getWritePtr(dst.get(), p), api->getStride(dst.get(), p),
                                   api->getFrameWidth(src.get(), p), api->getFrameHeight(src.get(), p),
                                   get_data_type(f.bytesPerSample, f.sampleType == stFloat), f.bitsPerSample,
                                   d.config.radius[p], d.config.threshold[p]);
      return dst.release();
    }
  } catch (const std::exception& e) {
    api->setFilterError(e.what(), ctx);
  } catch (...) {
    api->setFilterError("MiniDeen: frame processing failed.", ctx);
  }
  return nullptr;
}
void VS_CC free_filter(void* data, VSCore*, const VSAPI* api) {
  auto* d = static_cast<Instance*>(data);
  api->freeNode(d->node);
  delete d;
}
} // namespace
void VS_CC mini_deen_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* api) {
  const auto release = [api](VSNode* p) {
    if (p)
      api->freeNode(p);
  };
  std::unique_ptr<VSNode, decltype(release)> node(nullptr, release);
  try {
    int err = 0;
    node.reset(api->mapGetNode(in, "clip", 0, &err));
    require(!err && node, "MiniDeen: clip is required.");
    const auto& vi = *api->getVideoInfo(node.get());
    const auto& f = vi.format;
    require((f.colorFamily == cfGray && f.numPlanes == 1) ||
                ((f.colorFamily == cfYUV || f.colorFamily == cfRGB) && f.numPlanes == 3),
            "MiniDeen: constant GRAY, YUV or RGB format required.");
    require(vi.width > 0 && vi.height > 0 && vi.numFrames > 0, "MiniDeen: constant positive dimensions required.");
    require((f.sampleType == stInteger && ((f.bitsPerSample == 8 && f.bytesPerSample == 1) ||
                                           (f.bitsPerSample >= 9 && f.bitsPerSample <= 16 && f.bytesPerSample == 2))) ||
                (f.sampleType == stFloat && ((f.bitsPerSample == 16 && f.bytesPerSample == 2) ||
                                             (f.bitsPerSample == 32 && f.bytesPerSample == 4))),
            "MiniDeen: unsupported sample format");
    auto d = std::make_unique<Instance>();
    const FormatInfo format{f.colorFamily == cfGray  ? 1
                            : f.colorFamily == cfRGB ? 2
                                                     : 3,
                            f.sampleType == stFloat,
                            f.bitsPerSample,
                            f.bytesPerSample,
                            f.numPlanes,
                            f.subSamplingW,
                            f.subSamplingH};
    d->config = deen_config(deen_parameters(in, api, true), format, true);
    d->vi = vi;
    d->node = node.get();
    const VSFilterDependency dep{node.get(), rpStrictSpatial};
    api->createVideoFilter(out, "MiniDeen", &vi, get_frame, free_filter, fmParallel, &dep, 1, d.get(), core);
    d.release();
    node.release();
  } catch (const std::exception& e) {
    api->mapSetError(out, e.what());
  } catch (...) {
    api->mapSetError(out, "MiniDeen: creation failed.");
  }
}
} // namespace neo_smo::plugin
