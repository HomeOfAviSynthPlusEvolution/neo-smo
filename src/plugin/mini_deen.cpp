#include "plugin/mini_deen.hpp"
#include "algorithms/mini_deen.hpp"
#include "base/checked.hpp"
#include <array>
#include <memory>
namespace neo_smo::plugin {
namespace {
struct Instance {
  VSNode* node = nullptr;
  VSVideoInfo vi{};
  std::array<int, 3> radius{1, 1, 1}, threshold{10, 10, 10};
  std::array<bool, 3> process{true, true, true};
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
      const VSFrame* planes[3] = {d.process[0] ? nullptr : src.get(), d.process[1] ? nullptr : src.get(),
                                  d.process[2] ? nullptr : src.get()};
      const int indices[3] = {0, 1, 2};
      std::unique_ptr<VSFrame, FrameDeleter> dst(
          api->newVideoFrame2(&f, d.vi.width, d.vi.height, planes, indices, src.get(), core), {api});
      if (!dst)
        throw std::bad_alloc();
      for (int p = 0; p < f.numPlanes; ++p)
        if (d.process[p])
          mini_deen_process(api->getReadPtr(src.get(), p), api->getStride(src.get(), p), api->getWritePtr(dst.get(), p),
                            api->getStride(dst.get(), p), api->getFrameWidth(src.get(), p),
                            api->getFrameHeight(src.get(), p), f.bitsPerSample, d.radius[p], d.threshold[p]);
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
    require(f.sampleType == stInteger && ((f.bitsPerSample == 8 && f.bytesPerSample == 1) ||
                                          (f.bitsPerSample >= 9 && f.bitsPerSample <= 16 && f.bytesPerSample == 2)),
            "MiniDeen: only 8..16 bit integer clips are supported.");
    auto d = std::make_unique<Instance>();
    const auto read_array = [&](const char* key, std::array<int, 3>& values, int low, int high) {
      const int count = api->mapNumElements(in, key);
      require(count <= 3, "MiniDeen: at most three parameter values are allowed.");
      for (int p = 0; p < 3; ++p) {
        if (p < count) {
          const auto value = api->mapGetInt(in, key, p, &err);
          require(!err && value >= low && value <= high, std::string("MiniDeen: invalid ") + key + ".");
          values[p] = static_cast<int>(value);
        } else if (p > 0)
          values[p] = values[p - 1];
      }
    };
    read_array("radius", d->radius, 1, 7);
    read_array("threshold", d->threshold, 0, 255);
    const int count = api->mapNumElements(in, "planes");
    if (count > 0) {
      d->process.fill(false);
      for (int i = 0; i < count; ++i) {
        const auto p = api->mapGetInt(in, "planes", i, &err);
        require(!err && p >= 0 && p < f.numPlanes, "MiniDeen: invalid plane index.");
        require(!d->process[static_cast<std::size_t>(p)], "MiniDeen: duplicate plane index.");
        d->process[static_cast<std::size_t>(p)] = true;
      }
    }
    for (int p = 0; p < f.numPlanes; ++p)
      d->process[p] = d->process[p] && (d->threshold[p] * ((1u << f.bitsPerSample) - 1) / 255u > 1);
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
