#include "plugin/deen.hpp"
#include "algorithms/deen.hpp"
#include "base/checked.hpp"
#include <memory>

namespace neo_smo::plugin {
namespace {
struct Instance {
  explicit Instance(DeenOptions o) : filter(std::move(o)) {}
  Deen filter;
  VSNode* node = nullptr;
  VSVideoInfo vi{};
  std::array<bool, 3> process{true, true, true};
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
  const bool temporal = d.filter.temporal() && n > 0 && n < d.vi.numFrames - 1;
  const int count = temporal ? 3 : 1;
  const int indices[3] = {n, n - 1, n + 1};
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
      if (temporal && d.filter.options().scenechange) {
        // Evaluate both sides even if the first is a cut: all required samples must be valid.
        const bool left = deen_scene_cut(d.filter, input[0], input[1]);
        const bool right = deen_scene_cut(d.filter, input[0], input[2]);
        use_temporal = !left && !right;
      }
      const VSFrame* plane_src[3] = {d.process[0] ? nullptr : frames[0].get(), d.process[1] ? nullptr : frames[0].get(),
                                     d.process[2] ? nullptr : frames[0].get()};
      const int plane_indices[3] = {0, 1, 2};
      std::unique_ptr<VSFrame, FrameDeleter> dst(
          api->newVideoFrame2(&d.vi.format, d.vi.width, d.vi.height, plane_src, plane_indices, frames[0].get(), core),
          {api});
      if (!dst)
        throw std::bad_alloc();
      for (int p = 0; p < d.vi.format.numPlanes; ++p) {
        if (!d.process[p])
          continue;
        std::array<DeenPlane, 3> planes{};
        for (int i = 0; i < (use_temporal ? 3 : 1); ++i)
          planes[i] = input[i][p];
        deen_process(d.filter, d.vi.format.colorFamily == cfYUV && p > 0, use_temporal, planes,
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
    const auto& vi = *api->getVideoInfo(node.get());
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
    const auto number = [&](const char* key, double fallback) {
      const int count = api->mapNumElements(in, key);
      if (count < 0)
        return fallback;
      require(count == 1, "Deen: expected one numeric value.");
      const double v = api->mapGetFloat(in, key, 0, &err);
      require(!err, "Deen: invalid numeric parameter.");
      return v;
    };
    const auto integer = [&](const char* key, std::int64_t fallback) {
      const int count = api->mapNumElements(in, key);
      if (count < 0)
        return fallback;
      require(count == 1, "Deen: expected one integer value.");
      const auto v = api->mapGetInt(in, key, 0, &err);
      require(!err, "Deen: invalid integer parameter.");
      return v;
    };
    DeenOptions o;
    if (api->mapNumElements(in, "mode") >= 0) {
      require(api->mapNumElements(in, "mode") == 1, "Deen: expected one mode.");
      const auto* mode = api->mapGetData(in, "mode", 0, &err);
      require(!err && mode, "Deen: invalid mode.");
      const int size = api->mapGetDataSize(in, "mode", 0, &err);
      require(!err && size == 3, "Deen: invalid mode.");
      o.mode.assign(mode, static_cast<std::size_t>(size));
    }
    const auto r = integer("rad", 1);
    require(r >= 1 && r <= 7, "Deen: invalid radius.");
    o.radius = static_cast<int>(r);
    o.spatial_y = number("thrY", 7);
    o.spatial_uv = number("thrUV", 9);
    o.temporal_y = number("tthY", 4);
    o.temporal_uv = number("tthUV", 6);
    o.minimum = number("min", 0.5);
    o.scene_threshold = number("scd", 9);
    const auto scene = integer("scenechange", 1);
    require(scene == 0 || scene == 1, "Deen: scenechange must be boolean.");
    o.scenechange = scene != 0;
    auto d = std::make_unique<Instance>(std::move(o));
    const int count = api->mapNumElements(in, "planes");
    require(count != 0, "Deen: planes cannot be empty.");
    if (count > 0) {
      d->process.fill(false);
      for (int i = 0; i < count; ++i) {
        const auto p = api->mapGetInt(in, "planes", i, &err);
        require(!err && p >= 0 && p < f.numPlanes, "Deen: invalid plane index.");
        require(!d->process[static_cast<std::size_t>(p)], "Deen: duplicate plane index.");
        d->process[static_cast<std::size_t>(p)] = true;
      }
    }
    d->node = node.get();
    d->vi = vi;
    const VSFilterDependency dep{node.get(), d->filter.temporal() ? rpGeneral : rpStrictSpatial};
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
