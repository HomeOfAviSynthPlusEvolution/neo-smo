#include "algorithms/f3kdb_views.hpp"
#include "algorithms/f3kdb.hpp"
#include "plugin/f3kdb.hpp"
#include "plugin/f3kdb_params.hpp"
#include <cstring>
#include <memory>

namespace neo_smo::plugin {
namespace {
namespace f = neo_smo::f3kdb;
struct Reader {
  const VSMap* in;
  const VSAPI* api;
  bool present(std::size_t, const std::string& k) const {
    const int count = api->mapNumElements(in, k.c_str());
    require(count < 0 || count == 1, "F3kdb: expected a single parameter value");
    return count == 1;
  }
  std::int64_t integer(std::size_t, const std::string& k, bool) const {
    int e = 0;
    const auto v = api->mapGetInt(in, k.c_str(), 0, &e);
    require(!e, "F3kdb: invalid integer");
    return v;
  }
  double number(std::size_t, const std::string& k) const {
    int e = 0;
    const auto v = api->mapGetFloat(in, k.c_str(), 0, &e);
    require(!e, "F3kdb: invalid number");
    return v;
  }
  std::string string(std::size_t, const std::string& k) const {
    int e = 0;
    const char* v = api->mapGetData(in, k.c_str(), 0, &e);
    require(!e && v, "F3kdb: invalid string");
    const int size = api->mapGetDataSize(in, k.c_str(), 0, &e);
    require(!e && size >= 0, "F3kdb: invalid string length");
    return {v, static_cast<std::size_t>(size)};
  }
};
struct Instance {
  VSNode* node = nullptr;
  VSVideoInfo input{}, output{};
  std::unique_ptr<neo_smo::f3kdb::Filter> filter;
};
struct FrameDeleter {
  const VSAPI* api;
  void operator()(const VSFrame* f) const {
    if (f)
      api->freeFrame(f);
  }
};
const VSFrame* VS_CC get_frame(int n, int activation, void* data, void**, VSFrameContext* ctx, VSCore* core,
                               const VSAPI* api) {
  const auto& d = *static_cast<const Instance*>(data);
  try {
    if (activation == arInitial)
      api->requestFrameFilter(n, d.node, ctx);
    else if (activation == arAllFramesReady) {
      std::unique_ptr<const VSFrame, FrameDeleter> src(api->getFrameFilter(n, d.node, ctx), {api});
      if (!src)
        return nullptr;
      const auto& actual = *api->getVideoFrameFormat(src.get());
      const auto& fmt = d.input.format;
      require(actual.colorFamily == fmt.colorFamily && actual.sampleType == fmt.sampleType &&
                  actual.bitsPerSample == fmt.bitsPerSample && actual.bytesPerSample == fmt.bytesPerSample &&
                  actual.subSamplingW == fmt.subSamplingW && actual.subSamplingH == fmt.subSamplingH &&
                  actual.numPlanes == fmt.numPlanes,
              "F3kdb: input format changed");
      std::vector<f::PlaneView> input;
      for (int p = 0; p < fmt.numPlanes; ++p) {
        const int w = d.input.width >> (p ? fmt.subSamplingW : 0), h = d.input.height >> (p ? fmt.subSamplingH : 0);
        const auto stride = api->getStride(src.get(), p);
        require(api->getFrameWidth(src.get(), p) == w && api->getFrameHeight(src.get(), p) == h && stride > 0,
                "F3kdb: invalid input dimensions or stride");
        require(static_cast<std::size_t>(stride) <=
                    std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(h),
                "F3kdb: plane size overflow");
        input.push_back(
            {api->getReadPtr(src.get(), p), static_cast<std::size_t>(stride) * h, static_cast<std::size_t>(stride)});
      }
      std::unique_ptr<VSFrame, FrameDeleter> dst(
          api->newVideoFrame(&d.output.format, d.output.width, d.output.height, src.get(), core), {api});
      if (!dst)
        throw std::bad_alloc();
      std::vector<f::MutablePlaneView> output;
      for (int p = 0; p < d.output.format.numPlanes; ++p) {
        const auto stride = api->getStride(dst.get(), p);
        const auto height = d.output.height >> (p ? d.output.format.subSamplingH : 0);
        require(stride > 0 && static_cast<std::size_t>(stride) <= std::numeric_limits<std::size_t>::max() / height,
                "F3kdb: invalid output stride");
        output.push_back({api->getWritePtr(dst.get(), p), static_cast<std::size_t>(stride) * height,
                          static_cast<std::size_t>(stride)});
      }
      d.filter->process_into(static_cast<std::size_t>(n), input, output);
      return dst.release();
    }
  } catch (const std::exception& e) {
    api->setFilterError(e.what(), ctx);
  } catch (...) {
    api->setFilterError("F3kdb: processing failed", ctx);
  }
  return nullptr;
}
void VS_CC free_filter(void* data, VSCore*, const VSAPI* api) {
  auto* d = static_cast<Instance*>(data);
  api->freeNode(d->node);
  delete d;
}
} // namespace
void VS_CC f3kdb_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* api) {
  const auto release = [api](VSNode* p) {
    if (p)
      api->freeNode(p);
  };
  std::unique_ptr<VSNode, decltype(release)> node(nullptr, release);
  try {
    int e = 0;
    node.reset(api->mapGetNode(in, "clip", 0, &e));
    require(!e && node, "F3kdb: clip is required");
    const auto vi = *api->getVideoInfo(node.get());
    const auto& v = vi.format;
    require(vi.width > 0 && vi.height > 0 && vi.numFrames > 0, "F3kdb: constant positive dimensions required");
    require((v.colorFamily == cfGray && v.numPlanes == 1) || (v.colorFamily == cfYUV && v.numPlanes == 3),
            "F3kdb: Gray or YUV without alpha required");
    f::Format fmt;
    fmt.color = v.colorFamily == cfGray ? f::Color::gray : f::Color::yuv;
    fmt.sample = v.sampleType == stFloat ? f::Sample::floating : f::Sample::integer;
    fmt.depth = v.bitsPerSample;
    fmt.width = vi.width;
    fmt.height = vi.height;
    fmt.frames = vi.numFrames;
    fmt.sub_x = v.subSamplingW;
    fmt.sub_y = v.subSamplingH;
    auto d = std::make_unique<Instance>();
    d->input = vi;
    d->output = vi;
    d->filter = std::make_unique<neo_smo::f3kdb::Filter>(fmt, f3kdb_options(Reader{in, api}));
    const auto& output = d->filter->parameters().output;
    require(api->queryVideoFormat(&d->output.format, v.colorFamily,
                                  output.sample == f::Sample::floating ? stFloat : stInteger, output.depth,
                                  v.subSamplingW, v.subSamplingH, core) != 0,
            "F3kdb: output format unavailable in host");
    d->node = node.get();
    const VSFilterDependency dep{node.get(), rpStrictSpatial};
    api->createVideoFilter(out, "F3kdb", &d->output, get_frame, free_filter, fmParallel, &dep, 1, d.get(), core);
    d.release();
    node.release();
  } catch (const std::exception& e) {
    api->mapSetError(out, e.what());
  } catch (...) {
    api->mapSetError(out, "F3kdb: creation failed");
  }
}
} // namespace neo_smo::plugin
