#include "plugin/dct_filter.hpp"
#include "algorithms/dct_filter.hpp"
#include "base/checked.hpp"

namespace neo_smo::plugin {
namespace {
struct DctInstance {
  explicit DctInstance(const std::vector<double>& factors) : filter(factors) {}
  VSNode* node = nullptr;
  VSVideoInfo vi{};
  std::array<bool, 3> process{true, true, true};
  DctFilter filter;
};

const VSFrame* VS_CC get_frame(int n, int activation, void* data, void**, VSFrameContext* ctx, VSCore* core,
                               const VSAPI* api) {
  const auto& d = *static_cast<const DctInstance*>(data);
  try {
    if (activation == arInitial) {
      api->requestFrameFilter(n, d.node, ctx);
    } else if (activation == arAllFramesReady) {
      const auto free_frame = [api](const VSFrame* frame) {
        api->freeFrame(frame);
      };
      std::unique_ptr<const VSFrame, decltype(free_frame)> src(api->getFrameFilter(n, d.node, ctx), free_frame);
      if (!src)
        return nullptr;
      const VSFrame* plane_src[3] = {d.process[0] ? nullptr : src.get(), d.process[1] ? nullptr : src.get(),
                                     d.process[2] ? nullptr : src.get()};
      const int planes[3] = {0, 1, 2};
      std::unique_ptr<VSFrame, decltype(free_frame)> dst(
          api->newVideoFrame2(&d.vi.format, d.vi.width, d.vi.height, plane_src, planes, src.get(), core), free_frame);
      if (!dst)
        throw std::bad_alloc();
      DctScratch scratch(d.filter, static_cast<std::size_t>(d.vi.width));
      for (int p = 0; p < d.vi.format.numPlanes; ++p) {
        if (!d.process[static_cast<std::size_t>(p)])
          continue;
        const auto src_stride = api->getStride(src.get(), p);
        const auto dst_stride = api->getStride(dst.get(), p);
        require(src_stride > 0 && dst_stride > 0, "DCTFilter: unsupported plane stride.");
        scratch.process(d.filter, get_data_type(d.vi.format.bytesPerSample, d.vi.format.sampleType == stFloat),
                        d.vi.format.bitsPerSample, api->getReadPtr(src.get(), p), api->getWritePtr(dst.get(), p),
                        static_cast<std::size_t>(api->getFrameWidth(src.get(), p)),
                        static_cast<std::size_t>(api->getFrameHeight(src.get(), p)),
                        static_cast<std::size_t>(src_stride), static_cast<std::size_t>(dst_stride));
      }
      return dst.release();
    }
  } catch (const std::exception& e) {
    api->setFilterError(e.what(), ctx);
  } catch (...) {
    api->setFilterError("DCTFilter: frame processing failed.", ctx);
  }
  return nullptr;
}

void VS_CC free_filter(void* data, VSCore*, const VSAPI* api) {
  auto* d = static_cast<DctInstance*>(data);
  api->freeNode(d->node);
  delete d;
}
} // namespace

void VS_CC dct_filter_create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* api) {
  const auto free_node = [api](VSNode* node) {
    api->freeNode(node);
  };
  std::unique_ptr<VSNode, decltype(free_node)> node(nullptr, free_node);
  try {
    int err = 0;
    node.reset(api->mapGetNode(in, "clip", 0, &err));
    require(!err && node, "DCTFilter: clip is required.");
    const auto* vi = api->getVideoInfo(node.get());
    const auto& fmt = vi->format;
    require((fmt.colorFamily == cfGray || fmt.colorFamily == cfRGB || fmt.colorFamily == cfYUV) && vi->width > 0 &&
                vi->height > 0 && fmt.numPlanes >= 1 && fmt.numPlanes <= 3,
            "DCTFilter: only constant format and dimensions supported.");
    require((fmt.sampleType == stInteger &&
             ((fmt.bytesPerSample == 1 && fmt.bitsPerSample == 8) ||
              (fmt.bytesPerSample == 2 && fmt.bitsPerSample >= 9 && fmt.bitsPerSample <= 16))) ||
                (fmt.sampleType == stFloat && ((fmt.bytesPerSample == 2 && fmt.bitsPerSample == 16) ||
                                               (fmt.bytesPerSample == 4 && fmt.bitsPerSample == 32))),
            "DCTFilter: unsupported sample format.");
    require(api->mapNumElements(in, "factors") == 8, "DCTFilter: exactly eight factors are required.");
    std::vector<double> factors(8);
    for (int i = 0; i < 8; ++i) {
      factors[static_cast<std::size_t>(i)] = api->mapGetFloat(in, "factors", i, &err);
      require(!err, "DCTFilter: invalid factors.");
    }
    auto d = std::make_unique<DctInstance>(factors);
    const int count = api->mapNumElements(in, "planes");
    if (count >= 0) {
      d->process.fill(false);
      for (int i = 0; i < count; ++i) {
        const auto plane = api->mapGetInt(in, "planes", i, &err);
        require(!err && plane >= 0 && plane < fmt.numPlanes, "DCTFilter: plane index out of range.");
        const auto p = static_cast<std::size_t>(plane);
        require(!d->process[p], "DCTFilter: plane specified twice.");
        d->process[p] = true;
      }
    }
    d->node = node.get();
    d->vi = *vi;
    const VSFilterDependency dep{node.get(), rpStrictSpatial};
    api->createVideoFilter(out, "DCTFilter", vi, get_frame, free_filter, fmParallel, &dep, 1, d.get(), core);
    d.release();
    node.release();
  } catch (const std::exception& e) {
    api->mapSetError(out, e.what());
  } catch (...) {
    api->mapSetError(out, "DCTFilter: creation failed.");
  }
}
} // namespace neo_smo::plugin
