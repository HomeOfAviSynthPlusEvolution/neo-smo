#include "plugin/avs_filter.hpp"
#include "algorithms/scene_detection.hpp"
#include <algorithm>

namespace neo_smo::avs {
std::int64_t scene_property(const PVideoFrame& f, const char* key, IScriptEnvironment* env, bool required) {
  int err = 0;
  const auto value = env->propGetInt(env->getFramePropsRO(f), key, 0, &err);
  require(!required || !err, "scene change handling requires _SceneChangePrev and _SceneChangeNext frame properties");
  return err ? 0 : value;
}
namespace {
class SceneDetect final : public GenericVideoFilter {
public:
  SceneDetect(PClip clip, double threshold) : GenericVideoFilter(clip), threshold_(threshold) {
    require(!vi.IsRGB(), "automatic scene detection requires Gray or YUV; for RGB supply scene properties explicitly");
    require(vi.num_frames > 1, "automatic scene detection requires more than one frame");
  }
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override {
    const auto src = child->GetFrame(n, env);
    // Match SCDetect's shifted PlaneStats stream, including its first frame.
    const bool prev = difference(std::max(n - 1, 0), env) > threshold_;
    const bool next = difference(n, env) > threshold_;
    auto dst = env->NewVideoFrameP(vi, &src);
    for (int id : {PLANAR_Y, PLANAR_U, PLANAR_V}) {
      if (id != PLANAR_Y && vi.IsY())
        break;
      env->BitBlt(dst->GetWritePtr(id), dst->GetPitch(id), src->GetReadPtr(id), src->GetPitch(id), src->GetRowSize(id),
                  src->GetHeight(id));
    }
    auto* props = env->getFramePropsRW(dst);
    env->propSetInt(props, "_SceneChangePrev", prev, PROPAPPENDMODE_REPLACE);
    env->propSetInt(props, "_SceneChangeNext", next, PROPAPPENDMODE_REPLACE);
    return dst;
  }
  int __stdcall SetCacheHints(int hints, int) override { return hints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0; }

private:
  double threshold_;
  double difference(int n, IScriptEnvironment* env) const {
    if (n == vi.num_frames - 1)
      return 0;
    auto a = child->GetFrame(n, env), b = child->GetFrame(n + 1, env);
    const auto plane = [&](const PVideoFrame& f) {
      return DeenPlane{f->GetReadPtr(PLANAR_Y),
                       f->GetPitch(PLANAR_Y),
                       vi.width,
                       vi.height,
                       get_data_type(vi.ComponentSize(), vi.BitsPerComponent() == 32),
                       vi.BitsPerComponent()};
    };
    return scene_difference(plane(a), plane(b));
  }
};
} // namespace
PClip scene_detect(PClip clip, double threshold, IScriptEnvironment* env) {
  return env->Invoke("Cache", AVSValue(new SceneDetect(clip, threshold))).AsClip();
}
PClip chroma_luma(PClip clip, const FormatInfo& format, IScriptEnvironment* env) {
  auto y = env->Invoke("ExtractY", AVSValue(clip)).AsClip();
  if (format.subsampling_w || format.subsampling_h) {
    const AVSValue args[]{y, format.width >> format.subsampling_w, format.height >> format.subsampling_h};
    y = env->Invoke("BilinearResize", AVSValue(args, 3)).AsClip();
  }
  return y;
}
} // namespace neo_smo::avs
