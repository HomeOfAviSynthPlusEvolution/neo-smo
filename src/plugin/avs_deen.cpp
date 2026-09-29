#include "plugin/avs_deen.hpp"
#include "algorithms/deen.hpp"
#include "algorithms/mini_deen.hpp"
#include "base/checked.hpp"
#include "algorithms/deen_config.hpp"
#include "plugin/deen_descriptor.hpp"
#include "plugin/avs_filter.hpp"
#include <dualsynth/avisynth/video_bridge.hpp>
#include <memory>

namespace neo_smo::avs {
namespace {
DeenParameters parameters(const AVSValue& args, bool mini) {
  DeenParameters p;
  const auto descriptor = plugin::deen_descriptor(mini);
  for (std::size_t i = 1; i < descriptor.params.size(); ++i) {
    const auto value = args[static_cast<int>(i)];
    if (!value.Defined())
      continue;
    const auto& spec = descriptor.params[i];
    if (spec.is_array) {
      const int count = value.IsArray() ? value.ArraySize() : 1;
      if (spec.type == ds::ParamType::Integer) {
        auto& target = spec.name == "radius" ? p.radius : p.planes;
        target.emplace();
        for (int j = 0; j < count; ++j) {
          const auto v = value.IsArray() ? value[j] : value;
          require(v.IsInt(), spec.name + " requires integers");
          target->push_back(v.AsInt());
        }
      } else {
        auto& target = spec.name == "threshold" ? p.threshold
                       : spec.name == "minimum" ? p.minimum
                                                : p.temporal_threshold;
        target.emplace();
        for (int j = 0; j < count; ++j) {
          const auto v = value.IsArray() ? value[j] : value;
          require(v.IsFloat(), spec.name + " requires numbers");
          target->push_back(v.AsFloat());
        }
      }
    } else if (spec.name == "mode") {
      require(value.IsString(), "mode requires a string");
      p.mode = value.AsString();
    } else if (spec.name == "scalep") {
      require(value.IsBool(), "scalep requires a boolean");
      p.scalep = value.AsBool();
    } else {
      require(value.IsInt(), "scenechange requires an integer");
      p.scenechange = value.AsInt();
    }
  }
  return p;
}
PClip input_clip(const AVSValue& args) {
  require(args[0].IsClip(), "clip is required");
  return args[0].AsClip();
}
class DeenFilter final : public GenericVideoFilter {
public:
  DeenFilter(const AVSValue& args, bool mini, IScriptEnvironment* env)
      : GenericVideoFilter(input_clip(args)), mini_(mini) {
    require(vi.HasVideo() && vi.width > 0 && vi.height > 0 && vi.num_frames > 0 && vi.IsPlanar() &&
                (vi.NumComponents() == 1 || vi.NumComponents() == 3),
            "only planar Gray, RGB and YUV without alpha are supported");
    const int bits = vi.BitsPerComponent();
    require(bits == 8 || bits == 10 || bits == 12 || bits == 14 || bits == 16 || bits == 32,
            "unsupported sample format");
    const auto format = ds::avisynth::make_video_format(vi);
    require(format.has_value(), "unsupported pixel format");
    format_ = format.value();
    const FormatInfo fmt{
        vi.IsY() ? 1 : vi.IsRGB() ? 2 : 3, bits == 32, bits, vi.ComponentSize(), vi.NumComponents(), 0, 0};
    config_ = deen_config(parameters(args, mini), fmt, mini);
    if (!mini) {
      for (int p = 0; p < vi.NumComponents(); ++p)
        deen_[p] = std::make_unique<Deen>(deen_plane_filter(config_, p));
      if (config_.temporal && config_.scenechange > 0)
        child = scene_detect(child, config_.scenechange / 255.0, env);
    }
  }
  int __stdcall SetCacheHints(int hints, int) override { return hints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0; }
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override {
    try {
      const bool temporal = config_.temporal && n > 0 && n < vi.num_frames - 1;
      const int count = temporal ? 3 : 1;
      const int indices[3] = {n, n - 1, n + 1};
      std::array<PVideoFrame, 3> frames;
      std::array<std::vector<DeenPlane>, 3> inputs;
      for (int f = 0; f < count; ++f) {
        frames[f] = child->GetFrame(indices[f], env);
        require(!!frames[f], "failed to retrieve input frame");
        inputs[f] = views(frames[f]);
      }
      bool use_temporal = temporal;
      if (temporal && config_.scenechange != 0) {
        const bool left = scene_property(frames[0], "_SceneChangePrev", env) != 0;
        const bool right = scene_property(frames[0], "_SceneChangeNext", env) != 0;
        use_temporal = !left && !right;
      }
      auto dst = env->NewVideoFrameP(vi, &frames[0]);
      for (int p = 0; p < vi.NumComponents(); ++p) {
        const int id = ds::avisynth::plane_id(format_, p);
        const auto& src = inputs[0][p];
        if (!config_.process[p]) {
          env->BitBlt(dst->GetWritePtr(id), dst->GetPitch(id), frames[0]->GetReadPtr(id), frames[0]->GetPitch(id),
                      frames[0]->GetRowSize(id), frames[0]->GetHeight(id));
        } else if (mini_) {
          mini_deen_process_native(src.data, src.stride, dst->GetWritePtr(id), dst->GetPitch(id), src.width, src.height,
                                   src.type, src.bits, config_.radius[p], config_.threshold[p]);
        } else {
          std::array<DeenPlane, 3> planes{};
          for (int f = 0; f < (use_temporal ? 3 : 1); ++f)
            planes[f] = inputs[f][p];
          deen_process_native(*deen_[p], config_.threshold[p], config_.temporal_threshold[p], use_temporal, planes,
                              dst->GetWritePtr(id), dst->GetPitch(id));
        }
      }
      return dst;
    } catch (const AvisynthError&) {
      throw;
    } catch (const std::exception& e) {
      env->ThrowError("neo_smo_%s: %s", mini_ ? "MiniDeen" : "Deen", e.what());
    } catch (...) {
      env->ThrowError("neo_smo_%s: frame processing failed", mini_ ? "MiniDeen" : "Deen");
    }
    return {};
  }

private:
  std::vector<DeenPlane> views(const PVideoFrame& frame) const {
    std::vector<DeenPlane> out;
    for (int p = 0; p < vi.NumComponents(); ++p) {
      const int id = ds::avisynth::plane_id(format_, p);
      const bool uv = vi.IsYUV() && !vi.IsY() && p > 0;
      const int width = vi.width >> (uv ? vi.GetPlaneWidthSubsampling(id) : 0);
      const int height = vi.height >> (uv ? vi.GetPlaneHeightSubsampling(id) : 0);
      require(frame->GetRowSize(id) / vi.ComponentSize() == width && frame->GetHeight(id) == height,
              "input dimensions changed");
      out.push_back({frame->GetReadPtr(id), frame->GetPitch(id), width, height,
                     get_data_type(vi.ComponentSize(), vi.BitsPerComponent() == 32), vi.BitsPerComponent()});
    }
    return out;
  }
  bool mini_;
  ds::VideoFormat format_{};
  DeenConfig config_;
  std::array<std::unique_ptr<Deen>, 3> deen_;
};
template <bool Mini>
AVSValue __cdecl create(AVSValue args, void*, IScriptEnvironment* env) {
  try {
    env->CheckVersion(11);
    // PClip transfers ownership through the host SDK linkage table.
    // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
    PClip filter = new DeenFilter(args, Mini, env);
    return AVSValue(filter);
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& e) {
    env->ThrowError("neo_smo_%s: %s", Mini ? "MiniDeen" : "Deen", e.what());
  } catch (...) {
    env->ThrowError("neo_smo_%s: creation failed", Mini ? "MiniDeen" : "Deen");
  }
  return {};
}
} // namespace
void add_deen(IScriptEnvironment* env) {
  for (bool mini : {false, true}) {
    auto descriptor = plugin::deen_descriptor(mini);
    for (auto& p : descriptor.params)
      p.required = false;
    const auto result = ds::make_avisynth_signature(descriptor);
    require(result.has_value(), result.has_value() ? "" : result.error().message);
    const auto name = std::string("neo_smo_") + descriptor.name;
    env->AddFunction(env->SaveString(name.c_str()), env->SaveString(result.value().c_str()),
                     mini ? create<true> : create<false>, nullptr);
  }
}
} // namespace neo_smo::avs
