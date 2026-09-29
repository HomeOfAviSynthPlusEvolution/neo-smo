#include "plugin/avs_deen.hpp"
#include "algorithms/deen.hpp"
#include "algorithms/mini_deen.hpp"
#include "base/checked.hpp"
#include <dualsynth/avisynth/video_bridge.hpp>
#include <memory>

namespace neo_smo::avs {
namespace {
std::vector<int> integers(const AVSValue& value) {
  std::vector<int> out;
  if (!value.Defined())
    return out;
  const int count = value.IsArray() ? value.ArraySize() : 1;
  for (int i = 0; i < count; ++i) {
    const auto v = value.IsArray() ? value[i] : value;
    require(v.IsInt(), "parameters must contain integers");
    out.push_back(v.AsInt());
  }
  return out;
}
PClip input_clip(const AVSValue& args) {
  require(args[0].IsClip(), "clip is required");
  return args[0].AsClip();
}
class DeenFilter final : public GenericVideoFilter {
public:
  DeenFilter(const AVSValue& args, bool mini) : GenericVideoFilter(input_clip(args)), mini_(mini) {
    require(vi.HasVideo() && vi.width > 0 && vi.height > 0 && vi.num_frames > 0 && vi.IsPlanar() &&
                (vi.NumComponents() == 1 || vi.NumComponents() == 3),
            "only planar Gray, RGB and YUV without alpha are supported");
    const int bits = vi.BitsPerComponent();
    require(bits == 8 || bits == 10 || bits == 12 || bits == 14 || bits == 16 || (!mini && bits == 32),
            "unsupported sample format");
    const auto format = ds::avisynth::make_video_format(vi);
    require(format.has_value(), "unsupported pixel format");
    format_ = format.value();
    const auto& planes = args[mini ? 3 : 10];
    require(!planes.IsArray() || planes.ArraySize() > 0, "planes cannot be empty");
    const auto selected = integers(planes);
    process_.fill(selected.empty());
    for (int p : selected) {
      require(p >= 0 && p < vi.NumComponents() && !process_[p], "planes must contain distinct valid indices");
      process_[p] = true;
    }
    if (mini) {
      const auto parse = [&](int index, std::array<int, 3>& values, int lo, int hi) {
        const auto list = integers(args[index]);
        require(list.size() <= 3, "at most three parameter values are allowed");
        for (std::size_t i = 0; i < list.size(); ++i) {
          require(list[i] >= lo && list[i] <= hi, "invalid radius or threshold");
          values[i] = list[i];
        }
        for (std::size_t i = list.size(); i < 3; ++i)
          if (i > 0)
            values[i] = values[i - 1];
      };
      parse(1, radius_, 1, 7);
      parse(2, threshold_, 0, 255);
      for (int p = 0; p < vi.NumComponents(); ++p)
        process_[p] = process_[p] && threshold_[p] * ((1u << bits) - 1) / 255u > 1;
    } else {
      DeenOptions o;
      o.mode = args[1].AsString("c3d");
      o.radius = args[2].AsInt(1);
      o.spatial_y = args[3].AsFloat(7);
      o.spatial_uv = args[4].AsFloat(9);
      o.temporal_y = args[5].AsFloat(4);
      o.temporal_uv = args[6].AsFloat(6);
      o.minimum = args[7].AsFloat(0.5);
      o.scene_threshold = args[8].AsFloat(9);
      o.scenechange = args[9].AsBool(true);
      deen_ = std::make_unique<Deen>(std::move(o));
    }
  }
  int __stdcall SetCacheHints(int hints, int) override { return hints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0; }
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override {
    try {
      const bool temporal = !mini_ && deen_->temporal() && n > 0 && n < vi.num_frames - 1;
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
      if (temporal && deen_->options().scenechange) {
        const bool left = deen_scene_cut(*deen_, inputs[0], inputs[1]);
        const bool right = deen_scene_cut(*deen_, inputs[0], inputs[2]);
        use_temporal = !left && !right;
      }
      auto dst = env->NewVideoFrameP(vi, &frames[0]);
      for (int p = 0; p < vi.NumComponents(); ++p) {
        const int id = ds::avisynth::plane_id(format_, p);
        const auto& src = inputs[0][p];
        if (!process_[p]) {
          env->BitBlt(dst->GetWritePtr(id), dst->GetPitch(id), frames[0]->GetReadPtr(id), frames[0]->GetPitch(id),
                      frames[0]->GetRowSize(id), frames[0]->GetHeight(id));
        } else if (mini_) {
          mini_deen_process(src.data, src.stride, dst->GetWritePtr(id), dst->GetPitch(id), src.width, src.height,
                            src.bits, radius_[p], threshold_[p]);
        } else {
          std::array<DeenPlane, 3> planes{};
          for (int f = 0; f < (use_temporal ? 3 : 1); ++f)
            planes[f] = inputs[f][p];
          deen_process(*deen_, vi.IsYUV() && !vi.IsY() && p > 0, use_temporal, planes, dst->GetWritePtr(id),
                       dst->GetPitch(id));
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
  std::unique_ptr<Deen> deen_;
  std::array<bool, 3> process_{};
  std::array<int, 3> radius_{1, 1, 1}, threshold_{10, 10, 10};
};
template <bool Mini>
AVSValue __cdecl create(AVSValue args, void*, IScriptEnvironment* env) {
  try {
    env->CheckVersion(11);
    // PClip transfers ownership through the host SDK linkage table.
    // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
    PClip filter = new DeenFilter(args, Mini);
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
  env->AddFunction("neo_smo_MiniDeen", "[clip]c[radius]i*[threshold]i*[planes]i*", create<true>, nullptr);
  env->AddFunction("neo_smo_Deen",
                   "[clip]c[mode]s[rad]i[thrY]f[thrUV]f[tthY]f[tthUV]f[min]f[scd]f[scenechange]b[planes]i*",
                   create<false>, nullptr);
}
} // namespace neo_smo::avs
