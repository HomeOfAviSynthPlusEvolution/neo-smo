#include "algorithms/f3kdb.hpp"
#include "plugin/avs_f3kdb.hpp"
#include "plugin/f3kdb_params.hpp"
#include <dualsynth/avisynth/video_bridge.hpp>
#include <memory>

namespace neo_smo::avs {
namespace {
namespace f = neo_smo::f3kdb;
struct Reader {
  const AVSValue& args;
  bool present(std::size_t i, const std::string&) const { return args[static_cast<int>(i)].Defined(); }
  std::int64_t integer(std::size_t i, const std::string&, bool boolean) const {
    const auto v = args[static_cast<int>(i)];
    require(boolean ? v.IsBool() : v.IsInt(), "F3kdb: invalid integer/boolean");
    return boolean ? v.AsBool() : v.AsInt();
  }
  double number(std::size_t i, const std::string&) const {
    const auto v = args[static_cast<int>(i)];
    require(v.IsFloat(), "F3kdb: expected number");
    return v.AsFloat();
  }
  std::string string(std::size_t i, const std::string&) const {
    const auto v = args[static_cast<int>(i)];
    require(v.IsString(), "F3kdb: expected string");
    return v.AsString();
  }
};
PClip input_clip(const AVSValue& args) {
  require(args[0].IsClip(), "F3kdb: clip is required");
  return args[0].AsClip();
}
class F3kdbFilter final : public GenericVideoFilter {
public:
  explicit F3kdbFilter(const AVSValue& args) : GenericVideoFilter(input_clip(args)), input_(vi) {
    require(vi.HasVideo() && vi.width > 0 && vi.height > 0 && vi.num_frames > 0 && vi.IsPlanar() &&
                (vi.IsY() || vi.IsYUV()) && (vi.NumComponents() == 1 || vi.NumComponents() == 3),
            "F3kdb: planar Gray or YUV without alpha required");
    f::Format fmt;
    fmt.color = vi.IsY() ? f::Color::gray : f::Color::yuv;
    fmt.sample = vi.BitsPerComponent() == 32 ? f::Sample::floating : f::Sample::integer;
    fmt.depth = vi.BitsPerComponent();
    fmt.width = vi.width;
    fmt.height = vi.height;
    fmt.frames = vi.num_frames;
    if (!vi.IsY()) {
      fmt.sub_x = vi.GetPlaneWidthSubsampling(PLANAR_U);
      fmt.sub_y = vi.GetPlaneHeightSubsampling(PLANAR_U);
    }
    filter_ = std::make_unique<neo_smo::f3kdb::Filter>(fmt, plugin::f3kdb_options(Reader{args}));
    const auto converted = ds::avisynth::make_video_format(vi);
    require(converted.has_value(), "F3kdb: unsupported input format");
    auto output = converted.value();
    const auto& target = filter_->parameters().output;
    if (target.sample == f::Sample::floating) {
      require(target.depth == 32, "F3kdb: AviSynth does not support F16 output");
      output.sample_format = ds::SampleFormat::Float32;
    } else {
      switch (target.depth) {
        case 8:
          output.sample_format = ds::SampleFormat::UInt8;
          break;
        case 10:
          output.sample_format = ds::SampleFormat::UInt10;
          break;
        case 12:
          output.sample_format = ds::SampleFormat::UInt12;
          break;
        case 14:
          output.sample_format = ds::SampleFormat::UInt14;
          break;
        case 16:
          output.sample_format = ds::SampleFormat::UInt16;
          break;
        default:
          require(false, "F3kdb: output depth unavailable in AviSynth");
      }
    }
    vi.pixel_type = ds::avisynth::pixel_type(output);
    require(vi.pixel_type != VideoInfo::CS_UNKNOWN, "F3kdb: output format unavailable in AviSynth");
    require(static_cast<std::size_t>(vi.width) <=
                static_cast<std::size_t>(std::numeric_limits<int>::max()) / vi.ComponentSize(),
            "F3kdb: output row exceeds AviSynth limits");
  }
  int __stdcall SetCacheHints(int hints, int) override { return hints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0; }
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override {
    try {
      const auto src = child->GetFrame(n, env);
      require(!!src, "F3kdb: input frame unavailable");
      std::vector<f::PlaneView> views;
      for (int p = 0; p < input_.NumComponents(); ++p) {
        const int id = plane(p), stride = src->GetPitch(id), h = src->GetHeight(id);
        const int w = input_.width >> (p ? input_.GetPlaneWidthSubsampling(id) : 0);
        require(h == (input_.height >> (p ? input_.GetPlaneHeightSubsampling(id) : 0)) &&
                    src->GetRowSize(id) / input_.ComponentSize() == w && stride > 0,
                "F3kdb: invalid frame layout");
        require(static_cast<std::size_t>(stride) <=
                    std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(h),
                "F3kdb: plane size overflow");
        views.push_back({src->GetReadPtr(id), static_cast<std::size_t>(stride) * h, static_cast<std::size_t>(stride)});
      }
      const auto output = filter_->process(static_cast<std::size_t>(n), views);
      auto dst = env->NewVideoFrameP(vi, &src);
      for (int p = 0; p < vi.NumComponents(); ++p) {
        const auto& data = output.planes[p];
        env->BitBlt(dst->GetWritePtr(plane(p)), dst->GetPitch(plane(p)), data.data.data(),
                    static_cast<int>(data.stride), static_cast<int>(data.stride), static_cast<int>(data.height));
      }
      return dst;
    } catch (const AvisynthError&) {
      throw;
    } catch (const std::exception& e) {
      env->ThrowError("neo_smo_F3kdb: %s", e.what());
    } catch (...) {
      env->ThrowError("neo_smo_F3kdb: frame processing failed");
    }
    return {};
  }

private:
  static int plane(int p) { return p == 0 ? PLANAR_Y : p == 1 ? PLANAR_U : PLANAR_V; }
  VideoInfo input_;
  std::unique_ptr<neo_smo::f3kdb::Filter> filter_;
};
AVSValue __cdecl create(AVSValue args, void*, IScriptEnvironment* env) {
  try {
    env->CheckVersion(11);
    // PClip transfers ownership through the host SDK linkage table.
    // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
    PClip filter = new F3kdbFilter(args);
    return AVSValue(filter);
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& e) {
    env->ThrowError("neo_smo_F3kdb: %s", e.what());
  } catch (...) {
    env->ThrowError("neo_smo_F3kdb: creation failed");
  }
  return {};
}
} // namespace
void add_f3kdb(IScriptEnvironment* env) {
  auto descriptor = plugin::f3kdb_descriptor();
  descriptor.params.front().required = false;
  const auto signature = ds::make_avisynth_signature(descriptor);
  require(signature.has_value(), "F3kdb: signature generation failed");
  env->AddFunction("neo_smo_F3kdb", env->SaveString(signature.value().c_str()), create, nullptr);
}
} // namespace neo_smo::avs
