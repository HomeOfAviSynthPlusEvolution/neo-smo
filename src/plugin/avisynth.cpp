#include "plugin/descriptors.hpp"
#include "common/copy.hpp"
#include <dualsynth/avisynth/video_bridge.hpp>

const AVS_Linkage* AVS_linkage = nullptr;

namespace neo_smo::avs {
namespace av = ds::avisynth;

template <Algorithm Alg>
struct CoreFilter {
  inline static const char* name = algorithm_name(Alg);
  using State = FilterPlan;
  static constexpr int input_count = 1;
  static constexpr ds::OutputOrigin output_origin = ds::OutputOrigin::fresh(0);
  static constexpr ds::HostRequirements host_requirements{true, 0, 0};

  static ds::Result<ds::VideoInitStateResult<State>> init(ds::VideoInitContext& ctx) {
    try {
      const auto& in_info = ctx.input_info(0);
      FormatInfo fmt{};
      if (in_info.format.color_family == ds::ColorFamily::Gray) fmt.color_family = 1;
      else if (in_info.format.color_family == ds::ColorFamily::Rgb) fmt.color_family = 2;
      else if (in_info.format.color_family == ds::ColorFamily::Yuv) fmt.color_family = 3;
      fmt.is_float = (in_info.format.sample_format == ds::SampleFormat::Float32);
      fmt.bits_per_sample = ds::bits_per_sample(in_info.format.sample_format);
      fmt.bytes_per_sample = ds::bytes_per_sample(in_info.format.sample_format);
      fmt.num_planes = in_info.format.plane_count;
      fmt.subsampling_w = in_info.format.subsampling_w;
      fmt.subsampling_h = in_info.format.subsampling_h;
      fmt.width = in_info.width;
      fmt.height = in_info.height;

      const char* param_key = (Alg == Algorithm::Median) ? "radius" : "mode";
      std::vector<int> param_list;
      if (const auto* p = ctx.params.find(param_key)) {
        for (const auto& v : p->values) {
          param_list.push_back(static_cast<int>(std::get<std::int64_t>(v)));
        }
      }

      std::vector<int> planes_list;
      bool planes_specified = false;
      if (const auto* p = ctx.params.find("planes")) {
        planes_specified = true;
        for (const auto& v : p->values) {
          planes_list.push_back(static_cast<int>(std::get<std::int64_t>(v)));
        }
      }

      FilterPlan plan = build_plan(Alg, fmt, param_list, planes_list, planes_specified);
      return ds::Result<ds::VideoInitStateResult<State>>::success({in_info, std::move(plan)});
    } catch (const std::exception& e) {
      return ds::Result<ds::VideoInitStateResult<State>>::failure({ds::ErrorCode::InvalidArgument, e.what()});
    }
  }

  static ds::Result<ds::VideoRequestResult> request(ds::VideoRequestContext& ctx) {
    ctx.request_frame(0, ctx.output_frame);
    return ds::Result<ds::VideoRequestResult>::success({});
  }

  static ds::VideoRequestPattern request_pattern(int, const State&) {
    return ds::VideoRequestPattern::StrictSpatial;
  }

  static ds::Result<ds::VideoProcessResult> process(ds::VideoProcessContext& ctx) {
    const auto& plan = ctx.state<FilterPlan>();
    auto src_res = ctx.frames.get(0, ctx.output_frame);
    if (!src_res.has_value()) {
      return ds::Result<ds::VideoProcessResult>::failure(src_res.error());
    }
    const auto& src = src_res.value().view;
    auto& dst = ctx.output;

    for (int plane = 0; plane < plan.format.num_planes; ++plane) {
      const auto& sp = src.planes[static_cast<std::size_t>(plane)];
      auto& dp = dst.planes[static_cast<std::size_t>(plane)];
      const auto w = static_cast<std::size_t>(dp.width);
      const auto h = static_cast<std::size_t>(dp.height);
      const auto ss = static_cast<std::size_t>(sp.stride_bytes);
      const auto ds_bytes = static_cast<std::size_t>(dp.stride_bytes);

      if (plan.process[static_cast<std::size_t>(plane)]) {
        plugin::execute_plane(plan, plane, sp.data, dp.data, w, h, ss, ds_bytes);
      } else {
        copy_plane(dp.data, sp.data, w * static_cast<std::size_t>(plan.format.bytes_per_sample), h, ds_bytes, ss);
      }
    }
    return ds::Result<ds::VideoProcessResult>::success({});
  }
};

template <Algorithm Alg>
struct Bridge {
  using Core = CoreFilter<Alg>;
  static constexpr bool forward_audio = true;
  static constexpr std::size_t parity_source_index = 0;
  static constexpr av::MtMode avs_mt_mode = av::MtMode::NiceFilter;
  static constexpr const char* avs_format_error = "neo-smo: unsupported planar format";
  static ds::FilterDescriptor descriptor() {
    return plugin::descriptor(Alg);
  }
};

template <Algorithm Alg>
AVSValue __cdecl create(AVSValue args, void*, IScriptEnvironment* env) {
  try {
    env->CheckVersion(11);
    const auto d = Bridge<Alg>::descriptor();
    std::vector<AVSValue> values(d.params.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = args.IsArray() ? (i < static_cast<std::size_t>(args.ArraySize()) ? args[static_cast<int>(i)] : AVSValue())
                                 : (i == 0 ? args : AVSValue());
      if (d.params[i].is_array && values[i].Defined() && !values[i].IsArray()) {
        AVSValue scalar = values[i];
        values[i] = AVSValue(&scalar, 1);
      }
    }
    return av::create_video_filter_bridge<Bridge<Alg>>(AVSValue(values.data(), static_cast<int>(values.size())), env);
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& e) {
    env->ThrowError("neo-smo: %s", e.what());
  } catch (...) {
    env->ThrowError("neo-smo: AviSynth creation failed");
  }
  return {};
}

template <Algorithm Alg>
void add(IScriptEnvironment* env) {
  const auto d = Bridge<Alg>::descriptor();
  const std::string name = std::string("neo_smo_") + d.name;
  auto sig_d = d;
  for (auto& p : sig_d.params) {
    p.required = false;
  }
  const auto sig = ds::make_avisynth_signature(sig_d).value();
  env->AddFunction(env->SaveString(name.c_str()), env->SaveString(sig.c_str()), create<Alg>, nullptr);
}

} // namespace neo_smo::avs

#ifdef _WIN32
#define NEO_SMO_AVS_EXPORT extern "C" __declspec(dllexport)
#else
#define NEO_SMO_AVS_EXPORT extern "C" __attribute__((visibility("default")))
#endif

NEO_SMO_AVS_EXPORT const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* linkage) {
  AVS_linkage = linkage;
  try {
    env->CheckVersion(11);
    using namespace neo_smo;
    avs::add<Algorithm::Median>(env);
    avs::add<Algorithm::VerticalCleaner>(env);
    avs::add<Algorithm::RemoveGrain>(env);
    return "neo-smo AviSynth+ filters";
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& e) {
    env->ThrowError("neo-smo: %s", e.what());
  } catch (...) {
    env->ThrowError("neo-smo: AviSynth registration failed");
  }
  return nullptr;
}
