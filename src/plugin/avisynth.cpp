#include "plugin/avs_filter.hpp"
#include "plugin/avs_deen.hpp"

const AVS_Linkage* AVS_linkage = nullptr;
namespace neo_smo::avs {
template <Algorithm Alg>
AVSValue __cdecl create(AVSValue args, void*, IScriptEnvironment* env) {
  try {
    env->CheckVersion(11);
    // PClip transfers ownership through the SDK linkage table, opaque to the analyzer.
    // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
    PClip filter = new Filter(Alg, Params(Alg, args), env);
    return AVSValue(filter);
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& e) {
    env->ThrowError("neo_smo_%s: %s", algorithm_name(Alg), e.what());
  } catch (...) {
    env->ThrowError("neo_smo_%s: creation failed", algorithm_name(Alg));
  }
  return {};
}
template <Algorithm Alg>
void add(IScriptEnvironment* env) {
  auto d = plugin::descriptor(Alg);
  for (auto& p : d.params)
    p.required = false;
  const auto result = ds::make_avisynth_signature(d);
  require(result.has_value(), result.has_value() ? "" : result.error().message);
  const auto& signature = result.value();
  const auto name = std::string("neo_smo_") + d.name;
  env->AddFunction(env->SaveString(name.c_str()), env->SaveString(signature.c_str()), create<Alg>, nullptr);
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
#define ADD(name) avs::add<Algorithm::name>(env)
    ADD(Median);
    ADD(VerticalCleaner);
    ADD(RemoveGrain);
    ADD(Repair);
    ADD(Clense);
    ADD(ForwardClense);
    ADD(BackwardClense);
    ADD(InterQuartileMean);
    ADD(SmartMedian);
    ADD(TemporalMedian);
    ADD(TemporalSoften);
    ADD(TemporalRepair);
    ADD(DegrainMedian);
    ADD(FluxSmoothT);
    ADD(FluxSmoothST);
    ADD(TTempSmooth);
    ADD(CCD);
    ADD(Cnr4);
    ADD(DCTFilter);
#undef ADD
    avs::add_deen(env);
    return "neo-smo AviSynth+ filters";
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& e) {
    env->ThrowError("neo-smo: %s", e.what());
  } catch (...) {
    env->ThrowError("neo-smo: registration failed");
  }
  return nullptr;
}
