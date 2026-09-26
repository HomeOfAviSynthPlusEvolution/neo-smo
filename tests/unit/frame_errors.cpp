// Exercise the real VS callback with a fake VSAPI; no host or production fault hook.
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <new>
#include <cstdio>
static int fail_allocation = 0;
void* operator new(std::size_t n) {
  if (fail_allocation) {
    const int kind = fail_allocation;
    fail_allocation = 0;
    if (kind == 3)
      throw 42;
    throw std::bad_alloc();
  }
  if (void* p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc();
}
void operator delete(void* p) noexcept {
  std::free(p);
}
void operator delete(void* p, std::size_t) noexcept {
  std::free(p);
}
#include "plugin/vapoursynth.cpp"

namespace {
alignas(float) unsigned char src[49 * sizeof(float)]{}, dst[49 * sizeof(float)]{};
int bytes_per_sample = 1;
int freed = 0, errors = 0, scenario = 0;
int fetched = 0, null_at = -1;
bool bad_access = false;
const VSFrame* VS_CC get_frame(int, VSNode*, VSFrameContext*) noexcept {
  if (fetched++ == null_at) return nullptr;
  return reinterpret_cast<const VSFrame*>(src);
}
VSFrame* VS_CC new_frame(const VSVideoFormat*, int, int, const VSFrame**, const int*, const VSFrame*,
                         VSCore*) noexcept {
  return scenario == 2 ? nullptr : reinterpret_cast<VSFrame*>(dst);
}
void VS_CC free_frame(const VSFrame*) noexcept {
  ++freed;
}
void VS_CC set_error(const char*, VSFrameContext*) noexcept {
  ++errors;
}
int VS_CC dimension(const VSFrame*, int) noexcept {
  return 7;
}
ptrdiff_t VS_CC stride(const VSFrame*, int) noexcept {
  return 7 * bytes_per_sample;
}
const uint8_t* VS_CC read_ptr(const VSFrame* frame, int) noexcept {
  if (!frame) bad_access = true;
  return src;
}
uint8_t* VS_CC write_ptr(VSFrame*, int) noexcept {
  fail_allocation = scenario == 1 || scenario == 3 ? scenario : 0;
  return dst;
}

const VSFrame* VS_CC add_ref(const VSFrame* frame) noexcept {
  if (!frame) bad_access = true;
  ++fetched;
  return frame;
}
VSVideoInfo main_vi{}, ref_vi{};
VSNode* main_node = reinterpret_cast<VSNode*>(src);
VSNode* other_node = reinterpret_cast<VSNode*>(dst);
int dependency_case = 0, dependency_count = 0;
VSFilterDependency dependencies[3]{};
VSNode* VS_CC get_node(const VSMap*, const char* key, int, int* error) noexcept {
  *error = 0;
  if (std::strcmp(key, "clip") == 0)
    return main_node;
  return dependency_case == 0 ? main_node : other_node;
}
const VSVideoInfo* VS_CC video_info(VSNode* node) noexcept {
  return node == main_node ? &main_vi : &ref_vi;
}
int VS_CC num_elements(const VSMap*, const char*) noexcept { return -1; }
void VS_CC free_node(VSNode*) noexcept {}
void VS_CC create_filter(VSMap*, const char*, const VSVideoInfo*, VSFilterGetFrame, VSFilterFree free,
                        int, const VSFilterDependency* deps, int count, void* instance, VSCore* core) noexcept;
const VSAPI* creation_api = nullptr;
void VS_CC create_filter(VSMap*, const char*, const VSVideoInfo*, VSFilterGetFrame, VSFilterFree free,
                        int, const VSFilterDependency* deps, int count, void* instance, VSCore* core) noexcept {
  dependency_count = count;
  std::copy_n(deps, count, dependencies);
  free(instance, core, creation_api);
}
} // namespace
int main() {
  VSAPI api{};
  api.getFrameFilter = get_frame;
  api.newVideoFrame2 = new_frame;
  api.freeFrame = free_frame;
  api.setFilterError = set_error;
  api.getFrameWidth = dimension;
  api.getFrameHeight = dimension;
  api.getStride = stride;
  api.getReadPtr = read_ptr;
  api.getWritePtr = write_ptr;
  neo_smo::plugin::VSFilterInstance instance{};
  instance.vi.width = instance.vi.height = 7;
  instance.vi.format.numPlanes = 1;
  instance.plan.algorithm = neo_smo::Algorithm::Median;
  instance.plan.format.bytes_per_sample = 1;
  instance.plan.process[0] = true;
  instance.plan.params[0] = 1;
  for (scenario = 0; scenario <= 3; ++scenario) {
    freed = errors = 0;
    const auto* result =
        neo_smo::plugin::filter_get_frame(0, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api);
    if (scenario == 0) {
      if (!result || freed != 1 || errors != 0)
        return 1;
      api.freeFrame(result);
    } else if (result || errors != 1 || freed != (scenario == 2 ? 1 : 2)) {
      std::fprintf(stderr, "scenario=%d freed=%d errors=%d\n", scenario, freed, errors);
      return 1;
    }
  }

  neo_smo::plugin::VSRepairInstance rep_instance{};
  rep_instance.vi.width = rep_instance.vi.height = 7;
  rep_instance.vi.format.numPlanes = 1;
  rep_instance.plan.algorithm = neo_smo::Algorithm::Repair;
  rep_instance.plan.format.bytes_per_sample = 1;
  rep_instance.plan.process[0] = true;
  rep_instance.plan.params[0] = 1;
  for (scenario = 0; scenario <= 3; ++scenario) {
    freed = errors = 0;
    const auto* result =
        neo_smo::plugin::repair_get_frame(0, arAllFramesReady, &rep_instance, nullptr, nullptr, nullptr, &api);
    if (scenario == 0) {
      if (!result || freed != 2 || errors != 0)
        return 1;
      api.freeFrame(result);
    } else if (result || errors != 1 || freed != (scenario == 2 ? 2 : 3)) {
      std::fprintf(stderr, "repair scenario=%d freed=%d errors=%d\n", scenario, freed, errors);
      return 1;
    }
  }

  neo_smo::plugin::VSClenseInstance clense_instance{};
  clense_instance.vi.width = clense_instance.vi.height = 7;
  clense_instance.vi.numFrames = 10;
  clense_instance.vi.format.numPlanes = 1;
  clense_instance.plan.algorithm = neo_smo::Algorithm::Clense;
  clense_instance.plan.format.bytes_per_sample = 2;
  clense_instance.plan.format.is_float = true;
  bytes_per_sample = 2;
  clense_instance.plan.process[0] = true;
  for (scenario = 0; scenario <= 3; ++scenario) {
    freed = errors = 0;
    void* frame_data = reinterpret_cast<void*>(1);
    const auto* result =
        neo_smo::plugin::clense_get_frame<neo_smo::Algorithm::Clense>(1, arAllFramesReady, &clense_instance, &frame_data, nullptr, nullptr, &api);
    if (scenario == 0) {
      if (!result || freed != 3 || errors != 0)
        return 1;
      api.freeFrame(result);
    } else if (result || errors != 1 || freed != (scenario == 2 ? 3 : 4)) {
      std::fprintf(stderr, "clense scenario=%d freed=%d errors=%d\n", scenario, freed, errors);
      return 1;
    }
  }

  // New temporal callbacks must release every acquired frame on allocation
  // failure, unexpected exceptions, and a missing source frame.
  instance.vi.numFrames = 10;
  instance.plan.format.bytes_per_sample = 2;
  instance.plan.format.is_float = true;
  instance.plan.thresholds[0] = 0.1f;
  rep_instance.vi.numFrames = 10;
  rep_instance.plan.format.bytes_per_sample = 2;
  rep_instance.plan.format.is_float = true;
  const std::array<VSFilterGetFrame, 6> callbacks{
    neo_smo::plugin::temporal_filter_get_frame<neo_smo::Algorithm::TemporalMedian>,
    neo_smo::plugin::temporal_filter_get_frame<neo_smo::Algorithm::TemporalSoften>,
    neo_smo::plugin::trio_filter_get_frame<neo_smo::Algorithm::DegrainMedian>,
    neo_smo::plugin::trio_filter_get_frame<neo_smo::Algorithm::FluxSmoothT>,
    neo_smo::plugin::trio_filter_get_frame<neo_smo::Algorithm::FluxSmoothST>,
    neo_smo::plugin::temporal_repair_get_frame,
  };
  for (std::size_t filter = 0; filter < callbacks.size(); ++filter) {
    const int sources = filter == 5 ? 4 : 3;
    void* data = filter == 5 ? static_cast<void*>(&rep_instance) : static_cast<void*>(&instance);
    for (scenario = 0; scenario <= 3; ++scenario) {
      freed = errors = fetched = 0;
      const auto* result = callbacks[filter](2, arAllFramesReady, data, nullptr, nullptr, nullptr, &api);
      if (scenario == 0) {
        if (!result || freed != sources || errors) return 1;
        api.freeFrame(result);
      } else if (result || errors != 1 || freed != sources + (scenario != 2)) {
        std::fprintf(stderr, "phase4 filter=%zu scenario=%d freed=%d errors=%d\n", filter, scenario, freed, errors);
        return 1;
      }
    }
    scenario = 0;
    for (null_at = 0; null_at < sources; ++null_at) {
      freed = errors = fetched = 0;
      const auto* result = callbacks[filter](2, arAllFramesReady, data, nullptr, nullptr, nullptr, &api);
      if (result || freed != fetched - 1) return 1;
    }
    null_at = -1;
  }

  // Phase5 callbacks: missing inputs and allocation failure must release every
  // frame acquired, including Cnr4's duplicated current-frame references.
  bytes_per_sample = 1;
  api.addFrameRef = add_ref;
  neo_smo::plugin::VSTTempSmoothInstance tt{};
  neo_smo::plugin::VSCCDInstance ccd{};
  neo_smo::plugin::VSCnr4Instance cnr{};
  tt.vi = ccd.vi = cnr.vi = main_vi;
  for (auto* vi : {&tt.vi, &ccd.vi, &cnr.vi}) {
    vi->width = vi->height = 7; vi->numFrames = 10;
    vi->format.colorFamily = cfYUV;
    vi->format.numPlanes = 3; vi->format.bytesPerSample = 1;
    vi->format.bitsPerSample = 8;
  }
  tt.maxr = cnr.radius = 1;
  tt.scenechange = cnr.scenechange = false;
  tt.weight_mode = {1, 1, 1}; tt.center_weights = {0.5f, 0.5f, 0.5f};
  for (auto& weights : tt.temporal_weights) weights = {0.5f, 0.25f};
  ccd.weights = {1.0f}; ccd.points = {{0, 0}};
  const std::array<VSFilterGetFrame, 3> phase5_callbacks{
    neo_smo::plugin::ttempsmooth_get_frame, neo_smo::plugin::ccd_get_frame, neo_smo::plugin::cnr4_get_frame};
  void* phase5_instances[] = {&tt, &ccd, &cnr};
  for (int filter = 0; filter < 3; ++filter) {
    int normal_fetches = 0;
    for (scenario = 0; scenario <= 3; ++scenario) {
      // CCD's vector kernel performs no dynamic allocations after WritePtr.
      if (filter == 1 && (scenario == 1 || scenario == 3)) continue;
      freed = errors = fetched = 0; bad_access = false;
      const auto* result = phase5_callbacks[filter](4, arAllFramesReady, phase5_instances[filter], nullptr, nullptr, nullptr, &api);
      fail_allocation = 0;
      if (scenario == 0) {
        normal_fetches = fetched;
        if (!result || freed != fetched || errors || bad_access) return 1;
        api.freeFrame(result);
      } else if (result || errors != 1 || freed != fetched + (scenario != 2) || bad_access) {
        std::fprintf(stderr, "phase5 filter=%d scenario=%d fetched=%d freed=%d errors=%d\n", filter, scenario, fetched, freed, errors);
        return 1;
      }
    }
    scenario = 0;
    for (null_at = 0; null_at < normal_fetches; ++null_at) {
      freed = errors = fetched = 0; bad_access = false;
      const auto* result = phase5_callbacks[filter](4, arAllFramesReady, phase5_instances[filter], nullptr, nullptr, nullptr, &api);
      // The last two Cnr4 references are addFrameRef, not upstream fetches.
      if (result) { api.freeFrame(result); continue; }
      if (freed != fetched - 1 || bad_access) return 1;
    }
    null_at = -1;
  }

  // Observe the real create callback's cache contract, including aliased inputs.
  api.mapGetNode = get_node;
  api.getVideoInfo = video_info;
  api.mapNumElements = num_elements;
  api.freeNode = free_node;
  api.createVideoFilter = create_filter;
  creation_api = &api;
  main_vi.width = main_vi.height = 7;
  main_vi.numFrames = 10;
  main_vi.format.colorFamily = cfGray;
  main_vi.format.numPlanes = main_vi.format.bytesPerSample = 1;
  main_vi.format.bitsPerSample = 8;
  ref_vi = main_vi;
  for (dependency_case = 0; dependency_case < 2; ++dependency_case) {
    neo_smo::plugin::clense_create<neo_smo::Algorithm::Clense>(nullptr, nullptr, nullptr, nullptr, &api);
    if (dependency_count != dependency_case + 1 ||
        dependencies[dependency_case].requestPattern != rpGeneral)
      return 1;
  }
  neo_smo::plugin::clense_create<neo_smo::Algorithm::ForwardClense>(nullptr, nullptr, nullptr, nullptr, &api);
  if (dependency_count != 1 || dependencies[0].requestPattern != rpGeneral)
    return 1;
  neo_smo::plugin::clense_create<neo_smo::Algorithm::BackwardClense>(nullptr, nullptr, nullptr, nullptr, &api);
  if (dependency_count != 1 || dependencies[0].requestPattern != rpGeneral)
    return 1;
  for (dependency_case = 0; dependency_case < 2; ++dependency_case) {
    neo_smo::plugin::temporal_repair_create(nullptr, nullptr, nullptr, nullptr, &api);
    if (dependency_count != dependency_case + 1 || dependencies[dependency_case].requestPattern != rpGeneral)
      return 1;
  }
  std::puts("VS frame ownership, exceptions and temporal dependency contracts passed");
}
