// Exercise the real VS callback with a fake VSAPI; no host or production fault hook.
#include <cstdlib>
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
unsigned char src[49]{}, dst[49]{};
int freed = 0, errors = 0, scenario = 0;
const VSFrame* VS_CC get_frame(int, VSNode*, VSFrameContext*) noexcept {
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
  return 7;
}
const uint8_t* VS_CC read_ptr(const VSFrame*, int) noexcept {
  return src;
}
uint8_t* VS_CC write_ptr(VSFrame*, int) noexcept {
  fail_allocation = scenario == 1 || scenario == 3 ? scenario : 0;
  return dst;
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
  std::puts("VS success, allocation failure, null destination and unknown exception passed");
}
