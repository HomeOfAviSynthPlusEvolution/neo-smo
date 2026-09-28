// Exercise MiniDeen's real VS callback without requiring a host installation.
#include "plugin/mini_deen.cpp"
#include <cstdio>
#include <vector>
#include <cstring>
#include <limits>

namespace {
alignas(float) std::uint8_t source[36]{}, destination[36]{};
VSVideoFormat format{};
int fetched = 0, released = 0, errors = 0, missing = -1;
bool allocation_failure = false, invalid_format = false;
std::vector<int> requests;
const VSFrame* VS_CC fetch(int, VSNode*, VSFrameContext*) noexcept {
  return fetched++ == missing ? nullptr : reinterpret_cast<const VSFrame*>(source);
}
void VS_CC release(const VSFrame*) noexcept {
  ++released;
}
void VS_CC error(const char*, VSFrameContext*) noexcept {
  ++errors;
}
int VS_CC size(const VSFrame*, int) noexcept {
  return 3;
}
ptrdiff_t VS_CC stride(const VSFrame*, int) noexcept {
  return 3 * format.bytesPerSample;
}
const std::uint8_t* VS_CC read(const VSFrame*, int) noexcept {
  return source;
}
std::uint8_t* VS_CC write(VSFrame*, int) noexcept {
  return destination;
}
const VSVideoFormat* VS_CC frame_format(const VSFrame*) noexcept {
  static VSVideoFormat changed;
  changed = format;
  if (invalid_format)
    changed.bitsPerSample = 16;
  return &changed;
}
VSFrame* VS_CC allocate(const VSVideoFormat*, int, int, const VSFrame**, const int*, const VSFrame*, VSCore*) noexcept {
  return allocation_failure ? nullptr : reinterpret_cast<VSFrame*>(destination);
}
void VS_CC request(int n, VSNode*, VSFrameContext*) noexcept {
  requests.push_back(n);
}
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("MiniDeen VS error/request contract failed");
}
} // namespace
int main() {
  try {
    VSAPI api{};
    api.getFrameFilter = fetch;
    api.freeFrame = release;
    api.setFilterError = error;
    api.getFrameWidth = api.getFrameHeight = size;
    api.getStride = stride;
    api.getReadPtr = read;
    api.getWritePtr = write;
    api.getVideoFrameFormat = frame_format;
    api.newVideoFrame2 = allocate;
    api.requestFrameFilter = request;
    format.colorFamily = cfGray;
    format.sampleType = stInteger;
    format.bytesPerSample = 1;
    format.bitsPerSample = 8;
    format.numPlanes = 1;
    neo_smo::plugin::Instance instance;
    instance.vi.format = format;
    instance.vi.width = instance.vi.height = 3;
    instance.vi.numFrames = 3;
    instance.node = reinterpret_cast<VSNode*>(source);
    for (int n : {2, 0, 1}) {
      requests.clear();
      neo_smo::plugin::get_frame(n, arInitial, &instance, nullptr, nullptr, nullptr, &api);
      check(requests == std::vector<int>{n});
    }
    missing = 0;
    check(!neo_smo::plugin::get_frame(1, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api));
    check(released == 0 && errors == 0);
    missing = -1;
    allocation_failure = true;
    fetched = released = errors = 0;
    check(!neo_smo::plugin::get_frame(1, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api));
    check(released == 1 && errors == 1);
    allocation_failure = false;
    invalid_format = true;
    fetched = released = errors = 0;
    check(!neo_smo::plugin::get_frame(1, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api));
    check(released == 1 && errors == 1);
    invalid_format = false;
    fetched = released = errors = 0;
    const auto* result = neo_smo::plugin::get_frame(1, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api);
    check(result && released == 1 && errors == 0);
    api.freeFrame(result);
    check(released == 2);
    std::puts("MiniDeen VS failure ownership checks passed");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
