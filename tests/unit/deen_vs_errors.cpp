// Exercise Deen's real VS callback without requiring a host installation.
#include "plugin/deen.cpp"
#include <cstdio>
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
    throw std::runtime_error("Deen VS error/request contract failed");
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
    for (const char* mode : {"c2d", "c3d", "w2d", "w3d", "a2d", "a3d"}) {
      neo_smo::DeenParameters options;
      options.mode = mode;
      neo_smo::plugin::Instance instance(neo_smo::deen_config(options, {1, false, 8, 1, 1, 0, 0}, false));
      instance.vi.format = format;
      instance.vi.width = instance.vi.height = 3;
      instance.vi.numFrames = 3;
      instance.node = reinterpret_cast<VSNode*>(source);
      const int count = mode[1] == '3' ? 3 : 1;
      for (int n : {0, 1, 2}) {
        requests.clear();
        neo_smo::plugin::get_frame(n, arInitial, &instance, nullptr, nullptr, nullptr, &api);
        check(requests == (n == 1 && count == 3 ? std::vector<int>{1, 0, 2} : std::vector<int>{n}));
      }
      for (missing = 0; missing < count; ++missing) {
        fetched = released = errors = 0;
        check(!neo_smo::plugin::get_frame(1, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api));
        check(released == missing && fetched == missing + 1 && errors == 0);
      }
      missing = -1;
      allocation_failure = true;
      fetched = released = errors = 0;
      check(!neo_smo::plugin::get_frame(1, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api));
      check(released == count && errors == 1);
      allocation_failure = false;
      invalid_format = true;
      fetched = released = errors = 0;
      check(!neo_smo::plugin::get_frame(1, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api));
      check(released == 1 && errors == 1);
      invalid_format = false;
      format.sampleType = stFloat;
      format.bitsPerSample = 32;
      format.bytesPerSample = 4;
      instance.vi.format = format;
      const float nan = std::numeric_limits<float>::quiet_NaN();
      std::memcpy(source, &nan, 4);
      fetched = released = errors = 0;
      check(!neo_smo::plugin::get_frame(1, arAllFramesReady, &instance, nullptr, nullptr, nullptr, &api));
      check(released == (count + 1) && errors == 1);
      std::memset(source, 0, sizeof(source));
      format.sampleType = stInteger;
      format.bitsPerSample = 8;
      format.bytesPerSample = 1;
    }
    std::puts("Deen VS requests and failure ownership passed");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
