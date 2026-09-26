#pragma once

#include <vapoursynth/VapourSynth4.h>
#include <cstddef>
#include <utility>
#include <stdexcept>
#include <string>

namespace neo_smo {

struct TemporalWindowResult {
  std::size_t from_idx;
  std::size_t to_idx;
};

inline TemporalWindowResult resolve_temporal_window(
    const VSFrame* const* frames,
    std::size_t radius,
    std::size_t diameter,
    bool scenechange,
    const VSAPI* vsapi,
    const char* filter_name,
    bool require_properties) {
  std::size_t from_idx = 0;
  std::size_t to_idx = diameter - 1;

  if (!scenechange) {
    return {from_idx, to_idx};
  }

  // Check scene change properties on the first frame (frame 0)
  const VSMap* p0 = vsapi->getFramePropertiesRO(frames[0]);
  int err_prev = 0;
  int err_next = 0;
  vsapi->mapGetInt(p0, "_SceneChangePrev", 0, &err_prev);
  vsapi->mapGetInt(p0, "_SceneChangeNext", 0, &err_next);

  if (require_properties && (err_prev || err_next)) {
    throw std::runtime_error(std::string(filter_name) +
                             ": Scene change handling enabled, but input frame is missing scene change properties.");
  }

  // Walk backwards from radius (current frame)
  for (std::size_t i = radius; i > 0; --i) {
    const VSMap* p = vsapi->getFramePropertiesRO(frames[i]);
    int err = 0;
    const auto sc_prev = vsapi->mapGetInt(p, "_SceneChangePrev", 0, &err);
    if (!err && sc_prev != 0) {
      from_idx = i;
      break;
    }
  }

  // Walk forwards from radius (current frame)
  for (std::size_t i = radius; i < diameter - 1; ++i) {
    const VSMap* p = vsapi->getFramePropertiesRO(frames[i]);
    int err = 0;
    const auto sc_next = vsapi->mapGetInt(p, "_SceneChangeNext", 0, &err);
    if (!err && sc_next != 0) {
      to_idx = i;
      break;
    }
  }

  return {from_idx, to_idx};
}

} // namespace neo_smo
