#pragma once
#include "algorithms/deen_config.hpp"
#include <vapoursynth/VapourSynth4.h>

namespace neo_smo::plugin {
inline DeenParameters deen_parameters(const VSMap* in, const VSAPI* api, bool mini) {
  DeenParameters p;
  const auto integers = [&](const char* key, auto& output) {
    const int count = api->mapNumElements(in, key);
    if (count < 0)
      return;
    output.emplace();
    for (int i = 0; i < count; ++i) {
      int err = 0;
      const auto value = api->mapGetInt(in, key, i, &err);
      require(!err, std::string("invalid ") + key);
      output->push_back(value);
    }
  };
  const auto numbers = [&](const char* key, auto& output) {
    const int count = api->mapNumElements(in, key);
    if (count < 0)
      return;
    output.emplace();
    for (int i = 0; i < count; ++i) {
      int err = 0;
      const auto value = api->mapGetFloat(in, key, i, &err);
      require(!err, std::string("invalid ") + key);
      output->push_back(value);
    }
  };
  const auto integer = [&](const char* key, std::int64_t fallback) {
    if (api->mapNumElements(in, key) < 0)
      return fallback;
    require(api->mapNumElements(in, key) == 1, std::string("expected one ") + key);
    int err = 0;
    const auto value = api->mapGetInt(in, key, 0, &err);
    require(!err, std::string("invalid ") + key);
    return value;
  };
  integers("radius", p.radius);
  integers("planes", p.planes);
  numbers("threshold", p.threshold);
  const auto scalep = integer("scalep", 0);
  require(scalep == 0 || scalep == 1, "scalep must be boolean");
  p.scalep = scalep != 0;
  if (!mini) {
    numbers("temporal_threshold", p.temporal_threshold);
    numbers("minimum", p.minimum);
    p.scenechange = integer("scenechange", 0);
    if (api->mapNumElements(in, "mode") >= 0) {
      require(api->mapNumElements(in, "mode") == 1, "expected one mode");
      int err = 0;
      const char* value = api->mapGetData(in, "mode", 0, &err);
      require(!err && value, "invalid mode");
      const int size = api->mapGetDataSize(in, "mode", 0, &err);
      require(!err && size == 3, "invalid mode");
      p.mode.assign(value, static_cast<std::size_t>(size));
    }
  }
  return p;
}
} // namespace neo_smo::plugin
