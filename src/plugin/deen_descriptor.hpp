#pragma once
#include <dualsynth/video_bridge.hpp>

namespace neo_smo::plugin {
inline ds::FilterDescriptor deen_descriptor(bool mini) {
  using P = ds::ParamType;
  ds::FilterDescriptor d;
  d.name = mini ? "MiniDeen" : "Deen";
  const auto add = [&](const char* name, P type, bool array = false, bool required = false) {
    d.params.push_back({name,
                        type,
                        {},
                        required,
                        array,
                        true,
                        true,
                        array ? ds::AvisynthArrayBinding::Native : ds::AvisynthArrayBinding::Legacy});
  };
  add("clip", P::Clip, false, true);
  if (!mini)
    add("mode", P::String);
  add("radius", P::Integer, true);
  add("threshold", P::Float, true);
  if (!mini) {
    add("temporal_threshold", P::Float, true);
    add("minimum", P::Float, true);
    add("scenechange", P::Integer);
  }
  add("scalep", P::Boolean);
  add("planes", P::Integer, true);
  return d;
}
} // namespace neo_smo::plugin
