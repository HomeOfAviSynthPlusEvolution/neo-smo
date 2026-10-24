#pragma once
#include "algorithms/f3kdb.hpp"
#include "base/checked.hpp"
#include <dualsynth/video_bridge.hpp>
#include <limits>

namespace neo_smo::plugin {
inline ds::FilterDescriptor f3kdb_descriptor() {
  using P = ds::ParamType;
  ds::FilterDescriptor d;
  d.name = "F3kdb";
  const auto add = [&](const char* name, P type) {
    d.params.push_back({name, type, {}, false, false, true, true, ds::AvisynthArrayBinding::Legacy});
  };
  add("clip", P::Clip);
  d.params.back().required = true;
  add("range", P::Integer);
  for (auto name : {"y", "cb", "cr", "grainy", "grainc"})
    add(name, P::Float);
  add("sample_mode", P::Integer);
  add("seed", P::Integer);
  add("blur_first", P::Boolean);
  add("dynamic_grain", P::Boolean);
  add("mt", P::Boolean);
  add("dither_algo", P::Integer);
  add("keep_tv_range", P::Boolean);
  add("output_depth", P::Integer);
  add("random_algo_ref", P::Integer);
  add("random_algo_grain", P::Integer);
  add("random_param_ref", P::Float);
  add("random_param_grain", P::Float);
  add("preset", P::String);
  for (auto name : {"y_1", "cb_1", "cr_1", "y_2", "cb_2", "cr_2"})
    add(name, P::Float);
  add("scale", P::Boolean);
  add("angle_boost", P::Float);
  add("max_angle", P::Float);
  add("output_sample_type", P::String);
  return d;
}

// Both hosts use the same descriptor and option assignment, retaining unset fields.
template <class Reader>
neo_smo::f3kdb::Options f3kdb_options(Reader&& r) {
  neo_smo::f3kdb::Options o;
  const auto d = f3kdb_descriptor();
  for (std::size_t i = 1; i < d.params.size(); ++i) {
    const auto& p = d.params[i];
    const auto& key = p.name;
    if (!r.present(i, key))
      continue;
    if (p.type == ds::ParamType::String) {
      const auto v = r.string(i, key);
      if (key == "preset")
        o.presets.push_back(v);
      else {
        require(v == "same" || v == "integer" || v == "float", "F3kdb: invalid output_sample_type");
        o.output_sample_type = v == "same"      ? neo_smo::f3kdb::OutputSample::same
                               : v == "integer" ? neo_smo::f3kdb::OutputSample::integer
                                                : neo_smo::f3kdb::OutputSample::floating;
      }
    } else if (p.type == ds::ParamType::Integer || p.type == ds::ParamType::Boolean) {
      const auto v = r.integer(i, key, p.type == ds::ParamType::Boolean);
      require(v >= std::numeric_limits<int>::min() && v <= std::numeric_limits<int>::max(),
              "F3kdb: integer out of range");
      if (p.type == ds::ParamType::Boolean)
        require(v == 0 || v == 1, "F3kdb: expected boolean");
      const int n = static_cast<int>(v);
      if (key == "range")
        o.range = n;
      else if (key == "sample_mode")
        o.sample_mode = n;
      else if (key == "seed")
        o.seed = n;
      else if (key == "blur_first")
        o.blur_first = n != 0;
      else if (key == "dynamic_grain")
        o.dynamic_grain = n != 0;
      else if (key == "mt")
        o.mt = n != 0;
      else if (key == "dither_algo")
        o.dither_algo = n;
      else if (key == "keep_tv_range")
        o.keep_tv_range = n != 0;
      else if (key == "output_depth")
        o.output_depth = n;
      else if (key == "random_algo_ref")
        o.random_algo_ref = n;
      else if (key == "random_algo_grain")
        o.random_algo_grain = n;
      else if (key == "scale")
        o.scale = n != 0;
    } else {
      const double v = r.number(i, key);
      if (key == "grainy")
        o.grain_y = v;
      else if (key == "grainc")
        o.grain_c = v;
      else if (key == "random_param_ref")
        o.random_param_ref = v;
      else if (key == "random_param_grain")
        o.random_param_grain = v;
      else if (key == "angle_boost")
        o.angle_boost = v;
      else if (key == "max_angle")
        o.max_angle = v;
      else {
        const int plane = key[0] == 'y' ? 0 : key.compare(0, 2, "cb") == 0 ? 1 : 2;
        if (key.back() == '1')
          o.threshold1[plane] = v;
        else if (key.back() == '2')
          o.threshold2[plane] = v;
        else
          o.threshold[plane] = v;
      }
    }
  }
  return o;
}
} // namespace neo_smo::plugin
