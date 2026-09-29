#pragma once
#include "algorithms/deen.hpp"
#include "algorithms/plan.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <optional>
#include <type_traits>

namespace neo_smo {
struct DeenParameters {
  std::string mode = "c3d";
  std::optional<std::vector<std::int64_t>> radius, planes;
  std::optional<std::vector<double>> threshold, temporal_threshold, minimum;
  bool scalep = false;
  std::int64_t scenechange = 0;
};
struct DeenConfig {
  std::string mode;
  bool temporal = false;
  int scenechange = 0;
  std::array<int, 3> radius{};
  std::array<double, 3> threshold{}, temporal_threshold{}, minimum{};
  std::array<bool, 3> process{};
};
inline bool deen_finite(double value) {
  std::uint64_t bits;
  std::memcpy(&bits, &value, sizeof(bits));
  return (bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
}
inline DeenConfig deen_config(const DeenParameters& p, const FormatInfo& f, bool mini) {
  DeenConfig c;
  c.mode = p.mode;
  require(mini || p.mode == "c2d" || p.mode == "c3d" || p.mode == "a2d" || p.mode == "a3d" || p.mode == "w2d" ||
              p.mode == "w3d",
          "Deen: invalid mode");
  c.temporal = !mini && p.mode[1] == '3';
  require(p.scenechange >= -1 && p.scenechange <= 254, "scenechange must be in [-1,254]");
  c.scenechange = static_cast<int>(p.scenechange);
  require(!c.temporal || c.scenechange <= 0 || f.color_family != 2, "automatic scene detection requires Gray or YUV");
  const auto expand = [&](const auto& input, auto& output, double fallback, double lo, double hi, const char* name) {
    if (input)
      require(!input->empty() && input->size() <= static_cast<std::size_t>(f.num_planes),
              std::string(name) + " requires 1..num_planes values");
    for (int i = 0; i < f.num_planes; ++i) {
      const double value =
          input ? static_cast<double>((*input)[std::min<std::size_t>(i, input->size() - 1)]) : fallback;
      require(deen_finite(value) && value >= lo && value <= hi, std::string("invalid ") + name);
      using T = typename std::decay_t<decltype(output)>::value_type;
      output[i] = static_cast<T>(value);
    }
  };
  expand(p.radius, c.radius, 1, 0, c.temporal ? 4 : 7, "radius");
  expand(p.minimum, c.minimum, 0.5, 0, 1, "minimum");
  const double limit = p.scalep ? 255 : f.is_float ? 1 : (1u << f.bits_per_sample) - 1;
  expand(p.threshold, c.threshold, 0, 0, limit, "threshold");
  expand(p.temporal_threshold, c.temporal_threshold, 0, 0, limit, "temporal_threshold");
  // Use the same power-of-two / float conversion as all other filters.
  const double scale = f.is_float ? 1.0 / 255 : scale_to_format(f, 1.0f);
  for (int i = 0; i < f.num_planes; ++i) {
    const bool uv = f.color_family == 3 && i > 0;
    c.threshold[i] = p.threshold ? c.threshold[i] * (p.scalep ? scale : 1) : (mini ? 10 : uv ? 9 : 7) * scale;
    c.temporal_threshold[i] =
        p.temporal_threshold ? c.temporal_threshold[i] * (p.scalep ? scale : 1) : (uv ? 6 : 4) * scale;
    c.process[i] = !p.planes;
  }
  if (p.planes) {
    require(!p.planes->empty(), "planes cannot be empty");
    for (auto i : *p.planes) {
      require(i >= 0 && i < f.num_planes && !c.process[static_cast<std::size_t>(i)],
              "planes must contain distinct valid indices");
      c.process[static_cast<std::size_t>(i)] = true;
    }
  }
  for (int i = 0; i < f.num_planes; ++i)
    c.process[i] = c.process[i] && c.radius[i] != 0;
  return c;
}
inline Deen deen_plane_filter(const DeenConfig& c, int plane) {
  DeenOptions o;
  o.mode = c.mode;
  o.radius = std::max(1, c.radius[plane]);
  o.minimum = c.minimum[plane];
  o.scenechange = false;
  return Deen(std::move(o));
}
} // namespace neo_smo
