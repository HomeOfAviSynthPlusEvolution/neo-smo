#pragma once

#include "base/types.hpp"
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace neo_smo {
enum class DeenFamily { Constant, Weighted, Adaptive };
struct DeenOptions {
  std::string mode = "c3d";
  int radius = 1;
  double spatial_y = 7, spatial_uv = 9, temporal_y = 4, temporal_uv = 6;
  double minimum = 0.5, scene_threshold = 9;
  bool scenechange = true;
};
struct DeenByteWeights {
  std::array<std::uint16_t, 450> coefficients{};
  std::uint16_t rounding = 128;
};
class Deen {
public:
  explicit Deen(DeenOptions options);
  const DeenOptions& options() const { return options_; }
  DeenFamily family() const { return family_; }
  bool temporal() const { return temporal_; }
  const std::vector<double>& weights() const { return weights_; }
  const DeenByteWeights& byte_weights(bool temporal) const { return byte_weights_[temporal]; }

private:
  DeenOptions options_;
  DeenFamily family_{};
  bool temporal_{};
  std::vector<double> weights_;
  std::array<DeenByteWeights, 2> byte_weights_{};
};
struct DeenPlane {
  const std::uint8_t* data = nullptr;
  std::ptrdiff_t stride = 0;
  int width = 0, height = 0;
  DataType type = DataType::U8;
  int bits = 8;
};
// Frames are current, previous, next. Only current is consumed for spatial evaluation.
void deen_process(const Deen& filter, bool chroma, bool temporal, const std::array<DeenPlane, 3>& frames,
                  std::uint8_t* dst, std::ptrdiff_t stride);
bool deen_scene_cut(const Deen& filter, const std::vector<DeenPlane>& a, const std::vector<DeenPlane>& b);
} // namespace neo_smo
