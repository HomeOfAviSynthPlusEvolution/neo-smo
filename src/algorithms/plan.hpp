#pragma once

#include "base/checked.hpp"
#include "base/types.hpp"
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace neo_smo {

enum class Algorithm {
  Median,
  VerticalCleaner,
  RemoveGrain,
  Repair,
  Clense,
  ForwardClense,
  BackwardClense,
  InterQuartileMean,
  SmartMedian,
  TemporalMedian,
  TemporalSoften,
  TemporalRepair,
  DegrainMedian,
  FluxSmoothT,
  FluxSmoothST,
  TTempSmooth,
  CCD,
  Cnr4,
  DCTFilter,
};

struct FormatInfo {
  int color_family; // 1 = Gray, 2 = RGB, 3 = YUV
  bool is_float;
  int bits_per_sample;
  int bytes_per_sample;
  int num_planes;
  int subsampling_w;
  int subsampling_h;
  int width;
  int height;
};

struct FilterPlan {
  Algorithm algorithm;
  FormatInfo format;
  std::array<int, 3> params{0, 0, 0}; // radius for Median/IQM/SmartMedian/TemporalMedian/TemporalSoften, mode for others
  std::array<float, 3> thresholds{0.0f, 0.0f, 0.0f}; // primary threshold / limit
  std::array<float, 3> secondary_thresholds{0.0f, 0.0f, 0.0f}; // spatial threshold for FluxSmoothST
  std::array<bool, 3> process{false, false, false};
  bool scenechange = false; // for TemporalMedian / TemporalSoften
  bool interlaced = false;  // for DegrainMedian
  bool norow = false;       // for DegrainMedian
};

std::array<float, 3> flux_thresholds(const FormatInfo& fmt, const std::vector<float>& values, bool scalep);

const char* algorithm_name(Algorithm alg) noexcept;

inline float scale_to_format(const FormatInfo& fmt, float value) noexcept {
  if (fmt.is_float) {
    return value / 255.0f;
  }
  if (fmt.bits_per_sample > 8) {
    return value * static_cast<float>(1 << (fmt.bits_per_sample - 8));
  }
  return value;
}

FilterPlan build_plan(
  Algorithm alg,
  const FormatInfo& fmt,
  const std::vector<int>& param_list,
  const std::vector<int>& planes_list,
  bool planes_specified
);

FilterPlan build_plan(
  Algorithm alg,
  const FormatInfo& fmt,
  const std::vector<int>& param_list,
  const std::vector<float>& threshold_list,
  bool scalep,
  const std::vector<int>& planes_list,
  bool planes_specified
);

} // namespace neo_smo
