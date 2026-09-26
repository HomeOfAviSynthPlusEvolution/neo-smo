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
  std::array<int, 3> params{0, 0, 0}; // radius for Median, mode for VerticalCleaner/RemoveGrain
  std::array<bool, 3> process{false, false, false};
};

const char* algorithm_name(Algorithm alg) noexcept;

FilterPlan build_plan(
  Algorithm alg,
  const FormatInfo& fmt,
  const std::vector<int>& param_list,
  const std::vector<int>& planes_list,
  bool planes_specified
);

} // namespace neo_smo
