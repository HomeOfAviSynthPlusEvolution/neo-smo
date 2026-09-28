#pragma once

#include "base/types.hpp"
#include <neo_dct.h>
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace neo_smo {

using DctPlanPtr = std::unique_ptr<neo_dct_plan, decltype(&neo_dct_plan_destroy)>;
using DctWorkspacePtr = std::unique_ptr<neo_dct_workspace, decltype(&neo_dct_workspace_destroy)>;

// Immutable coefficients and transform plan can be shared by all frame requests.
class DctFilter {
public:
  explicit DctFilter(const std::vector<double>& factors);
  const neo_dct_plan* plan() const noexcept { return plan_.get(); }
  const std::array<float, 64>& weights() const noexcept { return weights_; }

private:
  DctPlanPtr plan_{nullptr, neo_dct_plan_destroy};
  std::array<float, 64> weights_{};
};

// One per frame request, reused across its planes. Storage is eight padded rows.
class DctScratch {
public:
  DctScratch(const DctFilter& filter, std::size_t max_width);
  // Source and destination planes must be disjoint, as in newVideoFrame2.
  void process(const DctFilter& filter, DataType type, int bits, const std::uint8_t* src, std::uint8_t* dst,
               std::size_t width, std::size_t height, std::size_t src_stride, std::size_t dst_stride);

private:
  DctWorkspacePtr workspace_{nullptr, neo_dct_workspace_destroy};
  std::vector<float> strip_;
};

} // namespace neo_smo
