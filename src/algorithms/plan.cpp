#include "algorithms/plan.hpp"

namespace neo_smo {

const char* algorithm_name(Algorithm alg) noexcept {
  switch (alg) {
    case Algorithm::Median:
      return "Median";
    case Algorithm::VerticalCleaner:
      return "VerticalCleaner";
    case Algorithm::RemoveGrain:
      return "RemoveGrain";
    case Algorithm::Repair:
      return "Repair";
    case Algorithm::Clense:
      return "Clense";
    case Algorithm::ForwardClense:
      return "ForwardClense";
    case Algorithm::BackwardClense:
      return "BackwardClense";
  }
  return "neo_smo";
}

namespace {

std::array<bool, 3> normalize_planes(int num_planes, const std::vector<int>& planes_list, bool planes_specified, const std::string& prefix) {
  if (!planes_specified) {
    return {true, true, true};
  }
  std::array<bool, 3> out{false, false, false};
  for (int p : planes_list) {
    require(p >= 0 && p < num_planes, prefix + ": Plane index out of range.");
    require(!out[static_cast<std::size_t>(p)], prefix + ": Plane specified twice.");
    out[static_cast<std::size_t>(p)] = true;
  }
  return out;
}

} // namespace

FilterPlan build_plan(
  Algorithm alg,
  const FormatInfo& fmt,
  const std::vector<int>& param_list,
  const std::vector<int>& planes_list,
  bool planes_specified
) {
  const std::string prefix = algorithm_name(alg);
  FilterPlan plan{};
  plan.algorithm = alg;
  plan.format = fmt;

  if (alg == Algorithm::Median) {
    const int count = static_cast<int>(param_list.size());
    require(count <= fmt.num_planes, "Median: Element count of radius must be less than or equal to the number of input planes.");
    if (count > 0) {
      for (int i = 0; i < 3; ++i) {
        if (i < count) {
          const int r = param_list[static_cast<std::size_t>(i)];
          require(r >= 0 && r <= 3, "Median: Invalid radius specified, only radius 0-3 supported.");
          plan.params[static_cast<std::size_t>(i)] = r;
        } else {
          plan.params[static_cast<std::size_t>(i)] = plan.params[static_cast<std::size_t>(i - 1)];
        }
      }
    } else {
      plan.params = {1, 1, 1};
    }
    const auto planes = normalize_planes(fmt.num_planes, planes_list, planes_specified, prefix);
    for (int i = 0; i < 3; ++i) {
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && planes[static_cast<std::size_t>(i)] && (plan.params[static_cast<std::size_t>(i)] > 0);
    }
  } else if (alg == Algorithm::VerticalCleaner) {
    const int count = static_cast<int>(param_list.size());
    require(count > 0, "VerticalCleaner: mode is required.");
    require(count <= fmt.num_planes, "VerticalCleaner: Number of modes must be equal or fewer than the number of input planes.");
    for (int i = 0; i < 3; ++i) {
      if (i < count) {
        const int m = param_list[static_cast<std::size_t>(i)];
        require(m >= 0 && m <= 2, "VerticalCleaner: Invalid mode specified, only modes 0-2 supported.");
        plan.params[static_cast<std::size_t>(i)] = m;
      } else {
        plan.params[static_cast<std::size_t>(i)] = plan.params[static_cast<std::size_t>(i - 1)];
      }
      const int plane_height = fmt.height >> (i > 0 ? fmt.subsampling_h : 0);
      if (plan.params[static_cast<std::size_t>(i)] == 1) {
        require(plane_height >= 3, "VerticalCleaner: corresponding plane's height must be greater than or equal to 3 for mode 1");
      } else if (plan.params[static_cast<std::size_t>(i)] == 2) {
        require(plane_height >= 5, "VerticalCleaner: corresponding plane's height must be greater than or equal to 5 for mode 2");
      }
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && (plan.params[static_cast<std::size_t>(i)] > 0);
    }
  } else if (alg == Algorithm::RemoveGrain) {
    const int count = static_cast<int>(param_list.size());
    require(count > 0, "RemoveGrain: mode is required.");
    require(count <= fmt.num_planes, "RemoveGrain: Number of modes must be equal or fewer than the number of input planes.");
    for (int i = 0; i < 3; ++i) {
      if (i < count) {
        const int m = param_list[static_cast<std::size_t>(i)];
        require(m >= 0 && m <= 24, "RemoveGrain: Invalid mode specified, only modes 0-24 supported.");
        plan.params[static_cast<std::size_t>(i)] = m;
      } else {
        plan.params[static_cast<std::size_t>(i)] = plan.params[static_cast<std::size_t>(i - 1)];
      }
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && (plan.params[static_cast<std::size_t>(i)] > 0);
    }
  } else if (alg == Algorithm::Repair) {
    const int count = static_cast<int>(param_list.size());
    require(count > 0, "Repair: mode is required.");
    require(count <= fmt.num_planes, "Repair: Number of modes must be equal or fewer than the number of input planes.");
    for (int i = 0; i < 3; ++i) {
      if (i < count) {
        const int m = param_list[static_cast<std::size_t>(i)];
        require(m >= 0 && m <= 24, "Repair: Invalid mode specified, only modes 0-24 supported.");
        plan.params[static_cast<std::size_t>(i)] = m;
      } else {
        plan.params[static_cast<std::size_t>(i)] = plan.params[static_cast<std::size_t>(i - 1)];
      }
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && (plan.params[static_cast<std::size_t>(i)] > 0);
    }
  } else if (alg == Algorithm::Clense || alg == Algorithm::ForwardClense || alg == Algorithm::BackwardClense) {
    const auto planes = normalize_planes(fmt.num_planes, planes_list, planes_specified, prefix);
    for (int i = 0; i < 3; ++i) {
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && planes[static_cast<std::size_t>(i)];
    }
  }

  return plan;
}

} // namespace neo_smo
