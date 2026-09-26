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
    case Algorithm::InterQuartileMean:
      return "InterQuartileMean";
    case Algorithm::SmartMedian:
      return "SmartMedian";
    case Algorithm::TemporalMedian:
      return "TemporalMedian";
    case Algorithm::TemporalSoften:
      return "TemporalSoften";
    case Algorithm::TemporalRepair:
      return "TemporalRepair";
    case Algorithm::DegrainMedian:
      return "DegrainMedian";
    case Algorithm::FluxSmoothT:
      return "FluxSmoothT";
    case Algorithm::FluxSmoothST:
      return "FluxSmoothST";
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

std::array<float, 3> flux_thresholds(const FormatInfo& fmt, const std::vector<float>& values, bool scalep) {
  std::array<float, 3> out{};
  for (std::size_t i = 0; i < out.size(); ++i) {
    if (i < values.size()) {
      const float value = values[i];
      require(!scalep || value < 0 || value <= 255, "FluxSmooth: scaled thresholds must be at most 255.");
      out[i] = scalep && value >= 0 ? scale_to_format(fmt, value) : value;
    } else {
      out[i] = i == 0 ? scale_to_format(fmt, 7.0f) : out[i - 1];
    }
  }
  return out;
}

FilterPlan build_plan(
  Algorithm alg,
  const FormatInfo& fmt,
  const std::vector<int>& param_list,
  const std::vector<float>& threshold_list,
  bool scalep,
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
  } else if (alg == Algorithm::InterQuartileMean) {
    const int count = static_cast<int>(param_list.size());
    require(count <= fmt.num_planes, "InterQuartileMean: Element count of radius must be less than or equal to the number of input planes.");
    if (count > 0) {
      for (int i = 0; i < 3; ++i) {
        if (i < count) {
          const int r = param_list[static_cast<std::size_t>(i)];
          require(r >= 0 && r <= 3, "InterQuartileMean: Invalid radius specified, only radius 0-3 supported.");
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
  } else if (alg == Algorithm::SmartMedian) {
    const int count = static_cast<int>(param_list.size());
    require(count <= fmt.num_planes, "SmartMedian: Element count of radius must be less than or equal to the number of input planes.");
    if (count > 0) {
      for (int i = 0; i < 3; ++i) {
        if (i < count) {
          const int r = param_list[static_cast<std::size_t>(i)];
          require(r >= 0 && r <= 3, "SmartMedian: Invalid radius specified, only radius 0-3 supported.");
          plan.params[static_cast<std::size_t>(i)] = r;
        } else {
          plan.params[static_cast<std::size_t>(i)] = plan.params[static_cast<std::size_t>(i - 1)];
        }
      }
    } else {
      plan.params = {1, 1, 1};
    }

    const int count_th = static_cast<int>(threshold_list.size());
    require(count_th <= fmt.num_planes, "SmartMedian: Element count of threshold must be less than or equal to the number of input planes.");
    const float format_max = scalep ? 255.0f : (fmt.is_float ? 1.0f : static_cast<float>((1 << fmt.bits_per_sample) - 1));
    if (count_th > 0) {
      for (int i = 0; i < 3; ++i) {
        if (i < count_th) {
          const float th = threshold_list[static_cast<std::size_t>(i)];
          require(th >= 0.0f && th <= format_max,
                  "SmartMedian: Invalid threshold, must be in the range of 0 - " +
                  (fmt.is_float && !scalep ? std::to_string(format_max) : std::to_string(static_cast<int>(format_max))) +
                  " with scalep = " + (scalep ? "true" : "false") + " for this bit depth");
          plan.thresholds[static_cast<std::size_t>(i)] = scalep ? scale_to_format(fmt, th) : th;
        } else {
          plan.thresholds[static_cast<std::size_t>(i)] = plan.thresholds[static_cast<std::size_t>(i - 1)];
        }
      }
    } else {
      const float fifty = scale_to_format(fmt, 50.0f);
      const float one_twenty_eight = scale_to_format(fmt, 128.0f);
      plan.thresholds = {
          plan.params[0] == 1 ? fifty : one_twenty_eight,
          plan.params[1] == 1 ? fifty : one_twenty_eight,
          plan.params[2] == 1 ? fifty : one_twenty_eight,
      };
    }
    const auto planes = normalize_planes(fmt.num_planes, planes_list, planes_specified, prefix);
    for (int i = 0; i < 3; ++i) {
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && planes[static_cast<std::size_t>(i)] && (plan.params[static_cast<std::size_t>(i)] > 0);
    }
  } else if (alg == Algorithm::TemporalMedian) {
    int r = 1;
    if (!param_list.empty()) {
      r = param_list[0];
    }
    require(r >= 1 && r <= 10, "TemporalMedian: Radius must be between 1 and 10 (inclusive)");
    plan.params = {r, r, r};
    const auto planes = normalize_planes(fmt.num_planes, planes_list, planes_specified, prefix);
    for (int i = 0; i < 3; ++i) {
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && planes[static_cast<std::size_t>(i)];
    }
  } else if (alg == Algorithm::TemporalSoften) {
    int r = 4;
    if (!param_list.empty()) {
      r = param_list[0];
    }
    require(r >= 1 && r <= 10, "TemporalSoften: Radius must be between 1 and 10 (inclusive)");
    plan.params = {r, r, r};

    const int count_th = static_cast<int>(threshold_list.size());
    if (count_th > 0) {
      for (int i = 0; i < 3; ++i) {
        if (i < count_th) {
          const float th = threshold_list[static_cast<std::size_t>(i)];
          require(!scalep || (th >= 0.0f && th <= 255.0f), "TemporalSoften: scaled threshold must be in 0-255.");
          const float value = scalep ? scale_to_format(fmt, th) : th;
          const bool chroma = fmt.is_float && fmt.color_family == 3 && i > 0;
          const float lo = chroma ? -0.5f : 0.0f;
          const float hi = fmt.is_float ? (chroma ? 0.5f : 1.0f) : static_cast<float>((1 << fmt.bits_per_sample) - 1);
          require(value >= lo && value <= hi, "TemporalSoften: Invalid threshold specified.");
          plan.thresholds[static_cast<std::size_t>(i)] = value;
        } else {
          plan.thresholds[static_cast<std::size_t>(i)] = plan.thresholds[static_cast<std::size_t>(i - 1)];
        }
      }
    } else {
      const float def_th = scale_to_format(fmt, 4.0f);
      plan.thresholds = {def_th, def_th, def_th};
    }

    require(fmt.color_family == 3 || plan.thresholds[0] != 0, "TemporalSoften: threshold 0 cannot be zero for RGB or Gray.");
    require(plan.thresholds[0] != 0 || plan.thresholds[1] != 0 || plan.thresholds[2] != 0, "TemporalSoften: All thresholds cannot be 0.");
    const auto planes = normalize_planes(fmt.num_planes, planes_list, planes_specified, prefix);
    for (int i = 0; i < 3; ++i) {
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && planes[static_cast<std::size_t>(i)] && plan.thresholds[static_cast<std::size_t>(i)] > 0;
    }
  } else if (alg == Algorithm::TemporalRepair) {
    const int count = static_cast<int>(param_list.size());
    require(count <= fmt.num_planes, "TemporalRepair: Number of modes must be equal or fewer than the number of input planes.");
    if (count > 0) {
      for (int i = 0; i < 3; ++i) {
        if (i < count) {
          const int m = param_list[static_cast<std::size_t>(i)];
          require(m >= 0 && m <= 4, "TemporalRepair: Invalid mode specified, only modes 0-4 supported.");
          plan.params[static_cast<std::size_t>(i)] = m;
        } else {
          plan.params[static_cast<std::size_t>(i)] = plan.params[static_cast<std::size_t>(i - 1)];
        }
      }
    } else {
      plan.params = {0, 0, 0};
    }
    const auto planes = normalize_planes(fmt.num_planes, planes_list, planes_specified, prefix);
    for (int i = 0; i < 3; ++i) {
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && planes[static_cast<std::size_t>(i)];
    }
  } else if (alg == Algorithm::DegrainMedian) {
    const int count = static_cast<int>(param_list.size());
    require(count <= fmt.num_planes, "DegrainMedian: Number of modes must be equal or fewer than the number of input planes.");
    if (count > 0) {
      for (int i = 0; i < 3; ++i) {
        if (i < count) {
          const int m = param_list[static_cast<std::size_t>(i)];
          require(m >= 0 && m <= 5, "DegrainMedian: Invalid mode specified, only modes 0-5 supported.");
          plan.params[static_cast<std::size_t>(i)] = m;
        } else {
          plan.params[static_cast<std::size_t>(i)] = plan.params[static_cast<std::size_t>(i - 1)];
        }
      }
    } else {
      plan.params = {1, 1, 1};
    }

    const int count_lim = static_cast<int>(threshold_list.size());
    require(count_lim <= fmt.num_planes, "DegrainMedian: limit has more elements than there are planes.");
    if (count_lim > 0) {
      for (int i = 0; i < 3; ++i) {
        if (i < count_lim) {
          const float lim = threshold_list[static_cast<std::size_t>(i)];
          const float formatMaximum = (fmt.color_family == 3 && i > 0) ? (fmt.is_float ? 0.5f : static_cast<float>((1 << fmt.bits_per_sample) - 1))
                                                                        : (fmt.is_float ? 1.0f : static_cast<float>((1 << fmt.bits_per_sample) - 1));
          require(!scalep || (lim >= 0.0f && lim <= 255.0f), "DegrainMedian: scaled limit must be in 0-255.");
          const float scaled_lim = scalep ? (formatMaximum * lim / 255.0f) : lim;
          const float formatMinimum = fmt.is_float && fmt.color_family == 3 && i > 0 ? -0.5f : 0.0f;
          require(scaled_lim >= formatMinimum && scaled_lim <= formatMaximum, "DegrainMedian: Invalid limit specified.");
          plan.thresholds[static_cast<std::size_t>(i)] = scaled_lim;
        } else {
          plan.thresholds[static_cast<std::size_t>(i)] = plan.thresholds[static_cast<std::size_t>(i - 1)];
        }
      }
    } else {
      plan.thresholds = {4.0f, 4.0f, 4.0f}; // Upstream defaults are native units, even with scalep.
    }

    require(plan.thresholds[0] != 0 || plan.thresholds[1] != 0 || plan.thresholds[2] != 0, "DegrainMedian: All limits cannot be 0.");

    const auto planes = normalize_planes(fmt.num_planes, planes_list, planes_specified, prefix);
    for (int i = 0; i < 3; ++i) {
      plan.process[static_cast<std::size_t>(i)] = (i < fmt.num_planes) && planes[static_cast<std::size_t>(i)] && (plan.thresholds[static_cast<std::size_t>(i)] > 0);
    }
  } else if (alg == Algorithm::FluxSmoothT || alg == Algorithm::FluxSmoothST) {
    plan.thresholds = flux_thresholds(fmt, threshold_list, scalep);
    const auto planes = normalize_planes(fmt.num_planes, planes_list, planes_specified, prefix);
    for (int i = 0; i < 3; ++i) {
      plan.process[i] = i < fmt.num_planes && planes[i];
      if (alg == Algorithm::FluxSmoothT) plan.process[i] = plan.process[i] && plan.thresholds[i] >= 0;
    }
  }

  return plan;
}

FilterPlan build_plan(
  Algorithm alg,
  const FormatInfo& fmt,
  const std::vector<int>& param_list,
  const std::vector<int>& planes_list,
  bool planes_specified
) {
  return build_plan(alg, fmt, param_list, {}, false, planes_list, planes_specified);
}

} // namespace neo_smo
