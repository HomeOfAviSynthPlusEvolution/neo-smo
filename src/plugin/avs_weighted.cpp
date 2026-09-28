#include "plugin/avs_filter.hpp"
#include <algorithm>
#include <limits>

namespace neo_smo::avs {
void Filter::init_weighted(const Params& p, IScriptEnvironment* env) {
  if (alg_ == Algorithm::TTempSmooth) {
    radius_ = p.integer("maxr", 3);
    require(radius_ >= 1 && radius_ <= 7, "maxr must be 1-7");
    auto th = p.integers("thresh"), md = p.integers("mdiff");
    require(th.size() <= 3 && md.size() <= 3, "thresh and mdiff support at most 3 elements");
    std::array<int, 3> mdiff{2, 3, 3};
    for (int i = 0; i < 3; ++i) {
      if (!th.empty())
        thresholds_[i] = th[std::min(static_cast<std::size_t>(i), th.size() - 1)];
      if (!md.empty())
        mdiff[i] = md[std::min(static_cast<std::size_t>(i), md.size() - 1)];
      require(thresholds_[i] >= 1 && thresholds_[i] <= 256, "thresh must be 1-256");
      require(mdiff[i] >= 0 && mdiff[i] <= 255, "mdiff must be 0-255");
    }
    const int strength = p.integer("strength", 2);
    require(strength >= 1 && strength <= 8, "strength must be 1-8");
    const float sc = static_cast<float>(p.number("scthresh", 12));
    require(sc >= -1 && sc <= 100, "scthresh must be -1 to 100");
    scene_ = sc != 0;
    fp_ = p.boolean("fp", true);
    ref_ = p.clip("pfclip");
    check_reference(ref_);
    if (sc > 0) {
      if (ref_)
        ref_ = scene_detect(ref_, static_cast<double>(sc) / 100, env);
      else
        child = scene_detect(child, static_cast<double>(sc) / 100, env);
    }
    auto calc_temp_weights = [](int r, int str, std::vector<float>& weights, float& cw) {
      weights.resize(static_cast<std::size_t>(r + 1));
      for (int i = 0; i <= r; ++i) {
        weights[static_cast<std::size_t>(i)] = (i < str) ? 1.0f : 1.0f / static_cast<float>(i - str + 2);
      }
      float sum = weights[0];
      for (int i = 1; i <= r; ++i) {
        sum += weights[static_cast<std::size_t>(i)] * 2.0f;
      }
      for (int i = 0; i <= r; ++i) {
        weights[static_cast<std::size_t>(i)] /= sum;
      }
      cw = weights[0];
    };

    auto calc_diff_weights = [](int thresh, int md, int r, int str, std::vector<float>& diff_weights, float& cw) {
      diff_weights.assign(static_cast<std::size_t>(r) * 256, 0.0f);
      std::vector<float> tw(static_cast<std::size_t>(r + 1), 0.0f);
      for (int i = 0; i <= r; ++i) {
        tw[static_cast<std::size_t>(i)] = (i < str) ? 1.0f : 1.0f / static_cast<float>(i - str + 2);
      }
      std::vector<float> dw(256, 0.0f);
      const float step = 256.0f / static_cast<float>(thresh - std::min(md, thresh - 1));
      float base = 256.0f;
      for (int diff = 0; diff < thresh; ++diff) {
        if (diff < md) {
          dw[static_cast<std::size_t>(diff)] = 256.0f;
        } else {
          if (base > 0.0f)
            dw[static_cast<std::size_t>(diff)] = base;
          else
            break;
          base -= step;
        }
      }
      float t_sum = tw[0];
      for (int rad = 1; rad <= r; ++rad) {
        t_sum += tw[static_cast<std::size_t>(rad)] * 2.0f;
        for (int diff = 0; diff < 256; ++diff) {
          diff_weights[static_cast<std::size_t>(rad - 1) * 256 + static_cast<std::size_t>(diff)] =
              tw[static_cast<std::size_t>(rad)] * dw[static_cast<std::size_t>(diff)] / 256.0f;
        }
      }
      for (int rad = 0; rad < r; ++rad) {
        for (int diff = 0; diff < 256; ++diff) {
          diff_weights[static_cast<std::size_t>(rad) * 256 + static_cast<std::size_t>(diff)] /= t_sum;
        }
      }
      cw = tw[0] / t_sum;
    };

    for (int i = 0; i < fmt_.num_planes; ++i) {
      if (thresholds_[i] > mdiff[i] + 1) {
        weight_mode_[i] = 0;
        calc_diff_weights(thresholds_[i], mdiff[i], radius_, strength, difference_weights_[i], center_weights_[i]);
      } else {
        weight_mode_[i] = 1;
        calc_temp_weights(radius_, strength, temporal_weights_[i], center_weights_[i]);
      }
    }
    return;
  }
  ref_ = p.clip("ref", child);
  check_reference(ref_);
  if (alg_ == Algorithm::CCD) {
    require(fmt_.color_family != 1, "CCD requires RGB or YUV");
    threshold_ = static_cast<float>(p.number("threshold", 4));
    require(std::isfinite(threshold_), "threshold is too large");
    threshold_ *= fmt_.is_float ? 1.f : static_cast<float>((1u << fmt_.bits_per_sample) - 1);
    threshold_ = threshold_ * threshold_ / (255.f * 255.f * 3.f);
    require(std::isfinite(threshold_), "threshold is too large");
    radius_ = p.integer("temporal_radius", 0);
    require(radius_ >= 0 && radius_ <= 10, "temporal_radius must be 0-10");
    weights_.assign(radius_ * 2 + 1, 1.f);
    for (int r = 0; r < radius_; ++r) {
      const float tr = static_cast<float>(radius_), fr = static_cast<float>(r);
      weights_[radius_ - 1 - r] = std::sqrt((tr + 1.f - fr) / ((tr + 1.f) * 2.f));
      weights_[radius_ + 1 + r] = std::sin((tr + 2.f - fr) / ((tr + 1.f) * 2.f));
    }
    scale_ = static_cast<float>(p.number("scale", static_cast<float>(fmt_.height >> fmt_.subsampling_h) / 240.f));
    require(std::isfinite(scale_) && scale_ >= 1 &&
                static_cast<double>(scale_) * 12 <= std::numeric_limits<int>::max() / 4,
            "scale is out of range");
    auto sets = p.integers("points");
    if (sets.empty())
      sets = {1, 1, 0};
    require(sets.size() == 3 && (sets[0] || sets[1] || sets[2]),
            "points requires 3 elements with at least one enabled");
    const std::vector<Point> groups[]{{{-4, -4}, {4, -4}, {-4, 4}, {4, 4}},
                                      {{-8, -8}, {0, -8}, {8, -8}, {-8, 0}, {8, 0}, {-8, 8}, {0, 8}, {8, 8}},
                                      {{-12, -12},
                                       {-4, -12},
                                       {4, -12},
                                       {12, -12},
                                       {-12, -4},
                                       {12, -4},
                                       {-12, 4},
                                       {12, 4},
                                       {-12, 12},
                                       {-4, 12},
                                       {4, 12},
                                       {12, 12}}};
    for (int i = 0; i < 3; ++i)
      if (sets[i])
        points_.insert(points_.end(), groups[i].begin(), groups[i].end());
    for (auto& pt : points_) {
      pt.x = static_cast<int>(std::round(static_cast<float>(pt.x >> fmt_.subsampling_w) * scale_));
      pt.y = static_cast<int>(std::round(static_cast<float>(pt.y >> fmt_.subsampling_h) * scale_));
      diameter_ = std::max({diameter_, std::abs(pt.x) * 2 + 1, std::abs(pt.y) * 2 + 1});
    }
    std::sort(points_.begin(), points_.end(), [](Point a, Point b) { return a.y != b.y ? a.y < b.y : a.x < b.x; });
    const double size = std::round(static_cast<double>(diameter_) * scale_);
    require(size <= fmt_.width && size <= fmt_.height, "scaled filter diameter exceeds the source dimensions");
    if (fmt_.color_family == 3)
      ref_luma_ = chroma_luma(ref_, fmt_, env);
    return;
  }
  require(fmt_.color_family == 3 && !fmt_.is_float, "Cnr4 requires 8-16 bit integer YUV");
  std::string mode = p.has("mode") ? p.get("mode").AsString() : "oxx";
  require(mode.size() == 3, "mode requires 3 characters");
  for (char c : mode)
    require(c == 'o' || c == 'x', "mode only accepts o and x");
  radius_ = p.integer("radius", 2);
  tmode_ = p.integer("tmode", 0);
  wmode_ = p.integer("wmode", 0);
  require(radius_ >= 1 && radius_ <= 10, "radius must be 1-10");
  require(tmode_ >= 0 && tmode_ <= 4, "tmode must be 0-4");
  require(wmode_ >= 0 && wmode_ <= 3, "wmode must be 0-3");
  scene_ = p.boolean("scenechange", true);
  auto sense = std::array<int, 3>{35, 47, 47}, strength = std::array<int, 3>{192, 255, 255};
  for (const auto& entry : {std::pair<const char*, std::array<int, 3>*>{"sense", &sense}, {"str", &strength}}) {
    auto values = p.integers(entry.first);
    if (!values.empty()) {
      require(values.size() == 3, std::string(entry.first) + " requires 3 elements");
      for (int i = 0; i < 3; ++i) {
        require(values[i] >= -1 && values[i] <= 255, std::string(entry.first) + " must be -1 to 255");
        if (values[i] != -1)
          (*entry.second)[i] = values[i];
      }
    }
  }
  auto power = p.floats("pow");
  if (power.empty())
    power = {1, 1, 1};
  require(power.size() == 3, "pow requires 3 elements");
  const float pi = 3.14159265358979323846f;
  for (int i = 0; i < 3; ++i) {
    require(power[i] >= 0, "pow must be nonnegative");
    auto& table = tables_[i];
    if (sense[i] == 0) {
      table[0] = static_cast<std::uint8_t>(strength[i]);
      continue;
    }
    const float inv = power[i] > 0 ? 1.f / power[i] : std::numeric_limits<float>::infinity();
    const float s = static_cast<float>(sense[i]);
    for (int l = 0; l <= strength[i]; ++l) {
      const float f = static_cast<float>(l);
      const float base =
          mode[i] == 'o' ? (1.f + std::cos(f * f * pi / (s * s))) / 2.f : (1.f + std::cos(f * pi / s)) / 2.f;
      const float shaped = std::isinf(inv) ? (base >= 1.f ? 1.f : 0.f) : std::pow(base, inv);
      table[l] = static_cast<std::uint8_t>(std::clamp(static_cast<float>(strength[i]) * shaped, 0.f, 255.f));
    }
  }
  luma_ = chroma_luma(child, fmt_, env);
  ref_luma_ = ref_ == child ? luma_ : chroma_luma(ref_, fmt_, env);
}

void Filter::process_weighted(int n, const PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env) {
  const auto dtype = get_data_type(fmt_.bytes_per_sample, fmt_.is_float);
  const int last = vi.num_frames - 1;
  if (alg_ == Algorithm::TTempSmooth) {
    std::vector<PVideoFrame> frames, refs;
    for (int i = -radius_; i <= radius_; ++i) {
      const int idx = static_cast<int>(std::clamp<std::int64_t>(static_cast<std::int64_t>(n) + i, 0, last));
      frames.push_back(frame(child, idx, env));
      refs.push_back(ref_ ? frame(ref_, idx, env) : frames.back());
    }
    int before = radius_, after = radius_;
    if (scene_) {
      for (int i = radius_; i > 0; --i)
        if (scene_property(refs[i], "_SceneChangePrev", env) == 1) {
          before = radius_ - i;
          break;
        }
      for (int i = radius_; i < 2 * radius_; ++i)
        if (scene_property(refs[i], "_SceneChangeNext", env) == 1) {
          after = i - radius_;
          break;
        }
    }
    for (int p = 0; p < fmt_.num_planes; ++p) {
      if (!plan_.process[p])
        continue;
      const int id = plane_id(p), h = src->GetHeight(id), w = src->GetRowSize(id) / fmt_.bytes_per_sample;
      const auto ss = src->GetPitch(id), rs = refs[radius_]->GetPitch(id);
      FrameStorage storage;
      std::vector<const std::uint8_t*> prev, prevref, next, nextref;
      for (int i = 1; i <= radius_; ++i) {
        prev.push_back(storage.plane(frames[radius_ - i], id, ss, h));
        prevref.push_back(storage.plane(refs[radius_ - i], id, rs, h));
        next.push_back(storage.plane(frames[radius_ + i], id, ss, h));
        nextref.push_back(storage.plane(refs[radius_ + i], id, rs, h));
      }
      process_ttempsmooth_plane(
          dtype, w, h, ss, rs, dst->GetPitch(id), src->GetReadPtr(id), refs[radius_]->GetReadPtr(id), prev.data(),
          prevref.data(), next.data(), nextref.data(), dst->GetWritePtr(id), radius_, before, after,
          scale_to_format(fmt_, static_cast<float>(thresholds_[p])), fp_, fmt_.bits_per_sample, weight_mode_[p],
          center_weights_[p], temporal_weights_[p].data(), difference_weights_[p].data());
    }
    return;
  }
  const bool rgb = fmt_.color_family == 2;
  const int output_id = plane_id(rgb ? 0 : 1);
  const int h = dst->GetHeight(output_id), w = dst->GetRowSize(output_id) / fmt_.bytes_per_sample;
  const std::size_t stride = dst->GetPitch(output_id);
  FrameStorage storage;
  // Keep host frames alive until the kernels finish reading their planes.
  std::vector<PVideoFrame> owners;
  auto get_planes = [&](PClip clip, PClip luma, int idx, bool ccd_source) {
    auto f = frame(clip, idx, env);
    owners.push_back(f);
    std::array<const std::uint8_t*, 3> result{};
    for (int p = 0; p < 3; ++p) {
      if (p == 0 && !rgb && ccd_source) {
        // YUV CCD does not read source luma; retain its identity for shared-guide dispatch.
        result[p] = f->GetReadPtr(plane_id(p));
        continue;
      }
      if (p == 0 && luma) {
        auto y = frame(luma, idx, env);
        owners.push_back(y);
        result[p] = storage.plane(y, PLANAR_Y, stride, h);
      } else
        result[p] = storage.plane(f, plane_id(p), stride, h);
    }
    return result;
  };
  if (alg_ == Algorithm::CCD) {
    std::vector<const std::uint8_t*> sources, references;
    for (int i = n - radius_; i <= n + radius_; ++i) {
      const auto a = get_planes(child, {}, i, true), b = get_planes(ref_, ref_luma_, i, false);
      sources.insert(sources.end(), a.begin(), a.end());
      references.insert(references.end(), b.begin(), b.end());
    }
    process_ccd_planes(dtype, rgb, w, h, stride, sources.data(), references.data(),
                       rgb ? dst->GetWritePtr(plane_id(0)) : nullptr, dst->GetWritePtr(plane_id(1)),
                       dst->GetWritePtr(plane_id(2)), threshold_, radius_, weights_.data(), points_.data(),
                       static_cast<int>(points_.size()), diameter_, scale_, fmt_.bits_per_sample);
    return;
  }
  std::vector<int> indices;
  std::vector<PVideoFrame> scene_frames;
  for (int i = -radius_; i <= radius_; ++i) {
    if (i == 0 && tmode_ == 0)
      continue;
    indices.push_back(static_cast<int>(std::clamp<std::int64_t>(static_cast<std::int64_t>(n) + i, 0, last)));
    if (scene_)
      scene_frames.push_back(frame(child, indices.back(), env));
  }
  const int count = static_cast<int>(indices.size());
  int from = 0, to = count - 1;
  if (scene_) {
    scene_property(scene_frames[0], "_SceneChangePrev", env, true);
    scene_property(scene_frames[0], "_SceneChangeNext", env, true);
    for (int i = tmode_ != 0 ? count / 2 : count / 2 - 1; i > 0; --i)
      if (scene_property(scene_frames[i], "_SceneChangePrev", env) == 1) {
        from = i;
        break;
      }
    for (int i = count / 2; i < count - 1; ++i)
      if (scene_property(scene_frames[i], "_SceneChangeNext", env) == 1) {
        to = i;
        break;
      }
    if (tmode_ == 0) {
      if (scene_property(src, "_SceneChangePrev", env) == 1)
        from = count / 2;
      if (scene_property(src, "_SceneChangeNext", env) == 1)
        to = count / 2 - 1;
    }
  }
  const auto current = get_planes(child, luma_, n, false), reference = get_planes(ref_, ref_luma_, n, false);
  std::vector<std::array<const std::uint8_t*, 3>> sources, references;
  for (int i = 0; i < count; ++i) {
    sources.push_back(i < from || i > to ? current : get_planes(child, luma_, indices[i], false));
    references.push_back(i < from || i > to ? reference : get_planes(ref_, ref_luma_, indices[i], false));
  }
  process_cnr4_frame(dtype, w, h, stride, fmt_.bits_per_sample, radius_, tmode_, wmode_, current.data(),
                     reference.data(), sources.data(), references.data(), sources.size(), dst->GetWritePtr(plane_id(1)),
                     dst->GetWritePtr(plane_id(2)), tables_[0].data(), tables_[1].data(), tables_[2].data());
}
} // namespace neo_smo::avs
