#include "plugin/avs_filter.hpp"
#include <algorithm>
#include <cstring>
#include <limits>

namespace neo_smo::avs {
Params::Params(Algorithm alg, const AVSValue& args) {
  const auto d = plugin::descriptor(alg);
  for (std::size_t i = 0; i < d.params.size(); ++i) {
    const auto& p = d.params[i];
    const auto& value = args[static_cast<int>(i)];
    require(!p.required || value.Defined(), p.name + " is required");
    if (value.IsArray())
      require(value.ArraySize() > 0, p.name + " cannot be empty");
    values_.emplace(p.name, value);
  }
}
AVSValue Params::get(const char* name) const {
  const auto it = values_.find(name);
  return it == values_.end() ? AVSValue() : it->second;
}
int Params::integer(const char* name, int fallback) const {
  const auto v = get(name);
  if (!v.Defined())
    return fallback;
  require(v.IsInt(), std::string(name) + " must be an integer");
  return v.AsInt();
}
double Params::number(const char* name, double fallback) const {
  const auto v = get(name);
  if (!v.Defined())
    return fallback;
  require(v.IsFloat(), std::string(name) + " must be numeric");
  const double result = v.AsFloat();
  require(std::isfinite(result), std::string(name) + " must be finite");
  return result;
}
bool Params::boolean(const char* name, bool fallback) const {
  const auto v = get(name);
  if (!v.Defined())
    return fallback;
  require(v.IsBool(), std::string(name) + " must be a bool");
  return v.AsBool();
}
std::vector<int> Params::integers(const char* name) const {
  const auto v = get(name);
  std::vector<int> result;
  if (!v.Defined())
    return result;
  const int count = v.IsArray() ? v.ArraySize() : 1;
  for (int i = 0; i < count; ++i) {
    const auto x = v.IsArray() ? v[i] : v;
    require(x.IsInt(), std::string(name) + " must contain integers");
    result.push_back(x.AsInt());
  }
  return result;
}
std::vector<double> Params::numbers(const char* name) const {
  const auto v = get(name);
  std::vector<double> result;
  if (!v.Defined())
    return result;
  const int count = v.IsArray() ? v.ArraySize() : 1;
  for (int i = 0; i < count; ++i) {
    const auto x = v.IsArray() ? v[i] : v;
    require(x.IsFloat(), std::string(name) + " must contain numbers");
    const double value = x.AsFloat();
    require(std::isfinite(value), std::string(name) + " must be finite");
    result.push_back(value);
  }
  return result;
}
std::vector<float> Params::floats(const char* name) const {
  std::vector<float> result;
  for (double x : numbers(name)) {
    require(std::abs(x) <= std::numeric_limits<float>::max(), std::string(name) + " is too large");
    result.push_back(static_cast<float>(x));
  }
  return result;
}
PClip Params::clip(const char* name, PClip fallback) const {
  const auto v = get(name);
  if (!v.Defined())
    return fallback;
  require(v.IsClip(), std::string(name) + " must be a clip");
  return v.AsClip();
}
const std::uint8_t* FrameStorage::plane(const PVideoFrame& f, int id, std::size_t stride, int height) {
  const int row = f->GetRowSize(id);
  require(height > 0 && height == f->GetHeight(id) && row > 0 && stride >= static_cast<std::size_t>(row),
          "invalid plane geometry");
  if (stride == static_cast<std::size_t>(f->GetPitch(id)))
    return f->GetReadPtr(id);
  require(stride <= std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(height), "plane is too large");
  auto& buf = rows.emplace_back(stride * height);
  for (int y = 0; y < height; ++y)
    std::memcpy(buf.data() + y * stride, f->GetReadPtr(id) + static_cast<std::ptrdiff_t>(y) * f->GetPitch(id), row);
  return buf.data();
}
void Filter::check_reference(PClip clip) const {
  if (!clip)
    return;
  const auto& r = clip->GetVideoInfo();
  require(r.pixel_type == vi.pixel_type && r.width == vi.width && r.height == vi.height,
          "reference format and dimensions must match the source");
}
PVideoFrame Filter::frame(PClip clip, int n, IScriptEnvironment* env) const {
  require(n >= 0 && n < clip->GetVideoInfo().num_frames, "reference clip does not contain the requested frame");
  return clip->GetFrame(n, env);
}
Filter::Filter(Algorithm alg, const Params& p, IScriptEnvironment* env)
    : GenericVideoFilter(p.clip("clip")), alg_(alg) {
  require(vi.HasVideo() && vi.width > 0 && vi.height > 0 && vi.num_frames > 0 && vi.IsPlanar() &&
              vi.NumComponents() <= 3,
          "only planar Gray, RGB and YUV without alpha are supported");
  const auto f = ds::avisynth::make_video_format(vi);
  require(f.has_value(), "unsupported pixel format");
  format_ = f.value();
  fmt_ = {vi.IsY()     ? 1
          : vi.IsRGB() ? 2
                       : 3,
          vi.BitsPerComponent() == 32,
          vi.BitsPerComponent(),
          vi.ComponentSize(),
          vi.NumComponents(),
          vi.IsYUV() && !vi.IsY() ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0,
          vi.IsYUV() && !vi.IsY() ? vi.GetPlaneHeightSubsampling(PLANAR_U) : 0,
          vi.width,
          vi.height};
  auto params = alg == Algorithm::Cnr4 ? std::vector<int>{} : p.integers(p.has("radius") ? "radius" : "mode");
  const char* th = alg == Algorithm::DegrainMedian                                   ? "limit"
                   : alg == Algorithm::FluxSmoothT || alg == Algorithm::FluxSmoothST ? "temporal_threshold"
                                                                                     : "threshold";
  plan_ = build_plan(alg, fmt_, params, p.floats(th), p.boolean("scalep"), p.integers("planes"), p.has("planes"));
  if (alg == Algorithm::DCTFilter || alg == Algorithm::TTempSmooth) {
    plan_.process.fill(false);
    if (!p.has("planes"))
      for (int i = 0; i < fmt_.num_planes; ++i)
        plan_.process[i] = true;
    else
      for (int i : p.integers("planes")) {
        require(i >= 0 && i < fmt_.num_planes && !plan_.process[i], "planes must contain distinct valid plane indices");
        plan_.process[i] = true;
      }
  }
  if (alg == Algorithm::DCTFilter)
    dct_ = std::make_unique<DctFilter>(p.numbers("factors"));
  ref_ = p.clip("repairclip");
  prev_ = p.clip("previous");
  next_ = p.clip("next");
  check_reference(ref_);
  check_reference(prev_);
  check_reference(next_);
  plan_.interlaced = p.boolean("interlaced");
  plan_.norow = p.boolean("norow");
  if (alg == Algorithm::FluxSmoothST) {
    plan_.secondary_thresholds = flux_thresholds(fmt_, p.floats("spatial_threshold"), p.boolean("scalep"));
    for (int i = 0; i < fmt_.num_planes; ++i)
      plan_.process[i] = plan_.process[i] && (plan_.thresholds[i] >= 0 || plan_.secondary_thresholds[i] >= 0);
  }
  if (alg == Algorithm::TemporalMedian)
    scene_ = p.boolean("scenechange");
  if (alg == Algorithm::TemporalSoften) {
    const int sc = p.integer("scenechange", 0);
    require(sc >= -1 && sc <= 254, "scenechange must be -1 to 254");
    scene_ = sc != 0;
    if (sc > 0)
      child = scene_detect(child, sc / 255.0, env);
  }
  if (alg == Algorithm::TTempSmooth || alg == Algorithm::CCD || alg == Algorithm::Cnr4)
    init_weighted(p, env);
}

PVideoFrame __stdcall Filter::GetFrame(int n, IScriptEnvironment* env) {
  try {
    auto src = frame(child, n, env);
    const int last = vi.num_frames - 1;
    const bool triple = alg_ == Algorithm::Clense || alg_ == Algorithm::TemporalRepair ||
                        alg_ == Algorithm::DegrainMedian || alg_ == Algorithm::FluxSmoothT ||
                        alg_ == Algorithm::FluxSmoothST;
    if ((triple && (n == 0 || n == last)) || (alg_ == Algorithm::ForwardClense && n >= last - 1) ||
        (alg_ == Algorithm::BackwardClense && n < 2) ||
        (alg_ == Algorithm::TemporalMedian && (n < plan_.params[0] || n > last - plan_.params[0])) ||
        (alg_ == Algorithm::CCD && (n < radius_ || n > last - radius_)))
      return src;
    auto dst = env->NewVideoFrameP(vi, &src);
    for (int p = 0; p < fmt_.num_planes; ++p) {
      const bool processed = alg_ == Algorithm::CCD    ? (fmt_.color_family == 2 || p > 0)
                             : alg_ == Algorithm::Cnr4 ? p > 0
                                                       : plan_.process[p];
      if (processed)
        continue;
      const int id = plane_id(p);
      env->BitBlt(dst->GetWritePtr(id), dst->GetPitch(id), src->GetReadPtr(id), src->GetPitch(id), src->GetRowSize(id),
                  src->GetHeight(id));
    }
    if (alg_ == Algorithm::TTempSmooth || alg_ == Algorithm::CCD || alg_ == Algorithm::Cnr4) {
      process_weighted(n, src, dst, env);
      return dst;
    }
    PVideoFrame a, b, c;
    if (alg_ == Algorithm::Repair)
      b = frame(ref_, n, env);
    if (triple) {
      a = frame(alg_ == Algorithm::TemporalRepair ? ref_ : prev_ ? prev_ : child, n - 1, env);
      c = frame(alg_ == Algorithm::TemporalRepair ? ref_ : next_ ? next_ : child, n + 1, env);
      if (alg_ == Algorithm::TemporalRepair)
        b = frame(ref_, n, env);
    }
    if (alg_ == Algorithm::ForwardClense || alg_ == Algorithm::BackwardClense) {
      const int direction = alg_ == Algorithm::ForwardClense ? 1 : -1;
      a = frame(child, n + direction, env);
      c = frame(child, n + 2 * direction, env);
    }
    std::vector<PVideoFrame> window;
    int center = 0, from = 0, to = 0;
    if (alg_ == Algorithm::TemporalMedian || alg_ == Algorithm::TemporalSoften) {
      const int r = plan_.params[0], start = std::max(0, n - r), end = n + std::min(r, last - n);
      for (int i = start; i <= end; ++i)
        window.push_back(frame(child, i, env));
      center = n - start;
      to = end - start;
      if (scene_) {
        const bool required = alg_ == Algorithm::TemporalMedian;
        scene_property(window[0], "_SceneChangePrev", env, required);
        scene_property(window[0], "_SceneChangeNext", env, required);
        for (int i = center; i > 0; --i)
          if (scene_property(window[i], "_SceneChangePrev", env)) {
            from = i;
            break;
          }
        for (int i = center; i < to; ++i)
          if (scene_property(window[i], "_SceneChangeNext", env)) {
            to = i;
            break;
          }
      }
    }
    const auto dtype = get_data_type(fmt_.bytes_per_sample, fmt_.is_float);
    std::unique_ptr<DctScratch> scratch;
    if (dct_)
      scratch = std::make_unique<DctScratch>(*dct_, vi.width);
    for (int p = 0; p < fmt_.num_planes; ++p) {
      if (!plan_.process[p])
        continue;
      const int id = plane_id(p), w = src->GetRowSize(id) / fmt_.bytes_per_sample, h = src->GetHeight(id);
      const auto* s = src->GetReadPtr(id);
      auto* d = dst->GetWritePtr(id);
      const auto ss = src->GetPitch(id), ds = dst->GetPitch(id);
      const bool chroma = fmt_.color_family == 3 && p > 0;
      switch (alg_) {
        case Algorithm::Repair:
          plugin::execute_repair_plane(plan_, p, s, b->GetReadPtr(id), d, w, h, ss, b->GetPitch(id), ds);
          break;
        case Algorithm::Clense:
          plugin::execute_clense_plane(plan_, s, a->GetReadPtr(id), c->GetReadPtr(id), d, w, h, ss, a->GetPitch(id),
                                       c->GetPitch(id), ds);
          break;
        case Algorithm::ForwardClense:
        case Algorithm::BackwardClense:
          plugin::execute_clense_forward_backward_plane(plan_, s, a->GetReadPtr(id), c->GetReadPtr(id), d, w, h, ss,
                                                        a->GetPitch(id), c->GetPitch(id), ds);
          break;
        case Algorithm::TemporalRepair:
          process_temporal_repair_plane(dtype, plan_.params[p], chroma, fmt_.bits_per_sample, s, a->GetReadPtr(id),
                                        b->GetReadPtr(id), c->GetReadPtr(id), d, w, h, ss, a->GetPitch(id),
                                        b->GetPitch(id), c->GetPitch(id), ds);
          break;
        case Algorithm::DegrainMedian:
          process_degrain_median_plane(dtype, plan_.params[p], plan_.thresholds[p], plan_.interlaced, plan_.norow,
                                       chroma, fmt_.bits_per_sample, a->GetReadPtr(id), s, c->GetReadPtr(id), d, w, h,
                                       a->GetPitch(id), ss, c->GetPitch(id), ds);
          break;
        case Algorithm::FluxSmoothT:
          process_fluxsmooth_t_plane(dtype, plan_.thresholds[p], a->GetReadPtr(id), s, c->GetReadPtr(id), d, w, h,
                                     a->GetPitch(id), ss, c->GetPitch(id), ds);
          break;
        case Algorithm::FluxSmoothST:
          process_fluxsmooth_st_plane(dtype, plan_.thresholds[p], plan_.secondary_thresholds[p], a->GetReadPtr(id), s,
                                      c->GetReadPtr(id), d, w, h, a->GetPitch(id), ss, c->GetPitch(id), ds);
          break;
        case Algorithm::TemporalMedian:
        case Algorithm::TemporalSoften: {
          FrameStorage storage;
          std::vector<const std::uint8_t*> ptrs;
          if (alg_ == Algorithm::TemporalMedian) {
            for (int i = from; i <= to; ++i)
              ptrs.push_back(storage.plane(window[i], id, ss, h));
            process_temporal_median_plane(dtype, static_cast<int>(ptrs.size()), ptrs.data(), d, w, h, ss, ds);
          } else {
            ptrs.push_back(s);
            for (int i = center - 1; i >= from; --i)
              ptrs.push_back(storage.plane(window[i], id, ss, h));
            for (int i = center + 1; i <= to; ++i)
              ptrs.push_back(storage.plane(window[i], id, ss, h));
            process_temporal_soften_plane(dtype, static_cast<int>(ptrs.size()), plan_.thresholds[p], ptrs.data(), d, w,
                                          h, ss, ds);
          }
          break;
        }
        case Algorithm::DCTFilter:
          scratch->process(*dct_, dtype, fmt_.bits_per_sample, s, d, w, h, ss, ds);
          break;
        default:
          plugin::execute_plane(plan_, p, s, d, w, h, ss, ds);
          break;
      }
    }
    return dst;
  } catch (const AvisynthError&) {
    throw;
  } catch (const std::exception& e) {
    env->ThrowError("neo_smo_%s: %s", algorithm_name(alg_), e.what());
  } catch (...) {
    env->ThrowError("neo_smo_%s: unexpected frame error", algorithm_name(alg_));
  }
  return {};
}
} // namespace neo_smo::avs
