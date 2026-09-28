#pragma once

#include "plugin/descriptors.hpp"
#include "algorithms/dct_filter.hpp"
#include "common/copy.hpp"
#include <dualsynth/avisynth/video_bridge.hpp>
#include <cmath>
#include <map>

namespace neo_smo::avs {
class Params {
public:
  Params(Algorithm alg, const AVSValue& args);
  AVSValue get(const char* name) const;
  bool has(const char* name) const { return get(name).Defined(); }
  int integer(const char* name, int fallback) const;
  double number(const char* name, double fallback) const;
  bool boolean(const char* name, bool fallback = false) const;
  std::vector<int> integers(const char* name) const;
  std::vector<double> numbers(const char* name) const;
  std::vector<float> floats(const char* name) const;
  PClip clip(const char* name, PClip fallback = {}) const;

private:
  std::map<std::string, AVSValue> values_;
};

// Normalize only mismatched pitches, including cropped and reference frames.
struct FrameStorage {
  std::vector<std::vector<std::uint8_t>> rows;
  const std::uint8_t* plane(const PVideoFrame& frame, int id, std::size_t stride, int height);
};

class Filter final : public GenericVideoFilter {
public:
  Filter(Algorithm algorithm, const Params& params, IScriptEnvironment* env);
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override;
  int __stdcall SetCacheHints(int hints, int) override { return hints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0; }

private:
  Algorithm alg_;
  FormatInfo fmt_{};
  ds::VideoFormat format_{};
  FilterPlan plan_{};
  PClip ref_, prev_, next_, luma_, ref_luma_;
  std::unique_ptr<DctFilter> dct_;
  int radius_ = 0;
  bool scene_ = false;
  bool fp_ = true;
  int tmode_ = 0, wmode_ = 0;
  float threshold_ = 0, scale_ = 1;
  int diameter_ = 0;
  std::vector<Point> points_;
  std::vector<float> weights_;
  std::array<int, 3> thresholds_{4, 5, 5}, weight_mode_{};
  std::array<float, 3> center_weights_{};
  std::array<std::vector<float>, 3> temporal_weights_, difference_weights_;
  std::array<std::array<std::uint8_t, 256>, 3> tables_{};
  int plane_id(int p) const { return ds::avisynth::plane_id(format_, p); }
  PVideoFrame frame(PClip clip, int n, IScriptEnvironment* env) const;
  void check_reference(PClip clip) const;
  void init_weighted(const Params& p, IScriptEnvironment* env);
  void process_weighted(int n, const PVideoFrame& src, PVideoFrame& dst, IScriptEnvironment* env);
};

PClip scene_detect(PClip clip, double threshold, IScriptEnvironment* env);
std::int64_t scene_property(const PVideoFrame& frame, const char* key, IScriptEnvironment* env, bool required = false);
PClip chroma_luma(PClip clip, const FormatInfo& format, IScriptEnvironment* env);
} // namespace neo_smo::avs
