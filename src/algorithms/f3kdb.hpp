#pragma once

// Production F3kdb configuration and immutable filter state.
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace neo_smo::f3kdb {
enum class Color { gray, yuv };
enum class Sample { integer, floating };
enum class OutputSample { same, integer, floating };
struct Format {
  Color color = Color::gray;
  Sample sample = Sample::integer;
  int depth = 8;
  std::size_t width = 1, height = 1, frames = 1;
  int sub_x = 0, sub_y = 0;
};

// Optional fields distinguish explicit overrides from defaults/presets.
struct Options {
  bool scale = false;
  OutputSample output_sample_type = OutputSample::same;
  int output_depth = -1;
  std::vector<std::string> presets;
  std::array<std::optional<double>, 3> threshold, threshold1, threshold2;
  std::optional<double> grain_y, grain_c;
  int range = 15, sample_mode = 2, dither_algo = 3;
  bool blur_first = true, dynamic_grain = false, keep_tv_range = false, mt = true;
  std::int32_t seed = 0;
  int random_algo_ref = 1, random_algo_grain = 1;
  double random_param_ref = 1, random_param_grain = 1;
  double angle_boost = 1.5, max_angle = 0.15;
};
struct Parameters {
  Format input, output;
  Options options;
  std::array<std::array<double, 3>, 3> thresholds{}; // plane, main/aux1/aux2; processing units
  std::array<double, 2> grain{};
};
Parameters normalize(const Format&, const Options&);

struct PlaneView {
  const void* data = nullptr;
  std::size_t bytes = 0, stride = 0;
};
struct MutablePlaneView {
  void* data = nullptr;
  std::size_t bytes = 0, stride = 0;
};
struct Plane {
  std::size_t width = 0, height = 0, stride = 0;
  std::vector<std::uint8_t> data;
};
struct Frame {
  Format format;
  std::vector<Plane> planes;
};
struct Distance {
  std::uint8_t first = 0, second = 0;
};

// Public pure helpers make mathematical boundary tests independent of RNG tables.
double round_even(double);
double response(double delta, double threshold);
double weight(const std::array<double, 4>&);
unsigned bayer16(std::size_t y, std::size_t x);
float decode_half(std::uint16_t);
std::uint16_t encode_half(double);
bool aligned(const std::array<double, 2>&, const std::array<double, 2>&, double angle);
double deband(double center, const std::array<double, 4>& refs, int mode, const std::array<double, 3>& thresholds,
              bool integer, bool blur_first, bool directions_aligned = false, double angle_boost = 1.5);

class Random {
public:
  explicit Random(std::uint32_t state) : state_(state) {}
  double next(int algorithm, double sigma);
  std::uint32_t state() const { return state_; }

private:
  std::uint32_t state_;
};

// Construction freezes all random state. process() is const and has no shared
// mutable scratch, so requests can run repeatedly, concurrently, or out of order.
class Filter {
public:
  Filter(const Format&, const Options& = {});
  const Parameters& parameters() const { return params_; }
  const std::vector<Distance>& distances() const { return distances_; }
  const std::vector<std::size_t>& offsets() const { return offsets_; }
  const std::vector<float>& grain_y() const { return grain_y_; }
  const std::vector<float>& grain_c() const { return grain_c_; }
  const std::vector<std::int16_t>& grain_y_integer() const { return grain_y_integer_; }
  const std::vector<std::int16_t>& grain_c_integer() const { return grain_c_integer_; }
  Frame process(std::size_t frame, const std::vector<PlaneView>& input) const;
  // Output storage must not overlap input or other output planes. Padding is untouched.
  // Storage is validated before writing; a later processing error may leave partial output.
  void process_into(std::size_t frame, const std::vector<PlaneView>& input,
                    const std::vector<MutablePlaneView>& output) const;

private:
  Parameters params_;
  std::vector<Distance> distances_;
  // Only the input format's representation is populated. Integer noise is
  // rounded and bounded by +/-16384; floating noise is already rounded to F32.
  std::vector<float> grain_y_, grain_c_;
  std::vector<std::int16_t> grain_y_integer_, grain_c_integer_;
  std::vector<std::size_t> offsets_;
};
} // namespace neo_smo::f3kdb
