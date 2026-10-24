#pragma once
#include "algorithms/f3kdb.hpp"
#include <cstddef>
#include <cstdint>
namespace neo_smo::f3kdb {
Frame process_frame(const Filter&, std::size_t frame, const std::vector<PlaneView>& input);
// Hardware targets use this pipeline; software targets call process_scalar.
Frame process_highway(const Filter&, std::size_t frame, const std::vector<PlaneView>& input);
// Full active row, unpacked to processing units. Gradients only used by mode 7.
struct Row {
  const double* center;
  const double* refs[4];
  const double* gx[5];
  const double* gy[5];
  const float* grain;
  const std::int16_t* integer_grain;
  double* output;
  std::size_t count;
};
// SIMD prefix of neighbour collection for modes 1-3. The caller handles the
// remaining pixels with the checked scalar collector.
struct SamplingRow {
  const double* center;
  const Distance* distances;
  const std::int64_t* vertical_offsets;
  double* samples[4];
  std::size_t width;
  std::size_t vertical_limit;
  int sub_x, sub_y;
};
struct CrossSamplingRow {
  SamplingRow sampling;
  const double* gx;
  const double* gy;
  double* grad_x[4];
  double* grad_y[4];
  std::size_t row_offset;
  // Caller-proven integral Sobel components in [-262140, 262140].
  // Their cross/dot products are exact in double and need no normalization.
  bool bounded_integer_gradients = false;
};
std::size_t collect_cross_row(const CrossSamplingRow&);
struct SamplingResult {
  std::size_t collected, processed;
};
SamplingResult process_cross_row(const CrossSamplingRow&, const Row&, const double*, bool integer, double angle,
                                 double boost);
// Fuse integer sampling/deband where profitable; otherwise collect references.
// The processed prefix leaves sample buffers untouched.
SamplingResult process_sampled_row(const SamplingRow&, const Row&, int mode, double threshold, bool integer, bool blur);
std::size_t collect_simple_row(const SamplingRow&, int mode);
void process_row(const Row&, int mode, const double* thresholds, bool integer, bool blur, double angle, double boost);
// Exact 16-bit processing units in 32-bit lanes. The caller handles the tail.
struct IntegerRow {
  const std::int32_t* center;
  const Distance* distances;
  const std::int16_t* grain;
  double* output;
  std::size_t width, vertical_limit;
  std::size_t row_offset; // center index from the plane origin; plane extent <= INT32_MAX
  int sub_x, sub_y;
  // Optional dither-1 integer output; otherwise write double processing units.
  std::uint8_t* packed_output = nullptr;
  int output_depth = 8;
  bool chroma = false, limited = false;
};
std::size_t process_integer_row(const IntegerRow&, int mode, const double* thresholds, bool blur);
struct IntegerInputPlane {
  const std::uint8_t* input;
  std::int32_t* output;
  std::size_t width, height, stride;
  int depth;
};
bool unpack_integer_plane(const IntegerInputPlane&);
// Exposed internally so tests can check the continuous weight before quantization.
void weights(const double* responses[4], double* output, std::size_t count);
// Contiguous rows of post-grain processing values. boundary contains raw F32
// residuals from the preceding row and is replaced with the last output row's
// residuals. Initialize it to zero once per plane, not once per stripe.
struct FsStripe {
  const double* values;
  std::uint8_t* output;
  float* boundary;
  std::size_t width;
  std::size_t height;
  std::size_t output_stride;
  int output_depth;
  bool integer_input;
  bool chroma;
  bool limited;
};
void diffuse_stripe(const FsStripe&);
struct InputPlane {
  const std::uint8_t* input;
  double* output;
  std::size_t width;
  std::size_t height;
  std::size_t stride;
  int depth;
  bool integer;
};
// False means an out-of-depth integer or a nonfinite float was encountered.
bool unpack_plane(const InputPlane&);
void sobel_plane(const double* source, double* gx, double* gy, std::size_t width, std::size_t height);
struct OutputRow {
  const double* input;
  std::uint8_t* output;
  std::size_t count;
  std::size_t y;
  int depth;
  int dither; // 1: ties-to-even; 2: Bayer. Floating outputs ignore dithering.
  bool integer_input;
  bool integer_output;
  bool chroma;
  bool limited;
};
void pack_row(const OutputRow&);
} // namespace neo_smo::f3kdb
