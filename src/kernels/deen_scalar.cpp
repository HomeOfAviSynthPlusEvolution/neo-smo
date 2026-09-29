#include "kernels/deen_scalar.hpp"
#include "base/fp16.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace neo_smo {
namespace {
template <class T>
T read(const std::uint8_t* p) {
  T v;
  std::memcpy(&v, p, sizeof(v));
  return v;
}
template <class T>
void write(std::uint8_t* p, T v) {
  std::memcpy(p, &v, sizeof(v));
}
int clamp_offset(int x, int dx, int size) {
  return static_cast<int>(std::clamp(std::int64_t{x} + dx, std::int64_t{0}, std::int64_t{size - 1}));
}
unsigned integer_sample(const DeenPlane& p, int x, int y) {
  const auto* row = p.data + y * p.stride;
  return p.type == DataType::U8 ? row[x] : read<std::uint16_t>(row + std::size_t(x) * 2);
}
float float_sample(const DeenPlane& p, int x, int y) {
  const auto* row = p.data + y * p.stride;
  return p.type == DataType::F16 ? fp16_to_fp32(read<std::uint16_t>(row + std::size_t(x) * 2))
                                 : read<float>(row + std::size_t(x) * 4);
}
void validate_integer(const DeenPlane& p) {
  if (p.type == DataType::U8 || p.bits == 16)
    return;
  const unsigned peak = (1u << p.bits) - 1;
  for (int y = 0; y < p.height; ++y)
    for (int x = 0; x < p.width; ++x)
      if (integer_sample(p, x, y) > peak)
        throw std::invalid_argument("Deen: sample exceeds bit depth.");
}
void validate_float(const DeenPlane& p) {
  const bool half = p.type == DataType::F16;
  const unsigned exponent = half ? 0x7c00u : 0x7f800000u;
  for (int y = 0; y < p.height; ++y)
    for (int x = 0; x < p.width; ++x) {
      const auto* ptr = p.data + y * p.stride + std::size_t(x) * (half ? 2 : 4);
      const unsigned bits = half ? read<std::uint16_t>(ptr) : read<std::uint32_t>(ptr);
      if ((bits & exponent) == exponent)
        throw std::invalid_argument("Deen: non-finite sample.");
    }
}
} // namespace

void deen_integer_scalar(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius,
                         double spatial, double temporal, const double* weights, const DeenByteWeights& table,
                         std::uint8_t* dst, std::ptrdiff_t stride) {
  const auto& p = frames[0];
  for (int f = 0; f < count; ++f)
    validate_integer(frames[f]);
  if (family == DeenFamily::Weighted && count == 1 && weights[0] == 1)
    family = DeenFamily::Constant;
  const bool byte = p.type == DataType::U8;
  const int side = 2 * radius + 1;
  std::array<std::array<unsigned, 225>, 2> limits{};
  for (int f = 0; f < 2; ++f)
    for (int i = 0; i < side * side; ++i)
      limits[f][i] = static_cast<unsigned>(
          std::floor((f ? temporal : spatial) * (family == DeenFamily::Adaptive ? weights[i] : 1)));
  for (int y = 0; y < p.height; ++y)
    for (int x = 0; x < p.width; ++x) {
      const unsigned center = integer_sample(p, x, y);
      unsigned sum = family == DeenFamily::Weighted ? (byte ? table.rounding : table.word_rounding) : 0;
      unsigned accepted = 0;
      for (int f = 0; f < count; ++f)
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const unsigned sample = integer_sample(frames[f], clamp_offset(x, dx - radius, p.width),
                                                   clamp_offset(y, dy - radius, p.height));
            const int tap = dy * side + dx;
            const bool pass = std::max(center, sample) - std::min(center, sample) <= limits[f != 0][tap];
            if (family == DeenFamily::Adaptive && !pass)
              continue;
            const unsigned value = pass ? sample : center;
            if (family == DeenFamily::Weighted) {
              const int index = (f != 0 ? 225 : 0) + tap;
              // Match the SIMD fixed-point products and centered rounding bias.
              sum += byte ? (value * table.coefficients[index]) >> 8
                          : static_cast<unsigned>((std::uint64_t(value) * table.word_coefficients[index]) >> 16);
            } else {
              sum += value;
              ++accepted;
            }
          }
      const unsigned result = family == DeenFamily::Weighted ? sum >> 8 : (sum + accepted / 2) / accepted;
      auto* out = dst + y * stride + std::size_t(x) * (byte ? 1 : 2);
      if (byte)
        *out = static_cast<std::uint8_t>(result);
      else
        write(out, static_cast<std::uint16_t>(result));
    }
}

void deen_float_scalar(DeenFamily family, const std::array<DeenPlane, 3>& frames, int count, int radius, double spatial,
                       double temporal, const double* weights, std::uint8_t* dst, std::ptrdiff_t stride) {
  const auto& p = frames[0];
  for (int f = 0; f < count; ++f)
    validate_float(frames[f]);
  const int side = 2 * radius + 1;
  double denominator = 0;
  for (int i = 0; i < side * side; ++i)
    denominator += weights[i];
  denominator = family == DeenFamily::Weighted ? denominator * (count == 3 ? 4 : 1) : side * side * count;
  const float reciprocal = static_cast<float>(1 / denominator);
  for (int y = 0; y < p.height; ++y)
    for (int x = 0; x < p.width; ++x) {
      const float center = float_sample(p, x, y);
      float sums[2]{}, accepted[2]{};
      float low = center, high = center;
      for (int f = 0; f < count; ++f)
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const double weight = weights[dy * side + dx];
            const float coefficient = static_cast<float>(weight * (count == 3 && f == 0 ? 2 : 1));
            if (family == DeenFamily::Weighted && coefficient == 0)
              continue;
            const float sample =
                float_sample(frames[f], clamp_offset(x, dx - radius, p.width), clamp_offset(y, dy - radius, p.height));
            const float difference = sample - center;
            const float limit =
                static_cast<float>((f ? temporal : spatial) * (family == DeenFamily::Adaptive ? weight : 1));
            if (std::abs(difference) > limit)
              continue;
            const int part = dx & 1;
            if (family == DeenFamily::Weighted)
              sums[part] += coefficient * difference;
            else
              sums[part] += difference;
            accepted[part] += 1;
            low = std::min(low, sample);
            high = std::max(high, sample);
          }
      const float residual = sums[0] + sums[1];
      const float mean = family == DeenFamily::Adaptive ? center + residual / (accepted[0] + accepted[1])
                                                        : center + residual * reciprocal;
      const float result = std::clamp(mean, low, high);
      auto* out = dst + y * stride + std::size_t(x) * (p.type == DataType::F16 ? 2 : 4);
      if (p.type == DataType::F16)
        write(out, fp32_to_fp16(result));
      else
        write(out, result);
    }
}

std::uint64_t deen_sad_scalar(const DeenPlane& a, const DeenPlane& b) {
  validate_integer(a);
  validate_integer(b);
  std::uint64_t sum = 0;
  for (int y = 0; y < a.height; ++y)
    for (int x = 0; x < a.width; ++x) {
      const unsigned av = integer_sample(a, x, y), bv = integer_sample(b, x, y);
      sum += std::max(av, bv) - std::min(av, bv);
    }
  return sum;
}

void deen_padded_scalar(DeenFamily family, const std::array<const double*, 3>& src, int count, std::size_t pitch,
                        int width, int height, int radius, double spatial, double temporal, const double* weights,
                        double* dst) {
  const int side = 2 * radius + 1;
  for (int y = 0; y < height; ++y)
    for (int x = 0; x < width; ++x) {
      const double center = src[0][(std::size_t(y) + radius) * pitch + x + radius];
      double sum = 0, denominator = 0, low = center, high = center;
      for (int f = 0; f < count; ++f)
        for (int dy = 0; dy < side; ++dy)
          for (int dx = 0; dx < side; ++dx) {
            const double factor = weights ? weights[dy * side + dx] : 1;
            const double weight = family == DeenFamily::Weighted ? factor * (count == 3 && f == 0 ? 2 : 1) : 1;
            if (weight == 0)
              continue;
            const double sample = src[f][(std::size_t(y) + dy) * pitch + x + dx];
            const double limit = (f ? temporal : spatial) * (family == DeenFamily::Adaptive ? factor : 1);
            const bool pass = std::abs(sample - center) <= limit;
            if (!pass && family == DeenFamily::Adaptive)
              continue;
            const double value = pass ? sample : center;
            sum += family == DeenFamily::Weighted ? weight * (value - center) : value;
            denominator += weight;
            low = std::min(low, value);
            high = std::max(high, value);
          }
      const double mean = family == DeenFamily::Weighted ? center + sum / denominator : sum / denominator;
      dst[std::size_t(y) * width + x] = family == DeenFamily::Constant ? mean : std::clamp(mean, low, high);
    }
}
} // namespace neo_smo
