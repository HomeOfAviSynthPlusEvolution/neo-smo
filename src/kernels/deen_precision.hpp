#pragma once

#include "common/deen_exact.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace neo_smo::deen_detail {
struct Threshold {
  double hi; // Fast estimate; the original numerator is retained for exact comparisons.
  double raw, peak;
};
inline Threshold scaled_threshold(double threshold, double peak) {
  return {threshold * peak / 255, threshold, peak};
}
inline ThresholdInteger double_magnitude(double value, bool& negative) {
  std::uint64_t bits;
  std::memcpy(&bits, &value, sizeof(bits));
  negative = (bits >> 63) != 0;
  const int exponent = static_cast<int>((bits >> 52) & 2047);
  const auto mantissa = (bits & 0xfffffffffffffu) | (exponent ? (std::uint64_t{1} << 52) : 0);
  return ThresholdInteger(mantissa).shifted(exponent ? exponent - 1 : 0);
}
inline ThresholdInteger exact_difference(double a, double b) {
  bool an, bn;
  auto ai = double_magnitude(a, an), bi = double_magnitude(b, bn);
  if (an != bn) {
    ai += bi;
    return ai;
  }
  return ThresholdInteger::compare(ai, 0, bi, 0) >= 0 ? ai - bi : bi - ai;
}
inline ThresholdInteger threshold_numerator(Threshold threshold) {
  bool sign;
  return double_magnitude(threshold.raw, sign).times(static_cast<std::uint64_t>(threshold.peak));
}
inline bool within_threshold(double sample, double center, Threshold threshold) {
  if (sample == center)
    return true;
  const double difference = std::abs(sample - center);
  const double uncertainty = 8 * std::numeric_limits<double>::epsilon() * std::max(difference, threshold.hi);
  if (std::abs(difference - threshold.hi) > uncertainty)
    return difference < threshold.hi;
  return ThresholdInteger::compare(exact_difference(sample, center).times(255), 0, threshold_numerator(threshold), 0) <=
         0;
}
inline bool adaptive_boundary(double sample, double center, Threshold threshold, double minimum, int dx, int dy,
                              int radius) {
  if (sample == center)
    return true;
  const auto numerator = threshold_numerator(threshold);
  const auto difference = exact_difference(sample, center).times(255);
  if (ThresholdInteger::compare(difference, 0, numerator, 0) > 0)
    return false;
  bool sign;
  const auto falloff = ThresholdInteger(1).shifted(1074) - double_magnitude(minimum, sign);
  const auto remaining = numerator - difference;
  const auto slope = numerator.times(falloff);
  // Compare (t*peak-255*d)^2*q >= (t*peak*(1-min))^2*k.
  // Integer scaling eliminates sqrt and division, preserving exact equality.
  const auto left = remaining.times(remaining).times(static_cast<std::uint64_t>(2 * radius * radius)).shifted(2148);
  const auto right = slope.times(slope).times(static_cast<std::uint64_t>(dx * dx + dy * dy));
  return ThresholdInteger::compare(left, 0, right, 0) >= 0;
}
} // namespace neo_smo::deen_detail
