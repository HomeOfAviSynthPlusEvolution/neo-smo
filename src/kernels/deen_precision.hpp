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
// Error-free subtraction for finite F16/F32/integer inputs promoted to double.
// Their range leaves ample headroom against binary64 overflow and underflow.
inline bool exact_sample_difference(double sample, double center, double& difference) {
  const double delta = sample - center;
  const double virtual_center = sample - delta;
  const double error = (sample - (delta + virtual_center)) + (virtual_center - center);
  difference = std::abs(delta);
  return error == 0;
}
inline bool exact_product(double a, double b, double& product) {
  product = a * b;
  // A 106-bit exact product needs room for its low bits, not only its rounded value.
  // Below 2^-969 a nonzero FMA residual could underflow to zero.
  return (a == 0 || b == 0 || std::abs(product) >= 0x1p-969) && std::fma(a, b, -product) == 0;
}
inline bool exact_linear_comparison(double sample, double center, Threshold threshold, double factor, bool& pass) {
  double difference, lhs, numerator, rhs;
  if (!exact_sample_difference(sample, center, difference) || !exact_product(difference, 255, lhs) ||
      !exact_product(threshold.raw, threshold.peak, numerator) || !exact_product(numerator, factor, rhs))
    return false;
  pass = lhs <= rhs;
  return true;
}
inline bool within_threshold(double sample, double center, Threshold threshold) {
  if (sample == center)
    return true;
  const double difference = std::abs(sample - center);
  const double uncertainty = 8 * std::numeric_limits<double>::epsilon() * std::max(difference, threshold.hi);
  if (std::abs(difference - threshold.hi) > uncertainty)
    return difference < threshold.hi;
  bool pass;
  if (exact_linear_comparison(sample, center, threshold, 1, pass))
    return pass;
  return ThresholdInteger::compare(exact_difference(sample, center).times(255), 0, threshold_numerator(threshold), 0) <=
         0;
}
inline bool adaptive_boundary(double sample, double center, Threshold threshold, double minimum, int dx, int dy,
                              int radius) {
  if (sample == center)
    return true;
  if ((dx == 0 && dy == 0) || minimum == 1)
    return within_threshold(sample, center, threshold);
  if (std::abs(dx) == radius && std::abs(dy) == radius) {
    bool pass;
    if (exact_linear_comparison(sample, center, threshold, minimum, pass))
      return pass;
  }
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
