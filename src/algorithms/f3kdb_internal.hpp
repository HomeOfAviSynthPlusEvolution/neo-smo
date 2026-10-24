#pragma once
#include "algorithms/f3kdb.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace neo_smo::f3kdb::detail {
[[noreturn]] inline void invalid(const char* message) {
  throw std::invalid_argument(message);
}
inline std::size_t add(std::size_t a, std::size_t b) {
  if (b > std::numeric_limits<std::size_t>::max() - a)
    throw std::overflow_error("size addition");
  return a + b;
}
inline std::size_t mul(std::size_t a, std::size_t b) {
  if (b && a > std::numeric_limits<std::size_t>::max() / b)
    throw std::overflow_error("size multiplication");
  return a * b;
}
template <class T>
void capacity(std::size_t count) {
  (void)mul(count, sizeof(T));
  if (count > std::vector<T>().max_size())
    throw std::overflow_error("vector capacity");
}
inline void bounded(double v, double lo, double hi, const char* name) {
  if (!std::isfinite(v) || v < lo || v > hi)
    invalid(name);
}
inline void valid_format(const Format& f) {
  if (f.color != Color::gray && f.color != Color::yuv)
    invalid("unsupported color family");
  if (f.sample != Sample::integer && f.sample != Sample::floating)
    invalid("unsupported sample type");
  if (f.sample == Sample::integer ? (f.depth < 8 || f.depth > 16) : (f.depth != 16 && f.depth != 32))
    invalid("unsupported depth");
  if (!f.width || !f.height || !f.frames)
    invalid("empty video");
  if (f.sub_x < 0 || f.sub_x > 4 || f.sub_y < 0 || f.sub_y > 4)
    invalid("subsampling exponent");
  if (f.color == Color::gray && (f.sub_x || f.sub_y))
    invalid("Gray cannot be subsampled");
  if ((f.width % (std::size_t{1} << f.sub_x)) || (f.height % (std::size_t{1} << f.sub_y)))
    invalid("nonintegral plane dimensions");
  (void)mul(f.width, f.height);
}
inline std::size_t bytes_per_sample(const Format& f) {
  return f.depth <= 8 ? 1 : f.depth <= 16 ? 2 : 4;
}
// Round each normal operation to F32. An out-of-range intermediate stays F64;
// the subsequent operation resumes F32 only when its result is representable.
inline double fp(double v) {
  return std::abs(v) <= std::numeric_limits<float>::max() ? static_cast<double>(static_cast<float>(v)) : v;
}
inline double away(double x) {
  return x > 0 ? std::floor(x + 0.5) : std::ceil(x - 0.5);
}
inline double uniform_value(std::uint32_t s) {
  const auto m = (std::uint64_t{s} << 20) + (s >> 12);
  return std::ldexp(static_cast<double>(m), -51) - 1;
}
} // namespace neo_smo::f3kdb::detail
