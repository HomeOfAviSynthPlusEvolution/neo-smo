#pragma once
#include "algorithms/deen.hpp"
#include "base/fp16.hpp"
#include <cmath>
#include <cstring>

namespace neo_smo {
namespace scene_detail {
template <class T>
T load(const std::uint8_t* row, int x) {
  T value;
  std::memcpy(&value, row + static_cast<std::size_t>(x) * sizeof(T), sizeof(T));
  return value;
}
template <class T, bool Half = false>
double difference(const DeenPlane& a, const DeenPlane& b) {
  double sum = 0;
  for (int y = 0; y < a.height; ++y) {
    const auto* ap = a.data + y * a.stride;
    const auto* bp = b.data + y * b.stride;
    for (int x = 0; x < a.width; ++x) {
      if constexpr (Half)
        sum += std::abs(static_cast<double>(fp16_to_fp32(load<T>(ap, x))) - fp16_to_fp32(load<T>(bp, x)));
      else
        sum += std::abs(static_cast<double>(load<T>(ap, x)) - load<T>(bp, x));
    }
  }
  const double peak = a.type == DataType::F16 || a.type == DataType::F32 ? 1.0 : (1u << a.bits) - 1u;
  return sum / (static_cast<double>(a.width) * a.height * peak);
}
} // namespace scene_detail
// Normalized luma MAD, shared by the host adapters. The caller supplies matching planes.
inline double scene_difference(const DeenPlane& a, const DeenPlane& b) {
  switch (a.type) {
    case DataType::U8:
      return scene_detail::difference<std::uint8_t>(a, b);
    case DataType::U16:
      return scene_detail::difference<std::uint16_t>(a, b);
    case DataType::F16:
      return scene_detail::difference<std::uint16_t, true>(a, b);
    case DataType::F32:
      return scene_detail::difference<float>(a, b);
  }
  return 0;
}
} // namespace neo_smo
