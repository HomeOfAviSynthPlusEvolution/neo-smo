#pragma once
#include "algorithms/f3kdb_internal.hpp"
#include <dualsynth/span2d.hpp>

namespace neo_smo::f3kdb::detail {
// span2d stores dimensions and element strides in signed 32-bit fields.
// Check before narrowing; byte views preserve odd host addresses and strides.
template <class T>
span2d::Plane<T> plane_view(T* data, std::size_t width, std::size_t height, std::size_t stride_bytes) {
  constexpr auto limit = static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
  if (width > limit || height > limit || stride_bytes % sizeof(T) || stride_bytes / sizeof(T) > limit)
    throw std::overflow_error("F3kdb: plane exceeds span2d limits");
  const auto row_bytes = mul(width, sizeof(T));
  if (stride_bytes < row_bytes)
    invalid("F3kdb: plane stride smaller than row");
  const auto extent = height ? add(mul(height - 1, stride_bytes), row_bytes) : 0;
  if (stride_bytes > static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()) ||
      extent > static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max()))
    throw std::overflow_error("F3kdb: plane address range");
  return {data, static_cast<std::int32_t>(width), static_cast<std::int32_t>(height),
          static_cast<std::ptrdiff_t>(stride_bytes)};
}
template <class T>
span2d::Plane<T> plane_view(T* data, std::size_t width, std::size_t height) {
  return plane_view(data, width, height, mul(width, sizeof(T)));
}
} // namespace neo_smo::f3kdb::detail
