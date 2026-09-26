#pragma once

#include <stdexcept>
#include <string>
#include <cstddef>
#include <limits>

namespace neo_smo {

inline std::size_t checked_product(std::size_t a, std::size_t b) {
  if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b) {
    throw std::length_error("neo-smo: workspace size overflow");
  }
  return a * b;
}

inline void require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::invalid_argument(message);
  }
}

} // namespace neo_smo
