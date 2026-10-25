#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

namespace neo_smo::f3kdb {
// For integral thresholds <=65535, the smallest positive response is at
// least 0.75/65535. Four such responses have a product above 2^-66.
// Index binary32 exponents [-66,0] with 256 linear segments per mantissa.
inline constexpr std::uint32_t integer_weight_first_bits = 61u << 23;
inline constexpr std::size_t integer_weight_table_size = 66 * 256 + 2;
inline const std::array<std::uint64_t, integer_weight_table_size>& integer_weight_table() {
  static const auto table = [] {
    std::array<std::uint64_t, integer_weight_table_size> result{};
    const auto root = [](std::size_t i) {
      const auto bits = integer_weight_first_bits + (static_cast<std::uint32_t>(i) << 15);
      float x;
      std::memcpy(&x, &bits, sizeof(x));
      return static_cast<float>(std::pow(static_cast<double>(x), 0.1));
    };
    float lower = root(0);
    for (std::size_t i = 0; i < result.size(); ++i) {
      const float pair[2] = {lower, root(i + 1)};
      // Preserve memory lane order on both little- and big-endian targets.
      std::memcpy(&result[i], pair, sizeof(pair));
      lower = pair[1];
    }
    return result;
  }();
  return table;
}
} // namespace neo_smo::f3kdb
