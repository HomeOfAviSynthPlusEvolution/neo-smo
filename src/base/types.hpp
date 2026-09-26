#pragma once

#include <cstdint>
#include <limits>
#include <type_traits>

namespace neo_smo {

enum class DataType {
  U8,
  U16,
  F16,
  F32,
};

inline DataType get_data_type(int bytes_per_sample, bool is_float) noexcept {
  if (bytes_per_sample == 1) {
    return DataType::U8;
  }
  if (bytes_per_sample == 2) {
    return is_float ? DataType::F16 : DataType::U16;
  }
  return DataType::F32;
}

template <typename T>
inline T get_format_maximum(int bits_per_sample, bool chroma) noexcept {
  if constexpr (std::is_floating_point_v<T>) {
    return static_cast<T>(chroma ? 0.5f : 1.0f);
  } else {
    return static_cast<T>((1u << bits_per_sample) - 1u);
  }
}

template <typename T>
inline T get_format_minimum(bool chroma) noexcept {
  if constexpr (std::is_floating_point_v<T>) {
    return static_cast<T>(chroma ? -0.5f : 0.0f);
  } else {
    return static_cast<T>(0);
  }
}

template <typename T>
inline T get_type_maximum(bool chroma) noexcept {
  if constexpr (std::is_floating_point_v<T>) {
    return static_cast<T>(chroma ? 0.5f : 1.0f);
  } else {
    return std::numeric_limits<T>::max();
  }
}

template <typename T>
inline T get_type_minimum(bool chroma) noexcept {
  if constexpr (std::is_floating_point_v<T>) {
    return static_cast<T>(chroma ? -0.5f : 0.0f);
  } else {
    return static_cast<T>(0);
  }
}

} // namespace neo_smo
