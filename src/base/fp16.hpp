#pragma once

#include <cstdint>
#include <cstring>

namespace neo_smo {

// Software IEEE-754 binary16 <-> binary32 conversion matching hardware F16C / Zig f16 semantics.
inline float fp16_to_fp32(std::uint16_t h) noexcept {
  const std::uint32_t sign = static_cast<std::uint32_t>(h & 0x8000u) << 16;
  const std::uint32_t exp = (h >> 10) & 0x1Fu;
  const std::uint32_t mant = h & 0x03FFu;

  std::uint32_t bits = 0;
  if (exp == 0) {
    if (mant == 0) {
      bits = sign;
    } else {
      // Subnormal fp16 -> normalized fp32
      std::uint32_t m = mant;
      std::uint32_t e = 127 - 14;
      while ((m & 0x0400u) == 0) {
        m <<= 1;
        --e;
      }
      m &= 0x03FFu;
      bits = sign | (e << 23) | (m << 13);
    }
  } else if (exp == 31) {
    bits = sign | 0x7F800000u | (mant << 13);
  } else {
    bits = sign | ((exp + (127 - 15)) << 23) | (mant << 13);
  }

  float out = 0.0f;
  std::memcpy(&out, &bits, sizeof(float));
  return out;
}

inline std::uint16_t fp32_to_fp16(float f) noexcept {
  std::uint32_t bits = 0;
  std::memcpy(&bits, &f, sizeof(float));

  const std::uint32_t sign = (bits >> 16) & 0x8000u;
  const std::int32_t exp = static_cast<std::int32_t>((bits >> 23) & 0xFFu) - 127 + 15;
  std::uint32_t mant = bits & 0x007FFFFFu;

  if (exp <= 0) {
    if (exp < -10) {
      return static_cast<std::uint16_t>(sign);
    }
    mant |= 0x00800000u;
    const std::uint32_t shift = static_cast<std::uint32_t>(14 - exp);
    const std::uint32_t round_bit = 1u << (shift - 1);
    const std::uint32_t sticky_mask = round_bit - 1u;
    std::uint32_t half_mant = mant >> shift;
    if ((mant & round_bit) != 0 && ((mant & sticky_mask) != 0 || (half_mant & 1u) != 0)) {
      ++half_mant;
    }
    return static_cast<std::uint16_t>(sign | half_mant);
  } else if (exp == 0xFF - 127 + 15) {
    if (mant == 0) {
      return static_cast<std::uint16_t>(sign | 0x7C00u);
    }
    return static_cast<std::uint16_t>(sign | 0x7C00u | (mant >> 13) | 1u);
  } else if (exp > 30) {
    return static_cast<std::uint16_t>(sign | 0x7C00u);
  }

  // Round to nearest, ties to even
  std::uint32_t half = sign | (static_cast<std::uint32_t>(exp) << 10) | (mant >> 13);
  const std::uint32_t remainder = mant & 0x1FFFu;
  if (remainder > 0x1000u || (remainder == 0x1000u && (half & 1u) != 0)) {
    ++half;
  }
  return static_cast<std::uint16_t>(half);
}

} // namespace neo_smo
