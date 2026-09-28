#pragma once

#include <array>
#include <cstdint>
#include <cstring>

namespace neo_smo::deen_detail {
// Bounded unsigned arithmetic for Deen decisions, not a general bigint API.
// A finite F32 sample in units of 2^-149 needs 277 bits. Including the maximum
// validated pixel count, two other plane sizes and the 255 scale needs < 480 bits.
// 640 bits also accommodates the denominator times a binary64 significand.
template <std::size_t Words>
class BoundedInteger {
public:
  explicit BoundedInteger(std::uint64_t value = 0) {
    words_[0] = static_cast<std::uint32_t>(value);
    words_[1] = static_cast<std::uint32_t>(value >> 32);
  }
  BoundedInteger& operator+=(const BoundedInteger& b) {
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < words_.size(); ++i) {
      const std::uint64_t v = static_cast<std::uint64_t>(words_[i]) + b.words_[i] + carry;
      words_[i] = static_cast<std::uint32_t>(v);
      carry = v >> 32;
    }
    return *this;
  }
  BoundedInteger operator-(const BoundedInteger& b) const {
    BoundedInteger out;
    std::uint64_t borrow = 0;
    for (std::size_t i = 0; i < words_.size(); ++i) {
      const std::uint64_t v = static_cast<std::uint64_t>(words_[i]) - b.words_[i] - borrow;
      out.words_[i] = static_cast<std::uint32_t>(v);
      borrow = v >> 63;
    }
    return out;
  }
  BoundedInteger times(std::uint64_t b) const {
    BoundedInteger out;
    for (std::size_t part = 0; part < 2; ++part) {
      const auto factor = static_cast<std::uint32_t>(b >> (32 * part));
      std::uint64_t carry = 0;
      for (std::size_t i = 0; i + part < words_.size(); ++i) {
        const std::uint64_t v = static_cast<std::uint64_t>(words_[i]) * factor + out.words_[i + part] + carry;
        out.words_[i + part] = static_cast<std::uint32_t>(v);
        carry = v >> 32;
      }
    }
    return out;
  }
  BoundedInteger times(const BoundedInteger& b) const {
    BoundedInteger out;
    for (std::size_t i = 0; i < words_.size(); ++i) {
      if (!words_[i])
        continue;
      std::uint64_t carry = 0;
      for (std::size_t j = 0; i + j < words_.size(); ++j) {
        if (!b.words_[j] && !carry)
          continue;
        const std::uint64_t v = static_cast<std::uint64_t>(words_[i]) * b.words_[j] + out.words_[i + j] + carry;
        out.words_[i + j] = static_cast<std::uint32_t>(v);
        carry = v >> 32;
      }
    }
    return out;
  }
  BoundedInteger shifted(int bits) const {
    BoundedInteger out;
    const auto whole = static_cast<std::size_t>(bits / 32);
    const int rest = bits % 32;
    for (std::size_t i = 0; i + whole < words_.size(); ++i) {
      out.words_[i + whole] |= words_[i] << rest;
      if (rest && i + whole + 1 < words_.size())
        out.words_[i + whole + 1] |= words_[i] >> (32 - rest);
    }
    return out;
  }
  int length() const {
    for (int i = static_cast<int>(words_.size()) - 1; i >= 0; --i)
      if (words_[i]) {
        int bits = 0;
        for (auto w = words_[i]; w; w >>= 1)
          ++bits;
        return i * 32 + bits;
      }
    return 0;
  }
  bool bit(int i) const {
    return i >= 0 && i < static_cast<int>(words_.size() * 32) && ((words_[i / 32] >> (i % 32)) & 1u);
  }
  // Compare logical shifts without materializing a potentially 1000-bit shift.
  static int compare(const BoundedInteger& a, int shift_a, const BoundedInteger& b, int shift_b) {
    const int al = a.length(), bl = b.length();
    const int as = al ? al + shift_a : 0, bs = bl ? bl + shift_b : 0;
    if (as != bs)
      return as > bs ? 1 : -1;
    for (int i = as - 1; i >= 0; --i) {
      const bool av = a.bit(i - shift_a), bv = b.bit(i - shift_b);
      if (av != bv)
        return av ? 1 : -1;
    }
    return 0;
  }

private:
  std::array<std::uint32_t, Words> words_{};
};

using SceneInteger = BoundedInteger<20>;
// Threshold squaring uses a 2^-1074 grid for exact binary64 parameters.
// After rejecting d > T, squared operands and scaling need fewer than 4400 bits.
using ThresholdInteger = BoundedInteger<144>;

inline SceneInteger float_magnitude(float value, bool& negative) {
  std::uint32_t bits;
  std::memcpy(&bits, &value, sizeof(bits));
  negative = (bits >> 31) != 0;
  const int exponent = static_cast<int>((bits >> 23) & 255);
  const auto mantissa = (bits & 0x7fffffu) | (exponent ? 0x800000u : 0);
  return SceneInteger(mantissa).shifted(exponent ? exponent - 1 : 0);
}
inline SceneInteger float_difference(float a, float b) {
  bool an, bn;
  auto ai = float_magnitude(a, an), bi = float_magnitude(b, bn);
  if (an != bn) {
    ai += bi;
    return ai;
  }
  return SceneInteger::compare(ai, 0, bi, 0) >= 0 ? ai - bi : bi - ai;
}
inline bool ratio_above(const SceneInteger& numerator, const SceneInteger& denominator, double threshold) {
  std::uint64_t bits;
  std::memcpy(&bits, &threshold, sizeof(bits));
  const int exponent = static_cast<int>((bits >> 52) & 2047);
  const auto mantissa = (bits & 0xfffffffffffffu) | (exponent ? (std::uint64_t{1} << 52) : 0);
  const int shift = exponent ? exponent - 1023 - 52 : -1074;
  const auto right = denominator.times(mantissa);
  return SceneInteger::compare(numerator, shift < 0 ? -shift : 0, right, shift > 0 ? shift : 0) > 0;
}
} // namespace neo_smo::deen_detail
