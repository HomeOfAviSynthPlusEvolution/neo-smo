#include "kernels/f3kdb.hpp"
#include "kernels/f3kdb_scalar.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#ifndef HWY_COMPILE_ALL_ATTAINABLE
#define HWY_COMPILE_ALL_ATTAINABLE
#endif
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/f3kdb_io.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo::f3kdb {
namespace HWY_NAMESPACE {
#if HWY_TARGET == HWY_SCALAR || HWY_TARGET == HWY_EMU128
bool unpack_plane_target(const InputPlane& plane) {
  return unpack_plane_scalar(plane);
}
void pack_row_target(const OutputRow& row) {
  pack_row_scalar(row);
}
void sobel_plane_target(const double* source, double* gx, double* gy, std::size_t width, std::size_t height) {
  sobel_plane_scalar(source, gx, gy, width, height);
}
#else
namespace hn = hwy::HWY_NAMESPACE;
using IoD = hn::ScalableTag<double>;

// Byte loads/stores are also valid for odd addresses/strides on fallback targets.
// LoadN/StoreN bound the physical access, not merely the lanes used afterwards.
template <class D>
HWY_INLINE auto load_bytes(D d, const std::uint8_t* p, std::size_t n) {
  const hn::Repartition<std::uint8_t, D> bytes;
  return hn::BitCast(d, n == hn::Lanes(d) ? hn::LoadU(bytes, p) : hn::LoadN(bytes, p, n * sizeof(hn::TFromD<D>)));
}
template <class D>
HWY_INLINE void store_bytes(D d, hn::VFromD<D> v, std::uint8_t* p, std::size_t n) {
  const hn::Repartition<std::uint8_t, D> bytes;
  if (n == hn::Lanes(d))
    hn::StoreU(hn::BitCast(bytes, v), bytes, p);
  else
    hn::StoreN(hn::BitCast(bytes, v), bytes, p, n * sizeof(hn::TFromD<D>));
}

bool unpack_plane_target(const InputPlane& plane) {
  const IoD d;
  const hn::Rebind<std::uint8_t, IoD> du8;
  const hn::Rebind<std::uint16_t, IoD> du16;
  const hn::Rebind<std::uint32_t, IoD> du32;
  const hn::Rebind<std::uint64_t, IoD> du64;
  const hn::Rebind<float, IoD> df;
  const auto scale = hn::Set(d, std::ldexp(1.0, 16 - plane.depth));
  const std::size_t bytes = plane.depth == 8 ? 1 : plane.depth == 16 || plane.integer ? 2 : 4;
  for (std::size_t y = 0; y < plane.height; ++y)
    for (std::size_t x = 0; x < plane.width; x += hn::Lanes(d)) {
      const auto n = std::min(hn::Lanes(d), plane.width - x);
      const auto* p = plane.input + y * plane.stride + x * bytes;
      auto value = hn::Zero(d);
      if (plane.integer) {
        const auto bits =
            plane.depth == 8 ? hn::PromoteTo(du32, load_bytes(du8, p, n)) : hn::PromoteTo(du32, load_bytes(du16, p, n));
        if (!hn::AllFalse(du32, hn::Gt(bits, hn::Set(du32, (std::uint32_t{1} << plane.depth) - 1))))
          return false;
        value = hn::Mul(hn::PromoteTo(d, bits), scale);
      } else if (plane.depth == 32) {
        const auto bits = load_bytes(du32, p, n);
        // Integer bit checks remain meaningful even under finite-math assumptions.
        if (!hn::AllFalse(du32, hn::Eq(hn::And(bits, hn::Set(du32, 0x7f800000u)), hn::Set(du32, 0x7f800000u))))
          return false;
        value = hn::PromoteTo(d, hn::BitCast(df, bits));
      } else {
        const auto raw = hn::PromoteTo(du64, load_bytes(du16, p, n));
        const auto exponent = hn::And(hn::ShiftRight<10>(raw), hn::Set(du64, 31));
        if (!hn::AllFalse(du64, hn::Eq(exponent, hn::Set(du64, 31))))
          return false;
        const auto mantissa = hn::And(raw, hn::Set(du64, 1023));
        const auto normal = hn::BitCast(
            d, hn::Or(hn::ShiftLeft<52>(hn::Add(exponent, hn::Set(du64, 1008))), hn::ShiftLeft<42>(mantissa)));
        const auto subnormal = hn::Mul(hn::PromoteTo(d, hn::TruncateTo(du32, mantissa)), hn::Set(d, 0x1p-24));
        value = hn::IfThenElse(hn::RebindMask(d, hn::Eq(exponent, hn::Zero(du64))), subnormal, normal);
        value =
            hn::BitCast(d, hn::Or(hn::BitCast(du64, value), hn::ShiftLeft<48>(hn::And(raw, hn::Set(du64, 0x8000)))));
      }
      if (n == hn::Lanes(d))
        hn::StoreU(value, d, plane.output + y * plane.width + x);
      else
        hn::StoreN(value, d, plane.output + y * plane.width + x, n);
    }
  return true;
}

// Convert directly from F64 to binary16. An intermediate F32 would double-round
// values immediately adjacent to half-way points. Input is finite; saturate to
// the largest finite half and preserve signed zero, including underflow to zero.
HWY_INLINE auto encode_half_vector(IoD d, hn::VFromD<IoD> value) {
  const hn::Rebind<std::uint64_t, IoD> du64;
  const hn::Rebind<std::int64_t, IoD> di64;
  const hn::Rebind<std::uint16_t, IoD> du16;
  const auto sign = hn::And(hn::ShiftRight<48>(hn::BitCast(du64, value)), hn::Set(du64, 0x8000));
  const auto magnitude = hn::Min(hn::Abs(value), hn::Set(d, 65504.0));
  const auto field = hn::Max(hn::ShiftRight<52>(hn::BitCast(du64, magnitude)), hn::Set(du64, 1009));
  const auto factor = hn::BitCast(d, hn::ShiftLeft<52>(hn::Sub(hn::Set(du64, 2056), field)));
  const auto significant = hn::BitCast(du64, hn::ConvertTo(di64, hn::Round(hn::Mul(magnitude, factor))));
  const auto normal = hn::Add(hn::ShiftLeft<10>(hn::Sub(field, hn::Set(du64, 1009))), significant);
  const auto subnormal = hn::BitCast(du64, hn::ConvertTo(di64, hn::Round(hn::Mul(magnitude, hn::Set(d, 0x1p24)))));
  const auto bits = hn::IfThenElse(hn::RebindMask(du64, hn::Lt(magnitude, hn::Set(d, 0x1p-14))), subnormal, normal);
  return hn::TruncateTo(du16, hn::Or(sign, bits));
}

void pack_row_target(const OutputRow& row) {
  const IoD d;
  const hn::Rebind<float, IoD> df;
  const hn::Rebind<std::int32_t, IoD> di32;
  const hn::Rebind<std::uint8_t, IoD> du8;
  const hn::Rebind<std::uint16_t, IoD> du16;
  const hn::Rebind<std::uint32_t, IoD> du32;
  const hn::Rebind<std::uint64_t, IoD> du64;
  const double offset = row.chroma ? 0.5 : 0;
  const std::size_t bytes = row.depth == 8 ? 1 : row.depth <= 16 ? 2 : 4;
  const double scale = std::ldexp(1.0, row.depth), tvscale = std::ldexp(1.0, row.depth - 8);
  const auto integer_scale = hn::Set(d, std::ldexp(1.0, row.depth - 16));
  const double lo = row.limited ? 16 * tvscale : 0;
  const double hi = row.limited ? (row.chroma ? 240 : 235) * tvscale : scale - 1;
  for (std::size_t x = 0; x < row.count; x += hn::Lanes(d)) {
    const auto n = std::min(hn::Lanes(d), row.count - x);
    auto value = n == hn::Lanes(d) ? hn::LoadU(d, row.input + x) : hn::LoadN(d, row.input + x, n);
    auto* p = row.output + x * bytes;
    if (!row.integer_output) {
      if (row.integer_input)
        value = hn::Sub(hn::Mul(value, hn::Set(d, 0x1p-16)), hn::Set(d, offset));
      if (row.limited)
        value = hn::Min(hn::Set(d, (row.chroma ? 240.0 : 235.0) / 256 - offset),
                        hn::Max(hn::Set(d, 16.0 / 256 - offset), value));
      if (row.depth == 16)
        store_bytes(du16, encode_half_vector(d, value), p, n);
      else {
        const auto maximum = hn::Set(d, double(std::numeric_limits<float>::max()));
        value = hn::Min(maximum, hn::Max(hn::Neg(maximum), value));
        store_bytes(df, hn::DemoteTo(df, value), p, n);
      }
    } else {
      if (row.integer_input)
        value = hn::Mul(value, integer_scale);
      else
        value = hn::Mul(hn::Add(value, hn::Set(d, offset)), hn::Set(d, scale));
      value = hn::Min(hn::Set(d, hi + 1), hn::Max(hn::Set(d, lo - 1), value));
      if (row.dither == 1)
        value = hn::Round(value);
      else {
        const auto positions = hn::Iota(du64, static_cast<std::uint64_t>(x));
        auto bayer = hn::Zero(du64);
        // b[ybit][xbit] = (2*xbit) XOR (3*ybit), most significant
        // Bayer digit comes from the least significant coordinate bit.
        for (int bit = 0; bit < 4; ++bit) {
          const auto xb = hn::And(hn::ShiftRightSame(positions, bit), hn::Set(du64, 1));
          const auto digit = hn::Xor(hn::ShiftLeft<1>(xb), hn::Set(du64, 3 * ((row.y >> bit) & 1)));
          bayer = hn::Or(bayer, hn::ShiftLeftSame(digit, 2 * (3 - bit)));
        }
        const auto threshold =
            hn::Mul(hn::Add(hn::PromoteTo(d, hn::TruncateTo(du32, bayer)), hn::Set(d, 0.5)), hn::Set(d, 1.0 / 256));
        value = hn::Floor(hn::Add(value, threshold));
      }
      value = hn::Min(hn::Set(d, hi), hn::Max(hn::Set(d, lo), value));
      const auto codes = hn::DemoteTo(du16, hn::DemoteTo(di32, value));
      if (row.depth == 8)
        store_bytes(du8, hn::DemoteTo(du8, codes), p, n);
      else
        store_bytes(du16, codes, p, n);
    }
  }
}

void sobel_plane_target(const double* source, double* gx, double* gy, std::size_t width, std::size_t height) {
  const IoD d;
  const hn::Rebind<std::int64_t, IoD> di;
  const auto zero = hn::Zero(di), last = hn::Set(di, static_cast<std::int64_t>(width - 1));
  const auto distance = hn::Set(di, 20);
  const auto twice = hn::Set(d, 2.0);
  for (std::size_t y = 0; y < height; ++y) {
    const auto ym = y < 20 ? 0 : y - 20, yp = height - 1 - y < 20 ? height - 1 : y + 20;
    const auto* top = source + ym * width;
    const auto* middle = source + y * width;
    const auto* bottom = source + yp * width;
    for (std::size_t x = 0; x < width; x += hn::Lanes(d)) {
      const auto n = std::min(hn::Lanes(d), width - x);
      const auto center = hn::Min(last, hn::Iota(di, static_cast<std::int64_t>(x)));
      const auto left = hn::Max(zero, hn::Sub(center, distance)), right = hn::Min(last, hn::Add(center, distance));
      const auto tl = hn::GatherIndex(d, top, left), tr = hn::GatherIndex(d, top, right);
      const auto ml = hn::GatherIndex(d, middle, left), mr = hn::GatherIndex(d, middle, right);
      const auto bl = hn::GatherIndex(d, bottom, left), br = hn::GatherIndex(d, bottom, right);
      const auto tc = hn::LoadN(d, top + x, n), bc = hn::LoadN(d, bottom + x, n);
      auto horizontal = hn::Add(tr, hn::Mul(twice, mr));
      horizontal = hn::Add(horizontal, br);
      horizontal = hn::Sub(horizontal, tl);
      horizontal = hn::Sub(horizontal, hn::Mul(twice, ml));
      horizontal = hn::Sub(horizontal, bl);
      auto vertical = hn::Add(tl, hn::Mul(twice, tc));
      vertical = hn::Add(vertical, tr);
      vertical = hn::Sub(vertical, bl);
      vertical = hn::Sub(vertical, hn::Mul(twice, bc));
      vertical = hn::Sub(vertical, br);
      hn::StoreN(horizontal, d, gx + y * width + x, n);
      hn::StoreN(vertical, d, gy + y * width + x, n);
    }
  }
}
#endif
} // namespace HWY_NAMESPACE
} // namespace neo_smo::f3kdb
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo::f3kdb {
HWY_EXPORT(unpack_plane_target);
HWY_EXPORT(pack_row_target);
HWY_EXPORT(sobel_plane_target);
bool unpack_plane(const InputPlane& plane) {
  return HWY_DYNAMIC_DISPATCH(unpack_plane_target)(plane);
}
void pack_row(const OutputRow& row) {
  HWY_DYNAMIC_DISPATCH(pack_row_target)(row);
}
void sobel_plane(const double* source, double* gx, double* gy, std::size_t width, std::size_t height) {
  HWY_DYNAMIC_DISPATCH(sobel_plane_target)(source, gx, gy, width, height);
}
} // namespace neo_smo::f3kdb
#endif
