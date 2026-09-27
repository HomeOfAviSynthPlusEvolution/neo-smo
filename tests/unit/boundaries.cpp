#include "kernels/dispatch.hpp"
#include "base/fp16.hpp"
#include "base/checked.hpp"
#include "common/padded_row.hpp"
#include "hwy/targets.h"
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>
#include "guarded.hpp"

// Independent safe extension of Zig's single reflection for tiny planes.
int reflect(int x, int n) {
  if (n == 1)
    return 0;
  if (x < 0)
    return std::min(-x, n - 1);
  if (x >= n)
    return std::max(0, 2 * n - 2 - x);
  return x;
}
template <class T>
void check(neo_smo::DataType type, int width, int height, bool end) {
  const size_t count = static_cast<size_t>(width) * height;
  Guarded source(count * sizeof(T), end), dest(count * sizeof(T), end);
  auto* src = reinterpret_cast<T*>(source.data);
  auto* dst = reinterpret_cast<T*>(dest.data);
  for (size_t i = 0; i < count; ++i) {
    const auto v = static_cast<unsigned>((i * 71 + i / width * 43) % 251);
    src[i] = type == neo_smo::DataType::F16 ? static_cast<T>(neo_smo::fp32_to_fp16(v / 256.0f)) : static_cast<T>(v);
  }
  const size_t pitch = static_cast<size_t>(width) * sizeof(T);
  for (int radius = 1; radius <= 3; ++radius) {
    neo_smo::process_median_plane(type, radius, source.data, dest.data, width, height, pitch, pitch);
    for (int y = 0; y < height; ++y)
      for (int x = 0; x < width; ++x) {
        std::vector<T> samples;
        for (int dy = -radius; dy <= radius; ++dy)
          for (int dx = -radius; dx <= radius; ++dx)
            samples.push_back(src[reflect(y + dy, height) * width + reflect(x + dx, width)]);
        std::sort(samples.begin(), samples.end());
        if (dst[y * width + x] != samples[samples.size() / 2])
          throw std::runtime_error("median mismatch");
      }
  }
  for (int mode = 1; mode <= 24; ++mode) {
    neo_smo::process_remove_grain_plane(type, mode, false, source.data, dest.data, width, height, pitch, pitch);
    if (mode >= 13 && mode <= 16)
      for (int y = 1; y < height - 1; ++y) {
        if ((y & 1) == ((mode == 13 || mode == 15) ? 1 : 0) && std::memcmp(src + y * width, dst + y * width, pitch))
          throw std::runtime_error("field parity");
      }
  }
  for (int mode = 1; mode <= 2; ++mode)
    if (height >= 2 * mode + 1) {
      neo_smo::process_vertical_cleaner_plane(type, mode, false, sizeof(T) == 1 ? 8 : 16, source.data, dest.data, width,
                                              height, pitch, pitch);
      if (std::memcmp(src, dst, pitch * mode) ||
          std::memcmp(src + (height - mode) * width, dst + (height - mode) * width, pitch * mode))
        throw std::runtime_error("vertical border copy");
    }

  Guarded repair_source(count * sizeof(T), end);
  auto* rep = reinterpret_cast<T*>(repair_source.data);
  for (size_t i = 0; i < count; ++i) {
    const auto v = static_cast<unsigned>((i * 53 + i / width * 29) % 251);
    rep[i] = type == neo_smo::DataType::F16 ? static_cast<T>(neo_smo::fp32_to_fp16(v / 256.0f)) : static_cast<T>(v);
  }
  for (int mode = 1; mode <= 24; ++mode) {
    neo_smo::process_repair_plane(type, mode, false, source.data, repair_source.data, dest.data, width, height, pitch, pitch, pitch);
  }
  Guarded next_source(count * sizeof(T), end);
  auto* next = reinterpret_cast<T*>(next_source.data);
  for (size_t i = 0; i < count; ++i)
    next[i] = rep[count - 1 - i];
  neo_smo::process_clense_plane(type, source.data, repair_source.data, next_source.data, dest.data, width, height, pitch, pitch, pitch, pitch);
  for (size_t i = 0; i < count; ++i)
    if (dst[i] != std::clamp(src[i], std::min(rep[i], next[i]), std::max(rep[i], next[i])))
      throw std::runtime_error("clense median mismatch");
  neo_smo::process_clense_forward_backward_plane(type, source.data, repair_source.data, next_source.data, dest.data, width, height, pitch, pitch, pitch, pitch);
  for (size_t i = 0; i < count; ++i) {
    const auto decode = [type](T value) -> double {
      return type == neo_smo::DataType::F16 ? neo_smo::fp16_to_fp32(static_cast<uint16_t>(value)) : value;
    };
    const double s = decode(src[i]), p = decode(rep[i]), n = decode(next[i]);
    double lo = 2 * std::min(p, n), hi = 2 * std::max(p, n);
    if (type == neo_smo::DataType::F16) {
      lo = neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(static_cast<float>(lo)));
      hi = neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(static_cast<float>(hi)));
    }
    lo -= n;
    hi -= n;
    if (type == neo_smo::DataType::U8 || type == neo_smo::DataType::U16) {
      lo = std::max(0.0, lo);
      hi = std::min(static_cast<double>(std::numeric_limits<T>::max()), hi);
    } else if (type == neo_smo::DataType::F16) {
      lo = neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(static_cast<float>(lo)));
      hi = neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(static_cast<float>(hi)));
    }
    if (decode(dst[i]) != std::clamp(s, lo, hi))
      throw std::runtime_error("clense forward/backward mismatch");
  }

  for (int radius = 1; radius <= 3; ++radius) {
    neo_smo::process_inter_quartile_mean_plane(type, radius, source.data, dest.data, width, height, pitch, pitch);
    neo_smo::process_smart_median_plane(type, radius, 50.0f, source.data, dest.data, width, height, pitch, pitch);
  }

  const uint8_t* t_planes[3] = {source.data, repair_source.data, next_source.data};
  neo_smo::process_temporal_median_plane(type, 3, t_planes, dest.data, width, height, pitch, pitch);
  neo_smo::process_temporal_soften_plane(type, 3, 50.0f, t_planes, dest.data, width, height, pitch, pitch);
  for (int mode = 0; mode <= 4; ++mode) {
    neo_smo::process_temporal_repair_plane(type, mode, false, sizeof(T) == 1 ? 8 : 16, source.data, repair_source.data, repair_source.data, next_source.data, dest.data, width, height, pitch, pitch, pitch, pitch, pitch);
  }
  for (int mode = 0; mode <= 5; ++mode) {
    for (bool interlaced : {false, true})
      for (bool norow : {false, true})
        neo_smo::process_degrain_median_plane(type, mode, 4.0f, interlaced, norow, false, sizeof(T) == 1 ? 8 : 16, repair_source.data, source.data, next_source.data, dest.data, width, height, pitch, pitch, pitch, pitch);
  }
  neo_smo::process_fluxsmooth_t_plane(type, 7.0f, repair_source.data, source.data, next_source.data, dest.data, width, height, pitch, pitch, pitch, pitch);
  neo_smo::process_fluxsmooth_st_plane(type, 7.0f, 7.0f, repair_source.data, source.data, next_source.data, dest.data, width, height, pitch, pitch, pitch, pitch);
}
// Exercise every possible 8/16-bit 3x3 sum through the dispatched mode-20 kernel.
template <class T>
void check_rg_mean(neo_smo::DataType type) {
  constexpr int peak = std::numeric_limits<T>::max();
  constexpr int batch = 4096, width = 3 * batch;
  std::vector<T> src(3 * width), dst(src.size());
  for (int first = 0; first <= 9 * peak; first += batch) {
    const int count = std::min(batch, 9 * peak + 1 - first);
    for (int i = 0; i < count; ++i) {
      const int sum = first + i;
      for (int j = 0; j < 9; ++j)
        src[(j / 3) * width + 3 * i + j % 3] = static_cast<T>(sum / 9 + (j < sum % 9));
    }
    neo_smo::process_remove_grain_plane(type, 20, false,
        reinterpret_cast<const uint8_t*>(src.data()), reinterpret_cast<uint8_t*>(dst.data()),
        width, 3, width * sizeof(T), width * sizeof(T));
    for (int i = 0; i < count; ++i)
      if (dst[width + 3 * i + 1] != (first + i + 4) / 9)
        throw std::runtime_error("RemoveGrain mode 20 integer mean mismatch");
  }
}
template <class T>
void check_rank_clamp(neo_smo::DataType type, int width, bool end) {
  constexpr int height = 3;
  const std::size_t src_pitch = width + 1, ref_pitch = width + 3, dst_pitch = width + 5;
  Guarded source(height * src_pitch * sizeof(T), end), reference(height * ref_pitch * sizeof(T), end), output(height * dst_pitch * sizeof(T), end);
  auto* src = reinterpret_cast<T*>(source.data);
  auto* ref = reinterpret_cast<T*>(reference.data);
  auto* dst = reinterpret_cast<T*>(output.data);
  for (int y = 0; y < height; ++y) for (int x = 0; x < width; ++x) {
    src[y * src_pitch + x] = T((x + y) % 3 == 0 ? std::numeric_limits<T>::max() : (x + y) % 3 == 1 ? 0 : x * 7919 + y * 11);
    ref[y * ref_pitch + x] = T(x * 17113 + y * 23);
  }
  for (bool repair : {false, true}) for (int mode : {1, 2, 3, 4, 11, 12, 13, 14}) {
    if (!repair && mode > 4) continue;
    if (repair) neo_smo::process_repair_plane(type, mode, false, source.data, reference.data, output.data,
        width, height, src_pitch * sizeof(T), ref_pitch * sizeof(T), dst_pitch * sizeof(T));
    else neo_smo::process_remove_grain_plane(type, mode, false, source.data, output.data,
        width, height, src_pitch * sizeof(T), dst_pitch * sizeof(T));
    for (int y = 0; y < height; ++y) for (int x = 0; x < width; ++x) {
      std::array<T, 9> values{};
      int count = 0;
      for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
        if ((!repair || mode >= 12) && dx == 0 && dy == 0) continue;
        const auto* samples = repair ? ref : src;
        const auto stride = repair ? ref_pitch : src_pitch;
        values[count++] = samples[reflect(y + dy, height) * stride + reflect(x + dx, width)];
      }
      std::sort(values.begin(), values.begin() + count);
      const int rank = mode == 11 ? 1 : mode >= 12 ? mode - 10 : mode;
      T lo = values[rank - 1], hi = values[count - rank];
      if (mode >= 12) {
        lo = std::min(lo, ref[y * ref_pitch + x]);
        hi = std::max(hi, ref[y * ref_pitch + x]);
      }
      const T expected = std::clamp(src[y * src_pitch + x], lo, hi);
      if (dst[y * dst_pitch + x] != expected) throw std::runtime_error("Rank clamp sorted oracle mismatch");
    }
  }
}


int main() {
  try {
    (void)neo_smo::checked_product(std::numeric_limits<size_t>::max(), 7);
    return 1;
  } catch (const std::length_error&) {
  }
  if (neo_smo::mirror_index(static_cast<int64_t>(INT32_MAX) + 1, INT32_MAX) != static_cast<size_t>(INT32_MAX - 3)) {
    return 1;
  }
  for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
    hwy::SetSupportedTargetsForTest(target);
    for (bool end : {false, true}) for (int width : {1, 15, 16, 17, 31, 32, 33, 65, 129}) {
      check_rank_clamp<uint8_t>(neo_smo::DataType::U8, width, end);
      check_rank_clamp<uint16_t>(neo_smo::DataType::U16, width, end);
      check_rank_clamp<float>(neo_smo::DataType::F32, width, end);
    }
    check_rg_mean<uint8_t>(neo_smo::DataType::U8);
    check_rg_mean<uint16_t>(neo_smo::DataType::U16);
    for (bool end : {false, true})
      for (int height : {1, 2, 3, 5, 7})
        for (int width : {1, 2, 3, 7, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129}) {
          check<uint8_t>(neo_smo::DataType::U8, width, height, end);
          check<uint16_t>(neo_smo::DataType::U16, width, height, end);
          check<uint16_t>(neo_smo::DataType::F16, width, height, end);
          check<float>(neo_smo::DataType::F32, width, height, end);
        }
    std::printf("%s: guarded boundaries passed\n", hwy::TargetName(target));
  }
  hwy::SetSupportedTargetsForTest(0);
}
