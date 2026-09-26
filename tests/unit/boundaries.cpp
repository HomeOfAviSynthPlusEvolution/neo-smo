#include "kernels/dispatch.hpp"
#include "base/fp16.hpp"
#include "base/checked.hpp"
#include "common/padded_row.hpp"
#include "hwy/targets.h"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

struct Guarded {
  unsigned char* base;
  unsigned char* data;
  size_t allocated;
  Guarded(size_t bytes, bool at_end) {
#ifdef _WIN32
    SYSTEM_INFO si{};
    GetSystemInfo(&si);
    const size_t page = si.dwPageSize;
#else
    const size_t page = static_cast<size_t>(sysconf(_SC_PAGESIZE));
#endif
    const size_t usable = (bytes + page - 1) / page * page;
    allocated = usable + 2 * page;
#ifdef _WIN32
    base = static_cast<unsigned char*>(VirtualAlloc(nullptr, allocated, MEM_RESERVE, PAGE_NOACCESS));
    if (!base || !VirtualAlloc(base + page, usable, MEM_COMMIT, PAGE_READWRITE))
      throw std::runtime_error("VirtualAlloc");
#else
    base = static_cast<unsigned char*>(mmap(nullptr, allocated, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
    if (base == MAP_FAILED || mprotect(base + page, usable, PROT_READ | PROT_WRITE))
      throw std::runtime_error("mmap");
#endif
    data = base + page + (at_end ? usable - bytes : 0);
  }
  ~Guarded() {
#ifdef _WIN32
    VirtualFree(base, 0, MEM_RELEASE);
#else
    munmap(base, allocated);
#endif
  }
};

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
