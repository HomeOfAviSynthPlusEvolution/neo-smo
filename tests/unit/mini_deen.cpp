#include "algorithms/mini_deen.hpp"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>
using namespace neo_smo;
namespace {
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("MiniDeen test failed");
}
unsigned load(const std::uint8_t* p, int bytes) {
  if (bytes == 1)
    return *p;
  std::uint16_t v;
  std::memcpy(&v, p, 2);
  return v;
}
void store(std::uint8_t* p, unsigned v, int bytes) {
  if (bytes == 1)
    *p = static_cast<std::uint8_t>(v);
  else {
    const auto t = static_cast<std::uint16_t>(v);
    std::memcpy(p, &t, 2);
  }
}
void run(int bits, int w, int h, int r, int t) {
  const int bytes = bits == 8 ? 1 : 2, ss = w * bytes + 3, ds = w * bytes + 7;
  const unsigned max = (1u << bits) - 1, threshold = t * max / 255;
  std::vector<std::uint8_t> input(ss * h + 2, 0xa5), output(ds * h + 2, 0xcc);
  auto* src = input.data() + 1;
  auto* dst = output.data() + 1;
  std::vector<unsigned> values(w * h);
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x) {
      // Small differences exercise t=1 in high depths; alternating maxima stress sums.
      unsigned v = (x + y) % 5 == 0 ? max : (max / 2 + (x * 7 + y * 3) % 19);
      values[y * w + x] = v;
      store(src + y * ss + x * bytes, v, bytes);
    }
  const auto original = input;
  mini_deen_process(src, ss, dst, ds, w, h, bits, r, t);
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x) {
      // Independent absolute-coordinate reference; reject positions outside the image.
      unsigned sum = 2 * values[y * w + x], count = 2;
      for (int yy = 0; yy < h; ++yy)
        for (int xx = 0; xx < w; ++xx) {
          if (std::abs(yy - y) > r || std::abs(xx - x) > r)
            continue;
          const unsigned v = values[yy * w + xx], c = values[y * w + x];
          if ((v > c ? v - c : c - v) < threshold) {
            sum += v;
            ++count;
          }
        }
      check(load(dst + y * ds + x * bytes, bytes) == (sum + count / 2) / count);
    }
  check(input == original && output.front() == 0xcc && output.back() == 0xcc);
  for (int y = 0; y < h; ++y)
    for (int x = w * bytes; x < ds; ++x)
      check(dst[y * ds + x] == 0xcc);
}
} // namespace
int main() {
  try {
    for (int bits = 8; bits <= 16; ++bits)
      for (int r = 1; r <= 7; ++r)
        for (int t : {0, 1, 2, 10, 255})
          for (int w : {1, 2, 7, 17})
            run(bits, w, w == 1 ? 9 : 3, r, t);
    // Explicit center weighting and strict equality: [100,110], threshold=10 rejects,
    // threshold=11 yields (3*100+110)/4=102.5 -> 103 at the first pixel.
    std::uint8_t src[]{100, 110}, dst[2]{};
    mini_deen_process(src, 2, dst, 2, 2, 1, 8, 1, 10);
    check(dst[0] == 100 && dst[1] == 110);
    mini_deen_process(src, 2, dst, 2, 2, 1, 8, 1, 11);
    check(dst[0] == 103 && dst[1] == 108);
    for (int bad : {0, 8}) {
      bool threw = false;
      try {
        mini_deen_process(src, 2, dst, 2, 2, 1, 8, bad, 10);
      } catch (const std::invalid_argument&) {
        threw = true;
      }
      check(threw);
    }
    bool threw = false;
    try {
      mini_deen_process(src, 1, dst, 2, 2, 1, 8, 1, 10);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    check(threw);
    std::puts("MiniDeen scalar checks passed");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
