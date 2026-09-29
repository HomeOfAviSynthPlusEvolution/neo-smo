#include "algorithms/mini_deen.hpp"
#include "base/checked.hpp"
#include "kernels/mini_deen.hpp"
#include <limits>
#include <cmath>
#include "algorithms/deen_config.hpp"
namespace neo_smo {
void mini_deen_process_native(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst,
                              std::ptrdiff_t dst_stride, int width, int height, DataType type, int bits, int radius,
                              double threshold) {
  require(src && dst && width > 0 && height > 0, "MiniDeen: invalid plane");
  const bool fp = type == DataType::F16 || type == DataType::F32;
  require((type == DataType::U8 && bits == 8) || (type == DataType::U16 && bits >= 9 && bits <= 16) ||
              (type == DataType::F16 && bits == 16) || (type == DataType::F32 && bits == 32),
          "MiniDeen: invalid format");
  require(radius >= 1 && radius <= 7, "MiniDeen: invalid radius");
  require(deen_finite(threshold) && threshold >= 0 && threshold <= (fp ? 1.0 : (1u << bits) - 1),
          "MiniDeen: invalid threshold");
  const auto bytes = type == DataType::U8 ? 1u : type == DataType::F32 ? 4u : 2u;
  const auto row_bytes = static_cast<std::size_t>(width) * bytes;
  const auto max = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max());
  for (auto stride : {src_stride, dst_stride})
    require(stride > 0 && static_cast<std::size_t>(stride) >= row_bytes &&
                static_cast<std::size_t>(stride) <= max / static_cast<std::size_t>(height),
            "MiniDeen: invalid stride");
  if (fp)
    mini_deen_float_kernel(src, src_stride, dst, dst_stride, width, height, type == DataType::F16, radius,
                           static_cast<float>(threshold));
  else
    mini_deen_kernel(src, src_stride, dst, dst_stride, width, height, bits == 8, radius,
                     static_cast<unsigned>(std::ceil(threshold)));
}
void mini_deen_process(const std::uint8_t* src, std::ptrdiff_t src_stride, std::uint8_t* dst, std::ptrdiff_t dst_stride,
                       int width, int height, int bits, int radius, int threshold) {
  require(src && dst && width > 0 && height > 0, "MiniDeen: invalid plane.");
  require(bits >= 8 && bits <= 16, "MiniDeen: expected 8..16 bit integer samples.");
  require(radius >= 1 && radius <= 7, "MiniDeen: radius must be in [1,7].");
  require(threshold >= 0 && threshold <= 255, "MiniDeen: threshold must be in [0,255].");
  const auto row_bytes = static_cast<std::size_t>(width) * (bits == 8 ? 1u : 2u);
  const auto max = static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max());
  for (auto stride : {src_stride, dst_stride})
    require(stride > 0 && static_cast<std::size_t>(stride) >= row_bytes &&
                static_cast<std::size_t>(stride) <= max / static_cast<std::size_t>(height),
            "MiniDeen: invalid stride or plane size.");
  const unsigned scaled = static_cast<unsigned>(threshold) * ((1u << bits) - 1u) / 255u;
  mini_deen_kernel(src, src_stride, dst, dst_stride, width, height, bits == 8, radius, scaled);
}
} // namespace neo_smo
