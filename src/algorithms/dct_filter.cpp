#include "algorithms/dct_filter.hpp"
#include "base/checked.hpp"
#include "kernels/dct_pixels.hpp"
#include <algorithm>
#include <cstring>
#include <limits>

namespace neo_smo {
namespace {
void check_dct(neo_dct_status status) {
  if (status != NEO_DCT_OK)
    throw std::runtime_error(std::string("DCTFilter: ") + neo_dct_status_string(status));
}

std::size_t padded_width(std::size_t width) {
  require(width > 0 && width <= std::numeric_limits<std::size_t>::max() - 7, "DCTFilter: invalid plane width.");
  return (width + 7) & ~std::size_t{7};
}
} // namespace

DctFilter::DctFilter(const std::vector<double>& factors) {
  require(factors.size() == 8, "DCTFilter: exactly eight factors are required.");
  for (double factor : factors) {
    // Keep validation effective even when the core is built with /fp:fast.
    std::uint64_t representation;
    std::memcpy(&representation, &factor, sizeof(factor));
    require((representation & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000) && factor >= 0.0 &&
                factor <= 1.0,
            "DCTFilter: factors must be finite and between 0 and 1.");
  }
  for (std::size_t y = 0; y < 8; ++y)
    for (std::size_t x = 0; x < 8; ++x)
      weights_[y * 8 + x] = static_cast<float>(factors[y] * factors[x]);
  neo_dct_plan* plan = nullptr;
  check_dct(neo_dct_plan_create(8, 8, NEO_DCT_HIGHWAY, NEO_DCT_FFTW, &plan));
  plan_.reset(plan);
}

DctScratch::DctScratch(const DctFilter& filter, std::size_t max_width)
    : strip_(checked_product(padded_width(max_width), 8)) {
  neo_dct_workspace* workspace = nullptr;
  check_dct(neo_dct_workspace_create(filter.plan(), &workspace));
  workspace_.reset(workspace);
}

void DctScratch::process(const DctFilter& filter, DataType type, int bits, const std::uint8_t* src, std::uint8_t* dst,
                         std::size_t width, std::size_t height, std::size_t src_stride, std::size_t dst_stride) {
  const auto stride = padded_width(width);
  require(height > 0 && checked_product(stride, 8) <= strip_.size(), "DCTFilter: invalid plane size.");
  const std::size_t bytes = type == DataType::U8 ? 1 : (type == DataType::F32 ? 4 : 2);
  const auto row_bytes = checked_product(width, bytes);
  require(src && dst && src != dst && src_stride >= row_bytes && dst_stride >= row_bytes,
          "DCTFilter: invalid or aliased plane buffers.");
  require((type == DataType::U8 && bits == 8) || (type == DataType::U16 && bits >= 9 && bits <= 16) ||
              (type == DataType::F16 && bits == 16) || (type == DataType::F32 && bits == 32),
          "DCTFilter: unsupported sample format.");
  checked_product(height, src_stride);
  checked_product(height, dst_stride);
  neo_dct_batch batch{strip_.data(), strip_.data(), stride, stride, 8, 8, stride / 8};
  for (std::size_t y = 0; y < height;) {
    const auto rows = std::min(std::size_t{8}, height - y);
    for (std::size_t j = 0; j < 8; ++j) {
      auto* row = strip_.data() + j * stride;
      // Match resize.Point's extension: reflect once, repeating the edge,
      // then clamp to the opposite edge for planes smaller than the padding.
      const auto source_y = j < rows ? y + j : height - 1 - std::min(j - rows, height - 1);
      dct_load_row(type, src + source_y * src_stride, row, width);
      for (std::size_t x = width; x < stride; ++x)
        row[x] = row[width - 1 - std::min(x - width, width - 1)];
    }
    check_dct(neo_dct_filter8(filter.plan(), workspace_.get(), &batch, filter.weights().data()));
    for (std::size_t j = 0; j < rows; ++j)
      dct_store_row(type, strip_.data() + j * stride, dst + (y + j) * dst_stride, width, bits);
    y += rows;
  }
}
} // namespace neo_smo
