#include "algorithms/deen.hpp"
#include "base/checked.hpp"
#include "base/fp16.hpp"
#include "kernels/deen.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <utility>

namespace neo_smo {
namespace {
// Largest-remainder quantization preserves each frame's total weight. Its
// coefficient L1 error is at most half a coefficient unit per spatial tap.
template <class T>
void quantize_coefficients(const std::vector<double>& weights, int total, T* dst) {
  const int taps = static_cast<int>(weights.size());
  std::array<int, 225> order{};
  std::array<double, 225> fractions{};
  double denominator = 0;
  for (double weight : weights)
    denominator += weight;
  int remainder = total;
  for (int i = 0; i < taps; ++i) {
    const double scaled = weights[i] * (total / denominator);
    const int whole = static_cast<int>(std::floor(scaled));
    dst[i] = static_cast<T>(whole);
    fractions[i] = scaled - whole;
    order[i] = i;
    remainder -= whole;
  }
  require(remainder >= 0 && remainder <= taps, "Deen: invalid weight normalization.");
  std::sort(order.begin(), order.begin() + taps,
            [&](int a, int b) { return fractions[a] != fractions[b] ? fractions[a] > fractions[b] : a < b; });
  for (int i = 0; i < remainder; ++i)
    ++dst[order[i]];
}
std::size_t product(std::size_t a, std::size_t b) {
  require(!b || a <= static_cast<std::size_t>(PTRDIFF_MAX) / b, "Deen: plane size overflow.");
  return a * b;
}
std::size_t sample_bytes(DataType t) {
  return t == DataType::U8 ? 1 : t == DataType::F32 ? 4 : 2;
}
void validate(const DeenPlane& p) {
  require(p.data && p.width > 0 && p.height > 0, "Deen: invalid plane.");
  require((p.type == DataType::U8 && p.bits == 8) || (p.type == DataType::U16 && p.bits >= 9 && p.bits <= 16) ||
              (p.type == DataType::F16 && p.bits == 16) || (p.type == DataType::F32 && p.bits == 32),
          "Deen: unsupported sample format.");
  require(p.stride > 0 && static_cast<std::size_t>(p.stride) >= product(p.width, sample_bytes(p.type)),
          "Deen: invalid stride.");
  product(static_cast<std::size_t>(p.stride), p.height);
}
void matching(const DeenPlane& a, const DeenPlane& b) {
  validate(b);
  require(a.width == b.width && a.height == b.height && a.type == b.type && a.bits == b.bits,
          "Deen: incompatible reference plane.");
}
double read(const DeenPlane& p, int x, int y) {
  const auto* v = p.data + y * p.stride + x * sample_bytes(p.type);
  if (p.type == DataType::U8)
    return *v;
  if (p.type == DataType::F32) {
    float f;
    std::memcpy(&f, v, sizeof(f));
    // Bit checks remain valid under the production compiler's fast FP assumptions.
    std::uint32_t u;
    std::memcpy(&u, v, sizeof(u));
    if ((u & 0x7f800000u) == 0x7f800000u)
      throw std::invalid_argument("Deen: non-finite sample.");
    return f;
  }
  std::uint16_t u;
  std::memcpy(&u, v, sizeof(u));
  if (p.type == DataType::F16) {
    if ((u & 0x7c00u) == 0x7c00u)
      throw std::invalid_argument("Deen: non-finite sample.");
    return fp16_to_fp32(u);
  }
  if (u > ((1u << p.bits) - 1))
    throw std::invalid_argument("Deen: sample exceeds bit depth.");
  return u;
}
void write(std::uint8_t* dst, DataType type, double value) {
  if (type == DataType::U8) {
    *dst = static_cast<std::uint8_t>(std::floor(value + 0.5));
    return;
  }
  if (type == DataType::U16) {
    const auto v = static_cast<std::uint16_t>(std::floor(value + 0.5));
    std::memcpy(dst, &v, sizeof(v));
  } else if (type == DataType::F16) {
    const auto v = fp32_to_fp16(static_cast<float>(value));
    std::memcpy(dst, &v, sizeof(v));
  } else {
    const auto v = static_cast<float>(value);
    std::memcpy(dst, &v, sizeof(v));
  }
}
} // namespace

Deen::Deen(DeenOptions options) : options_(std::move(options)) {
  const auto& p = options_;
  require(p.mode == "c2d" || p.mode == "c3d" || p.mode == "w2d" || p.mode == "w3d" || p.mode == "a2d" ||
              p.mode == "a3d",
          "Deen: invalid mode.");
  temporal_ = p.mode[1] == '3';
  family_ = p.mode[0] == 'c' ? DeenFamily::Constant : p.mode[0] == 'w' ? DeenFamily::Weighted : DeenFamily::Adaptive;
  require(p.radius >= 1 && p.radius <= (temporal_ ? 4 : 7), "Deen: invalid radius.");
  for (double t : {p.spatial_y, p.spatial_uv, p.temporal_y, p.temporal_uv})
    require(std::isfinite(t) && t >= 0 && t <= 255, "Deen: invalid threshold.");
  require(std::isfinite(p.minimum) && p.minimum >= 0 && p.minimum <= 1, "Deen: invalid min.");
  require(std::isfinite(p.scene_threshold) && p.scene_threshold >= 0, "Deen: invalid scd.");
  const int r = p.radius;
  for (int y = -r; y <= r; ++y)
    for (int x = -r; x <= r; ++x) {
      const double rho = std::sqrt(static_cast<double>(x * x + y * y)) / std::sqrt(2.0 * r * r);
      weights_.push_back(x == 0 && y == 0                       ? 1
                         : std::abs(x) == r && std::abs(y) == r ? p.minimum
                                                                : 1 - (1 - p.minimum) * rho);
    }
  if (family_ == DeenFamily::Weighted) {
    quantize_coefficients(weights_, 32768, byte_weights_[0].coefficients.data());
    if (temporal_) {
      quantize_coefficients(weights_, 16384, byte_weights_[1].coefficients.data());
      quantize_coefficients(weights_, 8192, byte_weights_[1].coefficients.data() + 225);
    }
    for (int temporal = 0; temporal <= static_cast<int>(temporal_); ++temporal) {
      auto& table = byte_weights_[temporal];
      quantize_coefficients(weights_, temporal ? (1 << 23) : (1 << 24), table.word_coefficients.data());
      if (temporal)
        quantize_coefficients(weights_, 1 << 22, table.word_coefficients.data() + 225);
      unsigned word_span = 0;
      for (int frame = 0; frame <= temporal; ++frame)
        for (std::size_t i = 0; i < weights_.size(); ++i)
          word_span += (65536u - std::gcd(table.word_coefficients[frame * 225 + i], 65536u)) * (frame == 1 ? 2 : 1);
      table.word_rounding = 128 + (word_span + 65536) / 131072;
      unsigned truncation_span = 0;
      for (int frame = 0; frame <= temporal; ++frame)
        for (std::size_t i = 0; i < weights_.size(); ++i) {
          auto& coefficient = table.coefficients[frame * 225 + i];
          // The fractional part of sample*coefficient/128 ranges from zero
          // to (128-gcd(coefficient,128))/128. Center that interval once per
          // output, then double Q15 coefficients for unsigned high multiply.
          truncation_span += (128u - std::gcd(unsigned(coefficient), 128u)) * (frame == 1 ? 2 : 1);
          coefficient *= 2;
        }
      table.rounding = static_cast<std::uint16_t>(128 + (truncation_span + 128) / 256);
    }
  }
}

void deen_process(const Deen& filter, bool chroma, bool temporal, const std::array<DeenPlane, 3>& frames,
                  std::uint8_t* dst, std::ptrdiff_t stride) {
  const auto& p = frames[0];
  validate(p);
  require(!temporal || filter.temporal(), "Deen: invalid temporal evaluation.");
  require(dst && stride > 0 && static_cast<std::size_t>(stride) >= product(p.width, sample_bytes(p.type)),
          "Deen: invalid output stride.");
  product(static_cast<std::size_t>(stride), p.height);
  const int count = temporal ? 3 : 1, r = filter.options().radius;
  if (p.type == DataType::U8 || p.type == DataType::U16) {
    for (int f = 0; f < count; ++f)
      matching(p, frames[f]);
    const auto& o = filter.options();
    const double peak = (1u << p.bits) - 1;
    deen_integer_kernel(filter.family(), frames, count, r, (chroma ? o.spatial_uv : o.spatial_y) * peak / 255,
                        (chroma ? o.temporal_uv : o.temporal_y) * peak / 255, filter.weights().data(),
                        filter.byte_weights(temporal), dst, stride);
    return;
  }
  if (filter.family() == DeenFamily::Weighted) {
    for (int f = 0; f < count; ++f)
      matching(p, frames[f]);
    const auto& o = filter.options();
    deen_weighted_float_kernel(frames, count, r, (chroma ? o.spatial_uv : o.spatial_y) / 255,
                               (chroma ? o.temporal_uv : o.temporal_y) / 255, filter.weights().data(), dst, stride);
    return;
  }
  const std::size_t pitch = static_cast<std::size_t>(p.width) + 2 * r;
  const auto length = product(pitch, static_cast<std::size_t>(p.height) + 2 * r);
  product(length, sizeof(double));
  std::array<std::vector<double>, 3> padded;
  std::array<const double*, 3> sources{};
  for (int f = 0; f < count; ++f) {
    matching(p, frames[f]);
    auto& buf = padded[f];
    buf.resize(length);
    for (int y = 0; y < p.height; ++y) {
      auto* row = buf.data() + (static_cast<std::size_t>(y) + r) * pitch;
      for (int x = 0; x < p.width; ++x)
        row[static_cast<std::size_t>(x) + r] = read(frames[f], x, y);
      std::fill_n(row, r, row[r]);
      std::fill_n(row + r + p.width, r, row[static_cast<std::size_t>(r) + p.width - 1]);
    }
    for (int y = 0; y < r; ++y) {
      std::copy_n(buf.data() + r * pitch, pitch, buf.data() + y * pitch);
      std::copy_n(buf.data() + (static_cast<std::size_t>(r) + p.height - 1) * pitch, pitch,
                  buf.data() + (static_cast<std::size_t>(r) + p.height + y) * pitch);
    }
    sources[f] = buf.data();
  }
  const auto out_size = product(p.width, p.height);
  product(out_size, sizeof(double));
  std::vector<double> output(out_size);
  const auto& o = filter.options();
  const double peak = p.type == DataType::U8 || p.type == DataType::U16 ? (1u << p.bits) - 1 : 1;
  deen_kernel(filter.family(), sources, count, pitch, p.width, p.height, r,
              (chroma ? o.spatial_uv : o.spatial_y) * peak / 255, (chroma ? o.temporal_uv : o.temporal_y) * peak / 255,
              filter.weights().data(), output.data());
  for (int y = 0; y < p.height; ++y)
    for (int x = 0; x < p.width; ++x)
      write(dst + y * stride + x * sample_bytes(p.type), p.type, output[static_cast<std::size_t>(y) * p.width + x]);
}

bool deen_scene_cut(const Deen& filter, const std::vector<DeenPlane>& a, const std::vector<DeenPlane>& b) {
  if (!filter.options().scenechange)
    return false;
  require((a.size() == 1 || a.size() == 3) && a.size() == b.size(), "Deen: invalid scene planes.");
  double metric = 0;
  for (std::size_t p = 0; p < a.size(); ++p) {
    validate(a[p]);
    matching(a[p], b[p]);
    require(a[p].type == a[0].type && a[p].bits == a[0].bits, "Deen: mixed scene formats.");
    const bool integer = a[p].type == DataType::U8 || a[p].type == DataType::U16;
    const auto samples = product(a[p].width, a[p].height);
    const double peak = integer ? (1u << a[p].bits) - 1 : 1;
    double sum = 0;
    if (integer) {
      require(samples <= UINT64_MAX / static_cast<std::uint64_t>(peak), "Deen: scene sum overflow.");
      sum = static_cast<double>(deen_integer_sad(a[p], b[p]));
    } else {
      double correction = 0;
      for (int y = 0; y < a[p].height; ++y)
        for (int x = 0; x < a[p].width; ++x) {
          const double difference = std::abs(read(a[p], x, y) - read(b[p], x, y));
          const double adjusted = difference - correction;
          const double next = sum + adjusted;
          correction = (next - sum) - adjusted;
          sum = next;
        }
    }
    metric += (sum / static_cast<double>(samples)) * (255 / peak);
  }
  metric /= static_cast<double>(a.size());
  return metric > filter.options().scene_threshold;
}
} // namespace neo_smo
