#include "kernels/dispatch.hpp"
#include "base/fp16.hpp"
#include "guarded.hpp"
#include "hwy/targets.h"
#include "hwy/per_target.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

int reflect(int x, int n) {
  if (n <= 1) return 0;
  if (x < 0) return std::min(-x, n-1);
  return x < n ? x : std::max(0, 2*n-2-x);
}

template<class T> void check(neo_smo::DataType type, int width, int height, bool end) {
  const bool half = type == neo_smo::DataType::F16;
  const bool integer = type == neo_smo::DataType::U8 || type == neo_smo::DataType::U16;
  const int bits = sizeof(T) == 1 ? 8 : 16;
  const size_t size = static_cast<size_t>(width) * height, pitch = width * sizeof(T);
  std::vector<std::unique_ptr<Guarded>> storage;
  std::array<const uint8_t*, 9> src{};
  std::array<uint8_t*, 3> dst{};
  auto decode = [&](T v) { return half ? neo_smo::fp16_to_fp32(static_cast<uint16_t>(v)) : double(v); };
  auto encode = [&](double v) { return half ? T(neo_smo::fp32_to_fp16(float(v))) : T(v); };
  auto round = [&](double v) { return half && hwy::HaveFloat16() ? double(neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(float(v)))) : integer ? v : double(float(v)); };
  for (int frame = 0; frame < 3; ++frame) for (int p = 0; p < 3; ++p) {
    storage.push_back(std::make_unique<Guarded>(size*sizeof(T), end));
    src[frame*3+p] = storage.back()->data;
    auto* data = reinterpret_cast<T*>(storage.back()->data);
    for (size_t i=0;i<size;++i) {
      const uint32_t value = uint32_t(i*17113 + frame*541 + p*211);
      data[i] = encode(integer ? double(value % (1u << bits)) : 0.48 + double(value % 129)/4096);
    }
  }
  for (int p=0;p<3;++p) {
    storage.push_back(std::make_unique<Guarded>(size*sizeof(T), end));
    dst[p] = storage.back()->data;
  }
  const neo_smo::Point points[] = {{-2,-1},{1,-1},{0,1},{2,1},
                                   {0,0},{-1,0},{1,0},{0,-1}};
  const float weights[] = {0.70710677f, 1.f, 0.68163878f};
  for (bool rgb : {false,true}) for (int radius : {0,1}) for (int num_points : {4,8}) {
    if (num_points == 8 && type != neo_smo::DataType::F32 && type != neo_smo::DataType::U8) continue;
    const float threshold = integer ? 123456.f : 0.0003f;
    const auto* sources = radius ? src.data() : src.data()+3;
    neo_smo::process_ccd_planes(type,rgb,width,height,pitch,sources,sources,dst[0],dst[1],dst[2],threshold,radius,
                               radius?weights:weights+1,points,num_points,5,1.f,bits);
    auto sample = [&](int f, int p, int x, int y) {return decode(reinterpret_cast<const T*>(src[f*3+p])[reflect(y,height)*width+reflect(x,width)]);};
    for(int y=0;y<height;++y) for(int x=0;x<width;++x) {
      double totals[3] = {sample(1,0,x,y),sample(1,1,x,y),sample(1,2,x,y)};
      int accepted=1;
      for(int point_index=0;point_index<num_points;++point_index) {
        const auto point=points[point_index];
        auto distance = [&](int f) {
          double squares[3];
          for(int p=0;p<3;++p) { const double diff=round(sample(f,p,x+point.x,y+point.y)-sample(1,p,x,y)); squares[p]=round(diff*diff); }
          if(!rgb) squares[0]=round(squares[0]*4);
          return round(round(squares[0]+squares[1])+squares[2]);
        };
        double ssd=distance(1);
        if(radius) {
          if(integer) ssd=std::floor((ssd+std::round(float(distance(0))*weights[0]+float(distance(2))*weights[2])+1)/3);
          else ssd=round(round(ssd+round(round(distance(0)*round(weights[0]))+round(distance(2)*round(weights[2]))))/3);
        }
        if(ssd<round(threshold)) {
          for(int p=0;p<3;++p) totals[p]=round(totals[p]+sample(1,p,x+point.x,y+point.y));
          ++accepted;
        }
      }
      for(int p=rgb?0:1;p<3;++p) {
        double expected = integer ? std::round(float(totals[p])/float(accepted)) : round(totals[p]/accepted);
        const double got=decode(reinterpret_cast<T*>(dst[p])[y*width+x]);
        // FP32 reciprocal/FMA rounding can cross a final binary16 midpoint.
        const double tolerance = integer || (half && hwy::HaveFloat16()) ? 0 :
            half ? std::max(std::abs(expected) * 0x1p-10, 0x1p-24) : 1e-6;
        if(std::abs(got-decode(encode(expected))) > tolerance) {
          std::fprintf(stderr,"CCD type=%d rgb=%d radius=%d width=%d xy=%d,%d got=%g expected=%g\n",int(type),rgb,radius,width,x,y,got,expected);
          throw std::runtime_error("CCD scalar oracle mismatch");
        }
      }
    }
  }
  if(integer) {
    std::array<uint8_t,256> table{}; for(int i=0;i<256;++i) table[i]=uint8_t(255-i);
    std::array<const uint8_t*,3> neighbors[20];
    for (int radius : {1, 2, 10}) {
      for (int i = 0; i < 2 * radius; ++i) {
        const int frame = i < radius ? 0 : 2;
        neighbors[i] = {src[frame * 3], src[frame * 3 + 1], src[frame * 3 + 2]};
      }
      for(int mode=0;mode<4;++mode) {
        neo_smo::process_cnr4_frame(type,width,height,pitch,bits,radius,0,mode,src.data()+3,src.data()+3,neighbors,neighbors,2*radius,dst[1],dst[2],table.data(),table.data(),table.data());
        const uint64_t max=uint64_t{1}<<(bits*2);
        for(size_t i=0;i<size;++i) for(int p=1;p<3;++p) {
          uint64_t sum=0;
          auto value=[&](int f,int plane) {return uint64_t(reinterpret_cast<const T*>(src[f*3+plane])[i]);};
          for(int j = 0; j < 2 * radius; ++j) {
            const int f = j < radius ? 0 : 2;
            const float weight_distance = float(j < radius ? j + 1 : 2 * radius - j);
            const float tw = mode == 0 ? 1.f : mode == 1 ? std::sqrt(weight_distance / float(2 * radius)) :
                mode == 2 ? std::sin((weight_distance + 1) / float(2 * radius)) :
                1.f / float(j < radius ? radius - j + 1 : j - radius + 2);
            const uint64_t dy=uint64_t(std::abs(int64_t(value(1,0))-int64_t(value(f,0))));
            const uint64_t dc=uint64_t(std::abs(int64_t(value(1,p))-int64_t(value(f,p))));
            uint64_t weight=(uint64_t(table[dy>>(bits-8)])*table[dc>>(bits-8)])<<((bits-8)*2);
            weight=uint64_t(std::round(float(weight)*tw));
            const uint64_t result=(weight*value(f,p)+(max-weight)*value(1,p)+max/2)>>(bits*2);
            sum+=(max-dy-dc)*result;
          }
          if(reinterpret_cast<T*>(dst[p])[i]!=(sum+max*radius)/(max*2*radius)) throw std::runtime_error("Cnr4 exact integer oracle mismatch");
        }
      }
    }
  }

  if(!half) {
    const float weights2[]={0.5f,0.25f};
    const uint8_t* prev[]={src[0]},*next[]={src[6]};
    neo_smo::process_ttempsmooth_plane(type,width,height,pitch,pitch,pitch,src[3],src[3],prev,prev,next,next,dst[0],1,1,1,0.01f,true,bits,1,0.5f,weights2,nullptr);
  }
}
// Exercise each packed comparison group, independent reference clips, and
// unequal pitches. Binary-exact weights isolate final division rounding.
template <class T>
void check_ttempsmooth_packed(int bits, int width, bool end) {
  constexpr int height = 2;
  const int scale = 1 << (bits - 8);
  const int source_pitch = width + 3, reference_pitch = width + 5, output_pitch = width + 7;
  std::vector<std::unique_ptr<Guarded>> storage;
  std::array<const uint8_t*, 7> source{}, reference{};
  for (int frame = 0; frame < 7; ++frame) {
    for (int is_ref = 0; is_ref < 2; ++is_ref) {
      const int pitch = is_ref ? reference_pitch : source_pitch;
      storage.push_back(std::make_unique<Guarded>((pitch + width) * sizeof(T), end));
      auto* data = reinterpret_cast<T*>(storage.back()->data);
      (is_ref ? reference : source)[frame] = storage.back()->data;
      for (int y = 0; y < height; ++y) for (int x = 0; x < width; ++x) {
        const int value = is_ref ? 100 + (x * 13 + x / 7 + frame * 3 +
            (frame % 2) * ((x >> 4) % 5) + y) % 9 : (x * 97 + frame * 33 + y * 9) % 256;
        data[y * pitch + x] = static_cast<T>(value * scale + (is_ref ? 0 : frame % scale));
      }
    }
  }
  storage.push_back(std::make_unique<Guarded>((output_pitch + width) * sizeof(T), end));
  uint8_t* output = storage.back()->data;
  const uint8_t* prev[] = {source[2], source[1], source[0]};
  const uint8_t* next[] = {source[4], source[5], source[6]};
  const uint8_t* prev_ref[] = {reference[2], reference[1], reference[0]};
  const uint8_t* next_ref[] = {reference[4], reference[5], reference[6]};
  const float weights[] = {0.25f, 0.125f, 0.0625f, 0.0625f};
  const auto value = [](const uint8_t* plane, int index) {
    return float(reinterpret_cast<const T*>(plane)[index]);
  };
  for (int threshold : {0, 1, 2, 256}) for (bool fp : {false, true}) {
    for (int num_prev : {0, 1, 3}) for (int num_next : {0, 3}) {
      neo_smo::process_ttempsmooth_plane(sizeof(T) == 1 ? neo_smo::DataType::U8 : neo_smo::DataType::U16,
          width, height, source_pitch * sizeof(T), reference_pitch * sizeof(T), output_pitch * sizeof(T),
          source[3], reference[3], prev, prev_ref, next, next_ref, output,
          3, num_prev, num_next, float(threshold * scale), fp, bits, 1, weights[0], weights, nullptr);
      for (int y = 0; y < height; ++y) for (int x = 0; x < width; ++x) {
        const float current = value(source[3], y * source_pitch + x);
        const float current_ref = value(reference[3], y * reference_pitch + x);
        float sum = current * weights[0], total = weights[0];
        for (int dir = 0; dir < 2; ++dir) {
          const int frames = dir == 0 ? num_prev : num_next;
          if (frames == 0) break;
          float previous = current_ref;
          for (int i = 0; i < frames; ++i) {
            const auto* refs = dir == 0 ? prev_ref : next_ref;
            const auto* sources = dir == 0 ? prev : next;
            const float neighbor = value(refs[i], y * reference_pitch + x);
            if (std::abs(current_ref - neighbor) >= threshold * scale ||
                (i != 0 && std::abs(previous - neighbor) >= threshold * scale)) break;
            total += weights[1 + i];
            sum += value(sources[i], y * source_pitch + x) * weights[1 + i];
            previous = neighbor;
          }
        }
        const double exact = fp ? double(current) * (1.0 - double(total)) + double(sum)
                                : double(sum) / double(total);
        const double expected = std::round(exact);
        const float got = value(output, y * output_pitch + x);
        // Fast division may land just below an exact half-way result, including
        // on the existing narrow-width path. Only allow the other adjacent
        // integer at that exact midpoint; all other values must match exactly.
        const bool lower_tie = exact - std::floor(exact) == 0.5 && got == std::floor(exact);
        if (got != expected && !lower_tie) {
          std::fprintf(stderr, "TTempSmooth packed bits=%d width=%d xy=%d,%d threshold=%d fp=%d prev=%d next=%d got=%g expected=%g\n",
              bits, width, x, y, threshold, fp, num_prev, num_next, got, expected);
          throw std::runtime_error("TTempSmooth packed scalar oracle mismatch");
        }
      }
    }
  }
}

void check_ccd_packed(int width, bool end) {
  constexpr int height = 5;
  const std::size_t pitch = width * sizeof(uint16_t);
  std::vector<std::unique_ptr<Guarded>> storage;
  std::array<const uint8_t*, 3> src{}, ref{};
  std::array<uint8_t*, 3> dst{};
  for (int p = 0; p < 3; ++p) {
    for (int kind = 0; kind < 3; ++kind) {
      storage.push_back(std::make_unique<Guarded>(height * pitch, end));
      auto* data = storage.back()->data;
      if (kind == 0) src[p] = data;
      else if (kind == 1) ref[p] = data;
      else dst[p] = data;
    }
  }
  std::vector<neo_smo::Point> points(127);
  for (int i = 0; i < 127; ++i) points[i] = {i % 5 - 2, i % 3 - 1};
  const float weight = 1;
  for (int bits : {10, 16}) {
    const int peak = (1 << bits) - 1;
    for (int p = 0; p < 3; ++p) for (int i = 0; i < width * height; ++i) {
      reinterpret_cast<uint16_t*>(const_cast<uint8_t*>(src[p]))[i] = uint16_t((i * 17113 + p * 113) & peak);
      reinterpret_cast<uint16_t*>(const_cast<uint8_t*>(ref[p]))[i] = uint16_t(i % 3 == 0 ? peak : (i % 3 == 1 ? 0 : (i * 7919 + p * 211) & peak));
    }
    auto sample = [&](const auto& planes, int p, int x, int y) -> int64_t {
      return reinterpret_cast<const uint16_t*>(planes[p])[reflect(y, height) * width + reflect(x, width)];
    };
    for (bool rgb : {false, true}) for (int count : {4, 127})
      for (float threshold : {0.f, 1.f, 123456.f, 4294836224.f, 4294967296.f, 1.e12f}) {
        neo_smo::process_ccd_planes(neo_smo::DataType::U16, rgb, width, height, pitch,
            src.data(), ref.data(), dst[0], dst[1], dst[2], threshold, 0, &weight,
            points.data(), count, 5, 1.f, bits);
        for (int y = 0; y < height; ++y) for (int x = 0; x < width; ++x) {
          int64_t sum[3] = {sample(src, 0, x, y), sample(src, 1, x, y), sample(src, 2, x, y)};
          int accepted = 1;
          for (int i = 0; i < count; ++i) {
            const int px = x + points[i].x, py = y + points[i].y;
            int64_t distance = 0;
            for (int p = 0; p < 3; ++p) {
              const auto diff = sample(ref, p, px, py) - sample(ref, p, x, y);
              distance += diff * diff * (!rgb && p == 0 ? 4 : 1);
            }
            if (double(distance) < std::floor(double(threshold))) {
              for (int p = 0; p < 3; ++p) sum[p] += sample(src, p, px, py);
              ++accepted;
            }
          }
          for (int p = rgb ? 0 : 1; p < 3; ++p) {
            const auto expected = (2 * sum[p] + accepted) / (2 * accepted);
            if (reinterpret_cast<const uint16_t*>(dst[p])[y * width + x] != expected)
              throw std::runtime_error("CCD packed exact integer oracle mismatch");
          }
        }
      }
  }
}


int main() {
  for(auto target:hwy::SupportedAndGeneratedTargets()) {
    hwy::SetSupportedTargetsForTest(target);
    for(bool end:{false,true}) for(int width:{1,2,7,17,33,65}) for(int height:{1,5}) {
      check<uint8_t>(neo_smo::DataType::U8,width,height,end);
      check<uint16_t>(neo_smo::DataType::U16,width,height,end);
      check<uint16_t>(neo_smo::DataType::F16,width,height,end);
      check<float>(neo_smo::DataType::F32,width,height,end);
    }
    for (bool end : {false, true}) for (int width : {17, 33, 65}) check_ccd_packed(width, end);
    for (bool end : {false, true}) for (int width : {1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127}) {
      check_ttempsmooth_packed<uint8_t>(8, width, end);
      check_ttempsmooth_packed<uint16_t>(10, width, end);
      check_ttempsmooth_packed<uint16_t>(16, width, end);
    }
    std::printf("%s: Phase5 oracles and guarded tails passed\n",hwy::TargetName(target));
  }
  hwy::SetSupportedTargetsForTest(0);
}
