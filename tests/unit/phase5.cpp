#include "kernels/dispatch.hpp"
#include "base/fp16.hpp"
#include "guarded.hpp"
#include "hwy/targets.h"
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
  const size_t size = width * height, pitch = width * sizeof(T);
  std::vector<std::unique_ptr<Guarded>> storage;
  std::array<const uint8_t*, 9> src{};
  std::array<uint8_t*, 3> dst{};
  auto decode = [&](T v) { return half ? neo_smo::fp16_to_fp32(static_cast<uint16_t>(v)) : double(v); };
  auto encode = [&](double v) { return half ? T(neo_smo::fp32_to_fp16(float(v))) : T(v); };
  auto round = [&](double v) { return half ? double(neo_smo::fp16_to_fp32(neo_smo::fp32_to_fp16(float(v)))) : integer ? v : double(float(v)); };
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
  const neo_smo::Point points[] = {{-2,-1},{1,-1},{0,1},{2,1}};
  const float weights[] = {0.70710677f, 1.f, 0.68163878f};
  for (bool rgb : {false,true}) for (int radius : {0,1}) {
    const float threshold = integer ? 123456.f : 0.0003f;
    const auto* sources = radius ? src.data() : src.data()+3;
    neo_smo::process_ccd_planes(type,rgb,width,height,pitch,sources,sources,dst[0],dst[1],dst[2],threshold,radius,
                               radius?weights:weights+1,points,4,5,1.f,bits);
    auto sample = [&](int f, int p, int x, int y) {return decode(reinterpret_cast<const T*>(src[f*3+p])[reflect(y,height)*width+reflect(x,width)]);};
    for(int y=0;y<height;++y) for(int x=0;x<width;++x) {
      double totals[3] = {sample(1,0,x,y),sample(1,1,x,y),sample(1,2,x,y)};
      int accepted=1;
      for(auto point:points) {
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
        if(std::abs(got-decode(encode(expected)))>(integer||half?0:1e-6)) {
          std::fprintf(stderr,"CCD type=%d rgb=%d radius=%d width=%d xy=%d,%d got=%g expected=%g\n",int(type),rgb,radius,width,x,y,got,expected);
          throw std::runtime_error("CCD scalar oracle mismatch");
        }
      }
    }
  }
  if(integer) {
    std::array<uint8_t,256> table{}; for(int i=0;i<256;++i) table[i]=uint8_t(255-i);
    std::array<const uint8_t*,3> neighbors[2] = {{{src[0],src[1],src[2]}},{{src[6],src[7],src[8]}}};
    for(int mode=0;mode<4;++mode) {
      neo_smo::process_cnr4_frame(type,width,height,pitch,bits,1,0,mode,src.data()+3,src.data()+3,neighbors,neighbors,2,dst[1],dst[2],table.data(),table.data(),table.data());
      const float tw=mode==0?1.f:mode==1?std::sqrt(0.5f):mode==2?std::sin(1.f):0.5f;
      const uint64_t max=uint64_t{1}<<(bits*2);
      for(size_t i=0;i<size;++i) for(int p=1;p<3;++p) {
        uint64_t sum=0;
        auto value=[&](int f,int plane) {return uint64_t(reinterpret_cast<const T*>(src[f*3+plane])[i]);};
        for(int f:{0,2}) {
          const uint64_t dy=uint64_t(std::abs(int64_t(value(1,0))-int64_t(value(f,0))));
          const uint64_t dc=uint64_t(std::abs(int64_t(value(1,p))-int64_t(value(f,p))));
          uint64_t weight=(uint64_t(table[dy>>(bits-8)])*table[dc>>(bits-8)])<<((bits-8)*2);
          weight=uint64_t(std::round(float(weight)*tw));
          const uint64_t result=(weight*value(f,p)+(max-weight)*value(1,p)+max/2)>>(bits*2);
          sum+=(max-dy-dc)*result;
        }
        if(reinterpret_cast<T*>(dst[p])[i]!=(sum+max)/(max*2)) throw std::runtime_error("Cnr4 exact integer oracle mismatch");
      }
    }
  }
  if(!half) {
    const float weights2[]={0.5f,0.25f};
    const uint8_t* prev[]={src[0]},*next[]={src[6]};
    neo_smo::process_ttempsmooth_plane(type,width,height,pitch,pitch,pitch,src[3],src[3],prev,prev,next,next,dst[0],1,1,1,0.01f,true,bits,1,0.5f,weights2,nullptr);
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
    std::printf("%s: Phase5 oracles and guarded tails passed\n",hwy::TargetName(target));
  }
  hwy::SetSupportedTargetsForTest(0);
}
