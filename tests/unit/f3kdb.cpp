#include "f3kdb_scalar_filter.hpp"
#include "algorithms/f3kdb.hpp"
#include "kernels/f3kdb.hpp"
#ifndef HWY_COMPILE_ALL_ATTAINABLE
#define HWY_COMPILE_ALL_ATTAINABLE
#endif
#include "hwy/targets.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <future>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
namespace ref=neo_smo::f3kdb;
namespace simd=neo_smo::f3kdb;
namespace {
std::size_t checks=0,pixels=0,different=0;
std::size_t fs_pixels=0;
double max_weight=0,max_integer=0,max_float=0,max_pixel=0;
void check(bool ok,const std::string& what) { ++checks; if (!ok) throw std::runtime_error(what); }
template<class T> T read(const void* p) { T v; std::memcpy(&v,p,sizeof(v)); return v; }
template<class T> void write(void* p,T v) { std::memcpy(p,&v,sizeof(v)); }
void fs_tests() {
  // Independent row-major diffusion: wavefront scheduling, vector lane count
  // and stripe height must not change either output codes or outgoing residuals.
  std::mt19937 rng(0x46535746);
  // Pair shapes instead of multiplying every dimension and format. Keep tiny
  // rows, partial vectors, multiple wavefront groups and the 32-row boundary.
  for (const auto shape : {std::array<std::size_t,2>{1,1}, {2,2}, {3,7}, {7,3}, {17,9}, {3,33}})
    for (int variant=0; variant<8; ++variant) {
        const auto width=shape[0], height=shape[1];
        const int depths[]={8,9,10,12,16};
        const int depth=depths[variant%5];
        const bool integer=(variant&1)!=0, chroma=(variant&2)!=0, limited=(variant&4)!=0;
        const std::size_t bytes=depth==8?1:2, stride=width*bytes+7;
        const double scale=std::ldexp(1.0,depth), tvscale=std::ldexp(1.0,depth-8);
        const double lo=limited?16*tvscale:0, hi=limited?(chroma?240:235)*tvscale:scale-1;
        std::vector<double> values(width*height);
        for (std::size_t i=0;i<values.size();++i) {
          double code=(rng()%static_cast<unsigned>(scale))+double(rng()%1024)/1024;
          if (i%11==0) code=17.5;
          if (i%11==1) code=18.5;
          if (i%11==2) code=std::nextafter(17.5,0.0);
          if (i%11==3) code=std::nextafter(18.5,32.0);
          if (i%11==4) code=-1e30;
          if (i%11==5) code=1e30;
          if (i%11==6) code=std::numeric_limits<float>::denorm_min();
          values[i]=integer?std::ldexp(code,16-depth):code/scale-(chroma?0.5:0);
        }
        std::vector<float> initial(width),previous(width),current(width);
        for (auto& e:initial) e=(static_cast<int>(rng()%1024)-512)/1024.0f;
        previous=initial;
        std::vector<std::uint8_t> expected(height*stride+2,0xa5);
        for (std::size_t y=0;y<height;++y) {
          for (std::size_t x=0;x<width;++x) {
            const double z=values[y*width+x];
            const double value=integer?std::ldexp(z,depth-16):(z+(chroma?0.5:0))*scale;
            float error=x?current[x-1]*(7.0f/16):0;
            if (x) error+=previous[x-1]*(1.0f/16);
            error+=previous[x]*(5.0f/16);
            if (x+1<width) error+=previous[x+1]*(3.0f/16);
            const double v=std::clamp(double(float(std::clamp(value,lo,hi)+error)),lo,hi);
            const double q=ref::round_even(v);
            current[x]=static_cast<float>(v-q);
            auto* ptr=expected.data()+1+y*stride+x*bytes;
            if (bytes==1) *ptr=static_cast<std::uint8_t>(q);
            else write(ptr,static_cast<std::uint16_t>(q));
          }
          previous.swap(current);
        }
        for (std::size_t stripe_height:{1,3,32}) {
          std::vector<std::uint8_t> actual(expected.size(),0xa5);
          auto boundary=initial;
          for (std::size_t y=0;y<height;y+=stripe_height) {
            simd::diffuse_stripe({values.data()+y*width,actual.data()+1+y*stride,boundary.data(),
                                 width,std::min(stripe_height,height-y),stride,depth,integer,chroma,limited});
          }
          fs_pixels+=values.size();
          check(actual==expected,"FS wavefront pixels/padding w="+std::to_string(width)+" h="+std::to_string(height)+
                                " depth="+std::to_string(depth)+" stripe="+std::to_string(stripe_height));
          check(boundary==previous,"FS outgoing signed residuals");
        }
      }
}
void weight_tests() {
  constexpr std::size_t n=65;
  std::array<std::vector<double>,4> c;
  for (auto& v:c) v.resize(n);
  std::vector<double> actual(n);
  // Simultaneously cover endpoints, values whose product underflows, and a
  // monotone sweep. Checking the weight directly avoids hiding error in pixels.
  for (int pattern=0;pattern<5;++pattern) {
    for (std::size_t i=0;i<n;++i) for (int j=0;j<4;++j)
      c[j][i]=pattern==0?double(i)/(n-1):pattern==1?(j==0?double(i)/(n-1):1.0):
              pattern==2?std::exp(-740.0+740.0*i/(n-1)): (i%2?1.0:0.0);
    if (pattern==4) {
      // Keep the original Highway Log regression explicitly when reducing sweeps.
      const double edges[]={0, std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::min(), 1.0/4098, 0x1p-12,
        std::nextafter(0x1p-12, 0.0), 0.5, std::nextafter(1.0, 0.0), 1};
      for (std::size_t i=0;i<n;++i) for (int j=0;j<4;++j) c[j][i]=edges[i%9];
    }
    const double* ptr[4]={c[0].data(),c[1].data(),c[2].data(),c[3].data()};
    simd::weights(ptr,actual.data(),n);
    for (std::size_t i=0;i<n;++i) {
      const double expected=ref::weight({c[0][i],c[1][i],c[2][i],c[3][i]});
      max_weight=std::max(max_weight,std::abs(actual[i]-expected));
      check(std::isfinite(actual[i]) && actual[i]>=0 && actual[i]<=1,"weight range");
      if (!(std::abs(actual[i]-expected)<=0x1p-20)) std::cerr<<"weight pattern="<<pattern<<" i="<<i<<" actual="<<actual[i]<<" expected="<<expected<<"\n";
      check(std::abs(actual[i]-expected)<=0x1p-20,"N05 weight accuracy");
      if (expected==0 || expected==1) check(actual[i]==expected,"weight endpoint");
      if (pattern<3 && i) check(actual[i]>=actual[i-1],"weight monotonicity pattern="+std::to_string(pattern)+" i="+std::to_string(i)+" previous="+std::to_string(actual[i-1])+" current="+std::to_string(actual[i]));
    }
  }
}
void row_tests() {
  constexpr std::size_t n=17;
  std::mt19937 rng(0x3f3db);
  std::array<std::vector<double>,16> storage;
  for (auto& v:storage) v.resize(n);
  std::vector<double> out(n),grain(n);
  simd::Row row{}; row.center=storage[0].data(); row.count=n; row.output=out.data(); row.grain=grain.data();
  for (int j=0;j<4;++j) row.refs[j]=storage[j+1].data();
  for (int j=0;j<5;++j) { row.gx[j]=storage[5+2*j].data(); row.gy[j]=storage[6+2*j].data(); }
  for (bool integer:{true,false}) for (int pattern=0;pattern<5;++pattern) {
    for (std::size_t i=0;i<n;++i) {
      for (int j=0;j<5;++j) {
        double v=32768+int(rng()%4097)-2048;
        if (pattern==1) v=32768+(int(i%7)-3)*64; // identical references and constant preservation
        if (!integer) {
          v=static_cast<float>((v-32768)/65536.0);
          if (pattern==2) v=std::numeric_limits<float>::denorm_min()*(1+int(rng()%4));
          if (pattern==3) v=(j%2?-1:1)*double(std::numeric_limits<float>::max());
          if (pattern==4) v=static_cast<float>(v*1e-25);
        }
        storage[j][i]=v;
        storage[5+2*j][i]=pattern==1?0:int(rng()%524281)-262140;
        storage[6+2*j][i]=pattern==1?0:int(rng()%524281)-262140;
        if (!integer && pattern>=2) { storage[5+2*j][i]*=1e-35; storage[6+2*j][i]*=1e-35; }
      }
      grain[i]=0;
    }
    for (int mode=1;mode<=7;++mode) for (bool blur:{true,false}) for (double angle:{0.0,0.15,0.25,0.5}) {
      for (double threshold:{0.0,1500.0,65535.0}) {
        const double t=integer?threshold:threshold/65536;
        const std::array<double,3> thresholds={t,t*0.75,t*0.5};
        simd::process_row(row,mode,thresholds.data(),integer,blur,angle,1.5);
        for (std::size_t i=0;i<n;++i) {
          bool aligned=true;
          for (int j=1;j<5;++j) aligned=aligned && ref::aligned({row.gx[0][i],row.gy[0][i]},{row.gx[j][i],row.gy[j][i]},angle);
          const double expected=ref::deband(row.center[i],{row.refs[0][i],row.refs[1][i],row.refs[2][i],row.refs[3][i]},mode,thresholds,integer,blur,aligned,1.5);
          const double error=std::abs(out[i]-expected);
          check(std::isfinite(out[i]),"finite row output");
          if (integer) { max_integer=std::max(max_integer,error); check(error<=(mode<=5?0:1),"integer row mode="+std::to_string(mode)); }
          else {
            max_float=std::max(max_float,error);
            const double scale=std::max(std::abs(expected),std::abs(row.center[i]));
            check(error<=8*std::numeric_limits<float>::epsilon()*scale+8*std::numeric_limits<float>::denorm_min(),
                  "float row mode="+std::to_string(mode)+" pattern="+std::to_string(pattern)+" error="+std::to_string(error));
          }
        }
      }
    }
  }
}
struct Input {
  std::vector<std::vector<std::uint8_t>> bytes;
  std::vector<ref::PlaneView> views;
};
Input input(const ref::Format& f,int pattern) {
  Input in;
  const int planes=f.color==ref::Color::gray?1:3;
  const std::size_t b=f.depth<=8?1:f.depth<=16?2:4;
  in.bytes.resize(planes);
  std::mt19937 rng(0x7f3db);
  for (int p=0;p<planes;++p) {
    const auto w=f.width>>(p?f.sub_x:0),h=f.height>>(p?f.sub_y:0),stride=w*b+7;
    in.bytes[p].resize(h*stride+1,0xa5);
    auto* data=in.bytes[p].data()+1; // deliberately unaligned, including rows
    for (std::size_t y=0;y<h;++y) for (std::size_t x=0;x<w;++x) {
      double v;
      if (pattern==0) v=(32768+int((x*17+y*29)%257)-128)/65536.0;
      else if (pattern==1) v=(rng()%65536)/65536.0;
      else if (pattern==2) v=0.25;
      else v=(x%3==0?0.0:x%3==1?1.0:0.5);
      auto* ptr=data+y*stride+x*b;
      if (f.sample==ref::Sample::integer) {
        const auto code=std::uint16_t(std::min(std::ldexp(1.0,f.depth)-1,std::floor(std::ldexp(v,f.depth))));
        if (b==1) *ptr=static_cast<std::uint8_t>(code); else write(ptr,code);
      } else {
        if (p) v-=0.5;
        if (b==2) write(ptr,ref::encode_half(v)); else write(ptr,static_cast<float>(v));
      }
    }
    in.views.push_back({data,h*stride,stride});
  }
  return in;
}
double sample(const ref::Frame& f,std::size_t p,std::size_t i) {
  const auto b=f.format.depth<=8?1:f.format.depth<=16?2:4;
  const auto* ptr=f.planes[p].data.data()+i*b;
  return f.format.sample==ref::Sample::integer?(b==1?*ptr:read<std::uint16_t>(ptr)):
         b==2?ref::decode_half(read<std::uint16_t>(ptr)):read<float>(ptr);
}
void compare(const ref::Frame& a,const ref::Frame& b,const ref::Options& o) {
  check(a.planes.size()==b.planes.size(),"output planes");
  for (std::size_t p=0;p<a.planes.size();++p) {
    check(a.planes[p].data.size()==b.planes[p].data.size(),"output storage");
    for (std::size_t i=0;i<a.planes[p].width*a.planes[p].height;++i) {
      const auto x=sample(a,p,i),y=sample(b,p,i),e=std::abs(x-y);
      ++pixels; different+=e!=0; max_pixel=std::max(max_pixel,e);
      // Exact discrete modes; blends allow one internal code. Quantization can
      // turn a tiny continuous error into one output code. FS residual feedback
      // adds at most one incoming code plus the final rounding boundary.
      const double limit=a.format.sample==ref::Sample::integer?(o.sample_mode<=5?0:o.dither_algo==3?3:1):
                         a.format.depth==16?0x1p-10:8*std::numeric_limits<float>::epsilon()*std::max(1.0,std::abs(x));
      check(std::isfinite(y) && e<=limit,"frame pixels mode="+std::to_string(o.sample_mode)+" dither="+std::to_string(o.dither_algo)+" delta="+std::to_string(e));
    }
  }
}
void frame_tests() {
  // Integer-output FS must be exercised explicitly; floating output bypasses
  // dithering even when dither_algo=3. These frames also cross the 32-row buffer.
  int case_index=0;
  for (int depth:{8,10,16,32,-16}) for (bool limited:{false,true}) {
      const int mode=case_index%7+1;
      const int outputs[]={8,10,16};
      const int output_depth=outputs[case_index++%3];
      ref::Format f;
      f.width=18; f.height=34; f.depth=std::abs(depth);
      f.sample=depth==32 || depth<0?ref::Sample::floating:ref::Sample::integer;
      f.color=ref::Color::yuv; f.sub_x=1; f.sub_y=1;
      ref::Options o; o.sample_mode=mode; o.dither_algo=3; o.keep_tv_range=limited;
      o.grain_y=3.5; o.grain_c=2.25; o.output_sample_type=ref::OutputSample::integer; o.output_depth=output_depth;
      const auto in=input(f,1);
      compare(ScalarFilter(f,o).process(0,in.views),simd::Filter(f,o).process(0,in.views),o);
    }
  // Kernel tests cover each operation. Integration uses rotating combinations
  // of formats, modes and options rather than their Cartesian product.
  for (int depth:{8,10,12,16,32,-16}) for (int mode=1;mode<=7;++mode) {
    const int shape=case_index%4, variant=(case_index/4+mode)%4;
    ++case_index;
    ref::Format f;
    f.depth=std::abs(depth); f.sample=(depth==32 || depth<0)?ref::Sample::floating:ref::Sample::integer;
    f.width=shape==0?17:shape==1?2:18; f.height=shape==1?2:shape==0?3:6; f.frames=3;
    if (shape) { f.color=ref::Color::yuv; f.sub_x=shape==2?1:0; f.sub_y=shape==2?1:shape==3?1:0; }
    {
      ref::Options o; o.sample_mode=mode; o.dither_algo=variant%3+1; o.dynamic_grain=variant%2!=0;
      o.grain_y=variant==0?0:5.25; o.grain_c=variant==0?0:3.5; o.keep_tv_range=variant==3;
      o.range=variant==2?255:15; o.blur_first=variant!=1; o.max_angle=variant==0?0:variant==1?0.25:variant==2?0.5:0.15;
      o.random_algo_ref=variant%3; o.random_algo_grain=(variant+1)%3;
      if (variant==1) { o.output_sample_type=ref::OutputSample::integer; o.output_depth=10; }
      if (variant==2) { o.output_sample_type=ref::OutputSample::floating; o.output_depth=32; }
      if (variant==3) { o.output_sample_type=ref::OutputSample::floating; o.output_depth=16; }
      const auto in=input(f,variant);
      const ScalarFilter oracle(f,o); const simd::Filter filter(f,o);
      compare(oracle.process(2,in.views),filter.process(2,in.views),o);
    }
  }
  // Range 255 must actually reach beyond signed-byte distances.
  ref::Options o; o.sample_mode=7; o.range=255; o.grain_y=0;
  for (bool vertical:{false,true}) {
    ref::Format large; large.width=vertical?3:519; large.height=vertical?519:3;
    const auto in=input(large,0);
    compare(ScalarFilter(large,o).process(0,in.views),simd::Filter(large,o).process(0,in.views),o);
  }
  ref::Format f; f.width=17; f.height=3; f.frames=3; f.sample=ref::Sample::floating; f.depth=32;
  auto small=input(f,1); o.dynamic_grain=true; const simd::Filter filter(f,o);
  auto a=std::async(std::launch::async,[&]{ return filter.process(2,small.views); });
  auto b=std::async(std::launch::async,[&]{ return filter.process(0,small.views); });
  check(a.get().planes[0].data==filter.process(2,small.views).planes[0].data,"concurrent repeat");
  check(b.get().planes[0].data==filter.process(0,small.views).planes[0].data,"concurrent order");
  for (float bad:{std::numeric_limits<float>::infinity(),std::numeric_limits<float>::quiet_NaN()}) {
    write(small.bytes[0].data()+1,bad);
    bool rejected=false; try { filter.process(0,small.views); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"nonfinite validation under production flags");
  }
  bool rejected=false; try { filter.process(3,small.views); } catch(const std::invalid_argument&) { rejected=true; }
  check(rejected,"frame index validation");
  small.views[0].stride=1; rejected=false;
  try { filter.process(0,small.views); } catch(const std::invalid_argument&) { rejected=true; }
  check(rejected,"stride validation");
}
} // namespace
int main() {
  try {
    const auto targets = hwy::SupportedAndGeneratedTargets();
    for (auto requested : {HWY_SSE4, HWY_AVX2, HWY_AVX3})
      if (std::find(targets.begin(), targets.end(), requested) == targets.end())
        std::cout << "Skipping unavailable target " << hwy::TargetName(requested) << '\n';
    for (auto target : targets) {
      // Normal x86 builds retain the reduced SIMD selection. Dedicated
      // builds exercise SCALAR and EMU128 without duplicating them here.
#if HWY_ARCH_X86 && !defined(HWY_COMPILE_ONLY_SCALAR) && !defined(HWY_COMPILE_ONLY_EMU128)
      if (target != HWY_SSE4 && target != HWY_AVX2 && target != HWY_AVX3)
        continue;
#endif
      hwy::SetSupportedTargetsForTest(target);
      std::cout<<"Testing "<<hwy::TargetName(target)<<std::endl;
      weight_tests(); row_tests(); fs_tests(); frame_tests();
      std::cout<<hwy::TargetName(target)<<": checks="<<checks<<" pixels="<<pixels<<" differing="<<different
               <<" fs_pixels="<<fs_pixels<<" max_weight="<<max_weight<<" max_integer_internal="<<max_integer<<" max_float_internal="<<max_float<<" max_pixel="<<max_pixel<<'\n';
    }
    if (!checks) throw std::runtime_error("No supported F3kdb test target executed");
    hwy::SetSupportedTargetsForTest(0);
  } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
