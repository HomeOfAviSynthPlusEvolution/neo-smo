#include "algorithms/f3kdb.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <future>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace f=neo_smo::f3kdb;
namespace {
int checks=0;
void check(bool ok,const char* description) { ++checks; if (!ok) throw std::runtime_error(description); }
void near(double a,double b,double error,const char* description) { check(std::abs(a-b)<=error,description); }
template<class F> void rejects(F&& action,const char* description) {
  bool rejected=false;
  try { action(); } catch (const std::exception&) { rejected=true; }
  check(rejected,description);
}
f::Options identity() { f::Options o; o.presets={"depth"}; o.dither_algo=1; return o; }
template<class T> f::PlaneView view(const std::vector<T>& values,std::size_t stride) { return {values.data(),values.size()*sizeof(T),stride*sizeof(T)}; }
template<class T> std::vector<T> values(const f::Plane& plane) {
  std::vector<T> v(plane.width*plane.height);
  for(std::size_t y=0;y<plane.height;++y) std::memcpy(v.data()+y*plane.width,plane.data.data()+y*plane.stride,plane.width*sizeof(T));
  return v;
}
void parameters() {
  f::Format fmt;
  for (bool scale:{false,true}) {
    auto o=identity(); o.scale=scale; o.presets={"/high//luma/"};
    auto p=f::normalize(fmt,o);
    check(p.thresholds[0]==std::array<double,3>{256,256,256},"preset scale invariance");
    check(p.thresholds[1]==std::array<double,3>{0,0,0} && p.grain[1]==0,"luma preset");
    check(p.grain[0]==256,"preset grain units");
    o.threshold[1]=1.25; o.threshold1[1]=-1;
    p=f::normalize(fmt,o);
    check(p.thresholds[1][0]==1.25*(scale?1:4) && p.thresholds[1][1]==p.thresholds[1][0],"explicit override then inheritance");
  }
  auto p=f::normalize(fmt,{}); check(p.thresholds[0]==std::array<double,3>{256,256,256},"default inheritance");
  auto o=identity(); o.output_sample_type=f::OutputSample::floating;
  check(f::normalize(fmt,o).output.depth==32,"implicit float depth");
  fmt.sample=f::Sample::floating; fmt.depth=16; o.output_sample_type=f::OutputSample::integer;
  check(f::normalize(fmt,o).output.depth==16,"implicit integer depth");
  o.output_depth=-2; rejects([&]{f::normalize(fmt,o);},"negative output depth");
  o=identity(); o.threshold1[2]=-2; rejects([&]{f::normalize(fmt,o);},"absent chroma still validated");
  o=identity(); o.random_param_ref=-1; rejects([&]{f::normalize(fmt,o);},"inactive sigma validated");
  o=identity(); o.max_angle=1.01; rejects([&]{f::normalize(fmt,o);},"angle upper bound");
  o=identity(); o.presets={"unknown"}; rejects([&]{f::normalize(fmt,o);},"unknown preset");
  for(double bad:{-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    o=identity(); o.grain_y=bad; rejects([&]{f::normalize(fmt,o);},"invalid grain");
    o=identity(); o.threshold[0]=bad; rejects([&]{f::normalize(fmt,o);},"invalid threshold");
  }
  o=identity(); o.dither_algo=4; rejects([&]{f::normalize(fmt,o);},"invalid float dither");
  fmt={}; fmt.width=std::numeric_limits<std::size_t>::max(); fmt.height=2;
  rejects([&]{f::Filter filter(fmt);},"size multiplication overflow before allocation");
  fmt.height=1; rejects([&]{f::Filter filter(fmt);},"sequence alignment overflow before allocation");
  fmt={}; fmt.color=f::Color::yuv; fmt.sub_x=1;
  rejects([&]{f::Filter filter(fmt);},"nonintegral chroma dimensions");
  fmt={}; fmt.color=static_cast<f::Color>(2); rejects([&]{f::Filter filter(fmt);},"RGB rejected");
}
void mathematics() {
  const std::array<double,3> high{100,100,100};
  std::array<double,4> permutation{0,0,0,1};
  do {
    check(f::deband(0,permutation,2,high,true,true)==0,"A4 permutation integer");
    check(f::deband(0,permutation,2,high,false,true)==0.25,"A4 permutation floating");
  } while(std::next_permutation(permutation.begin(),permutation.end()));
  check(f::deband(0,{0,0,1,1},2,high,true,true)==1,"A4 halfway");
  check(f::deband(0,{-0.25,0.5,0,0},1,high,false,true)==0.125,"signed float average");
  const double tiny=std::numeric_limits<float>::denorm_min();
  for(int mode=1;mode<=7;++mode) check(f::deband(tiny,{tiny,tiny,tiny,tiny},mode,high,false,true)==tiny,"subnormal constant preserved");
  for(int mode:{1,2,3,4}) {
    check(f::deband(0,{4,4,4,4},mode,{4,100,100},true,true)==0,"strict threshold equal rejects");
    check(f::deband(0,{4,4,4,4},mode,{4.01,100,100},true,true)==4,"fractional integer threshold");
  }
  check(f::deband(10,{8,8,20,20},4,{3,100,100},true,true)==9,"independent mode4 directions");
  check(f::deband(10,{0,20,10,10},1,{1,100,100},true,true)==10,"blur average acceptance");
  check(f::deband(11,{0,20,10,10},1,{2,100,100},true,false)==11,"individual reference rejection");
  check(f::deband(10,{8,8,8,8},5,{3,3,5},true,false)==8,"mode5 all conditions");
  for(auto t:{std::array<double,3>{2,3,5},std::array<double,3>{3,2,5},std::array<double,3>{3,3,4}})
    check(f::deband(10,{8,8,8,8},5,t,true,true)==10,"mode5 strict gates");
  check(f::response(0,0)==1 && f::response(1,0)==0,"zero threshold responses");
  check(f::response(1e-20,2e-20)==1,"tiny positive threshold");
  near(f::weight({1e-12,1e-12,1e-12,1e-12}),1.584893192461114e-5,1e-18,"underflow-safe weight");
  double last=0;
  for(int i=0;i<=256;++i) {
    const double c=i/256.0, w=f::weight({c,1,1,1});
    check(w>=last && w>=0 && w<=1,"weight monotonicity and range");
    near(w,static_cast<double>(std::pow(static_cast<long double>(c),0.1L)),0x1p-20,"independent long double weight oracle");
    last=w;
  }
  // Subnormal/normal transitions and representative logarithmic scales.
  for(int exp:{-1074,-1073,-1023,-1022,-1021,-512,-149,-148,-127,-126,-125,-64,-32,-2,-1}) {
    const double c=std::ldexp(1.0,exp);
    near(f::weight({c,c,c,c}),static_cast<double>(std::pow(static_cast<long double>(c),0.4L)),0x1p-20,"weight exponent sweep");
  }
  check(f::weight({0,1,1,1})==0 && f::weight({1,1,1,1})==1,"weight endpoints");
  check(f::aligned({0,0},{0,0},0),"both flat");
  check(!f::aligned({0,0},{1,0},1),"one flat even at max angle");
  check(f::aligned({1,1},{-1,-1},0),"undirected parallel");
  check(f::aligned({1,0},{1,1},0.25) && !f::aligned({1,0},{0,1},0.25),"45 degree endpoint");
  check(f::aligned({1,0},{0,1},0.5),"90 degree endpoint");
  check(f::aligned({1e-40,0},{1e-40,1e-40},0.25),"tiny gradient direction");
  check(!f::aligned({1,0},{1,1},0.2499),"angle outside tolerance");
  check(f::deband(10,{8,8,8,8},6,{100,100,100},true,true)==8,"continuous full acceptance");
  check(f::deband(10,{8,8,8,8},6,{0,0,0},true,true)==10,"continuous zero rejection");
  check(f::deband(10,{8,8,8,8},7,{3,3,3},true,true,false,2)==10,"mode7 unaligned no boost");
  check(f::deband(10,{8,8,8,8},7,{3,3,3},true,true,true,2)==8,"mode7 aligned boost");
  check(f::round_even(0.5)==0 && f::round_even(1.5)==2 && f::round_even(2.5)==2 && f::round_even(-1.5)==-2,"nearest even");
  std::array<bool,256> seen{};
  for(int y=0;y<16;++y) for(int x=0;x<16;++x) { auto b=f::bayer16(y,x); check(!seen[b],"Bayer permutation"); seen[b]=true; }
  check(f::bayer16(0,1)==128 && f::bayer16(1,0)==192 && f::bayer16(15,0)==255,"Bayer orientation");
  // The audit covers every finite exponent; keep focused boundaries here.
  for(unsigned h:{0u,1u,0x3ffu,0x400u,0x3bffu,0x3c00u,0x7bffu,0x8000u,0x8001u,0x83ffu,0x8400u,0xfbffu})
    check(f::encode_half(f::decode_half(static_cast<std::uint16_t>(h)))==h,"half boundary roundtrips");
  check(f::encode_half(1+0x1p-11)==0x3c00 && f::encode_half(1+3*0x1p-11)==0x3c02,"half ties even");
  check(f::encode_half(1e20)==0x7bff && f::encode_half(-1e20)==0xfbff,"half finite saturation");
}
void random_tables() {
  f::Random u(0); const double first=u.next(1,0);
  check(u.state()==1013904223u,"uniform recurrence");
  const auto m=(std::uint64_t{1013904223u}<<20)+(1013904223u>>12);
  check(first==2*static_cast<double>(m)/4503599627370496.0-1,"uniform mapping");
  f::Random zero(0); check(zero.next(0,0)==-1 && zero.state()==0,"legacy zero fixed point");
  f::Random gaussian(123), oracle(123); oracle.next(1,0);oracle.next(1,0);
  check(gaussian.next(2,0)==0 && gaussian.state()==oracle.state(),"zero sigma consumes two states");
  for(double sigma:{0.1,0.999,1.0,1e300}) for(int i=0;i<16;++i) { const double z=gaussian.next(2,sigma); check(z>-1 && z<1,"Gaussian open interval"); }
  rejects([&]{gaussian.next(2,-1);},"negative sigma");
  f::Format fmt; fmt.width=8;fmt.height=6;fmt.frames=3;
  f::Options o; o.dynamic_grain=true; o.range=255; o.sample_mode=2;
  f::Filter filter(fmt,o), repeat(fmt,o);
  check(filter.grain_y()==repeat.grain_y() && filter.offsets()==repeat.offsets(),"deterministic tables");
  const std::size_t length=128*((fmt.width+255)/128)*fmt.height;
  check(filter.grain_y().size()==3*length,"dynamic sequence capacity");
  for(auto offset:filter.offsets()) check(offset%16==0 && offset<=2*length,"wide aligned offsets");
  for(std::size_t y=0;y<fmt.height;++y) for(std::size_t x=0;x<fmt.width;++x) {
    const auto d=filter.distances()[y*fmt.width+x]; const auto r=std::min({x,y,fmt.width-1-x,fmt.height-1-y});
    check(d.first<=r && d.second<=r,"distance geometric bounds");
  }
  // Independent state-consumption replay includes discarded samples and zero
  // amplitude calls; it checks the first generated grain value exactly.
  f::Random replay((0x92D68CA2u-o.seed)^(8u<<16)^6u^(3u<<16)^3u);
  for(std::size_t y=0;y<6;++y) for(std::size_t x=0;x<8;++x) {
    replay.next(1,0); const auto r=std::min({x,y,7-x,5-y}); if(r){replay.next(1,0);replay.next(1,0);}
  }
  check(filter.grain_y()[0]==std::round(replay.next(1,0)*256),"RNG consumption order");
}
void output_cases() {
  f::Format fmt;fmt.width=3; auto o=identity();o.output_sample_type=f::OutputSample::floating;
  std::vector<std::uint8_t> bytes{128,255,0};
  auto result=f::Filter(fmt,o).process(0,{view(bytes,3)});
  check(values<float>(result.planes[0])==std::vector<float>{0.5f,255.0f/256,0},"integer to float exact scale");
  fmt.color=f::Color::yuv; std::vector<std::uint8_t> chroma{128,129,127};
  result=f::Filter(fmt,o).process(0,{view(bytes,3),view(chroma,3),view(chroma,3)});
  check(values<float>(result.planes[1])==std::vector<float>{0,1.0f/256,-1.0f/256},"zero centered chroma");
  fmt={};fmt.sample=f::Sample::floating;fmt.depth=32;fmt.width=3;
  std::vector<float> floats{1,1.25f,-0.25f}; o=identity();
  result=f::Filter(fmt,o).process(0,{view(floats,3)}); check(values<float>(result.planes[0])==floats,"float superrange identity");
  o.keep_tv_range=true; result=f::Filter(fmt,o).process(0,{view(floats,3)});
  check(values<float>(result.planes[0])==std::vector<float>{235.0f/256,235.0f/256,16.0f/256},"float TV clip");
  o.keep_tv_range=false; o.output_sample_type=f::OutputSample::integer; o.output_depth=8;
  result=f::Filter(fmt,o).process(0,{view(floats,3)});check(values<std::uint8_t>(result.planes[0])==std::vector<std::uint8_t>{255,255,0},"float integer saturation");
  floats={0.5f/256,1.5f/256,2.5f/256}; result=f::Filter(fmt,o).process(0,{view(floats,3)});
  check(values<std::uint8_t>(result.planes[0])==std::vector<std::uint8_t>{0,2,2},"output round even");
  fmt.width=16;fmt.height=16;floats.assign(256,0.5f/256);o.dither_algo=2;
  result=f::Filter(fmt,o).process(0,{view(floats,16)});auto v=values<std::uint8_t>(result.planes[0]);
  check(std::count(v.begin(),v.end(),0)==128 && std::count(v.begin(),v.end(),1)==128,"centered Bayer balance");
  fmt.width=2;fmt.height=1;floats={0.75f/256,0.75f/256};o.dither_algo=3;
  result=f::Filter(fmt,o).process(0,{view(floats,2)});check(values<std::uint8_t>(result.planes[0])==std::vector<std::uint8_t>{1,1},"signed FS residual");
  // F32 -> U16 still requires dither.
  fmt.width=16;fmt.height=16;floats.assign(256,0.5f/65536);o.output_depth=16;o.dither_algo=2;
  result=f::Filter(fmt,o).process(0,{view(floats,16)});auto words=values<std::uint16_t>(result.planes[0]);
  check(std::count(words.begin(),words.end(),1)==128,"16-bit float quantization remains active");
  fmt={};bytes={0};o=identity();o.keep_tv_range=true;
  result=f::Filter(fmt,o).process(0,{view(bytes,1)});check(result.planes[0].data[0]==16,"zero strength does not bypass TV");
  fmt.sample=f::Sample::floating;fmt.depth=32;floats={std::numeric_limits<float>::infinity()};o=identity();
  rejects([&]{f::Filter(fmt,o).process(0,{view(floats,1)});},"identity validates infinity");
  floats[0]=std::numeric_limits<float>::quiet_NaN();rejects([&]{f::Filter(fmt,o).process(0,{view(floats,1)});},"identity validates NaN");
  floats[0]=std::numeric_limits<float>::max();o.output_depth=16;
  result=f::Filter(fmt,o).process(0,{view(floats,1)});check(values<std::uint16_t>(result.planes[0])[0]==0x7bff,"float output finite half saturation");
  // Finite extremes must not produce a nonfinite intermediate comparison or
  // output; 20-pixel gradient replication also traverses both signs here.
  fmt.width=3;fmt.height=3;floats={-std::numeric_limits<float>::max(),0,std::numeric_limits<float>::max(),0,0,0,std::numeric_limits<float>::max(),0,-std::numeric_limits<float>::max()};
  o={};o.grain_y=0;o.output_depth=32;
  for(int mode=1;mode<=7;++mode) { o.sample_mode=mode;result=f::Filter(fmt,o).process(0,{view(floats,3)});for(float pixel:values<float>(result.planes[0]))check(std::isfinite(pixel),"finite F32 extremes stay finite"); }
  // Independent full-image residual transport oracle. Carefully chosen dyadic
  // input values and coefficients make the first rows exactly representable;
  // later rows explicitly reproduce F32 storage and prescribed incoming order.
  fmt.width=37;fmt.height=9;floats.resize(fmt.width*fmt.height);
  for(std::size_t i=0;i<floats.size();++i) floats[i]=static_cast<float>((i*19)%97+0.375)/256;
  o=identity();o.output_sample_type=f::OutputSample::integer;o.output_depth=8;o.dither_algo=3;
  result=f::Filter(fmt,o).process(0,{view(floats,fmt.width)});
  std::vector<float> residual(floats.size());std::vector<std::uint8_t> expected(floats.size());
  for(std::size_t y=0;y<fmt.height;++y)for(std::size_t x=0;x<fmt.width;++x){
    const auto i=y*fmt.width+x;float error=x?residual[i-1]*0.4375f:0;
    if(y){if(x)error+=residual[i-fmt.width-1]*0.0625f;error+=residual[i-fmt.width]*0.3125f;if(x+1<fmt.width)error+=residual[i-fmt.width+1]*0.1875f;}
    const float target=std::clamp(floats[i]*256+error,0.0f,255.0f);
    const float lower=std::floor(target),part=target-lower;
    const float rounded=part<0.5f?lower:part>0.5f?lower+1:(static_cast<int>(lower)&1)?lower+1:lower;
    expected[i]=static_cast<std::uint8_t>(rounded);residual[i]=target-rounded;
  }
  check(values<std::uint8_t>(result.planes[0])==expected,"FS preserves cross-row and potential tile-boundary errors");
}
void integration_matrix() {
  // Keep every mode/depth pair; rotate layout and dither instead of their
  // Cartesian product. Focused output tests cover dither-specific semantics.
  int case_index=0;
  for(int mode=1;mode<=7;++mode) for(int depth:{8,10,12,16,32,-16}) {
    const int layout=case_index%4, dither=(case_index/4+mode)%3+1;
    ++case_index;
    f::Format fmt;fmt.color=layout?f::Color::yuv:f::Color::gray;fmt.width=18;fmt.height=6;fmt.frames=3;
    fmt.sub_x=layout>=2?1:0;fmt.sub_y=layout==3?1:0;
    fmt.sample=(depth==32 || depth==-16)?f::Sample::floating:f::Sample::integer;fmt.depth=std::abs(depth);
    f::Options o;o.sample_mode=mode;o.dither_algo=dither;o.dynamic_grain=true;o.seed=-113;o.range=255;
    std::vector<std::vector<std::uint8_t>> buffers;std::vector<f::PlaneView> views;
    for(int p=0;p<(layout?3:1);++p) {
      const auto w=fmt.width>>(p?fmt.sub_x:0),h=fmt.height>>(p?fmt.sub_y:0),bytes=static_cast<std::size_t>(fmt.depth==8?1:fmt.depth==32?4:2),stride=(w+3)*bytes;
      buffers.emplace_back(stride*h,0xa5);auto& buffer=buffers.back();
      for(std::size_t y=0;y<h;++y) for(std::size_t x=0;x<w;++x) {
        const auto pattern=(x*17+y*13)%256;auto* ptr=buffer.data()+y*stride+x*bytes;
        if(depth==8) *ptr=static_cast<std::uint8_t>(pattern);
        else if(depth==32) { float value=static_cast<float>(pattern)/256-(p?0.5f:0);std::memcpy(ptr,&value,4); }
        else { const std::uint16_t value=depth==-16?f::encode_half(static_cast<double>(pattern)/256-(p?0.5:0)):static_cast<std::uint16_t>(pattern<<(depth-8));std::memcpy(ptr,&value,2); }
      }
      views.push_back({buffer.data(),buffer.size(),stride});
    }
    const auto untouched=buffers;f::Filter filter(fmt,o);
    auto frame2=filter.process(2,views);filter.process(0,views);auto again=filter.process(2,views);
    for(std::size_t p=0;p<again.planes.size();++p) check(frame2.planes[p].data==again.planes[p].data,"out-of-order determinism matrix");
    check(buffers==untouched,"input immutable including padded rows");
    // Cross-type output cannot change the random/deband processing domain.
    o.output_sample_type=fmt.sample==f::Sample::integer?f::OutputSample::floating:f::OutputSample::integer;
    auto cross=f::Filter(fmt,o).process(1,views);check(cross.planes.size()==frame2.planes.size(),"cross type matrix");
  }
  for(int mode=1;mode<=7;++mode) {
    f::Format fmt;f::Options o=identity();o.sample_mode=mode;
    std::vector<std::uint8_t> single{103};f::Filter filter(fmt,o);
    check(filter.process(0,{view(single,1)}).planes[0].data==single,"1x1 all modes");
    rejects([&]{filter.process(1,{view(single,1)});},"frame out of range");
    rejects([&]{filter.process(0,{{single.data(),0,1}});},"capacity rejected");
  }
  f::Format fmt;fmt.width=32;fmt.height=24;fmt.frames=4;f::Options o;o.dynamic_grain=true;
  std::vector<std::uint8_t> bytes(fmt.width*fmt.height,128);f::Filter filter(fmt,o);
  auto serial=filter.process(3,{view(bytes,32)});
  auto task=std::async(std::launch::async,[&]{return filter.process(3,{view(bytes,32)});});
  auto other=filter.process(1,{view(bytes,32)});(void)other;
  check(task.get().planes[0].data==serial.planes[0].data,"concurrent request determinism");
  o.mt=false;check(f::Filter(fmt,o).process(3,{view(bytes,32)}).planes[0].data==serial.planes[0].data,"mt does not alter results");
}
} // namespace
int main() {
  try { parameters(); mathematics(); random_tables(); output_cases(); integration_matrix(); std::cout<<"Passed "<<checks<<" checks\n"; return 0; }
  catch(const std::exception& e) { std::cerr<<"FAILED after "<<checks<<" checks: "<<e.what()<<'\n'; return 1; }
}
