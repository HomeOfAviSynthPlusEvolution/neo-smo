#include "algorithms/f3kdb.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <future>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
namespace f = neo_smo::f3kdb;
std::size_t checks = 0;
void check(bool b, const char* name) { ++checks; if (!b) throw std::runtime_error(name); }
template<class F> void rejects(F fn, const char* name) {
  bool threw = false; try { fn(); } catch (const std::exception&) { threw = true; }
  check(threw, name);
}
f::Options zero() {
  f::Options o;
  o.threshold.fill(0.0); o.threshold1.fill(0.0); o.threshold2.fill(0.0);
  o.grain_y = 0.0; o.grain_c = 0.0;
  return o;
}
template<class T> T pixel(const f::Plane& p, std::size_t x, std::size_t y=0) {
  T v; std::memcpy(&v,p.data.data()+y*p.stride+x*sizeof(T),sizeof(T)); return v;
}
double uniform(std::uint32_t& s) {
  s = s*1664525u+1013904223u;
  const std::uint64_t m=(std::uint64_t(s)<<20)|(s>>12);
  return std::ldexp(double(m),-51)-1.0;
}
void halves() {
  for (unsigned h=0;h<65536;++h) {
    const unsigned e=(h>>10)&31, m=h&1023;
    if (e==31) continue;
    const double sign=(h&32768)?-1.0:1.0;
    const double exact=sign*(e?std::ldexp(double(1024+m),int(e)-25):std::ldexp(double(m),-24));
    const double got=f::decode_half(std::uint16_t(h));
    check(got==exact,"half decode analytical oracle");
    check(f::encode_half(exact)==h,"half finite roundtrip");
  }
  check(f::encode_half(1.0+std::ldexp(1.0,-11))==0x3c00,"half tie even down");
  check(f::encode_half(1.0+3*std::ldexp(1.0,-11))==0x3c02,"half tie even up");
  check(f::encode_half(1e30)==0x7bff,"half finite saturation positive");
  check(f::encode_half(-1e30)==0xfbff,"half finite saturation negative");
}
void numerical() {
  check(f::response(0,0)==1 && f::response(1,0)==0,"response zero threshold");
  check(std::abs(f::response(5e-100,6e-100)-.5)<1e-14,"response no arbitrary epsilon");
  check(f::weight({0,1,1,1})==0 && f::weight({1,1,1,1})==1,"weight endpoints");
  double prev=0;
  for(int i=-3000;i<=0;++i) {
    const double c=std::pow(10.0,i/10.0);
    const double v=f::weight({c,c,c,c});
    const long double ref=std::exp(0.4L*std::log(static_cast<long double>(c)));
    check(v>=prev && v>=0 && v<=1,"weight monotone and bounded");
    check(std::abs(static_cast<long double>(v)-ref)<=std::ldexp(1.0L,-20),"weight independent logarithmic oracle");
    prev=v;
  }
  check(f::weight({1e-12,1e-12,1e-12,1e-12})>1e-5,"weight avoids product underflow");
  check(f::aligned({0,0},{0,0},0),"both zero align");
  check(!f::aligned({0,0},{1,0},1),"one zero rejects at max angle");
  check(f::aligned({1,0},{0,1},.5),"orthogonal max angle");
  check(!f::aligned({1,0},{0,1},.25),"orthogonal 45 rejects");
  check(f::aligned({1,0},{1,1},.25),"45 exact equality");
  check(f::aligned({100,1},{-100,1},.01),"direction periodicity");
  check(f::aligned({1e-300,0},{-1e-300,0},0),"tiny antiparallel vectors");
  check(!f::aligned({1e-300,0},{0,1e-300},.25),"tiny nonzero vectors are not flat");
  check(!f::aligned({1e300,0},{0,1e300},.25),"large vectors no infinity comparison");
  const double tiny=std::numeric_limits<float>::denorm_min();
  for(int mode=1;mode<=7;++mode)
    check(f::deband(tiny,{tiny,tiny,tiny,tiny},mode,{1,1,1},false,true,true)==tiny,"subnormal constant preserved");
  for (bool integer : {false,true}) {
    std::array<double,4> a{0,0,0,1};
    do { check(f::deband(.25,a,2,{100,100,100},integer,true)==(integer?0.0:.25),"four average permutation invariant"); }
    while(std::next_permutation(a.begin(),a.end()));
    check(f::deband(100,{90,90,130,130},4,{11,11,11},integer,true)==95,"mode4 separate directions");
    check(f::deband(100,{98,98,98,98},5,{3,3,4},integer,true)==100,"mode5 midpoint equality reject");
    check(f::deband(100,{98,98,98,98},5,{3,3,5},integer,false)==98,"mode5 midpoint passes");
    check(f::deband(100,{104,104,104,104},6,{100,100,100},integer,true)==104,"mode6 fully blended");
    check(f::deband(100,{104,104,104,104},6,{4,100,100},integer,true)==100,"mode6 hard zero response");
  }
}
void parameters() {
  f::Format fmt;
  {
    const auto p=f::normalize(fmt,{});
    check(p.thresholds[0]==std::array<double,3>{256,256,256},"integer auxiliary inheritance before scaling");
    auto floats=fmt;floats.sample=f::Sample::floating;floats.depth=32;
    const auto pf=f::normalize(floats,{});
    check(pf.thresholds[0]==std::array<double,3>{1.0/256,1.0/256,1.0/256},"floating auxiliary inheritance before scaling");
  }
  for(bool scale:{false,true}) {
    f::Options o; o.scale=scale; o.presets={"high"};
    auto p=f::normalize(fmt,o);
    check(p.thresholds[0][0]==256 && p.grain[0]==256,"preset semantic scale independence");
    o.presets={"high","luma"};
    p=f::normalize(fmt,o);
    check(p.thresholds[1]==std::array<double,3>{0,0,0} && p.thresholds[2]==std::array<double,3>{0,0,0} && p.grain[1]==0,"luma clears all chroma thresholds");
  }
  auto o=zero(); o.output_sample_type=f::OutputSample::floating;
  check(f::normalize(fmt,o).output.depth==32,"cross type default float32");
  fmt.sample=f::Sample::floating; fmt.depth=16; o.output_sample_type=f::OutputSample::integer;
  check(f::normalize(fmt,o).output.depth==16,"cross type default int16");
  o.output_depth=-2; rejects([&]{f::normalize(fmt,o);},"invalid negative depth");
  o=zero(); o.random_param_ref=-1; rejects([&]{f::normalize(fmt,o);},"negative sigma");
  o=zero(); o.random_param_ref=std::numeric_limits<double>::infinity(); rejects([&]{f::normalize(fmt,o);},"unused nonfinite sigma");
  o=zero(); o.dither_algo=4; rejects([&]{f::normalize(fmt,o);},"invalid dither even float output");
  o=zero(); o.range=256; rejects([&]{f::normalize(fmt,o);},"range upper bound");
  o=zero(); o.threshold[0]=std::numeric_limits<double>::quiet_NaN(); rejects([&]{f::normalize(fmt,o);},"nan threshold");
}
void random_tables() {
  std::uint32_t s=0;
  f::Random rng(0);
  for(int i=0;i<1000;++i) check(rng.next(1,1)==uniform(s) && rng.state()==s,"uniform state and mantissa mapping");
  f::Random g(12); std::uint32_t gs=12; uniform(gs);uniform(gs);
  check(g.next(2,0)==0 && g.state()==gs,"zero sigma consumes two uniform updates");
  for(double sigma:{1e-300,.1,.99,1.,100.,1e300}) {
    f::Random x(12);
    for(int i=0;i<100;++i) { const double v=x.next(2,sigma); check(std::isfinite(v) && v>-1 && v<1,"gaussian finite support extreme sigma"); }
  }
  f::Format fmt; fmt.width=600; fmt.height=600;
  auto o=zero();o.range=255; o.sample_mode=2;
  f::Filter filter(fmt,o);
  bool saw_large=false;
  const auto& d=filter.distances();
  check(d.size()==fmt.width*fmt.height,"distance table dimensions");
  for(std::size_t y=0;y<fmt.height;++y) for(std::size_t x=0;x<fmt.width;++x) {
    const unsigned bound=unsigned(std::min({std::size_t(255),x,y,fmt.width-1-x,fmt.height-1-y}));
    const auto v=d[y*fmt.width+x];
    check(v.first<=bound && v.second<=bound,"distance geometry bounds");
    saw_large=saw_large || v.first>127 || v.second>127;
  }
  check(saw_large,"full range survives uint8 storage");
}
void frames() {
  f::Format fmt;fmt.sample=f::Sample::floating;fmt.depth=32;fmt.width=2;
  auto o=zero();o.output_sample_type=f::OutputSample::integer;o.output_depth=8;o.dither_algo=3;
  const std::array<float,2> src{.75f/256,.75f/256};
  f::Filter fs(fmt,o);
  auto out=fs.process(0,{{src.data(),sizeof(src),sizeof(src)}});
  check(pixel<std::uint8_t>(out.planes[0],0)==1 && pixel<std::uint8_t>(out.planes[0],1)==1,"signed FS golden row");
  fmt.width=1;
  const float outside=1.25f;
  o=zero();f::Filter full(fmt,o);
  out=full.process(0,{{&outside,4,4}});
  check(pixel<float>(out.planes[0],0)==outside,"float overshoot not clipped");
  o.keep_tv_range=true;f::Filter tv(fmt,o);
  out=tv.process(0,{{&outside,4,4}});
  check(pixel<float>(out.planes[0],0)==235.f/256,"TV applies on identity");
  for(float bad:{std::numeric_limits<float>::infinity(),std::numeric_limits<float>::quiet_NaN()})
    rejects([&]{full.process(0,{{&bad,4,4}});},"nonfinite pixel rejected even identity");
  rejects([&]{full.process(1,{{&outside,4,4}});},"frame upper bound");
  rejects([&]{full.process(0,{{&outside,3,4}});},"short buffer");
  rejects([&]{full.process(0,{{nullptr,4,4}});},"null input");

  fmt.color=f::Color::yuv;fmt.width=2;fmt.height=2;fmt.sub_x=1;fmt.sub_y=1;
  o=zero();o.output_sample_type=f::OutputSample::integer;o.output_depth=8;o.dither_algo=1;
  f::Filter yuv(fmt,o);
  const std::array<float,4> yy{.5f,.5f,.5f,.5f}; const float cc=0;
  out=yuv.process(0,{{yy.data(),sizeof(yy),8},{&cc,4,4},{&cc,4,4}});
  check(pixel<std::uint8_t>(out.planes[0],0)==128 && pixel<std::uint8_t>(out.planes[1],0)==128,"float chroma offset only at integer output");

  fmt={};fmt.width=9;fmt.height=7;fmt.frames=3;
  o=zero();o.grain_y=8;o.dynamic_grain=true;o.dither_algo=1;
  f::Filter dynamic(fmt,o);
  std::vector<std::uint8_t> a(fmt.width*fmt.height,128);
  const std::vector<f::PlaneView> views{{a.data(),a.size(),fmt.width}};
  const auto first=dynamic.process(2,views).planes[0].data;
  dynamic.process(0,views);dynamic.process(1,views);
  check(dynamic.process(2,views).planes[0].data==first,"out of order deterministic");
  auto task=std::async(std::launch::async,[&]{return dynamic.process(2,views).planes[0].data;});
  check(dynamic.process(2,views).planes[0].data==first && task.get()==first,"concurrent deterministic");
}
int main() {
  try { halves();numerical();parameters();random_tables();frames(); }
  catch(const std::exception& e) { std::cerr<<"AUDIT FAIL after "<<checks<<": "<<e.what()<<'\n';return 1; }
  std::cout<<"Independent audit passed "<<checks<<" checks\n";
}
