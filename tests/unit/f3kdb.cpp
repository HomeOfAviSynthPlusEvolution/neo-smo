#include "f3kdb_scalar_filter.hpp"
#include "algorithms/f3kdb.hpp"
#include "kernels/f3kdb.hpp"
#include "kernels/f3kdb_weight_table.hpp"
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
  for (int pattern=0;pattern<9;++pattern) {
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
    if (pattern>=5)
      for (std::size_t i=0;i<n;++i) for (int j=0;j<4;++j)
        c[j][i]=i%2?1.0:(j==int(i%4)?0.0:0.375);
    if (pattern==6) c[0][n/2]=0.25; // an interior lane keeps the logarithmic path active
    if (pattern==7 || pattern==8) {
      // Unequal responses exercise the single-log product and the transition
      // to the logarithmic underflow fallback, including mixed endpoint lanes.
      std::mt19937 rng(0x57454947);
      for (std::size_t i=0;i<n;++i) for (int j=0;j<4;++j) {
        const double exponent = pattern==7 ? double(rng()%2000)/100 : 170+double(rng()%2000)/100;
        c[j][i]=std::exp(-exponent);
      }
      c[0][0]=0;
      for (int j=0;j<4;++j) c[j][1]=1;
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
void normal_weight_tests() {
  // Exercise every normal exponent, mantissa interiors, and adjacent doubles
  // at exponent joins. The reference uses libm, not the production polynomial.
  constexpr std::size_t n=260;
  std::array<std::vector<double>,4> c;
  for (auto& v:c) v.assign(n,1.0);
  std::vector<double> output(n);
  const double* ptr[4]={c[0].data(),c[1].data(),c[2].data(),c[3].data()};
  for (int e=-1022;e<=0;++e) {
    const double base=std::ldexp(1.0,e);
    c[0][0]=std::nextafter(base,0.0);
    c[0][1]=base;
    c[0][2]=std::min(1.0,std::nextafter(base,1.0));
    for (std::size_t i=3;i<n;++i)
      c[0][i]=std::min(1.0,std::ldexp(1.0+double(i-3)/256,e));
    std::sort(c[0].begin(),c[0].end());
    simd::weights(ptr,output.data(),n);
    for (std::size_t i=0;i<n;++i) {
      const auto expected=std::pow(c[0][i],0.1);
      const auto error=std::abs(output[i]-expected);
      max_weight=std::max(max_weight,error);
      check(std::isfinite(output[i]) && output[i]>=0 && output[i]<=1,"normal weight range");
      check(error<=0x1p-20,"normal weight accuracy");
      if (i) check(output[i]>=output[i-1],"normal weight exponent-join monotonicity");
      if (c[0][i]==1) check(output[i]==1,"normal weight exact one");
    }
  }
}
void sampling_tests() {
  // Unique source values expose wrong neighbours even when a later threshold
  // would reject them. Keep the distance buffer exact-sized for ASan.
  for (std::size_t width : {1, 3, 7, 17, 519}) {
    constexpr std::size_t height = 515;
    std::vector<double> source(width * height);
    for (std::size_t i = 0; i < source.size(); ++i) source[i] = double(i);
    std::vector<std::int32_t> packed_source(source.begin(), source.end());
    std::array<std::int64_t, 256> vertical{};
    for (std::size_t i = 1; i < vertical.size(); ++i) vertical[i] = vertical[i - 1] + width;
    for (int sx : {0, 1, 2}) for (int sy : {0, 1}) for (int mode : {1, 2, 3}) {
      for (std::size_t y : {0, 1, 257, 514}) {
        const auto step = std::size_t{1} << sx;
        std::vector<ref::Distance> distances(width * step, {255, 255});
        std::array<std::vector<double>, 4> samples;
        for (auto& sample : samples) sample.assign(width + 1, -1);
        const auto vlimit = std::min(y, height - 1 - y);
        for (std::size_t x = 0; x < width; ++x) {
          const auto hlimit = std::min(x, width - 1 - x);
          auto bound = std::size_t{255};
          if (mode != 3) bound = std::min(bound, vlimit << sy);
          if (mode != 1) bound = std::min(bound, hlimit << sx);
          distances[x * step] = {static_cast<std::uint8_t>(bound), static_cast<std::uint8_t>(bound / 2)};
        }
        simd::SamplingRow row{source.data() + y * width, distances.data(), vertical.data(),
                              {samples[0].data(), samples[1].data(), samples[2].data(), samples[3].data()},
                              width, vlimit, sx, sy};
        const auto count = simd::collect_simple_row(row, mode);
        check(count <= width && (sx < 2 || count == 0), "sampling prefix extent");
        for (std::size_t x = 0; x < count; ++x) {
          const auto first = distances[x * step].first, second = distances[x * step].second;
          const auto h1 = std::int64_t(first >> sx), v1 = std::int64_t(first >> sy);
          const auto h2 = std::int64_t(second >> sx), v2 = std::int64_t(second >> sy);
          std::array<std::int64_t, 4> offsets{};
          if (mode == 1) offsets = {v1 * std::int64_t(width), -v1 * std::int64_t(width), 0, 0};
          else if (mode == 3) offsets = {h1, -h1, 0, 0};
          else {
            const auto a = v2 * std::int64_t(width) + h1, b = h2 - v1 * std::int64_t(width);
            offsets = {a, b, -a, -b};
          }
          for (int j = 0; j < (mode == 2 ? 4 : 2); ++j)
            check(samples[j][x] == double(std::int64_t(y * width + x) + offsets[j]), "sampling neighbour");
        }
        for (int j = 0; j < 4; ++j)
          for (std::size_t x = j >= (mode == 2 ? 4 : 2) ? 0 : count; x <= width; ++x)
            check(samples[j][x] == -1, "sampling tail and unused outputs");
        for (bool blur : {false, true}) for (double threshold : {0.0, 0.5, 1.0, 255.25, 65535.0}) {
          std::vector<double> output(width + 1, -999);
          std::vector<std::int16_t> grain(width);
          for (std::size_t x = 0; x < width; ++x) grain[x] = double(int(x % 31) - 15);
          simd::Row batch{}; batch.center = row.center; batch.integer_grain = grain.data(); batch.output = output.data(); batch.count = width;
          const auto processed = simd::process_sampled_row(row, batch, mode, threshold, true, blur).processed;
          check(processed <= count, "fused sampling extent");
          for (std::size_t x = 0; x < processed; ++x) {
            std::array<double, 4> refs{samples[0][x], samples[1][x], mode == 2 ? samples[2][x] : row.center[x], mode == 2 ? samples[3][x] : row.center[x]};
            const auto expected = ref::deband(row.center[x], refs, mode, {threshold, threshold, threshold}, true, blur) + grain[x];
            check(output[x] == expected, "fused sampling value");
          }
          for (std::size_t x = processed; x <= width; ++x) check(output[x] == -999, "fused sampling tail");
          std::vector<float> floating_grain(width);
          for (std::size_t x = 0; x < width; ++x) floating_grain[x] = float(grain[x]) / 65536;
          batch.grain = floating_grain.data();
          std::fill(output.begin(), output.end(), -999);
          const auto floating_processed = simd::process_sampled_row(row, batch, mode, threshold, false, blur).processed;
          check(floating_processed <= count, "fused floating sampling extent");
          for (std::size_t x = 0; x < floating_processed; ++x) {
            std::array<double, 4> refs{samples[0][x], samples[1][x], mode == 2 ? samples[2][x] : row.center[x], mode == 2 ? samples[3][x] : row.center[x]};
            const auto expected = double(float(ref::deband(row.center[x], refs, mode, {threshold, threshold, threshold}, false, blur) + floating_grain[x]));
            check(output[x] == expected, "fused floating sampling value");
          }
          for (std::size_t x = floating_processed; x <= width; ++x)
            check(output[x] == -999, "fused floating sampling tail");
          std::fill(output.begin(), output.end(), -999);
          const simd::IntegerRow packed{packed_source.data() + y * width, distances.data(), grain.data(), output.data(), width, vlimit, y * width, sx, sy};
          const std::array<double, 3> thresholds{threshold, threshold, threshold};
          const auto prefix = simd::process_integer_row(packed, mode, thresholds.data(), blur);
          check(prefix <= width && (sx < 2 || prefix == 0), "packed sampling extent");
          for (std::size_t x = 0; x < prefix; ++x) {
            const auto d = distances[x * step];
            const auto h1 = std::int64_t(d.first >> sx), v1 = std::int64_t(d.first >> sy);
            const auto h2 = std::int64_t(d.second >> sx), v2 = std::int64_t(d.second >> sy), w = std::int64_t(width);
            const auto a = mode == 1 ? v1 * w : mode == 3 ? h1 : v2 * w + h1, b = h2 - v1 * w;
            const std::array<std::int64_t, 4> offsets = mode == 2 ? std::array<std::int64_t, 4>{a,b,-a,-b} : std::array<std::int64_t, 4>{a,-a,0,0};
            std::array<double, 4> refs{};
            for (int j = 0; j < 4; ++j) refs[j] = source[std::int64_t(y * width + x) + offsets[j]];
            check(output[x] == ref::deband(source[y * width + x], refs, mode, thresholds, true, blur) + grain[x], "packed sampling value");
          }
          for (std::size_t x = prefix; x <= width; ++x) check(output[x] == -999, "packed sampling tail");
        }
        if (count && y == 0) {
          distances[0].first = 255;
          bool rejected = false;
          try { simd::collect_simple_row(row, mode); } catch (const std::logic_error&) { rejected = true; }
          check(rejected, "sampling first distance bounds");
          if (mode == 2) {
            distances[0] = {0, 255}; rejected = false;
            try { simd::collect_simple_row(row, mode); } catch (const std::logic_error&) { rejected = true; }
            check(rejected, "sampling second distance bounds");
          }
        }
      }
    }
  }
}
void gradient_sampling_tests() {
  std::mt19937 rng(0x675ad);
  for (const auto shape : {std::array<std::size_t,2>{1,3},{3,7},{17,37},{519,515}}) {
    const auto w=shape[0],h=shape[1];
    std::vector<double> source(w*h),gx(w*h),gy(w*h);
    std::vector<std::int64_t> vertical(256);
    for (std::size_t y=0;y<256;++y) vertical[y]=std::int64_t(y*w);
    std::array<std::vector<double>,4> values,rx,ry;
    for (auto* array : {&values,&rx,&ry}) for (auto& v:*array) v.resize(w+1);
    std::vector<double> output(w+1);
    std::vector<std::int16_t> integer_grain(w);
    std::vector<float> grain(w);
    for (std::size_t x=0;x<w;++x) { integer_grain[x]=int(x%31)-15;grain[x]=float(int(x%31)-15)/65536; }
    for (int pattern=0;pattern<8;++pattern) {
      const bool integer=pattern==0 || pattern>=4;
      for (std::size_t i=0;i<source.size();++i) {
        source[i]=int(rng()%65536);gx[i]=int(rng()%524281)-262140;gy[i]=int(rng()%524281)-262140;
        if (pattern==4) source[i]=i%65536;
        if (pattern==5) source[i]=(rng()%256)*256;
        if (pattern==6) { gx[i]=(i%4)?262140:0;gy[i]=(i%5)?-262140:0; }
        if (pattern==7) { gx[i]=0;gy[i]=0; }
        if (pattern==1) { source[i]=float((source[i]-32768)/65536);gx[i]*=1e-25;gy[i]*=1e-25; }
        if (pattern==2) { source[i]=(i%2?-1:1)*double(std::numeric_limits<float>::max());gx[i]*=1e30;gy[i]*=1e30; }
        if (pattern==3) { source[i]=std::numeric_limits<float>::denorm_min()*(1+i%3);gx[i]=0;gy[i]=0; }
      }
      std::vector<std::int32_t> compact_source, compact_gx, compact_gy;
      if (integer) {
        compact_source.assign(source.begin(), source.end());
        compact_gx.assign(gx.begin(), gx.end());
        compact_gy.assign(gy.begin(), gy.end());
      }
      for (int sx : {0,1,2}) for (int sy : {0,1}) for (std::size_t y : {std::size_t{0},h/2,h-1}) {
        std::vector<ref::Distance> distances(w<<sx);
        const auto vlimit=std::min(y,h-1-y);
        for (std::size_t x=0;x<w;++x)
          distances[x<<sx]={std::uint8_t(std::min({std::size_t{255},vlimit<<sy,std::min(x,w-1-x)<<sx})),255};
        for (auto* array : {&values,&rx,&ry}) for (auto& v:*array) std::fill(v.begin(),v.end(),-999999);
        simd::CrossSamplingRow row{{source.data()+y*w,distances.data(),vertical.data(),
                                   {values[0].data(),values[1].data(),values[2].data(),values[3].data()},w,vlimit,sx,sy},
                                  gx.data()+y*w,gy.data()+y*w,{rx[0].data(),rx[1].data(),rx[2].data(),rx[3].data()},
                                  {ry[0].data(),ry[1].data(),ry[2].data(),ry[3].data()},y*w};
        const auto count=simd::collect_cross_row(row);
        check(count<=w && (sx<2 || count==0),"gradient prefix extent");
        for (std::size_t x=0;x<count;++x) {
          const auto dist=distances[x<<sx].first;const std::size_t dx=dist>>sx,dy=dist>>sy;
          const std::size_t indices[4]={(y+dy)*w+x,(y-dy)*w+x,y*w+x+dx,y*w+x-dx};
          for (int j=0;j<4;++j) {
            check(values[j][x]==source[indices[j]],"gradient sample coordinate");
            check(rx[j][x]==gx[indices[j]] && ry[j][x]==gy[indices[j]],"gradient coordinate");
          }
        }
        for (auto* array : {&values,&rx,&ry}) for (auto& v:*array)
          for (std::size_t x=count;x<=w;++x) check(v[x]==-999999,"gradient collection tail");
        simd::Row batch{};batch.center=row.sampling.center;batch.output=output.data();batch.count=w;
        batch.grain=grain.data();batch.integer_grain=integer_grain.data();
        // Exercise both wide-coordinate contracts and the bounded compact
        // planes against the same independent coordinate/direction oracle.
        for (int storage : {0,1,2}) {
          const bool known=storage!=0;
          if (known && !integer) continue;
          row.bounded_integer_gradients=known;
          row.compact_source=storage==2?compact_source.data():nullptr;
          row.compact_gx=storage==2?compact_gx.data():nullptr;
          row.compact_gy=storage==2?compact_gy.data():nullptr;
        for (double angle : {0.0,0.15,0.25,0.5}) for (double boost : {0.0,1.5,4.0}) {
          const double threshold=pattern==5?256:integer?65535:0.5;
          const std::array<double,3> t{threshold,threshold*0.75,threshold*0.5};
          std::fill(output.begin(),output.end(),-999999);
          const auto prefix=simd::process_cross_row(row,batch,t.data(),integer,angle,boost);
          check(prefix.collected==count && prefix.processed==count,"gradient fused extent");
          for (std::size_t x=0;x<count;++x) {
            bool aligned=true;
            for (int j=0;j<4;++j) aligned=aligned && ref::aligned({row.gx[x],row.gy[x]},{rx[j][x],ry[j][x]},angle);
            const std::array<double,4> refs{values[0][x],values[1][x],values[2][x],values[3][x]};
            double expected=ref::deband(batch.center[x],refs,7,t,integer,true,aligned,boost);
            expected=integer?expected+integer_grain[x]:double(float(expected+grain[x]));
            const double tolerance=integer?1:8*std::numeric_limits<float>::epsilon()*std::max(std::abs(expected),std::abs(batch.center[x]))+8*std::numeric_limits<float>::denorm_min();
            check(std::isfinite(output[x]) && std::abs(output[x]-expected)<=tolerance,"gradient fused value");
          }
          for (std::size_t x=count;x<=w;++x) check(output[x]==-999999,"gradient fused tail");
        }
        }
        if (count && y==0) {
          distances[0].first=255;bool rejected=false;
          try { simd::collect_cross_row(row); } catch(const std::logic_error&) { rejected=true; }
          check(rejected,"gradient bounds checked");
          rejected=false;const std::array<double,3> t{1,1,1};
          try { simd::process_cross_row(row,batch,t.data(),integer,0.15,1.5); } catch(const std::logic_error&) { rejected=true; }
          check(rejected,"gradient fused bounds checked");
        }
      }
    }
  }
}
void integer_cross_tests() {
  // Four distinct neighbours, full-range codes, 255-pixel distances and a
  // partial vector expose pair ordering and both stages of mode 4 rounding.
  constexpr std::size_t width=519,height=515;
  std::mt19937 rng(0x4f3db);
  std::vector<std::int32_t> source(width*height);
  for (auto& value : source) value=static_cast<std::int32_t>(rng()%65536);
  std::vector<std::int16_t> grain(width);
  for (std::size_t x=0;x<width;++x) grain[x]=x%2?16384:-16384;
  std::vector<double> output(width+1);
  for (int sx : {0,1}) for (int sy : {0,1}) for (int mode : {4,5,6})
    for (std::size_t y : {0,1,257,514}) {
      const auto vlimit=std::min(y,height-1-y);
      std::vector<ref::Distance> distances(width<<sx);
      for (std::size_t x=0;x<width;++x)
        distances[x<<sx]={static_cast<std::uint8_t>(std::min({std::size_t{255},vlimit<<sy,std::min(x,width-1-x)<<sx})),255};
      const simd::IntegerRow row{source.data()+y*width,distances.data(),grain.data(),output.data(),width,vlimit,y*width,sx,sy};
      for (bool blur : {false,true}) for (double threshold : {0.0,0.5,4.0,2048.0,1500.25,65532.0,65535.0,131070.5}) {
        const std::array<double,3> thresholds{threshold,threshold*0.75,threshold*0.5};
        std::fill(output.begin(),output.end(),-999999);
        const auto count=simd::process_integer_row(row,mode,thresholds.data(),blur);
        for (std::size_t x=0;x<count;++x) {
          const auto distance=distances[x<<sx].first;
          const std::size_t dx=distance>>sx,dy=distance>>sy;
          const std::array<double,4> refs{double(source[(y+dy)*width+x]),double(source[(y-dy)*width+x]),
                                        double(source[y*width+x+dx]),double(source[y*width+x-dx])};
          const auto expected=ref::deband(source[y*width+x],refs,mode,thresholds,true,blur)+grain[x];
          check(std::abs(output[x]-expected)<=(mode==6?1:0),"integer cross coordinate oracle");
        }
        for (std::size_t x=count;x<output.size();++x) check(output[x]==-999999,"integer cross tail untouched");
      }
    }
}
void integer_weight_table_tests() {
  const auto& table=simd::integer_weight_table();
  double worst=0;
  for(std::size_t i=0;i+1<table.size();++i) {
    float pair[2],next[2];std::memcpy(pair,&table[i],sizeof(pair));std::memcpy(next,&table[i+1],sizeof(next));
    const float upper=pair[0]+pair[1];
    check(upper==next[0],"integer root table shared endpoints");
    check(pair[1]>=0,"integer root table monotonic endpoints");
    float previous=pair[0];
    for(unsigned fraction : {0u,1u,8192u,16384u,24576u,32766u,32767u}) {
      const std::uint32_t bits=simd::integer_weight_first_bits+(std::uint32_t(i)<<15)+fraction;
      float x;std::memcpy(&x,&bits,sizeof(x));
      const float t=float(fraction)*0x1p-15f;
      const float actual=std::fma(t,pair[1],pair[0]);
      check(actual>=previous && actual<=upper,"integer root interpolation monotonic");
      previous=actual;
      worst=std::max(worst,std::abs(double(actual)-std::pow(double(x),0.1)));
    }
  }
  check(worst<2.4e-7,"integer root interpolation error");
  float endpoint[2];std::memcpy(endpoint,&table[66*256],sizeof(endpoint));
  check(endpoint[0]==1.0f,"integer root exact one");
  max_weight=std::max(max_weight,worst);
}
void integer_weight_boundary_tests() {
  // Smooth and full-range neighbours exercise interior weights, quarter-code
  // response distances, integral thresholds, and the fractional fallback.
  constexpr std::size_t width=1027;
  std::vector<std::int32_t> source(3*width);
  std::vector<ref::Distance> distances(width,{1,0});
  distances.front().first=distances.back().first=0;
  std::vector<std::int16_t> grain(width,0);
  std::vector<double> output(width+1);
  std::mt19937 rng(0x6f32);
  for(double threshold : {1.0,2.0,3.0,4.0,108.0,2048.0,32768.0,65535.0,
                         std::nextafter(108.0,109.0),2048.25}) {
    const std::array<double,3> thresholds{threshold,threshold,threshold};
    for(int pattern=0;pattern<12;++pattern) {
      for(auto& value:source) {
        if(pattern<8) {
          const int radius=std::min(32767,int(threshold));
          value=32768+int(rng()%unsigned(2*radius+1))-radius;
        } else value=int(rng()%65536);
      }
      const simd::IntegerRow row{source.data()+width,distances.data(),grain.data(),output.data(),width,1,width,0,0};
      std::fill(output.begin(),output.end(),-999999);
      const auto count=simd::process_integer_row(row,6,thresholds.data(),true);
      for(std::size_t x=0;x<count;++x) {
        const auto offset=distances[x].first;
        const std::array<double,4> refs{double(source[(1+offset)*width+x]),
          double(source[(1-offset)*width+x]),double(source[width+x+offset]),double(source[width+x-offset])};
        const double expected=ref::deband(source[width+x],refs,6,thresholds,true,true);
        const double error=std::abs(output[x]-expected);
        max_integer=std::max(max_integer,error);
        check(error<=1,"integer weighted boundary oracle");
      }
      for(std::size_t x=count;x<output.size();++x)
        check(output[x]==-999999,"integer weighted boundary tail");
    }
  }
}
void horizontal_table_tests() {
  std::mt19937 rng(0x7ab1e);
  for (std::size_t w : {31,32,47,48,49,63,64,65,519}) {
    std::vector<std::int32_t> source(w);
    for (auto& v:source) v=rng()%65536;
    std::vector<std::int16_t> grain(w);
    for(std::size_t x=0;x<w;++x) grain[x]=int(x%32769)-16384;
    std::vector<double> output(w+1);
    for (int sx : {0,1}) for (std::size_t range : {15,16,17,127}) {
      std::vector<ref::Distance> distances(w<<sx);
      for(std::size_t x=0;x<w;++x) distances[x<<sx]={std::uint8_t(std::min({range,x,w-1-x})<<sx),0};
      simd::IntegerRow row{source.data(),distances.data(),grain.data(),output.data(),w,0,0,sx,0};
      for(bool blur : {false,true}) for(double threshold : {0.0,0.5,100.25,65535.0}) {
        const std::array<double,3> t{threshold,threshold,threshold};
        std::fill(output.begin(),output.end(),-999999);
        const auto count=simd::process_integer_row(row,3,t.data(),blur);
        for(std::size_t x=0;x<count;++x) {
          const auto dx=distances[x<<sx].first>>sx;
          const std::array<double,4> refs{double(source[x+dx]),double(source[x-dx]),0,0};
          check(output[x]==ref::deband(source[x],refs,3,t,true,blur)+grain[x],"horizontal table distance/lane boundary");
        }
        for(std::size_t x=count;x<=w;++x) check(output[x]==-999999,"horizontal table tail");
      }
    }
  }
}
void integer_pack_tests() {
  // Exhaust the processing range including negative grain, half-way ties,
  // saturation and every output depth. An odd output address checks word stores.
  constexpr std::size_t width=98303;
  std::vector<std::int32_t> source(width);
  std::vector<std::int16_t> grain(width);
  std::vector<ref::Distance> distances(width);
  std::vector<double> intermediate(width,-999999);
  for (std::size_t x=0;x<width;++x) {
    const int value=int(x)-16384;
    source[x]=std::clamp(value,0,65535);
    grain[x]=static_cast<std::int16_t>(value-source[x]);
  }
  const std::array<double,3> thresholds{};
  for (int depth=8;depth<=16;++depth) for (bool chroma : {false,true}) for (bool limited : {false,true})
    for (int dither : {1,2}) {
    const std::size_t bytes=depth==8?1:2;
    std::vector<std::uint8_t> output(width*bytes+2,0xa5);
    simd::IntegerRow row{source.data(),distances.data(),grain.data(),intermediate.data(),width,0,0,0,0};
    row.packed_output=output.data()+1; row.output_depth=depth; row.chroma=chroma; row.limited=limited; row.dither=dither; row.y=depth+3;
    const auto count=simd::process_integer_row(row,1,thresholds.data(),true);
    const double scale=std::ldexp(1.0,depth-16);
    const double low=limited?std::ldexp(16.0,depth-8):0;
    const double high=limited?std::ldexp(chroma?240.0:235.0,depth-8):std::ldexp(1.0,depth)-1;
    for (std::size_t x=0;x<count;++x) {
      const double value=(int(x)-16384)*scale;
      const auto expected=std::clamp(dither==1?ref::round_even(value):
          std::floor(value+(ref::bayer16(row.y,x)+0.5)/256),low,high);
      const auto actual=depth==8?output[x+1]:read<std::uint16_t>(output.data()+1+x*2);
      check(actual==expected,"pack value");
    }
    check(output[0]==0xa5,"pack leading guard");
    for (std::size_t i=1+count*bytes;i<output.size();++i) check(output[i]==0xa5,"pack tail guard");
  }
  check(std::all_of(intermediate.begin(),intermediate.end(),[](double x){return x==-999999;}),"direct pack skips double output");
}
void row_tests() {
  constexpr std::size_t n=17;
  std::mt19937 rng(0x3f3db);
  std::array<std::vector<double>,16> storage;
  for (auto& v:storage) v.resize(n);
  std::vector<double> out(n);
  std::vector<float> grain(n);
  std::vector<std::int16_t> integer_grain(n);
  simd::Row row{}; row.center=storage[0].data(); row.count=n; row.output=out.data(); row.grain=grain.data(); row.integer_grain=integer_grain.data();
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
void identity_frame_tests() {
  // Test the exact input-code boundary in processing units. Wider-than-vector
  // rows exercise the ordinary scalar tail beside the identity SIMD prefix.
  int variant=0;
  for (int depth=8;depth<=16;++depth) for (int mode : {5,6}) {
    const double step=std::ldexp(1.0,16-depth);
    for (double threshold : {0.0,std::nextafter(step,0.0),step,std::nextafter(step,step*2),step*3})
      for (int dither : {1,2,3}) {
        ref::Format f; f.width=34; f.height=6; f.frames=3; f.depth=depth;
        f.color=ref::Color::yuv; f.sub_x=1; f.sub_y=1;
        auto in=input(f,0);
        // Adjacent codes produce both equal and one-step-different references,
        // rather than only large deltas rejected on either side of the boundary.
        for (std::size_t p=0;p<in.views.size();++p) {
          auto* data=in.bytes[p].data()+1;
          const auto w=f.width>>(p?1:0), h=f.height>>(p?1:0);
          const std::size_t bytes=depth==8?1:2;
          for (std::size_t y=0;y<h;++y) for (std::size_t x=0;x<w;++x) {
            const auto code=std::uint16_t((1<<(depth-1))+(x+3*y)%3);
            auto* out=data+y*in.views[p].stride+x*bytes;
            if (bytes==1) *out=static_cast<std::uint8_t>(code); else write(out,code);
          }
        }
        ref::Options o; o.sample_mode=mode; o.scale=true; o.threshold.fill(65535);
        o.threshold1.fill(threshold); o.threshold2.fill(65535); o.range=255;
        o.grain_y=3.75; o.grain_c=2.5; o.dynamic_grain=true; o.dither_algo=dither;
        o.blur_first=(variant&1)!=0; o.keep_tv_range=(variant&2)!=0;
        o.output_sample_type=(variant&4)?ref::OutputSample::floating:ref::OutputSample::integer;
        o.output_depth=(variant&4)?((variant&1)?16:32):((variant&1)?8:16);
        ++variant;
        compare(ScalarFilter(f,o).process(2,in.views),simd::Filter(f,o).process(2,in.views),o);
      }
  }
}
void external_output_tests() {
  int variant=0;
  for(int depth : {8,10,16,32,-16}) for(int mode=1;mode<=7;++mode) {
    ref::Format f;f.width=variant%2?194:70;f.height=38;f.frames=3;
    f.depth=std::abs(depth);f.sample=depth==32 || depth<0?ref::Sample::floating:ref::Sample::integer;
    f.color=ref::Color::yuv;f.sub_x=variant&1;f.sub_y=(variant>>1)&1;
    ref::Options o;o.sample_mode=mode;o.dither_algo=variant%3+1;
    o.dynamic_grain=true;o.grain_y=3.25;o.grain_c=2.5;
    o.blur_first=(variant&1)!=0;o.keep_tv_range=(variant&2)!=0;
    o.output_sample_type=variant%4==0?ref::OutputSample::floating:ref::OutputSample::integer;
    o.output_depth=variant%4==0?(variant%8?16:32):(variant%3==0?8:variant%3==1?10:16);
    const auto in=input(f,variant%4);const simd::Filter filter(f,o);
    const auto& out=filter.parameters().output;const auto b=out.depth==8?1:out.depth<=16?2:4;
    std::vector<std::vector<std::uint8_t>> buffers(3);
    std::vector<ref::MutablePlaneView> views;
    for(std::size_t plane=0;plane<3;++plane) {
      const auto w=f.width>>(plane?f.sub_x:0),h=f.height>>(plane?f.sub_y:0),row=w*b,stride=row+9;
      const auto extent=(h-1)*stride+row;
      buffers[plane].assign(extent+2,0xa5);
      views.push_back({buffers[plane].data()+1,extent,stride});
    }
    for(std::size_t frame : {std::size_t{2},std::size_t{0}}) {
      filter.process_into(frame,in.views,views);
      ref::Frame actual{out,{}};
      for(std::size_t plane=0;plane<3;++plane) {
        const auto w=f.width>>(plane?f.sub_x:0),h=f.height>>(plane?f.sub_y:0),row=w*b,stride=views[plane].stride;
        ref::Plane image{w,h,row,{}};image.data.resize(row*h);
        const auto* data=static_cast<const std::uint8_t*>(views[plane].data);
        for(std::size_t y=0;y<h;++y) {
          std::memcpy(image.data.data()+y*row,data+y*stride,row);
          if(y+1<h) for(std::size_t x=row;x<stride;++x) check(data[y*stride+x]==0xa5,"external output padding");
        }
        check(buffers[plane].front()==0xa5 && buffers[plane].back()==0xa5,"external output guards");
        actual.planes.push_back(std::move(image));
      }
      compare(ScalarFilter(f,o).process(frame,in.views),actual,o);
    }
    ++variant;
  }
  ref::Format f;f.width=67;f.height=5;f.color=ref::Color::yuv;f.frames=2;
  ref::Options o;o.dither_algo=1;const simd::Filter filter(f,o);const auto in=input(f,1);
  std::vector<std::vector<std::uint8_t>> buffers(3,std::vector<std::uint8_t>(67*5,0xa5));
  std::vector<ref::MutablePlaneView> views;
  for(auto& bytes:buffers) views.push_back({bytes.data(),bytes.size(),67});
  auto reject=[&](const std::vector<ref::MutablePlaneView>& output,const char* reason) {
    bool rejected=false;
    try { filter.process_into(0,in.views,output); }
    catch(const std::invalid_argument& e) { rejected=std::string(e.what()).find(reason)!=std::string::npos; }
    check(rejected,"external output rejects invalid storage/overlap");
    for(const auto& bytes:buffers) check(std::all_of(bytes.begin(),bytes.end(),[](auto b){return b==0xa5;}),"external output validates before writing");
  };
  auto bad=views;bad[0].stride=66;reject(bad,"storage");
  bad=views;--bad[2].bytes;reject(bad,"storage");
  bad=views;bad.pop_back();reject(bad,"plane count");
  bad=views;bad[0].data=nullptr;reject(bad,"storage");
  bad=views;bad[1]=bad[0];reject(bad,"overlap");
  bad=views;bad[0]={const_cast<void*>(in.views[0].data),in.views[0].bytes,in.views[0].stride};reject(bad,"overlap");
  auto bad_input=in.views;bad_input[2].data=nullptr;
  bool rejected=false;
  try { filter.process_into(0,bad_input,views); } catch(const std::invalid_argument&) { rejected=true; }
  check(rejected,"external output rejects invalid final input plane");
  rejected=false;
  try { filter.process_into(f.frames,in.views,views); } catch(const std::invalid_argument&) { rejected=true; }
  check(rejected,"external output rejects invalid frame");
  for(const auto& bytes:buffers) check(std::all_of(bytes.begin(),bytes.end(),[](auto b){return b==0xa5;}),"invalid input leaves external output intact");

}
void byte_minimum_extent_tests() {
  // The final gathered byte may be the final allocated byte. No trailing
  // padding is available for a four-byte load, even for unaligned input.
  for(int mode : {1,2,3,4}) for(std::size_t w : {1,2,3,7,16,17,31,32,33,65}) for(std::size_t h : {1,3,17})
    for(bool blur : {false,true}) {
      ref::Format f;f.width=w;f.height=h;f.frames=1;f.color=ref::Color::gray;
      ref::Options o;o.sample_mode=mode;o.dither_algo=1;o.blur_first=blur;
      o.range=255;
      o.threshold.fill(100);o.grain_y=3.75;
      const auto stride=w+3,extent=(h-1)*stride+w;
      std::vector<std::uint8_t> data(extent+1,0xa5);
      for(std::size_t y=0;y<h;++y) for(std::size_t x=0;x<w;++x)
        data[1+y*stride+x]=std::uint8_t((x*73+y*101)%256);
      const std::vector<ref::PlaneView> views{{data.data()+1,extent,stride}};
      compare(ScalarFilter(f,o).process(0,views),simd::Filter(f,o).process(0,views),o);
    }
}
void native_integer_frame_tests() {
  // Exact allocation extent, unaligned base and odd stride exercise the final
  // four-byte gather without granting readable padding beyond the last pixel.
  for (int mode=1;mode<=7;++mode) for (int depth : {8,10,16})
    for (std::size_t w : {1,2,3,7,17,33}) for (int dither : {1,2,3}) for (bool float_output : {false,true}) {
      const std::size_t h=35, bytes=depth==8?1:2, stride=w*bytes+3, extent=(h-1)*stride+w*bytes;
      ref::Format f;f.width=w;f.height=h;f.depth=depth;
      ref::Options o;o.sample_mode=mode;o.dither_algo=dither;o.range=255;o.scale=true;
      o.threshold.fill(2048.5);o.threshold1.fill(4096);o.threshold2.fill(8192);o.grain_y=3.75;
      o.blur_first=(w%2)!=0;o.keep_tv_range=dither==2;
      o.output_sample_type=float_output?ref::OutputSample::floating:ref::OutputSample::integer;
      o.output_depth=float_output?32:dither==2?10:depth==8?16:8;
      std::vector<std::uint8_t> data(extent+1,0xa5);
      for (std::size_t y=0;y<h;++y) for (std::size_t x=0;x<w;++x) {
        const auto code=std::uint16_t((x*73+y*101)&((1u<<depth)-1));
        auto* p=data.data()+1+y*stride+x*bytes;
        if(bytes==1) *p=std::uint8_t(code);else write(p,code);
      }
      const std::vector<ref::PlaneView> views{{data.data()+1,extent,stride}};
      compare(ScalarFilter(f,o).process(0,views),simd::Filter(f,o).process(0,views),o);
      if(depth==10) {
        write(data.data()+data.size()-2,std::uint16_t(1024));
        bool rejected=false;
        try {simd::Filter(f,o).process(0,views);} catch(const std::invalid_argument&) {rejected=true;}
        check(rejected,"native input validates last out-of-depth code");
      }
    }
}
void cached_reference_frame_tests() {
  for(int depth : {8,10,16}) for(int sx : {0,1}) for(int sy : {0,1}) for(int dither : {1,2,3}) {
    ref::Format f;f.width=66;f.height=38;f.frames=2;f.depth=depth;
    f.color=ref::Color::yuv;f.sub_x=sx;f.sub_y=sy;
    ref::Options o;o.sample_mode=6;o.scale=true;o.range=255;o.dynamic_grain=true;
    o.threshold.fill(2047);o.threshold1.fill(4093);o.threshold2.fill(8191);
    o.dither_algo=dither;o.output_depth=dither==2?8:depth;o.grain_y=o.grain_c=3.75;
    const auto padded=input(f,0);
    const simd::Filter filter(f,o);
    std::vector<std::vector<std::uint8_t>> storage(3);
    std::vector<ref::PlaneView> tight;
    for(std::size_t plane=0;plane<3;++plane) {
      const auto w=f.width>>(plane?sx:0),h=f.height>>(plane?sy:0),stride=w*(depth==8?1:2);
      check(filter.reference_offsets(plane).size()==2*w*h,"native reference cache extent");
      storage[plane].resize(h*stride+1);
      for(std::size_t y=0;y<h;++y)
        std::memcpy(storage[plane].data()+1+y*stride,
                    static_cast<const std::uint8_t*>(padded.views[plane].data)+y*padded.views[plane].stride,stride);
      tight.push_back({storage[plane].data()+1,h*stride,stride});
    }
    const auto cached=filter.process(1,tight),fallback=filter.process(1,padded.views);
    compare(ScalarFilter(f,o).process(1,tight),cached,o);
    for(std::size_t plane=0;plane<3;++plane)
      check(cached.planes[plane].data==fallback.planes[plane].data,"reference cache matches stride fallback");
    o.threshold.fill(2047.5);
    const simd::Filter fractional(f,o);
    for(std::size_t plane=0;plane<3;++plane)
      check(fractional.reference_offsets(plane).empty(),"F64 path does not allocate unused offsets");
  }
}
void native_pair_frame_tests() {
  // Exercise native word pairs independently of the existing byte-output path:
  // both blur rules, subsampling, saturation, threshold fallback, and exact
  // unaligned allocations. Width 130 covers 32-lane table interiors and tails.
  int variant = 0;
  for (int mode : {1, 3})
    for (int depth : {8, 10, 16})
      for (int sx : {0, 1})
        for (int sy : {0, 1})
          for (bool blur : {false, true}) {
            ref::Format f;
            f.width = 130;
            f.height = 38;
            f.depth = depth;
            f.color = ref::Color::yuv;
            f.sub_x = sx;
            f.sub_y = sy;
            ref::Options o;
            o.sample_mode = mode;
            o.blur_first = blur;
            o.scale = true;
            o.dither_algo = 1;
            o.output_depth = depth == 8 ? 16 : (variant & 1) ? 8 : 16;
            const double thresholds[] = {0.5, 2048.5, 65535.0, 0.0};
            o.threshold.fill(thresholds[variant % 4]);
            o.range = variant % 4 == 3 ? 255 : 15;
            o.grain_y = o.grain_c = 4096;
            o.keep_tv_range = (variant & 2) != 0;
            const auto generated = input(f, 0);
            std::vector<std::vector<std::uint8_t>> storage(3);
            std::vector<ref::PlaneView> views;
            for (std::size_t plane = 0; plane < 3; ++plane) {
              const auto w = f.width >> (plane ? sx : 0), h = f.height >> (plane ? sy : 0);
              const auto row_bytes = w * (depth == 8 ? 1 : 2), stride = row_bytes + variant % 2;
              const auto extent = (h - 1) * stride + row_bytes;
              storage[plane].resize(extent + 1);
              for (std::size_t y = 0; y < h; ++y)
                std::memcpy(storage[plane].data() + 1 + y * stride,
                            static_cast<const std::uint8_t*>(generated.views[plane].data) +
                                y * generated.views[plane].stride,
                            row_bytes);
              views.push_back({storage[plane].data() + 1, extent, stride});
            }
            compare(ScalarFilter(f, o).process(0, views), simd::Filter(f, o).process(0, views), o);
            ++variant;
          }
}
void native_weighted_gather_tests() {
  // Maximal valid distances at each coordinate include odd luma distances
  // at the right edge of subsampled planes. Use exact, unaligned allocations
  // to check both direct gathers and the final-row bounded fallback.
  constexpr std::size_t height=9;
  for(int depth : {8,10,16}) for(int sx : {0,1}) for(int sy : {0,1})
    for(std::size_t width : {32,33,64,65}) for(std::size_t padding : {0,1,3}) {
      const std::size_t bytes=depth==8?1:2,stride=width*bytes+padding,extent=(height-1)*stride+width*bytes;
      std::vector<std::uint8_t> data(extent+1,0xa5);
      std::vector<std::int32_t> source(width*height);
      for(std::size_t y=0;y<height;++y) for(std::size_t x=0;x<width;++x) {
        const auto code=std::uint16_t((x*37+y*113)&((1u<<depth)-1));
        auto* p=data.data()+1+y*stride+x*bytes;
        if(bytes==1)*p=std::uint8_t(code);else write(p,code);
        source[y*width+x]=std::int32_t(code)<<(16-depth);
      }
      std::vector<ref::Distance> distances(width<<sx);
      std::vector<std::int16_t> grain(width,0);
      std::vector<double> output(width+1);
      const std::array<double,3> thresholds{2047,4093,8191};
      std::vector<std::int32_t> vertical_offsets(width),horizontal_offsets(width);
      for(std::size_t y=0;y<height;++y) {
        const auto vlimit=std::min(y,height-1-y);
        for(std::size_t x=0;x<width;++x) {
          const auto hlimit=std::min(x,width-1-x);
          distances[x<<sx].first=std::uint8_t(std::min({std::size_t{255},((vlimit+1)<<sy)-1,((hlimit+1)<<sx)-1}));
          vertical_offsets[x]=std::int32_t((distances[x<<sx].first>>sy)*stride);
          horizontal_offsets[x]=std::int32_t((distances[x<<sx].first>>sx)*bytes);
        }
        for(int variant : {0,1,2}) {
        const bool validated=variant!=0;
        simd::IntegerRow row{nullptr,distances.data(),grain.data(),output.data(),width,vlimit,y*stride,sx,sy};
        row.validated_distances=validated;
        if(variant==2) {
          row.vertical_offsets=vertical_offsets.data();
          row.horizontal_offsets=horizontal_offsets.data();
        }
        row.native=data.data()+1;row.native_stride=stride;row.native_extent=extent;row.input_depth=depth;
        std::fill(output.begin(),output.end(),-999999);
        const auto count=simd::process_integer_row(row,6,thresholds.data(),true);
        for(std::size_t x=0;x<count;++x) {
          const auto distance=distances[x<<sx].first;
          const std::size_t dx=distance>>sx,dy=distance>>sy;
          const std::array<double,4> refs{double(source[(y+dy)*width+x]),double(source[(y-dy)*width+x]),
              double(source[y*width+x+dx]),double(source[y*width+x-dx])};
          const double expected=ref::deband(source[y*width+x],refs,6,thresholds,true,true);
          check(std::abs(output[x]-expected)<=1,"native weighted gather oracle");
        }
        for(std::size_t x=count;x<output.size();++x)check(output[x]==-999999,"native weighted gather tail");
        if(!validated && y==0 && count) {
          const auto original=distances[0].first;
          distances[0].first=255;
          bool rejected=false;
          try { simd::process_integer_row(row,6,thresholds.data(),true); }
          catch(const std::logic_error&) { rejected=true; }
          distances[0].first=original;
          check(rejected,"untrusted native weighted coordinates rejected");
        }
        }
      }
    }
}
void byte_frame_tests() {
  // Full byte vectors, odd tails, chroma layouts, large offsets, fractional
  // thresholds, signed grain saturation, and limited-range packing.
  int variant=0;
  for(int mode=1;mode<=4;++mode) for(int sx : {0,1}) for(int sy : {0,1})
    for(bool blur : {false,true}) {
      const int pattern=variant%4;
      const bool wide=(mode==2 || mode==3) && sx==0 && sy==0 && blur;
      ref::Format f;f.width=wide?518:130;f.height=wide?518:38;
      f.frames=3;f.color=ref::Color::yuv;f.sub_x=sx;f.sub_y=sy;
      ref::Options o;o.sample_mode=mode;o.blur_first=blur;o.dither_algo=1;o.scale=true;
      const double thresholds[]={0.0,0.5,257.25,65535.0};
      o.threshold.fill(thresholds[pattern]);o.range=wide?255:15;
      o.grain_y=pattern==3?4096:3.75;o.grain_c=pattern==3?4096:2.5;
      o.keep_tv_range=(variant&1)!=0;o.dynamic_grain=true;
      const auto in=input(f,pattern);const simd::Filter filter(f,o);
      compare(ScalarFilter(f,o).process(2,in.views),filter.process(2,in.views),o);
      std::vector<std::vector<std::uint8_t>> buffers(3);
      std::vector<ref::MutablePlaneView> views;
      for(std::size_t plane=0;plane<3;++plane) {
        const auto w=f.width>>(plane?sx:0),h=f.height>>(plane?sy:0),stride=w+7;
        const auto extent=(h-1)*stride+w;buffers[plane].assign(extent+2,0xa5);
        views.push_back({buffers[plane].data()+1,extent,stride});
      }
      filter.process_into(0,in.views,views);
      const auto expected=ScalarFilter(f,o).process(0,in.views);
      for(std::size_t plane=0;plane<3;++plane) {
        const auto& e=expected.planes[plane];const auto* data=static_cast<const std::uint8_t*>(views[plane].data);
        for(std::size_t y=0;y<e.height;++y) {
          check(std::memcmp(data+y*views[plane].stride,e.data.data()+y*e.stride,e.width)==0,"byte external exact pixels");
          if(y+1<e.height) for(std::size_t x=e.width;x<views[plane].stride;++x) check(data[y*views[plane].stride+x]==0xa5,"byte output padding");
        }
        check(buffers[plane].front()==0xa5 && buffers[plane].back()==0xa5,"byte output guards");
      }
      ++variant;
    }
}
void native_float_frame_tests() {
  // Aligned F32 sources exercise the borrowed-plane path. Existing integration
  // fixtures deliberately use odd addresses, so they cover its fallback only.
  for (int mode = 1; mode <= 7; ++mode)
    for (bool blur : {false, true})
      for (int variant = 0; variant < 7; ++variant) {
        ref::Format f;
        f.width = 130;
        f.height = 38;
        f.frames = 3;
        f.depth = 32;
        f.sample = ref::Sample::floating;
        f.color = ref::Color::yuv;
        f.sub_x = variant & 1;
        f.sub_y = (variant >> 1) & 1;
        ref::Options o;
        o.sample_mode = mode;
        o.blur_first = blur;
        o.range = 255;
        o.scale = true;
        o.threshold.fill(variant == 3   ? 0
                         : variant == 4 ? std::nextafter(2048.0, 0.0)
                         : variant == 5 ? std::nextafter(2048.0, 4096.0)
                                        : 2048.25);
        o.threshold1.fill(1234.5);
        o.threshold2.fill(3456.75);
        const double angles[] = {0.15, 0, 0.25, 0.5, std::nextafter(0.25, 0.0), std::nextafter(0.25, 1.0), 0.15};
        o.max_angle = angles[variant];
        if (variant == 6) {
          o.threshold.fill(65535);
          o.threshold1.fill(65535);
          o.threshold2.fill(65535);
        }
        o.grain_y = variant == 2 ? 0 : 3.25;
        o.grain_c = variant == 2 ? 0 : 2.5;
        o.dynamic_grain = true;
        o.keep_tv_range = variant == 3;
        std::array<std::vector<float>, 3> storage;
        std::vector<ref::PlaneView> views;
        for (std::size_t p = 0; p < 3; ++p) {
          const auto w = f.width >> (p ? f.sub_x : 0), h = f.height >> (p ? f.sub_y : 0), stride = w + 3;
          storage[p].assign((h - 1) * stride + w, std::numeric_limits<float>::quiet_NaN());
          for (std::size_t y = 0; y < h; ++y)
            for (std::size_t x = 0; x < w; ++x) {
              const auto hash = std::uint32_t(x * 1664525u + y * 1013904223u);
              float value = float(hash >> 8) / 16777216.0f * 0.08f - 0.04f + (p ? 0.0f : 0.5f);
              if (variant >= 4)
                value = float((x * 13 + y * 7) % 3) / 32.0f;
              if (variant == 6)
                value = float((x * 3 + y * 7) % 5) * 0.5f;
              if (variant == 1 && (x + y) % 11 == 0)
                value = (x % 2 ? -1.0f : 1.0f) * std::numeric_limits<float>::max();
              if (variant == 2)
                value = float(int((x * 3 + y) % 17) - 8) * std::numeric_limits<float>::denorm_min();
              storage[p][y * stride + x] = value;
            }
          views.push_back({storage[p].data(), storage[p].size() * 4, stride * 4});
        }
        const simd::Filter filter(f, o);
        const auto expected = ScalarFilter(f, o).process(2, views);
        compare(expected, filter.process(2, views), o);
        if (mode == 7) {
          auto unit_boost = o;
          unit_boost.angle_boost = 1;
          compare(ScalarFilter(f, unit_boost).process(2, views), simd::Filter(f, unit_boost).process(2, views),
                  unit_boost);
        }
        std::array<std::vector<std::uint8_t>, 3> outputs;
        std::vector<ref::MutablePlaneView> targets;
        for (std::size_t p = 0; p < 3; ++p) {
          const auto w = expected.planes[p].width, h = expected.planes[p].height, stride = w * 4 + 7,
                     extent = (h - 1) * stride + w * 4;
          outputs[p].assign(extent + 2, 0xa5);
          targets.push_back({outputs[p].data() + 1, extent, stride});
        }
        filter.process_into(2, views, targets);
        for (std::size_t p = 0; p < 3; ++p) {
          const auto& e = expected.planes[p];
          auto actual = e;
          for (std::size_t y = 0; y < e.height; ++y) {
            std::memcpy(actual.data.data() + y * e.stride, outputs[p].data() + 1 + y * targets[p].stride, e.stride);
            if (y + 1 < e.height)
              for (std::size_t x = e.stride; x < targets[p].stride; ++x)
                check(outputs[p][1 + y * targets[p].stride + x] == 0xa5, "native float output padding");
          }
          compare(ref::Frame{expected.format, {e}}, ref::Frame{expected.format, {actual}}, o);
          check(outputs[p].front() == 0xa5 && outputs[p].back() == 0xa5, "native float output guards");
        }
        if (variant == 0 || variant == 4 || variant == 5) {
          for (int output_depth : {8, 16})
            for (int dither : {1, 2, 3}) {
              auto converted = o;
              converted.output_sample_type = ref::OutputSample::integer;
              converted.output_depth = output_depth;
              converted.dither_algo = dither;
              compare(ScalarFilter(f, converted).process(2, views), simd::Filter(f, converted).process(2, views),
                      converted);
            }
          auto half = o;
          half.output_sample_type = ref::OutputSample::floating;
          half.output_depth = 16;
          compare(ScalarFilter(f, half).process(2, views), simd::Filter(f, half).process(2, views), half);
        }
        storage[2].back() = std::numeric_limits<float>::infinity();
        bool rejected = false;
        try {
          filter.process(0, views);
        } catch (const std::invalid_argument&) {
          rejected = true;
        }
        check(rejected, "native float nonfinite tail rejection");
      }
}
void native_half_frame_tests() {
  for (int mode = 1; mode <= 7; ++mode)
    for (bool blur : {false, true})
      for (int variant = 0; variant < 5; ++variant) {
        ref::Format f;
        f.width = variant == 0 ? 66 : 34;
        f.height = 34;
        f.frames = 3;
        f.depth = 16;
        f.sample = ref::Sample::floating;
        f.color = ref::Color::yuv;
        f.sub_x = variant & 1;
        f.sub_y = (variant >> 1) & 1;
        ref::Options o;
        o.sample_mode = mode;
        o.blur_first = blur;
        o.range = 255;
        o.scale = true;
        o.threshold.fill(variant == 3 ? 65535 : 2048.25);
        o.threshold1.fill(1234.5);
        o.threshold2.fill(3456.75);
        o.max_angle = variant == 3 ? 0.25 : 0.15;
        o.grain_y = 3.25;
        o.grain_c = 2.5;
        o.dynamic_grain = true;
        o.keep_tv_range = variant == 4;
        std::array<std::vector<std::uint8_t>, 3> storage;
        std::vector<ref::PlaneView> views;
        for (std::size_t p = 0; p < 3; ++p) {
          const auto w = f.width >> (p ? f.sub_x : 0), h = f.height >> (p ? f.sub_y : 0);
          const auto stride = w * 2 + 3, extent = (h - 1) * stride + w * 2;
          storage[p].assign(extent + 2, 0xff);
          for (std::size_t y = 0; y < h; ++y)
            for (std::size_t x = 0; x < w; ++x) {
              const auto hash = std::uint32_t(x * 1664525u + y * 1013904223u);
              double v = double(hash >> 8) / 16777216 * .08 - .04 + (p ? 0 : .5);
              if (variant == 1)
                v = ((x + y) % 2 ? -1 : 1) * 65504.0;
              if (variant == 2)
                v = (int((x * 3 + y) % 17) - 8) * 0x1p-24;
              if (variant == 3)
                v = ((x * 3 + y * 7) % 5) * .5;
              write(storage[p].data() + 1 + y * stride + 2 * x, ref::encode_half(v));
            }
          views.push_back({storage[p].data() + 1, extent, stride});
        }
        for (int output : {-16, 32, 8, 16}) {
          const int selected[] = {-16, 32, 16, -16, 8};
          if (variant != 0 && output != selected[variant])
            continue;
          o.output_depth = std::abs(output);
          o.output_sample_type = output < 0 || output == 32 ? ref::OutputSample::floating : ref::OutputSample::integer;
          o.dither_algo = variant % 3 + 1;
          compare(ScalarFilter(f, o).process(2, views), simd::Filter(f, o).process(2, views), o);
        }
        if (mode == 1 && !blur && variant == 0)
          for (std::uint16_t invalid :
               {std::uint16_t{0x7c00}, std::uint16_t{0xfc00}, std::uint16_t{0x7c01}, std::uint16_t{0x7e00}}) {
            write(storage[2].data() + storage[2].size() - 3, invalid);
            bool rejected = false;
            try {
              simd::Filter(f, o).process(0, views);
            } catch (const std::invalid_argument&) {
              rejected = true;
            }
            check(rejected, "native half nonfinite tail rejection");
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
  // Every luma and chroma row has a scalar tail after fused gradients.
  for (int depth : {16,32,-16}) for (double angle : {0.0,0.15,0.25}) {
    ref::Format f;f.width=38;f.height=46;f.frames=3;f.depth=std::abs(depth);
    f.sample=depth==16?ref::Sample::integer:ref::Sample::floating;
    f.color=ref::Color::yuv;f.sub_x=1;f.sub_y=1;
    ref::Options o;o.sample_mode=7;o.scale=true;o.threshold.fill(65535);o.max_angle=angle;o.angle_boost=2.5;
    o.grain_y=3.25;o.grain_c=2.5;o.dynamic_grain=true;o.dither_algo=1;o.range=255;
    const auto in=input(f,0);
    compare(ScalarFilter(f,o).process(2,in.views),simd::Filter(f,o).process(2,in.views),o);
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
      if (target != HWY_SSSE3 && target != HWY_SSE4 && target != HWY_AVX2 && target != HWY_AVX3 && target != HWY_AVX3_SPR)
        continue;
#endif
      hwy::SetSupportedTargetsForTest(target);
      std::cout<<"Testing "<<hwy::TargetName(target)<<std::endl;
      weight_tests();
      normal_weight_tests();
      sampling_tests();
      gradient_sampling_tests();
      integer_cross_tests();
      integer_weight_table_tests();
      integer_weight_boundary_tests();
      horizontal_table_tests();
      integer_pack_tests();
      row_tests();
      fs_tests();
      identity_frame_tests();
      external_output_tests();
      byte_minimum_extent_tests();
      native_integer_frame_tests();
      native_pair_frame_tests();
      cached_reference_frame_tests();
      native_weighted_gather_tests();
      byte_frame_tests();
      native_float_frame_tests();
      native_half_frame_tests();
      frame_tests();
      std::cout<<hwy::TargetName(target)<<": checks="<<checks<<" pixels="<<pixels<<" differing="<<different
               <<" fs_pixels="<<fs_pixels<<" max_weight="<<max_weight<<" max_integer_internal="<<max_integer<<" max_float_internal="<<max_float<<" max_pixel="<<max_pixel<<'\n';
    }
    if (!checks) throw std::runtime_error("No supported F3kdb test target executed");
    hwy::SetSupportedTargetsForTest(0);
  } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
