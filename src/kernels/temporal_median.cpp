#include "kernels/dispatch.hpp"
#include "common/copy.hpp"
#include "common/padded_row.hpp"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "kernels/temporal_median.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace neo_smo {
namespace HWY_NAMESPACE {

#include "common/sorting_networks.hpp"
#include "common/float_arithmetic.hpp"

template <class D, class V = hn::Vec<D>>
HWY_INLINE V avg2_int(D d, V a, V b) {
  const auto one = hn::Set(d, 1);
  return hn::Add(hn::Add(hn::ShiftRight<1>(a), hn::ShiftRight<1>(b)), hn::And(hn::And(a, b), one));
}

template <bool IsF16, class D, class V = hn::Vec<D>>
HWY_INLINE V avg2_float(D d, V a, V b) {
  const auto half = float_set(d, 0.5f);
  return float_mul<IsF16>(d, float_add<IsF16>(d, a, b), half);
}

template <class D, class V = hn::Vec<D>>
HWY_INLINE V median21(D d, V* v) {
  #define CS(i, j) compare_swap(d, v[i], v[j])
CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13); CS(14, 15); CS(16, 17); CS(18, 19);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7); CS(8, 10); CS(9, 11); CS(12, 14); CS(13, 15); CS(16, 18); CS(17, 19);
      CS(1, 5); CS(2, 6); CS(3, 15); CS(4, 16); CS(13, 17); CS(14, 18);
      CS(1, 14); CS(2, 13); CS(3, 7); CS(5, 18); CS(6, 17); CS(12, 16);
      CS(0, 16); CS(1, 2); CS(3, 19); CS(5, 13); CS(6, 14); CS(17, 18);
      CS(0, 4); CS(5, 14); CS(6, 10); CS(9, 13); CS(15, 19);
      CS(5, 8); CS(6, 12); CS(7, 13); CS(11, 14);
      CS(2, 12); CS(7, 17); CS(8, 9); CS(10, 11);
      CS(3, 9); CS(7, 11); CS(8, 12); CS(10, 16);
      CS(3, 10); CS(4, 12); CS(7, 15); CS(9, 16);
      CS(7, 10); CS(9, 12);
      CS(7, 9); CS(10, 12);
      CS(9, 10);
      CS(10, 20);
      CS(9, 10);
      return v[10];
  #undef CS
}

template <class D, class V = hn::Vec<D>>
HWY_NOINLINE V eval_temporal_median_int(D d, V* v, int diameter) {
  #define CS(i, j) compare_swap(d, v[i], v[j])
  switch (diameter) {
    case 1: return v[0];
    case 2: return avg2_int(d, v[0], v[1]);
    case 3: return median3(d, v[0], v[1], v[2]);
    case 4:
      CS(0, 1); CS(2, 3);
      CS(0, 2); CS(1, 3);
      return avg2_int(d, v[1], v[2]);
    case 5:
      CS(0, 1); CS(2, 3);
      CS(0, 2); CS(1, 3);
      CS(2, 4); CS(1, 2); CS(2, 4);
      return v[2];
    case 6:
      CS(0, 1); CS(4, 5);
      CS(0, 5); CS(1, 3); CS(2, 4);
      CS(0, 2); CS(1, 4); CS(3, 5);
      CS(1, 2); CS(3, 4);
      return avg2_int(d, v[2], v[3]);
    case 7:
      CS(0, 6); CS(1, 2); CS(3, 4);
      CS(0, 2); CS(1, 4); CS(3, 5);
      CS(0, 1); CS(2, 5); CS(4, 6);
      CS(1, 3); CS(2, 4);
      CS(3, 4); CS(2, 3);
      return v[3];
    case 8:
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7);
      CS(0, 4); CS(1, 5); CS(2, 6); CS(3, 7);
      CS(0, 1); CS(2, 4); CS(3, 5); CS(6, 7);
      CS(2, 3); CS(4, 5);
      CS(1, 4); CS(3, 6);
      return avg2_int(d, v[3], v[4]);
    case 9:
      return median9(d, v);
    case 10:
      CS(0, 1); CS(3, 5); CS(4, 6); CS(8, 9);
      CS(0, 3); CS(1, 5); CS(4, 8); CS(6, 9);
      CS(1, 3); CS(6, 8);
      CS(0, 6); CS(1, 4); CS(3, 9); CS(5, 8);
      CS(2, 6); CS(3, 7);
      CS(2, 3); CS(6, 7);
      CS(3, 4); CS(5, 6);
      CS(3, 5); CS(4, 6);
      return avg2_int(d, v[4], v[5]);
    case 11:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7);
      CS(1, 2); CS(5, 6);
      CS(0, 5); CS(1, 4); CS(2, 7); CS(3, 6);
      CS(2, 5); CS(3, 4);
      CS(2, 8); CS(4, 9);
      CS(8, 10);
      CS(3, 8); CS(5, 10);
      CS(4, 5);
      CS(4, 8);
      CS(5, 8);
      return v[5];
    case 12:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7); CS(8, 10); CS(9, 11);
      CS(0, 4); CS(1, 10); CS(2, 9); CS(5, 6); CS(7, 11);
      CS(2, 6); CS(3, 7); CS(4, 8); CS(5, 9);
      CS(1, 5); CS(2, 8); CS(3, 9); CS(6, 10);
      CS(3, 5); CS(6, 8);
      CS(3, 6); CS(5, 8);
      return avg2_int(d, v[5], v[6]);
    case 13:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11);
      CS(0, 8); CS(1, 9); CS(2, 4); CS(3, 5); CS(6, 10); CS(7, 11);
      CS(0, 6); CS(1, 7); CS(3, 4); CS(8, 10); CS(9, 11);
      CS(2, 6); CS(3, 10); CS(4, 8); CS(5, 9);
      CS(1, 3); CS(4, 6); CS(7, 8);
      CS(3, 12); CS(5, 6);
      CS(3, 5); CS(6, 12);
      CS(5, 7); CS(6, 10);
      CS(5, 6);
      CS(6, 7);
      return v[6];
    case 14:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13);
      CS(0, 9); CS(2, 10); CS(3, 11); CS(4, 13);
      CS(0, 12); CS(1, 13); CS(4, 6); CS(7, 9);
      CS(0, 4); CS(1, 10); CS(3, 12); CS(5, 7); CS(6, 8); CS(9, 13);
      CS(1, 8); CS(2, 4); CS(5, 12); CS(9, 11);
      CS(1, 5); CS(3, 4); CS(8, 12); CS(9, 10);
      CS(3, 5); CS(4, 9); CS(8, 10);
      CS(4, 6); CS(7, 9);
      CS(5, 7); CS(6, 8);
      CS(5, 6); CS(7, 8);
      return avg2_int(d, v[6], v[7]);
    case 15:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13);
      CS(0, 12); CS(1, 13); CS(2, 8); CS(3, 9); CS(10, 14);
      CS(3, 11); CS(10, 12);
      CS(2, 10); CS(3, 6); CS(5, 11); CS(7, 12);
      CS(0, 3); CS(4, 10); CS(5, 6); CS(7, 8); CS(9, 12);
      CS(1, 5); CS(3, 10); CS(6, 13); CS(8, 11);
      CS(1, 3); CS(5, 8); CS(6, 10);
      CS(3, 7); CS(5, 14); CS(6, 9);
      CS(4, 7); CS(5, 6); CS(8, 14); CS(9, 10);
      CS(6, 7); CS(8, 9);
      CS(6, 8);
      CS(7, 8);
      return v[7];
    case 16:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13); CS(14, 15);
      CS(0, 6); CS(2, 4); CS(9, 15); CS(11, 13);
      CS(4, 9); CS(6, 11);
      CS(1, 9); CS(3, 11); CS(4, 12); CS(6, 14);
      CS(0, 4); CS(1, 10); CS(2, 6); CS(3, 8); CS(5, 14); CS(7, 12); CS(9, 13); CS(11, 15);
      CS(1, 12); CS(3, 14); CS(4, 7); CS(5, 6); CS(8, 11); CS(9, 10);
      CS(1, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 14);
      CS(3, 5); CS(6, 8); CS(7, 9); CS(10, 12);
      CS(5, 8); CS(7, 10);
      CS(5, 7); CS(8, 10);
      return avg2_int(d, v[7], v[8]);
    case 17:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13); CS(14, 15);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7); CS(8, 10); CS(9, 11); CS(12, 14); CS(13, 15);
      CS(1, 2); CS(5, 13); CS(6, 14); CS(9, 10);
      CS(1, 9); CS(2, 10); CS(5, 6); CS(13, 14);
      CS(1, 4); CS(2, 13); CS(5, 8); CS(6, 9); CS(7, 10); CS(11, 14);
      CS(2, 6); CS(4, 7); CS(8, 11); CS(9, 13);
      CS(0, 8); CS(3, 11); CS(4, 12); CS(7, 15); CS(9, 16);
      CS(3, 8); CS(6, 9); CS(7, 12);
      CS(3, 7); CS(8, 16);
      CS(7, 9); CS(8, 12);
      CS(6, 7); CS(8, 9);
      CS(7, 8);
      return v[8];
    case 18:
      CS(1, 2); CS(3, 4); CS(5, 6); CS(7, 8); CS(9, 10); CS(11, 12); CS(13, 14); CS(15, 16);
      CS(0, 2); CS(3, 5); CS(4, 6); CS(7, 9); CS(8, 10); CS(11, 13); CS(12, 14); CS(15, 17);
      CS(0, 4); CS(1, 5); CS(2, 14); CS(3, 15); CS(12, 16); CS(13, 17);
      CS(0, 13); CS(1, 12); CS(2, 6); CS(4, 17); CS(5, 16); CS(11, 15);
      CS(0, 1); CS(4, 12); CS(5, 13); CS(16, 17);
      CS(4, 13); CS(5, 9); CS(8, 12);
      CS(4, 7); CS(5, 11); CS(6, 12); CS(10, 13);
      CS(1, 11); CS(6, 16); CS(7, 8); CS(9, 10);
      CS(2, 8); CS(6, 10); CS(7, 11); CS(9, 15);
      CS(2, 9); CS(3, 11); CS(6, 14); CS(8, 15);
      CS(6, 9); CS(8, 11);
      CS(6, 8); CS(9, 11);
      return avg2_int(d, v[8], v[9]);
    case 19:
      CS(1, 2); CS(3, 4); CS(5, 6); CS(7, 8); CS(9, 10); CS(11, 12); CS(13, 14); CS(15, 16);
      CS(0, 2); CS(3, 5); CS(4, 6); CS(7, 9); CS(8, 10); CS(11, 13); CS(12, 14); CS(15, 17);
      CS(0, 4); CS(1, 5); CS(2, 14); CS(3, 15); CS(12, 16); CS(13, 17);
      CS(0, 13); CS(1, 12); CS(2, 6); CS(4, 17); CS(5, 16); CS(11, 15);
      CS(0, 1); CS(4, 12); CS(5, 13); CS(16, 17);
      CS(4, 13); CS(5, 9); CS(8, 12);
      CS(4, 7); CS(5, 11); CS(6, 12); CS(10, 13);
      CS(1, 11); CS(6, 16); CS(7, 8); CS(9, 10);
      CS(2, 8); CS(6, 10); CS(7, 11); CS(9, 15);
      CS(2, 9); CS(3, 11); CS(6, 14); CS(8, 15);
      CS(6, 9); CS(8, 11);
      CS(6, 8); CS(9, 11);
      CS(8, 9);
      CS(9, 18);
      CS(8, 9);
      return v[9];
    case 20:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13); CS(14, 15); CS(16, 17); CS(18, 19);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7); CS(8, 10); CS(9, 11); CS(12, 14); CS(13, 15); CS(16, 18); CS(17, 19);
      CS(1, 5); CS(2, 6); CS(3, 15); CS(4, 16); CS(13, 17); CS(14, 18);
      CS(1, 14); CS(2, 13); CS(3, 7); CS(5, 18); CS(6, 17); CS(12, 16);
      CS(0, 16); CS(1, 2); CS(3, 19); CS(5, 13); CS(6, 14); CS(17, 18);
      CS(0, 4); CS(5, 14); CS(6, 10); CS(9, 13); CS(15, 19);
      CS(5, 8); CS(6, 12); CS(7, 13); CS(11, 14);
      CS(2, 12); CS(7, 17); CS(8, 9); CS(10, 11);
      CS(3, 9); CS(7, 11); CS(8, 12); CS(10, 16);
      CS(3, 10); CS(4, 12); CS(7, 15); CS(9, 16);
      CS(7, 10); CS(9, 12);
      CS(7, 9); CS(10, 12);
      return avg2_int(d, v[9], v[10]);
    case 21: return median21(d, v);
    default:
      return v[0];
  }
  #undef CS
}

template <bool IsF16, class D, class V = hn::Vec<D>>
HWY_NOINLINE V eval_temporal_median_float(D d, V* v, int diameter) {
  #define CS(i, j) compare_swap(d, v[i], v[j])
  switch (diameter) {
    case 1: return v[0];
    case 2: return avg2_float<IsF16>(d, v[0], v[1]);
    case 3: return median3(d, v[0], v[1], v[2]);
    case 4:
      CS(0, 1); CS(2, 3);
      CS(0, 2); CS(1, 3);
      return avg2_float<IsF16>(d, v[1], v[2]);
    case 5:
      CS(0, 1); CS(2, 3);
      CS(0, 2); CS(1, 3);
      CS(2, 4); CS(1, 2); CS(2, 4);
      return v[2];
    case 6:
      CS(0, 1); CS(4, 5);
      CS(0, 5); CS(1, 3); CS(2, 4);
      CS(0, 2); CS(1, 4); CS(3, 5);
      CS(1, 2); CS(3, 4);
      return avg2_float<IsF16>(d, v[2], v[3]);
    case 7:
      CS(0, 6); CS(1, 2); CS(3, 4);
      CS(0, 2); CS(1, 4); CS(3, 5);
      CS(0, 1); CS(2, 5); CS(4, 6);
      CS(1, 3); CS(2, 4);
      CS(3, 4); CS(2, 3);
      return v[3];
    case 8:
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7);
      CS(0, 4); CS(1, 5); CS(2, 6); CS(3, 7);
      CS(0, 1); CS(2, 4); CS(3, 5); CS(6, 7);
      CS(2, 3); CS(4, 5);
      CS(1, 4); CS(3, 6);
      return avg2_float<IsF16>(d, v[3], v[4]);
    case 9:
      return median9(d, v);
    case 10:
      CS(0, 1); CS(3, 5); CS(4, 6); CS(8, 9);
      CS(0, 3); CS(1, 5); CS(4, 8); CS(6, 9);
      CS(1, 3); CS(6, 8);
      CS(0, 6); CS(1, 4); CS(3, 9); CS(5, 8);
      CS(2, 6); CS(3, 7);
      CS(2, 3); CS(6, 7);
      CS(3, 4); CS(5, 6);
      CS(3, 5); CS(4, 6);
      return avg2_float<IsF16>(d, v[4], v[5]);
    case 11:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7);
      CS(1, 2); CS(5, 6);
      CS(0, 5); CS(1, 4); CS(2, 7); CS(3, 6);
      CS(2, 5); CS(3, 4);
      CS(2, 8); CS(4, 9);
      CS(8, 10);
      CS(3, 8); CS(5, 10);
      CS(4, 5);
      CS(4, 8);
      CS(5, 8);
      return v[5];
    case 12:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7); CS(8, 10); CS(9, 11);
      CS(0, 4); CS(1, 10); CS(2, 9); CS(5, 6); CS(7, 11);
      CS(2, 6); CS(3, 7); CS(4, 8); CS(5, 9);
      CS(1, 5); CS(2, 8); CS(3, 9); CS(6, 10);
      CS(3, 5); CS(6, 8);
      CS(3, 6); CS(5, 8);
      return avg2_float<IsF16>(d, v[5], v[6]);
    case 13:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11);
      CS(0, 8); CS(1, 9); CS(2, 4); CS(3, 5); CS(6, 10); CS(7, 11);
      CS(0, 6); CS(1, 7); CS(3, 4); CS(8, 10); CS(9, 11);
      CS(2, 6); CS(3, 10); CS(4, 8); CS(5, 9);
      CS(1, 3); CS(4, 6); CS(7, 8);
      CS(3, 12); CS(5, 6);
      CS(3, 5); CS(6, 12);
      CS(5, 7); CS(6, 10);
      CS(5, 6);
      CS(6, 7);
      return v[6];
    case 14:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13);
      CS(0, 9); CS(2, 10); CS(3, 11); CS(4, 13);
      CS(0, 12); CS(1, 13); CS(4, 6); CS(7, 9);
      CS(0, 4); CS(1, 10); CS(3, 12); CS(5, 7); CS(6, 8); CS(9, 13);
      CS(1, 8); CS(2, 4); CS(5, 12); CS(9, 11);
      CS(1, 5); CS(3, 4); CS(8, 12); CS(9, 10);
      CS(3, 5); CS(4, 9); CS(8, 10);
      CS(4, 6); CS(7, 9);
      CS(5, 7); CS(6, 8);
      CS(5, 6); CS(7, 8);
      return avg2_float<IsF16>(d, v[6], v[7]);
    case 15:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13);
      CS(0, 12); CS(1, 13); CS(2, 8); CS(3, 9); CS(10, 14);
      CS(3, 11); CS(10, 12);
      CS(2, 10); CS(3, 6); CS(5, 11); CS(7, 12);
      CS(0, 3); CS(4, 10); CS(5, 6); CS(7, 8); CS(9, 12);
      CS(1, 5); CS(3, 10); CS(6, 13); CS(8, 11);
      CS(1, 3); CS(5, 8); CS(6, 10);
      CS(3, 7); CS(5, 14); CS(6, 9);
      CS(4, 7); CS(5, 6); CS(8, 14); CS(9, 10);
      CS(6, 7); CS(8, 9);
      CS(6, 8);
      CS(7, 8);
      return v[7];
    case 16:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13); CS(14, 15);
      CS(0, 6); CS(2, 4); CS(9, 15); CS(11, 13);
      CS(4, 9); CS(6, 11);
      CS(1, 9); CS(3, 11); CS(4, 12); CS(6, 14);
      CS(0, 4); CS(1, 10); CS(2, 6); CS(3, 8); CS(5, 14); CS(7, 12); CS(9, 13); CS(11, 15);
      CS(1, 12); CS(3, 14); CS(4, 7); CS(5, 6); CS(8, 11); CS(9, 10);
      CS(1, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 14);
      CS(3, 5); CS(6, 8); CS(7, 9); CS(10, 12);
      CS(5, 8); CS(7, 10);
      CS(5, 7); CS(8, 10);
      return avg2_float<IsF16>(d, v[7], v[8]);
    case 17:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13); CS(14, 15);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7); CS(8, 10); CS(9, 11); CS(12, 14); CS(13, 15);
      CS(1, 2); CS(5, 13); CS(6, 14); CS(9, 10);
      CS(1, 9); CS(2, 10); CS(5, 6); CS(13, 14);
      CS(1, 4); CS(2, 13); CS(5, 8); CS(6, 9); CS(7, 10); CS(11, 14);
      CS(2, 6); CS(4, 7); CS(8, 11); CS(9, 13);
      CS(0, 8); CS(3, 11); CS(4, 12); CS(7, 15); CS(9, 16);
      CS(3, 8); CS(6, 9); CS(7, 12);
      CS(3, 7); CS(8, 16);
      CS(7, 9); CS(8, 12);
      CS(6, 7); CS(8, 9);
      CS(7, 8);
      return v[8];
    case 18:
      CS(1, 2); CS(3, 4); CS(5, 6); CS(7, 8); CS(9, 10); CS(11, 12); CS(13, 14); CS(15, 16);
      CS(0, 2); CS(3, 5); CS(4, 6); CS(7, 9); CS(8, 10); CS(11, 13); CS(12, 14); CS(15, 17);
      CS(0, 4); CS(1, 5); CS(2, 14); CS(3, 15); CS(12, 16); CS(13, 17);
      CS(0, 13); CS(1, 12); CS(2, 6); CS(4, 17); CS(5, 16); CS(11, 15);
      CS(0, 1); CS(4, 12); CS(5, 13); CS(16, 17);
      CS(4, 13); CS(5, 9); CS(8, 12);
      CS(4, 7); CS(5, 11); CS(6, 12); CS(10, 13);
      CS(1, 11); CS(6, 16); CS(7, 8); CS(9, 10);
      CS(2, 8); CS(6, 10); CS(7, 11); CS(9, 15);
      CS(2, 9); CS(3, 11); CS(6, 14); CS(8, 15);
      CS(6, 9); CS(8, 11);
      CS(6, 8); CS(9, 11);
      return avg2_float<IsF16>(d, v[8], v[9]);
    case 19:
      CS(1, 2); CS(3, 4); CS(5, 6); CS(7, 8); CS(9, 10); CS(11, 12); CS(13, 14); CS(15, 16);
      CS(0, 2); CS(3, 5); CS(4, 6); CS(7, 9); CS(8, 10); CS(11, 13); CS(12, 14); CS(15, 17);
      CS(0, 4); CS(1, 5); CS(2, 14); CS(3, 15); CS(12, 16); CS(13, 17);
      CS(0, 13); CS(1, 12); CS(2, 6); CS(4, 17); CS(5, 16); CS(11, 15);
      CS(0, 1); CS(4, 12); CS(5, 13); CS(16, 17);
      CS(4, 13); CS(5, 9); CS(8, 12);
      CS(4, 7); CS(5, 11); CS(6, 12); CS(10, 13);
      CS(1, 11); CS(6, 16); CS(7, 8); CS(9, 10);
      CS(2, 8); CS(6, 10); CS(7, 11); CS(9, 15);
      CS(2, 9); CS(3, 11); CS(6, 14); CS(8, 15);
      CS(6, 9); CS(8, 11);
      CS(6, 8); CS(9, 11);
      CS(8, 9);
      CS(9, 18);
      CS(8, 9);
      return v[9];
    case 20:
      CS(0, 1); CS(2, 3); CS(4, 5); CS(6, 7); CS(8, 9); CS(10, 11); CS(12, 13); CS(14, 15); CS(16, 17); CS(18, 19);
      CS(0, 2); CS(1, 3); CS(4, 6); CS(5, 7); CS(8, 10); CS(9, 11); CS(12, 14); CS(13, 15); CS(16, 18); CS(17, 19);
      CS(1, 5); CS(2, 6); CS(3, 15); CS(4, 16); CS(13, 17); CS(14, 18);
      CS(1, 14); CS(2, 13); CS(3, 7); CS(5, 18); CS(6, 17); CS(12, 16);
      CS(0, 16); CS(1, 2); CS(3, 19); CS(5, 13); CS(6, 14); CS(17, 18);
      CS(0, 4); CS(5, 14); CS(6, 10); CS(9, 13); CS(15, 19);
      CS(5, 8); CS(6, 12); CS(7, 13); CS(11, 14);
      CS(2, 12); CS(7, 17); CS(8, 9); CS(10, 11);
      CS(3, 9); CS(7, 11); CS(8, 12); CS(10, 16);
      CS(3, 10); CS(4, 12); CS(7, 15); CS(9, 16);
      CS(7, 10); CS(9, 12);
      CS(7, 9); CS(10, 12);
      return avg2_float<IsF16>(d, v[9], v[10]);
    case 21: return median21(d, v);
    default:
      return v[0];
  }
  #undef CS
}

template <typename T>
void temporal_median_int_impl(int diameter, const T* const* srcp_planes, T* dstp,
                              int width, int height, std::size_t src_stride, std::size_t dst_stride) {
  const hn::ScalableTag<T> d;
  const std::size_t lanes = hn::Lanes(d);

  for (int y = 0; y < height; ++y) {
    T* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const std::size_t rem = static_cast<std::size_t>(width - x);
      const std::size_t count = std::min(lanes, rem);
      hn::Vec<decltype(d)> vals[21];
      for (int i = 0; i < diameter; ++i) {
        vals[i] = hn::LoadN(d, srcp_planes[i] + static_cast<std::size_t>(y) * src_stride + x, count);
      }
      const auto med = eval_temporal_median_int(d, vals, diameter);
      if (rem >= lanes) {
        hn::StoreU(med, d, dst_row + x);
      } else {
        hn::StoreN(med, d, dst_row + x, rem);
      }
    }
  }
}

// The largest fixed window can keep its selection network in registers.
// Other diameters retain the shared evaluator and their existing code size.
template <bool IsF16, typename StorageT>
void temporal_median21_float(const StorageT* const* planes, StorageT* dstp,
    int width, int height, std::size_t src_stride, std::size_t dst_stride) {
  const hn::ScalableTag<FloatLane<IsF16>> d;
  const auto lanes = hn::Lanes(d);
  for (int y = 0; y < height; ++y) {
    for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
      const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
      hn::Vec<decltype(d)> values[21];
      for (int i = 0; i < 21; ++i) {
        const auto* src = planes[i] + static_cast<std::size_t>(y) * src_stride + x;
        if constexpr (IsF16) values[i] = load_f16(d, src, count);
        else values[i] = hn::LoadN(d, src, count);
      }
      const auto result = median21(d, values);
      auto* dst = dstp + static_cast<std::size_t>(y) * dst_stride + x;
      if constexpr (IsF16) store_f16(d, result, dst, count);
      else hn::StoreN(result, d, dst, count);
    }
  }
}

template <bool IsF16, typename StorageT>
void temporal_median_float_impl(int diameter, const StorageT* const* srcp_planes, StorageT* dstp,
                                int width, int height, std::size_t src_stride, std::size_t dst_stride) {
  if (diameter == 21) {
    temporal_median21_float<IsF16>(srcp_planes, dstp, width, height, src_stride, dst_stride);
    return;
  }
  using ComputeT = FloatLane<IsF16>;
  const hn::ScalableTag<ComputeT> d;
  const std::size_t lanes = hn::Lanes(d);

  if constexpr (IsF16) {
    for (int y = 0; y < height; ++y) {
      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const auto count = std::min(lanes, static_cast<std::size_t>(width) - x);
        hn::Vec<decltype(d)> vals[21];
        for (int i = 0; i < diameter; ++i)
          vals[i] = load_f16(d, srcp_planes[i] + static_cast<std::size_t>(y) * src_stride + x, count);
        store_f16(d, eval_temporal_median_float<true>(d, vals, diameter), dst_row + x, count);
      }
    }
  } else {
    for (int y = 0; y < height; ++y) {
      StorageT* dst_row = dstp + static_cast<std::size_t>(y) * dst_stride;
      for (std::size_t x = 0; x < static_cast<std::size_t>(width); x += lanes) {
        const std::size_t rem = static_cast<std::size_t>(width - x);
        const std::size_t count = std::min(lanes, rem);
        hn::Vec<decltype(d)> vals[21];
        for (int i = 0; i < diameter; ++i) {
          vals[i] = hn::LoadN(d, srcp_planes[i] + static_cast<std::size_t>(y) * src_stride + x, count);
        }
        const auto med = eval_temporal_median_float<false>(d, vals, diameter);
        if (rem >= lanes) {
          hn::StoreU(med, d, dst_row + x);
        } else {
          hn::StoreN(med, d, dst_row + x, rem);
        }
      }
    }
  }
}

void dispatch_temporal_median_target(DataType dtype, int diameter, const std::uint8_t* const* srcp_planes,
                                     std::uint8_t* dstp, std::size_t width, std::size_t height,
                                     std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  std::array<const std::uint16_t*, 21> u16_planes{};
  std::array<const float*, 21> f32_planes{};
  for (int i = 0; i < diameter; ++i) {
    u16_planes[i] = reinterpret_cast<const std::uint16_t*>(srcp_planes[i]);
    f32_planes[i] = reinterpret_cast<const float*>(srcp_planes[i]);
  }
  const int w = static_cast<int>(width);
  const int h = static_cast<int>(height);

  if (dtype == DataType::U8) {
    temporal_median_int_impl<std::uint8_t>(diameter, srcp_planes, dstp, w, h, src_stride_bytes, dst_stride_bytes);
  } else if (dtype == DataType::U16) {
    temporal_median_int_impl<std::uint16_t>(
        diameter, u16_planes.data(),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F16) {
    temporal_median_float_impl<true, std::uint16_t>(
        diameter, u16_planes.data(),
        reinterpret_cast<std::uint16_t*>(dstp), w, h, src_stride_bytes / 2, dst_stride_bytes / 2);
  } else if (dtype == DataType::F32) {
    temporal_median_float_impl<false, float>(
        diameter, f32_planes.data(),
        reinterpret_cast<float*>(dstp), w, h, src_stride_bytes / 4, dst_stride_bytes / 4);
  }
}

} // namespace HWY_NAMESPACE
} // namespace neo_smo
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace neo_smo {
HWY_EXPORT(dispatch_temporal_median_target);

void process_temporal_median_plane(DataType dtype, int diameter, const std::uint8_t* const* srcp_planes,
                                   std::uint8_t* dstp, std::size_t width, std::size_t height,
                                   std::size_t src_stride_bytes, std::size_t dst_stride_bytes) {
  HWY_DYNAMIC_DISPATCH(dispatch_temporal_median_target)(
      dtype, diameter, srcp_planes, dstp, width, height, src_stride_bytes, dst_stride_bytes);
}
} // namespace neo_smo
#endif
