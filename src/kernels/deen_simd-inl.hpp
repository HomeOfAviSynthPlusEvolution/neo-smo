// Included inside each Highway target namespace.
namespace hn = hwy::HWY_NAMESPACE;
template <class D>
auto deen_selection(D d, hn::VFromD<D> sample, hn::VFromD<D> center, double threshold) {
  return hn::Le(hn::Abs(hn::Sub(sample, center)), hn::Set(d, threshold));
}
