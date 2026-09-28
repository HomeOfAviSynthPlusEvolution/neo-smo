// Included inside each Highway target namespace.
namespace hn = hwy::HWY_NAMESPACE;
// Normal lanes stay vectorized; only ambiguous comparisons need exact arithmetic.
template <bool Adaptive, class D>
auto deen_selection(D d, hn::VFromD<D> sample, hn::VFromD<D> center, DeenThreshold threshold, double weight,
                    double minimum, int dx, int dy, int radius) {
  const auto difference = hn::Abs(hn::Sub(sample, center));
  const double limit = threshold.hi * (Adaptive ? weight : 1);
  auto pass = hn::Le(difference, hn::Set(d, limit));
  const auto equal = hn::Eq(sample, center);
  const auto uncertainty = hn::Mul(hn::Set(d, (Adaptive ? 16 : 8) * std::numeric_limits<double>::epsilon()),
                                   hn::Max(difference, hn::Set(d, threshold.hi)));
  const auto near = hn::AndNot(equal, hn::Le(hn::Abs(hn::Sub(difference, hn::Set(d, limit))), uncertainty));
  if (!hn::AllFalse(d, near)) {
    HWY_ALIGN double samples[hn::MaxLanes(d)], centers[hn::MaxLanes(d)], flags[hn::MaxLanes(d)], needs[hn::MaxLanes(d)];
    hn::Store(sample, d, samples);
    hn::Store(center, d, centers);
    hn::Store(hn::IfThenElse(pass, hn::Set(d, 1), hn::Zero(d)), d, flags);
    hn::Store(hn::IfThenElse(near, hn::Set(d, 1), hn::Zero(d)), d, needs);
    for (std::size_t i = 0; i < hn::Lanes(d); ++i)
      if (needs[i]) {
        if constexpr (Adaptive)
          flags[i] = deen_detail::adaptive_boundary(samples[i], centers[i], threshold, minimum, dx, dy, radius) ? 1 : 0;
        else
          flags[i] = deen_detail::within_threshold(samples[i], centers[i], threshold) ? 1 : 0;
      }
    pass = hn::Eq(hn::Load(d, flags), hn::Set(d, 1));
  }
  return hn::Or(pass, equal);
}
