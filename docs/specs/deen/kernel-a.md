# a2d / a3d: distance-threshold selection mean

Specification: `DEEN-A-001`. Uses the [common domain](kernel-domain.md) and [temporal selection](kernel-temporal.md).

Use F={0} for 2D and F={−1,0,+1} for 3D. T_f is Ts for the current frame and Tt for adjacent frames. Define:

```text
h_f(o) = T_f · g(o)
E = {(f,o) | f∈F, o∈Ω, |s_f(o)−c| ≤ h_f(o)}
z = Σ[(f,o)∈E] s_f(o) / |E|
```

Thresholds decrease linearly with two-dimensional Euclidean distance from T_f at the center to m·T_f at corners. Preserve fractional thresholds; do not round either the decrement or its result. Equality is accepted. Rejected taps contribute to neither numerator nor denominator. All accepted taps have weight 1.

The current-frame offset (0,0) is included exactly once and always passes, so |E|≥1. At borders, other offsets may map to the same center address; they remain separate taps and are not deduplicated. Centers from adjacent frames must pass Tt like any other sample. There is no w3d-style 2:1:1 temporal weighting. Fallback is exactly a2d. Setting m=1 removes distance attenuation of thresholds but does not turn an a mode into a c mode.

## Examples

For r=1, m=1, Ts=10, center 100, one neighbor 110, and seven other neighbors 130, E contains two samples and output is 105. The corresponding c2d output is 101.

For r=1, m=0.5, and Ts=20, the corner threshold is 10 and the axial threshold is 20·(1−0.5/√2), approximately 12.9289. An axial sample differing by 12 passes; a corner sample differing by 12 fails. A corner sample differing by exactly 10 passes.

With m=0, corner thresholds are zero. Corner samples equal to c still pass and count, unlike zero-weight corners in w modes. If all neighbors fail, output remains c.
