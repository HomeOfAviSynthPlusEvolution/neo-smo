# w2d / w3d: distance-weighted mean

Specification: `DEEN-W-001`. Uses the [common domain](kernel-domain.md) and [temporal selection](kernel-temporal.md).

Use the threshold-substituted values b_f(o) from [c modes](kernel-c.md), with spatial weights g(o). In 2D, α_0=1. In 3D, α_0=2 and α_−1=α_+1=1. Define:

```text
G = Σ[o∈Ω] g(o)
z = Σ[f∈F] α_f · Σ[o∈Ω] g(o)·b_f(o) / (G·Σ[f∈F] α_f)
```

Since G≥1, the denominator is strictly positive. Rejected taps retain their weights and substitute c. Weights depend on spatial distance only, not a three-dimensional distance including time. Form the weighted sum before normalization. Do not truncate individual tap products, transfer fixed-point residuals to the center, or add 1 to the final output.

In w3d, the current frame contributes half the total weight and each adjacent frame contributes one quarter. Fallback to w2d retains only the current frame and renormalizes. With m=0, corner weights vanish but the center still has weight 1. With m=1, w2d equals c2d; w3d does not equal c3d because their temporal weights differ.

## Examples

For r=1 and m=0, the center weight is 1, each axial neighbor has weight 1−1/√2, each corner has weight 0, and G=5−2√2. For center 100, axial neighbors 110, and sufficiently large thresholds, z=100+40·(1−1/√2)/(5−2√2), giving U8 output 105. Changing corner values does not change the mathematical output.

For constant previous/current/next frames 80/100/140, sufficiently large thresholds, and no scene fallback, w3d produces (80+2·100+140)/4=105 regardless of m or r.
