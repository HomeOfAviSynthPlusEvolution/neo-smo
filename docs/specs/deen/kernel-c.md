# c2d / c3d: fixed-denominator mean

Specification: `DEEN-C-001`. Uses the [common domain](kernel-domain.md) and [temporal selection](kernel-temporal.md).

The effective frame set is F={0} for 2D and F={−1,0,+1} for 3D. T_f is Ts for f=0 and Tt otherwise. Define:

```text
b_f(o) = s_f(o)    if |s_f(o)−c| ≤ T_f
         c         otherwise
z = Σ[f∈F] Σ[o∈Ω] b_f(o) / (|F|·N)
```

Threshold equality is accepted. A rejected tap is replaced by c and still occupies one denominator position; do not omit it and reduce the denominator. All frames and offsets have equal weight. The current center enters once through (f=0,o=(0,0)); do not seed an extra center value or count.

Parameter min has no effect on c modes. Current and adjacent frames use their respective thresholds but all comparisons reference the current output center c. Temporal fallback evaluates exactly c2d with N taps, without retaining the 3N denominator.

## Examples

For r=1, an interior center 100, one neighbor 110, seven other neighbors 130, and Ts=10, only the center and 110 pass. Replacing other values with 100 gives z=910/9 and U8 output 101. With Ts below 10 the output is 100.

For constant previous/current/next frames 80/100/140, sufficiently large thresholds, and no scene fallback, c3d produces round(320/3)=107. With Tt=20, 140 is replaced by 100 and the output is round(280/3)=93.
