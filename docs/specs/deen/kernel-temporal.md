# Temporal window and scene fallback

Specification: `DEEN-TEMPORAL-001`. Scene classification precedes pixel filtering and selects one effective 2D/3D mode for the entire frame.

## Frame set

A 2D mode uses only input n. A 3D mode uses n−1, n, and n+1; radius r changes only the spatial window. On the first or last frame, including a single-frame clip, use the corresponding 2D mode. Do not duplicate endpoint frames to fill the window or request out-of-range frames.

For n with both adjacent frames available, a scene cut on either side forces the entire frame to 2D. Do not continue with only one temporal reference or maintain a countdown inherited from earlier requests.

## Scene metric

Let P be the number of color planes, A/B adjacent input frames, and W_p/H_p each plane's visible dimensions. Define the plane mean absolute difference in 8-bit units:

```text
L = Q  for integer samples, or 1 for floating point
D_p(A,B) = 255/(L·W_p·H_p) · Σ[x,y] |A_p(x,y)−B_p(x,y)|
D(A,B) = (Σ[p=0..P−1] D_p(A,B)) / P
cut(A,B) = D(A,B) > scd
```

Use every color plane, including planes not selected for filtering. Planes have equal weight; do not pool their samples into a pixel-count-weighted mean. GRAY uses one plane. Read visible samples only, ignore padding, and use each frame's own stride. Do not derive the metric from output frames, frame properties, or previously filtered values.

Parameter scd is finite and nonnegative, in the same 8-bit difference units. Equality D=scd is not a cut. A binding may explicitly disable scene detection; endpoint fallback still applies. Do not retain the legacy negative-value debug-printing interpretation.

Integer SAD uses sufficiently wide exact accumulation with prior checks of dimensions and maximum sums. Floating differences and means use a stable reference calculation of at least binary64 precision. Production reductions may be optimized, but arbitrary rounding drift must not flip a whole-frame decision near scd. Re-evaluate in higher precision if the error bound straddles scd.

## Determinism and errors

Output n depends only on its input frames and immutable configuration. Forward, reverse, random, repeated, and concurrent requests must select the same effective mode. Caches may store pure calculation results but must not define subsequent frame state.

Failure to obtain a required reference, changed formats, or non-finite samples consumed by scene detection are frame errors, not implicit scene cuts. Neither 2D modes nor endpoint fallback require scene metrics. Non-finite values on an unprocessed plane are errors only when that plane is actually consumed by scene detection.

## Examples

For three constant GRAY frames 100/100/110 and scd=9, the current frame falls back to 2D because its right-hand difference is 10. With scd=10 it remains 3D, as it does with scene detection disabled.

For YUV plane mean differences 0/12/0, D=4 regardless of 4:4:4 or 4:2:0 subsampling. With scd=4 there is no cut. The two output frames adjacent to a cut can each identify it independently, without first requesting the preceding output frame.
