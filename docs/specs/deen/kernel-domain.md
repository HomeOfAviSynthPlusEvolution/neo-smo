# Sample domains, neighborhoods, and numerics

Specification: `DEEN-DOMAIN-001`. Applies to all six modes.

## Inputs and outputs

An operator receives the current plane, or the corresponding planes from three adjacent frames, plus an effective mode, radius r, spatial threshold Ts, temporal threshold Tt, and m. Input planes have identical sample types, bit depths, and positive dimensions W and H. Different planes and frames may have different row strides. Output dimensions and sample type match the current plane. Inputs are immutable.

Sample formats are unsigned 8- through 16-bit integers and IEEE binary16/binary32 floating point. Integer values lie in [0,Q], where Q=2^bits−1. Floating point uses unit amplitude: nominal luma/RGB range [0,1] and centered chroma range [−0.5,0.5]. Preserve finite values outside nominal ranges; do not pre-clamp or add a chroma offset. Consumed NaN/Inf samples are errors. Unprocessed planes may pass through bitwise without content validation unless scene detection reads them.

Radius r is an integer in [1,7] for 2D and [1,4] for 3D. Validate a 3D configuration against the 3D range even when an individual frame falls back to 2D. Parameter m must be finite and in [0,1]. Public threshold t must be finite and in [0,255]; fractional values are allowed. Convert to sample units as T=t·Q/255 for integers and T=t/255 for floating point. Thresholds are difference magnitudes, not pixel levels. Do not add a chroma offset, quantize thresholds to integers, or use the scale 2^(bits−8).

## Spatial coordinates

For output coordinate (x,y), let c be the current frame's center value. Define:

```text
Ω = {(dx,dy) | −r ≤ dx ≤ r, −r ≤ dy ≤ r}
N = (2r+1)²
s_f(dx,dy) = frame_f[clamp(y+dy,0,H−1), clamp(x+dx,0,W−1)]
ρ(dx,dy) = sqrt(dx²+dy²) / (r·sqrt(2))
g(dx,dy) = 1 − (1−m)·ρ(dx,dy)
```

Frame f=0 is current; f=−1/+1 denotes previous/next. Extend each frame and plane using its own visible dimensions. The same offset addresses the same spatial position in every frame; do not shift temporal references by r rows. Replicate the nearest visible sample at borders and filter every visible output pixel. Distinct offsets mapping to one address remain distinct taps. Radius is measured in samples of the current plane and is not scaled for chroma subsampling. This definition also covers 1×1, 1×H, and W×1 planes.

In exact arithmetic, ρ is in [0,1] and g is in [m,1]. The center has g=1; all four corners have g=m. Implementations should set these endpoints explicitly to avoid negative residual weights. Every difference comparison uses the current output center c, not the center of each reference frame.

## Numerical semantics

Equations describe the mathematical algorithm; implementations use ordinary floating-point arithmetic for thresholds, distances, and means. Comparisons use the computed values, including values rounded to equality. No particular SIMD accumulation order is prescribed. Final integer output is clamp(floor(z+0.5),0,Q), rounded only once at storage. Floating output uses round-to-nearest, ties-to-even in the destination format, without nominal-range clamping. Integer c/a sums and counts must be exact. Optimized division must preserve the specified final rounding. The maximum sum is 243·65535, so 16-bit accumulation is insufficient.

Do not truncate w weights or a distance thresholds into bytes or legacy fixed-point tables. Floating production paths may use FMA, reassociation, and validated reciprocal optimizations under the [acceptance criteria](acceptance.md). F16 describes input/output storage, not a requirement for F16 intermediates. Accumulation must provide at least F32 range and precision to avoid half-precision sum overflow. Instruction selection belongs to implementation design.

The mathematical result is a convex combination of participating values: constants must remain constant, and filtering must not introduce new extrema. Numerical recovery may constrain the computed result to the participating minimum/maximum, but must not conceal a wrong denominator, overflow, or non-finite intermediates. If finite inputs cause intermediate overflow, use stable evaluation or wider-precision recomputation, or report a frame error; never successfully return NaN/Inf.

Storage adaptation must not interpret row padding as samples. SIMD tails must not access outside valid allocations. Check arithmetic overflow in strides, dimensions, and workspace sizes. The core does not require in-place operation. Any integration offering it must preserve original input values for all subsequent tap reads.
