# Acceptance criteria

Specification: `DEEN-ACCEPT-001`. This document defines conformance criteria for implementations.

## Independent mathematical reference

Implement the defining loops directly in tests. Integer c/a reference paths use exact sums, counts, and final rounding. Distances, weights, and floating means use at least binary64. Do not generate expected results through production weight tables, border-mapping helpers, reciprocal tables, or the same SIMD kernel. The original DLL, MiniDeen, and supplied reconstructed code are not correctness oracles.

Check the hand-calculated examples in each kernel document before broader format, spatial, and temporal combinations. For three-frame mixing examples, disable scene detection or choose a sufficiently large scd; otherwise the test exercises 2D fallback instead.

## Error assessment

FMA and floating rounding differences are permitted. Error allowances do not excuse displaced neighborhoods, incorrect denominators, wrong threshold scaling, or incorrect scene decisions.

- Integer c/a: for identical accepted sets, output must follow exact integer rounding.
- Integer w: output may differ from the rounded binary64 reference by at most one integer code value. Persistent one-sided bias, shifted constants, or differences away from rounding boundaries still require investigation; a difference of one is not automatic acceptance.
- F32: absolute output error must not exceed `2e−6·max(1,M)`, where M is the maximum absolute participating value for that output.
- F16: compare the stored output, promoted exactly to binary64, against the unrounded reference result z. Allow the F32 budget above plus ULP16(q), where q is z rounded to binary16 with round-to-nearest, ties-to-even. Output must remain finite; half-precision accumulation overflow is not an acceptable explanation.

Define ULP16(q)=2^−24 for zero and subnormal q. For normal finite q, define ULP16(q)=2^(floor(log2(abs(q)))−10). At an exact power of two, use the spacing on the larger-magnitude side; the sign of q does not change the budget. This gives 2^−24 at zero and the smallest normal, 2^−10 at ±1, and 32 at ±65504 without referring to an infinite adjacent value. A finite F16 input convex combination has a finite rounded reference q.

These are initial acceptance ceilings, not algorithm parameters relaxing the equations. All formats must preserve representable constant values; floating signed-zero identity is not required. Unprocessed planes must remain bitwise identical. Results must not introduce extrema outside participating values.

Selection and scene classification are discrete decisions. Construct equality, just-below, and just-above fixtures to check inclusive sample thresholds and strict scene thresholds. Small output error does not justify a wrong branch. Recheck distance thresholds or scene means near comparison boundaries in higher precision. If a production optimization's error bound crosses a threshold, handle that boundary with reliable evaluation rather than fitting a legacy binary's decision.

## Required matrix

| Dimension | Minimum coverage |
| --- | --- |
| Modes | All six; normal 3D, endpoint fallback, and scene fallback |
| Samples | U8, 10/12/16-bit U16, F16, F32; negative chroma and finite values outside nominal float ranges |
| Radius | Every allowed value, especially 1, 2D maximum 7, and 3D maximum 4 |
| Geometry | 1×1, 1×H, W×1, smaller than the window; SIMD width minus one/exact/plus one and nonmultiples |
| Storage | Different per-frame strides, unaligned starts, independent plane padding, guard pages or equivalent bounds checking |
| Parameters | Threshold 0/255/fractional; min=0/0.5/1; scd=0/boundary/disabled |
| Patterns | Zero, maximum, arbitrary constants, impulses, steps, ramps, fixed-seed random images |
| Planes | GRAY, RGB, YUV 4:4:4/4:2:2/4:2:0; individual selections, all planes, and unselected passthrough |

Core tests cover the algorithm's sample types. Host tests cover only types representable by that host, checking explicit rejection otherwise rather than reporting unsupported formats as passed. Check NaN/Inf consumption: processed planes and planes actually used for scene detection must fail; untouched passthrough planes retain their original bits.

## Relations and integration

- With min=1, w2d equals c2d; a2d remains a selection mean and must not be conflated with c2d.
- When all three frames are samplewise identical and Ts=Tt, 3D and corresponding 2D results are mathematically equal. Do not assume this when Ts and Tt differ.
- With all sample thresholds zero, only values equal to the current center pass; every mode returns the center.
- A 3D fallback equals explicit matching 2D filtering, without shifted reference coordinates.
- Changing padding never changes output. Request order and concurrency never change effective modes or output values.
- Compare each runnable Highway target against the reference. Unexecuted targets are not passed targets.
- VS/AVS integration covers defaults, invalid arguments, dependency failures, properties, AVS audio passthrough, and resource release.

Correctness does not require original-binary performance or matching original output. Keep performance benchmarking separate; any later measurements must identify the platform, effective mode, and sample format.
