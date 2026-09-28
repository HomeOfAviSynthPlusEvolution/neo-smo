# Samples, thresholds, and precision

[Back to the knowledge index](../README.md) · [Common API conventions](../../../api/en/README.md)

## Planes and numerical ranges

The nominal integer range is `0 … 2^bits−1`. Floating-point luma and RGB typically use 0–1, while floating-point YUV chroma typically uses −0.5–0.5. Neutral floating-point chroma is 0; do not substitute the midpoint used by integer chroma.

These are nominal ranges used by parameters and algorithms, not a promise that every filter clips every output to them. Some order-statistic operations return input samples directly, while some repair operations have explicit bounds.

A radius is measured on each plane's own pixel grid. A 3×3 YUV420 chroma neighborhood covers a larger image area than a 3×3 luma neighborhood. CCD separately aligns reference luma to the chroma grid; ordinary per-plane median filters do not.

## Explicit and default thresholds

For the usual functions with scalep, explicitly setting `scalep=True` gives:

| Input format | An 8-bit-scale value t becomes |
|---|---|
| 8-bit integer | t |
| 10-bit integer | 4t |
| 16-bit integer | 256t |
| F16/F32 | t/255 |

For example, TemporalSoften's threshold 4 becomes 16 on 10-bit input when written as `threshold=[4], scalep=True`. Omitting threshold also uses that scaled default. Writing `threshold=[4]` alone uses the native value 4.

DegrainMedian's explicit limit is an exception: integer input scales by `(2^bits−1)/255`, floating-point luma/RGB by `1/255`, and floating-point chroma by `0.5/255`. Omitting limit always retains the native value 4. With floating-point input, omission can therefore behave substantially differently from explicitly using `[4], scalep=True`.

TTempSmooth's thresh/mdiff, Cnr4's sense/str, and CCD's threshold each have their own normalization rules and no scalep parameter. Do not apply one conversion function to all interfaces.

## F16 and F32

F16 is the half-precision input/output storage format, not another name for an F32 array. For kernels using half-precision arithmetic, compilation targets with native support can compute directly in F16; other targets compute in F32 and convert output back to F16. Some steps still use wider types as the algorithm requires, so this does not mean that every instruction in a filter uses 16-bit arithmetic.

DCTFilter always computes transforms in F32, including for F16 input. It converts half precision only at input/output and does not change arithmetic precision according to native F16 instruction support. Its factors do not use bit-depth threshold scaling, and floating-point output is not clipped to nominal ranges. See [Block DCT and frequency weighting](../dct-filter.md).

The native-half and F32-compute paths have different intermediate rounding. Summing several fractional values with half-precision rounding after each addition can produce a neighboring but different result from accumulating in F32 and converting only at the end. F16 input is not necessarily faster than F32; the benefit depends on the CPU, compilation target, and algorithm.

Production builds allow FMA. A fused multiply-add rounds once and can differ in low-order bits from separate multiplication and addition. Reference comparisons should distinguish numerical error from incorrect indexing, sampling, boundaries, or modes. Production paths need not add half-precision round trips solely for bitwise reference matching.

## SIMD and interpreting tests

Highway selects a SIMD path from the targets available in the build and on the running platform. The public interface has no SIMD OFF or opt parameter. Operation order can differ across targets, so floating-point output is not guaranteed to be bitwise identical.

When checking output, record the format, parameters, frame number, plane, and maximum difference. Near a threshold, one rounding difference can change whether a sample is accepted, so a count of passing tests alone does not describe the actual error. Conversely, not every difference should automatically be attributed to FMA.
