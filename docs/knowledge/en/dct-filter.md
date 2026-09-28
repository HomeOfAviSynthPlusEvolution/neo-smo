# Block DCT and frequency weighting

[Back to the knowledge index](README.md) · [DCTFilter API](../../api/en/dct-filter.md)

## Transform and weights

DCTFilter divides each plane into non-overlapping 8×8 blocks aligned to its top-left corner. Each block undergoes a two-dimensional DCT-II, element-wise coefficient weighting, and a normalized inverse transform. Using an orthonormal DCT convention, its mathematical action is:

```text
X[v,u] = DCT2(input)[v,u]
W[v,u] = factors[v] * factors[u]
output = IDCT2(X * W)          # * denotes element-wise multiplication
```

Here u is horizontal frequency and v is vertical frequency; 0 denotes DC and 7 the highest frequency along that axis. Weights are the outer product of the same eight factors. This interface cannot independently set arbitrary weights for all 64 coefficients or separate horizontal and vertical factor arrays.

For example, with `factors[0]=1` and `factors[1]=0.5`, `W[0,1]=W[1,0]=0.5`, but `W[1,1]=0.25`. The DC coefficient's weight is `factors[0]²`. Retaining only DC produces each block's mean; for edge blocks, this includes the added samples.

Attenuating high frequencies reduces both noise and detail. Blocks do not overlap, and neighboring blocks are not blended, so stronger attenuation can make block boundaries more visible. This is neither whole-frame FFT convolution nor a sliding per-pixel window filter.

## Edge extension

Only the last block at the right or bottom can be incomplete. The filter supplies the missing samples internally and preserves the original output dimensions; callers do not need to resize or crop the clip.

For a row or column of length N, sampling beyond its right or bottom edge uses:

```text
index(i) = i                      if 0 <= i < N
index(i) = max(2*N - 1 - i, 0)    if i >= N
```

For example, `[a,b,c]` extends to `[a,b,c,c,b,a,a,a,...]`. This neither repeats c indefinitely nor reflects periodically forever. Each plane has its own dimensions and block grid; an 8×8 YUV420 chroma block covers the image area of 16×16 luma pixels.

## Precision and normalization

User factors are multiplied in double precision before conversion to 64 F32 weights, matching Zig's construction order. All input pixels are converted to F32 for computation. F16 is used only for input/output storage, without intermediate half-precision rounding in the transform.

The forward and inverse transforms include the required normalization. Identity weights ideally reproduce the input, subject to actual floating-point rounding. Callers should not rescale the input or factors to compensate for transform normalization.

Production builds allow FMA. Integer output is clipped and rounded, with positive halfway values rounded upward. Floating-point output is not clipped; F16 output is finally converted back to half precision. These numerical effects must be distinguished from wrong coefficient indices, duplicated normalization, or incorrect edge extension.
