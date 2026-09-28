# DCTFilter

[Back to the API index](README.md)

```python
core.neo_smo.DCTFilter(clip, factors, planes=None)
```

Applies a DCT, frequency-coefficient weighting, and an inverse transform to fixed 8×8 blocks in each selected plane. It operates on the current frame only, without motion compensation or temporal averaging. Supports constant-format, constant-size Gray, RGB, and YUV with 8–16-bit integer, F16, or F32 samples (F16 is available only in VapourSynth). Output retains the input format, dimensions, and frame properties.

## Parameters

| Parameter | Requirements and meaning |
|---|---|
| `clip` | Required input video |
| `factors` | Required array of exactly eight floating-point values, each finite and in `[0,1]`; no default |
| `planes` | Omit to process all planes, or supply a nonempty array of plane indices such as `[0]` or `[1,2]`; duplicate, out-of-range, and empty arrays are errors |

`factors` runs from DC to the highest frequency. The two-dimensional weight at horizontal frequency u and vertical frequency v is `factors[v] * factors[u]`; the eight values do not set strengths for eight spatial pixels. Every selected plane uses the same factors. The array is not extended automatically. There is no `scalep`, adjustable block-size, or `opt` parameter.

- `[1] * 8` is the identity in ideal arithmetic, but still performs a DCT round trip and does not guarantee bitwise unchanged output.
- `[1,0,0,0,0,0,0,0]` retains only each block's mean.
- `[0] * 8` sets selected planes to numeric zero. Zero is not neutral chroma in integer YUV, so this is generally unsuitable for chroma planes.
- The DC weight is `factors[0] ** 2`. Use `factors[0]=1` to preserve the block mean.

## Boundaries and precision

The non-overlapping block grid starts at the top-left corner of each plane. Width and height need not be multiples of 8 or 16. Missing samples at the right and bottom follow the upstream `resize.Point` behavior: reflect once, repeating the edge pixel, then clamp to the opposite edge if the extension exceeds the entire original plane. Only the original output area is written. Unselected planes retain their original content.

Both integer and floating-point input use F32 DCT arithmetic. F16 input is promoted to F32 and converted back to F16 only on output; this filter does not switch to native half-precision DCT arithmetic. Production builds allow FMA.

Integer output is clamped to the legal bit-depth range and rounded to the nearest integer, with positive halfway values rounded upward. Floating-point output is not forcibly clipped to `[0,1]` for luma/RGB or `[-0.5,0.5]` for chroma.

## Examples

```python
# Attenuate higher frequencies on luma only.
out = core.neo_smo.DCTFilter(
    src, factors=[1, .9, .7, .5, .3, .2, .1, 0], planes=[0])

# Retain each 8×8 luma block's mean to make the block grid visible.
blocks = core.neo_smo.DCTFilter(
    src, factors=[1, 0, 0, 0, 0, 0, 0, 0], planes=[0])
```

See [Block DCT and frequency weighting](../../knowledge/en/dct-filter.md) for the algorithm and numerical conventions.

## AviSynth

```text
neo_smo_DCTFilter(clip, factors, planes=Undefined())
```

Parameter names and order follow the shared description above. See the [API index](README.md#avisynth-calls-and-builds) for host formats, array syntax, and audio/parity preservation.

`factors` must contain exactly eight numbers; scalars are not repeated. List all eight values explicitly in AVS rather than copying Python syntax such as `[1] * 8`. `planes` accepts an integer or array. AVS supports planar integer and F32 input, not F16; DCT computation remains F32.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_DCTFilter(src, factors=[1, 0.9, 0.7, 0.5, 0.3, 0.2, 0.1, 0], planes=[0])
```
