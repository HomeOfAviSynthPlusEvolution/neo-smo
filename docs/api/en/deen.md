# Deen, MiniDeen

[Back to the API index](README.md)

## Deen

```python
core.neo_smo.Deen(clip, mode="c3d", rad=1, thrY=7.0, thrUV=9.0,
                  tthY=4.0, tthUV=6.0, min=0.5, scd=9.0,
                  scenechange=True, planes=None)
```

Spatial or three-frame thresholded smoothing. Requires constant-format, constant-size planar GRAY/YUV/RGB: 8–16-bit integer, F16 or F32 in VapourSynth; 8/10/12/14/16-bit integer or F32 in AviSynth+. Alpha and packed formats are unsupported.

### Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `mode` | `"c3d"` | One of `c2d`, `c3d`, `a2d`, `a3d`, `w2d`, `w3d`; see below. |
| `rad` | `1` | Spatial radius in each plane's pixels. Integer 1–7 for 2d modes, 1–4 for 3d modes. Window size is `(2*rad+1)²`. |
| `thrY`, `thrUV` | `7`, `9` | Current-frame thresholds on a 0–255 scale. `thrUV` applies to YUV chroma; GRAY and all RGB planes use `thrY`. |
| `tthY`, `tthUV` | `4`, `6` | Previous/next-frame thresholds, with the same units and plane mapping. Used in 3d modes. |
| `min` | `0.5` | Distance factor at the window corners, from 0 to 1. Controls thresholds in a modes and weights in w modes; has no effect in c modes. |
| `scd` | `9` | Nonnegative scene-change threshold on the 0–255 difference scale. |
| `scenechange` | `True` | Enable built-in scene detection in 3d modes; no scene-change frame properties are required. |
| `planes` | all | Distinct valid plane indices; omitted selects all planes; `[]` is rejected. Unselected planes are copied. |

The four filtering thresholds must be finite and within 0–255. They always scale by `(2**bits-1)/255` for integer input or `1/255` for floating input. There is no `scalep` parameter. Floating samples must be finite; they are not restricted to 0–1.

### Modes

Every candidate is compared to the current frame's center pixel, with equality accepted. Spatial coordinates outside the plane are clamped to its edges.

| Family | Behavior |
|---|---|
| `c` | Accept samples within the threshold; replace rejected samples with the center, then average over the fixed window size. |
| `a` | Multiply the threshold by the distance factor; average only accepted samples, with an equal weight for each. The denominator is the accepted count. |
| `w` | Use the unmodified threshold; replace rejected samples with the center, then average with distance weights. In 3d, current-frame weights are doubled, giving frame totals in a 2:1:1 ratio. |

For spatial offset `(dx, dy)`, the distance factor is `1 - (1-min)*sqrt((dx²+dy²)/(2*rad²))`: 1 at the center and `min` at the corners. In a modes, a zero threshold still accepts equal samples. In w modes, zero weights contribute nothing.

2d modes use only the current frame. 3d modes use the current, previous and next frames, with the same spatial radius in each. The first and last frames use the corresponding spatial operation. If scene detection finds a cut on either side, both temporal neighbors are omitted and the current frame is filtered spatially.

Scene detection computes the mean absolute difference per plane, converts each mean to the 0–255 scale, and averages the plane means equally. It uses all input planes, including those excluded by `planes`. A cut requires a metric strictly greater than `scd`. Setting `scd=0` does not disable detection; use `scenechange=False`.

Integer c/a means use exact half-up rounding. Integer w output can differ from the rounded binary64 weighted equation by at most one code value; constant inputs remain unchanged and no new extrema are introduced. All floating modes read F16/F32 samples directly and calculate in F32. Threshold comparisons use the F32 difference and threshold, accepting rounded equality; a-mode thresholds include the distance factor before conversion to F32. Results can vary with arithmetic rounding and the SIMD target. `w2d` with `min=1` uses the exact c2d mean for integer input.

```python
# Spatial weighted smoothing on luma only.
out = core.neo_smo.Deen(src, mode="w2d", rad=2, thrY=7, min=0.5, planes=[0])
# Three-frame adaptive smoothing, with built-in scene detection.
out = core.neo_smo.Deen(src, mode="a3d", rad=2, thrY=7, tthY=4)
```

## MiniDeen

```python
core.neo_smo.MiniDeen(clip, radius=None, threshold=None, planes=None)
```

Spatial smoothing for planar GRAY/YUV/RGB integer input only: 8–16 bits in VapourSynth and 8/10/12/14/16 bits in AviSynth+. Floating input is unsupported. There is no temporal processing or scene detection.

- `radius`: up to three integers in 1–7, default `[1]`.
- `threshold`: up to three integers in 0–255, default `[10]`. Each value scales to `floor(threshold*(2**bits-1)/255)`.
- Short arrays repeat their last value across the remaining planes. Omitted or empty radius/threshold arrays select defaults.
- `planes`: distinct valid indices; omitted selects all planes; an empty array is rejected. A scaled threshold of 0 or 1 copies that plane unchanged.

Only samples whose difference from the center is **strictly less** than the scaled threshold participate. The square window includes the center, and two additional copies of the center are included in the sum and count. At image edges the window is clipped to valid coordinates; edge samples are not replicated. The resulting integer mean uses half-up rounding.

```python
out = core.neo_smo.MiniDeen(src, radius=[2, 1], threshold=[10, 6], planes=[0, 1, 2])
```

## AviSynth

```text
neo_smo_Deen(clip, mode="c3d", rad=1, thrY=7.0, thrUV=9.0, tthY=4.0, tthUV=6.0,
             min=0.5, scd=9.0, scenechange=true, planes=Undefined())
neo_smo_MiniDeen(clip, radius=Undefined(), threshold=Undefined(), planes=Undefined())
```

Names, order, defaults and behavior follow the descriptions above. Array parameters also accept a scalar as a one-element array. See the [AviSynth interface](README.md#avisynth-calls-and-builds) for loading and common host conventions.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_Deen(src, mode="w3d", rad=2, thrY=7, tthY=4, planes=[0])
```
