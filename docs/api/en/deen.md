# Deen, MiniDeen

[Back to the API index](README.md)

```python
core.neo_smo.Deen(clip, mode="c3d", radius=None, threshold=None,
                  temporal_threshold=None, minimum=None, scenechange=0,
                  scalep=False, planes=None)
core.neo_smo.MiniDeen(clip, radius=None, threshold=None,
                      scalep=False, planes=None)
```

Both accept constant-format, constant-size planar Gray/YUV/RGB without alpha. VapourSynth supports 8–16-bit integers, F16 and F32; AviSynth+ supports 8/10/12/14/16-bit integers and F32. Processed floating samples must be finite, but need not lie in 0–1.

## Shared parameters

| Parameter | Default | Meaning |
|---|---|---|
| `radius` | `[1]` | Per-plane spatial radius: 0–7 for MiniDeen and Deen 2d, 0–4 for Deen 3d. Zero copies the plane, including disabling its temporal processing. |
| `threshold` | see below | Per-plane floating thresholds for samples from the current frame. |
| `scalep` | `False` | Whether explicitly supplied thresholds use the 8-bit scale. |
| `planes` | all | Distinct valid plane indices. Unselected planes are copied. Empty, duplicate and out-of-range indices are errors. |

`radius`, `threshold`, `temporal_threshold` and `minimum` accept a scalar or an array with 1 to the actual number of planes. Short arrays repeat the last value; Gray accepts one value only. Indices refer to input planes, independently of `planes`. Omit a parameter to use its defaults; explicitly empty arrays are errors. Values for unselected planes are validated too. NaN and infinity are rejected in numeric parameters.

Pass integers to integer parameters. The VS Python binding may coerce a floating value to an integer before the plugin receives it; the plugin cannot recover the original type.

### Threshold units

- `scalep=False`: explicit values use native sample units, from 0 to `(2**bits-1)` for integers or 0 to 1 for floating input.
- `scalep=True`: explicit values must be in 0–255. Multiply by `2**(bits-8)` for integer input or divide by 255 for floating input.
- Omitted thresholds always scale from the defaults below, regardless of `scalep`.
- Differences on signed floating chroma use the same nonnegative threshold range, 0–1.

| Default in 8-bit units | Gray | RGB | YUV |
|---|---|---|---|
| Deen `threshold` | `[7]` | `[7,7,7]` | `[7,9,9]` |
| Deen `temporal_threshold` | `[4]` | `[4,4,4]` | `[4,6,6]` |
| MiniDeen `threshold` | `[10]` | `[10,10,10]` | `[10,10,10]` |

An explicit `threshold=7` repeats 7 across all planes. It does not generate a separate chroma default. For 16-bit input, `threshold=10, scalep=True` means 2560. Use `planes` or `radius=0` to disable a plane; zero threshold remains an algorithm threshold.

Deen accepts equality (`abs(sample-center) <= threshold`); MiniDeen uses strict comparison (`<`). Fractional integer thresholds preserve that distinction: a difference of 10 is accepted by MiniDeen at 10.5, but rejected at 10.

## Deen

Spatial or fixed three-frame smoothing, with a global `mode`:

| Family | Behavior |
|---|---|
| `c2d`, `c3d` | Replace rejected samples with the center and average over the fixed window size. |
| `a2d`, `a3d` | Multiply thresholds by the distance factor and average only accepted samples. |
| `w2d`, `w3d` | Replace rejected samples with the center and average with distance weights. In 3d, frame weights are 2:1:1 for current:previous:next. |

| Parameter | Default | Meaning |
|---|---|---|
| `mode` | `"c3d"` | One of the six modes above. |
| `temporal_threshold` | see table | Per-plane thresholds for both neighboring frames; units follow `scalep`. |
| `minimum` | `[0.5]` | Per-plane distance factor at window corners, in 0–1. Adjusts thresholds in a modes and weights in w modes. |
| `scenechange` | `0` | Scene handling described below. |

For offset `(dx,dy)`, the distance factor is `1-(1-minimum)*sqrt((dx²+dy²)/(2*radius²))`. Spatial coordinates outside the plane replicate its edges. Every candidate is compared to the current frame's center. `minimum` is unused in c modes; `temporal_threshold` and `scenechange` are unused in 2d modes, but supplied values must still be valid.

3d uses the previous, current and next frames. First/last frames use the matching 2d operation. If either neighboring boundary is a cut, both temporal neighbors are omitted.

### Scene detection

Uses the same convention as TemporalSoften, independently of `scalep`:

- `0`: disabled (default).
- `-1`: read `_SceneChangePrev` and `_SceneChangeNext` from the current frame. Missing properties mean no cut.
- `1..254`: automatic detection with threshold `scenechange/255`. VS uses `misc.SCDetect` when available and automatically falls back to the shared built-in luma detector when the plugin or function is absent. AVS uses that built-in detector. Both paths write `_SceneChangePrev` and `_SceneChangeNext`. Automatic detection supports Gray/YUV only and requires at least two frames. For RGB use 0 or supplied properties with -1.

Integer c/a means use exact half-up rounding. Integer w output can differ from the rounded binary64 weighted equation by at most one code value; constant input remains unchanged. `w2d` with `minimum=1` uses the exact c2d mean. Floating modes calculate differences and sums in F32; comparisons use F32 thresholds, including the distance factor in a modes. Arithmetic rounding can vary by SIMD target.

```python
out = core.neo_smo.Deen(src, mode="w2d", radius=[2,1], threshold=[7,9],
                       minimum=0.5, scalep=True, planes=[0])
out = core.neo_smo.Deen(src, mode="a3d", radius=2, threshold=7,
                       temporal_threshold=4, scenechange=12, scalep=True)
```

## MiniDeen

Spatial processing only. The square window is clipped to valid coordinates, includes the center, and adds two extra copies of the center to the sum and count. Only samples with a difference strictly below the threshold participate. Integer output uses exact half-up rounding. F16/F32 use F32 center-relative accumulation; F16 is converted when stored. Software targets use plain C++ implementations.

```python
out = core.neo_smo.MiniDeen(src, radius=[2,1], threshold=[10,6], scalep=True)
```

## AviSynth

Use `neo_smo_Deen` and `neo_smo_MiniDeen`, with the same names, parameter order, defaults and semantics as above. Use `Undefined()` to omit a value and `false`/`true` for booleans. Scalar shorthand and native numeric arrays are supported. See the [AviSynth interface](README.md#avisynth-calls-and-builds).

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_Deen(src, mode="w3d", radius=2, threshold=7,
                    temporal_threshold=4, scalep=true, planes=[0])
```

## Migration from the withdrawn interface

The old `rad`, `thrY`/`thrUV`, `tthY`/`tthUV`, `min` and `scd` names are removed. Use `radius`, per-plane `threshold`/`temporal_threshold`, `minimum` and the scene convention above. Scaling now follows the project's power-of-two rule rather than `(2**bits-1)/255`. MiniDeen thresholds accept fractional values; empty radius/threshold arrays are rejected. Scene detection defaults to off.
