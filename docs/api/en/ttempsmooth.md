# TTempSmooth

[Back to the API index](README.md)

```python
core.neo_smo.TTempSmooth(clip, maxr=3, thresh=[4, 5, 5], mdiff=[2, 3, 3],
                       strength=2, scthresh=12.0, fp=True, pfclip=None, planes=None)
```

Accepts Gray/RGB/YUV with 8–16-bit integer or F32 samples. **F16 is not supported.**

| Parameter | Range | Meaning |
|---|---|---|
| maxr | 1–7 | Maximum temporal radius |
| thresh | Each element 1–256 | Per-plane difference thresholds on the 8-bit scale |
| mdiff | Each element 0–255 | Per-plane interval before difference weights start to decay |
| strength | 1–8 | Range retaining high temporal weights |
| scthresh | −1–100 | Scene-change control; see below |
| fp | True/False | How unused weight is handled; **not a floating-point precision switch** |
| pfclip | Defaults to the main clip | Prefiltered reference for difference tests; output samples still come from the main clip |
| planes | All by default | Planes to process |

Short thresh/mdiff arrays repeat their last element. There is no `scalep` parameter: these values always use their specified 8-bit scale. `pfclip` must match the main clip's format and dimensions and be temporally aligned with it.

In VapourSynth, `scthresh>0` invokes `misc.SCDetect` with threshold `scthresh/100`, detecting on the reference clip when pfclip is supplied. RGB does not support this automatic detection path. Zero disables scene checks completely. Negative values read existing scene properties; use −1. Missing properties supply no cut point in this mode. In VapourSynth, the default 12.0 therefore requires the misc plugin providing SCDetect. AviSynth+ includes equivalent detection without an additional plugin.

`fp=True` assigns unused normalized weight to the center pixel. `False` renormalizes using the weight actually accumulated. This changes filter behavior; it is not a speed preset.

Frame indices outside the clip are clamped to its first or last frame. Reference samples are checked from near to far; when a sample fails the difference conditions in one direction, more distant samples in that direction are not treated as independent candidates.

```python
# A self-contained call without SCDetect:
out = core.neo_smo.TTempSmooth(src, maxr=3, scthresh=0)
# Guide temporal weights with a spatially filtered result:
guide = core.neo_smo.Median(src, radius=[1])
out = core.neo_smo.TTempSmooth(src, pfclip=guide, scthresh=0)
```

See [Temporal denoising and repair](../../knowledge/en/temporal.md).

## AviSynth

```text
neo_smo_TTempSmooth(clip, maxr=3, thresh=[4, 5, 5], mdiff=[2, 3, 3], strength=2, scthresh=12.0, fp=true, pfclip=Undefined(), planes=Undefined())
```

Parameter names and order follow the shared description above. See the [API index](README.md#avisynth-calls-and-builds) for host formats, array syntax, and audio/parity preservation.

`thresh`, `mdiff`, and `planes` accept integer arrays or scalars; `fp` takes a boolean and `pfclip` takes a clip. The default `scthresh=12.0` uses built-in detection on `pfclip` when supplied, without the misc plugin. For RGB, use 0 to disable detection or −1 to read previously supplied properties. F16 is also unsupported on AVS.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
guide = neo_smo_Median(src, radius=[1])
return neo_smo_TTempSmooth(src, pfclip=guide, scthresh=12.0, fp=true)
```
