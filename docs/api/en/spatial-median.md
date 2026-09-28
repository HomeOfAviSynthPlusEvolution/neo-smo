# Median, InterQuartileMean, SmartMedian

[Back to the API index](README.md)

```python
core.neo_smo.Median(clip, radius=[1], planes=None)
core.neo_smo.InterQuartileMean(clip, radius=[1], planes=None)
core.neo_smo.SmartMedian(clip, radius=[1], threshold=None, scalep=False, planes=None)
```

Supports the common formats. All three use a `(2r+1) × (2r+1)` spatial neighborhood. `radius` is a per-plane integer array with values 0–3, defaulting to 1; 0 copies the plane. `planes` further restricts which planes are processed.

| Function | Output |
|---|---|
| Median | Median of the neighborhood, including the center |
| InterQuartileMean | Weighted mean of the middle 50% of sorted samples, with fractional weights at quartile boundaries |
| SmartMedian | Find the two middle values after excluding the center; if the neighborhood passes a flatness test, clamp the center between them, otherwise retain it |

SmartMedian's `threshold` is a per-plane floating-point array. Omitted values use 50 on the 8-bit scale for radius 1, or 128 for radius 2/3, automatically scaled to the format. Explicit values with `scalep=True` must be in 0–255. With `False`, they must be in 0–the format maximum for integer input or 0–1 for floating-point input. A higher threshold makes correction more likely; it is not a "maximum difference between the center and the median."

Boundaries use mirrored sampling, with out-of-range coordinates on very small planes clamped to the valid range. Border pixels are not automatically left unchanged as an entire outer ring.

```python
# Process luma only.
out = core.neo_smo.Median(src, radius=[1], planes=[0])
# Apply the interquartile mean to every plane.
out = core.neo_smo.InterQuartileMean(src, radius=[2])
# Explicitly use threshold units on the 8-bit scale.
out = core.neo_smo.SmartMedian(src, radius=[1], threshold=[50], scalep=True)
```

See [Spatial order statistics and clamping](../../knowledge/en/spatial.md) for the algorithms and worked examples.

## AviSynth

```text
neo_smo_Median(clip, radius=[1], planes=Undefined())
neo_smo_InterQuartileMean(clip, radius=[1], planes=Undefined())
neo_smo_SmartMedian(clip, radius=[1], threshold=Undefined(), scalep=false, planes=Undefined())
```

Parameter names and order follow the shared description above. See the [API index](README.md#avisynth-calls-and-builds) for host formats, array syntax, and audio/parity preservation.

`radius`, `threshold`, and `planes` accept numeric arrays or scalars. SmartMedian takes a boolean `scalep`; omitting `threshold` retains the radius-dependent internal default.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_Median(src, radius=[1], planes=[0])
```

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_InterQuartileMean(src, radius=[2])
```

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_SmartMedian(src, radius=[1], threshold=[50], scalep=true)
```
