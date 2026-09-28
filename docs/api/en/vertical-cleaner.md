# VerticalCleaner

[Back to the API index](README.md)

```python
core.neo_smo.VerticalCleaner(clip, mode)
```

Supports the common formats. `mode` is a required per-plane integer array with values 0–2; its last element is repeated for remaining planes. There is no `planes` parameter.

| mode | Behavior | Dimensions and boundaries |
|---|---|---|
| 0 | Copy | No additional height requirement |
| 1 | Median of the previous, current, and next rows | Processed planes need at least 3 rows; copy the first and last row |
| 2 | Construct bounds from trends in two rows above and below, then clamp the current value | Processed planes need at least 5 rows; copy the first and last two rows |

The height check applies to each plane; YUV420 chroma, for example, has half the luma height. This filter does not sample horizontally or across frames.

```python
out = core.neo_smo.VerticalCleaner(src, mode=[1, 0, 0])
```

Mode 2 is not a five-point median. See [Spatial order statistics and clamping](../../knowledge/en/spatial.md).

## AviSynth

```text
neo_smo_VerticalCleaner(clip, mode)
```

Parameter names and order follow the shared description above. See the [API index](README.md#avisynth-calls-and-builds) for host formats, array syntax, and audio/parity preservation.

`mode` is required and accepts an integer or array. There are no `planes` or `y/u/v` arguments; `[1, 0, 0]` processes only Y in YUV or R in RGB.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_VerticalCleaner(src, mode=[1, 0, 0])
```
