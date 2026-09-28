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
