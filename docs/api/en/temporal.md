# TemporalMedian, TemporalSoften

[Back to the API index](README.md)

```python
core.neo_smo.TemporalMedian(clip, radius=1, planes=None, scenechange=False)
core.neo_smo.TemporalSoften(clip, radius=4, threshold=None, scenechange=0,
                          scalep=False, planes=None)
```

Supports the common formats. `radius` is a single integer in 1–10, giving a full window of `2r+1` frames. Both sample the same pixel coordinates across time, without motion compensation.

## TemporalMedian

Takes the median of the temporal window. The first and last `radius` frames are copied; this differs from clamping frame indices to the clip boundaries. `scenechange=True` shortens the window using existing `_SceneChangePrev` and `_SceneChangeNext` properties. It does not run scene detection and reports an error if a required property is missing. If the shortened window has an even number of samples, it averages the two middle values, rounding down for integer samples.

## TemporalSoften

`threshold` is a per-plane floating-point array. Omitted values use an automatically scaled value of 4 on the 8-bit scale. Explicit values are scaled from 8-bit units with `scalep=True`; otherwise they use native units. Only positive thresholds enable processing. Zero can disable a YUV plane; the first threshold cannot be zero for RGB/Gray, and all three thresholds cannot be zero.

Native integer thresholds range from 0 to the format maximum. Floating-point luma/RGB thresholds range from 0 to 1, and floating-point chroma thresholds from −0.5 to 0.5; nonpositive chroma thresholds disable processing. Scaled values must still pass the format-range check.

A neighboring frame's sample is used if its absolute difference from the center is no greater than the threshold; otherwise the center replaces it. The result averages the entire valid window. **The denominator is not just the number of samples that passed the threshold.** At clip boundaries, the window shrinks to the frames that exist; whole frames are not skipped.

| scenechange | Behavior |
|---|---|
| 0 (default) | Do not check scene properties |
| −1 | Use existing properties; missing properties provide no usable cut point |
| 1–254 | Invoke `misc.SCDetect` with this value / 255 as its threshold; automatic detection is unsupported for RGB |

Thus `scenechange=True` equals numeric 1 and enables automatic detection; it does not mean "read existing properties only."

```python
out = core.neo_smo.TemporalMedian(src, radius=2)
out = core.neo_smo.TemporalSoften(src, radius=3, threshold=[4, 6, 6], scalep=True)
# Requires the misc plugin providing SCDetect:
sc = core.misc.SCDetect(src, threshold=0.1)
out = core.neo_smo.TemporalMedian(sc, radius=2, scenechange=True)
```

See [Temporal denoising and repair](../../knowledge/en/temporal.md) and the [boundary behavior tables](../../knowledge/en/shared/temporal-boundaries.md).
