# CCD

[Back to the API index](README.md)

```python
core.neo_smo.CCD(clip, threshold=4.0, temporal_radius=0,
                points=[1, 1, 0], scale=None, ref=None)
```

Accepts 8–16-bit integer, F16, and F32 RGB/YUV. Gray is not supported. RGB processes all three planes; YUV processes U/V and retains the main clip's Y plane. There is no `planes` or `scalep` parameter.

| Parameter | Range and default | Meaning |
|---|---|---|
| threshold | Default 4.0; must be finite, with a finite square after scaling | Controls the multichannel distance test; used internally as a squared threshold |
| temporal_radius | 0–10; default 0 | Expands the time window used to estimate differences |
| points | Exactly three Boolean integers: 0 disables, nonzero enables; default `[1,1,0]`; at least one group must be enabled | Selects the near, middle, and far spatial sampling groups |
| scale | Finite and ≥1; defaults to chroma-plane height / 240, or image height / 240 for RGB | Scales the spacing between sampling points |
| ref | Defaults to the main clip | Reference for difference tests; must match the input format and dimensions |

The near group contains four corner points at distance 4, the middle group eight points at distance 8, and the far group twelve outer points at distance 12. See [Color denoising](../../knowledge/en/chroma.md) for the layout. These are not dense square convolution windows.

The default scale can be less than 1 and cause creation to fail, for example with YUV420 below 480 pixels high. You can specify `scale=1.0`, but the image must still accommodate the selected geometry. The current implementation also validates width and height against the largest sampling span and scale; larger scales are more likely to exceed those limits.

For subsampled YUV, reference luma is resized to chroma dimensions with Bilinear, requiring the standard `std`/`resize` plugins. The first and last `temporal_radius` frames are copied unchanged from the main clip, without spatial filtering. If no frame has a complete temporal window, the entire clip remains unchanged. There is no scene-change parameter.

Increasing `temporal_radius` adds temporal information to the acceptance test. Accepted spatial points contribute **pixels from the current frame of the main clip**, rather than directly averaging frames across time. A negative threshold is also squared and therefore behaves like its positive counterpart; use nonnegative values for normal calls.

```python
out = core.neo_smo.CCD(src, threshold=4.0, temporal_radius=0, scale=1.0)
```
