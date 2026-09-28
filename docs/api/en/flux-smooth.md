# FluxSmoothT, FluxSmoothST

[Back to the API index](README.md)

```python
core.neo_smo.FluxSmoothT(clip, temporal_threshold=None, planes=None, scalep=False)
core.neo_smo.FluxSmoothST(clip, temporal_threshold=None, spatial_threshold=None,
                        planes=None, scalep=False)
```

Supports the common formats. Thresholds are per-plane floating-point arrays, with the last element repeated. Each omitted threshold defaults to an automatically scaled value of 7 on the 8-bit scale. Explicit values use native units by default. With `scalep=True`, nonnegative values must be in 0–255 and are scaled. Negative thresholds disable the corresponding component; zero remains a valid threshold.

Smoothing is attempted only when the current pixel is strictly greater than both temporal neighbors or strictly less than both. FluxSmoothT checks those two neighbors; FluxSmoothST also checks the current frame's eight spatial neighbors. A neighbor participates in the average only if its absolute difference from the center is no greater than the corresponding threshold. The center always participates.

Even when a negative temporal threshold makes FluxSmoothST accumulate only spatial neighbors, the temporal-extremum condition still applies. It does not become a general spatial averaging filter. The first and last frames are copied; there is no scene-change detection.

```python
out = core.neo_smo.FluxSmoothT(src, temporal_threshold=[7], scalep=True)
out = core.neo_smo.FluxSmoothST(src, temporal_threshold=[7],
                              spatial_threshold=[5], planes=[0], scalep=True)
```

See [Temporal denoising and repair](../../knowledge/en/temporal.md) for how its denominator differs from TemporalSoften's.

## AviSynth

```text
neo_smo_FluxSmoothT(clip, temporal_threshold=Undefined(), planes=Undefined(), scalep=false)
neo_smo_FluxSmoothST(clip, temporal_threshold=Undefined(), spatial_threshold=Undefined(), planes=Undefined(), scalep=false)
```

Parameter names and order follow the shared description above. See the [API index](README.md#avisynth-calls-and-builds) for host formats, array syntax, and audio/parity preservation.

Thresholds and `planes` accept numeric arrays or scalars; `scalep` takes a boolean. Omitted thresholds use scaled defaults; explicit thresholds use native units unless scaling is enabled. Boundary copying and negative thresholds follow the rules above. Neither function has a scene-change argument.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_FluxSmoothT(src, temporal_threshold=[7], scalep=true)
```

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_FluxSmoothST(src, temporal_threshold=[7], spatial_threshold=[5], planes=[0], scalep=true)
```
