# DegrainMedian

[Back to the API index](README.md)

```python
core.neo_smo.DegrainMedian(clip, limit=None, mode=[1], interlaced=False,
                         norow=False, scalep=False)
```

Supports the common formats. Selects a correction from directional neighbors in the previous, current, and next frames, then limits its distance from the original value. The first and last frames are copied.

| Parameter | Default and range | Meaning |
|---|---|---|
| limit | Omitted: native value 4 on every plane | Maximum per-pixel correction; per-plane array with the last element repeated |
| mode | `[1]`; each element 0–5 | Direction-selection cost; 0 still processes pixels |
| interlaced | False | Changes vertical neighbor spacing from one row to two, without changing temporal frame spacing |
| norow | False | Excludes the left/right horizontal pair in the current frame |
| scalep | False | Whether to scale explicitly supplied limit values |

There is no `planes` parameter. Use `limit=0` to disable a plane; all three limits cannot be zero. Explicit native limits range from 0 to the format maximum for integer samples, 0–1 for floating-point luma/RGB, and −0.5–0.5 for floating-point chroma. Nonpositive chroma limits leave the plane unprocessed.

With `scalep=True`, explicit limits must be in 0–255 and are scaled by `limit × format maximum / 255`. The floating-point chroma maximum is 0.5; other floating-point planes use 1. This differs from the integer bit-shift scaling used by most filters.

**Omitting limit retains the native value 4, including for floating-point input; scalep=True does not automatically turn it into 4/255.** Specify both limit and scalep explicitly when reusing parameters across bit depths.

```python
out = core.neo_smo.DegrainMedian(src, limit=[4, 0, 0], mode=[1], scalep=True)
```

See [Temporal denoising and repair](../../knowledge/en/temporal.md) for direction selection and mode costs.

## AviSynth

```text
neo_smo_DegrainMedian(clip, limit=Undefined(), mode=[1], interlaced=false, norow=false, scalep=false)
```

Parameter names and order follow the shared description above. See the [API index](README.md#avisynth-calls-and-builds) for host formats, array syntax, and audio/parity preservation.

`limit` and `mode` accept numeric arrays or scalars; `interlaced`, `norow`, and `scalep` take booleans. There is no `planes` argument. The example processes only luma by setting chroma limits to zero. Omitting `limit` still uses the native value 4.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_DegrainMedian(src, limit=[4, 0, 0], mode=[1], interlaced=false, norow=false, scalep=true)
```
