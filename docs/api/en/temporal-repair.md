# TemporalRepair

[Back to the API index](README.md)

```python
core.neo_smo.TemporalRepair(clip, repairclip, mode=[0], planes=None)
```

Supports the common formats. `repairclip` is required and must match the main clip's format and dimensions. Uses the current main frame and the previous, current, and next reference frames. The first and last frames are copied from the main clip. There is no scene-change parameter.

`mode` is a per-plane integer array with values 0–4, repeating its last element. **Mode 0 is an active repair mode, not an off switch.** Use `planes` to select which planes to process.

| mode | Source of the constraints |
|---|---|
| 0 (default) | Minimum and maximum of the same-position pixels in three reference frames |
| 1 | Expand the current reference value using positive and negative temporal changes in its eight neighbors, then include the previous and next centers |
| 2 | Construct a symmetric range from the largest absolute temporal change in the reference 3×3 neighborhood |
| 3 | Find the maximum difference to the previous frame and to the next frame over the reference 3×3 neighborhood; use the smaller maximum for a symmetric range |
| 4 | Use only the three reference centers to construct bounds from the temporal trend |

Modes 1–3 mirror out-of-range spatial coordinates. Modes are not increasing strength levels.

```python
filtered = core.neo_smo.Median(src, radius=[1])
out = core.neo_smo.TemporalRepair(filtered, src, mode=[0], planes=[0])
```

See [Temporal denoising and repair](../../knowledge/en/temporal.md).

## AviSynth

```text
neo_smo_TemporalRepair(clip, repairclip, mode=[0], planes=Undefined())
```

Parameter names and order follow the shared description above. See the [API index](README.md#avisynth-calls-and-builds) for host formats, array syntax, and audio/parity preservation.

`repairclip` is required; its previous, current, and next frames are requested internally. `mode=0` still repairs pixels; use `planes` to select which planes to process. The reference must match the source format and dimensions and provide the requested frames.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
filtered = neo_smo_Median(src, radius=[1])
return neo_smo_TemporalRepair(filtered, src, mode=[0], planes=[0])
```
