# Repair

[Back to the API index](README.md)

```python
core.neo_smo.Repair(clip, repairclip, mode)
```

Supports the common formats. `clip` is the result to repair, and `repairclip` supplies the reference neighborhood from the same frame. Their formats and dimensions must match. `mode` is a required per-plane array with values 0–24, repeating its last element. Mode 0 copies the corresponding plane of `clip`. There is no `planes` parameter.

| mode | Behavior |
|---|---|
| 1, 11 | Clamp the main clip using the minimum and maximum of the reference 3×3 neighborhood |
| 2–4 | Sort the nine reference samples and use inner ranks as bounds |
| 5–9 | Select a directional range in the reference neighborhood using different costs |
| 10 | Select a reference value based on differences between the main value and reference samples |
| 12–14 | Sort the eight reference neighbors, select an inner range, then expand it to include the reference center |
| 15, 16 | Compute directional costs from the reference center, expand the selected range to include it, then clamp the main value |
| 17, 18 | Combine or select reference directional ranges and include the reference center |
| 19, 20 | Construct a symmetric range around the reference center using its smallest or second-smallest absolute difference from the eight neighbors |
| 21 | Construct a symmetric range around the reference center using the smallest of the opposing pairs' maximum differences |
| 22–24 | Use the main value as the distance origin in variants of modes 19–21, then clamp the reference center to that range |

Modes do not correspond one-to-one to RemoveGrain modes. In particular, Repair 11 is not a Gaussian weighted average. Spatial boundaries use mirrored sampling.

```python
filtered = core.neo_smo.RemoveGrain(src, mode=[2])
out = core.neo_smo.Repair(filtered, src, mode=[1])
```

This example limits the filtered result's deviation from the original clip; the reference is not simply another layer blended into the output. See [Spatial order statistics and clamping](../../knowledge/en/spatial.md).
