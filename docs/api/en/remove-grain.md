# RemoveGrain

[Back to the API index](README.md)

```python
core.neo_smo.RemoveGrain(clip, mode)
```

Supports the common formats. `mode` is required: a per-plane array with values 0–24, repeating its last element. There is no `planes` parameter; mode 0 preserves the corresponding plane. Other modes use a 3×3 neighborhood with mirrored boundaries.

| mode | Behavior |
|---|---|
| 0 | Copy |
| 1 | Clamp between the minimum and maximum of the eight surrounding samples |
| 2–4 | Sort the eight neighbors and progressively narrow the bounds; mode 4 uses the two middle neighbors |
| 5–9 | Choose the least-cost direction from four opposing neighbor pairs and clamp; see the knowledge page for costs |
| 10 | Select the neighbor closest to the center |
| 11, 12 | Identical 3×3 weighted average: center weight 4, axial neighbors 2, corners 1, total weight 16 |
| 13, 14 | Choose the vertical or diagonal pair with the smallest endpoint difference and take its mean; 13 processes even rows, 14 odd rows |
| 15, 16 | Compute a weighted average of the rows above and below, then clamp using the same best direction; 15 processes even rows, 16 odd rows |
| 17 | Construct bounds from the inner endpoints of the four opposing-pair ranges |
| 18 | Select a direction by the maximum difference between the center and either endpoint, then clamp |
| 19 | Average the eight neighbors, excluding the center |
| 20 | Average all nine samples, including the center |
| 21, 22 | Clamp the center to the range of the four opposing-pair means; integer rounding differs, while floating-point behavior is identical |
| 23, 24 | Limit upward and downward corrections using directional interval widths; mode 24 also weakens the correction as the distance beyond the interval increases |

Row numbers start at 0. Modes 13–16 copy unselected rows and do not read field-order properties. Mode numbers are not a common strength scale: mode 24 does not mean "24 times stronger than mode 1."

```python
out = core.neo_smo.RemoveGrain(src, mode=[2, 0, 0])
```

See [Spatial order statistics and clamping](../../knowledge/en/spatial.md) for formulas and direction selection.
