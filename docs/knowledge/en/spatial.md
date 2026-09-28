# Spatial order statistics and clamping

[Back to the knowledge index](README.md) · [Spatial median API](../../api/en/spatial-median.md)

Let c be the current value. `clamp(c, L, U)` limits c to `[L,U]`: unchanged inside the interval, L below it, and U above it.

## Median and InterQuartileMean

Median sorts the square neighborhood, including its center, and returns the middle value. It does not average every sample, so an isolated large value does not affect a sum as it would in a mean.

InterQuartileMean uses the middle 50% of the sorted samples with appropriate weights. For nine ascending 3×3 values x0…x8, its ideal real-valued formula is:

```text
(x3 + x4 + x5 + 0.75*(x2+x6)) / 4.5
```

For 5×5, the fully included samples are x7…x17, with fractional boundary samples x6/x18, divided by 12.5. For 7×7, they are x13…x35, with x12/x36 at the boundaries, divided by 24.5.

The integer implementation rounds the two 3/4-weight boundary terms before integer normalization. This is not equivalent to rounding the ideal formula only once at the end. Nor is it simply discarding the top and bottom quarters and equally averaging an integer number of remaining samples.

## SmartMedian

Exclude the center and sort the remaining 8, 24, or 48 neighbors. Let a and b be the two middle values, and m the mean of all those neighbors. The flatness measure is:

```text
v = 13 * sqrt((a-m)^2 + (b-m)^2)
out = clamp(c, min(a,b), max(a,b))   if v <= threshold
out = c                           otherwise
```

Integer arithmetic also rounds the mean and v as specified by the implementation. This v is not the standard deviation of all neighbors, and the threshold does not measure the difference between the center and the median.

For example, if all eight neighbors equal 100 and the center is 200, then a=b=m=100 and v=0, allowing correction to 100. A large center deviation does not prevent correction because the test examines the neighborhood itself.

## RemoveGrain directional costs

The four opposing neighbor pairs are horizontal, vertical, and the two diagonals. Each defines `[L,U]`. Let `p=clamp(c,L,U)`, `e=abs(c-p)`, and `d=U-L`.

| Mode | Direction-selection cost |
|---|---|
| 5 | e |
| 6 | 2e+d, saturated to the implementation's upper limit |
| 7 | e+d |
| 8 | e+2d, limited to the implementation's range |
| 9 | d |

The output is p from the selected direction. Equal costs have a fixed directional priority; changing this to whichever direction happens to be visited first need not preserve behavior. Mode 18 instead uses the maximum difference between the center and either endpoint of the pair.

Modes 23/24 do not simply select one direction and clamp. For each direction, let `h=c-U`, `l=L-c`, and `d=U-L`. The downward candidate is `min(h,d)` in mode 23 and `min(h,d-h)` in mode 24; upward candidates are analogous. Take the largest nonnegative correction in each sense across all directions, then output `c−downward correction+upward correction`.

See the [RemoveGrain API](../../api/en/remove-grain.md) for all modes. Modes 13–16 process selected even or odd rows rather than denoising every pixel.

## The role of Repair's reference

Repair takes the neighborhood from repairclip and the value to constrain from clip. Mode 1 constructs a range from the minimum and maximum of the nine reference samples. If that range is 90–110 and the main value is 120, the output is 110.

Modes 2–4 sort the nine reference samples and choose inner bounds. Modes 12–14 first sort the eight neighbors, then expand the range to include the reference center. The latter therefore always admit the reference center, while the former need not.

Modes 5–9 use directional ranges that include the reference center, with costs measured relative to the main value. Modes 15/16 select a direction from the relationship between the reference center and opposing neighbors, then expand its range to include the reference center. Modes 22–24 further exchange the roles of the distance origin and the value being clamped. Repair is therefore not just a blend of RemoveGrain's output with the original clip.

## VerticalCleaner and spatial boundaries

VerticalCleaner mode 1 is a vertical three-point median. Mode 2 estimates local trends from two rows above and below, constructs permissible bounds, then clamps the center; it is not a five-point median. The first and last one or two rows remain unprocessed.

Most other spatial neighborhoods use mirrored coordinates: −1 maps to 1, and coordinate N for a length-N plane maps to N−2. For length 1, or radii larger than very small planes, coordinates are additionally clamped to the valid range. This is not infinite periodic reflection. Internal memory padding serves SIMD access only; it does not add valid image pixels.
