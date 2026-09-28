# Temporal denoising and repair

[Back to the knowledge index](README.md)

## Medians and one-sided clamping

Clense uses the previous, current, and next values p, c, and n to output `clamp(c,min(p,n),max(p,n))`, equivalent to a three-point median. TemporalMedian extends this idea to as many as 21 frames. If scene boundaries shorten the window to an even number of samples, it averages the two middle values.

ForwardClense/BackwardClense use the temporally nearer reference a and farther reference b. They construct an interval from a and the extrapolated value `2a−b`, then clamp c, using saturation at type limits. For example, a=100 and b=90 give an ideal range of 100–110: c=105 stays unchanged, while c=130 becomes 110.

## TemporalSoften and FluxSmooth denominators

TemporalSoften tests each temporal sample with `abs(sample-c) <= threshold`. A passing sample is retained; a failing sample is replaced by c. The sum is divided by the valid window length.

FluxSmooth first requires c to be a strict temporal extremum, then accumulates only neighbors that pass their thresholds and divides by the actual selected count, including the center. FluxSmoothST adds eight spatial neighbors but still first checks for a temporal extremum.

For three temporal values `[80,100,95]` and threshold 10, TemporalSoften's real-valued mean is `(100+100+95)/3≈98.33`. FluxSmoothT processes because the center 100 exceeds both neighbors, and its mean is `(100+95)/2=97.5`. Actual integer results also follow each filter's rounding rules. Equal thresholds do not imply equal output.

TemporalSoften's integer implementation normalizes with a fixed-point reciprocal. Ordinary floating-point division followed by rounding is not a sufficient reference for bitwise agreement.

## TemporalRepair

The main clip supplies c, while the reference supplies a range based on temporal changes. Mode 0 clamps c between the minimum and maximum of the same-position samples in three reference frames. It does not mean "do not process."

Mode 1 estimates positive and negative temporal changes from the eight reference neighbors. Mode 2 uses the largest change magnitude in the reference 3×3 neighborhood for a symmetric permissible range. Mode 3 takes the smaller of the maximum difference to the previous frame and the maximum difference to the next frame. Mode 4 uses only the trend of the three reference centers. These are different constraints, not a linear strength scale.

Repair can limit changes made by a preceding spatial filter, but it does not guarantee recovery of information already discarded. There is no motion compensation or scene-change parameter; the calling processing chain must account for behavior near cuts.

## DegrainMedian

Selects a direction from opposing points across the previous and next frames and opposing spatial points in the current frame. Mode 0 chooses the pair with the smallest endpoint difference. Other modes consider both the change from clamping the center and the pair's width.

Let p be the center clamped to a candidate interval, `e=abs(c-p)` its change, and d the interval width. For nonnegative samples, ignoring saturation details, costs are:

| mode | Cost |
|---|---|
| 0 | d |
| 1 | e+4d |
| 2 | e+2d |
| 3 | e+d |
| 4 | 2e+d |
| 5 | e |

Actual integer arithmetic saturates, and floating-point chroma also follows its nominal negative lower-bound convention. This table explains direction selection; it is not a bitwise replacement formula for the kernel. Every mode finally uses limit to bound the correction relative to the current value.

`interlaced=True` changes the vertical spatial span to two rows, while still using adjacent full frames in time. `norow=True` removes only the current frame's left/right pair, not all horizontal information.

## TTempSmooth

Checks guide samples outward from the current frame in both temporal directions, considering differences from the center and temporal continuity. Accepted samples receive weights based on temporal distance and difference. A failing condition stops further accumulation in that direction.

At distance i, the base temporal weight is 1 when `i < strength`, otherwise `1/(i-strength+2)`, then normalized over the complete window. thresh and mdiff determine whether extra difference weighting applies and where it decays. Larger strength retains high weights for more nearby samples.

Let c be the current value, S the selected weighted sum, and W the used weight:

```text
fp=True:  S + c*(1-W)
fp=False: S/W
```

S and W already include the center. Thus fp=True assigns unused weight to the center; False renormalizes the remaining valid samples to fill the total weight. pfclip is used only for decisions and weights; accumulated output colors come from the main clip.

See [Temporal windows and scene changes](shared/temporal-boundaries.md) for boundaries and scene parameters.
