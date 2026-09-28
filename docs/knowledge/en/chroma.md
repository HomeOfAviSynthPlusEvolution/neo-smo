# Color denoising: CCD and Cnr4

[Back to the knowledge index](README.md) · [CCD API](../../api/en/ccd.md) · [Cnr4 API](../../api/en/cnr4.md)

## CCD: selecting spatial points by color differences

CCD checks sparse spatial sample points. RGB distance is the sum of squared differences in three channels. YUV uses `4*ΔY²+ΔU²+ΔV²`, giving luma changes greater influence over whether samples are blended.

The original point sets, expressed as `(x,y)` offsets from the center, are:

| Group | Points |
|---|---|
| Near, 4 points | `(±4,±4)` |
| Middle, 8 points | `(±8,±8)`, `(0,±8)`, `(±8,0)` |
| Far, 12 points | x=−12,−4,4,12 at y=±12; x=±12 at y=±4 |

YUV coordinates are reduced according to subsampling, then scaled by scale and rounded. Out-of-range positions use mirrored coordinates. Reference luma is aligned to chroma dimensions using Bilinear. RGB processes all three output planes; YUV retains main-clip luma and updates only chroma.

For format maximum M, or 1 for floating-point input, threshold parameter t gives the internal threshold:

```text
T = (t*M)^2 / (255^2*3)
```

A point is accepted only when its distance is strictly below the applicable threshold. Integer paths also round distance terms and thresholds as specified by the implementation. This is not a test that "each channel differs by less than t."

With temporal_radius=0, only the current reference frame is used. Increasing the radius includes distances from candidate positions in other temporal frames to the **center of the current reference frame**, with temporal weighting. Forward and backward weighting formulas are not identical; this cannot be rewritten as a simple symmetric temporal mean.

Regardless of how many frames contribute to the decision, an accepted point adds the pixel at that spatial position in the current main frame to the output average. The center always participates. Enlarging the decision's temporal radius therefore differs from averaging more output frames.

ref can supply more stable difference information from a preprocessed result, while output colors still come from the main clip. Excessive reference blurring can also make genuinely different regions more likely to pass the threshold.

## Cnr4: controlling temporal blending with Y/U/V differences

Cnr4 processes integer YUV only. It compares reference luma and corresponding chroma differences between the current frame and temporal neighbors, then combines lookup results into blending weights. Larger weights give neighboring frames more influence on output chroma. Luma affects this blending even though output luma itself is unchanged.

The three mode characters select Y/U/V curves; they are not processing switches. For 8-bit difference index l and sense s:

```text
o: base = (1 + cos(l*l*pi/(s*s))) / 2
x: base = (1 + cos(l*pi/s)) / 2
table[l] = floor(clamp(str * base^(1/pow), 0, 255))
```

Only the specified `0…str` index range is populated; the rest is zero. sense=0 and pow=0 have special handling and cannot be substituted directly into these division formulas. Because parameters affect both curve shape and valid index range, increasing sense or str does not necessarily increase every difference weight monotonically.

tmode 0 directly processes neighbors in one pass. Modes 1–4 preprocess inward from both ends of the temporal window. They differ in whether the local radius grows and whether source pixels used later are also updated. They do not use the previous output frame as global recursive state, but they incur multipass processing within each window.

For distance k, with k=1 for the nearest frame, and radius r, temporal weights are:

| wmode | Weight |
|---|---|
| 0 | 1 |
| 1 | `sqrt((r-k+1)/(2r))` |
| 2 | `sin((r-k+2)/(2r))`, in radians |
| 3 | `1/(k+1)` |

These temporal weights modify the difference weights; they do not multiply the entire result by a single luma coefficient. Main-clip scene properties determine which window samples are replaced by the center; see [Temporal windows and scene changes](shared/temporal-boundaries.md).
