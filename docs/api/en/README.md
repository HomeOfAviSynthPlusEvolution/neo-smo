# neo-smo API

neo-smo ports zsmooth's spatial and temporal denoising and repair algorithms to C++ / Highway. The public interface is VapourSynth, under `core.neo_smo`, with plugin ID `org.neofilters.neo_smo`. This documentation describes the current behavior, including DCTFilter. An AviSynth interface is not yet available.

## Quick start

```python
import vapoursynth as vs
core = vs.core
# If the plugin is not installed in an autoload directory, supply its actual path:
# core.std.LoadPlugin(path=r"C:\plugins\neo-smo.dll")

src = core.std.BlankClip(width=640, height=480, format=vs.YUV420P10, length=24)
out = core.neo_smo.Median(src, radius=[1], planes=[0])
out.set_output()
```

Examples on the following pages assume that `core` and `src` already exist. Functions return a `VideoNode`. Signatures use readable Python call notation; `None` means omitting an optional argument, not passing a numeric value to the plugin.

## Function index

| Category | Functions | Purpose |
|---|---|---|
| Spatial order statistics | [Median, InterQuartileMean, SmartMedian](spatial-median.md) | Median, interquartile mean, and conditional median correction |
| Vertical processing | [VerticalCleaner](vertical-cleaner.md) | Vertical median or extrapolation-based clamping |
| Spatial denoising | [RemoveGrain](remove-grain.md) | 24 neighborhood processing modes |
| Block frequency processing | [DCTFilter](dct-filter.md) | Coefficient weighting in fixed 8×8 DCT blocks |
| Spatial repair | [Repair](repair.md) | Constrain changes using a reference clip's spatial neighborhood |
| Three-frame processing | [Clense, ForwardClense, BackwardClense](clense.md) | Temporal clamping and one-sided extrapolation |
| Temporal order statistics and averaging | [TemporalMedian, TemporalSoften](temporal.md) | Temporal median or thresholded averaging |
| Temporal repair | [TemporalRepair](temporal-repair.md) | Temporal or spatiotemporal constraints from a reference clip |
| Spatiotemporal directional processing | [DegrainMedian](degrain-median.md) | Select a direction and limit each correction |
| Temporal extrema processing | [FluxSmoothT, FluxSmoothST](flux-smooth.md) | Thresholded averaging at temporal extrema |
| Weighted temporal smoothing | [TTempSmooth](ttempsmooth.md) | Weight samples by differences and temporal distance |
| Color denoising | [CCD](ccd.md) | Spatial sampling guided by multichannel differences |
| Temporal chroma denoising | [Cnr4](cnr4.md) | Blend using both luma and chroma differences |

## Common conventions

- Input must have a constant format and dimensions. Filters supporting the common formats accept Gray, RGB, and YUV with 8–16-bit integer, 16-bit floating-point (F16), or 32-bit floating-point (F32) samples. Exceptions: TTempSmooth does not accept F16; CCD does not accept Gray; Cnr4 accepts only 8–16-bit integer YUV.
- Omitting `planes` processes every plane. An explicitly supplied array must be nonempty: the host rejects `[]`. Plane indices 0/1/2 mean Y/U/V for YUV and R/G/B for RGB; Gray has only plane 0. Duplicate and out-of-range indices are errors.
- Per-plane arrays such as `radius` and `mode` generally repeat their last element for remaining planes. For example, `[1, 0]` becomes `[1, 0, 0]` on a three-plane clip. Each of Cnr4's three arrays must contain exactly three elements. DCTFilter's `factors` must contain exactly eight elements and describes frequencies, not planes. See the individual pages. A required `mode` array cannot be empty.
- Spatial radii use each plane's own pixel grid. A radius of 1 on a YUV420 chroma plane covers a different image area from a radius of 1 on the luma plane.
- Reference clips must match the main clip's dimensions, color family, sample type, bit depth, and chroma subsampling. The interface does not check for matching frame rates or frame counts; callers must ensure temporal alignment and valid frame ranges.
- Output retains the input format and dimensions. Unprocessed planes retain the main clip's content. There is no automatic bit-depth conversion, color-matrix conversion, or motion compensation.

## Threshold units

Functions with `scalep` default to `scalep=False`: **explicit values use the current format's units**. Setting it to `True` generally scales from an 8-bit reference: multiply by `2^(bits−8)` for integer samples, or divide by 255 for floating-point samples. DegrainMedian uses a different ratio; see its page.

Omitting a threshold can select an internally scaled default, so it need not equal explicitly passing the same number. For example, on 10-bit input, TemporalSoften uses 16 when `threshold` is omitted, but uses 4 for an explicit `[4]` without `scalep`.

See [Samples, thresholds, and precision](../../knowledge/en/shared/sample-and-precision.md). Scene-change parameters do not follow one common convention; see [Temporal windows and scene changes](../../knowledge/en/shared/temporal-boundaries.md).

## Migration and algorithms

- [Migrating from zsmooth](migration.md)
- [Algorithm knowledge index](../../knowledge/en/README.md)
