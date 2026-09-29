# neo-smo API

neo-smo ports zsmooth's spatial and temporal denoising and repair algorithms to C++ / Highway. Both VapourSynth and AviSynth+ are supported. VapourSynth uses `core.neo_smo` and plugin ID `org.neofilters.neo_smo`; all 21 AviSynth+ functions use the `neo_smo_` prefix.

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

VapourSynth examples on the following pages assume that `core` and `src` already exist. VapourSynth functions return a `VideoNode`; AviSynth functions return a clip. Signatures use readable Python call notation; `None` means omitting an optional argument, not passing a numeric value to the plugin.

See the [AviSynth+ interface](#avisynth-calls-and-builds) for loading, argument syntax, and host differences. Filter pages use Python notation for parameters shared by both hosts.

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
| Thresholded spatial and spatiotemporal smoothing | [Deen, MiniDeen](deen.md) | Fixed, adaptive and distance-weighted neighborhood means |
| Color denoising | [CCD](ccd.md) | Spatial sampling guided by multichannel differences |
| Temporal chroma denoising | [Cnr4](cnr4.md) | Blend using both luma and chroma differences |

## Common conventions

- Input must have a constant format and dimensions. Filters supporting the common formats accept Gray, RGB, and YUV with 8–16-bit integer, 16-bit floating-point (F16), or 32-bit floating-point (F32) samples. AviSynth+ does not provide F16 and accepts planar formats without alpha only. Exceptions: TTempSmooth does not accept F16; CCD does not accept Gray; Cnr4 accepts only 8–16-bit integer YUV.
- Omitting `planes` processes every plane. An explicitly supplied array must be nonempty; `planes=[]` is rejected. Plane indices 0/1/2 mean Y/U/V for YUV and R/G/B for RGB; Gray has only plane 0. Duplicate and out-of-range indices are errors.
- Per-plane arrays such as `radius` and `mode` generally repeat their last element for remaining planes. For example, `[1, 0]` becomes `[1, 0, 0]` on a three-plane clip. Each of Cnr4's three arrays must contain exactly three elements. DCTFilter's `factors` must contain exactly eight elements and describes frequencies, not planes. See the individual pages. A required `mode` array cannot be empty.
- Spatial radii use each plane's own pixel grid. A radius of 1 on a YUV420 chroma plane covers a different image area from a radius of 1 on the luma plane.
- Reference clips must match the main clip's dimensions, color family, sample type, bit depth, and chroma subsampling. The interface does not check for matching frame rates or frame counts; callers must ensure temporal alignment and valid frame ranges. AVS reports an error when a requested reference frame is out of range.
- Output retains the input format and dimensions. Unprocessed planes retain the main clip's content. There is no automatic bit-depth conversion, color-matrix conversion, or motion compensation.

## Threshold units

Functions with `scalep` default to `scalep=False`: **explicit values use the current format's units**. Setting it to `True` generally scales from an 8-bit reference: multiply by `2^(bits−8)` for integer samples, or divide by 255 for floating-point samples. DegrainMedian uses a different ratio; see its page.

Omitting a threshold can select an internally scaled default, so it need not equal explicitly passing the same number. For example, on 10-bit input, TemporalSoften uses 16 when `threshold` is omitted, but uses 4 for an explicit `[4]` without `scalep`.

See [Samples, thresholds, and precision](../../knowledge/en/shared/sample-and-precision.md). Scene-change parameters do not follow one common convention; see [Temporal windows and scene changes](../../knowledge/en/shared/temporal-boundaries.md).

## Migration and algorithms

- [Migrating from zsmooth](migration.md)
- [Algorithm knowledge index](../../knowledge/en/README.md)

## AviSynth calls and builds

In signatures, `Undefined()` means omitting an argument; normally leave that argument out of the call. It does not pass zero or an empty array. Each filter page provides AVS signatures and runnable examples.

AviSynth+ interface version 11 or later is required. Loading the same neo-smo plugin exposes all 21 functions with the `neo_smo_` prefix, including `neo_smo_Repair`, `neo_smo_TTempSmooth`, and `neo_smo_DCTFilter`.

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_Median(src, radius=[1], planes=[0])
```

### Arguments and formats

- After removing `neo_smo_`, names match the index above. Argument names, order, defaults, and algorithm behavior follow the individual filter pages.
- Write Python `True` / `False` as AVS `true` / `false`. `scalep`, `interlaced`, `norow`, `fp`, and boolean `scenechange` arguments accept booleans. TemporalSoften and Deen use integer `scenechange`.
- Numeric arrays use `[1, 2, 3]`; a scalar is also accepted as a one-element array. String arrays are not accepted. Empty arrays are normally rejected. CCD's `points` takes three integers: zero disables a set, nonzero enables it.
- `planes` uses Y/U/V or R/G/B indices 0/1/2; Gray has only plane 0. AVS storage order does not change this mapping. Only functions already exposing `planes` accept it.
- Planar Gray, YUV, and RGB support 8/10/12/14/16-bit integers and F32. The AVS interface has no F16 support and rejects packed formats and alpha. CCD rejects Gray; Cnr4 accepts integer YUV only.
- Reference clips must match the source format and dimensions. Frame rates and counts need not match, but callers must preserve temporal alignment. Requests beyond a reference clip's valid frame range report an error.
- Output preserves source audio, parity, and current-frame properties. Unprocessed planes retain source pixels. Each frame request owns its scratch space and supports `Prefetch`.

### Scene changes and resizing

Positive TemporalSoften/Deen `scenechange` and TTempSmooth `scthresh` values use built-in luma difference detection, without requiring VS's misc plugin. Threshold units remain `/255` and `/100`, respectively. TTempSmooth detects on `pfclip` when supplied. Automatic detection requires at least two frames and rejects RGB.

Existing scene information uses integer `_SceneChangePrev` and `_SceneChangeNext` properties. TemporalMedian with scene handling enabled, and Cnr4 with its default scene handling, report missing properties. Supply the properties or explicitly disable `scenechange`; TTempSmooth uses `scthresh=0` to disable detection.

For chroma-subsampled YUV, CCD and Cnr4 use AVS `ExtractY` / `BilinearResize` to obtain reference luma at chroma resolution. AVS and VS resizer rounding can differ, affecting neighbor acceptance and output pixels. These inputs therefore do not guarantee bitwise agreement between hosts.

### Building and testing

`NEO_SMO_BUILD_AVISYNTH` defaults to ON and is independent of `NEO_SMO_BUILD_VAPOURSYNTH`. Windows AVS builds require clang-cl or MSVC; MinGW builds must disable AVS. `NEO_SMO_AVS_SDK` selects local headers.

Enable `NEO_SMO_TEST_AVISYNTH` and set `NEO_SMO_AVISYNTH_RUNTIME` for host acceptance tests. `NEO_SMO_TEST_CROSS_HOST` additionally requires both interfaces, VapourSynth Python, and NumPy to compare output on identical inputs. Test subprocesses have time limits and suppress Windows crash dialogs.

Enable `NEO_SMO_TEST_VAPOURSYNTH` for real-host smoke tests covering all 21 functions. The selected Python environment must provide VapourSynth and NumPy; no zsmooth reference plugin is required.
