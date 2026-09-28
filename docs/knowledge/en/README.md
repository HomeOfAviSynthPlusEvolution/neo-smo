# neo-smo algorithm knowledge

This directory explains how the algorithms select samples, decide whether to modify pixels, and handle numerical and boundary conditions. For signatures, defaults, and ready-to-use examples, see the [API index](../../api/en/README.md).

| Topic | Contents |
|---|---|
| [Spatial order statistics and clamping](spatial.md) | Median, InterQuartileMean, SmartMedian, VerticalCleaner, RemoveGrain, Repair |
| [Temporal denoising and repair](temporal.md) | The Clense family, TemporalMedian/Soften/Repair, DegrainMedian, FluxSmooth, TTempSmooth |
| [Block DCT and frequency weighting](dct-filter.md) | DCTFilter's 8×8 blocks, outer-product weights, edge extension, and normalization |
| [Color denoising](chroma.md) | CCD's multichannel spatial sampling and Cnr4's temporal blending |
| [Samples, thresholds, and precision](shared/sample-and-precision.md) | Planes, bit depth, scalep, F16/F32, FMA |
| [Temporal windows and scene changes](shared/temporal-boundaries.md) | Clip boundaries, scene properties, and each filter's window rules |

## Understanding the algorithms

Spatial filters sample neighbors in the same frame; temporal filters sample the same coordinates in neighboring frames. Increasing the temporal radius does not expand a motion-search range: these neo-smo functions do not estimate motion. With strong motion, temporal samples can belong to different objects. Thresholds and scene boundaries help reduce unsuitable blending but do not replace motion compensation.

Medians mainly suppress isolated extreme values. Averaging can reduce random variation but also loses detail. Clamping only pulls values that exceed the permitted range back to its boundary. Filters such as Repair construct that range from a reference clip, serving a different purpose from directly blurring the main clip.

The formulas explain algorithm structure; integer saturation, integer rounding, and floating-point rounding still follow the implementation. Mode numbers, radii, and thresholds do not form a common "denoising strength" scale across filters.
