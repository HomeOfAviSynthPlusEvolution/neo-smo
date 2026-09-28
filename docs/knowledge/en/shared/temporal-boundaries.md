# Temporal windows and scene changes

[Back to the knowledge index](../README.md)

## Clip boundaries differ by filter

| Filter | Temporal range | When the full window is unavailable |
|---|---|---|
| Clense, TemporalRepair, DegrainMedian, FluxSmoothT/ST | n−1, n, n+1 | Copy the main clip's first and last frames |
| ForwardClense | n, n+1, n+2 | Copy the last two frames |
| BackwardClense | n−2, n−1, n | Copy the first two frames |
| TemporalMedian, CCD | n−r … n+r | Copy the first and last r frames; r is radius for TemporalMedian and temporal_radius for CCD |
| TemporalSoften | n−r … n+r | Shorten the window to frames that exist |
| TTempSmooth, Cnr4 | Parameter-defined temporal range | Clamp out-of-range frame indices to the first or last frame, potentially using boundary frames repeatedly |

For example, with TemporalMedian radius=3 on a five-frame clip, no frame has a complete window, so every frame remains unchanged. TemporalSoften can still process the same clip.

## Scene properties

`_SceneChangePrev` marks a cut between the current and previous frames; `_SceneChangeNext` marks a cut between the current and next frames. These are frame properties, not detection results automatically present on every clip. Whether a function generates, requires, or ignores them must be considered separately.

| Filter | Default | Enabling scene handling | Missing properties |
|---|---|---|---|
| TemporalMedian | Disabled | scenechange=True reads main-clip properties | Reports an error when a required property is missing |
| TemporalSoften | Disabled | −1 reads existing properties; 1–254 invokes SCDetect with threshold /255 | In −1 mode, missing properties provide no cut point |
| TTempSmooth | Automatic detection | scthresh>0 invokes SCDetect with threshold /100; negative values read existing properties | In read mode, missing properties provide no cut point |
| Cnr4 | Read existing properties | scenechange=True checks main-clip properties | Reports an error when a required property is missing |
| Other functions | No scene handling | No switch | Not read |

When pfclip is supplied, TTempSmooth uses scene properties from that guide path. Cnr4 uses properties from the main clip; setting ref does not change it to checking only ref.

TemporalMedian, TemporalSoften, and TTempSmooth restrict the candidate temporal range at cuts. Cnr4 keeps the window structure and replaces samples excluded by scene boundaries with the center. Different window shapes and replacement rules produce different output; these options are not interchangeable versions of one Boolean switch.

RGB cannot directly use TemporalSoften/TTempSmooth's automatic SCDetect path. To use existing properties on RGB, first detect scenes on suitable material, transfer the properties correctly, then select the mode that reads them.

Automatic detection invokes `misc.SCDetect` in VapourSynth. AviSynth+ includes normalized mean absolute luma difference detection with the same threshold units and property names. Automatic detection requires at least two frames; RGB still requires externally supplied scene properties.
