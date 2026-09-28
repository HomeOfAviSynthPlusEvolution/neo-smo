# Cnr4

[Back to the API index](README.md)

```python
core.neo_smo.Cnr4(clip, mode="oxx", radius=2, sense=[35, 47, 47],
                 str=[192, 255, 255], pow=[1.0, 1.0, 1.0],
                 tmode=0, wmode=0, scenechange=True, ref=None)
```

Accepts only 8–16-bit integer YUV. Only U/V are modified; Y retains the main clip's content. There is no `planes` parameter, and `mode` does not enable or disable planes.

| Parameter | Range | Meaning |
|---|---|---|
| mode | Exactly three lowercase `o` or `x` characters | Selects the lookup curve for Y/U/V differences |
| radius | 1–10 | Temporal radius |
| sense | Exactly three elements, each −1–255 | Difference-curve scale; −1 retains that position's default |
| str | Exactly three elements, each −1–255 | Lookup amplitude and valid table range; −1 retains that position's default |
| pow | Exactly three finite elements, each ≥0 | Difference-curve exponent parameter |
| tmode | 0–4 | Single-pass or multipass processing; see below |
| wmode | 0–3 | Temporal-distance weighting |
| scenechange | Default True | Uses existing scene-change properties; does not invoke a detector automatically |
| ref | Defaults to the main clip | Supplies Y/U/V difference information; must match the main clip's format and dimensions |

| tmode | Behavior |
|---|---|
| 0 | Single-pass blending directly from temporal neighbors |
| 1 | Preprocess inward from both ends of the window with local radius 1, updating reference results |
| 2 | Like 1, also updating source pixels used in subsequent processing |
| 3 | Gradually increase the preprocessing radius, updating reference results |
| 4 | Like 3, also updating source pixels |

wmode 0 uses equal weights; 1 uses square-root distance weights; 2 uses sine distance weights; 3 uses inverse-distance decay. This is independent of the curve selected by mode.

The default `scenechange=True` requires `_SceneChangePrev` and `_SceneChangeNext` properties on the main clip. Disable it explicitly if no scene detection has been performed. When enabled, main-clip properties determine the window boundaries, and excluded samples are replaced by the center. Clip boundaries repeat the first or last frame. Reference luma must be aligned to chroma dimensions; the reference path is constructed using `std`/`resize`.

```python
out = core.neo_smo.Cnr4(src, scenechange=False)
# With the misc plugin installed, generate scene properties first:
sc = core.misc.SCDetect(src, threshold=0.1)
out = core.neo_smo.Cnr4(sc)
```

See [Color denoising](../../knowledge/en/chroma.md) for the curves and the role of the reference clip.
