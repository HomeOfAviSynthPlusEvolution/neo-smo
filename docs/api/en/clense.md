# Clense, ForwardClense, BackwardClense

[Back to the API index](README.md)

```python
core.neo_smo.Clense(clip, previous=None, next=None, planes=None)
core.neo_smo.ForwardClense(clip, planes=None)
core.neo_smo.BackwardClense(clip, planes=None)
```

Supports the common formats. Omitting `planes` processes every plane.

| Function | References when processing frame n | Temporal boundaries left unchanged |
|---|---|---|
| Clense | `previous[n-1]` and `next[n+1]`; both clips default to `clip` | First and last frame |
| ForwardClense | `clip[n+1]` and `clip[n+2]` | Last two frames |
| BackwardClense | `clip[n-1]` and `clip[n-2]` | First two frames |

Optional reference clips must satisfy the common matching requirements. `previous` and `next` are video nodes; the interface already applies the ±1 frame offsets. Do not shift them again to obtain the same intended alignment.

Clense clamps the current value between the previous and next reference values, equivalent to a median of three. ForwardClense and BackwardClense use the near and far reference values for one-sided extrapolation and clamping; **they are not three-value medians**. None of the three checks scene-change properties.

```python
out = core.neo_smo.Clense(src, planes=[0])
out = core.neo_smo.ForwardClense(src)
```

See [Temporal denoising and repair](../../knowledge/en/temporal.md) for the formulas.
