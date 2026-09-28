# Migrating from zsmooth

[Back to the API index](README.md)

neo-smo currently provides the 19 functions listed in this directory. Start by changing the namespace to `core.neo_smo`, then check each signature and format requirement. Do not assume that every upstream function and auxiliary option is registered.

```python
# A spatial median call in an existing script can be written as:
out = core.neo_smo.Median(src, radius=[1])
```

## Behavior to check

- The second input to both Repair and TemporalRepair is `repairclip`. Repair mode 0 copies; TemporalRepair mode 0 processes.
- `scalep` is off by default. Explicit thresholds can differ from omitted thresholds. Specify units when reusing parameters across bit depths, especially with DegrainMedian.
- Scene-change parameters have different meanings and defaults in TemporalMedian, TemporalSoften, TTempSmooth, and Cnr4. Follow the individual pages.
- F16 input is not an alias for F32 input. Platforms supporting native half-precision arithmetic can use it; other platforms compute in F32 and store the result as F16. DCTFilter is an explicit exception: even F16 input always uses F32 transforms before conversion back to F16. Bitwise equality across platforms is not a compatibility guarantee.
- Production paths allow FMA and ordinary rounding differences. Near thresholds or ordering boundaries, small differences can change a branch decision; assess their magnitude and their effect on actual output.
- There is no public `opt`, SIMD OFF, or KernelInfo interface, and no AviSynth interface yet. DCTFilter is available and requires exactly eight finite factors in `[0,1]`; see its [API page](dct-filter.md).

These documents describe the current implementation, without claiming bitwise reproduction of every historical zsmooth version. See the [knowledge index](../../knowledge/en/README.md) for algorithms and numerical conventions.
