# Reference evidence and behavior choices

Specification: `DEEN-REF-001`. This specification family defines a new behavior contract, not an exact description of the original DLL.

## Materials and evidence limits

The documents were informed by seven user-supplied analyses under local `deen/out/spec/`: `c2d.md`, `c3d.md`, `w2d.md`, `w3d.md`, `a2d.md`, `a3d.md`, and `_aux.md`, plus the relevant formulas in `deen/src/deen.cpp`. They describe mode structure, threshold comparisons, distance falloff, temporal weights, and surrounding behavior.

Those materials report assembly or runtime verification. Their claims are secondary evidence, not independently established facts about the original binary. Conformance to these specifications does not require the supplied materials or the original DLL.

The mathematical definitions and fixtures in this specification family are the conformance authority. Reconstructed code and original-DLL output are not correctness oracles. Deen, MiniDeen, and different historical versions must not be treated as one pixel-exact contract.

## Mode structure and new rules

| Area | Structure or behavior reported by the supplied materials | This specification |
| --- | --- | --- |
| c | Replace rejected samples with the center; fixed tap count | Preserve structure; correct normalization and final rounding |
| w | Distance weights, rejection substitution, temporal weights 1:2:1 | Preserve structure and normalization; permit bounded U8 fixed-point rounding; no unconditional final +1 |
| a | Distance thresholds; average accepted samples | Preserve structure; count center once; no byte threshold quantization or approximate count reciprocals |
| Comparison | Inclusive threshold against the current center | Preserve; define high-bit-depth and floating difference units explicitly |
| Borders | Some paths read padding; borderfix has differing frame coordinates | Replicate borders consistently in every frame and filter all visible samples |
| Scene difference | Equally weighted plane means; some paths ignore stride | Preserve metric structure; use visible samples and actual strides |
| Temporal fallback | Spatial endpoints; scene behavior depends on request history | Spatial endpoints; either adjacent cut forces 2D without historical state |
| Numerics | Fixed-point reciprocals, truncation, saturation, platform floating details | Mathematical definitions; bounded FMA/rounding differences without systematic bias |
| Parameters | Byte narrowing and permissive or unsafe cases | Explicit domains and rejection; min restricted to [0,1] |
| fcf / borderfix | Legacy surrounding interface | Excluded from this contract; no claim of complete script compatibility |

Border replication, full-range integer threshold scaling, floating extensions, an explicit scene switch, and RGB/planes mappings are deliberate choices, not claims about the original author's implementation. Reusability comes from explicit contracts rather than emulation of the old binary.
