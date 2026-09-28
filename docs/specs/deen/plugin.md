# neo-smo integration contract

Specification: `DEEN-PLUGIN-001`. This document defines the neo-smo host-binding requirements. Other projects reusing the mathematics may define separate host bindings.

## Parameters

Function name: `Deen`. VS namespace: `neo_smo`. AVS name: `neo_smo_Deen`. Parameters are ordered as follows; all except clip are optional.

| Parameter | Default | Domain and behavior |
| --- | --- | --- |
| clip | Required | Fixed-format planar GRAY/YUV/RGB video with positive dimensions and frame count |
| mode | `"c3d"` | Exactly c2d/c3d/w2d/w3d/a2d/a3d, in lowercase |
| rad | 1 | Integer; 1–7 for 2D, 1–4 for 3D |
| thrY | 7 | Finite [0,255]; spatial threshold for Y/GRAY or all RGB planes |
| thrUV | 9 | Finite [0,255]; spatial threshold for YUV chroma |
| tthY | 4 | Finite [0,255]; temporal threshold for Y/GRAY or all RGB planes |
| tthUV | 6 | Finite [0,255]; temporal threshold for YUV chroma |
| min | 0.5 | Finite [0,1]; corner weight ratio for w, corner threshold ratio for a |
| scd | 9 | Finite and nonnegative; scene difference threshold in 8-bit units |
| scenechange | true | Enables scene detection; false does not disable endpoint fallback |
| planes | All | Array of actual color-plane indices; omission or an empty array selects all |

VS uses an integer array; AVS uses its supported array representation for planes. Explicit indices must exist and must not repeat. Interpret indices in the host's logical color-plane order, not physical memory order. Support one or three color planes only. Reject alpha-bearing, packed, variable-format, or variable-size inputs rather than converting implicitly.

Each host binding must support the following sample formats:

| Binding | Unsigned integer bit depths | Floating-point formats |
| --- | --- | --- |
| VapourSynth | Every integer bit depth from 8 through 16; 8-bit samples occupy one byte, 9–16-bit samples occupy two bytes | F16 and F32 |
| AviSynth+ | 8, 10, 12, 14, and 16 | F32; no F16 |

These rows apply to planar GRAY, RGB, and YUV formats representable by each host. Preserve the input's valid chroma subsampling and actual plane dimensions; do not resample. Reject sample formats outside the binding's row. In particular, AVS must not reinterpret 16-bit integers as F16 or silently convert unsupported input. Other projects may define a different binding subset while retaining the core's 8–16-bit integer and F16/F32 contracts.

Validate every parameter at creation, including parameters unused by the selected mode or color family. A valid min has no effect in c modes. Valid temporal thresholds, scd, and scenechange in 2D modes do not trigger temporal reads. Check host integer ranges before narrowing; never truncate negative or large values into bytes. Reject unknown modes, invalid argument types, and non-finite parameters.

There is no fcf, borderfix, legacy-rounding option, or opt/SIMD OFF parameter. Border replication is fixed by the algorithm. Do not silently accept and ignore legacy arguments. Any future per-frame configuration file requires a separate reviewed extension; parser state must not enter pixel kernels.

## Execution and lifetime

Convert host parameters into immutable configuration at creation. The core receives sample views, thresholds/weights, and the effective mode. Frame requests, error reporting, planes selection, and properties belong to the binding. The scene metric's mathematical definition belongs to the reusable algorithm layer.

Request only n for 2D modes and 3D endpoints. Other 3D frames request n−1, n, and n+1 and use the [temporal rules](kernel-temporal.md). Deduplicate identical node/frame dependencies. Keep inputs alive while used. Writable workspaces and output storage must not be shared between concurrent requests.

Output format, dimensions, frame count, and frame rate match clip. Copy all properties from clip[n] without adding diagnostics or automatically rewriting scene properties. Preserve AVS audio and other video information. Unselected planes pass through bitwise; every visible pixel of selected planes must be initialized. Padding must not affect the result.

Report allocation failures, required dependency failures, changed input descriptors, and numerical failures through the host's frame-error mechanism. Release all acquired resources and never publish partial output. Creation failures must not leak node/frame references. Production paths may use Highway dynamic dispatch and FMA. A test reference implementation is not an additional public scalar backend.
