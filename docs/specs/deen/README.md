# Deen specifications

Specification family: `DEEN-001`.

The goal is a mathematically consistent, memory-safe, request-order-independent filter based on the algorithmic structure of Deen's six modes. The original DLL, reconstructed code, and MiniDeen are not pixel-exact correctness oracles. Black-box matching is not required. Legacy integer approximations, out-of-bounds reads, and state-dependent behavior are not compatibility requirements.

Correctness means conformance to the explicit mathematical contracts below. It does not imply that border extension, scene metrics, or other design choices have only one reasonable definition. The distinction between reported legacy behavior and deliberate choices is recorded in the [reference mapping](reference-mapping.md).

## Organization

| Document | Responsibility |
| --- | --- |
| [kernel-domain.md](kernel-domain.md) | Sample domains, neighborhoods, threshold units, and numerics |
| [kernel-c.md](kernel-c.md) | c2d/c3d: threshold substitution with a fixed denominator |
| [kernel-w.md](kernel-w.md) | w2d/w3d: distance weights, substitution, and normalization |
| [kernel-a.md](kernel-a.md) | a2d/a3d: distance thresholds and accepted-sample averaging |
| [kernel-temporal.md](kernel-temporal.md) | Three-frame selection, scene detection, and spatial fallback |
| [plugin.md](plugin.md) | neo-smo interface contract, host adaptation, and errors |
| [acceptance.md](acceptance.md) | Independent reference calculations and acceptance criteria |
| [reference-mapping.md](reference-mapping.md) | Evidence, limitations, and deliberate differences |

The kernel specifications do not depend on AviSynth, VapourSynth, DualSynth2, or Highway. They describe values and operators without prescribing a C/C++ ABI, allocator, or SIMD instruction set. Other projects may reuse these definitions directly. Their interface names, parameter containers, and frame lifetimes belong in separate binding specifications, not in the core equations.

Scope includes all six modes, planar integer and floating-point samples, spatial boundaries, and temporal fallback. It excludes the legacy `fcf` configuration format, memory-layout emulation, legacy rounding emulation, and MiniDeen.
