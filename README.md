# neo-smo

**English** | [简体中文](README.zh-CN.md) | [日本語](README.ja.md)

neo-smo is an AviSynth+ / VapourSynth video filter plugin ported directly from [zsmooth](https://github.com/adworacz/zsmooth). Based on its Zig source, it is implemented and optimized in C++17 with Google Highway. It provides spatial and temporal denoising, neighborhood repair, color denoising, and block DCT filtering. VapourSynth uses the `core.neo_smo` namespace; AviSynth+ uses the `neo_smo_` function prefix.

## Origins

zsmooth was developed by Austin Dworaczyk Wiltshire (adworacz). neo-smo directly inherits its algorithms and interfaces. zsmooth itself builds on algorithms and implementations developed by the video filter community. Its documentation lists sources and references including:

- [AviSynth RgTools](https://github.com/pinterf/RgTools) and [VapourSynth RemoveGrain](https://github.com/vapoursynth/vs-removegrain): historical implementations and behavior references for the RemoveGrain family.
- [VapourSynth TemporalSoften](https://github.com/dubhater/vapoursynth-temporalsoften2), [VapourSynth TemporalMedian](https://github.com/dubhater/vapoursynth-temporalmedian), and [Neo Temporal Median](https://github.com/HomeOfAviSynthPlusEvolution/neo_TMedian): temporal smoothing and temporal median filters.
- [VapourSynth FluxSmooth](https://github.com/dubhater/vapoursynth-fluxsmooth): the FluxSmooth family.
- [Dogway's AviSynth Scripts](https://github.com/Dogway/Avisynth-Scripts): algorithm ideas for `ex_median`, interquartile mean, and SmartMedian.
- [End-of-Eternity's CCD](https://github.com/End-of-Eternity/vs-ccd) and [vs-jetpack's CCD](https://github.com/Jaded-Encoding-Thaumaturgy/vs-jetpack/blob/e0f47d86930150fd0bf92b0845ccc2b0491f7807/vsdenoise/ccd.py#L95): the main implementation references for zsmooth CCD. CCD was originally written by Sergey Stolyarevsky for VirtualDub.
- [Cnr2](http://avisynth.nl/index.php/Cnr2): one of the algorithm sources for zsmooth Cnr4.

Thanks to the authors and contributors of zsmooth and these upstream projects.

## Implementation

The computation core handles sampling, sorting, weighting, and transforms on image planes. The host layer manages parameters, frame requests, properties, and output allocation. Google Highway provides cross-platform SIMD and runtime instruction-set selection; DualSynth2 provides host integration infrastructure. DCTFilter uses neo-libdct's fixed 8×8 transforms and requires no FFTW runtime.

The port preserves zsmooth's public function names and parameter structure while fixing identified issues and optimizing computation paths. Production builds allow FMA and floating-point rounding differences; bitwise agreement across CPUs or historical zsmooth versions is not guaranteed. See the API documentation for formats, threshold units, and boundary behavior.

One plugin provides both VapourSynth and AviSynth+ interfaces; hosts can also be selected individually at build time. AviSynth+ requires interface version 11 or later.

## Supported operations

| Category | Functions and purpose |
|---|---|
| Spatial order statistics and cleanup | `Median`, `InterQuartileMean`, `SmartMedian`, `VerticalCleaner`: median, interquartile mean, conditional median correction, and vertical cleanup. |
| Neighborhood denoising and repair | `RemoveGrain`, `Repair`: spatial processing modes and clamping against a reference neighborhood. |
| Temporal clamping and repair | `Clense`, `ForwardClense`, `BackwardClense`, `TemporalRepair`: constrain current pixels using adjacent frames. |
| Temporal and spatiotemporal denoising | `TemporalMedian`, `TemporalSoften`, `DegrainMedian`, `FluxSmoothT`, `FluxSmoothST`, `TTempSmooth`: temporal order statistics, thresholded averaging, direction selection, and weighted smoothing. |
| Color denoising | `CCD`, `Cnr4`: spatial sampling guided by color differences and temporal chroma denoising. |
| Block frequency filtering | `DCTFilter`: DCT coefficient weighting in non-overlapping 8×8 blocks. |

There are 19 functions. Common support covers constant-format, constant-size planar Gray, YUV, and RGB with 8–16-bit integer, 16-bit floating-point (F16), and 32-bit floating-point (F32) samples. F16 is available only in VapourSynth; AviSynth+ accepts planar formats without alpha. TTempSmooth does not accept F16, CCD does not accept Gray, and Cnr4 accepts integer YUV only.

Output preserves the input format and dimensions. Functions exposing `planes` process all planes when it is omitted; `planes=[0]` processes only the first plane. Explicit empty arrays are rejected. Some functions use `mode` or thresholds to control plane processing, so `planes` cannot be passed to every function. These filters do not estimate motion vectors or perform motion compensation.

## Documentation and use

- [API reference](docs/api/en/README.md): function signatures, parameters, defaults, and examples.
- [Algorithm knowledge base](docs/knowledge/en/README.md): computation, numerical rules, and boundary behavior.
- [Migrating from zsmooth](docs/api/en/migration.md): namespaces, threshold units, scene changes, and precision conventions.

Place the plugin in VapourSynth's autoload directory or load it explicitly with `LoadPlugin`. The example uses a Windows filename; on other platforms, substitute the actual plugin path. The plugin ID is `org.neofilters.neo_smo`.

```python
import vapoursynth as vs

core = vs.core
core.std.LoadPlugin(path="/path/to/neo-smo.dll")

clip = core.std.BlankClip(width=640, height=480, format=vs.YUV420P10, length=24)
output = core.neo_smo.Median(clip, radius=[1], planes=[0])
output.set_output()
```

This minimal example uses a synthetic clip to demonstrate luma median filtering. Replace `clip` with your source for actual use. Functions with `scalep` interpret explicit thresholds in native units by default; omitted arguments may select internal defaults already scaled for bit depth. Check the specific function's parameter documentation when migrating scripts.

AviSynth+ example:

```avs
LoadPlugin("/path/to/neo-smo.dll")
clip = BlankClip(width=640, height=480, pixel_type="YUV420P10", length=24)
return neo_smo_Median(clip, radius=[1], planes=[0])
```

All 19 functions use the `neo_smo_` prefix, with argument names and order matching their API pages. Booleans use `true` / `false`; numeric arrays accept `[1, 2, 3]` or a single number. See the [AviSynth+ interface](docs/api/en/README.md#avisynth-calls-and-builds).

## SIMD and precision

Highway automatically selects a SIMD target included in the build and supported by the current CPU. There is no public `opt`, SIMD OFF, or `KernelInfo` interface.

Paths using half-precision arithmetic compute directly in F16 when the compilation target and CPU support it; otherwise, they compute in F32 and write back F16. Wider intermediate types are still used where required by the algorithm. DCTFilter always uses F32 transforms; F16 is used only for input and output storage.

FMA, operation order, and intermediate precision can affect rounding and threshold decisions. Wider SIMD or F16 storage does not guarantee that every filter runs faster. See [Samples, thresholds, and precision](docs/knowledge/en/shared/sample-and-precision.md).

## Building and testing

Requires CMake 3.24 or later, Git, and a C++17 compiler. Tests also require Python 3. Prefer clang-cl on Windows; Clang or GCC can be used on Linux. For native half-precision paths, Clang 22 or later with a supported CPU is recommended.

CMake fetches pinned versions of DualSynth2, Highway, and neo-libdct, and downloads the corresponding host headers when no local SDK is supplied. The default build includes both host interfaces and core tests. Core tests do not require an installed video host.

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/release --config Release --parallel 4
ctest --test-dir build/release -C Release --output-on-failure
```

For clang-cl on Windows, configure a Ninja build from a Visual Studio developer command prompt:

```sh
cmake -S . -B build/clang-cl -G Ninja -DCMAKE_CXX_COMPILER=clang-cl -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/clang-cl --parallel 4
ctest --test-dir build/clang-cl --output-on-failure
```

| Option | Purpose |
|---|---|
| `NEO_SMO_BUILD_AVISYNTH=OFF` | Disable the AviSynth+ interface. Required for Windows MinGW builds. |
| `NEO_SMO_BUILD_VAPOURSYNTH=OFF` | Disable the VapourSynth interface. Disabling both hosts builds only the computation core and its tests. |
| `BUILD_TESTING=OFF` | Disable tests. |
| `NEO_SMO_AVS_SDK=/path/to/sdk` | Use a local AviSynth+ SDK. |
| `NEO_SMO_TEST_AVISYNTH=ON` | Enable AVS host tests; also set `NEO_SMO_AVISYNTH_RUNTIME` to the runtime path. |
| `NEO_SMO_TEST_CROSS_HOST=ON` | Enable AVS/VS output comparisons; requires both interfaces, NumPy, and a VapourSynth Python environment. |
| `NEO_SMO_VS_SDK=/path/to/sdk` | Use a local VapourSynth SDK. |
| `Python3_EXECUTABLE=/path/to/python` | Select the Python interpreter for tests. |
| `NEO_SMO_TEST_BLACKBOX=ON` | Enable black-box comparisons against zsmooth, disabled by default; requires a VapourSynth Python environment, NumPy, and the reference plugin. |
| `NEO_SMO_REFERENCE_PLUGIN=/path/to/zsmooth` | Select the reference plugin file for black-box tests. |
| `FETCHCONTENT_SOURCE_DIR_DUALSYNTH2=/path/to/dualsynth2` | Use local DualSynth2 source. |
| `FETCHCONTENT_SOURCE_DIR_NEO_LIBDCT=/path/to/neo-libdct` | Use local neo-libdct source. |

The plugin build target is `neo_smo`, and its output basename is `neo-smo`. Tests cover sampling boundaries, integer and floating-point arithmetic, half-precision paths, DCT, and frame error handling. Comparisons against upstream are enabled separately.

## Performance

The table summarizes **neo-smo frame generation time / zsmooth frame generation time** by filter for tested VapourSynth configurations. **Values below 1 mean neo-smo is faster.**

| Filter | Time ratio range |
|---|---:|
| Median | 0.40–1.21× |
| InterQuartileMean | 0.31–1.23× |
| SmartMedian | 0.70–1.16× |
| VerticalCleaner | 0.50–1.00× |
| RemoveGrain | 0.42–1.14× |
| Repair | 0.48–1.02× |
| Clense | 0.50–1.00× |
| ForwardClense | 0.49–0.98× |
| BackwardClense | 0.49–0.98× |
| TemporalMedian | 0.49–1.02× |
| TemporalSoften | 0.66–1.19× |
| TemporalRepair | 0.27–1.20× |
| DegrainMedian | 0.29–0.68× |
| FluxSmoothT | 0.53–1.30× |
| FluxSmoothST | 0.76–1.18× |
| TTempSmooth | 0.55–0.58× |
| CCD | 0.69–1.37× |
| Cnr4 | 0.68–1.08× |
| DCTFilter | 0.23–0.28× |

F16 results use the corresponding upstream F32 path as the comparison baseline. Actual speed varies with input, parameters, hardware, and thread count.

## Development and contributions

Maintainers are responsible for technical direction, change review, and releases. Bug reports, suggestions, and contributions are welcome. Please discuss goals and approaches before changing numerical semantics, public interfaces, or major architecture.

This project uses AI assistance for implementation, testing, and review. Contributions should explain the problem, approach, validation, and how AI was involved. Reports should include the version, OS, CPU, compiler, build options, input/output formats, and a minimal reproducer. Numerical-difference reports should also identify the reference version, parameters, and request order; performance reports should describe the timing scope and thread configuration.

## Acknowledgments and license

The original zsmooth code's MIT copyright notice and full license text are preserved in the [zsmooth license file](LICENSES/zsmooth-MIT.txt).

neo-smo also uses the following libraries:

- [Google Highway](https://github.com/google/highway): cross-platform SIMD support.
- [DualSynth2](https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2): connects VapourSynth and AviSynth to the shared computation core.
- [neo-libdct](https://github.com/HomeOfAviSynthPlusEvolution/neo-libdct): block DCT, licensed under GNU GPL version 2 or later.

Thanks to the developers and users who contribute tests, reports, and improvements.

Thanks to [SB.SB](https://sb.sb) for sponsoring the LLM subscription used in this project's development.

neo-smo is licensed under the GNU General Public License, version 2 or later (`GPL-2.0-or-later`). See [LICENSE](LICENSE). Third-party components retain their own copyright notices and license terms.
