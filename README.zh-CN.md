# neo-smo

[English](README.md) | **简体中文** | [日本語](README.ja.md)

neo-smo 是直接从 [zsmooth](https://github.com/adworacz/zsmooth) 移植的 AviSynth+ / VapourSynth 视频滤镜插件，以其 Zig 源码为基础，改用 C++17 和 Google Highway 实现并优化。项目提供空间与时间降噪、邻域修复、色彩降噪及块 DCT 滤波，VapourSynth 使用 `core.neo_smo` 命名空间，AviSynth+ 使用 `neo_smo_` 函数前缀。

## 项目来源

zsmooth 由 Austin Dworaczyk Wiltshire（adworacz）开发。neo-smo 的算法与接口直接承接 zsmooth；zsmooth 又沿用了视频滤镜社区长期积累的算法与实现经验。其文档列出的来源和参考包括：

- [AviSynth RgTools](https://github.com/pinterf/RgTools) 与 [VapourSynth RemoveGrain](https://github.com/vapoursynth/vs-removegrain)：RemoveGrain 系列的历史实现与行为参考。
- [VapourSynth TemporalSoften](https://github.com/dubhater/vapoursynth-temporalsoften2)、[VapourSynth TemporalMedian](https://github.com/dubhater/vapoursynth-temporalmedian) 与 [Neo Temporal Median](https://github.com/HomeOfAviSynthPlusEvolution/neo_TMedian)：时间平滑与时间中值滤镜。
- [VapourSynth FluxSmooth](https://github.com/dubhater/vapoursynth-fluxsmooth)：FluxSmooth 系列。
- [Dogway 的 AviSynth Scripts](https://github.com/Dogway/Avisynth-Scripts)：`ex_median`、四分位均值与 SmartMedian 等算法思路。
- [End-of-Eternity 的 CCD](https://github.com/End-of-Eternity/vs-ccd) 与 [vs-jetpack 的 CCD](https://github.com/Jaded-Encoding-Thaumaturgy/vs-jetpack/blob/e0f47d86930150fd0bf92b0845ccc2b0491f7807/vsdenoise/ccd.py#L95)：zsmooth CCD 的主要实现参考；CCD 最初由 Sergey Stolyarevsky 为 VirtualDub 编写。
- [Cnr2](http://avisynth.nl/index.php/Cnr2)：zsmooth Cnr4 的算法来源之一。

感谢 zsmooth 及这些上游项目的作者与贡献者。

## 实现

计算核心负责图像平面的采样、排序、加权和变换，宿主层负责参数、帧请求、属性及输出分配。Google Highway 提供跨平台 SIMD 和运行时指令集选择，DualSynth2 提供宿主接入所需的基础设施。DCTFilter 使用 neo-libdct 的固定 8×8 变换，不依赖 FFTW 运行时。

移植保留 zsmooth 的公开函数名称和参数结构，同时修正已发现的问题并优化计算路径。生产构建允许 FMA 与浮点舍入差异，不保证不同 CPU 或历史 zsmooth 版本之间逐位一致。具体格式、阈值单位与边界行为见 API 文档。

同一个插件提供 VapourSynth 与 AviSynth+ 接口，也可在构建时单独选择宿主。AviSynth+ 需要接口版本 11 或更新版本。

## 支持的操作

| 类别 | 函数与用途 |
|---|---|
| 空间排序与清理 | `Median`、`InterQuartileMean`、`SmartMedian`、`VerticalCleaner`：中值、四分位均值、条件中值修正及垂直清理。 |
| 邻域去噪与修复 | `RemoveGrain`、`Repair`：多模式空间处理及参考邻域限幅。 |
| 时间限幅与修复 | `Clense`、`ForwardClense`、`BackwardClense`、`TemporalRepair`：利用相邻帧约束当前像素。 |
| 时间与时空降噪 | `TemporalMedian`、`TemporalSoften`、`DegrainMedian`、`FluxSmoothT`、`FluxSmoothST`、`TTempSmooth`：时间排序、阈值平均、方向选择及加权平滑。 |
| 色彩降噪 | `CCD`、`Cnr4`：颜色差异引导的空间采样与时间色度降噪。 |
| 块频率滤波 | `DCTFilter`：不重叠 8×8 块的 DCT 系数加权。 |

共 19 个函数。通常支持固定格式、固定尺寸的平面 Gray、YUV、RGB，以及 8–16 位整数、16 位浮点（F16）、32 位浮点（F32）样本。F16 仅适用于 VapourSynth；AviSynth+ 支持无 alpha 的平面格式。TTempSmooth 不接受 F16，CCD 不接受 Gray，Cnr4 仅接受整数 YUV。

输出保留输入格式和尺寸。支持 `planes` 的函数在省略该参数时处理全部平面，`planes=[0]` 只处理第一个平面；显式空数组不被接受。部分函数用 `mode` 或阈值控制平面处理，不能向所有函数统一传入 `planes`。这些滤镜不估计运动向量，也不做运动补偿。

## 文档与使用

- [API 使用参考](docs/api/zh-CN/README.md)：函数签名、参数、默认值与示例。
- [算法知识库](docs/knowledge/zh-CN/README.md)：计算过程、数值规则与边界行为。
- [从 zsmooth 迁移](docs/api/zh-CN/migration.md)：命名空间、阈值单位、场景切换与精度约定。

将插件放入 VapourSynth 的自动加载目录，或使用 `LoadPlugin` 显式加载。下面使用 Windows 文件名；其他平台请替换为实际插件路径。插件标识符为 `org.neofilters.neo_smo`。

```python
import vapoursynth as vs

core = vs.core
core.std.LoadPlugin(path="/path/to/neo-smo.dll")

clip = core.std.BlankClip(width=640, height=480, format=vs.YUV420P10, length=24)
output = core.neo_smo.Median(clip, radius=[1], planes=[0])
output.set_output()
```

这个最小示例使用合成剪辑展示亮度中值滤波。实际使用时替换 `clip` 即可。带 `scalep` 的函数默认使用原生单位解释显式阈值；省略参数所选的内部默认值可能已按位深缩放。迁移脚本时应核对具体函数的参数说明。

AviSynth+ 示例：

```avs
LoadPlugin("/path/to/neo-smo.dll")
clip = BlankClip(width=640, height=480, pixel_type="YUV420P10", length=24)
return neo_smo_Median(clip, radius=[1], planes=[0])
```

19 个函数均可使用 `neo_smo_` 前缀调用，参数顺序和名称与 API 各页一致；布尔值使用 `true` / `false`，数值数组支持 `[1, 2, 3]` 或单个数值。详见 [AviSynth+ 接口](docs/api/zh-CN/README.md#avisynth-调用与构建)。

## SIMD 与精度

Highway 自动选择构建中包含且当前 CPU 支持的 SIMD 目标。没有公开的 `opt`、SIMD OFF 或 `KernelInfo` 接口。

使用半精度算术的路径在编译目标和 CPU 支持时直接进行 F16 运算，否则使用 F32 计算后写回 F16。算法需要时仍会使用更宽的中间类型。DCTFilter 始终使用 F32 变换，F16 只用于输入输出存储。

FMA、运算顺序与中间精度可能影响舍入和阈值判断。更宽的 SIMD 或 F16 存储不保证每个滤镜都更快，详见[采样、阈值与精度](docs/knowledge/zh-CN/shared/sample-and-precision.md)。

## 构建与测试

需要 CMake 3.24 或更新版本、Git 和支持 C++17 的编译器；开启测试还需要 Python 3。Windows 优先使用 clang-cl，Linux 可使用 Clang 或 GCC。需要原生半精度路径时建议使用 Clang 22 或更新版本，并配合支持的 CPU。

CMake 获取固定版本的 DualSynth2、Highway 和 neo-libdct，并在没有本地 SDK 时下载相应宿主头文件。默认构建双宿主插件和核心测试，无需为核心测试安装视频宿主。

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/release --config Release --parallel 4
ctest --test-dir build/release -C Release --output-on-failure
```

Windows 使用 clang-cl 时，在 Visual Studio 开发者命令行中配置 Ninja 构建：

```sh
cmake -S . -B build/clang-cl -G Ninja -DCMAKE_CXX_COMPILER=clang-cl -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/clang-cl --parallel 4
ctest --test-dir build/clang-cl --output-on-failure
```

| 选项 | 用途 |
|---|---|
| `NEO_SMO_BUILD_AVISYNTH=OFF` | 不构建 AviSynth+ 接口。Windows MinGW 构建需关闭此项。 |
| `NEO_SMO_BUILD_VAPOURSYNTH=OFF` | 不构建 VapourSynth 接口；两个宿主均关闭时只构建计算核心及其测试。 |
| `BUILD_TESTING=OFF` | 不构建测试。 |
| `NEO_SMO_AVS_SDK=/path/to/sdk` | 指定本地 AviSynth+ SDK。 |
| `NEO_SMO_TEST_AVISYNTH=ON` | 启用 AVS 宿主测试，另需指定 `NEO_SMO_AVISYNTH_RUNTIME` 运行库路径。 |
| `NEO_SMO_TEST_CROSS_HOST=ON` | 启用 AVS/VS 输出对照，需两个宿主、NumPy 和 VapourSynth Python 环境。 |
| `NEO_SMO_VS_SDK=/path/to/sdk` | 指定本地 VapourSynth SDK。 |
| `Python3_EXECUTABLE=/path/to/python` | 指定测试使用的 Python。 |
| `NEO_SMO_TEST_BLACKBOX=ON` | 启用与 zsmooth 的黑盒对照，默认关闭；需要 VapourSynth Python 环境、NumPy 与参考插件。 |
| `NEO_SMO_REFERENCE_PLUGIN=/path/to/zsmooth` | 指定黑盒测试使用的参考插件文件。 |
| `FETCHCONTENT_SOURCE_DIR_DUALSYNTH2=/path/to/dualsynth2` | 使用本地 DualSynth2 源码。 |
| `FETCHCONTENT_SOURCE_DIR_NEO_LIBDCT=/path/to/neo-libdct` | 使用本地 neo-libdct 源码。 |

插件构建目标为 `neo_smo`，输出文件基本名为 `neo-smo`。测试覆盖采样边界、整数与浮点运算、半精度路径、DCT 和帧错误处理；与上游的对照测试单独启用。

## 性能

下表按滤镜汇总 VapourSynth 中已测配置的 **neo-smo 帧生成耗时 / zsmooth 帧生成耗时**；**小于 1 表示 neo-smo 更快**。

| 滤镜 | 耗时比范围 |
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

F16 结果以对应的上游 F32 为比较基准。实际速度随输入、参数、硬件和线程数变化。

## 开发与贡献

维护者负责技术方向、变更审核和发布。欢迎问题报告、建议与贡献；修改数值语义、公开接口或重要架构前，建议先讨论目标和方案。

本项目使用 AI 辅助实现、测试和审查。贡献应说明问题、方案、验证方法和 AI 参与方式。报告问题请提供版本、系统、CPU、编译器、构建选项、输入输出格式及最小复现；数值差异还应注明参考版本、参数和请求顺序，性能报告应注明计时范围及线程配置。

## 致谢与许可证

zsmooth 原始代码的 MIT 版权声明与完整许可文本保留在 [zsmooth 许可文件](LICENSES/zsmooth-MIT.txt) 中。

neo-smo 还使用了以下基础库：

- [Google Highway](https://github.com/google/highway)：提供跨平台 SIMD 支持。
- [DualSynth2](https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2)：连接 VapourSynth、AviSynth 与共享计算核心。
- [neo-libdct](https://github.com/HomeOfAviSynthPlusEvolution/neo-libdct)：块 DCT，采用 GNU GPL 第 2 版或更新版本。

感谢参与测试、报告问题和改进的开发者与用户。

感谢 [烧饼论坛](https://sb.sb) 赞助本项目开发使用的 LLM 订阅。

neo-smo 采用 GNU 通用公共许可证第 2 版或更新版本（`GPL-2.0-or-later`），完整条款见 [LICENSE](LICENSE)。第三方组件保留各自的版权声明和许可条款。
