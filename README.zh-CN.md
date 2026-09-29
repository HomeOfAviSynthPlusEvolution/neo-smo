# neo-smo

[English](README.md) | **简体中文** | [日本語](README.ja.md)

neo-smo 是直接从 [zsmooth](https://github.com/adworacz/zsmooth) 移植的 VapourSynth 和 AviSynth 视频滤镜插件，提供空间与时间降噪、邻域修复、色彩降噪及块 DCT 滤波。

以 zsmooth 的 Zig 源码为基础，内部使用 C++17 和 Google Highway 实现并优化，块 DCT 使用 neo-libdct。DualSynth2 将同一计算核心连接到两个宿主。VapourSynth 使用 `core.neo_smo`，AviSynth 使用带有 `neo_smo_` 前缀的函数。

## 项目来源

zsmooth 由 Austin Dworaczyk Wiltshire（adworacz）开发。neo-smo 的算法与接口直接承接 zsmooth；zsmooth 又沿用了视频滤镜社区长期积累的算法与实现经验。其文档列出的来源和参考包括：

- [AviSynth RgTools](https://github.com/pinterf/RgTools) 与 [VapourSynth RemoveGrain](https://github.com/vapoursynth/vs-removegrain)：RemoveGrain 系列的历史实现与行为参考。
- [VapourSynth TemporalSoften](https://github.com/dubhater/vapoursynth-temporalsoften2)、[VapourSynth TemporalMedian](https://github.com/dubhater/vapoursynth-temporalmedian) 与 [Neo Temporal Median](https://github.com/HomeOfAviSynthPlusEvolution/neo_TMedian)：时间平滑与时间中值滤镜。
- [VapourSynth FluxSmooth](https://github.com/dubhater/vapoursynth-fluxsmooth)：FluxSmooth 系列。
- [Dogway 的 AviSynth Scripts](https://github.com/Dogway/Avisynth-Scripts)：`ex_median`、四分位均值与 SmartMedian 等算法思路。
- [End-of-Eternity 的 CCD](https://github.com/End-of-Eternity/vs-ccd) 与 [vs-jetpack 的 CCD](https://github.com/Jaded-Encoding-Thaumaturgy/vs-jetpack/blob/e0f47d86930150fd0bf92b0845ccc2b0491f7807/vsdenoise/ccd.py#L95)：zsmooth CCD 的主要实现参考；CCD 最初由 Sergey Stolyarevsky 为 VirtualDub 编写。
- [Cnr2](http://avisynth.nl/index.php/Cnr2)：zsmooth Cnr4 的算法来源之一。

感谢 zsmooth 及这些上游项目的作者与贡献者。

## 设计

neo-smo 将滤镜计算与宿主帧管理分离。核心负责图像平面的采样、排序、加权和变换，宿主层负责参数、帧请求、属性及输出分配。核心可以独立构建和测试。DCTFilter 使用 neo-libdct 的固定 8×8 变换，不依赖 FFTW 运行时。

移植保留 zsmooth 的公开函数名称和参数结构，同时修正已发现的问题并优化计算路径。修正上游问题、FMA 和浮点舍入都可能影响结果，不保证与所有历史 zsmooth 版本逐像素一致。具体调用与数值约定见[迁移说明](docs/api/zh-CN/migration.md)。

滤镜和 DCT 内部均在宿主调用线程上计算，不创建工作线程或线程池。宿主仍可并发请求多个帧；SIMD 和批量 DCT 不表示内部多线程。

## 支持的操作

| 类别 | 函数与用途 |
|---|---|
| 空间排序与清理 | `Median`、`InterQuartileMean`、`SmartMedian`、`VerticalCleaner`：中值、四分位均值、条件中值修正及垂直清理。 |
| 邻域去噪与修复 | `RemoveGrain`、`Repair`：多模式空间处理及参考邻域限幅。 |
| 时间限幅与修复 | `Clense`、`ForwardClense`、`BackwardClense`、`TemporalRepair`：利用相邻帧约束当前像素。 |
| 时间与时空降噪 | `TemporalMedian`、`TemporalSoften`、`DegrainMedian`、`FluxSmoothT`、`FluxSmoothST`、`TTempSmooth`：时间排序、阈值平均、方向选择及加权平滑。 |
| 阈值邻域平滑 | [Deen、MiniDeen](docs/api/zh-CN/deen.md)：空间及三帧降噪，支持固定分母、自适应与距离加权平均。 |
| 色彩降噪 | `CCD`、`Cnr4`：颜色差异引导的空间采样与时间色度降噪。 |
| 块频率滤波 | `DCTFilter`：不重叠 8×8 块的 DCT 系数加权。 |

共 21 个函数。通常支持固定格式、固定尺寸的平面 GRAY/YUV/RGB，以及 8–16 位整数、16 位浮点（F16）、32 位浮点（F32）样本。F16 仅适用于 VapourSynth；AviSynth+ 支持无 alpha 的平面格式。TTempSmooth 不接受 F16，CCD 不接受 GRAY，Cnr4 仅接受整数 YUV。

输出保留输入格式、尺寸、帧数和帧率。支持 `planes` 的函数在省略该参数时处理全部平面，`planes=[0]` 只处理第一个平面；显式传入空平面数组会报错。部分函数用 `mode` 或阈值控制平面处理，不能向所有函数统一传入 `planes`。这些滤镜不估计运动向量，也不做运动补偿。

## 文档与使用

API 文档说明怎样调用函数，知识库说明输入怎样经过采样、排序、加权和变换成为输出。

- [API 使用参考](docs/api/zh-CN/README.md)：函数签名、参数、默认值及使用示例。
- [计算原理知识库](docs/knowledge/zh-CN/README.md)：数据表示、公式、计算顺序、边界及精度。
- [从 zsmooth 迁移](docs/api/zh-CN/migration.md)：命名空间、阈值单位、场景切换与精度约定。

可以显式加载构建出的插件，也可以将其放入 VapourSynth 的插件自动加载目录。下面使用 Windows 文件名；Linux 使用 `neo-smo.so`，其他平台请替换为实际插件路径。VapourSynth 插件标识符为 `org.neofilters.neo_smo`。

```python
import vapoursynth as vs

core = vs.core
core.std.LoadPlugin(path="/path/to/neo-smo.dll")

clip = core.std.BlankClip(width=640, height=360, format=vs.YUV420P8, length=24)
output = core.neo_smo.Median(clip, radius=[1], planes=[0])
output.set_output()
```

这个最小示例使用合成剪辑展示亮度中值滤波。实际使用时替换 `clip` 即可。带 `scalep` 的函数默认使用原生单位解释显式阈值；省略参数所选的内部默认值可能已按位深缩放。迁移脚本时应核对具体函数的参数说明。

同一个插件文件也提供 AviSynth C++ 接口，需要支持接口版本 11 的宿主，通过 `LoadPlugin` 加载：

```avs
LoadPlugin("/path/to/neo-smo.dll")
clip = BlankClip(width=640, height=360, length=24, pixel_type="YV12")
return neo_smo_Median(clip, radius=[1], planes=[0])
```

函数名和参数顺序对应 API 使用参考，函数名前加 `neo_smo_`。数组参数接受 `[0, 1]` 这样的原生数组，单值可作为一个元素的简写；DCTFilter 的 `factors` 仍需恰好八项；Cnr4 的 `sense`、`str`、`pow` 各须恰好三项，不能用单值替代。布尔参数使用 `true`/`false`。音频和场序从主输入剪辑传递，输出帧属性来自对应的源帧。详见 [AviSynth 接口](docs/api/zh-CN/README.md#avisynth-调用与构建)。

显式表达参数省略语义时，VapourSynth 使用 `None`，AviSynth 使用 `Undefined()`；也可以直接不写该参数。省略参数会使用相应默认行为，不等于传入零或空数组。例如 `planes=None` 或 `planes=Undefined()` 使用默认平面选择，滤镜不接受 `planes=[]`。

## SIMD 与 CPU 选择

Highway 自动选择构建中包含且当前 CPU 支持的 SIMD 目标。没有公开的 `opt`、SIMD OFF 或 `KernelInfo` 接口。

使用半精度算术的路径在编译目标和 CPU 支持时直接进行 F16 运算，否则使用 F32 计算后写回 F16。算法需要时仍会使用更宽的中间类型。DCTFilter 始终使用 F32 变换，F16 只用于输入输出存储。

FMA、运算顺序与中间精度可能影响舍入和阈值判断。更宽的 SIMD 或 F16 存储不保证每个滤镜都更快，详见[采样、阈值与精度](docs/knowledge/zh-CN/shared/sample-and-precision.md)。

## 构建与测试

需要 CMake 3.24 或更新版本、Git 及支持 C++17 的编译器；测试还需要 Python 3。CMake 获取固定版本的 DualSynth2、Highway 和 neo-libdct。两个宿主的 SDK 均可从本地发现或自动下载。

下面构建双宿主插件和核心测试，不要求本地已安装视频宿主。Windows 优先使用 clang-cl；MinGW 构建须关闭 AVS 入口。Linux 可使用 Clang 或 GCC；原生半精度路径需编译器和 CPU 同时支持，建议使用 Clang 22 或更新版本。

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
| `NEO_SMO_BUILD_VAPOURSYNTH=OFF` | 禁用 VapourSynth 入口。 |
| `NEO_SMO_BUILD_AVISYNTH=OFF` | 禁用 AviSynth 入口；两个宿主选项均为 OFF 时只构建核心及其测试。 |
| `BUILD_TESTING=OFF` | 不构建测试。 |
| `NEO_SMO_TEST_VAPOURSYNTH=ON` | 启用 VapourSynth 宿主测试，默认关闭；所选 Python 环境需安装 VapourSynth、NumPy 和匹配架构的运行时。 |
| `Python3_EXECUTABLE=/path/to/python` | 指定宿主测试所用 Python。 |
| `NEO_SMO_VS_SDK=/path/to/sdk` | 指定本地 VapourSynth SDK。 |
| `NEO_SMO_AVS_SDK=/path/to/sdk` | 指定本地 AviSynth SDK。 |
| `NEO_SMO_TEST_AVISYNTH=ON` | 启用 AviSynth 宿主测试，默认关闭；需要 Python 和匹配架构的运行时。 |
| `NEO_SMO_AVISYNTH_RUNTIME=/path/to/avisynth.dll` | 指定 AviSynth 宿主测试使用的运行时库。 |
| `NEO_SMO_TEST_CROSS_HOST=ON` | 启用 AVS/VS 输出对照，默认关闭；需两个宿主、NumPy 和 VapourSynth Python 环境。 |
| `NEO_SMO_TEST_BLACKBOX=ON` | 启用与 zsmooth 的黑盒对照，默认关闭；需要 VapourSynth Python 环境、NumPy 与参考插件。 |
| `NEO_SMO_REFERENCE_PLUGIN=/path/to/zsmooth` | 指定黑盒测试使用的参考插件文件。 |
| `FETCHCONTENT_SOURCE_DIR_DUALSYNTH2=/path/to/dualsynth2` | 使用本地 DualSynth2 源码代替固定版本下载。 |
| `FETCHCONTENT_SOURCE_DIR_NEO_LIBDCT=/path/to/neo-libdct` | 使用本地 neo-libdct 源码代替固定版本下载。 |

插件构建目标为 `neo_smo`，输出文件基本名为 `neo-smo`，默认包含两个宿主入口。测试覆盖采样边界、整数与浮点运算、半精度路径、DCT 及宿主行为，并包含标量参考计算核对和当前机器支持且已编译的不同 Highway 目标之间的数值对照。对照按各路径的精度规则判断，不要求所有浮点路径逐位一致。与 zsmooth 比较的黑盒测试另需参考插件。

CI 包含 Windows x64、Linux x64/ARM64、macOS ARM64 及 Linux ASan/UBSan 检查。Linux runner 使用 Ubuntu 26.04，GCC 使用系统默认版本，Clang 使用版本 22。发布工作流构建 Windows、Linux 和 macOS 的 x64/ARM64 产物，VapourSynth 宿主测试目前在 Windows x64 上运行，AviSynth 宿主测试通过上述选项单独启用。工作流产物的实际构建环境和测试范围记录在包内，不代表适用于所有 Linux 发行版。

## 性能

下表按滤镜汇总 VapourSynth 中已测配置的 **neo-smo 帧生成耗时 / zsmooth 帧生成耗时**；**小于 1 表示 neo-smo 更快**。下列范围跨越已测的 8/16 位整数、F16/F32、颜色格式及参数配置，不是平均值或置信区间；各滤镜仅包含其支持的已测配置。

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

这些数据来自单个滤镜的帧生成计时，不代表整条处理链的吞吐量。F16 结果以对应的上游 F32 为比较基准。实际速度随输入、参数、硬件和线程数变化。

## 开发与贡献

维护者负责技术方向、变更审核和发布。欢迎问题报告、建议与贡献；修改数值语义、公开接口或重要架构前，建议先讨论目标和方案。

本项目使用 AI 辅助实现、测试和审查。贡献应说明问题、方案、验证方法和 AI 参与方式。报告问题请提供版本、系统、CPU、编译器、构建选项、输入输出格式及最小复现；数值差异还应注明参考版本、参数和请求顺序，性能报告应注明计时范围及线程配置。

## 致谢与许可证

neo-smo 还使用了以下基础库：

- [Google Highway](https://github.com/google/highway)：提供跨平台 SIMD 支持。
- [DualSynth2](https://github.com/HomeOfAviSynthPlusEvolution/dualsynth2)：连接 VapourSynth、AviSynth 与共享计算核心。
- [neo-libdct](https://github.com/HomeOfAviSynthPlusEvolution/neo-libdct)：用于 DCTFilter 的固定 8×8 块变换。

感谢参与测试、报告问题和改进的开发者与用户。

感谢 [烧饼论坛](https://sb.sb) 赞助本项目开发使用的 LLM 订阅。

zsmooth 原始代码的 MIT 版权声明与完整许可文本保留在 [zsmooth 许可文件](LICENSES/zsmooth-MIT.txt) 中。

neo-smo 采用 GNU 通用公共许可证第 2 版或更新版本（`GPL-2.0-or-later`），完整条款见 [LICENSE](LICENSE)。第三方组件保留各自的版权声明和许可条款。
