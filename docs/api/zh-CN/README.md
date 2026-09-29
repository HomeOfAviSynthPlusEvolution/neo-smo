# neo-smo API

neo-smo 将 zsmooth 的空间、时域去噪与修复算法移植到 C++ / Highway。提供 VapourSynth 与 AviSynth+ 接口。VapourSynth 命名空间为 `core.neo_smo`，插件 ID 为 `org.neofilters.neo_smo`；AviSynth+ 的全部 21 个函数以 `neo_smo_` 为前缀。

## 快速开始

```python
import vapoursynth as vs
core = vs.core
# 未安装到自动加载目录时，填写实际插件路径：
# core.std.LoadPlugin(path=r"C:\plugins\neo-smo.dll")

src = core.std.BlankClip(width=640, height=480, format=vs.YUV420P10, length=24)
out = core.neo_smo.Median(src, radius=[1], planes=[0])
out.set_output()
```

以下各页的 VapourSynth 示例假定已有 `core` 与 `src`。VapourSynth 函数返回 `VideoNode`，AviSynth 函数返回剪辑；页面中的签名是便于阅读的 Python 调用形式，`None` 表示省略参数，不是向插件传入一个数值。

AviSynth+ 的加载方法、参数写法与宿主差异见 [AviSynth+ 接口](#avisynth-调用与构建)。各滤镜页面使用 Python 形式展示两端共用的参数。

## 接口目录

| 类别 | 接口 | 作用 |
|---|---|---|
| 空间排序 | [Median、InterQuartileMean、SmartMedian](spatial-median.md) | 中值、四分位均值、条件中值修正 |
| 垂直处理 | [VerticalCleaner](vertical-cleaner.md) | 垂直中值或外推限幅 |
| 空间去噪 | [RemoveGrain](remove-grain.md) | 24 种邻域处理模式 |
| 块频率处理 | [DCTFilter](dct-filter.md) | 固定 8×8 DCT 系数加权 |
| 空间修复 | [Repair](repair.md) | 根据参考片的空间邻域限制改动 |
| 三帧处理 | [Clense、ForwardClense、BackwardClense](clense.md) | 时域限幅和单向外推 |
| 时域排序与平均 | [TemporalMedian、TemporalSoften](temporal.md) | 时间窗口中值或阈值平均 |
| 时域修复 | [TemporalRepair](temporal-repair.md) | 参考片的时域／时空约束 |
| 时空方向处理 | [DegrainMedian](degrain-median.md) | 选择方向并限制单次改动 |
| 时域极值处理 | [FluxSmoothT、FluxSmoothST](flux-smooth.md) | 对时域极值做阈值平均 |
| 加权时域平滑 | [TTempSmooth](ttempsmooth.md) | 根据连续帧差异与距离加权 |
| 阈值空间与时空平滑 | [Deen、MiniDeen](deen.md) | 固定分母、自适应及距离加权邻域平均 |
| 色彩去噪 | [CCD](ccd.md) | 由多通道差异引导空间取样 |
| 色度时域去噪 | [Cnr4](cnr4.md) | 亮度与色度差异联合控制混合 |

## 通用约定

- 输入必须是固定格式、固定尺寸的视频。通用滤镜接受 Gray、RGB、YUV 的 8–16 位整数、16 位浮点（F16）与 32 位浮点（F32）。AviSynth+ 不提供 F16，只接受无 alpha 的平面格式。例外：TTempSmooth 不接受 F16；CCD 不接受 Gray；Cnr4 只接受 YUV 8–16 位整数。
- `planes` 省略时处理全部平面；显式传入须为非空数组，`planes=[]` 会报错。YUV 的 0/1/2 是 Y/U/V；RGB 是 R/G/B；Gray 只有 0。重复或越界编号会报错。
- `radius`、`mode` 等逐平面数组一般将最后一项重复到剩余平面，如 `[1, 0]` 对三平面等价于 `[1, 0, 0]`。Cnr4 的三个数组必须各有三项；DCTFilter 的 `factors` 必须恰好八项，按频率而非按平面解释，见各专页。必填的 `mode` 不能是空数组。
- 空间半径按各平面自己的像素计算。YUV420 的色度平面半径 1 并不等于亮度平面半径 1 的画面覆盖范围。
- 参考片必须与主片具有相同尺寸、颜色家族、采样类型、位深及色度抽样。接口不检查帧率和帧数是否一致；调用方应保证时间对齐和有效帧范围；AVS 在请求超出参考片帧数时报告错误。
- 输出保持输入的格式与尺寸，未处理平面沿用主片内容。没有自动位深转换、颜色矩阵转换或运动补偿。

## 阈值单位

含 `scalep` 的接口默认 `scalep=False`：**显式传入的数值按当前格式解释**。设为 `True` 时通常按 8 位标度缩放：整数乘 `2^(位深−8)`，浮点除以 255。DegrainMedian 使用另一套比例，见专页。

省略阈值可能触发内部缩放后的默认值，因此不一定等价于显式传入同一个数字。例如 10 位 TemporalSoften 省略 `threshold` 得到 16，而显式 `[4]` 且不设 `scalep` 得到 4。

详细说明见[采样、阈值与精度](../../knowledge/zh-CN/shared/sample-and-precision.md)。场景切换参数没有统一约定，见[时间窗口与场景切换](../../knowledge/zh-CN/shared/temporal-boundaries.md)。

## 迁移与原理

- [从 zsmooth 迁移](migration.md)
- [算法知识目录](../../knowledge/zh-CN/README.md)

## AviSynth 调用与构建

签名中的 `Undefined()` 表示省略参数；调用时通常直接不写该参数。它不会传入零或空数组。各滤镜页面分别提供 AVS 签名和可运行示例。

需要 AviSynth+ 接口版本 11 或更新版本。加载同一个 neo-smo 插件后，全部 21 个函数可通过 `neo_smo_` 前缀调用，例如 `neo_smo_Repair`、`neo_smo_TTempSmooth` 和 `neo_smo_DCTFilter`。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_Median(src, radius=[1], planes=[0])
```

### 参数与格式

- 函数名去掉 `neo_smo_` 后，与上方目录中的名称一致。参数名称、顺序、默认值及算法行为共用各滤镜页面的说明。
- Python 示例中的 `True` / `False` 在 AVS 中写作 `true` / `false`。`scalep`、`interlaced`、`norow`、`fp` 和布尔型 `scenechange` 接受布尔值；TemporalSoften 和 Deen 的 `scenechange` 为整数。
- 数值数组写作 `[1, 2, 3]`，也接受单个数值作为单元素数组。不解析字符串数组；空数组通常也不接受。CCD 的 `points` 使用三个整数，0 表示关闭，非零表示开启。
- `planes` 使用 Y/U/V 或 R/G/B 的 0/1/2 索引，Gray 只有 0；不受 AVS 内部平面存储顺序影响。只对原本提供 `planes` 的函数开放该参数。
- 支持 Gray、YUV、RGB 的平面 8/10/12/14/16 位整数与 F32。AVS 不提供 F16 接口，不接受 packed 或带 alpha 的格式。CCD 不支持 Gray，Cnr4 仅支持整数 YUV。
- 参考片必须具有相同格式和尺寸。帧率和帧数不强制相同，但调用方必须保证时间对齐；请求超过参考片有效范围时会报错。
- 输出保留主片音频、场序和当前帧属性，未处理平面复制主片内容。各帧请求使用独立工作空间，可配合 `Prefetch`。

### 场景切换与缩放

TemporalSoften/Deen 的正 `scenechange` 和 TTempSmooth 的正 `scthresh` 使用内置亮度帧差检测，无需安装 VS 的 misc 插件。阈值单位分别为 `/255` 和 `/100`；有 `pfclip` 时 TTempSmooth 在参考片上检测。自动检测要求至少两帧，不接受 RGB。

读取已有属性时，使用整数 `_SceneChangePrev` 与 `_SceneChangeNext`。TemporalMedian 启用检测、Cnr4 保持默认检测时，缺失属性会报错。可以先准备场景属性，或明确关闭 `scenechange`；TTempSmooth 则用 `scthresh=0` 关闭。

CCD 与 Cnr4 在 YUV 色度抽样输入上，使用 AVS 的 `ExtractY` / `BilinearResize` 将参考亮度缩放到色度尺寸。AVS 与 VS 的缩放取整可能不同，进而影响接受邻居的判断和最终像素；因此这类输入不保证跨宿主逐位一致。

### 构建与测试

`NEO_SMO_BUILD_AVISYNTH` 默认开启，可与 `NEO_SMO_BUILD_VAPOURSYNTH` 独立选择。Windows AVS 构建使用 clang-cl 或 MSVC；MinGW 必须关闭 AVS 接口。`NEO_SMO_AVS_SDK` 可指定本地头文件目录。

`NEO_SMO_TEST_AVISYNTH=ON` 配合 `NEO_SMO_AVISYNTH_RUNTIME` 启用独立宿主验收。`NEO_SMO_TEST_CROSS_HOST=ON` 额外需要两个宿主接口、VapourSynth Python 与 NumPy，进行同输入输出对照。测试子进程具有超时与 Windows 崩溃弹窗抑制。

设置 `NEO_SMO_TEST_VAPOURSYNTH=ON` 可启用覆盖全部 21 个函数的真实宿主冒烟测试。所选 Python 环境需安装 VapourSynth 和 NumPy，无需 zsmooth 参考插件。
