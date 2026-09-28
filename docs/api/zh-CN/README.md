# neo-smo API

neo-smo 将 zsmooth 的空间、时域去噪与修复算法移植到 C++ / Highway。目前公开接口为 VapourSynth，命名空间为 `core.neo_smo`，插件 ID 为 `org.neofilters.neo_smo`。本文按当前源码描述行为；已提供 DCTFilter，尚未提供 AviSynth 接口。

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

以下各页示例假定已有 `core` 与 `src`。函数返回 `VideoNode`；页面中的签名是便于阅读的 Python 调用形式，`None` 表示省略参数，不是向插件传入一个数值。

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
| 色彩去噪 | [CCD](ccd.md) | 由多通道差异引导空间取样 |
| 色度时域去噪 | [Cnr4](cnr4.md) | 亮度与色度差异联合控制混合 |

## 通用约定

- 输入必须是固定格式、固定尺寸的视频。通用滤镜接受 Gray、RGB、YUV 的 8–16 位整数、16 位浮点（F16）与 32 位浮点（F32）。例外：TTempSmooth 不接受 F16；CCD 不接受 Gray；Cnr4 只接受 YUV 8–16 位整数。
- `planes` 省略时处理全部平面；显式传入时须为非空数组，宿主不接受 `[]`。YUV 的 0/1/2 是 Y/U/V；RGB 是 R/G/B；Gray 只有 0。重复或越界编号会报错。
- `radius`、`mode` 等逐平面数组一般将最后一项重复到剩余平面，如 `[1, 0]` 对三平面等价于 `[1, 0, 0]`。Cnr4 的三个数组必须各有三项；DCTFilter 的 `factors` 必须恰好八项，按频率而非按平面解释，见各专页。必填的 `mode` 不能是空数组。
- 空间半径按各平面自己的像素计算。YUV420 的色度平面半径 1 并不等于亮度平面半径 1 的画面覆盖范围。
- 参考片必须与主片具有相同尺寸、颜色家族、采样类型、位深及色度抽样。接口不检查帧率和帧数是否一致；调用方应保证时间对齐和有效帧范围。
- 输出保持输入的格式与尺寸，未处理平面沿用主片内容。没有自动位深转换、颜色矩阵转换或运动补偿。

## 阈值单位

含 `scalep` 的接口默认 `scalep=False`：**显式传入的数值按当前格式解释**。设为 `True` 时通常按 8 位标度缩放：整数乘 `2^(位深−8)`，浮点除以 255。DegrainMedian 使用另一套比例，见专页。

省略阈值可能触发内部缩放后的默认值，因此不一定等价于显式传入同一个数字。例如 10 位 TemporalSoften 省略 `threshold` 得到 16，而显式 `[4]` 且不设 `scalep` 得到 4。

详细说明见[采样、阈值与精度](../../knowledge/zh-CN/shared/sample-and-precision.md)。场景切换参数没有统一约定，见[时间窗口与场景切换](../../knowledge/zh-CN/shared/temporal-boundaries.md)。

## 迁移与原理

- [从 zsmooth 迁移](migration.md)
- [算法知识目录](../../knowledge/zh-CN/README.md)
