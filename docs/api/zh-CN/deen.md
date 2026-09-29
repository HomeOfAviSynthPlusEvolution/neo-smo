# Deen、MiniDeen

[返回 API 索引](README.md)

```python
core.neo_smo.Deen(clip, mode="c3d", radius=None, threshold=None,
                  temporal_threshold=None, minimum=None, scenechange=0,
                  scalep=False, planes=None)
core.neo_smo.MiniDeen(clip, radius=None, threshold=None,
                      scalep=False, planes=None)
```

两者都接受固定格式、固定尺寸、无 alpha 的 Gray/YUV/RGB 平面输入。VapourSynth 支持 8–16 位整数、F16、F32；AviSynth+ 支持 8/10/12/14/16 位整数及 F32。参与处理的浮点样本必须有限，但不要求局限在 0–1。

## 公共参数

| 参数 | 默认值 | 含义 |
|---|---|---|
| `radius` | `[1]` | 逐平面空间半径。MiniDeen 和 Deen 2d 为 0–7；Deen 3d 为 0–4。0 完整复制该平面，包括关闭其时间处理。 |
| `threshold` | 见下表 | 当前帧的逐平面浮点阈值。 |
| `scalep` | `False` | 显式阈值是否使用 8 位标度。 |
| `planes` | 所有平面 | 不重复的有效平面索引。未选平面原样复制；空数组、重复或越界索引报错。 |

`radius`、`threshold`、`temporal_threshold`、`minimum` 接受单值或 1 到实际平面数个元素的数组。短数组重复最后一项；Gray 只接受一项。索引对应输入平面，与 `planes` 的顺序无关。省略参数使用默认值；显式空数组报错。未选平面对应的值也必须通过校验。数值参数拒绝 NaN 和 Inf。

整数参数应传整数。VS 的 Python 绑定可能在调用插件前把浮点数转成整数，插件无法恢复转换前的类型。

### 阈值单位

- `scalep=False`：显式值使用原生样本单位。整数范围为 0 到 `(2**bits-1)`；浮点范围为 0 到 1。
- `scalep=True`：显式值范围为 0–255。整数乘 `2**(bits-8)`，浮点除以 255。
- 省略阈值时，默认值始终按下表的 8 位标度换算，不受 `scalep` 影响。
- 有符号浮点色度的差值阈值同样使用非负的 0–1 范围。

| 8 位标度默认值 | Gray | RGB | YUV |
|---|---|---|---|
| Deen `threshold` | `[7]` | `[7,7,7]` | `[7,9,9]` |
| Deen `temporal_threshold` | `[4]` | `[4,4,4]` | `[4,6,6]` |
| MiniDeen `threshold` | `[10]` | `[10,10,10]` | `[10,10,10]` |

显式 `threshold=7` 会给所有平面重复 7，不再隐式生成色度特例。16 位输入下，`threshold=10, scalep=True` 表示 2560。关闭平面使用 `planes` 或 `radius=0`；零阈值仍按算法处理。

Deen 接受等号，即 `abs(sample-center) <= threshold`；MiniDeen 使用严格小于 `<`。整数输入的小数阈值也遵守此规则：差值 10 在 MiniDeen 阈值 10.5 下参与平均，阈值 10 下不参与。

## Deen

空间或固定三帧平滑，`mode` 控制整段剪辑：

| 模式族 | 行为 |
|---|---|
| `c2d`、`c3d` | 用中心替代超出阈值的样本，再按固定窗口大小平均。 |
| `a2d`、`a3d` | 阈值乘以距离因子，只对接受的样本等权平均。 |
| `w2d`、`w3d` | 用中心替代超出阈值的样本，再按距离权重平均。3d 的当前、前一、后一帧总权重比为 2:1:1。 |

| 参数 | 默认值 | 含义 |
|---|---|---|
| `mode` | `"c3d"` | 上述六种模式之一。 |
| `temporal_threshold` | 见上表 | 前后两帧共用的逐平面阈值；单位遵循 `scalep`。 |
| `minimum` | `[0.5]` | 窗口角点的逐平面距离因子，范围 0–1；a 模式调整阈值，w 模式调整权重。 |
| `scenechange` | `0` | 场景检测，见下文。 |

偏移 `(dx,dy)` 的距离因子为 `1-(1-minimum)*sqrt((dx²+dy²)/(2*radius²))`。空间越界坐标复制边缘像素。所有候选都与当前帧中心比较。c 模式不使用 `minimum`；2d 不使用 `temporal_threshold` 和 `scenechange`，但显式参数仍须合法。

3d 固定使用前、中、后三帧。首尾帧使用对应的 2d 操作；任一侧边界发生切换时，两侧时间邻居都不参与。

### 场景检测

沿用 TemporalSoften 的约定，不受 `scalep` 影响：

- `0`：关闭，默认值。
- `-1`：读取当前帧的 `_SceneChangePrev`、`_SceneChangeNext`，缺失属性视为没有切换。
- `1..254`：以 `scenechange/255` 为阈值自动检测。VS 优先使用 `misc.SCDetect`，插件或函数缺失时自动使用共享的内置亮度检测器；AVS 使用该内置检测器。两条路径都会写入 `_SceneChangePrev` 和 `_SceneChangeNext`。自动检测只支持至少两帧的 Gray/YUV；RGB 使用 0，或预先设置属性后使用 -1。

整数 c/a 平均采用精确的半值向上舍入。整数 w 与舍入后的双精度加权公式最多相差一个码值，常量保持不变；`w2d` 在 `minimum=1` 时采用精确 c2d 平均。浮点模式以 F32 计算差值与累加，阈值比较使用 F32，a 模式先乘距离因子。算术舍入可能随 SIMD target 有细微差异。

```python
out = core.neo_smo.Deen(src, mode="w2d", radius=[2,1], threshold=[7,9],
                       minimum=0.5, scalep=True, planes=[0])
out = core.neo_smo.Deen(src, mode="a3d", radius=2, threshold=7,
                       temporal_threshold=4, scenechange=12, scalep=True)
```

## MiniDeen

仅做空间处理。正方形窗口裁剪到有效坐标，窗口包含中心，并在总和与计数中额外加入两份中心。只接受差值严格小于阈值的样本。整数输出采用精确的半值向上舍入；F16/F32 使用 F32 中心差值累加，F16 在存储时转换。软件 target 使用普通 C++ 实现。

```python
out = core.neo_smo.MiniDeen(src, radius=[2,1], threshold=[10,6], scalep=True)
```

## AviSynth

使用 `neo_smo_Deen` 和 `neo_smo_MiniDeen`；参数名、顺序、默认值和语义同上。省略参数使用 `Undefined()`，布尔值使用 `false`/`true`。支持单值简写与原生数值数组。参见 [AviSynth 接口说明](README.md#avisynth-调用与构建)。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_Deen(src, mode="w3d", radius=2, threshold=7,
                    temporal_threshold=4, scalep=true, planes=[0])
```

## 从已撤回的接口迁移

旧参数 `rad`、`thrY`/`thrUV`、`tthY`/`tthUV`、`min`、`scd` 已移除。改用 `radius`、逐平面 `threshold`/`temporal_threshold`、`minimum` 和上述场景检测约定。缩放由 `(2**bits-1)/255` 改为项目通用的二次幂规则。MiniDeen 阈值支持小数，空 radius/threshold 数组报错。场景检测默认关闭。
