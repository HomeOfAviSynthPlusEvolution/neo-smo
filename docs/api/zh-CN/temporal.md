# TemporalMedian、TemporalSoften

[返回 API 目录](README.md)

```python
core.neo_smo.TemporalMedian(clip, radius=1, planes=None, scenechange=False)
core.neo_smo.TemporalSoften(clip, radius=4, threshold=None, scenechange=0,
                          scalep=False, planes=None)
```

支持通用格式，`radius` 是单个整数，范围 1–10，完整窗口为 `2r+1` 帧。二者均在相同像素坐标上取样，没有运动补偿。

## TemporalMedian

求时间窗口的中值。首尾各 `radius` 帧直接复制；这与把帧号限制到片头片尾不同。`scenechange=True` 按已有 `_SceneChangePrev`、`_SceneChangeNext` 缩短窗口，不自动检测场景，缺少所需属性时会报错。缩短后若样本数为偶数，取中间两值的平均，整数向下取整。

## TemporalSoften

`threshold` 为逐平面浮点数组，省略时使用自动缩放的 8 位标度 4。显式值在 `scalep=True` 时按 8 位标度缩放；否则使用原生单位。正阈值才会处理该平面，0 可用于停用 YUV 的某个平面；RGB/Gray 的第一个阈值不能为 0，三个阈值也不能全为 0。

原生整数阈值范围为 0–格式最大值；浮点亮度/RGB 为 0–1，浮点色度为 −0.5–0.5，非正色度阈值不处理。缩放后仍须通过格式范围检查。

与中心差的绝对值不大于阈值时使用邻帧值，否则用中心值替代，然后对整个有效窗口平均。**不是只除以通过阈值的样本数**。片头片尾缩短为实际存在的帧，不整帧跳过。

| scenechange | 行为 |
|---|---|
| 0（默认） | 不检查场景属性 |
| −1 | 使用已有属性；属性缺失时没有可用的切点 |
| 1–254 | 自动检测，阈值为此值 / 255；VS 调用 `misc.SCDetect`，AVS 使用内置检测；不支持 RGB 自动检测 |

在 VapourSynth 的 Python 调用中，传 `scenechange=True` 等价于数值 1，会开启自动检测，而不是表示“只读取已有属性”。

```python
out = core.neo_smo.TemporalMedian(src, radius=2)
out = core.neo_smo.TemporalSoften(src, radius=3, threshold=[4, 6, 6], scalep=True)
# 需要安装提供 SCDetect 的 misc 插件：
sc = core.misc.SCDetect(src, threshold=0.1)
out = core.neo_smo.TemporalMedian(sc, radius=2, scenechange=True)
```

参见[时域原理](../../knowledge/zh-CN/temporal.md)和[边界行为表](../../knowledge/zh-CN/shared/temporal-boundaries.md)。

## AviSynth

```text
neo_smo_TemporalMedian(clip, radius=1, planes=Undefined(), scenechange=false)
neo_smo_TemporalSoften(clip, radius=4, threshold=Undefined(), scenechange=0, scalep=false, planes=Undefined())
```

参数名称、顺序与上方共用说明一致；宿主格式、数组写法及音频/场序保留见 [API 目录](README.md#avisynth-调用与构建)。

TemporalMedian 的 `scenechange` 是布尔值：`true` 只读取已有整数 `_SceneChangePrev` / `_SceneChangeNext` 属性，不会自动检测。TemporalSoften 的同名参数是整数：0 关闭，−1 读已有属性，1–254 使用内置亮度帧差检测，阈值为参数 /255，无需 misc 插件。AVS 不接受用布尔值替代这个整数参数。自动检测只支持至少两帧的 Gray/YUV。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_TemporalMedian(src, radius=2, scenechange=false)
```

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_TemporalSoften(src, radius=3, threshold=[4, 6, 6], scenechange=12, scalep=true)
```
