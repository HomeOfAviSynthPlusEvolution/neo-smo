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
| 1–254 | 调用 `misc.SCDetect`，阈值为此值 / 255；不支持 RGB 自动检测 |

因此传 `scenechange=True` 等价于数值 1，会开启自动检测，而不是表示“只读取已有属性”。

```python
out = core.neo_smo.TemporalMedian(src, radius=2)
out = core.neo_smo.TemporalSoften(src, radius=3, threshold=[4, 6, 6], scalep=True)
# 需要安装提供 SCDetect 的 misc 插件：
sc = core.misc.SCDetect(src, threshold=0.1)
out = core.neo_smo.TemporalMedian(sc, radius=2, scenechange=True)
```

参见[时域原理](../../knowledge/zh-CN/temporal.md)和[边界行为表](../../knowledge/zh-CN/shared/temporal-boundaries.md)。
