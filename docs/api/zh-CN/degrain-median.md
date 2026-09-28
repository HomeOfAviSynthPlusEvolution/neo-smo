# DegrainMedian

[返回 API 目录](README.md)

```python
core.neo_smo.DegrainMedian(clip, limit=None, mode=[1], interlaced=False,
                         norow=False, scalep=False)
```

支持通用格式。使用前、当前、后三帧的方向邻点选择修正值，最后限制相对原值的改变量。首尾各一帧复制。

| 参数 | 默认与范围 | 含义 |
|---|---|---|
| limit | 省略为各平面原生数值 4 | 每像素最大修正量；逐平面数组，最后一项重复 |
| mode | `[1]`，每项 0–5 | 方向选择的代价；0 仍然处理 |
| interlaced | False | 垂直邻点间距由 1 行改为 2 行，不改变时间帧间距 |
| norow | False | 排除当前帧的左右水平邻点对 |
| scalep | False | 是否缩放显式传入的 limit |

没有 `planes` 参数，使用 `limit=0` 停用对应平面；不能把三个 limit 都设为 0。显式原生 limit 的范围为整数 0–格式最大值、浮点亮度/RGB 0–1、浮点色度 −0.5–0.5。非正色度 limit 不处理。

`scalep=True` 时显式 limit 范围为 0–255，按 `limit × 格式最大值 / 255` 缩放；浮点色度的“最大值”为 0.5，其他浮点平面为 1。这与多数滤镜的整数左移缩放不同。

**省略 limit 保留原生 4，包括浮点输入，也不因 scalep=True 自动变为 4/255。** 为便于跨位深复用参数，建议明确填写 limit 和 scalep。

```python
out = core.neo_smo.DegrainMedian(src, limit=[4, 0, 0], mode=[1], scalep=True)
```

模式代价和方向解释见[时域去噪与修复](../../knowledge/zh-CN/temporal.md)。

## AviSynth

```text
neo_smo_DegrainMedian(clip, limit=Undefined(), mode=[1], interlaced=false, norow=false, scalep=false)
```

参数名称、顺序与上方共用说明一致；宿主格式、数组写法及音频/场序保留见 [API 目录](README.md#avisynth-调用与构建)。

`limit`、`mode` 接受数值数组或单值，`interlaced`、`norow`、`scalep` 接受布尔值。没有 `planes` 参数；下例以非零亮度 limit、零色度 limit 仅处理亮度。省略 limit 仍使用原生数值 4。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_DegrainMedian(src, limit=[4, 0, 0], mode=[1], interlaced=false, norow=false, scalep=true)
```
