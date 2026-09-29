# Deen、MiniDeen

[返回 API 目录](README.md)

## Deen

```python
core.neo_smo.Deen(clip, mode="c3d", rad=1, thrY=7.0, thrUV=9.0,
                  tthY=4.0, tthUV=6.0, min=0.5, scd=9.0,
                  scenechange=True, planes=None)
```

空间或三帧阈值平滑。输入须为固定格式、固定尺寸的平面 GRAY/YUV/RGB。VapourSynth 支持 8–16 位整数、F16 和 F32；AviSynth+ 支持 8/10/12/14/16 位整数和 F32。不支持 alpha 或打包格式。

### 参数

| 参数 | 默认值 | 含义 |
|---|---|---|
| `mode` | `"c3d"` | `c2d`、`c3d`、`a2d`、`a3d`、`w2d` 或 `w3d`，见下文。 |
| `rad` | `1` | 空间半径，按各平面自身像素计算。2d 模式为整数 1–7，3d 模式为整数 1–4；窗口大小为 `(2*rad+1)²`。 |
| `thrY`、`thrUV` | `7`、`9` | 当前帧阈值，使用 0–255 标度。`thrUV` 用于 YUV 色度；GRAY 和 RGB 的所有平面使用 `thrY`。 |
| `tthY`、`tthUV` | `4`、`6` | 前后帧阈值，单位和平面映射同上，仅用于 3d 模式。 |
| `min` | `0.5` | 窗口角点的距离因子，范围 0–1。a 模式用它调整阈值，w 模式用它调整权重；c 模式不使用。 |
| `scd` | `9` | 非负场景切换阈值，使用 0–255 帧差标度。 |
| `scenechange` | `True` | 在 3d 模式启用内置场景检测，无需场景切换帧属性。 |
| `planes` | 全部 | 不重复的有效平面编号；省略时处理全部平面，`[]` 会报错。未选中的平面直接复制。 |

四个滤波阈值须为 0–255 内的有限数。整数输入始终按 `(2**bits-1)/255` 缩放，浮点输入按 `1/255` 缩放，没有 `scalep` 参数。浮点样本须为有限值，但不限于 0–1。

### 模式

所有候选样本都与当前帧的中心像素比较，差值等于阈值时接受。越出平面的空间坐标钳位到图像边缘。

| 模式族 | 行为 |
|---|---|
| `c` | 接受阈值内的样本，拒绝的样本替换成中心值，再按固定窗口大小求平均。 |
| `a` | 阈值乘以距离因子；仅对接受的样本等权平均，分母为接受数量。 |
| `w` | 使用未乘距离因子的阈值；拒绝的样本替换成中心值，再按距离权重平均。3d 模式将当前帧权重加倍，各帧总权重比为 2:1:1。 |

空间偏移 `(dx, dy)` 的距离因子为 `1 - (1-min)*sqrt((dx²+dy²)/(2*rad²))`，中心为 1，角点为 `min`。a 模式中，零阈值仍接受与中心相等的样本；w 模式中，零权重不参与计算。

2d 模式仅使用当前帧。3d 模式使用当前、前一和后一帧，每帧取相同空间半径。首尾帧使用对应的空间操作。若任一侧检测到场景切换，则同时舍弃两个时间邻帧，仅对当前帧做空间滤波。

场景检测先计算每个平面的平均绝对差，换算到 0–255 标度，再将各平面均值等权平均。检测使用全部输入平面，包括 `planes` 未选中的平面。指标严格大于 `scd` 才判为切换。`scd=0` 不会关闭检测；关闭请用 `scenechange=False`。

整数 c/a 均值使用精确的四舍五入（半值向上）。整数 w 结果相对双精度加权公式的舍入值最多相差 1 个码值，常量保持不变，不产生新极值。所有浮点模式直接读取 F16/F32 样本，以 F32 计算。阈值判断使用 F32 差值和阈值，接受舍入后的相等；a 模式先将阈值乘以距离因子，再转换到 F32。结果可能随算术舍入和 SIMD 目标略有变化。整数输入的 `w2d min=1` 使用精确的 c2d 均值。

```python
# 仅对亮度做空间加权平滑。
out = core.neo_smo.Deen(src, mode="w2d", rad=2, thrY=7, min=0.5, planes=[0])
# 三帧自适应平滑，使用内置场景检测。
out = core.neo_smo.Deen(src, mode="a3d", rad=2, thrY=7, tthY=4)
```

## MiniDeen

```python
core.neo_smo.MiniDeen(clip, radius=None, threshold=None, planes=None)
```

仅支持整数平面 GRAY/YUV/RGB 的空间平滑：VapourSynth 为 8–16 位，AviSynth+ 为 8/10/12/14/16 位。不接受浮点输入，不做时间处理或场景检测。

- `radius`：最多三个 1–7 内的整数，默认 `[1]`。
- `threshold`：最多三个 0–255 内的整数，默认 `[10]`。各值换算为 `floor(threshold*(2**bits-1)/255)`。
- 数组不足三项时，最后一项重复到剩余平面。省略或传入空的 radius/threshold 数组均使用默认值。
- `planes`：不重复的有效编号，省略时表示全部，空数组会报错。缩放后的阈值为 0 或 1 时，该平面直接复制。

只有与中心差值**严格小于**缩放阈值的样本参与平均。方形窗口包含中心；累加和与计数还额外包含两份中心值。图像边缘只取有效坐标，窗口截短，不复制边缘样本。最终整数均值按半值向上舍入。

```python
out = core.neo_smo.MiniDeen(src, radius=[2, 1], threshold=[10, 6], planes=[0, 1, 2])
```

## AviSynth

```text
neo_smo_Deen(clip, mode="c3d", rad=1, thrY=7.0, thrUV=9.0, tthY=4.0, tthUV=6.0,
             min=0.5, scd=9.0, scenechange=true, planes=Undefined())
neo_smo_MiniDeen(clip, radius=Undefined(), threshold=Undefined(), planes=Undefined())
```

名称、顺序、默认值和行为同上。数组参数也接受单值，视为单元素数组。加载方法与宿主约定见 [AviSynth 接口](README.md#avisynth-调用与构建)。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_Deen(src, mode="w3d", rad=2, thrY=7, tthY=4, planes=[0])
```
