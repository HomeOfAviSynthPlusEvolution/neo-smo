# DCTFilter

[返回 API 目录](README.md)

```python
core.neo_smo.DCTFilter(clip, factors, planes=None)
```

对每个选中平面的固定 8×8 块执行 DCT、频率系数加权和逆变换。它只处理当前帧，不做运动补偿或时间平均。支持固定格式、固定尺寸的 Gray、RGB、YUV，采样类型为 8–16 位整数、F16、F32（F16 仅限 VapourSynth）。输出格式、尺寸和帧属性保持不变。

## 参数

| 参数 | 要求与含义 |
|---|---|
| `clip` | 必填输入视频 |
| `factors` | 必填且恰好 8 个浮点数，每项必须有限且在 `[0,1]` 内；没有默认值 |
| `planes` | 省略时处理全部平面；显式传入非空平面编号数组，例如 `[0]` 或 `[1,2]`；重复、越界和空数组会报错 |

`factors` 从直流到最高频率排列。横向频率 u、纵向频率 v 对应的二维权重是 `factors[v] * factors[u]`，并非给 8 个空间像素分别设置强度；所有选中平面共用这一组系数。数组不会自动补齐，也没有 `scalep`、可调块尺寸或 `opt` 参数。

- `[1] * 8` 表示理想数学意义上的恒等变换，但仍执行 DCT 往返，不承诺逐位不变。
- `[1,0,0,0,0,0,0,0]` 只保留每块的平均值。
- `[0] * 8` 使选中平面的输出为数值 0；整数 YUV 色度的 0 不是中性色，通常不应对色度这样使用。
- 直流权重实际为 `factors[0] ** 2`；希望保留块平均值时设 `factors[0]=1`。

## 边界与精度

块网格从各平面左上角开始，块之间不重叠。宽高不必是 8 或 16 的倍数。右侧和底部缺少的样本按上游 `resize.Point` 的实际行为补齐：镜像一次，重复边缘像素；延伸超过整个原平面时再限制到另一侧边缘。输出只写回原有范围，未选择的平面保留原内容。

整数和浮点输入都以 F32 执行 DCT。F16 输入先转 F32，仅在写出时转回 F16；这里不会切换成原生半精度 DCT。生产允许 FMA。

整数输出限制到对应位深的合法范围后四舍五入，正数恰好半整数时向上取整。浮点输出不强制裁剪到亮度/RGB 的 `[0,1]` 或色度的 `[-0.5,0.5]`。

## 示例

```python
# 衰减较高频率，只处理亮度
out = core.neo_smo.DCTFilter(
    src, factors=[1, .9, .7, .5, .3, .2, .1, 0], planes=[0])

# 保留每个 8×8 亮度块的平均值，便于观察块网格
blocks = core.neo_smo.DCTFilter(
    src, factors=[1, 0, 0, 0, 0, 0, 0, 0], planes=[0])
```

原理与数值约定见[块 DCT 与频率加权](../../knowledge/zh-CN/dct-filter.md)。

## AviSynth

```text
neo_smo_DCTFilter(clip, factors, planes=Undefined())
```

参数名称、顺序与上方共用说明一致；宿主格式、数组写法及音频/场序保留见 [API 目录](README.md#avisynth-调用与构建)。

`factors` 必须为恰好八项的数值数组，不自动重复单值；AVS 示例须显式列出八项，不能照抄 Python 的 `[1] * 8` 写法。`planes` 接受整数或数组。AVS 支持平面整数与 F32，不支持 F16；DCT 仍以 F32 计算。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_DCTFilter(src, factors=[1, 0.9, 0.7, 0.5, 0.3, 0.2, 0.1, 0], planes=[0])
```
