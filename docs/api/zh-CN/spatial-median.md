# Median、InterQuartileMean、SmartMedian

[返回 API 目录](README.md)

```python
core.neo_smo.Median(clip, radius=[1], planes=None)
core.neo_smo.InterQuartileMean(clip, radius=[1], planes=None)
core.neo_smo.SmartMedian(clip, radius=[1], threshold=None, scalep=False, planes=None)
```

支持通用格式。三者使用 `(2r+1) × (2r+1)` 空间邻域；`radius` 是逐平面整数数组，每项范围 0–3，省略时为 1，0 表示复制该平面。`planes` 再限制实际处理的平面。

| 接口 | 输出含义 |
|---|---|
| Median | 包含中心像素的邻域中值 |
| InterQuartileMean | 排序后中间 50% 的加权均值，四分位边界按部分权重计入 |
| SmartMedian | 排除中心后求中间两个值；邻域满足平坦度条件时，将中心限制在两者之间，否则保持原值 |

SmartMedian 的 `threshold` 为逐平面浮点数组。省略时，半径 1 使用 8 位标度 50，半径 2/3 使用 128，均自动缩放到格式。显式指定时：`scalep=True` 要求 0–255；`False` 要求整数格式的 0–最大值，浮点格式的 0–1。阈值越高，越容易触发修正；它不是“中心与中值的最大差”。

边缘采用镜像取样，极小平面的越界坐标再限制到有效范围；不会因为处于边缘就整圈保留原像素。

```python
# 只处理亮度
out = core.neo_smo.Median(src, radius=[1], planes=[0])
# 各平面均做四分位均值
out = core.neo_smo.InterQuartileMean(src, radius=[2])
# 明确使用 8 位阈值标度
out = core.neo_smo.SmartMedian(src, radius=[1], threshold=[50], scalep=True)
```

原理与计算例子见[空间排序与限幅](../../knowledge/zh-CN/spatial.md)。
