# FluxSmoothT、FluxSmoothST

[返回 API 目录](README.md)

```python
core.neo_smo.FluxSmoothT(clip, temporal_threshold=None, planes=None, scalep=False)
core.neo_smo.FluxSmoothST(clip, temporal_threshold=None, spatial_threshold=None,
                        planes=None, scalep=False)
```

支持通用格式。阈值均为逐平面浮点数组，最后一项重复；省略值均为自动缩放的 8 位标度 7。显式值默认按原生单位解释，`scalep=True` 时非负值须在 0–255 内并缩放。负阈值停用相应分量，0 仍是有效阈值。

只在当前像素严格大于前后帧，或严格小于前后帧时尝试平滑。FluxSmoothT 检查两个时间邻点；FluxSmoothST 还检查当前帧的八个空间邻点。与中心的绝对差不大于对应阈值的邻点才参与平均，中心始终参与。

FluxSmoothST 即便将时间阈值设为负值、只累加空间邻点，仍保留上述时域极值门槛，并不变成普通空间均值滤镜。首尾各一帧复制；没有场景切换检测。

```python
out = core.neo_smo.FluxSmoothT(src, temporal_threshold=[7], scalep=True)
out = core.neo_smo.FluxSmoothST(src, temporal_threshold=[7],
                              spatial_threshold=[5], planes=[0], scalep=True)
```

与 TemporalSoften 的分母区别见[时域原理](../../knowledge/zh-CN/temporal.md)。
