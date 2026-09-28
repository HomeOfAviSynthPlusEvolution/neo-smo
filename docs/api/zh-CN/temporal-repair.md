# TemporalRepair

[返回 API 目录](README.md)

```python
core.neo_smo.TemporalRepair(clip, repairclip, mode=[0], planes=None)
```

支持通用格式。`repairclip` 必填且须匹配主片格式和尺寸。使用主片当前帧，以及参考片的前、当前、后三帧。首尾各一帧直接复制主片；没有场景切换参数。

`mode` 为逐平面整数数组，范围 0–4，最后一项重复。**mode 0 是有效修复模式，不是关闭。** 保留平面应通过 `planes` 选择。

| mode | 约束来源 |
|---|---|
| 0（默认） | 参考三帧同坐标像素的最小值与最大值 |
| 1 | 参考八邻点的时域正、负变化分别扩展当前参考值，再纳入前后帧中心 |
| 2 | 参考 3×3 内最大的时域绝对变化，构造对称范围 |
| 3 | 分别求参考 3×3 到前帧、后帧的最大差，取两者较小值构造对称范围 |
| 4 | 仅使用三帧中心值，按时域趋势构造限幅范围 |

模式 1–3 的空间越界坐标采用镜像。模式不是递增强度等级。

```python
filtered = core.neo_smo.Median(src, radius=[1])
out = core.neo_smo.TemporalRepair(filtered, src, mode=[0], planes=[0])
```

原理见[时域去噪与修复](../../knowledge/zh-CN/temporal.md)。
