# VerticalCleaner

[返回 API 目录](README.md)

```python
core.neo_smo.VerticalCleaner(clip, mode)
```

支持通用格式。`mode` 是必填逐平面整数数组，取值 0–2，最后一项重复到剩余平面。没有 `planes` 参数。

| mode | 行为 | 尺寸与边界 |
|---|---|---|
| 0 | 复制 | 无额外高度要求 |
| 1 | 上一行、当前行、下一行的中值 | 被处理平面至少 3 行；首尾各 1 行复制 |
| 2 | 根据上下各两行的趋势构造范围，再限制当前值 | 被处理平面至少 5 行；首尾各 2 行复制 |

检查的是平面高度；例如 YUV420 色度高度只有亮度的一半。该滤镜不沿水平方向采样，也不跨帧。

```python
out = core.neo_smo.VerticalCleaner(src, mode=[1, 0, 0])
```

模式 2 不是五点中值，见[空间排序与限幅](../../knowledge/zh-CN/spatial.md)。
