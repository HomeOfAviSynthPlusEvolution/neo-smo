# Repair

[返回 API 目录](README.md)

```python
core.neo_smo.Repair(clip, repairclip, mode)
```

支持通用格式。`clip` 是待修复结果，`repairclip` 提供同一帧的参考邻域；两者格式和尺寸须匹配。`mode` 是必填逐平面数组，0–24，最后一项重复。mode 0 复制 `clip` 对应平面。没有 `planes` 参数。

| mode | 行为概要 |
|---|---|
| 1、11 | 用参考 3×3 的最小／最大值限制主片 |
| 2–4 | 对参考九点排序，以内侧秩作为上下限 |
| 5–9 | 选择参考邻域的方向范围，以不同代价决定方向 |
| 10 | 根据主片与参考样本的差异选择参考值 |
| 12–14 | 排序参考八邻点，选内侧范围，再把参考中心纳入范围 |
| 15、16 | 按参考中心计算方向代价，将选中范围扩展到包含参考中心，再限制主片 |
| 17、18 | 参考方向范围的组合／选择，并纳入参考中心 |
| 19、20 | 以参考中心为基准，用其到八邻点的最小／第二小绝对差构造对称范围 |
| 21 | 以参考中心为基准，使用各对向组最大差中的最小值构造对称范围 |
| 22–24 | 将类似 19–21 的距离基准改为主片值，并把参考中心限制到该范围 |

模式含义与 RemoveGrain 并不逐项相同；尤其 Repair 11 不是高斯加权平均。空间边缘使用镜像取样。

```python
filtered = core.neo_smo.RemoveGrain(src, mode=[2])
out = core.neo_smo.Repair(filtered, src, mode=[1])
```

此例限制已处理结果相对原片的偏离；参考片不是输出中简单混合的一层。原理见[空间排序与限幅](../../knowledge/zh-CN/spatial.md)。

## AviSynth

```text
neo_smo_Repair(clip, repairclip, mode)
```

参数名称、顺序与上方共用说明一致；宿主格式、数组写法及音频/场序保留见 [API 目录](README.md#avisynth-调用与构建)。

`repairclip` 是必填剪辑，格式与尺寸须匹配主片；`mode` 是必填整数或数组。参考片按同一帧号读取，不自动对齐时间。没有 `planes` 参数，mode 0 复制主片对应平面。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
filtered = neo_smo_RemoveGrain(src, mode=[2])
return neo_smo_Repair(filtered, src, mode=[1])
```
