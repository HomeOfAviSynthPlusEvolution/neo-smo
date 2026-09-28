# CCD

[返回 API 目录](README.md)

```python
core.neo_smo.CCD(clip, threshold=4.0, temporal_radius=0,
                points=[1, 1, 0], scale=None, ref=None)
```

接受 RGB/YUV 的 8–16 位整数、F16、F32（F16 仅限 VapourSynth），不接受 Gray。RGB 处理三平面；YUV 只处理 U/V，Y 保持主片内容。没有 `planes` 或 `scalep` 参数。

| 参数 | 范围与默认 | 含义 |
|---|---|---|
| threshold | 默认 4.0，必须有限且缩放后平方有限 | 控制多通道距离判定；内部按平方使用 |
| temporal_radius | 0–10，默认 0 | 扩大差异估计的时间窗口 |
| points | 恰好三个布尔整数（0 关闭、非零开启），默认 `[1,1,0]`，至少一个启用 | 选择近、中、远三组空间采样点 |
| scale | 有限且 ≥1；默认色度平面高度 / 240（RGB 为图像高度 / 240） | 缩放采样点间距 |
| ref | 默认主片 | 提供差异判断的参考片，须匹配格式和尺寸 |

近组为距离 4 的四角点，中组为距离 8 的八点，远组为外圈距离 12 的十二点；具体布局见[色彩去噪原理](../../knowledge/zh-CN/chroma.md)。这些不是连续方形卷积窗口。

默认 scale 可能小于 1 并导致创建失败，例如高度不足 480 的 YUV420。可显式指定 `scale=1.0`，但图像仍须足够容纳选定几何范围。当前实现还会根据最大点跨度与 scale 检查宽高，较大的 scale 更容易超过尺寸限制。

YUV 抽样格式的参考亮度通过 Bilinear 缩放到色度尺寸，VapourSynth 使用标准 `std`/`resize` 插件，AVS 使用内置缩放。前后各 `temporal_radius` 帧直接保留主片原样，也不做空间滤波。如果所有帧都没有完整时间窗口，整段视频保持不变。没有场景切换参数。

`temporal_radius` 增大的是判定所用的时间信息；被接受的空间点贡献的是**当前帧主片像素**，不是把各时间帧直接平均输出。负 threshold 当前也会被平方，因此与相反数相同；正常调用使用非负值即可。

```python
out = core.neo_smo.CCD(src, threshold=4.0, temporal_radius=0, scale=1.0)
```

## AviSynth

```text
neo_smo_CCD(clip, threshold=4.0, temporal_radius=0, points=[1, 1, 0], scale=Undefined(), ref=Undefined())
```

参数名称、顺序与上方共用说明一致；宿主格式、数组写法及音频/场序保留见 [API 目录](README.md#avisynth-调用与构建)。

AVS 接受平面整数 RGB/YUV 与 F32，不支持 F16。`points` 必须为三个整数，不能写布尔数组；`ref` 接受剪辑。没有 `planes` 或 `scalep` 参数。色度抽样输入的参考亮度使用内置 `ExtractY` / `BilinearResize`；宿主缩放取整差异可能改变邻居选择，跨宿主不保证逐位一致。下例显式设置 scale，避免低分辨率输入的默认值小于 1。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_CCD(src, threshold=4.0, temporal_radius=0, points=[1, 1, 0], scale=1.0)
```
