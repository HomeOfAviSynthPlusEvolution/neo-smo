# Cnr4

[返回 API 目录](README.md)

```python
core.neo_smo.Cnr4(clip, mode="oxx", radius=2, sense=[35, 47, 47],
                 str=[192, 255, 255], pow=[1.0, 1.0, 1.0],
                 tmode=0, wmode=0, scenechange=True, ref=None)
```

仅接受 YUV 8–16 位整数。只修改 U/V，Y 沿用主片。没有 `planes` 参数；`mode` 也不是平面开关。

| 参数 | 范围 | 含义 |
|---|---|---|
| mode | 恰好三个小写 `o` 或 `x` | Y/U/V 差异查表曲线的类型 |
| radius | 1–10 | 时间半径 |
| sense | 恰好三项，每项 −1–255 | 差异曲线尺度；−1 保留该位置默认值 |
| str | 恰好三项，每项 −1–255 | 查表幅度与表的有效范围；−1 保留该位置默认值 |
| pow | 恰好三项，有限且 ≥0 | 差异曲线指数参数 |
| tmode | 0–4 | 单遍或多遍处理方式，见下表 |
| wmode | 0–3 | 时间距离权重 |
| scenechange | 默认 True | 使用已有场景切换属性，不自动调用检测器 |
| ref | 默认主片 | 提供 Y/U/V 差异信息，须匹配主片格式与尺寸 |

| tmode | 行为 |
|---|---|
| 0 | 直接使用时间邻帧做单遍混合 |
| 1 | 从窗口两端向内预处理，局部半径 1，更新参考结果 |
| 2 | 同 1，同时更新参与后续处理的源像素 |
| 3 | 预处理半径逐步扩大，更新参考结果 |
| 4 | 同 3，同时更新源像素 |

wmode 0 等权，1 使用平方根距离权重，2 使用正弦距离权重，3 按距离倒数衰减。它与 mode 的曲线选择是独立参数。

默认 scenechange=True 要求主片具有 `_SceneChangePrev` 和 `_SceneChangeNext` 属性。没有预先检测时可显式关闭；开启时按主片属性判定窗口边界，并用中心替代被排除的样本。片头片尾重复边界帧。参考亮度需与色度尺寸对齐，内部使用 `std`/`resize` 构造参考路径。

```python
out = core.neo_smo.Cnr4(src, scenechange=False)
# 有 misc 插件时也可先生成场景属性：
sc = core.misc.SCDetect(src, threshold=0.1)
out = core.neo_smo.Cnr4(sc)
```

曲线与参考片作用见[色彩去噪原理](../../knowledge/zh-CN/chroma.md)。
