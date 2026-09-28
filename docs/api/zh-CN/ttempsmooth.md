# TTempSmooth

[返回 API 目录](README.md)

```python
core.neo_smo.TTempSmooth(clip, maxr=3, thresh=[4, 5, 5], mdiff=[2, 3, 3],
                       strength=2, scthresh=12.0, fp=True, pfclip=None, planes=None)
```

接受 Gray/RGB/YUV 的 8–16 位整数与 F32，**不接受 F16**。

| 参数 | 范围 | 含义 |
|---|---|---|
| maxr | 1–7 | 最大时间半径 |
| thresh | 每项 1–256 | 8 位标度的差异阈值，逐平面数组 |
| mdiff | 每项 0–255 | 在差异权重开始衰减前的区间，逐平面数组 |
| strength | 1–8 | 时间权重中保持高权重的范围 |
| scthresh | −1–100 | 场景切换控制，见下文 |
| fp | True/False | 剩余权重的处理方式，**不是浮点精度开关** |
| pfclip | 默认主片 | 预滤参考片，用于差异判断，输出像素仍取自主片 |
| planes | 默认全部 | 处理平面 |

thresh/mdiff 的短数组重复最后一项；没有 `scalep` 参数，始终使用其规定的 8 位标度。`pfclip` 须匹配主片格式和尺寸，并与主片时间对齐。

`scthresh>0` 自动调用 `misc.SCDetect`，其阈值为 `scthresh/100`，有 pfclip 时在参考片上检测；RGB 不支持这条自动检测路径。0 完全关闭场景检查；负值读取已有场景属性，推荐写 −1，属性缺失时不会得到切点。默认值 12.0 因而需要安装提供 SCDetect 的 misc 插件。

`fp=True` 将没有使用的归一化权重补到中心像素；`False` 按实际累加权重重新归一化。它改变滤镜行为，不是速度预设。

时间越界帧号限制到首尾帧。参照帧从近到远检查，某方向遇到不满足差异条件的样本后，不继续把更远样本当成独立候选。

```python
# 自包含调用，不依赖 SCDetect
out = core.neo_smo.TTempSmooth(src, maxr=3, scthresh=0)
# 使用空间预滤结果引导时间权重
guide = core.neo_smo.Median(src, radius=[1])
out = core.neo_smo.TTempSmooth(src, pfclip=guide, scthresh=0)
```

原理见[时域去噪与修复](../../knowledge/zh-CN/temporal.md)。
