# 时间窗口与场景切换

[返回知识目录](../README.md)

## 片头片尾不是统一处理

| 滤镜 | 时间范围 | 不足窗口时 |
|---|---|---|
| Clense、TemporalRepair、DegrainMedian、FluxSmoothT/ST | n−1、n、n+1 | 首尾各一帧复制主片 |
| ForwardClense | n、n+1、n+2 | 最后两帧复制 |
| BackwardClense | n−2、n−1、n | 最前两帧复制 |
| TemporalMedian、CCD | n−r … n+r | 首尾各 r 帧复制；r 分别为 radius、temporal_radius |
| TemporalSoften | n−r … n+r | 缩短为实际存在的帧 |
| TTempSmooth、Cnr4 | 参数指定的时间范围 | 越界帧号限制到首尾帧，可能重复使用边界帧 |

例如五帧片段使用 TemporalMedian radius=3 时，每帧都没有完整窗口，所以都保持原片。TemporalSoften 在同样情况下仍可以处理。

## 场景属性

`_SceneChangePrev` 表示当前帧与前帧之间有切点；`_SceneChangeNext` 表示当前帧与后帧之间有切点。属性属于帧，不是自动附带的检测结果。接口是否生成、要求或忽略它们，必须分别看待。

| 滤镜 | 默认 | 开启方式 | 属性缺失 |
|---|---|---|---|
| TemporalMedian | 关闭 | scenechange=True 读取主片属性 | 检查到所需属性缺失时报错 |
| TemporalSoften | 关闭 | −1 读已有属性；1–254 自动 SCDetect，阈值 /255 | −1 模式下缺失属性不提供切点 |
| TTempSmooth | 自动检测 | scthresh>0 自动 SCDetect，阈值 /100；负值读已有属性 | 读取模式下缺失属性不提供切点 |
| Cnr4 | 读取已有属性 | scenechange=True，检查主片属性 | 检查到所需属性缺失时报错 |
| 其余接口 | 无场景处理 | 无开关 | 不读取 |

TTempSmooth 在有 pfclip 时使用该引导路径的场景属性。Cnr4 的场景属性来自主片，不因设置 ref 而改为只检查 ref。

TemporalMedian/Soften/TTempSmooth 会根据切点限制候选时间范围；Cnr4 保持窗口结构，将被边界判定排除的样本替换为中心。不同窗口形状与替代方式会造成输出差异，不能将这些选项视为同一个通用布尔开关。

RGB 不能直接走 TemporalSoften/TTempSmooth 的自动 SCDetect 路径。若要在 RGB 上使用已有属性，应先在合适的检测素材上生成并正确传递场景属性，再选读取模式。

自动检测在 VapourSynth 中调用 `misc.SCDetect`；AviSynth+ 内置按亮度平面归一化平均绝对帧差的检测，使用相同阈值单位和属性名。自动检测要求至少两帧，RGB 仍需外部提供场景属性。
