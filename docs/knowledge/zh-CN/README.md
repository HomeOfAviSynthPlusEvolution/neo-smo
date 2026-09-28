# neo-smo 算法知识

本目录说明算法如何选择样本、决定是否修改像素，以及数值与边界约定。参数签名、默认值和可直接使用的调用示例见 [API 目录](../../api/zh-CN/README.md)。

| 主题 | 内容 |
|---|---|
| [空间排序与限幅](spatial.md) | Median、InterQuartileMean、SmartMedian、VerticalCleaner、RemoveGrain、Repair |
| [时域去噪与修复](temporal.md) | Clense 家族、TemporalMedian/Soften/Repair、DegrainMedian、FluxSmooth、TTempSmooth |
| [块 DCT 与频率加权](dct-filter.md) | DCTFilter 的 8×8 块、系数外积、补边和归一化 |
| [色彩去噪](chroma.md) | CCD 的多通道空间取样、Cnr4 的时域混合 |
| [采样、阈值与精度](shared/sample-and-precision.md) | 平面、位深、scalep、F16/F32、FMA |
| [时间窗口与场景切换](shared/temporal-boundaries.md) | 片头片尾、场景属性、不同滤镜的窗口规则 |

## 如何理解这些算法

空间滤镜从同一帧取邻点；时间滤镜从相邻帧取同坐标像素。时间半径扩大并不意味着运动搜索范围扩大：neo-smo 这些接口不做运动估计。运动较强时，时间样本可能来自不同物体，阈值与场景边界负责减少不合适的混合，但不能替代运动补偿。

中值主要抑制孤立极值；平均可以降低随机波动，也会损失细节；限幅只将越过允许范围的值拉回边界。Repair 一类算法由参考片构造允许范围，作用与直接模糊主片不同。

本文公式先解释算法结构；整数的饱和、取整以及浮点舍入仍以实现为准。模式号、半径和阈值不是跨滤镜统一的“降噪强度”。
