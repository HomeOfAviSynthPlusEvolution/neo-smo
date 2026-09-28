# 从 zsmooth 迁移

[返回 API 目录](README.md)

neo-smo 当前提供本目录所列 19 个接口。迁移时先将调用命名空间改为 `core.neo_smo`，再逐项对照签名与格式要求；不要假定上游的全部接口和辅助选项已经注册。

```python
# 已有脚本的空间中值调用可以写为：
out = core.neo_smo.Median(src, radius=[1])
```

## 需要确认的行为

- Repair 与 TemporalRepair 的第二个输入都是 `repairclip`。Repair mode 0 复制；TemporalRepair mode 0 处理。
- `scalep` 默认关闭；显式阈值与省略阈值可能不同。跨位深复用时应明确指定单位，尤其是 DegrainMedian。
- TemporalMedian、TemporalSoften、TTempSmooth、Cnr4 的场景参数有不同含义和默认值，按各专页迁移。
- F16 输入不是 F32 输入的缩写；支持原生半精度的平台可使用半精度运算，不支持的平台使用 F32 运算后存回 F16。DCTFilter 是明确的例外：F16 输入也固定用 F32 做变换，再写回 F16。不能以跨平台逐位相同作为兼容承诺。
- 生产路径允许 FMA 和通常的舍入差异。接近阈值或排序边界时，微小差异可能改变分支选择，验收应观察差异幅度及实际输出影响。
- 当前没有公开的 `opt`、SIMD OFF 或 KernelInfo 接口，尚无 AviSynth 接口。DCTFilter 已接入，要求恰好八个有限且在 `[0,1]` 内的系数；参见 [API](dct-filter.md)。

这些文档描述当前实现，不宣称对所有历史 zsmooth 版本逐位复现。原理与数值约定见[知识目录](../../knowledge/zh-CN/README.md)。
