# Clense、ForwardClense、BackwardClense

[返回 API 目录](README.md)

```python
core.neo_smo.Clense(clip, previous=None, next=None, planes=None)
core.neo_smo.ForwardClense(clip, planes=None)
core.neo_smo.BackwardClense(clip, planes=None)
```

支持通用格式。`planes` 省略时处理全部平面。

| 接口 | 处理第 n 帧时的参考 | 不处理的时间边界 |
|---|---|---|
| Clense | `previous[n-1]`、`next[n+1]`；两片分别默认 `clip` | 首尾各 1 帧 |
| ForwardClense | `clip[n+1]`、`clip[n+2]` | 最后 2 帧 |
| BackwardClense | `clip[n-1]`、`clip[n-2]` | 最前 2 帧 |

可选参考片必须满足通用匹配约束。`previous`、`next` 是视频节点，接口内部已经施加 ±1 帧偏移；不要为了同一语义再手动偏移一次。

Clense 将当前值限制在前、后参考值之间，等价于三个值的中值。ForwardClense 与 BackwardClense 根据近、远两个参考值做单向外推限幅，**不是三个值的中值**。三者均不检查场景切换属性。

```python
out = core.neo_smo.Clense(src, planes=[0])
out = core.neo_smo.ForwardClense(src)
```

公式见[时域去噪与修复](../../knowledge/zh-CN/temporal.md)。
