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

## AviSynth

```text
neo_smo_Clense(clip, previous=Undefined(), next=Undefined(), planes=Undefined())
neo_smo_ForwardClense(clip, planes=Undefined())
neo_smo_BackwardClense(clip, planes=Undefined())
```

参数名称、顺序与上方共用说明一致；宿主格式、数组写法及音频/场序保留见 [API 目录](README.md#avisynth-调用与构建)。

Clense 的 `previous`、`next` 接受剪辑，省略时各自使用主片；内部读取的是 n−1 / n+1，调用方不应额外偏移。三者的 `planes` 接受整数或数组，均不读取场景切换属性。

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_Clense(src, planes=[0])
```

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_ForwardClense(src)
```

```avs
LoadPlugin("/path/to/neo-smo.dll")
src = BlankClip(width=640, height=480, length=24, pixel_type="YUV420P10")
return neo_smo_BackwardClense(src)
```
