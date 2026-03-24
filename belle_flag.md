# Belle 事件分类标志说明

## evtcls_flag[20] - 事件类型分类标志

用于标识事件属于哪种物理过程类型:

| 索引 | 事件类型 | 说明 |
|------|---------|------|
| `evtcls_flag[0]` | 强子事件 | Hadronic events |
| `evtcls_flag[1]` | Bhabha散射 | e⁺e⁻ → e⁺e⁻ |
| `evtcls_flag[2]` | 双光子事件 | γγ pair |
| `evtcls_flag[3]` | μ对产生 | μ⁺μ⁻ pair |
| `evtcls_flag[4]` | τ对产生 | τ⁺τ⁻ pair |
| `evtcls_flag[5]` | 双光子过程 | Two-photon process |
| `evtcls_flag[6]` | 宇宙线事件 | Cosmic ray |
| `evtcls_flag[7]` | 辐射Bhabha | Radiative Bhabha |
| `evtcls_flag[8]` | eeee末态 | Four electron final state |
| `evtcls_flag[9]` | 量能器Bhabha | Calorimeter Bhabha |
| `evtcls_flag[10]` | γ+φ事件 | Gamma + phi |
| `evtcls_flag[11]` | 辐射μ对 | Radiative mu pair |
| `evtcls_flag[17]` | 束流气体本底监测 | Beam gas monitoring |
| `evtcls_flag[18]` | 束流本底 | Beam background |
| `evtcls_flag[19]` | 未知/垃圾事件 | Unknown/junk events |

### 标志值含义

| 值 | 含义 |
|----|------|
| `10` | 满足基本选择条件 |
| `20` | 满足较松选择,但不满足严格条件 |
| `30` | 满足更严格的选择条件 |
| `40` | 满足最严格的选择条件 |

---

## hadronic_flag[20] - 强子事件子分类标志

专门用于强子事件的进一步细分(仅在 `evtcls_flag[0] >= 10` 时有效):

| 索引 | 子类别 | 说明 |
|------|--------|------|
| `hadronic_flag[0]` | 基本强子选择 | Basic hadronic selection |
| `hadronic_flag[1]` | HadronA选择 | 更严格的ECL和径迹要求 |
| `hadronic_flag[2]` | HadronB选择 | Brendan Casey的选择标准 |
| `hadronic_flag[3]` | HadronA + Ntrk≥5 | HadronA且至少5条径迹 |
| `hadronic_flag[4]` | HadronB + Ntrk≥5 | HadronB且至少5条径迹 |
| `hadronic_flag[5]` | 非常严格的强子选择 | Very tight hadronic selection |

---

## 使用示例

```cpp
// 检查是否为强子事件
if (evtcls_flag[0] >= 10) {
    // 是强子事件
    
    // 进一步检查强子子类别
    if (hadronic_flag[1] >= 10) {
        // 满足HadronA选择标准
    }
    
    if (hadronic_flag[2] >= 10) {
        // 满足HadronB选择标准
    }
}

// 检查是否为Bhabha事件
if (evtcls_flag[1] >= 10) {
    // 是Bhabha事件
}
```

