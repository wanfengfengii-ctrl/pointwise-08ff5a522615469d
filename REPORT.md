# 同批次重叠重采样趋势分析事故诊断报告

本报告只诊断，不修复。原始实现（`matrix.py`、`resample.py`、`summary.py`、`trend.py`、`pipeline.py`、`checks.py`）与 `data/`、`logs/` 均未改动。所有数值由 `analysis.py` 从 `data/incident.json` 的原始观测用精确有理数独立重建；`python3 analysis.py verify` 可重算、核对本报告锚点并重放原入口。

---

## 1. 结论摘要

1. **根因**：`pipeline.py:15` 按计划把窗口逐块交给 `summary.summarize`，而 `summary.summarize` 只在自己的块内选取线性无关行、只用块内窗口构造协方差并求逆（`summary.py:8-14`）；`summary.merge` 再把各块的 information、rhs 与 rank **直接相加**（`summary.py:30-36`）。这等于把同一批次的不同分块当成彼此独立的观测批次：块间共享源导出的协方差被丢掉，块内已消掉的相关在同一块内被算作独立信息，块间线性相关也被重复计秩。
2. **规范口径**：整批 8 个窗口的误差协方差秩为 **6**（不是 8），完整批次 GLS 的斜率为 `72028/36005`，斜率方差 `520/7201`，标准误 `0.26872326425340776`，95% 区间 `[1.47380233262852, 2.5271975285018784]`。任何极大线性无关窗口子集都给出同一解（17 个基逐一核对，全部相同）。
3. **点估计掩盖**：五个计划的斜率只在 `2.0004999…` 到 `2.0010334…` 之间变动（最大差 `10933009/20493109870 ≈ 5.335e-4`）。原因不是估计正确，而是窗口值几乎落在 `y=1+2t` 上，任何满足线性无偏约束（负载荷和 = 0、时间加权和 = 1）的权重都给出几乎相同的斜率。
4. **区间变化**：报告方差从 `all` 的 `520/7201`（`0.0722122`）最低降到 `single`/`striped` 的 `13930/284587`（`0.0489481`，仅为规范的 67.78%），报告半宽从 `0.5266976` 缩到 `0.4336348`（`0.82331` 倍）。但**现有估计量自身的真实方差**在同一计划下是 `6719310430/80989760569`（`0.0829649`），是它报告值的 `1.6949564` 倍，也已是规范 GLS 方差的 114.89%。报告区间与估计量本身都对不上，属于两种不同的量，不能互相代替。
5. **秩被高估**：`all` 报 6（碰巧正确），`cut` 报 7，`pairs`/`striped`/`single` 报 8；规范值一律为 6。高估量恰好等于该计划各块基大小之和减去全局秩。
6. **单因素对照**：仅换分块即可复现漂移；追加一个与 `w3` 共用 `b`、`c` 且权重相同的重复窗口在单窗口分块下把报告方差**再压低**到 `16010/334909`（`0.0478040`）并把秩抬到 9，而真实信息增量为 0；把每个窗口的引用源换成同数值的独立新 id 后，五种计划全部回到同一规范解（方差 `13930/284587`），且报告方差 = 实际方差。三组对照把原因锁定为“跨块共享源被当作新增独立测量 / 块间相关被丢弃”。

---

## 2. 规范口径下的完整批次（独立核定）

源观测 6 个（`a…f`），窗口 8 个（`w0…w7`），权重满足合同：`w0={a:1}`、`w1={a:1/2,b:1/2}`、`w2={b:1}`、`w3={b:1/2,c:1/2}`、`w4={c:1/2,d:1/2}`、`w5={d:1/2,e:1/2}`、`w6={e:1/2,f:1/2}`、`w7={f:1}`。窗口中共享源使权重行线性相关，例如整批的零空间含

- `w0 − 2·w1 + w2 = 0`
- `w0 − 2·w1 + 2·w3 − 2·w4 + 2·w5 − 2·w6 + w7 = 0`

因此 8×8 窗口协方差矩阵奇异，秩 = 6。按合同，`information_rank` 定义为完整窗口误差协方差的秩，即 6。

取任一极大线性无关子集（例如 `w0,w1,w3,w4,w5,w6`）做 GLS，窗口值为 `w0=1, w1=2, w2=3, w3=4, w4=121/20, w5=161/20, w6=10, w7=11`，窗口时间为 `0, 1/2, 1, 3/2, 5/2, 7/2, 9/2, 5`。结果：

| 量 | 精确值 | 浮点 |
|---|---|---|
| 规范 information `XᵀΣ⁻¹X` | `[[65/18, 103/12],[103/12, 137/4]]` | — |
| 规范 rhs `XᵀΣ⁻¹y` | `[[1871/90],[4627/60]]` | — |
| 协方差 `(XᵀΣ⁻¹X)⁻¹` | `[[4932/7201, −1236/7201],[−1236/7201, 520/7201]]` | — |
| 截距 | `36073/36005` | `1.0018886265796416` |
| 斜率 | `72028/36005` | `2.000499930565199` |
| 截距方差 | `4932/7201` | `0.6849048743230107` |
| 斜率方差 | `520/7201` | `0.0722121927510068` |
| 标准误 | — | `0.26872326425340776` |
| 95% 区间 | — | `[1.47380233262852, 2.5271975285018784]` |
| 信息秩 | 6 | — |

17 个极大无关基逐一核对，`beta`、协方差、源负载系数完全相同（`all_bases_equivalent = True`），说明“选哪个极大无关子集”不影响规范解，符合合同对伪逆实现不敏感的要求。

---

## 3. 源码根因：首次偏离业务口径的中间状态

按数据流定位：

1. `pipeline.py:15` 把每个计划块单独交给 `summarize(source, [lookup[wid] for wid in group])`。`source` 虽是整批，但 `summarize` 只使用本块窗口。
2. `summary.py:8-11` 在**块内**做贪心线性无关选择，得到 `basis_windows`。块内秩正确不代表并集秩正确。
3. `summary.py:12-13` 只用本块 `selected` 行构造协方差 `C_b = W_b diag(v) W_bᵀ`。跨块窗口之间由共享源产生的非零协方差在这里被截断。
4. `summary.py:14` 对本块 `C_b` 求逆，得到**块内边际精度**；`summary.py:17-20` 用该精度算块级 `information = X_bᵀ C_b⁻¹ X_b`、`rhs = X_bᵀ C_b⁻¹ y_b`。
5. `summary.py:30-36` 的 `merge` 把块级 `information`、`rhs` **逐元素相加**，并把 `information_rank` 相加为 `sum(b["information_rank"])`（第 34 行）。这正是业务口径被破坏的地方：它等价于假设块与块相互独立、且块内精度就是全批精度。
6. `trend.py:5-8` 直接对相加后的 `information` 求逆并取 `[1][1]` 作为斜率方差，因此报告方差 = `(Σ_b X_bᵀ C_b⁻¹ X_b)⁻¹[1][1]`，而不是规范 `(Xᵀ Σ_full⁻¹ X)⁻¹[1][1]`。

**首次偏离业务口径的中间状态**：`summary.py:12-14` 的“仅块内协方差 + 块内求逆”。`pipeline.py:15` 制造了分块输入，`summary.py:12` 起信息不再依赖整批误差协方差，`summary.py:32-35` 的相加把该偏差固化为最终结果，`summary.py:34` 的秩相加进一步把“块内独立”误当成“全局独立”。

核对：`checks.py:13-21` 的 `normal` 只要求三计划各自满足固定期望值且彼此相等，`checks.py:23-39` 的 `local` 只检查**块内**协方差对称、可逆和信息矩阵非负——它们都看不到跨块协方差被截断，所以 `logs/local.log` 会显示 `{"checked_block_count": 26, "local_covariance_checks": "passed"}`。这也解释了为什么局部检查通过却不构成完整批次正确性的保证。

---

## 4. 必要触发条件

必须同时满足：

1. **同一源 id 被不同块引用**（跨块共享）。只有块内有共享时，块内协方差仍能正确传播相关，块级估计是无偏的；一旦共享跨过块边界，`merge` 就把相关当独立。
2. **分块非平凡**，即至少两个块，且跨块窗口的权重行并集线性相关（本例整批秩 6 < 窗口数 8）。这使“块内各自满秩、并集降秩”，从而 `sum(rank_b) > rank_union`。
3. **报告方差取自 `inv(Σ_b J_b)`**，且不做跨块协方差修正（`summary.py:30-36` + `trend.py:5-8`）。这是区间漂移的直接通道。

本包中 `all` 只有一块，不触发条件 2，故恰好等于规范解；`cut`、`pairs`、`striped`、`single` 全部触发。具体地：

- `cut`：块1 `{w0,w1,w2,w3}`、块2 `{w4,w5,w6,w7}`，唯一的跨块共享是 `w3`/`w4` 共享 `c`。并集基为 `{w0,w1,w3,w4,w5,w6,w7}`，全局秩 6，块基大小 `3+4=7`，秩高估 1。
- `pairs`：块 `{w0,w1},{w2,w3},{w4,w5},{w6,w7}`，跨块的共享为 `w1`/`w2` 共享 `b`（协方差 `2`）、`w1`/`w3` 共享 `b`（协方差 `1`）、`w3`/`w4` 共享 `c`（协方差 `1/4`）、`w5`/`w6` 共享 `e`（协方差 `1`）。块基大小 `2+2+2+2=8`，全局秩 6，秩高估 2。
- `striped`：块 `{w0,w4},{w1,w5},{w2,w6},{w3,w7}`，块间存在大量共享；块基大小和 8，全局秩 6，高估 2。
- `single`：8 个单窗口块，全部共享都在块间；高估 2。

---

## 5. 本包影响：逐计划核定

下表中“报告”指原入口 `pipeline.py` 的 `result`；“实际”指现有估计量（同一线性估计量）的真实方差 `Σ_k v_k L_k²`，`L_k` 为该计划下斜率的源负载系数，并由窗口值+完整窗口协方差二次型独立复核（两者精确相等）。“规范”为第 2 节 GLS。

| 计划 | 块基 | 报告秩 / 并集秩 | 截距 | 斜率 | 报告斜率方差 | 实际斜率方差 | 实际/报告 | 报告标准误 | 实际标准误 | 报告 95% 区间 | 实际方差对应区间 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `all` | `w0,w1,w3,w4,w5,w6` | 6 / 6 | `36073/36005` | `72028/36005` | `520/7201` | `520/7201` | 1 | `0.26872326425340776` | `0.26872326425340776` | `[1.47380233262852, 2.5271975285018784]` | `[1.47380233262852, 2.5271975285018784]` |
| `cut` | `w0,w1,w3` + `w4,w5,w6,w7` | 7 / 6 | `3569/3565` | `7132/3565` | `664/9269` | `6231448/85914361` | `778931/769327` | `0.26765021624604396` | `0.26931565891598497` | `[1.4759665859754256, 2.525155433659918]` | `[1.4727023183423413, 2.5284197012930023]` |
| `pairs` | `w0,w1` + `w2,w3` + `w4,w5` + `w6,w7` | 8 / 6 | `123539/123085` | `49247/24617` | `1630/24617` | `46692208/605996689` | `23346104/20062855` | `0.25732159777152963` | `0.2775793024778584` | `[1.4961777587118728, 2.504878421976269]` | `[1.4564726574874687, 2.5445835232006733]` |
| `striped` | `w0,w4` + `w1,w5` + `w2,w6` + `w3,w7` | 8 / 6 | `2857789/2845870` | `5694681/2845870` | `13930/284587` | `6719310430/80989760569` | `95990149/56632813` | `0.22124223110993743` | `0.28803634766676006` | `[1.567398654412281, 2.434668200363236]` | `[1.4364821859609087, 2.5655846688146084]` |
| `single` | 8 个单窗口 | 8 / 6 | `2857789/2845870` | `5694681/2845870` | `13930/284587` | `6719310430/80989760569` | `95990149/56632813` | `0.22124223110993743` | `0.28803634766676006` | `[1.567398654412281, 2.434668200363236]` | `[1.4364821859609087, 2.5655846688146084]` |

相对规范值的比例（斜率方差）：

| 计划 | 报告/规范 | 实际/规范 | 报告半宽/规范半宽 | 实际半宽/规范半宽 |
|---|---|---|---|---|
| `all` | 1 | 1 | 1 | 1 |
| `cut` | `0.99203` | `1.00441` | `0.99601` | `1.00220` |
| `pairs` | `0.91694` | `1.06700` | `0.95757` | `1.03296` |
| `striped` / `single` | `0.67784` | `1.14890` | `0.82331` | `1.07187` |

截距同样受影响（规范截距方差 `4932/7201`）：`cut` 报告 `5508/9269`（规范 86.8%）、实际 `63156852/85914361`（规范 107.3%）；`pairs` 报告 `14612/24617`（规范 86.7%）、实际 `451341844/605996689`（规范 108.7%）；`striped`/`single` 报告 `129270/284587`（规范 66.3%）、实际 `76419039970/80989760569`（规范 137.8%）。截距的偏差幅度更大，但事故入口只用斜率方差与秩作签名（`checks.py:45-47`），所以本包以斜率区间漂移暴露。

---

## 6. 现有估计量的线性系数与实际方差（不拿规范方差代替）

对每个计划，现有斜率估计量可写成源观测的线性式 `Σ_k L_k y_k`。精确系数如下（分母为各计划公共分母）：

| 计划 | `a` | `b` | `c` | `d` | `e` | `f` |
|---|---|---|---|---|---|---|
| `all` | `−1236/7201` | `−179/7201` | `−196/7201` | `36/7201` | `211/7201` | `1364/7201` |
| `cut` | `−1524/9269` | `−215/9269` | `−392/9269` | `4/713` | `283/9269` | `1796/9269` |
| `pairs` | `−3764/24617` | `−1067/24617` | `−1052/24617` | `130/24617` | `1367/24617` | `4386/24617` |
| `striped` / `single` | `−42812/284587` | `−19572/284587` | `−4219/284587` | `2941/284587` | `14536/284587` | `49126/284587` |

所有计划均满足 `Σ L_k = 0`、`Σ L_k t_k = 1`（线性无偏约束）。实际方差 `Σ_k v_k L_k²`（`v=(1,4,1,9,4,1)`）即上表“实际斜率方差”列，并与“把估计量写成窗口值线性式后用完整 8×8 窗口协方差求二次型”的独立结果精确相等（这是关键交叉核对：块间的相关必须进入这个二次型）。

要点：

- **`all`**：现有估计量 = 规范 BLUE，实际方差 = 报告方差 = `520/7201`。
- **`cut`**：实际方差 `6231448/85914361 = 0.0725309`，比报告 `664/9269 = 0.0716366` 大 `1.0124836` 倍；也比规范略大 `0.441%`——它已不是 BLUE，只是差得很小。
- **`pairs`**：实际 `0.0770503`，报告 `0.0662144`，实际/报告 `1.1636481`。
- **`striped`/`single`**：实际 `0.0829649`，报告 `0.0489481`，实际/报告 `1.6949564`，报告值把真实方差低估了 41.0%（`1 − 1/1.6949564`）。

因此必须区分三个量：(i) 报告方差 `(Σ_b J_b)⁻¹[1][1]`；(ii) 现有估计量的真实方差 `Σ v_k L_k²`；(iii) 规范 GLS 方差 `520/7201`。它们只有在 `all` 与独立新来源对照下才重合。

---

## 7. 点估计为何掩盖问题

`cut` 的情形可以精确说明机制。把 `cut` 的合并信息与规范信息相减，得到

```
ΔJ = J_cut − J_canonical = [[1, 2], [2, 4]] = x xᵀ,   x = (1, 2)
Δrhs = r_cut − r_canonical = [[5], [10]] = 5 x
```

恰等于“额外加入一条 `time=2, value=5, variance=1` 的独立观测”所贡献的 `x xᵀ/v` 与 `x·y/v`。源 `c` 正是 `time=2, value=5, variance=1`。也就是说，`cut` 把 `c` 的同一测量经 `w3`、`w4` 两个块各用了一次，效果等于**把源 `c` 复制成两次独立测量**。所以 `cut` 的合并解精确等于“规范数据 + 一份重复的 `c`”的 GLS 解（已核对完全相等），而不是 `cut` 数据本身的 GLS。

复制带来的估计偏移由 `J_cut⁻¹` 作用于与规范直线正交的残差决定：残差 `Δr − ΔJ·β_canonical = (−104/36005, −208/36005)`，于是斜率仅移动 `1568/25671565 ≈ 6.11e-5`。因为 `c` 的读数 `5` 几乎正好落在 `y=1+2t`（`t=2`）上，副本携带的“噪声”只有 `−104/36005 ≈ −0.00289`，点估计几乎不动。

一般地，本包所有窗口值几乎都落在 `y=1+2t`：只有 `w4`（`121/20` vs 6）、`w5`（`161/20` vs 8）各偏 `0.05`，源 `d` 偏 `0.1`。任一满足无偏约束的权重向量都会把斜率放在 `2` 附近，因此五种计划的斜率最大差只有 `5.335e-4`。**点估计接近不是正确性证据**；真正的诊断量是方差、区间与秩。源负载系数则清楚暴露了差异：`d` 的斜率系数从规范的 `36/7201 ≈ 0.00500` 变到 `striped/single` 的 `2941/284587 ≈ 0.01033`（约 2.06 倍），`d` 正是唯一明显偏离直线的源。

---

## 8. 区间为何随分块改变

报告标准误直接由 `inv(Σ_b J_b)[1][1]` 开方（`trend.py:5-12`），而 `Σ_b J_b` 依赖分块：

- `single` 的每个块只有一个窗口，`C_b` 退化为该窗口标量方差，`J_b = x_w x_wᵀ / var(w)`。
- `striped` 的块虽各有两个窗口，但每一对块内窗口（`w0/w4`、`w1/w5`、`w2/w6`、`w3/w7`）都不共享源，块内协方差精确为 0，`C_b` 仍是对角阵，`J_b` 仍是两个边际项之和。
- 因此这两种计划相加后

  ```
  J_single = J_striped = Σ_w x_w x_wᵀ / var(w) = [[1393/260, 3257/260],[3257/260, 12927/260]]
  ```

  这正是“把 8 个窗口都当成相互独立的测量”的信息矩阵，与独立新来源对照的规范 information 完全相同（已用 `reported_information_equals_marginal_sum = True` 逐计划核对）。所以 `single`/`striped` 报告方差 `13930/284587` 就是**虚构的全独立模型**方差。
- `pairs` 的块内仍是两窗口联合协方差求逆，保留了块内相关，但丢掉了块间相关，故报告方差介于全独立与规范之间（`1630/24617`）。
- `cut` 只丢了一条跨块共享（`c`），信息增量最小且退化为秩 1 副本，报告方差离规范最近（`664/9269`）。
- `all` 只有一块，块内协方差即整批协方差，故正确。

同时 `merge` 的 rank 相加（`summary.py:34`）把这些“虚假独立”块计入：报告秩 6/7/8/8/8 对规范 6，`decomposition.<plan>.rank_deficit` 精确记录 0/1/2/2/2，且等于 `sum(块基大小) − 并集秩`。

因此区间变窄并非精度提高，而是把共享来源与跨块相关当成了额外信息，随着分块越细（`single`/`striped` 最细）虚增越多，区间越窄（半宽降到规范的 `0.82331` 倍）。而现有估计量本身（因为丢掉了块间相关、权重不再最优）实际方差反而**大于**规范（`1.14890` 倍），其应得区间比规范更宽（`1.07187` 倍）；同时报告标准误只有该估计量自身真实标准误的 `0.768105` 倍（真实标准差是报告值的 `1.301905` 倍）。这三者必须分开陈述。

---

## 9. 三组单因素对照（只改一个因素）

派生数据全部在 `analysis.py` 内构造，不写回 `data/`。

### 9.1 原包仅换分块（`partition_only`）

只把 `plans` 从 `all` 换成 `cut`/`pairs`/`striped`/`single`，源、窗口、权重完全不变。报告斜率方差 `0.0722122 → 0.0716366 → 0.0662144 → 0.0489481`，报告秩 `6 → 7 → 8 → 8`，而规范解恒为 `(72028/36005, 520/7201, 秩 6)`。**分块是处理安排，不能改变同一批次的规范结果**；此处改变说明实现把分块当成了独立批次。这是“仅换分块”的直接证据，也排除了数据或公式变动。

### 9.2 追加与 `w3` 共用来源及权重的重复窗口（`duplicate_shared_source`）

新增 `w8 = {b:1/2, c:1/2}`，与 `w3` 完全相同（同源、同权、同时间、同读数）。

- 单块计划（`all9`）：仍为 `72028/36005`、`520/7201`、秩 6，与规范完全一致——因为一块内贪心选择会剔除依赖的 `w8`，块内协方差求逆正确。
- 全单窗口分块（`singles9`）：报告方差 `16010/334909 = 0.0478040`，比原 `single` 的 `0.0489481` **再降 2.34%**；报告秩 `9`（全局秩仍 6，高估 3）；斜率漂到 `6702177/3349090`。其“实际方差” `9688337670/112164038281 = 0.0863765`，实际/报告 `968833767/536189309 = 1.80689`。

对照结论：一条与既有测量完全重复的窗口**没有增加任何信息**（`w8 = w3`，并集零空间增加向量 `w3 − w8`，全局秩仍 6），但在最细分块下被算成一条新的独立测量，于是报告方差下降、秩上升、区间更窄。这直接证明“共享来源的窗口被当作新增独立测量”。

### 9.3 改用独立新来源（`independent_sources`）

把每个窗口引用的 `{源:权重}` 替换为 `{窗口::源:同权重}`，并保留原 `time`、`value`、`variance`；13 个新 id 互不共享。此时数值不变但误差独立：

- 规范 information 变为 `[[1393/260, 3257/260],[3257/260, 12927/260]]`，秩 8，规范斜率 `5694681/2845870`，方差 `13930/284587`，标准误 `0.22124223110993743`，区间 `[1.567398654412281, 2.434668200363236]`。
- 原有五种计划（`all/cut/pairs/striped/single`）全部得到同一结果：点估计、报告方差、秩都等于规范值，且每计划“报告方差 = 实际方差”（`all_plans_equal_canonical = True`、`all_plans_ratio_one = True`）。

对照结论：实现的分块合并逻辑在**真独立**数据上是自洽的、分块不变；问题只在跨块共享源。这也表明原包 `single` 的报告方差 `13930/284587` 恰好等于“若这些窗口真的独立”的规范方差——它是把共享源误当独立源的结果，不是原数据的规范方差。

### 9.4 三组对照的区分

| 对照 | 改变的单一因素 | 规范方差 | 报告方差 | 报告秩 | 实际/报告 |
|---|---|---|---|---|---|
| 原包 | — | `520/7201` | 随分块 `520/7201…13930/284587` | 6…8 | 1…`95990149/56632813` |
| 仅换分块 | 只改 `plan` | `520/7201`（不变） | `520/7201…13930/284587` | 6…8 | 1…`95990149/56632813` |
| 追加重复窗口 | 增加同源同权 `w8` | `520/7201`（不变） | 单块 `520/7201`；单窗口 `16010/334909` | 6；9 | 1；`968833767/536189309` |
| 独立新来源 | 只改源 id 独立性 | `13930/284587`（数据语义变了） | 五计划皆 `13930/284587` | 8 | 1 |

注意 9.3 的规范方差变化是因为数据语义确实变了（源不再共享），不是 bug；它用来隔离“共享”这一因素。9.1 严格保持数据不变，9.2 严格只增加一条无信息窗口。

---

## 10. 原事故入口与回归

- `logs/normal.log`：`checks.py normal` 通过，三计划均为 `intercept=1`、`slope=2`、`slope_variance=1/5`、`rank=4`（正常包 4 个窗口互不共享来源，不触发条件）。
- `logs/local.log`：`{"checked_block_count": 26, "local_covariance_checks": "passed"}`；局部检查只看块内对象，无法发现跨块相关丢失。
- `logs/incident-1.log`、`logs/incident-2.log`：两次运行的业务 JSON 完全相同（sha256 一致），末行均为 `AssertionError: TREND_PARTITION_DRIFT`，退出码非零。`analysis.py verify` 会再次连续两次运行原 `checks.py incident`，要求两次都非零退出且 stderr 含 `TREND_PARTITION_DRIFT`，并要求两次 stdout/stderr 哈希一致（确定性复现）。
- 差异签名由 `checks.py:45-47` 定义为 `(slope_variance, information_rank)`；本包中 `all` 与 `cut`/`pairs`/`striped`/`single` 在这两项上都不同，故必然失败。若只比较斜率（最大差 `5.3e-4`）则会漏报——这也是不能把“斜率接近”当作正确的原因。

---

## 11. 运行交接与复现

`.pointwise-runtime.json` 已按合同接通：

- `verify` → `["python3", "analysis.py", "verify"]`：从原始数据重算全部证据、逐条核对本报告数值锚点、运行原 `normal`/`local` 回归、连续两次运行原 `incident` 并确认非零退出且含 `TREND_PARTITION_DRIFT`，并把独立模型与原 `pipeline.py` 五种计划逐项比对（精确数相等、标准误/区间误差 ≤ `1e-10`）。任一不符即失败。
- `demo` → `["python3", "analysis.py", "demo"]`：展示规范结果、现有 `cut` 结果（报告区间 vs 实际方差对应区间）以及“追加重复窗口”这一真实单因素差异。

两种模式均不写入 `data/`、`logs/` 或任何原实现文件；派生包只在内存中构造。

---

## 12. 数值锚点（供 `analysis.py verify` 逐条核对）

以下每行 `key = value` 均对应 `evidence.json` 的 `anchors`，由精确有理数或报告浮点生成。

```
    canonical.intercept = 36073/36005
    canonical.slope = 72028/36005
    canonical.slope_variance = 520/7201
    canonical.standard_error = 0.26872326425340776
    canonical.ci95_low = 1.47380233262852
    canonical.ci95_high = 2.5271975285018784
    canonical.information_rank = 6
    canonical.full_covariance_rank = 6
    canonical.intercept_variance = 4932/7201
    canonical.basis = w0 w1 w3 w4 w5 w6
    plans.all.slope = 72028/36005
    plans.all.reported_slope_variance = 520/7201
    plans.all.actual_slope_variance = 520/7201
    plans.all.variance_ratio = 1/1
    plans.all.information_rank = 6
    plans.all.union_window_rank = 6
    plans.all.reported_intercept_variance = 4932/7201
    plans.all.actual_intercept_variance = 4932/7201
    plans.all.ci95_low = 1.47380233262852
    plans.all.ci95_high = 2.5271975285018784
    plans.all.actual_ci95_low = 1.47380233262852
    plans.all.actual_ci95_high = 2.5271975285018784
    plans.cut.slope = 7132/3565
    plans.cut.reported_slope_variance = 664/9269
    plans.cut.actual_slope_variance = 6231448/85914361
    plans.cut.variance_ratio = 778931/769327
    plans.cut.information_rank = 7
    plans.cut.union_window_rank = 6
    plans.cut.reported_intercept_variance = 5508/9269
    plans.cut.actual_intercept_variance = 63156852/85914361
    plans.cut.ci95_low = 1.4759665859754256
    plans.cut.ci95_high = 2.525155433659918
    plans.cut.actual_ci95_low = 1.4727023183423413
    plans.cut.actual_ci95_high = 2.5284197012930023
    plans.pairs.slope = 49247/24617
    plans.pairs.reported_slope_variance = 1630/24617
    plans.pairs.actual_slope_variance = 46692208/605996689
    plans.pairs.variance_ratio = 23346104/20062855
    plans.pairs.information_rank = 8
    plans.pairs.union_window_rank = 6
    plans.pairs.reported_intercept_variance = 14612/24617
    plans.pairs.actual_intercept_variance = 451341844/605996689
    plans.pairs.ci95_low = 1.4961777587118728
    plans.pairs.ci95_high = 2.504878421976269
    plans.pairs.actual_ci95_low = 1.4564726574874687
    plans.pairs.actual_ci95_high = 2.5445835232006733
    plans.striped.slope = 5694681/2845870
    plans.striped.reported_slope_variance = 13930/284587
    plans.striped.actual_slope_variance = 6719310430/80989760569
    plans.striped.variance_ratio = 95990149/56632813
    plans.striped.information_rank = 8
    plans.striped.union_window_rank = 6
    plans.striped.reported_intercept_variance = 129270/284587
    plans.striped.actual_intercept_variance = 76419039970/80989760569
    plans.striped.ci95_low = 1.567398654412281
    plans.striped.ci95_high = 2.434668200363236
    plans.striped.actual_ci95_low = 1.4364821859609087
    plans.striped.actual_ci95_high = 2.5655846688146084
    plans.single.slope = 5694681/2845870
    plans.single.reported_slope_variance = 13930/284587
    plans.single.actual_slope_variance = 6719310430/80989760569
    plans.single.variance_ratio = 95990149/56632813
    plans.single.information_rank = 8
    plans.single.union_window_rank = 6
    plans.single.reported_intercept_variance = 129270/284587
    plans.single.actual_intercept_variance = 76419039970/80989760569
    plans.single.ci95_low = 1.567398654412281
    plans.single.ci95_high = 2.434668200363236
    plans.single.actual_ci95_low = 1.4364821859609087
    plans.single.actual_ci95_high = 2.5655846688146084
    contrast.duplicate.canonical_slope = 72028/36005
    contrast.duplicate.canonical_slope_variance = 520/7201
    contrast.duplicate.canonical_rank = 6
    contrast.duplicate.all_slope = 72028/36005
    contrast.duplicate.all_reported_variance = 520/7201
    contrast.duplicate.all_actual_variance = 520/7201
    contrast.duplicate.all_rank = 6
    contrast.duplicate.singles_slope = 6702177/3349090
    contrast.duplicate.singles_reported_variance = 16010/334909
    contrast.duplicate.singles_actual_variance = 9688337670/112164038281
    contrast.duplicate.singles_variance_ratio = 968833767/536189309
    contrast.duplicate.singles_rank = 9
    contrast.independent.canonical_slope = 5694681/2845870
    contrast.independent.canonical_slope_variance = 13930/284587
    contrast.independent.canonical_rank = 8
    contrast.independent.all_plans_equal_canonical = True
    contrast.independent.all_plans_ratio_one = True
    decomposition.all.rank_deficit = 0
    decomposition.all.information_difference_rank = 0
    decomposition.all.information_difference_trace = 0/1
    decomposition.all.information_difference_determinant = 0/1
    decomposition.all.information_difference_psd = True
    decomposition.all.shift_matches = True
    decomposition.cut.rank_deficit = 1
    decomposition.cut.information_difference_rank = 1
    decomposition.cut.information_difference_trace = 5/1
    decomposition.cut.information_difference_determinant = 0/1
    decomposition.cut.information_difference_psd = True
    decomposition.cut.shift_matches = True
    decomposition.pairs.rank_deficit = 2
    decomposition.pairs.information_difference_rank = 2
    decomposition.pairs.information_difference_trace = 3151/882
    decomposition.pairs.information_difference_determinant = 121/196
    decomposition.pairs.information_difference_psd = True
    decomposition.pairs.shift_matches = True
    decomposition.striped.rank_deficit = 2
    decomposition.striped.information_difference_rank = 2
    decomposition.striped.information_difference_trace = 8057/468
    decomposition.striped.information_difference_determinant = 268313/23400
    decomposition.striped.information_difference_psd = True
    decomposition.striped.shift_matches = True
    decomposition.single.rank_deficit = 2
    decomposition.single.information_difference_rank = 2
    decomposition.single.information_difference_trace = 8057/468
    decomposition.single.information_difference_determinant = 268313/23400
    decomposition.single.information_difference_psd = True
    decomposition.single.shift_matches = True
    decomposition.cut.phantom_is_source_c_squared = True
```
