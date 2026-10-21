# 冗余执行器与传感器辨识：近邻先例及主线收敛

检索/核验日期：2026-10-12。对象：4 个执行器、4 个传感器、3 个刚体自由度，目标估计前 4 个柔性模态，第 5 个保持未建模。本文仅做文献核验与研究判断，未修改或运行 MATLAB。

## 结论

**否决把“零空间激励 + RIV + 时延”直接当作创新；有条件保留一条主线：刚体动态模型失配时，利用输入、输出两侧刚体零空间形成交叉观测，并进行共享极点的结构 IV 柔性模态辨识。**

现有文献已经覆盖：冗余执行器与传感器、刚体解耦、零空间激励、约束最优试验设计、双侧 FRF 坐标变换、刚体先验约束下的频域 IV、加性连续时间闭环 IV、模态结构投影。单行/单列模态恢复及近节点问题也已有研究。因此候选问题必须落在“**刚体动态允许未知或失配，而所需柔性模态仍可恢复，且给出噪声和几何误差下的适用边界**”上。

本次定向检索尚未核实到完全相同的“左右刚体零空间交叉观测 + 完整非对称模态残量恢复 + 结构 IV”论文。这只是检索结果，不是新颖性证明。特别应与下列 S1、S2、S3、S5、S6 正面对比；简单拼接它们未必形成可发表增量。

## 九篇核心先例

阅读等级：A = 已读取作者稿/公开全文的相关推导段；B = 出版社或大学库摘要和元数据已核验；不将 B 等同于已审阅全部理论。

| 编号 | 已核验文献与原始来源 | 阅读等级与明确先例 | 本项目仍须区别之处 |
|---|---|---|---|
| S1 | Voorhoeve, de Rozario & Oomen (2016), *Identification for Motion Control: Incorporating Constraints and Numerical Considerations*. ACC, 6209–6214. DOI [10.1109/ACC.2016.7526645](https://doi.org/10.1109/ACC.2016.7526645)；[作者公开全文](https://gitlab.tue.nl/oomenpublic/pdfs/-/raw/main/pdf/VoorhoeveRozOom2016.pdf)；[TU/e 元数据](https://research.tue.nl/en/publications/identification-for-motion-control-incorporating-constraints-and-n) | **A，最直接的机械辨识近例。** §III-B 式(20)–(21)对输入、输出做刚体解耦并同步转换权重；§IV涵盖 SK/IV。实验恰为 4×4、3 个刚体自由度，第四压电输出仅测内部变形；§V采用闭环多正弦 FRF、时延预补偿。 | 已经做过“先解耦再辨识”，不能据此宣称创新。其刚体先验通过积分器结构嵌入全模型；尚未见以任意 (H_r(s)) 为 nuisance、只使用交叉子集并恢复完整柔性残量的构造。 |
| S2 | Häußler & Rixen (2017), *Optimal Transformation of Frequency Response Functions on Interface Deformation Modes*. IMAC, 225–237. DOI [10.1007/978-3-319-54930-9_20](https://doi.org/10.1007/978-3-319-54930-9_20)；[TUM 作者稿](https://mediatum.ub.tum.de/doc/1355455/397849.pdf)；[TUM 元数据](https://portal.fis.tum.de/en/publications/optimal-transformation-of-frequency-response-functions-on-interfa/) | **B+：大学摘要及作者稿检索可见段；全文读取接口未成功。** FRF 双侧变换、位移/力几何映射、加权最小二乘解释均已有。摘要明确允许按测量方差和机械能选择变换权重。 | 不是新提出投影或加权投影；其任务是接口变形模式与动态子结构耦合，不等同于本候选的模态参数 IV 和交叉残量恢复。提交前应逐式审阅完整作者稿。 |
| S3 | Pogačar, Bregar, Čepon & Boltežar (2021；2020 在线), *Near-to-node modal identification using multiple related response models*. Measurement 172,108793. DOI [10.1016/j.measurement.2020.108793](https://doi.org/10.1016/j.measurement.2020.108793)；[实验室作者稿](https://www.ladisk.si/data/pdf/Near-to-node%20modal%20identification%20using%20multiple%20related%20response%20models.pdf)；[实验室记录](https://www.ladisk.si/?ID=265&what=abstract) | **A，近节点与相关响应的直接先例。** 讨论单行/单列响应模型、非互补激励-响应的非对称模态常数；组合不同传感器响应并改造 LSFD，在近节点时提高弱模态常数的一致性。 | “少量行列就可恢复模态”和“弱参考点要融合其他响应”均非新。候选需说明几何零空间交叉观测与这些参考模型之间的区别，量化 pivot 很小时的病态性。 |
| S4 | González, Classens, Rojas, Welsh & Oomen (2025), *Identification of additive continuous-time systems in open and closed loop*. Automatica 173,112013. DOI [10.1016/j.automatica.2024.112013](https://doi.org/10.1016/j.automatica.2024.112013)；[出版社](https://www.sciencedirect.com/science/article/pii/S0005109824005077)；[作者预印本](https://arxiv.org/abs/2401.01263) | **A/B：已读公开摘要及相关方法段。** 加性连续时间模型的开闭环 RIV、一致性条件、边界稳定机械对象已有。 | 是算法基线，不能把加性闭环 RIV 本身列为贡献。模型截断或错误刚体结构不自动满足其正确模型的一致性前提。 |
| S5 | González, van der Hulst, Classens & Oomen (2025), *Statistically Optimal Structured Additive MIMO Continuous-time System Identification*. [arXiv:2505.14169](https://arxiv.org/abs/2505.14169)；[公开全文 v1](https://arxiv.org/html/2505.14169v1) | **A。** 结构化加性 MIMO 的两阶段 IV/投影、开闭环处理、参数协方差已有；开环统计效率与闭环 IV 类中的最优性应区分。闭环一般不达到所有估计器的 CRLB。 | 不能把“结构 IV + 模态投影 + 协方差”重复贡献。候选应改变刚体 nuisance 假设/可用观测并重新证明可辨识性及估计方程有效性。 |
| S6 | van der Hulst, González, Classens, Tacx, Dirkx, van de Wijdeven & Oomen (2026), *Structured identification of multivariable modal systems*. MSSP 247,113948. DOI [10.1016/j.ymssp.2026.113948](https://doi.org/10.1016/j.ymssp.2026.113948)；[作者预印本](https://arxiv.org/abs/2510.10820)；[全文](https://arxiv.org/html/2510.10820v1) | **A，用户论文所在直接技术线。** 频域加性 RIV 后投影到模态形式，处理一般/比例阻尼；13 执行器、4 传感器原型 wafer stage，3 刚体与17柔性。 | 冗余平台、完整模态残量、结构投影不是空白。可对比其全 FRF 联合拟合在刚体模型失配时的偏差，而不能只对比不加权、初始化较差的实现。 |
| S7 | Tacx, van de Vosse, Voorhoeve, Witvoet, Heertjes & Oomen (2025；2024 在线), *Spatio-temporal modeling for overactuated motion control*. Mechatronics 105,103270. DOI [10.1016/j.mechatronics.2024.103270](https://doi.org/10.1016/j.mechatronics.2024.103270)；[TNO 全文](https://publications.tno.nl/publication/34643448/VeZLCCzd/tacx-2025-spatio-temporal.pdf) | **A。** 冗余运动系统的时空模态建模、执行器/传感器空间信息融合及模型更新。§4.1式(12)明确要求 FRF 试验激励跨越完整输入空间。 | 只用一维输入零空间激励不能恢复全输入方向。这限制不会因时间信号足够丰富而消失；交叉观测的行部分仍需合适的其他输入实验。 |
| S8 | Dirkx, van de Wijdeven & Oomen (2020), *Frequency Response Function identification for multivariable motion control: Optimal experiment design with element-wise constraints*. Mechatronics 71,102440. DOI [10.1016/j.mechatronics.2020.102440](https://doi.org/10.1016/j.mechatronics.2020.102440)；[出版社](https://www.sciencedirect.com/science/article/pii/S0957415820301100)；[TU/e 全文](https://pure.tue.nl/ws/portalfiles/portal/173830510/1_s2.0_S0957415820301100_main.pdf) | **A。** 逐输入/输出功率约束下的频率相关多变量激励设计；7×8 冗余 wafer stage 验证；与 SIMO 和正交激励比较。 | 否决“用冗余空间改善 FRF，且限制刚体/输出运动”作为独立新主线。本项目本轮不再把 OED 扩成第二主线。 |
| S9 | Weiss, Leve, Kolmanovsky & Jah (2015), *Reaction Wheel Parameter Identification and Control through Receding Horizon-Based Null Motion Excitation*. *Advances in Estimation, Navigation, and Spacecraft Control*,477–494. DOI [10.1007/978-3-662-44785-7_25](https://doi.org/10.1007/978-3-662-44785-7_25)；[出版社](https://link.springer.com/chapter/10.1007/978-3-662-44785-7_25) | **B。** 滚动时域零空间激励，以参数误差协方差与激励代价为目标，在较小主体扰动下辨识反作用轮参数。 | 直接否决“借冗余零空间做辨识而不大幅扰动主运动”的概念新颖性；对象并非柔性模态，故不直接否定本文更窄的刚体失配问题。 |

另一个必须记录的公开技术先例：Kwan 与 Laro，*Method for controlling a motion of optical elements in lithography systems*，公开专利 [US20140204354A1](https://patents.google.com/patent/US20140204354A1/en)，公开日 2014-07-24，优先权日 2011-10-07。已核验全文 Eq.(9)–(10)：7 个传感器/6 个刚体自由度，由刚体传感器映射的零空间构造不含刚体信息的柔性输出，用于控制补偿。它强烈否定“输出零空间消刚体”本身的新颖性；不据此作专利权法律判断。

## 单一主线的精确定义

设

\[
G(s)=C_rH_r(s)B_r+\sum_{i=1}^{4}\frac{R_i}{a_i(s)}+G_5(s),\quad
B_r\in\mathbb R^{3\times4},\ C_r\in\mathbb R^{4\times3},
\]

其中 (B_r,C_r) 满秩，几何已知；(H_r(s)) 允许未知或建模错误，不局限于理想 (M^{-1}/s^2)。(G_5) 表示用户有意不拟合的第 5 阶及其余残差；**下述消除并不自动消除它**。选择单位向量

\[
B_rn=0,\qquad m^TC_r=0.
\]

两组观测为 (Gn) 和 (m^TG)。前者有 4 个元素，后者有 4 个元素，但共同包含 (m^TGn)，仅 7 个独立组合。将几何基补成两个正交基后，这就是变换后 4×4 FRF 的末列与末行，去掉左上 3×3 刚体块；不是仅保留标量 (m^TGn)。全部交叉观测均准确消去 (C_rH_rB_r)。

对一个满足 (R_i=l_ir_i^T) 的模态，记

\[
a_i^{\rm cross}=R_in,\quad b_i^{\rm cross}=m^TR_i,\quad c_i=m^TR_in.
\]

当 (c_i\ne0) 时有已知的 rank-one cross/skeleton 恒等式

\[
R_i=\frac{a_i^{\rm cross}b_i^{\rm cross}}{c_i}.
\]

**这条恒等式不是创新，也不能逐频率直接套在总 FRF 上**：多模态总和通常不是秩一。必须先用共同极点把各模态分开。比例阻尼的实二阶常数残量适用；一般阻尼时应按简单复极点残量分析，不能未经证明把整个二阶分子当作同一个秩一常矩阵。

本项目现有 `matlab/anlaysis/anlysis.m` 已包含 `Pbnull = null(Br)` 和 `Pcnull = null(Cr')'`。它们提供几何基础，不意味着当前数据已经满足交叉观测的可辨识与噪声要求。

### 可以研究的增量与不能研究成贡献的部分

- 可检验的研究问题：当刚体动态或刚体模型阶次错了，但刚体输入/输出子空间仍准确时，交叉观测能否保留目标柔性残量的可辨识性，并阻止刚体模型误差向柔性极点与阻尼传递？
- 潜在理论内容：交叉观测可恢复的充要条件、噪声加权的联合 IV 方程、几何子空间误差引起的泄漏界、弱交叉响应的拒绝/回退准则。新颖性仍待进一步对照已有辨识与 nuisance projection 理论。
- 已知基础：正交投影消去 nuisance、加权残差消元/Schur 补、rank-one 完成、最小二乘/IV 的坐标变换。把已知步骤串联，并不足以证明方法新。
- 若每个频率的 (H_r) 完全自由，白噪声且采用正交坐标时，对全 FRF 消去刚体 nuisance 的 profile least-squares 确实只剩 cross；一般协方差下应在白化空间消元。这个等价性应当用作基线/解释，不能自行把它称为新的统计原理。
- 7 个独立交叉残量恰对应非对称 4×4 秩一残量的 7 维自由度。它不额外提供秩约束降维带来的方差收益；收益假说是对刚体失配更稳健，而非无条件更精确。

## 必须先说清的失败边界

1. **小 pivot：** (c_i=(m^Tl_i)(r_i^Tn)) 小时残量恢复病态；任一因子为零便可能丢失该侧方向或整模态。相同现象与 S3 的近节点问题相连。4对3的冗余只有一维，不存在自由旋转 null 向量来修复。
2. **仅 input-null 不够：** 只施加 (u=nq) 得到 (Gn)，不能得到任意输入方向的 (m^TG)。行部分通常需要其他输入方向的外生激励；“7通道”并不等于 7 个独立物理传感器，也不自动节省试验。
3. **输出反馈不等于无刚体运动：** 名义 null 力不直接驱动刚体，但柔性传感器反馈、几何误差和通道动态可能在实际闭环中引起刚体运动。应在被辨识的实际输入/输出关系上推导。
4. **几何错误：** (B_r\) 或 (C_r\) 误差使 cross 不再精确消刚体；强刚体响应乘小角度误差也可显著。必须独立扫几何误差，不只测传感器白噪声。
5. **重复观测和协方差：** 直接堆叠4+4含一次重复的 pivot，协方差奇异，不能按8个独立误差算权重。应选正交 cross 基或显式处理冗余及相关性。
6. **第5阶仍在：** (G_5n) 与 (m^TG_5\) 一般非零。主方法仅免疫刚体子空间内失配，不能承诺消除任意遗漏柔性模态偏差。第四/第五模态靠近时必须单列结果。
7. **时延：** 共同标量时延保持几何消除，可作为同一模型中的次要参数；不同执行器/传感器时延一般破坏静态几何 null。不能从共同时延结论外推到所有通道时延。
8. **IV 有效性：** 闭环需使用独立于噪声的参考/外生激励及相应滤波构造工具变量；投影不自动赋予外生性。S4/S5 的正确模型一致性不能直接移植到仍遗漏第五模态的情况。

## 最小判别实验：只验证这一主线

在相同激励、数据长度、频率点和合理权重下，至少比较：

| 方法 | 用途 |
|---|---|
| 全16通道结构 RIV，正确刚体模型 | 信息充分的 oracle 参考；预计交叉法无理由普遍超过它。 |
| 全16通道结构 RIV，刚体模型有受控失配 | 检验目标问题是否实际存在，避免制造不合理弱基线。 |
| 全16通道中按频点消去任意刚体 nuisance，再估计柔性 | 强等价基线；若与 cross 一致是理论自检，不是额外性能胜利。 |
| 单侧 (Gn)、单侧 (m^TG)，以及两侧 cross 的联合结构 IV | 证明双侧恢复完整残量的必要性；单侧不能恢复全非对称残量须如实计分。 |

首先关掉第5模态，仅注入刚体动态失配，检验消除机制；随后恢复第5未建模模态、加入共同时延，作为同一机制的压力测试。最后独立扫 (c_i\) 大小和几何误差。报告第4模态频率/阻尼、完整 (R_i\) 的尺度一致误差、重构 FRF误差、失败比例与区间覆盖，不能仅报整体拟合百分比。

**继续条件：** 能清楚证明全模型刚体失配会污染目标柔性模态，而 cross 在理论允许区域内保持稳定，并能事先诊断弱 pivot / 几何泄漏而拒绝错误完成。

**否决条件：** 原项目的刚体模型已经足够准确、cross 各目标模态接近节点、或者与 profile nuisance 基线等价后只剩坐标改写且没有新的可辨识/稳健性结果。此时应把它作为工程预处理改进，不硬写成论文创新。

## 检索范围与剩余核验

本轮以作者/出版社/大学库为证据，检索了 overactuated/oversensed modal identification、null-motion excitation、rigid-body-free/bilateral/nullspace projection、FRF row/column reconstruction、near-node modal identification、González additive/structured IV 等词组。未用无搜索命中证明无人做。最需要进一步全文排查的是 S2 的双侧加权变换、S1 的刚体约束 IV 和 S3 的非对称相关响应恢复；专利是公开技术邻近证据，不作权利范围判断。

最终选择只有一个：**刚体模型不确定下，基于两侧刚体零空间交叉观测的柔性模态结构 IV。** 零空间激励最优设计不再独立展开；共同时延与第5阶未建模仅作为模型项或适用边界。
