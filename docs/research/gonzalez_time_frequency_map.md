# González 研究谱系：时域方法迁移到频域的真实增量

核验日期：2026-10-12。范围限定为 Rodrigo A. González 与合作者的连续时间 IV 辨识、additive/structured 模型与柔性运动系统。这里的“已读”只表示读到相应方法、假设与定理，不表示独立复现其全部证明和实验。研究结论尚不是新颖性证明。

## 1. 结论先行

**不建议把“把 SRIVC/RIVC 从时域改为频域”作为主创新。** 该团队已有频域 additive MIMO RIV、运动平台验证以及 2026 年正式发表的频域两阶段结构化模态辨识。2025 年的另一篇原文还在明确条件下证明了非参数 FRF 估计对于参数辨识的充分性及有限样本等价关系。这些结果会直接反驳“跳过非参数辨识就必然更准确”和“频域版必然是新方法”。

建议只选一个应用领域：**冗余执行器柔性运动平台中的闭环模态辨识**。最值得继续验证的问题是：当正常任务控制只激励低维执行器子空间时，怎样以有限任务扰动取得足够的模态信息，并在执行器相对时延和邻近漏模存在时可靠估计四个目标模态？这与单纯增加一个时延参数不同。

该方向的关键对象应是“实际输入分配后的可辨识参数”和“参考信号/独立探针生成的工具变量”。时延承担执行器通道校准及任务扰动解释；第五模态仍然作为未建模动态，不应偷偷并入四模态模型然后声称原问题已解决。

## 2. 阅读等级和论文谱系

- A：已阅读原文方法、关键假设和相应定理/讨论；下文说明具体边界。
- A-M：已阅读正式论文核心方法及实验解释；没有该论文所声称性质的独立完整证明。
- B：元数据和摘要核验，用于定位前后关系，不从它抽取未读定理。
- 发表状态以正式出版信息优先；arXiv 年份不能冒充期刊发表年份。作者个人网站只辅助核对当前清单。

| 核心论文 | 年份、状态与原文 | 阅读等级/定位 | 对我们可迁移的东西 | 已占据的贡献，不能重称创新 |
|---|---|---|---|---|
| Statistical Analysis of Block Coordinate Descent Algorithms for Linear Continuous-Time System Identification | IEEE Control Systems Letters 8，2024；DOI [10.1109/LCSYS.2024.3389939](https://doi.org/10.1109/LCSYS.2024.3389939)；[原文](https://arxiv.org/html/2404.09071v1) | A；Lemma IV.1、Assumptions IV.1–IV.4、Theorems IV.2/IV.4 | 分块更新的耦合项、相关输入对收敛极限的影响；应与执行器分配造成的相关性连接 | BCD、SRIVC 分块求解、把其交叉相关写成谱积分 |
| Consistency analysis of refined instrumental variable methods for continuous-time system identification in closed-loop | Automatica 166，2024，111697；DOI [10.1016/j.automatica.2024.111697](https://doi.org/10.1016/j.automatica.2024.111697)；[正式 PDF](https://pure.tue.nl/ws/portalfiles/portal/330610186/1-s2.0-S0005109824001912-main.pdf) | A；连续/离散控制器区分，Theorems 11/15、Corollary 13 | 用外部参考重建无噪声输入的 CLSRIVC 工具变量；采样保持和闭环相关性条件 | “闭环仍用原 SRIVC 即一致”不成立；参考工具变量本身也已知 |
| Identification of additive continuous-time systems in open and closed loop | Automatica 173，2025，112013；DOI [10.1016/j.automatica.2024.112013](https://doi.org/10.1016/j.automatica.2024.112013)；[正式 PDF](https://pure.tue.nl/ws/portalfiles/portal/348706385/1-s2.0-S0005109824005077-main.pdf) | A；Assumptions 1–4、Theorem 2 及闭环 IV 构造 | 四个目标模态分别参数化；开闭环 instruments 设计和同阶真模型下的分析框架 | additive 参数化与相应 RIV 迭代；2024 DOI 不代表应写成 2024 卷期 |
| Identification of Additive Multivariable Continuous-Time Systems | IEEE Control Systems Letters 9，2025，547–552；DOI [10.1109/LCSYS.2025.3575179](https://doi.org/10.1109/LCSYS.2025.3575179)；[原文 v2](https://arxiv.org/pdf/2504.01639) | A；MIMO 开闭环模型与工具变量；已核对论文自述其为已有频域方法的时域对应物 | MIMO 参考输入 instruments、不同输出噪声协方差与模态分块 | “首次将 additive RIV 推广 MIMO”或其直接 FFT |
| Frequency domain identification for multivariable motion control systems: Applied to a prototype wafer stage | IFAC-PapersOnLine 59(17)，2025，67–72；DOI [10.1016/j.ifacol.2025.10.141](https://doi.org/10.1016/j.ifacol.2025.10.141)；[正式 PDF](https://pure.tue.nl/ws/files/373699105/1-s2.0-S2405896325015344-main.pdf) | A-M；FRF 输入、频域 additive RIV 与运动平台实验 | 直接比较用的频域 additive 基线 | 多变量运动系统的频域 additive RIV 已存在 |
| Statistically Optimal Structured Additive MIMO Continuous-time System Identification | arXiv 2505.14169，2025；[原文 v1](https://arxiv.org/html/2505.14169v1)；[DOI](https://doi.org/10.48550/arXiv.2505.14169)；截至本次核验仍列为投稿稿件 | A；Section 4、Theorems 4.1/4.2、Remark 3、Section 5/Theorem 5.3 | 先估计冗余参数再按协方差投影；投影映射的可辨识性、满秩条件 | “两阶段结构化 MIMO + 协方差加权”已占据；闭环最优需限定 IV 类 |
| Recursive Identification of Structured Systems: An Instrumental-Variable Approach Applied to Mechanical Systems | European Journal of Control 84，2025，101238；DOI [10.1016/j.ejcon.2025.101238](https://doi.org/10.1016/j.ejcon.2025.101238)；[原文](https://arxiv.org/html/2504.17927v1) | A；Sections 2–4、Theorems 1/2、Remark 4.5 与实验装置 | 若未来确有跟踪需求，可用递归分块和遗忘；本轮不列入主线 | 递归 structured IV、冗余执行器柔性梁实验均已有；正文基础模型是 SISO，不能据装置名称认为解决全 MIMO |
| Finite Sample MIMO System Identification with Multisine Excitation: Nonparametric, Direct, and Two-step Parametric Estimators | arXiv 2510.26929，2025-10-30；[原文 v1](https://arxiv.org/html/2510.26929v1)；[DOI](https://doi.org/10.48550/arXiv.2510.26929)；截至核验仍为投稿稿件 | A；Assumptions 1–3、Theorem 4、Corollary 2 | 直接方法与两步 FRF 方法的严格公平比较；完整协方差的重要性 | 在其条件下 FRF 是充分统计，直接原始数据不提供额外统计信息 |
| Structured identification of multivariable modal systems | Mechanical Systems and Signal Processing 247，2026，113948；DOI [10.1016/j.ymssp.2026.113948](https://doi.org/10.1016/j.ymssp.2026.113948)；[正式记录及全文入口](https://repository.tudelft.nl/record/uuid%3Aab88a6ca-77d0-4286-99d1-636b72e24b2a) | A-M；本地正式 PDF Sections 3–6，式 (21)/(22)/(38)/(39)，Remarks 1–3 | 当前最接近的 FRF→additive RIV→模态结构投影基线 | 频域两阶段、秩一留数、一般/比例阻尼、空间模态运动平台验证均已存在 |

2020 年理论基础不可忽略：[Consistency analysis...，Automatica 113，108767](https://doi.org/10.1016/j.automatica.2019.108767) 为基础一致性工作（本次 B）；[Efficiency Analysis...，Automatica 121，109196](https://doi.org/10.1016/j.automatica.2020.109196) 的[原文](https://arxiv.org/pdf/2002.00518)已读 Assumptions 1–6、Theorems 1/2（A）。后者效率结论要求真模型阶次正确、稳定/互质、输入充分激励、输入与噪声独立、输入采样间行为准确已知，并假定估计一致及迭代收敛；其 CRLB 达到结果针对相应输出误差及白高斯噪声条件。它不是失配四模态模型的效率保证。

截至核验，不能把 2025 structured 预印本和 2026 MSSP 当成同一篇或同一条定理。前者是统计结构化框架，后者是正式发表的频域模态方法；两者相互相关但假设、数据对象和论证范围不同。[作者当前研究清单](https://rodrigoagv.github.io/research.html)

## 3. 关键定理究竟保证了什么

### 3.1 收敛点、迭代收敛、局部极小和统计一致不是一回事

2024 BCD 论文在内层迭代收敛、矩阵非奇异等条件下建立驻点关系。其一致性分析还需要样本数及外层迭代的极限、充分激励和采样条件；有限次分块更新不是该极限。输入相关时，其他子模型的当前误差会进入本子模型的极限方程。把相同输入送入多个 additive 分支，不能沿用独立 MISO 输入情形下交叉项消失的结论。[原文 Lemma IV.1、Theorems IV.2/IV.4](https://arxiv.org/html/2404.09071v1)

2026 MSSP 式 (38) 的固定点满足式 (25) 的一阶条件，只直接给出驻点。文中把驻点接着称为局部最优，数学上还需局部 Hessian 等条件；Section 6.6 又明确指出迭代一般不保证代价单调下降。这是推理边界，不等于整套方法无效。我们已用独立确定性例子核对完整 RIV 步可使代价上升，见同目录 chapter3_audit.md。加阻尼或线搜索是合理实现改进，但单独不足以形成强论文创新。[2026 正式论文](https://doi.org/10.1016/j.ymssp.2026.113948)

递归论文题为“global convergence”的 Theorem 1 需要协调子问题的唯一最小值及下降等条件，其含义是收敛子序列满足相应驻点性质，不能读成“全局最优”。Theorem 2 的一致性极限要求遗忘因子为 1、相应迭代收敛并且系统/模型最终趋于 LTI；固定遗忘因子小于 1 时不会消除持续噪声造成的随机波动。它没有证明任意持续时变系统可无偏跟踪。[原文 Section 4](https://arxiv.org/html/2504.17927v1)

### 3.2 闭环一致性：正确工具变量和采样间行为都是条件

2024 闭环论文区分连续控制器与离散控制器加保持器。输入采样间行为错误时，连续时间模型的采样辨识可能不一致；原 SRIVC 不能自动克服闭环输入噪声相关性。Theorem 11 的直接 SRIVC 结果额外依赖噪声和闭环严格真有理等条件。Theorem 15 的 CLSRIVC 使用外部参考生成 instruments，但仍有激励、稳定性、采样、扰动范数小于非奇异性裕度以及迭代收敛等要求。[原文 Theorems 11/15、Corollary 13](https://doi.org/10.1016/j.automatica.2024.111697)

2025 additive 原文 Theorem 2 和 structured 稿 Theorem 4.1 均出现“插值相关扰动小于理想矩阵最小奇异值”的条件，而不是对所有采样率/激励/初始化无条件成立。structured 稿 Remark 3 明确其统计分析不覆盖临界稳定系统；不要因为算法允许积分环节就宣称相应一致性定理已经覆盖自由刚体模态。[additive 原文](https://doi.org/10.1016/j.automatica.2024.112013)、[structured 原文](https://arxiv.org/html/2505.14169v1)

### 3.3 “统计最优”的范围

structured 稿 Theorem 4.2 在白高斯输出噪声、正确模型和收敛等条件下给出渐近协方差。闭环中其式 (41) 的 Fisher 信息比所构造 IV 信息多一项；因此闭环中是相应 IV 类中的最小协方差，不能泛称达到全部估计器的 ML/CRLB 最优。

Section 5 的二阶段投影还要求结构映射光滑、可辨识以及 Jacobian 满列秩。模态左/右向量的任意缩放会破坏一一性，必须先固定尺度规范；MSSP 的 Remark 1 也处理了这个问题。低秩输入分配、重复模态、输入模态形状不可见时，不能直接套用一个满秩投影协方差。[structured 原文 Sections 4–5](https://arxiv.org/html/2505.14169v1)

## 4. 为什么直接 FFT 通常不是新估计器

以下为独立推导，不能当成上述作者对全部情形的结论。

设时间域残差向量为 $e$，已知正定协方差为 $K$，$F$ 为单位酉 DFT。令 $E=Fe$，则

$$
e^H K^{-1}e
=E^H(FKF^H)^{-1}E.
$$

因此，对相同残差、全部频点、正确变换的完整协方差和边界处理，时域 GLS 与频域 GLS 具有相同目标函数。实数数据需正确处理正负频率配对；只用单边谱时相应常数及权重也要一致。把滤波后的数据 FFT、将和写为谱和，并不改变可辨识性或闭环噪声相关性。

有限记录下，卷积边界/初始状态不能凭 FFT 消失；任意非周期采样数据也不能无条件满足 $Y_k=G(j\omega_k)U_k+V_k$。若由离散输入序列经 ZOH 驱动连续对象，通常应使用包含采样保持和时延的离散等效映射 $G_d(e^{j\omega T_s};\rho,\tau)$，或准确处理连续输入及其谱。把 $G(j\omega)$ 直接乘采样 DFT 可能引入另一种模型失配。

所谓“可迁移”真正有价值的部分可能是：

1. 把有效闭环 instruments 与记录边界、频线选择和噪声协方差统一处理，并证明选频后仍满足所需矩条件。
2. 针对实际输入分配，分析目标模态 Jacobian 的秩，而不是假定全部执行器独立充分激励。
3. 把第五模态残差对目标参数的投影偏差写清，并设计激励降低它与相对时延方向的相关性。

这些仍需证明、实验和相关领域查重；“频域表达式更简洁”本身不是统计性能的提升。

## 5. 原始 DFT、FRF 拟合和 LRM 不能混成一个概念

### 5.1 两类输入数据

原始频域拟合保留 $U_k,Y_k$，必要时还保留独立参考 $R_k$、探针 $\Xi_k$ 与边界项。其残差是输出量纲，闭环中 $U_k$ 与输出噪声相关，需要合理的 IV 或噪声模型。

FRF 拟合保留 $\widehat G_k$ 及其协方差，以 $\widehat G_k-G_k(\rho)$ 为残差。输出/输入比值并不天然无偏，闭环 FRF 估计必须已有相应去偏机制。LRM 是一种非参数估计/瞬态处理选择，不等同于所有 FRF 方法，也不是原始数据不可辨识时的补信息工具。

### 5.2 已有精确等价结果，不能笼统声称原始数据优于 FRF

2025 多正弦有限样本论文要求稳定因果连续系统、已知连续多正弦输入、每个频率的多次实验输入矩阵满秩、采样频率之间无混叠重合、稳态记录与相应白噪声条件。Gaussian 条件用于充分统计结果；等价加权使用正确的完整协方差。Theorem 4 给出非参数 LS FRF 对参数辨识的充分性，Corollary 2 给出与时间域参数拟合相同的最小值解。

**该主结果不要求记录恰为整数周期。** 额外无泄漏/整数周期条件才使频点之间的表示和权重进一步分离。把“非整数周期”单独当成否定充分性的理由不正确。同样，这篇论文并没有使用“必须先做 LRM”的前提。[原文 Assumptions 1–3、Theorem 4、Corollary 2](https://arxiv.org/html/2510.26929v1)

我们能检验的真实差异应是：瞬态、输入分配秩不足、闭环输入噪声、频点协方差被忽略、未知相对时延与漏模同时出现时，哪些等价条件失效。不能先宣布原始 DFT 一定更好，再只挑 FRF 权重错误的基线证明它。

## 6. 与冗余执行器背景连接的独立推导

### 6.1 硬件输入多，不代表辨识信息多

设物理执行器数为 $m$，任务指令维数为 $r<m$，固定分配矩阵为 $T\in\mathbb R^{m\times r}$。正常任务输入

$$
u=T v,\qquad y=G u+e=(G T)v+e.
$$

只要 $\Delta G T=0$，$G$ 与 $G+\Delta G$ 对这些实验不可区分，且

$$
S_{uu}=T S_{vv}T^H,\qquad \operatorname{rank}(S_{uu})\le r.
$$

对模态展开

$$
G(s)=\sum_{i=1}^{4}\frac{\ell_i r_i^T}{d_i(s)}+G_{\rm miss}(s),
$$

任务数据只直接约束 $r_i^T T$，不自动约束完整执行器侧模态形状 $r_i$。**这不等于所有模态频率都不可辨识**：只要可见模态有足够输入/输出响应，部分频率和阻尼仍可辨识。严谨条件应写成所关心参数在 $\operatorname{vec}(G T)$ 中的 Jacobian 满列秩，消除尺度规范后再判断。

这也说明为什么 2025 多正弦定理不能原样用于固定低秩分配：它的逐频输入矩阵满秩假设未满足。但可针对被激励子空间建立降维充分统计，不能由此笼统宣称“FRF 路径失效”。

### 6.2 两条同一领域内的补信息方式

跨实验改变已知分配 $T_e$，得到 $G T_e$；若 $[T_1\ \cdots\ T_E]$ 满行秩且对象在实验间保持相同，就有机会恢复完整空间映射。若仅估计部分目标模态参数，可用更弱的参数 Jacobian 条件。

或者，在已有任务输入上加入外部探针：

$$
u=T v+N\xi,\qquad R_{\rm rb}N=0,
$$

其中 $R_{\rm rb}$ 是已知刚体力/力矩映射，$N$ 是其零空间基。探针可激励柔性动态而限制理想静态刚体任务扰动。这里的“限制”不等于实际闭环任务输出严格为零；柔性耦合、时延和建模误差仍需计入约束。

单纯“零空间激励”不应被直接申报新颖性，仍需和该领域已有动态控制分配、实验设计及执行器校准论文逐项对比。

### 6.3 闭环 instruments 应随分配结构构造

若控制律为

$$
u=T K(r-y)+N\xi,\qquad y=G_d u+e,
$$

则无噪声输入分量为

$$
u_0=(I_m+T K G_d)^{-1}[\,T K\quad N\,]
\begin{bmatrix}r\\\xi\end{bmatrix}.
$$

可以用当前模型替代 $G_d$，由与测量噪声独立的外部参考/探针构造 instruments。矩阵尺寸和稳定闭环逆必须正确；实际平台若有不同注入位置，表达式随之改变。这是现有 CLSRIVC 思路的结构化使用，真正待证明的是有限频域记录、分配和选频后工具变量的有效性及目标参数信息条件，不是外部参考的概念本身。

### 6.4 相对时延为什么与冗余执行器问题有关

执行器独立时延矩阵为

$$
D_\tau(s)=\operatorname{diag}(e^{-s\tau_1},\ldots,e^{-s\tau_m}).
$$

理想静态零空间探针满足 $R_{\rm rb}N=0$，但

$$
R_{\rm rb}D_\tau(s)N
\approx-s R_{\rm rb}\operatorname{diag}(\tau)N,
$$

一般不再为零。因此**通道相对时延**会把原本的刚体零空间探针泄漏为任务扰动，也会改变柔性模态输入形状的相位。若所有通道时延相同，则 $D_\tau=e^{-s\tau}I$，上述静态零空间关系仍保持；不能把公共时延和相对时延混为一谈。

这提供了把“时延作为其中一项”纳入主线的物理理由。它没有自动证明所有通道时延可分别估计；需要固定参考时延规范，检查延迟参数与模态、输入分配的局部 Jacobian 是否独立，并排除带宽不足导致的近似不可辨识。

### 6.5 固定四目标模态时，漏模偏差必须保留

若第五模态被有意留在模型外，即使无限样本和零测量噪声，估计通常也只趋于所选代价的伪真值。在线性化的加权实数残差坐标下，令 $A$ 为目标模态及其他允许参数的 Jacobian，$b$ 为一项时延的 Jacobian，$d$ 为第五模态等未建模残差，$Q=I-AA^\dagger$，则

$$
\delta\tau\approx\frac{b^TQd}{b^TQb}.
$$

这里已把权重的平方根吸收到 $A,b,d$。分母衡量排除其他参数后剩余的时延信息；分子衡量漏模残差伪装成时延的程度。增加信噪比不一定消除这个偏差。激励设计和干扰项建模的目标可同时约束这两个量，但要给出可检验的残差结构/上界，不能允许任意干扰函数后仍声称唯一辨识。

## 7. 一个收敛的研究题目与可证伪实验

建议工作题目：

**任务扰动约束下冗余执行器柔性平台的闭环模态辨识：输入分配、相对时延与漏模偏差。**

最小可交付的理论问题：

1. 给出在已有任务分配及附加探针下，四个目标模态参数和相对时延的局部可辨识条件；区分“完整空间模态形状”和“任务可见形状”。
2. 构造以外部参考/探针为来源的频域 IV，明确采样保持、瞬态和协方差假设。在真模型情形证明相应一致性；在保留第五模态的失配情形分析伪真值偏移，不冒称一致。
3. 在输入能量及实际任务响应约束下设计探针，降低目标模态与相对时延/漏模之间的混淆。该设计若只是标准 Fisher 信息优化，必须进一步说明和既有实验设计的差异。

需要与以下基线公平比较：

| 比较 | 必须保持一致 | 可反驳我们的结果 |
|---|---|---|
| 2025 时域 additive MIMO IV vs 同数据频域 IV | 输入、样本、模型类、边界、协方差和初始化 | 正确变换后完全等价，则不能宣称频域具有统计增益 |
| 2026 FRF 两阶段模态方法 vs 原始频域方法 | 有效数据、正确 FRF 协方差、相同目标模态数 | 等价条件成立时效果相同是合理结果；只能在明确失效条件下讨论差异 |
| 固定任务分配 vs 零空间/多分配实验 | 总输入能量、闭环任务误差约束、测量时长 | 激励增益仅来自额外能量，而非可辨识方向改善 |
| 已知时延/未知公共时延/未知通道相对时延 | 同一数据生成系统和优化容差 | 延迟信息近奇异或估计改善仅由多参数过拟合 |
| 四模态失配模型 vs 五模态 oracle | 四模态任务的评价只对四个真实目标模态 | 延迟估计吸收第五模态，训练拟合变好但验证参数更差 |

至少报告目标频率/阻尼/空间形状误差、延迟误差、任务输出扰动、输入能量、参数信息矩阵的最小奇异值及独立实验验证。仅报告 FRF 拟合误差不足以支撑“辨识成功”。

如果现有数据只有固定低秩分配且没有独立探针，首先报告可辨识的降维参数。算法优化不能恢复从未被激励、也未被已知结构唯一约束的空间自由度。实验权限或额外结构先验缺失时，这就是方向的边界。

## 8. 明确排除的路线

- 不以“去掉 LRM”直接宣称创新；它是数据处理选择，需要与正确非参数协方差基线比较。
- 不以“additive + MIMO + 频域 + 二阶段秩一投影”宣称新组合；2025–2026 谱系已覆盖主体。
- 不把闭环 IV 最小协方差写成无条件 ML 最优。
- 不把递归论文扩展为另一条在线跟踪主线；它会引入遗忘、时变与相关频窗的新问题，冲淡当前领域聚焦。
- 不把刚体临界稳定算法可运行当成统计定理已成立。
- 不把四模态失配的拟合最优称为四个物理模态的无偏估计。

本报告确认的是“哪些直接迁移已被占据，以及哪类问题值得进一步验证”。对所建议的冗余分配/探针/通道校准联合方向，仍须用该领域原始论文完成最终新颖性对照。
