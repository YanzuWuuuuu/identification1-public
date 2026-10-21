# 闭环参考矩、目标一致性与隔离性能证书：数学审查

日期：2026-10-17。范围：只推导，不运行新辨识或实验，不修改原论文。

## 0. 审查结论

**可以形成严格的条件闭环，但不是无条件证明现有数据已经足够。** 最稳妥的链条是

\[
\text{正确边界与外部参考}\Longrightarrow
\text{零均值矩条件}\Longrightarrow
\text{目标分离/等价类一致性}\Longrightarrow
\text{矩检验置信集}\Longrightarrow
\text{固定校准器下的最坏隔离性能界}.
\]

这里“目标分离”“观测模型正确”“频谱边界误差可忽略”是需要验证的条件，不能由优化收敛替代。弱辨识时应优先采用矩检验反演置信集；参数 Hessian 的伪逆不是自动有效的置信区间。校准器与置信集可以来自同一训练集，但必须对同一个**冻结校准器**遍历候选真实对象，不能每个候选对象重新校准后宣称泄漏都为零。

## 1. 数据、观测边界和重复实验

在输入外加信号的正确边界设

\[
y=G_0u+v,\qquad u=d-Ky.
\tag{1}
\]

\(d\) 与总加性扰动 \(v\) 外生，即其有关频谱满足 \(E[V_kD_k^H]=0\)。\(K\) 包含实际数字控制器、几何映射、增益、反馈滤波及延迟。总电机指令 \(u\) 通常与 \(v\) 相关，但这不破坏外部参考矩。

以下将 \(G_k(\theta)\) 定义为实际采样数据对应的频率算子。若使用连续模型，保持器、采样、抗混叠及传感链必须已正确纳入；不能在存在显著混叠/保持器差异时直接令其等于裸机械 \(G(j\omega_k)\)。

取预先确定的有限频率集合 \(\Omega=\{\omega_1,\ldots,\omega_K\}\)。重复索引 \(r=1,\ldots,R\) 表示独立的完整实验簇，每簇可以包含四组独立输入方向及多个周期，而不是把共享窗口的相邻频点算作独立重复。定义实数化算子 \(\mathcal R(z)=[\Re z^T,\Im z^T]^T\)，以及

\[
q_r(\theta)=\mathcal R\operatorname{col}_{k=1}^K
\operatorname{vec}\big[(Y_{rk}-G_k(\theta)U_{rk})D_{rk}^H\big],
\quad \bar q_R(\theta)=R^{-1}\sum_rq_r(\theta).
\tag{2}
\]

允许对式(2)左乘预先给定的实矩阵 \(F\) 压缩矩维数；下面的 \(q\) 均指实际保留的矩维数。压缩可能丢失辨识信息，必须重新核对目标分离。

有限窗口一般存在边界/瞬态项 \(Y_{rk}=G_{0k}U_{rk}+V_{rk}+B_{rk}\)。因此还需要

\[
E[B_{rk}D_{rk}^H]=0,
\]

或令记录长度 \(T\to\infty\) 时，该均值偏差 \(b_T\to0\)。若要使用以 \(\sqrt R\) 为速率的中心极限定理，还需 \(\sqrt R\,b_T\to0\)，或显式估计/计入该偏差。简单丢弃一段过渡后做 FFT，不自动证明此条件。

若日志满足真实控制律，可恢复 \(d=u+Ky\)。这是代数恢复，真实反馈测量噪声可以准确抵消；但控制器近似、日志边界错误、饱和、未记录前馈、不同滤波/时间戳和未处理状态初值都会使恢复误差与扰动相关。恢复有效性必须核对，不能只看重构波形平滑。

## 2. 矩条件及全局目标辨识

记 \(A_k=E[U_kD_k^H]\)，则

\[
g(\theta):=E[q_r(\theta)]
=\mathcal R\operatorname{col}_k\operatorname{vec}
\{[G_{0k}-G_k(\theta)]A_k\},\qquad g(\theta_0)=0.
\tag{3}
\]

如果每个 \(A_k\) 都有满行秩 \(n_u\)，则对应频点的零矩等价于 \(G_k(\theta)=G_{0k}\)。如果没有满行秩，也可能依靠跨频率结构辨识某些目标，但必须另证；总输入谱满秩并不等于外部激励矩 \(A_k\) 满秩。

定义零矩集合

\[
\Theta_0=\{\theta\in\Theta:g(\theta)=0\},
\]

以及所需的规范不变目标 \(\eta(\theta)\)，例如确定频带内的相对通道动态、支承块左右子空间投影矩阵、指定模态的投影响应。目标应使用投影矩阵或规范固定后的坐标，避免把零空间向量的符号/旋转不唯一当作误差。

**全局目标可辨识条件为**

\[
\eta(\theta)=\eta_0\quad\text{对所有 }\theta\in\Theta_0.
\tag{4}
\]

式(4)不要求所有物理参数唯一，却不能省略。仅有 \(g(\theta_0)=0\) 不够：例如 \(g(\theta)=\theta^2-1\)，在 \(\theta_0=1\) 处导数非零，但 \(-1\) 仍为全局零矩点，目标 \(\eta=\theta\) 不唯一。

### 定理1：目标一致性

假设：

1. \(\Theta\) 紧致，\(g\) 与 \(\eta\) 连续，真实参数包含于 \(\Theta\)；
2. \(\sup_{\theta\in\Theta}\|\bar q_R(\theta)-g(\theta)\|\to_p0\)；
3. \(W_R\to_pW\succ0\)，且准全局极小解满足
   \(Q_R(\hat\theta_R)\le\inf_\Theta Q_R+o_p(1)\)，其中 \(Q_R=\bar q_R^TW_R\bar q_R\)；
4. 式(4)成立。

则 \(\eta(\hat\theta_R)\to_p\eta_0\)。也有 \(\operatorname{dist}(\hat\theta_R,\Theta_0)\to_p0\)，但不必有唯一物理参数极限。

**证明。** 连续性、紧致性和假设2–3给出 \(\sup_\Theta|Q_R-Q|\to_p0\)，其中 \(Q=g^TWg\)。对任意 \(\epsilon>0\)，集合 \(B_\epsilon=\{\theta:\|\eta(\theta)-\eta_0\|\ge\epsilon\}\) 紧致，并与 \(\Theta_0\) 不相交，所以 \(c_\epsilon:=\min_{B_\epsilon}Q>0\)。准全局最小性与 \(Q_R(\theta_0)\to_p0\) 表明 \(Q(\hat\theta_R)\to_p0\)，从而进入 \(B_\epsilon\) 的概率趋零。参数到零矩集合的距离同理。证毕。

**具体的一致大数定律条件。** 固定有限频点，独立同分布实验簇，\(G_k(\theta)\) 连续且在紧致参数集上一致有界，并有

\[
E\sum_k(\|Y_k\|+\sup_\Theta\|G_k(\theta)\|\|U_k\|)\|D_k\|<\infty,
\]

即可由连续函数族的可积包络得到所需 ULLN。若还具有一致有界一阶导数，可以直接用紧致网格覆盖和 Lipschitz 包络给出初等证明。独立簇可放宽到满足相应 ULLN/CLT 的平稳混合序列，但协方差必须改用长程协方差。密集插值出来的频点不能充当新的独立实验簇。

### 连续频带唯一性不能直接变成固定网格唯一性

设 \(h(\omega,\theta)=[G_0(\omega)-G(\omega,\theta)]A(\omega)\) 在紧致频带×紧致参数集上连续，且连续频带零响应足以确定目标。对固定 \(\epsilon>0\)，每个 \(\theta\in B_\epsilon\) 至少存在一个 \(\omega_\theta\) 使 \(h(\omega_\theta,\theta)\ne0\)。由连续性，这一严格不等式在该参数的邻域仍成立。\(B_\epsilon\) 紧致，故可取有限子覆盖，得到有限频率集 \(\Omega_\epsilon\)，使

\[
\min_{\theta\in B_\epsilon}\sum_{\omega\in\Omega_\epsilon}\|h(\omega,\theta)\|_F^2>0.
\tag{5}
\]

这证明**任意预先给定的非零目标分辨率存在一个有限分离网格**。它不证明同一个任意固定网格对所有 \(\epsilon\) 都足够，更不证明有限信噪比达到该分辨率。若采用越来越密的网格做精确一致性，还需均匀数值积分与三角阵列 ULLN；后面的固定维数 \(\chi^2\) 置信结论不能不加说明地用于增长的矩维数。

## 3. 可执行目标函数、Jacobian 与 GN/LM

对给定一轮固定权矩阵 \(W\)，优化

\[
\min_{\theta\in\Theta}\frac12\bar q_R(\theta)^TW\bar q_R(\theta).
\tag{6}
\]

第 \(\ell\) 列 Jacobian 为

\[
J_{R,\ell}(\theta)=
-\mathcal R\operatorname{col}_k\operatorname{vec}
\{\partial_{\theta_\ell}G_k(\theta)\,\widehat S_{ud,k}\},
\quad \widehat S_{ud,k}=R^{-1}\sum_rU_{rk}D_{rk}^H.
\tag{7}
\]

观测 \(U,D\) 是固定数据，对参数求导时不能再次将 \(U\) 当作模型生成信号求导。若优化另一种闭环均值模型，则需使用那种模型对应的 Jacobian，不能混用。

设局部允许增量 \(\delta\theta=Z\delta\xi\)，其中 \(Z\) 是结构约束/规范选择的切空间基，则 LM 子问题为

\[
\{(J_RZ)^TW(J_RZ)+\lambda S^TS\}\delta\xi
=-(J_RZ)^TW\bar q_R.
\tag{8}
\]

\(S\) 用于尺度归一化，\(\lambda\) 由信赖域下降比调节。稳定性、模态频率分离和可实现通道模型应由参数化/约束维护。LM 正则项是求解手段，不是新增观测信息；不能直接将带正则的逆 Hessian 当作无偏统计协方差。

例如对行/列纯时延，

\[
\partial_{\tau_{y,j}}G=-sE_jG,
\qquad \partial_{\tau_{u,k}}G=-sGE_k.
\]

机械模态结构和通道模型的其余导数可用解析求导或可靠自动微分。自由 SRIVC 可作为初始化或结构优化的内部步骤，但并未因此证明其不动点等于式(6)的驻点；若要声称两者等价，必须针对实际仪器、预滤波与权重重新核对。

**权重处理。** 先用固定正定权得到一致初估，再估计矩协方差，并在第二阶段固定其极限一致估计，是可行 GMM。数据依赖权并非一律破坏一致性：一致收敛到正定确定极限的 GMM 权可以成立。问题在于任意由含噪 FRF 模倒数生成的频点权，是否改变了总体矩/极限准则，以及是否证明其极限性质。随 \(\theta\) 连续更新的权也可以研究，但其真实目标梯度含权导数，不能每轮冻结后仍宣称精确优化同一个连续更新目标。

## 4. 正则商空间中的渐近协方差

先选择可识别的光滑局部坐标 \(\xi\in\mathbb R^p\)。假设真实点内部、结构秩局部恒定、参数规范不引入奇点，\(J=\partial_\xi g(\xi_0)\) 满列秩，\(\sqrt R\bar q_R(\xi_0)\Rightarrow N(0,\Omega)\)，一阶导数一致收敛，且优化一阶残差为 \(o_p(R^{-1/2})\)。则

\[
\sqrt R(\hat\xi-\xi_0)\Rightarrow N(0,V_\xi),
\quad
V_\xi=(J^TWJ)^{-1}J^TW\Omega WJ(J^TWJ)^{-1}.
\tag{9}
\]

证明由 \(0=J^TW\{\bar q_R(\xi_0)+J(\hat\xi-\xi_0)\}+o_p(R^{-1/2})\) 展开得到。当 \(W=\Omega^{-1}\) 时，\(V_\xi=(J^T\Omega^{-1}J)^{-1}\)。这只是给定矩类中的效率，不能不加条件地等同完整似然的 CRLB。

对 \(H=\partial_\xi\eta(\xi_0)\)，

\[
\operatorname{Cov}(\hat\eta)\simeq R^{-1}HV_\xi H^T.
\tag{10}
\]

在冗余原坐标下，线性化信息 \(I=J_\theta^T\Omega^{-1}J_\theta\) 可以奇异。必要的一阶目标可估条件是

\[
\ker I\subseteq\ker J_\eta.
\tag{11}
\]

只有在奇异性来自局部恒秩的平滑等价类、目标在这些等价类上恒定、并存在上述正则商坐标时，才能将目标协方差写成

\[
R^{-1}J_\eta I^\dagger J_\eta^T.
\tag{12}
\]

**式(11)本身不够。** 在 \(g(\theta)=\theta^2\)、\(\eta(\theta)=\theta^2\)、\(\theta_0=0\) 处，两个导数都为零，式(11)成立，但噪声矩反演可产生 \(\hat\theta=O_p(R^{-1/4})\)、\(\sqrt R\hat\eta\) 的截断非高斯极限，式(12)给出的零方差错误。类似地，极点碰撞、秩变化和活动不等式约束可能非正则。不能用伪逆掩盖这些情况。

更换坐标可以改善优化条件数，却不能增加真实信息。若时延与滤波常数弱可分，适当坐标有助数值求解，但不会使有限数据中不存在的区分能力凭空出现。

## 5. 不依赖参数强辨识的置信集：推荐主证书

对每个候选 \(\theta\)，用完整实验簇估计矩协方差：

\[
\widehat\Omega_R(\theta)=\frac1{R-1}\sum_r
\{q_r(\theta)-\bar q_R(\theta)\}
\{q_r(\theta)-\bar q_R(\theta)\}^T.
\tag{13}
\]

若簇有相关性，改用与相关结构相容的长程协方差/HAC 或簇重采样。对真正同簇的跨频率相关性，不能随意丢弃。单记录 LRM 窗口共享数据尤其需要处理相关性。

假设真实参数处 \(\Omega(\theta_0)\succ0\)、固定矩维数为 \(q\)、协方差估计一致，定义

\[
\mathcal C_R=\left\{\theta\in\Theta:
R\bar q_R(\theta)^T\widehat\Omega_R(\theta)^{-1}\bar q_R(\theta)
\le\chi^2_{q,1-\alpha}\right\}.
\tag{14}
\]

则 \(P(\theta_0\in\mathcal C_R)\to1-\alpha\)。证明只需在真实参数处将 CLT 白化后平方，不要求 \(J_\theta\) 满秩，也不要求估计参数渐近正态。它是绝对矩检验反演，**自由度为矩维数 \(q\)，不是过度识别检验的 \(q-p\)**。

该覆盖为固定真实对象下的点态渐近覆盖。若希望对趋近弱辨识/秩变化的一族对象保持均匀覆盖，还需统一 CLT、统一矩协方差非退化等条件，不能自动推断。

若矩协方差有已知的结构性奇异性，应先投影到固定非退化矩空间，再用其维数。数据驱动 SVD 阈值和收缩逆矩阵不能未经校准直接代入 \(\chi^2\) 分位数。特别是样本式(13)秩至多为 \(R-1\)：几次重复不能估计数千维任意协方差。实际可行方案是增加独立实验簇、预定矩压缩、使用有依据的谱噪声协方差结构，或对完整实验簇重采样并校准相应统计量；每个方案都需明确假设。

有限样本的严格概率保证仅在附加分布条件下成立。例如独立同分布高斯矩、\(R>q\) 时，Hotelling 统计量可用相应 F 分布而不是 \(\chi^2\)；一般实验数据仍是渐近证书。若矩偏差有可信上界，也可扩大接受域；未知模型失配不能靠 sandwich 协方差自动覆盖。

目标置信集为 \(\eta(\mathcal C_R)\)。参数方向很弱时，\(\mathcal C_R\) 可以很长，但若目标沿该方向近似不变，目标集合仍可能较小；若目标集合也宽，就应诚实报告现有数据不能给出所需隔离保证。

## 6. 从对象置信集到冻结校准器的隔离证书

令训练数据生成实际可实现的校准器 \(\widehat C\)，对应左、右变换 \(L(s;\widehat C)\)、\(R(s;\widehat C)\)。在验证和证书计算中固定它们，候选真实对象记为 \(\theta\)：

\[
T(s;\theta,\widehat C)=L(s;\widehat C)G(s;\theta)R(s;\widehat C).
\tag{15}
\]

按刚体逻辑坐标 \(r\) 和冗余/柔性坐标 \(n\) 分块。若理论目标是消除支承刚体分量，应对可辨识支承块 \(G_s\) 定义

\[
X_{rn}=L_rG_sR_n,\quad X_{nr}=L_nG_sR_r,
\quad F_t=L_nG_{f,t}R_n.
\tag{16}
\]

如果实验只能验证总传递函数，就报告总交叉块 \(L_rGR_n,L_nGR_r\)，并明确它可能含真实柔性耦合；不能将其全部称为刚体泄漏。支承/指定柔性块的归属需要频谱分解的唯一性和物理依据。

频带范数、坐标单位和允许逻辑指令幅值应提前归一化，否则转角与位移、力矩与力的交叉块不能直接比较。定义带权频带范数后，取

\[
B_{rn}=\sup_{\theta\in\mathcal C_R}\|X_{rn}(\theta,\widehat C)\|_{\mathcal B},
\quad B_{nr}=\sup_{\theta\in\mathcal C_R}\|X_{nr}(\theta,\widehat C)\|_{\mathcal B},
\]

\[
A_t=\inf_{\theta\in\mathcal C_R}\|F_t(\theta,\widehat C)\|_{\mathcal B}.
\tag{17}
\]

在 \(\theta_0\in\mathcal C_R\) 的事件上，真实泄漏不超过 \(B_{rn},B_{nr}\)，目标柔性响应不低于 \(A_t\)。若 \(A_t>0\)，比例 \(B_{rn}/A_t,B_{nr}/A_t\) 即为相应最坏相对界。也可直接对归一化泄漏比取 supremum，通常更不保守。设定工程阈值后，只有整个候选集合均满足阈值才给出“已达到精度”的判定。

### 为什么同数据校准不必自动导致证书无效

式(17)对每个样本逐点成立，因此 \(\widehat C\) 即使由同一数据选择，只要 \(\mathcal C_R\) 覆盖真实对象，最坏界仍覆盖真实性能；不必为这一步强制独立分样。问题在于不能把式(15)改成 \(L(\theta)G(\theta)R(\theta)\)：每个候选对象都用它自己的精确零空间会造成自证。

例如模型定义 \(L_n(\theta)C_r(\theta)=0\)，则拟合内支承泄漏永远为零，其全导数也可能为零；这不代表真实对象上的误差方差为零。实际证书的对象导数必须是

\[
\partial_\theta T\big|_{\widehat C\text{ 固定}}
=L(\widehat C)(\partial_\theta G)R(\widehat C).
\tag{18}
\]

如果研究的是校准矩阵估计值自身的分布，才需要对 \(C(\theta)\) 求导。两类目标不能混用。独立验证仍十分有价值：它检查参数模型之外的失配、非线性、漂移、错误边界和选择偏差，尤其避免在同一个 FRF 估计上调到满意再宣称实机有效。

### 计算与实现上的边界

- 非线性置信集上的最坏值通常非凸。多起点优化仅找到一个较大的值，不保证得到上界；严格数值证书需全局上界、区间界、可信外包集合或其他保守方法。随机采样候选集合不是最坏性能证明。
- 正则局部情形可用椭球和 Taylor 上界：若 \(\delta\xi^TV^{-1}\delta\xi\le c\)，且 \(\|D^2f\|\le L_f\)，则
  \[
  \|f(\hat\xi+\delta\xi)\|\le\|f(\hat\xi)\|+
  \sqrt{c\lambda_{\max}(J_fVJ_f^T)}+
  \tfrac12L_fc\lambda_{\max}(V).
  \]
  省掉二阶余项只能称一阶近似证书。
- 若真实对象存在有界频率失配 \(\|\Delta_k\|_2\le\epsilon_k\)，相应变换块需额外加入 \(\|L_k\|_2\epsilon_k\|R_k\|_2\) 的界。协方差不能替代未知模型偏差。
- 补偿后的噪声增益必须按实际噪声注入位置计算；若反馈改变噪声传递，要包含闭环灵敏度，而不只计算左变换矩阵范数。
- 频率离散网格上的界不自动覆盖整个连续频带；可用有根据的频率 Lipschitz 上界补上网格间隔误差。
- 离线对象通道隔离不等于部署后的闭环稳定。若校准器改变反馈回路，还须独立检查整个置信对象集合下的闭环稳定性和控制输入限制。

## 7. 可以作为正式结论与尚不能声称的结论

**现阶段可严格证明：** 外部参考零矩；目标/等价类一致性定理；给定目标分辨率的有限网格存在性；正则商空间 GN 与 sandwich；无需参数强辨识的绝对矩置信集；冻结数据生成校准器后的集合性能证书。

**仍需特定模型证明或后续数据核实：** 实际参考有效性；所选频点和输入方向下的目标分离；物理参数规范与重复极点处的正则性；支承块与目标柔性块分解唯一性；噪声/协方差模型；矩偏差和机械模型失配边界；置信集合的可计算最坏上界；可实现补偿器与反馈稳定性。

所以合理的“数学闭环”结论是：**在明确且可核查的观测、可辨识和模型条件下，可从闭环数据得到有覆盖率的对象集合，再得到刚柔隔离与柔性保留的性能界；如果数据不足，界会变宽并拒绝给出通过结论。** 不能把“已经写出置信公式”说成“当前实验一定能达到所需精度”。

## 8. 来源与标准理论边界

本文件的目标分离证明、有限网格覆盖论证和冻结校准器集合推论为针对本项目直接推导；GMM、delta method、矩检验反演与置信集投影属于标准统计框架，不应单独包装为创新。

- Newey & McFadden, *Large Sample Estimation and Hypothesis Testing*, Handbook of Econometrics IV (1994)：极值估计的一致性、ULLN 和 GMM 渐近理论。[出版社](https://www.sciencedirect.com/science/article/pii/S1573441205800054)，[高校保存全文](https://www.ssc.wisc.edu/~xshi/econ715/chap36neweymacfadden.pdf)。
- Hansen, *Large Sample Properties of Generalized Method of Moments Estimators*, Econometrica 50 (1982),1029–1054：GMM 基础。[DOI](https://doi.org/10.2307/1912775)，[作者文献页](https://home.uchicago.edu/~lhansen/papers.html)。
- González et al., *Consistency analysis of refined instrumental variable methods for continuous-time system identification in closed-loop*, Automatica 166 (2024),111697：闭环 RIV 条件不能省略。[DOI](https://doi.org/10.1016/j.automatica.2024.111697)。
- González et al., *Statistically Optimal Structured Additive MIMO Continuous-time System Identification*, 作者预印本 (2025)：结构化加性与闭环 IV 已有基线。[原文](https://arxiv.org/abs/2505.14169)。
- van der Hulst et al., *Structured identification of multivariable modal systems*, MSSP 247 (2026),113948：低秩结构化模态细化与闭环实验基线。[DOI](https://doi.org/10.1016/j.ymssp.2026.113948)。
