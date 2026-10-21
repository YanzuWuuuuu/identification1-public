# 第 4 章独立推导审阅：IPEM、秩结构、协方差与时延

日期：2026-10-12。范围：只读论文与已有代码背景，不改算法，不运行共享 MATLAB 会话，不运行 Monte Carlo。

本文区分三种证据：**原文明确陈述**、**由公式推导出的适用边界**、**尚需文献与实验验证的研究候选**。发现论文公式问题不等于发现可发表创新。

## 1. 来源与页码核对

- Thesis：article/1486187-ModalIdentification_of_MultivariableSystems-42-76.pdf。节选 PDF 第 13–35 页对应印刷页 45–67，第 4 章。
- 更新论文：M. van der Hulst et al., *Structured identification of multivariable modal systems*, MSSP 247 (2026), 113948，DOI：10.1016/j.ymssp.2026.113948。文件位于 article/Structured identification of multivariable modal systems.pdf。
- 已渲染 thesis PDF 第 20、21、22、23、24、33 页，以及 MSSP 第 5、8、9 页；直接看图核对了主要公式，未只依赖提取文本。
- IPEM 原始文献：Söderström, Stoica, Friedlander (1991), *An indirect prediction error method for system identification*，https://doi.org/10.1016/0005-1098(91)90019-X 。两阶段与直接 PEM 的渐近统计性质一致是已有理论，不可作为新贡献。

## 2. 原文到底近似了什么

| 项目 | 原文位置 | 判定 | 2026 更新 |
|---|---|---|---|
| 将高阶未辨识模态视为静态增益 D | thesis p51，(4.20)–(4.23) 后文字 | 明确模型近似 | 仍采用，MSSP p3 文字及 p12 §6.3 |
| SVD 单模态秩一截断 | thesis p46；p55 (4.46)–(4.48) | 明确仅为 Frobenius 范数最优，带权或近邻共振不保证好 | 保留为初始化，再 IPEM；最终方法并非只做 SVD |
| 用参数空间二次距离替代原始数据代价 | thesis p52 §4.4.1，p53 (4.30) | 明确局部近似，在第一阶段最优参数附近最好 | MSSP p5 (22) 仍为该方法 |
| 以协方差逆作为二次度量 | thesis p53 (4.31) | 需要统计权重、噪声与正则条件；不应无条件等同任意工程权重的 Hessian | MSSP p8 (39) 明称协方差近似，详细推导不在文中 |
| 比例/一般黏性阻尼 | thesis p48–52，(4.6)、(4.17)、(4.24)、(4.27) | 已支持一般阻尼与复振型；不能说原方法只支持实振型 | 2026 两者均支持 |
| 映射 injective、Jacobian 满秩、参数可辨识 | thesis p52 最后一段 | 理论假设；无规范因子参数化不自动满足 | MSSP p5 Remark 1 补充单振型尺度规范，多刚体共同分母的混合自由度仍需区分 |
| 非凸优化与初始化 | thesis p56、p61 | 明确局部最优风险，仿真在真值附近初始化 | 2026 用 CMIF 初始化及实际晶圆台验证，仍无全局收敛保证 |

原文“一步秩约束方法的统计保证较难”是选择两阶段的理由（p46），不等于证明不存在一步一致/有效估计器。2026 的 RIV 推导解决加性模型驻点方程，再做结构映射；并未证明“在每轮 RIV 中加入任意秩约束后仍保留原统计结论”。

## 3. 推导 A：漏掉第五阶为何可能破坏两阶段等价

为避免复数转置及因子 2 歧义，将所有频率、通道的实部和虚部堆叠为实向量。记观测 z、加性模型 g(β)、残差 r=z−g(β)，固定对称正定权重 W：

\[
Q(\beta)=\tfrac12 r(\beta)^\top W r(\beta),\qquad
\widehat\beta=\arg\min_\beta Q(\beta).
\]

令 \(D=\partial g/\partial\beta^\top\)。在驻点附近：

\[
Q(f(\rho))
=Q(\widehat\beta)
+\tfrac12\Delta^\top H\Delta
+O(\|\Delta\|^3),\qquad
\Delta=f(\rho)-\widehat\beta,
\]

\[
H=D^\top WD-\sum_i(Wr)_i\,\nabla_\beta^2g_i.
\]

GN 只保留 \(H_{\rm GN}=D^\top WD\)。这里有两个近似：从原数据代价到参数空间二次代价的 Taylor 余项，以及从精确 Hessian 到 GN 时省去的“残差 × 二阶导数”项。

若真系统属于嵌套模态模型集、参数可辨识、第一阶段一致、统计权重正确，约束投影距离通常随样本数减小，余项渐近变小。这是两阶段理论的合理场景，有限样本也不意味着逐点完全相同。

若真实系统含第五阶、拟合固定四阶，真传递函数一般不在四阶模态集，通常也不在四阶加性集。第一阶段趋向伪真值 \(\beta_A^\star\)，它到四阶模态流形的距离不保证随样本数消失：

\[
\|\beta_A^\star-f(\rho^\star)\|=O(1)
\]

是可能的。所以增加数据量不自动消除局部替代误差。这是**由原文假设推出的边界**，不是原文给出的失败定理，也不是说四阶拟合总会失败。

可检验：

- 参数空间预测的代价增量与原数据实际增量的比值；
- 原数据代价在两阶段输出处的梯度；
- 两阶段输出与从该点启动直接结构化 PEM 的差异。

若结构投影代价下降、原数据代价却明显增加，应区分局部代价替代失真与优化未收敛。

### 3.1 静态 D 对紧邻频带的第五阶缺少小参数保证

thesis p51 明确用 D 近似高阶模态。单个漏模：

\[
\Delta_5(s)=\frac{B_5}{s^2+2\zeta_5\omega_5s+\omega_5^2}
=\frac{B_5}{\omega_5^2}
\left[
1-\frac{2\zeta_5}{\omega_5}s
+\frac{4\zeta_5^2-1}{\omega_5^2}s^2+\cdots
\right].
\]

只保留 D 要求工作频率相对漏模频率足够低，不能仅要求漏模频率略高于截止频率。当前 720 Hz 截止与名义 738.83 Hz 第五阶的比例约 0.975；按此前 ode3 等效第五阶约 755.50 Hz 也约 0.953。频带上端的静态近似没有小参数保证。

这支持研究“动态残差 + 时延 + 四个目标模态”。但残差复杂度不能无限增加，也不能暗中把第五阶完整辨识出来再宣称无需建模第五阶。

## 4. 推导 B：时延如何吸收漏模，如何定量区分

考虑共同时延：

\[
G_{\rm true}(s)=e^{-s\tau_0}
\{G_4(s,\theta_0)+\Delta_5(s)\}.
\]

拟合 \(G(s,\theta,\tau)=e^{-s\tau}G_4(s,\theta)\)。在基点附近实化并白化：

\[
\widetilde r_\Delta=W^{1/2}\operatorname{stack}
\{e^{-s\tau_0}\Delta_5(s)\},
\quad
A=W^{1/2}\frac{\partial g}{\partial\theta^\top},
\quad
b=W^{1/2}\frac{\partial g}{\partial\tau},
\]

\[
\frac{\partial G}{\partial\tau}
=-s\,e^{-s\tau}G_4(s,\theta).
\]

\(\theta\) 应使用去除冗余后的可辨识坐标。令 \(P_\perp=I-AA^\dagger\)，局部小失配线性化下：

\[
\boxed{\delta\tau\simeq
\frac{b^\top P_\perp\widetilde r_\Delta}
{b^\top P_\perp b}}
\]

推导：在线性模型
\(\widetilde r_\Delta\simeq A\delta\theta+b\delta\tau\)
中对 \(\delta\theta\) 最小二乘消元，再解余下标量问题。

- 分子：漏模与“模态参数解释之后剩余时延方向”的相关程度。
- 分母：时延在模态参数之外的独立信息。
- 分母很小意味着小 FRF 误差也无法保障时延可信；分子大意味着系统性漏模偏差。
- 这是局部诊断，不是全局偏差公式。一般工程 W 下分母只是条件曲率；正确噪声白化及正则条件下才可解释为 Fisher 信息。

若加残差模型 \(R(s,\gamma)\)，白化灵敏度为 B：

\[
\mathcal I_{\tau\mid\theta,\gamma}
=b^\top\{I-[A\ B][A\ B]^\dagger\}b.
\]

残差模型减轻偏差，却可能减少时延独立信息。因此可研究的核心是**残差表达能力与时延/目标模态可辨识性的权衡**。如果残差任意吸收时延方向，就不能声称二者同时可辨识。

一个具体候选：设计有限复杂度带外动态残差，按条件曲率/信息或最小奇异值选频带与复杂度，给出可辨识区间及局部偏差界。必须相对联合 PEM、时延网格搜索、静态 D 与常规动态残差证明优势。

### 4.1 时延应和第一阶段参数联合传递

第一阶段估计 \(\eta=(\beta,\tau)\) 后，第二阶段可写：

\[
\min_{\rho,\tau}
\left(
\begin{bmatrix}\widehat\beta\\\widehat\tau\end{bmatrix}
-
\begin{bmatrix}f(\rho)\\\tau\end{bmatrix}
\right)^\top C_\eta^{-1}
\left(
\begin{bmatrix}\widehat\beta\\\widehat\tau\end{bmatrix}
-
\begin{bmatrix}f(\rho)\\\tau\end{bmatrix}
\right).
\]

冻结时延、舍弃交叉协方差 \(C_{\beta\tau}\)，通常不是同一问题。正确似然度量下 nuisance 消元对应 Schur 补：

\[
\mathcal I_{\rho\mid\nu}
=\mathcal I_{\rho\rho}
-\mathcal I_{\rho\nu}
\mathcal I_{\nu\nu}^{-1}\mathcal I_{\nu\rho},
\qquad \nu=(\tau,\gamma).
\]

Schur 补与信息投影本身是标准工具，不是创新；候选贡献在于将其发展成此场景可验证的误差、选带、可辨识结论。

## 5. 推导 C：Hessian、协方差及确定性偏差不能混同

固定设计下，噪声协方差为 C，模型正确且局部线性化有效。记 \(H_{\rm GN}=D^\top WD\)：

\[
\operatorname{Cov}(\widehat\beta)
\simeq H_{\rm GN}^{-1}
D^\top WCWD\,H_{\rm GN}^{-1}.
\]

仅当 \(W=C^{-1}\)（或相应比例和正则条件）时化为：

\[
\operatorname{Cov}(\widehat\beta)
\simeq(D^\top C^{-1}D)^{-1}.
\]

工程权重 \(1/|G|\)、\(1/|G|^2\) 一般不满足上述条件。非零失配还影响精确 Hessian、伪真值及稳健方差；W 若由同一份带噪 FRF 得到，其随机性也不能自动忽略。

要区分：

1. 近似原工程加权代价，应保留原代价局部曲率。
2. 按第一阶段不确定度作最小距离结构估计，应使用有效参数协方差。

特殊统计加权下二者吻合，任意权重下不必吻合。

**Sandwich 修正的是随机波动的协方差，不会消除第五阶造成的确定性模型偏差。** 若置信区间围绕伪真值构造，即使协方差估计正确，也不保证覆盖物理真值。只有建模/约束/偏差界/稳健目标解决了相应偏差，才能讨论物理参数覆盖率。

### 5.1 已核对的公式和叙述疑点

**MSSP p8 式 (39)：**

\[
\widehat\Sigma_{\widehat\beta}
=\left[
\frac1N\sum_k\widehat\Phi_k
\widehat\Sigma_{\widehat G,k}\widehat\Phi_k^H
\right]^{-1}.
\]

正文将 \(\widehat\Sigma_{\widehat G}\) 定义为 FRF covariance，式 (26)–(28) 的 \(\widehat\Phi\) 是未含噪声白化的 Jacobian。按当前记号，该式可疑：把噪声协方差乘 c，参数协方差反而除 c。thesis p53 (4.31) 对应位置使用 \(C_{\operatorname{vec}G}^{-1}\)。可能是排版/记号问题，不能替作者断定原因，应从估计方程重推。

**N 的归一化：** 形如 \([\frac1N\sum D_k^\top C_k^{-1}D_k]^{-1}\) 常表示 \(\sqrt N(\widehat\beta-\beta_0)\) 的渐近协方差；估计量本身多一个 \(1/N\)。整体倍数不改变 IPEM 点估计，却改变置信区间。

**thesis p65 方差解释过强：** 本仿真分母系数方差几乎不变，不可普遍归因于分母参数数量不变。若极点与留数估计相关，约束留数也会改善极点方差。有效估计且局部线性化下，在加性坐标中：

\[
C_{\rm constrained}
=T(T^\top C^{-1}T)^{-1}T^\top,\qquad
T=\partial f/\partial\rho^\top.
\]

其极点子块不必等于 C 的极点子块。最小反例：两个参数 a,b 的协方差为
\(\left[\begin{smallmatrix}1&0.8\\0.8&1\end{smallmatrix}\right]\)，若真约束 b=0，最优投影为 \(a_c=\widehat a-0.8\widehat b\)，方差由 1 降到 0.36，尽管 a 的参数个数未变。

另外，单看权重矩阵不同参数的数值数量级会受单位和缩放影响，不能作为坐标不变的统计解释。

## 6. 秩结构与可辨识：已解决与仍需区分

### 6.1 单模态尺度：2026 已补充，不能当创新

\[
B=lr^\top=(cl)(r/c)^\top,\quad c\ne0.
\]

沿 \(\delta l=l,\delta r=-r\)，有 \(\delta B=0\)，未经规范的 Jacobian 必有零方向。复因子还有相位自由度。thesis p52 的 injective/full-rank 条件不能直接套到无规范的 (4.26)/(4.29)。

MSSP p5 Remark 1 已要求固定系数或范数，消除尺度不唯一。“发现需振型归一化”不是贡献。复数范数本身仍不固定相位，需要相位/系数规范；实数范数后离散符号歧义未必破坏局部可辨识。

### 6.2 多刚体混合自由度，比逐列尺度更大

\[
M=LR^\top,\qquad
(L,R)\mapsto(LT,RT^{-\top}),\quad T\in GL(r).
\]

同一 M 对应 \(r^2\) 维连续规范，固定秩矩阵流形维数为 \(r(p+q-r)\)，而非 \(r(p+q)\)。

当前 p=q=4，r=3：刚体因子 24 个数，只有 15 个有效方向。四个简单比例阻尼柔性模态每个 \(2+4+4-1=9\) 维，共 51 个有效方向；无 D 的 64 维 rho 至少有 13 维规范零方向。

SVD 选择数值基，不代表仅凭共同分母 FRF 唯一识别每条物理平移/转动方向。几何、质量等额外物理约束可以选基；仅逐列归一化不足以消除全部混合自由度。

2026 p8 已使用“刚体子空间”表述。稳妥输出是 M、可控/可观子空间，或明确物理规范下的振型。

### 6.3 重根与近邻模态

\[
\frac{B_1}{a(s)}+\frac{B_2}{a(s)}
=\frac{B_1+B_2}{a(s)}.
\]

分母完全相同，只能直接识别合并留数；独立 B1、B2 分解不唯一。近邻模态把精确不可辨识变成病态信息。CMIF 显示重数（thesis p57）不自动提供唯一分解或可靠协方差。

研究对象可改为模态簇/子空间并纳入共同时延。但相关研究多，优先级低于漏模-时延混淆。

### 6.4 已有研究的排重提醒

- 固定秩/商流形：Mishra et al., *Fixed-rank matrix factorizations and Riemannian low-rank optimization*，https://arxiv.org/abs/1209.0430 。
- 一般阻尼、复振型、近邻模态和 Hessian 不确定度：Rather & Bansal (2025), *Bayesian modal identification of non-classically damped systems using frequency domain data*，https://doi.org/10.1016/j.jsv.2025.119039 。
- 近邻模态两阶段 Bayesian FFT/JAD：Xie et al. (2025), *Efficient two-stage modal identification for structures with closely spaced modes by Bayesian FFT and joint approximate diagonalization*，https://doi.org/10.1016/j.jsv.2025.119365 。

后两者以 operational/ambient 模态识别为背景，未必完全覆盖强制激励 MIMO RIV；但“近邻模态+两阶段+不确定度”不能直接声称首创。

## 7. Thesis 两处复数公式，不能作为研究创新

由
\[
\frac{B_1s+B_0}{(s-\lambda)(s-\bar\lambda)}
=\frac{R}{s-\lambda}+\frac{\bar R}{s-\bar\lambda}
\]
可直接得到
\[
R=\frac{B_0+\lambda B_1}{\lambda-\bar\lambda}.
\]

thesis p55 (4.47) 印成
\(\lambda(B_0+B_1/\lambda)/(\lambda-\bar\lambda)\)，与 (4.22) 的系数约定不符；MSSP p8 (40) 已纠正。

若 \(R=U\Sigma V^H\)，又要写成 \(\psi_l\psi_r^\top\)，应取：
\[
\psi_l=\sqrt{\sigma_1}u_1,\qquad
\psi_r=\sqrt{\sigma_1}\bar v_1.
\]
thesis p55 (4.48) 直接取右奇异向量，按该记号需要共轭；实振型时区别消失。2026 未给此处显式向量式，不能据此判定其实现仍错。

## 8. 三个候选，按优先级排序

### 候选 1：邻近漏模下的时延-目标模态解耦

主问题：四阶模型无法完整表示第五阶，第四阶频率、阻尼和共同延时互相补偿。

候选贡献：有限复杂度动态残差；局部漏模偏差与条件信息推导；据此选频带/正则化；给出目标四阶与时延的误差及置信覆盖率。

比较基线：四阶无延时、四阶有延时、四阶+D+延时、四阶+常规动态残差+延时、五阶+延时参考。统一权重、数据、采样/时延假设、初始化及计算预算。

### 候选 2：失配时两阶段压缩有效性诊断与修正

主问题：加性最优点未必靠近结构流形；工程权重下曲率与协方差不是同一对象。

候选贡献：分开保存原代价曲率和统计协方差；建立可计算投影失真指标；必要时触发少量原数据结构化修正，保留时延交叉信息。

GN 精化、sandwich、delta method 本身不是创新；需有面向此失配场景的误差界/适用判据及成本收益。Sandwich 不消除确定性偏差。

### 候选 3：规范不变的模态簇估计

用留数矩阵/子空间为主要对象，处理多刚体、近重根、复振型规范，和共同时延联合估计。可作为前两者的理论支撑。固定秩几何和 Bayesian 近邻模态文献多，独立主线的新颖性风险较高。

## 9. 后续最小验证要求

1. 精确采样/ZOH 一致仿真，单独处理 ode3 极点偏移，避免把数值误差当贡献。
2. 真模型包含于模型集的五阶场景与故意四阶失配场景分开报告。
3. 对照使用同样原始/FRF 数据与同样目标函数；不能把权重变化的收益误归于优化方法。
4. 区分实际收敛、达到迭代上限、触边界和不可辨识。
5. 报告频率/阻尼/留数或子空间/时延的 bias、RMSE、覆盖率，以及结构代价、原代价和运算量。
6. 扫描第四/第五阶间距、带宽、时延、SNR、权重与残差复杂度，验证理论预计的混淆区域，不能只演示单次 100 µs 成功。

这些是有数学动机的研究候选，尚不是新颖性及性能已经确认的结论。
