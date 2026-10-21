# SRIVC 从数据到加性模态模型：与本项目代码逐项对应的学习笔记

日期：2026-10-17。本文核对了当前项目的 `frf_srivc_fit.m`、`raw_srivc_fit.m`、`raw_srivc_design.m`、新扫频预处理及拟合脚本，并对照 González 等的 2025 年 Automatica 论文和 van der Hulst 等的 2026 年 MSSP 论文。本文不新增实验结果，也不把后续结构优化称作已经完成的 SRIVC 理论创新。

## 1. 先把三个问题分开

1. **辨识问题：** 已知输入输出数据，估计一个能够预测频率响应或时间响应的模型。
2. **结构问题：** 判断这些二阶拟合分量是否能够表示真实的简单机械模态、共同支承空间及通道动态。
3. **解耦问题：** 利用模型和验证数据设计输入混合、输出混合或动态滤波器，使指定通道之间的影响减小，同时保留需要的运动响应。

SRIVC 主要解决第 1 个问题。它给解耦提供对象模型，不会自动保证解耦成功。特别是四个大响应相减形成的小零空间响应，可能比原始 16 路频响更难辨识。

## 2. 我们到底在估计什么模型

令 \(n_u\) 为输入数、\(n_y\) 为输出数、\(n_m\) 为二阶分量数。未包含通道时延的自由加性模型写为

\[
P(s)=D+\sum_{m=1}^{n_m}\frac{B_{m,0}+sB_{m,1}}{s^2+2\zeta_m\omega_m s+\omega_m^2},
\tag{1}
\]

其中 \(B_{m,0},B_{m,1},D\in\mathbb R^{n_y\times n_u}\)。同一个分量的频率、阻尼对所有输入输出通道相同，分子则描述该分量被不同执行器激发、被不同传感器观测到的程度和相位关系。

例如对第 \(j\) 个传感器、第 \(k\) 个电机，

\[
P_{jk}(s)=D_{jk}+\sum_m\frac{b_{m,0,jk}+sb_{m,1,jk}}{s^2+2\zeta_m\omega_m s+\omega_m^2}.
\]

这是“16 路共同辨识”，而不是互不相关地拟合 16 个 SISO 模型。共同极点使各路对同一机械共振的信息可以合并。

本项目使用归一化频率坐标改善数值尺度。对每个分量取固定参考频率 \(f_{m,\mathrm{ref}}>0\)，定义

\[
q_m=\frac{s}{\omega_{m,\mathrm{ref}}},\quad
A_m(q_m)=1+a_{m,1}q_m+a_{m,2}q_m^2,
\quad \omega_{m,\mathrm{ref}}=2\pi f_{m,\mathrm{ref}},
\tag{2}
\]

\[
P(s)=D+\sum_m\frac{C_{m,0}+q_mC_{m,1}}{A_m(q_m)}.
\tag{3}
\]

最终物理参数由

\[
f_m=\frac{f_{m,\mathrm{ref}}}{\sqrt{a_{m,2}}},\qquad
\zeta_m=\frac{a_{m,1}}{2\sqrt{a_{m,2}}},
\tag{4}
\]

\[
B_{m,0}=\frac{\omega_{m,\mathrm{ref}}^2}{a_{m,2}}C_{m,0},\qquad
B_{m,1}=\frac{\omega_{m,\mathrm{ref}}}{a_{m,2}}C_{m,1}
\tag{5}
\]

恢复。`coefficients` 内的 \(C_{m,0},C_{m,1}\) 不是式 (1) 的物理分子，使用时不能漏掉这一换算。

**参数规模。** 当前 4 输入、4 输出、15 个自由二阶分量，分母共有 \(2\times15=30\) 个实参数，分子有 \(2\times15\times16=480\) 个，加上 16 个实直通系数，合计 526 个实参数，尚未计入时延。自由矩阵分子不自动满足简单机械模态的秩一条件。因此“15 个拟合分量”不等于“已经确认 15 个独立机械模态”，也不能自动理解为一个仅 30 个状态的最小 MIMO 机械实现。

## 3. 为什么普通最小二乘不够：先看一个二阶 SISO 模型

用 \(q=p/\omega_{\mathrm{ref}}\)，\(p=d/dt\)，设

\[
y(t)=\frac{b_0+b_1q}{1+a_1q+a_2q^2}u(t)+v(t).
\tag{6}
\]

将分母移到左边，可以整理为关于参数的线性关系，但回归量中包含 \(qy,q^2y\)。测量输出 \(y\) 本身含噪，回归量也含噪；直接最小二乘通常会出现回归量与误差相关的问题。此外直接差分求导会放大高频噪声。

SRIVC 采用当前分母构造滤波器。第 \(r\) 次迭代令

\[
F^{(r)}(q)=\frac{1}{A^{(r)}(q)},\quad
\theta=[a_1,a_2,b_0,b_1]^T,
\]

构造被解释量及回归向量

\[
y_f^{(r)}=F^{(r)}y,
\]

\[
\phi_f^{(r)}=
\begin{bmatrix}
-qF^{(r)}y\\-q^2F^{(r)}y\\F^{(r)}u\\qF^{(r)}u
\end{bmatrix}.
\tag{7}
\]

这里的 \(qF^{(r)},q^2F^{(r)}\) 是滤波器，不是先对含噪数据做数值微分。对于当前二阶稳定分母，这些滤波器可以直接用状态空间实现。

由式 (6) 可验证

\[
y_f^{(r)}=(\phi_f^{(r)})^T\theta
+\frac{A(q)}{A^{(r)}(q)}v.
\tag{8}
\]

即使有了滤波，直接用 \(\phi_f\) 做最小二乘也不自动消除噪声相关性。因此先用当前模型生成无测量噪声的预测

\[
\widehat x^{(r)}=G(q,\theta^{(r)})u,
\]

再构造工具变量

\[
z_f^{(r)}=
\begin{bmatrix}
-qF^{(r)}\widehat x^{(r)}\\
-q^2F^{(r)}\widehat x^{(r)}\\
F^{(r)}u\\qF^{(r)}u
\end{bmatrix}.
\tag{9}
\]

在开环、输入与输出噪声独立、模型结构和采样处理正确等条件下，\(z_f\) 与回归量相关，却能避免直接使用含噪输出作为工具变量。新的参数由

\[
\theta^{(r+1)}=
\left[\sum_t z_f^{(r)}(t)(\phi_f^{(r)}(t))^T\right]^{-1}
\sum_t z_f^{(r)}(t)y_f^{(r)}(t)
\tag{10}
\]

更新。实际实现应采用分解求解而非显式求逆。更新模型后重新生成滤波器、预测输出和工具变量，循环至收敛。

这解释了名称中的 refined：工具变量和预滤波器会随当前模型改善。RIVC 通常还涉及噪声模型及相应预白化；SRIVC 使用简化的噪声设定。对于本项目 FRF 数据上的实现，最严谨的称法是“频域加性 RIV／SRIVC 型迭代”，不能仅凭名称赋予它任意有色噪声、任意闭环下的最优性。

## 4. 多个模态相加时：部分输出是关键

对

\[
y=Du+\sum_{m=1}^{n_m}G_m u+v,
\]

单独研究第 \(m\) 个分量时，先从测量输出减去其他分量的当前预测：

\[
\widetilde y_m^{(r)}
=y-D^{(r)}u-\sum_{\ell\ne m}G_\ell^{(r)}u.
\tag{11}
\]

令总残差 \(e^{(r)}=y-\widehat y^{(r)}\)，第 \(m\) 个分量的预测为 \(\widehat y_m^{(r)}=G_m^{(r)}u\)，则

\[
\widetilde y_m^{(r)}=\widehat y_m^{(r)}+e^{(r)}.
\tag{12}
\]

于是，第 \(m\) 个分量的测量回归量用 \(\widetilde y_m\)，工具变量用 \(\widehat y_m\)。注意 \(\widetilde y_m\) 是算法中的“扣除其他分量后的部分输出”，并不是实验上直接测得的纯模态输出；模型未收敛时，它还含有其他分量的估计误差。

这正对应当前代码的核心三行：

```matlab
part = B(:,ix)*C(ix,c);
z(:,ix) = -filt.*part;
x(:,ix) = -filt.*(part+res(:,c));
```

其中 `part` 是当前模态预测，`part+res` 是含测量残差的部分输出。若把两者都换成 `part`，就不是当前 SRIVC 更新了。

## 5. 搬到频域后，公式为何更简单

在频率 \(f_k\) 上，令 \(s_k=\mathrm j2\pi f_k\)，则 \(q_m=\mathrm jf_k/f_{m,\mathrm{ref}}\)。微分与滤波变为复数乘除。

### 5.1 以非参数 FRF 为输入

数据为 \(\widehat H_{jk}(\mathrm j\omega)\)。每条 FRF 可以理解为单位输入下的复数响应。对通道 \(c\) 的加性模型，定义

\[
\mathcal B_k(a)=
\begin{bmatrix}
A_1^{-1}&q_1A_1^{-1}&\cdots&A_{n_m}^{-1}&q_{n_m}A_{n_m}^{-1}&1
\end{bmatrix},
\tag{13}
\]

\[
\widehat P_c=\mathcal B(a)c_c,\qquad
e_c=\widehat H_c-\widehat P_c.
\tag{14}
\]

分母 \(a\) 对所有通道共用；每条通道有独立的实系数向量 \(c_c\)。先给定模态初始频率和阻尼，固定分母后，分子估计是一次实参数线性最小二乘。这只是初始化，不是整个 SRIVC。

### 5.2 以原始输入输出 DFT 为输入

无需先生成 FRF，也可以直接拟合

\[
Y_j(\omega_k)=\sum_{\ell=1}^{n_u}P_{j\ell}(\mathrm j\omega_k)U_\ell(\omega_k)
+T_j(\omega_k)+E_j(\omega_k),
\tag{15}
\]

其中 \(T_j\) 用于描述有限记录的初末状态、截断等项。此时基函数由 \(1/A_m,q_m/A_m\) 变成 \(U_\ell/A_m,q_mU_\ell/A_m\)。

项目的 `raw_srivc_design.m` 除这些输入相关列外，还放入了每个分量的 \(1/A_m,q_m/A_m\) 边界基函数，以及额外的常数和频率一次项。它们是有限记录中的干扰／边界描述，不应直接当作机械传递函数的分子或直通项导入 Simulink。

“不经过 LRM”与“无闭环偏差”不是同一个概念；原始总输入 \(U\) 仍可能和闭环噪声相关。采样保持方式、FFT 窗口和连续时间模型的关系也必须核对。直接把有限记录 DFT 的微分等同于 \(\mathrm j\omega\) 相乘而不处理边界，通常不是严格等式。

## 6. 复数数据如何求实参数

所有机械模型系数要求为实数，而 FRF 是复数。定义

\[
\mathscr R(v)=\begin{bmatrix}\Re v\\\Im v\end{bmatrix}.
\]

实参数 \(x\) 的复数线性回归 \(b=Ax\) 等价于

\[
\begin{bmatrix}\Re b\\\Im b\end{bmatrix}
=\begin{bmatrix}\Re A\\\Im A\end{bmatrix}x.
\tag{16}
\]

对于当前逐元素实权 \(W_{kc}\ge0\)，加权后使用 \(\mathscr R(\sqrt W\,v)\)。它满足

\[
\mathscr R(z)^T\mathscr R(x)=\Re(z^Hx),
\]

因此实数法正规方程与复数共轭内积的实部一致。不能省略共轭或直接取复数正规方程的某一半。

如果采用完整的通道相关协方差，权重一般是 Hermitian 矩阵 \(W_k=L_k^HL_k\)。应先用 \(L_k\) 对完整复残差白化，再堆叠实部、虚部；不能把相关协方差简单替换为 16 个独立标量方差。

## 7. 本项目的目标函数与“倒数权”到底是什么

当前 FRF 拟合目标为

\[
J(a,C)=\sum_{k,c}W_{kc}
\left|\widehat H_c(\omega_k)-P_c(\omega_k;a,C)\right|^2.
\tag{17}
\]

`frf_srivc_fit.m` 中 \(W\) 是**二次代价权重**，实际回归残差乘的是 \(\sqrt W\)。这与 `raw_srivc_fit.m` 不同：后者的 `opts.weight` 直接乘在残差上，其二次代价权重是 `opts.weight.^2`。

当前重测数据采用

\[
W_{kc}\propto
\frac{\Delta f_k\,\beta_k}
{\max(|\widehat H_c(\omega_k)|,\epsilon_c)},
\tag{18}
\]

其中 \(\beta_k=0.03\) 用于 25 Hz 以下，其余为 1，\(\epsilon_c\) 是每通道幅值下限。统一归一化 \(W\) 不改变最优点。

三种选择的意义不同：

| 二次代价权 \(W\) | 代价中的单项 | 含义 |
|---|---|---|
| \(1\) | \(\lvert e\rvert^2\) | 绝对复频响误差；容易由大幅值段主导 |
| \(1/\lvert H\rvert\) | \(\lvert e\rvert^2/\lvert H\rvert\) | 当前采用的中间方案，提升小幅值段的重要性 |
| \(1/\lvert H\rvert^2\) | \(\lvert e/H\rvert^2\) | 严格的相对复频响误差；凹口附近可能过度放大噪声 |

所以“倒数权”必须说明是乘残差的权，还是乘平方误差的权。文献可能以文字将两者都称作相对权，判断时应看公式定义。

幅值下限 \(\epsilon_c\) 改变凹口附近的拟合偏好。SRIVC 线性方程中的 ridge 则用于稳定更新方向。这是两个不同的正则化位置。**Schur 消元本身不是一个权函数，也不等于在凹口处加正则。**

另一个统计边界是：如果 \(W\) 从同一份含噪 \(\widehat H\) 计算，它和噪声一般相关，不能直接套用“确定性权重＋零均值误差”的一致性证明。工程上可以使用这种权重做频段平衡，但要将该做法与统计最优性的声明分开。

## 8. 当前 SRIVC 迭代的矩阵推导

下面给出与实际程序一致的增量形式。先不写列尺度归一化，只保留数学核心。

对每个通道 \(c\)，记 \(\mathcal B_c=\mathscr R(\sqrt{W_c}\mathcal B)\)、\(E_c=\mathscr R(\sqrt{W_c}e_c)\)。对于第 \(m\) 个分母的两个系数：

\[
Z_{c,m}=\mathscr R\!\left[
-\sqrt{W_c}\frac{[q_m,q_m^2]}{A_m}P_{m,c}
\right],
\tag{19}
\]

\[
X_{c,m}=\mathscr R\!\left[
-\sqrt{W_c}\frac{[q_m,q_m^2]}{A_m}(P_{m,c}+e_c)
\right].
\tag{20}
\]

将所有分量横向拼成 \(Z_c,X_c\)，它们分别是工具变量块和测量回归块。模型对分母的导数正是 \(Z_c\)：

\[
\frac{\partial P_{m,c}}{\partial a_{m,r}}
=-\frac{q_m^r}{A_m}P_{m,c},\qquad r=1,2.
\tag{21}
\]

增量 \(\Delta a,\Delta c_c\) 通过下面的工具变量矩方程求解：

\[
\begin{aligned}
\sum_c Z_c^T X_c\Delta a
+\sum_c Z_c^T\mathcal B_c\Delta c_c
&=\sum_c Z_c^T E_c,\\
\mathcal B_c^TX_c\Delta a
+\mathcal B_c^T\mathcal B_c\Delta c_c
&=\mathcal B_c^TE_c\quad\text{对每个 }c.
\end{aligned}
\tag{22}
\]

这里是 \(Z^TX\)，并非 \(Z^TZ\)。多通道共用分母体现在所有通道对第一行方程的累加。

### 8.1 Schur／QR 消元为什么有用

对每个 \(\mathcal B_c\) 做 QR 分解，令 \(Q_c\) 列正交、张成其列空间，定义

\[
\Pi_c=I-Q_cQ_c^T.
\]

将式 (22) 的分子增量消去，得到只有 \(2n_m\) 维的分母系统：

\[
\left[\sum_c Z_c^T\Pi_cX_c\right]\Delta a
=\sum_c Z_c^T\Pi_c E_c.
\tag{23}
\]

解出分母后再恢复

\[
\Delta c_c=\mathcal B_c^\dagger(E_c-X_c\Delta a).
\tag{24}
\]

这与求解完整矩方程代数等价，前提是相应线性块满列秩；代码还通过列归一化改善尺度。它主要减少求解规模、避免不必要的大矩阵病态，不创造新的实验信息，也没有改变原目标中的频率权重。

### 8.2 稳定与下降检查

程序使用尺度化系统，加 ridge 候选求更新，再按

\[
a^+=a+\alpha\Delta a,\qquad C^+=C+\alpha\Delta C,
\quad\alpha=1,\tfrac12,\tfrac14,\ldots
\]

回溯试探。只有满足频率范围、正阻尼等约束且式 (17) 下降的候选才接受。没有找到下降方向会返回 `no_descent`，达到迭代上限会返回 `maximum_iterations`。这两种状态都不能直接称为收敛。

代码的 score 同时检查分母工具变量相关性与线性分子相关性。仅仅“相邻两次参数变化小”，在病态问题中也可能只是步长太小；因此保留 score 是必要的。

## 9. SRIVC、SK、Gauss–Newton 的关系与区别

### 9.1 与 Gauss–Newton 的区别

对于同一个实堆叠输出误差，Gauss–Newton 使用模型 Jacobian \(J=[Z,\mathcal B]\)，求解

\[
J^TJ\Delta\theta=J^TE.
\tag{25}
\]

当前 SRIVC 则使用式 (22) 中的 \(Z^TX\)。两者的差别来自

\[
X_{c,m}-Z_{c,m}
=-\mathscr R\!\left[\sqrt{W_c}\frac{[q_m,q_m^2]}{A_m}e_c\right].
\tag{26}
\]

当残差趋于零时，\(X\to Z\)，两种更新会接近；残差不小时不能把二者混为一谈。代码保存了 `X-Z gap`，正是为了检查这个区别。

可以先用 SRIVC 得到可靠初值，再用 GN／LM 优化秩一复留数、相对时延、共同支承空间等非线性结构。但后续步骤应明确叫“结构化输出误差优化”或“闭环参考矩优化”，不能全部贴上 SRIVC 标签。

### 9.2 与 SK 的关系

对单个有理模型，将工具变量替换为测量回归量，可以得到 Sanathanan–Koerner 一类迭代加权方程误差求解。MSSP 论文 Remark 3 明确讨论了 \(K=1\) 时的关系。SK 与 RIV 可以有相近的滤波形式，但它们采用的左乘矩阵不同，固定点条件也不应无条件互换。

### 9.3 固定点为什么对应输出误差驻点

固定 \(W\)、固定给定时延，式 (17) 的实梯度为

\[
\nabla_aJ=-2\sum_c Z_c^TE_c,\qquad
\nabla_{c_c}J=-2\mathcal B_c^TE_c.
\tag{27}
\]

在内部点，若 SRIVC 的完整矩相关性已归零，则式 (27) 为零。Schur 形式与此一致，因为先有 \(\mathcal B_c^TE_c=0\)，再有 \(Z_c^T\Pi_cE_c=0\)，即可推出 \(Z_c^TE_c=0\)。

这是一个有限样本的代数关系，不需要渐近极限。但它仅说明**驻点**，没有自动证明：

- 驻点一定是局部极小，而不是鞍点；
- 当前程序一定找到全局最小值；
- 数据模型正确，估计一定一致或无偏；
- 任意权重下等价于极大似然；
- 约束边界上的解需要满足无约束梯度为零。

因此“SRIVC 固定点等于输出误差驻点”可以帮助统一算法推导，但不能替代闭环数据、噪声和可辨识性条件。

## 10. 从自由加性模型到机械模态模型还差一步

若一个简单非实极点 \(\lambda_m\) 来自有限维状态空间模型，其复留数为

\[
R_m=\frac{Cv_mw_m^TB}{w_m^Tv_m},\qquad\operatorname{rank}R_m\le1.
\tag{28}
\]

此处 \(v_m,w_m^T\) 分别按右、左特征向量定义；也可以采用共轭左特征向量的等价记号，但不能在同一个公式里混用。一般阻尼下，

\[
P_m(s)=\frac{R_m}{s-\lambda_m}+\frac{\overline R_m}{s-\overline\lambda_m},
\quad R_m=c_mb_m^T.
\tag{29}
\]

合并分母后，两个实分子分别为

\[
B_{m,1}=R_m+\overline R_m=2\Re R_m,
\]

\[
B_{m,0}=-(\overline\lambda_mR_m+\lambda_m\overline R_m)
=-2\Re(\overline\lambda_mR_m).
\tag{30}
\]

即使复留数 \(R_m\) 为秩一，两个实分子也未必分别秩一。强行分别将 \(B_{m,0},B_{m,1}\) 截为秩一可能破坏合法的一般阻尼模态。

MSSP 2026 论文的关键安排是：先用自由加性 RIV 得到模型，再用参数不确定度加权的投影将其映射到模态结构。后续 GN 是该结构投影的一部分，并不意味着第一阶段也是 GN。重复极点或未分离的近邻分量可能需要更高秩的聚合留数，不能机械地对所有分量强制单秩。

约 263 Hz 的附加分量也应先按“拟合分量”看待。它可以改善特定频段的频响重建，但单凭存在一个拟合极点，无法判定它是一个真实、独立的板模态。

## 11. 通道时延如何进入本模型

考虑每个执行器和每个传感器各自的通道时延：

\[
G(s)=H_y(s)P(s)H_u(s),
\quad
H_y=\operatorname{diag}(e^{-s\tau_{y,j}}),\quad
H_u=\operatorname{diag}(e^{-s\tau_{u,k}}).
\tag{31}
\]

于是

\[
G_{jk}(s)=e^{-s\ell_{jk}}P_{jk}(s),\quad
\ell_{jk}=\tau_{y,j}+\tau_{u,k}.
\tag{32}
\]

4×4 的通道总时延矩阵不是 16 个完全独立参数，而具有行时延加列时延结构。将输入侧全部增加同一常数、输出侧全部减去该常数，不改变 \(\ell\)。必须固定一个规范，才能唯一报告相对坐标。这个规范选择不等于真实硬件某路时延必然为零。

时延给定时，先对 FRF 进行

\[
\widetilde H_{jk}(\omega)=e^{+\mathrm j\omega\ell_{jk}}\widehat H_{jk}(\omega),
\tag{33}
\]

再用上述加性 SRIVC 拟合 \(P\)。物理模型由乘回 \(e^{-\mathrm j\omega\ell}\) 得到。

时延未知时，可以外层更新相对时延，内层更新加性模型；或在 SRIVC 初值之后进行联合 GN／LM。相应 Jacobian 例如

\[
\frac{\partial G_{jk}}{\partial\tau_{y,r}}
=-\mathrm j\omega\,G_{jk}\,\mathbf1_{j=r},\quad
\frac{\partial G_{jk}}{\partial\tau_{u,r}}
=-\mathrm j\omega\,G_{jk}\,\mathbf1_{k=r}.
\tag{34}
\]

模型自由度高、带宽不足或通道还含滤波时，留数相位、滤波一阶项和纯时延可能相互吸收。要报告拟合得到的**等效或相对动态**，不能在缺乏信息时声称准确测到了全部硬件物理延迟。

事实核对：新扫频的基线 `fit_retest_position.m` 使用旧 z=0 时延矩阵并将其固定。它不是新数据上联合重估全部时延的结果。后续结构验证中的时延更新属于另一轮结构细化，不能回写成基线一直在联合估计时延。

## 12. 为什么我们当前闭环数据还不能承诺无偏

平台测得的是闭环下的总电机输入 \(U\) 和输出 \(Y\)，一般有

\[
Y=GU+E.
\]

直接谱关系给出

\[
S_{YU}S_{UU}^{-1}=G+S_{EU}S_{UU}^{-1}.
\tag{35}
\]

反馈一般导致 \(S_{EU}\ne0\)。局部有理 FRF 可以改善频率分辨率、瞬态和局部曲线平滑，但不会仅靠名称就消除这个偏差。后端 SRIVC 拟合该 FRF，也不会自动产生一个与原始噪声独立的新参考信号。

这不表示闭环数据一定完全不可用；它表示当前实现没有建立闭环一致性保证，也不能仅凭闭环二字断言偏差必然很大。正确噪声模型下的直接闭环 PEM、参考工具变量以及已知闭环模型拟合都可以是合法路线。

最适合下一轮补齐的数据是已知独立外部激励 \(d\)、总输入 \(u\)、输出 \(y\) 以及实际控制器 \(K\)。例如采用边界约定

\[
u=d-Ky,\quad y=Gu+e,
\]

可使用参考矩条件

\[
\mathbb E\{[Y-G(\theta)U]D^H\}=0,
\tag{36}
\]

其中 \(D\) 为外部激励的 DFT，且与噪声独立。也可直接拟合

\[
Y=T(\theta)D+S_y(\theta)E,
\quad T=(I+GK)^{-1}G,\quad S_y=(I+GK)^{-1}.
\tag{37}
\]

若完整控制器、滤波、增益、采样延迟和求和边界完全已知，\(d=u+Ky\) 可以按实际代数关系恢复；不能笼统说重建参考必然含额外偏差。但若记录边界不明、存在未建模限幅或缺少控制路径，则这种恢复不能凭猜测完成。

## 13. 建议学习和复现时按这个顺序看

1. 看一条清楚的共振频响，理解二阶分母决定极点、分子决定通道幅相和反共振相消。
2. 固定所有极点，仅求线性分子；理解它为什么是一次加权最小二乘。
3. 手工构造一个分量的 `part`、`res`、`x` 和 `z`，理解部分输出和工具变量差别。
4. 对比 \(Z^TX\) 与 \(Z^TZ\)，再看 Schur／QR 如何消去线性分子。
5. 学会看四个诊断：训练代价、留出误差、score、稳定／边界约束；不要只看曲线是否光滑。
6. 再引入每路时延和简单极点的复留数秩一结构，明确它们属于模型结构细化。
7. 最后用同一批留出 FRF 检查逻辑轴、双向交叉、零空间目标响应与噪声放大，讨论是否能做动态解耦。

当前训练和检查来自同一段记录的两个半段，它们可以检查拟合外推和记录内稳定性，但不等于独立重测实验，也不能直接充当估计误差置信区间。

## 14. 可核对的代码与论文出处

本地代码行号以本笔记生成时版本为准：

- `experiment_data/srivc_no_lrm_20261015/frf_srivc_fit.m:25`：W 的平方根与分母初值。
- `frf_srivc_fit.m:31`：固定极点的加权线性初始化。
- `frf_srivc_fit.m:40`：预测、复残差与目标函数。
- `frf_srivc_fit.m:50`：部分模态输出；51–54 行构造 X 与 Z。
- `frf_srivc_fit.m:56`：QR 正交投影；58 行累加 Schur 矩方程。
- `frf_srivc_fit.m:64`：同时检查分母和分子相关性的 score。
- `frf_srivc_fit.m:77`：完整矩系统与 Schur 解的独立调试比较。
- `frf_srivc_fit.m:100`：工具变量与模型导数的有限差分核对。
- `frf_srivc_fit.m:109`：ridge 候选、步长回溯与稳定性限制。
- `frf_srivc_fit.m:144`：归一化加性基函数；151 行实虚部堆叠。
- `experiment_data/srivc_no_lrm_20261015/raw_srivc_fit.m:19`：原始 DFT 路径中的幅度权 w；34–42 行是 raw U/Y 的部分输出与 Schur 形式。
- `experiment_data/srivc_no_lrm_20261015/raw_srivc_design.m:11`：输入 DFT 基函数；14、21 行是边界／干扰项。
- `sweep/srivc_retest_20261016/prepare_retest_data.m:18`：输入倍率修正；24 行划分记录；33–34 行倒数幅值二次权；60–63 行局部有理估计。
- `sweep/srivc_retest_20261016/fit_retest_position.m:7`：旧时延读取；42 行调用真正的 SRIVC；55 行固定时延状态标记。

主要原始文献：

1. R. A. González et al., *Identification of additive continuous-time systems in open and closed loop*, Automatica 173 (2025), 112013. [原文 DOI](https://doi.org/10.1016/j.automatica.2024.112013)。重点读第 3 节的开闭环最优性与工具变量条件、第 4 节迭代。项目文本：`research/identification_innovation_20261017/additive_ct_text.txt`。
2. M. van der Hulst et al., *Structured identification of multivariable modal systems*, Mechanical Systems and Signal Processing 247 (2026), 113948. [原文 DOI](https://doi.org/10.1016/j.ymssp.2026.113948)。重点读式 (23) 固定极点初始化、式 (29)–(38) 的加性 RIV、Remark 3 的 SK 关系、第 5 节结构投影；式 (52) 是其权重定义。项目文本：`research/identification_innovation_20261017/mssp_structured_text.txt`。

阅读文献时也应检查每个结论的条件。例如，一阶驻点条件本身不足以证明局部最小；工程权重不自动等于噪声协方差逆；原论文证明的特定开闭环设定不能不加说明地转移到当前 U/Y 入口。
