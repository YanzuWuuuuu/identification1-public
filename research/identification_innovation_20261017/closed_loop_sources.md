# 闭环 SRIVC 与结构模态辨识：一手文献和可执行修复

核查日期：2026-10-17。只做资料与数学审计，没有运行 MATLAB、修改实验数据或改写论文。核对了用户提供的三篇本地 PDF、2025 年两份作者预印本及出版社/作者网站。本记录把文献结论与针对本项目的推论分开。

## 1. 最重要的判断

当前流程确实执行了频域加性 SRIVC，但“算法采用工具变量”不等于“当前闭环实验已无偏”。本轮非参数入口只包含闭环总电机指令和传感器输出，没有单独外部激励，训练权又由同一训练 FRF 的模倒数构造。数值 score 达标只能证明相应有限数据迭代收敛，不能独立排除第一阶段 FRF 偏差、错误输入边界、非线性和弱辨识。

用户提供的 MSSP 2026 论文与本项目的根本差别不仅是画图：它通过闭环间接法取得 FRF，采用 13 个单轴随机相位多正弦、10 周期、10 次实现；随后才执行加性模型 RIV 与结构投影。必须先建立可审查的数据入口，再讨论新算法收益。

## 2. 八篇直接相关一手来源

| 来源 | 已核实结论与对本项目的影响 | 文献状态与链接 |
|---|---|---|
| van der Hulst, González, Classens, Tacx, Dirkx, van de Wijdeven, Oomen, **Structured identification of multivariable modal systems** | 频域加性 RIV 后经 IPEM 投影至模态结构；兼容一般/比例阻尼。实验使用闭环间接 FRF，重复多正弦及 CMIF 选阶。不能把加性分解、秩约束或两阶段投影单独当作我们的新贡献。 | *MSSP* 247 (2026), 113948；[DOI](https://doi.org/10.1016/j.ymssp.2026.113948)；核对用户本地出版社 PDF 的第 6.2–6.5 节。 |
| González, van der Hulst, Classens, Oomen, **Statistically Optimal Structured Additive MIMO Continuous-time System Identification** | 已有时域结构加性 MIMO、闭环 IV、一致性、协方差最优投影。闭环数据明确要求 r,u,y；仪器由参考通过模型闭环灵敏度生成。其“最优”在闭环是 IV 类内最小协方差，通常不等于全 ML 的 CRLB。 | 2025-05-20，作者网站仍列 working manuscript，勿写成已发表 Automatica；[arXiv:2505.14169](https://arxiv.org/abs/2505.14169)。 |
| González, Pan, Rojas, Welsh, **Consistency analysis of refined instrumental variable methods for continuous-time system identification in closed-loop** | 证明连续控制器闭环与数字控制器闭环不能混为一谈。数字控制器场景的普通 SRIVC 一致性仍有白噪声及严格因果等限制；CLSRIVC 需正确构造参考仪器。它不证明所有含闭环数据的 SRIVC 均有偏，也不证明用了 SRIVC 就无偏。 | *Automatica* 166 (2024), 111697，2024-08；[DOI](https://doi.org/10.1016/j.automatica.2024.111697)，本地 PDF 定理 11、推论 13、定理 15。 |
| González, Classens, Rojas, Welsh, Oomen, **Identification of additive continuous-time systems in open and closed loop** | 加性 SISO 的开闭环统一 IV 迭代已有成熟推导。闭环参考与扰动独立、仪器与噪声不相关、交叉矩阵非奇异、正确保持器/采样和闭环稳定是关键条件。该文不能替本项目未知控制边界或相对动态模型失配兜底。 | *Automatica* 173 (2025), 112013，2025-03；DOI 年份为 2024；[DOI](https://doi.org/10.1016/j.automatica.2024.112013)。 |
| González, Classens, Rojas, Oomen, Hjalmarsson, **Finite Sample MIMO System Identification with Multisine Excitation: Nonparametric, Direct, and Two-step Parametric Estimators** | 已讨论任意 Schur 权影响统计性质，并在明确模型下证明 FRF 充分统计量和最优两步/时域估计等价。设定为已知连续多正弦、稳态、输出白噪声及多组线性独立输入，不能直接迁移到本项目闭环总 u 的 LRM 结果。 | 2025-10-30 作者预印本；[arXiv:2510.26929](https://arxiv.org/abs/2510.26929)。 |
| Voorhoeve, van der Maas, Oomen, **Non-parametric identification of multivariable systems: A local rational modeling approach with application to a vibration isolation benchmark** | MIMO LRM 的参数化、窗口、方向依赖及线性化误差都不是小细节。非迭代 Levy 式有额外分母加权；文中给出迭代方法。LRM 可以改善泄漏/轻阻尼共振估计，不能自动消除任何反馈相关噪声。 | *MSSP* 105 (2018),129–152，2018-05-15；[DOI](https://doi.org/10.1016/j.ymssp.2017.11.044)。 |
| Dirkx, van de Wijdeven, Oomen, **Frequency Response Function identification for multivariable motion control: Optimal experiment design with element-wise constraints** | 已针对 MIMO 方向性、输入输出逐元素功率约束进行 FRF 激励优化，并在 7×8 晶圆台验证。不能把普通多输入最优激励独立作为本研究创新；可以借其方法专门优化“校准目标”的信息。 | *Mechatronics* 71 (2020),102440，2020-11；[DOI](https://doi.org/10.1016/j.mechatronics.2020.102440)。 |
| Gilson, Garnier, Young, Van den Hof, **Optimal instrumental variable method for closed-loop identification** | 仪器、预滤波器与噪声模型决定闭环 IV 的偏差和方差；已区分控制器已知/未知情形。是参考矩路线的基础，不是本文新发明。 | *IET Control Theory & Applications* 5(10) (2011),1147–1154；[DOI](https://doi.org/10.1049/iet-cta.2009.0476)；[作者保存全文](https://publications.pvandenhof.nl/Paperfiles/Gilson%26etal_IET2011.pdf)。 |

补充边界：Classens, González, Oomen 的 **Recursive identification of structured systems: An instrumental-variable approach applied to mechanical systems**（*European Journal of Control* 84 (2025),101238，[DOI](https://doi.org/10.1016/j.ejcon.2025.101238)）已做递归加性 SRIVC 及随时间变化的机械共振。因此仅将五个 z 位置扩成在线跟踪也不够新。

## 3. 闭环偏差的最小数学说明（本项目推论）

在正确记录边界设

\[
y=Gu+v,\qquad u=d-Ky,
\]

其中 d 是外部激励，K 包含已实现的全部数字控制器、几何映射和已知反馈滤波。定义互谱 S_ab=E[AB^H]。

直接用输入输出谱矩估计有

\[
S_{yu}S_{uu}^{-1}=G+S_{vu}S_{uu}^{-1}.
\]

反馈通常使 S_vu 非零。窗口平均或局部有理逼近不能在没有额外条件时使该项自动消失。这里不是断言本次每个频点都有显著偏差；它指出当前尚缺怎样的排除证据。

若外部 d 与 v 独立，则

\[
S_{yd}=G S_{ud},\qquad
m_k(\theta)=\operatorname{vec}\{\widehat S_{yd,k}-G_k(\theta)\widehat S_{ud,k}\}.
\]

可先用外部参考估计 FRF 再 SRIVC，也可直接拟合 m_k，避免显式求逆病态 \(\widehat S_{ud}\)。但矩条件本身是标准 IV/GMM，不能因此声称新方法。有限样本 \(\widehat S_{yd}\widehat S_{ud}^{-1}\) 是随机量比值；一致性不等于有限样本精确无偏。

写成实虚堆叠后，合理的二阶段权为矩误差协方差的逆；在 model-dependent covariance 下应使用可行 GMM/固定每轮权并分析其极限，而不是无说明地把任意权更新命名为 ML。局部 LRM 的相邻频点共享样本，因此仅取各频点对角方差可能低估不确定性。重复实现估计与按完整实验分组的 bootstrap 更可信，不能把密集插值频点当独立样本数。

## 4. 没记录 d 是否还可挽救：可以，但条件严格

如果 u=d−K_eff y 正是硬件运行时的数字恒等式，并且记录了同一 u 和同一反馈 y，则

\[
\widehat d=u+K_{eff}y=d.
\]

这里即使反馈 y 含测量噪声，噪声也随真实控制律代数抵消；不能笼统说“重构参考天然有噪声”。必须满足：

- K_eff 是完整实际离散算法，包含一拍延迟、滤波、单位和静态几何矩阵；不是一个近似连续 PID。
- 输入在与辨识对象一致的加法节点记录；本次 7.74 倍校正正说明该检查有实际必要。
- y 日志与控制器使用的值相同且同步；若另经滤波、重采样或异步记录，要把差异建模。
- 饱和、限速、死区、抗积分饱和、前馈、重力补偿等需显式复现，或限于未触发非线性/已扣除已知项的数据段。
- 对有状态控制器，初始状态须可重放或排除其过渡段；跨片段直接 FFT 乘 K 还需处理有限记录边界/瞬态。

如果 \(u_{log}=u+\eta_u\)、\(y_{log}=y+\eta_y\)、\(\widehat K=K+\Delta K\)，则

\[
\widehat d-d=\eta_u+K\eta_y+\Delta K\,y+\Delta K\,\eta_y.
\]

此时以恢复信号作 IV 的有效性不能自动成立，应检验可重放的实际控制律并量化该误差。若参考只加在三根逻辑轴上，恢复得到的四维物理激励可能只有三个方向，仍不一定足以单独估计第四个冗余输入方向。

无法完整恢复时：现有数据仍可探索共振、参数变化和 FRF；可尝试带可靠噪声结构的闭环 PEM 作为敏感性对照，但“只有 u,y”在一般闭环系统下并不等价于已经有足够信息证明独立机械参数的一致性。最简可靠修复仍是下次同步记录外加 d。

## 5. 激励的秩：本地回归满秩不等于校准目标信息充分

每个频点单组 4×1 输入和 4×1 输出不能自由确定 4×4 FRF；需要多次独立方向实验、随机输入谱矩满秩，或跨频率的可信参数化。LRM 通过跨邻近频点借用平滑结构，可能产生满秩局部回归，但它不直接证明冗余方向受到足够外部激励。

应检查：外部参考谱/跨实验输入矩阵的最小奇异值，条件数随频率的变化，以及与噪声水平相比的有效第四方向。控制噪声可能让总输入谱矩看起来数值满秩，不等于外部参考提供了第四方向信息。

对现平台推荐四组独立物理激励方向；也可用三个刚体逻辑方向加一个名义零空间方向组成满秩基。必须保留所有四路原始输入输出。仅在名义零空间扫频会漏掉几何失配对其他方向的证据。

## 6. 对本项目创新的约束与机会

### 已有，不能再单独作为贡献

- 频域加性 RIV/SRIVC；低阶结构投影；模态残差秩约束。
- 闭环参考 IV、一般一致性与渐近协方差。
- 用完整协方差替换任意 Schur 权；两步 FRF 与参数辨识。
- 普通递归跟踪、普通多输入最优激励。

### 值得聚焦但尚待证明的增量

**把辨识精度的目标从主导物理 FRF 改为冗余系统的解耦交叉块，并给出可验证的误差和性能联系。** 例如同一整体物理 FRF 误差下，是否能更准确预测/降低逻辑刚体到冗余通道、冗余输入到刚体通道的泄漏，同时保留目标柔性响应。通道相对动态作为共同模型结构的一部分，避免把弱可辨识的独立纯时延/滤波时间常数当最终工程成果。

令目标响应为

\[
\psi(\theta)=\operatorname{stack}_{\omega\in\Omega}
\{L_rG(\theta)R_n,\ L_nG(\theta)R_r,\ L_nG(\theta)R_n\},
\]

可用 \(J_\psi P_\theta J_\psi^T\) 判断哪些参数弱方向真正损害目标，而不是只报原始参数 Hessian 条件数。若校准器本身依赖 θ，J_ψ 必须包含校准器的导数；不能将其当固定矩阵而漏项。该式是标准误差传播，潜在创新来自可实现的结构约束、目标可辨识性/误差界和工程实证的组合，而不是这一个公式。

当全模型正确、采用有效 GLS/IV 时，不应随意加目标权就承诺方差更小；目标权的价值主要在有限数据、有限阶模型失配和激励预算下。必须与已有最优结构 IV 在相同数据、阶数、约束、计算预算下比较，证明目标泄漏确有改善且不是牺牲主通道或把噪声拟合成“零空间”。

最低证明/验证要求：

1. 校准目标在通道公共动态/规范自由度下是否唯一或至少等价；
2. 从可辨识的外部参考矩到目标 ψ 的估计误差条件；
3. 一般刚体输入造成的一阶串扰必须纳入，不能用仅纯冗余激励成立的二阶界替代；
4. 独立重复实验上报告目标交叉块绝对值、置信区间、柔性保留率和噪声增益；
5. 在已知数字控制器和完整原始日志下验证数据入口，再去谈时延带来的校准收益。

这条路线可作为待验证研究假设；本轮文献核查不足以宣称“从未有人做过”或“已构成可发表新理论”。
