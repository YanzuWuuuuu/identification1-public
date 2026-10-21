## 11. 从本地论文可借鉴什么，以及本工作的边界

本轮已读用户 `article` 内相关全文，关键公式页已核对。以下编号与第 0–8 节对应，页码为 PDF 顺序页码。

| 编号及文献 | 核对位置 | 本版借鉴与限制 |
|---|---|---|
| R1：González, Pan, Rojas, Welsh, *Consistency analysis of refined instrumental variable methods for continuous-time system identification in closed-loop*, Automatica 166 (2024), 111697；[DOI](https://doi.org/10.1016/j.automatica.2024.111697) | 第 3–4 页 Lemma 3；第 7 页 Theorem 11；第 8–9 页 Theorem 15 | 解析域内的通用非奇异性、直接法一致性的条件、参考工具与正确采样间行为。其定理也以迭代已收敛等为前提，不能用来证明本算法全局收敛。 |
| R2：González, Classens, Rojas, Welsh, Oomen, *Identification of additive continuous-time systems in open and closed loop*, Automatica 173 (2025), 112013；[DOI](https://doi.org/10.1016/j.automatica.2024.112013) | 第 4 页 Lemma 1、Remark 3；第 4–5 页 Theorem 1；第 7 页 Theorem 2 | 使用工具—回归交叉矩阵、Taylor/CLT 与 sandwich；加性 SISO 理论不能不经推导直接改称 MIMO 定理。 |
| R3：Gilson, Garnier, Young, Van den Hof, *Optimal instrumental variable method for closed-loop identification*, IET Control Theory & Applications 5(10) (2011), 1147–1154；[DOI](https://doi.org/10.1049/iet-cta.2009.0476) | 第 6 页 §3.1–3.2；第 9–10 页 §4.2 | 一致性需要外生性及交叉矩阵满列秩；未知控制器可用外部参考辅助模型。最优性限于规定工具类及噪声条件。 |
| R4：van der Hulst et al., *Identification of additive multivariable continuous-time systems*, 本地 LCSS 已接受预印本，arXiv:2504.01639v2 (2025) | 第 3 页式 (11)、(14)；第 4 页式 (26) | 已有加性 MIMO、闭环参考 IV、完整输出协方差与残差更新。不能把这些基本要素整体称为本文首次提出。固定邻点谱权与其全记录池化协方差不同。 |
| R5：Blom, Van den Hof, *Multivariable Frequency Domain Identification using IV-based Linear Regression*, CDC 2010, 1148–1153 | 第 2–4 页 Propositions 4.1–4.4，式 (12)–(19) | 借鉴频域 IV 的驻点恒等式和加权形式；其 FRF 拟合结论不自动覆盖原始闭环 DFT 的噪声内生性。 |
| R6：Pintelon, Schoukens, Vandersteen, *Frequency Domain System Identification Using Arbitrary Signals*, IEEE TAC 42(12) (1997), 1717–1720 | 第 2 页式 (10)–(13) | 有限记录 DFT 的瞬态分子及联合估计已存在；本版细化逐块瞬态和两种渐近尺度。 |
| R7：Pintelon, Schoukens, *Box–Jenkins identification revisited—Part I: Theory*；以及 *Frequency domain maximum likelihood estimation of linear dynamic errors-in-variables models* | 本地全文及先前精读；本轮不将其所有估计器重新证明 | 频带选择、有限记录和 EIV/SML 都有成熟先例。不能声称 SML 一定需要正确控制器，也不能把未知控制器处理本身当作参考 IV 的独有优势。 |

本轮逐页记录：[literature_iv.md](C:/Users/Administrator/Desktop/my_article/identifiation1/research/closed_loop_consistency_20261019/gpt_runs_20261019/followup_DE/proof_repair_v4/literature_iv.md)。固定邻点随机权的完整证明：[random_weights_proof.md](C:/Users/Administrator/Desktop/my_article/identifiation1/research/closed_loop_consistency_20261019/gpt_runs_20261019/followup_DE/proof_repair_v4/random_weights_proof.md)。

**可讨论的论文内容**是：指定频带闭环频域加性模型的可计算偏差、无需 FRF 中间步的参考矩实现、正确 ZOH 与干扰瞬态处理、周期路线区分、以及实际随机加权的一致性条件。基本参考 IV 原理已有先例，固定邻点权证明的文献新颖性仍应作针对性检索，不能仅凭本地这批论文未写就宣称首创。

**不能主张：** 所有闭环直接法不一致；错误控制器只影响方差；仅少周期可用参考 IV；M3w 已达到最优效率；参考 IV 一定比 SML 准；短记录无偏或迭代全局收敛；由此解决零空间解耦。本次不扩展到解耦、柔性模态或间接 FRF。

## 12. 当前证明完成度与数值证据的关系

| 结论 | 本版状态 | 仍需满足的实际条件 |
|---|---|---|
| 确定性权的参考 IV 局部一致性 | 第 7.1 节给出交叉矩阵与 Gaussian 足够条件下的证明；路线 F 可用有限频点小噪声论证 | 正确模型、联合激励、局部根、工具滤波器有界 |
| 错误控制器 | 保留真根；局部一致性是条件性结论 | 交叉矩阵非退化，不能只检查工具 Gram |
| 固定周期、8/16 邻点估计权 | 第 7.3 节已补证明，允许同记录权重—残差相关；一或两次更新均覆盖 | 初步联合局部根、正谱、实际设计满秩、求解误差随样本缩小 |
| 固定少周期、8/16 邻点估计权 | 第 7.4 节在明确 Gaussian、平滑、局部正则条件下补证明 | 不能外推至所有平稳非 Gaussian 噪声或模型失配 |
| 确定性/oracle 权的结构参数协方差 | 第 7.2 节修正 cross/sandwich，并保留瞬态 Schur 补 | 足够精确的局部根、score CLT、相应设计极限 |
| M3w 的一般渐近正态及最优效率 | **未证明，也不作为本版结论** | 固定邻点权的分布/代入一阶影响需要另分析 |
| 当前有限预算程序必然一致 | **不能这样表述** | 固定极点步长阈值不等于完整 IV score 趋零，失败必须保留 |
| 实测谁更准确 | **没有真值，不能判定** | 只能报告方法差异、重复性、收敛与模型适切性 |

后续 D/E 报告已提供数值反证边界：预设规则 A、B 均未整体成立；P=2 的主要 M3w 设置三个阻尼 RMSE 均超过 M5n 的 1.1 倍，其他周期的相对表现随模态变化。因此即使一致性证明修好，也**不产生 SML 精度优势**。P=1 和 LM 失败问题属于有限样本与求解层面的限制，继续如实报告。

完整后续实验见：[后续实验D_E_报告.md](C:/Users/Administrator/Desktop/my_article/identifiation1/research/closed_loop_consistency_20261019/gpt_runs_20261019/followup_DE/后续实验D_E_报告.md)。本轮没有重新运行该蒙特卡洛，也未覆盖其任何结果文件。

## 13. v4 修订清单

| v3 的缺口/过强表述 | v4 修复 |
|---|---|
| 任意错误控制器使用同一 Gram 矩阵 | 改为 \(Z^H\mathcal WJ^d\) 交叉矩阵，正确控制器才化为 Gram |
| 错误控制器只影响效率 | 增加工具有界与交叉秩条件，指出弱工具/失秩可能 |
| 一切积分都按激励谱测度处理 | 激励测度与噪声频点计数测度分开 |
| 有界谱直接推出 a.s. 一致 LLN | 给出 Gaussian 创新下的方差、包络、有限覆盖论证；只声称依概率收敛 |
| 全单位圆排除极点，却又声称含刚体积分 | 只在固定辨识频带内要求有界，闭环边界另作假设 |
| 瞬态维数固定所以可直接删掉 | 区分物理边界与自由估计；两路线分别缩放；保留增广秩及 Schur 协方差 |
| 只写数据权重“渐近无关” | 给出固定邻点随机权的两条局部一致性证明，保留非消失共同瞬态平移 |
| 冻结输入等于求矩导数时也不用管工具变化 | 区分构造梯度与整体矩求导，说明零期望项的来源 |
| 离散极点与连续参数全局双射、逐频列空间相同足够 | 改为局部可逆、频率共用坐标变换的明确接口条件 |
| 固定周期重复与固定重复数延长混为一谈 | 单列路线 F/D，重写 \(1/P\) 推论与 P=1 异常解释 |
| LRM 等同全部频点直接法；SML 必须知道正确控制器 | 删除，改为按所估统计关系和噪声结构判断 |
| 有限样本几乎零偏等于证明成立 | 将局部统计定理、求解精度与数值证据分开 |

原 v3 完整留存于 [原稿_v3_修订前.md](C:/Users/Administrator/Desktop/my_article/identifiation1/research/closed_loop_consistency_20261019/gpt_runs_20261019/followup_DE/proof_repair_v4/原稿_v3_修订前.md)。旧版第 13 节的历史审查记录也保留在该备份；本版上表取代它作为当前状态。
