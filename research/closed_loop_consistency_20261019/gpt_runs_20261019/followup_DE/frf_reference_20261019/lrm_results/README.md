# 参考驱动局部有理 FRF 诊断

这是一条新增诊断路线，不替换预先指定的 LPM 主分析，不读取或拟合任何 SRIVC 参数、模态极点或时延。输入来自旧 `real_d/data/z{z}_P08_s01.npz` 的逐周期谱，四个位置为 z = −0.5/0/0.5/1 mm。

## 方法与固定选择规则

频率为每 6 s 周期的激励谱线，3–25 Hz、133 点。对八个输出 `[U−D, Y]` 分别拟合二阶局部有理模型。分子含四输入的二次多项式，共 12 个复系数；瞬态含 0/1/2 次项，共 3 个；分母常数固定为 1，另含 1/2 次项，共 17 个复系数。第四反馈输出数值精确为零，跳过该回归，重构时使用 `TU = I + TF`。最后 `G = TY / TU`，矩阵维数为频点×输出×输入。

六个候选为线性 LRM 或 SK3，局部窗 21/25/31 点。SK3 是初始线性解后，再进行三次冻结分母的 `1/abs(denominator)` 加权；数值分母下限仅为机器精度。未按与参数化模型的接近程度选择任何设置。

选择只用每位置前四周期。所有候选在共同中心集合（去除31点窗边界的首末15点）进行留一激励频线预测；中心频线在所有训练周期中都从回归中去除。八个输出分别以训练周期差分的全频带噪声均值标准误作固定缩放，第四反馈输出排除。缩放对全部候选相同，不使用参数模型。将四位置标准化预测 SSE 相加，除以共同预测项数，选择一个全局配置。

锁定结果为 **线性 LRM、31 点窗、0 次 SK 更新**。31点谱线窗跨度为5 Hz。锁定后，全数据、前半、后半均使用相同配置。分母拟合存在含噪回归项，本诊断并不声称其有限样本无偏或最优。

## 文件

- `configuration.json`：全局锁定设置与六组总 CV 分数。
- `z{z}_front_cv.json`：逐位置六配置CV、共同中心、固定尺度与失败计数。
- `z{z}_{full,front,back}_lrm.npz`：`f,G,TU,TY,se_db,se_phase_deg,se_complex`，另存17个复系数、局部回归秩与条件数、分母最小模、TU条件数、失败旗标、边界窗、预测、逐周期删除后的FRF。
- `z{z}_{split}_summary.json`：各记录摘要。
- `algebra_verification.json`：独立生成局部有理 MIMO 数据与瞬态，检查矩阵方向、17参数回归与留出频线预测；最大FRF绝对误差约4.1e−17。
- `run_all.log`：完整运行日志。
- `pilot_configuration.json`：仅 z=0 的初步结果，标记 provisional；正式结果以 `configuration.json` 为准。

## 不确定度与解释边界

对每个输出窗的 P 个周期，删除一个周期重新拟合并重构 FRF。复数重复性量为 `sqrt((P−1)/P * sum(abs(G_leave_one−mean(G_leave_one))**2))`；幅值和相位同理，分别使用 dB 与相对完整估计的主值相位。该 jackknife 量描述周期间随机重复性；相邻周期如果相关，其标准误解释还需额外条件。它不包含窗宽/局部模型偏差、激励恢复误差或硬件模型失配，不能视为“总误差棒”。

所有频点，包括边界窗、高条件数和很大的误差棒，均保留。全P8重构的最大幅值 jackknife 量约30.7 dB；不能仅用中位数掩盖局部弱响应。全P8模型对同P8 FRF属于同数据核查；前半模型对后半FRF属于分段验证。全局配置已用所有位置前半数据选择，因此后半仍是未用于配置选择的记录段。

## 复现

使用 Python 3.11 + NumPy，脚本自动限制BLAS为单线程，禁止写旧目录pycache。

```powershell
& 'C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe' 'C:/Users/Administrator/Desktop/my_article/identifiation1/research/closed_loop_consistency_20261019/gpt_runs_20261019/followup_DE/frf_reference_20261019/lrm_diagnostic.py' --all
```

`--pilot` 只计算 z=0 的六组CV及其暂优配置的完整记录。正式复现须运行 `--all`，冻结四位置全局配置。
