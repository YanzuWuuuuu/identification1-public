# 旧实测数据：非参数 FRF 参照

本目录只新增旧实验数据分析。主报告：[旧数据非参数FRF参照报告.md](旧数据非参数FRF参照报告.md)。

数据：z=-0.5/0/0.5/1 mm，5 kHz，6 s周期，第2–9周期，3–25 Hz。复用 `../real_d/data/` 的周期谱及既有 M1/M2/M3/M3w/M5/M5n 参数。全部16通道，无筛点。没有修改旧数据、旧拟合或理论文档。

## 运行

在本目录运行，Python 3.11 + numpy/scipy/matplotlib：

```powershell
$pythonExe = 'C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe'
& $pythonExe run_frf.py --z all
& $pythonExe lrm_diagnostic.py --all
& $pythonExe lrm_diagnostic.py --diagnostics
& $pythonExe verify_lpm.py
& $pythonExe summarize_frf.py
& $pythonExe plot_frf.py --z all --panels all --reference lpm
& $pythonExe plot_frf.py --z all --panels main --reference lrm
```

`run_frf.py` 把旧参数模型重新求值，不重新辨识；缺失的后半段M1保持不可用。LPM预定二阶21点，另有一阶15点、二阶31点敏感性检查。LRM是后续诊断；其窗口由前半数据留一频线预测选择，所选线性31点结果冻结到 `lrm_results/configuration.json`，不根据参数模型挑选。

## 内容

- `PROTOCOL.md`：范围、预定配置和后续诊断说明。
- `data_audit.md`、`method_audit.md`：数据边界、周期结构和独立方法检查。
- `frf_core.py`、`run_frf.py`：LPM及模型求值。
- `lrm_diagnostic.py`：LRM、训练段交叉验证、周期jackknife和频带/列分解。
- `verify_lpm.py`：无噪声代数核验、M5内部/外部模型与DFT尺度检查；**不是新的参数估计**。
- `summarize_frf.py`：报告、全指标CSV和输入哈希。
- `comparison_metrics.csv`：224组完整参照差、无逆残差和原求解状态，包括失败端点。
- `model_parameters_and_status.csv`：既有三支承模态参数、迭代/评估次数、原score与成功标记。
- `reference_repeatability.csv`：LPM/LRM各位置参照重复性。
- `results/`：LPM及参数模型FRF；`lrm_results/`：LRM、配置选择和诊断。
- `figures/`：LPM 48图、LRM 24图，均有PNG/PDF/SVG，源数值CSV和逐图manifest。
- `source_manifest.json`：实际读取的谱缓存及原拟合结果SHA256，不代表拥有硬件真值。

主预览：`figures/lrm_z0_full_vs_full_main_magnitude.png`；连续相位为同前缀 `_phase.png`；分段检查为 `lrm_z0_front_vs_back_main_magnitude.png`。LPM较差结果仍在不带 `lrm_` 的对应文件中。

## 解读边界

灰带是±2随机SE，**不是包含偏差的95%真值置信带**。LRM前后半机械FRF差约41%–58%，主要在8–18 Hz及第二输入列。不能据模型与此参照的接近程度给准确度排名。

M5在z=0、1虽旧优化器标记收敛，输出和FRF验证异常；FFT尺度与瞬态换算已排除，结果保留。所有方法名称仍对应原定义：M3为followup D的**谱线参考IV**。

复现过程只写本目录；各脚本禁用外部导入的字节码缓存，并限定数值库单线程。完整图批处理约几分钟，局部频响计算约数秒至数十秒。
