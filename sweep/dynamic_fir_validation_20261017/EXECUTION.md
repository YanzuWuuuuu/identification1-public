# 执行与验收记录

2026-10-17，本目录完成辅助开发验证，未改硬件和既有模型。

1. MATLAB只读导出5个位置的既有训练/检查FRF、SRIVC频响与坐标为本目录v7 MAT。
2. 初始pilot发现增益稀网格漏掉约0.00026的峰值；pilot另存，不进入正式筛选。正式网格上限1.0495，验收上限1.05。
3. 主程序 `python -u run_validation.py` 运行全部12个候选。三个parallel脚本执行互不重叠的SRIVC候选，主程序复用完成文件；没有依据检查集改变任何候选配置。
4. `selected_training_only.json` 在读取Hcheck之前生成。后半记录在过去开发中曾查看，不能称独立盲测；主程序docstring里的“Independent half-record”只对应分段，不能解释为统计独立。
5. `python summarize_results.py`、`python plot_results.py` 生成数值汇总与6张PNG/PDF文件。
6. `post_selection_repeatability.py` 在确认选择冻结后执行，并核对冻结JSON哈希不变。
7. MATLAB执行 `addpath('sweep/dynamic_fir_validation_20261017'); verify_fir_export();`。首次发现MIMO tf初始化不接受数值列向量，改为分子/分母cell；第二次发现MAT导出的延迟为整数类型，转为double。两项只影响导出接口，不改变优化系数和候选选择。
8. 最终MATLAB核验通过，最大误差4.54e-15：滤波器频响、纯延迟、完整YPU乘积、filter/lsim时域等价及FIR稳定性。
9. 三张PNG已目视核对。连续相位为原复频响独立展开，弱信号噪声可能导致多圈绕转；不能用圈数直接估时延。

冻结结论：非参数FIR路线有有限串扰收益，保护条件没有全部迁移到检查数据；SRIVC模型路线回退静态。没有为了通过验收放宽门槛或在检查集上重选参数。输出系数用于复现离线分析和Simulink探索，不等于已经通过硬件反馈验证。

正式优化配置及其执行代码哈希在manifest.json。所有模型/数据性能数值保存在CSV/NPZ中，说明性代码调整没有改写正式run_validation.py。
