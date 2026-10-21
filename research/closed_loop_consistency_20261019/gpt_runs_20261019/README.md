# 频域闭环 IV / 周期 SML / 间接 FRF 对照实验

主报告：[三组对照实验_完整报告.md](三组对照实验_完整报告.md)。所有本轮新输出集中在此目录；原理论、原代码、原实测记录、已有结果未修改。

## 完成范围

- A：5种周期数，180组配对仿真，M1–M4共720次拟合；SML140次主拟合（P2/4明确为收缩诊断）。
- A′：100组四次单轴实验；规范IV–GN160次主拟合，另160次同目标信赖域优化器对照。
- B：完整八模态，20组×3方法，60个最终拟合。
- C：96个CT/ZOH主拟合、4个第二半段参考估计、20个SML/不可用记录；4位置全覆盖。
- 保存完整参数、初值、随机种子、逐次状态、异常、迭代/评价次数和耗时；每完成一次即落盘。

周期SML和完整FRF协方差存在任务预算下的秩限制。原始不可用项与收缩版本分开列出，没有使用静默伪逆，也没有把不可用写成成功。

## 环境与运行

本次Python3.11.9、NumPy2.4.6、SciPy1.17.1；图使用Matplotlib。所有入口见本目录 `run_all.ps1`。在本workspace根目录执行：

```powershell
& 'research/closed_loop_consistency_20261019/gpt_runs_20261019/run_all.ps1'
```

运行器识别已有逐次检查点并跳过，不覆盖原始实验结果。上述命令可恢复/汇总本次运行。若要真正独立复跑，应创建新运行目录并保留原目录结构，不能删除本轮结果来伪造首次运行。拟合程序从原 `code/cl2.py` 只读导入；`baseline_snapshot/` 保存了本次源代码和文档副本，`source_inventory.json` 保存原件哈希。真实输入数据不复制，仍只读本机路径。

原始环境烟雾自检命令（仅复跑时将输出名设为新的文件，保留本轮chk.npz）：

```powershell
$env:NOFLEX='1'; $env:FC='60'; $env:AL='30'; $env:AH='100'
$env:FHI='100'; $env:FLO='3'; $env:SW='0.03'
$env:FSPR='0.01'; $env:ZSPR='0.4'; $env:NN='300000'; $env:NPER='10'
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:PYTHONDONTWRITEBYTECODE='1'
python research/closed_loop_consistency_20261019/code/run_mc5.py 0.2 3e-7 3 1002 NEW_SMOKE_OUTPUT.npz dir_band_W1,ref_band_W1
```

主运行的数值配置、种子及数据预算在 `PLAN.md` 和各子目录README；采用多进程时每个BLAS单线程。本轮A使用3个工作进程，B使用2个；SML/间接法分别使用不重叠的worker分片。子目录日志及单次JSON记录比控制台进度更详细。

## 数据和程序索引

|路径|内容|
|---|---|
|experiment_a/run_a.py|原IV模型的配对周期扫描与DFT缓存|
|experiment_a/results/|720个完整参数JSON；summary.csv/json、all_endpoints.npz为汇总|
|experiment_a/data/|180份共享DFT与逐周期DFT；SML完全复用|
|comparators/sml_core.py|参数相关协方差的完整SML残差及Jacobian|
|comparators/run_indirect_ivgn.py|用户所需梯度工具变量FRF拟合的无瞬态适配|
|comparators/README.md|SML文献核对、协方差秩、预算与算法区别|
|experiment_b/|完整模型、等价稀疏计算验证、20次主运行及全部异常日志|
|experiment_c/|原实测拟合复制件、ZOH扩展、120个结果、误差棒及完整表|
|experiment_c/diagnostics/|A少周期异常及对照方法的独立数值复核|
|plot_results.py|从保存的统计绘图|
|build_report.py|验证数量与源文件完整性，并生成主报告|

报告统计使用物理相对百分比 `100*expm1(delta_log)`。实测阻尼原值为百分数，旧/ZOH差以阻尼百分比点报告；仿真方法间差是相对误差的百分比点，两者不能混用。主表包含所有有限终值，同时报告失败；不得把删掉失败后的条件统计当作一般性能。

本轮结论：参考IV有减偏作用，少周期保证和全面优于SML的强叙事不受支持。详细数据与限制以主报告为准。
