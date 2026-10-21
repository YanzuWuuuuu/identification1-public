# P6 重参数化 / P7 控制器无关 jackknife 去偏

所有输入只读。新代码、端点、图表、日志和报告均在本目录。判据在开始前保存于 `PROTOCOL.md`；输入哈希在 `baseline/input_hashes.json`。

## 当前执行方式

Python：

```powershell
$py = 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe'
$run = 'C:\Users\Administrator\Desktop\my_article\identifiation1\research\closed_loop_consistency_20261019\gpt_runs_20261020\newdata\reparam_tests'
```

各脚本自动跳过已保存的端点。无需删除结果来恢复中断。

```powershell
& $py "$run\code\verify_profile.py"
& $py "$run\code\run_p6.py" --workers 4
& $py "$run\code\metrics.py"
& $py "$run\code\run_p7_real.py"
& $py "$run\code\run_p7_sim.py" --workers 4 --periods 4,8,16 --count 40
& $py "$run\code\run_p6.py" --jackknife --methods M3 --workers 4
& $py "$run\code\run_fixed_low.py"
& $py "$run\code\run_p6.py" --jackknife --methods M1 --warm-hybrid --workers 4
& $py "$run\code\analyze.py"
& $py "$run\code\report.py"
& $py "$run\code\audit.py"
```

`campaign.py` 按上述拟合顺序运行，避免两个四进程拟合池重叠。传入一个当前 P6 调度进程的 PID 时，先等待该进程结束。每个阶段的退出码保存在 `analysis/stage_*.json`。报告生成是独立步骤，不以启动了任务等同于完成任务。

## 文件

- `code/profile_core.py`：原 ZOH/瞬态模型的坐标变换，四种留数结构，剖面参考 IV 求根。
- `code/verify_profile.py`：原完整矩方程与剖面矩方程的根一致性；秩 1 坐标导数核查。
- `code/run_p6.py`：六记录、三窗口、四模型，M3 和 M1 完整拟合；公共时延 41 点网格；删周期拟合。
- `code/run_fixed_low.py`：低频公共时延固定到同记录两高频窗口 M3 估计的平均值；±25 μs 敏感性；配对 M1/M3 jackknife。
- `code/metrics.py`：冻结公共噪声标尺的 Q/E[Q]、局部分频段残差、时延剩余信息，以及同记录 M⁻¹b 预测。
- `code/run_p7_sim.py`：实验 A 原数据与种子的 P=4/8/16、各 40 次 MC；M2、M3、删除周期修正。
- `code/run_p7_real.py`：原 158 参数低频模型上的实测 P7；嵌套删二周期估计修正量的 SE。
- `code/analyze.py`：CSV、P6 对照和低频联合 jackknife z、MC 聚合。不会删掉失败有限端点。
- `code/report.py`：中文报告、P7 图（PNG 300 dpi / PDF）。
- `code/audit.py`：原输入哈希、完成数、失败数和矩阵结构的核对。
- `raw/`：每个拟合独立 JSON，包括全部初始网格、完整参数、求解历史与异常。
- `simdata/`：仅保存原实验 A 未运行的 P=16 后 20 个种子的数据；原数据文件不改。
- `analysis/`：结果表、SE、判定、阶段状态和验证。

## 复现与解释注意

参考工具依赖模型留数，故内层不是“一次闭式解即完成”，而是冻结工具闭式解与工具更新直到线性矩满足。秩 1 振型方向属于外层非线性参数。最终检验的是与原实现相同的 IV 矩方程，未改成输出误差最小二乘。

M3 的实测拟合使用激励谱线，M1 使用全频点；仿真 P7 的 M3 沿用实验 A 的全频点设置。无论非线性坐标是否收敛，原失败端点都保留。P7 MC 保留全部有限值统计并单列失败；不完整或含必要失败的判定标记为不充分。

P6 的主要时延误差棒针对 M3；M1 完整端点和时延 jackknife 也单列，低频核心的 M1/M3 都做删周期。`--warm-hybrid` 从完整数据端点直接解剖面矩根；与原 90 次 GN 再 hybrid 的两个 H300 对照根差小于 1.3e-12，核查文件见 `analysis/profile_hybrid_benchmark_*.json`。停止门槛和模型没有改变。

P7 实测为隔离去偏算法，继续使用原固定时延自由留数模型；不与 P6 改结构混杂。以上范围解释在报告的协议说明中明确列出。

低频锚定公共时延的 jackknife 是条件 SE；高频估计锚定值的随机误差没有自动包含在内。±25 μs 结果仅为敏感性，不能当作已传播全部不确定度。

## 最终完成状态

本轮拟合已结束，所有要求运行的端点均已落盘；原输入清单 4215 个文件的哈希全部不变。完成不等于判据通过。

- P6 完整端点 144/144，7 个未通过收敛判据。
- P6 主分支删周期：M1 465/465（4 失败），M3 472/472（0 失败）。失败父端点对应 M1 的 27 折、M3 的 20 折未运行，不能赋予有效 SE。
- 固定低频公共时延分支：M1 70/70（1 失败），M3 69/69（2 失败）；全部 72 个锚点及 ±25 μs 敏感性端点保存。
- P7 仿真各 P 均 40 次，共 120 次；240 个完整 M2/M3 拟合、1120 个删除周期拟合齐全。M3 在 P=4 的第 31 号重复失败，有限端点仍在统计中。
- P7 实测 6 条记录、215 个拟合全部完成并通过数值门槛。

`analysis/audit.json` 给出逐项清单、失败索引、哈希和缺项核对。报告生成顺序推荐：`analyze.py` → `audit.py` → `report.py`，使审计状态进入报告。

新增交付索引：

- `analysis/p6_cross_window_components.csv`、`p6_M1_cross_window_components.csv`：九个规范化时延分量的跨窗口有符号差和同折配对 SE。
- `analysis/p6_fixed_low_sensitivity.csv`：固定公共时延在锚点及 ±25 μs 时，各模态的频率、阻尼与失败状态。
- `analysis/p6_M1_dispersion.csv`：M1 六记录 SD 及 A/C 配对复现性。
- `figures/p6_model_comparison.png/.pdf`、`p7_damping_bias.png/.pdf`：300 dpi PNG 与矢量 PDF。

结论：P6 不支持统一删除 C1 或统一秩 1 替代；低频配对差 SE 未改善。P7 仿真阻尼 9/9 项符合冻结门槛。实测按预定平方和 SE 门槛阻尼 16/18 项通过，计入同记录相关性的配对 SE 诊断仅 5/18；不能把它写成两方法整体等价。
