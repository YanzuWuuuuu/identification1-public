# 时延理论的仿真与实测收尾

本目录保存新旧记录对齐、全八模态闭环仿真的稳定性/可辨识性审计、P2–P4 Monte Carlo 与 P5既有实测复核。所有新增输出在本目录；父级原数据、理论和历史结果只读。初始 SHA256 清单为 `input_hashes.json`。

## 计算安排与不变条件

- P20：M1，L=3–100、H=100–550、W=3–550 Hz，每组40次。
- P10/P20/P40：M3，L/H，每组40次。
- 总计360个原始拟合端点；八模态、自由C0/C1/D0/D2、逐模态瞬态和7时延，共375参数，真背景为零。
- 5kHz，每周期30000点，先稳定两周期；种子1002+r，初值扰动100200+r。
- 初始公共时延200μs恰好是真值，但相对时延不使用真值；已知局部分支，不执行实测的41点公共时延网格。
- 失败全部保留。正式主统计统一只取 `mc_qr_final/raw`，不存在时就是未完成，不回退至其他求解路线。

## 脚本

在本目录中运行：

```powershell
$py = 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe'
& $py -B code\validate_sim.py
& $py -B audit\stability_identifiability_audit.py
& $py -B code\run_qr_mc.py --periods 10,20,40 --count 40 --workers 4
& $py -B code\predict_bias.py
& $py -B p5\analyze_p5.py
& $py -B code\plot_audit.py
& $py -B code\summarize.py
```

对齐核查入口是 `alignment/check_alignment.py`。按当前已完成检查与任务调度复用结果，不需要为生成报告重复耗时计算。后台并发应服从主进程安排，勿在活动MC池之外再启动同一批任务。

## 报告汇总

`code/summarize.py` 仅读取结果，不拟合，也不会把未完成项填成零。正式估计来源唯一：

1. `mc_qr_final/raw/` 中同名最终端点，不论成功或失败；
2. 没有这一正式端点时明确标记未完成；
3. 不回退到正规方程、额外完整求根或中止试版 QR，也不在路线间挑选成功结果。

`mc/raw` 的正规方程路线与 `mc/refinements` 的额外375维矩求根只作数值诊断，状态和误差单列在端点清单，不纳入主统计。`mc_qr/raw` 是计算求和顺序校验修改前启动即中止的四条试版结果，仅列状态。未保存复算的原始失败不能称作“完整求根也失败”。

正式版把 QR/SVD 的 `1e-13` 截断施加在原设计矩阵，内层检查原始未压缩线性矩的 score `<1e-9`，完整原375维矩的 score `<1e-6` 及贴边标准不变。旧内层 `2e-12` 对 H 频带过严，这是明确记录的数值修复，不应声称“没有改任何求解容差”；它不是预先冻结的科学判据。正式外层依然调用完整剖面矩的 hybrid 求根。

正式 MC 的均值、标准误、SD与RMSE只对正式路线中的成功子集计算，并明确失败选择可能影响统计结论；其余端点在清单中保留。

主报告为 `时延理论仿真验证报告.md`。生成内容包括：

- `analysis/report_summary.json`：完成度、成功失败与P4状态。
- `analysis/mc_endpoint_inventory.csv`：360个计划样本、正式端点状态与误差、各诊断路线的单列状态与误差。
- `analysis/mc_groups.csv`：每组正式完成度及正规方程/额外root/中止QR试版覆盖率。
- `analysis/mc_statistics.csv`：仅正式路线的全部时延、八模态频率/阻尼成功子集统计。
- `analysis/p4_component_checks.csv`：两频带的P10/P20/P40时延趋势和P40的2MCSE检查；公共与零和相对分量并非独立。
- `analysis/prediction_diagnostics.csv`：总体偏差各积分/伪逆截断值，不选取一个有利截断代替稳定性核查。

`code/plot_audit.py` 只读取已保存的回路网格和可辨识性审计，不重新拟合。它输出 `figures/audit_high_band_loop_gain` 与 `figures/audit_band_information` 两组 PNG（300 dpi）和 PDF，以及对应 CSV 和 `analysis/diagnostic_figure_checks.json`。第二张图呈现数值雅可比与投影诊断，不能单凭它宣称精确 ZOH 模型严格不可辨识。

对应正式 MC 尚未全部保存时写“未完成，暂不能判定”，不能因未运行样本而填写“未通过”。全部运行完成后仍缺乏可靠数值端点，才写“未通过（数值不可检验）”，不称理论被证伪。P2/P3带外小回路前提在当前100–550Hz不成立，这是可独立报告的条件核查；尚未完成的MC比较仍标待完成。P5原文没有规定覆盖因子，报告1u/2u描述性比较而不追认为预设阈值。完整说明分别在 `alignment/`、`audit/` 与 `p5/` 的报告中。

## 自动进度更新与最终哈希核查

当前有限批次已配有 `code/watch_run.py` 伴随进程。它每新增8个正式端点自动重新生成报告，批次进程退出时再生成一次，并对冻结清单中的原文件逐个核验 SHA256，写入 `analysis/source_hash_verification.json`。该进程是当前已启动批次的有限监看，不是另建定时任务；不应重复启动相同监看或重复MC。

`RUN_STATUS.json` 保存实时的 `saved/successful/failed/remaining/running/complete`。`complete=false` 表示该批次尚未完成，报告也是相应时刻的进度稿。只有360个正式端点都存在、无无效JSON且批次正常完成标记存在时，监看才将 complete 置为 true；异常退出也会保留状态并执行最终报告和哈希核查，不会伪称完成。
