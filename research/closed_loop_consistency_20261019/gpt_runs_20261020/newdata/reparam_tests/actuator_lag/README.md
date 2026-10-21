# 输入端一阶滞后检验

此目录是 P6 留数模型比较之后的独立补充实验。判据在 `PROTOCOL.md` 冻结。已有原始数据、原模型和结果保持只读，新增拟合、核查和汇总均保存在本目录。

## 模型与任务

六条 A/C++++ 记录分别进行 H300（220–420 Hz）和 H493（450–550 Hz）的联合 M3 估计。两窗口共享七个结构化时延、四个输入端一阶滞后截止频率；保留各窗口原有极点、背景与输出瞬态。

- `a` / 模型 i：滞后＋自由 C0、C1。
- `c` / 模型 ii：滞后＋自由 C0、C1=0。
- 比较基线是上一级 `reparam_tests` 中原自由模型 a 的独立窗口拟合。
- 模型 ii 需要 12 个窗口残差比全部 ≤1.2、四路截止频率六记录 CV 全部 <20%，且六个端点全部有效，才建议采用。失败样本不删除。

截止频率从 2 kHz 初始化；技术搜索范围 20 Hz–200 kHz。公共时延沿用 ±200 μs、10 μs 网格。输出瞬态不乘输入端一阶滞后；不重新整体 ZOH 离散化，亦不再额外扣除半个采样周期。本次不启动新的 jackknife。

## 运行

以下 PowerShell 命令在本目录运行，使用本机已有 Python/numpy/scipy/matplotlib：

```powershell
$py = 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe'
& $py code\verify_lag.py
& $py code\run_lag.py --workers 4
& $py code\audit_lag.py
& $py code\summarize_lag.py
```

单个记录或模型可用已有 runner 的选项，例如：

```powershell
& $py code\run_lag.py --ids A_nominal --models a,c --workers 2
```

`run_lag.py` 遇到已保存的同名端点会读取并返回，因而可在中断后继续，不覆盖已完成端点。每个端点对应 `raw/{record}_{kind}_joint.json`，共同延迟初始化另有 `_initialization.json`。失败与异常同样落盘，不能手工删除后仅保留成功结果。

## 汇总不触发拟合

`code/summarize_lag.py` 仅读取保存的结果和原基线，不导入拟合模块，不运行优化。重复执行只更新本目录的报告、汇总表和图。端点未齐全时明确标注“尚未完成”，不会提前作模型选择结论。

```powershell
& $py code\summarize_lag.py
# 仅重建数字表和 Markdown 报告：
& $py code\summarize_lag.py --no-figures
```

最终交付前先完成数值/哈希审计，再重跑一次汇总；汇总器会读取 `analysis/*audit*.json` 与 `analysis/*verification*.json` 并写入报告。它不自行修改输入，也不声称尚不存在的审计已通过。

## 输出

- `执行器一阶滞后检验报告.md`：预定判据、结论、残差、截止频率、时延、模态差异及失败状态。
- `analysis/lag_summary.json`：完整性与固定判据的机器可读结论。
- `analysis/lag_endpoints.csv`：所有必需端点的有效性、score、边界、耗时、参数数。
- `analysis/lag_qe.csv`：每记录×模型×窗口的 Q、E[Q]、n、E[Q]/n、原 a 基线与相对比值。
- `analysis/lag_cutoffs.csv`、`lag_cutoff_stats.csv`：四个截止频率、六记录样本 SD/均值 CV、极差。
- `analysis/lag_delays.csv`：公共与输入/输出相对时延，单位 μs。
- `analysis/lag_modal_differences.csv`：282/333/359/493 Hz 的 f、ζ，及相对原 a 的绝对/相对变化。
- `analysis/lag_structural_ranges.csv`：三个结构选择间的模态范围，逐行标注是否所有端点有效；不把范围自动解释成概率不确定度。
- `figures/lag_residual_criterion.{png,pdf}`、`figures/lag_cutoff_reproducibility.{png,pdf}`：两张主判据图，PNG 为 300 dpi。

CSV 使用 UTF-8 BOM，便于在中文 Excel 中打开。阻尼绝对差采用“百分点”；频率和阻尼相对差采用 `100×expm1(log(新值)−log(基线值))`。

## 解释边界

本次新增滞后与共享跨窗口动态同时出现，而原基线各窗口时延独立；不能把所有改善单独归因于滞后。通过相对 Q/E[Q] 判据不等于绝对残差已经达到噪声水平。CV 只衡量六条记录间离散，不是单条记录内置信区间。模型 ii 只要任一预设条件不通过，主模型仍为原自由 C0、C1。

## 本轮已完成结果

十二个联合端点均已保存，九个通过原数值门槛。模型 i 的 A_lowgain、A_nominal、C_nominal++++ 三个端点未通过；A_highgain 和 C_highgain++++ 的第四路截止频率约 199.6 kHz，距技术上限不足 1%，并存在很大的缩放条件数。这项上限接近诊断不改变原成功状态，也不能作为硬件截止频率解释。

模型 ii 六个端点全部通过求解，四路截止频率 CV 依次约 2.12%、3.31%、19.70%、6.07%，数值均小于 20%。但是十二个窗口的残差比均超过 1.2：H300 为原自由模型 a 的约 2.31–3.40 倍，H493 为约 5.18–7.52 倍。因此按预先固定判据不采用模型 ii，主模型仍为自由 C0、C1。正式结论与全部精确数字由汇总脚本从保存结果生成。

原文件审计核对了 6221 个文件，均未改变。解析导数和矩方程核查通过；初始较大有限差分步长的失败核查仍保留。`analysis/residual_expectation_bootstrap.json` 对 A_nominal、模型 ii 使用 200 次抽样独立核对联合残差期望：两个窗口抽样均值与解析值之差分别约 −1.24 和 +0.87 个 Monte Carlo 标准误，独立迹公式相对差约 1e−14。
