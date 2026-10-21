# 周期 SML 与间接 FRF 对照的复现说明

日期：2026-10-19。只增加本目录的新代码和输出；原始数据、原理论和原结果不改动。

## 文件与运行顺序

- `sml_core.py`：通用周期 SML，包括均值协方差、模型相关协方差的精确 Cholesky 导数。
- `run_sml.py`：与实验 A 完全相同的 U/Y 记录、初始扰动；读取 `../experiment_a/data/Pxx_rxxx.npz`。
- `run_indirect.py`：四次单轴实验、间接 FRF、全协方差传播，以及信赖域最小二乘优化器对照。
- `run_indirect_ivgn.py`：主要 M6u/M6_reg 的 IV–Gauss–Newton 拟合；保留原 `cl2.iv_solve` 的正规方程、对角缩放、0.05 对数极点步长上限、80 次上限、1e-5 极点步长停止准则。
- `run_sensitivity.py`：有限配对样本的协方差收缩系数敏感性。
- `summarize.py`：保存六个模态参数的均值、SE、SD、RMSE，以及失败、耗时、条件数。
- `validate_sml.py`、`validation_*.json`：协方差导数、无噪声均值、间接 FRF 和快速仿真校验。
- `summary.csv/json`、`results_tables.md`：所有有限端点的统计。包含未收敛端点，不做异常值删除。
- `availability.json`：未正则化 M5/M6 不可计算的条件与样本数，未将不可计算计成成功。
- `sml_results/`、`indirect_results/`、`indirect_ivgn_results/`：完整 beta、beta0、truth、百分比误差、状态、求解次数、耗时；主运行每完成一次落盘。
- `indirect_data/`：四次独立实验的逐周期 DFT、FRF、协方差和条件数。
- `sources/`：本地用户论文的文本提取，仅为核对参考，不修改论文。

Python 3.11、NumPy 2.4、SciPy 1.17；使用原 `../../code/cl2.py`。参数 `NOFLEX=1 FC=60 AL=30 AH=100 FHI=100 SW=.03` 在运行文件中设置。每周期 30000 点，5000 Hz。所有拟合频率均为 3–100 Hz；保留 ZOH 精确求值。

在 workspace 根目录的 PowerShell 中：

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONDONTWRITEBYTECODE='1'
python research/closed_loop_consistency_20261019/gpt_runs_20261019/comparators/validate_sml.py
python research/closed_loop_consistency_20261019/gpt_runs_20261019/comparators/run_sml.py 0 1
python research/closed_loop_consistency_20261019/gpt_runs_20261019/comparators/run_indirect.py 0 1
python research/closed_loop_consistency_20261019/gpt_runs_20261019/comparators/run_indirect_ivgn.py 0 1
python research/closed_loop_consistency_20261019/gpt_runs_20261019/comparators/run_sensitivity.py
python research/closed_loop_consistency_20261019/gpt_runs_20261019/comparators/summarize.py
```

两个互不重叠的工作进程分别用参数 `0 2`、`1 2`。已有原始结果按文件名跳过，不覆盖。`run_sml.py` 和 `run_indirect_ivgn.py` 可等待上游逐次数据生成。SML 和间接法的原始结果为不同目录，不混合。额外执行的 P16 工作进程只复用/补齐已存在检查点；完整记录数见报告。

## 与需求中伪代码的必要区别

1. 周期确定性未知输入剖面 SML 的目标确为 sum e^H Ce(theta)^(-1)e；不额外添加 logdet(Ce)。但每轮仅冻结 Ce、忽略其导数的 IRLS 一般不求这个目标的驻点。本实现对白化残差整体求导，包含 Ce 导数。SML 使用信赖域 Gauss–Newton 最小二乘。
2. 2007 年一般随机输入 EIV-ML 的似然包含输入/噪声模型及行列式项；不能直接将其 logdet 移植到周期剖面 SML。周期公式直接核对 Hong 等 IFAC 2008 第 4 节式 (11)。
3. M5 原始四输出 Ce 在 P=2/4 秩至多 1/3，不能求逆。对应结果只作为 `M5_reg`，以固定 0.1 向联合 [Y,U] 协方差对角收缩。P=8/16 的 `M5` 不加正则。P=1 没有样本协方差。
4. 四次实验、每次 p=P/4 周期时，完整 vec(G) 的 16×16 样本协方差秩至多 4(p−1)=P−4。P=4/8/16 全部不足 16，因此原始 M6 均不可求逆。P=8/16 的 `M6_reg` 采用固定 0.1 对角收缩；P=4 连方差估计都不存在，只做 M6u。
5. 原 `cl2.iv_solve(free=...)` 会关闭非线性极点导数且只执行一次，并且只接受对角输出权；不能直接用于无瞬态、全协方差 FRF 拟合。因此新文件保留原 IV–GN 正规方程并做这两个必要适配。`*_IVGN` 为主要结果，信赖域优化器用于检查差结果是否只由求解器导致。

## 预算与随机数

SML 完全复用实验 A 数据：seed = 202610190 + 1000P + r，连初始极点扰动也相同。

间接法实验 a=0,1,2,3 的 seed = 202610190 + 1000P + r + 700000 + 100000a；四个实验彼此独立。每次只有第 a 个执行器激励，其他外加激励为零，所有实际电机总输入仍被反馈耦合。每次先丢两周期，记录 P/4 周期。四次记录总样本数等于同时四路实验的 P 周期；但预热成本为 8 周期而单次为 2 周期。保持每个活动执行器的激励幅值不变，因此间接法总注入能量为同时四路的约四分之一。这是按用户预算实施的对比，不能称为同时匹配记录时长和总激励能量。

## 来源

- Pintelon & Schoukens, *Frequency domain maximum likelihood estimation of linear dynamic errors-in-variables models*, Automatica 43 (2007), 621–630：本地全文已读；它讨论的一般随机输入 ML 不等于这里的周期 SML。
- Pintelon & Schoukens, *Box–Jenkins identification revisited—Part I: Theory*, Automatica 42 (2006)：本地全文核对；指定频带闭环 ML 的均值与噪声模型依赖控制器。
- Hong, Söderström, Schoukens & Pintelon, *Accuracy Analysis of Time-domain Maximum Likelihood Method and Sample Maximum Likelihood Method for Errors-in-Variables Identification*, IFAC 2008，第 4 节式 (10)–(12)：https://skoge.folk.ntnu.no/prost/proceedings/ifac2008/data/papers/0643.pdf 。其 P≥4 的明确论述针对 SISO，不能忽略 MIMO 协方差秩条件直接沿用。
- 2007 EIV 官方条目：https://www.sciencedirect.com/science/article/abs/pii/S0005109806004389
