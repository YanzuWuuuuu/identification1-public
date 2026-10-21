# 新采实测 A–E：复现说明

本目录为2026-10-20的新分析，原始MAT、旧结果和理论文档只读。入口报告是《新实测数据分析报告.md》。结论是部分支持，完整预设不变性/不确定度覆盖没有通过，不得只引用第二阶的有利结果。

## 数据与算法

- 原始数据：`C:/Users/Administrator/Desktop/my_article/identifiation1/sweep/sweep_data_new`。
- 27个MAT中D的1.5档与A nominal完全相同：26个独立记录，23个有激励、3个B噪声记录。
- 5kHz，30000点/周期；按公共时间轴、实际注入门控和预先固定的周期差规则截取。实际稳态P为5–8，不补造计划中的18周期。
- 模型保持既有实测结构：3支承模态、D0/D2、固定时延减100µs、精确ZOH和24个逐模态瞬态，共158参数；频带3–25Hz。
- M1全频点直接IV；M2谱线直接IV；M3参考IV；M3w一次16邻点噪声权更新；M5n邻点联合协方差SML。
- `baseline/`保存实际使用的旧代码副本及SHA256来源清单。旧文件不被导入后修改。
- 主计算94次，4次失败端点保留，26次未过严格score门槛。失败与严格检查是两套标志，不互相替代。

## 环境与命令

本机Python 3.11.9，numpy、scipy、matplotlib、Pillow；PDF检查另外使用PyMuPDF。准确版本见`analysis/delivery_checks.json`。无需MATLAB或新增硬件采集。所有脚本入口均相对自身目录定位，不依赖当前目录。

PowerShell运行顺序：

```powershell
Set-Location -LiteralPath 'C:/Users/Administrator/Desktop/my_article/identifiation1/research/closed_loop_consistency_20261019/gpt_runs_20261020/newdata'
$pythonExe = 'C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe'
& $pythonExe prepare_data.py
& $pythonExe identify.py --all --workers 3
& $pythonExe nonparametric.py
& $pythonExe uncertainty_bias.py
& $pythonExe sml_uncertainty.py
& $pythonExe evaluate_results.py
& $pythonExe verify_analysis.py
& $pythonExe plot_results.py
& $pythonExe build_report.py
& $pythonExe delivery_checks.py
```

单条调试：`& $pythonExe identify.py --id A_nominal --methods M1,M2,M3,M3w,M5n`。

现有拟合结果和谱缓存会被复用，失败记录也不会自动删除或重跑。仅在完全相同的预处理、设置与源数据下可直接续跑；改变任何设置时应复制到新的输出目录并使用新的prepared/spectra/fits/analysis目录，不能混用旧缓存。未实现自动依赖失效，不要把改设置后跳过缓存误当重新计算。各拟合完成后立即保存JSON/NPZ，重抽样也逐次保存。

## 脚本与产物

| 文件 | 作用与输出 |
| --- | --- |
| prepare_data.py | raw-deflate读取、时间对齐、门控、周期质检、去重；B前半ARX(8,8)、独立后半与各运行验证；prepared/、controllers/、data_inventory.json |
| identify.py | 原模型与求解器复用，94次主拟合；spectra/、fits/、fit_requests.json |
| nonparametric.py | C四次正交FRF、1000次成对周期bootstrap、跨控制器相容性、1–500Hz控制器/回路诊断 |
| uncertainty_bias.py | B谱、完整瞬态协方差、直接法M与M−C偏差、反馈噪声占比 |
| sml_uncertainty.py | 6组M5n各99次成对周期重估，保存594次状态及模态参数 |
| evaluate_results.py | 固定判据、独立A/C验证、条件Gaussian检验、幅值/位置/重复表和旧数据对照 |
| verify_analysis.py | 独立控制器恢复、正交符号/条件数、二次矩公式数值检查、原文件哈希 |
| plot_results.py | 20幅图，每幅300dpi PNG与文字可编辑PDF，来源manifest |
| build_report.py | 汇总报告与run_manifest.json，直接从落盘结果生成 |
| delivery_checks.py | 交付完整性、可解析性、参数行数、图片分辨率、PDF文字与源完整性检查 |

完整参数表：`analysis/all_parameter_estimates.csv`，564行；每个估计均有u、迭代数、失败及严格score状态。关键表还包括`controller_invariance_pairs.csv`、`bias_controller_drift.csv`、`A_C_independent_repeats.csv`、`frf_normalized_validation.csv`、`scans.csv`和`old_new_position_comparison.csv`。

## 统计口径与复用限制

1. 所有参数相对差用100·expm1(对数差)，频率Hz，阻尼表为百分数，阻尼u为百分点。图中误差棒是±2u，非不含条件的95%覆盖保证。
2. M3用定理2的完整158维sandwich代入；M3w冻结估计权；M1/M2用Gaussian线性+二次矩局部近似，协方差尚未含曲率C；M5n用99次条件周期bootstrap，不称其为定理2。
3. B噪声谱在不同控制器、A/B不同运行之间不完全一致；E噪声借用z=0，误差棒只能条件解释。
4. C参照采用未经正则化的YU⁻¹，nominal/highgain存在明显比值不稳定。灰带是点态周期bootstrap，不含运行间变化、系统偏差，也不是同时带。
5. C模型主验证采用独立A记录的参数；χ²仅给固定模型/已知PSD理想参照，另给共享参数扰动的条件模拟校准。M5n未提供全158参数u校正。
6. D是幅值扫描；同幅重复是A nominal加C nominal四次，共5次。不能把不同幅值当同分布重复样本。
7. 新E的−1mm记录虽通过门控检查，但3个方法未收敛；没有沿用旧“坏记录”标签删掉它。
8. 失败用红叉、严格数值检查未过用空心点。相位独立展开后只作每条曲线一个整数倍360°平移，所有复FRF和误差不变；低信噪比引起的局部相位解缠歧义没有人为消除。

## 日志与随机数

`logs/`保存各阶段运行记录。FRF bootstrap seed为20261020加控制器编号；B噪声区间seed=621020；条件Gaussian模型验证seed=421020；其他重抽样/矩检查种子写在对应脚本与JSON中。SML每次子估计的异常也单独落盘。

本目录是一份完整研究分析包；图形尺寸用于审阅，投稿时仍需按目标期刊栏宽缩排。原始文件哈希检查通过不等于物理模型/统计判据通过，二者已分别报告。
