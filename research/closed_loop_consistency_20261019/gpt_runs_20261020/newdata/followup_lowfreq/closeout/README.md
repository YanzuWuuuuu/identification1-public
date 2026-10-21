# 低频偏差预测收尾复现

原 `followup_lowfreq` 与其上级数据/结果只读；全部输出在本目录。先看 `收尾报告.md`；不能把新增背景的未收敛端点当作可靠效果。

环境：Python 3.11，NumPy/SciPy/Matplotlib；BLAS线程1，并发进程2。模型副本与原文件逐字节一致。`PLAN.md`记录本轮固定范围、诊断规则和背景对照约束。

在本目录运行：

```powershell
$py = 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe'
& $py -u code\curvature.py
& $py -u code\prediction_jackknife.py --workers 2
& $py -u code\sensitivity.py --workers 2
& $py -u code\present.py
& $py -u code\verify_closeout.py
& $py -u code\write_report.py
```

各命令可加 `*> logs\名称.log` 保存日志；本次日志已提供。前两类逐折/配置任务在已有结果时读取它们，支持中断恢复，不自动排除失败。

- `common.py`：隔离后的只读数据路径、谱与矩阵工具。
- `noise_bias.py`：从上一轮复制的预测/噪声函数，作为依赖使用，不需要单独运行其main。
- `residual_checks.py`：上一轮实现的保留副本；本轮敏感性残差由`sensitivity.py`计算。
- `curvature.py`：9种情形（6原始+3缩放），原/标准化范数、广义特征值、条件数、组贡献及6维Schur核对。
- `prediction_jackknife.py`：已有41个删除周期折，逐折重建DFT、0.5Hz噪声谱及两种预测；不重新筛选原82个端点。
- `sensitivity.py`：6记录×3频带×2背景=36配置，每配置M1/M3，共72端点。原背景158参数；额外二阶背景184参数。`--smoke`仅运行A_nominal两种3–25Hz模型。
- `present.py`：导出CSV和4套300dpi PNG/PDF图。
- `verify_closeout.py`：已有文件SHA-256、基线副本、新背景雅可比、退化模型、原结果及联合方差等式。
- `write_report.py`：从保存结果生成报告。

新增背景为全4×4实分子二阶项，初值40Hz/阻尼0.2；频率范围35–200Hz、阻尼范围0.005–2；增加8个对应瞬态系数。主模型和旧判据不改变。背景活动集只用于该敏感性模型，失败和边界均保存，不能以降低Q/E为由挑选失败端点。

误差棒为period jackknife 1SE，非模型无关置信区间。预测与实测同折相关，差值SE使用差值本身的jackknife，不作独立误差合成。新联合z是补充诊断，原点估计与判据仍有效地保持原样。
