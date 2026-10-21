# 实测 D 复现

本目录只写入新文件。`realfit_original_copy.py`、`realfit_extended.py` 来自上一轮 C 的逐字复制，`loader_copy.py` 仅更正目录层级。数据从原实验目录只读加载。

环境：Python 3.11，NumPy、SciPy。每个进程只使用一个 BLAS 线程。共同权矩阵/SML 实现从父目录 `weighted_core.py` 和 `sml_core_copy.py` 导入。

```powershell
$py = 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe'
$dir = 'C:\Users\Administrator\Desktop\my_article\identifiation1\research\closed_loop_consistency_20261019\gpt_runs_20261019\followup_DE\real_d'
& $py (Join-Path $dir 'run_real_d.py') --z 0 --p 2 --start 1
& $py (Join-Path $dir 'run_real_d.py') --all --workers 2
& $py (Join-Path $dir 'summarize_real_d.py')
```

已存在的单方法 JSON 自动跳过，便于断点续跑。方法顺序为 M2、M3、M3w、M5、M5n，M3w/M5n 用 M3 初步估计构建协方差。`--methods M2,M3` 可只跑未加权对照，但任何需要权重的方法均需要 M3。

M3w 采用 16 个总非激励频点、残差协方差对角收缩 0.1、一次权重更新。M5n 使用同一邻域的联合协方差，按联合矩阵收缩 0.1。二者均不使用激励谱线估计协方差；初估仍来自激励谱线，不能宣称权重在有限样本下严格独立。

所有方法只拟合谱线。前后半段独立初始化与求解；其中一周期半段的 M3w、M5、M5n 明确标为不适用。重复上下文（例如 P4 前半就是 P2 完整窗口）共享同一计算结果，并在报告中披露。

输出：`raw/` 逐估计 JSON/NPZ 和协方差；`data/` FFT 缓存；`summary.json` 180 行上下文、半段差异、失败计数、来源校验；`实测D_结果与半段差异.md` 可阅读表格。
