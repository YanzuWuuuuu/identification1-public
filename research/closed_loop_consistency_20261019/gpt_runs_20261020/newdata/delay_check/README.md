# 重现方式

所有命令在本目录执行，Python 3.11，依赖 numpy/scipy/matplotlib。
无需MATLAB；准备数据、控制器和旧拟合只读来自父目录。

```powershell
python -B code/validate.py
python -B code/run_real.py fixed --workers 4
python -B code/run_real.py joint --workers 4
python -B code/run_real.py folds --workers 4
python -B code/refine.py
python -B code/run_sim.py --validate
python -B code/run_sim.py --count 40 --workers 4
python -B code/refine_sim.py
python -B code/low_root_diagnostic.py
python -B code/analyze.py
python -B code/report.py
python -B code/audit.py
```

可把python换成本机
`C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe`。
单条试算如 `python -B code/run_real.py joint --id A_nominal --window H493 --method M3`。
脚本检测已有端点文件并跳过，精化结果另存refinements，不覆盖原失败。

- baseline：开始前复制的模型与cl2源文件。
- code/core.py：原模型的参数重排、7维零和时延、解析导数、IV方程及诊断。
- raw：252个固定偏移端点、36个全记录联合端点、246个删除周期端点、40个MC。
- refinements：原未通过端点的后续求根结果，含原失败标志；主统计优先读此目录，但同时报告原失败数。
- code/root_solver.py：保持相同IV矩条件、求完整矩函数的根；专门检查IV–GN失败。
- analysis：CSV、jackknife折数据、判据、数值验证和完整性检查。
- figures：3张PNG（300dpi）及PDF。
- logs：全过程与最初失稳仿真配置的预检失败。

每个文件中的微秒、阻尼百分比和百分点已有明确字段名。
差值使用`exp(log_new-log_old)-1`时为相对误差；`damping_percent`之差为百分点，
二者不得混为一谈。

M3的Q/E[Q]只在通过严格收敛的高频端点用于判定。低频失败时残差算子病态，
比值降低不是模型正确的证据。无噪声数值验证中归一化score的分母为舍入级，
应检查绝对残差及参数恢复误差，不能套含噪拟合的score停止阈值。
