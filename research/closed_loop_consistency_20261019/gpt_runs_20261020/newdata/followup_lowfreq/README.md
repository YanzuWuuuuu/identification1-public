# 3–25 Hz 低频补充分析复现

本目录为独立输出目录；上一级 `prepared`、`fits`、`spectra`、`controllers`、`analysis` 及原始采集文件只读。运行时不导入上一级带写入副作用的报告模块，所有模型导入都来自本目录的 `baseline` 副本。

环境：Python 3.11，NumPy、SciPy、Matplotlib、Pillow。BLAS线程固定1。主报告为《低频偏差预测与残差检验报告.md》。

## 执行

在本目录下，依次运行：

```powershell
$py = 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe'
& $py -u code\jackknife.py --workers 3
& $py -u code\noise_bias.py
& $py -u code\residual_checks.py
& $py -u code\summarize.py
& $py -u code\figures.py
& $py -u code\verify.py
& $py -u code\report.py
```

需要日志时，各命令后加 `*> logs\对应名称.log`。本次运行的日志已保存。不要改变工作目录为父目录后运行旧脚本。

- `common.py`：只读数据入口、模型副本加载、DFT、谱尺度和矩阵工具。
- `jackknife.py`：41个删除周期折×M1/M3=82次重估；每次立即存JSON。已有单折结果直接读取以支持中断恢复；不会悄悄筛选失败折。
- `noise_bias.py`：原B预测复现；0.5/1.0Hz同记录谱；两种偏差线性化；联合sandwich；裁切/瞬态诊断；条件激励缩放预测。可加 `--id A_nominal` 只重算一条。
- `residual_checks.py`：按原权构造158参数IV斜投影，实协方差尺度核对，200次proper Gaussian线性bootstrap，每2Hz分解，lowgain跨记录定位。
- `summarize.py`：固定判定、各CSV、C_highgain第一阶审计、D组响应线性核对。
- `figures.py`：11张图，每张300dpi PNG和矢量PDF；主图误差棒为周期jackknife 1SE。
- `verify.py`：检查受保护文件全部SHA-256、基线副本、原结果复现、投影代数、PNG/PDF导出。
- `report.py`：从JSON/CSV生成中文报告，不修改已有报告。

## 固定约定

- 频带3–25Hz；1–27Hz仅为主噪声谱窗口边界提供数据；同记录谱不额外收缩、不加正则。
- 模型为原158参数精确ZOH；时延已减100µs，不能再次减。
- DFT未归一化；`Phi=mean(v vᴴ)/N`，全长复协方差`N*Phi`。实协方差额外乘1/2。
- 差异用`100*expm1(delta_log)`；阻尼本身表为`100*zeta`。误差棒注明属于绝对阻尼百分数还是相对配对差百分数。
- 删除周期重估固定使用各方法在该记录保存的原始initial_beta，80步、tol=1e-5、cap=0.05；严格score<1e-6单列。
- 任务3使用M3/M3w原端点和原权。新的噪声谱不用于把M3w重新拟合成另一个方法。
- 图中的预测没有误差棒：尚未传播谱/模型/K估计不确定度；不把条件bootstrap解释成覆盖所有误差的置信区间。
- 原保护快照在本目录创建之前计算；`verify.py`只验证，不刷新快照。

全部失败、迭代上限、严格score异常和D组原有失败端点均保留。结果不满足预设判据，不调整模型/频带去凑成功。
