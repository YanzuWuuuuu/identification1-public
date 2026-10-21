# 高频偏差与复现性

所有结果只写入当前 `bias_check` 目录。原模型、数据、拟合结果及理论文档只读。
Python 3.11，numpy/scipy/matplotlib；无需MATLAB、LRM或新数据。

从本目录在PowerShell运行：

```powershell
$py = 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe'
& $py -u code\validate.py *> logs\validation.log
& $py -u code\bias.py --all --workers 2 *> logs\full.log
& $py -u code\bias.py --all --folds --workers 2 *> logs\jackknife.log
& $py -u code\analyze.py *> logs\analysis.log
& $py -u code\plot.py *> logs\plot.log
& $py -u code\report.py *> logs\report.log
& $py -u code\audit.py *> logs\audit.log
```

两个工作进程、各1个BLAS线程；按记录/折保存，可中断后继续。已存在的拟合/预测文件直接读取，
不会静默替换失败样本。代码改变若需重算，应另建命名结果目录，不删除本次审计记录。

- `baseline/modal_model.py`、`tau_seconds.npy`：原实现的逐字节副本。
- `code/common.py`：只读数据入口、原初始化/求解/精化流程和状态。
- `code/bias.py`：0.5Hz非激励点谱、两种偏差预测、条件sandwich。
- `code/validate.py`：原谱重现、分块/直接矩阵等价、谱估计与差分步长核查。
- `code/analyze.py`：全部jackknife误差棒、预设检验、复现性、带宽参照和D状态。
- `code/plot.py`、`report.py`：仅根据已保存结果生成PNG/PDF及中文报告。
- `code/audit.py`：原文件哈希、完整折数、联合方差恒等式、求解状态与图文件核查。

H300所有200参数、H493所有74参数均保留；仅不把H300第一宽背景分量作为机械模态报告。
主分析使用原diagnostics中已通过严格score的M1/M3端点；周期重估也执行同一精化程序。
低频图点与误差棒直接引用上一轮closeout结果，没有修改或精化低频端点。

主表：`analysis/high_frequency_tests.csv`（96行），
`bandwidth_inside_outside.csv`（低频第二阶+高频四模态，两种公式），
`repeat_pairs.csv`、`reproducibility.csv`、`D_amplitude.csv`。
完整每折模型、噪声谱和偏差矩阵在 `raw`；每个原始/精化状态都保留。

误差单位：频率/阻尼相对差均为100×expm1(log差)；阻尼参数本身为ζ×100；
复现性阻尼差的物理单位是“阻尼百分点”，不是相对百分比。
正式效应检测使用实测配对差SE，吻合z使用同折实测−预测的联合SE。

本次任务是高频旁证和不确定度审查，不以小绝对数值替代统计吻合，
不把三对运行的池化SD称为完整实验总不确定度。
