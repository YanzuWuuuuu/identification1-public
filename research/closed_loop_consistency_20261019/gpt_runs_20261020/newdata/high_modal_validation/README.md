# 高频实测辨识扩展

先读[中文报告](高频实测模态验证报告.md)，固定协议见[PLAN](PLAN.md)，未通过端点见[完整清单](失败与未通过端点清单.md)。本目录全部是新增实验；低频估计数组、原始MAT、理论文档没有改写。收尾时一次误调用重生成父报告及运行清单，详见[执行说明](implementation_note.md)。

## 内容

- `modal_model.py`：输出分块的加性ZOH模型、解析Jacobian及IV矩方程；结构与旧实现逐项数值核对。
- `run.py`：4固定窗口×6记录×5方法，120主估计。采用局部D0/D2背景，不是全3–800Hz模型。
- `analyze.py`：B噪声、M3/M3w条件sandwich、留数匹配、独立FRF/输出验证。
- `polish_diagnostic.py`：仅精化原已收敛且非边界IV解，61组；不覆盖原协议。
- `polish_uncertainty.py`：精化后条件u检查，避免停止阈值造成错误排除。
- `amplitude_check.py`：4条D记录×2窗口，仅M3，8组次要幅值诊断。1.5与A nominal相同，不另算。
- `plot_results.py`：41张图，各300dpi PNG与文字可编辑PDF，另存10张小预览（4参数、6FRF）及1张幅值预览。
- `build_report.py`：从保存端点自动生成报告和未通过清单。
- `qa.py`：统计完整性、PNG/图表/公式控制字符、源码语法、数值核查及输入/代码哈希。
- `baseline/`：两份旧依赖原样复制；运行还只读父目录`prepare_data.py`与其依赖。

## 运行（Windows PowerShell）

Python为既有Python3.11，依赖numpy/scipy/matplotlib/Pillow。先进入本目录，再依次执行。机器上可将`python`替换为`C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe`。

```powershell
python validate.py
python run.py --all --workers 3
python analyze.py
python polish_diagnostic.py
python polish_uncertainty.py
python amplitude_check.py
python plot_results.py
python build_report.py
python qa.py
```

单个主拟合例子：

```powershell
python run.py --id A_nominal --window H493 --methods M1,M2,M3,M3w,M5n
```

脚本按方法落盘。已有主端点会跳过，并核验窗口配置哈希；这不是对任意源码变更的缓存有效性保证。若改变模型/算法，应复制到另一个新实验目录，不能把旧缓存混作新运行。建议BLAS单线程，主任务最多3进程；配置已在入口设置。

## 输入与输出边界

输入来自父目录的`prepared/*.npz`、`controllers/*.npz`、`baseline/tau_seconds.npy`和`wideband_scope_audit/*_broadband.npz`，均对应新A–E记录；具体文件及SHA256见manifest。所有估计输出只在本目录。原MAT采集门控与原噪声模型的限制见父报告，本轮没有绕开原质量问题。

`raw`含120个端点JSON、120个模型/残差NPZ、24个局部协方差NPZ。JSON保留beta、初值、迭代/评估数、状态、score、边界、失败。`spectra`为24个独立窗口DFT缓存，模型与数据单位保持原记录。

`analysis/modal_parameters.csv`含360行；`controller_pairs.csv`360行（包括背景，主判据排除）；`independent_repeats.csv`180行；`FRF_validation.csv`60行。`summary.csv`20行。u仅按当前正确局部模型与B噪声条件计算，不是总实验不确定度。

`diagnostics`是61组精度续算及符合范围的协方差；`amplitude_diagnostic`是8组额外D结果。它们不写回原120组。

## 图索引

- `H100/H300/H493/H700_parameters`：五方法、三控制器、A/C两次记录。圆点A、三角C；叉号失败/贴边；淡色为原严格score未过；方框为模态对应未通过。
- `{lowgain,nominal,highgain}_{H100,H300,H493,H700}_{magnitude,phase,residual}`：独立C四实验FRF与A参考IV模型，4×4物理通道。残差为复数差的模。
- `amplitude_frequency_check`：同控制器D幅值诊断；不是误差棒，不声称唯一物理非线性机制。

相位分别unwrap，仅整条平移一个360°整数；未按模型逐点消除数据绕圈。灰色FRF是估计而非真值，绿色前后段差不是纯噪声。未收敛模型也保留图并标注。图中没有橙线。

本轮未做高频E位置扫描、完整宽频联合拟合、重新估计通道时延或新非线性模型；没有把这些内容列为已完成。
