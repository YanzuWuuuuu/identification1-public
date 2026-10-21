# 实验 C 复现

全部输入只读；根路径从本目录自动推导。`*_original_copy.py` 是原文件证据副本，不执行其云端路径；本机加载入口在 `run_c.py` 的 DATA、DELAY。

PowerShell（在本目录）：

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONDONTWRITEBYTECODE='1'
python verify_c.py
python run_c.py 0
python run_c.py -0.5
python run_c.py 0.5
python run_c.py 1
python run_c_sml.py 0
python run_c_sml.py -0.5
python run_c_sml.py 0.5
python run_c_sml.py 1
python summarize_c.py
```

结果存在时跳过，防止覆盖。单个结果有 JSON（完整参数和诊断）及 NPZ（参数数组）。如需完全重跑，请复制整个实现到另一个新目录或由用户明确指定新的输出位置，不能删除已有结果。

可选环境变量：C_P=8；C_MODE=ct,zoh；C_METHOD=M1,M2,M3；C_HALVES=0（省略两半）；C_SML_MAX=160。默认生成所有请求结果，IV最多80次，SML最多160次评价。SML依赖同级 `comparators/sml_core.py`。

注意：新旧模型比较采用相同初始模态参数和相同 M3 全频点规则；新 ZOH 版仅模态块精确采样，D0/D2 不变，时延减100µs。结果和失败定义详见 `experiment_c_report.md`。
