# 无噪声控制器预检查

这是正式 MC 之前的确定性频响和稳定性检查，不读取、生成或拟合随机实验数据。

PowerShell 运行方式（项目根目录）：

```powershell
& 'C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe' -B 'research/closed_loop_consistency_20261019/gpt_runs_20261020/newdata/delay_sim_v2/controller_check/check_controller.py'
& 'C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe' -B 'research/closed_loop_consistency_20261019/gpt_runs_20261020/newdata/delay_sim_v2/controller_check/write_report.py'
```

- `controller_check.json`：所有扫描、稳定性、幅相和实现核对原始数值。
- `控制器无噪声预检查.md`：可读结论。
- `six_mode_definition.npz`：本检查六模态的连续/离散状态矩阵、留数、极点、延迟。
- `controller_base_grid.npz`：原六模态对象的 0.02 Hz 频率网格与回路奇异值、谱半径。
- `dependency_hashes_before.json`：运行前只读依赖哈希；JSON 内另含运行后哈希。
- `artifact_hashes.json`：本目录交付文件哈希。
- `run.log`：运行日志。

不包含正式仿真、估计器或 MC 驱动器；该目录不能被误认为用户已批准新的实验方案。
