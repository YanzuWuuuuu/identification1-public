# MATLAB 辨识流程修复说明

本项目保留 **3 个刚体方向、4 阶柔性辨识模态**。仿真对象的第 5 阶柔性模态继续作为未建模动态，不加入辨识模型。

## 运行

将 MATLAB 当前文件夹切换到本目录，执行：

```matlab
result = run_identification(true, true);
```

第一个参数控制绘图；第二个参数控制是否从 `0/u1.mat`～`u4.mat`、`y1.mat`～`y4.mat` 重算 LRM。上述命令会重算频响，再顺序运行 1～4。它不会重新启动 Simulink 仿真，也不会覆盖原时域记录。

仅使用当前频响重跑：

```matlab
result = run_identification(false);
```

若需要重新仿真，运行 `0/Multiple_main_Control.m`；它仍包含五阶柔性模态，固定随机种子为 42，并在完成频响计算后导出 `0/LRM_data.mat`。本次验证使用已保存的时域记录，未重新运行 Simulink。

依赖：MATLAB；第 3 步的有界加权最小二乘细化需要 Optimization Toolbox。重新仿真还需要原项目使用的 Simulink、Control System Toolbox。

## 目录和数据流

| 编号 | 入口 | 作用 | 新输出 |
|---|---|---|---|
| 0 | `0/estimate_frf_from_saved_data.m` | 从已有四输入四输出记录重算 LRM | `0/LRM_data.mat` |
| 1 | `1/mian1.m` | 固定极点的加权线性初值 | `results/stage1.mat` |
| 2 | `2/main2.m` | 联合估计刚体与第一阶柔性，提取刚体矩阵 | `results/stage2.mat` |
| 3 | `3/main3.m` | 四阶加性模型 RIV，再做有界加权输出误差细化 | `results/stage3.mat` |
| 4 | `4/para_convert.m` | SVD 初始化和模态结构投影 | `results/stage4.mat` |

无需手工复制参数文件。原目录中的旧 `.mat` 和调试脚本保留供追溯，但新版主流程不再读取旧 `belta_init.mat`、`para_init.mat`、`vecM_est.mat` 或 `flexiblemodal_LRMfrf.mat`。优先使用 `0/LRM_data.mat`；它不存在时，使用原 `1/LRM_data.mat`。

统一设置在 `common/identification_config.m`。改变配置或源数据后，从第 1 步重跑；后续步骤会拒绝配置不一致的旧结果。

## 修复内容

- 去掉辨识路径中的真值调用；第 2 步残余频响只减去估计刚体，第 3 步图线和复数残差直接来自保存的估计参数。
- 第 1 步保存完整分母初值及 80 个分子参数；第 3 步正确读取前四阶分子，不再错位跳过第一阶。
- 第 4 步保留三个刚体奇异方向，输出 64 个模态参数；柔性残差矩阵保持秩一。
- 统一归一化参数 `[a1; a2; vec(B)]`，分母为 `1+a1*s+a2*s^2`；模态映射雅可比包括分母归一化对分子的导数。
- 加权线性求解采用列缩放和 QR/最小范数求解，去掉直接求逆和未缩放的固定 `1e5` 正则项。
- 投影使用有阻尼、只接受下降步的 Gauss–Newton，按梯度、步长或目标变化停止，单独保存收敛状态。
- LRM 只拟合一个共同输出瞬态，去掉重复列，修正瞬态索引。兼容输出 `Tr` 每行重复同一瞬态，不能解释为各输入独有的瞬态。
- 有界加权输出误差细化是对原 RIV 的额外数值优化，结果中单独记录 `refinement`；不将 RIV 停滞标记为收敛。

## 频段和权重

初始化为 40～720 Hz，刚体估计为 5～310 Hz，柔性辨识为 200～720 Hz。验证覆盖 5～800 Hz；720 Hz 以上保留第五阶及其他未建模影响。

四阶频率边界由初值 ±10% 构造，并限制上界为 720 Hz，防止遗漏的高频模态替代目标模态。边界活动会被报告，不能据此认定频率无偏。

频响权重为 `W_ii=1/max(abs(G_i),floor_i)`，每通道幅值下限为该通道最大幅值的 `1e-3`。它不是逆协方差。第 4 步采用加性参数处的频响最小二乘局部信息矩阵平方根，保存为 `information_metric`。现有单份 FRF 不足以支持论文中基于重复实验的统计协方差解释。

## 输出和验证

- `results/identification_result.mat`：四阶段结果、配置、频响、诊断及指标。
- `results/identified_modes.csv`：四阶频率、阻尼比。
- `results/modal_validation.png`：16 通道频响及复数残差，绘图开启时生成。
- `results/修复验证报告.md`：本次实际运行数据与限制。

```matlab
tests = runtests('tests');
assert(all([tests.Passed]));
```

测试包括多刚体重构、非方形尺寸、雅可比有限差分、加权投影和收敛、无噪声四阶 RIV 恢复、刚体相位、LRM 瞬态、频率边界、有效频点和完整数据传递。原来的 `test.m`、`test1.m`、`test2.m` 是历史实验脚本，不是这套回归测试的入口。

修改前完整备份位于项目根目录 `backup_before_fix_20261012.zip`。
