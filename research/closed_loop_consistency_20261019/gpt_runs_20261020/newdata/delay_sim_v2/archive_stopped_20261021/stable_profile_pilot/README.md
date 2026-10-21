# 外层数值修复的独立 pilot

这里只研究求解实现，不改变八模态、375 个参数、原边界和完整矩 score < 1e-6 的验收条件。所有结果与正式 `mc_qr_final` 隔离，不能挑选其中成功端点补入原统计。

## 修复内容

原外层再次形成线性正规块并截断，H 窗口的 352 个数值可分辨线性方向被削成 332 个；极大的窗外模态步长还将所有参数的共同步长压到约 5e-6。新的 `code/stable_profile.py` 使用

\[
C=Q_z^\top Q_x,\quad
\Pi=Q_xC^{-1}Q_z^\top,\quad
\bar J=(I-\Pi)J_n,\quad
\bar Z=(I-\Pi^\top)Z_n,\quad
S=\bar Z^\top\bar J.
\]

线性部分使用原高矩阵的 QR/SVD，截断仍为 1e-13；不对线性 Gram 块作伪逆。外层用有界 TRF，根据全部剩余矩选择步长，不再按最大单个分量缩放所有参数。

`projected` 模式的雅可比是冻结工具变量的 GN 近似。每次真正的残差评估仍重新估计线性参数并更新参考工具；该模式不是完整牛顿法。`central` 模式对重新剖面后的矩作中央差分，包含工具变量和线性参数更新。

最终检查调用原始未压缩 375 维矩函数。保持 score < 1e-6、内层线性 score < 1e-9 和原贴边判定，不以 TRF 的内部成功标志代替验收。

## 已执行的命令

在 `delay_sim_closeout` 目录下运行，Python 为本机 Python 3.11：

```text
python -B code/stable_profile.py --periods 10 --band H --method M3 --max-nfev 12 --max-profile-calls 24 --max-seconds 360 --tag v1
python -B code/stable_profile.py --periods 20 --band W --method M1 --max-nfev 12 --max-profile-calls 24 --max-seconds 360 --tag v1
python -B code/stable_profile.py --periods 10 --band H --method M3 --max-nfev 40 --max-profile-calls 40 --max-seconds 360 --tag v2
python -B code/central_profile_finish.py
```

前三次运行使用的原始数值代码已保存在 `code/stable_profile_v1_frozen.py`，各 `protocol.json` 记录代码哈希。随后 `stable_profile.py` 仅增加雅可比计算前后的时间检查，数值算法与阈值不变。已启动的中央差分进程也使用启动时加载的 v1 数值代码。

## 结果保存规则

- 每次实际 profile 调用写入 `eval_*.json`，包括失败和未接受的试探点。
- 每次投影雅可比写入 `jac_*.json`，保留逐输出秩、投影误差和 Schur 奇异值。
- `result.json` 是最终端点；不得把中间最小 score 点冒充最终点。
- `central_P10_r000_H_M3` 另外保存 `fdcheck_*.json` 与 `difference_step_check.json`。中央差分步长在外层缩放坐标中为 1e-3，并与 3e-3 作方向核查。
- 中央差分预算为 150 次实际求解 profile，加 12 次独立步长诊断；不把 scipy 的 `nfev` 冒充真实 profile 调用数。

时间预算采用安全计算边界检查，无法中断正在执行的 LAPACK。W/v1 在 360 秒后完成了在途 profile，又计算了一次投影，最终为 410 秒；该超时如实保留。后续版本增加雅可比前后检查，避免过期后再启动额外投影，但单次在途计算仍可能略超时。

## 数学与实现核查

独立审计保存在：

- `audit/schur_rank_check.json`：旧外层截断和共同缩步问题；
- `audit/stable_schur_equivalence.json`：新投影与稳定区间的完整 Schur 等价核查；
- `audit/stable_endpoint_check.json`：H/v2 剩余矩所在方向；
- `audit/正式QR批次独立实现审计.md`：原实现的独立审计。

W 的投影 Schur 与完整秩正规消元相对差约 3.10e-12，投影正交残差约 1.9e-15。该核查证明消元公式实现一致，不能证明所有参数实际可辨识，也不取代最终完整矩验收。
