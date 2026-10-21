# 新紧固 z=0 数据：动态 3+1 解耦入口审计

日期：2026-10-17。只读核对预处理、结构校准评价代码、CSV，以及主代理已读取的 MAT 字段。未修改原始 MAT、基线模型或评价脚本，未启动共享 MATLAB。Python 没有 h5py，不能直接读 v7.3；下列字段和尺寸由主代理 MATLAB 读取确认，建议用纯数组 -v7 导出后再交给 scipy。

## 1. 最小加载入口

~~~matlab
R = load('sweep/structured_validation_20261017/trial_z0.mat');
r = R.result;
f = r.f;                 % 1771 x 1, Hz
Htr = r.H_training;       % 1771 x 16, 已恢复记录边界单位
Hck = r.H_validation;     % 1771 x 16, 同记录后半开发检查 FRF
Hsrivc = r.baseline_H;    % 1771 x 16, 自由 SRIVC 预测，已含固定时延

E = load('sweep/structured_validation_20261017/evaluation_z0.mat');
e = E.evaluation;
Py = e.nominal_logical_output; % 3 x 4，每行单位范数
Pu = e.nominal_logical_input;  % 4 x 3，每列单位范数
ny = e.band_static_ny;         % 4 x 1，强静态基线
nu = e.band_static_nu;         % 4 x 1，强静态基线
Ty0 = [Py; ny.'];
Tu0 = [Pu, nu];
~~~

原始入口也可用 sweep/srivc_retest_20261016/data_z0.mat 的顶层 data，其字段为 f、trainingFRF、validationFRF、weights、training_indices、validation_indices、input_scale、output_scale、diagnostics。

model_z0.mat 的顶层 result 包含 selected.fit 和 selected.H。selected.H 是完整含时延预测，selected.fit.prediction 是给定时延补偿后的有理部分，不能混用。

## 2. 训练与检查划分

prepare_retest_data.m 丢弃首秒，然后把剩余连续记录均分：

| 项目 | z=0 数值 |
|---|---:|
| 原始样本数 | 303646 |
| 采样频率 | 5000 Hz |
| 原始持续时间 | 60.729 s |
| 训练 MATLAB 索引 | 5001:154323 |
| 检查 MATLAB 索引 | 154324:303646 |
| 每半样本数 | 149323 |
| 每半 FFT 周期长度 N/fs | 29.8646 s |
| 每半 FFT 频率间隔 | 0.0334844598622 Hz |
| LRM 91 点窗口首尾跨度 | 3.01360138760 Hz |

两半分别去均值、分别估计 LRM FRF。两半都用训练输入标准差做数值归一化，最后恢复单位。频率网格是 3:0.1:24.9 Hz 的 220 点，加上 25:0.5:800 Hz 的 1551 点；不能默认全带等间隔。

**后半已被多轮开发查看。** 本轮只能称“记录内留出开发检查”，不能称独立盲测或独立重测。LRM 邻近频率窗口重叠，1771 个频点不是 1771 次独立实验。

## 3. 输入倍率与 FRF 单位

新记录输入 SaturationU_0__0_ 至 SaturationU_0__3_，输出 y_sensor_0__0_ 至 y_sensor_0__3_，第一列时间、第二列信号。

用户确认多乘了

\[
k=\frac{10/6}{12.9}=0.129198966408269.
\]

预处理已经使用

\[
U_{\rm corrected}=U_{\rm recorded}/k=7.74U_{\rm recorded},
\qquad H_{\rm corrected}=kH_{\rm raw}.
\]

H_training/H_validation 已包含修正，不要重复乘除 7.74。LRM 内部 us=std(U_training)、ys=1e-6 仅用于数值缩放，存储时已按

~~~matlab
H(k,j+(u-1)*4) = theta(u)*ys/us(u);
~~~

还原。不要再次乘 ys/us。ys=1e-6 不是传感器单位已经证实为米的证据。当前单位是“传感器记录量／已修正电机指令记录量”，不能擅标 m/N。

## 4. 通道顺序和块标签

第 j 个输出、第 k 个输入列索引为 j+4(k-1)，即 output-fast。

~~~matlab
Gk = reshape(Htr(k,:),4,4);
Pk = Ty0*Gk*Tu0;
~~~

~~~python
Gk = Htr[k, :].reshape(4, 4, order="F")
Pk = Ty0 @ Gk @ Tu0
~~~

逻辑顺序 [Rx,Ry,Z,N]：

\[
P=\begin{bmatrix}A&b\\c&d\end{bmatrix},
\quad A\in\mathbb C^{3\times3},\
b\in\mathbb C^{3\times1},\
c\in\mathbb C^{1\times3},\
d\in\mathbb C.
\]

- A=P(1:3,1:3)：逻辑→逻辑。
- b=P(1:3,4)：N 输入→逻辑输出；output-fast 展平列 13:15。
- c=P(4,1:3)：逻辑输入→N 输出；展平列 [4 8 12]。
- d=P(4,4)：N→N；展平列 16。

## 5. 强静态基线准确重建

读取 sweep/srivc_retest_20261016/geometric_Ty_full.csv 和 geometric_Tu_full.csv。固定前三条逻辑方向，并逐行、逐列归一化：

\[
P_y(j,:)=T_y^{\rm geo}(j,:)/\|T_y^{\rm geo}(j,:)\|_2,\qquad
P_u(:,k)=T_u^{\rm geo}(:,k)/\|T_u^{\rm geo}(:,k)\|_2.
\]

只在 25–80 Hz **训练 FRF** 上，使用梯形积分权 q，构造：

\[
A_y=\sum_kq_k\Re\{(G_kP_u)(G_kP_u)^H\},
\qquad
A_u=\sum_kq_k\Re\{(P_yG_k)^H(P_yG_k)\}.
\]

对两矩阵取实对称部分后，分别取最小特征值对应的实单位特征向量 ny、nu，与名义 N 向量符号对齐。形成 Ty0=[Py;ny^T]、Tu0=[Pu,nu]，冻结后在检查数据评分。

这分别精确求解：

\[
\min_{\|n_y\|=1}\sum_kq_k\|n_y^TG_kP_u\|_2^2,\qquad
\min_{\|n_u\|=1}\sum_kq_k\|P_yG_kn_u\|_2^2.
\]

它没有使用结构模型，也没有时延优化。既有 z=0 的 25–80 Hz 检查中，相对名义几何矩阵，N→逻辑降低约 5.52 dB，逻辑→N 降低约 11.17 dB。**新动态方法应和该强基线比较，不能只和原几何矩阵比较。**

## 6. 尺度规范及实机坐标

之前评价是每端口单位范数。原几何矩阵的范数：

~~~text
Ty 行范数：[2.94901529, 5.27155939, 0.50000741, 1.0]
Tu 列范数：[5.11497293, 7.04998586, 0.50500670, 1.0]
~~~

名义单位端口完整矩阵条件数约为输出 1.05742、输入 1.15231。逻辑方向只是各自单位范数，并非严格正交，不能把 Py*Py^T 或 Pu^T*Pu 当作单位矩阵。更新 N 后仍须检查完整 Ty0、Tu0 的可逆性及条件数。

若单位坐标设计左右动态器 Lu、Ru，恢复原逻辑尺度时设

\[
T_y^{\rm phys}=D_yT_y^{\rm unit},\qquad
T_u^{\rm phys}=T_u^{\rm unit}D_u,
\]

则应使用

\[
L_{\rm phys}=D_yL_uD_y^{-1},\qquad
R_{\rm phys}=D_u^{-1}R_uD_u.
\]

规范化评价矩阵不能未经换算就代替保留 Rx/Ry/Z 单位及控制器增益的实机矩阵。

## 7. 动态验证关键检查

1. Htr、Hck、baseline_H 都已含实际观测边界的时延。直接使用完整频响设计 FIR，不要无说明先去一次时延、部署却只用 FIR。
2. 系数求解与候选选择只用训练数据，检查仅评分。既有后半已被开发查看的事实须保留。
3. 预先固定设计频带、检查频带。25–80 Hz 的静态最优不保证 100–170 Hz 或第一模态附近有效。
4. 保存 b、c 的复 RMS、A 的复响应变化、d 的幅相和共振窗口响应；不能只展示幅值变小。
5. 输入放大检查 Tu0*R(z)，输出噪声放大检查 L(z)*Ty0；只看 R、L 在逻辑坐标的范数未必代表实际执行器/传感器放大。
6. 没有实测噪声协方差时只能报告范数上界或白噪声假设下的量，不能声称实测噪声已证明不变。
7. 两侧理想消元通常将 d 变为 d-cA^{-1}b。Schur 补不等于原 NN 目标，交叉下降必须与柔性目标保留同时检查。
8. 既有频带 RMS 是

\[
{\rm RMS}_B(H)=
\sqrt{\frac{\int_{f_{\min}}^{f_{\max}}\|H(f)\|_F^2\,df}
{f_{\max}-f_{\min}}}.
\]

它没有再除以块内元素数；同块、同频带、同尺度的比值可比，不同大小块的数字不是每通道 RMS。积分要用真实频率网格及梯形权。

9. 原 CSV 的 change_dB_vs_nominal=20log10(new/nominal)，负值表示改善。若报告“降低量”则取负号。若与强静态比较，分母必须换成强静态响应。

## 8. 来源

- sweep/srivc_retest_20261016/prepare_retest_data.m：输入倍率、记录划分、LRM 和恢复单位。
- sweep/srivc_retest_20261016/raw_and_frf_audit.csv：N、fs、窗口跨度。
- sweep/structured_validation_20261017/evaluate_structured_trial.m：单位端口、静态特征向量基线、时延对齐、块选择和 RMS。
- sweep/structured_validation_20261017/实测检验结论.md：检查集使用历史及既有结果。
- 主代理本轮 MATLAB 读取：trial_z0.result 和 evaluation_z0.evaluation 的字段/尺寸。
