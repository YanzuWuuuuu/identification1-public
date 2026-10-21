from pathlib import Path
import json
import csv

P=Path(__file__).resolve().parent
read=lambda name:json.loads((P/name).read_text(encoding='utf-8'))
info=read('information_audit.json');bands=read('resonance_validation.json');coord=info['coordinates']
lines=['# 第二轮审阅建议核对与补充诊断（2026-10-14）','',
'结论：采用弱方向、控制器及共振频段诊断；不把换坐标写成消除弱辨识，也不预设带宽提高必然减少全频带有效信息。三项诊断已完成，未重跑主辨识、未设计新控制器、未新增硬件结论。','',
'## S1：等效斜率坐标、尺度与弱方向（C1/C2，正文第4章）','',
'同一主 N=32 所选解、数据、残差权重和 125 参数模型，仅作可逆变换 η=ℓ+a，保留 a、b。7 个 η、7 个 a、16 个 b 系数与原模型同维。不能只剩 (η,b) 后还声称只是换坐标。','',
'| 坐标 | 原残差尺度条件数 | 再列归一化的条件数 |','|---|---:|---:|',
f"| 原坐标 | {coord['old_condition']:.5g} | {coord['old_unit_column_condition']:.5g} |",
f"| η 坐标 | {coord['eta_condition']:.5g} | {coord['eta_unit_column_condition']:.5g} |",'',
f"同一相对响应目标的 C(JᵀJ)⁻¹C* 迹在两坐标下相对差仅 {coord['target_trace_relative_difference']:.3g}。该比较使用原固定拟合权重，只检验坐标不变性，不是置信覆盖验证。",'',
'展开 log[e^(−ℓs)/(1+as+bs²)]=−(ℓ+a)s+(a²/2−b)s²+(ab−a³/3)s³+O(s⁴)。最弱三个方向同时具有 δℓ≈−δa、δb≈aδa，即一阶斜率和二阶曲率都近乎不变。', '',
'| 弱方向 | σ | 斜率抵消比 ‖δa+δℓ‖/(‖δa‖+‖δℓ‖) | 曲率抵消比 |','|---|---:|---:|---:|']
for i,w in enumerate(coord['weak_directions']):
 lines.append(f"| {i+1} | {w['singular_value']:.4g} | {w['slope_cancellation_ratio']:.4g} | {w['curvature_cancellation_ratio']:.4g} |")
lines += ['',
'完整右奇异向量及逐通道 δa、δℓ、δb 保存在 coordinates.json。能量占比依赖参数尺度，不能把 log(b) 坐标的高占比解释成物理重要性排序。', '',
'因此，等效斜率坐标配合缩放值得作为后续求解预处理；数据本身的弱信息没有消失，不能只换尺度以通过绝对 score 阈值后宣称渐近正态已经验证。单个一阶链 t=35 μs、800 Hz 的幅值为 0.984875，相对同时间纯时延的相位差约 0.102106°；此例的数量级正确，但不是整个 125 参数模型条件数的唯一来源。', '',
'## S2：控制器扫描及有效信息（C1/C2，正文第4章、补充实现诊断）','',
'沿用原五阶对象、名义几何、5 kHz 采样、196 条 20:4:800 Hz 频线、每输入每线 0.1 N 峰值、100 nm 采样噪声、N=32。仅改变原控制律中的设计频率，微分滚降固定 120 Hz。设计频率不是实际交越频率。先检查完整离散闭环极点，失稳配置不计算可用实验的信息。','',
'信息使用真值处的均值 Jacobian、平稳频谱近似和完整输出谱协方差 Σ_y=S_yΣ_vS_y*。冻结真协方差后，可把白化均值灵敏度写为 Σ_v^(−1/2)δG S_u D。未计入协方差对参数的导数项，所以这里是均值信息/最优加权局部诊断，不是完整随机闭环似然的全部 Fisher 信息，也不是有限记录对角加权估计量的实测方差。','',
'相对目标取六个非参考输入/输出通道在 20、100、286.43、800 Hz 上的 δr/r，报告局部 RMS 标准差界；时延目标取输出与输入各三个相对时延差。两者均除以 N=32。所有 nuisance 列均保留，无截断弱方向；SVD 计算与 profile Schur 补核对相符。','',
'| 控制设计频率/Hz | 实际首次下降交越范围/Hz | 闭环谱半径 | 20 Hz 受控方向激励功率比 | 相对动态局部界/% | 相对时延局部界/μs |','|---|---|---:|---:|---:|---:|']
for r in info['scan']:
 cross=[v for v in r['first_downcross_hz'] if v is not None]
 cx=f'{min(cross):.1f}–{max(cross):.1f}' if cross else '未在扫描范围检出'
 if r['stable']:
  v=r['complex'];lines.append(f"| {r['controller_design_hz']} | {cx} | {r['closed_radius']:.7f} | {r['rigid_excitation_20hz_relative']:.5f} | {v['relative_response_rms_bound_percent']:.4f} | {v['delay_target_rms_bound_us']:.2f} |")
 else:lines.append(f"| {r['controller_design_hz']} | {cx} | {r['closed_radius']:.7f} | 失稳，不作有效实验评价 | — | — |")
lines += ['',
'20 Hz 受控方向的激励功率比从约 0.725 降到 0.00995，局部低频抑制确实存在。但完整频带的有效信息没有骤降：交越附近的灵敏度峰值、未被反馈压低的冗余输入方向，以及参数间的耦合都影响结果。固定注入 d 的总实际输入能量比随设计频率为 1.080、1.194、1.870，因此这不是实际执行器能量相同的比较，不能据此主张提高带宽改善辨识。100–200 Hz 设计频率的控制器已经失稳；不额外设计新控制律来凑齐扫描。','',
'按名义刚体右子空间投影，20–100 Hz 的平均激励功率比为 1.387、1.613、1.218；冗余方向投影的功率比均为 1。频带平均不能代替 20 Hz 单频的衰减。只有在正确完整协方差、该噪声入口及平稳近似下 S_y 才在均值信息内抵消；主求解器的对角权重、过程扰动或非稳态边界不能直接沿用这一简化。','',
'比例阻尼对照只在真实原模型上使用 90 参数局部切空间，保持机械极点、通道系数、数据和噪声完全相同。结果：','',
'| 设计频率/Hz | 复留数相对动态界/% | 实模态相对动态界/% | 复留数相对时延界/μs | 实模态相对时延界/μs |','|---|---:|---:|---:|---:|']
for r in info['scan']:
 if r['stable']:
  c=r['complex'];q=r['proportional'];lines.append(f"| {r['controller_design_hz']} | {c['relative_response_rms_bound_percent']:.4f} | {q['relative_response_rms_bound_percent']:.4f} | {c['delay_target_rms_bound_us']:.2f} | {q['delay_target_rms_bound_us']:.2f} |")
lines += ['',
'实模态先验在此真对象上减少 nuisance 自由度而增加信息，但没有检验一般阻尼失配，也不证明“吸收比例恒为一半”。相对时延局部界已达约 90–184 μs，和搜索域及真时延同量级或更大，意味着局部无约束正态近似不宜直接作实际置信区间。相对响应约 1% 的局部界与单独物理参数很弱可以同时成立。','',
'理想多输入频线整形应写为 D=(I+KG)U_target，而不是逐通道标量 1/|S_u|。使用真模型的理想诊断可精确恢复所选 U_target 的均值信息，但注入代价改变；15/30/60 Hz 设计下全频带 d 能量倍数为 0.951/0.984/1.882，最大单线峰值放大约 1.30/2.64/9.30。最大值只描述注入频线，不是总时域峰值或饱和保证。实际设计须使用先导估计并约束驱动力、位移、控制量及峰值。本轮不部署该整形。','',
'## S3：第一阶共振与整体模态频带（C3，正文第6章）','',
'复用全部已保存独立时域验证记录、所选模型、控制器及数字补偿，未重拟合。260–310 Hz 窗口实际包含 262:4:310 Hz，覆盖第一阶 286.43 Hz；此前 200–800 Hz 窗口实际包含 202:4:798 Hz。重新计算的原低频结果与旧日志一致。','',
'| 案例 | 第一阶频带 B1/B4 泄漏/nm | B4 柔性/真参数参考/nm | B4 噪声/nm | B1/B4 信噪比/dB |','|---|---:|---:|---:|---:|']
for e in bands:
 c=e['calibration'];s=c['static']['first_resonance'];j=c['joint']['first_resonance'];o=c['oracle']['first_resonance']
 lines.append(f"| {e['case']} | {s['rigid_rms_m']*1e9:.4f} / {j['rigid_rms_m']*1e9:.4f} | {j['flexible_rms_m']*1e9:.3f} / {o['flexible_rms_m']*1e9:.3f} | {j['noise_rms_m']*1e9:.4f} | {s['flexible_to_rigid_noise_db']:.3f} / {j['flexible_to_rigid_noise_db']:.3f} |")
lines += ['',
'主 N=32 在整个 202–798 Hz 模态频带，B1/B4 泄漏为 6.857/0.172 nm，B4 柔性约 730.043 nm、数字参考 730.179 nm；B4 噪声 25.710 nm，信噪比 28.770→29.065 dB（约 +0.294 dB）。第一阶频带约 +0.403 dB。它证明共振响应得到保留，但收益已经受噪声限制。','',
'22–98 Hz 的约 55 nm 柔性响应主要是模态尾部的准静态柔度，并非共振峰。该频段 +9.21 dB 的收益必须与共振频段约 +0.40 dB 分别报告。第一阶频带普通闭环输出接近数字参考，不替代先前双侧各阶投影的保留率；两者是不同任务指标。数字真参数参考也不是所有指标上的严格最小值，估计偏差可能偶然抵消数字实现误差。','',
'## 采用哪些意见，哪些不继续扩展','',
'- 采用：一阶/二阶分块误差结构作为 C3 表述；统计结论按条件性理论定位；补上共振验证；将随机种子与停止容差移到补充实现部分。','- 有条件采用：等效斜率与缩放是求解预处理候选，尚未重跑优化，不能声称已改善求解成功率。','- 不采用为事实：时延和滤波结构上完全分不开、闭环带宽提高必然使全频信息骤降、换坐标即可验证正态性、固定 1/2 或 98% 吸收比例。','- 暂不做：100–200 Hz 新控制器设计、完整激励最优化、一般摩擦/刚度建模、全套置信覆盖研究。它们会扩大本轮主线。','',
'证据：coordinates.json（含右奇异向量）、bandwidth_scan.json、information_audit.mat/json、resonance_validation.json；实现 run_information_audit.m、run_band_audit.m，父实现与旧结果保持不变。下一步返回论文整合与硬件入口验证。']
(P/'审阅建议核对_弱方向与闭环带宽.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
status={'status':'completed','route':'write','parent':'closed_loop_validation/outputs_structured',
 'S1':{'claim_links':['C1','C2'],'section_id':'4.6','paper_role':'supporting','takeaway':'eta plus column scaling reduces coordinate conditioning; target information unchanged.','comparability':'identical fit, data, weights and model; no refit'},
 'S2':{'claim_links':['C1','C2'],'section_id':'4.6 and supplementary implementation','paper_role':'auxiliary','takeaway':'low-frequency controlled directions are suppressed; full-band information not monotone; high-gain variants unstable','comparability':'fixed injected d, different actual force energy; full spectral covariance approximation differs from fit weights'},
 'S3':{'claim_links':['C3'],'section_id':'6.6','paper_role':'supporting','takeaway':'first resonance retained; SNR gain is noise-limited at about 0.40 dB','comparability':'same saved heldout records, no fit or selection changes'},
 'deferred':['reparameterized optimizer benchmark','new high-bandwidth controller','hardware validation','covariance coverage']}
(P/'campaign_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
print('Saved scoped review report and campaign status.')
