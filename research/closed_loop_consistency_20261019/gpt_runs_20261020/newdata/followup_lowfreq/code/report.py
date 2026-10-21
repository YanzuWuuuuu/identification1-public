from common import *
from summarize import csvout
def fmt(v,n=3):
    if isinstance(v,(float,np.floating)):return f'{v:.{n}f}'
    return str(v)
def table(head,rows):
    return '\n'.join(['| '+' | '.join(head)+' |','| '+' | '.join(['---']*len(head))+' |']+
                     ['| '+' | '.join(str(v) for v in r)+' |' for r in rows])+'\n'
def text(s):parts.append(s.strip()+'\n')
def tri(v,n=3):return '/'.join(fmt(x,n) for x in v)
def flag(v):return '是' if v else '否'
def main():
    global parts
    parts=[]
    bs={r:rd(HERE/'analysis'/f'{r}_noise_bias.json') for r in IDS}
    jk={r['id']:r for r in rd(HERE/'analysis/jackknife.json')}
    tests=rd(HERE/'analysis/paired_tests.json')
    res=rd(HERE/'analysis/residual_checks.json')
    dec=rd(HERE/'analysis/decisions.json')
    scale=rd(HERE/'analysis/scaling.json')
    ver=rd(HERE/'analysis/verification.json')
    text("""# 低频偏差预测与残差检验报告

日期：2026-10-20。数据：已有 A 与 C++++ 六条记录；仅 3–25 Hz；本次未采集新数据。

## 结论先行

固定模型、频带和 0.5 Hz 噪声估计下，**本轮三项预设验证均未整体通过**。这不是“所有计算失败”：82 次 jackknife 重估都达到原停止规则；问题是噪声层面的模型充分性和部分偏差预测仍与数据不相容。

- 第二阶阻尼的噪声替换预测，在四条 nominal/highgain 关键记录中通过 3 条；A_nominal 仍明显低估，不能用换噪声谱宣布问题解决。
- 参考 IV 同一记录残差的 Q/E[Q] 为 15.6–82.6；原加权参考 IV 为 7.7–19.2，全部超过预设门槛 2。残差存在当前条件噪声模型解释不了的成分。
- lowgain 跨运行检验中约 2600 倍的异常，96.95% 的 Q 位于 15–17 Hz，3–6 Hz 仅占 0.70%；“主要是低频漂移”的猜测不受这一检验支持。
- 非激励点残差中存在很强的“可由原瞬态列解释”的分量。把其总能量直接叫作加激励后的随机噪声增加，证据不足。新增瞬态诊断不替换主谱，也不改变判据。
- 线性 RMS 匹配的初始扰动系数为 lowgain≈0.856、nominal≈1.064、highgain=1.5；D 组显示线性缩放误差可达约 23%，不能保证这些设置一次就匹配响应。

这些结果不构成对条件一致性定理的反例；它们表明，不能把该定理的正确模型、噪声和求解条件未经检验地当成这六条有限记录已经满足的事实。

## 0. 复用设置与审计

参数顺序完全保留原实现：96 个模态 C0/C1 系数、32 个 D0/D2 系数、3 个 log(ω)、3 个 log(ζ)、24 个逐模态瞬态系数，共 158 个。时延使用原有效通道时延减 100 μs；支承模态精确 ZOH；瞬态分子为 τ0+τ1 z⁻¹。未新增柔性模态、未变频带、未重新拟合原端点。

M1 使用 3–25 Hz 全部 DFT 点；M3/M3w 使用激励谱线。M3 工具为 Û=(I+KĜ)⁻¹D。M3w 保留原来冻结的 16 邻点、α=0.1 收缩协方差权；本轮 0.5 Hz 原始谱仅用于噪声预测和残差检验，不能拿它重新定义原 M3w。
""")
    text(table(['记录','可用周期 P','全长点数 N','激励谱线数'],[
        [r,bs[r]['P'],bs[r]['N'],133] for r in IDS]))
    text(f"""写入前哈希快照与收尾比较：{ver['protected_files_checked']} 个已有文件，改变 {len(ver['changed'])}，丢失 {len(ver['missing'])}；四份模型/求解器/时延基线副本逐字节一致。新输出全部位于本目录。

复现原数据入口的最大误差为 {ver['numerical_checks']['source_spectra_max_error']:.1e}；复现原 B 谱偏差预测的最大绝对差为 {ver['numerical_checks']['original_B_prediction_max_abs_percent_error']:.2e} 个百分点。原报告、数据、拟合结果与理论文档未被改动。

## 1. 用同记录噪声谱重做偏差预测

### 1.1 DFT、谱与协方差尺度

采用未归一化 DFT：V_N(k)=Σₙv[n]exp(−j2πkn/N)。在 1–27 Hz 获取非激励点，选择规则 k mod P≠0；在各中心频率 f 上，使用 |f_j−f|≤W_s/2 的全部非激励点。0.5/1.0 Hz 是**总窗口宽度**。

v̂_j=Y_j−Ĝ_M3(f_j)U_j−T_M3(f_j)，

Φ̂_v(f)=(1/(N n_f))Σ_j v̂_j v̂_jᴴ。

这里 Φ 的单位与原 B 谱相同，即 N⁻¹ᐟ² DFT 的协方差；双边 PSD/Hz 为 Φ/fs，正频率单边 PSD 为 2Φ/fs。全长 DFT 残差的复协方差为 NΦ；单周期 DFT 为 P0Φ；P 周期平均 DFT 为 P0Φ/P。本轮没有再除或多乘一次周期数。

原 B 谱自身采用既有平滑/收缩；新谱按任务直接取局部二阶矩，**不增加对角加载或收缩，不减去局部频率均值**。所有新 4×4 谱矩阵 Cholesky 均成功。

偏差计算保留实际 D、恢复 K、全部 158 列和原 M1 频点，复用 M⁻¹b 与 (M−C_curv)⁻¹b 两种线性化。下表主结果使用后者，转换统一为 100[exp(Δlogζ)−1]；实测配对差是 100(ζ_M1/ζ_M3−1)，不是两个阻尼百分数的直接相减。
""")
    prows=[]
    for rid in IDS:
        b=bs[rid]
        for j in range(3):
            k=j+3
            prows.append([rid,f'ζ{j+1}',fmt(b['observed'][k]),fmt(b['old_B']['pred_MminusC_percent'][k]),
                 fmt(b['windows']['0.5']['pred_MminusC_percent'][k]),fmt(b['windows']['1.0']['pred_MminusC_percent'][k])])
    text(table(['记录','参数','实测配对差 %','原 B 预测 %','同记录 0.5 Hz %','同记录 1.0 Hz %'],prows))
    text("### 1.2 第二阶比值与窗口敏感性")
    text(table(['记录','实测/预测 (0.5 Hz)','预测/实测 (0.5 Hz)','1.0−0.5 Hz 预测差（百分点）','预设比值通过'],[
        [r,fmt(bs[r]['windows']['0.5']['observed_over_prediction'][4]),fmt(bs[r]['windows']['0.5']['prediction_over_observed'][4]),
         fmt(bs[r]['windows']['1.0']['pred_MminusC_percent'][4]-bs[r]['windows']['0.5']['pred_MminusC_percent'][4]),
         flag(.67<=bs[r]['windows']['0.5']['prediction_over_observed'][4]<=1.5) if r in ('A_nominal','A_highgain','C_nominal_pppp','C_highgain_pppp') else '不在任务1四条判定内'] for r in IDS]))
    text("""任务1整体不通过：A_nominal 的预测/实测为 0.386，另外三条关键记录分别为 1.156、0.878、0.890，落在 [0.67,1.5] 内。1.0 Hz 档的第二阶预测比0.5Hz低约11%–20%，不是完全可忽略的敏感性；但四条关键记录仍是同样的3条通过、A_nominal不通过。1.0 Hz 档仅作为敏感性分析，不据此选择“更好”窗口。

为了防止把曲率修正与谱替换混为一谈，以下同时保留 M⁻¹b 结果；斜线顺序始终为 ζ1/ζ2/ζ3。""")
    text(table(['记录','原 B：M⁻¹b %','0.5 Hz：M⁻¹b %','1.0 Hz：M⁻¹b %'],[
        [r,tri(bs[r]['old_B']['pred_M_percent'][3:]),tri(bs[r]['windows']['0.5']['pred_M_percent'][3:]),
         tri(bs[r]['windows']['1.0']['pred_M_percent'][3:])] for r in IDS]))
    text("### 1.3 多余谱能量、裁切周期和瞬态诊断")
    extra=[]
    for rid in IDS:
        a=np.load(HERE/'spectra'/f'{rid}_noise.npz');f=a['f']
        ex=np.maximum(np.trace(a['Phi_0.5']-a['B'],axis1=1,axis2=2).real,0)
        low=100*ex[f<6].sum()/ex.sum();modal=100*ex[(f>=9)&(f<18)].sum()/ex.sum()
        cuts=bs[rid]['extra_drop'];frozen=[x['frozen']['band'] for x in cuts];reb=[x['reboundary']['band'] for x in cuts]
        extra.append(dict(id=rid,excess_3_6_percent=low,excess_9_18_percent=modal,
                          frozen_ratios=frozen,transient_refit_ratios=reb))
    js(HERE/'analysis/noise_attribution.json',extra)
    text(table(['记录','全带残差谱/B','正的超额能量中 3–6 Hz %','9–18 Hz %','固定瞬态：多丢0/1/2周期','仅重估原瞬态：多丢0/1/2周期'],[
        [x['id'],fmt(x['frozen_ratios'][0]),fmt(x['excess_3_6_percent'],1),fmt(x['excess_9_18_percent'],1),
         tri(x['frozen_ratios']),tri(x['transient_refit_ratios'])] for x in extra]))
    text("""表中全带比为 ∫trΦ_new/∫trΦ_B；超额能量位置使用 max(tr(Φ_new−Φ_B),0)，不是把比值特别大的、但 B 本身能量很小的频点当作主导能量。

裁切后仍保持原 M3 的所有参数是主诊断。由于有限记录的瞬态系数取决于记录边界，另做一个明确分开的诊断：固定所有传递函数参数，仅在 3–25 Hz 非激励点对原有 24 个瞬态系数做线性最小二乘；维数、模型和频带均不改变。对“未额外裁切”也做同样诊断，避免把重新估计的收益误认成多丢周期的收益。此诊断**没有进入任何主预测、Q 或预设判定**。

逐条判断如下：

- **A_lowgain：**超额能量主要在支承频段；固定瞬态裁切后比值上升，重估瞬态后约 0.37–0.41，低于 B。不能断言噪声增加；原瞬态估计/有限记录边界分量明显。
- **A_nominal：**原约 6.14 倍降至仅重估瞬态后的约 1.16 倍，且多丢周期后仍约 1.13–1.17。所谓“六倍噪声”很大部分可由原瞬态列解释；不支持简单视为宽带随机噪声六倍增加。
- **A_highgain：**重估瞬态后仍约 1.59–1.67 倍，存在不能由这一瞬态诊断去掉的增加；现有数据不能区分有色噪声、模型误差和非平稳成分。
- **C_lowgain++++：**重估瞬态后约 0.80–0.86 倍，原 1.42 倍同样不能直接称为激励相关噪声增加。
- **C_nominal++++：**重估瞬态后约 1.55–1.71 倍，可能同时有边界分量与剩余有色残差；裁掉周期没有单调改善。
- **C_highgain++++：**重估瞬态后仍约 3.40–4.15 倍，并非低频 3–6 Hz 主导；变化不单调，无法仅归因于启动瞬态。

总体：89%–97% 的正超额能量位于 9–18 Hz；3–6 Hz 只占约 0.5%–1.1%。**不支持低频漂移是本项主因；支持存在“瞬态列可解释分量”，但不满足“多丢周期后明显衰减”的启动瞬态判据；剩余部分也不是已证实的白噪声/宽带噪声增加。** 六条主记录激励系数相同，不能据此建立噪声增加与激励幅值的因果关系。

每条记录的四个对角比值、迹比值及裁切对照均见对应 `figures/*_noise_ratios` PNG/PDF；原频率点、完整谱矩阵和裁切谱保留在 `spectra`。

## 2. 配对差的不确定度与吻合检验

### 2.1 周期 jackknife 实现及数值状态

对原可用的 P 个完整周期逐一删除一个，按原时间顺序拼接剩余周期，重建长度 (P−1)P0 的 DFT 和新的激励谱线索引。每种方法在该记录所有折中固定使用原文件保存的 initial_beta；M1 与 M3 保留各自在原流程中的初值，不用某一折端点给另一折热启动。停止规则仍为最多80步、模态步长阈值1e−5、步长上限0.05。

se_JK=sqrt[(P−1)/P Σᵢ(δ_(i)−δ̄)²]。未删除失败或异常折后重新汇总。

删除中间周期会形成一个拼接边界；保留原全记录瞬态模型而不添加新的边界参数，是本次按任务执行的约束。因此该 SE 是指定删除周期程序的重估离散；其理想 jackknife 解释仍要求周期足够稳定，不能自动覆盖运行间漂移或周期性非线性畸变。
""")
    text(table(['记录','P','M1/M3 原规则失败','M1/M3 严格score未通过','M1/M3 平均迭代数'],[
        [r,jk[r]['P'],f"{jk[r]['failed_M1']}/{jk[r]['failed_M3']}",
         f"{jk[r]['strict_failed_M1']}/{jk[r]['strict_failed_M3']}",
         f"{np.mean([v['M1']['solver_info']['iterations'] for v in jk[r]['folds']]):.1f}/{np.mean([v['M3']['solver_info']['iterations'] for v in jk[r]['folds']]):.1f}"] for r in IDS]))
    text("""82 次重估按原规则均收敛；22 次未达到更严格 score<1e−6。两种状态同时报告，不能把“未到迭代上限”写成“严格方程残差全部达标”。

### 2.2 检测门槛、z 与联合 sandwich

仅在 |预测|>2se_JK 时报告正式 z=(实测−预测)/se_JK；其他行标为“低于检测限”。完整 CSV 还保留未用于检验的代数 z，便于审计。此 z 不计预测谱估计误差、K 估计误差和线性化近似误差，不能解释为已校准的无条件标准正态检验。

联合 sandwich 在同一 M3 端点局部线性化：M1 保留 Gaussian 线性+中心化二次噪声矩及 M−C_curv，M3 保留线性参考 IV 矩，激励谱线上显式计算共享噪声交叉协方差；不能把 M1/M3 当成独立估计再相加。该近似的 PSD、模型和 K 均固定，不包括权重/谱估计和模型失配误差。
""")
    text(table(['记录','参数','实测 %','预测 %','JK SE','sandwich SE','正式 z','实测/预测','门槛/结论'],[
        [r['id'],r['parameter'],fmt(r['observed']),fmt(r['prediction']),fmt(r['jackknife_se']),fmt(r['sandwich_se']),
         fmt(r['formal_z'],2) if r['detectable'] else '—',fmt(r['observed_over_prediction'],2),
         ('|z|<2' if abs(r['formal_z'])<2 else '|z|≥2') if r['detectable'] else
         ('预测低于检测限；实测≤2SE' if r['below_detection_observed_within_2se'] else '预测低于检测限；实测>2SE')]
        for r in tests if r['parameter'].startswith('zeta')]))
    text("""第二阶有5条记录通过预测效应检测门槛，其中4条 |z|<2，A_lowgain 的 z≈−3.74 未通过。A_nominal 的预测 4.76% 小于2SE≈5.24%，按预设规则不作正式吻合检验；其实际差 12.35% 明显大于2SE，不能因此当作“吻合”。

第一、三阶低预测效应也并非全部满足实测≤2SE：A_nominal 第三阶、A_highgain 第一/三阶、C_nominal 第三阶、C_highgain 第一/三阶，共6项仍有超出。任务2不能整体通过。sandwich 与 JK 的差距也说明两者覆盖的误差来源和有限样本稳定性不同，不宜择小或择大选择有利的误差棒。

### 2.3 C_highgain++++ 第一阶 8.6% 配对差
""")
    ca=rd(HERE/'analysis/C_highgain_first_mode_audit.json')
    text(table(['方法','ζ1 (%)','迭代数','严格 score','原规则收敛'],[
        [met,fmt(ca['original'][met]['damping_percent'][0],6),ca['original'][met]['iterations'],
         f"{ca['original'][met]['score_normalized']:.3e}",flag(ca['original'][met]['converged'])] for met in ('M1','M3')]))
    text(table(['删除的周期（从0计）','第一阶配对差 %','M1 score','M3 score'],[
        [r['omit'],fmt(r['paired_zeta1']),f"{r['M1_score']:.2e}",f"{r['M3_score']:.2e}"] for r in ca['folds']]))
    text(f"""两端点都达到严格 score；求解器没有显式参数上下界，因此不存在“贴着设置的阻尼边界”这一解释。全部折的模态顺序保持不变。删除周期后的差为 {tri(ca['paired_range'])}%，SE={ca['paired_se']:.3f}%，新预测={ca['prediction']:.3f}%，未出现单折发散造成8.6%的证据；但其中3个 M3 折未通过更严格score，仍应保留数值限定。

可写成“本记录中稳定出现的方法间差异，未发现明显贴边或未收敛解释”。不能进一步把它断言为物理阻尼真实变化或纯反馈噪声偏差：二者来自同一记录，而且本记录残差充分性检验也未通过。

## 3. 同一记录内残差检验

### 3.1 统计量与线性化计算

每频率的实堆叠顺序为 [Re(y1:4),Im(y1:4)]，133条谱线给 n=1064。复协方差 NΦ 在 proper 复高斯约定下的实块为

C_R=(N/2) [[ReΦ,−ImΦ],[ImΦ,ReΦ]]。

采用任务给定 R=I−J(ZᵀWJ)⁻¹ZᵀW，保留全部158参数。以 C_R=LLᵀ 白化后计算 R_w=L⁻¹RL，Q=‖L⁻¹e‖²，E[Q]=‖R_w‖_F²。IV 的 R 通常为斜投影，**不能直接把 E[Q] 改为 n−p=906**。使用实测 U 的 J、参考工具的 Z 和原实际权；未将 Z 偷换成 J。

每个端点用200组 N(0,I) 实高斯向量生成 e*=RLξ，等价于所需 proper 复高斯噪声；Q* 用于核查条件线性化 E[Q]，不重新拟合非线性模型，也不包含谱估计不确定度。
""")
    text(table(['记录','方法','Q','E[Q]','Q/E[Q]','bootstrap均值','bootstrap均值MC-SE','Q* 95%分位'],[
        [r['id'],r['method'],fmt(r['Q'],1),fmt(r['EQ'],1),fmt(r['ratio'],2),
         fmt(r['bootstrap_mean'],1),fmt(r['bootstrap_mean_se'],1),fmt(r['bootstrap_q95'],1)] for r in res]))
    text(f"""12项检验均大于2。200次模拟均值与解析期望的最大距离为 {ver['numerical_checks']['bootstrap_mean_max_absolute_MC_SE']:.2f} 个 Monte Carlo 均值标准误；直接 Frobenius 计算与 QR 稳定化低秩迹公式的最大相对误差为 {ver['numerical_checks']['residual_trace_identity_max_relative_error']:.2e}，R 对 J 的消去误差最大 {ver['numerical_checks']['projector_annihilation_max_relative_error']:.2e}。

M3w 某些 Q/E[Q] 较低不等于拟合更准：其斜投影噪声放大使 E[Q] 明显增加，有些情况下原始 Q 本身反而更大，必须把 Q 与分母一起看。所有结论均保留“给定估计谱、线性化算子”的条件限定。

### 3.2 每2Hz频段分解
""")
    seg=[]
    for r in res:
        for q in r['frequency_bands']:seg.append(dict(id=r['id'],method=r['method'],**q,fraction_total_Q=q['Q']/r['Q']))
    csvout('residual_frequency_bands.csv',seg)
    text(table(['记录','方法','Q贡献最大的三个2Hz频段（占总Q）'],[
        [r['id'],r['method'],'；'.join(f"{x['lo']:g}–{x['hi']:g} Hz ({100*x['Q']/r['Q']:.1f}%)"
            for x in sorted(r['frequency_bands'],key=lambda q:q['Q'],reverse=True)[:3])] for r in res]))
    text("""全部频段 Q、E[Q] 和比值在 `analysis/residual_frequency_bands.csv`，对应图为 `figures/residual_frequency_bands`。未剔除高残差频段。不同记录和方法的超出频段不同，不支持用一个“3–6 Hz漂移”解释全部超出。

残差可同时包含结构/背景近似误差、周期性非线性畸变、瞬态估计误差、非平稳噪声和已估谱误差；仅凭本检验不能唯一分解。按用户约束，本轮不通过改模型、缩频带或另选谱来降低统计量。

### 3.3 lowgain 的跨运行异常
""")
    cr=rd(HERE/'analysis/lowgain_crossrecord.json')
    text(f"""精确复算原 A 模型预测四条 C 记录的统计量，Q/理想期望={cr['ratio']:.2f}。其中15–17Hz占96.95%，3–6Hz占{100*cr['fraction_3_6']:.2f}%。这项异常主要位于第三支承模态附近，不能再描述为“主要是3–6Hz漂移”。
""")
    text(table(['记录','四路3–6Hz RMS（原记录单位）','3–6Hz占3–25Hz功率比例 %'],[
        [r,tri(cr[r]['low_3_6_rms'],8),tri(100*np.array(cr[r]['low_3_6_power_fraction']),3)]
        for r in ('A_lowgain','C_lowgain_pppp')]))
    text("""两条 lowgain 记录的3–6Hz输出功率相近。注意：四条C组成的跨运行检验与仅C++++这条记录并不相同；A/C++++的低频功率相似也不能排除其他C方向记录发生局部模态变化。

本轮同记录检验采用该记录自己的拟合、瞬态、同记录谱和IV残差算子；不包含 A 与 C 运行之间的预测差。原跨运行检验使用 A 模型、四条 C 的周期均值和 B 谱，既有运行条件变化，又有谱来源不同，不能直接把二者比值当作单一误差成分的相减。

## 4. 补采的激励缩放系数

### 4.1 仅降低激励的线性初始值

用3–25Hz全带单边功率积分的传感器RMS，a_{c,j}=1.5 R_high,j/R_c,j。四路无法用一个标量同时精确匹配，统一建议取四路比例的最小值（且不超过1.5），保证在线性假设下每一路都不高于目标；另列总能量匹配值用于比较。
""")
    text(table(['控制器','输出1','输出2','输出3','输出4','统一建议','总能量匹配'],[
        [r['id']]+[fmt(v,4) for v in r['per_output_gain']]+[fmt(r['disturb_gain'],4),fmt(r['aggregate_gain'],4)]
        for r in scale['scalings']]))
    text(table(['记录','3–25Hz四路RMS（原记录单位）'],[
        [r['id'],tri(r['original_response_rms'],8)] for r in scale['scalings']]))
    text("""这些是补采初始设置，不是对硬件已执行的修改，也不是保证最终匹配的系数。

### 4.2 D组线性假设核对

以A_nominal的1.5工况线性外推，其余D组的实际RMS相对外推值如下。对实际D谱还做了相位/幅值核对：除以注入系数后，每频率 D Dᴴ 的相对差小于7×10⁻¹⁶；以固定 A_nominal 模型和实际 D/K 计算的纯激励图样效应约为机器精度。因此此处约5%–23%的变化不能归咎于四路随机相位换了一组。仍不能仅凭这些记录区分幅值非线性与随时间温漂。
""")
    text(table(['扰动系数','输出1偏离线性 %','输出2 %','输出3 %','输出4 %','原M3收敛'],[
        [r['gain']]+[fmt(v,2) for v in r['relative_linear_error_percent']]+[flag(r['M3_converged'])]
        for r in scale['amplitude_records']]))
    text("""0.6工况的原M3拟合失败，参数端点保留但不当作可信模态估计；其原始信号RMS仍可直接计算。即使只看已收敛的0.8–1.5，RMS偏离线性也可达约14.45%。D组只有nominal控制器，这些误差不能直接当作lowgain缩放误差的已测上界。

### 4.3 缩放后的条件偏差预测

按请求保持主0.5Hz噪声谱不变，把D及对应确定性瞬态幅值按α缩放；瞬态参数导数列仍保留。信号部分矩阵相应变换，噪声矩和曲率项继续计算，不能在有这些项时机械地把最终偏差除以α²。以下仍是原端点附近的局部线性化外推；较大偏差下近似误差可能不可忽略。
""")
    text(table(['控制器','扰动系数','M⁻¹b预测ζ2偏差 %','(M−C)⁻¹b预测ζ2偏差 %'],[
        [r['id'],fmt(r['disturb_gain'],4),fmt(r['pred_M_percent'][4]),fmt(r['pred_MminusC_percent'][4])]
        for r in scale['scalings']]))
    text(table(['控制器对','预测差绝对值（百分点）','原JK合成SE','差/SE','超过原SE的3倍'],[
        [r['a']+' vs '+r['b'],fmt(abs(r['prediction_difference'])),fmt(r['original_combined_jackknife_se']),
         fmt(r['difference_over_combined_se'],2),flag(r['exceeds_3_original_se'])] for r in scale['pairwise']]))
    text("""采用独立记录合成SE=sqrt(SE_a²+SE_b²)。highgain与另外两者相差超过当前合成SE的3倍；lowgain与nominal仅约1.26倍，**不能确认三个控制器两两都有充分差异**。

任务给定的“缩放后仍大于jackknife SE三倍”只能用现有记录SE做上述条件比较。降低激励后的SE尚未实测，噪声不变、系统不变的假设又已被D组提示为不充分，因此不能把该表写成补采后必然有3σ检出力；这项确认尚未成立。

## 5. 论文主图草稿

- `figures/main_damping_mode2.png/.pdf`：M1与M3第二阶阻尼，A/C++++用圆/方标记；误差棒为各方法按周期jackknife的1SE，单位是阻尼百分数的百分点。
- `figures/main_paired_prediction.png/.pdf`：实测相对配对差及其1SE，对比同记录0.5Hz谱预测；A/C++++分别用圆/方标记，预测为空心。预测点不附假造的置信区间，因为本轮未传播PSD估计、模型及K的不确定度。
- 另外6张逐记录谱诊断图、1张分频残差图、1张lowgain定位图、1张缩放图，共11张；均导出300dpi PNG和矢量PDF。

这些是记录内不确定度图，不是运行间复现性误差棒。不能用它们直接证明换控制器后绝对模态不变。

## 6. 预设判定、未解决事项和执行差异
""")
    text(table(['项目','预设要求','结果'],[
        ['任务1','nominal/highgain四记录预测/实测全部∈[0.67,1.5]','未通过；3/4，A_nominal不通过'],
        ['任务2第二阶','达到检测门槛后 |z|<2','未整体通过；5条可检，4条通过；A_lowgain失败；A_nominal低于预测检测限'],
        ['任务2第一/三阶','预测低于检测限时，实测≤2SE','未通过；部分第三阶仍显著超出'],
        ['任务3','M3/M3w Q/E[Q]≤2','未通过；0/12'],
        ['任务4','缩放后控制器间差异仍大于3SE','仅highgain与其余两者在当前SE下满足；lowgain/nominal不满足，未来SE未知']]))
    text("""执行差异与明确限制：

1. 主分析没有更改模型、频带、已有拟合、谱窗口或预设通过区间。0.5Hz是主结论，1Hz仅敏感性。
2. 增加了“固定传函，仅重估原24个瞬态系数”的裁切诊断，且包含drop=0；这是为分辨记录边界变化与噪声变化而作的附加分析，不用于替换任务1–3的主结果。
3. jackknife按请求拼接剩余周期，没有用连续四块替代；由此可能产生的内部边界是其适用限制，未通过增加瞬态模型去修饰结果。
4. 相同初值落实为各记录、各方法原保存的initial_beta在所有折中固定；保留原M1/M3不同初始化，而非强行将两方法改成同一向量。原停止规则和严格score分别报告。
5. 联合sandwich是局部Gaussian线性+二次近似，非新的一致性证明。Q的bootstrap仅是固定R和C的200次线性模拟，不是非线性重估bootstrap，也不构成精确χ²分布验证。
6. 缩放后的误差棒无法由现有不同激励的完整重复运行直接获得，未把当前jackknifeSE冒充未来SE；因此没有给出原任务所希望的“全部控制器差异仍有3SE”的确认。
7. 受保护文件哈希均一致。开发中修正过汇总脚本的NumPy布尔值打印，以及迹公式的浮点稳定性；后者改用QR核对，主Q/E计算和预设门槛未变。不存在删掉失败结果、改变权谱或调窗口凑判定。

论文可据此保留：控制器相关的方法配对差、部分记录上噪声谱修正后较好的定量对应、以及不能忽略的模型/运行条件误差。当前证据不足以声称“参考IV在这批实测上已达到噪声极限，且全部配对差由定理1解释”。下一步若补采，应控制响应幅值并做交错重复；本报告本身没有新增数据或开展下一轮模型优化。

## 7. 文件索引

`analysis/bias_predictions.csv`：18个阻尼参数×3种谱的两种偏差公式；`analysis/paired_tests.csv`：频率/阻尼配对差、JK和sandwich、门槛；`analysis/jackknife_status.csv`：失败、严格score和迭代；`analysis/residual_frequency_bands.csv`：全部2Hz分段；`analysis/amplitude_linearity.csv`：D组RMS。

`raw/*omit*.json`：82个重估端点、原始初值、score和求解状态；`raw/*prediction*.npz`：M、曲率、b和参数偏移；`raw/*paired_covariance.npz`：联合配对协方差；`raw/*residual.npz`：e、J、Z、C、W、分频Q和200个Q*；`spectra/*noise.npz`：完整谱、两窗口和裁切诊断。

`protected_hashes_before.json`与`analysis/verification.json`：只读保护审计；`code`、`baseline`、`logs`：复现代码、隔离副本和日志。运行步骤见README。
""")
    target=HERE/'低频偏差预测与残差检验报告.md'
    target.write_text('\n'.join(parts),encoding='utf-8')
    print('SAVED',target,'characters',len('\n'.join(parts)))
if __name__=='__main__':main()
