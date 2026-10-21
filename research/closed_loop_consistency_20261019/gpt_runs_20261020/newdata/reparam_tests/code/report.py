from profile_core import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

def fmt(x,n=3):
    if x is None:return '—'
    try:
        if not np.isfinite(float(x)):return '—'
        return f'{float(x):.{n}f}'
    except:return str(x)
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|',
                      *['| '+' | '.join(str(v) for v in row)+' |' for row in rows]])
def main():
    p6=rd(ROOT/'analysis/p6_summary.json') if (ROOT/'analysis/p6_summary.json').exists() else {}
    p7=rd(ROOT/'analysis/p7_sim_summary.json') if (ROOT/'analysis/p7_sim_summary.json').exists() else []
    endpoints=p6.get('endpoints',[]);gates=p6.get('residual_gates',[])
    p6raw=[rd(p) for p in (ROOT/'raw').glob('*_M[13]_full.json') if not p.name.startswith('p7')]
    text=['# 留数模型重参数化与控制器无关去偏检验报告',
          f'\n生成时间：{datetime.now():%Y-%m-%d %H:%M:%S}。原始数据、原结果、理论文档只读。判据采用理论 v1.1 第 8 节 P6、P7；未修改阈值。',
          '\n## 一、P6 模型选择结论',
          f'\n四模型 × 六记录 × 三窗口 × 两方法的完整数据拟合：已保存 {len(p6raw)}/144 个端点，其中严格通过 {sum(r.get("success",False) for r in p6raw)} 个。其余失败保留；下表数值不自动等同于有效结论。']
    g493=[g for g in gates if g['window']=='H493' and g['model'] in ('c','d') and g['valid']]
    fail493=sum(not g['increase_not_over_20pct'] for g in g493)
    if fail493:
        text+=['\n**当前证据要求保留 \(C_1\)**：在已有效求解的 493 Hz 对照中，'
               f'{fail493}/{len(g493)} 个 (c)/(d) 结果超过原定 20% 残差增长上限。不能以时延更稳定或参数更少为理由通过 P6。',
               '\n在一个模型覆盖全部窗口的前提下，自由 \(C_0,C_1\) 是保守选择。秩 1 共享振型 (b) 是否可在单个窗口使用，'
               '需结合下表判断；H300 包含原有宽背景项，该项也被施加秩 1 约束。因此，H300 拟合恶化不能单独归因为某一个真实模态不满足秩 1。']
    else:text+=['\n模型选择尚待完整有效端点，不提前宣布 P6 成立。']
    cross=p6.get('cross_window',[])
    crosssummary=[]
    for kind in 'abcd':
        xx=[r for r in cross if r['model']==kind and r['success']]
        if xx:crosssummary.append([kind,len(xx),fmt(np.mean([r['common_difference_us'] for r in xx])),
                                   fmt(np.mean([r['absolute_difference_us'] for r in xx]))])
    if crosssummary:
        text+=['\n'+table(['模型','有效配对数','平均 H493−H300 μs','平均绝对跨窗差 μs'],crosssummary),
               '\n(c)/(d) 的公共时延跨窗口差并未优于 (a)，因此 P6 的时延一致性预测本身也不获支持，'
               '并非仅因残差门槛未过。运行间离散同样并非所有模态都减小。'
               '模型 (b) 在单独的 H493 窗口保留了拟合质量，可作为局部模型的候选；目前证据不支持把它统一替换到三个窗口。']
    text+=['\n### 1.1 参数化与变量投影验证',
       '\n保持精确 ZOH、原模态数、原背景 \(D_0+D_2\)、原瞬态 \((t_0+t_1z^{-1})/a(z)\) 和原窗口。'
       '模型 (b)/(d) 采用单位输出振型，输入向量吸收幅值；(b) 进一步令 \((\\kappa_0,\\kappa_1)=(\\cos\\gamma,\\sin\\gamma)\)，消除剩余尺度自由度。',
       '\n参考输入依赖当前传函，所以不能把一次冻结工具的线性解误称为完整参考 IV 的闭式解。'
       '本实现固定非线性坐标后，反复更新参考输入并闭式解线性 IV 方程，直到线性矩为零；外层求剩余矩的根。'
       '秩 1 振型方向也是非线性坐标，外层并非只有极点和时延。若内层已收敛但外层矩未归零，则对剖面矩用 hybrid 求根；'
       '内层未收敛的情况保留为失败，不把未满足线性矩的结果当作有效剖面根。',
       '\n线性块满足 \(Z_l^H(Y-Hl)=0\)，冻结工具时 \(l=[\\Re(Z_l^HH)]^{-1}\\Re(Z_l^HY)\)。'
       '剖面外层的参考工具、极点导数和形状导数随当前参数重新计算。',
       '\n'+table(['窗口','(a) 自由 C0/C1','(b) 共享秩 1','(c) 自由 C0','(d) 秩 1 C0'],
                  [['3–25 Hz',165,93,117,90],['220–420 Hz',207,111,143,107],['450–550 Hz',81,57,65,56]])]
    vr=[]
    for p in (ROOT/'analysis').glob('verify_*.json'):
        r=rd(p);vr.append([r['window'],r['joint'],f"{r['max_absolute_difference']:.6e}",r['passed']])
    text+=['\n'+table(['验证窗口','联合时延','原根与剖面根最大绝对差','<1e−8'],vr),
           '\n秩 1 坐标方向导数的核查见 `analysis/coordinate_derivative.json`。仿真采用同一精确 ZOH 解析导数；'
           '原实验 A 使用前向差分导数，因此同时保留原 A 端点状态，数值细微变化不解释为统计改进。']
    text+=['\nM1/M2 在固定非线性坐标时工具等于输出梯度，故线性块可一次解正规方程；'
           '新增的直接线性剖面与原实现两个基准根的差分别为 9.05e−11、3.77e−13。'
           '这种参数根一致性以相应线性块可辨识为前提；若线性块秩亏，不同最小范数约定可能选择不同的零空间参数，'
           '不能扩展为“任何秩亏情况下参数根都唯一”。']
    text+=['\n### 1.2 同记录残差与预设阈值',
       '\n所有模型使用同一自由模型 M3 端点估计的 0.5 Hz 非激励点噪声协方差，避免换模型时改变残差标尺。'
       '对未收敛的低频自由模型，沿用原固定时延自由模型的协方差，并明确标记；失败端点的 Q/E 不用于模型通过判定。',
       '\nH300 的 333 Hz 和 359 Hz 局部评价段预先固定为 307.5–346 Hz、346–420 Hz，拟合仍使用整个 220–420 Hz。'
       '期望残差采用 IV 的非正交残差算子，不用简单的“样本数减参数数”。',
       '\n**弱辨识下须同时看 E[Q]/n。** 非正交 IV 可能极大放大噪声，使分母 E[Q] 膨胀。'
       '因此低频自由时延模型即便 Q/E[Q] 接近 1，也不能单凭该比值认定拟合达到了噪声水平或时延已可靠辨识。'
       '此处保留原 Q/E 判据，同时列原始 Q、E[Q] 和 E[Q]/n 揭示这种限制。',
       '\n'+table(['记录','模型','窗口','Q/E','333 Hz 局部 Q/E','359 Hz 局部 Q/E','收敛'],
         [[r['id'],r['model'],r['window'],fmt(r['QE']),fmt(r['QE333']),fmt(r['QE359']),r['success']]
          for r in endpoints if r['common_fixed_us'] is None]),
       '\n'+table(['记录','模型','窗口','Q','E[Q]','实数据维 n','E[Q]/n'],
          [[r['id'],r['model'],r['window'],fmt(r.get('Q'),1),fmt(r.get('EQ'),1),r.get('real_data_dimension'),fmt(r.get('expected_noise_amplification'))]
           for r in endpoints if r['common_fixed_us'] is None]),
       '\n'+table(['记录','模型','评价段','Q/E 相对 (a) 倍数','增长≤20%','端点有效'],
         [[g['id'],g['model'],g['window']+'/'+g['region'],fmt(g['QE_ratio_to_a']),g['increase_not_over_20pct'],g['valid']] for g in gates])]
    text+=['\n### 1.3 时延、跨窗口一致性及 jackknife',
       '\n完整的公共、四路输入和四路输出相对时延及 SE 见 `analysis/p6_delays.csv`。'
       '逐分量的 H493−H300 差值及同折配对 jackknife SE 见 `analysis/p6_cross_window_components.csv`；'
       '该表保留跨窗口的统计相关性，没有用两个 SE 简单相加替代。'
       '逐记录、逐模态的频率、阻尼及其 SE 见 `analysis/p6_modal_parameters.csv`，包含明确标记的宽背景分量。'
       '误差棒仅在完整删一周期折均有效时用于判据；少数成功折不用于替代失败折。',
       '\n'+table(['记录','窗口','模型','公共时延 μs','jackknife SE μs','已完成折','失败折','SE 有效'],
          [[r['id'],r['window'],r['model'],fmt(r['common_us']),fmt(r['common_se_us']),r['complete_folds'],r['failed_folds'],r['jackknife_valid']]
           for r in endpoints if r['common_fixed_us'] is None]),
       '\n'+table(['记录','模型','H493−H300 公共时延 μs','相对时延差 RMS μs','两端点有效'],
          [[r['id'],r['model'],fmt(r['common_difference_us']),fmt(r['relative_difference_rms_us']),r['success']]
           for r in p6.get('cross_window',[])])]
    text+=['\n### 1.4 运行间离散',
       '\n六条记录的 SD 同时含运行、控制器和状态差异；A/C 配对复现性单列在 CSV，不视为纯测量噪声。',
       '\n'+table(['窗口','模型','模态序号','平均 f Hz','f 的 SD Hz','平均阻尼 %','阻尼 SD 百分点','六端点有效'],
          [[r['window'],r['model'],r['mode_index']+1,fmt(r['frequency_mean_hz']),fmt(r['frequency_sd_hz'],6),fmt(r['damping_mean_percent'],6),fmt(r['damping_sd_pp'],6),r['all_success']]
           for r in p6.get('dispersion',[])])]
    direct=rd(ROOT/'analysis/p6_M1_summary.json') if (ROOT/'analysis/p6_M1_summary.json').exists() else []
    text+=['\n### 1.5 直接法 M1 的独立对照',
       '\nM1 使用全部频点，其端点和删周期时延 SE 单独报告；不把谱线拟合的 M3 残差算子用于 M1。'
       '全部九个规范化时延分量在 `analysis/p6_M1_delays.csv`，跨窗口分量差与配对 SE 在 '
       '`analysis/p6_M1_cross_window_components.csv`，六运行 SD 与 A/C 配对复现性在 '
       '`analysis/p6_M1_dispersion.csv`，模态及公共时延如下。',
       '\n'+table(['记录','窗口','模型','模态','f Hz','阻尼 %','公共时延 μs','SE μs','收敛','折数','SE有效'],
          [[r['id'],r['window'],r['model'],r['mode_index']+1,fmt(r['frequency_hz'],6),fmt(r['damping_percent'],6),
            fmt(r['common_us']),fmt(r['common_se_us']),r['success'],r['complete_folds'],r['jackknife_valid']] for r in direct])]
    text+=['\n## 二、低频核心检验重做',
       '\n自由估计与固定公共时延的结果分开。对低频公共时延漂移的 (a)/(b)，锚定值取同一记录两个高频窗口 M3 公共时延的等权平均，'
       '另外计算 ±25 μs 敏感性。此处仍估计六个相对时延，未为了收敛偷偷冻结它们。'
       '固定公共时延时的 jackknife SE 是条件 SE，不包括高频锚定值本身的不确定度。',
       '\n若 P6 不支持减少留数参数，则不能宣称“参数减少使低频配对差 SE 变小”。所需比较照样列出，失败也计入。']
    lowrows=[];lowse=[];lowresults=[]
    for rid in old.IDS:
        pp=ROOT/'analysis'/f'low_core_{rid}_a.json'
        if not pp.exists():continue
        r=rd(pp);lowresults.append(r);main=r['rows'][0];joint=r.get('joint_se',[np.nan]*6);z=r.get('z',[np.nan]*6)
        lowrows.append([rid,fmt(main.get('M1_damping',[np.nan]*3)[1],6),fmt(main.get('M3_damping',[np.nan]*3)[1],6),
                       fmt(main.get('observed',[np.nan]*6)[4]),fmt(main.get('prediction',[np.nan]*6)[4]),
                       fmt(joint[4]),fmt(z[4]),r['all_success'],r['complete']])
    text+=['\n'+table(['记录','M1 ζ2 %','M3 ζ2 %','配对差 %','M⁻¹b 预测 %','联合 SE','z','全成功','折数完整'],lowrows)]
    original=SRC/'followup_lowfreq/closeout/analysis/prediction_uncertainty_table.json'
    if original.exists():
        rr=[r for r in rd(original) if r['formula']=='M' and r['parameter']=='zeta2']
        text+=['\n原自由留数、固定原时延的主结果（只读引用）：',
               '\n'+table(['记录','实测配对差 %','M⁻¹b 预测 %','配对差 SE','联合 SE','z'],
                 [[r['id'],fmt(r['observed']),fmt(r['prediction']),fmt(r['observed_se']),fmt(r['joint_difference_se']),fmt(r['z_joint'])] for r in rr])]
        lookup={r['id']:r for r in rr}
        for r in lowresults:
            oldrow=lookup[r['id']]
            lowse.append([r['id'],fmt(oldrow['observed_se']),fmt(r.get('observed_se',[np.nan]*6)[4]),
                          fmt(oldrow['joint_difference_se']),fmt(r.get('joint_se',[np.nan]*6)[4]),r['all_success'] and r['complete']])
        text+=['\n'+table(['记录','原配对差 SE','新配对差 SE','原联合 SE','新联合 SE','新SE完整有效'],lowse)]
    if len(lowresults)==6:
        good=[r for r in lowresults if r['all_success'] and r['complete']]
        text+=[f'\n固定公共时延的六条核心对照中，只有 {len(good)}/6 条的配对父端点与全部折均有效。'
               'A_nominal 与 A_highgain 的配对差 SE 分别由约 2.620、6.822 增至 2.699、7.833 个相对百分比单位，'
               '**未发现配对差不确定度下降**。A_highgain 的预测与实测相容；A_nominal 的联合 z 约 2.18，仍超过原来的 2 门槛。'
               '其他记录的无效折不用于宣布改善。',
               '\n因此本次重参数化没有使“直接法随控制器变化、参考 IV 更稳定”这一主张获得更干净的整体实测证据。'
               '应保留原自由留数固定时延主结果，本次联合相对时延分支作为未通过的敏感性分析，而非替换主结果。']
    sp=ROOT/'analysis/p6_fixed_low_sensitivity.json'
    if sp.exists():
        ss=rd(sp)
        text+=['\n### 2.1 固定公共时延 ±25 μs 的敏感性',
               '\n以下列第二阶阻尼相对锚点的变化（绝对百分点），所有模态的频率、阻尼和失败标志保存在 '
               '`analysis/p6_fixed_low_sensitivity.csv`。锚点或扰动端点任一失败时，变化仅作数值记录，不用来证明稳健。',
               '\n'+table(['记录','模型','方法','锚点 μs','偏移 μs','ζ2 %','ζ2变化 百分点','锚点成功','偏移成功'],
                   [[r['id'],r['model'],r['method'],fmt(r['anchor_us']),fmt(r['shift_us'],0),
                     fmt(r['damping_percent'],6),fmt(r['damping_change_pp'],6),r['anchor_success'],r['success']]
                    for r in ss if r['mode_index']==1])]
    text+=['\n## 三、P7 控制器无关 jackknife 去偏',
       '\n主统计量在原 log 频率、log 阻尼坐标上计算：\n'
       '\\[\\hat\\theta_{JK}=P\\hat\\theta_P-(P-1)P^{-1}\\sum_i\\hat\\theta_{(-i)}.\\]\n'
       '随后统一转为 \(100\\operatorname{expm1}(\\hat\\theta-\\theta_0)\)。'
       'P7 的 \(O(P^{-2})\) 是阶数结论，并不保证在固定有限 P 下偏差小于两倍 MC 标准误；后者按原判据检验。',
       '\n附带诊断在 CSV 中列出 log 偏差及指数回变换的差额；它用于解释有限样本/Jensen 项，不替代指定的 expm1 百分比判据。',
       '\n**不能据此宣称物理阻尼偏差已验证为 \(O(P^{-2})\)。** 本次修正在 log 坐标上进行；'
       '指数回变换一般仍产生约为 log 估计方差一半的 Jensen 偏差，它可保持 \(O(P^{-1})\)。'
       '9 个有限样本门槛通过与验证渐近阶数是两件事。',
       '\n仿真沿用实验 A 的记录、设置与种子。已有 A 端点仅作剖面求根的热启动；未改变数据和目标方程。'
       '新增 P=16 的 20 条使用同一 seed 公式。M2 仅用激励谱线，M3 与 A 一样保留全部频点。'
       '删一周期后的激励线 DFT 是剩余各周期 DFT 的和，频率为原周期网格；瞬态项仍含常数项。',
       '\n### 3.1 Monte Carlo',
       '\n下表统计全部有限端点；失败不丢弃，并单列失败数。只有 40 条完成且无必要端点失败时，'
       '“偏差 <2 MC SE”才作为完整通过证据。方差比使用相同配对样本相对于 M2 的样本方差。',
       '\n'+table(['P','方法','参数','偏差 % ± MC SE','SD %','RMSE %','方差/M2','完成/40','失败','<2SE','判定有效'],
          [[r['P'],r['method'],r['parameter'],fmt(r['bias_percent'])+' ± '+fmt(r['standard_error']),
            fmt(r['sd_percent']),fmt(r['rmse_percent']),fmt(r['variance_ratio_to_M2']),r['completed'],r['failures'],
            r['below_two_MC_SE'],r['valid_primary_gate']] for r in p7])]
    jkrows=[r for r in p7 if r['method']=='JK' and r['parameter'].startswith('zeta')]
    if len(jkrows)==9 and all(r['valid_primary_gate'] for r in jkrows):
        count=sum(r['below_two_MC_SE'] for r in jkrows)
        text+=[f'\nP7 仿真阻尼判据：{count}/9 个 P×模态组合满足。'+('全部满足原判据。' if count==9 else '未全部满足，不能宣布 P7 仿真预测整体通过。')]
    else:text+=['\nP7 Monte Carlo 判定尚不完整；以表中完成数和失败数为准。']
    realrows=[]
    real_damping_pass=real_all_pass=real_all_count=real_damping_count=0
    paired_damping_pass=paired_all_pass=paired_second_pass=0
    for pp in (ROOT/'analysis').glob('p7_real_*.json'):
        r=rd(pp)
        for j in range(6):
            realrows.append([r['id'],['f1','f2','f3','ζ1','ζ2','ζ3'][j],fmt(r['difference_vs_M3_percent'][j],6),
                            fmt(r['combined_se_percent'][j],6),fmt(r['paired_se_percent'][j],6),
                            bool(r['passes_one_combined_se'][j]),r['failures'],r['total_fits']])
            real_all_count+=1;real_all_pass+=int(r['passes_one_combined_se'][j] and r['failures']==0)
            pairedpass=abs(r['difference_vs_M3_percent'][j])<r['paired_se_percent'][j] and r['failures']==0
            paired_all_pass+=int(pairedpass)
            if j==4:paired_second_pass+=int(pairedpass)
            if j>=3:
                real_damping_count+=1;real_damping_pass+=int(r['passes_one_combined_se'][j] and r['failures']==0)
                paired_damping_pass+=int(pairedpass)
    text+=['\n### 3.2 实测',
       '\n为隔离 P7 的去偏作用，此项保持原低频自由留数、固定原时延的 158 参数模型，不混入 P6 的结构和时延改变。'
       '两方法均由剖面 IV 求解。修正估计量的 SE 用嵌套删一/删二周期重估，不把未修正 M2 的 SE 冒充修正量的 SE。'
       '理论文档未具体规定“合成 SE”的协方差处理，本次运行前将主判据操作化为两者 SE 的平方和开根。'
       '这一约定忽略同记录估计之间的协方差，并不等于差值的联合 SE；同折相减得到的配对 SE 另列，作为相关性诊断。',
       '\n'+table(['记录','参数','JK(M2)−M3 %','合成 SE','配对 SE','差<1合成SE','失败拟合数','总拟合数'],realrows)]
    if realrows:
        text+=[f'\n实测原判据：阻尼 {real_damping_pass}/{real_damping_count} 项通过，频率与阻尼合计 {real_all_pass}/{real_all_count} 项通过。'
               '第二阶六条记录均相容，但 A_nominal 和 A_highgain 的第三阶阻尼未通过，不能宣布实测 P7 全部成立。',
               '\n这里的“相容”只指规定的误差棒检验，实测没有真值。合成 SE 与同记录配对 SE 不同，前者通常更宽。'
               '而 1 SE 也不是跨多个参数的联合 95% 检验；这些统计含义不改变原门槛及逐项结果。',
               f'\n**相关性诊断：**如果仅为诊断而使用同记录差值的配对 SE，阻尼只有 {paired_damping_pass}/{real_damping_count} 项，'
               f'全部参数 {paired_all_pass}/{real_all_count} 项、第二阶 {paired_second_pass}/6 条满足差小于 1 SE。'
               '这不替换事先冻结的平方和门槛，但说明其通过率不能解释为两方法已统计等价，'
               '也不能据此宣布 P7 在实测中消除了全部闭环偏差。']
    text+=['\n## 四、边界、失败与协议说明',
       '\n1. 原定 P6 的 20% 上限、运行间离散判据，以及 P7 的 2 MC SE / 1 合成 jackknife SE 门槛均未修改。',
       '\n2. “固定极点与时延即可一次线性闭式求出参考 IV 全部参数”对当前自适应参考工具不成立。'
       '本任务采用内层工具更新的条件线性求根；秩 1 的方向参数放入外层。这是必要的数学实现说明。',
       '\n3. H300 的宽背景拟合分量保留且同样受模型层级约束。不能把它当作已确认的真实板模态。',
       '\n4. Q/E[Q] 显著大于 1 时，噪声模型未覆盖全部误差；更小的条件 SE 不能证明硬件时延更准确。',
       '\n5. 成功端点的判据为归一化 score <1e−6 且未贴参数边界，求解目标通常为 2e−8；'
       '所有端点、网格初始化和失败信息保存在 raw 中。旧实验 A 的原失败标志另外保存，不被新求解覆盖。',
       '\n6. 未完成或失败的折不插值，不以少数成功折计算有效主误差棒。'
       '低频高频锚定公共时延的结果只给条件不确定度；±25 μs 为明确标注的敏感性分析。',
       '\n7. 计算分批落盘。失败父端点不继续赋予 jackknife 误差棒，其未运行折单列为“父端点失败而未运行”，'
       '不能把这些折计为成功。数值失败、未运行折和待完成任务分别统计。']
    ap=ROOT/'analysis/audit.json'
    if ap.exists():
        au=rd(ap)
        text+=['\n### 4.1 最终文件与任务审计',
               f"\n原输入文件哈希核查：{au['original_files']} 个文件；发生变化 {len(au['changed'])} 个，丢失 {len(au['missing'])} 个。",
               '\n'+table(['方法','P6主分支应完成删周期折','已完成','数值失败折','父端点失败而未运行折'],
                   [[m,r['expected'],r['completed'],r['failed'],r['skipped_failed_parent']] for m,r in au['p6_jackknife_counts'].items()]),
               '\n'+table(['方法','低频固定公共时延应完成折','已完成','数值失败折'],
                   [[m,r['expected'],r['completed'],r['failed']] for m,r in au['fixed_low_jackknife_counts'].items()]),
               f"\n其余必要折缺失数：{len(au['p6_missing_required_folds'])}。P7 仿真每个 P 的记录数：{au['p7_sim_counts']}。"
               f"实测 P7 汇总 {au['p7_real_summaries']}/6 条。所有要求运行的拟合已有结果记录：{au['all_required_fits_recorded']}。",
               '\n全部失败端点索引及原因见 `analysis/audit.json` 的 `failed_endpoint_index`；'
               'P7 仿真失败同时在每次 MC 结果中保存。完成计算不表示所有模型或判据通过。']
    (ROOT/'重参数化与去偏检验报告.md').write_text('\n'.join(text)+'\n',encoding='utf8')
    if jkrows:
        fig,axs=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
        fig.suptitle('40 paired runs per P; one M3 failure at P=4 retained among finite endpoints',fontsize=9)
        colors={'M2':'#737373','M3':'#1764ab','JK':'#39834b'}
        for j,ax in enumerate(axs):
            for k,met in enumerate(('M2','M3','JK')):
                rr=sorted([r for r in p7 if r['method']==met and r['parameter']==f'zeta{j+1}'],key=lambda x:x['P'])
                if rr:ax.errorbar(np.array([r['P'] for r in rr])+(k-1)*.22,[r['bias_percent'] for r in rr],
                                  yerr=[r['standard_error'] for r in rr],label=met,color=colors[met],marker='o',capsize=3)
            ax.axhline(0,color='black',linewidth=.7);ax.set_xticks([4,8,16]);ax.set_xlabel('Periods P');ax.set_title(f'Support mode {j+1}')
            ax.grid(alpha=.2)
        axs[0].set_ylabel('Damping bias (%) ± 1 MC SE');axs[-1].legend()
        fig.savefig(ROOT/'figures/p7_damping_bias.png',dpi=300);fig.savefig(ROOT/'figures/p7_damping_bias.pdf');plt.close(fig)
    if len(g493)==12:
        fig,axs=plt.subplots(1,2,figsize=(9,3.6),layout='constrained')
        colors={'a':'#777777','b':'#1764ab','c':'#7c52a1','d':'#39834b'}
        for k,kind in enumerate('bcd'):
            gg=[g for g in gates if g['window']=='H493' and g['model']==kind]
            axs[0].plot(np.full(len(gg),k)+np.linspace(-.12,.12,len(gg)),[g['QE_ratio_to_a'] for g in gg],
                        'o',color=colors[kind],label=kind)
        threshold=axs[0].axhline(1.2,color='black',ls='--',lw=1,label='20% increase threshold')
        axs[0].legend(handles=[threshold],loc='upper left',fontsize=8,frameon=False)
        axs[0].set_xticks(range(3),['b: shared rank 1','c: C1 = 0','d: rank 1, C1 = 0'])
        axs[0].set_ylabel('493 Hz: Q/E[Q] relative to model a');axs[0].grid(axis='y',alpha=.2)
        for k,kind in enumerate('abcd'):
            xx=[r['absolute_difference_us'] for r in cross if r['model']==kind and r['success']]
            axs[1].plot(np.full(len(xx),k)+np.linspace(-.12,.12,len(xx)),xx,'o',color=colors[kind])
        axs[1].set_xticks(range(4),list('abcd'));axs[1].set_xlabel('Residue model')
        axs[1].set_ylabel('Common delay gap (µs)\n|H493 − H300|');axs[1].grid(axis='y',alpha=.2)
        fig.savefig(ROOT/'figures/p6_model_comparison.png',dpi=300);fig.savefig(ROOT/'figures/p6_model_comparison.pdf');plt.close(fig)
if __name__=='__main__':main()
