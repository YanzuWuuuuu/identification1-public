from analyze import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def rows(name):
    with (ROOT/'analysis'/name).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def tab(head,body):
    return '| '+' | '.join(head)+' |\n| '+' | '.join(['---']*len(head))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in body)
def fmt(x,n=3):
    try:return f'{float(x):.{n}f}'
    except:return str(x)
def plots():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    dd=rows('joint_delays_and_se.csv');mm=rows('joint_modal_changes.csv');qq=rows('residual_check.csv')
    fig,ax=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
    for w,color,off in [('H300','#2459a6',-.09),('H493','#188477',.09)]:
        a=[next(r for r in dd if (r['id'],r['window'],r['method'],r['parameter'])==(rid,w,'M3','common')) for rid in IDS]
        ax[0].errorbar(np.arange(6)+off,[float(r['estimate_us']) for r in a],
            yerr=[float(r['jackknife_se_us']) for r in a],fmt='o',capsize=3,label=w,color=color)
    ax[0].axhline(TAU.mean()*1e6,color='.5',ls='--',label='Previously fixed')
    ax[0].set(xticks=range(6),xticklabels=['A low','A nom','A high','C low','C nom','C high'],
              ylabel='Effective common delay (us)',title='Reference IV: estimate and 1 jackknife SE')
    ax[0].tick_params(axis='x',rotation=35);ax[0].legend(fontsize=8)
    for k,(w,j,label) in enumerate([('H300',1,'282'),('H300',2,'333'),('H300',3,'359'),('H493',0,'493')]):
        rr=[r for r in mm if r['window']==w and int(r['mode'])==j and r['method']=='M3']
        ax[1].scatter(k+np.linspace(-.13,.13,len(rr)),[float(r['change_over_repro_sd']) for r in rr],s=28,color='#2459a6')
    ax[1].axhline(1,color='#a42946',ls='--');ax[1].set(xticks=range(4),xticklabels=['282','333','359','493'],
        xlabel='Mode frequency (Hz)',ylabel='Damping change / reproducibility SD',title='Effect of jointly estimating delay')
    for a in ax:a.grid(alpha=.18)
    for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'delay_and_damping.{ext}',dpi=300)
    plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11,3.8),layout='constrained')
    for a,w in zip(ax,['H300','H493']):
        rr=[next(r for r in qq if (r['id'],r['window'])==(rid,w)) for rid in IDS];x=np.arange(6)
        a.semilogy(x,[float(r['baseline_qe']) for r in rr],'o-',color='.5',label='Fixed delay')
        a.semilogy(x,[float(r['joint_qe_fixed_noise']) for r in rr],'s-',color='#2459a6',label='Joint delay')
        a.axhline(2,color='#a42946',ls='--',lw=1);a.set(xticks=x,xticklabels=['A low','A nom','A high','C low','C nom','C high'],
            ylabel='Q / E[Q]',title=w+'; identical frozen noise covariance')
        a.tick_params(axis='x',rotation=35);a.grid(alpha=.18);a.legend()
    for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'residual_comparison.{ext}',dpi=300)
    plt.close(fig)
    sim=rows('simulation_individual.csv');fig,ax=plt.subplots(figsize=(7.2,3.7),layout='constrained')
    for method,col,off in [('M1','#a42946',-.08),('M3','#2459a6',.08)]:
        rr=[r for r in sim if r['method']==method and r['parameter']=='common']
        ok=[r for r in rr if r['success']=='True'];bad=[r for r in rr if r['success']!='True']
        ax.scatter([int(r['replicate'])+off for r in ok],[float(r['error_us']) for r in ok],s=20,color=col,label=method)
        if bad:ax.scatter([int(r['replicate'])+off for r in bad],[float(r['error_us']) for r in bad],s=80,marker='x',color=col,label=method+' failed score')
    ax.axhline(0,color='.5',lw=1);ax.set(xlabel='Monte Carlo replicate',ylabel='Common-delay error (us)',title='All 40 trials retained; final root solver endpoints')
    ax.legend();ax.grid(alpha=.18)
    for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'simulation_delay_error.{ext}',dpi=300)
    plt.close(fig)
def build():
    counts=rd(ROOT/'analysis/counts.json');dec=rd(ROOT/'analysis/decision.json')
    dd=rows('joint_delays_and_se.csv');mm=rows('joint_modal_changes.csv');qq=rows('residual_check.csv')
    fx=rows('fixed_sensitivity.csv');cross=rows('delay_across_records.csv');tests=rows('rule_c_tests.csv')
    pp=rows('paired_delay_difference.csv');ss=rows('simulation_summary.csv');jk=rd(ROOT/'analysis/jackknife.json')
    v=rd(ROOT/'analysis/validation.json');sv=rd(ROOT/'analysis/simulation_validation.json')
    text=['# 时延联合估计可行性检查\n',
    '**结论：按预先给定的判据，建议论文保留时延估计这一项，至少必须报告其对高频阻尼测量的影响。判据(a)和(b)明确通过。当前证据不支持宣称“所有窗口均能可靠联合辨识物理时延”，也不支持把固定时延的影响写成低于测量复现性。**\n',
    '高频的24个全记录端点（6记录×2窗口×2方法）及164个删除周期端点最终通过严格score；低频原IV–GN的12个全记录端点全部未通过。追加完整求根后找到9个低频代数根，但仍出现毫秒级、包含大量负值的有效时延，并对初值敏感，不能当作可信物理时延标定。40次已知时延仿真中，原IV–GN的参考IV失败39次；改为同一IV方程的完整求根后仍有2次未通过。上述失败均保留。\n',
    '## 1. 范围、数据与计算口径\n',
    '仅使用A和C++++的6条既有记录。低频3–25 Hz保留158个原参数；H300保留200个原参数（约263 Hz背景分量仍是干扰参数，不当成板模态）；H493保留74个原参数。各自仅增加7个结构化时延参数，分别为165、207、81维。C0/C1、D0/D2、ZOH与逐模态瞬态没有替换。M1用窗口内全部DFT点，M3仅用外部周期激励谱线，参考输入随当前模型更新。噪声谱为同记录非激励点0.5 Hz平均，窗口两侧各外扩2 Hz；全长未归一化DFT的协方差为NΦ。\n',
    '公共时延网格以原有效公共值为中心，范围±200 μs、步长10 μs。每点在原极点处剖面求线性系数，随后从原端点和两个相隔至少50 μs的低残差网格点启动完整联合IV。网格用于初始化，并不是把时延约束在该区间内。输出/输入零和差异以Helmert正交基实现；计算单位100 μs。由exp(−sτ)给出的解析时延导数已核查。两个零和向量共8个分量但只有6个独立自由度，连公共项共7个。\n',
    '## 2. 原有效时延的加性结构\n',
    '以下均已包含此前减去100 μs的修正，本次未再减一次。约束最小二乘的公共项为 **496.655417 μs**。\n',
    tab(['侧','1','2','3','4'],[['输出 a (μs)']+[fmt(x,6) for x in v['structure']['output_us']],
                              ['输入 b (μs)']+[fmt(x,6) for x in v['structure']['input_us']]]),
    '\n有效矩阵τ（μs）：\n',
    tab(['输出\\输入','1','2','3','4'],[[j+1]+[fmt(x,6) for x in row] for j,row in enumerate(v['structure']['matrix_us'])]),
    f"\n最大绝对结构残差{np.max(abs(np.asarray(v['structure']['residual_us']))):.3g} μs，零和约束满足舍入精度。它说明**当前使用矩阵本来就是加性结构**，并非由新实验独立证明硬件时延可精确拆成两端纯时延。\n",
    '## 3. 统一偏移敏感性\n',
    '六个偏移档全部计算，并另做零偏移复算，共252个拟合；初始20个未通过端点在不改方程的半步松弛后全部通过。原端点与精化端点均保存。下表对6条记录、正负六档取最大绝对变化；零偏移不参与最大值。阻尼变化单位为百分点，避免与相对百分数混淆。\n']
    tb=[]
    for w,j,label in [('L',0,'低频1'),('L',1,'低频2'),('L',2,'低频3'),('H300',1,'282'),('H300',2,'333'),('H300',3,'359'),('H493',0,'493')]:
        for method in ['M1','M3']:
            rr=[r for r in fx if r['window']==w and int(r['mode'])==j and r['method']==method and float(r['shift_us'])!=0]
            ref=REPRO.get((w,j),np.nan);maxz=max(abs(float(r['change_zeta_pp'])) for r in rr)
            tb.append([label,method,fmt(max(abs(float(r['change_f_hz'])) for r in rr),6),fmt(maxz,6),
                       fmt(ref,6) if np.isfinite(ref) else '未给定',fmt(maxz/ref,2) if np.isfinite(ref) else '—'])
    text+=[tab(['模态/Hz','方法','max |Δf|/Hz','max |Δζ|/百分点','复现性SD/百分点','变化/SD'],tb),
    '\n全部逐记录、逐偏移、逐模态结果及收敛状态见 [fixed_sensitivity.csv](analysis/fixed_sensitivity.csv)。没有用最大值替代原始结果。低频第一、三阶缺少本任务指定的复现性SD，因此不补造比值。\n',
    '## 4. 联合估计与jackknife\n',
    '每次删除一个完整周期后，拼接剩余周期并重新构造DFT、频率网格与参考输入；共246折。SE为sqrt((P−1)/P·Σ(θ_(i)−θ̄)²)。这是同记录条件下的删除周期离散，不含运行漂移、模型失配和跨局部解风险。所有折都参与计算，失败没有删除；对应全记录或任一折失败时，SE标为无效诊断量，不能作置信区间。\n',
    '### 4.1 高频公共时延\n']
    tb=[]
    for rid in IDS:
        for w in ['H300','H493']:
            rr=[next(r for r in dd if (r['id'],r['window'],r['method'],r['parameter'])==(rid,w,m,'common')) for m in ['M1','M3']]
            tb.append([rid,w]+[f"{fmt(r['estimate_us'])} ± {fmt(r['jackknife_se_us'])}" for r in rr])
    text+=[tab(['记录','窗口','M1 τc±SE/μs','M3 τc±SE/μs'],tb),
    '\n每条记录的τc、a1…a4、b1…b4及其SE完整列于 [joint_delays_and_se.csv](analysis/joint_delays_and_se.csv)，16通道矩阵保存在对应raw或refinements的JSON中。没有把两端绝对硬件延时当成可独立求出的量。\n',
    '### 4.2 六次运行间的时延离散\n',
    '下表为M3的六记录平均±样本SD，含控制器与运行变化，不应叫作纯重复性。\n']
    tb=[]
    for p in PARAMS:
        rr=[next(r for r in cross if (r['window'],r['method'],r['parameter'])==(w,'M3',p)) for w in ['H300','H493']]
        tb.append([p]+[f"{fmt(r['mean_us'])} ± {fmt(r['sd_across_six_us'])}" for r in rr])
    text+=[tab(['时延分量','H300/μs','H493/μs'],tb),
    '\n两个窗口给出的通道相对时延并不完全相同，尤其输入侧分量。拟合所得是**给定窗口与背景模型下的有效时延**，可能吸收遗漏动态或残余相位误差；尚不能把它们直接当成跨频带不变的硬件延迟标定值。H300公共项平均300.915 μs，比原固定值小195.741 μs，接近5 kHz下的一个采样周期200 μs；这提示应核对旧时延的复用条件和信号记录时标，但本检查没有证据将差异直接归因为一拍记录错位或硬件改变。\n',
    '### 4.3 高频模态参数变化\n',
    '下表为M3联合时延相对原固定时延M3的变化。完整M1/M3频率、阻尼、SE与变化量见 [joint_modal_changes.csv](analysis/joint_modal_changes.csv)。\n']
    tb=[]
    for rid in IDS:
        for w,j in [('H300',1),('H300',2),('H300',3),('H493',0)]:
            r=next(r for r in mm if (r['id'],r['window'],r['method'],int(r['mode']))==(rid,w,'M3',j))
            tb.append([rid,fmt(r['frequency_hz'],3),fmt(r['damping_percent'],5),fmt(r['change_f_hz'],5),
                       fmt(r['change_zeta_pp'],6),fmt(r['change_over_repro_sd'],2)])
    text+=[tab(['记录','联合f/Hz','联合ζ/%','Δf/Hz','Δζ/百分点','|Δζ|/复现SD'],tb),
    '\n282 Hz变化为复现性SD的4.59–6.51倍，333 Hz为5.43–9.83倍；这两个模态的6条记录均超过门槛。359和493 Hz均未超过各自门槛。因此判据(a)通过，且不是仅一个异常记录触发。\n',
    '### 4.4 低频失败及初值敏感性\n']
    tb=[]
    for rid in IDS:
        for method in ['M1','M3']:
            r=load(ROOT/'raw'/f'joint_{rid}_L_{method}_full.json');kk=jk['|'.join((rid,'L',method))]
            tb.append([rid,method,fmt(r['delay']['common_us'],1),f"{r['solver']['score']:.3g}",
                       f"{kk['successful']}/{kk['P']}",'未通过'])
    text+=[tab(['记录','方法','失败端点τc/μs','严格score','通过的删除周期折','全记录'],tb),
    '\n低频原IV–GN时延端点漂移到毫秒尺度，部分有效通道时延变负；这12个端点均未达到score<1e−6。低频的Q/E[Q]有时表面降到1附近，但E[Q]因病态残差算子增大到约4×10⁵–1.8×10⁶，不能把这种比值下降解释成拟合已达到噪声水平。本任务不把这些端点计入(a)(b)(c)。\n',
    '高频初始24个全记录中有1个及164折中有2个未达门槛，均为C_lowgain、H493、M3。原迭代呈小幅交替振荡；追加半步松弛后通过。多个已收敛初值的最大时延矩阵差异见 [multistart_spread.csv](analysis/multistart_spread.csv)。这不构成全局唯一性证明。\n']
    diagnostics=[rd(p) for p in sorted((ROOT/'diagnostics').glob('low_fullroot*.json'))]
    dr=[]
    for r in diagnostics:
        for a in r['attempts']:
            d=a.get('delay',{});mat=np.asarray(d.get('matrix_us',[np.nan]));si=a.get('solver',{})
            dr.append([r['id'],r['method'],'原失败端点' if a['initial']=='previous_endpoint' else '原固定时延',
                       '通过' if si.get('success') else '未通过',fmt(si.get('score',np.nan),7),
                       fmt(d.get('common_us',np.nan),1),f"{fmt(mat.min(),1)}…{fmt(mat.max(),1)}"])
    text+=['\n追加低频完整矩函数求根诊断（每个初值最多2000次函数评估）：\n',
    tab(['记录','方法','初值','score门槛','score','τc/μs','16通道有效时延范围/μs'],dr),
    '\n从原失败端点重启得到9/12个代数根，从原固定时延重启则0/12通过。不能据原GN失败宣称低频没有根，但找到的这些根不能给出与高频一致、物理可信的纯延时。零和规范只移除两端分解的自由度，不能防止时延吸收模型失配。追加求根作为失败诊断保存，不替换主协议端点或其jackknife；这些根未另做删除周期验证，因此不提供伪装成可信置信区间的SE。低频联合时延的结论是未获可靠标定，而非已证明结构不可辨识。\n',
    '## 5. 同记录Q/E[Q]与判据(b)\n',
    '主比较冻结各记录原M3端点得到的同记录噪声谱；新旧模型各用自己的完整雅可比和参考工具计算残差算子R，E[Q]=tr(C⁻¹RCRᵀ)。额外列出联合估计后重新估谱的敏感性版本。以下只比较M3；M1端点文件中的ref-operator诊断不作为M1自身统计拟合优度使用。\n']
    tb=[]
    for rid in IDS:
        for w in ['H300','H493']:
            q=next(r for r in qq if (r['id'],r['window'])==(rid,w))
            tb.append([rid,w,fmt(q['baseline_qe'],1),fmt(q['joint_qe_fixed_noise'],1),
                fmt(q['reduction_percent'],1)+'%',fmt(q['joint_qe_reestimated_noise'],1)])
    text+=[tab(['记录','窗口','固定τ','联合τ/固定噪声','下降','联合τ/重估噪声'],tb),
    f"\n12组下降百分比的算术平均为 **{dec['b_mean_reduction_percent']:.2f}%**；先对Q/E[Q]取平均再算下降为{dec['b_ratio_of_means_reduction_percent']:.2f}%。重估噪声版本的平均下降为{dec['b_reestimated_mean_reduction_percent']:.2f}%。两种口径都远超过30%，判据(b)通过。联合后H300仍约307–708，H493约55–148，均远大于2，说明尚有明显模型失配或非线性/时变成分，时延不能解释全部误差。\n",
    '## 6. M1−M3时延差与控制器\n',
    '没有硬件时延真值，实测上只能把M1−M3称为方法间配对差，不能把它直接等同于M1的真实偏差。配对SE是在同一删除周期折内先相减再计算。完整9分量结果与所有逐分量检查均已导出，没有只保存显著结果。\n']
    tb=[]
    for t in tests:
        if t['pass_rule_c']!='True':continue
        arr=[]
        for ct in CTRS:
            rid=t['group']+'_'+ct+('_pppp' if t['group']=='C' else '')
            r=next(r for r in pp if (r['id'],r['window'],r['parameter'])==(rid,t['window'],t['parameter']))
            arr.append(f"{fmt(r['difference_us'],4)} ± {fmt(r['paired_se_us'],4)}")
        tb.append([t['window'],t['group'],t['parameter']]+arr)
    text+=[tab(['窗口','组','分量','lowgain/μs','nominal/μs','highgain/μs'],tb),
    '\n上列分量满足“样本均值随增益单调，三个增益下差值均超过2SE”的逐分量读法。但同一分量的该趋势未在A和C两组同时复现；也检查了多个相关分量，未作多重比较校正。b4的前两个均值尤其接近，不能把样本次序当成趋势的显著性检验。故(c)只能列为探索性、部分支持，不作为论文纳入时延的主要依据；决策已由(a)(b)充分触发。\n',
    '## 7. 已知结构化时延仿真：40次MC\n',
    '沿用cl2的3支承模态、精确ZOH、5 kHz采样、FC=60控制器、P0=30000、10周期、2周期settle、sd=0.2、sv=3e−7、输入扰动SD=0.03、AL=30/AH=100，频带3–100 Hz，种子1002…1041。仅在执行器和传感器链加入整数采样时延状态；原126个拟合参数加7个时延，不额外加入D0/D2。输入时延为[0,200,0,200] μs，输出为[200,0,200,0] μs，公共项200 μs，输出零和差异[100,−100,100,−100] μs、输入零和差异[−100,100,−100,100] μs。\n',
    f"预检时，较大的试设配置[0,200,400,200]与[400,200,0,200] μs使原闭环谱半径为1.001087而失稳，未用于MC；日志保留。最终稳定配置谱半径{sv['radius']:.9f}。无噪声输出与延迟ZOH模型的相对误差{sv['noiseless_relative_residual']:.3g}，与cl2解析频响对照误差{sv['response_vs_cl2']:.3g}；无噪声近真值恢复的最大时延误差{np.max(abs(np.asarray(sv['noiseless_delay_error_us']))):.3g} μs。无噪声归一化score在残差达舍入级时不适合作为停止判据，验证以输出及参数误差为准。\n",
    '原IV–GN：M1失败1/40，M3失败39/40。半步松弛的诊断没有解决参考IV的根本问题，因此在不改任何矩条件、数据或模型的情况下，追加完整IV矩函数的hybrid求根（包括工具和重建输入随参数变化的影响）。这与仅用ZᴴJ近似雅可比不同；解析时延输出导数仍使用−s exp(−sτ)。最终M1通过40/40，M3通过38/40。\n',
    '下表首先保留全部40个最终端点，失败端点也计入均值、SD和RMSE；另列通过样本的统计，二者不混用。单位μs。\n']
    tb=[]
    for r in ss:
        tb.append([r['method'],r['parameter'],r['failed']+'/40',
            fmt(r['mean_all_us'],2),fmt(r['sd_all_us'],2),fmt(r['rmse_all_us'],2),
            fmt(r['mean_converged_us'],2),fmt(r['sd_converged_us'],2)])
    text+=[tab(['方法','分量','最终失败','全体偏差','全体SD','全体RMSE','通过者偏差','通过者SD'],tb),
    '\n公共时延：M1全体偏差+89.00 μs、SD80.38 μs；M3全体偏差+8.53 μs、SD17.89 μs。参考IV明显减小本配置下的时延误差，但其38个通过样本仍有+7.97 μs平均偏差（SD17.66 μs）；不能称本组MC已经证明无偏。它只是单一N、单一控制器和有限SNR的结果。有限样本非线性、弱时延方向、根选择及离散延时寄存器的有限记录边界项可能共同影响，当前试验未分离这些原因。\n',
    '仿真使用真实离散延时寄存器，原逐模态瞬态类保持不变；有噪有限记录的延时寄存器边界项未另添FIR参数。此限制与原模型/范围约束一并保留，不能把有限样本残余偏差都归结为理论错误或都归结为数值失败。\n',
    '## 8. 预定判据与论文范围\n',
    tab(['判据','结果','可支持的结论'],[
        ['(a)任一高频阻尼变化超过复现SD','通过：282、333 Hz各6/6记录','固定时延会显著改变部分高频阻尼测量'],
        ['(b)高频Q/E[Q]平均下降>30%','通过：78.41%','联合有效时延能减少系统相位/拟合失配，但未达到噪声水平'],
        ['(c)时延配对差随增益单调且>2SE','逐分量部分命中；跨A/C未复现同分量','探索性支持，不主张已证实物理时延偏差规律']]),
    '\n**建议加入论文，但定位要收窄：把时延作为影响高频模态测量的待估干扰参数及误差来源，不把“低频闭环SRIVC一致性”和“所有窗口都能精确标定硬件时延”合并成一个已经完成的结论。** 原闭环偏差主线仍由支承模态检验支撑；高频时延结果可作为模型充分性和阻尼敏感性的扩展。当前不能写“时延固定的影响低于测量复现性”。\n',
    '## 9. 与任务细节的差异、失败和完整性\n',
    '1. 公共时延±200 μs网格解释为相对当前公共有效时延的偏移，不是绝对时延限制；每个网格点用于固定极点的线性剖面初始化，未在41个网格点都做完整非线性优化。联合指各窗口内部的模态参数和7个时延同时求解；未另做三个窗口共享同一套时延的跨窗口联合优化。\n'
    '2. 原迭代之外增加半步松弛和仿真完整矩函数求根；原失败文件不覆盖，精化文件另存refinements。仿真原GN失败39/40不能被最终38/40成功掩盖。\n'
    '3. 低频联合估计失败，SE及Q/E[Q]只作失败诊断，未因得到较低比值而宣告改善。原低频没有显式极点约束，本次加入宽于实际估计位置的模态窗口数值保护；记录未贴该保护边界。独立时延坐标±10 ms只作数值保护，不是物理先验。\n'
    '4. Q/E[Q]的主比较使用同一冻结噪声谱，另给重估谱敏感性；M1文件中的参考残差算子结果未拿来充当M1估计器的严格统计检验。\n'
    '5. 仿真新增时延先做闭环稳定性预检，失稳试设配置保留但不用于MC；实际40次使用稳定配置。M3按实测流程用激励谱线；未新增参考IV全频点方法或改变实测窗口。\n'
    '6. 判据(c)未指定公共项还是哪一个相对分量，也未预设多重检验。本报告列全分量和所有命中，不用挑出单一显著分量来替代跨运行验证。判据(a)(b)独立成立，结论不依赖这一歧义。\n'
    '7. 已有文件保护与哈希复核见analysis/audit.json；新增图为PNG 300dpi与矢量PDF，代码、原端点、所有折与精化端点全部保留。\n',
    '## 10. 图和逐项数据\n',
    '![时延与阻尼](figures/delay_and_damping.png)\n\n![残差对比](figures/residual_comparison.png)\n\n![40次MC](figures/simulation_delay_error.png)\n',
    '- [固定偏移全部结果](analysis/fixed_sensitivity.csv)\n'
    '- [联合时延及删除周期SE](analysis/joint_delays_and_se.csv)\n'
    '- [全部模态参数变化](analysis/joint_modal_changes.csv)\n'
    '- [时延配对差及联合SE](analysis/paired_delay_difference.csv)\n'
    '- [六记录间时延离散](analysis/delay_across_records.csv)\n'
    '- [Q/E[Q]完整结果含失败](analysis/residual_check.csv)\n'
    '- [MC逐次结果含失败](analysis/simulation_individual.csv)\n'
    '- [MC汇总](analysis/simulation_summary.csv)\n'
    '- [逐分量判据(c)](analysis/rule_c_tests.csv)\n'
    '- [运行方法](README.md)\n']
    (ROOT/'时延可行性报告.md').write_text('\n'.join(text),encoding='utf8')
if __name__=='__main__':plots();build()
