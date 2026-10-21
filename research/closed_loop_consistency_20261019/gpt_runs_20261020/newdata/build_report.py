from evaluate_results import *
import hashlib,platform

def table(h,rs):return '\n'.join(['| '+' | '.join(h)+' |','| '+' | '.join(['---']*len(h))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rs])
def fmt(v):return '/'.join(f'{x:.4g}' for x in v)
def pm(v,u):return '/'.join(f'{x:.4f} ± {y:.3g}' for x,y in zip(v,u))
def readjson(name):return json.loads((HERE/name).read_text(encoding='utf-8'))
def readcsv(name):return list(csv.DictReader((HERE/'analysis'/name).open(encoding='utf-8-sig')))

def interval_label(r):
    if r.get('duplicate_of'):return '同'+r['duplicate_of']
    if r['group']=='B':return f"记录 {r['common_start_s']:.3f}–{r['common_end_s']:.3f}"
    return f"注入 {r['injection_start_s']:.3f}–{r['injection_stop_s']:.3f}"

def main():
    inv=readjson('data_inventory.json');v=readjson('analysis/verification.json');core=readjson('analysis/core_decisions.json')
    no=readjson('analysis/noise_summary.json');fr=[readjson(f'analysis/{c}_orthogonal_frf.json') for c in CTRS]
    validation=readjson('analysis/frf_normalized_validation.json');sc=readjson('analysis/scans_summary.json')
    fits=[json.loads(p.read_text()) for p in sorted((HERE/'fits').glob('*.json'))];fails=[r for r in fits if r['failure']]
    out=['# 新实测数据分析报告：A–E五组闭环频域辨识','',
    '日期：2026-10-20。按提供的采集方案、任务清单和既有实测代码执行。输出限于本目录，原始27个MAT、旧算法与理论文档未修改。',
    '', '## 0. 结论与预设判据','',
    '**本批数据提供了有价值但不完整的支持：第二阶阻尼的直接法漂移清楚，参考IV削弱了这一变化；“三个模态全部在2倍不确定度内不变”没有通过。** '
    '不能把一个有利参数扩展为所有模态一致，不能把数值收敛扩展为模型正确。',
    '', 'A组三个控制器下的第二阶阻尼为：','',
    table(['控制器','M1全部频点直接法/%','M2谱线直接法/%','M3参考IV/%','M3w加权参考IV/%','M5n邻点SML/%'],
          [[c]+[f"{fit('A_'+c,m)['damping_percent'][1]:.4f}" for m in ('M1','M2','M3','M3w','M5n')] for c in CTRS]),
    '', '低增益→高增益时，M1第二阶阻尼变化+0.6167个百分点；M3变化−0.1098个百分点。'
    '用各控制器B噪声谱作局部偏差预测，预测M1偏差的变化为+0.3711个百分点，实测/预测为1.66。'
    '这支持反馈偏差机制的方向和量级，但不是全部参数的精确预测，也不能忽略M3自身的运行间变化。',
    '', table(['预设项目','结果','解释'],[
      ['M3三个模态f、ζ跨控制器均<2合成u','**未通过：8/18对满足**','第二阶阻尼3/3对满足，但其他参数未全部满足'],
      ['M3w同一要求','**未通过：9/18对满足**','还存在严格score未过的端点，冻结随机权u也不包含权估计影响'],
      ['M1的ζ1或ζ2跨控制器漂移>2u','满足局部模型近似下的数值条件','M1的u使用直接矩的Gaussian局部近似，不冒称定理2'],
      ['实测漂移与理论同向、同量级','第二阶阻尼部分支持；整体未通过','全部参数的方向和比值见表；原方案未给量级阈值，未事后新增'],
      ['C组不同控制器FRF相容','**未通过条件bootstrap检查**','nominal/highgain的比值反演仍不稳定，参照不能作真值'],
      ['D/重复实验不确定度覆盖','**不能认定通过**','D实际为幅值扫描；同幅A+C样本也有超出模型u的离散']]),
    '', f'完成{len(fits)}组主估计，{len(fails)}组按原迭代规则失败；全部端点与异常保留。'
    '共20幅图，每幅均有300dpi PNG和PDF；表格完整数值在analysis/*.csv。',
    '', '![跨控制器核心参数和预测漂移](figures/core_controller_zeta1_zeta2_f2.png)',
    '', '## 1. 数据整理、信号边界及质检','',
    'MAT内部压缩块存在已知zlib校验问题。按raw-deflate只读解压，检查流已结束、内部MAT可解析、16个变量完整；没有重写源MAT。'
    '按每路时间列对齐到5kHz公共采样格。时间打印精度误差约0.0005个采样，容许0.01个采样后取整数格，不对信号插值。',
    '', 'B组记录的disturb_input并不是实际注入：B中U=C到数值精度，实际参考为0。由B回归所得T是秩3控制输出子空间上的最小范数表示；'
    '由于C自身秩3，不能称4×4 T被唯一物理辨识出来。T在全部记录上用于核验U−TC=D的注入边界，且K=−TK_c在实际控制子空间内有效。',
    '', '注入区间以逐点‖U−TC−D‖<10⁻⁵max(‖D‖,10⁻⁴)且有效注入幅值非零检测。先丢1周期；'
    '若首个保留周期间去均值Y差大于后续差的中位数2倍，则继续丢首周期，至少留下3周期。这个规则在辨识前固定，完全不使用拟合结果。'
    'C_lowgain_ppmm共丢3周期，E的−0.5mm共丢2周期，其余激励记录丢1周期。全部原始相邻周期差与截取索引保存。',
    '', table(['记录','注入区间（B为记录区间）/s','丢弃周期','保留P','用途/备注'],[
        [r['id'],interval_label(r),
         r.get('discarded_initial_periods','—'),r.get('P','—'),('重复：'+r['duplicate_of']) if r.get('duplicate_of') else ('纯噪声，记录D未实际注入' if r['group']=='B' else ';'.join(r['warnings']) or '可用；不代表模型成立')]
        for r in inv]),
    '', '27个文件实际为26个独立记录，其中23个有激励、3个纯噪声。D_nominal_1p5与A_nominal字节完全相同，没有重复计数。'
    '没有看到电机极值长时间重复造成的明显削顶特征，但未提供硬件量程，不能完成绝对饱和/限位认证。传感器量化与物理记录单位保持原样。',
    '', '## 2. 控制器、噪声与实际激励','',
    '用B记录前半拟合完整ARX(8,8)的C←Y，后半与A/C/D/E独立数据验证。E的非零位置会增加常量设定值项：用训练半段恢复常量，独立后半验证；该常量不改变3–25Hz的传递函数。'
    '这消除了“位置非零就误报控制器改变”的问题，没有用近似低阶控制器替代。',
    '', table(['控制器','B后半ARX相对残差','独立A/B恢复K相对差','归一化C参考矩阵最大条件数'],[
        [c,f"{readjson('controllers/'+c+'.json')['validation_relative_rms']:.3g}",
         f"{r['K_A_vs_B_relative']:.3g}",f"{r['row_normalized_D_condition_max']:.6f}"] for c,r in zip(CTRS,v['independent_controller_and_excitation_checks'])]),
    '', '控制器恢复精度足以排除明显控制器符号/边界错误。所谓“精确”指记录算术精度范围内复现，不是代数系数零误差。'
    '控制器图用S_u=(I+KG)⁻¹、F=S_uK；1–500Hz的G来自C组非参数矩阵，未外推三支承模态模型。宽带回路/灵敏度曲线受该FRF不稳定影响，仅作诊断。',
    '', 'B噪声按每周期DFT形成V=Y−G_M3 U。以5邻频点平均、0.1对角收缩得到Φ_v，定义为N⁻¹/²归一化DFT的协方差；除以fs才是双边每Hz PSD。'
    '周期bootstrap区间只含给定G下的随机重复性。另保留用C非参数G提取的噪声作敏感性对照。',
    '', table(['控制器','B完整周期数','带内噪声能量指标','条件95%区间','A非激励残差谱/B谱'],[
        [r['controller'],r['P_noise'],f"{r['band_trace']:.4g}",fmt(r['band_trace_ci95']),
         f"{v['independent_controller_and_excitation_checks'][i]['A_nonexcited_residual_to_B_trace_ratio']:.3f}"] for i,r in enumerate(no['records'])]),
    '', table(['B噪声对比','能量比','条件95%比值区间','是否包含1'],[[r['a']+'/'+r['b'],f"{r['trace_ratio']:.3f}",fmt(r['ratio_ci95']),'是' if r['includes_one'] else '否'] for r in no['comparisons']]),
    '', '![B组输出等效噪声谱](figures/B_equivalent_noise_psd.png)',
    '', '**不能假定三次噪声工况等价。** lowgain的带内积分噪声较高；A非激励残差谱/B噪声谱分别约1.23、6.14、2.88倍。'
    '这可能包含运行噪声变化、未建模动态、非线性泄漏或瞬态残留，不能把这些残差全部解释为独立随机噪声。B代入u因此只是条件量，可能低估激励实验的不确定度。',
    '', table(['控制器','3–25Hz电机RMS（四路）','全部频点反馈噪声/观测U','谱线反馈噪声/观测U'],[
        [r['controller'],fmt(r['A_input_band_rms_per_channel'])]+[f"{q['feedback_noise_to_observed_input_percent']:.4g}%" for q in readjson(f"analysis/A_{r['controller']}_feedback_noise_fraction.json")]
        for r in v['independent_controller_and_excitation_checks']]),
    '', '三个控制器使用不同随机相位实现，未满足原计划的跨控制器完全相同激励；计算采用各自实际D，不假装设计相同。'
    '控制器改变也改变输入/位移幅值，可能改变弱非线性工作状态。C正交符号正确，但每条谱线上各通道幅值并不相同，原始D条件数可远大于1。',
    '', '## 3. 参数估计、标准不确定度与重复','',
    '模型完全沿用旧实测：3个支承模态+一般实系数C0/C1分子+D0+D2+固定通道时延减100µs+24个逐模态瞬态，精确ZOH，总158参数。'
    '初始频率10.6/11.7/16.1Hz，阻尼3%；IV最多80步，极点对数步长阈值10⁻⁵、截断0.05；M5n最多160函数评估。'
    'M1和M2使用实测U构造工具，M3/M3w使用精确恢复控制器和外部D构造工具。M3w和M5n共用16邻非激励频点、0.1收缩、一次更新。',
    '', 'M3的u由保留全部瞬态列的A⁻¹QA⁻ᵀ给出，Q含正频率实部的1/2因子及正确未归一化DFT尺度NΦ_v。'
    'M3w冻结实际估计权作同一代入计算，未包含权估计的一阶随机影响，不能冒称已证明其有限样本覆盖。'
    'M1/M2含输入噪声内生性，使用局部Gaussian线性+二次矩方差近似，并非直接套定理2；其标准差计算暂未包括曲率C校正。'
    'M5n采用99次成对周期重抽样并重新求解，冻结原联合噪声协方差；594次子拟合均收敛。',
    '', '表中u为一倍标准不确定度，阻尼以百分数表示，其u单位是百分点；频率单位Hz。数值成功与严格score<10⁻⁶分开。',
    '', table(['A控制器','方法','f1/f2/f3 ± u /Hz','ζ1/ζ2/ζ3 ± u /%','收敛/严格score'],[
        [c,method,pm(fit('A_'+c,method)['frequency_hz'],sd(unc('A_'+c,method))[:3]),pm(fit('A_'+c,method)['damping_percent'],sd(unc('A_'+c,method))[3:]),
         f"{fit('A_'+c,method)['converged']}/{fit('A_'+c,method)['strict_score_pass']}"] for c in CTRS for method in ('M1','M2','M3','M3w','M5n')]),
    '', '同一控制器A与C++++为两次独立运行，差值与合成u见 `A_C_independent_repeats.csv`。绝对两次差/√2也列出，但只有1个方差自由度，不能据此给精确覆盖率。'
    '全部C单次估计、D、E以及失败端点在 `all_parameter_estimates.csv`（564行）。',
    '', 'C组++++独立重复的完整结果如下，同样报告一倍条件标准不确定度：','',
    table(['C++++控制器','方法','f1/f2/f3 ± u /Hz','ζ1/ζ2/ζ3 ± u /%','收敛/严格score'],[
        [c,method,pm(fit('C_'+c+'_pppp',method)['frequency_hz'],sd(unc('C_'+c+'_pppp',method))[:3]),
         pm(fit('C_'+c+'_pppp',method)['damping_percent'],sd(unc('C_'+c+'_pppp',method))[3:]),
         f"{fit('C_'+c+'_pppp',method)['converged']}/{fit('C_'+c+'_pppp',method)['strict_score_pass']}"]
        for c in CTRS for method in ('M1','M2','M3','M3w','M5n')]),
    '', table(['失败记录','方法','终止','f/Hz','ζ/%'],[[r['id'],r['method'],str(r.get('iterations'))+'步',fmt(r['frequency_hz']),fmt(r['damping_percent'])] for r in fails]),
    '', f"此外，{sum(not r.get('strict_score_pass',False) for r in fits)}/94个端点未满足严格score门槛；报告未用很大的误差棒把这些端点算作理论验证通过。",
    '', '## 4. 直接法偏差与控制器漂移预测','',
    '对A和C++++各记录，以其M3端点作为局部G0、相应B谱作为噪声输入，使用实际参考谱和控制器计算M⁻¹b及(M−C)⁻¹b。'
    '完整158维矩阵包含瞬态列；未先删瞬态。C的外层中心差分步长2×10⁻⁴，避免嵌套过小步长。'
    '所有对数参数差均用100·expm1转换成相对百分比，而非把对数单位当百分比。',
    '', table(['记录','实测M1−M3阻尼相对差/%','M⁻¹b预测/%','(M−C)⁻¹b预测/%'],[
        [rid,fmt(r['observed_M1_minus_M3_percent'][3:]),fmt(r['pred_M_percent'][3:]),fmt(r['pred_MminusC_percent'][3:])]
        for rid in [f'{g}_{c}'+('_pppp' if g=='C' else '') for g in ('A','C') for c in CTRS]
        for r in [readjson(f'analysis/{rid}_bias.json')]]),
    '', '第二阶阻尼的跨控制器比较：','',
    table(['控制器变化','M1实测漂移/百分点','M3实测漂移/百分点','预测偏差漂移/百分点','M1漂移/预测','扣除M3漂移后/预测'],[
        [r['a']+'→'+r['b']]+[f"{float(r[k]):.4g}" for k in ('M1_observed_drift','M3_observed_drift','predicted_bias_drift','raw_drift_to_prediction','excess_drift_to_prediction')]
        for r in readcsv('bias_controller_drift.csv') if r['parameter']=='zeta2']),
    '', '**第二阶方向一致，量级接近但不精确；其他参数没有统一通过。** 例如低→高ζ1实测/预测约18.8，f1预测方向不一致，第三阶某些预测接近零而实测变化明显。'
    '不能只报ζ2的好结果。各控制器使用自身估计G0，属于局部代入诊断，不是独立已知真对象的预测；M3漂移明确单列，不把它归入直接法反馈偏差。',
    '', '## 5. C组正交非参数FRF与模型验证','',
    '每控制器四次实验按共同最小稳态周期数平均，逐频线拼D/U/Y矩阵，直接计算T_ud=UD⁻¹、T_yd=YD⁻¹、G=YU⁻¹。'
    '没有LRM、邻频平滑或逆矩阵正则化。每路D的行幅值归一化后条件数均为1到10⁻¹⁴精度，确认Hadamard方向正确。'
    '原始D的条件数较大来自各输入谱线幅值差，不能把正交等同于相同信噪比。',
    '', table(['控制器','共同P','max cond(D)','max cond(U)','前后段FRF差/%','随机SE/FRF范数/%','bootstrap均值偏移/%'],[
        [r['controller'],r['P_common']]+[f"{r[k]:.3g}" for k in ('D_condition_max','U_condition_max','G_halves_difference_percent','random_se_relative_percent','bootstrap_mean_shift_percent')] for r in fr]),
    '', 'C组解决了单线输入秩不足，但未自动解决闭环输入矩阵在共振处的近奇异。nominal/highgain出现明显有限比值曲率/不稳定，'
    'bootstrap均值偏移本身很大；其百分位区间不能视为无偏真实FRF置信带。lowgain参照相对稳定，但仍有运行间偏差。',
    '', '![低增益正交FRF与独立参数模型，8–18Hz](figures/C_lowgain_zoom_8_18.png)',
    '', '幅值/相位全通道图见figures/C_{lowgain,nominal,highgain}_{magnitude,phase}.png及同名PDF。'
    '相位各自独立unwrap，再对整条曲线作一个整数倍360°平移，使其与M3相位的中位差最小；跨控制器图以lowgain为共同显示分支。'
    '这只改变显示分支，不逐频跟随模型、不改变复FRF或误差统计。噪声导致的局部解缠歧义仍保留，不能把360°分支差当真实动力学误差。',
    '', table(['控制器对','复FRF差/%','点态95%带外比例','整体标准化统计/95%阈值','条件bootstrap相容'],[
        [r['controller_a']+'/'+r['controller_b'],f"{r['relative_difference_percent']:.2f}",f"{100*r['pointwise_outside_fraction']:.1f}%",
         f"{r['global_standardized_statistic']:.3g}/{r['centered_bootstrap_95_threshold']:.3g}",'是' if r['compatible_global_95'] else '否']
        for r in readjson('analysis/orthogonal_controller_comparison.json')]),
    '', '参数模型验证主比较为独立A模型对C四次实验，避免拿C++++拟合模型与包含自身的C参照做独立检验。'
    '统计量用Q=2Σ eᴴC_v⁻¹e，e=Ȳ−G_model Ū，C_v=P0Φ_B/P_common。若模型固定且准确、PSD已知、各频率proper Gaussian独立，则理想自由度为4256，期望4256，95%区间约[4077,4439]。'
    '实际PSD和模型均为估计，故这只是理想χ²参照。M1/M2/M3/M3w另加入A参数协方差并用499次条件Gaussian扰动校准，保留跨频率共同参数相关；仍未覆盖系统误差。',
    '', table(['控制器','方法','G对照差/%','求逆前残差/%','Q/理想期望','加入参数u后的Q/期望'],[
        [r['controller'],r['method'],f"{r['complex_G_difference_percent']:.2f}",f"{r['inverse_free_difference_percent']:.2f}",f"{r['statistic_over_expectation']:.2f}",
         f"{r['parameter_uncertainty_corrected_statistic']/4256:.2f}" if 'parameter_uncertainty_corrected_statistic' in r else '未提供全参数u校正'] for r in validation]),
    '', '没有任何方法达到“残差与纯随机噪声相当”的整体要求。加入参数u后仍有大失配，说明不能只用求逆放大解释；'
    '运行间对象变化、激励相关的未建模响应/非线性、B谱代表性不足都需要保留为解释，当前不能唯一归因。',
    '', '## 6. D幅值扫描、同幅重复性与E位置','',
    'D实际是0.6/0.8/1.0/1.2/1.5五档幅值，不是原计划的五次相同工况重复；1.5复用A nominal。'
    '因此D主要检验线性度/幅值相关变化，不能把其离散直接称为同分布重复覆盖率。nominal的A+四个C才提供五次同幅独立运行，激励方向不同。',
    '', table(['D增益','M3 f/Hz','M3 ζ/%','M3状态','M3w f/Hz','M3w ζ/%'],[
        [x,fmt(fit(rid,'M3')['frequency_hz']),fmt(fit(rid,'M3')['damping_percent']),fit(rid,'M3')['converged'],
         fmt(fit(rid,'M3w')['frequency_hz']),fmt(fit(rid,'M3w')['damping_percent'])]
        for rid,x in [('D_nominal_0p6',.6),('D_nominal_0p8',.8),('D_nominal_1p0',1),('D_nominal_1p2',1.2),('A_nominal',1.5)]]),
    '', '计算每次自己的u，中心为逆方差加权均值。因均值从同一五个样本估计，标准化残差RMS在理想条件下基准约√(4/5)=0.894，而不是硬性1。'
    '原方案的SD/u在0.7–1.4、约95%落入±2u仍按原样列出；五个样本无法精确验证95%覆盖。失败端点保留，相关组不作为通过证据。',
    '', table(['同幅A+C方法','参数','实测SD/RMS(u)','落入±2u/5','0.7–1.4内','数值条件均满足'],[
        [r['method'],r['parameter'],f"{r['observed_sd_to_rms_u']:.3g}",r['within_2u_count'],'是' if r['ratio_in_predeclared_0p7_to_1p4'] else '否',r['all_numerical_conditions']]
        for r in sc if r['kind']=='repeat']),
    '', table(['z/mm','M3 f/Hz','M3 ζ/%','M3收敛','M3w f/Hz','M3w ζ/%','M3w收敛'],[
        [x,fmt(fit(rid,'M3')['frequency_hz']),fmt(fit(rid,'M3')['damping_percent']),fit(rid,'M3')['converged'],
         fmt(fit(rid,'M3w')['frequency_hz']),fmt(fit(rid,'M3w')['damping_percent']),fit(rid,'M3w')['converged']]
        for rid,x in [('E_nominal_m1p0',-1),('E_nominal_m0p5',-.5),('A_nominal',0),('E_nominal_0p5',.5),('E_nominal_1p0',1)]]),
    '', '位置改变后的噪声只借用了z=0对应控制器B谱，E误差棒属于条件描述，不能用于位置专属覆盖认证。'
    '新z=−1记录本身周期和边界检查可用，但M1/M3/M3w求解未收敛，不能沿用旧数据“坏记录”标签直接删除，也不能将失败端点当准确模态。'
    '与旧z=−0.5/0/0.5/1的逐参数对比在 `old_new_position_comparison.csv`；旧z=−1不作为对照，装配/采集状态变化不归因于算法误差。',
    '', '## 7. 与任务计划不同之处、限制和后续判断','',
    '1. 数据实际短于原18周期计划，B也短于36周期；按实际公共完整周期分析，没有补造数据。部分记录按前述统一瞬态规则多丢周期。',
    '2. T只在秩3控制子空间上可定，报告为最小范数表示；非零位置增加常量项处理。控制器系数来自数据精确恢复，没有硬件导出系数，因此独立A/B恢复一致性另行核对。',
    '3. M1使用与本轮其他IV相同的80步、10⁻⁵、0.05设置；旧实验C的M1曾用10⁻⁶、0.1。模型没有变，停止规则差异已明确，没有用旧结果冒充新估计。',
    '4. 定理2只直接适用于符合条件的参考IV，不能给所有方法直接套同一个公式。M1/M2和M5n的u来源分别说明；M3w冻结随机权、B谱有限样本与控制器恢复误差等未被完整覆盖。',
    '5. C区间采用1000次成对周期bootstrap而非声称精确无偏；发现明显比值偏移则标注，不删频点。不按参数拟合好坏正则化FRF。',
    '6. 归一化残差给出理想χ²期望/区间及条件参数扰动校准；不把已估计协方差与共享模型误差下的统计量说成精确χ²。M5n全参数u校正未提供，仅给固定模型参照和周期重估u。',
    '7. D改为幅值扫描，不能按重复试验宣称覆盖率；A nominal没有重复计数。E的z=0复用A，噪声借用z=0只用于描述。',
    '8. 没有改模型结构、频带或固定时延，没有因未通过判据去加阶数/缩频带。失败和严格score未过结果全保留。',
    '', '**对论文的判断：可用第二阶阻尼结果支持“控制器相关偏差及参考IV减偏”的局部实验论点；暂不能写“全模态一致性及不确定度覆盖已在实机验证”。** '
    '下一步应先区分激励幅值/方向引起的实际响应变化与测量噪声，并改善支承峰附近的有效输入信噪比。'
    '不应先通过修改成功标准、只报第二阶或把不稳定C参照当真值来完成论文叙事。',
    '', '## 8. 交付与复现','',
    '- `README.md`：逐脚本运行命令及统计口径。',
    '- `prepared/`、`spectra/`、`controllers/`：统一时域记录、DFT、控制器；原始MAT只读。',
    '- `fits/`：94组完整估计与权协方差；`analysis/`：全部u、偏差、FRF/噪声/检验、bootstrap和CSV。',
    '- `figures/`：20幅PNG/PDF及来源manifest，重点图为 `core_controller_zeta1_zeta2_f2`、`all_modal_controller`、`C_*_zoom_8_18`。',
    '- `analysis/verification.json`：控制器独立恢复、正交结构、Gaussian二次得分矩代数检查、原文件哈希。独立数值检查的协方差相对差约1.25%，均值偏差在1.52个Monte Carlo标准误内。',
    '- `logs/`：每步日志；模型数据参数来源可追溯。本次没有新硬件实验或修改论文正文。','']
    (HERE/'新实测数据分析报告.md').write_text('\n'.join(out),encoding='utf-8')
    js(HERE/'run_manifest.json',dict(date='2026-10-20',python=sys.version,platform=platform.platform(),numpy=np.__version__,raw_files=len(inv),unique_records=26,
        main_fits=len(fits),failures=[dict(id=r['id'],method=r['method']) for r in fails],strict_score_failures=sum(not r.get('strict_score_pass',False) for r in fits),
        methods=['M1','M2','M3','M3w','M5n'],model_parameters=158,band_hz=[3,25],fs=FS,P0=P0,bootstrap_FRF=1000,bootstrap_SML_each=99,
        claim_update='partial support for second-mode damping bias; full predeclared invariance and uncertainty coverage not met',
        next_action='diagnose operating-condition dependence and effective resonance input SNR before claiming experimental consistency',original_integrity=v['source_integrity']))
    print('REPORT_SAVED',len(out),'blocks',flush=True)
if __name__=='__main__':main()
