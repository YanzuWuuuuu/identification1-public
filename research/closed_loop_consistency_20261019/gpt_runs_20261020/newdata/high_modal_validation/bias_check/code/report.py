from common import *
import csv
def readcsv(name):
    with (HERE/'analysis'/name).open(encoding='utf-8-sig') as h:return list(csv.DictReader(h))
def val(r,k):return float(r[k])
def fmt(x):return f'{float(x):.5g}'
def tab(head,rows):
    return '| '+' | '.join(head)+' |\n| '+' | '.join(['---']*len(head))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)+'\n'
def ident(r):return r['id'].replace('_pppp','++++')
def nm(r):return {'H300':{2:'282',3:'333',4:'359'},'H493':{1:'493'},'L3_25':{1:'L1',2:'L2',3:'L3'}}[r['window']][int(r['component'])]
def pm(r,k):return fmt(r[k+'_percent'])+' ± '+fmt(r[k+'_SE_percent'])
def main():
    high=readcsv('high_frequency_tests.csv');low=readcsv('low_frequency_all_modes.csv')
    rep=readcsv('reproducibility.csv');pairs=readcsv('repeat_pairs.csv');amp=readcsv('D_amplitude.csv')
    excluded=readcsv('excluded_windows_status.csv');summary=rd(HERE/'analysis/summary.json')
    validation=rd(HERE/'analysis/validation.json');parts=[];add=parts.append
    add("""# 高频模态的偏差预测与复现性报告

2026-10-20。原数据、原高频窗口模型和既有结果只读；本报告是低频闭环偏差主结果的高频补充，不能替代低频检验。

## 1. 结论：量级降低得到支持，噪声尺度下的吻合多数未通过

- H300与H493的12组完整记录计算及82个删除周期折已完成。164次新M1/M3拟合经原协议数值精化后均收敛、未贴边且通过严格score；没有删除失败端点。
- 两种公式给出的高频偏差都明显小于低频第二阶。曲率公式下，频率预测最大绝对值约0.000156%，阻尼约0.0694%；相应实测配对差最大约0.00574%、0.372%。后者虽小，但不是自动可忽略。
- 按预先固定的检测/吻合规则，两种公式均为：**频率4/24通过，阻尼4/24通过**。频率24项预测全部低于检测门槛，阻尼20/24项低于门槛；另外4项可检测的阻尼预测也未通过联合z检验。因此不能写“高频实测差异已在不确定度内等于零”。
- 约282/333/359/493Hz的跨运行频率复现性标准差分别约0.0483/0.0635/0.0260/0.0371Hz，比同记录谱的条件sandwich标准不确定度大约14.6/6.0/11.4/18.9倍。噪声不确定度仍不是实验总不确定度。
- D组四个模态的频率随激励增大下降，但每档仅一次运行；不能把幅值依赖和时间/温升漂移分开。

可用于论文的稳妥表述为：**高频的理论反馈偏差及实测方法差异在量级上均显著减弱；然而当前局部模型下，多数实测配对差仍超过周期重估给出的检测阈值，精确的噪声尺度吻合未获支持。运行间复现性是必须单列的不确定度来源。**

## 2. 固定范围、端点和数值实现

仅用A与C++++六条记录，H300为220–420Hz、原4分量/200参数；H493为450–550Hz、原1分量/74参数。H300约263Hz分量是宽背景，不列为机械模态，但其留数、极点和瞬态全部参与预测矩阵。其余约282、333、359、493Hz按原分量对应。未新增数据、未改变窗口、模型阶次、D0/D2、固定通道时延、精确ZOH或原参数边界。

完整记录引用原diagnostics的数值精化端点，核对原求解成功、精化score<10⁻⁶和非边界；24个M1/M3端点全部合格。没有把原协议的提前停止状态覆盖为精化状态。删除周期重估沿用每方法原OLS初始化、80步IV求解与原界限，随后执行既有最多25步数值精化（目标score<10⁻⁸）。每折单独保存原停止状态、原score、精化历史和最终score。

同一控制器的A/C分量频率对应保持原顺序；H300宽背景不与其余三个候选互换。这里“板模态候选”仍不构成其机械归属的独立证明。

### 2.1 同记录噪声谱和两种预测

在窗口两端外扩2Hz的频点上，用该完整记录/该折M3参数减去包含瞬态的模型输出。只取非激励频点j mod P≠0，0.5Hz全宽内求

\\[
\\widehat\\Phi_v(f)=\\frac{1}{N n_f}\\sum_{j\\in\\mathcal N_f}\\widehat v_j\\widehat v_j^H,
\\qquad C_k=N\\widehat\\Phi_v(f_k).
\\]

DFT未归一化；C_k是全长DFT残差的复协方差，不是单周期均值协方差。不加收缩或额外正则；所有局部协方差均通过Cholesky检查。

其余沿用低频流程。以M3端点代入G，令F=(I+KG)⁻¹K，U_d=(I+KG)⁻¹D：

\\[
M_0=\\operatorname{Re}\\sum J_d^H J_d,\\quad
M=M_0+\\mathbb E\\operatorname{Re}\\sum J_v^H J_v,\\quad
b_p=-\\operatorname{Re}\\sum_k\\operatorname{tr}(G_pF C_k).
\\]

其中J_d包含全部瞬态列，G_p是传函导数，J_v由−G_p F v构造。曲率C_curv由固定F和噪声谱的G二阶导数构造（避免与噪声协方差C_k混淆）。两种预测为

\\[
\\Delta_0=M^{-1}b,\\qquad
\\Delta_c=(M-C_{\\mathrm{curv}})^{-1}b.
\\]

外层差分步长沿用2×10⁻⁴。所有结果按100 expm1(Δlog f)或100 expm1(Δlog ζ)报告。实测配对差为100(参数_M1/参数_M3−1)，不是相对真值的误差。

### 2.2 周期jackknife与判定

六条记录的可用周期数分别为7、7、8、7、6、6；两个窗口共82个删除周期折。每折删一个整周期，将其余周期顺序拼接，重建DFT和噪声谱，再计算实测与预测。对x_i分别取实测δ_i、预测p_i及同折δ_i−p_i：

\\[
s_J(x)=\\sqrt{\\frac{P-1}{P}\\sum_i(x_i-\\bar x)^2},\\qquad
z=\\frac{\\delta_{\\rm full}-p_{\\rm full}}{s_J(\\delta-p)}.
\\]

图中点仍为完整记录估计，误差棒为1SE，未用折均值替代点估计。联合SE保留实测/预测的同折协方差，不能用两个边际SE独立相加代替。

任务中的“2SE”沿用低频主判据，预先明确为实测配对差的检测能力2s_J(δ)：若|p|<2s_J(δ)，检查|δ|<2s_J(δ)；否则检查联合|z|<2。预测自身SE另列，不用它事后切换检测分支。两种预测分别使用同一规则。

局限：只有6–8周期，±2SE是任务指定的近似判据，非精确小样本95%覆盖保证；48项检验也非独立，未做事后多重性修正。删除中间周期会增加人工拼接边界，原单段瞬态模型未增加新项，本轮沿用低频做法；因此jackknife本身也有近似性。K与固定时延未逐折重新辨识，模型类错误、运行间漂移和周期性非线性不在噪声协方差中。

### 2.3 独立实现核查
""")
    add(f"从prepared重建的f/U/Y/D/K与原高频谱逐元素相同。分块计算与直接稠密公式的M、b相对误差约{validation['chunk_vs_dense']['M']:.2g}/{validation['chunk_vs_dense']['b']:.2g}，曲率预测解相对差约{validation['chunk_vs_dense']['delta_MminusC']:.2g}；局部谱前缀和与逐点求平均相对差约{validation['local_spectrum_prefix_vs_direct']:.2g}。外层步长减半的曲率矩阵相对变化约{100*validation['curvature_step_2e4_vs_1e4_relative']:.4g}%，没有改变主步长。\n")
    add("## 3. 高频判定和低频对照\n\n### 3.1 总体计数\n")
    add(tab(['公式','参数','比较数','预测低于检测门槛','通过预定规则'],[
        ['M⁻¹b' if r['formula']=='M' else '(M−C)⁻¹b',r['parameter'],r['n'],r['nearzero'],r['formal_pass']] for r in summary['counts']]))
    agg=[]
    for par in ('f','zeta'):
        for name in ('282','333','359','493'):
            rr=[r for r in high if nm(r)==name and r['parameter']==par and r['formula']=='MminusC']
            agg.append([name,par,fmt(max(abs(val(r,'predicted_percent')) for r in rr)),
                 fmt(max(abs(val(r,'observed_percent')) for r in rr)),
                 str(sum(r['formal_rule_pass']=='True' for r in rr))+'/6'])
    add(tab(['模态/Hz','参数','最大绝对预测 %','最大绝对实测差 %','通过数'],agg))
    add("""预测小与实测不吻合不矛盾。两种方法使用同一记录，可消去跨运行漂移的混淆，但**不能保证模型失配在两种估计方程中抵消**：M1使用全部频点的实测输入工具，M3使用激励谱线参考工具，失配时它们可对应不同的近似参数。当前结果与这种可能性相容，也可能受有限记录边界、周期非线性、K/固定时延误差等影响；本轮没有做来源分离，不能把不吻合唯一归结到某一种原因。

这不是“高频一致性定理被否定”，因为实测局部模型正确等前提尚未成立；同样也不能把多数未通过的结果改写为定理得到实验证明。

### 3.2 带宽参考与低频第二阶对照

没有使用仿真60Hz替代实际带宽。利用已有C四符号FRF与恢复控制器的KG，取三个非零奇异值在3–90Hz内最后一个向下穿越1的位置：
""")
    add(tab(['控制器','三个交越位置/Hz','图中参考范围/Hz'],[
        [r['controller'],' / '.join(fmt(x) for x in r['last_downward_crossings_hz']),fmt(r['lower_hz'])+'–'+fmt(r['upper_hz'])] for r in summary['bandwidth']]))
    add("它是数据估计的主回路单位增益交越范围，不是每个逻辑轴的-3dB带宽或多变量稳定裕度。低频第一、二阶位于范围以下，第三阶在lowgain交越附近；282–493Hz在其上方。高频共振仍可能放大局部反馈，不能将此范围当作严格的偏差分界。\n\n下表按曲率主公式，所有数为百分比及1SE；低频原结果不改。\n")
    rr=[r for r in low if r['parameter']=='zeta' and r['component']=='2' and r['formula']=='MminusC']+[r for r in high if r['parameter']=='zeta' and r['formula']=='MminusC']
    rr.sort(key=lambda r:(IDS.index(r['id']),val(r,'frequency_hz')))
    add(tab(['记录','模态/Hz','实测 ± SE','预测 ± SE','联合SE','z'],[
        [ident(r),fmt(r['frequency_hz']),pm(r,'observed'),pm(r,'predicted'),fmt(r['joint_SE_percent']),fmt(r['z_joint'])] for r in rr]))
    add("""### 3.3 全部高频数值：不省略未通过结果

以下实测SE、预测SE和联合SE均为百分比单位。`零效应`分支检验实测是否小于2SE；此分支列出的z只作描述，不替代预设分支判定。
""")
    for par,lab in [('f','频率'),('zeta','阻尼')]:
        for form,title in [('M','M⁻¹b'),('MminusC','(M−C)⁻¹b')]:
            add(f"#### {lab}，{title}\n")
            rr=[r for r in high if r['parameter']==par and r['formula']==form]
            add(tab(['记录','模态','实测 ± SE','预测 ± SE','联合SE','z','分支','通过'],[
                [ident(r),nm(r),pm(r,'observed'),pm(r,'predicted'),fmt(r['joint_SE_percent']),fmt(r['z_joint']),
                 '零效应' if r['near_zero_prediction']=='True' else '可检测预测','是' if r['formal_rule_pass']=='True' else '否'] for r in rr]))
    add("""## 4. 复现性与条件统计不确定度

同一控制器的两个运行(A、C++++)做差d_g，再跨三个控制器池化：

\\[
s_{\\mathrm{rep}}=\\sqrt{\\frac{\\sum_{g=1}^3d_g^2}{2\\times3}}.
\\]

这是以每个控制器有各自均值、两次运行可交换且具有共同方差为前提的组内复现性SD估计，自由度只有3。未先扣除共同A→C位移，因而保留一致方向漂移；另外列平均差和中心化差的SD/√2。它是当前采集条件下的经验离散，不能仅凭三对记录建立完整可溯源总不确定度预算。

条件sandwich使用同记录谱，且保留全部干扰参数：

\\[
A=\\operatorname{Re}\\sum Z^H Z,\\quad
Q=\\tfrac12\\operatorname{Re}\\sum Z^H(N\\widehat\\Phi_v)Z,\\quad
\\mathrm{Cov}(\\widehat\\beta)=A^{-1}QA^{-T}.
\\]

这里沿用原高频实现的正确模型极限bread（ZᴴZ）；实际ZᴴJ版本也存入原始数组作核对，但未替代此统一口径。它不是模型失配稳健的总误差估计。log标准差按f或ζ乘回物理量。B谱的结果在同一精化端点重算，以避免混用端点。表中的u_RMS为六条记录各自标准不确定度的平方平均开方。
""")
    add(tab(['模态','参数/单位','复现性SD','同记录u_RMS','B谱u_RMS','SD/同记录u','平均C−A','中心化差SD/√2'],[
        [nm(r),'f/Hz' if r['parameter']=='f' else 'ζ/百分点',fmt(r['replication_SD']),fmt(r['conditional_same_record_u_RMS']),
         fmt(r['conditional_B_u_RMS']),fmt(r['replication_over_same_record']),fmt(r['mean_CminusA']),
         fmt(r['centered_pair_difference_SD_over_sqrt2'])] for r in rep]))
    add("### 4.1 每对运行差值\n")
    add(tab(['控制器','模态','参数','C−A (Hz或阻尼百分点)','相对差 %','同记录合成u'],[
        [r['controller'],nm(r),r['parameter'],fmt(r['difference_CminusA']),fmt(r['relative_difference_percent']),fmt(r['combined_same_record_u'])] for r in pairs]))
    add("""所有四个频率的复现性明显大于条件噪声u；阻尼的比值约1.7–9.5，并非每个模态都差两个数量级。约333Hz的同记录u与B谱u相近，其他多项同记录u更大。不能一概声称同记录谱总比B谱更保守。

同一记录M1−M3的联合SE可以远小于各自参数u，因为两种估计的随机误差高度相关；因此“运行间复现性大”不能被直接加到同一记录配对差SE中来使检验通过。本报告将两类比较分开。

## 5. D组幅值趋势：只作测量条件说明
""")
    arows=[]
    for gain in (.6,.8,1.,1.2,1.5):
        rr=sorted([r for r in amp if val(r,'gain')==gain],key=lambda r:val(r,'f'))
        arows.append([str(gain)]+[fmt(r['f'])+' ('+fmt(r['zeta_percent'])+')'+(' †' if r['strict_score_pass']!='True' else '') for r in rr])
    add(tab(['扰动系数','282Hz f(ζ%)','333Hz f(ζ%)','359Hz f(ζ%)','493Hz f(ζ%)'],arows))
    add("""†为原D估计尚未通过严格score，保留其状态，不混入主高频判定：H300的0.8/1.0/1.2三条原端点。所有D端点达到原停止规则且未贴边。1.5引用原A_nominal，未把同一文件重复当独立记录。

从0.6到1.5，四个频率约变化−0.330、−1.353、−0.293、−0.200Hz。每个幅值仅一次记录，采集时间/温升、幅值、激励实现等没有独立随机化；故不能单独辨识幅值非线性或认定片簧软化。不得将这些数据与A/C复现性合并后宣称样本量增加。

## 6. 排除窗口及失败/异常

H100和H700不参与本报告结论，下表只列原协议状态：
""")
    add(tab(['窗口','方法','次数','求解失败','贴边','严格score未过'],[
        [r['window'],r['method'],r['n'],r['failed'],r['bound'],r['strict_fail']] for r in excluded]))
    add("""H100的未加权方法原本可数值收敛，但其机械来源和本次范围不满足纳入要求；H700五方法全部存在求解失败。H100原有数值精化不改变其排除地位。

本次主窗口：24个完整记录端点合格，82折/164个新拟合端点全部合格；82折谱与预测无数值失败。原步长停止后有93/164个端点未过严格score，经事先规定的原数值精化流程后全部通过，最大最终score约9.98×10⁻⁹；不是把原停止成功直接等同于严格收敛。这里“合格”仅指数值求解，**不等于96项参数×公式检验通过**。每条未通过比较仍列在第3.3节和CSV中。任何原始提前停止、精化步数和谱条件数均留存，不只保存通过检验的记录。

汇总程序初次控制台打印遇到NumPy整型JSON序列化错误，结果文件已经正确写入；修复了该打印转换后完整重跑汇总，未重选数据或改变计算结果。

## 7. 论文图草稿及写作边界

- `paired_damping_frequency_MminusC`：原曲率主公式，A和C两面板，每个面板包含低频3模态+高频4模态；颜色为控制器，实心圆为实测、空心方块为预测，各带1SE。
- `paired_damping_frequency_M`：相同数据的M⁻¹b预测，不择优替换主图。
- `high_frequency_closeup`：高频频率与阻尼线性纵轴局部图，避免低频的大偏差把高频的不吻合遮住。
- `reproducibility_vs_noise`：复现性与两种噪声谱条件u并列。

主图横轴为实际M3模态频率（log轴，9–550Hz，无位置抖动），纵轴用有符号symlog，±0.3%内线性；同一模态的实测/预测锚定相同频率。低频仍引用原结果，不借高频数值精化修改低频score异常。误差棒是记录内jackknife，不能代表运行间复现性。

写作上保留两条不同层级的结论：

1. 支持量级对照：控制带宽附近的低频第二阶存在百分之几到几十的方法差异；高频实测阻尼差全部低于0.4%，预测更小。
2. 不支持精确吻合：高频多数配对差仍超出当前记录内误差棒，不能宣称全频段“预测与实验一致”或仅凭小百分比确认一致性。局部模型失配、周期性畸变、有限记录瞬态及固定动态误差仍是未分离的来源。

下一步可据此组织论文的“量级支持+测量不确定度局限”；若要声称噪声尺度下的定量吻合，需要另行解决这些未解释差异，不能在本任务中改模型或放宽判据凑通过。

## 8. 可复现文件

`README.md`列命令；`code`为新分析代码，`baseline`为未改动的原模型与时延副本；
`raw`保留全部拟合JSON、每折谱和预测矩阵NPZ；`analysis`含完整表、检验、复现性和低高频对照；
`figures`为PNG 300dpi/PDF。`protected_hashes_before.json`和`analysis/final_audit.json`记录已有文件的只读保护核查：2406个已有文件均未改变或丢失，原模型与时延副本逐字节一致。
""")
    target=HERE/'高频偏差与复现性报告.md';target.write_text('\n'.join(parts),encoding='utf-8');print('REPORT',target)
if __name__=='__main__':main()
