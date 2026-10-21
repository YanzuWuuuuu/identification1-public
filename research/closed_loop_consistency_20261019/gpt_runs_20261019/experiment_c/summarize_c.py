import json,csv,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
rows=[json.loads(p.read_text(encoding='utf-8')) for p in sorted((HERE/'raw').glob('*.json'))]
lookup={(str(x['z']),x['P'],x['start_cycle_zero_based'],x['sampling'],x['method']):x for x in rows}
ZS=['-0.5','0','0.5','1'];PS=[1,2,4,8]
def get(z,p,method,sampling='zoh',start=1):return lookup.get((z,p,start,sampling,method))
def status(x):
    if x is None:return 'NA'
    if x.get('status')=='covariance_rank_deficient':return '秩不足'
    if x.get('method','').startswith('M5'):
        if not x.get('success'):return '未收敛'
        return '终止/score通过' if x.get('score_normalized',9)<1e-6 else '终止/score未过'
    if x.get('exception'):return '异常'
    if not x.get('converged_pole_step'):return '迭代上限'
    return '收敛/score通过' if x.get('strict_score_pass') else '步长收敛/score未过'
def fmt(a,n=4):return '/'.join(f'{x:.{n}f}' for x in a)
def good(x):return x is not None and x.get('converged_pole_step') and x.get('strict_score_pass')
fields=['z','P','start_cycle_zero_based','sampling','method','status','iterations','score_normalized','elapsed_s','frequency_1_hz','frequency_2_hz','frequency_3_hz','damping_1_percent','damping_2_percent','damping_3_percent']
with (HERE/'experiment_c_all_fits.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
    for x in rows:
        y={k:x.get(k) for k in fields};y['status']=status(x)
        for i in range(3):y[f'frequency_{i+1}_hz']=x.get('frequency_hz',[np.nan]*3)[i];y[f'damping_{i+1}_percent']=x.get('damping_percent',[np.nan]*3)[i]
        w.writerow(y)

out=['# 实验 C：实测 ZOH 复算与周期数扫描','',
'## 范围和实现','',
'四个位置 z = −0.5、0、0.5、1 mm；P0 = 30000，fs = 5000 Hz。每个主结果取第 2 个周期起的前 P = 1、2、4、8 个完整周期。未使用 z = −1。',
'模型保留原代码的 3 个支承模态、D0、D2、固定逐通道时延及每模态每输出两个瞬态系数，共 158 个实参数。仅支承模态改为 ZOH 精确求值；瞬态分母仍为精确离散分母，分子仍为 tau0 + tau1 z^-1。固定时延矩阵统一减 100 µs。D0 和 D2 背景不作离散转换，因此这里的“ZOH 版”准确含义是模态块精确采样，不是整个背景加时延模型的精确离散化。减半采样周期是既有有效时延的重新分配假设，不是重新辨识出的物理时延。',
'M1 用 3–25 Hz 全部 DFT 点，M2 仅用激励谱线；M3 使用全部频点，工具输入为 (I+Khat G)^−1 D。旧 run_real.py 的 ref_K 实际只用谱线，本次已按任务更正为 ref_band；比较 CT/ZOH 时两者都采用相同新版频点规则。',
'每个窗口独立拟合 ARX(6,6) 控制器，K = −Q Kv；每方法/采样模型从相同模态频率 10.6/11.7/16.1 Hz、阻尼 3% 开始，再用各自模型线性初始化留数与瞬态。不以真值/其他方法结果初始化。IV-GN 保留原公式与 0.1 对数极点步长上限，最多 80 次，停止阈值 1e−6；额外记录完整归一化 score < 1e−6 和条件数。',
'实测没有真实频率/阻尼标签，下面的量是“参数估计和方法差异”，不是偏差或准确率。参数分量保持求解器标签；少周期时前两分量可能靠近或互换。',
'', '## 数值校验','']
v=json.loads((HERE/'verification.json').read_text())
out += [f"新 ZOH 函数与 cl2.zoh_2nd 最大绝对差 {v['zoh_max_absolute_vs_cl2']:.2e}；无噪声采样状态空间数据的频域恒等式（已估计边界瞬态）相对残差 {v['noiseless_DFT_identity_relative']:.3e}。模态雅可比与中心差分的最大相对差 {max(v['modal_jac_relative_error']):.3e}，对应沿用原代码 1e−6 前向差分的一阶误差。",'']
changes=[]
for p,h in v['original_SHA256'].items():
    if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h:changes.append(p)
out += [f'原代码、原数据及固定时延文件的 SHA256 复核：{len(v["original_SHA256"])} 个文件，变化数 {len(changes)}。','']
out += ['## C 主表：ZOH 版，全周期扫描','',
'频率和阻尼各三项以 / 分开；阻尼单位为 %。M5 在 P=2、4 的原始四输出样本残差协方差秩至多为 1、3，不能求逆，故明确标记不可计算。M5_reg 用 10% 联合对角收缩，是另一估计器，只作诊断。', '',
'|z (mm)|P|方法|频率 (Hz)|阻尼 (%)|迭代/评价次数|状态|',
'|---:|---:|---|---|---|---:|---|']
for z in ZS:
    for p in PS:
        methods=['M1','M2','M3']+(['M5','M5_reg'] if p in [2,4] else ['M5'] if p==8 else [])
        for method in methods:
            x=get(z,p,method)
            if not x:continue
            out += [f"|{z}|{p}|{method}|{fmt(x['frequency_hz'])}|{fmt(x['damping_percent'])}|{x.get('iterations')}|{status(x)}|"]
out += ['', '## C 采样对照：P=8，旧 CT 与新 ZOH','',
'Δ为新 ZOH 减旧 CT。阻尼差采用百分点（pp），不把百分数差误写成相对误差。所有 P 的完整结果在 experiment_c_all_fits.csv。','',
'|z|方法|CT频率|ZOH频率|最大绝对频率差 (Hz)|CT阻尼 (%)|ZOH阻尼 (%)|最大阻尼差 (pp)|状态 CT→ZOH|','|---:|---|---|---|---:|---|---|---:|---|']
maxdf=maxdz=0.
for z in ZS:
    for method in ['M1','M2','M3']:
        a=get(z,8,method,'ct');b=get(z,8,method)
        df=np.max(np.abs(np.array(b['frequency_hz'])-a['frequency_hz']));dz=np.max(np.abs(np.array(b['damping_percent'])-a['damping_percent']))
        maxdf=max(maxdf,df);maxdz=max(maxdz,dz)
        out += [f"|{z}|{method}|{fmt(a['frequency_hz'])}|{fmt(b['frequency_hz'])}|{df:.6f}|{fmt(a['damping_percent'])}|{fmt(b['damping_percent'])}|{dz:.6f}|{status(a)}→{status(b)}|"]
out += ['',f'P=8 所有 12 对端点的最大变化：频率 {maxdf:.6f} Hz，阻尼 {maxdz:.6f} 个百分点。接近/达到迭代上限的结果仍列出但不当成确定收敛解。少周期失败端点之间较大差异不能解释成 ZOH 的物理效应。','',
'## 方法差异是否随 P 缩小','',
'下表给出三个模态中最大的 |ζ_method−ζ_M3|，单位为百分点。星号表示该比较至少一端未通过“步长收敛且 score < 1e−6”，不宜用于检验渐近阶数。','',
'|z|P|M1−M3 最大绝对阻尼差|M2−M3 最大绝对阻尼差|','|---:|---:|---:|---:|']
for z in ZS:
    for p in PS:
        a,b,c=get(z,p,'M1'),get(z,p,'M2'),get(z,p,'M3')
        def gap(x):return f"{np.max(np.abs(np.array(x['damping_percent'])-c['damping_percent'])):.4f}"+('' if good(x) and good(c) else '*')
        out += [f'|{z}|{p}|{gap(a)}|{gap(b)}|']
out += ['', 'P=8 时 M2 通常比 M1 更靠近 M3，但个别位置/模态例外；窗口扫描不是独立重复实验，且前几项有求解失败，不能从这四点推出严格 1/P 规律，更不能把 M3 当实测真值。M1 与 M3 的差在 z=0、0.5 的第二模态仍明显；z=1 的差较小且方法排序不同。', '',
'## P=8 不确定度：两半重估','',
'将第2–9周期拆成第2–5周期与第6–9周期，各自重估 ZOH-M3。表中 ± 为两半估计绝对差的一半，只是重现性量级，不是置信区间或经验证的标准误。只有两半都通过收敛及 score 才画实心误差棒。','',
'|z|P8频率|两半半差 (Hz)|P8阻尼 (%)|两半半差 (pp)|两半状态|','|---:|---|---|---|---|---|']
halves=[]
for z in ZS:
    a=get(z,4,'M3');b=get(z,4,'M3',start=5);c=get(z,8,'M3')
    df=np.abs(np.array(a['frequency_hz'])-b['frequency_hz'])/2;dz=np.abs(np.array(a['damping_percent'])-b['damping_percent'])/2
    halves.append(dict(z=z,half_frequency_difference_hz=df.tolist(),half_damping_difference_pp=dz.tolist(),valid=good(a) and good(b)))
    out += [f"|{z}|{fmt(c['frequency_hz'])}|{fmt(df)}|{fmt(c['damping_percent'])}|{fmt(dz)}|{status(a)}；{status(b)}|"]
out += ['', 'z=−0.5 与 z=0.5 的第一半 M3 达到迭代上限，因此这两位置目前不能从该两半方法给出可信误差棒；表中差值仅暴露计算不确定性。z=0 与 z=1 两半均收敛，可作为描述性重复性检查。','',
'## SML 的结果及限制','',
'使用共同的 comparators/sml_core.py：周期均值协方差除以 P(P−1)，包含 C_YU，优化白化残差时计入协方差随参数变化的导数；没有把冻结协方差的一步 GN 误称为完整 SML。保留与其他方法相同模态及瞬态类。P=2、4 原始 SML 结构上不可计算；10%收缩仅是标明的辅助对照。',
'P=8 原始 SML 的前两模态有明显不同解，部分达到 160 次评价上限；不能据此声称 SML 无偏或参考 IV 优于正确实现的总体 SML。这是当前样本协方差、有限频点、非线性优化与近邻模态条件下的实际结果，应在论文保留。实测无真值，模型适配/协方差采样误差/局部极值均可能参与，现有结果不能分离其来源。', '',
'## 对任务已有结论的支持程度','',
'- “全频点直接法在控制带宽内存在与参考IV不同的估计”：本批实测支持这种方法差异，尤其 z=0、0.5 的第二模态；没有真值，不能单靠此表认定哪一个是偏差。',
'- “只用谱线的偏差按1/P下降”：P8总体更接近参考法，但实测序列不呈统一单调规律且小P失败，不能独立证实渐近阶。',
'- “参考IV无偏且少周期可用”：实测不能证明无偏；P1/P2多处不收敛，所以本实现不支持‘只有M3/M4可用’这种无条件工程说法。',
'- “SML在P≥4无偏”：四输出逐谱线协方差在P=4仍秩不足；原始SML需要至少P≥5才可能满秩，满秩仍不等于有限样本无偏。',
'- “ZOH替换会改变实测结论”：在P8收敛结果中改变很小，不是此前直接法/参考法差异的主要来源。',
'- 所有失败端点、完整 beta、迭代记录、score 与耗时均保留；没有剔除异常后只报告好结果。','']
(HERE/'experiment_c_report.md').write_text('\n'.join(out),encoding='utf-8')
(HERE/'halves_uncertainty.json').write_text(json.dumps(halves,indent=2),encoding='utf-8')

fig,axes=plt.subplots(4,3,figsize=(12,12),sharex=True)
colors={'M1':'#af3b3b','M2':'#55823d','M3':'#245b9b'}
for zi,z in enumerate(ZS):
    for mi in range(3):
        ax=axes[zi,mi]
        for method in ['M1','M2','M3']:
            xs=[get(z,p,method) for p in PS];yy=[x['damping_percent'][mi] for x in xs]
            ax.plot(PS,yy,color=colors[method],alpha=.75,lw=1.3,label=method)
            for p,y,x in zip(PS,yy,xs):ax.plot(p,y,'o' if good(x) else 'x',color=colors[method],ms=5)
        if halves[zi]['valid']:
            x=get(z,8,'M3');ax.errorbar(8,x['damping_percent'][mi],yerr=halves[zi]['half_damping_difference_pp'][mi],color=colors['M3'],capsize=3)
        ax.set_xscale('log',base=2);ax.set_xticks(PS,labels=PS);ax.grid(alpha=.2)
        ax.set_title(f'z = {z} mm, component {mi+1}')
        if mi==0:ax.set_ylabel('Estimated damping (%)')
        if zi==3:ax.set_xlabel('Number of periods P')
axes[0,0].legend(fontsize=9)
fig.suptitle('Measured data: exact-ZOH modal blocks\nx = solver/score criterion not passed; error bars = half-record spread (descriptive)',fontsize=12)
fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(HERE/'experiment_c_damping_scan.png',dpi=180);plt.close(fig)
summary=dict(n_fits=len(rows),max_P8_ct_zoh_frequency_change_hz=maxdf,max_P8_ct_zoh_damping_change_pp=maxdz,
             original_files_changed=changes,by_method={})
for method in ['M1','M2','M3','M5','M5_reg']:
    rr=[x for x in rows if x['sampling']=='zoh' and x['start_cycle_zero_based']==1 and x['method']==method]
    summary['by_method'][method]=dict(n=len(rr),status_counts={s:sum(status(x)==s for x in rr) for s in sorted(set(status(x) for x in rr))})
(HERE/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
