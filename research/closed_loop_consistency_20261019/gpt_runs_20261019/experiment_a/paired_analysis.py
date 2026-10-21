from pathlib import Path
import json
import numpy as np
R=Path(__file__).resolve().parent
d=[json.loads(p.read_text()) for p in (R/'results').glob('*.json')]
rows=[]
for P in [1,2,4,8,16]:
    by={(r['method'],r['replicate']):r for r in d if r['P']==P}
    out={'P':P}
    for method in ['dir_band_W1','dir_lines_W1','refD_band_W1']:
        pairs=[(a,by[('ref_band_W1',r)]) for (m,r),a in by.items() if m==method]
        diff=np.array([np.array(a['error_percent'])-np.array(b['error_percent']) for a,b in pairs])
        out[method+'_minus_ref']=dict(mean=diff.mean(0).tolist(),se=(diff.std(0,ddof=1)/np.sqrt(len(diff))).tolist(),n=len(diff))
    refs=np.array([r['error_percent'] for r in by.values() if r['method']=='ref_band_W1'])
    drefs=np.array([r['error_percent'] for r in by.values() if r['method']=='refD_band_W1'])
    out['M4_M3_sd_ratio']=(drefs.std(0,ddof=1)/refs.std(0,ddof=1)).tolist()
    out['M4_M3_variance_ratio']=(drefs.var(0,ddof=1)/refs.var(0,ddof=1)).tolist()
    rows.append(out)
(R/'paired_analysis.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
md=['# 配对分析：周期数阶数与效率','',
 '每个 P 内所有方法使用完全相同数据和初始扰动。以下保留所有有限端点，包括主表已标明的迭代失败，因此 P=1、2 的差异不能直接解释为已收敛估计器的渐近偏差。',
 '', '|P|M2−M3第一模态阻尼相对误差差值（均值±SE，百分比点）|P×均值|M4/M3阻尼标准差比（1/2/3）|方差比（1/2/3）|', '|---:|---|---:|---|---|']
for r in rows:
    a=r['dir_lines_W1_minus_ref'];fmt=lambda x:' / '.join(f'{v:.2f}' for v in x)
    md.append(f"|{r['P']}|{a['mean'][3]:.3f}±{a['se'][3]:.3f}|{r['P']*a['mean'][3]:.3f}|{fmt(r['M4_M3_sd_ratio'][3:])}|{fmt(r['M4_M3_variance_ratio'][3:])}|")
md+=['','P≥2 时第一模态配对差的 P×均值约保持同一数量级，与谱线直接法相对参考法的 1/P 衰减相容；小 P 求根失败和非线性有限样本偏差仍存在，因此这不是渐近阶的证明。',
'','“不用控制器时方差约两倍”与“标准差约两倍”不能混用：P=8 时标准差比约2.1–2.3，对应方差比约4.2–5.2。控制器缺失通常降低效率，但没有与系统无关的固定倍率。']
(R/'paired_analysis.md').write_text('\n'.join(md),encoding='utf-8')
print('Saved paired analysis for',len(rows),'period counts')
