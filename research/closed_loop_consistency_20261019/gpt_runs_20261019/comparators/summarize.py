import pathlib,json,csv,numpy as np
from scipy.stats import spearmanr
HERE=pathlib.Path(__file__).resolve().parent
groups={};files={}
for folder in ['sml_results','indirect_results','indirect_ivgn_results']:
    for path in (HERE/folder).glob('*.npz'):
        if 'sensitivity'in path.name:continue
        a=np.load(path);info=json.loads(str(a['info']));key=(info['method'],info['P'])
        groups.setdefault(key,[]).append((a['relative_error_percent'],info));files[(info['method'],info['P'],info['r'])]=a['relative_error_percent']
rows=[]
for (method,P),vals in sorted(groups.items()):
    e=np.array([v[0]for v in vals]);finite=np.all(np.isfinite(e),axis=1);ef=e[finite]
    row=dict(method=method,P=P,n=len(vals),finite=int(finite.sum()),solver_failures=sum(not v[1]['success']for v in vals),
         score_failures=sum(v[1].get('normalized_score',float('inf'))>1e-5 for v in vals),
         mean_seconds=float(np.mean([v[1].get('seconds',0.)for v in vals])),
         mean_iterations=float(np.mean([v[1].get('iterations',v[1].get('njev',0))for v in vals])))
    for metric,v in [('mean',np.mean(ef,0)),('se',np.std(ef,0,ddof=1)/np.sqrt(len(ef))),('sd',np.std(ef,0,ddof=1)),('rmse',np.sqrt(np.mean(ef**2,0)))]:
        row[metric]=v.tolist()
    row['all_endpoints_including_failed']=True;rows.append(row)
(HERE/'summary.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
with (HERE/'summary.csv').open('w',newline='',encoding='utf-8-sig') as f:
    names=['method','P','n','finite','solver_failures','score_failures','mean_seconds','mean_iterations']+[f'{metric}_{param}{i+1}' for metric in ['mean','se','sd','rmse'] for param in ['f','zeta'] for i in range(3)]
    wr=csv.DictWriter(f,fieldnames=names);wr.writeheader()
    for row in rows:
        d={k:row[k]for k in names[:8]}
        for metric in ['mean','se','sd','rmse']:
            for j,key in enumerate([f'{metric}_{param}{i+1}' for param in ['f','zeta']for i in range(3)]):d[key]=row[metric][j]
        wr.writerow(d)
conds=[]
for P in [4,8,16]:
    vv=[];qq=[];errs=[]
    for p in (HERE/'indirect_data').glob(f'P{P:02d}_r*.npz'):
        a=np.load(p);r=int(p.stem.split('r')[1]);cc=a['condition_Tud'];vv.extend(cc.tolist());qq.append(float(np.max(cc)))
        er=files[('M6u_IVGN',P,r)];errs.append(float(np.linalg.norm(er[3:])))
    co=dict(P=P,all_line_percentiles=dict(zip(['min','p50','p95','p99','max'],np.percentile(vv,[0,50,95,99,100]).tolist())),
        replicate_max_percentiles=np.percentile(qq,[0,50,95,100]).tolist(),
        spearman_maxcond_vs_damping_error=float(spearmanr(qq,errs).statistic),
        number_lines_condition_over_1e3=int(np.sum(np.array(vv)>1e3)),number_lines=len(vv));conds.append(co)
(HERE/'conditioning.json').write_text(json.dumps(conds,indent=2),encoding='utf-8')
availability=[]
for P in [1,2,4,8,16]:
    n=20 if P==16 else 40
    if P<=4:availability.append(dict(method='M5',P=P,n=n,status='unavailable',reason='No sample covariance at P=1; residual covariance rank <= P-1 < 4 at P=2,4',rank_upper_bound=max(0,P-1)))
    if P in [4,8,16]:availability.append(dict(method='M6',P=P,n=n,status='unavailable',reason='Full 16x16 FRF sample covariance rank <= P-4 < 16',rank_upper_bound=P-4))
(HERE/'availability.json').write_text(json.dumps(availability,indent=2),encoding='utf-8')
text=['# SML 与间接 FRF 对照：数值结果','', '下表为所有有限端点，包括到达迭代上限的端点；未删除异常值。误差均为 `100*expm1(log参数误差)`。均值±SE、SD、RMSE 均按三个阻尼分量列出。`M6*_IVGN` 为用户要求的梯度工具变量 Gauss–Newton 主结果；`M6*` 为同一目标函数的 LM/TRF 优化器对照。','',
'| 方法 | P | n | 阻尼均值±SE (%) | SD (%) | RMSE (%) | 求解失败/score未过 |','|---|---:|---:|---|---|---|---|']
for row in rows:
    fmt=lambda x:' / '.join(f'{v:.2f}'for v in x[3:])
    mm=' / '.join(f'{m:.2f}±{s:.2f}' for m,s in zip(row['mean'][3:],row['se'][3:]))
    text.append(f"|{row['method']}|{row['P']}|{row['n']}|{mm}|{fmt(row['sd'])}|{fmt(row['rmse'])}|{row['solver_failures']} / {row['score_failures']}|")
text+=['','完整六参数统计、耗时及迭代次数见 summary.csv/json。score 门槛为归一化全梯度最大分量 ≤1e-5；这与单纯的求解器成功标记分开。','', '## 间接逆矩阵的条件数','', '|P|逐线中位数|95%|99%|最大值|最大条件数与阻尼误差的 Spearman 相关|','|---:|---:|---:|---:|---:|---:|']
for co in conds:
    q=co['all_line_percentiles'];text.append(f"|{co['P']}|{q['p50']:.2f}|{q['p95']:.2f}|{q['p99']:.2f}|{q['max']:.2f}|{co['spearman_maxcond_vs_damping_error']:.3f}|")
(HERE/'results_tables.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
print(json.dumps(conds,indent=2));print('Saved tables',len(rows))
