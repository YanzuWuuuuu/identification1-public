from pathlib import Path
import json,csv
import numpy as np
HERE=Path(__file__).resolve().parent
METHODS={'dir_band_W1':'M1','dir_lines_W1':'M2','ref_band_W1':'M3','refD_band_W1':'M4'}
def metrics(arr):
    return dict(mean=np.mean(arr,0).tolist(),se=(np.std(arr,0,ddof=1)/np.sqrt(len(arr))).tolist(),sd=np.std(arr,0,ddof=1).tolist(),rmse=np.sqrt(np.mean(arr**2,0)).tolist()) if len(arr)>1 else {}
def main():
    records=[json.loads(p.read_text()) for p in sorted((HERE/'results').glob('*.json'))]
    rows=[]; summary=[]
    for P in [1,2,4,8,16]:
        for method,label in METHODS.items():
            rr=[r for r in records if r['P']==P and r['method']==method]
            if not rr:continue
            finite=[r for r in rr if r['error_percent'] is not None and np.all(np.isfinite(r['error_percent']))]
            conv=[r for r in finite if r['status']=='pole_step_converged']
            out=dict(P=P,method=method,label=label,n=len(rr),finite=len(finite),failures=sum(r['status']!='pole_step_converged' for r in rr),
                all_finite=metrics(np.array([r['error_percent'] for r in finite])),converged_only=metrics(np.array([r['error_percent'] for r in conv])),
                iterations_mean=float(np.mean([r['iterations'] for r in rr if r['iterations'] is not None])),
                seconds_median=float(np.median([r['elapsed_seconds'] for r in rr])),
                score_above_1e4=sum(r['score_max'] is not None and r['score_max']>1e-4 for r in rr))
            summary.append(out)
            for q in range(6):
                rows.append(dict(P=P,method=label,quantity='frequency' if q<3 else 'damping',mode=q%3+1,
                    **{key:out['all_finite'][key][q] for key in ('mean','se','sd','rmse')},failures=out['failures'],n=out['n'],iterations_mean=out['iterations_mean']))
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    with (HERE/'summary.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    md=['# 实验 A：原始 IV 迭代的周期数扫描','',
        '主表包含所有有限的最终迭代端点（包括达到迭代上限者），避免剔除失败后只呈现好结果。失败按原算法 80 次上限及非有限/异常计数；“收敛”仅为原算法模态步长门槛，不等于统计一致性。收敛子集另见 summary.json。百分比=100 expm1(log误差)。',
        '', '|P|方法|阻尼偏差均值±SE（1/2/3，%）|SD（%）|RMSE（%）|失败/总数|平均迭代|', '|---:|---|---|---|---|---:|---:|']
    for r in summary:
        a=r['all_finite'];fmt=lambda xs:' / '.join(f'{x:.2f}' for x in xs)
        ms=' / '.join(f'{m:.2f}±{se:.2f}' for m,se in zip(a['mean'][3:],a['se'][3:]))
        md.append(f"|{r['P']}|{r['label']}|{ms}|{fmt(a['sd'][3:])}|{fmt(a['rmse'][3:])}|{r['failures']}/{r['n']}|{r['iterations_mean']:.1f}|")
    md+=['','## 频率相对误差','', '|P|方法|频率偏差均值±SE（1/2/3，%）|SD（%）|RMSE（%）|','|---:|---|---|---|---|']
    for r in summary:
        a=r['all_finite']; ms=' / '.join(f'{m:.4f}±{se:.4f}' for m,se in zip(a['mean'][:3],a['se'][:3]));fmt=lambda xs:' / '.join(f'{x:.4f}' for x in xs)
        md.append(f"|{r['P']}|{r['label']}|{ms}|{fmt(a['sd'][:3])}|{fmt(a['rmse'][:3])}|")
    p1={m:[r for r in records if r['P']==1 and r['method']==m] for m in list(METHODS)[:2]}
    eq=max(np.max(abs(np.array(a['beta'])-np.array(b['beta']))) for a,b in zip(*p1.values()))
    md+=['',f'P=1 时 M1/M2 最终参数最大绝对差：{eq:.3g}；所有频点均为激励谱线，符合恒等检查。','',
      '## 原始异常清单','', '|P|重复|方法|状态|相对矩分数|阻尼相对误差（%）|','|---:|---:|---|---|---:|---|']
    for r in records:
        if r['status']!='pole_step_converged':
            md.append(f"|{r['P']}|{r['replicate']}|{METHODS[r['method']]}|{r['status']}|{r['score_max']}|{r['error_percent'][3:] if r['error_percent'] else r.get('exception')}|")
    (HERE/'report_a.md').write_text('\n'.join(md),encoding='utf-8')
    np.savez_compressed(HERE/'all_endpoints.npz',P=np.array([r['P'] for r in records]),method=np.array([r['method'] for r in records]),replicate=np.array([r['replicate'] for r in records]),beta=np.array([r['beta'] if r['beta'] else [np.nan]*126 for r in records]),error_percent=np.array([r['error_percent'] if r['error_percent'] else [np.nan]*6 for r in records]),status=np.array([r['status'] for r in records]))
    print('summarized',len(records),'fits')
if __name__=='__main__':main()
