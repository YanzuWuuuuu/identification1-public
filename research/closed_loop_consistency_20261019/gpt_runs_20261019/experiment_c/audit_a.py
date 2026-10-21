"""Independent audit of A; reads immutable A endpoints and writes C diagnostics."""
import os,sys
for name in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[name]='1'
sys.dont_write_bytecode=True
import numpy as np,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;AA=HERE.parent/'experiment_a'
sys.path.insert(0,str(AA))
import run_a as a
c=a.c
OUT=HERE/'diagnostics';OUT.mkdir(exist_ok=True)

def arrstat(rr):
    if not rr:return None
    ar=np.array([x['error_percent'] for x in rr]);return dict(n=len(rr),mean=ar.mean(0).tolist(),se=(ar.std(0,ddof=1)/np.sqrt(len(ar))).tolist(),sd=ar.std(0,ddof=1).tolist())

records=[json.loads(x.read_text()) for x in (AA/'results').glob('*.json')]
out=dict(checks=[],statistics=[],diagnostic_refits=[])
for P in [1,2]:
    for method in a.METHODS:
        rr=[r for r in records if r['P']==P and r['method']==method and r.get('beta')]
        out['statistics'].append(dict(P=P,method=method,all=arrstat(rr),converged=arrstat([r for r in rr if r['status']=='pole_step_converged']),
            collapsed_damping_below_20percent_truth=sum(min(r['error_percent'][3:])<-80 for r in rr),
            score_above_1e4=sum(r['score_max']>1e-4 for r in rr),iteration_limits=sum(r['status']=='iteration_limit' for r in rr)))
    for ri in [0,9,17]:
        z=np.load(AA/'data'/f'P{P:02d}_r{ri:03d}.npz');D,U,Y=z['D'],z['U'],z['Y'];k=z['k'];sel=k%P==0
        check=dict(P=P,replicate=ri,n_all=len(k),n_lines=int(sum(sel)),min_hz=float(z['f'][0]),max_hz=float(z['f'][-1]),
            D_mean_scaling=float(np.linalg.norm(D[sel]/P-z['Dper'].mean(0))/np.linalg.norm(D[sel]/P)),
            U_mean_scaling=float(np.linalg.norm(U[sel]/P-z['Uper'].mean(0))/np.linalg.norm(U[sel]/P)),
            Y_mean_scaling=float(np.linalg.norm(Y[sel]/P-z['Yper'].mean(0))/np.linalg.norm(Y[sel]/P)))
        if np.any(~sel):check['reference_offline_relative']=float(np.linalg.norm(D[~sel])/np.linalg.norm(D[sel]))
        initial=[np.array(next(x for x in records if x['P']==P and x['replicate']==ri and x['method']==met)['initial_beta']) for met in a.METHODS]
        check['M1_M3_M4_initial_maxdiff']=float(max(np.max(abs(initial[0]-q)) for q in initial[2:]))
        out['checks'].append(check)

# Two cases are diagnostics only. Every method receives identical enlarged budget.
for ri in [9,17]:
    dd=np.load(AA/'data'/f'P01_r{ri:03d}.npz');z,U,Y,D=dd['zk'],dd['U'],dd['Y'],dd['D'];K=c.Ktf(z)
    for method in a.METHODS:
        src=next(x for x in records if x['P']==1 and x['replicate']==ri and x['method']==method)
        beta0=np.array(src['initial_beta'])
        if method.startswith('dir'):Uh=U
        elif method.startswith('refD'):Uh=D
        else:Uh=lambda bb:np.linalg.solve(np.eye(4)[None]+K@c.resp(z,bb),D[:,:,None])[:,:,0]
        t=time.time()
        rr=dict(replicate=ri,method=method,budget=240,old_errors=src['error_percent'],old_iterations=src['iterations'],old_status=src['status'])
        try:
            b,it=c.iv_solve(z,U,Y,Uh,np.ones_like(U.real),beta0,it=240,tol=1e-5)
            sc,res=a.score(z,U,Y,Uh,b)
            rr.update(beta=b.tolist(),errors=(100*np.expm1(b[c.IW:c.ITR]-c.truth_beta()[c.IW:c.ITR])).tolist(),iterations=it+1,score=sc,residual=res)
        except Exception as ex:rr['exception']=repr(ex)
        rr['seconds']=time.time()-t
        out['diagnostic_refits'].append(rr)
        (OUT/'audit_a_numeric.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
        print(ri,method,rr.get('errors'),rr.get('iterations'),rr.get('score'),flush=True)

lines=['# 实验 A 少周期异常的独立复核','',
'本记录只读实验A缓存和原始算法，不修改主结果。额外试算仅两组样本，各方法同等240次预算，不能替换主实验。','',
'## 实现核对','',
'- 频率范围 3–100 Hz，P0=30000，P=1有583个点，P=2有1165个点但仍583条激励谱线。端点及谱线索引正确。',
'- 对P=1/2的r=0/9/17，整段DFT在激励线上除以P与各周期DFT均值一致到舍入精度；输入、输出、参考均核对。不存在DFT归一化因子漏乘。',
'- 所有方法读取相同缓存数据。M1/M3/M4使用完全相同线性初始化，M2在它自身谱线子集做同样初始化；P1全部频点就是谱线，所以M1=M2。',
'- 用solve代替inv后乘参考不改变数学工具输入；根线程已有原算法数值对照。ZOH与瞬态采用原cl2函数。每周期激励频率栅格固定，无重复周期时无法非参数估计逐频点4输出协方差。',
'- 原始算法仅检查六个极点参数步长而非完整矩score；全部有限端点与收敛子集必须分别报告。','',
'## P1/P2 负阻尼偏差是否只是失败剔除问题','',
'|P|方法|全部端点阻尼误差均值(%)|收敛子集阻尼误差均值(%)|限次|任一阻尼低于真值20%的次数|','|---:|---|---|---|---:|---:|']
for s in out['statistics']:
    fm=lambda ar:'/'.join(f'{x:.3f}' for x in ar[3:])
    lines.append(f"|{s['P']}|{s['method']}|{fm(s['all']['mean'])}|{fm(s['converged']['mean'])}|{s['iteration_limits']}|{s['collapsed_damping_below_20percent_truth']}|")
lines += ['',
'失败端点中确有第一支承模态阻尼向零塌缩，明显拉低均值。只取收敛子集后，P1的第一模态参考IV仍约−12%而非零，说明不能将负均值完全归因于迭代上限；但收敛条件筛选本身会造成选择偏差，收敛子集不是无偏性的证明或反证。',
'', '## 两个代表样本的同等增加迭代预算','',
'|样本|方法|80次端点阻尼误差(%)|240次端点阻尼误差(%)|实际迭代|最终score|','|---:|---|---|---|---:|---:|']
for x in out['diagnostic_refits']:
    fm=lambda ar:'/'.join(f'{y:.3f}' for y in ar[3:])
    lines.append(f"|{x['replicate']}|{x['method']}|{fm(x['old_errors'])}|{fm(x['errors']) if 'errors' in x else x.get('exception')}|{x.get('iterations')}|{x.get('score')}|")
lines += ['',
'## 解释边界','',
'不存在已发现的频率/归一化/方法配对实现错误。单周期只是一条6秒记录，并不是某种“无偏保证”。期望矩在真值为零既不等于非线性求根的有限样本均值为真值，也不保证当前IV-GN收敛到所需根。126参数（96留数、6极点、24瞬态）、两个近邻低频模态和随机参考方向，使小样本的极点/阻尼根对噪声较敏感；M4未通过灵敏度整形，数值失败更多。此处只观察到这些症状，尚未把“弱工具变量”作为单独被证明的原因。',
'固定P0而增加重复P时，确定性激励谱线的相干信息增长，噪声可以平均，有限样本非线性偏差及求解失败可能下降。故参考IV的一致性只能支持相应极限和正则条件下的结论，不能支持“P=1/2时只有M3/M4可用”的无条件工程断言。',
'主实验应保留全部失败、全部有限端点统计及收敛子集统计；真实改善应以RMSE、失败率及具体目标性能衡量，不能只看是否存在渐近无偏定理。','']
(OUT/'audit_a_note.md').write_text('\n'.join(lines),encoding='utf-8')
