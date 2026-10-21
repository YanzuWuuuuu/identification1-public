"""Diagnostic reference-driven local rational FRFs; no parametric fit is read.

Frozen protocol: numerator/order 2, transient/order 2, denominator/order 2;
linear or three SK reweighting updates, windows 21/25/31. Front-half
leave-one-excitation-line-out prediction selects ONE global configuration.
This diagnostic does not replace the pre-specified primary LPM analysis.
"""
import os
for _n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_n]='1'
import sys
sys.dont_write_bytecode=True
if hasattr(sys.stdout,'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
import argparse,json,time
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
DATA=HERE.parent/'real_d'/'data'
OUT=HERE/'lrm_results'
ZS=('-0.5','0','0.5','1')
CONFIGS=[dict(method=m,window=w,sk_updates=0 if m=='linear' else 3)
         for m in ('linear','SK3') for w in (21,25,31)]
NPAR=17


def atomic_json(path,data):
    path.parent.mkdir(exist_ok=True,parents=True)
    tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=True),encoding='utf-8')
    tmp.replace(path)


def load(z):
    with np.load(DATA/f'z{z}_P08_s01.npz') as a:
        f=a['f'][a['line_indices']]
        D,U,Y=(a[k].copy() for k in ('Dper','Uper','Yper'))
    return f,D,U,Y


def output_data(D,U,Y):
    O=np.concatenate((U-D,Y),axis=-1)
    # Recorded Q has exactly zero fourth row. Check before enforcing identity.
    zero=np.max(np.abs(O),axis=(0,1)) <= 1e-13*max(np.max(np.abs(D)),1e-300)
    return O,zero


def fit_one(f,D,y,k,window,sk_updates=0,heldout=False):
    nf=len(f);start=min(max(k-window//2,0),nf-window)
    idx=np.arange(start,start+window)
    scale_f=max(np.max(np.abs(f[idx]-f[k])),1e-100)
    if heldout:idx=idx[idx!=k]
    x=(f[idx]-f[k])/scale_f
    yy=y[idx]
    X=np.column_stack([D[idx]*x[:,None]**p for p in range(3)]+
                      [x[:,None]**p for p in range(3)]+[-yy*x,-yy*x*x])
    weights=np.ones(len(idx));coef=None;cond=np.nan;rank=0;denmin=np.nan
    error=None
    for update in range(sk_updates+1):
        A=X*weights[:,None];rhs=yy*weights
        sc=np.linalg.norm(A,axis=0)
        if np.any(sc<=np.finfo(float).tiny):
            error='zero_column';break
        As=A/sc
        s=np.linalg.svd(As,compute_uv=False)
        cond=float(s[0]/s[-1]);rank=int(np.count_nonzero(s>s[0]*1e-12))
        if rank<NPAR:
            error='rank_deficient';break
        coef=np.linalg.lstsq(As,rhs,rcond=1e-12)[0]/sc
        den=1+coef[15]*x+coef[16]*x*x
        denmin=float(np.min(np.abs(den)))
        if not np.all(np.isfinite(coef)) or not np.all(np.isfinite(den)):
            error='nonfinite';break
        if update<sk_updates:
            # Purely numerical floor; no fit-quality based denominator rejection.
            weights=1/np.maximum(np.abs(den),np.finfo(float).eps)
    if error is not None:
        coef=np.full(NPAR,np.nan+1j*np.nan)
    prediction=D[k]@coef[:4]+coef[12]
    return coef,dict(condition=cond,rank=rank,denominator_min=denmin,
                     failure=error is not None,error=error,prediction=prediction,
                     first_index=int(idx[0]),last_index=int(idx[-1]))


def front_cv(z):
    f,Dp,Up,Yp=load(z);Dp,Up,Yp=Dp[:4],Up[:4],Yp[:4]
    Op,zero=output_data(Dp,Up,Yp)
    D=Dp.mean(0);O=Op.mean(0)
    # Same fixed weights for all six candidates. No model-dependent weights.
    noise=np.sqrt(np.mean(np.abs(Op-Op.mean(0))**2,axis=(0,1))*4/3)/np.sqrt(4)
    amplitude=np.sqrt(np.mean(np.abs(O)**2,axis=0))
    noise=np.maximum(noise,amplitude*1e-8)
    noise[zero]=1.
    centers=np.arange(15,len(f)-15)
    rows=[]
    for cfg in CONFIGS:
        err=np.full((len(centers),8),np.nan+1j*np.nan)
        cond=[];bad=[]
        for i,k in enumerate(centers):
            for o in np.flatnonzero(~zero):
                _,info=fit_one(f,D,O[:,o],k,cfg['window'],cfg['sk_updates'],True)
                err[i,o]=(info['prediction']-O[k,o])/noise[o]
                cond.append(info['condition'])
                if info['failure']:bad.append(dict(center=int(k),output=int(o),error=info['error']))
        vals=err[:,~zero]
        score=float(np.mean(np.abs(vals)**2)) if np.all(np.isfinite(vals)) else float('inf')
        rows.append(dict(**cfg,score_mse=score,n_terms=int(vals.size),sse=score*vals.size,
                         failures=bad,condition_median=float(np.nanmedian(cond)),
                         condition_max=float(np.nanmax(cond)),
                         channel_mse=[None if zero[o] else float(np.mean(np.abs(err[:,o])**2)) for o in range(8)]))
    result=dict(z=z,split='front',period_indices_zero_based=[0,1,2,3],
                centers=centers.tolist(),center_frequency_hz=f[centers].tolist(),
                zero_channels=np.flatnonzero(zero).tolist(),noise_mean_standard_error=noise.tolist(),
                weights='fixed per-output global period-difference standard error, with 1e-8 amplitude floor',
                rows=rows)
    atomic_json(OUT/f'z{z}_front_cv.json',result)
    return result


def local_mean(f,Dp,Up,Yp,cfg):
    Op,zero=output_data(Dp,Up,Yp)
    D=Dp.mean(0);O=Op.mean(0);nf=len(f)
    coef=np.zeros((nf,8,NPAR),complex)
    condition=np.full((nf,8),np.nan);rank=np.zeros((nf,8),int)
    denmin=np.full((nf,8),np.nan);failed=np.zeros((nf,8),bool)
    prediction=np.zeros((nf,8),complex);bounds=np.zeros((nf,2),int)
    for k in range(nf):
        for o in np.flatnonzero(~zero):
            c,i=fit_one(f,D,O[:,o],k,cfg['window'],cfg['sk_updates'])
            coef[k,o]=c;condition[k,o]=i['condition'];rank[k,o]=i['rank']
            denmin[k,o]=i['denominator_min'];failed[k,o]=i['failure']
            prediction[k,o]=i['prediction'];bounds[k]=[i['first_index'],i['last_index']]
    TU=coef[:,:4,:4]+np.eye(4)[None]
    TY=coef[:,4:,:4]
    G=np.full((nf,4,4),np.nan+1j*np.nan);tucond=np.full(nf,np.nan)
    for k in range(nf):
        if failed[k].any():continue
        try:
            tucond[k]=np.linalg.cond(TU[k])
            G[k]=np.linalg.solve(TU[k].T,TY[k].T).T
        except np.linalg.LinAlgError:
            failed[k,:]=True
    edge=(np.arange(nf)<cfg['window']//2)|(np.arange(nf)>=nf-cfg['window']//2)
    return dict(f=f,G=G,TU=TU,TY=TY,coefficients=coef,
                design_condition=condition,design_rank=rank,
                denominator_min=denmin,failed_output=failed,tud_condition=tucond,
                edge=edge,bounds=bounds,prediction=prediction,
                zero_channels=zero,nperiod=np.array(len(Dp)),window=np.array(cfg['window']),
                order=np.array(2),n_coefficients=np.array(NPAR),sk_updates=np.array(cfg['sk_updates']))


def estimate(z,split,cfg):
    f,D,U,Y=load(z)
    sl={'full':slice(None),'front':slice(0,4),'back':slice(4,8)}[split]
    D,U,Y=D[sl],U[sl],Y[sl]
    result=local_mean(f,D,U,Y,cfg);P=len(D);jack=[]
    jfailure=[]
    for j in range(P):
        ix=np.arange(P)!=j
        jj=local_mean(f,D[ix],U[ix],Y[ix],cfg)
        jack.append(jj['G']);jfailure.append(jj['failed_output'])
    jack=np.asarray(jack);jm=jack.mean(0)
    result['se_complex']=np.sqrt((P-1)/P*np.sum(np.abs(jack-jm)**2,axis=0))
    db=20*np.log10(np.maximum(np.abs(jack),1e-300))
    phase=np.angle(jack/np.where(np.abs(result['G'])>1e-300,result['G'],np.nan))*180/np.pi
    result['se_db']=np.sqrt((P-1)/P*np.sum((db-db.mean(0))**2,axis=0))
    result['se_phase_deg']=np.sqrt((P-1)/P*np.sum((phase-phase.mean(0))**2,axis=0))
    result['jackknife_G']=jack;result['jackknife_failed_output']=np.asarray(jfailure)
    result['jackknife_valid']=np.all(np.isfinite(jack),axis=(0,2,3))
    np.savez_compressed(OUT/f'z{z}_{split}_lrm.npz',**result)
    finite=np.isfinite(result['G']).all(axis=(1,2))
    summary=dict(z=z,split=split,configuration=cfg,frequency_count=len(f),
                 failed_frequency_count=int(np.count_nonzero(~finite)),
                 failed_output_count=int(result['failed_output'].sum()),
                 jackknife_failed_outputs=int(np.sum(jfailure)),
                 median_Tud_condition=float(np.nanmedian(result['tud_condition'])),
                 max_Tud_condition=float(np.nanmax(result['tud_condition'])),
                 median_se_db=float(np.nanmedian(result['se_db'])),
                 max_se_db=float(np.nanmax(result['se_db'])),
                 median_se_phase_deg=float(np.nanmedian(result['se_phase_deg'])))
    atomic_json(OUT/f'z{z}_{split}_summary.json',summary)
    return summary


def select(allcv):
    scores=[]
    for cfg in CONFIGS:
        rr=[next(r for r in c['rows'] if r['method']==cfg['method'] and r['window']==cfg['window']) for c in allcv]
        sse=sum(r['sse'] for r in rr);n=sum(r['n_terms'] for r in rr)
        scores.append(dict(**cfg,score_mse=sse/n,sse=sse,n_terms=n,
                           failed_fits=sum(len(r['failures']) for r in rr)))
    best=min(scores,key=lambda r:r['score_mse'])
    chosen={k:best[k] for k in ('method','window','sk_updates')}
    atomic_json(OUT/'configuration.json',dict(
        role='diagnostic; does not replace pre-specified LPM primary',
        selection='minimum aggregate front-half leave-one-excitation-line-out standardized prediction MSE',
        positions=[c['z'] for c in allcv],selected=chosen,all_scores=scores,
        heldout_line='removed from all periods and the local rational fit',
        period_jackknife='random repeatability only; excludes smoothing/systematic/model bias',
        parameter_model_used=False))
    return chosen


def verify_local_rational():
    """Algebra check on independently generated local-rational transfer data."""
    rng=np.random.default_rng(418217)
    f=np.linspace(3,25,61);x=(f-14)/11
    D=rng.normal(size=(len(f),4))+1j*rng.normal(size=(len(f),4))
    H=np.zeros((len(f),8,4),complex);O=np.zeros((len(f),8),complex)
    for o in (0,1,2,4,5,6,7):
        a=.2*(rng.normal(size=2)+1j*rng.normal(size=2))
        den=1+a[0]*x+a[1]*x*x
        B=(rng.normal(size=(3,4))+1j*rng.normal(size=(3,4)))*(.02 if o<4 else .001)
        t=(rng.normal(size=3)+1j*rng.normal(size=3))*.01
        H[:,o]=(B[0]+x[:,None]*B[1]+x[:,None]**2*B[2])/den[:,None]
        O[:,o]=np.einsum('ki,ki->k',H[:,o],D)+(t[0]+t[1]*x+t[2]*x*x)/den
    U=D+O[:,:4];Y=O[:,4:]
    G=np.array([np.linalg.solve((np.eye(4)+H[k,:4]).T,H[k,4:].T).T for k in range(len(f))])
    rows=[]
    for updates in (0,3):
        r=local_mean(f,D[None],U[None],Y[None],dict(window=25,sk_updates=updates))
        ge=float(np.max(np.abs(r['G']-G)))
        pe=[]
        for k in (15,30,45):
            for o in (0,1,2,4,5,6,7):
                _,info=fit_one(f,D,O[:,o],k,25,updates,True)
                pe.append(abs(info['prediction']-O[k,o]))
        rows.append(dict(sk_updates=updates,frf_max_absolute_error=ge,
                         heldout_prediction_max_absolute_error=float(max(pe)),
                         failed_outputs=int(r['failed_output'].sum())))
    assert all(r['frf_max_absolute_error']<1e-10 and r['heldout_prediction_max_absolute_error']<1e-10 and not r['failed_outputs'] for r in rows),rows
    atomic_json(OUT/'algebra_verification.json',dict(passed=True,seed=418217,rows=rows))


def diagnostic_bands():
    """Read-only comparison of already saved endpoints; no fitting or selection."""
    bands=[('3-8',3.,8.),('8-18',8.,18.),('18-25',18.,25.)]
    output=dict(role='post-fit descriptive diagnostic; no frequency exclusions or re-selection',
                band_convention='left inclusive, right exclusive; 25 Hz included in last band',
                norm='Euclidean norm over all complex entries; uniform frequency bins',
                G_columns='physical input directions',TU_TY_columns='external reference injection directions',
                repeatability=[],model_moments=[])
    def make_rows(z,kind,pair,delta,baseline,f,totaldelta=None,extra=None):
        rows=[];all_den=float(np.sum(np.abs(baseline)**2))
        all_delta=float(np.sum(np.abs(delta)**2)) if totaldelta is None else totaldelta
        for name,lo,hi in [('all',3.,25.)]+bands:
            mask=(f>=lo)&((f<=hi) if hi==25 else (f<hi))
            for column in [None,0,1,2,3]:
                de=delta[mask] if column is None else delta[mask,:,column]
                ba=baseline[mask] if column is None else baseline[mask,:,column]
                dn=float(np.sum(np.abs(de)**2));bn=float(np.sum(np.abs(ba)**2))
                rows.append(dict(z=z,quantity=kind,pair=pair,band=name,n_frequencies=int(mask.sum()),
                    input_column=None if column is None else column+1,
                    complex_relative_percent=100*np.sqrt(dn/bn) if bn>0 else None,
                    difference_energy_fraction_percent=100*dn/all_delta if all_delta>0 else None,
                    baseline_energy_fraction_percent=100*bn/all_den if all_den>0 else None,
                    delta_squared_norm=dn,baseline_squared_norm=bn,
                    finite=bool(np.isfinite(de).all() and np.isfinite(ba).all()),**(extra or {})))
        return rows
    for z in ZS:
        refs={s:dict(np.load(OUT/f'z{z}_{s}_lrm.npz')) for s in ('full','front','back')}
        f=refs['full']['f']
        for kind in ('G','TU','TY','TF'):
            vals={s:(r['TU']-np.eye(4)[None] if kind=='TF' else r[kind]) for s,r in refs.items()}
            for first,second in [('front','back'),('front','full'),('back','full')]:
                output['repeatability']+=make_rows(z,kind,first+'-'+second,vals[first]-vals[second],vals['full'],f)
        for fit,ref in [('full','full'),('front','front'),('back','back'),('front','back'),('back','front')]:
            mpath=HERE/'results'/f'z{z}_{fit}_M3_model.npz'
            if not mpath.exists():continue
            model=dict(np.load(mpath));assert np.array_equal(model['f'],f)
            p,s=(8,1) if fit=='full' else (4,1) if fit=='front' else (4,5)
            meta=json.loads((HERE.parent/'real_d/raw'/f'z{z}_P{p:02d}_s{s:02d}_M3.json').read_text())
            rr=refs[ref]['TY']-model['G']@refs[ref]['TU']
            output['model_moments']+=make_rows(z,'TY-Gmodel*TU',fit+'->'+ref,rr,refs[ref]['TY'],f,
                extra=dict(model_failure=bool(meta['failure']),model_iterations=meta['iterations'],
                           model_score=meta['score_normalized']))
    atomic_json(OUT/'diagnostic_bands.json',output)
    def select(rows,z,kind,band='all',column=None,pair='front-back'):
        return next(x for x in rows if x['z']==z and x['quantity']==kind and x['band']==band and
                    x['input_column']==column and x['pair']==pair)
    lines=['# LRM 重复性分解与无逆闭环矩条件诊断','',
           '只读取已保存LRM和M3端点，没有重新拟合或按结果挑选频带。全部133条频线均保留。分带为[3,8)、[8,18)、[18,25] Hz。',
           '重复性定义：100×‖前半−后半‖/‖完整记录‖，每一行使用本行频带及输入列的完整记录范数。百分比大也可能是分母响应弱，因此同时保留差异能量占比。', '',
           '|z/mm|G差异/%|TU差异/%|TY差异/%|反馈TF=TU−I差异/%|',
           '|---:|---:|---:|---:|---:|']
    for z in ZS:
        vals=[select(output['repeatability'],z,k)['complex_relative_percent'] for k in ('G','TU','TY','TF')]
        lines.append('|'+z+'|'+'|'.join(f'{v:.3f}' for v in vals)+'|')
    lines+=['','## 按频带分解：全部4输入列','',
            '|z|频带/Hz|G差异/%|占G总差异能量/%|TU差异/%|TY差异/%|',
            '|---:|---|---:|---:|---:|---:|']
    for z in ZS:
        for band,_,_ in bands:
            g=select(output['repeatability'],z,'G',band)
            tu=select(output['repeatability'],z,'TU',band)['complex_relative_percent']
            ty=select(output['repeatability'],z,'TY',band)['complex_relative_percent']
            lines.append(f"|{z}|{band}|{g['complex_relative_percent']:.3f}|{g['difference_energy_fraction_percent']:.3f}|{tu:.3f}|{ty:.3f}|")
    lines+=['','## 按输入方向分解：完整3–25 Hz','',
            '|z|物理输入列|G差异/%|占G总差异能量/%|完整G响应能量占比/%|',
            '|---:|---:|---:|---:|---:|']
    for z in ZS:
        for col in range(1,5):
            r=select(output['repeatability'],z,'G',column=col)
            lines.append(f"|{z}|{col}|{r['complex_relative_percent']:.3f}|{r['difference_energy_fraction_percent']:.3f}|{r['baseline_energy_fraction_percent']:.3f}|")
    lines+=['','## M3的无逆矩条件残差','',
            '量为100×‖TY−G_M3 TU‖/‖TY‖；不用TU逆，不把LRM对象FRF当真值。比较不同数据段只替换稳定FRF，不转移瞬态系数。', '',
            '|z|模型段→LRM段|全部频带/%|3–8/%|8–18/%|18–25/%|模型是否原失败端点|',
            '|---:|---|---:|---:|---:|---:|---|']
    for z in ZS:
        for pair in ('full->full','front->back','back->front'):
            rr=[select(output['model_moments'],z,'TY-Gmodel*TU',band=b,pair=pair) for b in ('all','3-8','8-18','18-25')]
            lines.append('|'+z+'|'+pair+'|'+'|'.join(f"{r['complex_relative_percent']:.3f}" for r in rr)+'|'+('是' if rr[0]['model_failure'] else '否')+'|')
    lines+=['','JSON还保存每一频带×每一输入列的分解、front−full/back−full差异、同半段和交叉半段矩残差，以及全部失败旗标。TF用于避免TU中的精确单位直通掩盖反馈部分的变化。',
            '', '## 诊断结论','']
    for z in ZS:
        bands_r=[select(output['repeatability'],z,'G',b) for b,_,_ in bands]
        cols_r=[select(output['repeatability'],z,'G',column=c) for c in range(1,5)]
        bb=max(bands_r,key=lambda r:r['difference_energy_fraction_percent'])
        cc=max(cols_r,key=lambda r:r['difference_energy_fraction_percent'])
        lines.append(f"- z={z}：G的前后差异能量最多位于{bb['band']} Hz（{bb['difference_energy_fraction_percent']:.1f}%），物理输入列{cc['input_column']}贡献最大（{cc['difference_energy_fraction_percent']:.1f}%）。这些是完整数据的事后分解，未据此删除点或重选模型。")
    lines+=['','TU/TY的变化、TU求逆的放大，以及真实对象可能的时间变化/未建模部分都可能影响G重复性；本分解不能单独将问题归因于某一种。即使无逆残差较小，也只说明闭环传递矩的一致程度，不能证明参数模型的真实精度。']
    (OUT/'diagnostic_bands.md').write_text('\n'.join(lines),encoding='utf-8')
    print('\n'.join(lines),flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--all',action='store_true');ap.add_argument('--diagnostics',action='store_true')
    args=ap.parse_args();OUT.mkdir(exist_ok=True)
    if args.diagnostics:
        diagnostic_bands();return
    tic=time.time()
    verify_local_rational()
    if args.pilot:
        cv=front_cv('0')
        print(json.dumps(cv['rows'],ensure_ascii=False,indent=2),flush=True)
        cfg={k:min(cv['rows'],key=lambda r:r['score_mse'])[k] for k in ('method','window','sk_updates')}
        atomic_json(OUT/'pilot_configuration.json',dict(selected_z0_only=cfg,provisional=True))
        print(json.dumps(estimate('0','full',cfg),indent=2),flush=True)
    else:
        cvs=[]
        for z in ZS:
            cv=front_cv(z);cvs.append(cv)
            print('CV',z,[(r['method'],r['window'],r['score_mse']) for r in cv['rows']],flush=True)
        cfg=select(cvs);print('GLOBAL_FROZEN',cfg,flush=True)
        for z in ZS:
            for split in ('full','front','back'):
                print(json.dumps(estimate(z,split,cfg),ensure_ascii=False),flush=True)
    print('elapsed_s',time.time()-tic,flush=True)


if __name__=='__main__':main()
