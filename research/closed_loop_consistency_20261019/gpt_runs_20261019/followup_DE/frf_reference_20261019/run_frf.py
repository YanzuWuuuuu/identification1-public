import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np
from frf_core import local_reference_frf,metrics

HERE=Path(__file__).resolve().parent
REAL=HERE.parent/'real_d'
EXP_C=HERE.parents[1]/'experiment_c'
sys.path.insert(0,str(REAL))
import realfit_extended as m
import loader_copy as loader

CONFIGS=((2,21),(1,15),(2,31))
METHODS=('M1','M2','M3','M3w','M5','M5n')

def writejson(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=True),encoding='utf-8')

def fit_path(z,P,start,method):
    if method=='M1':
        return EXP_C/'raw'/f'z{z}_P{P}_s{start}_zoh_M1.json'
    return REAL/'raw'/f'z{z}_P{P:02d}_s{start:02d}_{method}.json'

def run(z):
    dest=HERE/'results';dest.mkdir(parents=True,exist_ok=True)
    # Reuse immutable spectra cache. Parameter model not used by nonparametric fit.
    cache=REAL/'data'/f'z{z}_P08_s01.npz'
    a=np.load(cache); f=a['f'][a['line_indices']]
    arrays=[a[x] for x in ('Dper','Uper','Yper')]
    periods={'full':slice(0,8),'front':slice(0,4),'back':slice(4,8)}
    refs={};summary=dict(z_mm=z,source=str(cache),sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
                        fmin=float(f[0]),fmax=float(f[-1]),n_lines=len(f),configs=[],fits=[],comparisons=[])
    summary['reference_repeat_relative']=float(np.linalg.norm(arrays[0]-arrays[0].mean(0))/np.linalg.norm(arrays[0]))
    for name,sl in periods.items():
        for order,win in CONFIGS:
            result=local_reference_frf(f,*[x[sl] for x in arrays],order,win)
            np.savez_compressed(dest/f'z{z}_{name}_q{order}_w{win}.npz',**result)
            refs[name,order,win]=result
            summary['configs'].append(dict(split=name,order=order,window=win,
               design_cond_median=float(np.median(result['design_condition'])),
               design_cond_max=float(np.max(result['design_condition'])),
               rank_min=int(np.min(result['design_rank'])),
               tud_cond_median=float(np.median(result['tud_condition'])),
               tud_cond_max=float(np.max(result['tud_condition'])),
               nonfinite_frf_bins=int(np.count_nonzero(~np.isfinite(result['G']).all(axis=(1,2)))),
               se_db_median=float(np.nanmedian(result['se_db'])),
               rms_random_relative_percent=float(100*np.linalg.norm(result['se_complex'])/np.linalg.norm(result['G']))))
    for name in periods:
        ref=refs[name,2,21]
        for q,w in CONFIGS[1:]:
            summary.setdefault('window_sensitivity',[]).append(dict(split=name,order=q,window=w,
                **metrics(refs[name,q,w]['G'],ref['G'])))
    summary['halves_difference']=metrics(refs['front',2,21]['G'],refs['back',2,21]['G'])
    summary['half_repeatability_relative_to_full_percent']=float(100*np.linalg.norm(refs['front',2,21]['G']-refs['back',2,21]['G'])/np.linalg.norm(refs['full',2,21]['G']))
    ctx=m.Ctx(f,loader.TAU-m.Ts/2,'zoh')
    for split,P,start in (('full',8,1),('front',4,1),('back',4,5)):
        for method in METHODS:
            path=fit_path(z,P,start,method)
            meta=dict(split=split,method=method,P=P,start=start,source=str(path))
            if not path.exists():
                meta.update(available=False,reason='prior_fit_missing');summary['fits'].append(meta);continue
            obj=json.loads(path.read_text(encoding='utf-8'))
            if 'beta' not in obj:
                meta.update(available=False,reason='no_parameter_endpoint',status=obj.get('status'));summary['fits'].append(meta);continue
            success=bool(obj.get('converged',obj.get('converged_pole_step',False)))
            meta.update(available=True,converged=success,iterations=obj.get('iterations'),
                        score=obj.get('score_normalized'),frequency_hz=obj.get('frequency_hz'),
                        damping_percent=obj.get('damping_percent'))
            summary['fits'].append(meta)
            g=m.resp(ctx,np.asarray(obj['beta']))
            np.savez_compressed(dest/f'z{z}_{split}_{method}_model.npz',f=f,G=g,beta=obj['beta'])
            targets=[split] if split=='full' else [split,'back' if split=='front' else 'front']
            for target in targets:
                for q,w in CONFIGS:
                    ref=refs[target,q,w]
                    mask=np.ones(len(f),bool)
                    summary['comparisons'].append(dict(method=method,fit_split=split,reference_split=target,
                        order=q,window=w,converged=success,comparison='same_segment' if split==target else 'opposite_half',
                        **metrics(g[mask],ref['G'][mask])))
    writejson(dest/f'z{z}_summary.json',summary)
    print(json.dumps({k:summary[k] for k in ('z_mm','n_lines','configs','half_repeatability_relative_to_full_percent')},ensure_ascii=True),flush=True)
    for row in summary['comparisons']:
        if row['fit_split']=='full' and row['window']==21:
            print(row['method'],row['converged'],row['complex_relative_percent'],row['magnitude_rmse_db'],row['phase_rmse_deg'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--z',default='0');args=parser.parse_args()
    for z in (('-0.5','0','0.5','1') if args.z=='all' else [args.z]):run(z)
