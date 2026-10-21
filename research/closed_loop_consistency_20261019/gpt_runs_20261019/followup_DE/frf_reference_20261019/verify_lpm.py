"""Read-only independent algebra / saved-result checks; no estimator re-fitting."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parent/'real_d'))
sys.path.insert(0,str(HERE.parent))
from frf_core import local_reference_frf
import realfit_extended as model
import loader_copy as loader
from sml_core_copy import gmodel_from_output, SMLResidual

def rel(a,b): return float(np.linalg.norm(a-b)/np.linalg.norm(b))

def toy():
    rng=np.random.default_rng(1933)
    nf,P=133,4
    f=np.linspace(3,25,nf); x=(f-14)/11
    D=np.exp(2j*np.pi*rng.random((nf,4)))
    A=rng.normal(size=(3,4,4))+1j*rng.normal(size=(3,4,4))
    B=rng.normal(size=(3,4,4))+1j*rng.normal(size=(3,4,4))
    v=np.array([x**r for r in range(3)])
    TU=np.eye(4)[None]+.03*np.einsum('ri,roj->ioj',v,A)
    TY=np.einsum('ri,roj->ioj',v,B)
    U=np.repeat(np.einsum('koi,ki->ko',TU,D)[None],P,axis=0)
    Y=np.repeat(np.einsum('koi,ki->ko',TY,D)[None],P,axis=0)
    for pp in range(P):
        U[pp]+=(pp+1)*.005*(1+x[:,None]+2j*x[:,None]**2)
        Y[pp]+=(pp+1)*.01*(1+2*x[:,None]+1j*x[:,None]**2)
    r=local_reference_frf(f,np.repeat(D[None],P,axis=0),U,Y)
    G=np.linalg.solve(TU.transpose(0,2,1),TY.transpose(0,2,1)).transpose(0,2,1)
    return dict(TU_relative_error=rel(r['TU'],TU),TY_relative_error=rel(r['TY'],TY),G_relative_error=rel(r['G'],G))

def saved_lpm():
    a=np.load(HERE/'results/z0_full_q2_w21.npz'); g=a['G']; edge=a['edge']
    result=dict(edge_count=int(edge.sum()),methods={})
    for met in ('M1','M2','M3','M3w','M5','M5n'):
        m=np.load(HERE/f'results/z0_full_{met}_model.npz')['G']; delta=m-g
        rr=a['TY']-m@a['TU']
        result['methods'][met]=dict(relative_G_percent=100*rel(m,g),interior_relative_G_percent=100*rel(m[~edge],g[~edge]),
            inverse_free_relative_percent=100*float(np.linalg.norm(rr)/np.linalg.norm(a['TY'])),
            error_energy_edge_fraction=float(np.sum(abs(delta[edge])**2)/np.sum(abs(delta)**2)))
    for name in ('TU','TY'):
        aa=np.load(HERE/'results/z0_front_q2_w21.npz')[name]; bb=np.load(HERE/'results/z0_back_q2_w21.npz')[name]
        result['halves_'+name+'_relative_percent']=100*float(np.linalg.norm(aa-bb)/np.linalg.norm((aa+bb)/2))
    return result

def sml_saved(z):
    d=np.load(HERE.parent/f'real_d/data/z{z}_P08_s01.npz')
    f=d['f'][d['line_indices']]
    c=model.Ctx(f,loader.TAU-model.Ts/2,'zoh')
    out=lambda b,u:model.out_jac(c,b,u)
    gm=gmodel_from_output(out,model.NB,len(f))
    j=json.loads((HERE.parent/f'real_d/raw/z{z}_P08_s01_M5.json').read_text(encoding='utf-8'))
    b=np.asarray(j['beta']); b0=np.asarray(j['initial_beta'])
    gg=gm(b)[0]; g=model.resp(c,b)
    um=d['Uper'].mean(0); ym=d['Yper'].mean(0)
    y1=out(b,um)[0]; scaled=b.copy(); scaled[model.ITR:]*=8
    y8=out(scaled,um*8)[0]
    obj=SMLResidual(d['Uper'],d['Yper'],out,gm,shrink=0.)
    val=obj.fun(b); ini=obj.fun(b0)
    def cov_stats(bb):
        gx=model.resp(c,bb); gh=gx.conj().transpose(0,2,1)
        ce=obj.CY-obj.CYU@gh-gx@obj.CYU.conj().transpose(0,2,1)+gx@obj.CU@gh
        return dict(median_trace=float(np.median(np.trace(ce,axis1=1,axis2=2).real)),
                    median_condition=float(np.median(np.linalg.cond(ce))))
    lrmp=HERE/f'lrm_results/z{z}_full_lrm.npz'
    result=dict(gmodel_resp_relative_error=rel(gg,g),period_to_full_prediction_scaling_relative_error=rel(y8,8*y1),
        model_FRF_unchanged_after_transient_scaling=rel(model.resp(c,scaled),g),
        SML_cost_initial=float(ini@ini),SML_cost_final=float(val@val),
        raw_output_initial_relative_percent=100*float(np.linalg.norm(ym-out(b0,um)[0])/np.linalg.norm(ym)),
        raw_output_relative_percent=100*float(np.linalg.norm(ym-y1)/np.linalg.norm(ym)),
        covariance_initial=cov_stats(b0),covariance_final=cov_stats(b),
        status=j.get('status'),converged=j.get('converged'),iterations=j.get('iterations'),frequency_hz=j.get('frequency_hz'),damping_percent=j.get('damping_percent'))
    if lrmp.exists():
        r=np.load(lrmp)
        result['LRM_keys']=r.files
        if 'G' in r.files:
            result['LRM_FRF_difference_percent']=100*rel(g,r['G'])
        if 'TU' in r.files and 'TY' in r.files:
            result['LRM_inverse_free_relative_percent']=100*float(np.linalg.norm(r['TY']-g@r['TU'])/np.linalg.norm(r['TY']))
    return result

if __name__=='__main__':
    ans=dict(noiseless_polynomial_test=toy(),saved_LPM=saved_lpm(),saved_M5={z:sml_saved(z) for z in ('0','1')})
    print(json.dumps(ans,indent=2,ensure_ascii=False))
    (HERE/'verify_lpm_result.json').write_text(json.dumps(ans,indent=2,ensure_ascii=False),encoding='utf-8')
