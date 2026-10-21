"""Read-only check of LM normal equations versus stable augmented least squares."""
import os
for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[n]='1'
os.environ.update(NOFLEX='1',FC='60',AL='30',AH='100',FHI='100',FLO='3',SW='.03',FSPR='.01',ZSPR='.4')
import sys
sys.dont_write_bytecode=True
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'code'))
import cl2 as c
from iv_lm import moments,exact_moment_jacobian

dat=np.load(HERE.parent/'experiment_a/data/P02_r000.npz')
k,z,U,Y,D=[dat[q] for q in ('k','zk','U','Y','D')]
line=k%2==0;zz=z[line];uu=U[line];yy=Y[line];dd=D[line];K=c.Ktf(zz)
model=lambda b,u:c.out_jac(zz,b,u)
inst=lambda b:np.linalg.solve(np.eye(4)[None]+K@c.resp(zz,b),dd[:,:,None])[:,:,0]
record=json.loads((HERE/'short_e/P02_r000_M3w_E3.json').read_text())
initial=np.array(record['initial_beta']);final=np.array(record['beta'])
residual=Y-c.out_jac(z,initial,U,need_jac=False)[0]
non=np.flatnonzero(~line);idx=np.flatnonzero(line)
nearest=np.array([non[np.argsort(np.abs(k[non]-k[q]),kind='stable')[:16]] for q in idx])
vv=residual[nearest];C=np.einsum('kli,klj->kij',vv,vv.conj())/16
C=.9*C+.1*np.eye(4)[None]*np.real(np.diagonal(C,axis1=1,axis2=2))[:,None,:]
W=np.linalg.inv(C)
out=[]
for name,b in [('initial',initial),('final',final)]:
    A,g,score=exact_moment_jacobian(b,uu,yy,model,inst,W,c.IW,c.ITR,K,lambda bb:c.resp(zz,bb))
    scale=1/np.sqrt(np.maximum(np.abs(np.diag(A)),1e-300))
    As=A*scale[:,None]*scale[None,:];gs=g*scale
    normal=As.T@As;rhs=As.T@gs;diag=np.maximum(np.diag(normal),1e-24)
    state=dict(point=name,score=score,condition_scaled_jacobian=float(np.linalg.cond(As)),
               condition_normal=float(np.linalg.cond(normal)),last_damping=record['history'][-1]['damping'],steps=[])
    for damp in (1e-3,1e-6,1e-9,1e-12,record['history'][-1]['damping']):
        N=normal+damp*np.diag(diag)
        dx_normal=np.linalg.solve(N,rhs)
        augmented=np.vstack([As,np.diag(np.sqrt(damp*diag))])
        target=np.r_[gs,np.zeros(len(gs))]
        dx_aug=np.linalg.lstsq(augmented,target,rcond=None)[0]
        step_n=scale*dx_normal;step_a=scale*dx_aug
        step=dict(damping=damp,condition_damped_normal=float(np.linalg.cond(N)),
                  relative_step_difference=float(np.linalg.norm(dx_normal-dx_aug)/np.linalg.norm(dx_aug)),
                  max_modal_normal_step=float(np.max(np.abs(step_n[c.IW:c.ITR]))),
                  max_modal_augmented_step=float(np.max(np.abs(step_a[c.IW:c.ITR]))),
                  max_modal_step_difference=float(np.max(np.abs((step_n-step_a)[c.IW:c.ITR]))),
                  predicted_normal_merit=float(np.linalg.norm(As@dx_normal-gs)**2),
                  predicted_augmented_merit=float(np.linalg.norm(As@dx_aug-gs)**2))
        for label,delta in [('normal',step_n),('augmented',step_a)]:
            try:
                _,gg,ss=moments(b+delta,uu,yy,model,inst,W)
                step[label+'_actual_merit']=float(np.linalg.norm(gg*scale)**2)
                step[label+'_score']=ss
            except Exception as exc:step[label+'_exception']=repr(exc)
        state['steps'].append(step)
    out.append(state)
(HERE/'audit_e_lm_condition.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
