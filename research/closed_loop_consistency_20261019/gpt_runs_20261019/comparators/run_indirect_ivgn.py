"""Canonical IV-Gauss-Newton FRF fit, adapted from cl2.iv_solve.

Removes transient parameters, keeps nonlinear pole derivatives, and supports the
full 16x16 FRF precision via Cholesky whitening. Exactly the original normal-
equation step, 0.05 pole-step cap and 1e-5 pole stopping rule; no LM globalization.
"""
import sys,time,json,pathlib,numpy as np
from run_indirect import c,HERE,atomic_npz
from sml_core import gmodel_from_output

def solve(data,pert,shrink=None,maxiter=80):
    zs=data['zk'];G=data['G'];bt=c.truth_beta();b=bt.copy()
    b[c.IW:c.IZ]+=np.log(1+.01*pert[:c.NM]);b[c.IZ:c.ITR]+=.4*pert[c.NM:];b[:c.NLIN]=0;b[c.ITR:]=0
    h0,h1,_=c.blocks(zs,b[c.IW:c.IZ],b[c.IZ:c.ITR]);H=np.stack([h for pair in zip(h0,h1) for h in pair],axis=1)
    coeff=np.linalg.lstsq(np.r_[H.real,H.imag],np.r_[G.real.reshape(len(zs),16),G.imag.reshape(len(zs),16)],rcond=None)[0]
    b[:c.NLIN]=coeff.ravel();b0=b.copy()
    Li=None
    if shrink is not None:
        C=data['cov_frf'];di=np.real(np.diagonal(C,axis1=1,axis2=2));C=(1-shrink)*C+shrink*np.eye(16)[None]*di[:,None,:]
        Li=np.linalg.inv(np.linalg.cholesky((C+C.conj().transpose(0,2,1))/2))
    model=lambda b,U:c.out_jac(zs,b,U);gm=gmodel_from_output(model,c.NB,len(zs));t=time.perf_counter();steps=[]
    def normal(beta):
        gh,jh=gm(beta);e=(G-gh).transpose(0,2,1).reshape(len(zs),16);J=jh[:,:,:,:c.ITR].transpose(0,2,1,3).reshape(len(zs),16,c.ITR)
        if Li is not None:e=np.einsum('kij,kj->ki',Li,e);J=np.einsum('kij,kjm->kim',Li,J)
        J=J.reshape(-1,c.ITR);e=e.ravel();A=np.real(J.conj().T@J);g=np.real(J.conj().T@e)
        return A,g,float(np.vdot(e,e).real)
    try:
        success=False
        for it in range(maxiter):
            A,g,cost=normal(b);sc=1/np.sqrt(np.abs(np.diag(A))+1e-300)
            step=sc*np.linalg.solve(A*sc[:,None]*sc[None,:],g*sc)
            mx=np.max(np.abs(step[c.IW:c.ITR]));
            if mx>.05:step*=.05/mx
            b[:c.ITR]+=step;steps.append(float(np.max(np.abs(step[c.IW:c.ITR]))))
            if steps[-1]<1e-5:success=True;break
        A,g,cost=normal(b);score=float(np.max(np.abs(g)/(np.sqrt(np.abs(np.diag(A)))*max(np.sqrt(cost),1e-30)+1e-30)))
        info=dict(success=success,status='pole_step_tolerance' if success else 'iteration_limit',iterations=it+1,
                  cost=cost,normalized_score=score,seconds=time.perf_counter()-t,shrink=shrink,
                  scaled_condition=float(np.linalg.cond(A/np.sqrt(np.outer(np.diag(A),np.diag(A))))))
    except Exception as e:
        info=dict(success=False,status=repr(e),iterations=len(steps),seconds=time.perf_counter()-t,shrink=shrink)
    return b,b0,info,np.array(steps)

def one(path):
    data=dict(np.load(path));P=int(path.stem[1:3]);r=int(path.stem.split('r')[1]);od=HERE/'indirect_ivgn_results';od.mkdir(exist_ok=True)
    a=np.load(HERE.parent/'experiment_a'/'data'/path.name);pert=a['pert']
    for method,shr in [('M6u_IVGN',None)]+([('M6_reg_IVGN',.1)] if P>4 else []):
        target=od/f'{path.stem}_{method}.npz'
        if target.exists():continue
        b,b0,info,steps=solve(data,pert,shr);bt=c.truth_beta();info.update(P=P,r=r,method=method,covariance_rank_upper_bound=P-4)
        atomic_npz(target,beta=b,beta0=b0,truth=bt,steps=steps,relative_error_percent=100*np.expm1(b[c.IW:c.ITR]-bt[c.IW:c.ITR]),info=json.dumps(info))
        print(target.name,info,flush=True)

if __name__=='__main__':
    worker=int(sys.argv[1]) if len(sys.argv)>1 else 0;nworkers=int(sys.argv[2]) if len(sys.argv)>2 else 1
    jobs=[(P,r) for P in [4,8,16] for r in range(20 if P==16 else 40) if r%nworkers==worker]
    while jobs:
        done=[]
        for P,r in jobs:
            path=HERE/'indirect_data'/f'P{P:02d}_r{r:03d}.npz'
            if path.exists():one(path);done.append((P,r))
        jobs=[j for j in jobs if j not in done]
        if jobs:time.sleep(10)
