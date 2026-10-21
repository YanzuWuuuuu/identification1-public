"""Frozen-covariance full-output-matrix IV and nonexcited-bin statistics.

This module copies the algebra/stopping rule of prior A, without mutating it.
All covariances returned by local_covariances use the supplied full-record DFT
units. No hidden spectral normalization or diagonal loading is applied.
"""
import numpy as np
import time

def precision(C):
    C=(C+C.conj().transpose(0,2,1))*.5
    L=np.linalg.cholesky(C)
    Li=np.linalg.inv(L)
    return Li.conj().transpose(0,2,1)@Li

def whitener(W,nf=None,ny=None):
    if W is None:
        return np.broadcast_to(np.eye(ny),(nf,ny,ny))
    if W.ndim==2:
        return np.eye(W.shape[1])[None]*np.sqrt(W)[:,:,None]
    # np.linalg.cholesky W=R R^H, hence L=R^H and W=L^H L.
    return np.linalg.cholesky((W+W.conj().transpose(0,2,1))*.5).conj().transpose(0,2,1)

def iv_matrices(model,U,Y,Uh,beta,L,free=None):
    ub=Uh(beta) if callable(Uh) else Uh
    yh,J=model(beta,U)
    Z=J if ub is U else model(beta,ub)[1]
    if free is not None: J=J[:,:,free]; Z=Z[:,:,free]
    jw=np.einsum('koi,kip->kop',L,J)
    zw=np.einsum('koi,kip->kop',L,Z)
    ew=np.einsum('koi,ki->ko',L,Y-yh)
    j=jw.reshape(-1,jw.shape[-1]); z=zw.reshape(-1,zw.shape[-1]); e=ew.ravel()
    A=np.real(z.conj().T@j); g=np.real(z.conj().T@e)
    return A,g,z,e

def iv_score(model,U,Y,Uh,W,beta):
    L=whitener(W,len(U),U.shape[1]); A,g,z,e=iv_matrices(model,U,Y,Uh,beta,L)
    scale=np.linalg.norm(z,axis=0)*np.linalg.norm(e)
    return dict(score_max=float(np.max(np.abs(g)/np.maximum(scale,1e-300))),
                weighted_residual_norm=float(np.linalg.norm(e)))

def iv_solve(model,U,Y,Uh,W,beta0,iw,itr,maxiter=80,tol=1e-5,cap=.05,free=None):
    t=time.perf_counter(); b=beta0.copy(); L=whitener(W,len(U),U.shape[1])
    fr=np.arange(len(b)) if free is None else np.asarray(free)
    info=dict(status='exception',success=False,iterations=0)
    try:
        for it in range(maxiter):
            A,g,_,_=iv_matrices(model,U,Y,Uh,b,L,fr)
            sc=1/np.sqrt(np.abs(np.diag(A))+1e-300)
            step=sc*np.linalg.solve(A*sc[:,None]*sc[None,:],g*sc)
            if free is not None:
                b[fr]+=step; info.update(status='linear_step',success=True,iterations=1); break
            mx=np.max(np.abs(step[iw:itr]))
            if mx>cap: step*=cap/mx
            b+=step
            info['iterations']=it+1
            if np.max(np.abs(step[iw:itr]))<tol: break
        if free is None:
            finite=bool(np.all(np.isfinite(b)))
            ok=finite and info['iterations']<maxiter
            info.update(success=ok,status='pole_step_converged' if ok else ('iteration_limit' if finite else 'nonfinite'))
        if np.all(np.isfinite(b)): info.update(iv_score(model,U,Y,Uh,W,b))
    except Exception as exc:
        info.update(status='exception',success=False,exception=repr(exc))
        # Keep last finite endpoint for audit, but never count it as converged.
    info['seconds']=time.perf_counter()-t
    return b,info

def shrink_cov(C,alpha):
    C=(C+C.conj().transpose(0,2,1))*.5
    return (1-alpha)*C+alpha*np.eye(C.shape[-1])[None]*np.real(np.diagonal(C,axis1=1,axis2=2))[:,None,:]

def local_covariances(k,U,Y,beta,model_full,nper,nl=16,alpha=.1):
    if nper<2: raise ValueError('P=1 has no nonexcited bins')
    k=np.asarray(k); lines=np.flatnonzero(k%nper==0); noise=np.flatnonzero(k%nper!=0)
    if len(noise)<nl: raise ValueError('Insufficient nonexcited bins')
    # The nearest nl TOTAL bins, not nl on each side. Ties lower-frequency first.
    ix=np.array([noise[np.argsort(np.abs(k[noise]-k[q]),kind='stable')[:nl]] for q in lines])
    yh,_=model_full(beta,U); trans,_=model_full(beta,np.zeros_like(U))
    residual=Y-yh; X=np.concatenate([Y-trans,U],axis=1)
    rv=residual[ix]; xx=X[ix]
    Cv=np.einsum('kni,knj->kij',rv,rv.conj())/nl
    Cx=np.einsum('kni,knj->kij',xx,xx.conj())/nl
    Cvs=shrink_cov(Cv,alpha); Cxs=shrink_cov(Cx,alpha)
    ny=Y.shape[1]
    return dict(covariance_residual=Cvs,CY=Cxs[:,:ny,:ny],CU=Cxs[:,ny:,ny:],CYU=Cxs[:,:ny,ny:],
                covariance_residual_raw=Cv,covariance_joint_raw=Cx,indices=ix,
                line_indices=lines,distances=np.abs(k[ix]-k[lines,None]),nl=nl,alpha=alpha)

def statistics(records):
    out={}
    for name,selected in [('all_finite',[r for r in records if r.get('error_percent') is not None and np.all(np.isfinite(r['error_percent']))]),
                          ('converged_only',[r for r in records if r.get('success') and r.get('error_percent') is not None and np.all(np.isfinite(r['error_percent']))])]:
        a=np.array([r['error_percent'] for r in selected]); n=len(a)
        if n:
            sd=np.std(a,axis=0,ddof=1) if n>1 else np.full(a.shape[1],np.nan)
            out[name]=dict(n=n,mean=a.mean(0).tolist(),se=(sd/np.sqrt(n)).tolist(),sd=sd.tolist(),rmse=np.sqrt(np.mean(a*a,axis=0)).tolist())
        else: out[name]=dict(n=0)
    it=[r['iterations'] for r in records if r.get('iterations') is not None]
    out.update(n_total=len(records),failures=sum(not r.get('success',False) for r in records),mean_iterations=float(np.mean(it)) if it else None)
    return out
