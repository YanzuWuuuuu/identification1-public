"""Periodic profile-SML with the derivative of its model-dependent covariance.

Uper/Yper: (period, frequency, channel), unnormalized period DFTs.
model(beta,U) -> (Yhat,J), gmodel(beta) -> (G,Gjac).
All parameters are real. Gjac dimensions: (frequency,out,in,parameter).
For P<=4, the unregularized four-output residual covariance is singular.
The optional joint diagonal shrinkage is an explicitly DIFFERENT estimator.
"""
import time
import numpy as np
from scipy.optimize import least_squares

def sample_covariances(Uper,Yper,shrink=0.):
    p=len(Uper)
    if p<2: raise ValueError('At least two periods are required')
    X=np.concatenate([Yper,Uper],axis=2)
    xc=X-X.mean(axis=0)
    C=np.einsum('pki,pkj->kij',xc,xc.conj())/(p*(p-1))
    if shrink:
        diag=np.real(np.diagonal(C,axis1=1,axis2=2))
        C=(1-shrink)*C+shrink*np.eye(C.shape[1])[None]*diag[:,None,:]
    ny=Yper.shape[2]
    return C[:,:ny,:ny],C[:,ny:,ny:],C[:,:ny,ny:]

class SMLResidual:
    def __init__(self,Uper,Yper,model,gmodel,shrink=0.):
        self.U=np.mean(Uper,axis=0); self.Y=np.mean(Yper,axis=0)
        self.CY,self.CU,self.CYU=sample_covariances(Uper,Yper,shrink)
        self.model=model; self.gmodel=gmodel; self.cache=None
    def calc(self,beta):
        if self.cache is not None and np.array_equal(beta,self.cache[0]): return self.cache[1:]
        yh,J=self.model(beta,self.U); G,Gj=self.gmodel(beta)
        GH=G.conj().transpose(0,2,1)
        Ce=self.CY-self.CYU@GH-G@self.CYU.conj().transpose(0,2,1)+G@self.CU@GH
        Ce=(Ce+Ce.conj().transpose(0,2,1))/2
        L=np.linalg.cholesky(Ce); Li=np.linalg.inv(L)
        r=np.einsum('kij,kj->ki',Li,self.Y-yh)
        # dCe=dG(CU G^H-CYU^H)+its Hermitian transpose.
        X=self.CU@GH-self.CYU.conj().transpose(0,2,1)
        dC=np.einsum('kipm,kpj->kijm',Gj,X)
        dC=dC+dC.conj().transpose(0,2,1,3)
        A=np.einsum('kai,kijm,kbj->kabm',Li,dC,Li.conj(),optimize=True)
        tri=np.tril(np.ones((L.shape[1],L.shape[1])))-.5*np.eye(L.shape[1])
        A*=tri[None,:,:,None]
        jr=-np.einsum('kij,kjm->kim',Li,J)-np.einsum('kijm,kj->kim',A,r)
        rr=np.r_[r.real.ravel(),r.imag.ravel()]
        jj=np.concatenate([jr.real.reshape(-1,len(beta)),jr.imag.reshape(-1,len(beta))],axis=0)
        self.cache=(beta.copy(),rr,jj)
        return rr,jj
    def fun(self,beta): return self.calc(beta)[0]
    def jac(self,beta): return self.calc(beta)[1]

class SMLCovResidual(SMLResidual):
    """Same profile objective with externally supplied covariance in U/Y units."""
    def __init__(self,U,Y,CY,CU,CYU,model,gmodel):
        self.U=U; self.Y=Y; self.CY=CY; self.CU=CU; self.CYU=CYU
        self.model=model; self.gmodel=gmodel; self.cache=None

def solve_sml_cov(U,Y,CY,CU,CYU,beta0,model,gmodel,maxiter=120,free=None):
    t0=time.perf_counter(); obj=SMLCovResidual(U,Y,CY,CU,CYU,model,gmodel)
    fr=np.arange(len(beta0)) if free is None else np.asarray(free)
    def expand(x):
        b=beta0.copy(); b[fr]=x; return b
    try:
        sol=least_squares(lambda x:obj.fun(expand(x)),beta0[fr],jac=lambda x:obj.jac(expand(x))[:,fr],
                          x_scale='jac',method='trf',ftol=1e-9,xtol=1e-9,gtol=1e-7,max_nfev=maxiter)
        b=expand(sol.x); rr,jj=obj.calc(b); g=jj.T@rr
        score=float(np.max(np.abs(g[fr])/(np.linalg.norm(jj[:,fr],axis=0)*max(np.linalg.norm(rr),1e-30)+1e-30)))
        info=dict(success=bool(sol.success),status=str(sol.message),nfev=int(sol.nfev),njev=int(sol.njev or 0),
                  cost=float(2*sol.cost),optimality=float(sol.optimality),normalized_score=score,seconds=time.perf_counter()-t0)
    except Exception as exc:
        b=np.full_like(beta0,np.nan); info=dict(success=False,status=repr(exc),nfev=0,njev=0,seconds=time.perf_counter()-t0)
    return b,info

def solve_sml(Uper,Yper,beta0,model,gmodel,maxiter=120,shrink=0.,free=None):
    """Returns beta, diagnostic dict. No silent pseudoinverse or covariance ridge."""
    t0=time.perf_counter()
    if len(Uper)<=Yper.shape[2] and shrink==0:
        return np.full_like(beta0,np.nan),dict(success=False,status='covariance_rank_deficient',nfev=0,njev=0,seconds=0.,shrink=shrink)
    obj=SMLResidual(Uper,Yper,model,gmodel,shrink)
    fr=np.arange(len(beta0)) if free is None else np.asarray(free)
    def expand(x):
        b=beta0.copy(); b[fr]=x; return b
    def fun(x): return obj.fun(expand(x))
    def jac(x): return obj.jac(expand(x))[:,fr]
    try:
        sol=least_squares(fun,beta0[fr],jac=jac,x_scale='jac',method='trf',
                          ftol=1e-9,xtol=1e-9,gtol=1e-7,max_nfev=maxiter)
        b=expand(sol.x); rr,jj=obj.calc(b); g=jj.T@rr
        score=float(np.max(np.abs(g[fr])/(np.linalg.norm(jj[:,fr],axis=0)*max(np.linalg.norm(rr),1e-30)+1e-30)))
        info=dict(success=bool(sol.success),status=str(sol.message),nfev=int(sol.nfev),njev=int(sol.njev or 0),cost=float(2*sol.cost),
                  optimality=float(sol.optimality),normalized_score=score,seconds=time.perf_counter()-t0,shrink=shrink)
    except Exception as e:
        b=np.full_like(beta0,np.nan);info=dict(success=False,status=repr(e),nfev=0,njev=0,seconds=time.perf_counter()-t0,shrink=shrink)
    return b,info

def gmodel_from_output(model,nb,nfreq,nu=4):
    """Fallback adapter: subtract U=0 response/derivative to remove transients."""
    def gg(beta):
        zero=np.zeros((nfreq,nu),complex); y0,j0=model(beta,zero)
        G=np.empty((nfreq,y0.shape[1],nu),complex); J=np.empty(G.shape+(nb,),complex)
        for q in range(nu):
            uu=zero.copy(); uu[:,q]=1.; yy,jj=model(beta,uu)
            G[:,:,q]=yy-y0; J[:,:,q,:]=jj-j0
        return G,J
    return gg
