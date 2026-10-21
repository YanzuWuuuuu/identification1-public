"""Experiment C extension. Original realfit.py is retained byte-for-byte nearby.
Only modal transfer functions change under zoh; D0,D2 remain the original
band-limited background. Fixed effective delays decrease by Ts/2 externally.
"""
import numpy as np
import scipy.linalg as sla
import realfit_original_copy as rf
from realfit_original_copy import *

class Ctx(rf.Ctx):
    def __init__(self, f, tau, sampling='ct'):
        super().__init__(f,tau)
        self.sampling=sampling
    def sub(self,sl):
        c=Ctx.__new__(Ctx)
        for a in ['f','s','zi','E']: setattr(c,a,getattr(self,a)[sl])
        c.sampling=self.sampling
        return c

def zoh_2nd(w,z,zk):
    A=np.array([[0.,1.],[-w*w,-2*z*w]])
    B=np.array([[0.],[1.]])
    E=sla.expm(np.block([[A,B],[np.zeros((1,3))]])*Ts)
    Ad,Bd=E[:2,:2],E[:2,2]
    a=np.array([1.,-np.trace(Ad),np.linalg.det(Ad)])
    den=zk*zk+a[1]*zk+a[2]
    h0=((zk-Ad[1,1])*Bd[0]+Ad[0,1]*Bd[1])/den
    h1=(Ad[1,0]*Bd[0]+(zk-Ad[0,0])*Bd[1])/den
    return w*w*h0,w*h1,a

def blocks(c,lw,lz):
    out=[]
    for i in range(NMODE):
        w,z=np.exp(lw[i]),np.exp(lz[i])
        if c.sampling=='zoh': h0,h1,ad=zoh_2nd(w,z,1/c.zi)
        else:
            den=c.s*c.s+2*z*w*c.s+w*w
            h0,h1,ad=w*w/den,w*c.s/den,ad_zoh(w,z)
        out.append((h0,h1,ad[0]+ad[1]*c.zi+ad[2]*c.zi*c.zi))
    return out

# Original routines resolve their global blocks function at runtime.
rf.blocks=blocks

def solve(c,U,Y,Uh,W,b0,it=80,tol=1e-6):
    """Same IV-GN step as original; exposes full final score and failures."""
    b=b0.copy(); hist=[]; numerical_failure=None
    def moments(beta):
        ub=Uh(beta) if callable(Uh) else Uh
        yh,J=out_jac(c,beta,U)
        Z=J if ub is U else out_jac(c,beta,ub)[1]
        zw=(Z.conj()*W[:,:,None]).reshape(-1,NB)
        A=np.real(zw.T@J.reshape(-1,NB))
        g=np.real(zw.T@(Y-yh).ravel())
        return A,g,yh,Z
    for n in range(it):
        try:
            A,g,yh,Z=moments(b)
            sc=1/np.sqrt(np.abs(np.diag(A))+1e-300)
            step=sc*np.linalg.solve(A*sc[:,None]*sc[None,:],g*sc)
            mx=np.max(np.abs(step[IW:ITR]))
            if mx>.1: step*=.1/mx
            b+=step
            hist.append(dict(iteration=n+1,max_log_modal_step=float(np.max(np.abs(step[IW:ITR])))))
            if not np.all(np.isfinite(b)): raise FloatingPointError('nonfinite beta')
            if max(np.exp(b[IW:IZ]))>2*np.pi*200 or min(np.exp(b[IW:IZ]))<2*np.pi*.1 or max(np.exp(b[IZ:ITR]))>3:
                raise FloatingPointError('modal bounds exceeded; no clipping')
            if hist[-1]['max_log_modal_step']<tol: break
        except Exception as e:
            numerical_failure=repr(e);break
    try:
        A,g,yh,Z=moments(b)
        den=np.linalg.norm(Z.reshape(-1,NB),axis=0)*np.linalg.norm(Y-yh)
        score=float(np.max(np.abs(g)/np.maximum(den,1e-300)))
        sc=1/np.sqrt(np.abs(np.diag(A))+1e-300)
        cond=float(np.linalg.cond(A*sc[:,None]*sc[None,:]))
        resid=float(np.linalg.norm(Y-yh)/np.linalg.norm(Y))
    except Exception: score=cond=resid=float('nan')
    small=bool(hist and hist[-1]['max_log_modal_step']<tol and numerical_failure is None)
    return b,dict(iterations=len(hist),converged_pole_step=small,score_normalized=score,
                  strict_score_pass=bool(np.isfinite(score) and score<1e-6),scaled_condition=cond,
                  rel_residual=resid,exception=numerical_failure,history=hist)
