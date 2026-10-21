"""Variable-projection IV with an inner reference-instrument fixed point.

The original output, ZOH, transient and IV moment functions are imported read-only.
Only coordinates and the numerical solution of their moment equations change.
"""
import os, sys
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
sys.dont_write_bytecode=True
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT.parent
sys.path.insert(0,str(SRC/'delay_check/code'))
import core as old
import numpy as np
import scipy.linalg as la
from scipy.sparse import csc_matrix
from scipy.optimize import root
import time

js,rd=old.js,old.rd

def scaled_solve(A,b,cut=1e-13):
    sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
    return sc*la.lstsq(A*sc[:,None]*sc[None,:],sc*b,cond=cut,lapack_driver='gelsy')[0]

class ResidueModel:
    """Coordinates for models a/b/c/d; common delay can optionally be fixed."""
    def __init__(self,w,kind='a',joint=True,common_fixed_us=None,pivots=None):
        self.w=w;self.kind=kind;self.nm=len(old.SPECS[w]['f'])
        self.full=old.setup(w,joint=joint)[1];self.joint=joint
        self.common_fixed_us=common_fixed_us
        self.pivots=list(pivots or [0]*self.nm)
        m=self.nm
        self.nres=(32 if kind=='a' else 16 if kind=='c' else 4)*m
        self.nl=self.nres+8*m+32
        self.nd=(7-int(common_fixed_us is not None)) if joint else 0
        self.nshape=(4 if kind=='b' else 3 if kind=='d' else 0)*m
        self.nq=2*m+self.nd+self.nshape
        self.nb=self.nl+self.nq
        self.iw=self.nl;self.iz=self.nl+m
        self.delay_slice=slice(self.nl+2*m,self.nl+2*m+self.nd)
        self.shape_start=self.nl+2*m+self.nd
        s=old.SPECS[w]
        self.lb=np.r_[np.log(2*np.pi*np.array(s['lo'])),np.log(np.full(m,1e-5)),np.full(self.nd,-100.),np.full(self.nshape,-np.inf)]
        self.ub=np.r_[np.log(2*np.pi*np.array(s['hi'])),np.log(np.full(m,.3)),np.full(self.nd,100.),np.full(self.nshape,np.inf)]
    def shape(self,b,j):
        ns=4 if self.kind=='b' else 3
        x=b[self.shape_start+ns*j:self.shape_start+ns*(j+1)]
        pivot=self.pivots[j];ii=[i for i in range(4) if i!=pivot]
        v=np.zeros(4);v[pivot]=1.;v[ii]=x[:3]
        nv=np.linalg.norm(v);psi=v/nv
        dp=(np.eye(4)-np.outer(psi,psi))[:,ii]/nv
        gamma=x[3] if self.kind=='b' else 0.
        return psi,dp,gamma
    def expand(self,b,deriv=True):
        f=self.full;m=self.nm
        fb=np.zeros(f.nb);T=np.zeros((f.nb,self.nb)) if deriv else None
        for j in range(m):
            rows0=np.array([[o*f.L+10*j+i for i in range(4)] for o in range(4)])
            rows1=rows0+4
            if self.kind in ('a','c'):
                width=32 if self.kind=='a' else 16;st=width*j
                fb[rows0.ravel()]=b[st:st+16]
                if deriv:T[rows0.ravel(),np.arange(st,st+16)]=1
                if self.kind=='a':
                    fb[rows1.ravel()]=b[st+16:st+32]
                    if deriv:T[rows1.ravel(),np.arange(st+16,st+32)]=1
            else:
                psi,dp,gamma=self.shape(b,j);phi=b[4*j:4*j+4]
                ca,sa=np.cos(gamma),np.sin(gamma)
                outer=np.outer(psi,phi);fb[rows0]=ca*outer;fb[rows1]=sa*outer
                if deriv:
                    q=self.shape_start+(4 if self.kind=='b' else 3)*j
                    for o in range(4):
                        T[rows0[o],np.arange(4*j,4*j+4)]=ca*psi[o]
                        T[rows1[o],np.arange(4*j,4*j+4)]=sa*psi[o]
                    for h in range(3):
                        T[rows0.ravel(),q+h]=(ca*np.outer(dp[:,h],phi)).ravel()
                        T[rows1.ravel(),q+h]=(sa*np.outer(dp[:,h],phi)).ravel()
                    if self.kind=='b':
                        T[rows0.ravel(),q+3]=(-sa*outer).ravel()
                        T[rows1.ravel(),q+3]=(ca*outer).ravel()
            for o in range(4):
                rows=np.array([o*f.L+10*j+8,o*f.L+10*j+9])
                cols=np.arange(self.nres+8*j+2*o,self.nres+8*j+2*o+2)
                fb[rows]=b[cols]
                if deriv:T[rows,cols]=1
        for o in range(4):
            rows=np.arange(o*f.L+10*m,o*f.L+10*m+8)
            cols=np.r_[np.arange(self.nres+8*m+4*o,self.nres+8*m+4*o+4),
                       np.arange(self.nres+8*m+16+4*o,self.nres+8*m+20+4*o)]
            fb[rows]=b[cols]
            if deriv:T[rows,cols]=1
        rows=np.arange(f.iw,f.iw+2*m);cols=np.arange(self.nl,self.nl+2*m)
        fb[rows]=b[cols]
        if deriv:T[rows,cols]=1
        if self.joint:
            ds=f.base_nb
            if self.common_fixed_us is not None:
                fb[ds]=(self.common_fixed_us-old.TAU.mean()*1e6)/100.
                rows=np.arange(ds+1,ds+7)
            else:rows=np.arange(ds,ds+7)
            cols=np.arange(self.delay_slice.start,self.delay_slice.stop)
            fb[rows]=b[cols]
            if deriv:T[rows,cols]=1
        return fb,T
    def from_full(self,fb):
        f=self.full;m=self.nm;b=np.zeros(self.nb)
        for j in range(m):
            C0=np.array([fb[o*f.L+10*j:o*f.L+10*j+4] for o in range(4)])
            C1=np.array([fb[o*f.L+10*j+4:o*f.L+10*j+8] for o in range(4)])
            if self.kind=='a':b[32*j:32*j+32]=np.r_[C0.ravel(),C1.ravel()]
            elif self.kind=='c':b[16*j:16*j+16]=C0.ravel()
            else:
                U,_,_=la.svd(np.c_[C0,C1] if self.kind=='b' else C0,full_matrices=False)
                psi=U[:,0];pivot=int(np.argmax(abs(psi)));self.pivots[j]=pivot
                if psi[pivot]<0:psi=-psi
                if self.kind=='b':
                    U2,sv,V=la.svd(np.vstack((psi@C0,psi@C1)),full_matrices=False)
                    gamma=np.arctan2(U2[1,0],U2[0,0]);phi=sv[0]*V[0]
                else:gamma=0.;phi=psi@C0
                b[4*j:4*j+4]=phi
                st=self.shape_start+(4 if self.kind=='b' else 3)*j
                b[st:st+3]=psi[[i for i in range(4) if i!=pivot]]/psi[pivot]
                if self.kind=='b':b[st+3]=gamma
            for o in range(4):b[self.nres+8*j+2*o:self.nres+8*j+2*o+2]=fb[o*f.L+10*j+8:o*f.L+10*j+10]
        for o in range(4):
            b[self.nres+8*m+4*o:self.nres+8*m+4*o+4]=fb[o*f.L+10*m:o*f.L+10*m+4]
            b[self.nres+8*m+16+4*o:self.nres+8*m+20+4*o]=fb[o*f.L+10*m+4:o*f.L+10*m+8]
        b[self.nl:self.nl+2*m]=fb[f.iw:f.iw+2*m]
        if self.joint:b[self.delay_slice]=fb[f.base_nb+(self.common_fixed_us is not None):]
        return b
    def output(self,c,b,U,deriv=True):
        fb,T=self.expand(b,deriv);y,J=self.full.output(c,fb,U,deriv)
        return y,(J.reshape(-1,self.full.nb)@csc_matrix(T)).reshape(len(U),4,self.nb) if deriv else None
    def resp(self,c,b):return self.full.resp(c,self.expand(b,False)[0])
    def modal(self,b):return np.exp(b[self.nl:self.nl+self.nm])/2/np.pi,100*np.exp(b[self.nl+self.nm:self.nl+2*self.nm])
    def delay(self,b):return self.full.delay(self.expand(b,False)[0]) if self.joint else self.full.tau
    def gmodel(self,c,b):
        fb,T=self.expand(b);G,J=self.full.gmodel(c,fb)
        return G,(J.reshape(-1,self.full.nb)@csc_matrix(T)).reshape(len(c['f']),4,4,self.nb)

class ProfileIV:
    def __init__(self,m,a,method):
        self.m=m;self.method=method
        ix=slice(None) if method=='M1' else a['line_indices']
        self.c=m.full.context(a['f'][ix])
        self.U,self.Y,self.D,self.K=[a[k][ix] for k in ('U','Y','D','K')]
        self.calls=0;self.inner_failures=0;self.history=[]
    def uh(self,b):
        return self.U if self.method in ('M1','M2') else np.linalg.solve(np.eye(4)[None]+self.K@self.m.resp(self.c,b),self.D[:,:,None])[:,:,0]
    def eq(self,b):
        fb,T=self.m.expand(b);u=self.uh(b)
        Af,gf,_,ee=old.equations(self.m.full,self.c,self.U,self.Y,u,fb)
        A=T.T@Af@T;g=T.T@gf
        # Exact instrument column norms including shared shape columns.
        _,Z=self.m.output(self.c,b,u)
        zz=np.sum(abs(Z)**2,axis=(0,1))
        score=float(np.max(abs(g)/np.sqrt(np.maximum(zz*ee,1e-300))))
        return A,g,score,ee,zz
    def profile(self,q,lin,maxinner=60):
        b=np.r_[lin,q];nl=self.m.nl
        if self.method in ('M1','M2'):
            # With measured-input instruments the linear block is ordinary LS.
            # Form only its nonzero output blocks, not the full pole/delay Jacobian.
            fb,T=self.m.expand(b)
            _,LL,_=self.m.full.compact(self.c,fb,self.U,False)
            A=np.zeros((nl,nl));rhs=np.zeros(nl);blocks=[]
            for o in range(4):
                Tblock=T[o*self.m.full.L:(o+1)*self.m.full.L,:nl]
                cols=np.flatnonzero(np.any(Tblock!=0,axis=0))
                X=LL[:,o]@Tblock[:,cols]
                A[np.ix_(cols,cols)]+=(X.conj().T@X).real
                rhs[cols]+=(X.conj().T@self.Y[:,o]).real
                blocks.append((cols,X))
            b[:nl]=scaled_solve(A,rhs)
            gl=np.zeros(nl);ee=0.
            for o,(cols,X) in enumerate(blocks):
                e=self.Y[:,o]-X@b[cols];gl[cols]+=(X.conj().T@e).real;ee+=np.sum(abs(e)**2)
            sc=float(np.max(abs(gl)/np.sqrt(np.maximum(np.diag(A)*ee,1e-300))))
            ok=bool(sc<1e-9 and np.isfinite(b).all());self.calls+=1
            if not ok:self.inner_failures+=1
            return b,ok,None,1
        hist=[];ok=False
        for k in range(maxinner):
            A,g,score,ee,zz=self.eq(b)
            sc=float(np.max(abs(g[:nl])/np.sqrt(np.maximum(zz[:nl]*ee,1e-300))))
            if sc<2e-12:ok=True;break
            step=scaled_solve(A[:nl,:nl],g[:nl])
            b[:nl]+=step
            hist.append(sc)
            if not np.isfinite(b).all():break
            if self.method in ('M1','M2') and k>=2 and sc<1e-9:ok=True;break
            if len(hist)>8 and hist[-1]>hist[-5]*1.5:break
        if not ok:
            self.inner_failures+=1
        self.calls+=1
        return b,ok,(A,g,score,ee,zz),len(hist)
    def solve(self,b0,maxouter=70,tol=2e-8,hybrid=True):
        start=time.perf_counter();b=b0.copy();hist=[];status='iteration_limit';nl=self.m.nl
        for it in range(maxouter):
            b,ok,eq,inner=self.profile(b[nl:],b[:nl])
            A,g,sc,ee,zz=self.eq(b)
            hist.append(dict(iteration=it,score=sc,inner_iterations=inner,inner_ok=ok))
            if sc<tol and ok:status='score_converged';break
            if not ok:status='inner_linear_moment_not_converged';break
            # Schur step; nuisance coordinates already satisfy their equations.
            sca=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
            As=A*sca[:,None]*sca[None,:];gs=g*sca
            x=la.lstsq(As[:nl,:nl],As[:nl,nl:],cond=1e-13,lapack_driver='gelsy')[0]
            S=As[nl:,nl:]-As[nl:,:nl]@x
            st=sca[nl:]*la.lstsq(S,gs[nl:],cond=1e-12,lapack_driver='gelsy')[0]
            cap=np.r_[np.full(2*self.m.nm,.04),np.full(self.m.nd,.4),np.full(self.m.nshape,.15)]
            fac=min(1.,1/max(np.max(abs(st)/cap),1e-300))
            pos=st>0;neg=st<0
            if pos.any():fac=min(fac,float(np.min(.99*(self.m.ub[pos]-b[nl:][pos])/st[pos])))
            if neg.any():fac=min(fac,float(np.min(.99*(self.m.lb[neg]-b[nl:][neg])/st[neg])))
            if fac<1e-9:status='boundary_or_stalled';break
            b[nl:]+=fac*st
        b,ok,_,_=self.profile(b[nl:],b[:nl]);A,g,sc,ee,zz=self.eq(b)
        hybrid_info=None
        if hybrid and sc>=tol and ok:
            # Differentiate the profiled moments, including instrument updates.
            startb=b.copy();sc0=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
            As=A*sc0[:,None]*sc0[None,:]
            S=As[nl:,nl:]-As[nl:,:nl]@la.lstsq(As[:nl,:nl],As[:nl,nl:],cond=1e-13)[0]
            qscale=np.minimum(sc0[nl:]/np.sqrt(np.maximum(abs(np.diag(S)),1e-300))*np.sqrt(ee),
                              np.r_[np.full(2*self.m.nm,.05),np.full(self.m.nd,1.),np.full(self.m.nshape,.3)])
            qscale=np.maximum(qscale,1e-6)
            denom=np.sqrt(np.maximum(zz[nl:]*ee,1e-300))
            cache={}
            def fun(x):
                q=startb[nl:]+qscale*x
                if np.any(q<self.m.lb) or np.any(q>self.m.ub):
                    return np.full(self.m.nq,1e3)+x
                bb,innerok,_,_=self.profile(q,startb[:nl])
                AA,gg,ss,eee,zzz=self.eq(bb)
                cache['b']=bb
                return gg[nl:]/denom
            opt=root(fun,np.zeros(self.m.nq),method='hybr',options=dict(xtol=2e-9,maxfev=600,factor=.2,eps=1e-9))
            b,ok,_,_=self.profile(startb[nl:]+qscale*opt.x,startb[:nl])
            A,g,sc,ee,zz=self.eq(b)
            hybrid_info=dict(message=str(opt.message),evaluations=opt.nfev,success=bool(opt.success))
            status='profile_hybrid_root' if sc<tol else 'profile_hybrid_failed'
        near=bool(np.min(np.minimum(b[nl:]-self.m.lb,self.m.ub-b[nl:]))<1e-5)
        return b,dict(success=bool(sc<1e-6 and ok and not near),status=status,score=sc,
                     iterations=len(hist),inner_failures=self.inner_failures,profile_calls=self.calls,
                     seconds=time.perf_counter()-start,near_bound=near,history=hist,hybrid=hybrid_info)

def load_joint(rid,w,method):
    name=f'joint_{rid}_{w}_{method}_full.json'
    paths=[SRC/'delay_check'/s/name for s in ('refinements','raw')]
    for p in paths:
        if p.exists():
            r=rd(p);return np.asarray(r['beta']),r
    raise FileNotFoundError(name)
