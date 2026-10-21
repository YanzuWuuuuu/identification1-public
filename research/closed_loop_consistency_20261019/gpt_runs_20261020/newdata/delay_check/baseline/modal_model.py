"""Output-block permutation of the existing additive ZOH model.

Per output: (mode C0[4], C1[4], transient[2])*m, D0[4], D2[4].
Then common log angular frequencies and log damping ratios. No structural
constraint on the 4x4 residues. Analytic Frechet derivatives of the ZOH map.
"""
import numpy as np
from scipy.linalg import expm,expm_frechet

class ModalModel:
    def __init__(self,nm,tau,wr,fs=5000.):
        self.nm=nm;self.L=10*nm+8;self.iw=4*self.L;self.iz=self.iw+nm;self.nb=self.iz+nm
        self.tau=tau;self.wr=2*np.pi*wr;self.fs=fs
    def context(self,f):
        f=np.asarray(f);s=2j*np.pi*f;zi=np.exp(-s/self.fs)
        return dict(f=f,s=s,zi=zi,z=1/zi,E=np.exp(-s[:,None,None]*self.tau[None]))
    def bounds(self,lower,upper):
        lb=np.full(self.nb,-np.inf);ub=-lb
        lb[self.iw:self.iz]=np.log(2*np.pi*np.asarray(lower));ub[self.iw:self.iz]=np.log(2*np.pi*np.asarray(upper))
        lb[self.iz:]=np.log(1e-5);ub[self.iz:]=np.log(.3)
        return lb,ub
    def blocks(self,c,b,deriv=True):
        hs=[]
        for j in range(self.nm):
            w,ze=np.exp(b[self.iw+j]),np.exp(b[self.iz+j])
            A=np.array([[0.,1.,0.],[-w*w,-2*ze*w,1.],[0.,0.,0.]])/self.fs
            E=expm(A);Ad=E[:2,:2];Bd=E[:2,2];det=np.linalg.det(Ad)
            tr=np.trace(Ad);z=c['z'];zi=c['zi'];den=z*z-tr*z+det;adz=1-tr*zi+det*zi*zi
            n0=(z-Ad[1,1])*Bd[0]+Ad[0,1]*Bd[1]
            n1=Ad[1,0]*Bd[0]+(z-Ad[0,0])*Bd[1]
            h0=w*w*n0/den;h1=w*n1/den;ds=[]
            if deriv:
                for q in range(2):
                    dA=np.zeros((3,3));dA[1,0]=-2*w*w/self.fs if q==0 else 0.;dA[1,1]=-2*ze*w/self.fs
                    dE=expm_frechet(A,dA,compute_expm=False);dAd=dE[:2,:2];dBd=dE[:2,2]
                    dtr=np.trace(dAd);ddet=det*np.trace(np.linalg.solve(Ad,dAd));dden=-dtr*z+ddet
                    dn0=-dAd[1,1]*Bd[0]+(z-Ad[1,1])*dBd[0]+dAd[0,1]*Bd[1]+Ad[0,1]*dBd[1]
                    dn1=dAd[1,0]*Bd[0]+Ad[1,0]*dBd[0]-dAd[0,0]*Bd[1]+(z-Ad[0,0])*dBd[1]
                    dh0=(2*h0 if q==0 else 0.)+w*w*(dn0/den-n0*dden/den**2)
                    dh1=(h1 if q==0 else 0.)+w*(dn1/den-n1*dden/den**2)
                    ds.append((dh0,dh1,-dtr*zi+ddet*zi*zi))
            hs.append((h0,h1,adz,ds))
        return hs
    def compact(self,c,b,U,deriv=True):
        n=len(U);lin=b[:self.iw].reshape(4,self.L);EU=c['E']*U[:,None,:]
        LL=np.empty((n,4,self.L),complex);N=np.zeros((n,4,2*self.nm),complex)
        for j,(h0,h1,den,ds) in enumerate(self.blocks(c,b,deriv)):
            a=10*j;LL[:,:,a:a+4]=h0[:,None,None]*EU;LL[:,:,a+4:a+8]=h1[:,None,None]*EU
            LL[:,:,a+8]=1/den[:,None];LL[:,:,a+9]=c['zi'][:,None]/den[:,None]
            if deriv:
                c0=np.einsum('koi,oi->ko',EU,lin[:,a:a+4]);c1=np.einsum('koi,oi->ko',EU,lin[:,a+4:a+8])
                tn=lin[:,a+8][None]+c['zi'][:,None]*lin[:,a+9][None]
                for q,(dh0,dh1,dden) in enumerate(ds):
                    N[:,:,q*self.nm+j]=dh0[:,None]*c0+dh1[:,None]*c1-tn*dden[:,None]/den[:,None]**2
        a=10*self.nm;LL[:,:,a:a+4]=EU;LL[:,:,a+4:a+8]=(c['s']/self.wr)[:,None,None]**2*EU
        y=np.einsum('kol,ol->ko',LL,lin)
        return y,LL,N
    def output(self,c,b,U,need_jac=True):
        y,LL,N=self.compact(c,b,U,need_jac)
        if not need_jac:return y,None
        J=np.zeros((len(U),4,self.nb),complex)
        for o in range(4):J[:,o,o*self.L:(o+1)*self.L]=LL[:,o]
        J[:,:,self.iw:]=N
        return y,J
    def gmodel(self,c,b,deriv=True):
        n=len(c['f']);lin=b[:self.iw].reshape(4,self.L);G=np.zeros((n,4,4),complex)
        Gj=np.zeros((n,4,4,self.nb),complex) if deriv else None
        for j,(h0,h1,_,ds) in enumerate(self.blocks(c,b,deriv)):
            a=10*j;G+=h0[:,None,None]*lin[:,a:a+4][None]+h1[:,None,None]*lin[:,a+4:a+8][None]
            if deriv:
                for o in range(4):
                    for i in range(4):
                        Gj[:,o,i,o*self.L+a+i]=h0;Gj[:,o,i,o*self.L+a+4+i]=h1
                for q,(dh0,dh1,_) in enumerate(ds):Gj[:,:,:,self.iw+q*self.nm+j]=dh0[:,None,None]*lin[:,a:a+4][None]+dh1[:,None,None]*lin[:,a+4:a+8][None]
        a=10*self.nm;h=(c['s']/self.wr)**2;G+=lin[:,a:a+4][None]+h[:,None,None]*lin[:,a+4:a+8][None]
        if deriv:
            for o in range(4):
                for i in range(4):Gj[:,o,i,o*self.L+a+i]=1;Gj[:,o,i,o*self.L+a+4+i]=h
            Gj*=c['E'][:,:,:,None]
        return G*c['E'],Gj
    def resp(self,c,b):return self.gmodel(c,b,False)[0]
    def modal(self,b):return np.exp(b[self.iw:self.iz])/2/np.pi,100*np.exp(b[self.iz:])
    def initial(self,c,U,Y,f,zeta):
        b=np.zeros(self.nb);b[self.iw:self.iz]=np.log(2*np.pi*np.asarray(f));b[self.iz:]=np.log(zeta)
        _,ll,_=self.compact(c,b,U,False);ranks=[]
        for o in range(4):
            x=np.vstack((ll[:,o].real,ll[:,o].imag));y=np.r_[Y[:,o].real,Y[:,o].imag]
            sc=1/np.maximum(np.linalg.norm(x,axis=0),1e-300)
            sol,_,rank,_=np.linalg.lstsq(x*sc,y,rcond=1e-12)
            b[o*self.L:(o+1)*self.L]=sol*sc;ranks.append(int(rank))
        return b,ranks

def equations(m,c,U,Y,Uh,b,W=None,chunk=1800):
    """Exactly Re sum Z^H W J and Z^H W e, exploiting output sparsity."""
    A=np.zeros((m.nb,m.nb));g=np.zeros(m.nb);zz=np.zeros(m.nb);ee=0.
    for start in range(0,len(U),chunk):
        sl=slice(start,min(start+chunk,len(U)));cc={k:v[sl] for k,v in c.items()}
        yh,jl,jn=m.compact(cc,b,U[sl]);_,zl,zn=(yh,jl,jn) if Uh is U else m.compact(cc,b,Uh[sl])
        e=Y[sl]-yh
        if W is None:we=e;wjn=jn;wzn=zn;ee+=np.sum(abs(e)**2)
        else:
            w=W[sl];we=np.einsum('koi,ki->ko',w,e);wjn=np.einsum('koi,kip->kop',w,jn);wzn=np.einsum('koi,kip->kop',w,zn)
            ee+=np.einsum('ko,ko->',e.conj(),we).real
        nn=slice(m.iw,m.nb)
        A[nn,nn]+=np.einsum('kop,koq->pq',zn.conj(),wjn,optimize=True).real
        g[nn]+=np.einsum('kop,ko->p',zn.conj(),we).real
        zz[nn]+=np.einsum('kop,kop->p',zn.conj(),wzn).real
        for o in range(4):
            oi=slice(o*m.L,(o+1)*m.L)
            g[oi]+=(zl[:,o].conj().T@we[:,o]).real
            zz[oi]+=np.sum(abs(zl[:,o])**2*(1 if W is None else w[:,o,o].real[:,None]),axis=0)
            A[oi,nn]+=(zl[:,o].conj().T@wjn[:,o]).real
            for i in (range(4) if W is not None else (o,)):
                ii=slice(i*m.L,(i+1)*m.L);v=jl[:,i]*(1 if W is None else w[:,o,i,None])
                A[oi,ii]+=(zl[:,o].conj().T@v).real
            # Sum conjugate Z_n W_{*,o}; Hermitian W makes this wzn_o^H.
            A[nn,oi]+=(wzn[:,o].conj().T@jl[:,o]).real
    score=float(np.max(abs(g)/np.sqrt(np.maximum(zz*ee,1e-300))))
    return A,g,score,float(ee)

def iv_solve(m,c,U,Y,Uh,W,b0,lb,ub,maxit=80):
    b=b0.copy();hist=[];status='iteration_limit';exc=None
    try:
        for it in range(maxit):
            u=Uh(b) if callable(Uh) else Uh;A,g,score,ee=equations(m,c,U,Y,u,b,W)
            sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300));st=sc*np.linalg.solve(A*sc[:,None]*sc[None,:],g*sc)
            raw=float(np.max(abs(st[m.iw:])));fac=min(1.,.05/max(raw,1e-300))
            pos=st>0;neg=st<0
            if np.any(pos):fac=min(fac,float(np.min(.995*(ub[pos]-b[pos])/st[pos])))
            if np.any(neg):fac=min(fac,float(np.min(.995*(lb[neg]-b[neg])/st[neg])))
            if fac<1e-7 and raw>1e-5:status='boundary_blocked';break
            b+=fac*st;hist.append(dict(iteration=it+1,raw_pole_step=raw,step_factor=fac,score=score,residual_squared=ee))
            if not np.isfinite(b).all():status='nonfinite';break
            if raw<1e-5:status='pole_step_converged';break
    except Exception as e:status='exception';exc=repr(e)
    try:
        A,g,score,ee=equations(m,c,U,Y,Uh(b) if callable(Uh) else Uh,b,W)
        sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300));cond=float(np.linalg.cond(A*sc[:,None]*sc[None,:]))
    except Exception:score=ee=cond=float('nan')
    return b,dict(success=status=='pole_step_converged',status=status,iterations=len(hist),score_max=score,
        scaled_condition=cond,weighted_residual_squared=ee,exception=exc,history=hist)
