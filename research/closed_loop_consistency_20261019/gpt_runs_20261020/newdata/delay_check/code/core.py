"""Isolated delay feasibility campaign. Existing files are read only."""
import os,sys
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
from pathlib import Path
import json,time,traceback
import numpy as np
import scipy.linalg as la
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT.parent
HIGH=SRC/'high_modal_validation'
sys.path.insert(0,str(ROOT/'baseline'))
from modal_model import ModalModel,equations
FS=5000.;P0=30000
TAU=np.load(ROOT/'baseline/tau_seconds.npy')
CTRS=['lowgain','nominal','highgain']
IDS=[g+'_'+c+('_pppp' if g=='C' else '') for g in ['A','C'] for c in CTRS]
SPECS={
 'L':dict(band=[3,25],f=[10.6,11.7,16.1],z=[.03]*3,lo=[8,10.8,14],hi=[11.5,14,19],report=[0,1,2]),
 'H300':dict(band=[220,420],f=[263,282,340,360],z=[.059,.023,.0045,.003],lo=[220,250,310,345],hi=[290,310,355,400],report=[1,2,3]),
 'H493':dict(band=[450,550],f=[493],z=[.0035],lo=[460],hi=[530],report=[0])}
HELM=la.helmert(4,full=False).T
DESIGN=np.zeros((4,4,7));DESIGN[:,:,0]=1
DESIGN[:,:,1:4]=HELM[:,None,:];DESIGN[:,:,4:7]=HELM[None,:,:]
UNIT=1e-4
def rd(p):return json.loads(Path(p).read_text(encoding='utf8'))
def js(p,x):
    def enc(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,Path):return str(v)
        raise TypeError(type(v).__name__)
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=enc),encoding='utf8');q.replace(p)
def delay_values(tau):
    c=tau.mean();a=tau.mean(1)-c;b=tau.mean(0)-c
    return dict(common_us=c*1e6,output_us=a*1e6,input_us=b*1e6,matrix_us=tau*1e6,
                residual_us=(tau-c-a[:,None]-b[None,:])*1e6)
def evalK(f,rid):
    ct=next(c for c in CTRS if c in rid)
    a=np.load(SRC/'controllers'/f'{ct}.npz');A,B,T=(a[k] for k in ('A','B','T'))
    zi=np.exp(-2j*np.pi*np.asarray(f)/FS)
    aa=np.eye(4)[None]-sum(A[i][None]*zi[:,None,None]**(i+1) for i in range(8))
    bb=sum(B[i][None]*zi[:,None,None]**i for i in range(9))
    return -T[None]@np.linalg.solve(aa,bb)
def spectrum(rid,w,omit=None,extended=False):
    a=np.load(SRC/'prepared'/f'{rid}.npz');p=int(a['P']);keep=np.arange(p)
    if omit is not None:keep=keep[keep!=omit]
    P=len(keep);N=P*P0;lo,hi=SPECS[w]['band']
    if extended:lo-=2;hi+=2
    k=np.arange(int(np.ceil(lo*N/FS)),int(np.floor(hi*N/FS))+1);f=k*FS/N
    out=dict(k=k,f=f,N=N,P=P,line_indices=np.flatnonzero(k%P==0),kept_periods=keep)
    for name in ('U','Y','D'):
        xx=a[name].reshape(p,P0,4)[keep].reshape(N,4);out[name]=np.fft.rfft(xx,axis=0)[k]
    out['K']=evalK(f,rid);return out
def old_order():
    order=[]
    for o in range(4):
        for j in range(3):order+=list(range(32*j+4*o,32*j+4*o+4))+list(range(32*j+16+4*o,32*j+16+4*o+4))+[134+8*j+2*o,135+8*j+2*o]
        order+=list(range(96+4*o,100+4*o))+list(range(112+4*o,116+4*o))
    return order+list(range(128,134))
def basefit(rid,w,method):
    if w=='L':
        r=rd(SRC/'fits'/f'{rid}_{method}.json');return np.asarray(r['beta'])[old_order()],r
    r=rd(HIGH/'diagnostics'/f'{rid}_{w}_{method}.json');return np.asarray(r['beta']),r

class DelayModel(ModalModel):
    """Seven gauge-fixed, jointly estimated delays; all original nuisance terms retained."""
    def __init__(self,nm,tau,wr,fs=5000.):
        super().__init__(nm,tau,wr,fs);self.base_nb=self.nb;self.nb+=7
    def ctx(self,c,b):
        cc=c.copy();cc['E']=np.exp(-c['s'][:,None,None]*self.delay(b)[None]);return cc
    def delay(self,b):return self.tau+UNIT*np.einsum('oip,p->oi',DESIGN,b[self.base_nb:])
    def compact(self,c,b,U,deriv=True):
        cc=self.ctx(c,b);y,L,N=super().compact(cc,b,U,deriv)
        if deriv:
            G=super().gmodel(cc,b,False)[0]
            Jt=-c['s'][:,None,None]*UNIT*np.einsum('koi,ki,oip->kop',G,U,DESIGN)
        else:Jt=np.zeros((len(U),4,7),complex)
        return y,L,np.concatenate((N,Jt),axis=2)
    def gmodel(self,c,b,deriv=True):
        cc=self.ctx(c,b);G,J=super().gmodel(cc,b,deriv)
        if deriv:J[:,:,:,self.base_nb:]=-c['s'][:,None,None,None]*UNIT*G[:,:,:,None]*DESIGN[None]
        return G,J
    def modal(self,b):return np.exp(b[self.iw:self.iz])/2/np.pi,100*np.exp(b[self.iz:self.base_nb])
    def bounds(self,lower,upper):
        lb,ub=super().bounds(lower,upper);lb[self.base_nb:]=-100.;ub[self.base_nb:]=100.
        return lb,ub

def setup(w,shift=0.,joint=False):
    s=SPECS[w];cls=DelayModel if joint else ModalModel
    m=cls(len(s['f']),TAU+shift*1e-6,s['band'][1]);lb,ub=m.bounds(s['lo'],s['hi'])
    return s,m,lb,ub
def inputs(m,a,method):
    ix=slice(None) if method=='M1' else a['line_indices'];c=m.context(a['f'][ix])
    U,Y,D,K=(a[k][ix] for k in ('U','Y','D','K'))
    Uh=U if method=='M1' else lambda b:np.linalg.solve(np.eye(4)[None]+K@m.resp(c,b),D[:,:,None])[:,:,0]
    return c,U,Y,Uh
def linsolve(m,c,U,Y,Uh,b):
    """One exact frozen-pole IV solve of original linear nuisance parameters."""
    u=Uh(b) if callable(Uh) else Uh;A,g,_,_=equations(m,c,U,Y,u,b)
    aa=A[:m.iw,:m.iw];sc=1/np.sqrt(np.maximum(abs(np.diag(aa)),1e-300))
    bb=b.copy();bb[:m.iw]+=sc*la.lstsq(aa*sc[:,None]*sc[None,:],g[:m.iw]*sc,cond=1e-13,lapack_driver='gelsy')[0]
    return bb
def solve(m,c,U,Y,Uh,b0,lb,ub,maxit=100,relaxation=1.):
    b=b0.copy();hist=[];status='iteration_limit';best=(np.inf,b.copy())
    for it in range(maxit):
        try:
            u=Uh(b) if callable(Uh) else Uh;A,g,score,ee=equations(m,c,U,Y,u,b)
            if score<best[0]:best=score,b.copy()
            if score<1e-8:status='score_converged';break
            sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300));As=A*sc[:,None]*sc[None,:]
            st=sc*la.lstsq(As,g*sc,cond=1e-13,lapack_driver='gelsy')[0]
            end=m.base_nb if isinstance(m,DelayModel) else m.nb
            fac=min(1.,.05/max(abs(st[m.iw:end]).max(),1e-300))
            if end<m.nb:fac=min(fac,.5/max(abs(st[end:]).max(),1e-300))
            pos=st>0;neg=st<0
            if pos.any():fac=min(fac,float(np.min(.995*(ub[pos]-b[pos])/st[pos])))
            if neg.any():fac=min(fac,float(np.min(.995*(lb[neg]-b[neg])/st[neg])))
            if fac<1e-8:status='boundary_or_stalled';break
            fac*=relaxation
            b+=fac*st
            hist.append(dict(iteration=it+1,score=score,ee=ee,step_factor=fac,delay_step_us=(st[end:]*UNIT*1e6).tolist()))
            if not np.isfinite(b).all():status='nonfinite';break
        except Exception as e:status='exception:'+repr(e);break
    # Retain final endpoint and score; never silently pick a favourable failed iterate.
    try:
        A,g,score,ee=equations(m,c,U,Y,Uh(b) if callable(Uh) else Uh,b)
        sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300));sv=la.svdvals(A*sc[:,None]*sc[None,:])
        end=m.base_nb if isinstance(m,DelayModel) else m.nb
        bound=bool(np.min(np.minimum(b[m.iw:]-lb[m.iw:],ub[m.iw:]-b[m.iw:]))<1e-4)
        good=bool(score<1e-6 and not bound)
        return b,dict(success=good,status=status,iterations=len(hist),score=score,ee=ee,near_bound=bound,
                      condition=sv[0]/sv[-1],rank_1e12=int(sum(sv>sv[0]*1e-12)),history=hist,best_score=best[0])
    except Exception as e:return b,dict(success=False,status=status,exception=repr(e),history=hist,iterations=len(hist))
def realvec(x):return np.concatenate((x.real,x.imag),axis=1)
def realmat(C):
    return np.concatenate((np.concatenate((C.real,-C.imag),2),np.concatenate((C.imag,C.real),2)),1)
def local_phi(m,a,b,f,width=.5):
    yh=m.output(m.context(a['f']),b,a['U'],False)[0]
    ix=a['k']%a['P']!=0;v=(a['Y']-yh)[ix];ff=a['f'][ix]
    cs=np.concatenate((np.zeros((1,4,4),complex),np.cumsum(v[:,:,None]*v[:,None,:].conj(),axis=0)))
    l=np.searchsorted(ff,f-width/2-1e-12);r=np.searchsorted(ff,f+width/2+1e-12,side='right')
    Phi=(cs[r]-cs[l])/(r-l)[:,None,None]/a['N'];return (Phi+Phi.conj().transpose(0,2,1))/2
def qe(m,a,b,Phi):
    c,U,Y,Uh=inputs(m,a,'M3');li=a['line_indices']
    yh,J=m.output(c,b,U);_,Z=m.output(c,b,Uh(b));e=realvec(Y-yh)
    Jr=realvec(J);Zr=realvec(Z);Cr=.5*realmat(a['N']*Phi[li]);L=np.linalg.cholesky(Cr)
    ee=np.linalg.solve(L,e[:,:,None])[:,:,0].reshape(-1)
    B=np.linalg.solve(L,Jr).reshape(-1,m.nb)
    A=np.einsum('kop,koq->pq',Zr,Jr,optimize=True)
    ZL=np.einsum('kop,kol->pkl',Zr,L,optimize=True).reshape(m.nb,-1)
    sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
    H=sc[:,None]*la.lstsq(A*sc[:,None]*sc[None,:],sc[:,None]*ZL,cond=1e-13,lapack_driver='gelsy')[0]
    scale=np.maximum(np.linalg.norm(B,axis=0),1e-300);qb,rb=np.linalg.qr(B/scale,mode='reduced')
    rh=rb@(H*scale[:,None])
    eq=len(ee)-2*np.einsum('ip,pi->',qb,rh)+np.sum(rh**2)
    return dict(Q=float(ee@ee),EQ=float(eq),ratio=float((ee@ee)/eq),n=len(ee))
def jk_se(x):
    x=np.asarray(x);P=len(x);return np.sqrt((P-1)/P*np.sum((x-x.mean(0))**2,axis=0))
