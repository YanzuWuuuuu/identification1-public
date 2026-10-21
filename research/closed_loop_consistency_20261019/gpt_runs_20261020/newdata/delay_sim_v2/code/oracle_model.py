"""Approved six-mode oracle model, shared by fits and population calculations."""
import os, sys, json, time
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np
import scipy.linalg as la
from scipy import signal

ROOT=Path(__file__).resolve().parents[1]
NEW=ROOT.parent
sys.path.insert(0,str(NEW/'delay_check/baseline'))
from modal_model import ModalModel
sys.path.insert(0,str(ROOT/'controller_check'))
import check_controller as ctl
FS=5000.; P0=30000; UNIT=1e-4
FM=np.array([10.6,11.7,16.1,286.573,344.945,490.390])
ZM=np.array([.03,.04,.05,.0119,.0132,.0186])
NI=np.array([0,1,0,1]); NO=np.array([1,0,1,0])
TAU=(NO[:,None]+NI[None,:])/FS
HELM=la.helmert(4,full=False).T
DESIGN=np.zeros((4,4,7)); DESIGN[:,:,0]=1
DESIGN[:,:,1:4]=HELM[:,None,:]; DESIGN[:,:,4:]=HELM[None,:,:]
BANDS={'L':(3.,100.,[0,1,2]),'H':(100.,550.,[3,4,5]),'W':(3.,550.,list(range(6)))}
with np.load(ROOT/'controller_check/six_mode_definition.npz') as fh:
    C0=fh['C0'].copy(); AD=fh['Ad'].copy(); BD=fh['Bd'].copy(); CD=fh['Cd'].copy()

def enc(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,np.generic): return x.item()
    if isinstance(x,Path): return str(x)
    raise TypeError(type(x).__name__)

def save(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=enc),encoding='utf8')
    # Windows readers/indexers may briefly deny an atomic replacement.
    # Keep the complete temporary file and retry, never truncate the target.
    for attempt in range(30):
        try:tmp.replace(path);return
        except PermissionError:
            if attempt==29:raise
            time.sleep(.2)

def rr(x): return np.concatenate([x.real,x.imag],axis=0)

class OracleModel:
    def __init__(self,band):
        self.band=band; self.lo,self.hi,active=BANDS[band]
        self.active=np.array(active); self.outside=np.setdiff1d(np.arange(6),self.active)
        self.nm=len(active); self.L=4*self.nm+23; self.nl=4*self.L
        self.nq=2*self.nm+7; self.nb=self.nl+self.nq
        self.id0=4*self.nm; self.id2=self.id0+4; self.it=self.id2+4; self.ip=self.it+12
        self.base=ModalModel(6,np.zeros((4,4)),self.hi,fs=FS)
        self.lower=np.r_[np.log(2*np.pi*np.array([8,11.1,13,240,310,420])[self.active]),
                         np.full(self.nm,np.log(1e-5)),np.full(7,-100.)]
        self.upper=np.r_[np.log(2*np.pi*np.array([11.1,13,20,310,400,600])[self.active]),
                         np.full(self.nm,np.log(.3)),np.full(7,100.)]
    def truth(self):
        l=np.zeros((4,self.L))
        for j,i in enumerate(self.active): l[:,4*j:4*j+4]=C0[i]
        q=np.r_[np.log(2*np.pi*FM[self.active]),np.log(ZM[self.active]),
                2.,HELM.T@np.array([1.,-1.,1.,-1.]),HELM.T@np.array([-1.,1.,-1.,1.])]
        return np.r_[l.ravel(),q]
    def initial(self,r):
        b=self.truth(); q=b[self.nl:]
        rng=np.random.default_rng(100200+r)
        dw=rng.normal(size=6)*.01; dz=rng.normal(size=6)*.4
        q[:self.nm]+=dw[self.active]; q[self.nm:2*self.nm]+=dz[self.active]
        q[-7:]=0.; q[-7]=2.
        return b
    def context(self,f): return self.base.context(np.asarray(f))
    def tau(self,q): return UNIT*np.einsum('oip,p->oi',DESIGN,q[-7:])
    def blocks(self,c,q,deriv=True):
        b=np.zeros(self.base.nb); b[self.base.iw:self.base.iz]=np.log(2*np.pi*FM)
        b[self.base.iz:]=np.log(ZM)
        b[self.base.iw+self.active]=q[:self.nm]
        b[self.base.iz+self.active]=q[self.nm:2*self.nm]
        return self.base.blocks(c,b,deriv)
    def design(self,c,q,U,deriv=True,transient_mask=None):
        """Y=offset+X*l. dX/dq and offset derivative include all fixed tails."""
        n=len(U); h=self.blocks(c,q,deriv); s=c['s']; zi=c['zi']
        E=np.exp(-s[:,None,None]*self.tau(q)); EU=E*U[:,None,:]
        X=np.zeros((n,4,self.L),complex)
        dX=np.zeros((n,4,self.L,self.nq),complex) if deriv else None
        Go=sum((h[i][0][:,None,None]*C0[i] for i in self.outside),start=np.zeros((n,4,4),complex))*E
        offset=np.einsum('koi,ki->ko',Go,U)
        dO=np.zeros((n,4,self.nq),complex) if deriv else None
        for j,i in enumerate(self.active):
            X[:,:,4*j:4*j+4]=h[i][0][:,None,None]*EU
            if deriv:
                for p in range(2): dX[:,:,4*j:4*j+4,p*self.nm+j]=h[i][3][p][0][:,None,None]*EU
        X[:,:,self.id0:self.id0+4]=EU
        X[:,:,self.id2:self.id2+4]=(s/self.base.wr)[:,None,None]**2*EU
        for i in range(6):
            den=h[i][2]; at=self.it+2*i
            X[:,:,at]=1/den[:,None]; X[:,:,at+1]=zi[:,None]/den[:,None]
            if deriv and i in self.active:
                j=int(np.where(self.active==i)[0][0])
                for p in range(2):
                    dd=h[i][3][p][2]
                    dX[:,:,at,p*self.nm+j]=-dd[:,None]/den[:,None]**2
                    dX[:,:,at+1,p*self.nm+j]=-zi[:,None]*dd[:,None]/den[:,None]**2
        X[:,:,self.ip:]=np.stack([zi**p for p in range(3)],axis=-1)[:,None,:]
        if transient_mask is not None:
            X[:,:,self.it:]*=np.asarray(transient_mask)[:,None,None]
            if deriv: dX[:,:,self.it:,:]*=np.asarray(transient_mask)[:,None,None,None]
        if deriv:
            # Delay affects only transfer terms, never finite-record nuisance.
            for p in range(7):
                hdel=-s[:,None,None]*UNIT*DESIGN[:,:,p]
                for j in range(self.nm): dX[:,:,4*j:4*j+4,2*self.nm+p]=X[:,:,4*j:4*j+4]*hdel
                dX[:,:,self.id0:self.id0+4,2*self.nm+p]=X[:,:,self.id0:self.id0+4]*hdel
                dX[:,:,self.id2:self.id2+4,2*self.nm+p]=X[:,:,self.id2:self.id2+4]*hdel
                dO[:,:,2*self.nm+p]=np.einsum('koi,ki->ko',Go*hdel,U)
        return offset,X,dO,dX
    def compact(self,c,b,U,deriv=True,transient_mask=None):
        l=b[:self.nl].reshape(4,self.L); q=b[self.nl:]
        off,X,_,_=self.design(c,q,U,False,transient_mask)
        y=off+np.einsum('kol,ol->ko',X,l)
        Jn=None
        if deriv:
            h=self.blocks(c,q,True);s=c['s'];zi=c['zi']
            E=np.exp(-s[:,None,None]*self.tau(q));EU=E*U[:,None,:]
            Jn=np.zeros((len(U),4,self.nq),complex)
            for j,i in enumerate(self.active):
                cu=np.einsum('koi,oi->ko',EU,l[:,4*j:4*j+4])
                tn=l[:,self.it+2*i][None]+zi[:,None]*l[:,self.it+2*i+1][None]
                if transient_mask is not None:tn*=np.asarray(transient_mask)[:,None]
                for p in range(2):
                    Jn[:,:,p*self.nm+j]=h[i][3][p][0][:,None]*cu-tn*h[i][3][p][2][:,None]/h[i][2][:,None]**2
            G=self.resp(c,b)
            for p in range(7):Jn[:,:,2*self.nm+p]=np.einsum('koi,ki->ko',-s[:,None,None]*UNIT*G*DESIGN[:,:,p],U)
        return y,X,Jn
    def resp(self,c,b):
        q=b[self.nl:]; l=b[:self.nl].reshape(4,self.L); h=self.blocks(c,q,False)
        G=sum((h[i][0][:,None,None]*C0[i] for i in self.outside),start=np.zeros((len(c['f']),4,4),complex))
        for j,i in enumerate(self.active): G+=h[i][0][:,None,None]*l[:,4*j:4*j+4]
        G+=l[:,self.id0:self.id0+4]+(c['s']/self.base.wr)[:,None,None]**2*l[:,self.id2:self.id2+4]
        return G*np.exp(-c['s'][:,None,None]*self.tau(q))
    def jac(self,c,b,U,transient_mask=None):
        y,X,N=self.compact(c,b,U,True,transient_mask)
        J=np.zeros((len(U),4,self.nb),complex)
        for o in range(4): J[:,o,o*self.L:(o+1)*self.L]=X[:,o]
        J[:,:,self.nl:]=N
        return y,J

def plant(fc,open_loop=False):
    """Exact ZOH plant+integer delay registers and controller; state recursion."""
    ak,bk,ck,dk=ctl.cascaded_controller(fc,None)
    np_,nk=len(AD),len(ak); ii=[];oo=[];off=np_
    for nn in NI: ii.append(np.arange(off,off+nn)); off+=nn
    for nn in NO: oo.append(np.arange(off,off+nn)); off+=nn
    op=off; n=op if open_loop else op+nk
    def step(x,e):
        raw=CD@x[:np_]
        yp=np.array([x[ix[-1]] if len(ix) else raw[j] for j,ix in enumerate(oo)])
        y=yp+e[4:8]
        u=e[:4] if open_loop else e[:4]-ck@x[op:]-dk@y
        uw=u+e[8:]; up=np.array([x[ix[-1]] if len(ix) else uw[j] for j,ix in enumerate(ii)])
        xn=np.zeros(n);xn[:np_]=AD@x[:np_]+BD@up
        if not open_loop: xn[op:]=ak@x[op:]+bk@y
        for j,ix in enumerate(ii):
            if len(ix):xn[ix]=np.r_[uw[j],x[ix[:-1]]]
        for j,ix in enumerate(oo):
            if len(ix):xn[ix]=np.r_[raw[j],x[ix[:-1]]]
        return xn,np.r_[y,u]
    ac=np.zeros((n,n));bc=np.zeros((n,12));cc=np.zeros((8,n));dc=np.zeros((8,12))
    for i in range(n): ac[:,i],cc[:,i]=step(np.eye(n)[i],np.zeros(12))
    for i in range(12):bc[:,i],dc[:,i]=step(np.zeros(n),np.eye(12)[i])
    return ac,bc,cc,dc

def excitation(r):
    rng=np.random.default_rng(1002+r);k=np.arange(1,P0//2)
    dp=np.zeros((P0//2+1,4),complex)
    shape=ctl.cl2.dshape(np.exp(2j*np.pi*k/P0))
    dp[k]=.2*np.sqrt(P0*shape)[:,None]*np.exp(2j*np.pi*rng.random((len(k),4)))
    return dp,rng

def deterministic_spectra(P,r,band,fc,method):
    N=P*P0;lo,hi,_=BANDS[band];k=np.arange(int(lo*N/FS),int(hi*N/FS)+1)
    if method=='M3': k=k[k%P==0]
    f=k*FS/N;dp,_=excitation(r);D=np.zeros((len(k),4),complex);ix=k%P==0
    D[ix]=P*dp[k[ix]//P]
    m=OracleModel(band);c=m.context(f);G=m.resp(c,m.truth());K=ctl.kresp(f,fc)
    S=np.linalg.inv(np.eye(4)+K@G);F=S@K;Sy=np.eye(4)-G@F
    U=np.einsum('koi,ki->ko',S,D);Y=np.einsum('koi,ki->ko',G,U)
    Pv=N*(3e-7**2*np.eye(4)+.03**2*(G@G.conj().transpose(0,2,1)))
    return dict(P=P,r=r,fc=fc,band=band,method=method,N=N,k=k,f=f,D=D,U=U,Y=Y,G=G,K=K,S=S,F=F,Sy=Sy,Pv=Pv)

def generate(P,r,fc):
    path=ROOT/'pilot/data'/f'FC{fc}_P{P}_r{r:03d}.npz'
    if path.exists():
        with np.load(path) as a:return {k:a[k] for k in a.files}
    dp,rng=excitation(r);d=np.tile(np.fft.irfft(dp,n=P0,axis=0),(P+2,1));nn=len(d)
    v=rng.normal(size=(nn,4))*3e-7;w=rng.normal(size=(nn,4))*.03
    a,b,c,e=plant(fc);rad=float(max(abs(la.eigvals(a))))
    assert rad<1
    _,out,_=signal.dlsim((a,b,c,e,1/FS),np.c_[d,v,w])
    N=P*P0;k=np.arange(int(3*N/FS),int(550*N/FS)+1)
    obj=dict(P=P,r=r,fc=fc,N=N,k=k,f=k*FS/N,radius=rad)
    for key,x in [('D',d),('U',out[:,4:]),('Y',out[:,:4])]:obj[key]=np.fft.rfft(x[2*P0:],axis=0)[k]
    path.parent.mkdir(exist_ok=True,parents=True);np.savez_compressed(path,**obj)
    return obj

def select(data,band,method):
    lo,hi,_=BANDS[band];ix=(data['f']>=lo-1e-9)&(data['f']<=hi+1e-9)
    if method=='M3':ix&=data['k']%int(data['P'])==0
    out={k:data[k][ix] for k in ('f','k','U','Y','D')};out['K']=ctl.kresp(out['f'],int(data['fc']))
    return out
