"""Shared delays and actuator lags over two existing additive modal windows."""
import os,sys
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent
sys.path.insert(0,str(BASE/'code'))
import profile_core as pc
import numpy as np
import scipy.linalg as la
import time
old=pc.old
js,rd=pc.js,pc.rd
WINDOWS=('H300','H493')

class LagFull:
    def __init__(self,w):
        self.base=old.setup(w,joint=True)[1]
        self.L=self.base.L;self.iw=self.base.iw;self.nm=self.base.nm
        self.nb=self.base.nb+4
    def context(self,f):return self.base.context(f)
    def factors(self,c,b):
        x=c['s'][:,None]/np.exp(b[-4:])[None]
        L=1/(1+x)
        return L,L*(1-L)
    def compact(self,c,b,U,deriv=True):
        L,dL=self.factors(c,b)
        y,LL,N=self.base.compact(c,b[:-4],U*L,deriv)
        dn=self.base.resp(c,b[:-4])*(U*dL)[:,None,:] if deriv else np.zeros((len(U),4,4),complex)
        return y,LL,np.concatenate((N,dn),axis=2)
    def output(self,c,b,U,deriv=True):
        y,LL,N=self.compact(c,b,U,deriv)
        if not deriv:return y,None
        J=np.zeros((len(U),4,self.nb),complex)
        for o in range(4):J[:,o,o*self.L:(o+1)*self.L]=LL[:,o]
        J[:,:,self.iw:]=N
        return y,J
    def resp(self,c,b):
        L,_=self.factors(c,b)
        return self.base.resp(c,b[:-4])*L[:,None,:]

class JointModel:
    def __init__(self,kind):
        assert kind in ('a','c')
        self.kind=kind
        self.res={w:pc.ResidueModel(w,kind) for w in WINDOWS}
        self.full={w:LagFull(w) for w in WINDOWS}
        self.nl=sum(m.nl for m in self.res.values())
        self.nm=5;self.nd=7;self.nshape=4 # reuse solver step caps for four log-cutoffs
        self.nq=21;self.nb=self.nl+self.nq
        self.maps={};self.active={};self.linear_blocks={}
        self.delay_slice=slice(self.nl+10,self.nl+17)
        self.lag_slice=slice(self.nl+17,self.nb)
        self.lb=[];self.ub=[]
        lstart=0;qstart=self.nl
        for w,m in self.res.items():
            _,T=m.expand(np.zeros(m.nb))
            local_to_global=np.r_[np.arange(lstart,lstart+m.nl),np.arange(qstart,qstart+2*m.nm),
                                  np.arange(self.delay_slice.start,self.delay_slice.stop)]
            rr,cc=np.nonzero(T)
            assert np.all(T[rr,cc]==1) and len(set(rr))==len(rr)
            fullrows=np.r_[rr,np.arange(m.full.nb,m.full.nb+4)]
            globalcols=np.r_[local_to_global[cc],np.arange(self.lag_slice.start,self.lag_slice.stop)]
            self.maps[w]=(fullrows,globalcols)
            self.active[w]=np.unique(globalcols)
            blocks=[]
            for o in range(4):
                mask=(fullrows>=o*m.full.L)&(fullrows<(o+1)*m.full.L)
                blocks.append((fullrows[mask]-o*m.full.L,globalcols[mask]))
            self.linear_blocks[w]=blocks
            self.lb.extend(m.lb[:2*m.nm]);self.ub.extend(m.ub[:2*m.nm])
            lstart+=m.nl;qstart+=2*m.nm
        self.lb=np.r_[self.lb,np.full(7,-100.),np.full(4,np.log(2*np.pi*20))]
        self.ub=np.r_[self.ub,np.full(7,100.),np.full(4,np.log(2*np.pi*200000))]
    def local(self,w,b):
        r,c=self.maps[w];fb=np.zeros(self.full[w].nb);fb[r]=b[c]
        return fb
    def output(self,w,c,b,U,deriv=True):
        y,J=self.full[w].output(c,self.local(w,b),U,deriv)
        if not deriv:return y,None
        rows,cols=self.maps[w]
        out=np.zeros((len(U),4,self.nb),complex);out[:,:,cols]=J[:,:,rows]
        return y,out
    def initial(self,rid):
        b=np.zeros(self.nb);lstart=0;qstart=self.nl;delays=[]
        for w,m in self.res.items():
            r=rd(BASE/'raw'/f'{rid}_{w}_{self.kind}_M3_full.json')
            loc=np.asarray(r['beta'])
            b[lstart:lstart+m.nl]=loc[:m.nl]
            b[qstart:qstart+2*m.nm]=loc[m.nl:m.nl+2*m.nm]
            delays.append(loc[m.delay_slice])
            lstart+=m.nl;qstart+=2*m.nm
        b[self.delay_slice]=np.mean(delays,axis=0)
        b[self.lag_slice]=np.log(2*np.pi*2000.)
        return b
    def unpack(self,b):
        result={};offset=self.nl
        for w,m in self.res.items():
            result[w]=dict(frequency_hz=np.exp(b[offset:offset+m.nm])/2/np.pi,
                           damping_percent=100*np.exp(b[offset+m.nm:offset+2*m.nm]))
            offset+=2*m.nm
        tau=old.TAU+old.UNIT*np.einsum('oip,p->oi',old.DESIGN,b[self.delay_slice])
        return dict(windows=result,delay=old.delay_values(tau),cutoff_hz=np.exp(b[self.lag_slice])/2/np.pi)

class JointSolver(pc.ProfileIV):
    def __init__(self,m,data):
        self.m=m;self.data=data;self.method='M3'
        self.calls=0;self.inner_failures=0;self.history=[]
        self.ctx={w:m.full[w].context(a['f']) for w,a in data.items()}
    def uh_window(self,w,b):
        a=self.data[w]
        return np.linalg.solve(np.eye(4)[None]+a['K']@self.m.full[w].resp(self.ctx[w],self.m.local(w,b)),
                               a['D'][:,:,None])[:,:,0]
    def eq(self,b):
        m=self.m;A=np.zeros((m.nb,m.nb));g=np.zeros(m.nb);zz=np.zeros(m.nb);ee=0.
        for w,a in self.data.items():
            fb=m.local(w,b);u=self.uh_window(w,b);fm=m.full[w];c=self.ctx[w]
            af,gf,_,e=old.equations(fm,c,a['U'],a['Y'],u,fb)
            _,zl,zn=fm.compact(c,fb,u)
            zfull=np.r_[np.sum(abs(zl)**2,axis=0).ravel(),np.sum(abs(zn)**2,axis=(0,1))]
            rows,cols=m.maps[w]
            A[np.ix_(cols,cols)]+=af[np.ix_(rows,rows)]
            g[cols]+=gf[rows];zz[cols]+=zfull[rows];ee+=e
        score=float(np.max(abs(g)/np.sqrt(np.maximum(zz*ee,1e-300))))
        return A,g,score,float(ee),zz
    def profile(self,q,lin,maxinner=60):
        b=np.r_[lin,q];cache={}
        for w,a in self.data.items():
            _,LL,_=self.m.full[w].compact(self.ctx[w],self.m.local(w,b),a['U'],False)
            cache[w]=LL
        ok=False;history=[]
        for it in range(maxinner):
            gl=np.zeros(self.m.nl);zl2=np.zeros(self.m.nl);ee=0.;steps=[]
            for w,a in self.data.items():
                fm=self.m.full[w];u=self.uh_window(w,b)
                _,ZL,_=fm.compact(self.ctx[w],self.m.local(w,b),u,False)
                for o,(lr,gc) in enumerate(self.m.linear_blocks[w]):
                    X=cache[w][:,o,lr];Z=ZL[:,o,lr]
                    e=a['Y'][:,o]-X@b[gc]
                    gl[gc]=(Z.conj().T@e).real
                    zl2[gc]=np.sum(abs(Z)**2,axis=0);ee+=np.sum(abs(e)**2)
                    steps.append((gc,pc.scaled_solve((Z.conj().T@X).real,gl[gc])))
            score=float(np.max(abs(gl)/np.sqrt(np.maximum(zl2*ee,1e-300))))
            history.append(score)
            if score<2e-12:ok=True;break
            for gc,step in steps:b[gc]+=step
            if not np.isfinite(b).all() or (it>8 and history[-1]>history[-5]*1.5):break
        self.calls+=1
        if not ok:self.inner_failures+=1
        return b,ok,None,len(history)

def get_data(rid):
    result={}
    for w in WINDOWS:
        a=old.spectrum(rid,w);li=a['line_indices']
        result[w]={k:a[k][li] for k in ('f','U','Y','D','K')}
        result[w]['N']=a['N']
        cov=np.load(BASE/'analysis'/f'noise_{rid}_{w}.npz')
        assert np.allclose(cov['f'],a['f'],rtol=0,atol=1e-10)
        result[w]['Phi']=cov['Phi'][li]
    return result
