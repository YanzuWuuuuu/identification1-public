from core import *
for k,v in dict(NOFLEX='1',FC='60',AL='30',AH='100',SW='.03',FSPR='.01',ZSPR='.4').items():os.environ[k]=v
import cl2
from scipy import signal
from concurrent.futures import ProcessPoolExecutor,as_completed
NI=np.array([0,1,0,1]);NO=np.array([1,0,1,0])
TRUTH_TAU=(NO[:,None]+NI[None,:])/FS

class SimModel(DelayModel):
    """Original cl2 126 parameters, plus seven delays; no extra background."""
    def __init__(self):
        self.full=DelayModel(3,np.full((4,4),200e-6),100)
        self.nm=3;self.L=30;self.iw=120;self.iz=123;self.base_nb=126;self.nb=133
        self.fs=FS;self.tau=self.full.tau;self.wr=self.full.wr
        self.map=np.r_[np.concatenate([np.arange(o*38,o*38+30) for o in range(4)]),np.arange(152,165)]
    def expand(self,b):
        full=np.zeros(165);full[self.map]=b;return full
    def compact(self,c,b,U,deriv=True):
        y,L,N=self.full.compact(c,self.expand(b),U,deriv);return y,L[:,:,:30],N
    def gmodel(self,c,b,deriv=True):
        G,J=self.full.gmodel(c,self.expand(b),deriv);return G,None if J is None else J[:,:,:,self.map]
    def delay(self,b):return self.full.delay(self.expand(b))
    def truth(self):
        b=np.zeros(self.nb)
        for o in range(4):
            for j in range(3):b[o*30+10*j:o*30+10*j+4]=cl2.RES[j][o]
        b[self.iw:self.iz]=np.log(2*np.pi*cl2.FM);b[self.iz:self.base_nb]=np.log(cl2.ZM)
        dif=TRUTH_TAU-self.tau;b[self.base_nb]=dif.mean()/UNIT
        b[self.base_nb+1:self.base_nb+4]=HELM.T@dif.mean(1)/UNIT
        b[self.base_nb+4:]=HELM.T@dif.mean(0)/UNIT
        return b

def plant_with_delays():
    A,B,C=cl2.ct_ss();Ad,Bd,Cd,_,_=signal.cont2discrete((A,B,C,np.zeros((4,4))),1/FS)
    AK,BK,CK,DK=cl2.controller()
    np_=len(A);nk=len(AK);n=np_+sum(NI)+sum(NO)+nk
    ii=[];oo=[];off=np_
    for nn in NI:ii.append(np.arange(off,off+nn));off+=nn
    for nn in NO:oo.append(np.arange(off,off+nn));off+=nn
    def update(x,e):
        xp=x[:np_];xk=x[-nk:];raw=Cd@xp
        y=np.array([x[ix[-1]] if len(ix) else raw[j] for j,ix in enumerate(oo)])+e[4:8]
        u=e[:4]-CK@xk-DK@y;uw=u+e[8:]
        up=np.array([x[ix[-1]] if len(ix) else uw[j] for j,ix in enumerate(ii)])
        xn=np.zeros(n);xn[:np_]=Ad@xp+Bd@up;xn[-nk:]=AK@xk+BK@y
        for j,ix in enumerate(ii):
            if len(ix):xn[ix]=np.r_[uw[j],x[ix[:-1]]]
        for j,ix in enumerate(oo):
            if len(ix):xn[ix]=np.r_[raw[j],x[ix[:-1]]]
        return xn,np.r_[y,u]
    AA=np.zeros((n,n));BB=np.zeros((n,12));CC=np.zeros((8,n));DD=np.zeros((8,12))
    for i in range(n):
        x=np.eye(n)[i];AA[:,i],CC[:,i]=update(x,np.zeros(12))
    for i in range(12):BB[:,i],DD[:,i]=update(np.zeros(n),np.eye(12)[i])
    return AA,BB,CC,DD

def generate(seed,noise=True):
    rng=np.random.default_rng(seed);N=300000;P=10
    kk=np.arange(1,P0//2);z=np.exp(2j*np.pi*kk/P0)
    dp=np.zeros((P0//2+1,4),complex)
    dp[kk]=.2*np.sqrt(P0*cl2.dshape(z))[:,None]*np.exp(2j*np.pi*rng.random((len(kk),4)))
    d0=np.fft.irfft(dp,n=P0,axis=0);d=np.tile(d0,(P+2,1));nn=len(d)
    v=rng.normal(size=(nn,4))*3e-7*noise;w=rng.normal(size=(nn,4))*.03*noise
    AA,BB,CC,DD=plant_with_delays();radius=max(abs(np.linalg.eigvals(AA)))
    if radius>=1:raise ValueError(('unstable closed loop',radius))
    _,out,_=signal.dlsim((AA,BB,CC,DD,1/FS),np.c_[d,v,w])
    d=d[2*P0:];out=out[2*P0:]
    k=np.arange(180,6001);f=k*FS/N
    a=dict(k=k,f=f,N=N,P=P,line_indices=np.flatnonzero(k%P==0),K=cl2.Ktf(np.exp(2j*np.pi*f/FS)))
    for name,x in [('D',d),('Y',out[:,:4]),('U',out[:,4:])]:a[name]=np.fft.rfft(x,axis=0)[k]
    return a,radius

def one(i):
    path=ROOT/'raw'/f'sim_{i:02d}.json'
    if path.exists():return rd(path)
    t=time.perf_counter();r=dict(replicate=i,seed=1002+i)
    try:
        a,radius=generate(1002+i);m=SimModel();bt=m.truth();lb,ub=m.bounds([8,9,13],[13,15,20])
        r.update(spectral_radius=radius,truth_delay=delay_values(TRUTH_TAU),methods={})
        for method in ['M1','M3']:
            c,U,Y,Uh=inputs(m,a,method);b0=bt.copy();b0[m.base_nb:]=0
            rng=np.random.default_rng(100200+i);b0[m.iw:m.iz]+=rng.normal(size=3)*.01
            b0[m.iz:m.base_nb]+=rng.normal(size=3)*.4
            b0=linsolve(m,c,U,Y,Uh,b0)
            grid=[]
            for us in range(-200,201,10):
                b=b0.copy();b[m.base_nb]=us/100;b=linsolve(m,c,U,Y,Uh,b)
                ee=equations(m,c,U,Y,Uh(b) if callable(Uh) else Uh,b)[3];grid.append((ee,us,b))
            initial=min(grid,key=lambda x:x[0]);b,info=solve(m,c,U,Y,Uh,initial[2],lb,ub,maxit=160)
            r['methods'][method]=dict(beta=b,solver=info,grid_offset_us=initial[1],delay=delay_values(m.delay(b)),
                delay_error_us=(m.delay(b)-TRUTH_TAU)*1e6,frequency_hz=m.modal(b)[0],damping_percent=m.modal(b)[1],
                modal_relative_error_percent=100*np.expm1(b[m.iw:m.base_nb]-bt[m.iw:m.base_nb]))
    except Exception as e:r.update(exception=repr(e),traceback=traceback.format_exc())
    r['seconds']=time.perf_counter()-t;js(path,r);print('SIM',i,{k:v['solver']['success'] for k,v in r.get('methods',{}).items()},round(r['seconds'],1),flush=True);return r

def validate():
    a,radius=generate(19001,False);m=SimModel();bt=m.truth();c,U,Y,Uh=inputs(m,a,'M3')
    rel=np.linalg.norm(Y-m.output(c,bt,U,False)[0])/np.linalg.norm(Y)
    f=np.linspace(3,100,100);co=m.context(f);G0=m.resp(co,bt)
    old=cl2.resp(np.exp(2j*np.pi*f/FS),cl2.truth_beta())*np.exp(-2j*np.pi*f[:,None,None]*TRUTH_TAU)
    b0=bt.copy();b0[m.base_nb:]+=.001
    lb,ub=m.bounds([8,9,13],[13,15,20]);b,info=solve(m,c,U,Y,Uh,b0,lb,ub,maxit=160)
    r=dict(radius=radius,noiseless_relative_residual=rel,response_vs_cl2=np.linalg.norm(G0-old)/np.linalg.norm(old),
           noiseless_solver=info,noiseless_delay_error_us=(m.delay(b)-TRUTH_TAU)*1e6,
           true_delay=delay_values(TRUTH_TAU),n=300000,P=10,settle=2,fs=FS,
           boundary_note='Integer input/output delay states simulated exactly. Modal transient family is retained; finite noise boundary terms from delay registers are not extra fit parameters.')
    js(ROOT/'analysis/simulation_validation.json',r);print('SIMVALID',r['noiseless_relative_residual'],r['response_vs_cl2'],info['success'])
    assert r['response_vs_cl2']<1e-10 and rel<1e-7
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--validate',action='store_true');p.add_argument('--count',type=int,default=40);p.add_argument('--workers',type=int,default=4)
    args=p.parse_args()
    if args.validate:validate()
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for f in as_completed([pool.submit(one,i) for i in range(args.count)]):f.result()
