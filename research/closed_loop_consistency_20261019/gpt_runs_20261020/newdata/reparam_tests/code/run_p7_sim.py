from profile_core import *
os.environ.update(NOFLEX='1',FC='60',AL='30',AH='100',FHI='100',FLO='3',SW='0.03',FSPR='0.01',ZSPR='0.4')
RESEARCH=SRC.parents[1]
sys.path.insert(0,str(RESEARCH/'code'))
import cl2 as cl
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,traceback
AOLD=RESEARCH/'gpt_runs_20261019/experiment_a'

class SimModel(ResidueModel):
    def __init__(self):
        super().__init__('L','a',joint=False)
        self.full=old.ModalModel(3,np.zeros((4,4)),100.)
        self.nl=120;self.nq=6;self.nb=126;self.iw=120;self.iz=123
        self.lb=np.r_[np.log(2*np.pi*cl.FM*.7),np.log(np.full(3,1e-5))]
        self.ub=np.r_[np.log(2*np.pi*cl.FM*1.3),np.log(np.full(3,.3))]
    def expand(self,b,deriv=True):
        # Reuse the exact same modal ZOH/transient implementation, with D0/D2 fixed at zero.
        temp=np.r_[b[:120],np.zeros(32),b[120:]]
        oldnl,oldnb=self.nl,self.nb;self.nl=152;self.nb=158
        try:fb,T=super().expand(temp,deriv)
        finally:self.nl=oldnl;self.nb=oldnb
        if deriv:T=T[:,np.r_[0:120,152:158]]
        return fb,T
    def from_cl(self,b):
        return np.r_[b[:cl.NLIN],b[cl.ITR:],b[cl.IW:cl.ITR]]
    def to_cl(self,b):
        return np.r_[b[:96],b[120:],b[96:120]]

def dataset(P,r):
    source=AOLD/'data'/f'P{P:02d}_r{r:03d}.npz'
    target=ROOT/'simdata'/source.name
    if source.exists():return np.load(source),str(source)
    if target.exists():return np.load(target),str(target)
    target.parent.mkdir(exist_ok=True)
    seed=202610190+1000*P+r;rng=np.random.default_rng(seed)
    d,u,y=cl.simulate_periodic(30000,P,.2,3e-7,rng,settle=2)
    pert=rng.standard_normal(6);N=30000*P
    k=np.arange(int(3*N/cl.fs),int(100*N/cl.fs)+1);f=k*cl.fs/N
    kp=np.arange(18,601);fp=kp*cl.fs/30000
    arrays=dict(k=k,f=f,zk=np.exp(2j*np.pi*k/N),pert=pert,seed=seed,P=P,
                fper=fp,zkper=np.exp(2j*np.pi*kp/30000))
    for key,xx in zip(('D','U','Y'),(d,u,y)):
        arrays[key]=np.fft.fft(xx,axis=0)[k]
        arrays[key+'per']=np.fft.fft(xx.reshape(P,30000,4),axis=1)[:,kp]
    np.savez_compressed(target,**arrays)
    return np.load(target),str(target)

def pack(a,P,omit=None,full=False):
    if full:
        f=a['f'];k=a['k'];out={s:a[s] for s in ('D','U','Y')}
        out.update(f=f,k=k,line_indices=np.arange(len(f)),K=cl.Ktf(a['zk']))
    else:
        keep=[i for i in range(P) if i!=omit];f=a['fper']
        out={s:a[s+'per'][keep].sum(0) for s in ('D','U','Y')}
        out.update(f=f,line_indices=np.arange(len(f)),K=np.zeros((len(f),4,4),complex))
    return out

def solve_sim(a,P,method,b0,omit=None,full=False):
    m=SimModel();data=pack(a,P,omit,full);pr=ProfileIV(m,data,method)
    b,info=pr.solve(b0,maxouter=80,tol=2e-8)
    return dict(beta=b,solver=info,success=info['success'],frequency_hz=m.modal(b)[0],damping_percent=m.modal(b)[1],
                log_parameters=b[m.nl:],error_percent=100*np.expm1(b[m.nl:]-np.r_[np.log(2*np.pi*cl.FM),np.log(cl.ZM)]))

def run(P,r):
    dest=ROOT/'raw'/f'p7_sim_P{P:02d}_r{r:03d}.json'
    if dest.exists():return
    t=time.perf_counter();out=dict(P=P,replicate=r,seed=202610190+1000*P+r,status='exception',folds=[])
    try:
        a,source=dataset(P,r);m=SimModel();bt=cl.truth_beta();b0=m.from_cl(bt)
        b0[120:123]+=np.log(1+.01*a['pert'][:3]);b0[123:]+=.4*a['pert'][3:];b0[:120]=0.
        initial={}
        for method,oldname in [('M2','dir_lines_W1'),('M3','ref_band_W1')]:
            pp=AOLD/'results'/f'P{P:02d}_r{r:03d}_{oldname}.json'
            previous=rd(pp) if pp.exists() else None
            initial[method]=m.from_cl(np.asarray(previous['beta'])) if previous and previous.get('beta') else b0.copy()
            out[method+'_previous_status']=previous.get('status') if previous else 'not previously run'
        d=solve_sim(a,P,'M2',initial['M2'])
        # Reference uses all A frequencies, as the original reference baseline.
        # ProfileIV M3 chooses its line_indices: here explicitly all frequencies.
        ref=solve_sim(a,P,'M3',initial['M3'],full=True)
        out.update(data_source=source,M2=d,M3=ref)
        js(dest.with_name(dest.stem+'_progress.json'),out)
        for i in range(P):
            fold=solve_sim(a,P,'M2',np.asarray(d['beta']),omit=i)
            fold['omit']=i;out['folds'].append(fold)
            js(dest.with_name(dest.stem+'_progress.json'),out)
        lp=P*np.asarray(d['log_parameters'])-(P-1)*np.mean([x['log_parameters'] for x in out['folds']],axis=0)
        jk=dict(log_parameters=lp,error_percent=100*np.expm1(lp-np.r_[np.log(2*np.pi*cl.FM),np.log(cl.ZM)]),
                success=bool(d['success'] and all(x['success'] for x in out['folds'])))
        out.update(JK=jk,status='complete')
    except Exception as e:out.update(exception=repr(e),traceback=traceback.format_exc())
    out['seconds']=time.perf_counter()-t;js(dest,out)
    print('P7',P,r,out['status'],out.get('JK',{}).get('success'),round(out['seconds'],1),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--periods',default='4,8,16')
    ap.add_argument('--count',type=int,default=40);args=ap.parse_args()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for fu in as_completed([pool.submit(run,P,r) for P in map(int,args.periods.split(',')) for r in range(args.count)]):fu.result()
