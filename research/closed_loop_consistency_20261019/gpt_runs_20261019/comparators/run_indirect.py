"""Four independent single-axis experiments; full propagated FRF covariance.
The raw covariance cannot be inverted for requested P<=16. M6_reg is distinct.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
for key,val in dict(NOFLEX='1',FC='60',AL='30',AH='100',FHI='100',SW='.03').items():os.environ[key]=val
import sys,time,json,pathlib,numpy as np,scipy.signal as ss,scipy.linalg as sla
sys.dont_write_bytecode=True
from scipy.optimize import least_squares
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parents[1]/'code'))
import cl2 as c
from sml_core import gmodel_from_output
from run_sml import atomic_npz

_CL=c.closed_loop(); ev,V=sla.eig(_CL[0]);VB=sla.solve(V,_CL[1]);CV=_CL[2]@V
def fast_sim(signals):
    drive=signals@VB.T
    for q in range(len(ev)):drive[:,q]=ss.lfilter([0.,1.],[1.,-ev[q]],drive[:,q])
    return np.real(drive@CV.T)+signals@_CL[3].T

def simulate_single(P0,p,axis,seed,settle=2):
    rng=np.random.default_rng(seed);kk=np.arange(1,P0//2);zs=np.exp(2j*np.pi*kk/P0)
    amp=.2*np.sqrt(P0*c.dshape(zs));dl=np.zeros(P0,complex)
    dl[kk]=amp*np.exp(2j*np.pi*rng.random(len(kk)));dl[P0-kk]=dl[kk].conj()
    d=np.zeros(((p+settle)*P0,4));d[:,axis]=np.tile(np.fft.ifft(dl).real,p+settle)
    v=3e-7*rng.standard_normal(d.shape);w=.03*rng.standard_normal(d.shape)
    out=fast_sim(np.hstack([d,v,w]));cut=settle*P0
    return d[cut:],out[cut:,4:8],out[cut:,:4]

def make_frf(P,r):
    p=P//4;P0=30000;k=np.arange(18,601);zs=np.exp(2j*np.pi*k/P0)
    Ups=[];Yps=[];Dps=[]
    for a in range(4):
        d,u,y=simulate_single(P0,p,a,202610190+P*1000+r+700000+100000*a)
        Dps.append(np.fft.fft(d.reshape(p,P0,4),axis=1)[:,k,:])
        Ups.append(np.fft.fft(u.reshape(p,P0,4),axis=1)[:,k,:])
        Yps.append(np.fft.fft(y.reshape(p,P0,4),axis=1)[:,k,:])
    Uper=np.array(Ups);Yper=np.array(Yps);Dper=np.array(Dps)
    # Matrices columns are independent experiments, rows are physical channels.
    U=Uper.mean(1).transpose(1,2,0);Y=Yper.mean(1).transpose(1,2,0);D=Dper.mean(1).transpose(1,2,0)
    tud=U@np.linalg.inv(D);tyd=Y@np.linalg.inv(D);G=tyd@np.linalg.inv(tud)
    out=dict(G=G,Uper=Uper,Yper=Yper,Dper=Dper,zk=zs,f=k*c.fs/P0,
             condition_U=np.linalg.cond(U),condition_Tud=np.linalg.cond(tud))
    if p>1:
        C=np.zeros((len(k),16,16),complex);Uinv=np.linalg.inv(U)
        for a in range(4):
            uc=Uper[a]-Uper[a].mean(0);yc=Yper[a]-Yper[a].mean(0)
            e=yc-np.einsum('kij,pkj->pki',G,uc)
            Ce=np.einsum('pki,pkj->kij',e,e.conj())/(p*(p-1))
            L=np.einsum('kj,oi->kjoi',Uinv[:,a,:],np.eye(4)).reshape(len(k),16,4)
            C+=L@Ce@L.conj().transpose(0,2,1)
        out['cov_frf']=C
    return out

def fit_frf(data,pert,shrink=None):
    zs=data['zk'];G=data['G'];bt=c.truth_beta();b0=bt.copy()
    b0[c.IW:c.IZ]+=np.log(1+.01*pert[:c.NM]);b0[c.IZ:c.ITR]+=.4*pert[c.NM:];b0[:c.NLIN]=0;b0[c.ITR:]=0
    h0,h1,_=c.blocks(zs,b0[c.IW:c.IZ],b0[c.IZ:c.ITR]);H=np.stack([h for pair in zip(h0,h1) for h in pair],axis=1)
    coeff=np.linalg.lstsq(np.r_[H.real,H.imag],np.r_[G.real.reshape(len(zs),16),G.imag.reshape(len(zs),16)],rcond=None)[0]
    b0[:c.NLIN]=coeff.reshape(-1)
    model=lambda b,U:c.out_jac(zs,b,U)
    gm=gmodel_from_output(model,c.NB,len(zs))
    Li=None
    if shrink is not None:
        C=data['cov_frf'];di=np.real(np.diagonal(C,axis1=1,axis2=2))
        C=(1-shrink)*C+shrink*np.eye(16)[None]*di[:,None,:]
        Li=np.linalg.inv(np.linalg.cholesky((C+C.conj().transpose(0,2,1))/2))
    cache=None
    def calc(x):
        nonlocal cache
        if cache is not None and np.array_equal(cache[0],x):return cache[1:]
        b=np.r_[x,np.zeros(c.NB-c.ITR)];gh,jh=gm(b)
        e=(gh-G).transpose(0,2,1).reshape(len(zs),16);J=jh[:,:,:,:c.ITR].transpose(0,2,1,3).reshape(len(zs),16,c.ITR)
        if Li is not None:e=np.einsum('kij,kj->ki',Li,e);J=np.einsum('kij,kjm->kim',Li,J)
        rr=np.r_[e.real.ravel(),e.imag.ravel()];jj=np.r_[J.real.reshape(-1,c.ITR),J.imag.reshape(-1,c.ITR)]
        cache=(x.copy(),rr,jj);return rr,jj
    t=time.perf_counter()
    try:
        sol=least_squares(lambda x:calc(x)[0],b0[:c.ITR],jac=lambda x:calc(x)[1],x_scale='jac',max_nfev=100,
                           ftol=1e-9,xtol=1e-9,gtol=1e-12)
        b=np.r_[sol.x,np.zeros(c.NB-c.ITR)];rr,jj=calc(sol.x)
        score=np.max(np.abs(jj.T@rr)/(np.linalg.norm(jj,axis=0)*max(np.linalg.norm(rr),1e-30)+1e-30))
        info=dict(success=bool(sol.success),status=str(sol.message),nfev=int(sol.nfev),njev=int(sol.njev or 0),cost=float(2*sol.cost),normalized_score=float(score),seconds=time.perf_counter()-t,shrink=shrink)
    except Exception as e:b=np.full(c.NB,np.nan);info=dict(success=False,status=repr(e),seconds=time.perf_counter()-t)
    return b,b0,info

def one(P,r):
    od=HERE/'indirect_results';od.mkdir(exist_ok=True);dd=HERE/'indirect_data';dd.mkdir(exist_ok=True)
    dp=dd/f'P{P:02d}_r{r:03d}.npz'
    if dp.exists():data=dict(np.load(dp))
    else:
        t=time.perf_counter();data=make_frf(P,r);data['simulation_seconds']=time.perf_counter()-t;atomic_npz(dp,**data)
    ap=HERE.parent/'experiment_a'/'data'/dp.name
    if ap.exists():pert=np.load(ap)['pert']
    else:pert=np.random.default_rng(202610190+P*1000+r).standard_normal(2*c.NM)
    for name,shr in [('M6u',None)]+([('M6_reg',.1)] if P>4 else []):
        path=od/f'P{P:02d}_r{r:03d}_{name}.npz'
        if path.exists():continue
        b,b0,info=fit_frf(data,pert,shr);bt=c.truth_beta();info.update(P=P,r=r,method=name,periods_per_experiment=P//4,
              covariance_rank_upper_bound=P-4,simulation_seconds=float(data['simulation_seconds']))
        atomic_npz(path,beta=b,beta0=b0,truth=bt,relative_error_percent=100*np.expm1(b[c.IW:c.ITR]-bt[c.IW:c.ITR]),info=json.dumps(info))
        print(path.name,info,flush=True)

if __name__=='__main__':
    worker=int(sys.argv[1]) if len(sys.argv)>1 else 0;nworkers=int(sys.argv[2]) if len(sys.argv)>2 else 1
    for P in [4,8,16]:
        for r in range(20 if P==16 else 40):
            if r%nworkers==worker:one(P,r)
