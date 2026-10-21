import os,sys
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
from pathlib import Path
import json,time,traceback
import numpy as np
import scipy.linalg as sla
HERE=Path(__file__).resolve().parents[1]
HIGH=HERE.parent
SRC=HIGH.parent
LOW=SRC/'followup_lowfreq'
sys.path.insert(0,str(HERE/'baseline'))
from modal_model import ModalModel,equations,iv_solve
FS=5000.;P0=30000
TAU=np.load(HERE/'baseline/tau_seconds.npy')
CTRS=('lowgain','nominal','highgain')
IDS=[g+'_'+c+('_pppp' if g=='C' else '') for g in ('A','C') for c in CTRS]
SPECS={
 'H300':dict(band=[220,420],f=[263.,282.,340.,360.],z=[.059,.023,.0045,.003],lo=[220,250,310,345],hi=[290,310,355,400],report=[1,2,3]),
 'H493':dict(band=[450,550],f=[493.],z=[.0035],lo=[460],hi=[530],report=[0])}
def rd(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def js(p,data):
    def enc(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,Path):return str(v)
        raise TypeError(type(v).__name__)
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=enc),encoding='utf-8');q.replace(p)
def ctrl(rid):return next(c for c in CTRS if c in rid)
def setup(w):
    s=SPECS[w];m=ModalModel(len(s['f']),TAU,s['band'][1]);lb,ub=m.bounds(s['lo'],s['hi']);return s,m,lb,ub
def evalK(f,ct):
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
        xx=a[name].reshape(p,P0,4)[keep].reshape(N,4)
        out[name]=np.fft.rfft(xx,axis=0)[k]
    out['K']=evalK(f,ctrl(rid));return out
def basefit(rid,w,method):
    r=rd(HIGH/'raw'/f'{rid}_{w}_{method}.json')
    q=rd(HIGH/'diagnostics'/f'{rid}_{w}_{method}.json')
    _,m,lb,ub=setup(w);b=np.asarray(q['beta'])
    margin=np.minimum(b[m.iw:]-lb[m.iw:],ub[m.iw:]-b[m.iw:])
    return dict(id=rid,window=w,method=method,beta=b,frequency_hz=q['frequency_hz'],damping_percent=q['damping_percent'],
       original_converged=r['converged'],original_strict=r['strict_score_pass'],
       converged=r['converged'],strict_score_pass=q['polished_strict_pass'],near_bound=bool(np.any(margin<1e-4)),
       score=q['polished_score'],original_score=r['score'],extra_iterations=q['extra_iterations'],
       source=str(HIGH/'diagnostics'/f'{rid}_{w}_{method}.json'))
def eligible(r):return bool(r.get('converged',False) and r.get('strict_score_pass',False) and not r.get('near_bound',True))
def scaled_solve(A,B):
    sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300));S=A*sc[:,None]*sc[None,:]
    return sc[:,None]*np.linalg.solve(S,B*sc[:,None]) if B.ndim==2 else sc*np.linalg.solve(S,B*sc)
def diagnostics(A):
    sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300));sv=np.linalg.svd(A*sc[:,None]*sc[None,:],compute_uv=False)
    return dict(condition=sv[0]/sv[-1],rank_1e12=int(np.sum(sv>sv[0]*1e-12)),dimension=len(A))
def jk_se(a):
    a=np.asarray(a);P=len(a);return np.sqrt((P-1)/P*np.sum((a-a.mean(0))**2,axis=0))
def polish(m,c,U,Y,Uh,b,lb,ub):
    b=b.copy();hist=[]
    for _ in range(25):
        AA,g,score,ee=equations(m,c,U,Y,Uh(b) if callable(Uh) else Uh,b)
        if score<1e-8:break
        step=scaled_solve(AA,g);fac=min(1.,.05/max(np.max(abs(step[m.iw:])),1e-300))
        pp=step>0;nn=step<0
        if pp.any():fac=min(fac,float(np.min(.995*(ub[pp]-b[pp])/step[pp])))
        if nn.any():fac=min(fac,float(np.min(.995*(lb[nn]-b[nn])/step[nn])))
        hist.append(dict(score=score,factor=fac))
        if fac<1e-7:break
        b+=fac*step
    score=equations(m,c,U,Y,Uh(b) if callable(Uh) else Uh,b)[2]
    return b,dict(score=score,strict_score_pass=bool(score<1e-6),polish_iterations=len(hist),polish_history=hist)
def foldfit(rid,w,a,omit,method):
    dest=HERE/'raw'/f'{rid}_{w}_omit{omit}_{method}.json'
    if dest.exists():return rd(dest)
    t=time.perf_counter();s,m,lb,ub=setup(w);r=dict(id=rid,window=w,omit=omit,method=method)
    try:
        ix=slice(None) if method=='M1' else a['line_indices']
        c=m.context(a['f'][ix]);U,Y,D,K=(a[k][ix] for k in ('U','Y','D','K'))
        Uh=U if method=='M1' else lambda bb:np.linalg.solve(np.eye(4)[None]+K@m.resp(c,bb),D[:,:,None])[:,:,0]
        b0,ranks=m.initial(c,U,Y,s['f'],s['z'])
        b,info=iv_solve(m,c,U,Y,Uh,None,b0,lb,ub)
        r.update(converged=info['success'],solver=info,initial_ranks=ranks,initial_beta=b0)
        margin=np.minimum(b[m.iw:]-lb[m.iw:],ub[m.iw:]-b[m.iw:])
        if info['success'] and min(margin)>=1e-4:
            b,pi=polish(m,c,U,Y,Uh,b,lb,ub);r.update(pi)
        else:r.update(score=info['score_max'],strict_score_pass=bool(info['score_max']<1e-6),polish_iterations=0)
        f,z=m.modal(b);margin=np.minimum(b[m.iw:]-lb[m.iw:],ub[m.iw:]-b[m.iw:])
        r.update(beta=b,frequency_hz=f,damping_percent=z,near_bound=bool(min(margin)<1e-4),iterations=info['iterations'])
    except Exception as e:r.update(exception=repr(e),traceback=traceback.format_exc(),converged=False,strict_score_pass=False,near_bound=True)
    r['seconds']=time.perf_counter()-t;js(dest,r);return r
