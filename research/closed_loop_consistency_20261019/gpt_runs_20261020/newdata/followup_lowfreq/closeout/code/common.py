"""Read-only access to original results; all computations use copied model files."""
import os, sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
sys.dont_write_bytecode=True
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import json, time, traceback
import numpy as np
import scipy.linalg as sla
HERE=Path(__file__).resolve().parents[1]
BASE=HERE.parent
SRC=BASE.parent
sys.path.insert(0,str(HERE/'baseline'))
import realfit_extended as m
import weighted_core as wc
FS=5000.; P0=30000
TAU=np.load(HERE/'baseline/tau_seconds.npy')
CTRS=('lowgain','nominal','highgain')
IDS=tuple('A_'+c for c in CTRS)+tuple('C_'+c+'_pppp' for c in CTRS)
MOD=np.arange(m.IW,m.ITR)
def js(path,data):
    def enc(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=enc),encoding='utf-8')
    tmp.replace(path)
def rd(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def fit(rid,method):return rd(SRC/'fits'/f'{rid}_{method}.json')
def meta(rid):return rd(SRC/'prepared'/f'{rid}.json')
def evalK(f,ctrl):
    a=np.load(SRC/'controllers'/f'{ctrl}.npz');A,B,T=(a[q] for q in ('A','B','T'))
    zi=np.exp(-2j*np.pi*np.asarray(f)/FS)
    aa=np.eye(4)[None]-sum(A[i][None]*zi[:,None,None]**(i+1) for i in range(8))
    bb=sum(B[j][None]*zi[:,None,None]**j for j in range(9))
    return -T[None]@np.linalg.solve(aa,bb)
def spectrum(rid,drop=0,omit=None,extended=False):
    a=np.load(SRC/'prepared'/f'{rid}.npz');p=int(a['P'])
    keep=np.arange(p)[drop:]
    if omit is not None:keep=keep[keep!=omit]
    P=len(keep);N=P*P0;lo,hi=(1,27) if extended else (3,25)
    k=np.arange(int(np.ceil(lo*N/FS)),int(np.floor(hi*N/FS))+1);f=k*FS/N
    out=dict(k=k,f=f,P=P,N=N,line_indices=np.flatnonzero(k%P==0),kept_periods=keep)
    for name in ('U','Y','D'):
        x=a[name].reshape(p,P0,4)[keep].reshape(N,4)
        out[name]=np.fft.rfft(x,axis=0)[k]
    out['K']=evalK(f,meta(rid)['controller'])
    return out
def oldspec(rid):
    with np.load(SRC/'spectra'/f'{rid}.npz') as a:return {k:a[k] for k in a.files}
def interp_cov(f,fp,C):
    out=np.empty((len(f),4,4),complex)
    for i in range(4):
        for j in range(4):
            out[:,i,j]=np.interp(f,fp,C[:,i,j].real)+1j*np.interp(f,fp,C[:,i,j].imag)
    return (out+out.conj().transpose(0,2,1))/2
def g_jac(c,b):
    # Difference removes transient terms, exactly as original gmodel_from_output.
    z=np.zeros((len(c.f),4)); y0,j0=m.out_jac(c,b,z)
    G=np.empty((len(c.f),4,4),complex);Gj=np.empty((len(c.f),4,4,m.NB),complex)
    for j in range(4):
        u=z.copy();u[:,j]=1
        y,jj=m.out_jac(c,b,u);G[:,:,j]=y-y0;Gj[:,:,j,:]=jj-j0
    return G,Gj
def scaled_solve(A,B):
    sc=1/np.sqrt(np.maximum(np.abs(np.diag(A)),1e-300))
    return sc[:,None]*np.linalg.solve(A*sc[:,None]*sc[None,:],B*sc[:,None]) if B.ndim==2 else sc*np.linalg.solve(A*sc[:,None]*sc[None,:],B*sc)
def diagnostics(A):
    sc=1/np.sqrt(np.maximum(np.abs(np.diag(A)),1e-300))
    sv=np.linalg.svd(A*sc[:,None]*sc[None,:],compute_uv=False)
    return dict(condition=sv[0]/sv[-1],rank_1e12=int(np.sum(sv>sv[0]*1e-12)),dimension=len(A))
def realvec(a):return np.concatenate((a.real,a.imag),axis=1)
def realmat(C):
    return np.concatenate((np.concatenate((C.real,-C.imag),axis=2),
                           np.concatenate((C.imag,C.real),axis=2)),axis=1)
def bands(f):
    return [(float(lo),float(min(lo+2,25)),(f>=lo)&((f<lo+2) if lo<23 else (f<=25))) for lo in np.arange(3,25,2)]
def rms_record(rid):
    a=oldspec(rid);N=int(a['N'])
    return np.sqrt(2*np.sum(abs(a['Y'])**2,axis=0))/N
