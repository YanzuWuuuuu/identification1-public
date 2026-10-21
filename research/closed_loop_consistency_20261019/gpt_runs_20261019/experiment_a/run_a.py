"""Paired Monte Carlo extension. Imports original cl2 read-only; no theory changes."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
os.environ.update(NOFLEX='1',FC='60',AL='30',AH='100',FHI='100',FLO='3',SW='0.03',FSPR='0.01',ZSPR='0.4')
import sys, json, time, hashlib, argparse
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
sys.path.insert(0,str(BASE/'code'))
import numpy as np
import cl2 as c
from concurrent.futures import ProcessPoolExecutor, as_completed

METHODS=['dir_band_W1','dir_lines_W1','ref_band_W1','refD_band_W1']
GRID={}
def atomic_json(path,obj):
    temp=path.with_suffix('.tmp'); temp.write_text(json.dumps(obj,indent=2,allow_nan=True),encoding='utf-8'); temp.replace(path)

def grid(P):
    if P not in GRID:
        N=30000*P; k=np.arange(int(np.ceil(3*N/c.fs)),int(np.floor(100*N/c.fs))+1)
        f=k*c.fs/N; z=np.exp(2j*np.pi*k/N)
        GRID[P]=(k,f,z,c.Ktf(z))
    return GRID[P]

def score(z,U,Y,Uh,b):
    uh=Uh(b) if callable(Uh) else Uh
    g=np.zeros(c.NB); zn=np.zeros(c.NB); en=0.
    for start in range(0,len(z),1200):
        sl=slice(start,min(start+1200,len(z)))
        yh,J=c.out_jac(z[sl],b,U[sl]); Z=J if uh is U else c.out_jac(z[sl],b,uh[sl])[1]
        e=Y[sl]-yh
        g+=np.real(np.einsum('kop,ko->p',Z.conj(),e))
        zn+=np.sum(np.abs(Z)**2,axis=(0,1)); en+=np.sum(np.abs(e)**2)
    return float(np.max(np.abs(g)/np.maximum(np.sqrt(zn*en),1e-300))),float(np.sqrt(en/np.sum(np.abs(Y)**2)))

def one(P,r):
    t0=time.perf_counter(); seed=202610190+1000*P+r
    root=HERE/'results'; root.mkdir(exist_ok=True); datadir=HERE/'data'; datadir.mkdir(exist_ok=True)
    k,f,z,Kz=grid(P); bt=c.truth_beta(); rng=np.random.default_rng(seed)
    datafile=datadir/f'P{P:02d}_r{r:03d}.npz'
    if datafile.exists():
        a=np.load(datafile); D,U,Y,pert=[a[x] for x in ('D','U','Y','pert')]
    else:
        d,u,y=c.simulate_periodic(30000,P,.2,3e-7,rng,settle=2)
        D=np.fft.fft(d,axis=0)[k]; U=np.fft.fft(u,axis=0)[k]; Y=np.fft.fft(y,axis=0)[k]
        pert=rng.standard_normal(2*c.NM)
        kp=np.arange(18,601); fp=kp*c.fs/30000; zp=np.exp(2j*np.pi*kp/30000)
        Dper=np.fft.fft(d.reshape(P,30000,4),axis=1)[:,kp]
        Uper=np.fft.fft(u.reshape(P,30000,4),axis=1)[:,kp]
        Yper=np.fft.fft(y.reshape(P,30000,4),axis=1)[:,kp]
        tmp=datafile.with_suffix('.tmp.npz')
        np.savez_compressed(tmp,D=D,U=U,Y=Y,k=k,f=f,zk=z,pert=pert,Dper=Dper,Uper=Uper,Yper=Yper,fper=fp,zkper=zp,seed=seed,P=P)
        tmp.replace(datafile)
    for method in METHODS:
        dest=root/f'P{P:02d}_r{r:03d}_{method}.json'
        if dest.exists(): continue
        start=time.perf_counter()
        result=dict(P=P,replicate=r,seed=seed,method=method,status='exception',iterations=None,beta=None,log_error=None,error_percent=None,score_max=None,residual=None)
        try:
            sel=k%P==0 if method=='dir_lines_W1' else np.ones(len(k),bool)
            zs,Us,Ys,Ds,Ks=z[sel],U[sel],Y[sel],D[sel],Kz[sel]
            W=np.ones_like(Us.real)
            b0=bt.copy(); b0[c.IW:c.IZ]+=np.log(1+.01*pert[:c.NM]); b0[c.IZ:c.ITR]+=.4*pert[c.NM:]
            b0[:c.NLIN]=0; b0[c.ITR:]=0
            b0=c.init_lin(zs,Us,Ys,W,b0)
            if method.startswith('dir'): Uh=Us
            elif method.startswith('refD'): Uh=Ds
            else:
                def Uh(b):
                    return np.linalg.solve(np.eye(4)[None]+Ks@c.resp(zs,b),Ds[:,:,None])[:,:,0]
            b,it=c.iv_solve(zs,Us,Ys,Uh,W,b0,it=80,tol=1e-5)
            finite=bool(np.all(np.isfinite(b)))
            sc,rr=score(zs,Us,Ys,Uh,b) if finite else (np.nan,np.nan)
            result.update(status='pole_step_converged' if finite and it<79 else ('iteration_limit' if finite else 'nonfinite'),iterations=int(it+1),beta=b.tolist(),initial_beta=b0.tolist(),log_error=(b[c.IW:c.ITR]-bt[c.IW:c.ITR]).tolist(),error_percent=(100*np.expm1(b[c.IW:c.ITR]-bt[c.IW:c.ITR])).tolist(),score_max=sc,residual=rr)
        except Exception as exc:
            result['exception']=repr(exc)
        result['elapsed_seconds']=time.perf_counter()-start
        atomic_json(dest,result)
    return dict(P=P,r=r,seconds=time.perf_counter()-t0)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--workers',type=int,default=3); ap.add_argument('--periods',default='1,2,4,8,16'); ap.add_argument('--count',type=int)
    args=ap.parse_args(); jobs=[(p,r) for p in map(int,args.periods.split(',')) for r in range(args.count or (20 if p==16 else 40))]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(one,p,r) for p,r in jobs]
        for task in as_completed(futures): print(json.dumps(task.result()),flush=True)
if __name__=='__main__': main()
