"""Task E: same line-only IV equations, initialisation and damped-root ablation."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ.update(NOFLEX='1',FC='60',AL='30',AH='100',FHI='100',FLO='3',SW='.03',FSPR='.01',ZSPR='.4')
import sys,json,time,argparse
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
sys.path.insert(0,str(BASE/'code'))
import numpy as np
import cl2 as c
from iv_lm import solve as lm_solve,moments
from concurrent.futures import ProcessPoolExecutor,as_completed
OUT=HERE/'short_e';OUT.mkdir(exist_ok=True)

def save(path,obj):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=True),encoding='utf-8');tmp.replace(path)

def basic_solve(beta0,U,Y,model,instrument,W):
    b=beta0.copy()
    for it in range(80):
        A,g,score=moments(b,U,Y,model,instrument,W)
        scale=1/np.sqrt(np.abs(np.diag(A))+1e-300)
        step=scale*np.linalg.solve(A*scale[:,None]*scale[None,:],g*scale)
        mx=np.max(np.abs(step[c.IW:c.ITR]))
        if mx>.05:step*=.05/mx
        b+=step
        if np.max(np.abs(step[c.IW:c.ITR]))<1e-5:break
    _,_,score=moments(b,U,Y,model,instrument,W)
    return b,dict(success=it<79 and bool(np.isfinite(b).all()),status='pole_step_converged' if it<79 else 'iteration_limit',iterations=it+1,score_max=score)

def one(P,r,strategies=('E1','E2','E3'),weighted=True):
    a=np.load(HERE.parent/'experiment_a'/'data'/f'P{P:02d}_r{r:03d}.npz')
    k,z,D,U,Y,pert=[a[x] for x in ('k','zk','D','U','Y','pert')]
    sel=k%P==0;zs,Ds,Us,Ys=z[sel],D[sel],U[sel],Y[sel];Ks=c.Ktf(zs)
    model=lambda b,u:c.out_jac(zs,b,u)
    instrument=lambda b:np.linalg.solve(np.eye(4)[None]+Ks@c.resp(zs,b),Ds[:,:,None])[:,:,0]
    truth=c.truth_beta()
    initial=truth.copy();initial[c.IW:c.IZ]+=np.log(1+.01*pert[:3]);initial[c.IZ:c.ITR]+=.4*pert[3:]
    initial[:c.NLIN]=0;initial[c.ITR:]=0
    initial=c.init_lin(zs,Us,Ys,np.ones_like(Us.real),initial)
    direct=json.loads((HERE.parent/'experiment_a'/'results'/f'P{P:02d}_r{r:03d}_dir_lines_W1.json').read_text())
    for strategy in strategies:
        b0=initial.copy() if strategy=='E1' else np.asarray(direct['beta'])
        solver=lm_solve if strategy=='E3' else basic_solve
        path=OUT/f'P{P:02d}_r{r:03d}_M3_{strategy}.json'
        if path.exists():
            record=json.loads(path.read_text());b=np.asarray(record['beta']) if record.get('beta') is not None else None
        else:
            started=time.perf_counter();record=dict(P=P,replicate=r,seed=int(a['seed']),method='M3',strategy=strategy,status='exception',success=False,beta=None,iterations=0,initial_beta=b0.tolist())
            try:
                if strategy=='E3':b,info=solver(b0,Us,Ys,model,instrument,None,c.IW,c.ITR,K=Ks,response=lambda bb:c.resp(zs,bb))
                else:b,info=solver(b0,Us,Ys,model,instrument,None)
                record.update(info,beta=b.tolist(),error_percent=(100*np.expm1(b[c.IW:c.ITR]-truth[c.IW:c.ITR])).tolist())
            except Exception as ex:record['exception']=repr(ex);b=None
            record['seconds']=time.perf_counter()-started;save(path,record)
        if P!=2 or not weighted or b is None:continue
        path=OUT/f'P{P:02d}_r{r:03d}_M3w_{strategy}.json'
        if path.exists():continue
        started=time.perf_counter();record=dict(P=P,replicate=r,seed=int(a['seed']),method='M3w',strategy=strategy,n_l=16,outer=1,status='exception',success=False,beta=None,iterations=0)
        try:
            non=np.flatnonzero(~sel);idx=np.flatnonzero(sel)
            residual=Y-c.out_jac(z,b,U,need_jac=False)[0]
            nearest=np.array([non[np.argsort(np.abs(k[non]-k[q]),kind='stable')[:16]] for q in idx])
            vv=residual[nearest]
            C=np.einsum('kli,klj->kij',vv,vv.conj())/16
            C=.9*C+.1*np.eye(4)[None]*np.real(np.diagonal(C,axis1=1,axis2=2))[:,None,:]
            W=np.linalg.inv(C)
            if strategy=='E3':bw,info=solver(b,Us,Ys,model,instrument,W,c.IW,c.ITR,K=Ks,response=lambda bb:c.resp(zs,bb))
            else:bw,info=solver(b,Us,Ys,model,instrument,W)
            record.update(info,beta=bw.tolist(),initial_beta=b.tolist(),error_percent=(100*np.expm1(bw[c.IW:c.ITR]-truth[c.IW:c.ITR])).tolist())
        except Exception as ex:record['exception']=repr(ex)
        record['seconds']=time.perf_counter()-started;save(path,record)
    return dict(P=P,r=r)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--periods',default='1,2');ap.add_argument('--workers',type=int,default=2);ap.add_argument('--count',type=int,default=40);ap.add_argument('--start',type=int,default=0);ap.add_argument('--strategies',default='E1,E2,E3');ap.add_argument('--no-weighted',action='store_true');args=ap.parse_args()
    jobs=[(p,r) for p in map(int,args.periods.split(',')) for r in range(args.start,args.start+args.count)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        fut=[pool.submit(one,p,r,tuple(args.strategies.split(',')),not args.no_weighted) for p,r in jobs]
        for f in as_completed(fut):print(json.dumps(f.result()),flush=True)
if __name__=='__main__':main()
