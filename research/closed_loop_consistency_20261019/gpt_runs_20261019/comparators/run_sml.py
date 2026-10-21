import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']: os.environ[key]='1'
for key,val in dict(NOFLEX='1',FC='60',AL='30',AH='100',FHI='100',SW='.03',FSPR='.01',ZSPR='.4').items():os.environ[key]=val
import sys,time,json,pathlib,numpy as np
sys.dont_write_bytecode=True
HERE=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'code'))
import cl2 as c
from sml_core import solve_sml,gmodel_from_output

def atomic_npz(path,**kw):
    tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**kw);tmp.replace(path)

def one(path,shrink=None,suffix=''):
    a=np.load(path);P=len(a['Uper']);r=int(path.stem.split('r')[1]);zs=a['zkper'];pert=a['pert']
    shrink=(.1 if P<=4 else 0.) if shrink is None else shrink
    name='M5_reg' if shrink else 'M5'; out=HERE/'sml_results';out.mkdir(exist_ok=True)
    target=out/f'{path.stem}_{name}{suffix}.npz'
    if target.exists():return
    bt=c.truth_beta();b0=bt.copy();b0[c.IW:c.IZ]+=np.log(1+.01*pert[:c.NM]);b0[c.IZ:c.ITR]+=.4*pert[c.NM:];b0[:c.NLIN]=0;b0[c.ITR:]=0
    t=time.perf_counter()
    try:
        um=a['Uper'].mean(0);ym=a['Yper'].mean(0)
        b0=c.init_lin(zs,um,ym,np.ones_like(um.real),b0)
        model=lambda b,U:c.out_jac(zs,b,U)
        gm=gmodel_from_output(model,c.NB,len(zs))
        b,info=solve_sml(a['Uper'],a['Yper'],b0,model,gm,maxiter=120,shrink=shrink)
    except Exception as e:
        b=np.full(c.NB,np.nan);info=dict(success=False,status=repr(e))
    info.update(P=P,r=r,method=name,shrink=shrink,elapsed_total=time.perf_counter()-t,seed=int(a['seed']) if 'seed'in a else 202610190+P*1000+r)
    atomic_npz(target,beta=b,beta0=b0,truth=bt,relative_error_percent=100*np.expm1(b[c.IW:c.ITR]-bt[c.IW:c.ITR]),info=json.dumps(info))
    print(target.name,info,flush=True)

if __name__=='__main__':
    worker=int(sys.argv[1]) if len(sys.argv)>1 else 0;nworkers=int(sys.argv[2]) if len(sys.argv)>2 else 1
    jobs=[(P,r) for P in [2,4,8,16] for r in range(20 if P==16 else 40) if r%nworkers==worker]
    dat=HERE.parent/'experiment_a'/'data'
    while jobs:
        done=[]
        for P,r in jobs:
            p=dat/f'P{P:02d}_r{r:03d}.npz'
            if p.exists():one(p);done.append((P,r))
        jobs=[j for j in jobs if j not in done]
        if jobs:time.sleep(10)
