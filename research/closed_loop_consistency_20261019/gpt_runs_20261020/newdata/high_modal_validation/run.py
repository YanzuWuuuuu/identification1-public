import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys,json,time,hashlib,traceback,argparse
sys.dont_write_bytecode=True
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent;PARENT=HERE.parent
sys.path.insert(0,str(PARENT));sys.path.insert(0,str(HERE/'baseline'));sys.path.insert(0,str(HERE))
from modal_model import ModalModel,equations,iv_solve
import weighted_core as wc
from sml_core_copy import SMLCovResidual
from prepare_data import evalK,FS,P0
TAU=np.load(PARENT/'baseline/tau_seconds.npy')
SPECS={
 'H100':dict(band=[90,190],f=[129.,136.],z=[.007,.011],lo=[100,100],hi=[180,180],background=[False,False]),
 'H300':dict(band=[220,420],f=[263.,282.,340.,360.],z=[.059,.023,.0045,.003],lo=[220,250,310,345],hi=[290,310,355,400],background=[True,False,False,False]),
 'H493':dict(band=[450,550],f=[493.],z=[.0035],lo=[460],hi=[530],background=[False]),
 'H700':dict(band=[620,800],f=[669.,679.,682.,706.,778.],z=[.0013,.001,.0013,.0018,.0005],lo=[640,672,678,697,745],hi=[675,684,697,735,799],background=[False]*5)}
METHODS=('M1','M2','M3','M3w','M5n')
IDS=[g+'_'+c+('_pppp' if g=='C' else '') for g in ('A','C') for c in ('lowgain','nominal','highgain')]
CONFIG_HASH=hashlib.sha256(json.dumps(SPECS,sort_keys=True).encode()).hexdigest()

def js(p,a):
    t=p.with_suffix('.tmp.json');t.write_text(json.dumps(a,ensure_ascii=False,indent=2),encoding='utf-8');os.replace(t,p)
def get_data(rid,window):
    path=HERE/'spectra'/f'{rid}_{window}.npz'
    if path.exists():return np.load(path)
    a=np.load(PARENT/'prepared'/f'{rid}.npz');p=int(a['P']);N=p*P0;lo,hi=SPECS[window]['band']
    k=np.arange(int(np.ceil(lo*N/FS)),int(np.floor(hi*N/FS))+1);f=k*FS/N;li=np.flatnonzero(k%p==0)
    ctrl=next(c for c in ('lowgain','nominal','highgain') if c in rid);co=np.load(PARENT/'controllers'/f'{ctrl}.npz')
    d=dict(k=k,f=f,line_indices=li,N=N,P=p,K=evalK(f,co['A'],co['B'],co['T']))
    for x in ('U','Y','D'):d[x]=np.fft.rfft(a[x],axis=0)[k]
    np.savez_compressed(path,**d);return np.load(path)

def setup(window):
    s=SPECS[window];m=ModalModel(len(s['f']),TAU,s['band'][1]);lb,ub=m.bounds(s['lo'],s['hi']);return s,m,lb,ub

def sml(m,c,U,Y,cv,b0,lb,ub):
    obj=SMLCovResidual(U,Y,cv['CY'],cv['CU'],cv['CYU'],lambda b,u:m.output(c,b,u),lambda b:m.gmodel(c,b))
    so=least_squares(obj.fun,b0,jac=obj.jac,bounds=(lb,ub),x_scale='jac',ftol=1e-9,xtol=1e-9,gtol=1e-7,max_nfev=160)
    e,j=obj.calc(so.x);score=np.max(abs(j.T@e)/np.maximum(np.linalg.norm(j,axis=0)*np.linalg.norm(e),1e-300))
    return so.x,dict(success=bool(so.success),status=str(so.message),iterations=int(so.nfev),iteration_unit='function_evaluations',score_max=float(score),cost=float(2*so.cost))

def run_one(rid,window,methods=METHODS):
    s,m,lb,ub=setup(window);a=get_data(rid,window);li=a['line_indices'];cf=m.context(a['f']);cl=m.context(a['f'][li])
    U,Y,D=(a[k][li] for k in ('U','Y','D'));K=a['K'][li]
    instr=lambda b:np.linalg.solve(np.eye(4)[None]+K@m.resp(cl,b),D[:,:,None])[:,:,0]
    b0,rank=m.initial(cl,U,Y,s['f'],s['z']);cache={};cv=None
    for method in METHODS:
        if method not in methods and not (method=='M3' and any(x in methods for x in ('M3w','M5n'))):continue
        path=HERE/'raw'/f'{rid}_{window}_{method}.json'
        if path.exists():
            cache[method]=json.loads(path.read_text(encoding='utf-8'));assert cache[method]['config_hash']==CONFIG_HASH;continue
        start=time.perf_counter();b=b0.copy()
        r=dict(id=rid,window=window,method=method,config_hash=CONFIG_HASH,spec=s,P=int(a['P']),N=int(a['N']),nb=m.nb,
               initial_linear_ranks=rank,initial_beta=b0.tolist(),controller=next(c for c in ('lowgain','nominal','highgain') if c in rid),
               model='local additive ZOH C0/C1 + D0/D2 + fixed delays + full modal transients',sampling='zoh',frequency_lines=len(li),
               full_frequency_bins=len(a['f']),comparison='same model within window; local background differs from global full-band model')
        try:
            if method=='M1':
                bstart,rr=m.initial(cf,a['U'],a['Y'],s['f'],s['z']);r['initial_beta']=bstart.tolist();r['initial_linear_ranks']=rr
                b,info=iv_solve(m,cf,a['U'],a['Y'],a['U'],None,bstart,lb,ub)
            elif method in ('M2','M3'):
                b,info=iv_solve(m,cl,U,Y,U if method=='M2' else instr,None,b0,lb,ub)
            else:
                bm=np.asarray(cache['M3']['beta']);r['initial_M3_success']=cache['M3']['converged']
                if cv is None:
                    cv=wc.local_covariances(a['k'],a['U'],a['Y'],bm,lambda bb,uu:m.output(cf,bb,uu,False),int(a['P']),nl=16,alpha=.1)
                    np.savez_compressed(HERE/'raw'/f'{rid}_{window}_local_cov.npz',**cv)
                if method=='M3w':b,info=iv_solve(m,cl,U,Y,instr,wc.precision(cv['covariance_residual']),bm,lb,ub)
                else:b,info=sml(m,cl,U,Y,cv,b0,lb,ub)
            ff,zz=m.modal(b);score=info.get('score_max',float('nan'))
            margin=np.minimum(b[m.iw:]-lb[m.iw:],ub[m.iw:]-b[m.iw:]);bound=bool(np.any(margin<1e-4))
            valid=np.isfinite(b).all();r.update(beta=b.tolist(),frequency_hz=ff.tolist(),damping_percent=zz.tolist(),solver=info,
                converged=bool(info['success'] and valid),strict_score_pass=bool(np.isfinite(score) and score<1e-6),near_bound=bound,
                bound_margins_log=margin.tolist(),iterations=info['iterations'],score=score,finite=bool(valid))
            yh,_=m.output(cl,b,U,False);G=m.resp(cl,b)
            r['line_residual_percent']=float(100*np.linalg.norm(Y-yh)/np.linalg.norm(Y))
            r['minimum_frequency_separation_hz']=float(np.min(np.diff(np.sort(ff)))) if m.nm>1 else None
            np.savez_compressed(path.with_suffix('.npz'),beta=b,f=cl['f'],G=G,residual=Y-yh)
        except Exception as e:
            r.update(beta=b.tolist(),converged=False,strict_score_pass=False,finite=bool(np.isfinite(b).all()),exception=repr(e),traceback=traceback.format_exc())
            ff,zz=m.modal(b);r.update(frequency_hz=ff.tolist(),damping_percent=zz.tolist())
        r['seconds']=time.perf_counter()-start;js(path,r);cache[method]=r
        print(rid,window,method,r['converged'],r.get('near_bound'),r.get('iterations'),np.round(r['frequency_hz'],3),np.round(r['damping_percent'],3),round(r['seconds'],1),flush=True)
    return rid,window

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--all',action='store_true');p.add_argument('--id',default='A_nominal');p.add_argument('--window',default='H493');p.add_argument('--methods',default=','.join(METHODS));p.add_argument('--workers',type=int,default=3);args=p.parse_args()
    methods=args.methods.split(',')
    if args.all:
        jobs=[(rid,w) for w in ('H493','H300','H100','H700') for rid in IDS];js(HERE/'run_requests.json',[dict(id=i,window=w,methods=methods,config_hash=CONFIG_HASH) for i,w in jobs])
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for f in as_completed([pool.submit(run_one,i,w,methods) for i,w in jobs]):print('FINISHED',f.result(),flush=True)
    else:run_one(args.id,args.window,methods)
