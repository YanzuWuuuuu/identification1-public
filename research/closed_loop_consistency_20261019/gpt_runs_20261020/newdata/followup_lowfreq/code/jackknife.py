from common import *
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse
def run_one(rid,omit,method):
    dest=HERE/'raw'/f'{rid}_omit{omit}_{method}.json'
    if dest.exists():return rd(dest)
    t=time.perf_counter();old=fit(rid,method)
    r=dict(id=rid,omit=omit,method=method,initial_rule='fixed original saved initial_beta',old_iterations=old['iterations'])
    try:
        a=spectrum(rid,omit=omit)
        sel=np.arange(len(a['f'])) if method=='M1' else a['line_indices']
        c=m.Ctx(a['f'][sel],TAU,'zoh');U,Y,D,K=(a[x][sel] for x in ('U','Y','D','K'))
        model=lambda b,u:m.out_jac(c,b,u)
        Uh=U if method=='M1' else lambda b:np.linalg.solve(np.eye(4)[None]+K@m.resp(c,b),D[:,:,None])[:,:,0]
        b0=np.array(old['initial_beta'])
        b,si=wc.iv_solve(model,U,Y,Uh,np.ones_like(U.real),b0,m.IW,m.ITR,maxiter=80,tol=1e-5,cap=.05)
        ff,zz=m.modal(b);finite=bool(np.all(np.isfinite(b)))
        r.update(P=int(a['P']),N=int(a['N']),kept_periods=a['kept_periods'],beta=b,initial_beta=b0,solver_info=si,
                 converged=bool(si.get('success') and finite and si['iterations']<80),finite_endpoint=finite,
                 frequency_hz=ff,damping_percent=zz,strict_score_pass=bool(si.get('score_max',np.inf)<1e-6),
                 modal_order=np.argsort(ff),seconds=time.perf_counter()-t)
    except Exception as exc:r.update(converged=False,exception=repr(exc),traceback=traceback.format_exc())
    js(dest,r);print('JK',rid,omit,method,r['converged'],r.get('solver_info',{}).get('iterations'),round(time.perf_counter()-t,2),flush=True)
    return r
def summarize():
    results=[]
    for rid in IDS:
        P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P']);folds=[]
        for i in range(P):
            a,b=(rd(HERE/'raw'/f'{rid}_omit{i}_{method}.json') for method in ('M1','M3'))
            row=dict(omit=i,M1=a,M3=b)
            if 'beta'in a and 'beta'in b:
                row['paired_percent']=100*np.expm1(np.array(a['beta'])[MOD]-np.array(b['beta'])[MOD])
            folds.append(row)
        valid=all(x['M1']['converged'] and x['M3']['converged'] for x in folds)
        finite=all('paired_percent' in x and np.all(np.isfinite(x['paired_percent'])) for x in folds)
        row=dict(id=rid,P=P,all_folds_converged=valid,all_folds_finite=finite,
                 failed_M1=sum(not x['M1']['converged'] for x in folds),failed_M3=sum(not x['M3']['converged'] for x in folds),
                 strict_failed_M1=sum(not x['M1'].get('strict_score_pass',False) for x in folds),
                 strict_failed_M3=sum(not x['M3'].get('strict_score_pass',False) for x in folds),folds=folds)
        if finite:
            v=np.array([x['paired_percent'] for x in folds])
            row.update(paired_mean=v.mean(0),paired_se=np.sqrt((P-1)/P*np.sum((v-v.mean(0))**2,axis=0)),
                       se_interpretation='all endpoints, no failed samples deleted; inferentially qualified only if all folds converge')
            for method in ('M1','M3'):
                v=np.array([x[method]['frequency_hz']+x[method]['damping_percent'] for x in folds])
                row[method+'_se']=np.sqrt((P-1)/P*np.sum((v-v.mean(0))**2,axis=0))
        results.append(row)
    js(HERE/'analysis/jackknife.json',results)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=3);args=p.parse_args()
    jobs=[(r,i,method) for r in IDS for i in range(int(np.load(SRC/'prepared'/f'{r}.npz')['P'])) for method in ('M1','M3')]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        fut=[pool.submit(run_one,*j) for j in jobs]
        for f in as_completed(fut):f.result()
    summarize()
