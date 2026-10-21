from common import *
from noise_bias import local_phi,cache_moments,predict
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse
def one(rid,i):
    dest=HERE/'raw'/f'{rid}_omit{i}_prediction.json'
    if dest.exists():return rd(dest)
    t=time.perf_counter()
    r=dict(id=rid,omit=i)
    try:
        f1=rd(BASE/'raw'/f'{rid}_omit{i}_M1.json');f3=rd(BASE/'raw'/f'{rid}_omit{i}_M3.json')
        r.update(M1_converged=f1['converged'],M3_converged=f3['converged'],
                 M1_strict=f1['strict_score_pass'],M3_strict=f3['strict_score_pass'])
        a=spectrum(rid,omit=i);ext=spectrum(rid,omit=i,extended=True)
        b=np.array(f3['beta']);C,nc,_,_=local_phi(ext,b,a['f'],.5)
        cache=cache_moments(a,b);pred=predict(cache,int(a['N'])*C)
        obs=100*np.expm1(np.array(f1['beta'])[MOD]-b[MOD])
        r.update(P=a['P'],N=a['N'],observed=obs,pred_M=pred['pred_M_percent'],pred_MminusC=pred['pred_MminusC_percent'],
                 difference_M=obs-pred['pred_M_percent'],difference_MminusC=obs-pred['pred_MminusC_percent'],
                 numerical_M=pred['numerical_M'],numerical_MminusC=pred['numerical_MminusC'],
                 minimum_neighbors=nc.min(),maximum_neighbors=nc.max(),success=True,seconds=time.perf_counter()-t)
        np.savez_compressed(HERE/'raw'/f'{rid}_omit{i}_prediction.npz',f=a['f'],Phi=C,counts=nc,
            M=pred['M'],Ccurvature=pred['Ccurvature'],b=pred['b'],delta_M=pred['delta_M'],delta_MminusC=pred['delta_MminusC'])
    except Exception as e:r.update(success=False,exception=repr(e),traceback=traceback.format_exc())
    js(dest,r);print('PRED_JK',rid,i,r['success'],round(time.perf_counter()-t,2),flush=True)
    return r
def summarize():
    out=[]
    for rid in IDS:
        old=rd(BASE/'analysis'/f'{rid}_noise_bias.json');P=old['P'];folds=[rd(HERE/'raw'/f'{rid}_omit{i}_prediction.json') for i in range(P)]
        row=dict(id=rid,P=P,folds=folds,prediction_failures=sum(not x['success'] for x in folds),
                 original_solver_failures=sum(not(x['M1_converged'] and x['M3_converged']) for x in folds),
                 original_strict_failures=sum(not(x['M1_strict'] and x['M3_strict']) for x in folds),
                 observed=old['observed'])
        if all(x['success'] for x in folds):
            for key in ('observed','pred_M','pred_MminusC','difference_M','difference_MminusC'):
                a=np.array([x[key] for x in folds]);se=np.sqrt((P-1)/P*np.sum((a-a.mean(0))**2,axis=0))
                row[key+'_se']=se;row[key+'_fold_mean']=a.mean(0)
            for key in ('M','MminusC'):
                pred=np.array(old['windows']['0.5']['pred_'+key+'_percent'])
                row['pred_'+key]=pred
                row['z_joint_'+key]=(np.array(old['observed'])-pred)/row['difference_'+key+'_se']
                row['detectable_using_joint_'+key]=abs(pred)>2*row['difference_'+key+'_se']
                aa=np.array([x['observed'] for x in folds]);bb=np.array([x['pred_'+key] for x in folds])
                row['observed_prediction_JK_covariance_'+key]=(P-1)/P*np.sum((aa-aa.mean(0))*(bb-bb.mean(0)),axis=0)
        out.append(row)
    js(HERE/'analysis/prediction_jackknife.json',out)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=2);args=p.parse_args()
    jobs=[(r,i) for r in IDS for i in range(rd(BASE/'analysis'/f'{r}_noise_bias.json')['P'])]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        fut=[pool.submit(one,*j) for j in jobs]
        for f in as_completed(fut):f.result()
    summarize()
