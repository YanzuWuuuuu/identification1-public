from profile_core import *
from concurrent.futures import ProcessPoolExecutor,as_completed
import traceback

def line_data(rid,omit=(),reference=True):
    x=np.load(SRC/'prepared'/f'{rid}.npz');P=int(x['P']);keep=[i for i in range(P) if i not in omit]
    kp=np.arange(18,151);f=kp/6.
    a=dict(f=f,line_indices=np.arange(len(f)),K=old.evalK(f,rid) if reference else np.zeros((len(f),4,4),complex),P=len(keep))
    for key in ('U','Y','D'):
        per=np.fft.rfft(x[key].reshape(P,30000,4),axis=1)[:,kp]
        a[key]=per[keep].sum(0)
    return a

def solve_record(rid,method,omit,b0):
    tag='full' if not omit else 'omit'+'_'.join(map(str,omit))
    dest=ROOT/'raw'/f'p7_real_{rid}_{method}_{tag}.json'
    if dest.exists():return rd(dest)
    m=ResidueModel('L','a',joint=False)
    r=dict(id=rid,method=method,omitted=list(omit),success=False)
    try:
        a=line_data(rid,omit,reference=method=='M3');pr=ProfileIV(m,a,method);b,info=pr.solve(b0,maxouter=80)
        r.update(beta=b,solver=info,success=info['success'],log_parameters=b[m.nl:],
                 frequency_hz=m.modal(b)[0],damping_percent=m.modal(b)[1])
    except Exception as e:r.update(exception=repr(e),traceback=traceback.format_exc())
    js(dest,r);return r

def run(rid):
    t=time.perf_counter();m=ResidueModel('L','a',joint=False);P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])
    dest=ROOT/'analysis'/f'p7_real_{rid}.json'
    if dest.exists():return
    full={};folds={}
    for method in ('M2','M3'):
        fb,_=old.basefit(rid,'L',method);full[method]=solve_record(rid,method,(),m.from_full(fb))
        folds[method]=[solve_record(rid,method,(i,),np.asarray(full[method]['beta'])) for i in range(P)]
    pairs={}
    for i in range(P):
        for j in range(i+1,P):
            pairs[i,j]=solve_record(rid,'M2',(i,j),np.asarray(full['M2']['beta']))
    lp2=np.asarray(full['M2']['log_parameters']);lp3=np.asarray(full['M3']['log_parameters'])
    lpf=np.array([x['log_parameters'] for x in folds['M2']])
    lpr=np.array([x['log_parameters'] for x in folds['M3']])
    jk=P*lp2-(P-1)*lpf.mean(0)
    nested=[]
    for i in range(P):
        other=np.array([pairs[min(i,j),max(i,j)]['log_parameters'] for j in range(P) if j!=i])
        nested.append((P-1)*lpf[i]-(P-2)*other.mean(0))
    nested=np.asarray(nested)
    # Compute SE in a common fixed relative scale for each modal coordinate.
    main=100*np.expm1(jk-lp3)
    nestedrel=100*np.expm1(nested-lp3);refrel=100*np.expm1(lpr-lp3)
    sejk=old.jk_se(nestedrel);seref=old.jk_se(refrel);combined=np.hypot(sejk,seref)
    paired=old.jk_se(100*np.expm1(nested-lpr))
    entries=list(full.values())+folds['M2']+folds['M3']+list(pairs.values())
    result=dict(id=rid,P=P,model='Original free C0/C1 + fixed original delays, exact ZOH, original 158 parameters',
                jackknife_log_parameters=jk,difference_vs_M3_percent=main,
                corrected_estimator_se_percent=sejk,M3_se_percent=seref,combined_se_percent=combined,
                paired_se_percent=paired,passes_one_combined_se=np.abs(main)<combined,
                failures=sum(not x['success'] for x in entries),total_fits=len(entries),
                full_M2=full['M2'],full_M3=full['M3'],seconds=time.perf_counter()-t,
                qualification='Nested delete-two fits estimate the SE of the corrected estimator; primary gate uses quadrature combined SE, paired SE is also reported.')
    js(dest,result);print('P7 real',rid,result['failures'],np.round(main[3:],3),np.round(combined[3:],3),flush=True)

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=4) as pool:
        for fu in as_completed([pool.submit(run,rid) for rid in old.IDS]):fu.result()
