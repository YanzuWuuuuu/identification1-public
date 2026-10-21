from core import *
import csv
REPRO={('L',1):.10,('H300',1):.01091974928,('H300',2):.006012159718,
       ('H300',3):.001069751086,('H493',0):.001288859477}
PARAMS=['common','a1','a2','a3','a4','b1','b2','b3','b4']
def delayvec(r):
    d=r['delay'];return np.r_[d['common_us'],d['output_us'],d['input_us']]
def load(p):
    q=ROOT/'refinements'/p.name
    return rd(q if q.exists() else p)
def csvout(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (ROOT/'analysis'/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def run():
    fixed=[load(p) for p in (ROOT/'raw').glob('fixed*.json')]
    joint=[load(p) for p in (ROOT/'raw').glob('joint*full.json')]
    folds=[load(p) for p in (ROOT/'raw').glob('joint*omit*.json')]
    jr={(r['id'],r['window'],r['method']):r for r in joint}
    jk={}
    for r in joint:
        key=(r['id'],r['window'],r['method'])
        ff=sorted([q for q in folds if (q['id'],q['window'],q['method'])==key],key=lambda q:q['omit'])
        dv=np.array([delayvec(q) for q in ff])
        jk['|'.join(key)]=dict(P=len(ff),successful=sum(q['success'] for q in ff),
            full_success=r['success'],valid=bool(r['success'] and all(q['success'] for q in ff)),
            delay_se_us=jk_se(dv),delay_fold_values_us=dv,
            frequency_se_hz=jk_se([q['frequency_hz'] for q in ff]),
            damping_se_pp=jk_se([q['damping_percent'] for q in ff]))
    paired=[]
    for rid in IDS:
        for w in SPECS:
            r1=jr[rid,w,'M1'];r3=jr[rid,w,'M3'];P=r1['P']
            f1=sorted([q for q in folds if (q['id'],q['window'],q['method'])==(rid,w,'M1')],key=lambda q:q['omit'])
            f3=sorted([q for q in folds if (q['id'],q['window'],q['method'])==(rid,w,'M3')],key=lambda q:q['omit'])
            se=jk_se(np.array([delayvec(x)-delayvec(y) for x,y in zip(f1,f3)]));dv=delayvec(r1)-delayvec(r3)
            valid=jk['|'.join((rid,w,'M1'))]['valid'] and jk['|'.join((rid,w,'M3'))]['valid']
            for i,param in enumerate(PARAMS):
                paired.append(dict(id=rid,window=w,parameter=param,difference_us=dv[i],paired_se_us=se[i],
                    exceeds_2se=bool(abs(dv[i])>2*se[i]),valid=valid))
    tests=[]
    for w in SPECS:
        for g in ['A','C']:
            ids=[g+'_'+ct+('_pppp' if g=='C' else '') for ct in CTRS]
            for p in PARAMS:
                rr=[next(r for r in paired if (r['id'],r['window'],r['parameter'])==(rid,w,p)) for rid in ids]
                x=np.array([r['difference_us'] for r in rr]);dx=np.diff(x)
                monotonic=bool(np.all(dx>0) or np.all(dx<0))
                tests.append(dict(window=w,group=g,parameter=p,monotonic=monotonic,
                    all_above_2se=all(r['exceeds_2se'] for r in rr),valid=all(r['valid'] for r in rr),
                    pass_rule_c=bool(monotonic and all(r['exceeds_2se'] and r['valid'] for r in rr))))
    fr=[]
    for r in fixed:
        for i in SPECS[r['window']]['report']:
            ref=REPRO.get((r['window'],i),float('nan'));dd=r['delta_damping_pp'][i]
            fr.append(dict(id=r['id'],window=r['window'],method=r['method'],shift_us=r['shift_us'],
                mode=i,frequency_hz=r['frequency_hz'][i],damping_percent=r['damping_percent'][i],
                change_f_hz=r['delta_frequency_hz'][i],change_zeta_pp=dd,repro_zeta_sd_pp=ref,
                change_over_repro_sd=abs(dd)/ref,success=r['success'],score=r['solver'].get('score'),
                original_failure=r.get('original_failure',False)))
    jc=[];dp=[];qrows=[];startrows=[]
    for r in joint:
        key='|'.join((r['id'],r['window'],r['method']));kk=jk[key]
        _,m,_,_=setup(r['window']);base,_=basefit(r['id'],r['window'],r['method']);f0,z0=m.modal(base)
        for i in SPECS[r['window']]['report']:
            ref=REPRO.get((r['window'],i),np.nan);dd=r['damping_percent'][i]-z0[i]
            jc.append(dict(id=r['id'],window=r['window'],method=r['method'],mode=i,
                frequency_hz=r['frequency_hz'][i],damping_percent=r['damping_percent'][i],
                change_f_hz=r['frequency_hz'][i]-f0[i],change_zeta_pp=dd,
                frequency_jk_se_hz=kk['frequency_se_hz'][i],damping_jk_se_pp=kk['damping_se_pp'][i],
                repro_zeta_sd_pp=ref,change_over_repro_sd=abs(dd)/ref,success=r['success'],jackknife_valid=kk['valid'],
                qualifies_a=bool(r['window']!='L' and r['method']=='M3' and r['success'] and abs(dd)>ref)))
        for i,p in enumerate(PARAMS):
            dp.append(dict(id=r['id'],window=r['window'],method=r['method'],parameter=p,
                estimate_us=delayvec(r)[i],jackknife_se_us=kk['delay_se_us'][i],
                valid=kk['valid'],full_success=r['success'],successful_folds=kk['successful'],P=kk['P']))
        if r['method']=='M3':
            q0=r['qe_baseline']['ratio'];q1=r['qe_fixed_noise']['ratio'];q2=r['qe_reestimated_noise']['ratio']
            qrows.append(dict(id=r['id'],window=r['window'],baseline_qe=q0,joint_qe_fixed_noise=q1,
                joint_qe_reestimated_noise=q2,reduction_percent=100*(1-q1/q0),
                reduction_reestimated_percent=100*(1-q2/q0),success=r['success'],
                joint_eq=r['qe_fixed_noise']['EQ'],baseline_eq=r['qe_baseline']['EQ']))
        ss=r.get('all_starts',[])
        good=[x for x in ss if x['solver']['success']]
        sm=setup(r['window'],joint=True)[1]
        mats=[sm.delay(np.asarray(x['beta']))*1e6 for x in good]
        startrows.append(dict(id=r['id'],window=r['window'],method=r['method'],original_successful_starts=len(good),
            max_delay_spread_us=float(np.max(np.ptp(np.asarray(mats),axis=0))) if mats else np.nan))
    cross=[]
    for w in SPECS:
        for method in ['M1','M3']:
            rr=[jr[rid,w,method] for rid in IDS];dv=np.array([delayvec(r) for r in rr])
            for i,p in enumerate(PARAMS):cross.append(dict(window=w,method=method,parameter=p,
                mean_us=dv[:,i].mean(),sd_across_six_us=dv[:,i].std(ddof=1),all_success=all(r['success'] for r in rr)))
    sim=[];sim_summary=[]
    for p in sorted((ROOT/'raw').glob('sim_*.json')):
        rr=load(p);raw=rd(p)
        for method,r in rr.get('methods',{}).items():
            truth=rr['truth_delay'];tv=np.r_[truth['common_us'],truth['output_us'],truth['input_us']]
            err=delayvec(r)-tv
            for i,param in enumerate(PARAMS):sim.append(dict(replicate=rr['replicate'],method=method,parameter=param,
                error_us=err[i],success=r['solver']['success'],original_success=raw['methods'][method]['solver']['success'],
                score=r['solver'].get('score')))
    for method in ['M1','M3']:
        for param in PARAMS:
            v=[r for r in sim if r['method']==method and r['parameter']==param];x=np.array([r['error_us'] for r in v])
            ok=np.array([r['success'] for r in v])
            sim_summary.append(dict(method=method,parameter=param,N=len(x),failed=int(sum(~ok)),
                original_failed=sum(not r['original_success'] for r in v),mean_all_us=x.mean(),
                sd_all_us=x.std(ddof=1),rmse_all_us=np.sqrt(np.mean(x*x)),
                mean_converged_us=x[ok].mean() if ok.any() else np.nan,
                sd_converged_us=x[ok].std(ddof=1) if sum(ok)>1 else np.nan))
    highq=[q for q in qrows if q['window']!='L' and q['success']]
    decision=dict(a_any=any(r['qualifies_a'] for r in jc),
        a_count=sum(r['qualifies_a'] for r in jc),a_total=sum(r['window']!='L' and r['method']=='M3' for r in jc),
        b_mean_reduction_percent=float(np.mean([q['reduction_percent'] for q in highq])),
        b_ratio_of_means_reduction_percent=float(100*(1-np.mean([q['joint_qe_fixed_noise'] for q in highq])/np.mean([q['baseline_qe'] for q in highq]))),
        b_reestimated_mean_reduction_percent=float(np.mean([q['reduction_reestimated_percent'] for q in highq])),
        b_count=len(highq),c_any=any(r['pass_rule_c'] for r in tests))
    decision['b_pass']=decision['b_mean_reduction_percent']>30
    decision['recommend_include']=decision['a_any'] or decision['b_pass'] or decision['c_any']
    counts=dict(fixed_total=len(fixed),fixed_original_fail=sum(not rd(p)['success'] for p in (ROOT/'raw').glob('fixed*.json')),
        fixed_final_fail=sum(not r['success'] for r in fixed),joint_full=len(joint),
        joint_high_success=sum(r['success'] and r['window']!='L' for r in joint),
        joint_low_success=sum(r['success'] and r['window']=='L' for r in joint),folds=len(folds),
        high_fold_success=sum(r['success'] and r['window']!='L' for r in folds),
        low_fold_success=sum(r['success'] and r['window']=='L' for r in folds))
    for name,rows in [('fixed_sensitivity.csv',fr),('joint_modal_changes.csv',jc),('joint_delays_and_se.csv',dp),
        ('paired_delay_difference.csv',paired),('rule_c_tests.csv',tests),('residual_check.csv',qrows),
        ('multistart_spread.csv',startrows),('delay_across_records.csv',cross),
        ('simulation_individual.csv',sim),('simulation_summary.csv',sim_summary)]:csvout(name,rows)
    js(ROOT/'analysis/jackknife.json',jk);js(ROOT/'analysis/decision.json',decision);js(ROOT/'analysis/counts.json',counts)
    print(json.dumps(counts,indent=2));print(json.dumps(decision,indent=2))
if __name__=='__main__':run()
