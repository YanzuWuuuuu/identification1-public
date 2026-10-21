from uncertainty_bias import *
from scipy.stats import chi2
import csv,itertools

def fit(rid,method):return json.loads((HERE/'fits'/f'{rid}_{method}.json').read_text(encoding='utf-8'))
def unc(rid,method):
    p=HERE/'analysis'/f'{rid}_{method}_uncertainty.json'
    return json.loads(p.read_text(encoding='utf-8')) if p.exists() else {}
def vec(r):return np.r_[r['frequency_hz'],r['damping_percent']]
def sd(r):return np.r_[r.get('frequency_u_hz',[np.nan]*3),r.get('damping_u_percentage_points',[np.nan]*3)]
def savecsv(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/'analysis'/name).open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(rows)

def normalized_validation():
    rows=[];rng=np.random.default_rng(421020)
    for ctrl in CTRS:
        a=np.load(HERE/'analysis'/f'{ctrl}_orthogonal_frf.npz');f=a['f'];Cref=np.load(HERE/'analysis'/f'{ctrl}_noise.npz')['Phi']*P0/int(a['P_common'])
        iv=np.linalg.inv(Cref);ci=chi2.ppf([.025,.975],2*len(f)*16);nf=len(f)
        for method in ('M1','M2','M3','M3w','M5n'):
            rr=fit('A_'+ctrl,method);b=np.asarray(rr['beta']);G,Gj=g_jac(m.Ctx(f,TAU,'zoh'),b)
            E=a['Y']-G@a['U'];t=2*np.einsum('kor,koi,kir->',E.conj(),iv,E,optimize=True).real
            dof=2*len(f)*16
            row=dict(controller=ctrl,method=method,model_source='A_'+ctrl,reference_source='C four independent runs',
                complex_G_difference_percent=float(100*np.linalg.norm(G-a['G'])/np.linalg.norm(a['G'])),
                inverse_free_difference_percent=float(100*np.linalg.norm(a['TY']-G@a['TU'])/np.linalg.norm(a['TY'])),
                chi_square_type_statistic=float(t),known_fixed_model_expectation=dof,known_fixed_model_ci95_low=float(ci[0]),known_fixed_model_ci95_high=float(ci[1]),
                statistic_over_expectation=float(t/dof),ideal_fixed_model_exceeds95=bool(t>ci[1]),
                scope='fixed exact model + known Gaussian PSD + independent frequencies ideal reference only; B PSD estimated; A parameters estimated; no exact chi-square claim',
                solver_converged=rr['converged'],strict_score_pass=rr['strict_score_pass'])
            cp=HERE/'analysis'/f'A_{ctrl}_{method}_covariance.npz'
            if cp.exists():
                P=np.load(cp)['covariance'];ev,V=np.linalg.eigh((P+P.T)/2)
                if ev.min()>=-1e-9*max(ev.max(),1e-300):
                    Lp=V*np.sqrt(np.maximum(ev,0))[None,:]
                    JJ=np.einsum('koip,kir->korp',Gj,a['U'],optimize=True)
                    Cp=np.einsum('korp,pq,kirq->kroi',JJ,P,JJ.conj(),optimize=True)
                    Ct=Cp+Cref[:,None];It=np.linalg.inv(Ct)
                    et=E.transpose(0,2,1)
                    stat=2*np.einsum('kro,kroi,kri->',et.conj(),It,et,optimize=True).real
                    L=np.linalg.cholesky(Cref);tb=[]
                    # Retain cross-frequency and cross-experiment dependence from
                    # the common real parameter perturbation in A.
                    for j in range(499):
                        db=Lp@rng.normal(size=m.NB)
                        eps=(rng.normal(size=(nf,4,4))+1j*rng.normal(size=(nf,4,4)))/np.sqrt(2)
                        pred=np.einsum('korp,p->kor',JJ,db)+L@eps
                        e=pred.transpose(0,2,1)
                        tb.append(float(2*np.einsum('kro,kroi,kri->',e.conj(),It,e,optimize=True).real))
                    row.update(parameter_uncertainty_corrected_statistic=float(stat),conditional_gaussian_mc_mean=float(np.mean(tb)),
                        conditional_gaussian_mc_ci95=np.percentile(tb,[2.5,97.5]).tolist(),conditional_gaussian_mc_p=float((1+np.sum(np.asarray(tb)>=stat))/500),
                        conditional_gaussian_mc_note='linearized A-parameter covariance + independent B-noise; conditions and estimated PSD limitations retained')
            rows.append(row)
    js(HERE/'analysis/frf_normalized_validation.json',rows)
    savecsv('frf_normalized_validation.csv',[{k:v for k,v in r.items() if not isinstance(v,list)} for r in rows])

def core_decisions():
    pairs=[]
    for method in ('M1','M2','M3','M3w','M5n'):
        for a,b in itertools.combinations(CTRS,2):
            ra,rb=fit('A_'+a,method),fit('A_'+b,method);ua,ub=unc('A_'+a,method),unc('A_'+b,method)
            diff=vec(rb)-vec(ra);comb=np.sqrt(sd(ua)**2+sd(ub)**2)
            for j in range(6):
                pairs.append(dict(method=method,controller_a=a,controller_b=b,parameter=['f1','f2','f3','zeta1','zeta2','zeta3'][j],
                    observed_difference=float(diff[j]),combined_standard_uncertainty=float(comb[j]),z_score=float(diff[j]/comb[j]),
                    less_than_2u=bool(abs(diff[j])<2*comb[j]),greater_than_2u=bool(abs(diff[j])>2*comb[j]),
                    numerical_conditions_ok=bool(ua.get('locally_valid_numerics',False) and ub.get('locally_valid_numerics',False))))
    savecsv('controller_invariance_pairs.csv',pairs)
    drift=[]
    for a,b in itertools.combinations(CTRS,2):
        ba=json.loads((HERE/'analysis'/f'A_{a}_bias.json').read_text());bb=json.loads((HERE/'analysis'/f'A_{b}_bias.json').read_text())
        va,vb=vec(fit('A_'+a,'M3')),vec(fit('A_'+b,'M3'))
        pa=va*np.array(ba['pred_MminusC_percent'])/100;pb=vb*np.array(bb['pred_MminusC_percent'])/100
        observed=vec(fit('A_'+b,'M1'))-vec(fit('A_'+a,'M1'));reference=vb-va;pred=pb-pa
        excess=observed-reference
        for j in range(6):drift.append(dict(a=a,b=b,parameter=['f1','f2','f3','zeta1','zeta2','zeta3'][j],
            M1_observed_drift=float(observed[j]),M3_observed_drift=float(reference[j]),predicted_bias_drift=float(pred[j]),
            observed_excess_drift=float(excess[j]),raw_drift_same_sign=bool(observed[j]*pred[j]>0),excess_drift_same_sign=bool(excess[j]*pred[j]>0),
            raw_drift_to_prediction=float(observed[j]/pred[j]),excess_drift_to_prediction=float(excess[j]/pred[j])))
    savecsv('bias_controller_drift.csv',drift)
    core={}
    for method in ('M3','M3w'):
        rr=[r for r in pairs if r['method']==method]
        core[method]=dict(passing_pairs=sum(r['less_than_2u'] for r in rr),total_pairs=len(rr),all_below_2u=all(r['less_than_2u'] for r in rr),
                          all_numerical_conditions=all(r['numerical_conditions_ok'] for r in rr))
    core['M1_zeta1_or2_drift_above_2u']=any(r['greater_than_2u'] for r in pairs if r['method']=='M1' and r['parameter'] in ('zeta1','zeta2'))
    core['prediction_note']='No predeclared numerical meaning of same order: report signed ratios; do not invent a pass cutoff. Show M3 drift separately because it is not zero.'
    js(HERE/'analysis/core_decisions.json',core)
    rep=[]
    for ctrl in CTRS:
        for method in ('M1','M2','M3','M3w','M5n'):
            ra,rb=fit('A_'+ctrl,method),fit('C_'+ctrl+'_pppp',method)
            du=vec(rb)-vec(ra);sc=np.sqrt(sd(unc('A_'+ctrl,method))**2+sd(unc('C_'+ctrl+'_pppp',method))**2)
            for j in range(6):rep.append(dict(controller=ctrl,method=method,parameter=['f1','f2','f3','zeta1','zeta2','zeta3'][j],
                A_C_repeat_difference=float(du[j]),two_run_sample_sd=abs(float(du[j]))/np.sqrt(2),combined_u=float(sc[j]),normalized_difference=float(du[j]/sc[j]),
                note='two independent runs only; absolute difference/sqrt2 has one degree of freedom, not precise uncertainty validation'))
    savecsv('A_C_independent_repeats.csv',rep)
    print('DECISIONS',core,flush=True)

def scans():
    groups=dict(amplitude=[('D_nominal_0p6',.6),('D_nominal_0p8',.8),('D_nominal_1p0',1.),('D_nominal_1p2',1.2),('A_nominal',1.5)],
                repeat=[('A_nominal',0)]+[('C_nominal_'+s,j+1) for j,s in enumerate(('pppp','pmpm','ppmm','pmmp'))],
                position=[('E_nominal_m1p0',-1.),('E_nominal_m0p5',-.5),('A_nominal',0),('E_nominal_0p5',.5),('E_nominal_1p0',1.)])
    rows=[];scores=[]
    for kind,rs in groups.items():
        methods=['M2','M3','M3w'] if kind in ('amplitude','repeat') else ['M1','M2','M3','M3w']
        for method in methods:
            rr=[fit(rid,method) for rid,x in rs];us=[unc(rid,method) for rid,x in rs]
            xx=np.array([vec(r) for r in rr]);ss=np.array([sd(r) for r in us]);ww=1/ss**2;mu=(ww*xx).sum(0)/ww.sum(0);z=(xx-mu)/ss
            for k,((rid,x),r,u) in enumerate(zip(rs,rr,us)):
                for j in range(6):rows.append(dict(kind=kind,id=rid,x=x,method=method,parameter=['f1','f2','f3','zeta1','zeta2','zeta3'][j],
                    estimate=float(xx[k,j]),standard_uncertainty=float(ss[k,j]),z_score=float(z[k,j]),solver_converged=r['converged'],
                    numerical_conditions_ok=u.get('locally_valid_numerics',False)))
            if kind!='position':
                for j in range(6):
                    valid=all(r['converged'] for r in rr) and all(u.get('locally_valid_numerics',False) for u in us)
                    ratio=np.std(xx[:,j],ddof=1)/np.sqrt(np.mean(ss[:,j]**2))
                    scores.append(dict(kind=kind,method=method,parameter=['f1','f2','f3','zeta1','zeta2','zeta3'][j],
                        observed_sd_to_rms_u=float(ratio),within_2u_count=int(np.sum(abs(z[:,j])<=2)),total=5,
                        raw_z_rms=float(np.sqrt(np.mean(z[:,j]**2))),expected_z_rms_after_mean_fit=float(np.sqrt(4/5)),
                        ratio_in_predeclared_0p7_to_1p4=bool(.7<=ratio<=1.4),all_numerical_conditions=valid,
                        interpretation='amplitude changes: linearity diagnostic, not IID repetitions' if kind=='amplitude' else 'same amplitude independent records; 5 samples, weak coverage evidence'))
    savecsv('scans.csv',rows);js(HERE/'analysis/scans_summary.json',scores)
    old=HERE.parents[1]/'gpt_runs_20261019/followup_DE/real_d/raw';oldrows=[]
    for rid,x in groups['position']:
        if x==-1:continue
        zs=str(int(x)) if float(x).is_integer() else str(x)
        for method in ('M2','M3','M3w'):
            p=old/f'z{zs}_P08_s01_{method}.json'
            if not p.exists():continue
            before=json.loads(p.read_text());after=fit(rid,method)
            for j in range(6):oldrows.append(dict(z_mm=x,method=method,parameter=['f1','f2','f3','zeta1','zeta2','zeta3'][j],
                old=float(vec(before)[j]),new=float(vec(after)[j]),change_percent=float(100*(vec(after)[j]/vec(before)[j]-1)),
                note='different acquisition/assembly conditions; descriptive, not estimator error'))
    savecsv('old_new_position_comparison.csv',oldrows)

def compile_table():
    rows=[]
    for p in sorted((HERE/'fits').glob('*.json')):
        r=json.loads(p.read_text());u=unc(r['id'],r['method']);v=vec(r);ss=sd(u)
        for j in range(6):rows.append(dict(id=r['id'],group=r['group'],controller=r['controller'],z_mm=r['z_mm'],gain=r['gain'],P=r['P'],method=r['method'],
            parameter=['f1','f2','f3','zeta1','zeta2','zeta3'][j],estimate=float(v[j]),standard_uncertainty=float(ss[j]),converged=r['converged'],
            iterations=r['iterations'],score=r['score_normalized'],strict_score_pass=r['strict_score_pass'],locally_valid_uncertainty=u.get('locally_valid_numerics',False),
            uncertainty_interpretation=u.get('interpretation','unavailable')))
    savecsv('all_parameter_estimates.csv',rows)

if __name__=='__main__':
    normalized_validation();core_decisions();scans();compile_table()
