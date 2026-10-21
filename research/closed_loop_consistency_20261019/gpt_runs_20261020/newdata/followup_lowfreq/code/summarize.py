from common import *
import csv,itertools

def csvout(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/'analysis'/name).open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=keys);w.writeheader();w.writerows(rows)

def main():
    jk={r['id']:r for r in rd(HERE/'analysis/jackknife.json')}
    br={r:rd(HERE/'analysis'/f'{r}_noise_bias.json') for r in IDS}
    rows=[];predrows=[];summary=[]
    for rid in IDS:
        b=br[rid];j=jk[rid]
        obs=np.array(b['observed']);pr=np.array(b['windows']['0.5']['pred_MminusC_percent'])
        p1=np.array(b['windows']['1.0']['pred_MminusC_percent']);se=np.array(j['paired_se'])
        for k in range(6):
            name=('f' if k<3 else 'zeta')+str(k%3+1)
            detect=bool(abs(pr[k])>2*se[k])
            row=dict(id=rid,parameter=name,observed=obs[k],prediction=pr[k],jackknife_se=se[k],
                     sandwich_se=b['paired_sandwich']['paired_se'][k],prediction_1Hz=p1[k],
                     abs_prediction_over_2se=abs(pr[k])/(2*se[k]),detectable=detect,
                     z_diagnostic=(obs[k]-pr[k])/se[k],formal_z=(obs[k]-pr[k])/se[k] if detect else '',
                     observed_over_prediction=obs[k]/pr[k],prediction_over_observed=pr[k]/obs[k],
                     below_detection_observed_within_2se=bool(abs(obs[k])<=2*se[k]),
                     all_folds_converged=j['all_folds_converged'],
                     all_folds_strict=not(j['strict_failed_M1'] or j['strict_failed_M3']))
            rows.append(row)
            if k>=3:
                for key,label in [('old_B','B'),('0.5','W0.5'),('1.0','W1.0')]:
                    v=b['old_B'] if key=='old_B' else b['windows'][key]
                    predrows.append(dict(id=rid,parameter=name,spectrum=label,observed=obs[k],
                                         pred_M=v['pred_M_percent'][k],pred_MminusC=v['pred_MminusC_percent'][k]))
        summary.append(dict(id=rid,P=j['P'],M1_failed=j['failed_M1'],M3_failed=j['failed_M3'],
                            M1_strict_failed=j['strict_failed_M1'],M3_strict_failed=j['strict_failed_M3'],
                            mean_M1_iterations=np.mean([x['M1']['solver_info']['iterations'] for x in j['folds']]),
                            mean_M3_iterations=np.mean([x['M3']['solver_info']['iterations'] for x in j['folds']])))
    csvout('paired_tests.csv',rows);csvout('bias_predictions.csv',predrows);csvout('jackknife_status.csv',summary)
    js(HERE/'analysis/paired_tests.json',rows)
    # Special-case investigation: no optimizer bounds exist in the retained solver.
    rid='C_highgain_pppp';j=jk[rid]
    audit=dict(id=rid,explicit_bounds_used=False,original={method:{q:fit(rid,method).get(q) for q in
                 ('converged','iterations','score_normalized','strict_score_pass','frequency_hz','damping_percent')}
                 for method in ('M1','M3')},folds=[])
    for r in j['folds']:
        audit['folds'].append(dict(omit=r['omit'],paired_zeta1=r['paired_percent'][3],
            M1_zeta1=r['M1']['damping_percent'][0],M3_zeta1=r['M3']['damping_percent'][0],
            M1_score=r['M1']['solver_info'].get('score_max'),M3_score=r['M3']['solver_info'].get('score_max'),
            M1_iterations=r['M1']['solver_info']['iterations'],M3_iterations=r['M3']['solver_info']['iterations'],
            M1_converged=r['M1']['converged'],M3_converged=r['M3']['converged'],
            M1_order=r['M1']['modal_order'],M3_order=r['M3']['modal_order']))
    vv=np.array([r['paired_zeta1'] for r in audit['folds']])
    audit.update(paired_range=[vv.min(),vv.max()],paired_mean=vv.mean(),paired_se=j['paired_se'][3],
                 prediction=br[rid]['windows']['0.5']['pred_MminusC_percent'][3])
    js(HERE/'analysis/C_highgain_first_mode_audit.json',audit)
    # D-group amplitude check uses original preparation and all endpoints, including failures.
    amp=[]
    ref=rms_record('A_nominal')
    asp=oldspec('A_nominal');refli=asp['line_indices']
    refd=asp['D'][refli]/int(asp['N'])/1.5
    refcov=refd[:,:,None]*refd[:,None,:].conj()
    af=asp['f'][refli];ag=m.resp(m.Ctx(af,TAU,'zoh'),np.array(fit('A_nominal','M3')['beta']))
    ak=asp['K'][refli]
    refyd=np.linalg.solve(np.eye(4)[None]+ag@ak,(ag@(asp['D'][refli]/int(asp['N']))[:,:,None]))[:,:,0]
    refyrms=np.sqrt(2*np.sum(abs(refyd)**2,axis=0))
    for rid,gain in [('D_nominal_0p6',.6),('D_nominal_0p8',.8),('D_nominal_1p0',1.),('D_nominal_1p2',1.2),('A_nominal',1.5)]:
        r=rms_record(rid);linear=ref*gain/1.5;fr=fit(rid,'M3')
        sp=oldspec(rid);d=sp['D'][sp['line_indices']]/int(sp['N'])
        dcov=(d/gain)[:,:,None]*(d/gain)[:,None,:].conj()
        yd=np.linalg.solve(np.eye(4)[None]+ag@ak,(ag@d[:,:,None]))[:,:,0]
        yrms=np.sqrt(2*np.sum(abs(yd)**2,axis=0))
        amp.append(dict(id=rid,gain=gain,response_rms=r,linear_prediction=linear,
                        response_over_linear=r/linear,relative_linear_error_percent=100*(r/linear-1),
                        normalized_excitation_covariance_relative_difference=np.linalg.norm(dcov-refcov)/np.linalg.norm(refcov),
                        normalized_injection_rms=np.sqrt(2*np.sum(abs(d/gain)**2,axis=0)),
                        same_fixed_model_response_rms=yrms,
                        excitation_pattern_only_rms_change_percent=100*(yrms/(refyrms*gain/1.5)-1),
                        observed_over_fixed_model_percent=100*(r/yrms-1),
                        M3_converged=fr['converged'],M3_strict_score=fr['strict_score_pass'],
                        frequency_hz=fr['frequency_hz'],damping_percent=fr['damping_percent']))
    rms_table=[]
    for a in amp:
        for k in range(4):
            rms_table.append(dict(id=a['id'],gain=a['gain'],output=k+1,rms=a['response_rms'][k],
                rms_linear=a['linear_prediction'][k],relative_linear_error_percent=a['relative_linear_error_percent'][k],
                M3_converged=a['M3_converged']))
    csvout('amplitude_linearity.csv',rms_table)
    scales=[dict(id='A_'+c,**br['A_'+c]['scaled_excitation']) for c in CTRS]
    pairs=[]
    for a,b in itertools.combinations(scales,2):
        sea=jk[a['id']]['paired_se'][4];seb=jk[b['id']]['paired_se'][4];combined=np.hypot(sea,seb)
        diff=a['pred_MminusC_percent'][4]-b['pred_MminusC_percent'][4]
        pairs.append(dict(a=a['id'],b=b['id'],prediction_difference=diff,original_combined_jackknife_se=combined,
                          difference_over_combined_se=abs(diff)/combined,exceeds_3_original_se=bool(abs(diff)>3*combined),
                          prospective_uncertainty_verified=False))
    scaleout=dict(amplitude_records=amp,scalings=scales,pairwise=pairs,
                  max_abs_linear_error_percent=float(np.max(np.abs([r['relative_linear_error_percent'] for r in amp]))),
                  limitation='D group has nominal controller only; RMS linearity does not prove fixed modal dynamics or constant noise under rescaling. Current jackknife SE is not future scaled-run SE.')
    js(HERE/'analysis/scaling.json',scaleout)
    # Pre-registered decisions, without choosing a better window/model.
    keyrecords=['A_nominal','A_highgain','C_nominal_pppp','C_highgain_pppp']
    task1=[r for r in rows if r['id'] in keyrecords and r['parameter']=='zeta2']
    task2mode2=[r for r in rows if r['parameter']=='zeta2']
    task2others=[r for r in rows if r['parameter'] in ('zeta1','zeta3') and not r['detectable']]
    rr=rd(HERE/'analysis/residual_checks.json')
    decisions=dict(
        task1_success=all(.67<=r['prediction_over_observed']<=1.5 for r in task1),
        task1_records=[dict(id=r['id'],prediction_over_observed=r['prediction_over_observed'],
                           passed=.67<=r['prediction_over_observed']<=1.5) for r in task1],
        task2_detectable_zeta2_success=all(abs(r['formal_z'])<2 for r in task2mode2 if r['detectable']),
        task2_mode2_nondetectable=[r['id'] for r in task2mode2 if not r['detectable']],
        task2_low_effect_other_modes_success=all(r['below_detection_observed_within_2se'] for r in task2others),
        task2_not_detectable_other_modes=task2others,
        task3_success=all(r['ratio']<=2 for r in rr),
        task3_pass_count=sum(r['ratio']<=2 for r in rr),task3_total=len(rr),
        numerical_caveat='all original stopping-rule JK fits retained; stricter scores reported separately')
    js(HERE/'analysis/decisions.json',decisions)
    print((HERE/'analysis/decisions.json').read_text(encoding='utf-8'))
if __name__=='__main__':main()
