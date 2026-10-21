from common import *
import csv

def csvout(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/'analysis'/name).open('w',newline='',encoding='utf-8-sig') as h:
        d=csv.DictWriter(h,fieldnames=keys);d.writeheader();d.writerows(rows)

def get_replicates():
    rows=[];summaries=[]
    for w in SPECS:
        s,m,_,_=setup(w)
        for rid in IDS:
            r=rd(HERE/'raw'/f'{rid}_{w}_full_prediction.json')
            P=r['P'];folds=[rd(HERE/'raw'/f'{rid}_{w}_omit{i}_prediction.json') for i in range(P)]
            q=dict(id=rid,window=w,P=P,folds=folds,full=r,
                 failed_predictions=sum(not x['success'] for x in folds),
                 unqualified_folds=sum(not x.get('base_eligible',False) for x in folds))
            if all(x['success'] for x in folds):
                for key in ('observed','pred_M','pred_MminusC','difference_M','difference_MminusC'):
                    q[key+'_se']=jk_se([x[key] for x in folds])
                for j in s['report']:
                    for par,index in [('f',j),('zeta',m.nm+j)]:
                        for formula in ('M','MminusC'):
                            ob=r['observed'][index];pr=r['pred_'+formula][index]
                            so=q['observed_se'][index];sp=q['pred_'+formula+'_se'][index];sj=q['difference_'+formula+'_se'][index]
                            z=(ob-pr)/sj
                            nearzero=abs(pr)<2*so
                            passed=abs(ob)<2*so if nearzero else abs(z)<2
                            valid=bool(r['base_eligible'] and q['unqualified_folds']==0)
                            rows.append(dict(id=rid,controller=ctrl(rid),group=rid[0],window=w,component=j+1,
                                frequency_hz=r['frequency_hz'][j],parameter=par,formula=formula,
                                observed_percent=ob,predicted_percent=pr,observed_SE_percent=so,predicted_SE_percent=sp,
                                joint_SE_percent=sj,z_joint=z,near_zero_prediction=nearzero,raw_rule_pass=passed,
                                eligible=valid,formal_rule_pass=bool(passed and valid),
                                observed_below_2SE=abs(ob)<2*so,P=P,unqualified_folds=q['unqualified_folds'],
                                detection_rule='abs(pred)<2*observed_SE: check abs(obs)<2*observed_SE; otherwise abs(z_joint)<2'))
            summaries.append(q)
    js(HERE/'analysis/jackknife_summary.json',summaries);csvout('high_frequency_tests.csv',rows)
    return rows,summaries

def low_table():
    out=[]
    for q in rd(LOW/'closeout/analysis/prediction_jackknife.json'):
        rid=q['id'];bf=rd(SRC/'fits'/f'{rid}_M3.json');freq=np.exp(np.array(bf['beta'])[128:131])/2/np.pi
        for j in range(3):
            for par,k in [('f',j),('zeta',3+j)]:
                for formula in ('M','MminusC'):
                    ob=q['observed'][k];pr=q['pred_'+formula][k];so=q['observed_se'][k]
                    sp=q['pred_'+formula+'_se'][k];sj=q['difference_'+formula+'_se'][k];z=(ob-pr)/sj
                    nearzero=abs(pr)<2*so;passed=abs(ob)<2*so if nearzero else abs(z)<2
                    out.append(dict(id=rid,controller=ctrl(rid),group=rid[0],window='L3_25',component=j+1,
                         frequency_hz=freq[j],parameter=par,formula=formula,observed_percent=ob,predicted_percent=pr,
                         observed_SE_percent=so,predicted_SE_percent=sp,joint_SE_percent=sj,z_joint=z,
                         near_zero_prediction=nearzero,raw_rule_pass=passed,
                         prior_unqualified_folds=q['original_strict_failures'],scope='existing low-frequency result, unchanged'))
    csvout('low_frequency_all_modes.csv',out);return out

def bread_noise_uncertainty(rid,w):
    # Reevaluate old B-spectrum at polished endpoint for an apples-to-apples noise-source comparison.
    s,m,_,_=setup(w);r=basefit(rid,w,'M3');b=np.array(r['beta']);a=spectrum(rid,w);li=a['line_indices'];f=a['f'][li]
    noise=np.load(HIGH/'analysis'/f'noise_{ctrl(rid)}_{w}.npz')
    Phi=np.empty((len(f),4,4),complex)
    for o in range(4):
        for i in range(4):Phi[:,o,i]=np.interp(f,noise['f'],noise['Phi'][:,o,i].real)+1j*np.interp(f,noise['f'],noise['Phi'][:,o,i].imag)
    c=m.context(f);G=m.resp(c,b);Ud=np.linalg.solve(np.eye(4)[None]+a['K'][li]@G,a['D'][li,:,None])[:,:,0]
    Z=m.output(c,b,Ud)[1];A=np.einsum('kop,koq->pq',Z.conj(),Z,optimize=True).real
    Q=.5*np.einsum('kop,koi,kiq->pq',Z.conj(),a['N']*Phi,Z,optimize=True).real
    iv=scaled_solve(A,np.eye(m.nb));cov=iv@Q@iv.T;d=np.diag(cov)[m.iw:];sd=np.sqrt(np.where(d>=0,d,np.nan))
    ff,zz=m.modal(b)
    return dict(f_u=ff*sd[:m.nm],z_u=zz*sd[m.nm:])

def reproducibility():
    pairs=[];pooled=[]
    for w,s in SPECS.items():
        data={};bs={}
        for rid in IDS:
            data[rid]=rd(HERE/'raw'/f'{rid}_{w}_full_prediction.json');bs[rid]=bread_noise_uncertainty(rid,w)
        for j in s['report']:
            for parameter,value,unit,bunit in [('f','frequency_hz','frequency_u_hz','f_u'),('zeta','damping_percent','damping_u_pp','z_u')]:
                diff=[];unc=[];bu=[]
                for ct in CTRS:
                    ar,cr='A_'+ct,'C_'+ct+'_pppp';a,c=data[ar],data[cr];d=c[value][j]-a[value][j]
                    u1,u2=a[unit][j],c[unit][j];b1,b2=bs[ar][bunit][j],bs[cr][bunit][j]
                    diff.append(d);unc.extend([u1,u2]);bu.extend([b1,b2])
                    pairs.append(dict(window=w,component=j+1,parameter=parameter,controller=ct,
                        A_value=a[value][j],C_value=c[value][j],difference_CminusA=d,
                        relative_difference_percent=100*np.expm1(np.log(c[value][j]/a[value][j])),
                        A_sandwich_u=u1,C_sandwich_u=u2,A_B_spectrum_u=b1,C_B_spectrum_u=b2,
                        combined_same_record_u=np.hypot(u1,u2),
                        eligible=bool(a['base_eligible'] and c['base_eligible'])))
                diff=np.array(diff);sr=np.sqrt(np.sum(diff**2)/6)
                su=np.sqrt(np.mean(np.array(unc)**2));sb=np.sqrt(np.mean(np.array(bu)**2))
                pooled.append(dict(window=w,component=j+1,frequency_hz=np.mean([data[r]['frequency_hz'][j] for r in IDS]),
                    parameter=parameter,unit='Hz' if parameter=='f' else 'damping percentage point',
                    controller_pairs=3,degrees_of_freedom=3,replication_SD=sr,mean_CminusA=np.mean(diff),
                    centered_pair_difference_SD_over_sqrt2=np.std(diff,ddof=1)/np.sqrt(2),
                    conditional_same_record_u_RMS=su,conditional_B_u_RMS=sb,replication_over_same_record=sr/su,
                    replication_over_B=sr/sb))
    csvout('repeat_pairs.csv',pairs);csvout('reproducibility.csv',pooled)
    js(HERE/'analysis/reproducibility.json',dict(pairs=pairs,pooled=pooled))
    return pairs,pooled

def amplitude():
    rows=[]
    for gain in (.6,.8,1.,1.2,1.5):
        for w,s in SPECS.items():
            if gain==1.5:r=rd(HIGH/'raw'/f'A_nominal_{w}_M3.json')
            else:r=rd(HIGH/'amplitude_diagnostic/raw'/f"D_nominal_{str(gain).replace('.','p')}_{w}_M3.json")
            for j in s['report']:
                rows.append(dict(gain=gain,window=w,component=j+1,f=r['frequency_hz'][j],zeta_percent=r['damping_percent'][j],
                      converged=r['converged'],strict_score_pass=r['strict_score_pass'],near_bound=r['near_bound'],
                      role='measurement condition only; single run per amplitude, order/thermal drift confounded'))
    csvout('D_amplitude.csv',rows);return rows

def bandwidth():
    rows=[]
    for ct in CTRS:
        a=np.load(SRC/'wideband_scope_audit'/f'{ct}_broadband.npz');f=a['f'];sv=a['loop_singular_values']
        allcross=[];last=[]
        for j in range(3):
            ix=np.flatnonzero((sv[:-1,j]>=1)&(sv[1:,j]<1)&(f[:-1]>=3)&(f[1:]<=90))
            cr=[]
            for i in ix:
                t=-np.log(sv[i,j])/(np.log(sv[i+1,j])-np.log(sv[i,j]))
                cr.append(float(np.exp(np.log(f[i])+t*np.log(f[i+1]/f[i]))))
            allcross.append(cr);last.append(cr[-1] if cr else np.nan)
        rows.append(dict(controller=ct,last_downward_crossings_hz=last,all_downward_crossings_hz=allcross,
                lower_hz=float(np.nanmin(last)),upper_hz=float(np.nanmax(last)),
                convention='three nonzero singular values of K G_C; last descending unity crossing in 3-90 Hz; not per-axis eigenloop margins or closed-loop -3dB bandwidth'))
    js(HERE/'analysis/bandwidth.json',rows);return rows

def statuses():
    rows=[]
    for w in ('H100','H700'):
        for method in ('M1','M2','M3','M3w','M5n'):
            rr=[rd(HIGH/'raw'/f'{rid}_{w}_{method}.json') for rid in IDS]
            rows.append(dict(window=w,method=method,n=len(rr),failed=sum(not x['converged'] for x in rr),
                bound=sum(x.get('near_bound',True) for x in rr),strict_fail=sum(not x['strict_score_pass'] for x in rr),
                used_for_conclusions=False))
    csvout('excluded_windows_status.csv',rows);return rows

def main():
    high,summaries=get_replicates();low=low_table();pairs,repro=reproducibility();amp=amplitude();bw=bandwidth();ex=statuses()
    csvout('bandwidth_inside_outside.csv',[r for r in low if r['parameter']=='zeta' and r['component']==2]+[r for r in high if r['parameter']=='zeta'])
    counts=[]
    for formula in ('M','MminusC'):
        for par in ('f','zeta'):
            rr=[r for r in high if r['formula']==formula and r['parameter']==par]
            counts.append(dict(formula=formula,parameter=par,n=len(rr),eligible=sum(r['eligible'] for r in rr),
              nearzero=sum(r['near_zero_prediction'] for r in rr),formal_pass=sum(r['formal_rule_pass'] for r in rr),
              raw_pass=sum(r['raw_rule_pass'] for r in rr),
              max_abs_prediction=max(abs(r['predicted_percent']) for r in rr),
              max_abs_observation=max(abs(r['observed_percent']) for r in rr)))
    js(HERE/'analysis/summary.json',dict(counts=counts,bandwidth=bw,reproducibility=repro,
       cases=len(summaries),folds=sum(x['P'] for x in summaries),
       prediction_failures=sum(x['failed_predictions'] for x in summaries),
       unqualified_folds=sum(x['unqualified_folds'] for x in summaries)))
    print(json.dumps(counts,indent=2,default=lambda x:x.item()));print('ANALYSIS COMPLETE')
if __name__=='__main__':main()
