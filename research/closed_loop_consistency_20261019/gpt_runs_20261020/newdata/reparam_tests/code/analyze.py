from metrics import *
import csv,traceback

def write_csv(name,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (ROOT/'analysis'/name).open('w',newline='',encoding='utf-8-sig') as f:
        wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(rows)

def p6():
    rows=[];delays=[];comparisons=[];dispersion=[];modalvalues=[];cross_components=[]
    for p in (ROOT/'raw').glob('*_M3_full*.json'):
        if p.name.startswith('p7'):continue
        r=rd(p)
        if 'beta' not in r:continue
        metric=endpoint_metric(p)
        prefix=p.name[:-5]
        if '_fixed' in prefix:
            stem,suffix=prefix.split('_full');foldstem=stem+'_omit';suffix=suffix
        else:foldstem=prefix.replace('_full','_omit');suffix=''
        folds=[]
        P=int(np.load(SRC/'prepared'/f"{r['id']}.npz")['P'])
        for i in range(P):
            pp=ROOT/'raw'/f'{foldstem}{i}{suffix}.json'
            if pp.exists():folds.append(rd(pp))
        usable=len(folds)==P and all(x.get('success') for x in folds) and r['success']
        d=lambda x:np.r_[x['delay']['common_us'],x['delay']['output_us'],x['delay']['input_us']]
        se=old.jk_se([d(x) for x in folds]) if len(folds)==P and all('delay' in x for x in folds) else np.full(9,np.nan)
        row={k:r.get(k) for k in ('id','window','model','success','common_fixed_us','parameter_count','linear_count','nonlinear_count')}
        row.update(score=r['solver']['score'],common_us=r['delay']['common_us'],common_se_us=se[0],
                   complete_folds=len(folds),failed_folds=sum(not x.get('success') for x in folds),
                   jackknife_valid=usable,QE=metric.get('qe',{}).get('all',{}).get('ratio'),
                   QE333=metric.get('qe',{}).get('mode333',{}).get('ratio'),QE359=metric.get('qe',{}).get('mode359',{}).get('ratio'))
        total=metric.get('qe',{}).get('all',{})
        nreal=8*(int(round((old.SPECS[r['window']]['band'][1]-old.SPECS[r['window']]['band'][0])*6))+1)
        row.update(Q=total.get('Q'),EQ=total.get('EQ'),real_data_dimension=nreal,
                   expected_noise_amplification=total.get('EQ',np.nan)/nreal)
        rows.append(row)
        nf=len(r['frequency_hz'])
        fse=old.jk_se([x['frequency_hz'] for x in folds]) if len(folds)==P and all('frequency_hz' in x for x in folds) else np.full(nf,np.nan)
        zse=old.jk_se([x['damping_percent'] for x in folds]) if len(folds)==P and all('damping_percent' in x for x in folds) else np.full(nf,np.nan)
        for j in range(nf):
            modalvalues.append(dict(id=r['id'],window=r['window'],model=r['model'],mode_index=j,
                role='broad background' if r['window']=='H300' and j==0 else 'reported modal component',
                common_fixed_us=r.get('common_fixed_us'),frequency_hz=r['frequency_hz'][j],damping_percent=r['damping_percent'][j],
                frequency_jk_se_hz=fse[j],damping_jk_se_pp=zse[j],success=r['success'],jackknife_valid=usable))
        for j,name in enumerate(['common','output1','output2','output3','output4','input1','input2','input3','input4']):
            delays.append(dict(id=r['id'],window=r['window'],model=r['model'],fixed=r.get('common_fixed_us'),
                               component=name,delay_us=d(r)[j],se_us=se[j],valid=usable))
    for rid in old.IDS:
        for kind in 'abcd':
            rr=[r for r in rows if r['id']==rid and r['model']==kind and r['window'] in ('H300','H493') and r['common_fixed_us'] is None]
            if len(rr)==2:
                rr=sorted(rr,key=lambda x:x['window'])
                comparisons.append(dict(id=rid,model=kind,common_difference_us=rr[1]['common_us']-rr[0]['common_us'],
                                        absolute_difference_us=abs(rr[1]['common_us']-rr[0]['common_us']),
                                        success=all(x['success'] for x in rr)))
                raw=[rd(ROOT/'raw'/f'{rid}_{w}_{kind}_M3_full.json') for w in ('H300','H493')]
                rel=np.r_[raw[1]['delay']['output_us'],raw[1]['delay']['input_us']]-np.r_[raw[0]['delay']['output_us'],raw[0]['delay']['input_us']]
                comparisons[-1].update(relative_difference_rms_us=float(np.sqrt(np.mean(rel**2))))
                P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])
                fold_pairs=[]
                for i in range(P):
                    paths=[ROOT/'raw'/f'{rid}_{w}_{kind}_M3_omit{i}.json' for w in ('H300','H493')]
                    if all(p.exists() for p in paths):fold_pairs.append([rd(p) for p in paths])
                get_delay=lambda r:np.r_[r['delay']['common_us'],r['delay']['output_us'],r['delay']['input_us']]
                difference=get_delay(raw[1])-get_delay(raw[0])
                valid=len(fold_pairs)==P and all(x['success'] for pair in fold_pairs for x in pair) and all(x['success'] for x in raw)
                se=old.jk_se([get_delay(pair[1])-get_delay(pair[0]) for pair in fold_pairs]) if len(fold_pairs)==P and all('delay' in x for pair in fold_pairs for x in pair) else np.full(9,np.nan)
                for j,name in enumerate(['common','output1','output2','output3','output4','input1','input2','input3','input4']):
                    cross_components.append(dict(id=rid,model=kind,component=name,
                        H493_minus_H300_us=difference[j],paired_jk_se_us=se[j],
                        complete_paired_folds=len(fold_pairs),jackknife_valid=valid))
    for w in old.SPECS:
        for kind in 'abcd':
            rr=[rd(ROOT/'raw'/f'{rid}_{w}_{kind}_M3_full.json') for rid in old.IDS if (ROOT/'raw'/f'{rid}_{w}_{kind}_M3_full.json').exists()]
            if len(rr)!=6 or not all('frequency_hz' in x for x in rr):continue
            for j in old.SPECS[w]['report']:
                f=np.array([x['frequency_hz'][j] for x in rr]);z=np.array([x['damping_percent'][j] for x in rr])
                dispersion.append(dict(window=w,model=kind,mode_index=j,frequency_mean_hz=f.mean(),frequency_sd_hz=f.std(ddof=1),
                     damping_mean_percent=z.mean(),damping_sd_pp=z.std(ddof=1),
                     paired_frequency_repro_sd_hz=np.sqrt(np.mean((f[3:]-f[:3])**2)/2),
                     paired_damping_repro_sd_pp=np.sqrt(np.mean((z[3:]-z[:3])**2)/2),all_success=all(x['success'] for x in rr)))
    gates=[]
    for rid in old.IDS:
        for w,key in [('H300','QE333'),('H300','QE359'),('H493','QE')]:
            aa=next((r for r in rows if r['id']==rid and r['window']==w and r['model']=='a'),None)
            for kind in ('b','c','d'):
                bb=next((r for r in rows if r['id']==rid and r['window']==w and r['model']==kind),None)
                if aa and bb and aa.get(key) and bb.get(key):
                    ratio=bb[key]/aa[key]
                    gates.append(dict(id=rid,window=w,region=key,model=kind,QE_ratio_to_a=ratio,increase_not_over_20pct=ratio<=1.2,
                                      valid=aa['success'] and bb['success']))
    js(ROOT/'analysis/p6_summary.json',dict(endpoints=rows,delays=delays,cross_window=comparisons,cross_window_components=cross_components,dispersion=dispersion,residual_gates=gates))
    for name,x in [('p6_endpoints.csv',rows),('p6_delays.csv',delays),('p6_cross_window.csv',comparisons),('p6_cross_window_components.csv',cross_components),('p6_dispersion.csv',dispersion),('p6_residual_gates.csv',gates),('p6_modal_parameters.csv',modalvalues)]:write_csv(name,x)
    direct=[];direct_delays=[]
    for p in (ROOT/'raw').glob('*_M1_full.json'):
        if p.name.startswith('p7'):continue
        r=rd(p)
        if 'beta' not in r:continue
        P=int(np.load(SRC/'prepared'/f"{r['id']}.npz")['P'])
        folds=[rd(pp) for i in range(P) if (pp:=ROOT/'raw'/p.name.replace('_full',f'_omit{i}')).exists()]
        getd=lambda x:np.r_[x['delay']['common_us'],x['delay']['output_us'],x['delay']['input_us']]
        se=old.jk_se([getd(x) for x in folds]) if len(folds)==P and all('delay' in x for x in folds) else np.full(9,np.nan)
        valid=len(folds)==P and r['success'] and all(x.get('success') for x in folds)
        for j in old.SPECS[r['window']]['report']:
            direct.append(dict(id=r['id'],window=r['window'],model=r['model'],mode_index=j,
                  frequency_hz=r['frequency_hz'][j],damping_percent=r['damping_percent'][j],
                  common_us=r['delay']['common_us'],common_se_us=se[0],success=r['success'],
                  score=r['solver']['score'],complete_folds=len(folds),failed_folds=sum(not x.get('success') for x in folds),jackknife_valid=valid))
        for j,name in enumerate(['common','output1','output2','output3','output4','input1','input2','input3','input4']):
            direct_delays.append(dict(id=r['id'],window=r['window'],model=r['model'],component=name,delay_us=getd(r)[j],se_us=se[j],valid=valid))
    js(ROOT/'analysis/p6_M1_summary.json',direct);write_csv('p6_M1_endpoints.csv',direct);write_csv('p6_M1_delays.csv',direct_delays)
    direct_cross=[]
    for rid in old.IDS:
        P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])
        for kind in 'abcd':
            parents=[rd(ROOT/'raw'/f'{rid}_{w}_{kind}_M1_full.json') for w in ('H300','H493')]
            pairs=[]
            for i in range(P):
                paths=[ROOT/'raw'/f'{rid}_{w}_{kind}_M1_omit{i}.json' for w in ('H300','H493')]
                if all(p.exists() for p in paths):pairs.append([rd(p) for p in paths])
            getd=lambda r:np.r_[r['delay']['common_us'],r['delay']['output_us'],r['delay']['input_us']]
            delta=getd(parents[1])-getd(parents[0])
            se=old.jk_se([getd(pair[1])-getd(pair[0]) for pair in pairs]) if len(pairs)==P and all('delay' in r for pair in pairs for r in pair) else np.full(9,np.nan)
            valid=len(pairs)==P and all(r['success'] for pair in pairs for r in pair) and all(r['success'] for r in parents)
            for j,name in enumerate(['common','output1','output2','output3','output4','input1','input2','input3','input4']):
                direct_cross.append(dict(id=rid,model=kind,component=name,H493_minus_H300_us=delta[j],
                    paired_jk_se_us=se[j],complete_paired_folds=len(pairs),jackknife_valid=valid))
    write_csv('p6_M1_cross_window_components.csv',direct_cross)
    direct_dispersion=[]
    for w in old.SPECS:
        for kind in 'abcd':
            rr=[rd(ROOT/'raw'/f'{rid}_{w}_{kind}_M1_full.json') for rid in old.IDS]
            for j in old.SPECS[w]['report']:
                f=np.array([r['frequency_hz'][j] for r in rr]);z=np.array([r['damping_percent'][j] for r in rr])
                direct_dispersion.append(dict(window=w,model=kind,mode_index=j,frequency_mean_hz=f.mean(),
                    frequency_sd_hz=f.std(ddof=1),damping_mean_percent=z.mean(),damping_sd_pp=z.std(ddof=1),
                    paired_frequency_repro_sd_hz=np.sqrt(np.mean((f[3:]-f[:3])**2)/2),
                    paired_damping_repro_sd_pp=np.sqrt(np.mean((z[3:]-z[:3])**2)/2),
                    all_success=all(r['success'] for r in rr)))
    write_csv('p6_M1_dispersion.csv',direct_dispersion)
    sensitivity=[]
    for ap in (ROOT/'analysis').glob('fixed_low_*_anchor.json'):
        stem=ap.stem[len('fixed_low_'):-len('_anchor')]
        rid,kind=stem.rsplit('_',1);anchor=rd(ap)['common_fixed_us']
        for method in ('M1','M3'):
            parent=rd(ROOT/'raw'/f'{rid}_L_{kind}_{method}_full_fixed{anchor:.6f}.json')
            for shift in (-25.,0.,25.):
                p=ROOT/'raw'/f'{rid}_L_{kind}_{method}_full_fixed{anchor+shift:.6f}.json'
                if not p.exists():continue
                r=rd(p)
                for j in range(3):
                    sensitivity.append(dict(id=rid,model=kind,method=method,mode_index=j,
                        anchor_us=anchor,shift_us=shift,success=r.get('success',False),
                        anchor_success=parent.get('success',False),
                        frequency_hz=r.get('frequency_hz',[np.nan]*3)[j],
                        damping_percent=r.get('damping_percent',[np.nan]*3)[j],
                        frequency_change_hz=r.get('frequency_hz',[np.nan]*3)[j]-parent.get('frequency_hz',[np.nan]*3)[j],
                        damping_change_pp=r.get('damping_percent',[np.nan]*3)[j]-parent.get('damping_percent',[np.nan]*3)[j]))
    js(ROOT/'analysis/p6_fixed_low_sensitivity.json',sensitivity)
    write_csv('p6_fixed_low_sensitivity.csv',sensitivity)

def low_one(rid,kind='a'):
    anchor=ROOT/'analysis'/f'fixed_low_{rid}_{kind}_anchor.json'
    if not anchor.exists():return
    fixed=rd(anchor)['common_fixed_us'];base=f'{rid}_L_{kind}'
    dest=ROOT/'analysis'/f'low_core_{rid}_{kind}.json'
    if dest.exists():return
    rows=[]
    P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])
    for omit in [None,*range(P)]:
        tag='full' if omit is None else f'omit{omit}'
        pp=[ROOT/'raw'/f'{base}_{met}_{tag}_fixed{fixed:.6f}.json' for met in ('M1','M3')]
        if not all(x.exists() for x in pp):continue
        rr=[rd(x) for x in pp]
        row=dict(omit=omit,success=all(x['success'] for x in rr))
        try:
            m,b=model_from(rr[1]);a=old.spectrum(rid,'L',omit);ext=old.spectrum(rid,'L',omit,extended=True)
            Phi=old.local_phi(m.full,ext,m.expand(b,False)[0],a['f'])
            pred=prediction(m,a,b,Phi)
            ob=100*np.expm1(np.asarray(rr[0]['beta'])[m.nl:m.nl+6]-b[m.nl:m.nl+6])
            row.update(observed=ob,prediction=pred['prediction_percent'],condition=pred['condition_scaled'],
                       M1_damping=rr[0]['damping_percent'],M3_damping=rr[1]['damping_percent'])
        except Exception as e:row.update(exception=repr(e),success=False)
        rows.append(row)
    if not rows:return
    out=dict(id=rid,model=kind,common_fixed_us=fixed,rows=rows,complete=len(rows)==P+1,all_success=all(x['success'] for x in rows))
    if len(rows)==P+1 and all('prediction' in x for x in rows):
        obs=np.array([x['observed'] for x in rows[1:]]);pred=np.array([x['prediction'] for x in rows[1:]])
        se=old.jk_se(obs-pred)
        out.update(observed_se=old.jk_se(obs),prediction_se=old.jk_se(pred),joint_se=se,
                   z=(np.asarray(rows[0]['observed'])-rows[0]['prediction'])/se)
    js(dest,out)

def p7():
    rows=[];files=list((ROOT/'raw').glob('p7_sim_P*_r[0-9][0-9][0-9].json'))
    for P in (4,8,16):
        rr=[rd(p) for p in files if rd(p).get('P')==P]
        for method in ('M2','M3','JK'):
            valid=[x for x in rr if method in x and np.isfinite(x[method]['error_percent']).all()]
            x=np.array([r[method]['error_percent'] for r in valid])
            if not len(x):continue
            for j,name in enumerate(['f1','f2','f3','zeta1','zeta2','zeta3']):
                sd=x[:,j].std(ddof=1) if len(x)>1 else np.nan;se=sd/np.sqrt(len(x))
                paired=[r for r in valid if 'M2' in r and np.isfinite(r['M2']['error_percent']).all()]
                ratio=np.var([r[method]['error_percent'][j] for r in paired],ddof=1)/np.var([r['M2']['error_percent'][j] for r in paired],ddof=1) if len(paired)>1 else np.nan
                failures=sum(not r[method]['success'] for r in valid)+len(rr)-len(valid)
                rows.append(dict(P=P,method=method,parameter=name,completed=len(rr),finite=len(x),required=40,bias_percent=x[:,j].mean(),
                     standard_error=se,sd_percent=sd,rmse_percent=np.sqrt(np.mean(x[:,j]**2)),failures=failures,
                     log_bias_times100=float(np.mean(100*np.log1p(x[:,j]/100))),
                     exp_transform_contribution=float(np.mean(x[:,j]-100*np.log1p(x[:,j]/100))),
                     variance_ratio_to_M2=ratio,below_two_MC_SE=bool(abs(x[:,j].mean())<2*se),
                     valid_primary_gate=bool(len(rr)==40 and failures==0),
                     original_A_failures=sum(r.get(method+'_previous_status') not in (None,'not previously run','pole_step_converged') for r in rr)))
    js(ROOT/'analysis/p7_sim_summary.json',rows);write_csv('p7_sim_summary.csv',rows)

if __name__=='__main__':
    p6()
    for rid in old.IDS:
        try:low_one(rid)
        except Exception as e:js(ROOT/'analysis'/f'low_core_{rid}_ERROR.json',dict(exception=repr(e),traceback=traceback.format_exc()))
    p7()
