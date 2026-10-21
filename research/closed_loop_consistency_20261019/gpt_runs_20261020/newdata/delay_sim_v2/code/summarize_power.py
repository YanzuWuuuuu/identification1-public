from oracle_model import *
import csv

def main():
    rows=[];p3=[]
    for folder in sorted((ROOT/'power').glob('FC*')):
        if not (folder/'summary.json').exists():continue
        _,band,pp,method=folder.name.split('_');fc=int(folder.name.split('_')[0][2:]);P=int(pp[1:]);m=OracleModel(band)
        src=json.loads((folder/'summary.json').read_text(encoding='utf8'))['rows'];good=[r for r in src if r['success']]
        if len(good)!=30:
            rows.append(dict(configuration=folder.name,status='incomplete',n=len(good)));continue
        b=np.array([r['bias'] for r in good]);sd=np.array([r['standard_deviation'] for r in good])
        items=[(i,f'f{FM[a]:g}','%',False) for i,a in zip(range(m.nl,m.nl+m.nm),m.active)]
        items +=[(i,f'zeta{FM[a]:g}','%',True) for i,a in zip(range(m.nl+m.nm,m.nl+2*m.nm),m.active)]
        items +=[(m.nb-7+i,('tau_c','aH1','aH2','aH3','bH1','bH2','bH3')[i],'us',False) for i in range(7)]
        for idx,name,unit,is_zeta in items:
            if unit=='us':pred=b[:,idx]*100;std=sd[:,idx]*100
            else:pred=100*np.expm1(b[:,idx]);std=100*np.exp(b[:,idx])*sd[:,idx]
            mu=float(np.mean(pred));se=float(la.norm(std)/30)
            ratio=abs(mu)/(3*se) if se else np.nan
            rows.append(dict(configuration=folder.name,fc=fc,P=P,band=band,method=method,parameter=name,unit=unit,
                predicted=mu,SE_MC=se,two_SE=2*se,ratio_to_3SE=ratio,
                N_required=(int(np.ceil(9*np.mean(std**2)/mu**2)) if mu else None),
                status=('ratio_power_sufficient' if ratio>=1 else 'ratio_underpowered') if method=='M1' else 'zero_bias_detection_scale',n=30))
        if band=='W' and method=='M1':
            # Eq (5.4) keeps all SIX high-mode pole coordinates. Eliminating
            # high-mode frequencies prematurely changes the propagation split.
            estimates=[];stds=[];audits=[]
            keep=np.r_[m.nl+np.array([3,4,5]),m.nl+m.nm+np.array([3,4,5]),np.arange(m.nb-7,m.nb)]
            rest=np.setdiff1d(np.arange(m.nb),keep)
            for r in range(30):
                with np.load(folder/f'r{r:03d}.npz') as a:
                    M=a['M'];bv=a['b'];scale=a['scales'];delta=a['bias'];cov=a['covariance']
                sn=scale[rest];A=M[np.ix_(rest,rest)]*sn[:,None]*sn[None,:];cross=M[np.ix_(rest,keep)]
                sol=sn[:,None]*la.solve(A,sn[:,None]*cross,assume_a='pos')
                ss=sn*la.solve(A,sn*bv[rest],assume_a='pos')
                H=M[np.ix_(keep,keep)]-cross.T@sol;bb=bv[keep]-cross.T@ss
                direct=la.solve(H[:6,:6],bb[:6]);prop=-la.solve(H[:6,:6],H[:6,6:]@delta[-7:])
                audits.append(dict(r=r,high_parameters=keep[:6],direct_plus_nuisance_log=direct,delay_propagation_log=prop,
                    total_log=delta[keep[:6]],sum_error=float(la.norm(direct+prop-delta[keep[:6]]))))
                estimates.append(100*np.expm1(prop[3:]));stds.append(100*np.exp(delta[keep[3:6]])*np.sqrt(np.diag(cov)[keep[3:6]]))
            save(ROOT/'power/p3_schur_audit'/f'FC{fc}.json',audits)
            estimates=np.array(estimates);stds=np.array(stds)
            for j,f in enumerate(FM[3:]):
                mu=float(np.mean(estimates[:,j]));se=float(la.norm(stds[:,j])/30)
                p3.append(dict(fc=fc,frequency=float(f),predicted_percent=mu,SE_MC_percent=se,ratio_to_3SE=abs(mu)/(3*se),
                    N_required=int(np.ceil(9*np.mean(stds[:,j]**2)/mu**2)) if mu else None,
                    status='ratio_power_sufficient' if abs(mu)>=3*se else 'ratio_underpowered',
                    applicability='primary_Ksmall' if fc==15 else 'extended_control_not_small_loop'))
    save(ROOT/'power/summary_all.json',dict(rows=rows,P3=p3))
    for name,table in [('power_table',rows),('p3_power',p3)]:
        keys=list(dict.fromkeys(k for r in table for k in r))
        with (ROOT/'power'/f'{name}.csv').open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,keys);w.writeheader();w.writerows(table)
    print(json.dumps(p3,default=enc),flush=True)

if __name__=='__main__':main()
