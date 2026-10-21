from common import *
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42})
COL=['#245A81','#258780','#78599A','#666666']
def csvout(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/'analysis'/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def save(fig,name):
    fig.savefig(HERE/'figures'/f'{name}.png',dpi=300,bbox_inches='tight')
    fig.savefig(HERE/'figures'/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)
def main():
    cv=rd(HERE/'analysis/curvature.json');jk=rd(HERE/'analysis/prediction_jackknife.json');ss=rd(HERE/'analysis/sensitivity.json')
    ctab=[];atab=[]
    for r in cv:
        row={k:v for k,v in r.items() if not isinstance(v,(dict,list))}
        row.update(zeta2_pred_M=r['pred_M'][4],zeta2_pred_MminusC=r['pred_MminusC'][4])
        ctab.append(row)
        for g in r['correction_by_source']:
            atab.append(dict(id=r['id'],scenario=r['scenario'],**{k:v for k,v in g.items() if k!='correction_vector'}))
    csvout('curvature_metrics.csv',ctab);csvout('curvature_contributions.csv',atab)
    jt=[]
    for r in jk:
        for key in ('M','MminusC'):
            obs=np.array(r['observed']);pred=np.array(r['pred_'+key])
            so=np.array(r['observed_se']);sp=np.array(r['pred_'+key+'_se']);sd=np.array(r['difference_'+key+'_se'])
            cov=np.array(r['observed_prediction_JK_covariance_'+key])
            for j in range(3):
                k=j+3
                jt.append(dict(id=r['id'],formula=key,parameter='zeta'+str(j+1),observed=obs[k],prediction=pred[k],
                    observed_se=so[k],prediction_se=sp[k],joint_difference_se=sd[k],
                    covariance=cov[k],correlation=cov[k]/(so[k]*sp[k]),z_joint=r['z_joint_'+key][k],
                    original_effect_gate=bool(abs(pred[k])>2*so[k]),
                    prediction_failures=r['prediction_failures'],original_strict_failed_folds=r['original_strict_failures']))
    csvout('prediction_uncertainty.csv',jt);js(HERE/'analysis/prediction_uncertainty_table.json',jt)
    stab=[];summary=[];edge=[]
    for r in ss:
        row=dict(id=r['id'],hi=r['hi'],background=r['background'],parameters=r['nb'])
        for method in ('M1','M3'):
            ff=r['fits'][method]
            row.update({method+'_zeta2':ff.get('damping_percent',[None,None])[1],
                method+'_converged':ff['converged'],method+'_strict':ff.get('strict_score_pass',False),
                method+'_iterations':ff.get('solver_info',{}).get('iterations'),
                method+'_score':ff.get('solver_info',{}).get('score_max'),
                method+'_background_f':ff.get('background_f_hz'),
                method+'_background_zeta':ff.get('background_damping'),
                method+'_bound':str(ff.get('background_bound_active',[]))})
        row.update(paired_zeta2=r.get('paired_percent',[None]*5)[4],
            Q=r.get('residual',{}).get('Q'),EQ=r.get('residual',{}).get('EQ'),Q_over_EQ=r.get('residual',{}).get('ratio'),
            scaled_condition=r.get('residual',{}).get('A_numerics',{}).get('condition'),
            edge_Q=r.get('edge_23_25',{}).get('Q'),edge_raw_residual=r.get('edge_23_25',{}).get('raw_residual_norm'))
        stab.append(row)
        if r['hi']==25:
            b=next(x for x in ss if x['id']==r['id'] and x['hi']==25 and x['background']=='D0D2')
            edge.append(dict(id=r['id'],background=r['background'],relative_common_covariance_Q=r['edge_23_25']['Q']/b['edge_23_25']['Q'],
                relative_raw_residual=r['edge_23_25']['raw_residual_norm']/b['edge_23_25']['raw_residual_norm'],
                M3_converged=r['fits']['M3']['converged'],M3_strict=r['fits']['M3']['strict_score_pass']))
    csvout('sensitivity_all.csv',stab);csvout('common_edge_23_25.csv',edge)
    for bg in ('D0D2','extra40'):
        for hi in (20,25,30):
            for group in ('A','C'):
                rows=[next(x for x in ss if x['hi']==hi and x['background']==bg and x['id']==(group+'_'+c+('_pppp' if group=='C' else ''))) for c in CTRS]
                z1=np.array([r['fits']['M1']['damping_percent'][1] for r in rows])
                z3=np.array([r['fits']['M3']['damping_percent'][1] for r in rows])
                dp=np.array([r['paired_percent'][4] for r in rows])
                summary.append(dict(background=bg,hi=hi,group=group,M1_range_percentage_points=np.ptp(z1),M3_range_percentage_points=np.ptp(z3),
                     M1_monotone=bool(np.all(np.diff(z1)>0)),paired_monotone=bool(np.all(np.diff(dp)>0)),
                     all_converged=all(v['converged'] for x in rows for v in x['fits'].values()),
                     all_strict=all(v['strict_score_pass'] for x in rows for v in x['fits'].values()),
                     M1_values=z1,M3_values=z3,paired_values=dp))
    js(HERE/'analysis/sensitivity_summary.json',summary)
    # Two prediction formulas are both visible; old plot remains untouched.
    fig,axs=plt.subplots(1,2,figsize=(10.8,4.1),sharey=True,layout='constrained')
    jj={r['id']:r for r in jk}
    for ax,key,title in zip(axs,('M','MminusC'),(r'$M^{-1}b$',r'$(M-C)^{-1}b$ (original main prediction)')):
        for group,marker,color,off in [('A','o',COL[0],-.09),('C','s',COL[1],.09)]:
            ids=[group+'_'+c+('_pppp' if group=='C' else '') for c in CTRS];xx=np.arange(3)+off
            yy=[jj[r]['observed'][4] for r in ids];yp=[jj[r]['pred_'+key][4] for r in ids]
            ey=[jj[r]['observed_se'][4] for r in ids];ep=[jj[r]['pred_'+key+'_se'][4] for r in ids]
            ax.errorbar(xx-.025,yy,yerr=ey,marker=marker,color=color,ls='-',ms=4,capsize=2,label=group+' observed')
            ax.errorbar(xx+.025,yp,yerr=ep,marker=marker,color=color,ls='--',mfc='white',ms=5,capsize=2,label=group+' predicted')
        ax.axhline(0,color='.7',lw=.6);ax.set(xticks=range(3),xticklabels=['Low','Nominal','High'],xlabel='Controller gain',title=title)
        ax.grid(axis='y',alpha=.15)
    axs[0].set_ylabel(r'Mode 2 paired damping difference (%)')
    axs[0].legend(ncol=2,fontsize=7)
    fig.suptitle('Both observed and predicted values: period jackknife ±1 SE',fontsize=11)
    save(fig,'main_paired_prediction')
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    labs=[r['id'].replace('_pppp','').replace('gain','')+(' scaled' if r['scenario']=='scaled' else '') for r in cv]
    xx=np.arange(len(cv))
    axs[0].semilogy(xx,[r['raw_H_norm2'] for r in cv],'o-',color=COL[3],label='Raw Euclidean norm')
    axs[0].semilogy(xx,[r['max_abs_eigenvalue'] for r in cv],'s-',color=COL[0],label='M-metric norm q')
    axs[0].axhline(1,color='.5',ls='--');axs[0].legend(fontsize=7);axs[0].set_ylabel('Curvature norm')
    axs[1].plot(xx,[r['relative_correction_in_M_metric']*100 for r in cv],'o-',color=COL[1],label='Whole vector in M metric')
    axs[1].plot(xx,[r['zeta2_relative_change_between_predictions']*100 for r in cv],'s-',color=COL[0],label='Mode 2 percentage prediction')
    axs[1].set_ylabel('Relative change after curvature correction (%)');axs[1].legend(fontsize=7)
    for ax in axs:
        ax.set_xticks(xx,labs,rotation=60,ha='right',fontsize=7);ax.grid(axis='y',alpha=.15)
    save(fig,'curvature_scale_and_component_sensitivity')
    fig,axs=plt.subplots(2,2,figsize=(9.4,6.5),layout='constrained')
    for row,group in enumerate(('A','C')):
        for hi,color in zip((20,25,30),COL):
            x=next(r for r in summary if r['background']=='D0D2' and r['hi']==hi and r['group']==group)
            axs[row,0].plot(range(3),x['M1_values'],'o--',color=color,ms=4,label=f'M1, 3–{hi} Hz')
            axs[row,0].plot(range(3),x['M3_values'],'s-',color=color,ms=4,label=f'M3, 3–{hi} Hz')
            axs[row,1].plot(range(3),x['paired_values'],'o-',color=color,ms=4,label=f'3–{hi} Hz')
        for ax in axs[row]:
            ax.set_xticks(range(3),['Low','Nominal','High']);ax.grid(alpha=.15)
        axs[row,0].set(title=group+' group: D0+D2 background',ylabel=r'Mode 2 damping (%)')
        axs[row,1].set(title=group+' group: paired difference',ylabel='Relative paired difference (%)')
    axs[0,0].legend(fontsize=6.5,ncol=2);axs[0,1].legend(fontsize=7)
    save(fig,'frequency_band_sensitivity')
    fig,axs=plt.subplots(1,2,figsize=(10.5,3.8),layout='constrained')
    for i,rid in enumerate(IDS):
        b=next(r for r in ss if r['id']==rid and r['hi']==25 and r['background']=='D0D2')
        e=next(r for r in ss if r['id']==rid and r['hi']==25 and r['background']=='extra40')
        axs[0].semilogy(i-.1,b['residual']['ratio'],'o',color=COL[0])
        axs[0].semilogy(i+.1,e['residual']['ratio'],'x',color=COL[3],ms=7)
        er=next(r for r in edge if r['id']==rid and r['background']=='extra40')
        axs[1].plot(i,er['relative_common_covariance_Q'],'x',color=COL[3],ms=7)
    axs[0].plot([],[],'o',color=COL[0],label='Original background')
    axs[0].plot([],[],'x',color=COL[3],label='Extra background: unconverged')
    axs[0].axhline(2,color='.6',ls='--');axs[0].legend(fontsize=7);axs[0].set_ylabel('Within-record Q / E[Q]')
    axs[1].axhline(1,color='.5',ls='--');axs[1].set_ylabel('23–25 Hz Q: extra / original (fixed covariance)')
    for ax in axs:
        ax.set_xticks(range(6),[r.replace('_pppp','').replace('gain','') for r in IDS],rotation=35,ha='right',fontsize=8)
        ax.grid(alpha=.15)
    fig.suptitle('Background sensitivity at 3–25 Hz: failed endpoints are diagnostic only',fontsize=10)
    save(fig,'background_failure_diagnostic')
if __name__=='__main__':main()
