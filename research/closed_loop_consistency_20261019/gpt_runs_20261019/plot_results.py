"""Figures from saved numerical summaries; no fitting or hidden selection."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
data=json.loads((R/'experiment_a/summary.json').read_text())
names={'M1':'Direct: all bins','M2':'Direct: excited lines','M3':'Reference IV: correct K','M4':'Reference IV: input D'}
colors={'M1':'#b43d40','M2':'#497b39','M3':'#2366a1','M4':'#7a4b9f'}
fig,ax=plt.subplots(2,3,figsize=(13.6,7.3),sharex=True)
for m in names:
    rr=[r for r in data if r['label']==m];pp=np.array([r['P'] for r in rr])
    for q in range(3):
        y=[r['all_finite']['mean'][3+q] for r in rr];e=[r['all_finite']['se'][3+q] for r in rr]
        ax[0,q].errorbar(pp,y,yerr=e,color=colors[m],marker='o',capsize=3,lw=1.5,label=names[m])
        ax[1,q].plot(pp,[r['all_finite']['rmse'][3+q] for r in rr],color=colors[m],marker='o',lw=1.5)
for q in range(3):
    ax[0,q].set_title(f'Support mode {q+1} ({[10.6,11.7,16.1][q]} Hz)')
    ax[0,q].axhline(0,color='.5',lw=.8,ls='--')
    ax[1,q].set_xlabel('Recorded periods P')
for a in ax.flat:
    a.set_xscale('log',base=2);a.set_xticks([1,2,4,8,16],labels=['1','2','4','8','16']);a.grid(alpha=.2)
ax[0,0].set_ylabel('Mean damping error (%)');ax[1,0].set_ylabel('Damping RMSE (%)')
fig.legend(*ax[0,0].get_legend_handles_labels(),loc='upper center',ncol=4,bbox_to_anchor=(.5,1.01),frameon=False)
fig.suptitle('Paired Monte Carlo: 40 trials per P; 20 at P=16\nAll finite endpoints retained; error bars = standard error (not standard deviation)',y=.94,fontsize=12)
fig.tight_layout(rect=(0,0,1,.88))
for ext in ['png','pdf','svg']:fig.savefig(R/f'experiment_A_period_scan.{ext}',bbox_inches='tight')
plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4))
ps=[1,2,4,8,16];x=np.arange(5)
for j,m in enumerate(names):
    rr=[r for r in data if r['label']==m];v=[r['failures']/r['n']*100 for r in rr]
    bars=ax.bar(x+(j-1.5)*.19,v,width=.18,color=colors[m],label=names[m])
    for bar,r in zip(bars,rr):
        if r['failures']:ax.annotate(f"{r['failures']}/{r['n']}",(bar.get_x()+bar.get_width()/2,bar.get_height()),ha='center',va='bottom',fontsize=8)
ax.set_xticks(x,labels=ps);ax.set_xlabel('Recorded periods P');ax.set_ylabel('Iteration-limit / numerical failure (%)');ax.set_ylim(0,60)
ax.legend(frameon=False,ncol=2);ax.grid(axis='y',alpha=.2);fig.tight_layout()
fig.savefig(R/'experiment_A_failure_rates.png',bbox_inches='tight');plt.close(fig)
print('saved A figures')
