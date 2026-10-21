from common import *
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':9,'axes.spines.right':False,'axes.spines.top':False,'legend.frameon':False,
 'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.8})
COLORS={'lowgain':'#28618A','nominal':'#248D87','highgain':'#7C609E'}
LABEL={'lowgain':'Low','nominal':'Nominal','highgain':'High'}
def readcsv(name):
    with (HERE/'analysis'/name).open(encoding='utf-8-sig') as h:return list(csv.DictReader(h))
def num(r,key):return float(r[key])
def save(fig,name):
    fig.savefig(HERE/'figures'/f'{name}.png',dpi=300,bbox_inches='tight')
    fig.savefig(HERE/'figures'/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)
def main():
    high=readcsv('high_frequency_tests.csv');low=readcsv('low_frequency_all_modes.csv')
    bw=rd(HERE/'analysis/bandwidth.json')
    for formula,label in [('M',r'$M^{-1}b$'),('MminusC',r'$(M-C)^{-1}b$')]:
        fig,axs=plt.subplots(2,1,figsize=(9.2,6.6),sharex=True)
        for gi,g in enumerate(('A','C')):
            ax=axs[gi]
            for bt in bw:
                color=COLORS[bt['controller']]
                ax.axvspan(bt['lower_hz'],bt['upper_hz'],color=color,alpha=.10,lw=0)
                ax.axvline(bt['upper_hz'],color=color,alpha=.6,lw=.7,ls=':')
            for ct in CTRS:
                rr=[r for r in low+high if r['formula']==formula and r['parameter']=='zeta' and r['group']==g and r['controller']==ct]
                rr.sort(key=lambda x:num(x,'frequency_hz'));x=np.array([num(r,'frequency_hz') for r in rr])
                for mode,mk,fac in [('observed','o',1.),('predicted','s',1.)]:
                    y=[num(r,mode+'_percent') for r in rr];se=[num(r,mode+'_SE_percent') for r in rr]
                    ax.errorbar(x*fac,y,yerr=se,color=COLORS[ct],marker=mk,ms=4,linestyle='none',
                        mfc=COLORS[ct] if mode=='observed' else 'white',mew=1,capsize=2,lw=.9)
                for r in rr:
                    if r.get('eligible')=='False':
                        ax.plot(num(r,'frequency_hz'),num(r,'observed_percent'),marker='x',color='black',ms=6,linestyle='none')
            ax.axhline(0,color='#777777',lw=.7)
            ax.set_xscale('log');ax.set_yscale('symlog',linthresh=.3,linscale=.7)
            ax.set_xlim(9,550);ax.set_ylabel('Paired damping difference (%)\nSigned symlog; linear within ±0.3%')
            ax.grid(axis='y',alpha=.16);ax.set_title(f'{chr(97+gi)}  '+('A records' if g=='A' else 'C++++ records'),loc='left',fontweight='bold')
            ax.text(.99,.95,'Shading: loop unity-crossing reference',transform=ax.transAxes,ha='right',va='top',fontsize=8,color='#555555')
        handles=[Line2D([],[],color=COLORS[c],lw=2,label=LABEL[c]+' gain') for c in CTRS]
        handles+=[Line2D([],[],color='#444444',marker='o',linestyle='none',label='Observed'),
                  Line2D([],[],color='#444444',marker='s',mfc='white',linestyle='none',label='Predicted')]
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.935),ncol=5)
        axs[-1].set_xticks([10,15,20,30,50,100,282,359,493],labels=['10','15','20','30','50','100','282','359','493'])
        axs[-1].set_xlabel('Mode frequency (Hz; log axis)')
        fig.suptitle('Measured and predicted closed-loop differences: '+label+'\nAll error bars: period jackknife ±1 SE',y=.995,fontsize=11)
        fig.tight_layout(rect=(0,0,1,.89))
        save(fig,'paired_damping_frequency_'+formula)
    fig,axs=plt.subplots(2,2,figsize=(10.2,6.6))
    for i,g in enumerate(('A','C')):
        for j,par in enumerate(('f','zeta')):
            ax=axs[i,j]
            for ci,ct in enumerate(CTRS):
                rr=[r for r in high if r['group']==g and r['controller']==ct and r['parameter']==par and r['formula']=='MminusC']
                rr.sort(key=lambda x:num(x,'frequency_hz'))
                x=np.arange(4)+(ci-1)*.22
                for typ,mk,dx in [('observed','o',-.035),('predicted','s',.035)]:
                    ax.errorbar(x+dx,[num(r,typ+'_percent') for r in rr],yerr=[num(r,typ+'_SE_percent') for r in rr],
                        color=COLORS[ct],marker=mk,mfc=COLORS[ct] if typ=='observed' else 'white',
                        linestyle='none',ms=4,capsize=2,lw=.9)
            ax.axhline(0,c='gray',lw=.7);ax.grid(axis='y',alpha=.15)
            ax.set_xticks(range(4),['282','333','359','493']);ax.set_xlabel('Mode (approximate frequency, Hz)')
            ax.set_ylabel('Paired '+('frequency' if par=='f' else 'damping')+' difference (%)')
            ax.set_title(('A' if g=='A' else 'C++++')+' — '+('frequency' if par=='f' else 'damping'),loc='left')
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.945),ncol=5)
    fig.suptitle('High-frequency close-up: curvature prediction, linear vertical axes; ±1 SE',y=.995,fontsize=11)
    fig.tight_layout(rect=(0,0,1,.89))
    save(fig,'high_frequency_closeup')
    rr=readcsv('reproducibility.csv')
    fig,axs=plt.subplots(1,2,figsize=(9,3.6),constrained_layout=True)
    for ax,par in zip(axs,('f','zeta')):
        rows=[r for r in rr if r['parameter']==par]
        x=np.arange(len(rows))
        for key,mk,color,lab in [('replication_SD','o','#28618A','A/C pooled repeat SD'),
              ('conditional_same_record_u_RMS','s','#248D87','Same-record noise u'),
              ('conditional_B_u_RMS','^','#999999','B-record noise u')]:
            ax.plot(x,[num(r,key) for r in rows],marker=mk,color=color,label=lab,lw=1)
        ax.set_yscale('log');ax.set_xticks(x,['282','333','359','493']);ax.set_xlabel('Mode (Hz)')
        ax.set_ylabel('Frequency SD / u (Hz)' if par=='f' else 'Damping SD / u (percentage points)')
        ax.grid(axis='y',alpha=.15)
    axs[0].legend(fontsize=8);fig.suptitle('Between-run reproducibility vs conditional noise uncertainty\nThree paired runs across controllers; pooled repeat SD has 3 degrees of freedom',fontsize=10)
    save(fig,'reproducibility_vs_noise')
if __name__=='__main__':main()
