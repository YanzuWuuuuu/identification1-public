from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,
    'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,
    'axes.linewidth':.7,'savefig.facecolor':'white'})
COL=['#245A81','#258780','#78599A','#666666']
audit=[]
def save(fig,name):
    fig.savefig(HERE/'figures'/f'{name}.png',dpi=300,bbox_inches='tight')
    fig.savefig(HERE/'figures'/f'{name}.pdf',bbox_inches='tight')
    fig.canvas.draw()
    audit.append(dict(name=name,figure_size_inches=list(fig.get_size_inches()),dpi=300,
                      axes=len(fig.axes),png_bytes=(HERE/'figures'/f'{name}.png').stat().st_size))
    plt.close(fig)
def short(rid):return rid.replace('gain','').replace('_pppp',' ++++').replace('_',' ')
def main():
    jk={r['id']:r for r in rd(HERE/'analysis/jackknife.json')}
    br={r:rd(HERE/'analysis'/f'{r}_noise_bias.json') for r in IDS}
    fig,ax=plt.subplots(figsize=(7.4,4.1),layout='constrained')
    for g,marker,off in [('A','o',-.07),('C','s',.07)]:
        ids=['A_'+c if g=='A' else 'C_'+c+'_pppp' for c in CTRS]
        for method,color,line in [('M1',COL[3],'--'),('M3',COL[0],'-')]:
            yy=[fit(r,method)['damping_percent'][1] for r in ids]
            se=[jk[r][method+'_se'][4] for r in ids]
            ax.errorbar(np.arange(3)+off,yy,yerr=se,marker=marker,color=color,ls=line,
                        lw=1.3,ms=5,capsize=3,label=f'{method}, {g}'+(' ++++' if g=='C' else ''))
    ax.set(xticks=range(3),xticklabels=['Low gain','Nominal','High gain'],ylabel=r'Mode 2 damping, $\zeta_2$ (%)')
    ax.grid(axis='y',alpha=.17);ax.legend(ncol=2,loc='upper left')
    fig.suptitle('Within-record identification; error bars: period jackknife SE',fontsize=10)
    save(fig,'main_damping_mode2')
    fig,ax=plt.subplots(figsize=(7.4,4.1),layout='constrained')
    for g,marker,color,off in [('A','o',COL[0],-.07),('C','s',COL[1],.07)]:
        ids=['A_'+c if g=='A' else 'C_'+c+'_pppp' for c in CTRS];xx=np.arange(3)+off
        obs=[br[r]['observed'][4] for r in ids];se=[jk[r]['paired_se'][4] for r in ids]
        pred=[br[r]['windows']['0.5']['pred_MminusC_percent'][4] for r in ids]
        ax.errorbar(xx,obs,yerr=se,marker=marker,color=color,ls='-',ms=5,capsize=3,label=f'{g}: observed, ±1 SE')
        ax.plot(xx,pred,marker=marker,color=color,ls='--',mfc='white',ms=6,label=f'{g}: predicted')
    ax.axhline(0,color='.65',lw=.7);ax.set(xticks=range(3),xticklabels=['Low gain','Nominal','High gain'],
         ylabel=r'Paired damping difference, $100(\zeta_{M1}/\zeta_{M3}-1)$ (%)')
    ax.grid(axis='y',alpha=.17);ax.legend(ncol=2,loc='upper left')
    fig.suptitle('Same-record noise prediction; fixed 0.5 Hz smoothing window',fontsize=10)
    save(fig,'main_paired_prediction')
    for rid in IDS:
        n=np.load(HERE/'spectra'/f'{rid}_noise.npz');f=n['f']
        fig,axs=plt.subplots(3,2,figsize=(8.4,7.2),sharex=True,layout='constrained')
        for i,ax in enumerate(axs.flat):
            ax.axhline(1,color='.6',lw=.7,ls=':')
            if i<5:
                for key,lab,color in [('0.5','0.5 Hz',COL[0]),('1.0','1.0 Hz',COL[1])]:
                    rat=n['ratio_diag_'+key][:,i] if i<4 else n['ratio_trace_'+key]
                    ax.semilogy(f,rat,color=color,lw=1,label=lab)
                ax.set_title(f'Output {i+1}: diagonal ratio' if i<4 else 'Trace ratio')
            else:
                for d,color in zip(range(3),COL):
                    ax.semilogy(f,n[f'ratio_trace_drop{d}_frozen'],color=color,lw=.85,label=f'Drop {d}: fixed transient')
                    ax.semilogy(f,n[f'ratio_trace_drop{d}_reboundary'],color=color,lw=.9,ls='--',label=f'Drop {d}: refit transient')
                ax.set_title('Trace: extra cycle removal')
            ax.set_xlim(3,25);ax.grid(alpha=.13);ax.set_ylabel('Same-record / B noise')
        axs[0,0].legend();axs[2,1].legend(fontsize=6.5,ncol=1)
        for ax in axs[-1]:ax.set_xlabel('Frequency (Hz)')
        fig.suptitle(short(rid)+' — nonexcited-bin residual spectra',fontsize=11)
        save(fig,rid+'_noise_ratios')
    rr=rd(HERE/'analysis/residual_checks.json')
    fig,axs=plt.subplots(2,3,figsize=(10.4,5.9),sharex=True,layout='constrained')
    for rid,ax in zip(IDS,axs.flat):
        for method,color in [('M3',COL[0]),('M3w',COL[1])]:
            r=next(x for x in rr if x['id']==rid and x['method']==method)
            bins=r['frequency_bands'];xx=[(x['lo']+x['hi'])/2 for x in bins]
            ax.semilogy(xx,[x['ratio'] for x in bins],marker='o' if method=='M3' else 's',ms=3,
                        color=color,label=method+f" (total {r['ratio']:.1f})")
        ax.axhline(2,color='.45',ls='--',lw=.8);ax.set_title(short(rid));ax.legend(fontsize=7)
        ax.set(ylabel=r'$Q_{\rm band}/E[Q_{\rm band}]$',xlabel='Frequency (Hz)',xlim=(3,25))
        ax.grid(alpha=.15)
    fig.suptitle('Within-record residual checks — 2 Hz bins; dashed line: threshold 2',fontsize=11)
    save(fig,'residual_frequency_bands')
    cross=np.load(HERE/'raw/lowgain_crossrecord.npz');f=cross['f']
    fig,axs=plt.subplots(1,3,figsize=(10.5,3.4),layout='constrained')
    axs[0].semilogy(f,cross['q']/32,color=COL[0],lw=1)
    axs[0].axvspan(15,17,color=COL[1],alpha=.12)
    axs[0].set(xlabel='Frequency (Hz)',ylabel='Cross-record Q / expected per bin',title='A model to four C runs')
    for rid,lab,color in [('A_lowgain','A',COL[0]),('C_lowgain_pppp','C ++++',COL[1])]:
        fq=cross[rid+'_f'];psd=cross[rid+'_psd'].sum(1)
        axs[1].semilogy(fq,psd,color=color,lw=.8,label=lab)
        mask=fq<=6;axs[2].semilogy(fq[mask],psd[mask],color=color,lw=.8,label=lab)
    axs[1].set(xlabel='Frequency (Hz)',ylabel='Sum output PSD (recorded units²/Hz)',title='Low-gain output spectra')
    axs[2].set(xlabel='Frequency (Hz)',ylabel='Sum output PSD (recorded units²/Hz)',title='3–6 Hz detail')
    for ax in axs:ax.grid(alpha=.15)
    axs[1].legend();axs[2].legend()
    save(fig,'lowgain_crossrecord_diagnostic')
    sc=rd(HERE/'analysis/scaling.json');amp=sc['amplitude_records']
    fig,axs=plt.subplots(1,2,figsize=(8.5,3.6),layout='constrained')
    for i,color in enumerate(COL):
        axs[0].plot([r['gain'] for r in amp],[r['response_over_linear'][i] for r in amp],
                    'o-',color=color,ms=3,label=f'Output {i+1}')
    axs[0].axhline(1,color='.5',ls='--');axs[0].set(xlabel='Disturbance coefficient',ylabel='Measured RMS / linear prediction',title='D-group response scaling')
    axs[0].legend(ncol=2)
    for k,name in enumerate(('f1','f2','f3')):
        base=amp[-1]['frequency_hz'][k]
        vals=[r['frequency_hz'][k]-base for r in amp]
        # Keep failed estimate visibly separate instead of connecting it to valid results.
        axs[1].plot([r['gain'] for r in amp[1:]],vals[1:],'o-',color=COL[k],ms=4,label=name)
        axs[1].plot(amp[0]['gain'],vals[0],marker='x',color=COL[k],ms=6)
    axs[1].set(xlabel='Disturbance coefficient',ylabel='Frequency shift from gain 1.5 (Hz)',title='Support modes (crosses: failed M3 fit)')
    axs[1].legend()
    for ax in axs:ax.grid(alpha=.15)
    save(fig,'excitation_scaling_check')
    js(HERE/'analysis/figure_export_audit.json',audit)
if __name__=='__main__':main()
