from run import *
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter,MaxNLocator
from analyze import loadfit

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':9,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False})
COLORS={'M1':'#81708e','M2':'#444444','M3':'#1769aa','M3w':'#208074','M5n':'#b43b77'}
LABELS={'M1':'M1: direct, all bins','M2':'M2: direct, lines','M3':'M3: reference IV','M3w':'M3w: weighted IV','M5n':'M5n: local-cov SML'}
CTRL=('lowgain','nominal','highgain')
def save(fig,name,preview=False):
    fig.savefig(HERE/'figures'/f'{name}.png',dpi=300,bbox_inches='tight',facecolor='white')
    fig.savefig(HERE/'figures'/f'{name}.pdf',bbox_inches='tight',facecolor='white')
    if preview:fig.savefig(HERE/'figures'/f'{name}_preview.png',dpi=120,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def parameters():
    rows=list(csv.DictReader((HERE/'analysis/modal_parameters.csv').open(encoding='utf-8-sig')))
    for w,s in SPECS.items():
        n=len(s['f']);fig,axes=plt.subplots(n,2,figsize=(10,2.3*n+1.0),squeeze=False)
        for j in range(n):
            for co,par in enumerate(('frequency_hz','damping_percent')):
                ax=axes[j,co]
                for mi,method in enumerate(METHODS):
                    for ci,ctrl in enumerate(CTRL):
                        for gi,group in enumerate(('A','C')):
                            r=next(r for r in rows if r['window']==w and r['controller']==ctrl and r['method']==method and r['record_group']==group and int(r['anchor_component'])==j+1)
                            bad=r['converged']!='True' or r['near_bound']=='True';match=r['candidate_match']=='True'
                            marker='x' if bad else ('o' if gi==0 else '^')
                            alpha=1 if r['strict_score_pass']=='True' else .45
                            ax.plot(ci+(mi-2)*.13+(gi-.5)*.04,float(r[par]),marker=marker,ls='none',ms=5,
                                color=COLORS[method],alpha=alpha,mfc=COLORS[method] if gi==0 or bad else 'white',mew=1)
                            if not match:ax.plot(ci+(mi-2)*.13+(gi-.5)*.04,float(r[par]),marker='s',ls='none',ms=8,mfc='none',mec='#888888',mew=.5)
                ax.set_xticks(range(3),['Low gain','Nominal','High gain']);ax.grid(axis='y',alpha=.2)
                ax.set_ylabel('Frequency (Hz)' if co==0 else r'Damping $\zeta$ (%)');ax.ticklabel_format(useOffset=False,style='plain',axis='y')
                ax.set_title(f"Component {j+1}"+(' (background only)' if s['background'][j] else ''),loc='left',fontsize=10)
                ax.spines[['top','right']].set_visible(False)
        handles=[Line2D([],[],color=COLORS[a],marker='o',ls='none',label=LABELS[a]) for a in METHODS]
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,1.02),ncol=3,fontsize=9)
        fig.suptitle(f"{s['band'][0]}–{s['band'][1]} Hz: parameters from independent A / C records",y=1.075,fontsize=13)
        fig.text(.5,.005,'Circle: A; triangle: C++++; cross: solver failure or bound. Faint: strict score not passed.\nSquare outline: residue/frequency correspondence rejected. Points are individual estimates, not confidence intervals.',ha='center',fontsize=9)
        fig.tight_layout(rect=[0,.10 if n==1 else .045,1,1]);save(fig,f'{w}_parameters',True)

def bode():
    for ctrl in CTRL:
        dat=np.load(PARENT/'wideband_scope_audit'/f'{ctrl}_broadband.npz')
        for w,s in SPECS.items():
            _,m,_,_=setup(w);mask=(dat['f']>=s['band'][0])&(dat['f']<=s['band'][1]);f=dat['f'][mask]
            r=loadfit('A_'+ctrl,w,'M3');GG=m.resp(m.context(f),np.asarray(r['beta']));R=dat['G'][mask];halfdiff=dat['G_front'][mask]-dat['G_back'][mask]
            status=('solver converged' if r['converged'] else 'FAILED solver')+('; boundary endpoint' if r['near_bound'] else '')+('; strict score passed' if r['strict_score_pass'] else '; strict score NOT passed')
            for kind in ('magnitude','phase','residual'):
                fig,axs=plt.subplots(4,4,figsize=(13.2,9.7),sharex=True)
                for o in range(4):
                    for i in range(4):
                        ax=axs[o,i];rr=R[:,o,i];gg=GG[:,o,i]
                        if kind=='phase':
                            pp=np.rad2deg(np.unwrap(np.angle(rr)));qq=np.rad2deg(np.unwrap(np.angle(gg)))
                            pp+=360*np.round(np.median(qq-pp)/360)
                            ax.plot(f,pp,color='#a0a0a0',lw=.7);ax.plot(f,qq,color=COLORS['M3'],lw=1.1)
                        else:
                            db=lambda x:20*np.log10(np.maximum(abs(x),1e-300))
                            ax.plot(f,db(rr),color='#a0a0a0',lw=.8)
                            ax.plot(f,db(gg if kind=='magnitude' else gg-rr),color=COLORS['M3'],lw=1.1)
                            if kind=='residual':ax.plot(f,db(halfdiff[:,o,i]),color='#3c8a6d',lw=.6,alpha=.8)
                        ax.set_title(rf'$u_{i+1}\to y_{o+1}$',fontsize=10);ax.grid(which='both',alpha=.20,lw=.5)
                        ax.set_xscale('log');ax.set_xlim(*s['band']);ticks=np.linspace(*s['band'],5);ax.set_xticks(ticks);ax.xaxis.set_major_formatter(ScalarFormatter());ax.minorticks_off()
                        if o==3:ax.set_xlabel('Frequency (Hz)')
                        if i==0:ax.set_ylabel('Phase (deg)' if kind=='phase' else 'Magnitude (dB)')
                labels=['Independent C FRF estimate','A reference-IV model' if kind!='residual' else '|Model − C FRF|']
                colors=['#a0a0a0',COLORS['M3']]
                if kind=='residual':labels+=['|C front − C back| (repeatability diagnostic)'];colors+=['#3c8a6d']
                fig.legend([Line2D([],[],color=c,lw=1.4) for c in colors],labels,loc='upper center',bbox_to_anchor=(.5,.96),ncol=3)
                fig.suptitle(f"{ctrl}: {s['band'][0]}–{s['band'][1]} Hz | {status}",y=.997,fontsize=12)
                foot='FRF is an independent estimate, not ground truth. All channels and frequencies retained; record units.'
                if kind=='phase':foot+=' Independently unwrapped; one constant 360° shift per curve.'
                fig.text(.5,.005,foot,ha='center',fontsize=8)
                fig.tight_layout(rect=[0,.02,1,.94]);save(fig,f'{ctrl}_{w}_{kind}',ctrl=='nominal' and w in ('H493','H300'))
                print(ctrl,w,kind,flush=True)

def amplitude():
    rows=list(csv.DictReader((HERE/'analysis/amplitude_diagnostic.csv').open(encoding='utf-8-sig')))
    fig,axs=plt.subplots(2,2,figsize=(9,6))
    for ax,(w,j) in zip(axs.flat,[('H300',2),('H300',3),('H300',4),('H493',1)]):
        rr=sorted([r for r in rows if r['window']==w and int(r['component'])==j],key=lambda r:float(r['disturb_gain']))
        ar=loadfit('A_nominal',w,'M3');cr=loadfit('C_nominal_pppp',w,'M3')
        xx=[float(r['disturb_gain']) for r in rr]+[1.5];ff=[float(r['f_hz']) for r in rr]+[ar['frequency_hz'][j-1]]
        ax.plot(xx,ff,'o-',color=COLORS['M3'],lw=1.1,label='D scan + A at gain 1.5')
        ax.plot(1.5,cr['frequency_hz'][j-1],'^',mfc='white',mec='#208074',ms=7,label='Independent C++++ at gain 1.5')
        ax.set_xlabel('Injected disturbance gain');ax.set_ylabel('Identified frequency (Hz)');ax.set_xticks([.6,.8,1.,1.2,1.5]);ax.grid(alpha=.2)
        ax.ticklabel_format(useOffset=False,style='plain',axis='y');ax.set_title(f'Candidate near {ar["frequency_hz"][j-1]:.1f} Hz')
        ax.spines[['top','right']].set_visible(False)
    handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.95),ncol=2)
    fig.suptitle('Same controller and position: amplitude scan (reference IV)',fontsize=13,y=.99)
    fig.text(.5,.008,'One record per gain; lines guide the eye. Trend does not by itself establish a physical nonlinear mechanism.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.03,1,.90]);save(fig,'amplitude_frequency_check',True)

if __name__=='__main__':parameters();bode();amplitude()
