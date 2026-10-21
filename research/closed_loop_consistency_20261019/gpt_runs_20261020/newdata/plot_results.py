from evaluate_results import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':10,'axes.titlesize':10,'legend.fontsize':8,
                     'pdf.fonttype':42,'savefig.dpi':300,'axes.linewidth':.7,'legend.frameon':False})
COL={'M1':'#8648a3','M2':'#ad9cbb','M3':'#0862b9','M3w':'#268366','M5n':'#4d5963'}
CC={'lowgain':'#0862b9','nominal':'#268366','highgain':'#8648a3'}
NAMES={'M1':'Direct, all bins','M2':'Direct, lines','M3':'Reference IV','M3w':'Weighted reference IV','M5n':'Local-noise SML'}
manifest=[]
def save(fig,name,source):
    fig.savefig(HERE/'figures'/f'{name}.png',dpi=300,bbox_inches='tight')
    fig.savefig(HERE/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
    manifest.append(dict(name=name,source=source,png_dpi=300,editable_pdf=True))
def linegrid(ax):ax.grid(True,which='both',color='#d9d9d9',alpha=.6,lw=.45)

def phase_branch(g, anchor):
    """Independent unwrap, then ONE integer 360-degree translation per channel.

    No frequency-dependent branch adjustment or model-based smoothing is used.
    The median avoids choosing a global branch from one noisy endpoint.
    """
    v=np.unwrap(np.angle(g),axis=0)*180/np.pi
    return v+360*np.round(np.median(anchor-v,axis=0)/360)

def controller_figs():
    for band in ((1,500),(3,25)):
        fig,axs=plt.subplots(2,2,figsize=(8.2,5.3),layout='constrained')
        specs=[('singular_K','Controller K (3 nonzero singular values)',3),('singular_loop','Input loop K G (3 nonzero)',3),
               ('singular_Su','Input sensitivity S_u (4 singular values)',4),('singular_F','Feedback path S_u K (3 nonzero)',3)]
        for ctrl in CTRS:
            a=np.load(HERE/'analysis'/f'{ctrl}_wide_frf.npz');sl=(a['f']>=band[0])&(a['f']<=band[1])
            for ax,(key,title,nsv) in zip(axs.flat,specs):
                for j in range(nsv):ax.semilogx(a['f'][sl],20*np.log10(np.maximum(a[key][sl,j],1e-300)),color=CC[ctrl],ls=('-', '--',':','-.')[j],lw=.9,label=ctrl if j==0 else None)
                ax.set_title(title);ax.set_xlabel('Frequency (Hz)');ax.set_ylabel('Singular value (dB)');linegrid(ax)
        axs[0,0].legend(ncol=3);fig.suptitle('Measured controller / loop diagnostics; G from orthogonal C records',fontsize=11)
        save(fig,f'controller_singular_{band[0]}_{band[1]}','controllers/*.npz; analysis/*_wide_frf.npz')
    fig,axs=plt.subplots(2,2,figsize=(8,5),layout='constrained')
    for ctrl in CTRS:
        a=np.load(HERE/'analysis'/f'{ctrl}_noise.npz')
        for j,ax in enumerate(axs.flat):
            ax.semilogy(a['f'],a['Phi'][:,j,j].real/FS,color=CC[ctrl],label=ctrl);ax.set_title(f'Equivalent output noise v{j+1}');ax.set_xlabel('Frequency (Hz)');ax.set_ylabel('Two-sided PSD (recorded units²/Hz)');linegrid(ax)
    axs[0,0].legend();save(fig,'B_equivalent_noise_psd','analysis/*_noise.npz')

def core():
    for indices,name in (([3,4,1],'core_controller_zeta1_zeta2_f2'),(list(range(6)),'all_modal_controller')):
        n=len(indices);fig,axs=plt.subplots(1 if n==3 else 2,3,figsize=(9,3 if n==3 else 5.6),layout='constrained')
        for ax,j in zip(axs.flat,indices):
            for mi,method in enumerate(('M1','M3','M3w')):
                r=[fit('A_'+c,method) for c in CTRS];u=[unc('A_'+c,method) for c in CTRS];x=np.arange(3)+(mi-1)*.05
                y=[vec(v)[j] for v in r];err=[2*sd(v)[j] for v in u]
                ax.errorbar(x,y,yerr=err,color=COL[method],marker='o',ms=4,lw=1,capsize=2,label=NAMES[method])
                for k,(v,w) in enumerate(zip(r,u)):
                    if not v['converged']:ax.plot(x[k],y[k],'x',color='red',ms=7,zorder=6)
                    elif not w.get('locally_valid_numerics'):ax.plot(x[k],y[k],'o',mfc='white',mec=COL[method],ms=4,zorder=6)
            biasv=[]
            for c in CTRS:
                b=json.loads((HERE/'analysis'/f'A_{c}_bias.json').read_text());biasv.append(vec(fit('A_'+c,'M3'))[j]*b['pred_MminusC_percent'][j]/100)
            pred=vec(fit('A_nominal','M1'))[j]+np.asarray(biasv)-biasv[1]
            ax.plot(np.arange(3),pred,'--',color=COL['M1'],lw=.8,label='Predicted M1 drift, nominal anchor')
            ax.set_xticks(range(3),['Low','Nominal','High']);ax.set_xlabel('Controller')
            ax.set_ylabel(f'f{j+1} (Hz)' if j<3 else f'Damping ζ{j-2} (%)');linegrid(ax)
        h,l=axs.flat[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=2)
        fig.suptitle('Independent A records; bars = ±2 conditional standard uncertainty\nOpen circles: strict numerical check not met',fontsize=10)
        save(fig,name,'fits/A_*; analysis/A_*_uncertainty.json; analysis/A_*_bias.json')

def frf():
    for ctrl in CTRS:
        a=np.load(HERE/'analysis'/f'{ctrl}_orthogonal_frf.npz');f=a['f'];g=a['G']
        anchor=np.unwrap(np.angle(np.load(HERE/'fits'/f'A_{ctrl}_M3.npz')['G']),axis=0)*180/np.pi
        phase=phase_branch(g,anchor)
        for kind in ('magnitude','phase','zoom_8_18'):
            fig,axs=plt.subplots(4,4,figsize=(10.5,9),sharex=True,layout='constrained')
            idx=(f>=8)&(f<=18) if kind=='zoom_8_18' else np.ones(len(f),bool)
            for o in range(4):
                for i in range(4):
                    ax=axs[o,i]
                    if kind=='phase':ref=phase[:,o,i];lo=ref+a['phase_ci_relative'][0,:,o,i];hi=ref+a['phase_ci_relative'][1,:,o,i]
                    else:ref=20*np.log10(np.maximum(abs(g[:,o,i]),1e-300));lo=a['mag_ci'][0,:,o,i];hi=a['mag_ci'][1,:,o,i]
                    ax.fill_between(f[idx],lo[idx],hi[idx],color='#aaaaaa',alpha=.28,lw=0);ax.plot(f[idx],ref[idx],color='#777777',lw=1.1)
                    for method in ('M1','M3','M3w','M5n'):
                        gm=np.load(HERE/'fits'/f'A_{ctrl}_{method}.npz')['G'][:,o,i]
                        if kind=='phase':v=phase_branch(gm,anchor[:,o,i])
                        else:v=20*np.log10(np.maximum(abs(gm),1e-300))
                        ax.plot(f[idx],v[idx],color=COL[method],lw=.8,ls={'M1':'-','M3':'-','M3w':'--','M5n':':'}[method])
                    ax.set_title(f'u{i+1} → y{o+1}');linegrid(ax)
                    if o==3:ax.set_xlabel('Frequency (Hz)')
                    if i==0:ax.set_ylabel('Phase (deg)' if kind=='phase' else 'Magnitude (dB)')
            handles=[Line2D([0],[0],color='#777777',label='C orthogonal FRF (pointwise 95% bootstrap)')]+[Line2D([0],[0],color=COL[x],ls={'M3w':'--','M5n':':'}.get(x,'-'),label=NAMES[x]) for x in ('M1','M3','M3w','M5n')]
            fig.legend(handles=handles,loc='outside lower center',ncol=3);fig.suptitle(f'{ctrl}: independent A models vs C orthogonal FRF\nAll bins retained; finite-ratio and between-run biases not covered by bands',fontsize=11)
            save(fig,f'C_{ctrl}_{kind}','analysis/*_orthogonal_frf.npz; fits/A_*_model NPZ')
    for kind in ('magnitude','phase'):
        fig,axs=plt.subplots(4,4,figsize=(10.5,9),sharex=True,layout='constrained')
        anchor=np.unwrap(np.angle(np.load(HERE/'analysis/lowgain_orthogonal_frf.npz')['G']),axis=0)*180/np.pi
        for ctrl in CTRS:
            a=np.load(HERE/'analysis'/f'{ctrl}_orthogonal_frf.npz');f=a['f'];g=a['G']
            val=phase_branch(g,anchor) if kind=='phase' else 20*np.log10(np.maximum(abs(g),1e-300))
            for o in range(4):
                for i in range(4):
                    axs[o,i].plot(f,val[:,o,i],lw=.9,color=CC[ctrl],label=ctrl);axs[o,i].set_title(f'u{i+1} → y{o+1}');linegrid(axs[o,i])
        for ax in axs[-1]:ax.set_xlabel('Frequency (Hz)')
        fig.legend(*axs[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=3)
        fig.suptitle('C orthogonal reference across controllers; no inverse regularization',fontsize=11)
        save(fig,'C_controller_comparison_'+kind,'analysis/*_orthogonal_frf.npz')

def scans_fig():
    rows=list(csv.DictReader((HERE/'analysis/scans.csv').open(encoding='utf-8-sig')))
    for kind in ('amplitude','repeat','position'):
        fig,axs=plt.subplots(2,3,figsize=(9,5.6),layout='constrained')
        for j,ax in enumerate(axs.flat):
            parameter=['f1','f2','f3','zeta1','zeta2','zeta3'][j]
            for method in (('M1','M2','M3','M3w') if kind=='position' else ('M2','M3','M3w')):
                rr=[r for r in rows if r['kind']==kind and r['parameter']==parameter and r['method']==method]
                x=np.array([float(r['x']) for r in rr]);y=np.array([float(r['z_score'] if kind=='repeat' else r['estimate']) for r in rr])
                err=None if kind=='repeat' else np.array([2*float(r['standard_uncertainty']) for r in rr])
                ax.errorbar(x,y,yerr=err,color=COL[method],marker='o',ms=3,lw=.8,capsize=2,label=NAMES[method])
                for k,r in enumerate(rr):
                    if r['solver_converged']=='False':ax.plot(x[k],y[k],'x',color='red',ms=7,zorder=6)
                    elif r['numerical_conditions_ok']=='False':ax.plot(x[k],y[k],'o',mfc='white',mec=COL[method],ms=3,zorder=6)
            if kind=='repeat':ax.axhline(2,c='#aaaaaa',ls='--',lw=.7);ax.axhline(-2,c='#aaaaaa',ls='--',lw=.7)
            ax.set_title(parameter);ax.set_xlabel({'amplitude':'Disturbance amplitude gain','repeat':'Independent run (A + four C)','position':'z position (mm)'}[kind])
            ax.set_ylabel('Centered estimate / standard uncertainty' if kind=='repeat' else 'Hz' if j<3 else 'Damping (%)');linegrid(ax)
        fig.legend(*axs[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=3)
        fig.suptitle({'amplitude':'D amplitude scan: linearity diagnostic, not IID repetitions','repeat':'Nominal same-amplitude repetitions: conditional uncertainty check','position':'E position scan (A nominal reused at z=0)'}[kind]+'\nOpen circles: numerical check not met; red crosses: solver failure',fontsize=10)
        save(fig,kind+'_modal_scan','analysis/scans.csv')

def residual():
    rows=json.loads((HERE/'analysis/frf_normalized_validation.json').read_text());fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    methods=('M1','M2','M3','M3w','M5n')
    for j,method in enumerate(methods):
        vals=[next(r['statistic_over_expectation'] for r in rows if r['controller']==c and r['method']==method) for c in CTRS]
        ax.bar(np.arange(3)+(j-2)*.15,vals,width=.14,color=COL[method],label=NAMES[method])
    ax.axhline(1,color='k',ls='--',lw=.8);ax.set_yscale('log');ax.set_xticks(range(3),['Low','Nominal','High']);ax.set_ylabel('Q / ideal fixed-model expectation');ax.set_title('Independent A model vs four C records; B-noise normalization')
    ax.legend(ncol=3,loc='upper center');linegrid(ax);save(fig,'C_normalized_residual','analysis/frf_normalized_validation.json')

if __name__=='__main__':
    controller_figs();core();frf();scans_fig();residual()
    js(HERE/'figures/manifest.json',manifest)
    print('FIGURES',len(manifest),flush=True)
