from common import *
import argparse

def local_phi(a,b,fquery,width,refit_transient=False):
    noise=a['k']%int(a['P'])!=0
    c=m.Ctx(a['f'],TAU,'zoh');bb=b.copy()
    if refit_transient:
        yy,J=m.out_jac(c,bb,a['U'])
        train=noise&(a['f']>=3)&(a['f']<=25)
        X=realvec(J[train,:,m.ITR:]).reshape(-1,m.NB-m.ITR)
        y=realvec(a['Y'][train]-yy[train]).ravel()
        sc=1/np.maximum(np.linalg.norm(X,axis=0),1e-300)
        bb[m.ITR:]+=sc*np.linalg.lstsq(X*sc,y,rcond=1e-12)[0]
    yy,_=m.out_jac(c,bb,a['U'],need_jac=False)
    v=(a['Y']-yy)[noise];fn=a['f'][noise]
    C=[];counts=[]
    for fq in fquery:
        use=np.abs(fn-fq)<=width/2+1e-12
        vv=v[use];counts.append(len(vv))
        if len(vv)<4:raise ValueError('Less than four nonexcited points')
        C.append(vv.T@vv.conj()/len(vv)/int(a['N']))
    C=np.asarray(C)
    np.linalg.cholesky(C) # Do not add a loading or change smoothing if this fails.
    return C,np.array(counts),v,bb

def cache_moments(a,b):
    c=m.Ctx(a['f'],TAU,'zoh')
    G,Gj=g_jac(c,b);K=a['K']
    F=np.linalg.solve(np.eye(4)[None]+K@G,K)
    Ud=np.linalg.solve(np.eye(4)[None]+K@G,a['D'][:,:,None])[:,:,0]
    _,J=m.out_jac(c,b,Ud);flat=J.reshape(-1,m.NB);M0=(flat.conj().T@flat).real
    # Cache only six mixed transfer-function derivatives. Outer step is fixed.
    derivatives=[]
    for j in MOD:
        bp=b.copy();bm=b.copy();bp[j]+=2e-4;bm[j]-=2e-4
        derivatives.append((g_jac(c,bp)[1]-g_jac(c,bm)[1])/(4e-4))
    return dict(c=c,G=G,Gj=Gj,F=F,Ud=Ud,J=J,M0=M0,derivatives=derivatives)

def predict(cache,C,signal_scale=1.):
    Gj,F=cache['Gj'],cache['F'];L=np.linalg.cholesky(C)
    nj=-np.einsum('koip,kil->kolp',Gj,F@L,optimize=True).reshape(-1,m.NB)
    Mn=(nj.conj().T@nj).real
    bv=-np.einsum('koip,kil,klo->p',Gj,F,C,optimize=True).real
    CC=np.zeros((m.NB,m.NB))
    for j,der in zip(MOD,cache['derivatives']):
        v=-np.einsum('koip,kil,klo->p',der,F,C,optimize=True).real
        CC[:,j]=v;CC[j,:]=v
    CC=(CC+CC.T)/2
    # Scaling deterministic input and boundary transient amplitudes by alpha:
    # all physical derivative columns scale alpha; transient-coordinate columns do not.
    scaling=np.ones(m.NB);scaling[:m.ITR]=signal_scale
    M0=cache['M0']*scaling[:,None]*scaling[None,:]
    M=M0+Mn;d1=scaled_solve(M,bv);d2=scaled_solve(M-CC,bv)
    return dict(M=M,M0=M0,Ccurvature=CC,b=bv,delta_M=d1,delta_MminusC=d2,
                pred_M_percent=100*np.expm1(d1[MOD]),pred_MminusC_percent=100*np.expm1(d2[MOD]),
                numerical_M=diagnostics(M),numerical_MminusC=diagnostics(M-CC))

def paired_sandwich(cache,C,li,pred,observed):
    J,Gj,F=cache['J'],cache['Gj'],cache['F']
    A1=pred['M']-pred['Ccurvature']
    JL=J[li];A3=np.einsum('kop,koq->pq',JL.conj(),JL,optimize=True).real
    Q1=.5*np.einsum('kop,koi,kiq->pq',J.conj(),C,J,optimize=True).real
    L=np.linalg.cholesky(C);H=-np.einsum('koip,kij->kojp',Gj,F,optimize=True)
    H=(H+H.conj().transpose(0,2,1,3))/2
    HL=np.einsum('kai,kabp,kbj->kijp',L.conj(),H,L,optimize=True).reshape(-1,m.NB)
    Q1+=(HL.conj().T@HL).real
    Q3=.5*np.einsum('kop,koi,kiq->pq',JL.conj(),C[li],JL,optimize=True).real
    I1=scaled_solve(A1,np.eye(m.NB));I3=scaled_solve(A3,np.eye(m.NB))
    P1=I1@Q1@I1.T;P3=I3@Q3@I3.T
    # Shared linear noise at the lines; centered Gaussian quadratic moment has zero cross-covariance with linear moment.
    P13=I1@Q3@I3.T
    Pd=P1+P3-P13-P13.T;Pd=(Pd+Pd.T)/2
    se=(100+observed)*np.sqrt(np.maximum(np.diag(Pd)[MOD],0))
    return dict(paired_se=se,paired_log_variance=np.diag(Pd)[MOD],
                M1=diagnostics(A1),M3=diagnostics(A3),negative_diagonal=bool(np.any(np.diag(Pd)[MOD]<0)),
                interpretation='conditional local Gaussian linear+quadratic approximation at M3; includes M-C curvature and shared noise; PSD/model/K estimation and non-Gaussian errors excluded'),Pd

def ratios(C,B,f):
    d=np.diagonal(C,axis1=1,axis2=2).real/np.diagonal(B,axis1=1,axis2=2).real
    tr=np.trace(C,axis1=1,axis2=2).real/np.trace(B,axis1=1,axis2=2).real
    out=dict(diagonal=d,trace=tr)
    for name,mask in [('band',np.ones(len(f),bool)),('low_3_6',(f>=3)&(f<6)),('modal_9_18',(f>=9)&(f<18)),('upper_18_25',f>=18)]:
        out[name]=np.trace(C[mask],axis1=1,axis2=2).real.sum()/np.trace(B[mask],axis1=1,axis2=2).real.sum()
    return out

def run(rid):
    t=time.perf_counter();a=oldspec(rid);b=np.array(fit(rid,'M3')['beta'])
    observed=100*np.expm1(np.array(fit(rid,'M1')['beta'])[MOD]-b[MOD])
    f=a['f'];N=int(a['N']);extended=spectrum(rid,extended=True)
    # Check copied data entry and controller agree exactly with original spectra.
    check=spectrum(rid)
    selfcheck={key:float(np.max(abs(check[key]-a[key]))) for key in ('f','U','Y','D','K')}
    if max(selfcheck.values())>1e-10:raise ValueError(selfcheck)
    no=np.load(SRC/'analysis'/f"{meta(rid)['controller']}_noise.npz");B=interp_cov(f,no['f'],no['Phi'])
    cache=cache_moments(a,b)
    old=predict(cache,N*B);oldj=rd(SRC/'analysis'/f'{rid}_bias.json')
    reproduction=float(np.max(abs(old['pred_MminusC_percent']-np.array(oldj['pred_MminusC_percent']))))
    if reproduction>1e-5:raise ValueError(('old bias mismatch',reproduction))
    rec=dict(id=rid,P=int(a['P']),N=N,selfcheck=selfcheck,old_prediction_reproduction_max_abs=reproduction,
             observed=observed,old_B={key:old[key] for key in ('pred_M_percent','pred_MminusC_percent')},windows={},
             extra_drop=[],primary_covariance='raw centered-frequency window second moments; no shrinkage/no diagonal loading')
    arrays=dict(f=f,B=B)
    for width in (.5,1.):
        C,count,v,_=local_phi(extended,b,f,width)
        pred=predict(cache,N*C);rat=ratios(C,B,f)
        key=str(width)
        rec['windows'][key]={q:pred[q] for q in ('pred_M_percent','pred_MminusC_percent','numerical_M','numerical_MminusC')}
        rec['windows'][key].update(count_min=count.min(),count_max=count.max(),covariance_condition_median=np.median(np.linalg.cond(C)),
                   ratios={k:v for k,v in rat.items() if np.ndim(v)==0},
                   observed_over_prediction=observed/pred['pred_MminusC_percent'],
                   prediction_over_observed=pred['pred_MminusC_percent']/observed)
        arrays.update({f'Phi_{key}':C,f'count_{key}':count,f'ratio_diag_{key}':rat['diagonal'],f'ratio_trace_{key}':rat['trace']})
        np.savez_compressed(HERE/'raw'/f'{rid}_prediction_{key}.npz',**{q:v for q,v in pred.items() if isinstance(v,np.ndarray)})
        if width==.5:
            ss,Pd=paired_sandwich(cache,N*C,a['line_indices'],pred,observed)
            rec['paired_sandwich']=ss;np.savez_compressed(HERE/'raw'/f'{rid}_paired_covariance.npz',covariance=Pd)
            if rid.startswith('A_'):
                rh=rms_record('A_highgain');rr=rms_record(rid)
                alpha=float(min(1,np.min(rh/rr)))
                sc=predict(cache,N*C,signal_scale=alpha)
                rec['scaled_excitation']=dict(alpha=alpha,disturb_gain=1.5*alpha,
                    per_output_gain=1.5*rh/rr,aggregate_gain=1.5*np.linalg.norm(rh)/np.linalg.norm(rr),
                    predicted_response_rms=rr*alpha,original_response_rms=rr,target_response_rms=rh,
                    pred_M_percent=sc['pred_M_percent'],pred_MminusC_percent=sc['pred_MminusC_percent'],
                    convention='common alpha is minimum channel target/current RMS, clipped at one; deterministic transients scale alpha, noise held fixed')
    for drop in (0,1,2):
        cut=spectrum(rid,drop=drop,extended=True)
        row=dict(extra_periods_dropped=drop,P=cut['P'])
        for refit in (False,True):
            C,count,v,bc=local_phi(cut,b,f,.5,refit_transient=refit)
            rat=ratios(C,B,f);key=f'drop{drop}_'+('reboundary' if refit else 'frozen')
            arrays['Phi_'+key]=C;arrays['ratio_trace_'+key]=rat['trace']
            arrays['ratio_diag_'+key]=rat['diagonal']
            row['reboundary' if refit else 'frozen']={k:v for k,v in rat.items() if np.ndim(v)==0}
            if refit:row['refitted_transient_beta']=bc[m.ITR:]
        rec['extra_drop'].append(row)
    np.savez_compressed(HERE/'spectra'/f'{rid}_noise.npz',**arrays)
    rec['seconds']=time.perf_counter()-t;js(HERE/'analysis'/f'{rid}_noise_bias.json',rec)
    print('BIAS',rid,'obs',np.round(observed[3:],3),'new',np.round(rec['windows']['0.5']['pred_MminusC_percent'][3:],3),'sec',round(rec['seconds'],1),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--id');args=p.parse_args()
    for rid in [args.id] if args.id else IDS:
        try:run(rid)
        except Exception as e:
            js(HERE/'analysis'/f'{rid}_noise_bias_ERROR.json',dict(exception=repr(e),traceback=traceback.format_exc()))
            raise
