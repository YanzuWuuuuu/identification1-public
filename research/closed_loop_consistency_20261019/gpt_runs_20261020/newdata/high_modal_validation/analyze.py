from run import *
import csv,itertools
from scipy.optimize import linear_sum_assignment

def loadfit(rid,w,method):return json.loads((HERE/'raw'/f'{rid}_{w}_{method}.json').read_text(encoding='utf-8'))
def csvout(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/'analysis'/name).open('w',encoding='utf-8-sig',newline='') as h:
        writer=csv.DictWriter(h,fieldnames=keys);writer.writeheader();writer.writerows(rows)
def residual_cov(V):
    out=[];n=V.shape[1]
    for j in range(n):
        start=min(max(j-2,0),n-5);r=V[:,start:start+5].reshape(-1,4)
        out.append(r.T@r.conj()/len(r)/P0)
    return wc.shrink_cov(np.asarray(out),.1)

def noise_spectra():
    rows=[]
    for ctrl in ('lowgain','nominal','highgain'):
        a=np.load(PARENT/'prepared'/f'B_{ctrl}.npz');P=int(a['P'])
        Uf=np.fft.rfft(a['U'].reshape(P,P0,4),axis=1);Yf=np.fft.rfft(a['Y'].reshape(P,P0,4),axis=1)
        for w,s in SPECS.items():
            _,m,_,_=setup(w);kp=np.arange(round(s['band'][0]*P0/FS),round(s['band'][1]*P0/FS)+1);f=kp*FS/P0
            r=loadfit('A_'+ctrl,w,'M3');G=m.resp(m.context(f),np.asarray(r['beta']))
            V=Yf[:,kp]-np.einsum('koi,pki->pko',G,Uf[:,kp]);Phi=residual_cov(V)
            np.savez_compressed(HERE/'analysis'/f'noise_{ctrl}_{w}.npz',f=f,Phi=Phi)
            rows.append(dict(controller=ctrl,window=w,reference_solver_converged=r['converged'],
                reference_strict_score=r['strict_score_pass'],reference_near_bound=r.get('near_bound',True),
                noise_periods=P,noise_trace_sum=float(np.trace(Phi,axis1=1,axis2=2).real.sum()),
                median_noise_covariance_condition=float(np.median(np.linalg.cond(Phi)))))
    csvout('noise_summary.csv',rows)

def uncertainty(rid,w,method,fit_override=None,output_dir=None):
    s,m,_,_=setup(w);r=loadfit(rid,w,method) if fit_override is None else fit_override;a=get_data(rid,w);li=a['line_indices'];f=a['f'][li];c=m.context(f);b=np.asarray(r['beta'])
    G=m.resp(c,b);Ud=np.linalg.solve(np.eye(4)[None]+a['K'][li]@G,a['D'][li,:,None])[:,:,0]
    Z=m.output(c,b,Ud)[1]
    W=np.broadcast_to(np.eye(4),(len(f),4,4))
    if method=='M3w':W=wc.precision(np.load(HERE/'raw'/f'{rid}_{w}_local_cov.npz')['covariance_residual'])
    wz=np.einsum('koi,kip->kop',W,Z);AA=np.einsum('kop,koq->pq',Z.conj(),wz,optimize=True).real
    Phi=np.load(HERE/'analysis'/f"noise_{r['controller']}_{w}.npz")['Phi'];C=int(a['N'])*Phi
    QQ=.5*np.einsum('kop,koi,kiq->pq',wz.conj(),C,wz,optimize=True).real
    sc=1/np.sqrt(np.maximum(abs(np.diag(AA)),1e-300));As=AA*sc[:,None]*sc[None,:]
    sv=np.linalg.svd(As,compute_uv=False);rank=int(np.count_nonzero(sv>sv[0]*1e-12))
    inv=sc[:,None]*np.linalg.solve(As,np.diag(sc));cov=inv@QQ@inv.T;cov=(cov+cov.T)/2
    var=np.diag(cov)[m.iw:];ss=np.sqrt(np.where(var>=0,var,np.nan));ff,zz=m.modal(b)
    valid=r['converged'] and r['strict_score_pass'] and not r.get('near_bound',True) and rank==m.nb and np.isfinite(ss).all()
    out=dict(id=rid,window=w,method=method,frequency_u_hz=(ff*ss[:m.nm]).tolist(),damping_u_pp=(zz*ss[m.nm:]).tolist(),
        local_numerical_conditions=bool(valid),scaled_rank=rank,dimension=m.nb,scaled_condition=float(sv[0]/sv[-1]),
        log_modal_sd=ss.tolist(),scope='B noise at z=0; full transient nuisance block; local correct-model sandwich; not model/operating-condition uncertainty',
        estimated_weight_frozen=method=='M3w')
    dest=HERE/'analysis' if output_dir is None else output_dir
    js(dest/f'{rid}_{w}_{method}_u.json',out)
    np.savez_compressed(dest/f'{rid}_{w}_{method}_cov.npz',cov=cov)
    return out

def residues(m,b):
    lin=b[:m.iw].reshape(4,m.L);f,z=m.modal(b);out=[]
    for j in range(m.nm):
        ww=2*np.pi*f[j];ze=z[j]/100;pole=ww*(-ze+1j*np.sqrt(1-ze*ze))
        out.append((ww*ww*lin[:,10*j:10*j+4]+ww*pole*lin[:,10*j+4:10*j+8])/(2*pole+2*ze*ww))
    return np.asarray(out)

def match(m,base,b,s):
    f0,_=m.modal(base);f,_=m.modal(b);R0=residues(m,base).reshape(m.nm,-1);R=residues(m,b).reshape(m.nm,-1)
    mac=abs(R0.conj()@R.T)**2/np.maximum(np.sum(abs(R0)**2,axis=1)[:,None]*np.sum(abs(R)**2,axis=1)[None,:],1e-300)
    freq=abs(np.log(f[None]/f0[:,None]));cost=1-mac+freq/np.log(1.1)
    bg=np.asarray(s['background']);cost[bg[:,None]!=bg[None,:]]=1e6
    ii,jj=linear_sum_assignment(cost);perm=np.empty(m.nm,int);perm[ii]=jj
    return [dict(anchor=j,fit_index=int(perm[j]),MAC=float(mac[j,perm[j]]),relative_frequency_difference=float(abs(f[perm[j]]/f0[j]-1)),
        candidate_match=bool(mac[j,perm[j]]>=.7 and abs(f[perm[j]]/f0[j]-1)<.10),background=bool(bg[j])) for j in range(m.nm)]

def collect():
    rows=[];fits=[];allmatches={}
    for w,s in SPECS.items():
        _,m,_,_=setup(w);anchor=np.asarray(loadfit('A_nominal',w,'M3')['beta'])
        for rid in IDS:
            for method in METHODS:
                r=loadfit(rid,w,method);r['matched']=match(m,anchor,np.asarray(r['beta']),s);fits.append(r)
                up=HERE/'analysis'/f'{rid}_{w}_{method}_u.json';u=json.loads(up.read_text()) if up.exists() else {}
                for mt in r['matched']:
                    j=mt['fit_index'];row=dict(id=rid,controller=r['controller'],record_group=rid[0],window=w,method=method,anchor_component=mt['anchor']+1,
                        estimated_component=j+1,frequency_hz=r['frequency_hz'][j],damping_percent=r['damping_percent'][j],
                        frequency_u_hz=u.get('frequency_u_hz',[None]*m.nm)[j],damping_u_pp=u.get('damping_u_pp',[None]*m.nm)[j],
                        candidate_match=mt['candidate_match'],residue_MAC=mt['MAC'],background=mt['background'],converged=r['converged'],
                        strict_score_pass=r['strict_score_pass'],near_bound=r.get('near_bound',True),iterations=r.get('iterations'),score=r.get('score'),
                        line_residual_percent=r.get('line_residual_percent'),numerical_u_valid=u.get('local_numerical_conditions',False))
                    rows.append(row)
    csvout('modal_parameters.csv',rows);js(HERE/'analysis/matched_fits.json',fits)
    lookup={(r['id'],r['window'],r['method'],r['anchor_component']):r for r in rows}
    pairs=[];repeats=[]
    for w,s in SPECS.items():
        for method in METHODS:
            for comp in range(1,len(s['f'])+1):
                for ca,cb in itertools.combinations(('lowgain','nominal','highgain'),2):
                    aa=lookup[('A_'+ca,w,method,comp)];bb=lookup[('A_'+cb,w,method,comp)]
                    for par,un in (('frequency_hz','frequency_u_hz'),('damping_percent','damping_u_pp')):
                        delta=bb[par]-aa[par];cu=np.hypot(aa[un],bb[un]) if aa[un] is not None and bb[un] is not None else np.nan
                        ok=aa['numerical_u_valid'] and bb['numerical_u_valid'] and aa['candidate_match'] and bb['candidate_match']
                        pairs.append(dict(window=w,method=method,component=comp,background=aa['background'],parameter=par,a=ca,b=cb,
                            difference=delta,relative_difference_percent=100*np.expm1(np.log(bb[par]/aa[par])),combined_u=float(cu),
                            abs_z=float(abs(delta)/cu) if np.isfinite(cu) and cu>0 else None,
                            below_2u=bool(abs(delta)<2*cu) if np.isfinite(cu) else None,eligible_for_check=bool(ok)))
                for ctrl in ('lowgain','nominal','highgain'):
                    aa=lookup[('A_'+ctrl,w,method,comp)];bb=lookup[('C_'+ctrl+'_pppp',w,method,comp)]
                    repeats.append(dict(window=w,controller=ctrl,method=method,component=comp,background=aa['background'],
                        f_difference_hz=bb['frequency_hz']-aa['frequency_hz'],zeta_difference_pp=bb['damping_percent']-aa['damping_percent'],
                        f_relative_difference_percent=100*np.expm1(np.log(bb['frequency_hz']/aa['frequency_hz'])),
                        zeta_relative_difference_percent=100*np.expm1(np.log(bb['damping_percent']/aa['damping_percent'])),
                        both_solver_converged=aa['converged'] and bb['converged'],both_candidate_match=aa['candidate_match'] and bb['candidate_match']))
    csvout('controller_pairs.csv',pairs);csvout('independent_repeats.csv',repeats)
    summary=[]
    for w,s in SPECS.items():
        for method in METHODS:
            rr=[r for r in fits if r['window']==w and r['method']==method];pp=[r for r in pairs if r['window']==w and r['method']==method and not r['background']]
            summary.append(dict(window=w,method=method,total_fits=len(rr),solver_failures=sum(not r['converged'] for r in rr),
                near_bounds=sum(r.get('near_bound',True) for r in rr),strict_score_failures=sum(not r['strict_score_pass'] for r in rr),
                median_seconds=float(np.median([r['seconds'] for r in rr])),
                nonbackground_parameter_pairs=len(pp),eligible_pairs=sum(r['eligible_for_check'] for r in pp),
                eligible_below_2u=sum(r['eligible_for_check'] and r['below_2u'] for r in pp),
                raw_below_2u=sum(bool(r['below_2u']) for r in pp)))
    js(HERE/'analysis/summary.json',summary);csvout('summary.csv',summary)
    return rows

def validation():
    out=[]
    for ctrl in ('lowgain','nominal','highgain'):
        a=np.load(PARENT/'wideband_scope_audit'/f'{ctrl}_broadband.npz')
        for w,s in SPECS.items():
            _,m,_,_=setup(w);mask=(a['f']>=s['band'][0])&(a['f']<=s['band'][1]);f=a['f'][mask];c=m.context(f)
            for method in METHODS:
                r=loadfit('A_'+ctrl,w,method);G=m.resp(c,np.asarray(r['beta']));Gref=a['G'][mask];err=a['Y'][mask]-G@a['U'][mask]
                b=get_data('C_'+ctrl+'_pppp',w);li=b['line_indices'];YY=b['Y'][li];UU=b['U'][li]
                # Ignore A-specific transients for transfer prediction on C.
                Cpred=G@UU[:,:,None];Cerr=YY-Cpred[:,:,0]
                out.append(dict(window=w,controller=ctrl,method=method,model_source='A',reference_source='four independent C runs',
                    complex_FRF_relative_difference_percent=float(100*np.linalg.norm(G-Gref)/np.linalg.norm(Gref)),
                    inverse_free_C_four_output_residual_percent=float(100*np.linalg.norm(err)/np.linalg.norm(a['Y'][mask])),
                    independent_C_pppp_raw_output_residual_percent=float(100*np.linalg.norm(Cerr)/np.linalg.norm(YY)),
                    includes_test_transients='no fitted test transient; may include small finite-record transient residual',
                    converged=r['converged'],near_bound=r.get('near_bound',True)))
    csvout('FRF_validation.csv',out)

def main():
    noise_spectra()
    for w in SPECS:
        for rid in IDS:
            for method in ('M3','M3w'):
                try:uncertainty(rid,w,method)
                except Exception as e:js(HERE/'analysis'/f'{rid}_{w}_{method}_u.json',dict(error=repr(e),local_numerical_conditions=False))
    collect();validation();print('ANALYSIS_COMPLETE',flush=True)

if __name__=='__main__':main()
