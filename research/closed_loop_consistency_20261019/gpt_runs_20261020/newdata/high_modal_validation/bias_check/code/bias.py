from common import *

def local_phi(m,a,b,fquery,width=.5):
    yh,_=m.output(m.context(a['f']),b,a['U'],False)
    ix=a['k']%int(a['P'])!=0;v=(a['Y']-yh)[ix];ff=a['f'][ix]
    outer=v[:,:,None]*v[:,None,:].conj()
    cs=np.concatenate((np.zeros((1,4,4),complex),np.cumsum(outer,axis=0)))
    left=np.searchsorted(ff,fquery-width/2-1e-12,side='left')
    right=np.searchsorted(ff,fquery+width/2+1e-12,side='right')
    count=right-left
    if min(count)<4:raise ValueError('Fewer than 4 nonexcited neighbors')
    Phi=(cs[right]-cs[left])/count[:,None,None]/int(a['N'])
    Phi=(Phi+Phi.conj().transpose(0,2,1))/2
    np.linalg.cholesky(Phi)
    return Phi,count

def moments(m,a,b,C,chunk=256,step=2e-4):
    """Low-frequency prediction equations, summed in bounded-memory chunks.

    C=N*Phi is covariance of unnormalised full-record DFT. G derivatives
    exclude deterministic boundary terms; Z derivatives include them.
    The feedback F and C are frozen in the curvature derivative.
    """
    nb=m.nb;M0=np.zeros((nb,nb));Mn=M0.copy();CC=M0.copy();bv=np.zeros(nb)
    A3=M0.copy();Q3=M0.copy();A3_observed=M0.copy()
    for start in range(0,len(a['f']),chunk):
        sl=slice(start,min(start+chunk,len(a['f'])))
        c=m.context(a['f'][sl]);G,Gj=m.gmodel(c,b)
        K=a['K'][sl];D=a['D'][sl];cov=C[sl]
        F=np.linalg.solve(np.eye(4)[None]+K@G,K)
        Ud=np.linalg.solve(np.eye(4)[None]+K@G,D[:,:,None])[:,:,0]
        _,J=m.output(c,b,Ud);flat=J.reshape(-1,nb);M0+=(flat.conj().T@flat).real
        L=np.linalg.cholesky(cov)
        nj=-np.einsum('koip,kil->kolp',Gj,F@L,optimize=True).reshape(-1,nb)
        Mn+=(nj.conj().T@nj).real
        bv-=np.einsum('koip,kil,klo->p',Gj,F,cov,optimize=True).real
        # Derive b's Jacobian with F and C fixed, exploiting real linearity.
        for j in range(m.iw,m.nb):
            bp=b.copy();bm=b.copy();bp[j]+=step;bm[j]-=step
            der=(m.gmodel(c,bp)[1]-m.gmodel(c,bm)[1])/(2*step)
            v=-np.einsum('koip,kil,klo->p',der,F,cov,optimize=True).real
            CC[:,j]+=v
        il=a['k'][sl]%int(a['P'])==0
        Z=J[il];fl=Z.reshape(-1,nb);A3+=(fl.conj().T@fl).real
        Q3+=.5*np.einsum('kop,koi,kiq->pq',Z.conj(),cov[il],Z,optimize=True).real
        cl={k:v[il] for k,v in c.items()};_,Jactual=m.output(cl,b,a['U'][sl][il])
        A3_observed+=np.einsum('kop,koq->pq',Z.conj(),Jactual,optimize=True).real
    # Linear-linear second derivatives are exactly zero. Complete mixed block.
    CC[m.iw:,:m.iw]=CC[:m.iw,m.iw:].T
    CC=(CC+CC.T)/2
    M=M0+Mn;d1=scaled_solve(M,bv);d2=scaled_solve(M-CC,bv)
    I3=scaled_solve(A3,np.eye(nb));P3=I3@Q3@I3.T;P3=(P3+P3.T)/2
    Iobs=scaled_solve(A3_observed,np.eye(nb));Pobs=Iobs@Q3@Iobs.T;Pobs=(Pobs+Pobs.T)/2
    return dict(M=M,M0=M0,Mn=Mn,Ccurvature=CC,b=bv,delta_M=d1,delta_MminusC=d2,
        pred_M=100*np.expm1(d1[m.iw:]),pred_MminusC=100*np.expm1(d2[m.iw:]),
        reference_covariance=P3,reference_observed_bread_covariance=Pobs,A3=A3,Q3=Q3,A3_observed=A3_observed,
        numerical_M=diagnostics(M),numerical_MminusC=diagnostics(M-CC),numerical_A3=diagnostics(A3))

def run_case(rid,w,omit=None,step=2e-4,save=True):
    tag=f'{rid}_{w}_'+('full' if omit is None else f'omit{omit}')
    dest=HERE/'raw'/f'{tag}_prediction.json'
    if save and dest.exists():return rd(dest)
    t=time.perf_counter();r=dict(id=rid,window=w,omit=omit)
    try:
        s,m,_,_=setup(w);a=spectrum(rid,w,omit);ext=spectrum(rid,w,omit,True)
        fits={method:(basefit(rid,w,method) if omit is None else foldfit(rid,w,a,omit,method)) for method in ('M1','M3')}
        b=np.asarray(fits['M3']['beta']);b1=np.asarray(fits['M1']['beta'])
        Phi,count=local_phi(m,ext,b,a['f'])
        pr=moments(m,a,b,int(a['N'])*Phi,step=step)
        obs=100*np.expm1(b1[m.iw:]-b[m.iw:])
        f,z=m.modal(b);diag=np.diag(pr['reference_covariance'])[m.iw:]
        sd=np.sqrt(np.where(diag>=0,diag,np.nan))
        r.update(success=True,P=a['P'],N=a['N'],nb=m.nb,report_components=s['report'],
             base_eligible=all(eligible(fits[k]) for k in fits),
             fits={k:{q:v for q,v in x.items() if q not in ('beta','solver','initial_beta','polish_history')} for k,x in fits.items()},
             observed=obs,pred_M=pr['pred_M'],pred_MminusC=pr['pred_MminusC'],
             difference_M=obs-pr['pred_M'],difference_MminusC=obs-pr['pred_MminusC'],
             frequency_hz=f,damping_percent=z,reference_log_sd=sd,frequency_u_hz=f*sd[:m.nm],
             damping_u_pp=z*sd[m.nm:],minimum_neighbors=count.min(),maximum_neighbors=count.max(),
             covariance_condition_median=np.median(np.linalg.cond(Phi)),
             numerical_M=pr['numerical_M'],numerical_MminusC=pr['numerical_MminusC'],
             numerical_A3=pr['numerical_A3'])
        if save:np.savez_compressed(HERE/'raw'/f'{tag}_prediction.npz',f=a['f'],Phi=Phi,counts=count,
              **{k:v for k,v in pr.items() if isinstance(v,np.ndarray)})
    except Exception as e:r.update(success=False,exception=repr(e),traceback=traceback.format_exc())
    r['seconds']=time.perf_counter()-t
    if save:js(dest,r)
    print('PRED',tag,r['success'],round(r['seconds'],2),flush=True)
    return r

if __name__=='__main__':
    import argparse
    from concurrent.futures import ProcessPoolExecutor,as_completed
    p=argparse.ArgumentParser();p.add_argument('--id');p.add_argument('--window',default='H493')
    p.add_argument('--folds',action='store_true');p.add_argument('--all',action='store_true');p.add_argument('--workers',type=int,default=2)
    args=p.parse_args()
    if args.all:
        jobs=[]
        for w in SPECS:
            for rid in IDS:
                if args.folds:
                    P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])
                    jobs.extend((rid,w,i) for i in range(P))
                else:jobs.append((rid,w,None))
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for f in as_completed([pool.submit(run_case,*job) for job in jobs]):f.result()
    else:run_case(args.id or 'A_nominal',args.window,0 if args.folds else None)
