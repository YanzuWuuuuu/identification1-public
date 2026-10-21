from common import *

def run_one(rid,method,seed):
    a=oldspec(rid);li=a['line_indices'];f=a['f'][li];N=int(a['N'])
    b=np.array(fit(rid,method)['beta']);c=m.Ctx(f,TAU,'zoh')
    U,Y,D,K=(a[q][li] for q in ('U','Y','D','K'))
    Uh=np.linalg.solve(np.eye(4)[None]+K@m.resp(c,b),D[:,:,None])[:,:,0]
    yh,J=m.out_jac(c,b,U);Z=m.out_jac(c,b,Uh)[1];e=Y-yh
    ns=np.load(HERE/'spectra'/f'{rid}_noise.npz')
    C=N*ns['Phi_0.5'][li]
    W=np.broadcast_to(np.eye(4,dtype=complex),(len(f),4,4))
    if method=='M3w':
        cv=np.load(SRC/'fits'/f'{rid}_local_covariance.npz')
        W=wc.precision(cv['covariance_residual'])
    Cr=.5*realmat(C);Wr=realmat(W);L=np.linalg.cholesky(Cr)
    Jr,Zr,er=realvec(J),realvec(Z),realvec(e)
    A=np.einsum('kop,koi,kiq->pq',Zr,Wr,Jr,optimize=True)
    B=np.linalg.solve(L,Jr).reshape(-1,m.NB)
    ZWL=np.einsum('kop,koi,kij->kpj',Zr,Wr,L,optimize=True).transpose(1,0,2).reshape(m.NB,-1)
    H=scaled_solve(A,ZWL);Rw=np.eye(len(B))-B@H
    ew=np.linalg.solve(L,er[:,:,None])[:,:,0]
    qf=np.sum(ew**2,axis=1);ef=np.sum(Rw**2,axis=1).reshape(len(f),8).sum(1)
    Q=qf.sum();EQ=ef.sum()
    # QR-stabilized low-rank trace identity: avoid squaring ill-scaled Gram matrices.
    bscale=np.linalg.norm(B,axis=0)
    bq,br=np.linalg.qr(B/bscale,mode='reduced')
    rh=br@(H*bscale[:,None])
    EQ_lowrank=len(B)-2*np.einsum('ij,ji',bq,rh)+np.sum(rh**2)
    rng=np.random.default_rng(seed);v=rng.standard_normal((len(B),200));estar=Rw@v
    qb=np.sum(estar**2,axis=0)
    sf=[]
    for lo,hi,mask in bands(f):
        sf.append(dict(lo=lo,hi=hi,Q=qf[mask].sum(),EQ=ef[mask].sum(),ratio=qf[mask].sum()/ef[mask].sum()))
    score=wc.iv_score(lambda bb,u:m.out_jac(c,bb,u),U,Y,lambda bb:np.linalg.solve(np.eye(4)[None]+K@m.resp(c,bb),D[:,:,None])[:,:,0],W,b)
    row=dict(id=rid,method=method,lines=len(f),real_dimension=len(B),n_parameters=m.NB,
             Q=Q,EQ=EQ,ratio=Q/EQ,nominal_n_minus_p=len(B)-m.NB,numerics=diagnostics(A),
             EQ_lowrank=EQ_lowrank,trace_identity_relative_error=abs(EQ-EQ_lowrank)/EQ,
             projector_relative_error=np.linalg.norm(Rw@B)/np.linalg.norm(B),
             bootstrap_n=200,bootstrap_seed=seed,bootstrap_mean=qb.mean(),bootstrap_sd=qb.std(ddof=1),
             bootstrap_mean_se=qb.std(ddof=1)/np.sqrt(200),bootstrap_q95=np.percentile(qb,95),
             bootstrap_interval95=np.percentile(qb,[2.5,97.5]),
             bootstrap_mean_standardized=(qb.mean()-EQ)/(qb.std(ddof=1)/np.sqrt(200)),
             frequency_bands=sf,score_recomputed=score,original_score=fit(rid,method)['score_normalized'],
             covariance='same-record M3 endpoint, Ws=0.5 Hz; full raw 4x4 spectral matrices; proper complex convention',
             weights='identity' if method=='M3' else 'original frozen 16-neighbor covariance inverse (alpha=0.1)',
             interpretation='conditional linearized IV projection; not exact chi-square when PSD estimated, instruments depend on fit, or model misspecified')
    np.savez_compressed(HERE/'raw'/f'{rid}_{method}_residual.npz',
                        f=f,e=e,C=C,W=W,J=J,Z=Z,q_frequency=qf,expected_frequency=ef,
                        bootstrap_Q=qb,A=A,beta=b)
    js(HERE/'analysis'/f'{rid}_{method}_residual.json',row)
    print('RESIDUAL',rid,method,'Q/E',round(Q/EQ,3),'E',round(EQ,1),'bootz',round(row['bootstrap_mean_standardized'],2),flush=True)
    return row

def lowgain_check():
    ctrl='lowgain'
    a=np.load(SRC/'analysis'/f'{ctrl}_orthogonal_frf.npz');f=a['f']
    nb=np.load(SRC/'analysis'/f'{ctrl}_noise.npz');C=P0/int(a['P_common'])*nb['Phi']
    b=np.array(fit('A_lowgain','M3')['beta']);G=m.resp(m.Ctx(f,TAU,'zoh'),b)
    e=a['Y']-G@a['U']
    q=2*np.einsum('kor,koi,kir->k',e.conj(),np.linalg.inv(C),e,optimize=True).real
    tab=[]
    for lo,hi,mask in bands(f):
        tab.append(dict(lo=lo,hi=hi,Q=q[mask].sum(),fraction=q[mask].sum()/q.sum(),ratio=q[mask].sum()/(32*np.sum(mask))))
    sel=f<6
    row=dict(id='lowgain_A_to_four_C',Q=q.sum(),dof=32*len(f),ratio=q.sum()/(32*len(f)),fraction_3_6=q[sel].sum()/q.sum(),bands=tab,
             P_common=int(a['P_common']),note='original independent validation: A model, four C records, B covariance, no transient correction to periodic means')
    arrays=dict(f=f,q=q)
    for rid in ('A_lowgain','C_lowgain_pppp'):
        sp=oldspec(rid);N=int(sp['N'])
        # Single-record one-sided power spectral density, includes excited lines.
        psd=2*abs(sp['Y'])**2/(FS*N)
        arrays[rid+'_f']=sp['f'];arrays[rid+'_psd']=psd
        df=FS/N;low=sp['f']<6
        row[rid]=dict(total_band_rms=np.sqrt(psd.sum(0)*df),low_3_6_rms=np.sqrt(psd[low].sum(0)*df),
                      low_3_6_power_fraction=psd[low].sum(0)/psd.sum(0))
    np.savez_compressed(HERE/'raw/lowgain_crossrecord.npz',**arrays)
    js(HERE/'analysis/lowgain_crossrecord.json',row)
    print('LOWGAIN_CROSS',row['ratio'],'fraction3to6',row['fraction_3_6'],flush=True)
if __name__=='__main__':
    rows=[]
    for j,rid in enumerate(IDS):
        for k,method in enumerate(('M3','M3w')):
            rows.append(run_one(rid,method,102020+2*j+k))
    js(HERE/'analysis/residual_checks.json',rows)
    lowgain_check()
