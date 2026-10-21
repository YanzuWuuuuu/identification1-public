from identify import *
from scipy.stats import chi2

CTRS=('lowgain','nominal','highgain')
SIGNS=('pppp','pmpm','ppmm','pmmp')
def right(a,b):return np.linalg.solve(b.swapaxes(-1,-2),a.swapaxes(-1,-2)).swapaxes(-1,-2)

def estimate(ctrl,B=1000):
    rng=np.random.default_rng(20261020+CTRS.index(ctrl));blocks=[]
    for sign in SIGNS:
        rid=f'C_{ctrl}_{sign}';a=np.load(HERE/'prepared'/f'{rid}.npz');p=int(a['P'])
        blocks.append({x:np.fft.rfft(a[x].reshape(p,P0,4),axis=1)[:,6:3001] for x in ('U','Y','D')})
    pp=min(x['U'].shape[0] for x in blocks);f=np.arange(6,3001)*FS/P0
    mu={x:np.stack([a[x][:pp].mean(0) for a in blocks],axis=2) for x in ('U','Y','D')}
    TU=right(mu['U'],mu['D']);TY=right(mu['Y'],mu['D']);G=right(mu['Y'],mu['U'])
    co=np.load(HERE/'controllers'/f'{ctrl}.npz');K=evalK(f,co['A'],co['B'],co['T'])
    Su=np.linalg.inv(np.eye(4)[None]+K@G);Sy=np.linalg.inv(np.eye(4)[None]+G@K);F=Su@K
    np.savez_compressed(HERE/'analysis'/f'{ctrl}_wide_frf.npz',f=f,G=G,TU=TU,TY=TY,D=mu['D'],U=mu['U'],Y=mu['Y'],K=K,Su=Su,Sy=Sy,F=F,
                        singular_K=np.linalg.svd(K,compute_uv=False),singular_loop=np.linalg.svd(K@G,compute_uv=False),
                        singular_Su=np.linalg.svd(Su,compute_uv=False),singular_Sy=np.linalg.svd(Sy,compute_uv=False),singular_F=np.linalg.svd(F,compute_uv=False))
    sl=(f>=3)&(f<=25);ff=f[sl];gg=G[sl]
    # Each run's U/Y pair is resampled together; all frequencies share the same
    # draw to preserve within-period cross-frequency correlation.
    u=np.stack([a['U'][:pp,sl] for a in blocks]);y=np.stack([a['Y'][:pp,sl] for a in blocks]);d=np.stack([a['D'][:pp,sl] for a in blocks])
    draws=[];bad=0
    for j in range(B):
        idx=rng.integers(pp,size=(4,pp))
        ub=np.stack([u[r,idx[r]].mean(0) for r in range(4)],axis=2)
        yb=np.stack([y[r,idx[r]].mean(0) for r in range(4)],axis=2)
        try:gb=right(yb,ub)
        except np.linalg.LinAlgError:gb=np.full_like(gg,np.nan);bad+=1
        draws.append(gb)
    draws=np.asarray(draws);delta=draws-draws.mean(0)
    vc=delta.reshape(B,len(ff),16)
    cov=np.einsum('bki,bkj->kij',vc,vc.conj())/(B-1)
    se=np.sqrt(np.real(np.diagonal(cov,axis1=1,axis2=2))).reshape(gg.shape)
    mag=20*np.log10(np.maximum(np.abs(draws),1e-300));ph=np.angle(draws/gg[None])*180/np.pi
    n1=pp//2
    g1=right(np.stack([v[:n1].mean(0) for v in y],axis=2),np.stack([v[:n1].mean(0) for v in u],axis=2))
    g2=right(np.stack([v[n1:].mean(0) for v in y],axis=2),np.stack([v[n1:].mean(0) for v in u],axis=2))
    # Linear delta-method covariance in real/imag form preserves non-proper
    # sample behavior and separates finite-ratio curvature from sampling spread.
    U=mu['U'][sl];Y=mu['Y'][sl];inv=np.linalg.inv(U);CV=np.zeros((len(ff),16,16),complex)
    CR=np.zeros((len(ff),32,32))
    for r in range(4):
        er=y[r]-y[r].mean(0)-np.einsum('koi,pki->pko',gg,u[r]-u[r].mean(0))
        dg=np.einsum('pko,ki->pkoi',er,inv[:,r,:]).reshape(pp,len(ff),16)
        CV+=np.einsum('pki,pkj->kij',dg,dg.conj())/(pp*(pp-1))
        dr=np.concatenate([dg.real,dg.imag],axis=2)
        CR+=np.einsum('pki,pkj->kij',dr,dr)/(pp*(pp-1))
    summary=dict(controller=ctrl,P_common=pp,P_each=[v['U'].shape[0] for v in blocks],n_draws=B,bootstrap_failed=bad,
        D_condition_max=float(np.max(np.linalg.cond(mu['D'][sl]))),U_condition_max=float(np.max(np.linalg.cond(U))),
        G_halves_difference_percent=float(100*np.linalg.norm(g1-g2)/np.linalg.norm(gg)),
        random_se_relative_percent=float(100*np.linalg.norm(se)/np.linalg.norm(gg)),
        bootstrap_mean_shift_percent=float(100*np.linalg.norm(draws.mean(0)-gg)/np.linalg.norm(gg)),
        instantaneous_gated_reference=True,interval='pointwise paired-period bootstrap percentile; not simultaneous, excludes systematic/model error')
    np.savez_compressed(HERE/'analysis'/f'{ctrl}_orthogonal_frf.npz',f=ff,G=gg,U=U,Y=Y,D=mu['D'][sl],TU=TU[sl],TY=TY[sl],
        Uper=u,Yper=y,Dper=d,G_boot=draws,se_complex=se,cov_complex=cov,cov_delta=CV,cov_real_delta=CR,
        mag_ci=np.percentile(mag,[2.5,97.5],axis=0),phase_ci_relative=np.percentile(ph,[2.5,97.5],axis=0),G_front=g1,G_back=g2,P_common=pp)
    js(HERE/'analysis'/f'{ctrl}_orthogonal_frf.json',summary);print('FRF',summary,flush=True)
    return summary

def compare():
    rows=[]
    for i,a in enumerate(CTRS):
        for b in CTRS[i+1:]:
            aa=np.load(HERE/'analysis'/f'{a}_orthogonal_frf.npz');bb=np.load(HERE/'analysis'/f'{b}_orthogonal_frf.npz')
            diff=aa['G']-bb['G'];se2=aa['se_complex']**2+bb['se_complex']**2
            # Global standardized squared discrepancy, calibrated by centered
            # bootstrap curves rather than pretending frequencies are independent.
            T=float(np.mean(abs(diff)**2/np.maximum(se2,1e-300)))
            dev=(aa['G_boot']-aa['G_boot'].mean(0))-(bb['G_boot']-bb['G_boot'].mean(0))
            Tboot=np.mean(abs(dev)**2/np.maximum(se2,1e-300),axis=(1,2,3))
            crit=float(np.percentile(Tboot,95));p=(1+np.count_nonzero(Tboot>=T))/(len(Tboot)+1)
            # Marginal complex-radius bootstrap is a pointwise check only.
            rad=np.percentile(abs(dev),95,axis=0);outside=np.abs(diff)>rad
            rows.append(dict(controller_a=a,controller_b=b,relative_difference_percent=float(100*np.linalg.norm(diff)/np.linalg.norm((aa['G']+bb['G'])/2)),
                global_standardized_statistic=T,centered_bootstrap_95_threshold=crit,bootstrap_p_value=p,compatible_global_95=bool(T<=crit),
                pointwise_outside_fraction=float(outside.mean()),interpretation='conditional repeatability calibration; not exact finite-sample coverage'))
    js(HERE/'analysis/orthogonal_controller_comparison.json',rows)
    return rows
if __name__=='__main__':
    for c in CTRS:estimate(c)
    print(compare(),flush=True)
