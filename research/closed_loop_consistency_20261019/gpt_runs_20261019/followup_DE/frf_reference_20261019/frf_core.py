"""Reference-driven MIMO local-polynomial FRF, with no parametric model input."""
import numpy as np


def local_reference_frf(f, Dper, Uper, Yper, order=2, window=21):
    """Arrays per: periods x frequencies x channels; f in Hz.

    Fits 8 closed-loop channels jointly with exactly the same complex design.
    All regressors depend only on measured periodic reference D and frequency.
    Each period is fitted separately with the fixed mean reference design,
    then the closed-loop transfer estimates are averaged before inversion.
    """
    P, nf, nu = Dper.shape
    assert nu == 4 and Uper.shape == Yper.shape == Dper.shape
    assert window <= nf and window % 2 == 1
    D = Dper.mean(0)
    O = np.concatenate([Uper, Yper], axis=-1)
    nc = (nu+1)*(order+1)
    assert window > nc
    TU = np.empty((P,nf,4,4), complex)
    TY = np.empty_like(TU)
    G = np.full((nf,4,4), np.nan+1j*np.nan)
    influence = np.full((P,nf,4,4), np.nan+1j*np.nan)
    design_condition=np.empty(nf); design_rank=np.empty(nf,int)
    tud_condition=np.empty(nf); edge=np.empty(nf,bool)
    pred = np.empty_like(O); residual_scale=np.empty((nf,8))
    bounds=np.empty((nf,2),int)
    for k in range(nf):
        start=min(max(k-window//2,0),nf-window)
        idx=np.arange(start,start+window); bounds[k]=[idx[0],idx[-1]]
        edge[k]=k-window//2<0 or k+window//2>=nf
        x=(f[idx]-f[k])/max(np.max(np.abs(f[idx]-f[k])),1e-100)
        X=np.column_stack([D[idx]*x[:,None]**p for p in range(order+1)]+
                          [x[:,None]**p for p in range(order+1)])
        scale=np.linalg.norm(X,axis=0)
        Xs=X/scale
        s=np.linalg.svd(Xs,compute_uv=False)
        design_condition[k]=s[0]/s[-1]
        design_rank[k]=np.count_nonzero(s>s[0]*1e-12)
        # Column stacking preserves period/output layout.
        target=O[:,idx,:].transpose(1,0,2).reshape(window,P*8)
        coef=np.linalg.lstsq(Xs,target,rcond=1e-12)[0]/scale[:,None]
        coeff=coef.reshape(nc,P,8)
        TU[:,k]=coeff[:4,:,:4].transpose(1,2,0)
        TY[:,k]=coeff[:4,:,4:].transpose(1,2,0)
        xc=np.r_[D[k],np.zeros(nu*order),1.,np.zeros(order)]
        pred[:,k]=np.einsum('a,apo->po',xc,coeff)
        resid=target-X@coef
        residual_scale[k]=np.sqrt(np.mean(np.abs(resid.reshape(window,P,8))**2,axis=(0,1)))
        tu=TU[:,k].mean(0); ty=TY[:,k].mean(0)
        tud_condition[k]=np.linalg.cond(tu)
        if design_rank[k]<nc or not np.isfinite(tud_condition[k]):
            continue
        try:
            G[k]=np.linalg.solve(tu.T,ty.T).T
            for p in range(P):
                delta=(TY[p,k]-ty)-G[k]@(TU[p,k]-tu)
                influence[p,k]=np.linalg.solve(tu.T,delta.T).T
        except np.linalg.LinAlgError:
            pass
    se=np.sqrt(np.sum(np.abs(influence)**2,axis=0)/(P*(P-1))) if P>1 else np.full(G.shape,np.nan)
    relative=influence/np.where(np.abs(G)>1e-300,G,complex(np.nan))
    semag=20/np.log(10)*np.sqrt(np.sum(relative.real**2,axis=0)/(P*(P-1))) if P>1 else np.full(G.shape,np.nan)
    sephase=180/np.pi*np.sqrt(np.sum(relative.imag**2,axis=0)/(P*(P-1))) if P>1 else np.full(G.shape,np.nan)
    return dict(f=f,G=G,TU_period=TU,TY_period=TY,TU=TU.mean(0),TY=TY.mean(0),
                influence=influence,se_complex=se,se_db=semag,se_phase_deg=sephase,
                design_condition=design_condition,design_rank=design_rank,
                tud_condition=tud_condition,edge=edge,bounds=bounds,
                prediction=pred,closed_loop_residual_scale=residual_scale,
                order=np.array(order),window=np.array(window),nperiod=np.array(P),
                n_coefficients=np.array(nc))


def metrics(Gmodel, Gref):
    delta=Gmodel-Gref
    norm=lambda x: np.sqrt(np.sum(np.abs(x)**2))
    mag=20*np.log10(np.maximum(np.abs(Gmodel),1e-300)/np.maximum(np.abs(Gref),1e-300))
    phase=np.angle(Gmodel*np.conj(Gref))*180/np.pi
    return dict(complex_relative_percent=100*norm(delta)/norm(Gref),
                channel_relative_percent=(100*np.sqrt(np.sum(np.abs(delta)**2,axis=0)/
                                                       np.sum(np.abs(Gref)**2,axis=0))).tolist(),
                magnitude_rmse_db=float(np.sqrt(np.mean(mag**2))),
                phase_rmse_deg=float(np.sqrt(np.mean(phase**2))))
