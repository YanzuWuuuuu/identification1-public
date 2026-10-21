from profile_core import *

def model_from(r):
    m=ResidueModel(r['window'],r['model'],common_fixed_us=r.get('common_fixed_us'),pivots=r['pivots'])
    return m,np.asarray(r['beta'])

def covariance(rid,w,a):
    """A fixed common covariance for the four residue models."""
    name=ROOT/'analysis'/f'noise_{rid}_{w}.npz'
    if name.exists():
        z=np.load(name);return z['Phi'],str(z['source'])
    path=ROOT/'raw'/f'{rid}_{w}_a_M3_full.json'
    r=rd(path);ext=old.spectrum(rid,w,extended=True)
    if r['success']:
        m,b=model_from(r);fb=m.expand(b,False)[0]
        Phi=old.local_phi(m.full,ext,fb,a['f']);source='profiled model (a), joint delays'
    else:
        m=old.setup(w)[1];b,_=old.basefit(rid,w,'M3')
        Phi=old.local_phi(m,ext,b,a['f']);source='original fixed-delay model (a); joint-delay parent failed'
    np.savez_compressed(name,Phi=Phi,f=a['f'],source=source)
    return Phi,source

def residual_test(m,a,b,Phi,method='M3'):
    # Evaluate on the excited lines, for all models, with common frozen covariance.
    if method!='M3':
        raise ValueError('This residual expectation is for line-fitted M3 only; all-bin M1 needs a different residual operator.')
    li=a['line_indices'];c=m.full.context(a['f'][li]);U,Y,D,K=[a[k][li] for k in ('U','Y','D','K')]
    Uh=U if method=='M1' else np.linalg.solve(np.eye(4)[None]+K@m.resp(c,b),D[:,:,None])[:,:,0]
    yh,J=m.output(c,b,U);_,Z=m.output(c,b,Uh)
    Jr=old.realvec(J);Zr=old.realvec(Z);e=old.realvec(Y-yh)
    L=np.linalg.cholesky(.5*old.realmat(a['N']*Phi[li]))
    ev=np.linalg.solve(L,e[:,:,None])[:,:,0]
    B=np.linalg.solve(L,Jr).reshape(-1,m.nb)
    A=np.einsum('kop,koq->pq',Zr,Jr,optimize=True)
    ZL=np.einsum('kop,kol->pkl',Zr,L,optimize=True).reshape(m.nb,-1)
    sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
    H=sc[:,None]*la.lstsq(A*sc[:,None]*sc[None,:],sc[:,None]*ZL,cond=1e-13,lapack_driver='gelsy')[0]
    scale=np.maximum(np.linalg.norm(B,axis=0),1e-300)
    qb,rb=np.linalg.qr(B/scale,mode='reduced');rh=rb@(H*scale[:,None])
    # Diagonal of Rw Rw^T, not the diagonal of an orthogonal LS projector.
    ed=1-2*np.einsum('ip,pi->i',qb,rh)+np.einsum('ip,pq,iq->i',qb,rh@rh.T,qb,optimize=True)
    ed=ed.reshape(len(li),8).sum(1);qline=(ev**2).sum(1)
    bands={'all':np.ones(len(li),bool)}
    if m.w=='H300':
        bands.update(mode333=(a['f'][li]>=307.5)&(a['f'][li]<346),mode359=a['f'][li]>=346)
    out={k:dict(Q=float(qline[v].sum()),EQ=float(ed[v].sum()),ratio=float(qline[v].sum()/ed[v].sum())) for k,v in bands.items()}
    out['negative_expected_entries']=int(sum(ed<0))
    return out

def prediction(m,a,b,Phi):
    c=m.full.context(a['f']);G,Gj=m.gmodel(c,b);C=a['N']*Phi
    F=np.linalg.solve(np.eye(4)[None]+a['K']@G,a['K'])
    Ud=np.linalg.solve(np.eye(4)[None]+a['K']@G,a['D'][:,:,None])[:,:,0]
    _,J=m.output(c,b,Ud)
    M0=np.einsum('kop,koq->pq',J.conj(),J,optimize=True).real
    chol=np.linalg.cholesky(C)
    nj=-np.einsum('koip,kil->kolp',Gj,F@chol,optimize=True).reshape(-1,m.nb)
    M=M0+(nj.conj().T@nj).real
    bv=-np.einsum('koip,kil,klo->p',Gj,F,C,optimize=True).real
    d=scaled_solve(M,bv)
    return dict(delta_log_parameters=d[m.nl:m.nl+2*m.nm],prediction_percent=100*np.expm1(d[m.nl:m.nl+2*m.nm]),
                condition_scaled=float(np.linalg.cond(M/np.sqrt(abs(np.diag(M)))[:,None]/np.sqrt(abs(np.diag(M)))[None,:])))

def endpoint_metric(path):
    dest=ROOT/'analysis'/('metric_'+path.name)
    if dest.exists():return rd(dest)
    r=rd(path);out={k:r.get(k) for k in ('id','window','model','method','success','omit')}
    try:
        m,b=model_from(r);a=old.spectrum(r['id'],r['window'])
        Phi,source=covariance(r['id'],r['window'],a)
        out.update(qe=residual_test(m,a,b,Phi,r['method']),covariance_source=source)
        # Instrument-whitened common delay direction remaining after other parameters.
        pr=ProfileIV(m,a,r['method']);_,J=m.output(pr.c,b,pr.uh(b))
        X=old.realvec(J).reshape(-1,m.nb)
        delayidx=list(range(m.delay_slice.start,m.delay_slice.stop))
        scales=np.linalg.norm(X,axis=0);Xs=X/np.maximum(scales,1e-300)
        if delayidx:
            rest=[i for i in range(m.nb) if i not in delayidx]
            proj=Xs[:,delayidx]-Xs[:,rest]@la.lstsq(Xs[:,rest],Xs[:,delayidx],cond=1e-12)[0]
            out['delay_residual_fraction']=np.sum(proj**2,axis=0)
            out['delay_conditional_singular_values']=la.svdvals(proj)
    except Exception as e:out['exception']=repr(e)
    js(dest,out);return out

if __name__=='__main__':
    for p in (ROOT/'raw').glob('*_M3_full.json'):
        if not p.name.startswith('p7_'):
            r=endpoint_metric(p)
            print(p.name,r.get('qe',{}).get('all',{}).get('ratio'),r.get('exception'),flush=True)
