from lag_core import *
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,traceback

def residual(pr,b):
    m=pr.m;pieces=[];ranges={};start=0
    for w,a in pr.data.items():
        y,J=m.output(w,pr.ctx[w],b,a['U']);_,Z=m.output(w,pr.ctx[w],b,pr.uh_window(w,b))
        Jr=old.realvec(J);Zr=old.realvec(Z);e=old.realvec(a['Y']-y)
        L=np.linalg.cholesky(.5*old.realmat(a['N']*a['Phi']))
        B=np.linalg.solve(L,Jr).reshape(-1,m.nb)
        ev=np.linalg.solve(L,e[:,:,None])[:,:,0].ravel()
        ZL=np.einsum('kop,kol->pkl',Zr,L,optimize=True).reshape(m.nb,-1)
        A=np.einsum('kop,koq->pq',Zr,Jr,optimize=True)
        pieces.append((B,ev,ZL,A));ranges[w]=slice(start,start+len(ev));start+=len(ev)
    B=np.concatenate([p[0] for p in pieces]);e=np.concatenate([p[1] for p in pieces])
    ZL=np.concatenate([p[2] for p in pieces],axis=1);A=sum(p[3] for p in pieces)
    sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
    H=sc[:,None]*la.lstsq(A*sc[:,None]*sc[None,:],sc[:,None]*ZL,cond=1e-13,lapack_driver='gelsy')[0]
    # Rw = I - B H; row-block expectations include all windows through H H^T.
    scale=np.maximum(np.linalg.norm(B,axis=0),1e-300)
    qb,rb=la.qr(B/scale,mode='economic');rh=rb@(H*scale[:,None])
    ed=1-2*np.einsum('ip,pi->i',qb,rh)+np.einsum('ip,pq,iq->i',qb,rh@rh.T,qb,optimize=True)
    result={}
    for w,s in ranges.items():
        Q=float(e[s]@e[s]);EQ=float(ed[s].sum())
        result[w]=dict(Q=Q,EQ=EQ,ratio=Q/EQ,n=len(e[s]),expected_noise_amplification=EQ/len(e[s]))
    result['joint']=dict(Q=float(e@e),EQ=float(ed.sum()),ratio=float(e@e/ed.sum()),n=len(e))
    # Joint normalized instrument information after all modal/linear columns.
    xn=np.concatenate([old.realvec(m.output(w,pr.ctx[w],b,pr.uh_window(w,b))[1]).reshape(-1,m.nb)
                       for w in WINDOWS])
    xn/=np.maximum(np.linalg.norm(xn,axis=0),1e-300)
    j=np.arange(m.delay_slice.start,m.nb);rest=np.arange(m.delay_slice.start)
    proj=xn[:,j]-xn[:,rest]@la.lstsq(xn[:,rest],xn[:,j],cond=1e-12,lapack_driver='gelsy')[0]
    result['joint_delay_lag_singular_values']=la.svdvals(proj)
    result['negative_expected_rows']=int(np.sum(ed<0))
    return result

def fit(rid,kind):
    dest=ROOT/'raw'/f'{rid}_{kind}_joint.json'
    if dest.exists():return rd(dest)
    begin=time.perf_counter()
    r=dict(id=rid,model=kind,success=False,grid=[],status='started')
    try:
        m=JointModel(kind);b=m.initial(rid);pr=JointSolver(m,get_data(rid))
        for shift in np.arange(-200.,201.,10.):
            candidate=b.copy();candidate[m.delay_slice.start]+=shift/100.
            candidate,ok,_,nit=pr.profile(candidate[m.nl:],candidate[:m.nl],maxinner=30)
            _,_,score,ee,_=pr.eq(candidate)
            r['grid'].append(dict(shift_us=shift,inner_success=ok,inner_iterations=nit,score=score,
                                  residual_squared=ee,beta=candidate))
        good=[g for g in r['grid'] if g['inner_success'] and np.isfinite(g['residual_squared'])]
        if not good:raise RuntimeError('All common-delay initialization grid points failed the inner moment equation.')
        initial=np.asarray(min(good,key=lambda x:x['residual_squared'])['beta'])
        js(dest.with_name(dest.stem+'_initialization.json'),r)
        best,info=pr.solve(initial,maxouter=55,tol=2e-8)
        r.update(beta=best,solver=info,success=info['success'],parameter_count=m.nb,
                 linear_count=m.nl,nonlinear_count=m.nq,**m.unpack(best))
        r['qe']=residual(pr,best)
        r['status']='complete'
    except Exception as exc:r.update(status='exception',exception=repr(exc),traceback=traceback.format_exc())
    r['seconds']=time.perf_counter()-begin;js(dest,r)
    print(rid,kind,r['success'],r.get('solver',{}).get('score'),r.get('cutoff_hz'),round(r['seconds'],1),flush=True)
    return r
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--ids',default=','.join(old.IDS))
    ap.add_argument('--models',default='a,c');ap.add_argument('--workers',type=int,default=4);args=ap.parse_args()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for fu in as_completed([pool.submit(fit,rid,kind) for rid in args.ids.split(',') for kind in args.models.split(',')]):fu.result()
