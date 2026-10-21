from profile_core import *
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,traceback

def fit_one(rid,w,kind,method,omit=None,fixed=None):
    name=f'{rid}_{w}_{kind}_{method}_'+('full' if omit is None else f'omit{omit}')
    if fixed is not None:name+=f'_fixed{fixed:.6f}'
    dest=ROOT/'raw'/f'{name}.json'
    if dest.exists():return rd(dest)
    start=time.perf_counter();r=dict(id=rid,window=w,model=kind,method=method,omit=omit,common_fixed_us=fixed,success=False)
    try:
        m=ResidueModel(w,kind,common_fixed_us=fixed)
        fb,orig=load_joint(rid,w,method)
        # Failed low-window delay endpoints are not initialization evidence.
        if w=='L':
            fb0,_=old.basefit(rid,w,method);fb=np.r_[fb0,np.zeros(7)]
        b=m.from_full(fb)
        a=old.spectrum(rid,w,omit);pr=ProfileIV(m,a,method)
        grid=[]
        if omit is not None:
            parent=ROOT/'raw'/f'{rid}_{w}_{kind}_{method}_full.json'
            if fixed is not None:parent=ROOT/'raw'/f'{rid}_{w}_{kind}_{method}_full_fixed{fixed:.6f}.json'
            p=rd(parent);m.pivots=p['pivots'];b=np.asarray(p['beta'])
        elif fixed is None:
            for delta in np.arange(-200.,201.,10.):
                bg=b.copy();bg[m.delay_slice.start]=delta/100.
                bg,ok,eq,ninner=pr.profile(bg[m.nl:],bg[:m.nl],maxinner=25)
                A,g,sc,ee,zz=pr.eq(bg)
                grid.append(dict(delta_us=delta,ee=ee,score=sc,inner_ok=ok,inner_iterations=ninner,beta=bg.tolist()))
            valid=[v for v in grid if v['inner_ok'] and np.isfinite(v['ee'])]
            if valid:b=np.asarray(min(valid,key=lambda x:x['ee'])['beta'])
        warm=omit is not None and os.environ.get('REPARAM_WARM_HYBRID')=='1'
        maxouter=(20 if kind in ('a','b') else 0) if warm else 90
        bout,info=pr.solve(b,maxouter=maxouter,tol=1e-10 if warm and kind in ('a','b') else 2e-8)
        f,z=m.modal(bout);fbout,_=m.expand(bout)
        r.update(beta=bout,full_beta=fbout,pivots=m.pivots,solver=info,success=info['success'],
                 frequency_hz=f,damping_percent=z,delay=old.delay_values(m.delay(bout)),
                 parameter_count=m.nb,linear_count=m.nl,nonlinear_count=m.nq,grid=grid,
                 numerical_route=('profile hybrid from full-data endpoint' if maxouter==0 else f'profile GN ({maxouter} max) then hybrid'))
    except Exception as exc:r.update(exception=repr(exc),traceback=traceback.format_exc())
    r['seconds']=time.perf_counter()-start;js(dest,r)
    print(name,r['success'],r.get('solver',{}).get('score'),round(r['seconds'],2),flush=True)
    return r

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--windows',default='L,H300,H493');ap.add_argument('--models',default='a,b,c,d')
    ap.add_argument('--ids',default=','.join(old.IDS));ap.add_argument('--methods',default='M3,M1');ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--jackknife',action='store_true');ap.add_argument('--warm-hybrid',action='store_true');args=ap.parse_args()
    if args.warm_hybrid:os.environ['REPARAM_WARM_HYBRID']='1'
    jobs=[]
    for rid in args.ids.split(','):
        for w in args.windows.split(','):
            for kind in args.models.split(','):
                for method in args.methods.split(','):
                    if args.jackknife:
                        full=ROOT/'raw'/f'{rid}_{w}_{kind}_{method}_full.json'
                        if not full.exists():continue
                        P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])
                        # Keep failed parent endpoints visible; do not claim their uncertainty.
                        if not rd(full).get('success'):continue
                        jobs.extend((rid,w,kind,method,i) for i in range(P))
                    else:jobs.append((rid,w,kind,method))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for fu in as_completed([pool.submit(fit_one,*j) for j in jobs]):fu.result()
