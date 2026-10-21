from core import *
from concurrent.futures import ProcessPoolExecutor,as_completed

def fixed(rid,w,shift,method):
    path=ROOT/'raw'/f'fixed_{rid}_{w}_{shift:+g}_{method}.json'
    if path.exists():return rd(path)
    t=time.perf_counter();r=dict(id=rid,window=w,shift_us=shift,method=method)
    try:
        s,m,lb,ub=setup(w,shift);a=spectrum(rid,w);b0,_=basefit(rid,w,method)
        c,U,Y,Uh=inputs(m,a,method);b,info=solve(m,c,U,Y,Uh,b0,lb,ub)
        f,z=m.modal(b);f0,z0=m.modal(b0)
        r.update(beta=b,solver=info,success=info['success'],frequency_hz=f,damping_percent=z,
                 delta_frequency_hz=f-f0,delta_damping_pp=z-z0)
    except Exception as e:r.update(success=False,exception=repr(e),traceback=traceback.format_exc())
    r['seconds']=time.perf_counter()-t;js(path,r);print('FIX',rid,w,shift,method,r['success'],round(r['seconds'],1),flush=True);return r

def identifiability(m,c,b,U,Uh):
    _,J=m.output(c,b,U);_,Z=m.output(c,b,Uh(b) if callable(Uh) else Uh)
    Jr=realvec(J).reshape(-1,m.nb);Zr=realvec(Z).reshape(-1,m.nb)
    sc=1/np.maximum(np.linalg.norm(Jr,axis=0),1e-300)
    sv=la.svdvals(Jr*sc);bn=m.base_nb
    Q,_=np.linalg.qr(Jr[:,:bn]*sc[:bn],mode='reduced')
    D=Jr[:,bn:];Dp=D-Q@(Q.T@D)
    sdelay=la.svdvals(Dp);raw=la.svdvals(D)
    return dict(jacobian_scaled_condition=sv[0]/sv[-1],joint_rank_1e12=int(sum(sv>sv[0]*1e-12)),
                delay_projected_singular_values=sdelay,delay_raw_singular_values=raw,
                delay_remaining_column_energy=np.sum(Dp**2,axis=0)/np.sum(D**2,axis=0))

def joint(rid,w,method,omit=None):
    tag='full' if omit is None else f'omit{omit}'
    path=ROOT/'raw'/f'joint_{rid}_{w}_{method}_{tag}.json'
    if path.exists():return rd(path)
    t=time.perf_counter();r=dict(id=rid,window=w,method=method,omit=omit)
    try:
        s,m,lb,ub=setup(w,joint=True);a=spectrum(rid,w,omit);c,U,Y,Uh=inputs(m,a,method)
        base,_=basefit(rid,w,method);zero=np.r_[base,np.zeros(7)]
        if omit is None:
            # Evaluate every requested common-delay grid start, profiling linear
            # coefficients at fixed modal poles. The grid is initialization only.
            grid=[];starts=[]
            for us in np.arange(-200,201,10):
                b=zero.copy();b[m.base_nb]=us/100
                b=linsolve(m,c,U,Y,Uh,b)
                _,_,score,ee=equations(m,c,U,Y,Uh(b) if callable(Uh) else Uh,b)
                grid.append(dict(offset_us=int(us),score=score,ee=ee,beta=b))
            # Residual energy selects starting points, not the final IV estimate.
            # Include baseline plus two separated profile minima for robustness.
            starts=[zero];chosen=[0]
            for g in sorted(grid,key=lambda x:x['ee']):
                if min(abs(g['offset_us']-v) for v in chosen)>=50:
                    starts.append(np.asarray(g['beta']));chosen.append(g['offset_us'])
                if len(starts)==3:break
            r['grid']=grid;r['start_offsets_us']=chosen
        else:
            full=rd(ROOT/'raw'/f'joint_{rid}_{w}_{method}_full.json')
            starts=[np.asarray(full['beta'])];r['full_endpoint_success']=full['success']
        runs=[]
        for ix,b0 in enumerate(starts):
            b,info=solve(m,c,U,Y,Uh,b0,lb,ub,maxit=120)
            runs.append(dict(start=ix,beta=b,solver=info))
        good=[q for q in runs if q['solver']['success']]
        selected=min(good or runs,key=lambda q:q['solver'].get('score',float('inf')))
        b=np.asarray(selected['beta']);info=selected['solver'];f,z=m.modal(b)
        r.update(beta=b,success=info['success'],solver=info,all_starts=runs,selected_start=selected['start'],
                 frequency_hz=f,damping_percent=z,delay=delay_values(m.delay(b)),P=a['P'])
        r['identifiability']=identifiability(m,c,b,U,Uh)
        if omit is None:
            _,mf,_,_=setup(w);ext=spectrum(rid,w,extended=True);bref,_=basefit(rid,w,'M3')
            Phi=local_phi(mf,ext,bref,a['f'])
            r['qe_fixed_noise']=qe(m,a,b,Phi)
            r['qe_baseline']=qe(mf,a,bref,Phi)
            Phi_new=local_phi(m,ext,b,a['f'])
            r['qe_reestimated_noise']=qe(m,a,b,Phi_new)
            np.savez_compressed(ROOT/'raw'/f'noise_{rid}_{w}_{method}.npz',f=a['f'],Phi_baseline=Phi,Phi_joint=Phi_new)
    except Exception as e:r.update(success=False,exception=repr(e),traceback=traceback.format_exc())
    r['seconds']=time.perf_counter()-t;js(path,r)
    print('JOINT',rid,w,method,tag,r['success'],r.get('solver',{}).get('score'),round(r['seconds'],1),flush=True);return r

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('kind',choices=['fixed','joint','folds']);p.add_argument('--id');p.add_argument('--window')
    p.add_argument('--method');p.add_argument('--workers',type=int,default=4)
    args=p.parse_args();jobs=[];ids=[args.id] if args.id else IDS;ws=[args.window] if args.window else SPECS
    for rid in ids:
        for w in ws:
            for method in ([args.method] if args.method else ['M1','M3']):
                if args.kind=='fixed':
                    jobs.extend((rid,w,sh,method) for sh in [-100,-50,-25,0,25,50,100])
                elif args.kind=='joint':jobs.append((rid,w,method,None))
                else:jobs.extend((rid,w,method,i) for i in range(int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])))
    fn=fixed if args.kind=='fixed' else joint
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(fn,*job) for job in jobs]):future.result()
