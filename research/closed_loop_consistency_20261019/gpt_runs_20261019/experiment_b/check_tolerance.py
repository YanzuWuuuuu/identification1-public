"""Additional diagnostics; does not replace any primary MC result."""
from run_b import *

results=[]
for rep in (0,19):
    rng=np.random.default_rng(52061019+rep)
    d,u,y=c.simulate_periodic(30000,10,.2,3e-7,rng,settle=2)
    k=np.arange(180,48001); z=np.exp(2j*np.pi*k/300000)
    U=np.fft.fft(u,axis=0)[k];Y=np.fft.fft(y,axis=0)[k];D=np.fft.fft(d,axis=0)[k];K=c.Ktf(z)
    bt=c.truth_beta()
    for method in METHODS:
        fn=HERE/'raw'/f'rep{rep:03d}_{method}.npz'
        retry=HERE/'raw'/f'rep{rep:03d}_{method}_retry1.npz'
        if retry.exists(): fn=retry
        with np.load(fn) as a: beta=a['beta']
        sel=(k%10==0) if '_lines_' in method else np.ones(len(k),bool)
        zs=z[sel];us=U[sel];ys=Y[sel];ds=D[sel];ks=K[sel]
        uh=(lambda b:np.linalg.solve(np.eye(4)[None]+ks@c.resp(zs,b),ds[...,None])[...,0]) if method.startswith('ref') else us
        t0=time.perf_counter(); b,status,hist,sc,obj=solve(zs,us,ys,uh,beta,tol=1e-8,maxit=20)
        diff=100*(np.exp(b[c.IW:c.ITR]-bt[c.IW:c.ITR])-np.exp(beta[c.IW:c.ITR]-bt[c.IW:c.ITR]))
        meta=dict(rep=rep,method=method,status=status,score=sc,iterations=len(hist),seconds=time.perf_counter()-t0,
                  parameter_percent_point_change=diff.tolist(),max_abs_change=float(max(abs(diff))))
        results.append(meta)
        atomic_npz(HERE/f'tolerance_rep{rep:03d}_{method}.npz',beta=b,primary_beta=beta,history=hist,metadata_json=np.array(json.dumps(meta)))
        atomic_json(HERE/'tolerance_check.json',results)
        print(json.dumps(meta),flush=True)
