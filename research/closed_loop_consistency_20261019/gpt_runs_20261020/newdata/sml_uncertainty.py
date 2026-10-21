from identify import *

def one(rid,B=99):
    a=prepare_spectra(rid);li=a['line_indices'];P=int(a['P']);f=a['f'][li]
    fit=json.loads((HERE/'fits'/f'{rid}_M5n.json').read_text());b0=np.asarray(fit['beta'])
    cv=np.load(HERE/'fits'/f'{rid}_local_covariance.npz');c=m.Ctx(f,TAU,'zoh')
    model=lambda b,u:m.out_jac(c,b,u);gm=gmodel_from_output(model,m.NB,len(f))
    rng=np.random.default_rng(20261020+sum(map(ord,rid)));out=HERE/'analysis/sml_bootstrap';out.mkdir(exist_ok=True)
    rows=[]
    for j in range(B):
        path=out/f'{rid}_{j:03d}.json'
        ix=rng.integers(P,size=P)
        if path.exists():rows.append(json.loads(path.read_text()));continue
        U=P*a['Uper'][ix].mean(0);Y=P*a['Yper'][ix].mean(0)
        b,si=solve_sml_cov(U,Y,cv['CY'],cv['CU'],cv['CYU'],b0,model,gm,maxiter=160)
        fz,zz=m.modal(b);r=dict(draw=j,success=bool(si['success']),frequency_hz=fz.tolist(),damping_percent=zz.tolist(),
                             modal_log=b[m.IW:m.ITR].tolist(),solver=si,indices=ix.tolist())
        js(path,r);rows.append(r)
    good=[r for r in rows if r['success'] and np.isfinite(r['modal_log']).all()]
    vals=np.array([r['modal_log'] for r in good]);p=np.cov(vals,rowvar=False);sd=vals.std(0,ddof=1)
    f,z=m.modal(b0);real=np.exp(vals);real[:,:3]/=2*np.pi;real[:,3:]*=100
    result=dict(id=rid,method='M5n',frequency_u_hz=real[:,:3].std(0,ddof=1).tolist(),damping_u_percentage_points=real[:,3:].std(0,ddof=1).tolist(),
        log_parameter_sd=sd.tolist(),n_bootstrap=B,failed_bootstraps=B-len(good),locally_valid_numerics=fit['converged'] and B==len(good),
        interpretation='conditional paired-period bootstrap; nonlinear SML refit; original joint noise covariance frozen; not theorem 2; no independent known-noise coverage guarantee')
    js(HERE/'analysis'/f'{rid}_M5n_uncertainty.json',result);print('SML_BOOT',rid,'failed',B-len(good),flush=True)
    return rid
if __name__=='__main__':
    ids=[p.stem[:-4] for p in (HERE/'fits').glob('*_M5n.json')]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(one,rid) for rid in ids]):f.result()
