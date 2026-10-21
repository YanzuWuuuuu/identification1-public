"""Task D paired simulation. All source artifacts are read-only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
os.environ.update(NOFLEX='1',FC='60',AL='30',AH='100',FHI='100',FLO='3',SW='.03',FSPR='.01',ZSPR='.4')
import sys,time,json,argparse
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent; OLD=HERE.parent; BASE=OLD.parent
sys.path.insert(0,str(BASE/'code'))
import numpy as np
import cl2 as c
from weighted_core import iv_solve,iv_score,precision,local_covariances,statistics,whitener,iv_matrices
from sml_core_copy import solve_sml,solve_sml_cov,gmodel_from_output,sample_covariances,SMLCovResidual,SMLResidual
from concurrent.futures import ProcessPoolExecutor,as_completed
OUT=HERE/'sim_d'; OUT.mkdir(exist_ok=True)
for sub in ('results','data','diagnostics'): (OUT/sub).mkdir(exist_ok=True)

def write_json(path,obj):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=True),encoding='utf-8');tmp.replace(path)

def load_data(P,r):
    name=f'P{P:02d}_r{r:03d}.npz'; source=OLD/'experiment_a'/'data'/name
    if source.exists(): return np.load(source)
    dest=OUT/'data'/name
    if dest.exists(): return np.load(dest)
    seed=202610190+1000*P+r; rng=np.random.default_rng(seed)
    d,u,y=c.simulate_periodic(30000,P,.2,3e-7,rng,settle=2)
    N=30000*P;k=np.arange(int(np.ceil(3*N/c.fs)),int(np.floor(100*N/c.fs))+1)
    kp=np.arange(18,601);pert=rng.standard_normal(6)
    tmp=dest.with_suffix('.tmp.npz')
    np.savez_compressed(tmp,D=np.fft.fft(d,axis=0)[k],U=np.fft.fft(u,axis=0)[k],Y=np.fft.fft(y,axis=0)[k],
        k=k,f=k*c.fs/N,zk=np.exp(2j*np.pi*k/N),pert=pert,
        Dper=np.fft.fft(d.reshape(P,30000,4),axis=1)[:,kp],Uper=np.fft.fft(u.reshape(P,30000,4),axis=1)[:,kp],Yper=np.fft.fft(y.reshape(P,30000,4),axis=1)[:,kp],
        fper=kp*c.fs/30000,zkper=np.exp(2j*np.pi*kp/30000),seed=seed,P=P)
    tmp.replace(dest);return np.load(dest)

def context(P,r):
    a=load_data(P,r);k=a['k'];sel=k%P==0;z=a['zk'][sel];U=a['U'][sel];Y=a['Y'][sel];D=a['D'][sel]
    K=c.Ktf(z);bt=c.truth_beta();pert=a['pert']
    b0=bt.copy();b0[c.IW:c.IZ]+=np.log(1+.01*pert[:3]);b0[c.IZ:c.ITR]+=.4*pert[3:];b0[:c.NLIN]=0;b0[c.ITR:]=0
    b0=c.init_lin(z,U,Y,np.ones_like(U.real),b0)
    model=lambda b,u:c.out_jac(z,b,u)
    model_full=lambda b,u:c.out_jac(a['zk'],b,u,need_jac=False)
    Uh=lambda b:np.linalg.solve(np.eye(4)[None]+K@c.resp(z,b),D[:,:,None])[:,:,0]
    return a,z,U,Y,D,K,bt,b0,model,model_full,Uh

def one(P,r,methods=None):
    start=time.perf_counter(); a,z,U,Y,D,K,bt,b0,model,model_full,Uh=context(P,r)
    stem=f'P{P:02d}_r{r:03d}'; base=dict(P=P,replicate=r,seed=int(a['seed']))
    def save(name,b,info,extra=None):
        row={**base,'method':name,**info,'beta':b.tolist(),'initial_beta':b0.tolist(),
             'error_percent':(100*np.expm1(b[c.IW:c.ITR]-bt[c.IW:c.ITR])).tolist(),
             'log_error':(b[c.IW:c.ITR]-bt[c.IW:c.ITR]).tolist()}
        if extra:row.update(extra)
        write_json(OUT/'results'/f'{stem}_{name}.json',row);return row
    def fit(name,instr,W,initial):
        path=OUT/'results'/f'{stem}_{name}.json'
        if path.exists():row=json.loads(path.read_text(encoding='utf-8'));return np.array(row['beta']),row
        b,info=iv_solve(model,U,Y,instr,W,initial,c.IW,c.ITR)
        return b,save(name,b,info,{'initial_beta':initial.tolist()})
    unity=np.ones_like(U.real)
    b2,r2=fit('M2',U,unity,b0)
    b3,r3=fit('M3',Uh,unity,b0)
    G0=c.resp(z,bt);N=30000*P;Pv=N*(3e-7**2*np.eye(4)[None]+c.SW**2*G0@G0.conj().transpose(0,2,1));Wo=precision(Pv)
    bo,ro=fit('M3o',Uh,Wo,b0)
    # Conditional sandwich given the realised excitation. Retain transient nuisance.
    theo=OUT/'diagnostics'/f'{stem}_oracle_theory.npz'
    if not theo.exists():
        utrue=Uh(bt); _,Z=model(bt,utrue);WZ=np.einsum('kij,kjp->kip',Wo,Z)
        A=np.real(Z.reshape(-1,c.NB).conj().T@WZ.reshape(-1,c.NB))
        Q=.5*np.real(np.einsum('kip,kij,kjq->pq',WZ.conj(),Pv,WZ,optimize=True))
        sc=1/np.sqrt(np.maximum(np.abs(np.diag(A)),1e-300));As=A*sc[:,None]*sc[None,:];Qs=Q*sc[:,None]*sc[None,:]
        ss=np.linalg.solve(As,Qs);C=np.linalg.solve(As,ss.T).T*sc[:,None]*sc[None,:]
        np.savez_compressed(theo,covariance=C,percent_std=100*np.sqrt(np.maximum(np.diag(C)[c.IW:c.ITR],0)),condition_scaled=np.linalg.cond(As))
    gm=gmodel_from_output(model,c.NB,len(z))
    for nl in (8,16):
        cov=local_covariances(a['k'],a['U'],a['Y'],b3,model_full,P,nl,.1)
        bw1,rw1=fit(f'M3w_n{nl}_o1',Uh,precision(cov['covariance_residual']),b3)
        cov2=local_covariances(a['k'],a['U'],a['Y'],bw1,model_full,P,nl,.1)
        bw2,rw2=fit(f'M3w_n{nl}_o2',Uh,precision(cov2['covariance_residual']),bw1)
        name=f'M5n_n{nl}';path=OUT/'results'/f'{stem}_{name}.json'
        if not path.exists():
            b,info=solve_sml_cov(U,Y,cov['CY'],cov['CU'],cov['CYU'],b3,model,gm,maxiter=120)
            info['iterations']=info.get('nfev',0)
            save(name,b,info,{'initial_beta':b3.tolist(),'covariance_source':'nonexcited_full_DFT','nl':nl,'alpha':.1})
        write_json(OUT/'diagnostics'/f'{stem}_n{nl}_neighbors.json',dict(P=P,nl=nl,median_span_Hz=float(np.median(cov['distances'].max(1))*c.fs/N),max_span_Hz=float(cov['distances'].max()*c.fs/N),counts=len(cov['indices'])))
    name='M5';path=OUT/'results'/f'{stem}_{name}.json'
    if not path.exists():
        bper=b0.copy();bper[c.ITR:]/=P
        shrink=.1 if P<5 else 0.
        b,info=solve_sml(a['Uper'],a['Yper'],bper,model,gm,maxiter=120,shrink=shrink)
        info['iterations']=info.get('nfev',0)
        save(name,b,info,{'initial_beta':bper.tolist(),'covariance_source':'period_mean','transient_beta_units':'period_mean_DFT'})
    return dict(P=P,r=r,seconds=time.perf_counter()-start)

def gate():
    P=16;r=0;a,z,U,Y,D,K,bt,b0,model,model_full,Uh=context(P,r)
    unity=np.ones_like(U.real);bg,ig=iv_solve(model,U,Y,Uh,unity,b0,c.IW,c.ITR)
    bold,it=c.iv_solve(z,U,Y,Uh,unity,b0,it=80,tol=1e-5)
    gm=gmodel_from_output(model,c.NB,len(z));bp=b0.copy();bp[c.ITR:]/=P
    CY,CU,CYU=sample_covariances(a['Uper'],a['Yper'],0.)
    op=SMLCovResidual(U/P,Y/P,CY,CU,CYU,model,gm)
    of=SMLCovResidual(U,Y,CY*P**2,CU*P**2,CYU*P**2,model,gm)
    rp,jp=op.calc(bp);rf,jf=of.calc(b0)
    dscale=np.ones(c.NB);dscale[c.ITR:]=P
    relres=np.linalg.norm(rp-rf)/np.linalg.norm(rp);reljac=np.linalg.norm(jp-jf*dscale)/np.linalg.norm(jp)
    bfitp,ip=solve_sml_cov(U/P,Y/P,CY,CU,CYU,bp,model,gm)
    bfitf,iff=solve_sml_cov(U,Y,CY*P**2,CU*P**2,CYU*P**2,b0,model,gm)
    out=dict(P=P,r=r,full_DFT_equals_P_periodmean_relative=float(np.linalg.norm(U-P*a['Uper'].mean(0))/np.linalg.norm(U)),
             iv_parity_max_modal_log_difference=float(np.max(abs(bg[c.IW:c.ITR]-bold[c.IW:c.ITR]))),iv_iterations_new=ig['iterations'],iv_iterations_original=it+1,
             covariance_scale_factor=P**2,residual_scale_relative=float(relres),jacobian_scale_relative=float(reljac),
             estimate_scale_modal_log_difference=float(np.max(abs(bfitp[c.IW:c.ITR]-bfitf[c.IW:c.ITR]))),periodmean_info=ip,fullDFT_info=iff)
    out['pass']=bool(out['full_DFT_equals_P_periodmean_relative']<1e-12 and out['iv_parity_max_modal_log_difference']<2e-4 and relres<1e-9 and reljac<1e-7 and out['estimate_scale_modal_log_difference']<2e-3)
    write_json(OUT/'diagnostics'/'scale_parity_gate.json',out);print(json.dumps(out),flush=True)
    return out

def summarize():
    rows=[json.loads(p.read_text(encoding='utf-8')) for p in sorted((OUT/'results').glob('*.json'))]
    groups={}
    for row in rows:groups.setdefault((row['P'],row['method']),[]).append(row)
    summary=[dict(P=p,method=m,**statistics(rs)) for (p,m),rs in sorted(groups.items())]
    write_json(OUT/'summary.json',summary)
    theory=[]
    for P in (2,4,8,16):
        cfs=[np.load(p)['covariance'][c.IW:c.ITR,c.IW:c.ITR] for p in (OUT/'diagnostics').glob(f'P{P:02d}_*_oracle_theory.npz')]
        if cfs:
            variances=np.maximum(np.diagonal(np.array(cfs),axis1=1,axis2=2),0)
            th_delta=100*np.sqrt(np.mean(variances,axis=0))
            # Conditional log-normal transform of percentage exp(error)-1;
            # include the between-design change in the conditional mean.
            mu=np.exp(.5*variances)
            th=100*np.sqrt(np.maximum(np.mean(np.exp(2*variances),axis=0)-np.mean(mu,axis=0)**2,0))
            within=100*np.sqrt(np.mean(np.exp(variances)*np.expm1(variances),axis=0))
            rs=groups.get((P,'M3o'),[]);st=statistics(rs)
            theory.append(dict(P=P,n=len(cfs),theory_percent_sd=th.tolist(),theory_delta_percent_sd=th_delta.tolist(),
                               theory_within_design_exp_percent_sd=within.tolist(),theory_exp_bias_percent=(100*(np.mean(mu,axis=0)-1)).tolist(),
                               theory_assumption='zero-mean Gaussian modal log errors conditional on realised excitation; stationary Fourier covariance; nuisance projected',**st))
    write_json(OUT/'oracle_theory_summary.json',theory)
    print('summaries',len(rows),'records',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--gate',action='store_true');ap.add_argument('--summary',action='store_true');ap.add_argument('--periods',default='2,4,8,16');ap.add_argument('--count',type=int,default=40);ap.add_argument('--workers',type=int,default=3);ap.add_argument('--start',type=int,default=0);args=ap.parse_args()
    if args.gate:gate()
    elif args.summary:summarize()
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            fs=[pool.submit(one,P,r) for P in map(int,args.periods.split(',')) for r in range(args.start,args.count)]
            for fut in as_completed(fs):
                try:print(json.dumps(fut.result()),flush=True)
                except Exception as exc:print(json.dumps({'worker_exception':repr(exc)}),flush=True)
        summarize()
