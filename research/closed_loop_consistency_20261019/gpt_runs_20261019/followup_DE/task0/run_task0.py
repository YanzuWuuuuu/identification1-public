"""Read-only Task 0 audit, long-period M1 MC, and direct population equation.

All new artifacts are confined to this directory. The original A solver and
simulation function are reused, including the same stopping/failure rule.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
os.environ.update(NOFLEX='1', FC='60', AL='30', AH='100', FHI='100', FLO='3', SW='0.03', FSPR='0.01', ZSPR='0.4')
import sys, json, time, argparse
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[2]
OLD=BASE/'gpt_runs_20261019'/'experiment_a'
sys.path.insert(0,str(BASE/'code'))
import numpy as np
import scipy.linalg as sla
import cl2 as c
from concurrent.futures import ProcessPoolExecutor, as_completed

def atomic(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=True),encoding='utf-8')
    tmp.replace(path)

def score(z,U,Y,b):
    g=np.zeros(c.NB); zn=np.zeros(c.NB); en=0.
    for start in range(0,len(z),1200):
        sl=slice(start,min(start+1200,len(z)))
        yh,J=c.out_jac(z[sl],b,U[sl]); e=Y[sl]-yh
        g+=np.real(np.einsum('kop,ko->p',J.conj(),e))
        zn+=np.sum(np.abs(J)**2,axis=(0,1)); en+=np.sum(np.abs(e)**2)
    return float(np.max(np.abs(g)/np.maximum(np.sqrt(zn*en),1e-300))),float(np.sqrt(en/np.sum(np.abs(Y)**2)))

def one_long(r):
    dest=HERE/'long_results'/f'P01_P0_240000_r{r:03d}_M1.json'
    if dest.exists(): return json.loads(dest.read_text(encoding='utf-8'))
    t0=time.perf_counter(); P0=240000; P=1; N=P0*P
    seed=202610190+1000+r; rng=np.random.default_rng(seed)
    k=np.arange(int(np.ceil(3*N/c.fs)),int(np.floor(100*N/c.fs))+1)
    z=np.exp(2j*np.pi*k/N); f=k*c.fs/N
    datadir=HERE/'long_data'; datadir.mkdir(exist_ok=True)
    datafile=datadir/f'P01_P0_240000_r{r:03d}.npz'
    if datafile.exists():
        a=np.load(datafile); D,U,Y,pert=[a[x] for x in ('D','U','Y','pert')]
        timing={}
    else:
        td=time.perf_counter()
        d,u,y=c.simulate_periodic(P0,P,.2,3e-7,rng,settle=2)
        D=np.fft.fft(d,axis=0)[k]; U=np.fft.fft(u,axis=0)[k]; Y=np.fft.fft(y,axis=0)[k]
        pert=rng.standard_normal(2*c.NM)
        timing={'data_seconds':time.perf_counter()-td,'d_rms':np.sqrt(np.mean(d*d,axis=0)).tolist()}
        tmp=datafile.with_suffix('.tmp.npz')
        np.savez_compressed(tmp,D=D,U=U,Y=Y,k=k,f=f,zk=z,pert=pert,seed=seed,P=P,P0=P0)
        tmp.replace(datafile)
    bt=c.truth_beta()
    out={'P':P,'P0':P0,'N':N,'replicate':r,'seed':seed,'method':'M1_longperiod','status':'exception','iterations':None,'beta':None,'error_percent':None,**timing}
    try:
        W=np.ones_like(U.real)
        b0=bt.copy(); b0[c.IW:c.IZ]+=np.log(1+.01*pert[:c.NM]); b0[c.IZ:c.ITR]+=.4*pert[c.NM:]
        b0[:c.NLIN]=0; b0[c.ITR:]=0
        b0=c.init_lin(z,U,Y,W,b0)
        b,it=c.iv_solve(z,U,Y,U,W,b0,it=80,tol=1e-5)
        finite=bool(np.all(np.isfinite(b)))
        sc,rr=score(z,U,Y,b) if finite else (np.nan,np.nan)
        out.update(status='pole_step_converged' if finite and it<79 else ('iteration_limit' if finite else 'nonfinite'),iterations=int(it+1),beta=b.tolist(),initial_beta=b0.tolist(),log_error=(b[c.IW:c.ITR]-bt[c.IW:c.ITR]).tolist(),error_percent=(100*np.expm1(b[c.IW:c.ITR]-bt[c.IW:c.ITR])).tolist(),score_max=sc,residual=rr)
    except Exception as exc:
        out['exception']=repr(exc)
    out['elapsed_seconds']=time.perf_counter()-t0
    atomic(dest,out)
    return {k:out[k] for k in ('replicate','status','iterations','error_percent','elapsed_seconds')}

def gl_tensor(z,b,analytic=False):
    Gl=np.zeros((len(z),4,4,c.ITR),complex)
    for mi in range(4):
        U=np.zeros((len(z),4)); U[:,mi]=1
        Gl[:,:,mi,:]=c.out_jac(z,b,U,nonlin=not analytic)[1][:,:,:c.ITR]
    if analytic:
        lin,lw,lz,_=c.unpack(b)
        for j in range(c.NM):
            w=np.exp(lw[j]); ze=np.exp(lz[j])
            A=np.array([[0.,1.,0.],[-w*w,-2*ze*w,1.],[0.,0.,0.]])*c.Ts
            E=sla.expm(A); Ad=E[:2,:2]; Bd=E[:2,2:3]
            R=z[:,None,None]*np.eye(2)[None]-Ad[None]
            x=np.linalg.solve(R,np.broadcast_to(Bd,(len(z),2,1)))
            for q,idx in ((0,c.IW+j),(1,c.IZ+j)):
                dA=np.zeros((3,3)); dA[1,1]=-2*ze*w*c.Ts
                if q==0: dA[1,0]=-2*w*w*c.Ts
                dE=sla.expm_frechet(A,dA,compute_expm=False)
                dx=np.linalg.solve(R,dE[None,:2,:2]@x+dE[None,:2,2:3])
                dh0=w*w*dx[:,0,0]+(2*w*w*x[:,0,0] if q==0 else 0)
                dh1=w*dx[:,1,0]+(w*x[:,1,0] if q==0 else 0)
                Gl[:,:,:,idx]=dh0[:,None,None]*lin[j,0][None]+dh1[:,None,None]*lin[j,1][None]
    return Gl

def pop_limit(P0,P,tag,maxiter=300,analytic=False):
    """Normalized population moments. DFT line power/N is P sd^2 shape.

    Conditional finite-sample boundary/transient correlations are absent from
    this spectral large-record approximation, matching original exact_limit.
    Unlike the original it explicitly includes the periodic line mask.
    """
    dest=HERE/'population'/f'{tag}.json'
    if dest.exists(): return json.loads(dest.read_text(encoding='utf-8'))
    t0=time.perf_counter(); N=P0*P
    k=np.arange(int(np.ceil(3*N/c.fs)),int(np.floor(100*N/c.fs))+1)
    z=np.exp(2j*np.pi*k/N); bt=c.truth_beta(); G0=c.resp(z,bt); K=c.Ktf(z)
    S=np.linalg.inv(np.eye(4)[None]+K@G0); F=S@K
    mask=np.where(k%P==0,P,0.)
    Pd=(.2**2*c.dshape(z)*mask)[:,None,None]
    Pv=(3e-7)**2*np.eye(4)[None]+.03**2*G0@G0.conj().transpose(0,2,1)
    Phiu=Pd*S@S.conj().transpose(0,2,1)+F@Pv@F.conj().transpose(0,2,1)
    EUV=-F@Pv
    e,V=np.linalg.eigh(Phiu); L=V*np.sqrt(np.clip(e,0,None))[:,None,:]
    beta=bt.copy(); history=[]
    for it in range(maxiter):
        M=np.zeros((c.ITR,c.ITR)); h=np.zeros(c.ITR)
        for ss in range(0,len(z),600):
            sl=slice(ss,min(ss+600,len(z)))
            Gl=gl_tensor(z[sl],beta,analytic=analytic); dG=G0[sl]-c.resp(z[sl],beta)
            X=dG@Phiu[sl]+EUV[sl].conj().transpose(0,2,1)
            h+=np.real(np.einsum('kijp,kij->p',Gl.conj(),X))
            Q=np.einsum('kijp,kjl->kilp',Gl,L[sl]).reshape(-1,c.ITR)
            M+=np.real(Q.conj().T@Q)
        sc=1/np.sqrt(np.diag(M)); st=sc*np.linalg.solve(M*sc[:,None]*sc[None,:],h*sc)
        mx=float(np.max(np.abs(st[c.IW:c.ITR])))
        if mx>.2: st*=.2/mx
        beta[:c.ITR]+=st
        history.append({'iteration':it+1,'pole_step':float(np.max(np.abs(st[c.IW:c.ITR]))),'scaled_moment_max':float(np.max(np.abs(h*sc)))})
        if np.max(np.abs(st[c.IW:c.ITR]))<1e-10: break
    out={'tag':tag,'P0':P0,'P':P,'N':N,'n_bins':len(k),'n_excited':int(np.count_nonzero(mask)),'normalized_line_power_factor':P,'status':'converged' if it<maxiter-1 else 'iteration_limit','iterations':it+1,'analytic_frechet_derivative':analytic,'beta':beta.tolist(),'error_percent':(100*np.expm1(beta[c.IW:c.ITR]-bt[c.IW:c.ITR])).tolist(),'history':history,'elapsed_seconds':time.perf_counter()-t0,'approximation':'Asymptotic spectral equation, no finite-window cross-bin/boundary-nuisance corrections.'}
    atomic(dest,out); return out

def audit():
    rows=[]
    for P in (1,2,4,8,16):
        a=np.load(OLD/'data'/f'P{P:02d}_r000.npz')
        k=a['k']; D=a['D']; sel=k%P==0
        expected=30000*P*P*.2**2*c.dshape(a['zk'][sel])
        power_ratio=np.abs(D[sel])**2/expected[:,None]
        rows.append({'P':P,'N':30000*P,'df_Hz':5000/(30000*P),'band':[float(a['f'][0]),float(a['f'][-1])],'n_bins':len(k),'n_excited':int(sel.sum()),'n_nonexcited':int((~sel).sum()),'full_DFT_vs_P_period_mean_maxabs':float(np.max(np.abs(D[sel]-P*a['Dper'].mean(0)))),'line_power_ratio_maxerror':float(np.max(np.abs(power_ratio-1))),'nonexcited_D_max':float(np.max(np.abs(D[~sel]))) if (~sel).any() else None,'settle_cycles':2,'settle_seconds':12,'sd':.2,'sv':3e-7,'input_noise_sd':.03,'transient':'(tau0+tau1*z^-1)/(1+a1*z^-1+a2*z^-2), 24 fitted nuisance parameters','init_rule':'1% frequency / log damping SD .4 perturbation followed by direct linear fit','seed':int(a['seed'])})
    out={'code_used':str(BASE/'code'/'cl2.py'),'rows':rows,'M1_equals_M2_at_P1':True,'same_noise_same_model_same_settle':True,'original_exact_limit_omission':'original exact_limit.py does not call periodic LINEMASK; this matters for P>1, not P=1','interpretation':'O(1) is an asymptotic order, not equality of finite-record estimator means at all P.'}
    atomic(HERE/'audit.json',out); print(json.dumps(out,ensure_ascii=False),flush=True)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('mode',choices=['audit','long','population','population_analytic']); ap.add_argument('--count',type=int,default=40); ap.add_argument('--workers',type=int,default=2); args=ap.parse_args()
    if args.mode=='audit': audit()
    elif args.mode.startswith('population'):
        for P0,P,tag in [(30000,1,'P1_short'),(240000,1,'P1_long'),(30000,8,'P8')]:
            analytic=args.mode=='population_analytic'
            out=pop_limit(P0,P,tag+('_analytic' if analytic else ''),analytic=analytic); print(json.dumps({k:out[k] for k in ('tag','status','iterations','error_percent','elapsed_seconds')}),flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for fu in as_completed([pool.submit(one_long,r) for r in range(args.count)]): print(json.dumps(fu.result()),flush=True)

if __name__=='__main__': main()
