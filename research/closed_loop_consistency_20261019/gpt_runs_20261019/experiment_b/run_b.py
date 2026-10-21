"""Experiment B: unchanged cl2 model, output-sparse algebra for its IV-GN equations."""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
ENV = dict(NOFLEX='0', FC='60', AL='30', AH='100', FHI='800', FLO='3',
           SW='0.03', FSPR='0.01', ZSPR='0.4', NN='300000', NPER='10')
os.environ.update(ENV)
import sys, json, time, traceback, hashlib
sys.dont_write_bytecode = True
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
HERE = Path(__file__).resolve().parent
ORIGINAL = HERE.parent.parent / 'code'
sys.path.insert(0, str(ORIGINAL))
import cl2 as c

METHODS = ('dir_band_W1', 'dir_lines_W1', 'ref_band_W1')
Q = 2*c.NM
L = 10*c.NM
IDS = [np.array([p for j in range(c.NM) for p in
       list(range(32*j+4*o,32*j+4*o+4)) + list(range(32*j+16+4*o,32*j+16+4*o+4)) +
       [c.ITR+8*j+2*o,c.ITR+8*j+2*o+1]]) for o in range(4)]
NID = np.arange(c.IW,c.ITR)

def atomic_json(path, obj):
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=True), encoding='utf-8')
    os.replace(temp,path)

def atomic_npz(path, **arr):
    tmp = path.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, **arr)
    os.replace(tmp,path)

def compact(z, beta, U, nonlinear=True):
    """Same Y, linear and nonlinear derivatives as cl2.out_jac, without zero columns."""
    lin,lw,lz,tr = c.unpack(beta)
    h0,h1,ad = c.blocks(z,lw,lz)
    zi=1/z; n=len(z)
    LF=np.empty((n,L),complex)
    N=np.zeros((n,4,Q),complex)
    Y=np.zeros((n,4),complex)
    for j in range(c.NM):
        den=ad[j][0]+ad[j][1]*zi+ad[j][2]*zi**2
        l0=U@lin[j,0].T; l1=U@lin[j,1].T
        tn=tr[j,:,0][None]+tr[j,:,1][None]*zi[:,None]
        Y += h0[j][:,None]*l0+h1[j][:,None]*l1+tn/den[:,None]
        off=10*j
        LF[:,off:off+4]=h0[j][:,None]*U
        LF[:,off+4:off+8]=h1[j][:,None]*U
        LF[:,off+8]=1/den; LF[:,off+9]=zi/den
        if nonlinear:
            for which,idx in ((0,j),(1,c.NM+j)):
                a0,a1,a2=c.zoh_2nd(np.exp(lw[j]+(1e-6 if which==0 else 0)),
                                  np.exp(lz[j]+(1e-6 if which==1 else 0)),z)
                den2=a2[0]+a2[1]*zi+a2[2]*zi**2
                N[:,:,idx]=((a0-h0[j])[:,None]*l0+(a1-h1[j])[:,None]*l1+
                            tn/den2[:,None]-tn/den[:,None])/1e-6
    return Y,LF,N

def equations(z,U,Y,Uh,beta,nonlinear=True,chunk=3000):
    A=[np.zeros((L,L)) for o in range(4)]
    B=[np.zeros((L,Q)) for o in range(4)]
    C=[np.zeros((Q,L)) for o in range(4)]
    D=np.zeros((Q,Q)); gl=[np.zeros(L) for o in range(4)]; gn=np.zeros(Q)
    zzl=[np.zeros(L) for o in range(4)]; zzn=np.zeros(Q); errsq=0.
    for lo in range(0,len(z),chunk):
        sl=slice(lo,min(lo+chunk,len(z)))
        yh,jl,jn=compact(z[sl],beta,U[sl],nonlinear)
        if Uh is U: zl,zn=jl,jn
        else: _,zl,zn=compact(z[sl],beta,Uh[sl],nonlinear)
        err=Y[sl]-yh; errsq+=np.sum(abs(err)**2)
        aa=np.real(zl.conj().T@jl)
        for o in range(4):
            A[o]+=aa
            gl[o]+=np.real(zl.conj().T@err[:,o])
            zzl[o]+=np.sum(abs(zl)**2,axis=0)
            if nonlinear:
                B[o]+=np.real(zl.conj().T@jn[:,o,:])
                C[o]+=np.real(zn[:,o,:].conj().T@jl)
                D+=np.real(zn[:,o,:].conj().T@jn[:,o,:])
                gn+=np.real(zn[:,o,:].conj().T@err[:,o])
                zzn+=np.sum(abs(zn[:,o,:])**2,axis=0)
    return A,B,C,D,gl,gn,zzl,zzn,errsq

def getstep(eq,nonlinear=True):
    A,B,C,D,gl,gn,*_=eq
    step=np.zeros(c.NB)
    if nonlinear:
        sn=1/np.sqrt(abs(np.diag(D))+1e-300)
        ds=D*sn[:,None]*sn[None,:]; gs=gn*sn
    sols=[]; scl=[]
    for o in range(4):
        sc=1/np.sqrt(abs(np.diag(A[o]))+1e-300)
        ao=A[o]*sc[:,None]*sc[None,:]
        if nonlinear:
            bo=B[o]*sc[:,None]*sn[None,:]
            co=C[o]*sn[:,None]*sc[None,:]
            sol=np.linalg.solve(ao,np.column_stack((gl[o]*sc,bo)))
            ds-=co@sol[:,1:]; gs-=co@sol[:,0]
        else:
            sol=np.linalg.solve(ao,gl[o]*sc)
        sols.append(sol); scl.append(sc)
    if nonlinear:
        dn=np.linalg.solve(ds,gs); step[NID]=sn*dn
    for o in range(4):
        step[IDS[o]]=scl[o]*(sols[o][:,0]-sols[o][:,1:]@dn) if nonlinear else scl[o]*sols[o]
    return step

def score(eq):
    *_,gl,gn,zzl,zzn,errsq=eq
    g=np.concatenate(gl+[gn]); zz=np.concatenate(zzl+[zzn])
    return float(np.max(abs(g)/np.sqrt(np.maximum(zz*errsq,1e-300))))

def solve(z,U,Y,Uh,b0,maxit=80,tol=1e-5):
    beta=b0.copy(); hist=[]; status='iteration_limit'
    for it in range(maxit):
        t0=time.perf_counter(); ub=Uh(beta) if callable(Uh) else Uh
        eq=equations(z,U,Y,ub,beta)
        step=getstep(eq); raw=float(np.max(abs(step[NID])))
        if raw>.05: step*=.05/raw
        used=float(np.max(abs(step[NID])))
        beta+=step
        hist.append([it+1,raw,used,score(eq),eq[-1],time.perf_counter()-t0])
        if not np.all(np.isfinite(beta)): status='nonfinite'; break
        if used<tol: status='pole_step_converged'; break
    # Independent full moment diagnostic at returned parameter, not just pole stopping rule.
    ub=Uh(beta) if callable(Uh) else Uh
    eq=equations(z,U,Y,ub,beta)
    return beta,status,np.array(hist),score(eq),float(eq[-1])

def run(rep):
    (HERE/'raw').mkdir(exist_ok=True)
    def final_path(m):
        base=HERE/'raw'/f'rep{rep:03d}_{m}.npz'
        if base.exists():
            with np.load(base) as a:
                meta=json.loads(str(a['metadata_json']))
            if meta['status']=='exception':
                return HERE/'raw'/f'rep{rep:03d}_{m}_retry1.npz'
        return base
    if all(final_path(m).exists() for m in METHODS): return dict(rep=rep,cached=True)
    seed=52061019+rep
    rng=np.random.default_rng(seed); start=time.perf_counter()
    N=300000; P=10
    d,u,y=c.simulate_periodic(30000,P,.2,3e-7,rng,settle=2)
    simsec=time.perf_counter()-start
    k=np.arange(int(np.ceil(3*N/c.fs)),int(np.floor(800*N/c.fs))+1)
    z=np.exp(2j*np.pi*k/N)
    D=np.fft.fft(d,axis=0)[k]; U=np.fft.fft(u,axis=0)[k]; Y=np.fft.fft(y,axis=0)[k]
    del d,u,y
    kt=c.Ktf(z); bt=c.truth_beta(); pert=rng.standard_normal(Q)
    initial=bt.copy(); initial[c.IW:c.IW+c.NM]+=np.log(1+.01*pert[:c.NM])
    initial[c.IZ:c.ITR]+=.4*pert[c.NM:]
    initial[:c.NLIN]=0.; initial[c.ITR:]=0.
    results=[]
    for method in METHODS:
        outfile=final_path(method)
        if outfile.exists(): continue
        t0=time.perf_counter(); sel=(k%P==0) if '_lines_' in method else np.ones(len(k),bool)
        zs=z[sel]; us=U[sel]; ys=Y[sel]; ds=D[sel]; ks=kt[sel]
        err=''; beta=np.full(c.NB,np.nan); b0=initial.copy(); status='exception'; hist=np.empty((0,6)); sc=np.nan; obj=np.nan
        try:
            b0+=getstep(equations(zs,us,ys,us,b0,nonlinear=False),nonlinear=False)
            if method.startswith('ref'):
                def instr(b):
                    return np.linalg.solve(np.eye(4)[None]+ks@c.resp(zs,b),ds[...,None])[...,0]
                uh=instr
            else: uh=us
            beta,status,hist,sc,obj=solve(zs,us,ys,uh,b0)
        except Exception:
            err=traceback.format_exc()
        elapsed=time.perf_counter()-t0
        errors=100*np.expm1(beta[c.IW:c.ITR]-bt[c.IW:c.ITR])
        meta=dict(rep=rep,seed=seed,method=method,status=status,iterations=len(hist),
                  duration_seconds=elapsed,simulation_seconds=simsec,nbeta=c.NB,nfreq=int(sel.sum()),
                  score_relative=sc,objective=obj,error=err,pole_step_tolerance=1e-5,max_iterations=80)
        atomic_npz(outfile,beta=beta,beta0=b0,truth_beta=bt,error_percent=errors,history=hist,
                   perturbation=pert,metadata_json=np.array(json.dumps(meta)))
        atomic_json(outfile.with_suffix('.json'),meta)
        print(json.dumps(meta),flush=True); results.append(meta)
    return dict(rep=rep,elapsed=time.perf_counter()-start,results=results)

def validate():
    rng=np.random.default_rng(83); z=np.exp(2j*np.pi*np.geomspace(3,800,1000)/c.fs)
    beta=c.truth_beta(); beta[c.IW:c.IW+c.NM]+=.001
    u=rng.normal(size=(len(z),4))+1j*rng.normal(size=(len(z),4))
    uh=u*.93+(.1+.1j)*rng.normal(size=u.shape)
    y=c.out_jac(z,c.truth_beta(),u,need_jac=False)[0]
    yh,j=c.out_jac(z,beta,u); zz=c.out_jac(z,beta,uh)[1]
    yy,ll,nn=compact(z,beta,u)
    jc=np.zeros_like(j)
    for o in range(4): jc[:,o,IDS[o]]=ll
    jc[:,:,NID]=nn
    eq=equations(z,u,y,uh,beta)
    A=np.real(zz.reshape(-1,c.NB).conj().T@j.reshape(-1,c.NB))
    g=np.real(zz.reshape(-1,c.NB).conj().T@(y-yh).ravel())
    Ac=np.zeros_like(A); gc=np.zeros_like(g)
    for o in range(4):
        Ac[np.ix_(IDS[o],IDS[o])]=eq[0][o]
        Ac[np.ix_(IDS[o],NID)]=eq[1][o]; Ac[np.ix_(NID,IDS[o])]=eq[2][o]
        gc[IDS[o]]=eq[4][o]
    Ac[np.ix_(NID,NID)]=eq[3]; gc[NID]=eq[5]
    scale=1/np.sqrt(abs(np.diag(A))+1e-300)
    d0=scale*np.linalg.solve(A*scale[:,None]*scale[None,:],g*scale)
    d1=getstep(eq)
    checks=dict(nm=c.NM,nb=c.NB,closed_loop_spectral_radius=float(c.closed_loop()[-1]),
                output_relerr=float(np.linalg.norm(yy-yh)/np.linalg.norm(yh)),
                jacobian_relerr=float(np.linalg.norm(jc-j)/np.linalg.norm(j)),
                scaled_normal_matrix_relerr=float(np.linalg.norm((Ac-A)*scale[:,None]*scale[None,:])/np.linalg.norm(A*scale[:,None]*scale[None,:])),
                scaled_gradient_relerr=float(np.linalg.norm((gc-g)*scale)/np.linalg.norm(g*scale)),
                pole_step_relerr=float(np.linalg.norm(d1[NID]-d0[NID])/np.linalg.norm(d0[NID])),
                original_sha256=hashlib.sha256((ORIGINAL/'cl2.py').read_bytes()).hexdigest(),environment=ENV)
    atomic_json(HERE/'validation.json',checks); print(json.dumps(checks,indent=2),flush=True)
    assert checks['closed_loop_spectral_radius']<1
    assert checks['jacobian_relerr']<1e-12
    assert checks['scaled_normal_matrix_relerr']<1e-10
    assert checks['pole_step_relerr']<1e-4

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(); p.add_argument('--validate',action='store_true');p.add_argument('--reps',type=int,default=20);p.add_argument('--workers',type=int,default=2);p.add_argument('--start',type=int,default=0)
    a=p.parse_args()
    if a.validate: validate()
    else:
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            fut=[pool.submit(run,r) for r in range(a.start,a.reps)]
            for f in as_completed(fut): print('REP_DONE '+json.dumps(f.result()),flush=True)
