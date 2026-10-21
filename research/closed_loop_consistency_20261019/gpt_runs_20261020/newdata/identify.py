import os
for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[n]='1'
import sys
sys.dont_write_bytecode=True
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import json,time,argparse,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'baseline'))
import realfit_extended as m
import weighted_core as wc
from sml_core_copy import solve_sml_cov,gmodel_from_output
from prepare_data import evalK,js,FS,P0
TAU=np.load(HERE/'baseline/tau_seconds.npy')

def prepare_spectra(rid):
    p=HERE/'spectra'/f'{rid}.npz'
    if p.exists():return np.load(p)
    a=np.load(HERE/'prepared'/f'{rid}.npz');P=int(a['P']);N=P*P0
    k=np.arange(1,N//2);f=k*FS/N;sl=(f>=3)&(f<=25);k=k[sl];f=f[sl]
    kp=np.arange(18,151);fp=kp*FS/P0
    data=dict(k=k,f=f,fp=fp,P=P,N=N,line_indices=np.flatnonzero(k%P==0))
    for x in ('U','Y','D'):
        data[x]=np.fft.rfft(a[x],axis=0)[k]
        data[x+'per']=np.fft.rfft(a[x].reshape(P,P0,4),axis=1)[:,kp]
    meta=json.loads((HERE/'prepared'/f'{rid}.json').read_text(encoding='utf-8'))
    c=np.load(HERE/'controllers'/f"{meta['controller']}.npz")
    data['K']=evalK(f,c['A'],c['B'],c['T']);data['K_lines']=data['K'][data['line_indices']]
    np.savez_compressed(p,**data);return np.load(p)

def initial(c,U,Y):
    b=np.zeros(m.NB);b[m.IW:m.IZ]=np.log(2*np.pi*np.array([10.6,11.7,16.1]));b[m.IZ:m.ITR]=np.log(.03)
    return m.init_lin(c,U,Y,np.ones_like(U.real),b)

def default_methods(meta):
    if meta['group']=='D':return ['M2','M3','M3w']
    if meta['group']=='E':return ['M1','M2','M3','M3w']
    if meta['group']=='A' or meta['signs']==[1,1,1,1]:return ['M1','M2','M3','M3w','M5n']
    return ['M1','M2','M3','M3w']

def run(rid,methods=None):
    meta=json.loads((HERE/'prepared'/f'{rid}.json').read_text(encoding='utf-8'))
    req=methods or default_methods(meta);a=prepare_spectra(rid);li=a['line_indices']
    ca=m.Ctx(a['f'],TAU,'zoh');cl=m.Ctx(a['f'][li],TAU,'zoh')
    U,Y,D=(a[k][li] for k in ('U','Y','D'));model=lambda b,u:m.out_jac(cl,b,u);modela=lambda b,u:m.out_jac(ca,b,u)
    K=a['K_lines'];Uh=lambda b:np.linalg.solve(np.eye(4)[None]+K@m.resp(cl,b),D[:,:,None])[:,:,0]
    cache={};cv=None
    for method in ['M1','M2','M3','M3w','M5n']:
        if method not in req and not (method=='M3' and ('M3w' in req or 'M5n' in req)):continue
        path=HERE/'fits'/f'{rid}_{method}.json'
        if path.exists():cache[method]=json.loads(path.read_text(encoding='utf-8'));continue
        t=time.perf_counter();info=dict(id=rid,method=method,group=meta['group'],controller=meta['controller'],z_mm=meta['z_mm'],
            gain=meta['disturb_gain'],P=int(a['P']),N=int(a['N']),band_hz=[3,25],model_parameters=m.NB,
            sampling='zoh',delay='old_effective_delay_minus_100us',reference_lines_only=method!='M1',max_iterations=80,tolerance=1e-5,cap=.05)
        try:
            b0=initial(cl,U,Y)
            if method=='M1':
                ib=initial(ca,a['U'],a['Y']);b,si=wc.iv_solve(modela,a['U'],a['Y'],a['U'],np.ones_like(a['U'].real),ib,m.IW,m.ITR,maxiter=80,tol=1e-5,cap=.05)
            elif method in ('M2','M3'):
                ib=b0;b,si=wc.iv_solve(model,U,Y,U if method=='M2' else Uh,np.ones_like(U.real),ib,m.IW,m.ITR,maxiter=80,tol=1e-5,cap=.05)
            else:
                bm=np.asarray(cache['M3']['beta']);info['initial_M3_failure']=cache['M3']['failure']
                if cv is None:
                    cv=wc.local_covariances(a['k'],a['U'],a['Y'],bm,modela,int(a['P']),nl=16,alpha=.1)
                    np.savez_compressed(HERE/'fits'/f'{rid}_local_covariance.npz',**cv)
                if method=='M3w':
                    ib=bm;b,si=wc.iv_solve(model,U,Y,Uh,wc.precision(cv['covariance_residual']),ib,m.IW,m.ITR,maxiter=80,tol=1e-5,cap=.05)
                else:
                    ib=b0;gm=gmodel_from_output(model,m.NB,len(U))
                    b,si=solve_sml_cov(U,Y,cv['CY'],cv['CU'],cv['CYU'],ib,model,gm,maxiter=160)
            finite=bool(np.isfinite(b).all());conv=bool(si.get('success',False)) and finite
            if method!='M5n':conv=conv and si.get('iterations',80)<80
            score=si.get('score_max',si.get('normalized_score'));f,z=m.modal(b)
            info.update(beta=b.tolist(),initial_beta=ib.tolist(),solver_info=si,converged=conv,failure=not conv,
                iterations=si.get('iterations',si.get('nfev')),iteration_unit='function_evaluations' if method=='M5n' else 'IV_steps',
                score_normalized=score,strict_score_pass=bool(score is not None and np.isfinite(score) and score<1e-6),
                frequency_hz=f.tolist(),damping_percent=z.tolist(),modal_order=np.argsort(f).tolist(),finite_endpoint=finite)
            if finite:
                yy,_=model(b,U);info['line_raw_output_residual_percent']=float(100*np.linalg.norm(Y-yy)/np.linalg.norm(Y))
                G=m.resp(cl,b);np.savez_compressed(HERE/'fits'/f'{rid}_{method}.npz',beta=b,f=cl.f,G=G)
        except Exception as e:
            info.update(converged=False,failure=True,status='exception',exception=repr(e),traceback=traceback.format_exc())
        info['seconds']=time.perf_counter()-t;js(path,info);cache[method]=info
        print(rid,method,info.get('converged'),info.get('iterations'),np.round(info.get('frequency_hz',[]),4),np.round(info.get('damping_percent',[]),4),round(info['seconds'],1),flush=True)
    return rid

def main():
    p=argparse.ArgumentParser();p.add_argument('--id');p.add_argument('--all',action='store_true');p.add_argument('--workers',type=int,default=3);p.add_argument('--methods')
    args=p.parse_args();methods=args.methods.split(',') if args.methods else None
    if args.all:
        inv=json.loads((HERE/'data_inventory.json').read_text(encoding='utf-8'))
        rows=[r for r in inv if r['group']!='B' and not r.get('duplicate_of')]
        rows.sort(key=lambda r:(0 if r['group']=='A' else 1 if r['group']=='C' and r['signs']==[1]*4 else 2,r['id']))
        js(HERE/'fit_requests.json',[dict(id=r['id'],methods=methods or default_methods(r)) for r in rows])
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            fut=[pool.submit(run,r['id'],methods) for r in rows]
            for f in as_completed(fut):print('RECORD_COMPLETE',f.result(),flush=True)
    else:run(args.id or 'A_nominal',methods)
if __name__=='__main__':main()
