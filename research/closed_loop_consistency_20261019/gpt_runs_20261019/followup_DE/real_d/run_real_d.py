"""Task D real data: unchanged ZOH model, copied loader, line-only estimators.

Writes only beneath real_d. Original observations and earlier results are read-only.
Run `python run_real_d.py --z 0 --p 2 --start 1` for a pilot.
Run `python run_real_d.py --all --workers 2` for all unique windows.
"""
import os
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import json
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
FOLLOW = HERE.parent
sys.path.insert(0, str(FOLLOW))
import loader_copy as load
import realfit_extended as m
import weighted_core as wc
from sml_core_copy import solve_sml, solve_sml_cov, gmodel_from_output

METHODS = ('M2', 'M3', 'M3w', 'M5', 'M5n')
ZPOSITIONS = ('-0.5', '0', '0.5', '1')


def atomic_json(path, obj):
    path.parent.mkdir(exist_ok=True, parents=True)
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=True), encoding='utf-8')
    tmp.replace(path)


def tag(z, p, start, method):
    return f'z{z}_P{p:02d}_s{start:02d}_{method}'


def initial(ctx, U, Y):
    b = np.zeros(m.NB)
    b[m.IW:m.IZ] = np.log(2*np.pi*np.array([10.6,11.7,16.1]))
    b[m.IZ:m.ITR] = np.log(.03)
    return m.init_lin(ctx, U, Y, np.ones_like(U.real), b)


def get_beta(result):
    return np.asarray(result['beta'], float)


def base_metadata(d, method):
    return dict(z_mm=d['z'], P=d['nper'], start_cycle_zero_based=d['start'], method=method,
                sampling='zoh', band_hz=[3,25], n_parameters=m.NB,
                fixed_delay_correction_s=-m.Ts/2, fitting_selection='excited_lines_only',
                K_arx_relative_max_residual=d['kres'],
                reference_repeat_relative=d['reference_repeat_relative'],
                boundary_identity_relative=d['boundary_identity_relative'],
                n_frequency_bins=int(np.count_nonzero(d['k'] % d['nper'] == 0)),
                n_nonexcited_bins=int(np.count_nonzero(d['k'] % d['nper'] != 0)),
                covariance_nl=16 if method in ('M3w','M5n') else None,
                covariance_outer_updates=1 if method in ('M3w','M5n') else None,
                covariance_shrinkage=.1 if method in ('M3w','M5n') or (method=='M5' and d['nper']<5) else 0.,
                real_truth_available=False)


def normalize_info(beta, info, method):
    """One explicit campaign rule; preserve the solver's original fields."""
    out = dict(info)
    out['solver_info'] = dict(info)
    finite = bool(np.all(np.isfinite(beta)))
    if method in ('M5','M5n'):
        iters = int(info.get('nfev', info.get('iterations', 0)))
        success = bool(info.get('success', False)) and finite
        # The shared SML core exposes 'status' and 'success' consistently;
        # retain status separately rather than reclassifying score as failure.
        if 'success' not in info:
            success = bool(info.get('converged', False)) and finite
        out['iterations'] = iters
        out['iteration_unit'] = 'objective_function_evaluations'
        out['score_normalized'] = info.get('normalized_score', info.get('score_normalized'))
        out['max_evaluations'] = 160
        out['converged'] = success
        out['failure'] = not success
    else:
        iters = int(info.get('iterations', 0))
        original_status = str(info.get('status', ''))
        success = finite and iters < 80 and original_status in ('pole_step_converged','converged')
        if 'converged' in info:
            success = bool(info['converged']) and finite and iters < 80
        out['iterations'] = iters
        out['iteration_unit'] = 'IV_Gauss_Newton_steps'
        out['score_normalized'] = info.get('score_normalized', info.get('score_max', info.get('normalized_score')))
        out['max_iterations'] = 80
        out['pole_step_tolerance'] = 1e-5
        out['pole_step_cap'] = .05
        out['converged'] = success
        out['failure'] = not success
    f, zeta = m.modal(beta)
    out.update(beta=beta.tolist(), finite_endpoint=finite,
               frequency_hz=f.tolist(), damping_percent=zeta.tolist())
    return out


def run_window(z, p, start, methods=METHODS):
    d = load.prepare(z,p,start)
    line = d['k'] % p == 0
    ctx = m.Ctx(d['f'][line], load.TAU-m.Ts/2, 'zoh')
    ctx_all = m.Ctx(d['f'], load.TAU-m.Ts/2, 'zoh')
    U,Y,D = [d['ff'][key][line] for key in ('U','Y','D')]
    model = lambda b,u: m.out_jac(ctx,b,u)
    model_all = lambda b,u: m.out_jac(ctx_all,b,u)
    K = -d['Q'][None] @ d['kv'](ctx.zi)
    Uh = lambda b: np.linalg.solve(np.eye(4)[None]+K@m.resp(ctx,b),D[:,:,None])[:,:,0]
    b0 = initial(ctx,U,Y)
    cache = {}
    cov = None
    cache_dir = HERE/'data'
    cache_dir.mkdir(exist_ok=True)
    npzpath = cache_dir/f'z{z}_P{p:02d}_s{start:02d}.npz'
    if not npzpath.exists():
        np.savez_compressed(npzpath,k=d['k'],f=d['f'],U=d['ff']['U'],Y=d['ff']['Y'],D=d['ff']['D'],
                            Uper=d['per']['U'],Yper=d['per']['Y'],Dper=d['per']['D'],K_lines=K,
                            initial_beta=b0,line_indices=np.flatnonzero(line))
    for method in METHODS:
        if method not in methods and method!='M3':
            continue
        path = HERE/'raw'/(tag(z,p,start,method)+'.json')
        if path.exists():
            cache[method] = json.loads(path.read_text(encoding='utf-8'))
            continue
        info = base_metadata(d,method)
        t0 = time.perf_counter()
        if p==1 and method in ('M3w','M5','M5n'):
            info.update(status='not_applicable',converged=None,failure=False,
                        unavailable=True,reason='One-period half: no non-excited bins and no periodic sample covariance.')
            atomic_json(path,info)
            cache[method]=info
            continue
        try:
            if method in ('M2','M3'):
                b,si = wc.iv_solve(model,U,Y,U if method=='M2' else Uh,np.ones((len(U),4)),b0,
                                   m.IW,m.ITR,maxiter=80,tol=1e-5,cap=.05)
                ib = b0
            elif method=='M3w':
                bm3 = get_beta(cache['M3'])
                cov = wc.local_covariances(d['k'],d['ff']['U'],d['ff']['Y'],bm3,model_all,p,nl=16,alpha=.1)
                C = cov['covariance_residual']
                W = wc.precision(C)
                b,si = wc.iv_solve(model,U,Y,Uh,W,bm3,m.IW,m.ITR,maxiter=80,tol=1e-5,cap=.05)
                ib = bm3
                info['initial_M3_failure']=cache['M3'].get('failure')
                info['weight_condition_median']=float(np.median(np.linalg.cond(C)))
                info['weight_condition_max']=float(np.max(np.linalg.cond(C)))
                np.savez_compressed(HERE/'raw'/(tag(z,p,start,method)+'_cov.npz'),**cov)
            elif method=='M5':
                up,yp = d['per']['U'],d['per']['Y']
                ib = initial(ctx,up.mean(0),yp.mean(0))
                gm = gmodel_from_output(model,m.NB,len(U))
                b,si = solve_sml(up,yp,ib,model,gm,maxiter=160,shrink=.1 if p<5 else 0.)
                info['DFT_units']='mean_of_unnormalized_period_DFTs'
            else:
                bm3 = get_beta(cache['M3'])
                if cov is None:
                    cov = wc.local_covariances(d['k'],d['ff']['U'],d['ff']['Y'],bm3,model_all,p,nl=16,alpha=.1)
                ib = b0
                gm = gmodel_from_output(model,m.NB,len(U))
                b,si = solve_sml_cov(U,Y,cov['CY'],cov['CU'],cov['CYU'],ib,model,gm,maxiter=160)
                info['DFT_units']='unnormalized_full_record_DFTs'
                info['initial_M3_failure']=cache['M3'].get('failure')
            info.update(normalize_info(b,si,method))
            info['initial_beta']=ib.tolist()
            info['status']='converged' if info['converged'] else 'solver_failure'
        except Exception as exc:
            info.update(status='exception',exception=repr(exc),traceback=traceback.format_exc(),
                        converged=False,failure=True,iterations=0)
        info['elapsed_seconds']=time.perf_counter()-t0
        atomic_json(path,info)
        if 'beta' in info:
            np.savez_compressed(path.with_suffix('.npz'),beta=np.asarray(info['beta']),initial_beta=np.asarray(info['initial_beta']))
        cache[method]=info
        print(tag(z,p,start,method),info.get('status'),info.get('iterations'),
              info.get('frequency_hz'),info.get('damping_percent'),flush=True)
    return z,p,start


def window_requests():
    rows=[]
    for z in ZPOSITIONS:
        for p in (2,4,8):
            for role,nper,start in [('full',p,1),('front_half',p//2,1),('back_half',p//2,1+p//2)]:
                rows.append(dict(z=z,parent_P=p,role=role,nper=nper,start=start))
    return rows


def audit_sources():
    files=[load.DATA/f'cl_data_z{z}.mat' for z in ZPOSITIONS]+[load.DELAY]
    original=HERE.parents[1]/'experiment_c'
    files += [original/f for f in ('realfit_original_copy.py','realfit_extended.py','run_c.py')]
    records=[]
    for path in files:
        records.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size))
    atomic_json(HERE/'source_inventory.json',records)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--all',action='store_true');ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--z',default='0');ap.add_argument('--p',type=int,default=2);ap.add_argument('--start',type=int,default=1)
    ap.add_argument('--methods',default=','.join(METHODS))
    args=ap.parse_args()
    (HERE/'raw').mkdir(exist_ok=True);(HERE/'logs').mkdir(exist_ok=True)
    audit_sources()
    atomic_json(HERE/'window_requests.json',window_requests())
    if args.all:
        todo=sorted({(r['z'],r['nper'],r['start']) for r in window_requests()})
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(run_window,z,p,s,args.methods.split(',')) for z,p,s in todo]
            for fut in as_completed(futures):
                print('WINDOW_COMPLETE',fut.result(),flush=True)
    else:
        run_window(args.z,args.p,args.start,args.methods.split(','))


if __name__=='__main__':
    main()
