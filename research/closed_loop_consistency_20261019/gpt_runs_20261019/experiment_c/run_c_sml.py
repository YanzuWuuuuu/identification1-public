import os
for n in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']: os.environ.setdefault(n,'1')
import sys, json, time
from pathlib import Path
import numpy as np
import run_c as r
import realfit_extended as m
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'comparators'))
from sml_core import solve_sml,gmodel_from_output,sample_covariances

def run(z,p,shrink):
    d=r.prepare(z,p);c=m.Ctx(d['fp'],r.TAU-m.Ts/2,'zoh')
    up,yp=d['per']['U'],d['per']['Y'];um=up.mean(0);ym=yp.mean(0)
    b=np.zeros(m.NB);b[m.IW:m.IZ]=np.log(2*np.pi*np.array([10.6,11.7,16.1]));b[m.IZ:m.ITR]=np.log(.03)
    b=m.init_lin(c,um,ym,np.ones((len(c.f),4)),b)
    model=lambda bb,uu:m.out_jac(c,bb,uu)
    gmodel=gmodel_from_output(model,m.NB,len(c.f))
    tic=time.time();bb,info=solve_sml(up,yp,b,model,gmodel,maxiter=int(os.environ.get('C_SML_MAX','160')),shrink=shrink)
    f,zeta=m.modal(bb)
    info.update(z=z,P=p,sampling='zoh',method='M5' if shrink==0 else 'M5_reg',start_cycle_zero_based=1,
                frequency_hz=f.tolist(),damping_percent=zeta.tolist(),elapsed_s=time.time()-tic,
                K_arx_relative_max_residual=d['kres'],reference_repeat_relative=d['reference_repeat_relative'],
                n_frequency_bins=len(c.f),n_parameters=m.NB,beta=bb.tolist(),initial_beta=b.tolist(),
                covariance_rank_upper_bound=min(p-1,4),period_average_covariance_divisor=p*(p-1),
                iterations=info.get('nfev'),score_normalized=info.get('normalized_score'),
                converged_pole_step=None)
    return info

if __name__=='__main__':
    z=sys.argv[1]
    for p in [2,4,8]:
        for shrink in ([0.,.1] if p<8 else [0.]):
            meth='M5' if shrink==0 else 'M5_reg'
            if (r.HERE/'raw'/f'z{z}_P{p}_s1_zoh_{meth}.json').exists():continue
            r.save_result(run(z,p,shrink))
