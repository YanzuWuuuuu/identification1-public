import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
import sys,json,numpy as np
sys.dont_write_bytecode=True
from run_sml import one as sml_one,HERE,atomic_npz
from run_indirect_ivgn import solve,c

for P in [2,4]:
    for r in [0,1,2]:
        path=HERE.parent/'experiment_a'/'data'/f'P{P:02d}_r{r:03d}.npz'
        for shr in [.01,.3]:sml_one(path,shr,suffix=f'_sensitivity{shr}')
out=HERE/'sensitivity_results';out.mkdir(exist_ok=True)
for P in [8,16]:
    for r in [0,1]:
        data=dict(np.load(HERE/'indirect_data'/f'P{P:02d}_r{r:03d}.npz'))
        pert=np.load(HERE.parent/'experiment_a'/'data'/f'P{P:02d}_r{r:03d}.npz')['pert']
        for shr in [.01,.3]:
            b,b0,info,steps=solve(data,pert,shr);bt=c.truth_beta();info.update(P=P,r=r,method='M6_reg_IVGN_sensitivity')
            atomic_npz(out/f'P{P:02d}_r{r:03d}_M6reg_shrink{shr}.npz',beta=b,beta0=b0,truth=bt,
              relative_error_percent=100*np.expm1(b[c.IW:c.ITR]-bt[c.IW:c.ITR]),info=json.dumps(info),steps=steps)
            print(P,r,shr,info,flush=True)
