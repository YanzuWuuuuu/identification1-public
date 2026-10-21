"""Independent DFT/covariance/transient unit check on real z=0,P=2."""
import os
for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[n]='1'
import sys
sys.dont_write_bytecode=True
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import loader_copy as load
import realfit_extended as m
from sml_core_copy import SMLCovResidual,gmodel_from_output

d=load.prepare('0',2,1);p=2;sl=d['k']%p==0
ctx=m.Ctx(d['f'][sl],load.TAU-m.Ts/2,'zoh')
model=lambda b,u:m.out_jac(ctx,b,u)
gm=gmodel_from_output(model,m.NB,len(ctx.f))
cov=np.load(HERE/'raw/z0_P02_s01_M3w_cov.npz')
b=np.asarray(json.loads((HERE/'raw/z0_P02_s01_M2.json').read_text(encoding='utf-8'))['beta'])
U,Y=d['ff']['U'][sl],d['ff']['Y'][sl]
objfull=SMLCovResidual(U,Y,cov['CY'],cov['CU'],cov['CYU'],model,gm)
objmean=SMLCovResidual(U/p,Y/p,cov['CY']/p**2,cov['CU']/p**2,cov['CYU']/p**2,model,gm)
bm=b.copy();bm[m.ITR:]/=p
rf,jf=objfull.calc(b);rm,jm=objmean.calc(bm);jm[:,m.ITR:]/=p
result=dict(P=p,z_mm=0,
            DFT_U_relative=float(np.linalg.norm(U/p-d['per']['U'].mean(0))/np.linalg.norm(U/p)),
            DFT_Y_relative=float(np.linalg.norm(Y/p-d['per']['Y'].mean(0))/np.linalg.norm(Y/p)),
            whitened_residual_relative=float(np.linalg.norm(rf-rm)/np.linalg.norm(rf)),
            transformed_jacobian_relative=float(np.linalg.norm(jf-jm)/np.linalg.norm(jf)),
            mean_covariance_divisor=p**2,mean_transient_divisor=p)
result['pass']=all(result[k]<1e-8 for k in ('DFT_U_relative','DFT_Y_relative','whitened_residual_relative','transformed_jacobian_relative'))
(HERE/'diagnostics').mkdir(exist_ok=True)
(HERE/'diagnostics/real_scale_check.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
