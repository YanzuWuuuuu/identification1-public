"""Meaningful numerical checks for covariance differentiation and exact mean data."""
import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
os.environ['NOFLEX']='1'
import numpy as np,json,pathlib,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parents[1]/'code'))
import cl2 as c
from sml_core import SMLResidual,solve_sml,gmodel_from_output
rng=np.random.default_rng(302);zs=np.exp(2j*np.pi*np.linspace(3,100,180)/5000);bt=c.truth_beta()
model=lambda b,U:c.out_jac(zs,b,U)
gm=gmodel_from_output(model,c.NB,len(zs))
U0=(rng.normal(size=(len(zs),4))+1j*rng.normal(size=(len(zs),4)))*100
Y0=np.einsum('kij,kj->ki',c.resp(zs,bt),U0)
du=(rng.normal(size=(16,len(zs),4))+1j*rng.normal(size=(16,len(zs),4)))*.01;du-=du.mean(0)
dy=(rng.normal(size=du.shape)+1j*rng.normal(size=du.shape))*1e-6;dy-=dy.mean(0)
U=U0[None]+du;Y=Y0[None]+dy
obj=SMLResidual(U,Y,model,gm,.1);btest=bt.copy();btest[c.IZ:c.ITR]+=.03
rr,jj=obj.calc(btest);direc=rng.normal(size=c.NB);direc[:c.NLIN]*=1e-5;direc[c.ITR:]*=1e-8
h=1e-5;fd=(obj.fun(btest+h*direc)-obj.fun(btest-h*direc))/(2*h)
err=np.linalg.norm(fd-jj@direc)/np.linalg.norm(fd)
assert err<1e-4,err
b0=bt.copy();b0[c.IW:c.IZ]+=.001;b0[c.IZ:c.ITR]+=.05;b0[:c.NLIN]=0
b0=c.init_lin(zs,U.mean(0),Y.mean(0),np.ones((len(zs),4)),b0)
b,info=solve_sml(U,Y,b0,model,gm,maxiter=100,shrink=0.)
rel=np.linalg.norm(c.resp(zs,b)-c.resp(zs,bt))/np.linalg.norm(c.resp(zs,bt))
pole=np.max(np.abs(np.expm1(b[c.IW:c.ITR]-bt[c.IW:c.ITR])))
assert rel<1e-9,(rel,pole,info)
out=dict(gradient_direction_relative_error=float(err),zero_noise_mean_response_relative_error=float(rel),
         zero_noise_mean_pole_relative_error=float(pole),solver=info,
         note='Mean U,Y exactly noiseless; zero-mean replicate deviations supply nonsingular covariance. Raw zero covariance has no inverse.')
(HERE/'validation_sml.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
