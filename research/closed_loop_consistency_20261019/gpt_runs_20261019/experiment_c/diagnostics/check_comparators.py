"""Read-only independent algebra and numerical checks of comparator helpers."""
import os,sys
for n in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[n]='1'
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np,json
import scipy.signal as ss
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'comparators'))
import run_indirect as ri
from sml_core import sample_covariances,SMLResidual
rng=np.random.default_rng(17293);K=17;ny=nu=4;nb=3;P=12
cplx=lambda s:rng.normal(size=s)+1j*rng.normal(size=s)
A=cplx((K,4,4));B=cplx((K,4,4,nb))*.2
C=cplx((K,4))*.1;Uper=cplx((P,K,4));Yper=cplx((P,K,4))
def gm(b):
    h=b+.2*b*b;G=A+np.einsum('kijm,m->kij',B,h)
    Gj=B*(1+.4*b)[None,None,None,:]
    return G,Gj
def model(b,U):
    G,Gj=gm(b);t=b[0]*b[1]*C
    y=np.einsum('kij,kj->ki',G,U)+t
    J=np.einsum('kijm,kj->kim',Gj,U);J[:,:,0]+=b[1]*C;J[:,:,1]+=b[0]*C
    return y,J
b=np.array([.2,-.4,.1]);obj=SMLResidual(Uper,Yper,model,gm)
rr,jj=obj.calc(b);err=[]
for j in range(nb):
    bp=b.copy();bm=b.copy();bp[j]+=1e-6;bm[j]-=1e-6
    fd=(obj.fun(bp)-obj.fun(bm))/(2e-6)
    err.append(float(np.linalg.norm(fd-jj[:,j])/np.linalg.norm(fd)))
cy,cu,cyu=sample_covariances(Uper,Yper)
G,_=gm(b);GH=G.conj().transpose(0,2,1)
ce=cy-cyu@GH-G@cyu.conj().transpose(0,2,1)+G@cu@GH
ep=Yper-np.einsum('kij,pkj->pki',G,Uper);ep-=ep.mean(0)
ce2=np.einsum('pki,pkj->kij',ep,ep.conj())/(P*(P-1))
covdiff=float(np.linalg.norm(ce-ce2)/np.linalg.norm(ce2))
costdiff=float(abs(rr@rr-np.einsum('ki,kij,kj->', (Yper.mean(0)-model(b,Uper.mean(0))[0]).conj(),np.linalg.inv(ce),Yper.mean(0)-model(b,Uper.mean(0))[0]).real)/(rr@rr))

# FRF covariance vec ordering: dG=(dY-G dU) U^-1. Sum independent experiments.
U=cplx((K,4,4))+np.eye(4)[None]*5;G=cplx((K,4,4));Y=G@U;Uinv=np.linalg.inv(U)
dU=cplx((K,4,4));dY=cplx((K,4,4));e=dY-G@dU;full=(e@Uinv).transpose(0,2,1).reshape(K,16)
via=np.zeros_like(full)
for a in range(4):
    L=np.einsum('kj,oi->kjoi',Uinv[:,a,:],np.eye(4)).reshape(K,16,4)
    via+=np.einsum('kij,kj->ki',L,e[:,:,a])
vecdiff=float(np.linalg.norm(via-full)/np.linalg.norm(full))
eps=1e-7;new=(Y+eps*dY)@np.linalg.inv(U+eps*dU)
fdf=((new-G)/eps).transpose(0,2,1).reshape(K,16)
propdiff=float(np.linalg.norm(fdf-full)/np.linalg.norm(full))

# Compare new eigendecomposed simulation directly to original state-space dlsim.
signals=rng.normal(size=(1000,12));fast=ri.fast_sim(signals)
_,slow,_=ss.dlsim(ss.dlti(*ri._CL[:4],dt=1/5000),signals)
simdiff=float(np.linalg.norm(fast-slow)/np.linalg.norm(slow))
out=dict(sml_jacobian_relative_errors=err,mean_covariance_residual_identity_relative=covdiff,
         whitened_cost_vs_quadratic_relative=costdiff,FRF_vecF_operator_relative=vecdiff,
         FRF_linearization_finite_difference_relative=propdiff,fast_vs_dlsim_relative=simdiff)
assert max(err)<1e-7 and covdiff<1e-12 and costdiff<1e-12 and vecdiff<1e-12 and propdiff<1e-6 and simdiff<1e-10
(HERE/'comparator_independent_checks.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
