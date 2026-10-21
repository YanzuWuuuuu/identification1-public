import run_e as e
import numpy as np,json
from iv_lm import moments,exact_moment_jacobian
from weighted_core import iv_matrices,whitener,local_covariances,precision
c=e.c
a=np.load(e.HERE.parent/'experiment_a'/'data'/'P02_r000.npz')
k,z,D,U,Y=[a[x] for x in ('k','zk','D','U','Y')];sel=k%2==0
zs,Ds,Us,Ys=z[sel],D[sel],U[sel],Y[sel];Ks=c.Ktf(zs)
b=np.array(json.loads((e.HERE.parent/'experiment_a'/'results'/'P02_r000_dir_lines_W1.json').read_text())['beta'])
model=lambda b,u:c.out_jac(zs,b,u)
instrument=lambda b:np.linalg.solve(np.eye(4)[None]+Ks@c.resp(zs,b),Ds[:,:,None])[:,:,0]
cv=local_covariances(k,U,Y,b,lambda b,u:c.out_jac(z,b,u),2)
W=precision(cv['covariance_residual']);out={}
A,g,sc=moments(b,Us,Ys,model,instrument,W)
A2,g2,_,_=iv_matrices(model,Us,Ys,instrument,b,whitener(W))
out['matrix_weight_parity_A']=float(np.linalg.norm(A-A2)/np.linalg.norm(A2))
out['matrix_weight_parity_g']=float(np.linalg.norm(g-g2)/np.linalg.norm(g2))
J,g,sc=exact_moment_jacobian(b,Us,Ys,model,instrument,W,c.IW,c.ITR,Ks,lambda b:c.resp(zs,b))
rng=np.random.default_rng(20261019)
directions=[]
for q in range(3):
    v=rng.normal(size=c.NB)*np.maximum(abs(b),1e-10);v[c.IW:c.ITR]=rng.normal(size=6)*.1
    eps=1e-4
    fd=(moments(b+eps*v,Us,Ys,model,instrument,W)[1]-moments(b-eps*v,Us,Ys,model,instrument,W)[1])/(2*eps)
    # Compare scaled rows to avoid units allowing large residue rows to dominate.
    scale=1/np.maximum(np.sqrt(abs(np.diag(A))),1e-30)
    directions.append(float(np.linalg.norm((fd+J@v)*scale)/np.linalg.norm(fd*scale)))
out['moment_derivative_relative_errors']=directions
old=json.loads((e.HERE.parent/'experiment_a'/'results'/'P01_r000_ref_band_W1.json').read_text())
new=json.loads((e.OUT/'P01_r000_M3_E1.json').read_text())
out['E1_prior_log_pole_max_diff']=float(np.max(np.abs(np.array(old['beta'])[c.IW:c.ITR]-np.array(new['beta'])[c.IW:c.ITR])))
out['passed']=out['matrix_weight_parity_A']<1e-10 and out['matrix_weight_parity_g']<1e-10 and max(directions)<.005 and out['E1_prior_log_pole_max_diff']<1e-5
(e.HERE/'validation_e.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
assert out['passed']
