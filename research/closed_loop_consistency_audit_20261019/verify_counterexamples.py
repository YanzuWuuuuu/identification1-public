"""Small deterministic checks of the review; no experiment/model modification."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import root_scalar

out=Path(__file__).resolve().parent
w=np.linspace(0,np.pi/2,100001);q=np.exp(-1j*w)
Su=1-.5*q;Phiu=np.abs(Su)**2+.25;Phivu=-.5
integ=lambda x:float(np.trapezoid(np.real(x),w)/(2*np.pi))
rows=[]
for eps in [1e-2,1e-3,1e-4]:
    J=q*q+eps*q
    b=integ(J.conj()*Phivu);M=integ(np.abs(J)**2*Phiu)
    C=integ(q.conj()*Phivu)
    def score(theta):
        delta=theta*J+.5*theta*theta*q
        deriv=J+theta*q
        return integ(deriv.conj()*(-delta*Phiu+Phivu))
    sol=root_scalar(score,bracket=(-eps,eps),xtol=1e-15).root
    rows.append(dict(epsilon=eps,b=b,M=M,C=C,exact_local_root=sol,
                     paper_b_over_M=b/M,correct_b_over_MminusC=b/(M-C)))
wf=np.linspace(0,np.pi,100001);qf=np.exp(-1j*wf)
# Phi_v=I, W=diag(1,2), and stable strictly causal G_l F.
tr=qf/(1-.3*qf)+2*qf/(1+.2*qf)
non_scalar_bias=-float(np.trapezoid(tr.real,wf)/(2*np.pi))
result=dict(description='Analytic counterexamples independently quadrature-checked.',
            nonlinear_bias=rows,constant_nonscalar_W_bias=non_scalar_bias,
            assumptions='Stable SISO true G=q^-1, K=.5/(1-.5q^-1), independent unit white d/v; B=(0,pi/2).')
(out/'counterexample_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
