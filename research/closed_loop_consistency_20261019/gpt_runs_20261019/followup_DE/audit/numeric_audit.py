"""Independent read-only algebra checks; no Monte Carlo or estimator changes."""
import sys, os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
sys.dont_write_bytecode=True
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from sim_d_run import context,c,OUT
from weighted_core import local_covariances,precision,whitener,iv_matrices,shrink_cov
from sml_core_copy import SMLCovResidual
from iv_lm import moments,exact_moment_jacobian
import numpy as np, json

def main():
    rows=[]
    for P in (2,4,8,16):
        a,z,U,Y,D,K,bt,b0,model,model_full,Uh=context(P,0)
        b=np.array(json.loads((OUT/'results'/f'P{P:02d}_r000_M3.json').read_text(encoding='utf-8'))['beta'])
        G=c.resp(z,b); T=model_full(b,np.zeros_like(a['U']))[0]
        for nl in (8,16):
            cv=local_covariances(a['k'],a['U'],a['Y'],b,model_full,P,nl,.1)
            A=np.concatenate([np.broadcast_to(np.eye(4),(len(z),4,4)),-G],axis=2)
            projected=A@cv['covariance_joint_raw']@A.conj().transpose(0,2,1)
            ix=cv['indices']; frozen=a['Y'][ix]-T[ix]-np.einsum('koi,kni->kno',G,a['U'][ix])
            frozenC=np.einsum('kni,knj->kij',frozen,frozen.conj())/nl
            GH=G.conj().transpose(0,2,1)
            Ce=cv['CY']-cv['CYU']@GH-G@cv['CYU'].conj().transpose(0,2,1)+G@cv['CU']@GH
            Li=np.linalg.inv(np.linalg.cholesky(cv['covariance_residual']))
            ratios=np.linalg.eigvalsh(Li@Ce@Li.conj().transpose(0,2,1))
            W=precision(cv['covariance_residual']);L=whitener(W)
            yh,J=model(b,U);Z=model(b,Uh(b))[1];E=Y-yh
            Ad=np.real(np.einsum('kop,kol,klq->pq',Z.conj(),W,J,optimize=True))
            gd=np.real(np.einsum('kop,kol,kl->p',Z.conj(),W,E,optimize=True))
            Aw,gw,_,_=iv_matrices(model,U,Y,Uh,b,L)
            Al,gl,_=moments(b,U,Y,model,Uh,W)
            gap=np.linalg.norm(projected-cv['covariance_residual_raw'],axis=(1,2))/np.linalg.norm(cv['covariance_residual_raw'],axis=(1,2))
            f=a['f'][cv['line_indices']]; near=np.any(abs(f[:,None]-c.fsup[None,:])<1,axis=1)
            row={'P':P,'replicate':0,'nl':nl,
              'projection_identity_relative':float(np.linalg.norm(projected-frozenC)/np.linalg.norm(projected)),
              'raw_local_Gj_vs_center_Gk_global_relative':float(np.linalg.norm(projected-cv['covariance_residual_raw'])/np.linalg.norm(cv['covariance_residual_raw'])),
              'raw_local_Gj_vs_center_Gk_median_near_support':float(np.median(gap[near])),
              'raw_local_Gj_vs_center_Gk_max':float(gap.max()),
              'raw_gap_peak_frequency_Hz':float(f[np.argmax(gap)]),
              'joint_vs_residual_shrink_map_alone_relative':float(np.linalg.norm(Ce-shrink_cov(projected,.1))/np.linalg.norm(shrink_cov(projected,.1))),
              'total_M5n_Ce_vs_M3w_Cv_relative':float(np.linalg.norm(Ce-cv['covariance_residual'])/np.linalg.norm(cv['covariance_residual'])),
              'M5n_Ce_vs_M3w_Cv_generalized_eigenratio_global_min':float(ratios.min()),
              'M5n_Ce_vs_M3w_Cv_generalized_eigenratio_global_max':float(ratios.max()),
              'M5n_Ce_vs_M3w_Cv_generalized_eigenratio_median_min':float(np.median(ratios[:,0])),
              'M5n_Ce_vs_M3w_Cv_generalized_eigenratio_median_max':float(np.median(ratios[:,-1])),
              'median_halfspan_Hz':float(np.median(cv['distances'].max(1))*c.fs/(30000*P)),
              'max_halfspan_Hz':float(cv['distances'].max()*c.fs/(30000*P)),
              'whitened_A_relative':float(np.linalg.norm(Ad-Aw)/np.linalg.norm(Ad)),
              'whitened_g_relative':float(np.linalg.norm(gd-gw)/np.linalg.norm(gd)),
              'lm_A_relative':float(np.linalg.norm(Al-Aw)/np.linalg.norm(Aw)),
              'lm_g_relative':float(np.linalg.norm(gl-gw)/np.linalg.norm(gw)),
              'excited_neighbor_count':int(np.count_nonzero(a['k'][ix]%P==0))}
            rows.append(row)
    # Independent finite-difference check of complex Cholesky SML derivative.
    rng=np.random.default_rng(801923);nf=9;ny=4;nu=4;nb=3
    G0=rng.normal(size=(nf,ny,nu))+1j*rng.normal(size=(nf,ny,nu))
    Gj=.1*(rng.normal(size=(nf,ny,nu,nb))+1j*rng.normal(size=(nf,ny,nu,nb)))
    U=rng.normal(size=(nf,nu))+1j*rng.normal(size=(nf,nu))
    Y=rng.normal(size=(nf,ny))+1j*rng.normal(size=(nf,ny))
    H=rng.normal(size=(nf,8,8))+1j*rng.normal(size=(nf,8,8))
    C=H@H.conj().transpose(0,2,1)+np.eye(8)[None]
    def gm(b): return G0+np.einsum('koip,p->koi',Gj,b),Gj
    def md(b,u):
        G,_=gm(b);return np.einsum('koi,ki->ko',G,u),np.einsum('koip,ki->kop',Gj,u)
    obj=SMLCovResidual(U,Y,C[:,:4,:4],C[:,4:,4:],C[:,:4,4:],md,gm)
    b=np.array([.3,-.2,.1]);r,J=obj.calc(b); h=1e-5
    Jfd=np.column_stack([(obj.fun(b+np.eye(nb)[j]*h)-obj.fun(b-np.eye(nb)[j]*h))/(2*h) for j in range(nb)])
    sml={'residual_jacobian_central_diff_relative':float(np.linalg.norm(J-Jfd)/np.linalg.norm(Jfd)),
         'step':h,'includes_complex_offdiagonal_covariance':True}
    # Oracle covariance equals .5*A^-1 under proper-complex independent-bin
    # noise, and its modal block equals the nuisance-profiled Schur result.
    a,z,U,Y,D,K,bt,b0,model,model_full,Uh=context(16,0)
    G0=c.resp(z,bt);Pv=30000*16*((3e-7)**2*np.eye(4)[None]+c.SW**2*G0@G0.conj().transpose(0,2,1))
    Wo=precision(Pv);_,Z=model(bt,Uh(bt));WZ=np.einsum('kij,kjp->kip',Wo,Z)
    Ao=np.real(Z.reshape(-1,c.NB).conj().T@WZ.reshape(-1,c.NB))
    Qo=.5*np.real(np.einsum('kip,kij,kjq->pq',WZ.conj(),Pv,WZ,optimize=True))
    sc=1/np.sqrt(np.diag(Ao));As=Ao*sc[:,None]*sc[None,:];Qs=Qo*sc[:,None]*sc[None,:]
    C=np.linalg.solve(As,np.linalg.solve(As,Qs).T).T*sc[:,None]*sc[None,:]
    idx=np.arange(c.IW,c.ITR);nuis=np.r_[0:c.IW,c.ITR:c.NB]
    S=As[np.ix_(idx,idx)]-As[np.ix_(idx,nuis)]@np.linalg.solve(As[np.ix_(nuis,nuis)],As[np.ix_(nuis,idx)])
    Cp=.5*np.linalg.inv(S)*sc[idx,None]*sc[None,idx]
    oracle={'Q_equals_half_A_relative':float(np.linalg.norm(Qs-.5*As)/np.linalg.norm(Qs)),
            'full_inverse_vs_nuisance_Schur_modal_covariance_relative':float(np.linalg.norm(C[np.ix_(idx,idx)]-Cp)/np.linalg.norm(Cp)),
            'minimum_scaled_information_eigenvalue':float(np.linalg.eigvalsh(As).min()),
            'condition_scaled':float(np.linalg.cond(As)),
            'n_structural_plus_transient_parameters':c.NB,'n_interest_parameters':len(idx),'n_nuisance_parameters':len(nuis)}
    out={'rows':rows,'sml_derivative':sml,'oracle_sandwich':oracle}
    (HERE/'numeric_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
