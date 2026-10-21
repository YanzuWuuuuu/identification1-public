"""Independent algebra/scale diagnostics for full-matrix weighting and local SML."""
from sim_d_run import *

def main():
    a,z,U,Y,D,K,bt,b0,model,model_full,Uh=context(16,0)
    p=OUT/'results'/'P16_r000_M3.json'; b=np.array(json.loads(p.read_text(encoding='utf-8'))['beta'])
    cv=local_covariances(a['k'],a['U'],a['Y'],b,model_full,16,16,.1)
    G=c.resp(z,b);ny=4;Cx=cv['covariance_joint_raw'];A=np.concatenate([np.broadcast_to(np.eye(4),(len(z),4,4)),-G],axis=2)
    projected=A@Cx@A.conj().transpose(0,2,1)
    trans,_=model_full(b,np.zeros_like(a['U']))
    ix=cv['indices']; frozen=a['Y'][ix]-trans[ix]-np.einsum('koi,kni->kno',G,a['U'][ix])
    frozen_cov=np.einsum('kni,knj->kij',frozen,frozen.conj())/16
    identity=np.linalg.norm(projected-frozen_cov)/np.linalg.norm(projected)
    W=precision(cv['covariance_residual']);L=whitener(W)
    yhat,J=model(b,U);Z=model(b,Uh(b))[1];E=Y-yhat
    AW=np.real(np.einsum('kop,kol,klq->pq',Z.conj(),W,J,optimize=True))
    gw=np.real(np.einsum('kop,kol,kl->p',Z.conj(),W,E,optimize=True))
    Aw,gg,zw,ew=iv_matrices(model,U,Y,Uh,b,L)
    Ce=cv['CY']-cv['CYU']@G.conj().transpose(0,2,1)-G@cv['CYU'].conj().transpose(0,2,1)+G@cv['CU']@G.conj().transpose(0,2,1)
    out=dict(raw_joint_projection_equals_frozen_center_residual_covariance_relative=float(identity),
      local_Gj_versus_center_Gk_covariance_relative=float(np.linalg.norm(projected-cv['covariance_residual_raw'])/np.linalg.norm(projected)),
      weighted_normal_matrix_whitening_relative=float(np.linalg.norm(AW-Aw)/np.linalg.norm(AW)),
      weighted_moment_whitening_relative=float(np.linalg.norm(gw-gg)/np.linalg.norm(gw)),
      covariance_times_precision_relative=float(np.linalg.norm(cv['covariance_residual']@W-np.eye(4))/np.sqrt(4*len(z))),
      shrink_maps_not_identical_relative=float(np.linalg.norm(Ce-cv['covariance_residual'])/np.linalg.norm(cv['covariance_residual'])),
      max_excited_neighbor_count=int(np.sum(a['k'][cv['indices']]%16==0)))
    out['pass']=bool(identity<1e-8 and out['weighted_normal_matrix_whitening_relative']<1e-10 and out['weighted_moment_whitening_relative']<1e-10 and out['max_excited_neighbor_count']==0)
    write_json(OUT/'diagnostics'/'matrix_covariance_validation.json',out);print(json.dumps(out))
if __name__=='__main__':main()
