"""Deterministic controller/finite-record nuisance supplement; no fits or RNG."""
import check_jacobian as c
import numpy as np
import scipy.linalg as la
import json, hashlib
from pathlib import Path

def complete_rank_and_direct_check(X,T,cut=1e-12):
    X=X/np.maximum(la.norm(X,axis=0),1e-300)
    T=T/np.maximum(la.norm(T,axis=0),1e-300)
    ss=la.svdvals(np.c_[X,T],check_finite=False)
    direct=[]
    for j in range(7):
        Z=np.c_[X,T[:,[k for k in range(7) if k!=j]]]
        u,s,_=la.svd(Z,full_matrices=False,check_finite=False)
        q=u[:,s>cut*s[0]]
        v=T[:,j]-q@(q.T@T[:,j])
        direct.append(float(v@v))
    return dict(joint_column_normalized_singular_values=ss,
        joint_rank=int(sum(ss>cut*ss[0])),total_parameters=X.shape[1]+7,
        joint_column_normalized_condition=float(ss[0]/ss[-1]),
        direct_all_other_columns_delay_energy=direct)

def main():
    before=json.loads((c.OUT/'source_hashes_before.json').read_text(encoding='utf8'))
    assert all(hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest().lower()==r['SHA256'].lower() for r in before)
    res=c.fixed_residues()
    with np.load(c.OLD/'mc/data/P20_r000.npz') as a:
        keep=a['k']%int(a['P'])==0
        f=a['f'][keep];D=a['D'][keep]/(int(a['P'])*30000)
    results=[]
    for fc in (60.,15.):
        for band,(lo,hi,active) in c.BANDS.items():
            ix=(f>=lo-1e-10)&(f<=hi+1e-10);ff=f[ix];dd=D[ix]
            gg=c.full_response(ff,res);K=c.controller(ff,fc=fc)
            U=np.linalg.solve(np.eye(4)[None]+K@gg,dd[:,:,None])[:,:,0]
            # Both window-tail transients are free even though their poles/residues are known.
            checks=[('exact_tail_free_bg',list(range(6)))]
            if fc==15.:checks+=[('local_free_bg',active)]
            for case,ta in checks:
                X,T,meta=c.matrices(ff,active,res,case,'reference_output',U,transient_active=ta)
                info=c.project(X,T)
                verification=complete_rank_and_direct_check(X,T)
                main=next(t for t in info['thresholds'] if t['cut']==1e-12)
                verification['sequential_vs_direct_max_abs_error']=float(np.max(abs(np.asarray(main['all_other_parameters_delay_energy'])-verification['direct_all_other_columns_delay_energy'])))
                r=dict(controller_fc=fc,band=band,case=case,observation='reference_output',
                    all_six_modal_transients=ta==list(range(6)),frequency_points=len(ff),
                    max_sigma_GK_on_excitation_grid=float(np.linalg.svd(gg@K,compute_uv=False)[:,0].max()),
                    **meta,**info,**verification)
                results.append(r)
                c.save(c.OUT/f'FC{int(fc)}_{band}_{case}_T{len(ta)}.json',r)
                print(fc,band,case,'nparam',r['total_count'],'rank',r['joint_rank'],
                    'minr',min(main['all_other_parameters_delay_energy']),
                    'joint',main['joint_weakest_energy_ratio'],'check',verification['sequential_vs_direct_max_abs_error'],flush=True)
    c.save(c.OUT/'selected_design_results.json',dict(no_fit=True,no_simulation=True,no_random_generation=True,results=results))
    assert all(hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest().lower()==r['SHA256'].lower() for r in before)

if __name__=='__main__':main()
