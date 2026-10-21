from common import *
from noise_bias import cache_moments,predict
GROUPS={'modal_residues':np.arange(96),'modal_poles':MOD,'D0D2':np.arange(96,128),'transients':np.arange(134,158)}

def inspect(M,C,b):
    H=scaled_solve(M,C);delta=scaled_solve(M,b);dc=scaled_solve(M-C,b)
    s=1/np.sqrt(np.abs(np.diag(M)));Ms=M*s[:,None]*s[None,:];Cs=C*s[:,None]*s[None,:]
    A=Ms-Cs;Hs=np.linalg.solve(Ms,Cs)
    eig,v=sla.eigh(Cs,Ms);q=np.max(abs(eig))
    lead=np.argsort(abs(eig))[-6:][::-1]
    out=dict(raw_H_norm2=np.linalg.norm(H,2),diagonal_scaled_H_norm2=np.linalg.norm(Hs,2),
        eigen_min=eig.min(),eigen_max=eig.max(),max_abs_eigenvalue=q,nearest_to_one=eig[np.argmin(abs(1-eig))],
        distance_to_one=np.min(abs(1-eig)),raw_condition_M=np.linalg.cond(M),raw_condition_MminusC=np.linalg.cond(M-C),
        scaled_condition_M=np.linalg.cond(Ms),scaled_condition_MminusC=np.linalg.cond(A),
        M_positive_definite=True,energy_inverse_gain=np.max(1/abs(1-eig)),
        energy_relative_correction_bound=q/(1-q) if q<1 else None,
        diagnostic='no_Neumann_guarantee' if q>=1 else 'large_correction_warning' if q>=.5 else 'moderate_correction' if q>.1 else 'small_curvature',
        pred_M=100*np.expm1(delta[MOD]),pred_MminusC=100*np.expm1(dc[MOD]))
    direct=[]
    normC=np.linalg.norm(Cs)
    for name,idx in GROUPS.items():
        direct.append(dict(group=name,standardized_C_row_frobenius_fraction=np.linalg.norm(Cs[idx,:])**2/normC**2,
                           raw_C_row_norm=np.linalg.norm(C[idx,:]),raw_C_column_norm=np.linalg.norm(C[:,idx])))
    contribution=[]
    denom=dc[132]-delta[132]
    factor=100*(np.exp(dc[132])-np.exp(delta[132]))/denom if abs(denom)>1e-15 else 100*np.exp(delta[132])
    for name,idx in GROUPS.items():
        d=scaled_solve(M-C,C[:,idx]@delta[idx])
        contribution.append(dict(group=name,zeta2_log_correction=d[132],zeta2_percent_point_correction=d[132]*factor,
                                 correction_vector=d))
    eigrows=[]
    for j in lead:
        x=v[:,j];frac={name:float(np.sum(abs(x[idx])**2)/np.sum(abs(x)**2)) for name,idx in GROUPS.items()}
        eigrows.append(dict(eigenvalue=eig[j],diagonal_scaled_eigenvector_group_fraction=frac))
    out.update(direct_C_groups=direct,correction_by_source=contribution,leading_eigendirections=eigrows,
        decomposition_error=float(abs(sum(x['zeta2_log_correction'] for x in contribution)-(dc[132]-delta[132]))),
        resolvent_relative_error=float(np.linalg.norm(dc-delta-scaled_solve(M-C,C@delta))/np.linalg.norm(dc-delta)))
    u=delta/s;du=(dc-delta)/s
    out['relative_correction_in_M_metric']=np.sqrt(du@Ms@du/(u@Ms@u))
    out['zeta2_relative_change_between_predictions']=abs(out['pred_MminusC'][4]-out['pred_M'][4])/abs(out['pred_M'][4])
    # Profile out all linear nuisance directions: both the Schur matrix and RHS change.
    nuisance=np.setdiff1d(np.arange(m.NB),MOD)
    def profile(A):
        Ann=A[np.ix_(nuisance,nuisance)];Apn=A[np.ix_(MOD,nuisance)]
        S=A[np.ix_(MOD,MOD)]-Apn@scaled_solve(Ann,A[np.ix_(nuisance,MOD)])
        rhs=b[MOD]-Apn@scaled_solve(Ann,b[nuisance])
        return S,rhs
    S0,r0=profile(M);S1,r1=profile(M-C)
    out['profile_poles']=dict(
        delta_M=scaled_solve(S0,r0),delta_MminusC=scaled_solve(S1,r1),
        nuisance_C_nn_norm=np.linalg.norm(C[np.ix_(nuisance,nuisance)]),
        rhs_relative_change=np.linalg.norm(r1-r0)/np.linalg.norm(r0),
        schur_condition_M=diagnostics(S0),schur_condition_MminusC=diagnostics(S1))
    out['profile_max_abs_reproduction_error']=float(max(np.max(abs(out['profile_poles']['delta_M']-delta[MOD])),
        np.max(abs(out['profile_poles']['delta_MminusC']-dc[MOD]))))
    return out,dict(M=M,C=C,b=b,H=H,eigenvalues=eig,eigenvectors_standardized=v,delta_M=delta,delta_MminusC=dc)

def main():
    rows=[]
    for rid in IDS:
        a=np.load(BASE/'raw'/f'{rid}_prediction_0.5.npz')
        out,arrays=inspect(a['M'],a['Ccurvature'],a['b'])
        out.update(id=rid,scenario='original',alpha=1.)
        js(HERE/'analysis'/f'{rid}_curvature.json',out)
        np.savez_compressed(HERE/'raw'/f'{rid}_curvature.npz',**arrays);rows.append(out)
        if rid.startswith('A_'):
            alpha=rd(BASE/'analysis'/f'{rid}_noise_bias.json')['scaled_excitation']['alpha']
            S=np.ones(m.NB);S[:m.ITR]=alpha
            Ms=a['M']-a['M0']+a['M0']*S[:,None]*S[None,:]
            so,sa=inspect(Ms,a['Ccurvature'],a['b'])
            so.update(id=rid,scenario='scaled',alpha=alpha)
            js(HERE/'analysis'/f'{rid}_scaled_curvature.json',so)
            np.savez_compressed(HERE/'raw'/f'{rid}_scaled_curvature.npz',**sa);rows.append(so)
        print('CURVATURE',rid,'q',out['max_abs_eigenvalue'],'range',out['eigen_min'],out['eigen_max'],
              'zeta2',out['pred_M'][4],out['pred_MminusC'][4],flush=True)
    js(HERE/'analysis/curvature.json',rows)
if __name__=='__main__':main()
