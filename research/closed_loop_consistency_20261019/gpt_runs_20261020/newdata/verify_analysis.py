from uncertainty_bias import *
from prepare_data import fit_controller
import hashlib

def main():
    checks={};rows=[]
    for ctrl in CTRS:
        x=np.load(HERE/'prepared'/f'A_{ctrl}.npz');d={k:x[k] for k in ('U','Y','C','D')}
        A,B,coef,info=fit_controller(d);old=np.load(HERE/'controllers'/f'{ctrl}.npz');f=np.geomspace(3,25,150)
        Ka=evalK(f,A,B,old['T']);Kb=evalK(f,old['A'],old['B'],old['T'])
        r=dict(controller=ctrl,A_recovery_validation=info['validation_relative_rms'],K_A_vs_B_relative=float(np.linalg.norm(Ka-Kb)/np.linalg.norm(Kb)),
               max_frequency_K_relative=float(np.max(np.linalg.norm(Ka-Kb,axis=(1,2))/np.linalg.norm(Kb,axis=(1,2)))))
        fr=np.load(HERE/'analysis'/f'{ctrl}_orthogonal_frf.npz');Dn=fr['D']/np.linalg.norm(fr['D'],axis=2,keepdims=True)
        r['row_normalized_D_condition_max']=float(np.max(np.linalg.cond(Dn)))
        cov=np.load(HERE/'fits'/f'A_{ctrl}_local_covariance.npz')['covariance_residual_raw']/len(x['U'])
        noise=np.load(HERE/'analysis'/f'{ctrl}_noise.npz')['Phi']
        r['A_nonexcited_residual_to_B_trace_ratio']=float(np.trace(cov,axis1=1,axis2=2).real.sum()/np.trace(noise,axis1=1,axis2=2).real.sum())
        sp=prepare_spectra('A_'+ctrl);N=int(sp['N']);r['A_input_band_rms_per_channel']=np.sqrt(2*np.sum(abs(sp['U'])**2,axis=0)/N**2).tolist()
        r['A_output_band_rms_per_channel']=np.sqrt(2*np.sum(abs(sp['Y'])**2,axis=0)/N**2).tolist()
        rows.append(r)
    checks['independent_controller_and_excitation_checks']=rows
    # Exact covariance algebra tested on synthetic proper Gaussian vectors;
    # this is a small implementation check, not additional evidence for the paper.
    rng=np.random.default_rng(90020);n=40000;ny=4;np_=7
    L=(rng.normal(size=(ny,ny))+1j*rng.normal(size=(ny,ny)))*.2;C=L@L.conj().T
    J=rng.normal(size=(ny,np_))+1j*rng.normal(size=(ny,np_));B=rng.normal(size=(ny,ny,np_))+1j*rng.normal(size=(ny,ny,np_))
    H=(B+B.conj().transpose(1,0,2))/2
    pred_mean=np.einsum('ijp,ji->p',H,C).real
    HL=np.einsum('ai,abp,bj->ijp',L.conj(),H,L)
    pred_cov=.5*(J.conj().T@C@J).real+(HL.reshape(-1,np_).conj().T@HL.reshape(-1,np_)).real
    v=(rng.normal(size=(n,ny))+1j*rng.normal(size=(n,ny)))@L.T/np.sqrt(2)
    scores=(v@J.conj()).real+np.einsum('ni,ijp,nj->np',v.conj(),H,v).real
    chk=dict(covariance_relative_error=float(np.linalg.norm(np.cov(scores,rowvar=False)-pred_cov)/np.linalg.norm(pred_cov)),
             mean_max_standard_errors=float(np.max(abs(scores.mean(0)-pred_mean)/np.sqrt(np.diag(pred_cov)/n))),seed=90020,samples=n)
    assert chk['covariance_relative_error']<.06 and chk['mean_max_standard_errors']<5,chk
    checks['Gaussian_quadratic_score_algebra']=chk
    # Original baseline and raw inputs remain byte-identical.
    inv=json.loads((HERE/'data_inventory.json').read_text());bad=[]
    for r in inv:
        if hashlib.sha256(Path(r['source']).read_bytes()).hexdigest()!=r['sha256']:bad.append(r['id'])
    base=json.loads((HERE/'baseline/sources.json').read_text())
    for r in base:
        if hashlib.sha256(Path(r['source']).read_bytes()).hexdigest()!=r['sha256']:bad.append(r['source'])
    checks['source_integrity']=dict(unchanged=not bad,raw_records=len(inv),baseline_sources=len(base),changed=bad)
    js(HERE/'analysis/verification.json',checks);print(json.dumps(checks,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
