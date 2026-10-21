from common import *
from sensitivity import Model,band_spectrum
import hashlib,datetime
def main():
    before=rd(HERE/'protected_hashes_before.json');changed=[];missing=[]
    for name,item in before['files'].items():
        p=Path(name)
        if not p.exists():missing.append(name);continue
        h=hashlib.sha256()
        with p.open('rb') as f:
            for c in iter(lambda:f.read(8*1024*1024),b''):h.update(c)
        if h.hexdigest()!=item['sha256']:changed.append(name)
    # Independent central-difference check of the extra background block.
    a=band_spectrum('A_nominal',25);li=a['line_indices'];c=m.Ctx(a['f'][li],TAU,'zoh');U=a['U'][li]
    mo=Model(True);rng=np.random.default_rng(1921)
    b=np.r_[np.array(fit('A_nominal','M3')['beta']),np.zeros(26)]
    b[174]=np.log(2*np.pi*40);b[175]=np.log(.2)
    b[158:174]=rng.normal(size=16)*1e-7;b[176:]=rng.normal(size=8)*1e-8
    y,J=mo.out(c,b,U);checks=[]
    for q in (158,164,173,174,175,176,183):
        h=1e-5 if q in (174,175) else 1e-7
        bp=b.copy();bm=b.copy();bp[q]+=h;bm[q]-=h
        fd=(mo.out(c,bp,U,False)[0]-mo.out(c,bm,U,False)[0])/(2*h)
        checks.append(dict(parameter=q,relative_error=np.linalg.norm(fd-J[:,:,q])/np.linalg.norm(fd)))
    z=b.copy();z[158:174]=0;z[176:]=0
    nested=np.max(abs(mo.out(c,z,U,False)[0]-m.out_jac(c,z[:158],U,False)[0]))
    cv=rd(HERE/'analysis/curvature.json');jk=rd(HERE/'analysis/prediction_jackknife.json');ss=rd(HERE/'analysis/sensitivity.json')
    oldjk={r['id']:r for r in rd(BASE/'analysis/jackknife.json')}
    seerr=max(np.max(abs(np.array(r['observed_se'])-np.array(oldjk[r['id']]['paired_se']))) for r in jk)
    covidentity=[]
    for r in jk:
        for key in ('M','MminusC'):
            lhs=np.array(r['difference_'+key+'_se'])**2
            rhs=np.array(r['observed_se'])**2+np.array(r['pred_'+key+'_se'])**2-2*np.array(r['observed_prediction_JK_covariance_'+key])
            covidentity.append(np.max(abs(lhs-rhs))/max(np.max(lhs),1e-30))
    sources_equal=all((HERE/'baseline'/n).read_bytes()==(BASE/'baseline'/n).read_bytes()
        for n in ('realfit_extended.py','realfit_original_copy.py','weighted_core.py','tau_seconds.npy'))
    out=dict(time=datetime.datetime.now().isoformat(),protected_files=len(before['files']),changed=changed,missing=missing,
        all_protected_unchanged=not(changed or missing),model_copies_identical=sources_equal,
        curvature_cases=len(cv),max_resolvent_error=max(r['resolvent_relative_error'] for r in cv),
        max_profile_reproduction_error=max(r['profile_max_abs_reproduction_error'] for r in cv),
        prediction_folds=sum(r['P'] for r in jk),prediction_fold_failures=sum(r['prediction_failures'] for r in jk),
        observed_jackknife_SE_reproduction_error=seerr,joint_variance_identity_relative_error=max(covidentity),
        sensitivity_configurations=len(ss),sensitivity_endpoints=sum(len(r['fits']) for r in ss),
        sensitivity_solver_failures=sum(not f['converged'] for r in ss for f in r['fits'].values()),
        sensitivity_strict_failures=sum(not f['strict_score_pass'] for r in ss for f in r['fits'].values()),
        baseline_reproduction=[dict(id=r['id'],**r['baseline_reproduction']) for r in ss if 'baseline_reproduction' in r],
        extra_block_jacobian=checks,nested_model_max_error=nested,
        figures=[p.name for p in (HERE/'figures').glob('*.png')],
        code_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'code').glob('*.py')})
    js(HERE/'analysis/verification.json',out)
    assert out['all_protected_unchanged'] and sources_equal
    assert nested==0 and max(x['relative_error'] for x in checks)<1e-5
    assert seerr<1e-12 and max(covidentity)<1e-10
    assert len(ss)==36 and out['prediction_folds']==41
    print('VERIFIED',out['protected_files'],'protected files;',out['prediction_folds'],'prediction folds;',out['sensitivity_endpoints'],'sensitivity endpoints')
if __name__=='__main__':main()
