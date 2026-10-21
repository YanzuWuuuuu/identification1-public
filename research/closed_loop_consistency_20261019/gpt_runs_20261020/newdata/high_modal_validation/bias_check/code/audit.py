from common import *
import hashlib
from PIL import Image
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for a in iter(lambda:f.read(8*1024*1024),b''):h.update(a)
    return h.hexdigest()
def main():
    protected=rd(HERE/'protected_hashes_before.json')
    changed=[];missing=[]
    for r in protected:
        p=Path(r['path'])
        if not p.exists():missing.append(str(p))
        elif digest(p)!=r['sha256']:changed.append(str(p))
    rows=rd(HERE/'analysis/jackknife_summary.json')
    max_identity=0.;fits=[];min_neighbors=100000
    for r in rows:
        fold=r['folds'];P=len(fold);ob=np.array([x['observed'] for x in fold])
        min_neighbors=min(min_neighbors,*[x['minimum_neighbors'] for x in fold])
        for form in ('M','MminusC'):
            pr=np.array([x['pred_'+form] for x in fold])
            cov=(P-1)/P*np.sum((ob-ob.mean(0))*(pr-pr.mean(0)),axis=0)
            identity=jk_se(ob)**2+jk_se(pr)**2-2*cov
            err=np.max(abs(jk_se(ob-pr)**2-identity)/np.maximum(jk_se(ob-pr)**2,1e-30))
            max_identity=max(max_identity,float(err))
        for i in range(P):
            for method in ('M1','M3'):
                rr=rd(HERE/'raw'/f"{r['id']}_{r['window']}_omit{i}_{method}.json")
                fits.append(rr)
    figures=[]
    for p in sorted((HERE/'figures').glob('*.png')):
        with Image.open(p) as im:
            figures.append(dict(file=p.name,size=list(im.size),dpi=im.info.get('dpi'),pdf=p.with_suffix('.pdf').exists()))
    result=dict(protected_count=len(protected),changed=changed,missing=missing,
        model_copy_unchanged=digest(HERE/'baseline/modal_model.py')==digest(HIGH/'modal_model.py'),
        tau_copy_unchanged=digest(HERE/'baseline/tau_seconds.npy')==digest(SRC/'baseline/tau_seconds.npy'),
        full_cases=len(rows),jackknife_folds=sum(len(r['folds']) for r in rows),refit_endpoints=len(fits),
        refit_failures=sum(not r['converged'] for r in fits),
        refit_boundary=sum(r['near_bound'] for r in fits),
        refit_strict_failures=sum(not r['strict_score_pass'] for r in fits),
        before_polish_strict_failures=sum(r['solver']['score_max']>=1e-6 for r in fits),
        max_final_score=max(r['score'] for r in fits),
        joint_variance_identity_max_relative_error=max_identity,minimum_neighbors=min_neighbors,
        figures=figures,model_validation=rd(HERE/'analysis/validation.json')['passed'],
        code_hashes={str(p.relative_to(HERE)):digest(p) for p in sorted((HERE/'code').glob('*.py'))})
    result['passed']=bool(not changed and not missing and result['model_copy_unchanged'] and result['tau_copy_unchanged']
           and result['refit_endpoints']==164 and result['jackknife_folds']==82 and max_identity<1e-10
           and result['model_validation'] and len(figures)==4 and all(x['pdf'] for x in figures))
    js(HERE/'analysis/final_audit.json',result)
    print('AUDIT',result['passed'],'protected',len(protected),'changed',len(changed),'folds',result['jackknife_folds'],
          'fits',len(fits),'before-polish strict failures',result['before_polish_strict_failures'])
    if not result['passed']:raise RuntimeError('Audit failed')
if __name__=='__main__':main()
