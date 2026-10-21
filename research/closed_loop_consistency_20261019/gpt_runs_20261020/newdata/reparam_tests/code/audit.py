from profile_core import *
import hashlib

def main():
    hashes=rd(ROOT/'baseline/input_hashes.json');changed=[];missing=[]
    for name,expected in hashes.items():
        p=Path(name)
        if not p.exists():missing.append(name)
        elif hashlib.sha256(p.read_bytes()).hexdigest()!=expected:changed.append(name)
    full=[rd(p) for p in (ROOT/'raw').glob('*_M[13]_full.json') if not p.name.startswith('p7')]
    sims=[rd(p) for p in (ROOT/'raw').glob('p7_sim_P*_r[0-9][0-9][0-9].json')]
    inventory=[];failed=[];missing_folds=[]
    parents=list((ROOT/'raw').glob('*_M[13]_full*.json'))
    parents=[p for p in parents if not p.name.startswith('p7')]
    for p in parents:
        r=rd(p);P=int(np.load(SRC/'prepared'/f"{r['id']}.npz")['P'])
        if r.get('common_fixed_us') is not None:
            anchor=rd(ROOT/'analysis'/f"fixed_low_{r['id']}_{r['model']}_anchor.json")['common_fixed_us']
            # Offset anchor fits are sensitivity endpoints, not jackknife parents.
            if abs(r['common_fixed_us']-anchor)>1e-4:continue
        expected=P if r.get('success') else 0
        folds=[p.with_name(p.name.replace('_full',f'_omit{i}')) for i in range(P)]
        saved=[q for q in folds if q.exists()]
        pending=[str(q.relative_to(ROOT)) for q in folds if not q.exists()] if expected else []
        missing_folds+=pending
        inventory.append(dict(parent=str(p.relative_to(ROOT)),method=r['method'],fixed=r.get('common_fixed_us'),
            success=r.get('success',False),periods=P,expected_folds=expected,completed_folds=len(saved),
            failed_folds=sum(not rd(q).get('success') for q in saved),
            unrun_due_to_failed_parent=P if not expected else 0,missing_folds=pending))
    for p in (ROOT/'raw').glob('*.json'):
        if p.name.startswith('p7_sim'):continue
        r=rd(p)
        if 'success' in r and not r['success']:
            failed.append(dict(file=str(p.relative_to(ROOT)),method=r.get('method'),
                reason=r.get('exception') or r.get('solver',{}).get('status'),
                score=r.get('solver',{}).get('score')))
    sim_incomplete=[];sim_failures=[];sim_folds=0;sim_full=0
    for r in sims:
        P=r['P'];folds=r.get('folds',[])
        sim_folds+=len(folds);sim_full+=sum(k in r for k in ('M2','M3'))
        if sorted(x.get('omit',-1) for x in folds)!=list(range(P)) or any(k not in r for k in ('M2','M3','JK')):
            sim_incomplete.append(dict(P=P,replicate=r['replicate']))
        for name,x in [('M2',r.get('M2',{})),('M3',r.get('M3',{})),('JK',r.get('JK',{})),
                       *[(f"M2_omit{fold.get('omit')}",fold) for fold in folds]]:
            if not x.get('success'):
                sim_failures.append(dict(P=P,replicate=r['replicate'],method=name,
                    score=x.get('solver',{}).get('score'),solver=x.get('solver',{}).get('status')))
    real_expected=[];real_missing=[];real_failures=[]
    for rid in old.IDS:
        P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])
        for method in ('M2','M3'):
            tags=['full']+[f'omit{i}' for i in range(P)]
            if method=='M2':tags += [f'omit{i}_{j}' for i in range(P) for j in range(i+1,P)]
            for tag in tags:
                p=ROOT/'raw'/f'p7_real_{rid}_{method}_{tag}.json'
                real_expected.append(str(p.relative_to(ROOT)))
                if not p.exists():real_missing.append(str(p.relative_to(ROOT)))
                elif not rd(p).get('success'):real_failures.append(str(p.relative_to(ROOT)))
    shape_checks=[]
    for r in full:
        if 'beta' not in r:continue
        from metrics import model_from
        m,b=model_from(r);fb,T=m.expand(b)
        for j in range(m.nm):
            C0=np.array([fb[o*m.full.L+10*j:o*m.full.L+10*j+4] for o in range(4)])
            C1=np.array([fb[o*m.full.L+10*j+4:o*m.full.L+10*j+8] for o in range(4)])
            if m.kind in ('b','d'):
                sv=la.svdvals(C0);sv1=la.svdvals(np.c_[C0,C1])
                shape_checks.append(dict(id=r['id'],window=r['window'],model=m.kind,method=r['method'],mode=j,
                         rank1_fraction=float(sv[1]/max(sv[0],1e-300)),shared_fraction=float(sv1[1]/max(sv1[0],1e-300)),
                         unit_shape_error=abs(np.linalg.norm(m.shape(b,j)[0])-1)))
            if m.kind in ('c','d') and np.max(abs(C1))>1e-25:
                raise ValueError('A C1=0 model has nonzero C1')
    out=dict(original_files=len(hashes),changed=changed,missing=missing,original_unchanged=not(changed or missing),
             p6_full_endpoints=len(full),p6_full_required=144,p6_full_failures=sum(not r.get('success') for r in full),
             p7_sim_counts={str(P):sum(r.get('P')==P for r in sims) for P in (4,8,16)},
             p7_sim_required_per_P=40,
             p7_real_summaries=len(list((ROOT/'analysis').glob('p7_real_*.json'))),
             p6_jackknife_inventory=inventory,p6_missing_required_folds=missing_folds,
             p6_jackknife_counts={method:dict(
                 expected=sum(x['expected_folds'] for x in inventory if x['method']==method and x['fixed'] is None),
                 completed=sum(x['completed_folds'] for x in inventory if x['method']==method and x['fixed'] is None),
                 failed=sum(x['failed_folds'] for x in inventory if x['method']==method and x['fixed'] is None),
                 skipped_failed_parent=sum(x['unrun_due_to_failed_parent'] for x in inventory if x['method']==method and x['fixed'] is None))
                 for method in ('M1','M3')},
             fixed_low_jackknife_counts={method:dict(
                 expected=sum(x['expected_folds'] for x in inventory if x['method']==method and x['fixed'] is not None),
                 completed=sum(x['completed_folds'] for x in inventory if x['method']==method and x['fixed'] is not None),
                 failed=sum(x['failed_folds'] for x in inventory if x['method']==method and x['fixed'] is not None))
                 for method in ('M1','M3')},
             p7_sim_fold_count=sim_folds,p7_sim_fold_expected=1120,p7_sim_full_count=sim_full,p7_sim_full_expected=240,
             p7_sim_incomplete=sim_incomplete,p7_sim_failed_endpoints=sim_failures,
             p7_real_fits_expected=len(real_expected),p7_real_fits_completed=len(real_expected)-len(real_missing),
             p7_real_missing=real_missing,p7_real_failed=real_failures,
             all_required_fits_recorded=len(full)==144 and not missing_folds and not real_missing and not sim_incomplete
                 and len(list((ROOT/'raw').glob('*_M[13]_full_fixed*.json')))==72
                 and all(sum(r.get('P')==P for r in sims)==40 for P in (4,8,16)),
             failed_endpoint_index=failed,
             rank1_checks=shape_checks,
             all_rank1_checks_pass=all(x['rank1_fraction']<1e-12 and x['shared_fraction']<1e-12 and x['unit_shape_error']<1e-12 for x in shape_checks))
    js(ROOT/'analysis/audit.json',out)
    print({k:v for k,v in out.items() if k not in ('rank1_checks','p6_jackknife_inventory','failed_endpoint_index')},flush=True)
if __name__=='__main__':main()
