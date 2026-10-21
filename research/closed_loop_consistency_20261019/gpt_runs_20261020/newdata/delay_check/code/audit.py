from core import *
import hashlib,csv
from root_solver import moments_only
def main():
    before=rd(ROOT/'hash_before.json');changed=[];missing=[]
    for p,h in before.items():
        path=Path(p)
        if not path.exists():missing.append(p)
        elif hashlib.sha256(path.read_bytes()).hexdigest()!=h:changed.append(p)
    results=dict(protected_files=len(before),changed=changed,missing=missing,
                 fixed_files=len(list((ROOT/'raw').glob('fixed*.json'))),
                 joint_full_files=len(list((ROOT/'raw').glob('joint*full.json'))),
                 fold_files=len(list((ROOT/'raw').glob('joint*omit*.json'))),
                 simulation_files=len(list((ROOT/'raw').glob('sim_*.json'))))
    # Independent dense check of the sparse IV equations and the moment-only
    # function used by the simulation root solver.
    _,m,_,_=setup('H493',joint=True)
    r=rd(ROOT/'raw/joint_A_nominal_H493_M3_full.json');b=np.asarray(r['beta'])
    a=spectrum('A_nominal','H493');c,U,Y,Uh=inputs(m,a,'M3')
    y,J=m.output(c,b,U);_,Z=m.output(c,b,Uh(b));e=Y-y
    Ad=np.einsum('kop,koq->pq',Z.conj(),J,optimize=True).real
    gd=np.einsum('kop,ko->p',Z.conj(),e,optimize=True).real
    A,g,score,_=equations(m,c,U,Y,Uh(b),b)
    gm,sm=moments_only(m,c,U,Y,Uh,b)
    scale=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
    results['sparse_vs_dense_A_relative']=float(np.linalg.norm((Ad-A)*scale[:,None]*scale[None,:])/np.linalg.norm(A*scale[:,None]*scale[None,:]))
    results['moment_only_vs_equations_absolute_scaled']=float(np.max(abs((gm-g)*scale)))
    results['dense_vs_sparse_g_absolute_scaled']=float(np.max(abs((gd-g)*scale)))
    results['score_difference']=float(abs(sm-score))
    # All constraint gauges must remain exact, including failed endpoints.
    maxzero=0.;maxadd=0.
    for p in (ROOT/'raw').glob('joint*.json'):
        q=rd(p);d=q.get('delay')
        if not d:continue
        maxzero=max(maxzero,abs(sum(d['output_us'])),abs(sum(d['input_us'])))
        pred=d['common_us']+np.asarray(d['output_us'])[:,None]+np.asarray(d['input_us'])[None,:]
        maxadd=max(maxadd,np.max(abs(np.asarray(d['matrix_us'])-pred)))
    results['max_gauge_sum_us']=maxzero;results['max_additive_matrix_error_us']=maxadd
    files=[ROOT/'时延可行性报告.md',ROOT/'README.md']
    files +=[ROOT/'figures'/f'{name}.{ext}' for name in ['delay_and_damping','residual_comparison','simulation_delay_error'] for ext in ['png','pdf']]
    results['deliverables_present']=all(p.exists() and p.stat().st_size>100 for p in files)
    results['counts_match']=results['fixed_files']==252 and results['joint_full_files']==36 and results['fold_files']==246 and results['simulation_files']==40
    results['passed']=not changed and not missing and results['counts_match'] and results['deliverables_present'] and results['sparse_vs_dense_A_relative']<1e-10 and maxadd<1e-8
    js(ROOT/'analysis/audit.json',results)
    print((ROOT/'analysis/audit.json').read_text(encoding='utf8'))
    assert results['passed']
if __name__=='__main__':main()
