from run_p6 import *

def one_record(rid,kind):
    os.environ['REPARAM_WARM_HYBRID']='1'
    hh=[rd(ROOT/'raw'/f'{rid}_{w}_{kind}_M3_full.json') for w in ('H300','H493')]
    if not all(h['success'] for h in hh):
        js(ROOT/'analysis'/f'fixed_low_{rid}_{kind}_unavailable.json',
           dict(reason='At least one matching high-window M3 parent did not converge. No substitute selected.'))
        return
    common=float(np.mean([h['delay']['common_us'] for h in hh]))
    js(ROOT/'analysis'/f'fixed_low_{rid}_{kind}_anchor.json',
       dict(common_fixed_us=common,rule='Equal mean of the same-record two high-window M3 common delays; not selected by low-window residual.'))
    P=int(np.load(SRC/'prepared'/f'{rid}.npz')['P'])
    for method in ('M3','M1'):
        parent=fit_one(rid,'L',kind,method,fixed=common)
        if parent['success']:
            for i in range(P):fit_one(rid,'L',kind,method,omit=i,fixed=common)
        # A non-converged parent is kept, not assigned a spurious SE.
        for delta in (-25.,25.):fit_one(rid,'L',kind,method,fixed=common+delta)

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=4) as pool:
        for fu in as_completed([pool.submit(one_record,rid,kind) for rid in old.IDS for kind in ('a','b')]):fu.result()
