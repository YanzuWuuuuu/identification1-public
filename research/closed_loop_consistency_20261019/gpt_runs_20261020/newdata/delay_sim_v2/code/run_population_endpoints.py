"""Load-balanced deterministic roots, four workers; no simulated data."""
from population_root import *
from concurrent.futures import ProcessPoolExecutor,as_completed
from functools import lru_cache
import shutil

@lru_cache(maxsize=2)
def setup(P,band,fc,method):
    base=ROOT/'power'/f'FC{fc}_{band}_P{P}_{method}'
    with np.load(base/'common.npz') as a:common={k:a[k] for k in a.files}
    return Population(P,band,fc,method,common),common

def job(case):
    P,band,fc,method,r=case;folder=ROOT/'population'/f'FC{fc}_{band}_P{P}_{method}'
    folder.mkdir(exist_ok=True,parents=True);path=folder/f'r{r:03d}.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old['success']:return dict(case=case,state='existing',success=True)
        ar=folder/'pre_broyden';ar.mkdir(exist_ok=True)
        if not (ar/path.name).exists():shutil.copy2(path,ar/path.name)
    pop,common=setup(P,band,fc,method)
    base=ROOT/'power'/f'FC{fc}_{band}_P{P}_{method}'
    with np.load(base/f'r{r:03d}.npz') as a:power={k:a[k] for k in a.files}
    try:
        row=run_one(pop,folder,r,common,power)
        return dict(case=case,state=row['status'],success=row['success'],seconds=row['seconds'])
    except Exception as exc:
        row=dict(r=r,success=False,error=repr(exc));save(path,row)
        return dict(case=case,state='exception',success=False,error=repr(exc))

if __name__=='__main__':
    configs=[(20,b,fc,'M1') for fc in (15,60) for b in BANDS]
    configs +=[(P,b,fc,'M3') for fc in (15,60) for P in (10,20,40) for b in ('L','H')]
    jobs=[(*c,r) for c in configs for r in range(30)];rows=[]
    with ProcessPoolExecutor(max_workers=4) as ex:
        for future in as_completed([ex.submit(job,c) for c in jobs]):
            row=future.result();rows.append(row)
            save(ROOT/'population/endpoint_status.json',rows)
            if row['state']!='existing':print(row,flush=True)
    save(ROOT/'population/COMPLETE.json',dict(endpoints=len(rows),successes=sum(r['success'] for r in rows),
        stochastic_data_generated=False))
