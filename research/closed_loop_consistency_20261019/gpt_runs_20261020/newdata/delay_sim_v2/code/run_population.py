"""Deterministic auxiliary roots. No pilot or stochastic simulation is called."""
from oracle_model import *
from concurrent.futures import ProcessPoolExecutor,as_completed
from population_root import run

def job(x):
    P,band,fc,method=x
    base=ROOT/'power'/f'FC{fc}_{band}_P{P}_{method}'
    while not (base/'summary.json').exists():
        if (base/'error.json').exists():return dict(configuration=x,state='power_failed')
        time.sleep(5)
    rows=json.loads((base/'summary.json').read_text(encoding='utf8'))['rows']
    if not all(r['success'] for r in rows):return dict(configuration=x,state='power_incomplete')
    try:run(*x);return dict(configuration=x,state='finished')
    except Exception as exc:return dict(configuration=x,state='failed',error=repr(exc))

if __name__=='__main__':
    acceptance=json.loads((ROOT/'validation/population_acceptance.json').read_text())
    assert acceptance['passed']
    configs=[(20,b,fc,'M1') for fc in (15,60) for b in BANDS]
    configs +=[(P,b,fc,'M3') for fc in (15,60) for P in (10,20,40) for b in ('L','H')]
    rows=[]
    with ProcessPoolExecutor(max_workers=4) as ex:
        for fut in as_completed([ex.submit(job,x) for x in configs]):
            row=fut.result();rows.append(row);save(ROOT/'population/run_status.json',rows);print(row,flush=True)
