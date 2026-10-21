"""Deterministic predictions only; never starts simulated MC endpoints."""
from oracle_model import *
from concurrent.futures import ProcessPoolExecutor,as_completed
from population_math import run

def job(x):
    P,band,fc,method=x
    folder=ROOT/'power'/f'FC{fc}_{band}_P{P}_{method}'
    if (folder/'summary.json').exists():return {'configuration':x,'state':'existing'}
    try:run(*x);return {'configuration':x,'state':'finished'}
    except Exception as exc:
        row={'configuration':x,'state':'failed','error':repr(exc)}
        save(folder/'error.json',row);return row

if __name__=='__main__':
    configs=[(20,b,fc,'M1') for fc in (15,60) for b in BANDS]
    configs +=[(P,b,fc,'M3') for fc in (15,60) for P in (10,20,40) for b in ('L','H')]
    rows=[]
    with ProcessPoolExecutor(max_workers=2) as ex:
        for fut in as_completed([ex.submit(job,x) for x in configs]):
            row=fut.result();rows.append(row);save(ROOT/'power/run_status.json',rows);print(row,flush=True)
