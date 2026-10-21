"""ONLY r=0..4. Requires completed deterministic/power/population acceptance.

Never dispatches r>=5. Exit after the 90 planned endpoints, regardless of results.
Each worker owns a controller/P/r data record, avoiding simultaneous file writes.
"""
from oracle_model import *
from profile_solver import Profile
from concurrent.futures import ProcessPoolExecutor,as_completed

def acceptance():
    checks={}
    checks['model']=json.loads((ROOT/'validation/model_acceptance.json').read_text())['acceptance_pass']
    checks['derivatives']=all(x['passed'] for x in json.loads((ROOT/'validation/profile_derivative_acceptance.json').read_text()))
    checks['implicit_derivatives']=all(x['passed'] for x in json.loads((ROOT/'validation/implicit_derivative_acceptance.json').read_text()))
    checks['population_derivative']=json.loads((ROOT/'validation/population_acceptance.json').read_text())['passed']
    checks['power']=len(list((ROOT/'power').glob('FC*/summary.json')))==18
    paths=list((ROOT/'population').glob('FC*/r[0-9][0-9][0-9].json'))
    checks['population_complete']=len(paths)==540
    checks['population_roots_pass']=all(json.loads(p.read_text())['success'] for p in paths) and len(paths)==540
    save(ROOT/'pilot/start_gate.json',checks)
    if not all(checks.values()):raise RuntimeError('Pilot NOT STARTED: prerequisite incomplete or failed: '+str(checks))

def job(fc,P,r):
    assert 0<=r<5,'Only approved pilot realizations are allowed'
    t0=time.perf_counter();rows=[]
    configs=[('M3','L'),('M3','H')]
    if P==20:configs=[('M1','L'),('M1','H'),('M1','W')]+configs
    try:data=generate(P,r,fc)
    except Exception as exc:
        for method,band in configs:
            folder=ROOT/'pilot'/f'FC{fc}_{band}_P{P}_{method}'/f'r{r:03d}'
            row=dict(fc=fc,P=P,r=r,method=method,band=band,success=False,status='generation_failed',error=repr(exc))
            save(folder/'result.json',row);rows.append(row)
        return rows
    for method,band in configs:
        folder=ROOT/'pilot'/f'FC{fc}_{band}_P{P}_{method}'/f'r{r:03d}'
        if (folder/'result.json').exists():
            rows.append(json.loads((folder/'result.json').read_text()));continue
        m=OracleModel(band);b0=m.initial(r);sub=select(data,band,method)
        save(folder/'initial.json',dict(beta=b0,data_seed=1002+r,initial_seed=100200+r,parameters=m.nb,
            f_min=float(min(sub['f'])),f_max=float(max(sub['f'])),bins=len(sub['f']),fc=fc,P=P,r=r,method=method,band=band))
        try:
            pr=Profile(m,sub,method);result=pr.solve(b0,folder)
            result.update(fc=fc,P=P,r=r,band=band,method=method)
            delta=np.asarray(result['beta'])[m.nl:]-m.truth()[m.nl:]
            result.update(frequency_error_percent=100*np.expm1(delta[:m.nm]),
                damping_error_percent=100*np.expm1(delta[m.nm:2*m.nm]),delay_error_us=delta[-7:]*100)
        except Exception as exc:
            result=dict(fc=fc,P=P,r=r,band=band,method=method,success=False,status='uncaught_exception',error=repr(exc))
        save(folder/'result.json',result);rows.append(result)
        print('PILOT',fc,P,r,method,band,result['success'],result.get('score'),result.get('seconds'),flush=True)
    return rows

if __name__=='__main__':
    acceptance();rows=[]
    configurations=[(fc,P,r) for r in range(5) for fc in (15,60) for P in (10,20,40)]
    with ProcessPoolExecutor(max_workers=4) as ex:
        for future in as_completed([ex.submit(job,*x) for x in configurations]):
            rows+=future.result();save(ROOT/'pilot/completed_status.json',rows)
    save(ROOT/'pilot/HARD_STOP.json',dict(completed=True,maximum_r=4,number_endpoints=len(rows),
        remaining_authorized=False,message='STOP. User confirmation required for r=5..29.'))
    print('PILOT FINISHED. HARD STOP. No r>=5 scheduled.',flush=True)
