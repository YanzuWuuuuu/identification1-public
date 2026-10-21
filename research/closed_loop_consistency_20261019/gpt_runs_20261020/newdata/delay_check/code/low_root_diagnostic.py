from core import *
from root_solver import exact_root
from concurrent.futures import ProcessPoolExecutor,as_completed
def one(rid,method):
    dest=ROOT/'diagnostics'/f'low_fullroot_{rid}_{method}.json'
    if dest.exists():return rd(dest)
    a=spectrum(rid,'L');_,m,lb,ub=setup('L',joint=True);c,U,Y,Uh=inputs(m,a,method)
    old=rd(ROOT/'raw'/f'joint_{rid}_L_{method}_full.json');base,_=basefit(rid,'L',method)
    out=dict(id=rid,method=method,attempts=[])
    for name,b0 in [('previous_endpoint',np.asarray(old['beta'])),('zero_delay_increment',np.r_[base,np.zeros(7)])]:
        t=time.perf_counter()
        try:
            b,info=exact_root(m,c,U,Y,Uh,b0,maxfev=2000)
            f,z=m.modal(b)
            out['attempts'].append(dict(initial=name,beta=b,solver=info,delay=delay_values(m.delay(b)),
                 frequency_hz=f,damping_percent=z,seconds=time.perf_counter()-t))
        except Exception as e:out['attempts'].append(dict(initial=name,exception=repr(e),seconds=time.perf_counter()-t))
    js(dest,out)
    print(rid,method,[(x['initial'],x.get('solver',{}).get('success'),x.get('solver',{}).get('score')) for x in out['attempts']],flush=True)
    return out
if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(one,rid,m) for rid in IDS for m in ['M1','M3']]):f.result()
