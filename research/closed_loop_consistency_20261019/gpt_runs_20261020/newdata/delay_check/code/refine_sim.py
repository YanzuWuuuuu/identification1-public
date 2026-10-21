from run_sim import *
from root_solver import exact_root
def one_ref(i):
    dest=ROOT/'refinements'/f'sim_{i:02d}.json'
    if dest.exists():return rd(dest)
    r=rd(ROOT/'raw'/f'sim_{i:02d}.json');r['original_methods']=r['methods'];r['methods']={}
    a,_=generate(1002+i);m=SimModel()
    for method,old in r['original_methods'].items():
        if old['solver']['success']:r['methods'][method]=old;continue
        c,U,Y,Uh=inputs(m,a,method);b,info=exact_root(m,c,U,Y,Uh,np.asarray(old['beta']))
        bt=m.truth();r['methods'][method]=dict(beta=b,solver=info,delay=delay_values(m.delay(b)),
            delay_error_us=(m.delay(b)-TRUTH_TAU)*1e6,frequency_hz=m.modal(b)[0],damping_percent=m.modal(b)[1],
            modal_relative_error_percent=100*np.expm1(b[m.iw:m.base_nb]-bt[m.iw:m.base_nb]))
    js(dest,r);print('SIMROOT',i,{k:v['solver']['success'] for k,v in r['methods'].items()},flush=True);return r
if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(one_ref,i) for i in range(40)]):f.result()
