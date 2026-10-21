from core import *
from run_real import identifiability
from concurrent.futures import ProcessPoolExecutor,as_completed
def refine(p):
    dest=ROOT/'refinements'/p.name
    if dest.exists():return rd(dest)
    r=rd(p)
    if r['success'] or r['window']=='L':return
    t=time.perf_counter();r['original_failure']=True;r['original_solver']=r.get('solver');r['parent_file']=str(p)
    rid,w,method=r['id'],r['window'],r['method'];joint='joint_' in p.name;omit=r.get('omit')
    try:
        s,m,lb,ub=setup(w,r.get('shift_us',0),joint=joint);a=spectrum(rid,w,omit)
        c,U,Y,Uh=inputs(m,a,method);b,info=solve(m,c,U,Y,Uh,np.asarray(r['beta']),lb,ub,maxit=240,relaxation=.5)
        f,z=m.modal(b);r.update(beta=b,solver=info,success=info['success'],frequency_hz=f,damping_percent=z,
                              refinement='same IV equation, half-step relaxation, at most 240 additional iterations')
        if joint:
            r['delay']=delay_values(m.delay(b));r['identifiability']=identifiability(m,c,b,U,Uh)
            if omit is None and method=='M3':
                _,mf,_,_=setup(w);ext=spectrum(rid,w,extended=True);base,_=basefit(rid,w,'M3')
                Phi=local_phi(mf,ext,base,a['f']);r['qe_fixed_noise']=qe(m,a,b,Phi)
                r['qe_reestimated_noise']=qe(m,a,b,local_phi(m,ext,b,a['f']))
        else:
            base,_=basefit(rid,w,method);f0,z0=m.modal(base)
            r.update(delta_frequency_hz=f-f0,delta_damping_pp=z-z0)
    except Exception as e:r.update(success=False,refinement_exception=repr(e))
    r['refinement_seconds']=time.perf_counter()-t;js(dest,r)
    print('REFINE',p.name,r['success'],r.get('solver',{}).get('score'),flush=True)
    return r
if __name__=='__main__':
    paths=[p for p in (ROOT/'raw').glob('*.json') if p.name.startswith(('joint_','fixed_'))]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(refine,p) for p in paths]):f.result()
