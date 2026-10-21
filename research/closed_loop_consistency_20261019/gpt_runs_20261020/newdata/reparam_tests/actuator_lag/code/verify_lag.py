from lag_core import *

def main():
    results=[]
    for kind in ('a','c'):
        m=JointModel(kind);b=m.initial('A_nominal');data=get_data('A_nominal')
        errors=[];rng=np.random.default_rng(716)
        for w,a in data.items():
            select=np.arange(0,len(a['f']),31);c=m.full[w].context(a['f'][select]);U=a['U'][select]
            y,J=m.output(w,c,b,U)
            for _ in range(3):
                scale=np.maximum(abs(b),1e-10);scale[m.nl:]=1
                d=rng.normal(size=m.nb)*scale;eps=1e-6
                yp=m.output(w,c,b+eps*d,U,False)[0];ym=m.output(w,c,b-eps*d,U,False)[0]
                fd=(yp-ym)/(2*eps);an=np.einsum('kop,p->ko',J,d)
                errors.append(float(np.linalg.norm(fd-an)/np.linalg.norm(an)))
            lim=b.copy();lim[m.lag_slice]=np.log(2*np.pi*1e16)
            full=m.local(w,lim);y0=m.full[w].base.output(c,full[:-4],U,False)[0]
            limit=float(np.linalg.norm(m.output(w,c,lim,U,False)[0]-y0)/np.linalg.norm(y0))
            assert limit<1e-10
        pr=JointSolver(m,data);start=time.perf_counter()
        bp,ok,_,iterations=pr.profile(b[m.nl:],b[:m.nl])
        A,g,score,ee,zz=pr.eq(bp)
        directA=np.zeros_like(A);directg=np.zeros_like(g);directz=np.zeros_like(zz)
        for w,a in data.items():
            y,J=m.output(w,pr.ctx[w],bp,a['U'])
            _,Z=m.output(w,pr.ctx[w],bp,pr.uh_window(w,bp))
            directA+=np.einsum('kop,koq->pq',Z.conj(),J,optimize=True).real
            directg+=np.einsum('kop,ko->p',Z.conj(),a['Y']-y,optimize=True).real
            directz+=np.sum(abs(Z)**2,axis=(0,1))
        errA=float(np.linalg.norm(A-directA)/np.linalg.norm(directA))
        errg=float(np.linalg.norm(g-directg)/max(np.linalg.norm(directg),1e-300))
        errz=float(np.linalg.norm(zz-directz)/np.linalg.norm(directz))
        result=dict(kind=kind,parameter_count=m.nb,finite_difference_step=eps,derivative_relative_errors=errors,
            zero_lag_limit_relative_error=limit,moment_matrix_relative_error=errA,
            moment_vector_relative_error=errg,instrument_norm_relative_error=errz,
            profile_inner_success=ok,profile_inner_iterations=iterations,
            profile_linear_score=float(np.max(abs(g[:m.nl])/np.sqrt(zz[:m.nl]*ee))),
            full_score=score,seconds=time.perf_counter()-start)
        result['passed']=max(errors)<1e-6 and limit<1e-10 and errA<1e-10 and errz<1e-10 and ok
        results.append(result)
    js(ROOT/'analysis/verification.json',results)
    print(results,flush=True)
    assert all(r['passed'] for r in results)
if __name__=='__main__':main()
