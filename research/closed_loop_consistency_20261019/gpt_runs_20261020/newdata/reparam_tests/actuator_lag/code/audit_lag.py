from lag_core import *
import hashlib

def main():
    hashes=rd(ROOT/'baseline/input_hashes.json')
    changed=[];missing=[]
    for p,h in hashes.items():
        path=Path(p)
        if not path.exists():missing.append(p)
        elif hashlib.sha256(path.read_bytes()).hexdigest()!=h:changed.append(p)
    rows=[]
    for rid in old.IDS:
        for kind in ('a','c'):
            p=ROOT/'raw'/f'{rid}_{kind}_joint.json'
            if not p.exists():rows.append(dict(id=rid,model=kind,present=False));continue
            r=rd(p);row=dict(id=rid,model=kind,present=True,success=r.get('success',False))
            if 'beta' in r:
                m=JointModel(kind);b=np.asarray(r['beta']);pr=JointSolver(m,get_data(rid))
                A,g,sc,ee,zz=pr.eq(b);scale=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
                sv=la.svdvals(A*scale[:,None]*scale[None,:])
                row.update(score_recomputed=sc,score_saved=r['solver']['score'],scaled_condition=float(sv[0]/sv[-1]),
                           scaled_rank_1e12=int(sum(sv>sv[0]*1e-12)),parameters=m.nb,
                           cutoff_near_technical_upper_1percent=(np.asarray(r['cutoff_hz'])>198000).tolist(),
                           joint_Q_sum_error=float(sum(r['qe'][w]['Q'] for w in WINDOWS)-r['qe']['joint']['Q']),
                           joint_EQ_sum_error=float(sum(r['qe'][w]['EQ'] for w in WINDOWS)-r['qe']['joint']['EQ']))
            rows.append(row)
    result=dict(original_files=len(hashes),changed=changed,missing=missing,original_unchanged=not(changed or missing),
                endpoints=rows,required=12,completed=sum(x['present'] for x in rows),
                successful=sum(x.get('success',False) for x in rows),
                all_results_present=all(x['present'] for x in rows),
                verification_passed=all(r['passed'] for r in rd(ROOT/'analysis/verification.json')))
    js(ROOT/'analysis/audit.json',result)
    print({k:v for k,v in result.items() if k!='endpoints'},flush=True)
    # Independent parametric check of the joint residual expectation; no fitted data resampled.
    rid='A_nominal';kind='c';r=rd(ROOT/'raw'/f'{rid}_{kind}_joint.json')
    m=JointModel(kind);b=np.asarray(r['beta']);pr=JointSolver(m,get_data(rid))
    pieces=[];ranges={};start=0
    for w,a in pr.data.items():
        _,J=m.output(w,pr.ctx[w],b,a['U']);_,Z=m.output(w,pr.ctx[w],b,pr.uh_window(w,b))
        Jr=old.realvec(J);Zr=old.realvec(Z)
        L=np.linalg.cholesky(.5*old.realmat(a['N']*a['Phi']))
        B=np.linalg.solve(L,Jr).reshape(-1,m.nb)
        ZL=np.einsum('kop,kol->pkl',Zr,L,optimize=True).reshape(m.nb,-1)
        A=np.einsum('kop,koq->pq',Zr,Jr,optimize=True)
        pieces.append((B,ZL,A));ranges[w]=slice(start,start+B.shape[0]);start+=B.shape[0]
    B=np.concatenate([p[0] for p in pieces]);ZL=np.concatenate([p[1] for p in pieces],axis=1);A=sum(p[2] for p in pieces)
    sc=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
    H=sc[:,None]*la.lstsq(A*sc[:,None]*sc[None,:],sc[:,None]*ZL,cond=1e-13,lapack_driver='gelsy')[0]
    rng=np.random.default_rng(2026102071);v=rng.standard_normal((B.shape[0],200));res=v-B@(H@v)
    check={}
    for w,s in ranges.items():
        q=np.sum(res[s]**2,axis=0);mean=q.mean();se=q.std(ddof=1)/np.sqrt(200);expected=r['qe'][w]['EQ']
        # Direct low-rank trace check uses all H, including the other window.
        eq2=(s.stop-s.start)-2*np.einsum('ip,pi->',B[s],H[:,s])+np.einsum('ij,ji->',B[s].T@B[s],H@H.T)
        check[w]=dict(analytic=expected,independent_low_rank_trace=float(eq2),bootstrap_mean=float(mean),
                      bootstrap_MC_SE=float(se),difference_in_MC_SE=float((mean-expected)/se),
                      trace_relative_error=float(abs(eq2-expected)/expected))
    js(ROOT/'analysis/residual_expectation_bootstrap.json',dict(id=rid,model=kind,draws=200,windows=check))
    print(check,flush=True)
if __name__=='__main__':main()
