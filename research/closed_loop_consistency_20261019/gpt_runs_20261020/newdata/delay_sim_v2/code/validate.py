"""Deterministic acceptance before power prediction or formal pilot endpoints."""
from oracle_model import *
from profile_solver import factor,coefficient,Profile

def main():
    out={};a,bb,cc,dd=plant(60,True)
    # Exact arbitrary endpoint boundary: no simulation, no random numbers.
    dx=np.sin(np.arange(len(a))+.3)*1e-4
    for band in BANDS:
        m=OracleModel(band);b=m.truth();f=np.linspace(m.lo,m.hi,241)
        c=m.context(f);k=np.arange(len(f))[:,None];i=np.arange(4)[None,:]
        U=np.exp(1j*(.071*k*k+.3*k*i+.4*i))*np.array([1,1.1,.9,1.2])
        y,J=m.jac(c,b,U);G=m.resp(c,b)
        baseG=ctl.plant(f)
        response=la.norm(G-baseG)/la.norm(baseG)
        rngdir=np.cos(np.arange(m.nb)+.27)
        rngdir[:m.nl]*=1e-10;rngdir[m.nl:2*m.nm+m.nl]*=.01;rngdir[-7:]*=.1
        h=1e-3
        fd=(m.compact(c,b+h*rngdir,U,False)[0]-m.compact(c,b-h*rngdir,U,False)[0])/(2*h)
        deriv=la.norm(fd-np.einsum('kop,p->ko',J,rngdir))/max(la.norm(fd),1e-300)
        boundary=np.array([cc[:4]@la.solve(z*np.eye(len(a))-a,z*dx) for z in c['z']])
        X=m.design(c,b[m.nl:],U,False)[1][:,:,m.it:]
        errs=[]
        for o in range(4):
            ff=factor(rr(X[:,o]));pred=ff['Q']@(ff['Q'].T@rr(boundary[:,o]))
            errs.append(la.norm(pred-rr(boundary[:,o]))/la.norm(rr(boundary[:,o])))
        out[band]=dict(parameters=m.nb,response_vs_independent_state= float(response),
            analytic_jacobian_relative=float(deriv),boundary_relative=max(errs),
            steady_output_relative=float(la.norm(y-np.einsum('koi,ki->ko',G,U))/la.norm(y)))
    out['acceptance_pass']=all(r['response_vs_independent_state']<1e-9 and r['analytic_jacobian_relative']<1e-4 and r['boundary_relative']<1e-8 for r in out.values())
    save(ROOT/'validation/model_acceptance.json',out);print(json.dumps(out,default=enc,indent=2),flush=True)
    if not out['acceptance_pass']:raise SystemExit(2)
    # Full profile derivative stability: use deterministic nonzero perturbation,
    # not a formal MC record. Same reference-IV map; no stochastic draws.
    rows=[]
    for fc in (15,60):
      for band in BANDS:
        m=OracleModel(band);b=m.truth();data=deterministic_spectra(10,0,band,fc,'M3')
        f=data['f'];data['Y']=data['Y']+1e-7*np.exp(1j*(f[:,None]*np.arange(1,5)[None,:]+.3))
        pr=Profile(m,data,'M3');s=pr.profile(b[m.nl:],b[:m.nl]);q=b[m.nl:]
        scale=np.r_[np.full(2*m.nm,.04),np.full(7,.4)];h=1e-3
        d=np.sin(np.arange(m.nq)+.9);d/=la.norm(d)
        vals=[];passed=s['inner_ok']
        for step in (h,h/2):
            sp=pr.profile(q+scale*d*step,s['beta'][:m.nl]);sm=pr.profile(q-scale*d*step,s['beta'][:m.nl])
            passed &= sp['inner_ok'] and sm['inner_ok']
            denom=np.sqrt(np.maximum(s['zz'][m.nl:]*s['ee'],1e-300))
            vals.append((sp['g'][m.nl:]-sm['g'][m.nl:])/(2*step)/denom)
        rel=la.norm(vals[0]-vals[1])/max(la.norm(vals[1]),1e-300)
        rows.append(dict(fc=fc,band=band,inner_pass=bool(passed),derivative_step_relative=float(rel),passed=bool(passed and rel<1e-4)))
        save(ROOT/'validation/profile_derivative_acceptance.json',rows)
        print(rows[-1],flush=True)
    if not all(r['passed'] for r in rows):raise SystemExit(3)

if __name__=='__main__':main()
