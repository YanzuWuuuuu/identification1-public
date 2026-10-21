"""Complete profiled moment derivative, including changing reference tools.

Hessian-vector contractions use central partial derivatives at fixed linear
coefficients and fixed Uh. The implicit linear solve adds dl/dq, and the
explicit instrument derivative adds dUh/dbeta. No frozen-tool approximation.
"""
from oracle_model import *
from profile_solver import coefficient
from population_math import transfer_jac

def derivative(pr,state):
    m=pr.m;b=state['beta'];n=len(pr.U);p=m.nb;nl=m.nl;nq=m.nq
    Uh=pr.uh(b);er=np.array(state['es']);e=(er[:,:n]+1j*er[:,n:]).T
    Rinv=[np.column_stack([coefficient(f,np.eye(m.L)[:,j]) for j in range(m.L)]) for f in state['fx']]
    H=np.zeros((p,p));h=1e-4
    for j in range(nq):
        bp=b.copy();bm=b.copy();bp[nl+j]+=h;bm[nl+j]-=h
        _,Xp,Np=m.compact(pr.c,bp,Uh,True,pr.mask)
        _,Xm,Nm=m.compact(pr.c,bm,Uh,True,pr.mask)
        for o in range(4):
            hh=(rr((Xp[:,o]-Xm[:,o])/(2*h)).T@state['es'][o])
            H[o*m.L:(o+1)*m.L,nl+j]=Rinv[o].T@hh
        H[nl:,nl+j]=sum((rr((Np[:,o]-Nm[:,o])/(2*h)).T@state['es'][o] for o in range(4)),start=np.zeros(nq))
    H[nl:,:nl]=H[:nl,nl:].T
    H[nl:,nl:]=(H[nl:,nl:]+H[nl:,nl:].T)/2
    A=H.copy()
    for o in range(4):
        J=np.zeros((2*n,p));Z=J.copy();oi=slice(o*m.L,(o+1)*m.L)
        J[:,oi]=state['fx'][o]['Q'];J[:,nl:]=rr(state['JN'][:,o])
        Z[:,oi]=rr(state['ZL'][:,o])@Rinv[o];Z[:,nl:]=rr(state['ZN'][:,o])
        A-=Z.T@J
    if pr.method!='M1':
        G=m.resp(pr.c,b);F=np.linalg.solve(np.eye(4)+pr.K@G,pr.K)
        # Sparse output blocks avoid a dense p-by-p change of coordinates.
        Gj=transfer_jac(m,pr.c,b)
        for o in range(4):
            oi=slice(o*m.L,(o+1)*m.L)
            Gj[:,o,:,oi]=Gj[:,o,:,oi]@Rinv[o]
        Vu=np.einsum('koip,ko->kip',Gj.conj(),e,optimize=True)
        Gu=np.einsum('koip,ki->kop',Gj,Uh,optimize=True)
        dUh=-np.einsum('koi,kip->kop',F,Gu,optimize=True)
        A+=(Vu.reshape(-1,p).T@dUh.reshape(-1,p).conj()).real
    # Differentiating g_l(l(q),q)=0 gives dl/dq=-A_ll^-1 A_lq.
    return A[nl:,nl:]-A[nl:,:nl]@la.solve(A[:nl,:nl],A[:nl,nl:],check_finite=False)

def validate():
    from profile_solver import Profile
    rows=[]
    for method in ('M1','M3'):
      for fc in (15,60):
       for band in BANDS:
        m=OracleModel(band);b=m.truth()
        # Validation uses exact excited-line spectra, no noise simulation.
        d=deterministic_spectra(10,0,band,fc,'M3');f=d['f']
        d['Y']=d['Y']+1e-7*np.exp(1j*(f[:,None]*np.arange(1,5)[None,:]+.3))
        pr=Profile(m,d,method);s=pr.profile(b[m.nl:],b[:m.nl]);A=derivative(pr,s)
        direction=np.sin(np.arange(m.nq)+.9);direction/=la.norm(direction)
        scale=np.r_[np.full(2*m.nm,.04),np.full(7,.4)];delta=scale*direction;h=5e-4
        sp=pr.profile(b[m.nl:]+h*delta,s['beta'][:m.nl]);sm=pr.profile(b[m.nl:]-h*delta,s['beta'][:m.nl])
        fd=(sp['g'][m.nl:]-sm['g'][m.nl:])/(2*h)
        norm=np.sqrt(s['zz'][m.nl:]*s['ee']);err=la.norm((A@delta-fd)/norm)/la.norm(fd/norm)
        row=dict(method=method,fc=fc,band=band,relative=float(err),passed=bool(err<1e-4 and s['inner_ok'] and sp['inner_ok'] and sm['inner_ok']))
        rows.append(row);save(ROOT/'validation/implicit_derivative_acceptance.json',rows);print(row,flush=True)
    if not all(r['passed'] for r in rows):raise SystemExit(2)

if __name__=='__main__':validate()
