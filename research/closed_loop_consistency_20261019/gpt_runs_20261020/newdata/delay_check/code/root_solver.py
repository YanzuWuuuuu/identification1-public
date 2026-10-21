"""Full IV moment root: checks GN failures without changing the moment equation."""
from core import *
from scipy.optimize import root
def moments_only(m,c,U,Y,Uh,b):
    yh,_,_=m.compact(c,b,U);uh=Uh(b) if callable(Uh) else Uh
    _,zl,zn=m.compact(c,b,uh)
    e=Y-yh;g=np.zeros(m.nb);zz=np.zeros(m.nb)
    for o in range(4):
        oi=slice(o*m.L,(o+1)*m.L);g[oi]=(zl[:,o].conj().T@e[:,o]).real
        zz[oi]=np.sum(abs(zl[:,o])**2,axis=0)
    g[m.iw:]=np.einsum('kop,ko->p',zn.conj(),e).real
    zz[m.iw:]=np.sum(abs(zn)**2,axis=(0,1));ee=np.sum(abs(e)**2)
    return g,float(np.max(abs(g)/np.sqrt(np.maximum(zz*ee,1e-300))))
def exact_root(m,c,U,Y,Uh,b0,maxfev=4000):
    b0=b0.copy();A,g,score,ee=equations(m,c,U,Y,Uh(b0) if callable(Uh) else Uh,b0)
    scale=1/np.sqrt(np.maximum(abs(np.diag(A)),1e-300))
    # Freeze variable and equation scaling. Numerical Jacobian differentiates
    # Z(beta), recovered input, and residual, not just Z^H J.
    evaluations=0;exceptions=[]
    def fun(x):
        nonlocal evaluations
        evaluations+=1;b=b0+scale*x
        if np.any(abs(b[m.iw:m.base_nb])>30) or not np.isfinite(b).all():
            return np.ones(m.nb)*1e10
        try:return scale*moments_only(m,c,U,Y,Uh,b)[0]/max(np.sqrt(ee),1e-300)
        except Exception as e:
            exceptions.append(repr(e));return np.ones(m.nb)*1e10
    out=root(fun,np.zeros(m.nb),method='hybr',options=dict(xtol=1e-9,maxfev=maxfev,factor=.1,eps=1e-10))
    b=b0+scale*out.x;_,score=moments_only(m,c,U,Y,Uh,b)
    return b,dict(success=bool(score<1e-6),score=score,status=str(out.message),scipy_success=bool(out.success),
                  evaluations=evaluations,exceptions=exceptions,solver='full IV moment hybrid root, frozen scaling')
