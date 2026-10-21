"""Conditional stationary-DFT population moments; no random noise generation.

Coordinates retain every approved parameter. Ill-conditioned endpoint columns
are orthogonalized BEFORE forming a Gram matrix. Q includes the quadratic
Gaussian term for direct tools, not only the usual deterministic-tool term.
"""
from oracle_model import *
from profile_solver import factor,coefficient

def herm(a): return (a+a.conj().swapaxes(-1,-2))/2
def gram(a):
    a=a.reshape(-1,a.shape[-1]);return (a.conj().T@a).real

def transfer_jac(m,c,b):
    q=b[m.nl:];l=b[:m.nl].reshape(4,m.L);hs=m.blocks(c,q,True)
    E=np.exp(-c['s'][:,None,None]*m.tau(q));n=len(c['f'])
    J=np.zeros((n,4,4,m.nb),complex)
    for j,i in enumerate(m.active):
        for o in range(4):
            for a in range(4):J[:,o,a,o*m.L+4*j+a]=hs[i][0]*E[:,o,a]
        for p in range(2):J[:,:,:,m.nl+p*m.nm+j]=hs[i][3][p][0][:,None,None]*l[:,4*j:4*j+4]*E
    for o in range(4):
        for a in range(4):
            J[:,o,a,o*m.L+m.id0+a]=E[:,o,a]
            J[:,o,a,o*m.L+m.id2+a]=(c['s']/m.base.wr)**2*E[:,o,a]
    G=m.resp(c,b)
    for p in range(7):J[:,:,:,m.nb-7+p]=-c['s'][:,None,None]*UNIT*DESIGN[:,:,p]*G
    return J

def indices(m):
    ti=np.concatenate([np.arange(o*m.L+m.it,(o+1)*m.L) for o in range(4)])
    return np.setdiff1d(np.arange(m.nb),ti),ti

def noise_bias(J,F,Pv):
    return -np.einsum('koip,kia,kao->p',J,F,Pv,optimize=True).real

def stable_solve(M,b):
    scale=1/np.sqrt(np.maximum(np.diag(M),1e-300));A=M*scale[:,None]*scale[None,:]
    w=la.eigvalsh(A);condition=float(w[-1]/w[0])
    if w[0]<=1e-13*w[-1]:raise la.LinAlgError('equilibrated population matrix lacks numerical rank')
    fac=la.cho_factor(A,check_finite=False)
    x=scale*la.cho_solve(fac,scale*b,check_finite=False)
    # Independent full-rank symmetric factorization, with no truncation.
    x2=scale*la.solve(A,scale*b,assume_a='sym',check_finite=False)
    residual=la.norm(A@(x/scale)-scale*b)/max(la.norm(scale*b),1e-300)
    difference=la.norm((x-x2)/scale)/max(la.norm(x/scale),1e-300)
    return x,dict(condition=condition,residual=float(residual),factorization_difference=float(difference)),scale

def solve_general(A,b,scale):
    B=A*scale[:,None]*scale[None,:];sv=la.svdvals(B)
    if sv[-1]<=sv[0]*1e-13:raise la.LinAlgError('curvature-corrected population matrix lacks numerical rank')
    x=scale*la.solve(B,scale*b,check_finite=False)
    return x,dict(condition=float(sv[0]/sv[-1]),residual=float(la.norm(B@(x/scale)-scale*b)/max(la.norm(scale*b),1e-300)))

def sandwich(A,Q,scale):
    B=A*scale[:,None]*scale[None,:];Q=Q*scale[:,None]*scale[None,:]
    sol=la.solve(B,Q,check_finite=False)
    V=la.solve(B,sol.T,check_finite=False).T
    return (V+V.T)/2*scale[:,None]*scale[None,:]

def prepare(P,band,fc,method,chunk=256):
    """Noise-only and endpoint blocks, identical for all 30 phase realizations."""
    t0=time.perf_counter();m=OracleModel(band);b=m.truth()
    d=deterministic_spectra(P,0,band,fc,method);c=m.context(d['f']);nf=len(d['f']);p=m.nb
    tf,ti=indices(m);X=m.design(c,b[m.nl:],np.zeros_like(d['U']),False)[1][:,0,m.it:]
    ff=factor(rr(X));qt=ff['Q'][:nf]+1j*ff['Q'][nf:]
    transform=np.eye(p);Rinv=np.column_stack([coefficient(ff,np.eye(15)[:,i]) for i in range(15)])
    for o in range(4):transform[np.ix_(ti[o*15:(o+1)*15],ti[o*15:(o+1)*15])]=Rinv
    mn=np.zeros((p,p));q2=np.zeros((p,p));qtt=np.zeros((p,p));bv=np.zeros(p);C=np.zeros((p,p))
    bt=m.truth();h=1e-4;bc=[(bt.copy(),bt.copy()) for _ in range(m.nq)]
    for j,(bp,bm) in enumerate(bc):bp[m.nl+j]+=h;bm[m.nl+j]-=h
    common_delay_formula=0.
    for a in range(0,nf,chunk):
        sl=slice(a,min(a+chunk,nf));cc=m.context(d['f'][sl]);F=d['F'][sl];Pv=d['Pv'][sl]
        Lv=np.linalg.cholesky(Pv);Jn=transfer_jac(m,cc,bt)
        vj=-np.einsum('koip,kia,kab->kobp',Jn,F,Lv,optimize=True)
        mn+=gram(vj);bv+=noise_bias(Jn,F,Pv)
        # H_i=Herm(Lv^H [-G_i F] Lv), Cov(v^H H_i v,v^H H_j v)=tr(H_i H_j).
        B=np.einsum('koa,kobp->kabp',Lv.conj(),vj,optimize=True)
        B=(B+B.conj().swapaxes(1,2))/2;q2+=gram(B)
        jt=np.zeros((len(cc['f']),4,p),complex)
        for o in range(4):jt[:,o,ti[o*15:(o+1)*15]]=qt[sl]
        qtt+=.5*gram(np.einsum('koa,kop->kap',Lv.conj(),jt,optimize=True))
        for j,(bp,bm) in enumerate(bc):
            C[:,m.nl+j]+=(noise_bias(transfer_jac(m,cc,bp),F,Pv)-noise_bias(transfer_jac(m,cc,bm),F,Pv))/(2*h)
        common_delay_formula+=float(np.sum(-2*np.pi*cc['f']*UNIT*np.einsum('koi,kij,kjo->k',d['G'][sl],F,Pv,optimize=True).imag))
    C[m.nl:,:m.nl]=C[:m.nl,m.nl:].T
    asym=la.norm(C-C.T)/max(la.norm(C),1e-300);C=(C+C.T)/2
    delaycheck=abs(common_delay_formula-bv[-7])/max(abs(bv[-7]),1e-300)
    assert delaycheck<1e-8
    line=d['k']%P==0
    # The conditional mean is nonzero only on these fixed excited lines.
    Jline=transfer_jac(m,m.context(d['f'][line]),bt)
    obj=dict(P=P,band=band,fc=fc,method=method,f=d['f'],k=d['k'],S=d['S'][line],
        Lv=np.linalg.cholesky(d['Pv'][line]),Jline=Jline,qt=qt[line],
        mn=mn,q2=q2,qtt=qtt,b=bv,C=C,transform=transform,
        tf=tf,ti=ti,noise_curvature_asymmetry=asym,delay_identity_relative=delaycheck,
        seconds=time.perf_counter()-t0)
    return m,obj

def realization(m,o,r):
    P=o['P'];dp,_=excitation(r);kl=o['k'][o['k']%P==0]//P;D=P*dp[kl]
    U=np.einsum('koi,ki->ko',o['S'],D);J=np.einsum('koip,ki->kop',o['Jline'],U,optimize=True)
    for a in range(4):J[:,a,o['ti'][a*15:(a+1)*15]]=o['qt']
    M=gram(J);Q=.5*gram(np.einsum('koa,kop->kap',o['Lv'].conj(),J,optimize=True))
    ti=o['ti'];M[np.ix_(ti,ti)]=np.eye(len(ti));Q[np.ix_(ti,ti)]=o['qtt'][np.ix_(ti,ti)]
    if o['method']=='M1':M+=o['mn'];Q+=o['q2'];bv=o['b'];C=o['C']
    else:bv=np.zeros(m.nb);C=np.zeros_like(M)
    x,diag,sc=stable_solve(M,bv)
    xc,dc=solve_general(M-C,bv,sc)
    V=sandwich(M-C,Q,sc)
    # Retain the complete decomposition of the high-frequency modal block.
    propagation=np.full(3,np.nan);direct=propagation.copy();schur={}
    if m.band=='W':
        hi=m.nl+m.nm+np.array([3,4,5]);keep=np.r_[hi,np.arange(m.nb-7,m.nb)]
        nuisance=np.setdiff1d(np.arange(m.nb),keep)
        an=M[np.ix_(nuisance,nuisance)];cross=M[np.ix_(nuisance,keep)]
        sn=sc[nuisance];aq=an*sn[:,None]*sn[None,:]
        sol=sn[:,None]*la.solve(aq,sn[:,None]*cross,assume_a='pos')
        sb=sn*la.solve(aq,sn*bv[nuisance],assume_a='pos')
        A=M[np.ix_(keep,keep)]-cross.T@sol;bb=bv[keep]-cross.T@sb
        direct=la.solve(A[:3,:3],bb[:3]);propagation=-la.solve(A[:3,:3],A[:3,3:]@x[-7:])
        schur=dict(high_direct_log=direct,high_delay_log=propagation,
            high_total_log=x[hi],decomposition_error=float(la.norm(direct+propagation-x[hi])))
    return dict(r=r,bias=x,bias_curvature=xc,covariance=V,diagnostic=diag,curvature_diagnostic=dc,
                schur=schur,M=M,Q=Q,C=C,b=bv,scales=sc)

def run(P,band,fc,method):
    folder=ROOT/'power'/f'FC{fc}_{band}_P{P}_{method}';folder.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();m,o=prepare(P,band,fc,method)
    np.savez_compressed(folder/'common.npz',**{k:v for k,v in o.items() if k!='Jline'})
    save(folder/'common.json',{k:o[k] for k in ('noise_curvature_asymmetry','delay_identity_relative','seconds')})
    print('prepared',fc,band,P,method,o['seconds'],flush=True)
    rows=[]
    for r in range(30):
        try:
            a=realization(m,o,r)
            np.savez_compressed(folder/f'r{r:03d}.npz',**{k:a[k] for k in ('bias','bias_curvature','covariance','M','Q','C','b','scales')})
            summary={k:v for k,v in a.items() if k not in ('covariance','M','Q','C','b','scales')}
            summary.update(success=True,standard_deviation=np.sqrt(np.maximum(np.diag(a['covariance']),0)))
        except Exception as exc:summary=dict(r=r,success=False,error=repr(exc))
        rows.append(summary);save(folder/f'r{r:03d}.json',summary)
        print(fc,band,P,method,r,summary['success'],flush=True)
    save(folder/'summary.json',dict(rows=rows,parameters=m.nb,seconds=time.perf_counter()-start))

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--P',type=int,required=True);ap.add_argument('--band',required=True)
    ap.add_argument('--fc',type=int,required=True);ap.add_argument('--method',required=True)
    a=ap.parse_args();run(a.P,a.band,a.fc,a.method)
