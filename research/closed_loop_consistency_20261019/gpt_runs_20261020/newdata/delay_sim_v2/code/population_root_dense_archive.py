"""Auxiliary root of the complete direct-method population moment.

This integrates the specified stationary DFT second moments exactly on the
fitting grid. No noise realization is drawn. It is not a claim that finite-record
DFT noise bins are exactly independent. All nuisance columns remain present.
The chord/Newton linear algebra is in the SAME nonsingular coordinates as power
calculations; acceptance uses the original full moment, not its reduced block.
"""
from population_math import *

class Population:
    def __init__(self,P,band,fc,grid_method,common):
        self.m=OracleModel(band);self.bt=self.m.truth();self.d=deterministic_spectra(P,0,band,fc,grid_method)
        self.Lv=np.linalg.cholesky(self.d['Pv']);self.Un=-self.d['F']@self.Lv
        self.transform=common['transform'];self.tf,self.ti=indices(self.m)
        self.calls=0;self.start=time.perf_counter()
    def set_phase(self,r):
        dp,_=excitation(r);d=self.d;P=d['P'];D=np.zeros_like(d['D']);ix=d['k']%P==0
        D[ix]=P*dp[d['k'][ix]//P];self.U=np.einsum('koi,ki->ko',d['S'],D)
    def equation(self,eta,chunk=640):
        m=self.m;d=self.d;b=self.bt+self.transform@eta
        gt=np.zeros(m.nb);gp=gt.copy();zz=gt.copy();ee=0.;self.calls+=1
        for a in range(0,len(d['f']),chunk):
            sl=slice(a,min(a+chunk,len(d['f'])));c=m.context(d['f'][sl]);U=self.U[sl]
            J=transfer_jac(m,c,b);G=m.resp(c,b);diff=d['G'][sl]-G
            Yn,X,N=m.compact(c,b,np.zeros_like(U),True)
            Em=np.einsum('koi,ki->ko',diff,U)-Yn
            En=self.Lv[sl]+diff@self.Un[sl]
            Jm=np.einsum('koip,ki->kop',J,U,optimize=True);Jm[:,:,m.nl:]+=N
            for o in range(4):Jm[:,o,o*m.L+m.it:(o+1)*m.L]=X[:,o,m.it:]
            Jn=np.einsum('koip,kia->koap',J,self.Un[sl],optimize=True)
            gp+=np.einsum('kop,ko->p',Jm.conj(),Em,optimize=True).real
            ng=np.einsum('koap,koa->p',Jn.conj(),En,optimize=True).real
            gp+=ng;zz+=np.sum(abs(Jm)**2,axis=(0,1))+np.sum(abs(Jn)**2,axis=(0,1,2))
            for o in range(4):
                ids=self.ti[o*15:(o+1)*15];Jm[:,o,ids]=Jm[:,o,ids]@self.transform[np.ix_(ids,ids)]
            gt+=np.einsum('kop,ko->p',Jm.conj(),Em,optimize=True).real+ng
            ee+=float(np.sum(abs(Em)**2)+np.sum(abs(En)**2))
        scores=abs(gp)/np.sqrt(np.maximum(zz*ee,1e-300))
        return dict(g=gt,original_g=gp,score=float(max(scores)),inner_score=float(max(scores[:m.nl])),
            ee=ee,zz=zz,beta=b)

def run_one(pop,folder,r,common,power):
    m=pop.m;pop.set_phase(r);t0=time.perf_counter();startcalls=pop.calls
    # Direct moments on the same grid are also reported for reference-IV designs.
    M=power['M'].copy();C=common['C'];bv=common['b']
    if pop.d['method']=='M3':M+=common['mn']
    x,diag,scale=stable_solve(M,bv);xc,dc=solve_general(M-C,bv,scale)
    A=(M-C)*scale[:,None]*scale[None,:];lu=la.lu_factor(A)
    eta=xc.copy();history=[];status='iteration_limit';state=pop.equation(eta)
    # An auxiliary deterministic population-root solver. Not an MC trial and
    # not a replacement for the approved profile trust-region MC solver.
    for it in range(30):
        if time.perf_counter()-t0>(2700 if m.band=='W' else 900):status='wall_time_limit';break
        merit=la.norm(scale*state['g'])
        row=dict(iteration=it,score=state['score'],inner_score=state['inner_score'],scaled_moment_norm=float(merit),
            beta=state['beta'],seconds=time.perf_counter()-t0);history.append(row)
        save(folder/f'r{r:03d}_progress.json',row)
        if state['score']<1e-6 and state['inner_score']<1e-9:status='converged';break
        step=scale*la.lu_solve(lu,scale*state['g'])
        accepted=False
        for alpha in (1.,.5,.25,.125,.0625,.03125):
            proposal=eta+alpha*step;q=pop.bt[m.nl:]+proposal[m.nl:]
            if np.any(q<=m.lower) or np.any(q>=m.upper):continue
            trial=pop.equation(proposal)
            if la.norm(scale*trial['g'])<merit:
                eta=proposal;state=trial;accepted=True;break
        if not accepted:status='no_decreasing_population_step';break
    q=state['beta'][m.nl:];bound=bool(np.min(np.minimum(q-m.lower,m.upper-q))<1e-5)
    success=state['score']<1e-6 and state['inner_score']<1e-9 and not bound
    row=dict(r=r,success=bool(success),status=status,score=state['score'],inner_score=state['inner_score'],near_bound=bound,
        beta=state['beta'],bias=eta,first_order=x,curvature_corrected=xc,
        seconds=time.perf_counter()-t0,calls=pop.calls-startcalls,history=history,
        first_order_diagnostic=diag,curvature_diagnostic=dc)
    save(folder/f'r{r:03d}.json',row);return row

def run(P,band,fc,grid_method,rstart=0,rstop=30):
    base=ROOT/'power'/f'FC{fc}_{band}_P{P}_{grid_method}'
    folder=ROOT/'population'/f'FC{fc}_{band}_P{P}_{grid_method}';folder.mkdir(parents=True,exist_ok=True)
    with np.load(base/'common.npz') as a:common={k:a[k] for k in a.files}
    pop=Population(P,band,fc,grid_method,common)
    for r in range(rstart,rstop):
        if (folder/f'r{r:03d}.json').exists():continue
        with np.load(base/f'r{r:03d}.npz') as a:power={k:a[k] for k in a.files}
        try:
            result=run_one(pop,folder,r,common,power)
            print(fc,band,P,grid_method,r,result['status'],result['score'],result['seconds'],flush=True)
        except Exception as exc:
            save(folder/f'r{r:03d}.json',dict(r=r,success=False,error=repr(exc)))
            print(fc,band,P,grid_method,r,repr(exc),flush=True)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--P',type=int,required=True);ap.add_argument('--band',required=True)
    ap.add_argument('--fc',type=int,required=True);ap.add_argument('--grid-method',required=True)
    ap.add_argument('--rstart',type=int,default=0);ap.add_argument('--rstop',type=int,default=30)
    a=ap.parse_args();run(a.P,a.band,a.fc,a.grid_method,a.rstart,a.rstop)
