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
    def equation(self,eta,chunk=1600):
        m=self.m;d=self.d;b=self.bt+self.transform@eta
        gt=np.zeros(m.nb);gp=gt.copy();zz=gt.copy();ee=0.;self.calls+=1
        for a in range(0,len(d['f']),chunk):
            sl=slice(a,min(a+chunk,len(d['f'])));c=m.context(d['f'][sl]);U=self.U[sl]
            Yh,X,N=m.compact(c,b,U,True);G=m.resp(c,b);diff=d['G'][sl]-G
            Em=np.einsum('koi,ki->ko',d['G'][sl],U)-Yh
            En=self.Lv[sl]+diff@self.Un[sl]
            q=b[m.nl:];l=b[:m.nl].reshape(4,m.L);hs=m.blocks(c,q,True)
            E=np.exp(-c['s'][:,None,None]*m.tau(q))
            basis=np.stack([hs[i][0] for i in m.active]+[np.ones(len(U)),(c['s']/m.base.wr)**2],axis=1)
            EU=E[:,:,:,None]*self.Un[sl][:,None,:,:]
            Xnoise=(basis[:,None,:,None,None]*EU[:,:,None,:,:]).transpose(0,1,4,2,3).reshape(len(U),4,4,m.it)
            Gq=np.zeros((len(U),4,4,m.nq),complex)
            for j,i in enumerate(m.active):
                for p in range(2):Gq[:,:,:,p*m.nm+j]=hs[i][3][p][0][:,None,None]*l[:,4*j:4*j+4]*E
            for p in range(7):Gq[:,:,:,2*m.nm+p]=-c['s'][:,None,None]*UNIT*G*DESIGN[:,:,p]
            Nnoise=np.einsum('koip,kia->koap',Gq,self.Un[sl],optimize=True)
            gn=np.einsum('kop,ko->p',N.conj(),Em,optimize=True).real+np.einsum('koap,koa->p',Nnoise.conj(),En,optimize=True).real
            gp[m.nl:]+=gn;gt[m.nl:]+=gn
            zz[m.nl:]+=np.sum(abs(N)**2,axis=(0,1))+np.sum(abs(Nnoise)**2,axis=(0,1,2))
            for o in range(4):
                oi=slice(o*m.L,(o+1)*m.L);ids=self.ti[o*15:(o+1)*15]
                g=(X[:,o].conj().T@Em[:,o]).real
                g[:m.it]+=np.einsum('kal,ka->l',Xnoise[:,o].conj(),En[:,o],optimize=True).real
                gp[oi]+=g;gt[oi]+=g
                gt[ids]-=g[m.it:]
                gt[ids]+=((X[:,o,m.it:]@self.transform[np.ix_(ids,ids)]).conj().T@Em[:,o]).real
                zz[oi]+=np.sum(abs(X[:,o])**2,axis=0)
                zz[o*m.L:o*m.L+m.it]+=np.sum(abs(Xnoise[:,o])**2,axis=(0,1))
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
    A=(M-C)*scale[:,None]*scale[None,:]
    eta=xc.copy();history=[];status='iteration_limit';state=pop.equation(eta)
    # An auxiliary deterministic population-root solver. Not an MC trial and
    # not a replacement for the approved profile trust-region MC solver.
    for it in range(100):
        if time.perf_counter()-t0>(2700 if m.band=='W' else 900):status='wall_time_limit';break
        merit=la.norm(scale*state['g'])
        row=dict(iteration=it,score=state['score'],inner_score=state['inner_score'],scaled_moment_norm=float(merit),
            beta=state['beta'],seconds=time.perf_counter()-t0);history.append(row)
        save(folder/f'r{r:03d}_progress.json',row)
        if state['score']<1e-6 and state['inner_score']<1e-9:status='converged';break
        step=scale*la.solve(A,scale*state['g'])
        accepted=False
        for alpha in (1.,.5,.25,.125,.0625,.03125):
            proposal=eta+alpha*step;q=pop.bt[m.nl:]+proposal[m.nl:]
            if np.any(q<=m.lower) or np.any(q>=m.upper):continue
            trial=pop.equation(proposal)
            if la.norm(scale*trial['g'])<merit:
                dx=(proposal-eta)/scale;dy=scale*(state['g']-trial['g'])
                if dx@dx>1e-30:A+=np.outer(dy-A@dx,dx)/(dx@dx)
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
        if (folder/f'r{r:03d}.json').exists():
            old=json.loads((folder/f'r{r:03d}.json').read_text())
            if old['success']:continue
            # Preserve unsuccessful auxiliary solves; never discard them.
            archive=folder/'pre_broyden';archive.mkdir(exist_ok=True)
            import shutil
            shutil.copy2(folder/f'r{r:03d}.json',archive/f'r{r:03d}.json')
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
