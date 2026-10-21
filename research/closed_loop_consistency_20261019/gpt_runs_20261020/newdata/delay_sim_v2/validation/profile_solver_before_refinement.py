"""Scaled bounded trust-region search of the complete profiled IV moment."""
from oracle_model import *
from scipy.optimize import least_squares

RCOND=1e-13

def factor(X):
    norms=np.maximum(la.norm(X,axis=0),1e-300)
    Q,R=la.qr(X/norms,mode='economic',check_finite=False)
    u,s,v=la.svd(R,full_matrices=False,check_finite=False)
    if s[-1]<=RCOND*s[0]:raise la.LinAlgError('tall design rank deficient')
    return dict(Q=Q@u,R=s[:,None]*v*norms[None,:],s=s,norms=norms,v=v)

def coefficient(f,gamma):return (f['v'].T@(gamma/f['s']))/f['norms']

class Limit(Exception):pass
class InnerFailure(Exception):pass

class Profile:
    def __init__(self,m,data,method):
        self.m=m;self.data=data;self.c=m.context(data['f']);self.U=data['U'];self.Y=data['Y']
        self.D=data.get('D');self.K=data.get('K');self.method=method;self.mask=data.get('transient_mask')
        self.calls=0;self.inner_failures=0
    def uh(self,b):
        if self.method=='M1':return self.U
        G=self.m.resp(self.c,b)
        return np.linalg.solve(np.eye(4)+self.K@G,self.D[:,:,None])[:,:,0]
    def profile(self,q,lin,maxinner=50):
        self.calls+=1;m=self.m;b=np.r_[lin.copy(),q]
        off,XL,_,_=m.design(self.c,q,self.U,False,self.mask)
        fx=[factor(rr(XL[:,o])) for o in range(4)]
        targets=[rr(self.Y[:,o]-off[:,o]) for o in range(4)]
        history=[];gamma=[];es=[]
        for it in range(1 if self.method=='M1' else maxinner):
            Uh=self.uh(b)
            ZL=XL if self.method=='M1' else m.design(self.c,q,Uh,False,self.mask)[1]
            gamma=[];es=[]
            for o in range(4):
                Qx=fx[o]['Q']
                if self.method=='M1':g=Qx.T@targets[o]
                else:
                    Qz=factor(rr(ZL[:,o]))['Q'];C=Qz.T@Qx
                    g=la.solve(C,Qz.T@targets[o],assume_a='gen',check_finite=False)
                gamma.append(g);es.append(targets[o]-Qx@g)
                b[o*m.L:(o+1)*m.L]=coefficient(fx[o],g)
            # Re-evaluate tools at the returned coefficients (actual IV root).
            Uh=self.uh(b)
            if self.method!='M1':ZL=m.design(self.c,q,Uh,False,self.mask)[1]
            ee=sum(float(e@e) for e in es)
            gl=np.concatenate([rr(ZL[:,o]).T@es[o] for o in range(4)])
            zz=np.concatenate([np.sum(abs(ZL[:,o])**2,axis=0) for o in range(4)])
            sc=float(np.max(abs(gl)/np.sqrt(np.maximum(zz*ee,1e-300))))
            history.append(sc)
            if sc<1e-9 or ee<1e-36:break
        ok=bool(sc<1e-9 or ee<1e-36)
        if not ok:self.inner_failures+=1
        # Nonlinear moments use exact original columns, evaluated at final l,q.
        JN=m.compact(self.c,b,self.U,True,self.mask)[2]
        ZN=JN if self.method=='M1' else m.compact(self.c,b,Uh,True,self.mask)[2]
        gn=sum((rr(ZN[:,o]).T@es[o] for o in range(4)),start=np.zeros(m.nq))
        zzn=np.sum(abs(ZN)**2,axis=(0,1));g=np.r_[gl,gn];norm=np.r_[zz,zzn]
        fullscore=float(np.max(abs(g)/np.sqrt(np.maximum(norm*ee,1e-300))))
        return dict(beta=b,inner_ok=ok,inner_iterations=len(history),inner_score=sc,
            score=fullscore,g=g,zz=norm,ee=ee,es=es,fx=fx,JN=JN,ZN=ZN,ZL=ZL,XL=XL,
            linear_history=history,minimum_tall_singular_ratio=min(f['s'][-1]/f['s'][0] for f in fx))
    def solve(self,b0,out,maxseconds=None):
        m=self.m;out=Path(out);out.mkdir(parents=True,exist_ok=True)
        start=time.perf_counter();maxseconds=maxseconds or (2700 if m.band=='W' else 900)
        q0=b0[m.nl:].copy();scale=np.r_[np.full(2*m.nm,.04),np.full(7,.4)]
        lower=(m.lower-q0)/scale;upper=(m.upper-q0)/scale
        history=[];jac_history=[];cache={};last=None;den=None;status='not_started';opt=None
        def budget():
            if self.calls>=2000:raise Limit('profile_call_limit')
            if time.perf_counter()-start>maxseconds:raise Limit('wall_time_limit')
        def evaluate(x):
            nonlocal den,last
            key=np.asarray(x).tobytes()
            if key in cache:return cache[key]
            budget();state=self.profile(q0+scale*x,b0[:m.nl] if last is None else last['beta'][:m.nl])
            if den is None:den=np.sqrt(np.maximum(state['zz'][m.nl:]*state['ee'],1e-300))
            state['fun']=state['g'][m.nl:]/den
            record={k:state[k] for k in ('beta','inner_ok','inner_iterations','inner_score','score','ee','minimum_tall_singular_ratio')}
            record.update(evaluation=self.calls,seconds=time.perf_counter()-start,x=x)
            save(out/f'eval_{self.calls:04d}.json',record);history.append(record)
            for heavy in ('es','fx','JN','ZN','ZL','XL'):state.pop(heavy,None)
            cache[key]=state
            if not state['inner_ok']:raise InnerFailure('unresolved_inner_moments')
            return state
        def fun(x):return evaluate(x)['fun']
        def jac(x):
            nonlocal last
            last=evaluate(x)
            # Complete profile derivative, including updated instruments and l(q).
            h=1e-3;J=np.empty((m.nq,m.nq))
            for j in range(m.nq):
                xp=x.copy();xm=x.copy();xp[j]+=h;xm[j]-=h
                J[:,j]=(fun(xp)-fun(xm))/(2*h)
            jac_history.append(dict(evaluation=self.calls,score=last['score']))
            save(out/'progress.json',dict(calls=self.calls,iterations=len(jac_history),last_score=last['score'],seconds=time.perf_counter()-start))
            if len(jac_history)>100:raise Limit('outer_iteration_limit')
            return J
        try:
            opt=least_squares(fun,np.zeros(m.nq),jac=jac,method='trf',bounds=(lower,upper),
                ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=100,x_scale=np.ones(m.nq))
            last=evaluate(opt.x);status=str(opt.message)
        except (Limit,InnerFailure,la.LinAlgError,ValueError) as exc:status=type(exc).__name__+':'+str(exc)
        if last is None and history:
            b=np.asarray(history[0]['beta']);last=self.profile(b[m.nl:],b[:m.nl])
        if last is None:raise RuntimeError(status)
        b=last['beta'];near=bool(np.min(np.minimum(b[m.nl:]-m.lower,m.upper-b[m.nl:]))<1e-5)
        result=dict(beta=b,success=bool(last['score']<1e-6 and last['inner_ok'] and not near and np.isfinite(b).all()),
            score=last['score'],inner_score=last['inner_score'],inner_ok=last['inner_ok'],near_bound=near,
            status=status,profile_calls=self.calls,iterations=len(jac_history),seconds=time.perf_counter()-start,
            nq=m.nq,parameters=m.nb,full_moment_derivative='central profile including instrument updates',
            history=[{k:v for k,v in a.items() if k not in ('beta','x')} for a in history])
        save(out/'result.json',result);return result
