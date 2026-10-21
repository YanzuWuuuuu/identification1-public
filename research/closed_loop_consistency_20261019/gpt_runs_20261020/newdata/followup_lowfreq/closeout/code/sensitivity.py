from common import *
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse

def band_spectrum(rid,hi,extended=False):
    a=np.load(SRC/'prepared'/f'{rid}.npz');P=int(a['P']);N=P*P0
    lo,upper=(1,hi+2) if extended else (3,hi)
    k=np.arange(int(np.ceil(lo*N/FS)),int(np.floor(upper*N/FS))+1);f=k*FS/N
    out=dict(k=k,f=f,P=P,N=N,line_indices=np.flatnonzero(k%P==0))
    for q in ('U','Y','D'):out[q]=np.fft.rfft(a[q],axis=0)[k]
    out['K']=evalK(f,meta(rid)['controller'])
    return out

class Model:
    def __init__(self,extra):
        self.extra=extra;self.nb=184 if extra else 158
        self.poles=np.r_[MOD,174,175] if extra else MOD
        self.linear=np.setdiff1d(np.arange(self.nb),self.poles)
        self.bounds={174:(np.log(2*np.pi*35),np.log(2*np.pi*200)),175:(np.log(.005),np.log(2.))} if extra else {}
    def block(self,c,b):
        h,_,ad=m.zoh_2nd(np.exp(b[174]),np.exp(b[175]),1/c.zi)
        return h,ad[0]+ad[1]*c.zi+ad[2]*c.zi**2
    def resp(self,c,b):
        G=m.resp(c,b[:158])
        if self.extra:
            h,_=self.block(c,b);G=G+h[:,None,None]*b[158:174].reshape(4,4)[None]*c.E
        return G
    def out(self,c,b,U,need_jac=True):
        y,j=m.out_jac(c,b[:158],U,need_jac=need_jac)
        if not self.extra:return y,j
        h,ad=self.block(c,b);C=b[158:174].reshape(4,4);tr=b[176:].reshape(4,2)
        tn=tr[:,0][None]+tr[:,1][None]*c.zi[:,None]
        EU=c.E*U[:,None,:]
        y+=h[:,None]*np.einsum('oi,koi->ko',C,EU)+tn/ad[:,None]
        if not need_jac:return y,None
        J=np.zeros((len(c.f),4,self.nb),complex);J[:,:,:158]=j
        for o in range(4):
            J[:,o,158+o*4:162+o*4]=h[:,None]*EU[:,o,:]
            J[:,o,176+2*o]=1/ad;J[:,o,177+2*o]=c.zi/ad
        for p in (174,175):
            bp=b.copy();bp[p]+=1e-6;hp,ap=self.block(c,bp)
            J[:,:,p]=((hp-h)[:,None]*np.einsum('oi,koi->ko',C,EU)+tn/ap[:,None]-tn/ad[:,None])/1e-6
        return y,J
    def initial(self,c,U,Y):
        b=np.zeros(self.nb);b[128:131]=np.log(2*np.pi*np.array([10.6,11.7,16.1]));b[131:134]=np.log(.03)
        if not self.extra:return m.init_lin(c,U,Y,np.ones_like(U.real),b)
        b[174]=np.log(2*np.pi*40);b[175]=np.log(.2)
        yh,J=self.out(c,b,U);X=J[:,:,self.linear].reshape(-1,len(self.linear))
        A=(X.conj().T@X).real;g=(X.conj().T@(Y-yh).ravel()).real
        b[self.linear]+=scaled_solve(A,g)
        return b
    def solve(self,c,U,Y,Uh,b):
        fun=lambda beta,u:self.out(c,beta,u)
        if not self.extra:
            return wc.iv_solve(fun,U,Y,Uh,np.ones_like(U.real),b,m.IW,m.ITR,maxiter=80,tol=1e-5,cap=.05)
        b=b.copy();hist=[];info=dict(success=False,iterations=0)
        L=wc.whitener(None,len(U),4)
        try:
            for it in range(80):
                A,g,_,_=wc.iv_matrices(fun,U,Y,Uh,b,L)
                step=scaled_solve(A,g);active=[]
                # Remove outward-pointing variables that are at an explicit background bound.
                for _ in range(3):
                    new=[q for q,(lo,hi) in self.bounds.items() if
                         (b[q]<=lo+1e-10 and step[q]<0) or (b[q]>=hi-1e-10 and step[q]>0)]
                    if set(new).issubset(active):break
                    active=sorted(set(active+new));free=np.setdiff1d(np.arange(self.nb),active)
                    step=np.zeros(self.nb);step[free]=scaled_solve(A[np.ix_(free,free)],g[free])
                mx=np.max(abs(step[self.poles]))
                if mx>.05:step*=.05/mx
                bn=b+step
                for q,(lo,hi) in self.bounds.items():bn[q]=np.clip(bn[q],lo,hi)
                mx=np.max(abs(bn[self.poles]-b[self.poles]));b=bn
                if not np.all(np.isfinite(b)):raise FloatingPointError('nonfinite endpoint')
                hist.append(dict(iteration=it+1,pole_step=mx,active=active,
                    background_f_hz=float(np.exp(b[174])/2/np.pi),background_damping=float(np.exp(b[175]))))
                if mx<1e-5:break
            info.update(iterations=len(hist),success=bool(len(hist)<80),status='converged' if len(hist)<80 else 'iteration_limit')
            info.update(wc.iv_score(fun,U,Y,Uh,np.ones_like(U.real),b))
        except Exception as e:info.update(success=False,status='exception',exception=repr(e),iterations=len(hist))
        info.update(history=hist,bound_active=[q for q,(lo,hi) in self.bounds.items() if min(abs(b[q]-lo),abs(b[q]-hi))<1e-9])
        return b,info

def noise_phi(a,b,model,query):
    c=m.Ctx(a['f'],TAU,'zoh');Yh,_=model.out(c,b,a['U'],False)
    ni=a['k']%a['P']!=0;v=(a['Y']-Yh)[ni];fn=a['f'][ni];C=[]
    for f in query:
        r=v[abs(fn-f)<=.25+1e-12]
        if len(r)<4:raise ValueError('too few noise points')
        C.append(r.T@r.conj()/len(r)/a['N'])
    C=np.array(C);np.linalg.cholesky(C)
    return C

def residual(model,b,a,phi):
    li=a['line_indices'];f=a['f'][li];c=m.Ctx(f,TAU,'zoh')
    U,Y,D,K=(a[q][li] for q in ('U','Y','D','K'))
    Ud=np.linalg.solve(np.eye(4)[None]+K@model.resp(c,b),D[:,:,None])[:,:,0]
    yh,J=model.out(c,b,U);Z=model.out(c,b,Ud)[1];e=Y-yh
    L=np.linalg.cholesky(.5*realmat(a['N']*phi))
    Jr,Zr=realvec(J),realvec(Z);A=np.einsum('kop,koq->pq',Zr,Jr,optimize=True)
    Bw=np.linalg.solve(L,Jr).reshape(-1,model.nb)
    zL=np.einsum('kop,koi->kpi',Zr,L,optimize=True).transpose(1,0,2).reshape(model.nb,-1)
    H=scaled_solve(A,zL);Rw=np.eye(len(Bw))-Bw@H
    ew=np.linalg.solve(L,realvec(e)[:,:,None])[:,:,0];qf=np.sum(ew**2,axis=1)
    ef=np.sum(Rw**2,axis=1).reshape(len(f),8).sum(1)
    return dict(Q=qf.sum(),EQ=ef.sum(),ratio=qf.sum()/ef.sum(),A_numerics=diagnostics(A),
                projector_error=np.linalg.norm(Rw@Bw)/np.linalg.norm(Bw)),dict(f=f,e=e,q_frequency=qf,expected_frequency=ef)

def run(rid,hi,extra):
    label=f'{rid}_hi{hi}_'+('extra40' if extra else 'base')
    dest=HERE/'analysis'/f'{label}.json'
    if dest.exists():return rd(dest)
    started=time.perf_counter();model=Model(extra)
    a=band_spectrum(rid,hi);ca=m.Ctx(a['f'],TAU,'zoh');li=a['line_indices'];cl=ca.sub(li)
    info=dict(id=rid,hi=hi,background='extra40' if extra else 'D0D2',nb=model.nb,fits={})
    for method in ('M1','M3'):
        t=time.perf_counter()
        try:
            sel=np.arange(len(a['f'])) if method=='M1' else li;c=ca if method=='M1' else cl
            U,Y,D,K=(a[q][sel] for q in ('U','Y','D','K'))
            init=model.initial(c,U,Y)
            Uh=U if method=='M1' else lambda bb:np.linalg.solve(np.eye(4)[None]+K@model.resp(c,bb),D[:,:,None])[:,:,0]
            b,si=model.solve(c,U,Y,Uh,init)
            rr=dict(beta=b,initial_beta=init,solver_info=si,converged=bool(si.get('success',False) and np.isfinite(b).all()),
                    strict_score_pass=bool(si.get('score_max',np.inf)<1e-6),frequency_hz=np.exp(b[128:131])/2/np.pi,
                    damping_percent=100*np.exp(b[131:134]),seconds=time.perf_counter()-t)
            if extra:rr.update(background_f_hz=np.exp(b[174])/2/np.pi,background_damping=np.exp(b[175]),background_bound_active=si.get('bound_active',[]))
        except Exception as e:rr=dict(converged=False,exception=repr(e),traceback=traceback.format_exc())
        info['fits'][method]=rr;js(HERE/'raw'/f'{label}_{method}.json',rr)
    try:
        b1=np.array(info['fits']['M1']['beta']);b=np.array(info['fits']['M3']['beta'])
        info['paired_percent']=100*np.expm1(b1[MOD]-b[MOD])
        ext=band_spectrum(rid,hi,True);phi=noise_phi(ext,b,model,a['f'][li])
        r,arr=residual(model,b,a,phi);info['residual']=r
        np.savez_compressed(HERE/'raw'/f'{label}_residual.npz',Phi=phi,**arr)
        # All configurations are compared at 23--25 Hz with the SAME original covariance.
        ae=band_spectrum(rid,25);eidx=ae['line_indices'];f=ae['f'][eidx];sel=f>=23
        fc=f[sel];c=m.Ctx(fc,TAU,'zoh')
        U,Y=ae['U'][eidx][sel],ae['Y'][eidx][sel]
        yh,_=model.out(c,b,U,False);e=Y-yh
        oldnoise=np.load(BASE/'spectra'/f'{rid}_noise.npz')
        C=ae['N']*interp_cov(fc,oldnoise['f'],oldnoise['Phi_0.5'])
        q=2*np.einsum('ko,koi,ki->k',e.conj(),np.linalg.inv(C),e,optimize=True).real
        info['edge_23_25']=dict(Q=q.sum(),raw_residual_norm=np.linalg.norm(e),points=len(fc),
                               common_noise=True,extrapolation=hi<25)
        np.savez_compressed(HERE/'raw'/f'{label}_edge.npz',f=fc,e=e,Q=q,C=C)
        if hi==25 and not extra:
            old=rd(BASE/'analysis'/f'{rid}_M3_residual.json')
            info['baseline_reproduction']=dict(
                max_abs_beta_M1=np.max(abs(b1-np.array(fit(rid,'M1')['beta']))),
                max_abs_beta_M3=np.max(abs(b-np.array(fit(rid,'M3')['beta']))),
                Q_ratio_relative=abs(r['ratio']-old['ratio'])/old['ratio'])
    except Exception as e:info['analysis_failure']=dict(exception=repr(e),traceback=traceback.format_exc())
    info['seconds']=time.perf_counter()-started;js(dest,info)
    print('SENS',label,[(k,x.get('converged'),x.get('solver_info',{}).get('iterations')) for k,x in info['fits'].items()],
          'Q/E',info.get('residual',{}).get('ratio'),round(info['seconds'],1),flush=True)
    return info
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=2);p.add_argument('--smoke',action='store_true');args=p.parse_args()
    if args.smoke:
        run('A_nominal',25,False);run('A_nominal',25,True)
    else:
        jobs=[(r,h,e) for e in (False,True) for h in (20,25,30) for r in IDS]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            fut=[pool.submit(run,*j) for j in jobs]
            for f in as_completed(fut):f.result()
        js(HERE/'analysis/sensitivity.json',[rd(HERE/'analysis'/f"{r}_hi{h}_{'extra40' if e else 'base'}.json") for r,h,e in jobs])
