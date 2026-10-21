from core import *
def main():
    r=dict(structure=delay_values(TAU),cases=[])
    rng=np.random.default_rng(11020)
    for w in SPECS:
        _,m,lb,ub=setup(w,joint=True);b0,_=basefit('A_nominal',w,'M3');b=np.r_[b0,rng.normal(size=7)*.1]
        f=np.linspace(*SPECS[w]['band'],211);c=m.context(f);U=rng.normal(size=(211,4))+1j*rng.normal(size=(211,4))
        y,J=m.output(c,b,U);G,Gj=m.gmodel(c,b);err=[];errg=[]
        for j in range(m.base_nb,m.nb):
            step=1e-4;bp=b.copy();bm=b.copy();bp[j]+=step;bm[j]-=step
            fd=(m.output(c,bp,U,False)[0]-m.output(c,bm,U,False)[0])/(2*step)
            fdg=(m.resp(c,bp)-m.resp(c,bm))/(2*step)
            err.append(np.linalg.norm(fd-J[:,:,j])/np.linalg.norm(fd))
            errg.append(np.linalg.norm(fdg-Gj[:,:,:,j])/np.linalg.norm(fdg))
        _,base,_,_=setup(w);bzero=np.r_[b0,np.zeros(7)]
        identical=np.linalg.norm(m.output(c,bzero,U,False)[0]-base.output(c,b0,U,False)[0])/np.linalg.norm(y)
        q=dict(window=w,delay_derivative_max=max(err),transfer_derivative_max=max(errg),zero_delay_vs_base=identical)
        assert max(err)<1e-5 and max(errg)<1e-5 and identical<1e-12
        r['cases'].append(q)
    js(ROOT/'analysis/validation.json',r);print(json.dumps(r,indent=2,default=lambda x:x.tolist()))
if __name__=='__main__':main()
