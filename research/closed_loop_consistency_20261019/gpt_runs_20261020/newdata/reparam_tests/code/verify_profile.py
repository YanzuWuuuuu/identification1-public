from profile_core import *
from root_solver import exact_root

def verify(w,joint):
    rid='A_nominal';method='M3'
    m=ResidueModel(w,'a',joint=joint)
    fb,orig=load_joint(rid,w,method) if joint else old.basefit(rid,w,method)
    a=old.spectrum(rid,w);pr=ProfileIV(m,a,method)
    b0=m.from_full(fb);rebuilt,T=m.expand(b0)
    r=dict(window=w,joint=joint,source_score=orig.get('solver',{}).get('score'),
           coordinate_roundtrip=float(np.max(abs(fb-rebuilt))),parameters=m.nb)
    # Original full moment equation refined to a genuinely tight root.
    full=m.full;c,U,Y,Uh=old.inputs(full,a,method)
    if not joint:full.base_nb=full.nb
    froot,fi=exact_root(full,c,U,Y,Uh,fb,maxfev=2500)
    bseed=m.from_full(froot)
    # Test from a nearby perturbed nonlinear state, not just re-evaluate the root.
    bseed[m.nl:m.nl+2*m.nm]+=np.linspace(-1,1,2*m.nm)*1e-6
    bout,info=pr.solve(bseed,maxouter=25,tol=2e-12)
    mapped,_=m.expand(bout)
    # Absolute root difference in stored coordinates; separately scale linear terms.
    dd=mapped-froot
    scale=np.maximum(abs(froot),1e-8)
    r.update(original_full_refinement=fi,profile=info,
             max_absolute_difference=float(np.max(abs(dd))),
             max_scaled_difference=float(np.max(abs(dd)/scale)),
             nonlinear_max_difference=float(np.max(abs(dd[full.iw:]))),
             passed=bool(np.max(abs(dd))<1e-8 and fi['score']<1e-8 and info['score']<1e-8),
             original_beta=froot,profile_beta=mapped)
    js(ROOT/'analysis'/f'verify_{w}_{joint}.json',r)
    print(w,joint,{k:r[k] for k in ('passed','max_absolute_difference','max_scaled_difference','nonlinear_max_difference')},fi['score'],info['score'],flush=True)

def derivative_test():
    m=ResidueModel('H493','b');fb,_=load_joint('A_nominal','H493','M3');b=m.from_full(fb)
    f,T=m.expand(b);rng=np.random.default_rng(731);d=rng.normal(size=m.nb)
    # Linear directions use their own physical scale.
    d[:m.nl]*=np.maximum(abs(b[:m.nl]),1e-10)
    eps=1e-5
    numerical=(m.expand(b+eps*d,False)[0]-m.expand(b-eps*d,False)[0])/(2*eps)
    analytical=T@d
    err=np.linalg.norm(numerical-analytical)/np.linalg.norm(analytical)
    js(ROOT/'analysis/coordinate_derivative.json',dict(relative_error=err,passed=bool(err<1e-7)))
    print('rank1 coordinate derivative',err,flush=True)
if __name__=='__main__':
    derivative_test()
    verify('L',False)
    verify('H493',True)
