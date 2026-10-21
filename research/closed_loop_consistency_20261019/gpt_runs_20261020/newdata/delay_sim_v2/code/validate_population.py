from population_root import *

if __name__=='__main__':
    m=OracleModel('L');base=ROOT/'power/FC15_L_P20_M1'
    with np.load(base/'common.npz') as a:o={k:a[k] for k in a.files}
    with np.load(base/'r000.npz') as a:p={k:a[k] for k in a.files}
    pop=Population(20,'L',15,'M1',o);pop.set_phase(0)
    s=pop.equation(np.zeros(m.nb));scale=p['scales'];gerr=la.norm(scale*(s['g']-o['b']))/la.norm(scale*o['b'])
    # Bounded, scaled direction with a nonzero nuisance component.
    v=np.cos(np.arange(m.nb)+.31)*scale;v*=.01/max(abs(v[m.nl:]))
    h=1e-3
    fd=(pop.equation(h*v)['g']-pop.equation(-h*v)['g'])/(2*h)
    exact=-(p['M']-o['C'])@v
    derr=la.norm(scale*(fd-exact))/la.norm(scale*exact)
    row=dict(bias_moment_identity_relative=float(gerr),full_curvature_derivative_relative=float(derr),
        passed=bool(gerr<1e-8 and derr<1e-4),includes_all_parameters=m.nb)
    save(ROOT/'validation/population_acceptance.json',row);print(row,flush=True)
    if not row['passed']:raise SystemExit(2)
