from identify import *
from scipy.stats import chi2

CTRS=('lowgain','nominal','highgain')
def interp_cov(f,fp,C):
    out=np.empty((len(f),4,4),complex)
    for i in range(4):
        for j in range(4):out[:,i,j]=np.interp(f,fp,C[:,i,j].real)+1j*np.interp(f,fp,C[:,i,j].imag)
    return (out+out.conj().transpose(0,2,1))/2
def smooth_cov(V,nl=5,alpha=.1):
    P,nf,ny=V.shape;out=[]
    for k in range(nf):
        a=min(max(k-nl//2,0),nf-nl);r=V[:,a:a+nl].reshape(-1,ny)
        out.append(r.T@r.conj()/len(r)/P0)
    return wc.shrink_cov(np.asarray(out),alpha)
def noise():
    summ=[];rng=np.random.default_rng(621020)
    for ctrl in CTRS:
        a=np.load(HERE/'prepared'/f'B_{ctrl}.npz');p=int(a['P']);kp=np.arange(18,151);f=kp*FS/P0
        U=np.fft.rfft(a['U'].reshape(p,P0,4),axis=1)[:,kp];Y=np.fft.rfft(a['Y'].reshape(p,P0,4),axis=1)[:,kp]
        b=np.array(json.loads((HERE/'fits'/f'A_{ctrl}_M3.json').read_text())['beta']);G=m.resp(m.Ctx(f,TAU,'zoh'),b)
        V=Y-np.einsum('koi,pki->pko',G,U);Phi=smooth_cov(V)
        # Conditional period-block bootstrap of band-integrated trace, fixed G.
        energy=np.sum(abs(V)**2,axis=(1,2))/P0
        draws=np.mean(energy[rng.integers(p,size=(1000,p))],axis=1)
        co=np.load(HERE/'controllers'/f'{ctrl}.npz');K=evalK(f,co['A'],co['B'],co['T'])
        Su=np.linalg.inv(np.eye(4)[None]+K@G);Sy=np.linalg.inv(np.eye(4)[None]+G@K);F=Su@K
        # A second, explicitly diagnostic extraction uses C orthogonal G.
        Gc=np.load(HERE/'analysis'/f'{ctrl}_orthogonal_frf.npz')['G']
        Vc=Y-np.einsum('koi,pki->pko',Gc,U);Pc=smooth_cov(Vc)
        np.savez_compressed(HERE/'analysis'/f'{ctrl}_noise.npz',f=f,Phi=Phi,Phi_C_diagnostic=Pc,Vperiod=V,Uperiod=U,Yperiod=Y,
                            energy_boot=draws,K=K,Su=Su,Sy=Sy,F=F,G_A_reference=G)
        summ.append(dict(controller=ctrl,P_noise=p,band_trace=float(np.mean(energy)),band_trace_ci95=np.percentile(draws,[2.5,97.5]).tolist(),
                         covariance_condition_median=float(np.median(np.linalg.cond(Phi))),alpha=.1,neighbor_bins=5,
                         phi_units='variance of N^-1/2 DFT; divide by fs for two-sided PSD per Hz',
                         C_diagnostic_band_trace=float(np.sum(np.trace(Pc,axis1=1,axis2=2).real))))
    pairs=[]
    for i,a in enumerate(CTRS):
        for b in CTRS[i+1:]:
            aa=np.load(HERE/'analysis'/f'{a}_noise.npz');bb=np.load(HERE/'analysis'/f'{b}_noise.npz')
            rat=aa['energy_boot']/bb['energy_boot'];ci=np.percentile(rat,[2.5,97.5])
            pairs.append(dict(a=a,b=b,trace_ratio=float(aa['energy_boot'].mean()/bb['energy_boot'].mean()),ratio_ci95=ci.tolist(),
                              includes_one=bool(ci[0]<=1<=ci[1]),median_frequency_trace_ratio=float(np.median(np.trace(aa['Phi'],axis1=1,axis2=2).real/np.trace(bb['Phi'],axis1=1,axis2=2).real))))
    js(HERE/'analysis/noise_summary.json',dict(records=summ,comparisons=pairs))

def solve_cov(A,Q):
    sc=1/np.sqrt(np.maximum(np.abs(np.diag(A)),1e-300));As=A*sc[:,None]*sc[None,:]
    sv=np.linalg.svd(As,compute_uv=False);rank=int(np.count_nonzero(sv>sv[0]*1e-12))
    inv=sc[:,None]*np.linalg.solve(As,np.diag(sc));P=inv@Q@inv.T;P=(P+P.T)/2
    return P,dict(scaled_condition=float(sv[0]/sv[-1]),scaled_rank=rank,dimension=len(A),min_variance=float(np.diag(P).min()))

def g_jac(c,b):
    return gmodel_from_output(lambda bb,u:m.out_jac(c,bb,u),m.NB,len(c.f))(b)

def direct_moments(c,b,Ud,K,C,need_quad=True):
    G,Gj=g_jac(c,b);F=np.linalg.solve(np.eye(4)[None]+K@G,K)
    _,J=m.out_jac(c,b,Ud);flat=J.reshape(-1,m.NB);M0=(flat.conj().T@flat).real
    L=np.linalg.cholesky(C);LF=F@L
    noisej=-np.einsum('koip,kil->kolp',Gj,LF,optimize=True).reshape(-1,m.NB)
    Mn=(noisej.conj().T@noisej).real
    bvec=-np.einsum('koip,kil,klo->p',Gj,F,C,optimize=True).real
    qlin=.5*np.einsum('kop,koi,kiq->pq',J.conj(),C,J,optimize=True).real
    Q=qlin
    if need_quad:
        # Re[-v^H G_p F v] = v^H Herm(-G_p F) v.
        H=-np.einsum('koip,kij->kojp',Gj,F,optimize=True)
        H=(H+H.conj().transpose(0,2,1,3))/2
        HL=np.einsum('kai,kabp,kbj->kijp',L.conj(),H,L,optimize=True)
        h=HL.reshape(-1,m.NB);Q=Q+(h.conj().T@h).real
    return M0+Mn,bvec,Q,Gj,F,M0

def covariance_one(path):
    info=json.loads(path.read_text(encoding='utf-8'));rid=info['id'];method=info['method']
    if not info.get('finite_endpoint'):return
    # M5n uses a separate conditional period bootstrap, not theorem 2.
    if method=='M5n':return
    a=prepare_spectra(rid);sel=np.arange(len(a['f'])) if method=='M1' else a['line_indices'];f=a['f'][sel];b=np.asarray(info['beta'])
    c=m.Ctx(f,TAU,'zoh');U=a['U'][sel];Y=a['Y'][sel];D=a['D'][sel];K=a['K'][sel]
    G=m.resp(c,b);Ud=np.linalg.solve(np.eye(4)[None]+K@G,D[:,:,None])[:,:,0]
    no=np.load(HERE/'analysis'/f"{info['controller']}_noise.npz")
    C=int(a['N'])*interp_cov(f,no['f'],no['Phi'])
    note='B-noise spectrum plug-in; fixed reference design; complete transient nuisance block retained'
    if method in ('M1','M2'):
        A,bb,Q,_,_,_=direct_moments(c,b,Ud,K,C)
        note+='; direct method local Gaussian linear+quadratic moment approximation (M curvature omitted here), not theorem 2'
    else:
        _,J=m.out_jac(c,b,Ud);Z=J
        W=np.broadcast_to(np.eye(4),(len(f),4,4))
        if method=='M3w':
            cv=np.load(HERE/'fits'/f'{rid}_local_covariance.npz');W=wc.precision(cv['covariance_residual']);note+='; estimated W frozen, weight-estimation uncertainty omitted'
        wz=np.einsum('koi,kip->kop',W,Z);A=np.einsum('kop,koq->pq',Z.conj(),wz,optimize=True).real
        Q=.5*np.einsum('kop,koi,kiq->pq',wz.conj(),C,wz,optimize=True).real
    P,diag=solve_cov(A,Q);v=np.diag(P)[m.IW:m.ITR];s=np.sqrt(np.maximum(v,0));ff,zz=m.modal(b)
    good=info['converged'] and info.get('strict_score_pass',False) and diag['scaled_rank']==m.NB and np.all(v>=0)
    result=dict(id=rid,method=method,frequency_u_hz=(ff*s[:3]).tolist(),damping_u_percentage_points=(zz*s[3:]).tolist(),
        log_parameter_sd=s.tolist(),relative_u_percent=(100*np.sqrt(np.expm1(v)*np.exp(v))).tolist(),
        locally_valid_numerics=bool(good),solver_converged=info['converged'],strict_score_pass=info.get('strict_score_pass'),
        covariance_conditions=diag,interpretation=note,theta_order='f1 f2 f3 zeta1 zeta2 zeta3',noise_borrowed_from_z0=info['z_mm']!=0)
    # The E spectrum is borrowed, so descriptive uncertainty cannot certify
    # location-specific coverage without independent noise at that location.
    np.savez_compressed(HERE/'analysis'/f'{rid}_{method}_covariance.npz',covariance=P,A=A,Q=Q)
    js(HERE/'analysis'/f'{rid}_{method}_uncertainty.json',result)

def bias(rid):
    info=json.loads((HERE/'fits'/f'{rid}_M3.json').read_text());b=np.asarray(info['beta']);a=prepare_spectra(rid)
    c=m.Ctx(a['f'],TAU,'zoh');K=a['K'];G=m.resp(c,b);Ud=np.linalg.solve(np.eye(4)[None]+K@G,a['D'][:,:,None])[:,:,0]
    no=np.load(HERE/'analysis'/f"{info['controller']}_noise.npz");C=int(a['N'])*interp_cov(a['f'],no['f'],no['Phi'])
    M,bv,Q,Gj,F,M0=direct_moments(c,b,Ud,K,C)
    # Only log modal coordinates affect the transfer-function basis; linear-
    # linear and transient second derivatives vanish. Central outer step>=1e-4.
    CC=np.zeros_like(M);h=2e-4
    for j in range(m.IW,m.ITR):
        bp=b.copy();bm=b.copy();bp[j]+=h;bm[j]-=h
        gp=g_jac(c,bp)[1];gm=g_jac(c,bm)[1]
        derivative=-np.einsum('koip,kil,klo->p',(gp-gm)/(2*h),F,C,optimize=True).real
        CC[:,j]=derivative;CC[j,:]=derivative
    CC=(CC+CC.T)/2
    def sol(A):
        sc=1/np.sqrt(np.maximum(np.abs(np.diag(A)),1e-300));return sc*np.linalg.solve(A*sc[:,None]*sc[None,:],bv*sc)
    d0=sol(M0);d1=sol(M);d2=sol(M-CC)
    direct=json.loads((HERE/'fits'/f'{rid}_M1.json').read_text());db=np.asarray(direct['beta'])-b
    loc=np.r_[m.IW:m.ITR]
    out=dict(id=rid,controller=info['controller'],pred_M_percent=(100*np.expm1(d1[loc])).tolist(),pred_MminusC_percent=(100*np.expm1(d2[loc])).tolist(),
             pred_M0_percent=(100*np.expm1(d0[loc])).tolist(),observed_M1_minus_M3_percent=(100*np.expm1(db[loc])).tolist(),
             ratio_observed_to_MminusC=(np.expm1(db[loc])/np.expm1(d2[loc])).tolist(),reference_converged=info['converged'],
             scope='local first-order finite-design plug-in using A reference IV and same-controller B spectrum; not hardware truth',outer_difference_step=h,
             transient_block_retained=True)
    np.savez_compressed(HERE/'analysis'/f'{rid}_bias.npz',M=M,Ccurvature=CC,b=bv,Q=Q,delta_M=d1,delta_MminusC=d2,delta_M0=d0)
    js(HERE/'analysis'/f'{rid}_bias.json',out)
    # Feedback-noise fractions use all-band physical input energy at the same
    # full-record FFT scaling; report ratios even if >100% rather than clipping.
    cu=F@C@F.conj().transpose(0,2,1);fb=np.trace(cu,axis1=1,axis2=2).real
    rows=[]
    for name,mask in [('all',np.ones(len(c.f),bool)),('lines',a['k']%int(a['P'])==0)]:
        observed=np.sum(abs(a['U'][mask])**2);pred=fb[mask].sum();sig=np.sum(abs(Ud[mask])**2)
        rows.append(dict(selection=name,feedback_noise_to_observed_input_percent=float(100*pred/observed),feedback_noise_fraction_predicted_total_percent=float(100*pred/(pred+sig))))
    js(HERE/'analysis'/f'{rid}_feedback_noise_fraction.json',rows)
    print('BIAS',rid,out['observed_M1_minus_M3_percent'],out['pred_MminusC_percent'],flush=True)

def main():
    noise()
    for p in sorted((HERE/'fits').glob('*.json')):
        try:covariance_one(p)
        except Exception as e:js(HERE/'analysis'/f'{p.stem}_uncertainty.json',dict(failure=True,exception=repr(e),source=str(p)))
    for ctrl in CTRS:
        bias(f'A_{ctrl}');bias(f'C_{ctrl}_pppp')
if __name__=='__main__':main()
