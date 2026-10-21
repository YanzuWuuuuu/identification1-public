from bias import *

def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))
def main():
    rid='A_nominal';w='H300';_,m,_,_=setup(w);full=spectrum(rid,w)
    with np.load(HIGH/'spectra'/f'{rid}_{w}.npz') as old:
        reconstruction={k:float(np.max(abs(full[k]-old[k]))) for k in ('f','U','Y','D','K')}
    b=np.asarray(basefit(rid,w,'M3')['beta'])
    ix=np.linspace(0,len(full['f'])-1,650).astype(int)
    a={k:(v[ix] if isinstance(v,np.ndarray) and len(v)==len(full['f']) else v) for k,v in full.items()}
    # Use positive, cross-channel covariance for the independent dense check.
    rng=np.random.default_rng(48215);L=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    cov=np.broadcast_to(1e-10*(L@L.conj().T+np.eye(4)),(len(ix),4,4)).copy()
    pr=moments(m,a,b,cov,chunk=113)
    c=m.context(a['f']);G,Gj=m.gmodel(c,b);F=np.linalg.solve(np.eye(4)[None]+a['K']@G,a['K'])
    Ud=np.linalg.solve(np.eye(4)[None]+a['K']@G,a['D'][:,:,None])[:,:,0];_,J=m.output(c,b,Ud)
    flat=J.reshape(-1,m.nb);M0=(flat.conj().T@flat).real
    NJ=-np.einsum('koip,kil->kolp',Gj,F@np.linalg.cholesky(cov)).reshape(-1,m.nb)
    M=M0+(NJ.conj().T@NJ).real
    bv=-np.einsum('koip,kil,klo->p',Gj,F,cov).real
    CC=np.zeros((m.nb,m.nb))
    for j in range(m.iw,m.nb):
        bp=b.copy();bm=b.copy();bp[j]+=2e-4;bm[j]-=2e-4
        der=(m.gmodel(c,bp)[1]-m.gmodel(c,bm)[1])/4e-4
        v=-np.einsum('koip,kil,klo->p',der,F,cov).real
        CC[:,j]=v;CC[j,:]=v
    CC=(CC+CC.T)/2
    errors=dict(M0=rel(pr['M0'],M0),M=rel(pr['M'],M),b=rel(pr['b'],bv),C=rel(pr['Ccurvature'],CC),
        delta_M=rel(pr['delta_M'],scaled_solve(M,bv)),
        delta_MminusC=rel(pr['delta_MminusC'],scaled_solve(M-CC,bv)))
    ext=spectrum(rid,w,extended=True);Phi,nc=local_phi(m,ext,b,full['f'])
    vv=(ext['Y']-m.output(m.context(ext['f']),b,ext['U'],False)[0])[ext['k']%ext['P']!=0]
    fn=ext['f'][ext['k']%ext['P']!=0];direct=[]
    for fq in full['f'][ix]:
        q=vv[np.abs(fn-fq)<=.25+1e-12];direct.append(q.T@q.conj()/len(q)/ext['N'])
    spectrum_error=rel(Phi[ix],np.array(direct))
    # The step is shared with low frequency. Report sensitivity, do not tune.
    pr_half=moments(m,a,b,cov,chunk=113,step=1e-4)
    step_difference=rel(pr_half['Ccurvature'],pr['Ccurvature'])
    selected=[]
    for ww in SPECS:
        for rr in IDS:
            for method in ('M1','M3'):
                z=basefit(rr,ww,method);selected.append(dict(id=rr,window=ww,method=method,
                    eligible=eligible(z),original_score=z['original_score'],score=z['score']))
    out=dict(reconstruction=reconstruction,chunk_vs_dense=errors,local_spectrum_prefix_vs_direct=spectrum_error,
        curvature_step_2e4_vs_1e4_relative=step_difference,selected_endpoints=selected,
        passed=bool(max(reconstruction.values())<1e-10 and max(errors.values())<1e-5 and spectrum_error<1e-10 and all(x['eligible'] for x in selected)))
    js(HERE/'analysis/validation.json',out);print(json.dumps(out,indent=2))
    if not out['passed']:raise RuntimeError('Validation failed')
if __name__=='__main__':main()
