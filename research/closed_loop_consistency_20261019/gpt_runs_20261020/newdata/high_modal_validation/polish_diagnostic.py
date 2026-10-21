"""Post-run numerical diagnostic; never replaces fixed-protocol primary results."""
from run import *
from analyze import loadfit,csvout

def main():
    dest=HERE/'diagnostics';dest.mkdir(exist_ok=True);rows=[]
    for w in ('H100','H300','H493'):
        _,m,lb,ub=setup(w)
        for rid in IDS:
            a=get_data(rid,w);li=a['line_indices']
            for method in ('M1','M2','M3','M3w'):
                r=loadfit(rid,w,method)
                if not r['converged'] or r['near_bound']:continue
                ii=slice(None) if method=='M1' else li
                c=m.context(a['f'][ii]);U,Y,D,K=(a[k][ii] for k in ('U','Y','D','K'))
                instr=(lambda b:np.linalg.solve(np.eye(4)[None]+K@m.resp(c,b),D[:,:,None])[:,:,0]) if method in ('M3','M3w') else U
                W=None if method!='M3w' else wc.precision(np.load(HERE/'raw'/f'{rid}_{w}_local_cov.npz')['covariance_residual'])
                b=np.asarray(r['beta']);f0,z0=m.modal(b);history=[]
                for it in range(25):
                    uh=instr(b) if callable(instr) else instr
                    AA,g,score,ee=equations(m,c,U,Y,uh,b,W)
                    if score<1e-8:break
                    sc=1/np.sqrt(np.maximum(abs(np.diag(AA)),1e-300));step=sc*np.linalg.solve(AA*sc[:,None]*sc[None,:],g*sc)
                    factor=min(1.,.05/max(np.max(abs(step[m.iw:])),1e-300))
                    pp=step>0;nn=step<0
                    if pp.any():factor=min(factor,float(np.min(.995*(ub[pp]-b[pp])/step[pp])))
                    if nn.any():factor=min(factor,float(np.min(.995*(lb[nn]-b[nn])/step[nn])))
                    history.append(dict(score=float(score),factor=float(factor)))
                    if factor<1e-7:break
                    b=b+factor*step
                f,z=m.modal(b);score=equations(m,c,U,Y,instr(b) if callable(instr) else instr,b,W)[2]
                row=dict(id=rid,window=w,method=method,primary_score=r['score'],polished_score=score,
                    primary_strict_pass=r['strict_score_pass'],polished_strict_pass=score<1e-6,
                    extra_iterations=len(history),max_frequency_change_hz=float(np.max(abs(f-f0))),
                    max_damping_change_pp=float(np.max(abs(z-z0))))
                rows.append(row);js(dest/f'{rid}_{w}_{method}.json',dict(**row,beta=b.tolist(),frequency_hz=f.tolist(),damping_percent=z.tolist(),history=history))
    csvout('polish_diagnostic.csv',rows);print('POLISH_DIAGNOSTIC',len(rows),flush=True)

if __name__=='__main__':main()
