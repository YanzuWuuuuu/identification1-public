"""Independent, short read-only audit of saved QR endpoints; no refitting."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"code"))
from qr_profile import *

def label(i, m):
    if i < m.nl:
        o, r = divmod(i, m.full.L)
        if r < 80:
            j, h = divmod(r, 10)
            return f"y{o+1}:mode{j+1}:" + (f"C0_u{h+1}" if h<4 else f"C1_u{h-3}" if h<8 else f"transient{h-8}")
        return f"y{o+1}:" + (f"D0_u{r-79}" if r<84 else f"D2_u{r-83}")
    if i<m.iz: return f"log_frequency_{i-m.nl+1}"
    if i<m.delay_slice.start:return f"log_damping_{i-m.iz+1}"
    return f"delay_coordinate_{i-m.delay_slice.start}"

def check(band, which):
    a=dict(np.load(ROOT/"mc/data/P10_r000.npz"))
    m=SimModel(BANDS[band][1])
    rr=rd(ROOT/"mc_qr_final/raw"/f"P10_r000_{band}_M3.json")
    b=np.array(rr["endpoint"]["beta"]) if which=="endpoint" else m.truth() if which=="truth" else np.array(rr["initial_beta"])
    pr=QRProfileIV(m,select(a,band,"M3"),"M3")
    uh=pr.uh(b)
    y,J=m.output(pr.c,b,pr.U)
    _,Z=m.output(pr.c,b,uh)
    e=pr.Y-y
    jf=J.reshape(-1,m.nb);zf=Z.reshape(-1,m.nb);ef=e.ravel()
    Ad=(zf.conj().T@jf).real
    gd=(zf.conj().T@ef).real
    zz=np.sum(abs(zf)**2,axis=0);ee=float(np.sum(abs(ef)**2))
    scores=abs(gd)/np.sqrt(zz*ee)
    A,g,sc,ec,zc=pr.eq(b)
    sort=np.argsort(scores)[-12:][::-1]
    _,XL,N=m.full.compact(pr.c,b,pr.U)
    comp=XL*b[:m.nl].reshape(1,4,m.full.L)
    cancellation=float(la.norm(np.sum(abs(comp),axis=2))/la.norm(y))
    # Analytic derivative checks in all nonlinear coordinates; no fit.
    fd={}
    for i in range(m.nl,m.nb):
        h=1e-5
        bp=b.copy();bm=b.copy();bp[i]+=h;bm[i]-=h
        fy=(m.output(pr.c,bp,pr.U,False)[0]-m.output(pr.c,bm,pr.U,False)[0])/(2*h)
        fd[label(i,m)]=float(la.norm(fy-J[:,:,i])/max(la.norm(J[:,:,i]),1e-300))
    result=dict(band=band,at=which,
        dense_chunk_A_relative=float(la.norm(A-Ad)/la.norm(Ad)),
        dense_chunk_g_relative=float(la.norm(g-gd)/la.norm(gd)),
        dense_chunk_zz_relative=float(la.norm(zc-zz)/la.norm(zz)),
        dense_score=float(max(scores)),chunk_score=sc,linear_score=float(max(scores[:m.nl])),
        highest_score=[dict(parameter=label(int(i),m),score=float(scores[i])) for i in sort],
        output_relative_residual=float(la.norm(e)/la.norm(pr.Y)),
        cancellation_factor=cancellation,
        max_coefficient=float(max(abs(b[:m.nl]))),
        nonlinear_fd_relative=fd,linear_profile=rr["linear_profile_diagnostics"]["inner_success"] if which=="endpoint" else None,
        feedback_solve_identity_relative=float(la.norm(uh+np.einsum('koi,ki->ko',pr.K@m.resp(pr.c,b),uh)-pr.D)/la.norm(pr.D)))
    print(band,which, result["dense_score"], result["linear_score"],result["highest_score"][:2], "cancel",cancellation,flush=True)
    return result

if __name__=="__main__":
    rows=[]
    for band in ("L","H"):
        for point in ("truth","initial","endpoint"):
            rows.append(check(band,point))
            js(ROOT/"audit/implementation_check.json",rows)
