"""Separate exact known stochastic forcing from finite-record boundary terms.

No object fitting: only projection of an exactly known boundary term onto
existing transient columns and optional diagnostic polynomial columns.
"""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from sim_core import *
from qr_profile import real_rows

def plant_open():
    A,B,C=cl.ct_ss()
    Ad,Bd,Cd,_,_=signal.cont2discrete((A,B,C,np.zeros((4,4))),1/FS)
    nm=len(A);n=nm+sum(NI)+sum(NO);ii=[];oo=[];off=nm
    for nn in NI:ii.append(np.arange(off,off+nn));off+=nn
    for nn in NO:oo.append(np.arange(off,off+nn));off+=nn
    def up(x,q):
        raw=Cd@x[:nm]
        y=np.array([x[ix[-1]] if len(ix) else raw[j] for j,ix in enumerate(oo)])
        qd=np.array([x[ix[-1]] if len(ix) else q[j] for j,ix in enumerate(ii)])
        xn=np.zeros(n);xn[:nm]=Ad@x[:nm]+Bd@qd
        for j,ix in enumerate(ii):
            if len(ix):xn[ix]=np.r_[q[j],x[ix[:-1]]]
        for j,ix in enumerate(oo):
            if len(ix):xn[ix]=np.r_[raw[j],x[ix[:-1]]]
        return xn,y
    Ao=np.zeros((n,n));Bo=np.zeros((n,4));Co=np.zeros((4,n));Do=np.zeros((4,4))
    for i in range(n):Ao[:,i],Co[:,i]=up(np.eye(n)[i],np.zeros(4))
    for i in range(4):Bo[:,i],Do[:,i]=up(np.zeros(n),np.eye(4)[i])
    return Ao,Bo,Co,Do

P=10;seed=1002
rng=np.random.default_rng(seed)
kk=np.arange(1,P0//2);dp=np.zeros((P0//2+1,4),complex)
dp[kk]=.2*np.sqrt(P0*cl.dshape(np.exp(2j*np.pi*kk/P0)))[:,None]*np.exp(2j*np.pi*rng.random((len(kk),4)))
d=np.tile(np.fft.irfft(dp,n=P0,axis=0),(P+2,1));nn=len(d)
v=rng.normal(size=(nn,4))*3e-7;w=rng.normal(size=(nn,4))*.03
AA,BB,CC,DD=plant_with_delays();forcing=np.c_[d,v,w]
_,out,states=signal.dlsim((AA,BB,CC,DD,1/FS),forcing)
saved=dict(np.load(ROOT/"mc/data/P10_r000.npz"));k=saved["k"]
Uf=np.fft.rfft(out[2*P0:,4:],axis=0)[k]
Yf=np.fft.rfft(out[2*P0:,:4],axis=0)[k]
Vf=np.fft.rfft(v[2*P0:],axis=0)[k]
Wf=np.fft.rfft(w[2*P0:],axis=0)[k]
m=SimModel(550.);b=m.truth();c=m.full.context(saved["f"]);G=m.resp(c,b)
noise=Vf+np.einsum("koi,ki->ko",G,Wf)
boundary=Yf-np.einsum("koi,ki->ko",G,Uf)-noise
Ao,Bo,Co,Do=plant_open();ns=len(Ao)
dx=states[2*P0,:ns]-(AA@states[-1]+BB@forcing[-1])[:ns]
exact=np.empty_like(boundary)
for start in range(0,len(k),400):
    sl=slice(start,start+400);z=c["z"][sl]
    exact[sl]=np.linalg.solve(z[:,None,None]*np.eye(ns)-Ao,z[:,None,None]*dx[None,:,None])[:,:,0]@Co.T
summary=dict(seed=seed,periods=P,
    regenerated_U_relative=float(la.norm(Uf-saved["U"])/la.norm(Uf)),
    regenerated_Y_relative=float(la.norm(Yf-saved["Y"])/la.norm(Yf)),
    boundary_state_formula_relative=float(la.norm(exact-boundary)/la.norm(boundary)),
    formula="boundary = C (z I - A_open)^(-1) z (x_start - x_end); noise = V + G W",
    rows=[])
_,LL,_=m.full.compact(c,b,Uf,False)
cols=np.ravel([[10*j+8,10*j+9] for j in range(8)])
T=LL[:,0,cols]
for band in ("L","H","W"):
    for select_type in ("all","excited_lines"):
        ix=(saved["f"]>=BANDS[band][0])&(saved["f"]<=BANDS[band][1])
        if select_type=="excited_lines":ix &= k%P==0
        X=real_rows(T[ix]);B=real_rows(boundary[ix])
        scale=1/la.norm(X,axis=0);theta,_,rank,_=la.lstsq(X*scale,B,cond=1e-13)
        residual=B-(X*scale)@theta
        pp=np.column_stack([np.ones(sum(ix)),c["zi"][ix],c["zi"][ix]**2])
        XX=real_rows(np.column_stack([T[ix],pp]));sc=1/la.norm(XX,axis=0)
        theta2,_,rank2,_=la.lstsq(XX*sc,B,cond=1e-13)
        residual2=B-(XX*sc)@theta2
        eta=noise[ix]
        row=dict(band=band,selection=select_type,points=int(sum(ix)),
            noise_norm=float(la.norm(eta)),boundary_norm=float(la.norm(B)),
            boundary_over_noise=float(la.norm(B)/la.norm(eta)),
            modal_transient_rank=int(rank),
            after_modal_transients_norm=float(la.norm(residual)),
            after_modal_transients_over_noise=float(la.norm(residual)/la.norm(eta)),
            after_modal_and_register_polynomials_norm=float(la.norm(residual2)),
            after_modal_and_register_polynomials_over_noise=float(la.norm(residual2)/la.norm(eta)),
            augmented_rank=int(rank2))
        summary["rows"].append(row)
        print(row,flush=True)
js(ROOT/"audit/boundary_noise_check.json",summary)
np.savez_compressed(ROOT/"audit/boundary_noise_components.npz",f=saved["f"],k=k,
                    known_noise=noise,boundary=boundary,exact_state_boundary=exact)
