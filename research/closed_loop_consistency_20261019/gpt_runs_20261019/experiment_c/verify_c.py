import os
os.environ['NOFLEX']='1'
import sys,json,hashlib
sys.dont_write_bytecode=True
import numpy as np
import scipy.signal as ss
import realfit_extended as m
import run_c as r
sys.path.insert(0,str(r.HERE.parents[1]/'code'))
import cl2

def check():
    rng=np.random.default_rng(1002);f=np.geomspace(3,800,300);zk=np.exp(2j*np.pi*f/m.fs)
    maxdiff=0.
    for w,z in zip(2*np.pi*np.array([10.6,11.7,16.1]),[.03,.04,.05]):
        h=m.zoh_2nd(w,z,zk);orig=cl2.zoh_2nd(w,z,zk)
        maxdiff=max(maxdiff,max(float(np.max(np.abs(a-b))) for a,b in zip(h,orig)))
    # Forward DFT identity: noiseless exact sampled dynamics and unknown endpoint
    # states are reconstructed with each modal (tau0+tau1*z^-1)/a_d term.
    N=30000;U=rng.normal(size=(N,4));A,B,C=cl2.ct_ss()
    Ad,Bd,Cd,Dd,_=ss.cont2discrete((A,B,C,np.zeros((4,4))),1/m.fs)
    _,Y,_=ss.dlsim((Ad,Bd,Cd,Dd,1/m.fs),U)
    k=np.arange(18,151);f=k*m.fs/N;Uq=np.fft.fft(U,axis=0)[k];Yq=np.fft.fft(Y,axis=0)[k]
    c=m.Ctx(f,np.zeros((4,4)),'zoh');b=np.zeros(m.NB)
    b[m.IW:m.IZ]=np.log(2*np.pi*cl2.FM);b[m.IZ:m.ITR]=np.log(cl2.ZM)
    for i in range(3):b[32*i:32*i+16]=cl2.RES[i].ravel()
    yh,J=m.out_jac(c,b,Uq);Jt=J[:,:,m.ITR:].reshape(-1,24)
    rhs=(Yq-yh).ravel()
    bt=np.linalg.lstsq(np.r_[Jt.real,Jt.imag],np.r_[rhs.real,rhs.imag],rcond=None)[0]
    b[m.ITR:]=bt;yh,J=m.out_jac(c,b,Uq)
    noiserel=float(np.linalg.norm(Yq-yh)/np.linalg.norm(Yq))
    # Central differences check the implemented forward-difference analytic blocks.
    errs=[]
    for j in range(m.IW,m.ITR):
        h=1e-5;bp=b.copy();bm=b.copy();bp[j]+=h;bm[j]-=h
        fd=(m.out_jac(c,bp,Uq,False)[0]-m.out_jac(c,bm,Uq,False)[0])/(2*h)
        errs.append(float(np.linalg.norm(fd-J[:,:,j])/np.linalg.norm(fd)))
    paths=list((r.HERE.parents[1]/'real_data').glob('*.py'))+list(r.DATA.glob('cl_data_z*.mat'))+[r.DELAY]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    out=dict(zoh_max_absolute_vs_cl2=maxdiff,noiseless_DFT_identity_relative=noiserel,
             modal_jac_relative_error=errs,original_SHA256=hashes)
    assert maxdiff<1e-12
    assert noiserel<1e-10
    assert max(errs)<1e-3
    (r.HERE/'verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='original_SHA256'},indent=2))
check()
