import os
for n in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ.setdefault(n,'1')
import numpy as np, scipy.io as sio, json, sys, time, hashlib
from pathlib import Path
import realfit_extended as m
ROOT=Path(__file__).resolve().parents[5]
HERE=Path(__file__).resolve().parent
DATA=ROOT/'experiment_data/srivc_no_lrm_20261015/closed_loop_consistency_20261019'
DELAY=ROOT/'experiment_data/srivc_no_lrm_20261015/delay_identification_20261016/effective_delay_matrix_us.csv'
TAU=np.loadtxt(DELAY,delimiter=',')*1e-6

def controller(y,v,na=6,nb=6):
    L=max(na,nb)
    Phi=np.hstack([v[L-i:len(v)-i] for i in range(1,na+1)]+[y[L-j:len(y)-j] for j in range(nb+1)])
    th=np.linalg.lstsq(Phi,v[L:],rcond=None)[0]
    rr=v[L:]-Phi@th
    A=[th[3*i:3*(i+1)].T for i in range(na)]
    B=[th[3*na+4*j:3*na+4*(j+1)].T for j in range(nb+1)]
    def kv(zi):
        az=np.eye(3)[None]-sum(A[i][None]*zi[:,None,None]**(i+1) for i in range(na))
        bz=sum(B[j][None]*zi[:,None,None]**j for j in range(nb+1))
        return np.linalg.solve(az,bz)
    return kv,float(np.max(np.abs(rr))/np.max(np.abs(v)))

def prepare(z,nper,start=1):
    raw=sio.loadmat(DATA/f'cl_data_z{z}.mat')
    p=int(raw['P'][0,0]);fs=float(raw['fs'][0,0]);N=p*nper
    assert fs==m.fs and (start+nper)*p<=len(raw['U'])
    seg=slice(start*p,(start+nper)*p)
    vals={a:raw[a][seg] for a in ['U','Y','D','V']}
    kv,kres=controller(vals['Y'],vals['V'])
    k=np.arange(1,N//2);f=k*fs/N;sel=(f>=3)&(f<=25)
    ff=f[sel];kk=k[sel]
    ffvals={a:np.fft.fft(v,axis=0)[kk] for a,v in vals.items()}
    # Period means use each unnormalised P0-point DFT; covariance conventions
    # are consistent with this scaling and mean transient nuisance coefficients.
    kp=np.arange(1,p//2);fp=kp*fs/p;sp=(fp>=3)&(fp<=25);kp=kp[sp];fp=fp[sp]
    per={a:np.fft.fft(v.reshape(nper,p,-1),axis=1)[:,kp] for a,v in vals.items()}
    repeat=np.linalg.norm(vals['D'].reshape(nper,p,4)-vals['D'][:p][None])/max(np.linalg.norm(vals['D']),1e-300)
    iderr=np.linalg.norm(vals['U']-vals['D']-vals['V']@raw['Q'].T)/np.linalg.norm(vals['U'])
    return dict(z=z,nper=nper,start=start,p=p,fs=fs,f=ff,k=kk,ff=ffvals,kv=kv,Q=raw['Q'],kres=kres,
                fp=fp,per=per,reference_repeat_relative=float(repeat),boundary_identity_relative=float(iderr))

def fit_one(d,sampling,method):
    sl=np.ones(len(d['f']),bool) if method in ['M1','M3'] else d['k']%d['nper']==0
    f=d['f'][sl]; U,Y,D=[d['ff'][x][sl] for x in ['U','Y','D']]
    tau=TAU-(m.Ts/2 if sampling=='zoh' else 0.)
    c=m.Ctx(f,tau,sampling)
    W=np.ones((len(f),4));b=np.zeros(m.NB)
    b[m.IW:m.IZ]=np.log(2*np.pi*np.array([10.6,11.7,16.1]));b[m.IZ:m.ITR]=np.log(.03)
    b=m.init_lin(c,U,Y,W,b)
    if method=='M3':
        K=-d['Q'][None]@d['kv'](c.zi)
        Uh=lambda bb:np.linalg.solve(np.eye(4)[None]+K@m.resp(c,bb),D[:,:,None])[:,:,0]
    else:Uh=U
    tic=time.time();beta,info=m.solve(c,U,Y,Uh,W,b)
    fz,zz=m.modal(beta)
    info.update(z=d['z'],P=d['nper'],start_cycle_zero_based=d['start'],sampling=sampling,method=method,
                frequency_hz=fz.tolist(),damping_percent=zz.tolist(),elapsed_s=time.time()-tic,
                K_arx_relative_max_residual=d['kres'],reference_repeat_relative=d['reference_repeat_relative'],
                boundary_identity_relative=d['boundary_identity_relative'],n_frequency_bins=len(f),n_parameters=m.NB,
                beta=beta.tolist(),initial_beta=b.tolist())
    return info

def save_result(r):
    out=HERE/'raw';out.mkdir(exist_ok=True)
    tag=f"z{r['z']}_P{r['P']}_s{r['start_cycle_zero_based']}_{r['sampling']}_{r['method']}"
    p=out/(tag+'.json');p.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
    np.savez_compressed(out/(tag+'.npz'),beta=r['beta'],initial_beta=r.get('initial_beta',[]))
    print(tag, np.round(r['frequency_hz'],5), np.round(r['damping_percent'],4),r.get('iterations'),r.get('converged_pole_step'),r.get('score_normalized'),flush=True)

def main():
    z=sys.argv[1] if len(sys.argv)>1 else '0'
    ps=[int(x) for x in os.environ.get('C_P','1,2,4,8').split(',')]
    modes=os.environ.get('C_MODE','ct,zoh').split(',')
    methods=os.environ.get('C_METHOD','M1,M2,M3').split(',')
    for p in ps:
        d=prepare(z,p)
        for mode in modes:
            for method in methods:
                tag=f'z{z}_P{p}_s1_{mode}_{method}.json'
                if (HERE/'raw'/tag).exists():continue
                try:r=fit_one(d,mode,method);save_result(r)
                except Exception as e:
                    r=dict(z=z,P=p,sampling=mode,method=method,start_cycle_zero_based=1,status='exception',exception=repr(e))
                    (HERE/'raw').mkdir(exist_ok=True);(HERE/'raw'/tag).write_text(json.dumps(r,indent=2),encoding='utf-8');print(r,flush=True)
    if os.environ.get('C_HALVES','1')=='1':
        for start in [1,5]:
            tag=f'z{z}_P4_s{start}_zoh_M3.json'
            if (HERE/'raw'/tag).exists():continue
            save_result(fit_one(prepare(z,4,start),'zoh','M3'))

if __name__=='__main__': main()
