"""Causal FIR block decoupling. Design/selection see ONLY training arrays.

No hardware operations. Independent half-record checks are loaded for final
evaluation after the candidates have been selected and serialized.
"""
from pathlib import Path
import argparse, json, hashlib, time, platform
import numpy as np
from scipy.io import loadmat, savemat
from scipy.linalg import eigh
import cvxpy as cp

OUT = Path(__file__).resolve().parent
CFG = dict(fs=5000., taps=[40,100], delays=[0,2,5], main=[25,80],
           guards=[[3,25],[80,170]], targets=[[260,310],[320,380],[450,800]],
           gain_max=1.05, gain_design_margin=.0005, guard_max=1.05, target_max=.10, ridge=1e-5,
           max_sweeps=4, feasibility_tolerance=2e-4)

def read_coords():
    c=loadmat(OUT/'coordinates.mat', squeeze_me=True)
    return {k:np.asarray(c[k]) for k in ['Py','Pu','ny','nu','nyNom','nuNom']}

def physical(H):
    return H.reshape(len(H),4,4).transpose(0,2,1)

def project(H,c,nom=False):
    ny=c['nyNom'] if nom else c['ny']; nu=c['nuNom'] if nom else c['nu']
    Ty=np.vstack([c['Py'],ny]); Tu=np.column_stack([c['Pu'],nu])
    return np.einsum('ab,kbc,cd->kad',Ty,physical(H),Tu)

def idx(f,band): return (f>=band[0]-1e-8)&(f<=band[1]+1e-8)
def quad(f):
    df=np.diff(f)
    return np.r_[df[0],df[:-1]+df[1:],df[-1]]/2
def energy(V,f):
    return float(np.sum(quad(f)*np.sum(np.abs(V.reshape(len(f),-1))**2,axis=1)))
def rms(V,f): return np.sqrt(energy(V,f)/(f[-1]-f[0]))
def basis(f,taps): return np.exp(-2j*np.pi*f[:,None]*np.arange(taps)/CFG['fs'])
def filt(f,h): return basis(f,len(h))@h

def apply(P,f,hu,hy,delay):
    Fu=filt(f,hu); Fy=filt(f,hy); q=np.exp(-2j*np.pi*f*delay/CFG['fs'])
    A=P[:,:3,:3]; b=P[:,:3,3]; c=P[:,3,:3]; d=P[:,3,3]
    Eu=q[:,None]*b-np.einsum('kij,kj->ki',A,Fu)
    Ey=q[:,None]*c-np.einsum('ki,kij->kj',Fy,A)
    dn=q*q*d-q*np.sum(c*Fu,axis=1)-np.sum(Fy*Eu,axis=1)
    new=P.copy();new[:,:3,3]=Eu;new[:,3,:3]=Ey;new[:,3,3]=dn
    return new,Fu,Fy

def cross_design(P,phi,side):
    A=P[:,:3,:3]
    if side=='y': A=A.transpose(0,2,1)
    return np.einsum('kij,kl->kilj',A,phi).reshape(len(P)*3,-1)

def factor(M,t,weights):
    """Exact real Gram compression of ||sqrt(w)(Mh-t)|| for real h."""
    Z=np.column_stack([M,-t])*np.sqrt(weights)[:,None]
    gram=np.real(Z.conj().T@Z);gram=(gram+gram.T)/2
    ev,V=eigh(gram);B=np.sqrt(np.maximum(ev,0))[:,None]*V.T
    err=np.linalg.norm(B.T@B-gram)/max(np.linalg.norm(gram),1e-30)
    assert err<1e-10
    return B

def cross_factor(P,f,phi,q,side,band):
    m=idx(f,band); ff=f[m]; block=P[m,3,:3] if side=='y' else P[m,:3,3]
    denom=energy(block,ff)
    return factor(cross_design(P[m],phi[m],side),(q[m,None]*block).ravel(),
                  np.repeat(quad(ff)/denom,3))

def metrics(P,f,hu,hy,delay,coords,gain_dense=4097):
    new,Fu,Fy=apply(P,f,hu,hy,delay)
    res={}
    for band in [CFG['main']]+CFG['guards']:
        m=idx(f,band); key=f'{band[0]}_{band[1]}'
        for side in ['u','y']:
            B=new[m,:3,3] if side=='u' else new[m,3,:3]
            R=P[m,:3,3] if side=='u' else P[m,3,:3]
            res[f'{side}_{key}']=np.sqrt(energy(B,f[m])/energy(R,f[m]))
    q=np.exp(-2j*np.pi*f*delay/CFG['fs'])
    for band in CFG['targets']:
        m=idx(f,band); ref=q[m]**2*P[m,3,3]
        res[f'nn_{band[0]}_{band[1]}']=np.sqrt(energy(new[m,3,3]-ref,f[m])/energy(ref,f[m]))
    fg=np.linspace(0,CFG['fs']/2,gain_dense); Q=np.exp(-2j*np.pi*fg*delay/CFG['fs'])
    up=Q[:,None]*coords['nu']-filt(fg,hu)@coords['Pu'].T
    yp=Q[:,None]*coords['ny']-filt(fg,hy)@coords['Py']
    res['gain_u']=float(np.max(np.linalg.norm(up,axis=1)))
    res['gain_y']=float(np.max(np.linalg.norm(yp,axis=1)))
    res['score_max']=max(res['u_25_80'],res['y_25_80'])
    res['score_sum']=res['u_25_80']**2+res['y_25_80']**2
    tol=CFG['feasibility_tolerance']
    res['feasible']=bool(all(res[f'{s}_{a}_{b}']<=CFG['guard_max']+tol for s in ['u','y'] for a,b in CFG['guards'])
       and all(res[f'nn_{a}_{b}']<=CFG['target_max']+tol for a,b in CFG['targets'])
       and max(res['gain_u'],res['gain_y'])<=CFG['gain_max']+tol)
    return res

def solve_side(P,f,phi,q,coords,side,other,x0,delay):
    taps=phi.shape[1];n=taps*3
    x=cp.Variable(n);h=cp.reshape(x,(taps,3),order='C')
    Bmain=cross_factor(P,f,phi,q,side,CFG['main'])
    objective=cp.sum_squares(Bmain[:,:n]@x+Bmain[:,n])+CFG['ridge']*cp.sum_squares(x)
    constraints=[]
    for band in CFG['guards']:
        B=cross_factor(P,f,phi,q,side,band)
        constraints.append(cp.norm(B[:,:n]@x+B[:,n])<=CFG['guard_max'])
    otherF=phi@other
    for band in CFG['targets']:
        m=idx(f,band);A=P[m,:3,:3];bb=P[m,:3,3];cc=P[m,3,:3]
        of=otherF[m]
        if side=='y':
            v=q[m,None]*bb-np.einsum('kij,kj->ki',A,of)
            constant=-q[m]*np.sum(cc*of,axis=1)
        else:
            v=q[m,None]*cc-np.einsum('ki,kij->kj',of,A)
            constant=-q[m]*np.sum(of*bb,axis=1)
        M=np.einsum('kl,kj->klj',phi[m],v).reshape(sum(m),n)
        B=factor(M,constant,quad(f[m])/energy(P[m,3,3],f[m]))
        constraints.append(cp.norm(B[:,:n]@x+B[:,n])<=CFG['target_max'])
    # Gain grid covers DC through Nyquist, not only the measured band.
    fg=np.linspace(0,CFG['fs']/2,max(321,6*taps+1));pg=basis(fg,taps)
    qg=np.exp(-2j*np.pi*fg*delay/CFG['fs'])
    directions=coords['Py'] if side=='y' else coords['Pu'].T
    anchor=coords['ny'] if side=='y' else coords['nu']
    port=qg[:,None]*anchor-(pg@h)@directions
    constraints.append(cp.norm(cp.hstack([cp.real(port),cp.imag(port)]),axis=1)<=CFG['gain_max']-CFG['gain_design_margin'])
    problem=cp.Problem(cp.Minimize(objective),constraints)
    x.value=x0.ravel()
    started=time.perf_counter()
    problem.solve(solver='CLARABEL',max_iter=150,tol_gap_abs=2e-8,tol_gap_rel=2e-8,tol_feas=2e-8)
    status=problem.status
    if x.value is None: raise RuntimeError(f'No solution: {status}')
    info=dict(side=side,status=status,objective=float(problem.value),seconds=time.perf_counter()-started)
    return np.asarray(x.value).reshape(taps,3),info

def run_candidate(source,taps,delay,f,Pdesign,Ptrain,coords):
    name=f'{source}_n{taps}_d{delay}';path=OUT/f'{name}.npz'
    if path.exists():
        with np.load(path) as d:return json.loads(str(d['summary']))
    hu=np.zeros((taps,3));hy=hu.copy();phi=basis(f,taps);q=np.exp(-2j*np.pi*f*delay/CFG['fs'])
    assert metrics(Pdesign,f,hu,hy,delay,coords)['feasible']
    history=[];previous=np.inf
    for sweep in range(CFG['max_sweeps']):
        hy,iy=solve_side(Pdesign,f,phi,q,coords,'y',hu,hy,delay)
        hu,iu=solve_side(Pdesign,f,phi,q,coords,'u',hy,hu,delay)
        sm=metrics(Pdesign,f,hu,hy,delay,coords)
        history.append(dict(sweep=sweep+1,y=iy,u=iu,design_metrics=sm))
        print(name,'sweep',sweep+1,'design u/y',round(sm['u_25_80'],5),round(sm['y_25_80'],5),
              'gain',round(sm['gain_u'],5),round(sm['gain_y'],5),flush=True)
        total=iy['objective']+iu['objective']
        if abs(previous-total)<1e-5:break
        previous=total
    tr=metrics(Ptrain,f,hu,hy,delay,coords)
    # Diagnostic same-basis unconstrained fit, never treated as deployable.
    free={}
    for side in ['u','y']:
        B=cross_factor(Pdesign,f,phi,q,side,CFG['main']);n=taps*3
        sol=np.linalg.lstsq(np.vstack([B[:,:n],np.sqrt(CFG['ridge'])*np.eye(n)]),
                            np.r_[-B[:,n],np.zeros(n)],rcond=None)[0]
        free[side]=float(np.linalg.norm(B[:,:n]@sol+B[:,n]))
    summary=dict(name=name,source=source,taps=taps,delay_samples=delay,train=tr,
                 design=sm,history=history,unconstrained_regularized_design_rms=free)
    np.savez_compressed(path,hu=hu,hy=hy,summary=json.dumps(summary))
    print('SAVED',name,'train u/y',tr['u_25_80'],tr['y_25_80'],'feasible',tr['feasible'],flush=True)
    return summary

def write_csv(path,rows):
    import csv
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with open(path,'w',newline='',encoding='utf-8-sig') as h:
        writer=csv.DictWriter(h,fieldnames=keys);writer.writeheader();writer.writerows(rows)

def evaluate(selected,f,coords):
    # All choices serialized before check data are read below.
    results={};frozen={}
    for source,choice in selected.items():
        if choice['name']=='static_fallback':hu=np.zeros((1,3));hy=hu.copy();delay=0
        else:
            data=np.load(OUT/(choice['name']+'.npz'));hu=data['hu'];hy=data['hy'];delay=choice['delay_samples']
        frozen[source]=(hu,hy,delay)
    rows=[];npz={'f':f,'fs':CFG['fs']}
    for tag in ['z0','zm0p6','zm0p3','z0p5','z1']:
        raw=loadmat(OUT/f'input_{tag}.mat',squeeze_me=True);ff=raw['f']
        P=project(raw['Hcheck'],coords);Pnom=project(raw['Hcheck'],coords,nom=True)
        if tag=='z0':npz.update(Pstatic=P,Pnom=Pnom)
        for source,(hu,hy,delay) in frozen.items():
            met=metrics(P,ff,hu,hy,delay,coords)
            rows.append(dict(position=tag,source=source,**met))
            new,Fu,Fy=apply(P,ff,hu,hy,delay)
            # Independent full matrix multiplication sanity check.
            for k in [0,len(ff)//2,len(ff)-1]:
                q=np.exp(-2j*np.pi*ff[k]*delay/CFG['fs'])
                U=np.eye(4,dtype=complex);Y=U.copy();U[:3,3]=-Fu[k];U[3,3]=q
                Y[3,:3]=-Fy[k];Y[3,3]=q
                assert np.linalg.norm(Y@P[k]@U-new[k])<1e-12*max(np.linalg.norm(P[k]),1e-30)
            if tag=='z0':
                npz.update({f'P_{source}':new,f'Fu_{source}':Fu,f'Fy_{source}':Fy,f'delay_{source}':delay})
    np.savez_compressed(OUT/'selected_results.npz',**npz)
    write_csv(OUT/'heldout_metrics.csv',rows)
    payload={'Ts':1/CFG['fs'],'Ty':np.vstack([coords['Py'],coords['ny']]),
             'Tu':np.column_stack([coords['Pu'],coords['nu']])}
    for source,(hu,hy,delay) in frozen.items():
        payload.update({f'hu_{source}':hu,f'hy_{source}':hy,f'delay_{source}':delay})
    savemat(OUT/'dynamic_fir_coefficients.mat',payload)
    print('HELDOUT z0',json.dumps(rows[:2],indent=2),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    coords=read_coords();raw=loadmat(OUT/'input_z0.mat',squeeze_me=True,variable_names=['f','Htrain','Hmodel'])
    f=raw['f'];Ptrain=project(raw['Htrain'],coords);Pmodel=project(raw['Hmodel'],coords)
    # Do not retain or read Hcheck for selection.
    del raw
    manifest=dict(config=CFG,python=platform.python_version(),numpy=np.__version__,cvxpy=cp.__version__,
      command='python -u run_validation.py',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      data_sha256=hashlib.sha256((OUT/'input_z0.mat').read_bytes()).hexdigest(),scope='offline/dev, no hardware')
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    allrows=[];summaries=[]
    for source,P in [('nonparam',Ptrain),('srivc',Pmodel)]:
        for taps in CFG['taps']:
            for delay in CFG['delays']:
                if args.pilot and (source,taps,delay)!=('nonparam',40,0):continue
                sm=run_candidate(source,taps,delay,f,P,Ptrain,coords);summaries.append(sm)
                allrows.append(dict(name=sm['name'],source=source,taps=taps,delay_samples=delay,**sm['train']))
                write_csv(OUT/'candidate_training_metrics.csv',allrows)
    if args.pilot:return
    selected={}
    for source in ['nonparam','srivc']:
        choices=[s for s in summaries if s['source']==source and s['train']['feasible']]
        if choices:
            best=min(choices,key=lambda s:(s['train']['score_max'],s['train']['score_sum']))
            if best['train']['score_max']<1:
                selected[source]={k:best[k] for k in ['name','source','taps','delay_samples','train']};continue
        selected[source]=dict(name='static_fallback',source=source,taps=1,delay_samples=0)
    (OUT/'selected_training_only.json').write_text(json.dumps(selected,indent=2),encoding='utf-8')
    evaluate(selected,f,coords)

if __name__=='__main__':main()
