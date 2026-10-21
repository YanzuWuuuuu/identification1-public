"""Deterministic no-noise Jacobian audit only: no fit, MC, RNG or simulation."""
import os, sys, ast, json, hashlib
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
sys.dont_write_bytecode=True
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import scipy.linalg as la
from scipy import signal

OUT=Path(__file__).resolve().parent
NEW=OUT.parents[1]
RESEARCH=NEW.parents[1]
OLD=NEW/'delay_sim_closeout'
sys.path.insert(0,str(NEW/'delay_check/baseline'))
from modal_model import ModalModel
FS=5000.; UNIT=1e-4
FM=np.array([10.6,11.7,16.1,286.573,344.945,490.39])
ZM=np.array([.03,.04,.05,.0119,.0132,.0186])
NI=np.array([0,1,0,1]); NO=np.array([1,0,1,0])
TAU=(NO[:,None]+NI[None,:])/FS
HELM=la.helmert(4,full=False).T
DESIGN=np.zeros((4,4,7));DESIGN[:,:,0]=1
DESIGN[:,:,1:4]=HELM[:,None,:];DESIGN[:,:,4:]=HELM[None,:,:]
BANDS={'L':(3.,100.,[0,1,2]),'H':(100.,550.,[3,4,5]),'W':(3.,550.,list(range(6)))}
CUTS=(1e-8,1e-10,1e-12,1e-14)

def enc(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x)
    raise TypeError(type(x).__name__)

def save(p,x):
    Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=enc),encoding='utf8')

def fixed_residues():
    p=OLD/'stable_profile_pilot/v1_P10_r000_H_M3/protocol.json'
    b=np.asarray(json.loads(p.read_text(encoding='utf8'))['initial_beta'])
    # This saved initializer only perturbed log poles and zeroed delay offsets.
    # Its linear C0 entries are exactly the original fixed plant residues.
    return np.array([[b[o*88+10*j:o*88+10*j+4] for o in range(4)] for j in range(6)])

def controller(f,fc=60.,flp=120.):
    # Extract only the deterministic original controller function. Importing
    # cl2 itself would execute its top-level random residue generation.
    tree=ast.parse((RESEARCH/'code/cl2.py').read_text(encoding='utf8'))
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='controller')
    xa=np.array([.2,-.2,-.2,.2]);ya=np.array([.15,.15,-.15,-.15])
    xs=np.array([.18,-.18,-.18,.18]);ys=np.array([.13,.13,-.13,-.13])
    Br=np.vstack([np.ones(4)/5.,ya/.05,-xa/.08])
    Cr=np.vstack([np.ones(4),ys,-xs]).T
    scope=dict(np=np,sla=la,ss=signal,os=os,Ts=1/FS,Br=Br,Cr=Cr,FLP=flp)
    exec(compile(ast.Module(body=[node],type_ignores=[]),'original_controller_only','exec'),scope)
    A,B,C,D=scope['controller'](fc=fc)
    z=np.exp(2j*np.pi*f/FS)
    return np.array([C@la.solve(zz*np.eye(len(A))-A,B)+D for zz in z])

def local_model(f,active,res):
    m=ModalModel(len(active),TAU,float(f[-1]),fs=FS)
    b=np.zeros(m.nb)
    for o in range(4):
        for j,i in enumerate(active):b[o*m.L+10*j:o*m.L+10*j+4]=res[i,o]
    b[m.iw:m.iz]=np.log(2*np.pi*FM[active]);b[m.iz:]=np.log(ZM[active])
    c=m.context(f)
    return m,b,c,m.blocks(c,b,True)

def full_response(f,res):
    m,b,c,_=local_model(f,list(range(6)),res)
    return m.resp(c,b)

def matrices(f,active,res,case,kind,U,transient_active=None):
    m,b,c,hs=local_model(f,active,res)
    bg=case!='exact_tail_no_bg'
    cols=16*len(active)+(32 if bg else 0)+2*len(active)
    J=np.zeros((len(f),4,4,cols),complex); labels=[];p=0
    for j,i in enumerate(active):
        h0=hs[j][0]
        for o in range(4):
            for a in range(4):
                J[:,o,a,p]=c['E'][:,o,a]*h0
                labels.append(f'C0_mode{i+1}_y{o+1}_u{a+1}');p+=1
    if bg:
        for name,h in [('D0',np.ones(len(f))),('D2',(c['s']/m.wr)**2)]:
            for o in range(4):
                for a in range(4):
                    J[:,o,a,p]=c['E'][:,o,a]*h
                    labels.append(f'{name}_y{o+1}_u{a+1}');p+=1
    for q,qn in enumerate(('logw','logz')):
        for j,i in enumerate(active):
            J[:,:,:,p]=c['E']*hs[j][3][q][0][:,None,None]*res[i][None]
            labels.append(f'{qn}_mode{i+1}');p+=1
    G_local=m.resp(c,b);G_true=full_response(f,res)
    G_model=G_local if case=='local_free_bg' else G_true
    Z=-c['s'][:,None,None,None]*UNIT*G_model[:,:,:,None]*DESIGN[None]
    if kind=='FRF':
        X=J.reshape(-1,cols);T=Z.reshape(-1,7)
        missing=la.norm(G_true-G_model)/max(la.norm(G_true),1e-300)
    else:
        X=np.einsum('koip,ki->kop',J,U).reshape(-1,cols)
        T=np.einsum('koip,ki->kop',Z,U).reshape(-1,7)
        ta=active if transient_active is None else transient_active
        ths=hs if ta==active else local_model(f,ta,res)[3]
        tr=np.zeros((len(f),4,8*len(ta)),complex);it=0
        for j,i in enumerate(ta):
            den=ths[j][2]
            for o in range(4):
                tr[:,o,it]=1/den;tr[:,o,it+1]=c['zi']/den;it+=2
                labels.extend([f'transient0_mode{i+1}_y{o+1}',f'transient1_mode{i+1}_y{o+1}'])
        X=np.c_[X,tr.reshape(-1,tr.shape[2])]
        missing=la.norm(np.einsum('koi,ki->ko',G_true-G_model,U))/max(
            la.norm(np.einsum('koi,ki->ko',G_true,U)),1e-300)
    # Analytic delay derivative is independently checked without fitting.
    direction=np.array([1.,.3,-.2,.1,-.1,.2,-.3])
    dt=UNIT*np.einsum('oip,p->oi',DESIGN,direction);h=1e-4
    fd=(np.exp(-c['s'][:,None,None]*h*dt)*G_model-
        np.exp(c['s'][:,None,None]*h*dt)*G_model)/(2*h)
    analytic=np.einsum('koip,p->koi',Z,direction)
    check=float(la.norm(fd-analytic)/max(la.norm(analytic),1e-300))
    return np.r_[X.real,X.imag],np.r_[T.real,T.imag],dict(
        labels=labels,nuisance_count=len(labels),total_count=len(labels)+7,
        exact_six_mode_steady_mean_represented=case!='local_free_bg' or len(active)==6,
        missing_model_relative_norm=float(missing),delay_derivative_relative_error=check,
        background_at_evaluation='D0=D2=0',free_background=bg,
        transient_columns_included=kind!='FRF',
        transient_modes=(active if transient_active is None else transient_active) if kind!='FRF' else [],
        FRF_note='Transfer function has no finite-record transient parameters; all modal transients are included in reference-output check.')

def project(X,T):
    xn=la.norm(X,axis=0);tn=la.norm(T,axis=0)
    X=X/np.maximum(xn,1e-300);T=T/np.maximum(tn,1e-300)
    U,s,_=la.svd(X,full_matrices=False,check_finite=False)
    Qt,Rt=la.qr(T,mode='economic',check_finite=False)
    results=[]
    for cut in CUTS:
        keep=s>cut*s[0];Q=U[:,keep]
        residual=T-Q@(Q.T@T)
        ratios=np.sum(residual**2,axis=0)
        # Condition on the other six delay columns as requested, after the
        # nuisance projection; use an absolute tolerance on original unit columns.
        conditional=[];other_ranks=[]
        for j in range(7):
            ii=[k for k in range(7) if k!=j]
            Ur,sr,_=la.svd(residual[:,ii],full_matrices=False,check_finite=False)
            kk=sr>cut*max(s[0],la.norm(T,2))
            v=residual[:,j]-Ur[:,kk]@(Ur[:,kk].T@residual[:,j])
            conditional.append(float(v@v));other_ranks.append(int(sum(kk)))
        Rq=Qt-Q@(Q.T@Qt)
        _,sq,Vh=la.svd(Rq,full_matrices=False,check_finite=False)
        v=la.solve_triangular(Rt,Vh[-1],check_finite=False)/tn
        v=v/la.norm(v)
        matrix=np.einsum('oip,p->oi',DESIGN,v)
        matrix=matrix/max(abs(matrix).max(),1e-300)
        results.append(dict(cut=cut,nuisance_rank=int(sum(keep)),
            nuisance_only_delay_energy=ratios,
            all_other_parameters_delay_energy=conditional,
            other_delay_projected_rank=other_ranks,
            joint_weakest_energy_ratio=float(sq[-1]**2),
            joint_projected_orthonormal_delay_singular_values=sq,
            weakest_parameter_direction=v,
            weakest_delay_matrix_normalized_max_one=matrix,
            all_seven_conditional_at_least_1e3=bool(min(conditional)>=1e-3)))
    return dict(column_norm_nuisance=xn,column_norm_delay=tn,
                nuisance_normalized_singular_values=s,thresholds=results)

def main():
    before=json.loads((OUT/'source_hashes_before.json').read_text(encoding='utf8'))
    for r in before:
        assert hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest().lower()==r['SHA256'].lower()
    res=fixed_residues()
    data_path=OLD/'mc/data/P20_r000.npz'
    with np.load(data_path) as a:
        keep=a['k']%int(a['P'])==0
        f=a['f'][keep];D=a['D'][keep]/(int(a['P'])*30000)
    result=dict(date_utc=datetime.now(timezone.utc).isoformat(),
        no_fit=True,no_simulation=True,no_random_generation=True,
        fixed_residues=res,frequencies_hz=FM,damping=ZM,
        delay_seconds=TAU,reference_source=data_path,
        controller=dict(FC=60,FLP=120,FSC=3),threshold=1e-3,
        main_projection_cut=1e-12,
        note='Known out-of-band mechanical terms remain inside the same unknown delay E(tau); their residues and poles alone are fixed.',
        cases=[])
    for band,(lo,hi,active) in BANDS.items():
        ix=(f>=lo-1e-10)&(f<=hi+1e-10);ff=f[ix];dd=D[ix]
        gg=full_response(ff,res);K=controller(ff)
        U=np.linalg.solve(np.eye(4)[None]+K@gg,dd[:,:,None])[:,:,0]
        for case in ('local_free_bg','exact_tail_free_bg','exact_tail_no_bg'):
            for kind in ('FRF','reference_output'):
                X,T,meta=matrices(ff,active,res,case,kind,U)
                info=project(X,T)
                r=dict(band=band,band_hz=[lo,hi],active_modes=[i+1 for i in active],
                    case=case,observation=kind,frequency_points=len(ff),**meta,**info)
                result['cases'].append(r)
                save(OUT/f'{band}_{case}_{kind}.json',r)
                main=next(t for t in r['thresholds'] if t['cut']==1e-12)
                print(band,case,kind,'nuisance',np.min(main['nuisance_only_delay_energy']),
                      'conditional',min(main['all_other_parameters_delay_energy']),
                      'joint',main['joint_weakest_energy_ratio'],flush=True)
    save(OUT/'jacobian_results.json',result)
    after=[]
    for r in before:
        h=hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest()
        after.append(dict(path=r['Path'],unchanged=h.lower()==r['SHA256'].lower(),sha256=h))
    save(OUT/'source_hash_verification.json',after)
    assert all(r['unchanged'] for r in after)

if __name__=='__main__':main()
