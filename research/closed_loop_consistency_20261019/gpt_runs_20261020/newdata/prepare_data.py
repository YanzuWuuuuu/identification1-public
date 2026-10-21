import os
for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[n]='1'
import sys
sys.dont_write_bytecode=True
sys.stdout.reconfigure(encoding='utf-8')
from pathlib import Path
import struct,zlib,io,json,re,hashlib,time,shutil
import numpy as np
import scipy.io as sio
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
DATA=ROOT/'sweep/sweep_data_new'
OLD=HERE.parents[1]/'gpt_runs_20261019/followup_DE'
FS=5000.;P0=30000
def js(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=True),encoding='utf-8');tmp.replace(path)
def load_mat(path):
    b=path.read_bytes();p=128;out={};checks=[]
    if b[126:128]!=b'IM':raise ValueError('Unexpected MATLAB endian/header')
    while p+8<=len(b):
        typ,n=struct.unpack('<II',b[p:p+8]);p+=8;ch=b[p:p+n];p+=n
        if len(ch)!=n:raise ValueError('Truncated MAT block')
        if typ==15:
            de=zlib.decompressobj(-15);raw=de.decompress(ch[2:]);raw+=de.flush()
            if not de.eof:raise ValueError('Truncated raw deflate payload')
            try:zlib.decompress(ch);valid=True
            except zlib.error:valid=False
            a=sio.loadmat(io.BytesIO(b[:128]+raw));out.update({k:v for k,v in a.items() if not k.startswith('__')})
            checks.append(valid)
        else:raise ValueError(f'Unexpected block {typ}; refusing silent skip')
    if len(out)!=16:raise ValueError(f'Expected16 vars, got{len(out)}')
    return out,dict(blocks=len(checks),zlib_checksum_pass=sum(checks),sha256=hashlib.sha256(b).hexdigest())
def align(raw):
    names=[f'y_sensor_0__{i}_' for i in range(4)]+[f'controller_u_0__{i}_' for i in range(4)]+[f'disturb_input{i+1}' for i in range(4)]+[f'u_Motor{i+1}' for i in range(4)]
    ticks=[np.rint(raw[k][:,0]*FS).astype(np.int64) for k in names]
    for k,t in zip(names,ticks):
        if np.any(np.diff(t)!=1):raise ValueError(f'Nonuniform or duplicate time:{k}')
        if np.max(abs(raw[k][:,0]*FS-t))>1e-2:raise ValueError('Clock not on 5kHz grid within 0.01 sample')
    lo=max(t[0] for t in ticks);hi=min(t[-1] for t in ticks)
    vals=np.column_stack([raw[k][lo-t[0]:hi-t[0]+1,1] for k,t in zip(names,ticks)])
    if not np.isfinite(vals).all():raise ValueError('nonfinite samples')
    return dict(Y=vals[:,:4],C=vals[:,4:8],D=vals[:,8:12],U=vals[:,12:16],t=(lo+np.arange(len(vals)))/FS),dict(lengths={k:len(raw[k]) for k in names},common_length=len(vals),common_start_s=lo/FS,common_end_s=hi/FS)
def metadata(p):
    s=p.stem;ctrl=next(c for c in ('highgain','nominal','lowgain') if '_'+c+'_' in s);g=p.parent.name
    pos=s.split('_'+ctrl)[0][2:];z=0. if pos=='0' else float(pos.replace('plus','').replace('minus','-').replace('_','.'))
    gain=None if g=='B' else float(s.split('disturbgain_')[1].split('_motor')[0].replace('_','.'))
    signs=[1]*4
    if '_motor_' in s:signs=[1 if a=='Plus' else -1 for a in re.findall('Plus|Minus',s.split('1234_')[1])]
    rid=g+'_'+ctrl
    if g=='C':rid+='_'+''.join('p' if a>0 else 'm' for a in signs)
    if g=='D':rid+='_'+str(gain).replace('.','p')
    if g=='E':rid+='_'+str(z).replace('-','m').replace('.','p')
    return dict(id=rid,group=g,controller=ctrl,z_mm=z,disturb_gain=gain,signs=signs,source=str(p))
def design(Y,C,ix,na=8,nb=8):
    return np.column_stack([C[ix-i] for i in range(1,na+1)]+[Y[ix-j] for j in range(nb+1)])
def fit_controller(d):
    n=len(d['Y']);ix=np.arange(8,n//2)
    X=design(d['Y'],d['C'],ix);sc=np.linalg.norm(X,axis=0)
    coef,res,rank,s=np.linalg.lstsq(X/sc,d['C'][ix],rcond=1e-13);coef/=sc[:,None]
    A=np.array([coef[i*4:(i+1)*4].T for i in range(8)])
    B=np.array([coef[32+j*4:32+(j+1)*4].T for j in range(9)])
    ix2=np.arange(n//2,n);err=d['C'][ix2]-design(d['Y'],d['C'],ix2)@coef
    return A,B,coef,dict(order=[8,8],design_rank=int(rank),n_columns=X.shape[1],scaled_singular_values=s.tolist(),
                        validation_relative_rms=float(np.linalg.norm(err)/np.linalg.norm(d['C'][ix2])),validation_max_absolute=float(np.max(abs(err))))
def evalK(f,A,B,T):
    zi=np.exp(-2j*np.pi*np.asarray(f)/FS)
    aa=np.eye(4)[None]-sum(A[i][None]*zi[:,None,None]**(i+1) for i in range(8))
    bb=sum(B[j][None]*zi[:,None,None]**j for j in range(9))
    return -T[None]@np.linalg.solve(aa,bb)
def main():
    for name in ('prepared','controllers','logs','baseline','fits','spectra','analysis','figures'): (HERE/name).mkdir(exist_ok=True)
    srcs=[(OLD/'real_d'/f,f) for f in ('realfit_original_copy.py','realfit_extended.py')]+[(OLD/f,f) for f in ('weighted_core.py','sml_core_copy.py')]
    manifest=[]
    for p,name in srcs:
        target=HERE/'baseline'/name
        if target.exists() and target.read_bytes()!=p.read_bytes():raise ValueError('Baseline copy mismatch')
        if not target.exists():shutil.copyfile(p,target)
        manifest.append(dict(source=str(p),copy=str(target),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    delay=ROOT/'experiment_data/srivc_no_lrm_20261015/delay_identification_20261016/effective_delay_matrix_us.csv'
    np.save(HERE/'baseline/tau_seconds.npy',np.loadtxt(delay,delimiter=',')*1e-6-1/FS/2)
    manifest.append(dict(source=str(delay),sha256=hashlib.sha256(delay.read_bytes()).hexdigest()))
    js(HERE/'baseline/sources.json',manifest)
    Ts=[];ctrls={};bmeta={}
    for p in sorted((DATA/'B').glob('*.mat')):
        raw,ch=load_mat(p);d,al=align(raw);meta=metadata(p)
        T,_,rank,s=np.linalg.lstsq(d['C'],d['U'],rcond=1e-12);T=T.T
        err=d['U']-d['C']@T.T;Ts.append(T)
        A,B,coef,info=fit_controller(d)
        info.update(T=T.tolist(),input_C_rank=int(rank),T_fit_relative=float(np.linalg.norm(err)/np.linalg.norm(d['U'])),
                    U_equals_C_relative=float(np.linalg.norm(d['U']-d['C'])/np.linalg.norm(d['U'])),
                    T_note='minimum-norm map; if C rank<4 this is not a uniquely identified physical matrix')
        np.savez_compressed(HERE/'controllers'/f"{meta['controller']}.npz",A=A,B=B,T=T,coef=coef)
        js(HERE/'controllers'/f"{meta['controller']}.json",info);ctrls[meta['controller']]=(T,coef)
        bmeta[meta['controller']]=info
        print('CONTROLLER',meta['controller'],info['validation_relative_rms'],'rankC',rank,flush=True)
    rows=[];seen={}
    for p in sorted(DATA.glob('*/*.mat')):
        raw,ch=load_mat(p);d,al=align(raw);meta=metadata(p);meta.update(ch);meta.update(al)
        rid=meta['id'];ctrl=meta['controller'];T,coef=ctrls[ctrl]
        if ch['sha256'] in seen:
            meta['duplicate_of']=seen[ch['sha256']];rows.append(meta);print('DUPLICATE',rid,meta['duplicate_of'],flush=True);continue
        seen[ch['sha256']]=rid
        ix=np.arange(8,len(d['Y']));er=d['C'][ix]-design(d['Y'],d['C'],ix)@coef
        meta['controller_record_relative_rms']=float(np.linalg.norm(er)/np.linalg.norm(d['C'][ix]))
        offset=er[:len(er)//2].mean(0)
        meta['controller_constant_offset']=offset.tolist()
        meta['controller_validation_after_train_offset_relative']=float(np.linalg.norm(er[len(er)//2:]-offset)/np.linalg.norm(d['C'][ix[len(er)//2:]]))
        effective=d['U']-d['C']@T.T
        e1=np.linalg.norm(effective-d['D'],axis=1);e0=np.linalg.norm(effective,axis=1);dd=np.linalg.norm(d['D'],axis=1)
        active=(e1<1e-5*np.maximum(dd,1e-4))&(e0>.1*dd)
        edges=np.flatnonzero(np.diff(np.r_[False,active,False]));runs=list(zip(edges[::2],edges[1::2]))
        runs=[(int(a),int(b)) for a,b in runs if b-a>100]
        meta['injection_runs']=runs
        warnings=[]
        if meta['group']=='B':
            start=0;P=len(d['U'])//P0;stop=P*P0;meta['injection_status']='noise_only'
            meta['stored_D_is_applied']=False
        else:
            if not runs:raise ValueError(f'No injection {rid}')
            a,b=max(runs,key=lambda ab:ab[1]-ab[0]);start=a+P0;P=(b-start)//P0;stop=start+P*P0
            drops=1;drift_history=[]
            while P>=4:
                yy=d['Y'][start:stop].reshape(P,P0,4);yy=yy-yy.mean(1,keepdims=True)
                change=np.linalg.norm(np.diff(yy,axis=0),axis=(1,2))/np.maximum(np.linalg.norm(yy[:-1],axis=(1,2)),1e-300)
                drift_history.append(dict(discarded=drops,adjacent_relative=change.tolist()))
                if change[0]<=2*np.median(change[1:]):break
                start+=P0;P-=1;drops+=1
            meta.update(injection_start_s=float(d['t'][a]),injection_stop_s=float(d['t'][b-1]+1/FS),injection_samples=b-a,
                        injection_relative_rms=float(np.linalg.norm(effective[a:b]-d['D'][a:b])/np.linalg.norm(d['D'][a:b])),discarded_initial_periods=drops,transient_check_history=drift_history)
            if P<2:warnings.append('fewer_than_2_periods')
        seg={k:v[start:stop] for k,v in d.items()}
        if meta['group']=='B':seg['D_stored']=seg['D'].copy();seg['D']=np.zeros_like(seg['D'])
        yp=seg['Y'].reshape(P,P0,4);dp=seg['D'].reshape(P,P0,4)
        dy=np.linalg.norm(np.diff(yp,axis=0),axis=(1,2));ys=np.linalg.norm(yp,axis=(1,2))
        # Also remove each period's mean so position offsets cannot mask drift.
        yac=yp-yp.mean(1,keepdims=True)
        dac=np.linalg.norm(np.diff(yac,axis=0),axis=(1,2))/np.maximum(np.linalg.norm(yac[:-1],axis=(1,2)),1e-300)
        meta.update(P=P,P0=P0,fs=FS,start_index=start,stop_index=stop,usable_start_s=float(seg['t'][0]),usable_stop_s=float(seg['t'][-1]+1/FS),
                    D_repeat_relative=float(np.linalg.norm(dp-dp.mean(0))/max(np.linalg.norm(dp),1e-300)),
                    adjacent_Y_ac_relative=dac.tolist(),Y_period_mean=yp.mean(1).tolist(),
                    U_min=seg['U'].min(0).tolist(),U_max=seg['U'].max(0).tolist(),Y_min=seg['Y'].min(0).tolist(),Y_max=seg['Y'].max(0).tolist())
        # Clipping flag is diagnostic only; true hardware limits are not supplied.
        meta['U_extreme_repeat_fraction']=[float(max(np.mean(seg['U'][:,j]==seg['U'][:,j].min()),np.mean(seg['U'][:,j]==seg['U'][:,j].max()))) for j in range(4)]
        if len(dac)>2 and dac[0]>2*np.median(dac[1:]):warnings.append('first_retained_period_difference_exceeds_twice_later_median')
        if meta['controller_validation_after_train_offset_relative']>1e-8:warnings.append('controller_validation_after_offset_above_1e-8')
        if meta['D_repeat_relative']>1e-8:warnings.append('reference_not_exactly_periodic')
        meta['warnings']=warnings;meta['hardware_range_verified']=False
        np.savez_compressed(HERE/'prepared'/f'{rid}.npz',**seg,P=P,P0=P0,fs=FS,T=T,start_index=start,stop_index=stop)
        js(HERE/'prepared'/f'{rid}.json',meta);rows.append(meta)
        print('PREP',rid,'P',P,'ctrl',meta['controller_record_relative_rms'],'warn',warnings,flush=True)
    js(HERE/'data_inventory.json',rows)
    js(HERE/'controllers/comparison.json',dict(max_T_difference=float(np.max(np.abs(np.asarray(Ts)-Ts[0]))),controllers=bmeta))
if __name__=='__main__':main()
