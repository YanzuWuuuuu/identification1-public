"""Read-only signal/timestamp and controller recurrence audit; no plant fits."""
import os
for _key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_key] = "1"
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import io, json, struct, zlib, hashlib, csv, time
import numpy as np
import scipy.io as sio
import scipy.linalg as sla

OUT = Path(__file__).resolve().parent
NEW = OUT.parents[1]
ROOT = NEW.parents[3]
sys.path.insert(0, str(NEW))
import prepare_data as prep

def save(name, obj):
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    (OUT/name).write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=convert), encoding="utf8")

def read_old(path):
    source = path.read_bytes()
    pos, arrays, checks = 128, {}, []
    while pos + 8 <= len(source):
        typ, n = struct.unpack("<II", source[pos:pos+8])
        assert typ == 15
        data = source[pos+8:pos+8+n]
        decoder = zlib.decompressobj(-15)
        raw = decoder.decompress(data[2:]) + decoder.flush()
        assert decoder.eof
        arrays.update({k:v for k,v in sio.loadmat(io.BytesIO(source[:128]+raw)).items() if not k.startswith("__")})
        try:
            zlib.decompress(data); checks.append(True)
        except zlib.error:
            checks.append(False)
        pos += n+8
    assert pos == len(source)
    return arrays, dict(source=str(path), sha256=hashlib.sha256(source).hexdigest(),
                        variables=len(arrays), valid_zlib_checksums=sum(checks))

def lag_comparison(U, C, T, lags=(-2,-1,0,1,2)):
    ix=np.arange(2,len(U)-2)
    rows=[]
    for lag in lags:
        e=U[ix]-C[ix-lag]@T.T
        rows.append(dict(lag_samples=lag, definition="U[t] - T C[t-lag]",
                         relative_rms=np.linalg.norm(e)/np.linalg.norm(U[ix]),
                         per_channel_relative_rms=np.linalg.norm(e,axis=0)/np.linalg.norm(U[ix],axis=0)))
    return rows

def arx(Y, C, nk, na=8, nb=8):
    """C[t]=sum A_i C[t-i] + sum B_j Y[t-(nk+j)], no intercept."""
    n=len(Y); margin=max(na,nk+nb,0)+2
    train=np.arange(margin,n//2)
    val=np.arange(n//2,n-max(0,-nk)-2)
    def design(ix):
        return np.column_stack([C[ix-i] for i in range(1,na+1)]+[Y[ix-nk-j] for j in range(nb+1)])
    X=design(train); norms=np.linalg.norm(X,axis=0)
    coef,_,rank,sv=sla.lstsq(X/norms,C[train],cond=1e-13,lapack_driver="gelsd")
    coef=coef/norms[:,None]
    e=C[val]-design(val)@coef
    nout=C.shape[1]
    B=np.array([coef[nout*na+4*j:nout*na+4*(j+1)].T for j in range(nb+1)])
    centered=C[val]-C[val].mean(0)
    result=dict(nk=nk,na=na,nb=nb,output_dimension=nout,input_dimension=4,
                input_lags=list(range(nk,nk+nb+1)),training_samples=len(train),validation_samples=len(val),
                rank=int(rank),columns=X.shape[1],singular_values=sv,
                validation_relative_rms=float(np.linalg.norm(e)/np.linalg.norm(C[val])),
                validation_relative_to_centered_output=float(np.linalg.norm(e)/np.linalg.norm(centered)),
                validation_per_channel_relative_rms=np.linalg.norm(e,axis=0)/np.linalg.norm(C[val],axis=0),
                validation_max_abs=float(abs(e).max()),
                input_tap_fro_norms=np.linalg.norm(B,axis=(1,2)),
                exact_to_1e8=bool(np.linalg.norm(e)/np.linalg.norm(C[val])<1e-8))
    return result,coef

def main():
    started=time.perf_counter(); allout=[]; flattened=[]
    for ctrl in ("lowgain","nominal","highgain"):
        path=ROOT/"sweep/sweep_data_new/B"/f"z_0_{ctrl}_B_noise.mat"
        raw,source=prep.load_mat(path); dat,alignment=prep.align(raw)
        saved=np.load(NEW/"controllers"/f"{ctrl}.npz")
        rec=dict(id="new_B_"+ctrl,kind="new",source=source,alignment=alignment,
                 raw_U_C_exact_equal=bool(np.array_equal(dat["U"],dat["C"])),
                 raw_U_C_max_abs=float(abs(dat["U"]-dat["C"]).max()),
                 C_rank=int(np.linalg.matrix_rank(dat["C"])),
                 T=saved["T"],
                 lag_comparison=lag_comparison(dat["U"],dat["C"],saved["T"]),
                 time_all_identical=all(np.array_equal(next(iter(raw.values()))[:,0],v[:,0]) for v in raw.values()),
                 timestamp_roundoff_max_samples=max(float(abs(v[:,0]*5000-np.rint(v[:,0]*5000)).max()) for v in raw.values()),
                 arx=[])
        for nk in (-1,0,1,2,3):
            fit,coef=arx(dat["Y"],dat["C"],nk)
            rec["arx"].append(fit)
            np.savez_compressed(OUT/f'{rec["id"]}_nk{nk}_arx.npz',coef=coef)
            flattened.append(dict(id=rec["id"],nk=nk,na=8,nb=8,relative_rms=fit["validation_relative_rms"],
                                  ac_relative_rms=fit["validation_relative_to_centered_output"],rank=fit["rank"],
                                  columns=fit["columns"],exact=fit["exact_to_1e8"]))
            print(rec["id"],nk,fit["validation_relative_rms"],flush=True)
        rec["largest_input_lag_retaining_exactness"]=max(v["nk"] for v in rec["arx"] if v["exact_to_1e8"])
        save(rec["id"]+".json",rec);allout.append(rec)
    oldraw,oldmeta=read_old(ROOT/"experiment_data/auto_identify_motor1234_z0.mat")
    oldcl=sio.loadmat(ROOT/"experiment_data/srivc_no_lrm_20261015/closed_loop_consistency_20261019/cl_data_z0.mat")
    time0=oldraw["u_Motor1"][:,0]
    oldU=np.column_stack([oldraw[f"u_Motor{i}"][:,1] for i in range(1,5)])
    oldY=np.column_stack([oldraw[f"pre_y_{i}__0_"][:,1] for i in range(4)])
    oldV=np.column_stack([oldraw[k][:,1] for k in ("u_z","u_Rx","u_Ry")])
    quiet=oldU[:,3]==0
    Q=oldcl["Q"];qerr=oldU[quiet]-oldV[quiet]@Q.T
    rec=dict(id="old_z0",kind="old",source=oldmeta,
             time_all_identical=all(np.array_equal(time0,v[:,0]) for v in oldraw.values()),
             timestamp_roundoff_max_samples=float(abs(time0*5000-np.rint(time0*5000)).max()),
             processed_U_equals_raw=bool(np.array_equal(oldU,oldcl["U"])),
             processed_Y_equals_raw=bool(np.array_equal(oldY,oldcl["Y"])),
             processed_V_equals_raw=bool(np.array_equal(oldV,oldcl["V"])),
             processed_U_raw_max_abs=float(abs(oldU-oldcl["U"]).max()),
             processed_Y_raw_max_abs=float(abs(oldY-oldcl["Y"]).max()),
             processed_V_raw_max_abs=float(abs(oldV-oldcl["V"]).max()),
             quiet_count=int(quiet.sum()),quiet_U_QV_relative=float(np.linalg.norm(qerr)/np.linalg.norm(oldU[quiet])),
             quiet_U_QV_max_abs=float(abs(qerr).max()),arx=[])
    for na,nb in ((6,6),(8,8)):
        for nk in (-1,0,1,2,3):
            fit,coef=arx(oldY,oldV,nk,na,nb)
            rec["arx"].append(fit)
            np.savez_compressed(OUT/f'old_z0_nk{nk}_na{na}_arx.npz',coef=coef)
            flattened.append(dict(id=rec["id"],nk=nk,na=na,nb=nb,relative_rms=fit["validation_relative_rms"],
                                  ac_relative_rms=fit["validation_relative_to_centered_output"],rank=fit["rank"],
                                  columns=fit["columns"],exact=fit["exact_to_1e8"]))
            print(rec["id"],na,nk,fit["validation_relative_rms"],flush=True)
    rec["largest_input_lag_retaining_exactness_by_order"]={str(na):max((v["nk"] for v in rec["arx"] if v["na"]==na and v["exact_to_1e8"]),default=None) for na in (6,8)}
    save("old_z0.json",rec);allout.append(rec)
    with (OUT/"arx_delay_scan.csv").open("w",encoding="utf-8-sig",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(flattened[0]));wr.writeheader();wr.writerows(flattened)
    save("alignment_summary.json",dict(records=allout,seconds=time.perf_counter()-started,
        distinction="These are identities/recurrences between recorded signals. No plant model was fit; ADC/DAC and hardware block locations are not identified."))

if __name__=="__main__": main()
