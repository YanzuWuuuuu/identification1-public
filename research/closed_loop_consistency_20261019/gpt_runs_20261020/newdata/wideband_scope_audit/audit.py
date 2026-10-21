"""Read-only feasibility check for extending real-data validation to 100–800 Hz.
No modal refitting, change of original model, or replacement of A–E results.
"""
import os,sys,json,csv
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
sys.path.insert(0,str(PARENT))
from prepare_data import evalK,FS,P0

CTRS=('lowgain','nominal','highgain')
SIGNS=('pppp','pmpm','ppmm','pmmp')
BANDS=((3,25),(100,170),(250,400),(450,800),(100,800))

def right(a,b):return np.linalg.solve(b.swapaxes(-1,-2),a.swapaxes(-1,-2)).swapaxes(-1,-2)
def writecsv(name,rows):
    with (HERE/name).open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    kp=np.arange(18,4801);f=kp*FS/P0
    rows=[];peakrows=[];allg={}
    for ctrl in CTRS:
        blocks=[]
        for s in SIGNS:
            a=np.load(PARENT/'prepared'/f'C_{ctrl}_{s}.npz');p=int(a['P'])
            blocks.append({k:np.fft.rfft(a[k].reshape(p,P0,4),axis=1)[:,kp] for k in ('U','Y','D')})
        p=min(b['U'].shape[0] for b in blocks);n=p//2
        def means(key,sl):return np.stack([b[key][sl].mean(0) for b in blocks],axis=2)
        U=means('U',slice(0,p));Y=means('Y',slice(0,p));D=means('D',slice(0,p));G=right(Y,U)
        Uf=means('U',slice(0,n));Yf=means('Y',slice(0,n));Ub=means('U',slice(n,p));Yb=means('Y',slice(n,p))
        Gf=right(Yf,Uf);Gb=right(Yb,Ub)
        c=np.load(PARENT/'controllers'/f'{ctrl}.npz');K=evalK(f,c['A'],c['B'],c['T'])
        Su=np.linalg.inv(np.eye(4)[None]+K@G);svloop=np.linalg.svd(K@G,compute_uv=False)
        ucond=np.linalg.cond(U);dcond=np.linalg.cond(D)
        allg[ctrl]=G
        np.savez_compressed(HERE/f'{ctrl}_broadband.npz',f=f,G=G,G_front=Gf,G_back=Gb,U=U,Y=Y,D=D,
            loop_singular_values=svloop,Su=Su,P_common=p)
        for lo,hi in BANDS:
            ix=(f>=lo)&(f<=hi)
            # Entire raw complex matrices, uniform grid, no bin exclusions.
            ef=Yb-Gf@Ub;eb=Yf-Gb@Uf
            rows.append(dict(controller=ctrl,lo_hz=lo,hi_hz=hi,P_common=p,P_front=n,P_back=p-n,
                frequency_lines=int(ix.sum()),FRF_half_difference_percent=float(100*np.linalg.norm((Gf-Gb)[ix])/np.linalg.norm(G[ix])),
                inverse_free_cross_half_residual_percent=float(100*np.sqrt((np.linalg.norm(ef[ix])**2+np.linalg.norm(eb[ix])**2)/(np.linalg.norm(Yb[ix])**2+np.linalg.norm(Yf[ix])**2))),
                cond_U_median=float(np.median(ucond[ix])),cond_U_max=float(ucond[ix].max()),cond_D_max=float(dcond[ix].max()),
                loop_largest_singular_median=float(np.median(svloop[ix,0])),loop_largest_singular_p95=float(np.percentile(svloop[ix,0],95)),
                Su_minus_I_norm_median=float(np.median(np.linalg.norm(Su[ix]-np.eye(4),axis=(1,2)))),
                finite_all_bins=bool(np.isfinite(G[ix]).all())))
        for lo,hi in ((110,155),(270,300),(325,350),(350,380),(480,510),(650,700),(700,730),(760,800)):
            ix=np.flatnonzero((f>=lo)&(f<=hi));j=ix[np.argmax(np.linalg.norm(G[ix],axis=(1,2)))]
            peakrows.append(dict(controller=ctrl,lo_hz=lo,hi_hz=hi,FRF_norm_maximum_hz=float(f[j]),
                interpretation='window maximum of raw matrix FRF; not an estimated natural frequency or mode assignment'))
    cross=[]
    for ia,a in enumerate(CTRS):
        for b in CTRS[ia+1:]:
            for lo,hi in BANDS:
                ix=(f>=lo)&(f<=hi);ga=allg[a][ix];gb=allg[b][ix]
                cross.append(dict(controller_a=a,controller_b=b,lo_hz=lo,hi_hz=hi,
                    complex_FRF_difference_percent=float(200*np.linalg.norm(ga-gb)/np.linalg.norm(ga+gb)),
                    interpretation='different independent runs and random-phase excitations; descriptive, no uncertainty or invariance pass claim'))
    writecsv('band_readiness.csv',rows);writecsv('cross_controller.csv',cross);writecsv('window_peaks.csv',peakrows)
    # Exact check against the original completed 3–25 Hz calculation.
    check=[]
    for ctrl in CTRS:
        old=np.load(PARENT/'analysis'/f'{ctrl}_orthogonal_frf.npz');now=allg[ctrl][:len(old['f'])]
        err=float(np.linalg.norm(now-old['G'])/np.linalg.norm(old['G']));assert err<1e-12
        check.append(dict(controller=ctrl,old_band_relative_difference=err))
    manifest=dict(parent_claim='closed-loop modal measurement across controllers',question='Should validation include high-frequency flexible dynamics?',
        design='same prepared C records, same gate/discard rules and common P; raw 3–800 Hz orthogonal FRF only',
        resource_envelope='CPU numpy/scipy, one process, sequential FFT, no Monte Carlo refit, bounded memory <1 GB',
        paper_role='supporting',section_id='experimental validation',item_id='wideband-feasibility',
        fixed_conditions=['raw data read only','original A–E fit results unchanged','all raw bins retained','no LRM or inverse regularization'],
        comparability='half-record repeatability comparable by band; frequency-dependent signal-to-noise differs; no modal consistency claims',
        validation=check,status='completed descriptive audit',next_action='separate extended modal-model experiment with frozen model order and fair M1/M2/M3/M3w/SML comparisons')
    (HERE/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(rows,ensure_ascii=False,indent=2));print('CROSS',json.dumps(cross,ensure_ascii=False));print('PEAKS',json.dumps(peakrows,ensure_ascii=False))

if __name__=='__main__':main()
