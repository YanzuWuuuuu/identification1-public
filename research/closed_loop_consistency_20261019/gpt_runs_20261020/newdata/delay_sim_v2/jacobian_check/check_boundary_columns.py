"""Single deterministic Jacobian check including delay-register boundary columns."""
import check_jacobian as c
from check_selected_design import complete_rank_and_direct_check
import numpy as np
import scipy.linalg as la
import json, hashlib
from pathlib import Path

def main():
    before=json.loads((c.OUT/'source_hashes_before.json').read_text(encoding='utf8'))
    assert all(hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest().lower()==r['SHA256'].lower() for r in before)
    res=c.fixed_residues()
    with np.load(c.OLD/'mc/data/P20_r000.npz') as a:
        keep=a['k']%int(a['P'])==0
        f=a['f'][keep];D=a['D'][keep]/(int(a['P'])*30000)
    results=[]
    for fc in (60.,15.):
        for band,(lo,hi,active) in c.BANDS.items():
            ix=(f>=lo-1e-10)&(f<=hi+1e-10);ff=f[ix];dd=D[ix]
            gg=c.full_response(ff,res);K=c.controller(ff,fc=fc)
            U=np.linalg.solve(np.eye(4)[None]+K@gg,dd[:,:,None])[:,:,0]
            X,T,meta=c.matrices(ff,active,res,'exact_tail_free_bg','reference_output',U,transient_active=list(range(6)))
            poly=np.zeros((len(ff),4,12),complex);labels=[]
            zi=np.exp(-2j*np.pi*ff/c.FS)
            for o in range(4):
                for p in range(3):
                    poly[:,o,3*o+p]=zi**p;labels.append(f'delay_boundary_y{o+1}_zi{p}')
            poly=poly.reshape(-1,12)
            X=np.c_[X,np.r_[poly.real,poly.imag]]
            meta['labels']+=labels;meta['nuisance_count']+=12;meta['total_count']+=12
            info=c.project(X,T);verify=complete_rank_and_direct_check(X,T)
            main=next(t for t in info['thresholds'] if t['cut']==1e-12)
            verify['sequential_vs_direct_max_abs_error']=float(np.max(abs(np.asarray(main['all_other_parameters_delay_energy'])-verify['direct_all_other_columns_delay_energy'])))
            r=dict(controller_fc=fc,band=band,case='exact_tail_free_bg',observation='reference_output',
                all_six_modal_transients=True,delay_boundary_powers=[0,1,2],frequency_points=len(ff),
                **meta,**info,**verify)
            results.append(r);c.save(c.OUT/f'FC{int(fc)}_{band}_exact_tail_T6_boundary3.json',r)
            print(fc,band,'nparam',r['total_count'],'rank',r['joint_rank'],
                'minr',min(main['all_other_parameters_delay_energy']),
                'joint',main['joint_weakest_energy_ratio'],'cond',r['joint_column_normalized_condition'],flush=True)
    c.save(c.OUT/'boundary_design_results.json',dict(no_fit=True,no_simulation=True,no_random_generation=True,results=results))
    assert all(hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest().lower()==r['SHA256'].lower() for r in before)

if __name__=='__main__':main()
