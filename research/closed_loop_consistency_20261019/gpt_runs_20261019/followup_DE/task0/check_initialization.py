"""Supplementary root check on the original P1 data; primary MC is unchanged."""
import sys
sys.dont_write_bytecode=True
from run_task0 import *

def main():
    bpop=np.array(json.loads((HERE/'population'/'P1_short_analytic.json').read_text(encoding='utf-8'))['beta'])
    rows=[]
    for r in range(40):
        a=np.load(OLD/'data'/f'P01_r{r:03d}.npz'); z,U,Y=a['zk'],a['U'],a['Y']
        old=json.loads((OLD/'results'/f'P01_r{r:03d}_dir_band_W1.json').read_text(encoding='utf-8'))
        row={'replicate':r,'original_error_percent':old['error_percent'],'status':'exception'}
        try:
            W=np.ones_like(U.real); b0=bpop.copy(); b0[:c.NLIN]=0; b0[c.ITR:]=0
            b0=c.init_lin(z,U,Y,W,b0)
            b,it=c.iv_solve(z,U,Y,U,W,b0,it=80,tol=1e-6)
            e=100*np.expm1(b[c.IW:c.ITR]-c.truth_beta()[c.IW:c.ITR])
            sc,rr=score(z,U,Y,b)
            row.update(status='pole_step_converged' if it<79 else 'iteration_limit',iterations=it+1,error_percent=e.tolist(),delta_percentage_points=(e-np.array(old['error_percent'])).tolist(),score_max=sc,residual=rr)
        except Exception as exc: row['exception']=repr(exc)
        rows.append(row); atomic(HERE/'initialization_results'/f'r{r:03d}.json',row)
    e=np.array([a['error_percent'] for a in rows if 'error_percent' in a]); dif=np.array([a['delta_percentage_points'] for a in rows if 'error_percent' in a])
    out={'description':'Same cached P1 data, initial poles at direct population pseudo-true values, linear coefficients reinitialized, tolerance tightened 1e-5 to 1e-6; diagnostic only.','n':len(rows),'failures':sum(a['status']!='pole_step_converged' for a in rows),'mean_error_percent':e.mean(0).tolist(),'max_abs_change_percentage_points':np.abs(dif).max(0).tolist(),'mean_change_percentage_points':dif.mean(0).tolist(),'rows':rows}
    atomic(HERE/'initialization_check.json',out); print(json.dumps({k:v for k,v in out.items() if k!='rows'}),flush=True)
if __name__=='__main__': main()
