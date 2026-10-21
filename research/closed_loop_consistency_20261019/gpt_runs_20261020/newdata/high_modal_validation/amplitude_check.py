"""Secondary D-amplitude check; same frozen local model, reference IV only."""
import run as r
from pathlib import Path
import json,csv

def main():
    root=r.HERE;r.HERE=root/'amplitude_diagnostic'
    for sub in ('raw','spectra'): (r.HERE/sub).mkdir(parents=True,exist_ok=True)
    ids=['D_nominal_0p6','D_nominal_0p8','D_nominal_1p0','D_nominal_1p2']
    rows=[]
    for w in ('H300','H493'):
        for rid in ids:
            r.run_one(rid,w,('M3',))
            a=json.loads((r.HERE/'raw'/f'{rid}_{w}_M3.json').read_text())
            for j,(f,z) in enumerate(zip(a['frequency_hz'],a['damping_percent'])):
                rows.append(dict(id=rid,window=w,component=j+1,disturb_gain=float(rid.split('_')[-1].replace('p','.')),
                    f_hz=f,zeta_percent=z,converged=a['converged'],near_bound=a['near_bound'],strict_score_pass=a['strict_score_pass'],score=a['score']))
    with (root/'analysis/amplitude_diagnostic.csv').open('w',encoding='utf-8-sig',newline='') as h:
        dw=csv.DictWriter(h,fieldnames=list(rows[0]));dw.writeheader();dw.writerows(rows)
    r.js(r.HERE/'scope.json',dict(purpose='Same-controller amplitude check, secondary to fixed 120 comparisons',methods=['M3'],windows=['H300','H493'],records=ids,primary_results_replaced=False))

if __name__=='__main__':main()
