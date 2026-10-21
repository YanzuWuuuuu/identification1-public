"""Final independent bookkeeping and percentage-error verification."""
import os,sys,json,hashlib
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
HERE=Path(__file__).resolve().parent
bt=np.r_[np.log(2*np.pi*np.array([10.6,11.7,16.1])),np.log([.03,.04,.05])]
checks={};maxerr=0.;counts={};failures={}
for name,folder,pattern,nexpected in [('D',HERE/'sim_d/results','*.json',1600),('E',HERE/'short_e','P*.json',360)]:
    files=sorted(folder.glob(pattern));counts[name]=len(files);checks[name+'_count']=len(files)==nexpected
    failures[name]=0
    for p in files:
        r=json.loads(p.read_text(encoding='utf-8'));failures[name]+=not r['success']
        if r.get('beta') is not None and r.get('error_percent') is not None:
            b=np.array(r['beta']);err=np.array(r['error_percent']);expected=100*np.expm1(b[96:102]-bt)
            if np.isfinite(err).all():maxerr=max(maxerr,float(np.max(abs(err-expected))))
checks['exp_minus_one_error_formula']=maxerr<1e-10
checks['task0_long_40']=len(list((HERE/'task0/long_results').glob('*.json')))==40
real=json.loads((HERE/'real_d/summary.json').read_text(encoding='utf-8'))
checks['real_180_contexts']=real['context_rows']==180
checks['real_140_unique']=real['unique_raw']==140
checks['original_1219_unchanged']=json.loads((HERE/'original_integrity_check.json').read_text())['unchanged']
checks['scale_gate']=json.loads((HERE/'sim_d/diagnostics/scale_parity_gate.json').read_text())['pass']
checks['matrix_gate']=json.loads((HERE/'sim_d/diagnostics/matrix_covariance_validation.json').read_text())['pass']
checks['moment_derivative_gate']=json.loads((HERE/'validation_e.json').read_text())['passed']
checks['report_exists']=(HERE/'后续实验D_E_报告.md').stat().st_size>10000
out={'checks':checks,'all_passed':all(checks.values()),'formal_counts':counts,'retained_numerical_failures':failures,'max_percentage_error_discrepancy':maxerr}
(HERE/'final_verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
manifest={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.rglob('*.py') if '__pycache__' not in p.parts}
(HERE/'new_code_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2));assert out['all_passed']
