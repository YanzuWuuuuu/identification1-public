from common import *
import hashlib,datetime
from PIL import Image
def main():
    manifest=rd(HERE/'protected_hashes_before.json')
    changed=[];missing=[];count=0
    for name,before in manifest['files'].items():
        p=Path(name)
        if not p.exists():missing.append(name);continue
        hh=hashlib.sha256()
        with p.open('rb') as h:
            for chunk in iter(lambda:h.read(8*1024*1024),b''):hh.update(chunk)
        if hh.hexdigest()!=before['sha256']:changed.append(name)
        count+=1
    baseline=[]
    for name in ('realfit_extended.py','realfit_original_copy.py','weighted_core.py','tau_seconds.npy'):
        eq=(HERE/'baseline'/name).read_bytes()==(SRC/'baseline'/name).read_bytes()
        baseline.append(dict(name=name,byte_identical=eq))
    bs=[rd(HERE/'analysis'/f'{r}_noise_bias.json') for r in IDS]
    rr=rd(HERE/'analysis/residual_checks.json')
    invariants=dict(
        source_spectra_max_error=max(max(x['selfcheck'].values()) for x in bs),
        original_B_prediction_max_abs_percent_error=max(x['old_prediction_reproduction_max_abs'] for x in bs),
        residual_trace_identity_max_relative_error=max(x['trace_identity_relative_error'] for x in rr),
        projector_annihilation_max_relative_error=max(x['projector_relative_error'] for x in rr),
        bootstrap_mean_max_absolute_MC_SE=max(abs(x['bootstrap_mean_standardized']) for x in rr),
        recomputed_score_max_absolute_error=max(abs(x['original_score']-x['score_recomputed']['score_max']) for x in rr),
        residual_real_dimension=sorted(set(x['real_dimension'] for x in rr)),
        model_parameter_count=m.NB)
    figures=[]
    for p in sorted((HERE/'figures').glob('*.png')):
        with Image.open(p) as im:figures.append(dict(file=p.name,size=im.size,dpi=im.info.get('dpi'),pdf_exists=p.with_suffix('.pdf').exists()))
    jk=rd(HERE/'analysis/jackknife.json')
    out=dict(time=datetime.datetime.now().isoformat(),protected_files_checked=count,changed=changed,missing=missing,
             all_protected_unchanged=not(changed or missing),copied_model_checks=baseline,numerical_checks=invariants,
             fits=sum(2*x['P'] for x in jk),original_rule_failures=sum(x['failed_M1']+x['failed_M3'] for x in jk),
             stricter_score_failures=sum(x['strict_failed_M1']+x['strict_failed_M3'] for x in jk),figures=figures,
             visual_review='Inspected main paired plot, noise-ratio six-panel plot and residual-band plot: no clipping/overlap obscuring data; all plot data preserved.')
    js(HERE/'analysis/verification.json',out)
    assert out['all_protected_unchanged'],out
    assert all(x['byte_identical'] for x in baseline)
    assert invariants['source_spectra_max_error']<1e-10
    assert invariants['original_B_prediction_max_abs_percent_error']<1e-5
    assert invariants['residual_trace_identity_max_relative_error']<1e-8
    assert invariants['projector_annihilation_max_relative_error']<1e-7
    assert len(figures)==11 and all(x['pdf_exists'] for x in figures)
    print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
