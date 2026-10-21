from run import *
import ast,csv,platform,importlib.util
from PIL import Image
import scipy,matplotlib

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    for module in ('run','analyze','plot_results','build_report','qa'):
        assert Path(importlib.util.find_spec(module).origin).resolve()==HERE/(module+'.py'),module
    rows={}
    expected={'modal_parameters.csv':360,'controller_pairs.csv':360,'independent_repeats.csv':180,'FRF_validation.csv':60,'summary.csv':20,'polish_diagnostic.csv':61,'amplitude_diagnostic.csv':20}
    for name,n in expected.items():
        a=list(csv.DictReader((HERE/'analysis'/name).open(encoding='utf-8-sig')));rows[name]=len(a);assert len(a)==n,(name,len(a))
    rr=[json.loads((HERE/'raw'/f'{rid}_{w}_{method}.json').read_text()) for w in SPECS for rid in IDS for method in METHODS]
    assert len(rr)==120 and len(list((HERE/'raw').glob('*.json')))==120
    assert sum('exception' in a for a in rr)==0
    for a in rr:
        assert a['config_hash']==CONFIG_HASH and a['finite']
        assert np.isfinite(a['frequency_hz']).all() and np.isfinite(a['damping_percent']).all()
        p=HERE/'raw'/f"{a['id']}_{a['window']}_{a['method']}.npz";b=np.load(p)
        assert np.array_equal(b['beta'],a['beta']) and np.isfinite(b['G']).all()
    figure_stems=[f'{w}_parameters' for w in SPECS]+[f'{c}_{w}_{k}' for c in ('lowgain','nominal','highgain') for w in SPECS for k in ('magnitude','phase','residual')]+['amplitude_frequency_check']
    for stem in figure_stems:
        png=HERE/'figures'/f'{stem}.png';pdf=HERE/'figures'/f'{stem}.pdf'
        with Image.open(png) as im:im.verify()
        with Image.open(png) as im:assert min(im.info['dpi'])>299
        assert pdf.read_bytes()[:5]==b'%PDF-'
    for p in HERE.glob('*.py'):ast.parse(p.read_text(encoding='utf-8'))
    report=(HERE/'高频实测模态验证报告.md').read_text(encoding='utf-8');assert not any(ord(c)<32 and c not in '\n\t' for c in report)
    assert '\\tau' in report and '\\right' in report and '\t' not in report
    for name in ('weighted_core.py','sml_core_copy.py'):assert digest(HERE/'baseline'/name)==digest(PARENT/'baseline'/name)
    val=json.loads((HERE/'validation.json').read_text());assert val['noiseless_max_log_modal_error']<1e-10 and val['output_vs_old_relative']<1e-12 and val['scaled_moment_matrix_relative']<1e-12
    input_ids=IDS+[f'C_{c}_{s}' for c in ('lowgain','nominal','highgain') for s in ('ppmm','pmpm','pmmp')]+['B_'+c for c in ('lowgain','nominal','highgain')]+['D_nominal_'+v for v in ('0p6','0p8','1p0','1p2')]
    inputs=[PARENT/'prepared'/f'{i}.npz' for i in input_ids]+list((PARENT/'controllers').glob('*.npz'))+[PARENT/'baseline/tau_seconds.npy',PARENT/'prepare_data.py']+list((PARENT/'wideband_scope_audit').glob('*_broadband.npz'))
    for p in inputs:assert p.exists(),p
    manifest=dict(date='2026-10-20',task='High-frequency real-data modal validation',config_hash=CONFIG_HASH,
        runtime=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        primary_fits=120,polished_endpoints=61,secondary_amplitude_fits=8,figures=len(figure_stems),
        input_sha256={str(p):digest(p) for p in inputs},code_sha256={str(p.relative_to(HERE)):digest(p) for p in list(HERE.glob('*.py'))+list((HERE/'baseline').glob('*.py'))},
        scope='four frozen local windows, not a global 3–800 Hz fit; original inputs/results read only',
        status='completed with failures retained; high-frequency invariance not established',
        incident='During final combined import, parent build_report was resolved and regenerated its report and run_manifest. Original measurements/fits/theory were not written. Import precedence repaired and guarded; see implementation_note.md.',
        limits=['high-frequency dense-window solver failure','weighted local models can hit bounds','B-based u omits model/operating variation','no high-frequency E position scan'])
    js(HERE/'manifest.json',manifest)
    js(HERE/'qa.json',dict(completed=True,table_rows=rows,primary_count=120,solver_failures=sum(not a['converged'] for a in rr),near_bounds=sum(a['near_bound'] for a in rr),strict_score_failures=sum(not a['strict_score_pass'] for a in rr),program_exceptions=0,figures_verified=len(figure_stems),baseline_copies_equal=True,numerical_validation=val))
    print('QA_PASS',len(rr),len(figure_stems),rows)

if __name__=='__main__':main()
