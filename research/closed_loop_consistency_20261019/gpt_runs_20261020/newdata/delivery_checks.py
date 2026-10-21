"""Checks artifact integrity, not whether the scientific hypothesis passed."""
import ast,csv,hashlib,json,platform,sys
from pathlib import Path
import numpy as np
import scipy,matplotlib,PIL,fitz
from PIL import Image

HERE=Path(__file__).resolve().parent

def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    inv=read(HERE/'data_inventory.json');requests=read(HERE/'fit_requests.json')
    expected={(r['id'],m) for r in requests for m in r['methods']}
    fs=[read(p) for p in (HERE/'fits').glob('*.json')]
    assert {(r['id'],r['method']) for r in fs}==expected and len(fs)==94
    assert len(inv)==27 and sum(bool(r.get('duplicate_of')) for r in inv)==1
    uncertainty=[]
    for r in fs:
        assert len(r['beta'])==158 and np.isfinite(r['beta']).all()
        assert len(r['frequency_hz'])==3 and len(r['damping_percent'])==3
        a=np.load(HERE/'fits'/f"{r['id']}_{r['method']}.npz")
        assert a['G'].shape==(133,4,4) and np.isfinite(a['G']).all()
        u=read(HERE/'analysis'/f"{r['id']}_{r['method']}_uncertainty.json")
        assert np.isfinite(u['frequency_u_hz']).all() and np.isfinite(u['damping_u_percentage_points']).all()
        uncertainty.append(u)
    pars=list(csv.DictReader((HERE/'analysis/all_parameter_estimates.csv').open(encoding='utf-8-sig')))
    assert len(pars)==564
    draws=list((HERE/'analysis/sml_bootstrap').glob('*.json'))
    assert len(draws)==594
    for p in HERE.glob('*.py'):ast.parse(p.read_text(encoding='utf-8'),filename=p.name)
    figures=[]
    for item in read(HERE/'figures/manifest.json'):
        png=HERE/'figures'/f"{item['name']}.png";pdf=png.with_suffix('.pdf')
        with Image.open(png) as im:
            dpi=im.info.get('dpi',[0,0]);assert min(dpi)>299.9
            dimensions=im.size;im.verify()
        with fitz.open(pdf) as doc:
            assert len(doc)==1 and len(doc[0].get_text().strip())>50
            text_chars=len(doc[0].get_text())
        figures.append(dict(name=item['name'],size_px=dimensions,dpi=dpi,pdf_text_characters=text_chars))
    assert len(figures)==20
    # Independent unwrap and a global branch translation preserve the complex
    # phase. This does not test away local unwrap ambiguity in noisy data.
    phase_error=0.
    for ctrl in ('lowgain','nominal','highgain'):
        g=np.load(HERE/'analysis'/f'{ctrl}_orthogonal_frf.npz')['G']
        phase=np.unwrap(np.angle(g),axis=0)
        phase_error=max(phase_error,float(np.max(abs(np.exp(1j*phase)-g/abs(g)))))
    assert phase_error<1e-12
    changed=[]
    for r in inv+read(HERE/'baseline/sources.json'):
        if sha(Path(r['source']))!=r['sha256']:changed.append(r['source'])
    assert not changed
    required=('README.md','CHECKLIST.md','QA_FIGURES.md','新实测数据分析报告.md')
    assert all((HERE/p).is_file() for p in required)
    result=dict(artifact_checks_passed=True,scientific_hypotheses_all_passed=False,
        environment=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,
                         matplotlib=matplotlib.__version__,Pillow=PIL.__version__,PyMuPDF=fitz.VersionBind),
        fits=len(fs),parameter_rows=len(pars),finite_uncertainty_files=len(uncertainty),SML_bootstrap_records=len(draws),
        solver_failures=[dict(id=r['id'],method=r['method']) for r in fs if r['failure']],
        strict_score_failures=sum(not r['strict_score_pass'] for r in fs),figures=figures,
        source_files_unchanged=True,unwrap_complex_invariance_max_error=phase_error,
        syntax_checked=[p.name for p in HERE.glob('*.py')])
    (HERE/'analysis/delivery_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    checksums={str(p.relative_to(HERE)):sha(p) for p in list(HERE.glob('*.py'))+list(HERE.glob('*.md'))+list((HERE/'figures').glob('*.pdf'))}
    (HERE/'analysis/delivery_checksums.json').write_text(json.dumps(checksums,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('figures','syntax_checked','environment')},ensure_ascii=False))

if __name__=='__main__':main()
