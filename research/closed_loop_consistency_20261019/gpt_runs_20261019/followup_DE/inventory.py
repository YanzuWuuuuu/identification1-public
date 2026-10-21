"""Hash original source files and supplied data without modifying them."""
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(4*1024*1024),b''):h.update(x)
    return h.hexdigest()
def main():
    target=HERE/'original_inventory.json'
    if '--verify' in sys.argv:
        old=json.loads(target.read_text(encoding='utf-8'))
        changed=[q for q,v in old.items() if not Path(q).exists() or digest(Path(q))!=v]
        out={'count':len(old),'unchanged':not changed,'changed':changed}
        (HERE/'original_integrity_check.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(out);return
    paths=[]
    for top in [BASE/'code',BASE/'real_data',HERE.parent]:
        paths.extend(p for p in top.rglob('*') if p.is_file() and HERE not in p.parents and '__pycache__' not in p.parts and p.suffix.lower() in {'.py','.md','.json','.ps1'})
    paths.extend(BASE.glob('*.md'))
    workspace=BASE.parents[1]
    paths.extend((workspace/'experiment_data'/'srivc_no_lrm_20261015'/'closed_loop_consistency_20261019').glob('cl_data_z*.mat'))
    paths.extend((HERE.parent/'experiment_a'/'data').glob('*.npz'))
    out={str(p):digest(p) for p in sorted(set(paths))}
    target.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8');print({'files':len(out)})
if __name__=='__main__':main()
