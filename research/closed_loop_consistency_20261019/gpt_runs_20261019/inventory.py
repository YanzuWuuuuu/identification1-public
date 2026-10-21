"""Immutable source snapshot and hash verification; only writes within this run."""
from pathlib import Path
import json,hashlib,shutil,sys,platform
import numpy,scipy
RUN=Path(__file__).resolve().parent; BASE=RUN.parent; WORK=BASE.parents[1]
def files():
    return sorted(list((BASE/'code').glob('*.*'))+list((BASE/'real_data').glob('*.py'))+list(BASE.glob('*.md'))+list((WORK/'experiment_data/srivc_no_lrm_20261015/closed_loop_consistency_20261019').glob('cl_data_z*.mat')))
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()
if __name__=='__main__':
    dest=RUN/'source_inventory.json'
    if '--verify' in sys.argv:
        old=json.loads(dest.read_text(encoding='utf-8'))
        changed=[x['path'] for x in old['files'] if not Path(x['path']).exists() or digest(Path(x['path']))!=x['sha256']]
        (RUN/'source_integrity_check.json').write_text(json.dumps(dict(unchanged=not changed,changed=changed,checked=len(old['files'])),indent=2),encoding='utf-8')
        print('unchanged',not changed,'checked',len(old['files']))
    else:
        entries=[]
        for p in files():
            entries.append(dict(path=str(p),bytes=p.stat().st_size,sha256=digest(p)))
            if p.suffix in ('.py','.md','.txt'):
                target=RUN/'baseline_snapshot'/p.relative_to(BASE);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
        dest.write_text(json.dumps(dict(python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,platform=platform.platform(),files=entries),indent=2,ensure_ascii=False),encoding='utf-8')
        print('saved',len(entries),'source hashes')
