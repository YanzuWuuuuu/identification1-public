"""Copy saved v1 endpoints/logs without modifying or deleting originals."""
from pathlib import Path
import hashlib, json, shutil
from datetime import datetime

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'delay_sim_closeout'
DEST = ROOT / 'archive_stopped_20261021'
assert not DEST.exists(), 'Refuse to overwrite an existing archive'
DEST.mkdir()
manifest = []
for path in sorted(OLD.rglob('*')):
    if not path.is_file():
        continue
    rel = path.relative_to(OLD)
    # Spectral data remain intact in place; endpoints, logs and audit files copied.
    if rel.parts[:2] == ('mc', 'data'):
        continue
    if '__pycache__' in rel.parts:
        continue
    content = path.read_bytes()
    item = dict(source=str(path), relative=str(rel), size=len(content),
                sha256=hashlib.sha256(content).hexdigest())
    target = DEST / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, target)
    assert hashlib.sha256(target.read_bytes()).hexdigest() == item['sha256']
    manifest.append(item)
data_files = [dict(path=str(p), size=p.stat().st_size) for p in (OLD/'mc/data').glob('*') if p.is_file()]
status = json.loads((OLD/'RUN_STATUS.json').read_text(encoding='utf-8-sig'))
result = dict(archived_at=datetime.now().isoformat(), source=str(OLD),
              copied_files=manifest, data_retained_in_place=data_files,
              stopped_status=status, note='No old endpoint or log deleted. Formal MC is stopped by user request.')
(DEST/'ARCHIVE_MANIFEST.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(dict(copied_files=len(manifest), copied_bytes=sum(x['size'] for x in manifest),
                     data_retained_in_place=len(data_files), status=status),ensure_ascii=False))
