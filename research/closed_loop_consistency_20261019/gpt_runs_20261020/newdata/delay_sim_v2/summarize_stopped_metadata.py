"""Bookkeeping only: summarize existing saved runtimes and interrupted jobs."""
from pathlib import Path
import json, statistics

root=Path(__file__).resolve().parent
raw=root.parent/'delay_sim_closeout/mc_qr_final/raw'
groups={}
for p in raw.glob('*.json'):
    if p.name.endswith('.started.json'): continue
    q=json.loads(p.read_text(encoding='utf8'))
    groups.setdefault(f"{q['band']}_{q['method']}",[]).append(q.get('seconds',0))
interrupted=[str(p) for p in raw.glob('*.started.json') if not p.with_name(p.name.replace('.started.json','.json')).exists()]
out=dict(completed_runtimes={k:dict(n=len(v),median_seconds=statistics.median(v),mean_seconds=statistics.mean(v),max_seconds=max(v)) for k,v in groups.items()},
         interrupted_checkpoints=interrupted,
         note='Saved v1 times only; not a v2 benchmark. No model/data generation or fitting.')
(root/'stopped_metadata.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(out,ensure_ascii=False,indent=2))
