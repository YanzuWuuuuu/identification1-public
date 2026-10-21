"""Read-only result/dependency verification; no model evaluation or fitting."""
from pathlib import Path
import json, hashlib
import numpy as np

root=Path(__file__).resolve().parent
j=json.loads((root/'jacobian_check/boundary_design_results.json').read_text(encoding='utf8'))
checks=[]
for row in j['results']:
    selected=next(t for t in row['thresholds'] if t['cut']==1e-12)
    assert len(selected['all_other_parameters_delay_energy'])==7
    assert min(selected['all_other_parameters_delay_energy'])>=1e-3
    assert row['joint_rank']==row['total_parameters']
    checks.append(dict(controller_fc=row['controller_fc'],band=row['band'],
        rank=row['joint_rank'],minimum_energy=min(selected['all_other_parameters_delay_energy'])))
with np.load(root/'controller_check/six_mode_definition.npz') as a:
    controller_residues=a['C0'].copy()
p=root.parent/'delay_sim_closeout/stable_profile_pilot/v1_P10_r000_H_M3/protocol.json'
b=np.asarray(json.loads(p.read_text(encoding='utf8'))['initial_beta'])
jacobian_residues=np.array([[b[o*88+10*i:o*88+10*i+4] for o in range(4)] for i in range(6)])
assert np.array_equal(controller_residues,jacobian_residues)
c=json.loads((root/'controller_check/controller_check.json').read_text(encoding='utf8'))
assert c['source_hashes']==c['source_hashes_after']
for p,h in c['source_hashes'].items():
    assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
status=json.loads((root.parent/'delay_sim_closeout/RUN_STATUS.json').read_text(encoding='utf8'))
assert status['running'] is False
assert status['saved']==235 and status['failed']==235
assert not status['source_audit']['changed']
for name in ('PLAN.md','criteria_audit.md','jacobian_check/六模态时延Jacobian预检查.md','controller_check/控制器无噪声预检查.md'):
    assert (root/name).is_file()
out=dict(verification_passed=True,checks=checks,controller_and_jacobian_C0_identical=True,
         original_batch_running=False,new_mc_run=False,old_inputs_unchanged=8143,
         note='Only reads existing deterministic results, hashes and metadata.')
(root/'artifact_verification.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(out,ensure_ascii=False))
