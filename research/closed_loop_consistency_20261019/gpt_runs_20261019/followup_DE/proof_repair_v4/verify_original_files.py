from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
inventory = json.loads((HERE.parent / "original_inventory.json").read_text(encoding="utf-8"))
allowed = HERE.parents[2] / "闭环一致性_频域加性SRIVC_数学推导与验证.md"
changed, missing = [], []
for name, expected in inventory.items():
    path = Path(name)
    if not path.is_file():
        missing.append(name)
        continue
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for part in iter(lambda: handle.read(2**20), b""):
            digest.update(part)
    if digest.hexdigest() != expected:
        changed.append({"path": name, "authorized_proof_revision": path.resolve() == allowed.resolve()})
result = {"checked": len(inventory), "missing": missing, "changed": changed,
          "passed": not missing and all(row["authorized_proof_revision"] for row in changed)}
(HERE / "original_files_recheck.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(result, ensure_ascii=True, indent=2))
if not result["passed"]:
    raise SystemExit(1)
