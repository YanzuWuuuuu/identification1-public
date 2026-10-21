"""Apply the authorized proof repair; preserve all old experiment tables.

Reads the archived v3, never old result/code files. Re-running is deterministic.
"""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = ROOT / "闭环一致性_频域加性SRIVC_数学推导与验证.md"
old = (HERE / "原稿_v3_修订前.md").read_text(encoding="utf-8-sig")
history = old[old.index("## 9. 数值验证"):old.index("## 11. 新颖性边界")]
history = history.replace(
    "## 9. 数值验证\n",
    "## 9. 数值验证（保留历史实验）\n\n"
    "> v4 范围注：本节数字和表格原样保留，未重新计算。它们是特定有限样本结果；"
    "模型正确性、局部根和新定理条件须分别核查。旧实验中的初始模型权与本轮 M3w "
    "邻点谱权不是同一种权重。第 9.4 节仅保留原有解耦诊断，不新增解耦主张。\n",
)
replacements = {
    '只有输入侧在长数据下有可分辨的差别': '在这些历史仿真中，只有输入侧在长数据下有可分辨的差别',
    '闭环偏差主要进入阻尼和频率，即时间特性；空间方向的误差主要由权重决定。所以**闭环一致性修正本身不会改善零空间方向，也就不能直接解释或解决此前的解耦问题。**':
    '在这些仿真中，闭环偏差主要表现于阻尼和频率，空间方向误差随所用权重变化；这些结果**没有提供足以解决此前解耦问题的改善证据**。不能把该有限组实验扩展为所有模型下空间方向必定不变的定理。',
    '直接法·全部频点：用全部 1057 个 DFT 频点，相当于 LRM 或非周期数据的做法；':
    '直接法·全部频点：用全部 1057 个 DFT 频点构造本文的直接 IV 矩；这不等同于所有 LRM 或其他非参数流程；',
    '4. **实践建议。** 闭环下估计低频支承模态时，不要使用无激励频点，即不要对 \\(u,y\\) 直接做跨频点平滑的 LRM 再拟合。应只用周期谱线，或者用参考 IV；参考 IV 不依赖激励的周期性。':
    '4. **实践建议（v4 修正）。** 本组数据支持避免把非激励频点的反馈噪声当作对象信息，直接套入未经闭环校正的输入梯度矩。可比较周期谱线直接法、外部参考 IV 或适当的闭环 EIV/非参数方法；不能据此笼统禁止 LRM。参考 IV 是否适用于非周期激励，仍需该设计下的外生性、信息与随机极限条件。',
    '**推论 3 得到实测支持。**': '**这些实测差异与推论 3 的机制相容。**',
    '即控制器被精确恢复。': '即记录信号之间的控制器拟合残差很小；这不单独证明所有未记录环节均已准确建模。',
}
for a, b in replacements.items():
    if a not in history:
        raise RuntimeError(f"Expected historical qualification text missing: {a[:70]}")
    history = history.replace(a, b)

head = (HERE / "theory_head_v4.md").read_text(encoding="utf-8")
tail = (HERE / "theory_tail_v4.md").read_text(encoding="utf-8")
new = head.rstrip() + "\n\n" + history.strip() + "\n\n" + tail.strip() + "\n"
table_rows = lambda s: [x.rstrip() for x in s.splitlines() if x.startswith("|")]
old_history = old[old.index("## 9. 数值验证"):old.index("## 11. 新颖性边界")]
assert table_rows(old_history) == table_rows(history), "Historical numeric tables changed"
assert new.count("\\[") == new.count("\\]")
assert new.count("\\(") == new.count("\\)")
assert "（v4）" in new.splitlines()[0]
TARGET.write_text(new, encoding="utf-8")
audit = {
    "target": str(TARGET),
    "source_backup": str(HERE / "原稿_v3_修订前.md"),
    "backup_sha256": hashlib.sha256((HERE / "原稿_v3_修订前.md").read_bytes()).hexdigest(),
    "revised_sha256": hashlib.sha256(TARGET.read_bytes()).hexdigest(),
    "historical_table_rows_preserved": len(table_rows(history)),
    "history_tables_identical": True,
    "display_math_pairs": new.count("\\["),
    "inline_math_pairs": new.count("\\("),
    "model_data_and_experiment_code_changed": False,
    "scope": "Authorized proof document repair, no new Monte Carlo runs",
}
(HERE / "revision_verification.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(audit, ensure_ascii=False, indent=2))
