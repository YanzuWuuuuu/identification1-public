# v4 证明修订文件说明

主文件在研究目录：`闭环一致性_频域加性SRIVC_数学推导与验证.md`。本目录保存修订备份、依据、独立核查和重建脚本；不改实验实现和数据。

- `原稿_v3_修订前.md`：修改前的完整版本。
- `theory_head_v4.md`、`theory_tail_v4.md`：修订后的第 0–8.5 节和第 11–13 节。
- `assemble_revision.py`：从备份保留第 9–10 节表格，加入明确的范围说明，原位重建主文件。会验证历史表格完全一致和公式定界符成对。
- `修订说明与核查结果.md`：读者优先阅读的修改摘要和限制。
- `random_weights_proof.md`：两条渐近路线、固定邻点权和逆矩包络的完整笔记。
- `deterministic_proof_audit.md`、`literature_iv.md`、`literature_draft_audit.md`：独立数学与原文核查。
- `check_cross_schur.py`、`cross_schur_checks.json`：总体矩导数和 Schur 协方差代数检查。
- `source_inventory.json`、`pdf_text/`：本地 PDF 的哈希与带页码的提取文本。p01–p18 对应索引中的原始文件；渲染页用于核对排版公式。
- `revision_verification.json`：主文件、备份及表格保留检查。
- `verify_original_files.py`、`original_files_recheck.json`：相对前次实验清单的只读哈希核查，只允许本次授权的主理论文档变化。

在本目录打开 PowerShell，执行：

```powershell
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe' '.\check_cross_schur.py'
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe' '.\assemble_revision.py'
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe' '.\verify_original_files.py'
```

这里只重建文档及核查公式，不启动任何新的 MC。`assemble_revision.py` 会覆盖主理论文档；人工继续编辑主文后，不应无意重跑覆盖那些后续编辑。
