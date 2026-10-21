# P5 既有结果的补充核查

只读复用 `newdata/reparam_tests` 中六条 A/C++++ 记录、H300/H493、M3、模型 a/b 的全样本与配对删周期结果。不新拟合，不分析公共时延。

运行：

```powershell
& 'C:\Users\Administrator\AppData\Local\Programs\Python\Python311\python.exe' -B analyze_p5.py
```

脚本独立从 raw 重建跨窗口差及其配对 jackknife 协方差，并与既有两个 CSV 核对；使用父目录已冻结的 `input_hashes.json` 核对读取过的文件。输出只写本目录。

理论 P5 没有指定覆盖因子，所以报告并列 1u/2u 描述性比较，不称其为新的预设阈值或 95% 置信区间。运行间复现性保留两个窗口的相关性，按三个控制器内部的 A/C 重复差估计；原始复现性加 JK 的合成结果与扣除已估统计噪声的敏感性结果都保留。八个相对分量有两个和为零约束，因此 48 个比较不是 48 个独立检验。

主报告：`P5_相对时延跨窗口核查.md`。机器可读结论：`p5_summary.json`。完整原始配对折、跨分量协方差、逐分量明细和审计见本目录 JSON/CSV。
