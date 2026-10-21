"""Assemble the final report from completed, immutable per-fit checkpoints."""
from pathlib import Path
import json,re,hashlib
import numpy as np
R=Path(__file__).resolve().parent
def read(p):return (R/p).read_text(encoding='utf-8')
def load(p):return json.loads(read(p))
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def validate():
    a=[json.loads(p.read_text(encoding='utf-8')) for p in (R/'experiment_a/results').glob('*.json')]
    assert len(a)==720
    for p in [1,2,4,8,16]:
        for m in ['dir_band_W1','dir_lines_W1','ref_band_W1','refD_band_W1']:
            assert sum(x['P']==p and x['method']==m for x in a)==(20 if p==16 else 40)
    a1={(x['replicate'],x['method']):x for x in a if x['P']==1}
    assert all(np.array_equal(a1[r,'dir_band_W1']['beta'],a1[r,'dir_lines_W1']['beta']) for r in range(40))
    cs=load('comparators/summary.json')
    assert sum(x['n'] for x in cs if x['method'] in ('M5','M5_reg'))==140
    assert sum(x['n'] for x in cs if x['method'].endswith('_IVGN'))==160
    assert sum(x['n'] for x in cs if x['method'] in ('M6u','M6_reg'))==160
    assert len(list((R/'experiment_c/raw').glob('*.json')))==120
    b=load('experiment_b/summary.json')
    assert b['n_mc_requested']==20 and len(b['methods'])==3
    sources=load('source_inventory.json')['files']
    changed=[p['path'] for p in sources if digest(p['path'])!=p['sha256']]
    assert not changed,changed
    out=dict(A_fits=720,A_replicates=180,B_final_fits=60,C_cases_including_unavailable=120,
             SML_main_fits=140,indirect_IVGN_main_fits=160,indirect_LM_optimizer_checks=160,
             original_files_checked=len(sources),original_files_changed=changed,P1_M1_equals_M2=True,
             complete_requested_replicate_counts=True)
    (R/'package_validation.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    return out

def make_report():
    val=validate()
    head='''# 频域闭环 SRIVC／参考 IV：三组对照实验报告

日期：2026-10-19。所有新增代码、结果、图和记录均在本目录；原理论、原代码、原始数据和已有结果未修改。以下是实验结果，不重写理论，也不将预期当成结论。

## 1. 结论先行

**参考 IV 能明显减少本平台低频支承模态的闭环偏差，但本次结果不支持“单周期只有参考 IV 可用”，也不支持“参考 IV 全面优于经典周期 SML”。**

- 8 周期、40 次配对仿真中，第一支承模态阻尼的均值相对误差：全部频点直接法 M1 为 **+45.67%±2.34%**，谱线直接法 M2 为 **+2.80%±1.68%**，正确控制器参考 IV M3 为 **−1.59%±1.63%**。M3 有明确的减偏价值。
- 单周期时 M3 有 **6/40** 次、无控制器参考 M4 有 **19/40** 次达到迭代上限。M3 第一模态平均误差为 −22.65%，即使只看收敛子集仍约 −12.39%。这不是归一化或数据配对错误；两个代表样本增加至 240 步也没有修复。渐近一致性不保证单周期的求根可靠性。
- 周期 SML 在 P=8/16 的主运行均达到求解及 score 门槛，阻尼 RMSE 小于本次均匀权 M3。例如 P=8 第一模态 **7.72% 对 10.28%**，P=16 为 **3.04% 对 9.91%**。SML 利用了噪声协方差，M3 此处为均匀权，因此不能把差异只归于 IV 与 ML 两种框架本身。
- 带五个柔性模态时，M3 的柔性参数均值与零误差相容；M1 在部分柔性参数上仍有小而可分辨的偏差。“大偏差主要集中在控制带宽内”有支持，“带宽外严格零偏差”不成立。
- 实测 P=8 时，按要求替换 ZOH 模态块并重分配 100 µs 延迟后，新旧估计最大差仅 **0.000718 Hz、0.005913 个阻尼百分点**。ZOH 修正不是此前方法差异的主要来源。

## 2. 运行范围、统计口径和数值核查

|组别|完成范围|原始结果位置|
|---|---|---|
|A，M1–M4|P=1/2/4/8 各40次，P=16为20次，共720次拟合|experiment_a/results；data|
|A，SML|140次主拟合；P=2/4为明确标注的收缩诊断，P=8/16为原始SML|comparators/sml_results|
|A′，间接FRF|100组四次单轴实验；160次规范IV–GN拟合；另160次信赖域优化器对照|comparators/indirect_data；indirect_ivgn_results；indirect_results|
|B|20组配对数据×3方法，共60个最终结果|experiment_b/raw|
|C|96个新旧模型主结果、4个第二半段结果、20个SML/不可用结果记录|experiment_c/raw|

所有仿真参数误差使用 `100*expm1(log估计−log真值)`。SD 为重复间样本标准差，SE=SD/√n，RMSE 为所有有限端点相对真值的均方根误差。主表不删除失败端点或异常值；数值异常、迭代上限、协方差不可用和完整 score 未通过分别记录。A/B 另提供仅收敛子集统计，不能把条件子集当作无偏证明。实测没有真值，只报告估计与方法差异。多参数的均值/SE是描述性诊断，没有进行事后挑选的显著性宣称。

原环境自检按指定原脚本完成：3次阻尼平均误差 M1=51.43/17.54/5.35%，M3=2.28/−3.96/1.04%，数量级一致。Python3.11.9、NumPy2.4.6、SciPy1.17.1。采样/噪声/初始化沿用任务；对不同方法使用同一记录及相同极点初始扰动。

主要正确性检查：P=1 的 M1/M2 参数逐项相同；新 A 与原脚本固定样本的直接法及M4完全相同，M3模态对数参数差约1e−6；B稀疏重排的正规矩阵误差5.78e−16；ZOH无噪声DFT恒等式残差3.83e−12；SML无噪声均值恢复模型误差3.53e−12；间接FRF无噪声模型误差7.69e−12。独立人工模型核对SML白化雅可比与中心差分误差约1e−9，FRF协方差向量化算子误差1.29e−16。

## 3. 经典对照的必要实现修正

### 3.1 周期 SML 的目标和导数

周期确定性潜在输入的剖面目标采用任务中的 `Σ eᴴ Ce(θ)⁻¹ e`，均值协方差除以 P(P−1)，保留输入输出交叉协方差。**Ce 依赖参数，故优化时对白化残差整体求导，包含 Cholesky 因子的导数；不把冻结Ce的IRLS误称为完整SML。** 此周期剖面目标不额外加入 logdet(Ce)。2007年一般随机输入 EIV-ML 的目标与此不同，不能混用。公式核对见 [Hong等 IFAC 2008，§4式(10)–(12)](https://skoge.folk.ntnu.no/prost/proceedings/ifac2008/data/papers/0643.pdf)；同时核对了用户本地的2007 EIV与BJ Part I全文。来源及实现细节在 comparators/README.md。

### 3.2 少周期协方差的秩限制

这些限制针对本任务指定的“逐频点、只用周期间样本协方差、没有频率平滑/已知噪声结构”的实现，不是所有EIV或FRF方法都不能用少周期。

|原始方法|P|协方差维数|秩上界|处理|
|---|---|---:|---:|---|
|M5 SML|1|4×4|无样本协方差|不可用|
|M5 SML|2/4|4×4|1/3|原始逆不存在；另跑M5_reg|
|M5 SML|8/16|4×4|最多4|跑原始M5，不加收缩|
|M6间接FRF|4/8/16|16×16|P−4=0/4/12|原始完整逆不存在；M6u可用，P8/16另跑M6_reg|

M5_reg 将联合[Y,U]样本协方差向其对角收缩10%；M6_reg 将完整FRF协方差向其对角收缩10%。它们是明确标注的不同估计器，不冒充原始SML或最优权。另在固定样本上做1%/30%收缩敏感性；完整结果保留。仅按秩条件，四输出SML至少需5周期才可能满秩；四次单轴FRF的完整协方差至少需每次5周期、总20周期才可能满秩，满秩仍不保证估计稳定。

### 3.3 间接 FRF 的拟合与预算

每次单轴实验先丢弃2周期，以四次记录的均值组成 T_ud=UD⁻¹、T_yd=YD⁻¹，再取 G=T_yd T_ud⁻¹。协方差按 δG=(δY−GδU)U⁻¹ 一阶传播，包括通道间相关；该有限样本比值不是严格无偏的线性观测。

主拟合保留原 cl2 的梯度工具变量 IV–GN 正规方程、对角缩放、极点步长上限及80次停止预算，移除全部瞬态自由参数。原函数的free开关会关闭非线性更新、且只支持对角输出权，因此需要新建适配器；没有直接误用该开关。额外LM/TRF求同一目标的结果单列，不能与主IV–GN混合。

四次实验总保留样本数等于单次四路同时激励，但预热成本为8周期而不是2周期；同一单路幅值下，总注入能量约为同时四路激励的1/4。因此本对照是等保留数据量，**不是同时等总采集时间与等激励能量**。不把间接法本轮表现差解释成该类方法普遍不适合低频。

## 4. 实验 A：周期数扫描

M1=全频点直接法；M2=仅激励谱线直接法；M3=正确K参考IV；M4=直接用D构造工具输入。下表均为主运行全部有限端点。
'''
    a=read('experiment_a/report_a.md')
    at=a[a.index('|P|方法|阻尼'):a.index('## 频率相对误差')]
    cp=read('comparators/results_tables.md')
    cpt=cp[cp.index('| 方法 | P'):cp.index('完整六参数统计')]
    paired=read('experiment_a/paired_analysis.md').split('\n',2)[2]
    between='''
完整频率误差、收敛子集、逐次异常清单见 experiment_a/report_a.md 与 summary.csv/json。每次完整参数和分数均保留。

### A 的少周期异常核查

没有发现DFT归一化、频带索引、方法数据配对或初值实现错误。P1的M3/M4不仅有失败端点，收敛子集第一模态阻尼均值仍约−12%；两个固定样本给所有方法同等增加到240步，部分IV阻尼继续趋零，不能靠增加迭代次数解决。当前证据指向有限样本非线性求根、近邻模态与敏感方向共同影响；尚未通过因素消融唯一归因于“弱工具变量”。详见 experiment_c/diagnostics/audit_a_note.md。

### A 的配对趋势与效率

'''
    comps='''
## 5. SML 与间接 FRF 的主对照及优化器检查

下表 `M6*_IVGN` 是任务要求的规范 IV–GN 主结果；没有后缀的 `M6*` 是同一目标的信赖域优化器检查。`M5_reg` 与 `M6_reg` 是收缩版本，原始不可用条件见第3节。失败数与完整score未过数分别列出（可能重叠，不相加）。

'''
    cpa=cp[cp.index('## 间接逆矩阵的条件数'):]
    endcomp='''
SML在P8/16的RMSE优于当前均匀权参考IV，是本次必须保留的正面对照。P16第一模态均值1.52±0.60%仍约2.5SE，不能把“有限样本无偏”当作已证明结论。M5_reg少周期运行也不差，但它依赖额外收缩假设，不能据此代表未正则化经典SML。

间接FRF的M6u在少周期出现较大散布和求解失败；信赖域优化减少部分失败，但没有消除均值偏移和大散布。M6_reg收敛较好，却仍有明显偏差，例如P8第一模态+31.68±4.79%。该结果同时受有限周期比值、plug-in协方差、收缩和较少激励能量影响，不能单独归因于间接辨识原理。条件数与异常有相关性，但本轮没有做因果隔离。

## 6. 实验 B：完整支承与柔性模型

'''
    b=read('experiment_b/report_b.md').split('\n',2)[2]
    c=read('experiment_c/experiment_c_report.md').split('\n',2)[2]
    conclusion='''
## 8. 逐条对照已有结论

|已有结论/预期|本轮判定|依据和限制|
|---|---|---|
|M1在带宽内存在不随P消失的大偏差|支持其大样本趋势；小P并非常数|P4/8/16第一阻尼均值45.20/45.67/55.10%；P1受有限样本影响为−0.17%，不能要求全部P都约45%|
|M2只用谱线，偏差按1/P下降|与趋势相容，未单凭有限点证明阶数|M2−M3第一阻尼配对差，P×均值在P2/4/8/16约31.7/30.1/35.1/31.1；小P求根失败需保留|
|M3无偏|支持较长记录减偏；不支持有限样本无条件无偏|P8/16接近真值；P1有负偏移、近零阻尼和6/40次失败|
|M4不需要正确控制器|有限数据中仍能减偏，但效率及求解更脆弱|P1失败19/40；P8全收敛，第一阻尼−7.18±3.53%；不能把任意工具的外生性等同于足够相关性|
|M4方差约M3两倍|不支持这一固定倍率|P8标准差比约2.1–2.3，方差比约4.2–5.2；“标准差”与“方差”应区分|
|P1/2只有M3/M4可用|本实现不支持|M1/M2也可求解；M3/M4存在失败和有限样本误差。原始SML不可用只说明协方差不足，不说明IV必定可靠|
|M5在P≥4无偏|需要修正|四输出P4仍秩不足；P8/16原始SML表现较好但不是严格有限样本无偏；P2/4只能报告收缩诊断|
|M6无偏但方差较大|本次小周期实验不支持“无偏”|原始完整协方差逆不存在；M6u/M6_reg存在均值偏移、散布或失败；等保留N未匹配总能量|
|柔性模态三法均无偏，偏差只在带宽内|强表述不支持|M1部分柔性偏差虽小仍可分辨；M3柔性结果与零误差相容|
|实测ZOH精确求值应核对|已完成，P8影响很小|最大0.000718Hz/0.005913阻尼百分点；只模态块精确，背景和有效延迟重分配假设仍保留|
|实测M1−M3差不随P消失，M2更靠近M3|部分支持，不能确认普遍阶数|P8的z0/.5第二模态有明显差异；少周期失败，z1排序有例外；实测没有真值|

## 9. SML 与间接 FRF 对比参考 IV 的优劣与适用场景

|方面|参考IV（M3/M4）|周期SML（M5）|间接FRF再拟合（M6）|
|---|---|---|---|
|数据入口|单次实验的参考、总输入、输出；直接DFT|重复周期的输入输出均值与联合协方差|四次独立单轴实验或其他足够秩的独立激励设计|
|周期数|方程可在1周期建立；本次1/2周期可靠性不足|本任务四输出逐频点协方差至少5周期才可能满秩；8/16周期有效|本任务每轴至少5周期、总20周期才可能估计可逆完整16维协方差；均匀/结构权另论|
|是否需要控制器|M3使用；M4不用但方差/失败率更高|这里的周期EIV版本不需要控制器；不能与依赖正确K的BJ-ML混称|已记录外激励与总输入输出时，FRF间接构造不需要控制器|
|是否需要噪声模型|均匀权不用参数噪声模型，也不用样本协方差|不用参数化噪声传函，但需要可靠的非参数联合协方差|均匀权不用；协方差权需要足够数据/噪声统计估计|
|本次优势|P8/16明显减小直接法大偏差，单次实验、不做FRF求逆|P8/16阻尼RMSE优于当前均匀权IV|可得到完整可检查的FRF，用于模态结构检查和多种后续参数模型|
|本次代价|少周期可能根不稳、阻尼塌缩；不用K效率下降|协方差秩与有限样本误差、非线性求解；实测P8本次未通过严格score|多次实验与更多预热；弱闭环输入方向求逆放大噪声；本次少周期原始协方差不可逆|

**本轮可保留的工程结论**：当有外生参考、希望单次采集并避免FRF矩阵求逆时，参考IV是有用入口，尤其能修正全频点直接法在低频支承模态的大偏差；但它的价值不是“单周期保证成功”或“必然优于SML”。已有足够重复周期且协方差可可靠估计时，周期SML是很强的基线，本次仿真甚至更有效率。间接FRF适合需要完整非参数矩阵、能安排足够独立实验的场景；本次有限预算和正则化版本的差结果不能否定成熟的多周期、多实现流程。

## 10. 交付与复核

- `README.md`、`run_all.ps1`：环境、命令、顺序、恢复运行方式。
- `source_inventory.json`、`source_integrity_check.json`：26个原始代码/文档/数据文件SHA256与复核；C另核对固定时延文件。
- `package_validation.json`：请求的主重复次数、文件数、P1恒等检查和原文件完整性。
- `experiment_A_period_scan.png/pdf/svg`、`experiment_A_failure_rates.png`：仿真均值、RMSE及失败率。
- `experiment_c/experiment_c_damping_scan.png`：实测扫描；失败点与描述性两半误差棒已区分。
- 各子目录的原始NPZ/JSON、完整参数、种子、初始化、停止状态与耗时是统计来源。报告未替换原有论文或理论文档。

原始不可计算方法没有虚构数值。少周期失败、初期NumPy接口错误及同种子修复重算、严格门槛补查未通过、SML实测不稳定、间接法较差表现均保留。新的图表和文字仅据本次实际运行结果形成。
'''
    text=head+at+between+paired+comps+cpt+cpa+endcomp+b+'\n\n## 7. 实验 C：实测 ZOH 与周期数扫描\n\n'+c+conclusion
    (R/'三组对照实验_完整报告.md').write_text(text,encoding='utf-8')
    print('Report saved; validation:',val)

if __name__=='__main__':make_report()
