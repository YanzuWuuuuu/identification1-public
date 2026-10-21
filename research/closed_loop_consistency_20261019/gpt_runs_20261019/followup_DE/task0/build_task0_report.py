"""Build Task 0 report from retained per-replicate results, including failures."""
import sys
sys.dont_write_bytecode=True
from run_task0 import *

def readmany(pattern):
    return [json.loads(p.read_text(encoding='utf-8')) for p in sorted(pattern.parent.glob(pattern.name))]

def stats(rows):
    def ss(sub):
        e=np.array([r['error_percent'] for r in sub if r.get('error_percent') is not None and np.all(np.isfinite(r['error_percent']))])
        return {'n':len(e),'mean':e.mean(0).tolist(),'se':(e.std(0,ddof=1)/np.sqrt(len(e))).tolist(),'sd':e.std(0,ddof=1).tolist(),'rmse':np.sqrt(np.mean(e*e,axis=0)).tolist()}
    its=[x['iterations'] for x in rows if x.get('iterations') is not None]
    return {'n_total':len(rows),'failures':sum(x['status']!='pole_step_converged' for x in rows),'mean_iterations':float(np.mean(its)),'all_finite':ss(rows),'converged':ss([x for x in rows if x['status']=='pole_step_converged']),'statuses':{s:sum(x['status']==s for x in rows) for s in sorted(set(x['status'] for x in rows))},'max_score':max(x.get('score_max') or 0 for x in rows)}

def main():
    data={
        'P1_short_M1':readmany(OLD/'results'/'P01_r*_dir_band_W1.json'),
        'P1_long_M1':readmany(HERE/'long_results'/'*.json'),
        'P8_M1':readmany(OLD/'results'/'P08_r*_dir_band_W1.json'),
        'P1_short_M3':readmany(OLD/'results'/'P01_r*_ref_band_W1.json')}
    assert all(len(x)==40 for x in data.values())
    summaries={k:stats(v) for k,v in data.items()}
    pops={x:json.loads((HERE/'population'/f'{x}_analytic.json').read_text(encoding='utf-8')) for x in ['P1_short','P1_long','P8']}
    audit=json.loads((HERE/'audit.json').read_text(encoding='utf-8'))
    init=json.loads((HERE/'initialization_check.json').read_text(encoding='utf-8'))
    out={'audit':audit,'summary':summaries,'population':pops,'initialization_diagnostic':init,'conclusion':'No P1-specific implementation mismatch found. Lengthening one period restores a large positive direct-method bias. Short-record nonlinear/nuisance estimation effects mask the population bias; no claim of exact finite-sample agreement or decomposition of those effects.','original_P1_M1_M3_rerun_required':False}
    atomic(HERE/'task0_summary.json',out)
    a=summaries['P1_short_M1']['all_finite']; b=summaries['P1_long_M1']['all_finite']; d=summaries['P8_M1']['all_finite']
    lines=['# 任务 0：单周期直接法偏差异常排查','','结论：**未发现 P=1 特有的激励、噪声、频带、稳态截取或瞬态实现错误。单周期加长八倍后，第一支承模态的直接法阻尼偏差由约 −0.17% 变为 +38.13%，接近谱矩方程预测的 +41.83%。** 原短记录的近零均值不能用于宣称直接法在 P=1 无偏。它体现了短记录估计的有限样本效应与渐近闭环偏差叠加；本次没有将这种叠加进一步唯一分解为瞬态参数估计、非线性小样本偏差和随机谱线方向各自的贡献。','','## 0a. 逐项核对','','| 项目 | P=1 | P≥2 | 核查结论 |','|---|---|---|---|',
    '| 周期长度 | 30,000 点，6 s | 同左 | 一致 |',
    '| 总记录长度 | 30,000 点 | 30,000 P 点 | 唯一预定的长度变化 |',
    '| 每周期激励 | 四路随机相位多正弦，每个正频内部谱线均激励 | 同一周期重复 P 次 | 不是 P=1 漏加噪声或错误留谱线 |',
    '| 频带 | 3–100 Hz，583 个点，全为激励线 | 同频带；激励线仍为 583 个，其他点不激励 | M1/M2 在 P=1 必然相同 |',
    '| 稳态处理 | 先运行并丢掉 2 周期，即 12 s | 同左 | 一致 |',
    '| 噪声 | sd=.2, sv=3e−7, SW=.03 | 同左 | 一致 |',
    '| 控制器/对象 | FC=60 Hz；仅 3 支承模态；5 kHz 采样，精确 ZOH | 同左 | 一致 |',
    '| 瞬态项 | 每模态/输出 (τ₀+τ₁z⁻¹)/(1+a₁z⁻¹+a₂z⁻²)，24 个未知系数 | 同左 | 未写成 z⁻¹、z⁻² |',
    '| 初始化/停止 | 频率 1% 扰动、log 阻尼 .4 扰动、线性初始化；80 次上限、极点步 1e−5 | 同左 | 一致 |',
    '| 失败规则 | 异常、非有限、达到第 80 次 | 同左 | 本报告全部保留 |','','实测缓存一致性检查：全长 DFT 激励线等于 P 倍单周期 DFT 均值，最大绝对差 <8×10⁻¹¹；激励线功率与 N·P·sd²·shape 的相对误差 <7×10⁻¹⁵；P>1 的未激励线 DFT 为数值零。详情见 `audit.json`。','','## 0b. 延长单周期的 40 次蒙特卡洛','','新配置 P=1、P₀=240,000，总有效记录 48 s，与普通 P=8 相同。继续保持 settle=2，因此新配置丢弃 96 s；普通配置丢弃 12 s。两者均先进入稳态。这项长度变化按用户要求保留；不能把它误写成相同的丢弃秒数。种子沿用原 P=1 的 202611190+r（r=0,…,39），随机生成顺序保持原脚本：生成完整时域数据后才抽初始化扰动，因此初始化扰动向量和旧短记录不逐项相同。','','所有误差为 100·expm1(log 估计−log 真值)。表中“±”为均值标准误 SE，SD 为单次估计散布，RMSE 为对真值的均方根误差。下表的三组 M1 均为 40/40 有限且收敛，无删样。','','| 配置 | 参数 | 平均误差 ± SE (%) | SD (%) | RMSE (%) | 失败/总数 | 平均迭代数 |','|---|---|---:|---:|---:|---:|---:|']
    labels={'P1_short_M1':'P=1，P₀=30,000','P1_long_M1':'P=1，P₀=240,000','P8_M1':'P=8，P₀=30,000'}
    names=['f₁','f₂','f₃','ζ₁','ζ₂','ζ₃']
    for name,label in labels.items():
        sm=summaries[name]; ss=sm['all_finite']
        for j,nm in enumerate(names):
            lines.append(f"| {label} | {nm} | {ss['mean'][j]:+.5f} ± {ss['se'][j]:.5f} | {ss['sd'][j]:.5f} | {ss['rmse'][j]:.5f} | {sm['failures']}/{sm['n_total']} | {sm['mean_iterations']:.2f} |")
    lines+=['','## 0c. 总体谱矩方程与仿真对照','','复用 `code/exact_limit.py` 的直接法总体谱矩方程，但在复制版中显式加入周期激励掩码。对未归一化全长 DFT 除以 N 后，激励协方差为：激励线 P·sd²·shape，非激励线 0。输出等效噪声谱为 sv²I+SW²G₀G₀ᴴ；输入–输出等效噪声互谱为 −S K Φᵥ，其中 S=(I+KG₀)⁻¹。原始 `exact_limit.py` 没有调用 `predict2.LINEMASK`，故它在 P>1 时没有逐点表示真实激励谱。本次只修正诊断副本，未修改原文件。','','该方程仍是总体/渐近谱近似：采用结构参数 102 维，不包含 24 个有限窗口瞬态参数的拟合及其边界相关效应，也不表示单条有限随机相位记录的条件均值。对有限 N 的 MC 均值不应要求严格数值相等。','','| 配置 | 总体 ζ₁偏差 (%) | 总体 ζ₂偏差 (%) | 总体 ζ₃偏差 (%) | MC ζ₁±SE (%) | MC ζ₂±SE (%) | MC ζ₃±SE (%) |','|---|---:|---:|---:|---:|---:|---:|']
    for key,smkey,label in [('P1_short','P1_short_M1','短 P=1'),('P1_long','P1_long_M1','长 P=1'),('P8','P8_M1','普通 P=8')]:
        pe=pops[key]['error_percent']; se=summaries[smkey]['all_finite']; lines.append('| '+label+' | '+' | '.join(f'{pe[j]:+.5f}' for j in range(3,6))+' | '+' | '.join(f"{se['mean'][j]:+.5f} ± {se['se'][j]:.5f}" for j in range(3,6))+' |')
    gap=(b['mean'][3]-pops['P1_long']['error_percent'][3])/b['se'][3]
    lines+=['',f'长 P=1 的第一阻尼均值仍比总体预测低 {abs(b["mean"][3]-pops["P1_long"]["error_percent"][3]):.3f} 个百分点（{abs(gap):.2f} 个 MC 均值标准误），因此报告不能写成“完全吻合”或人为取预期的 40–50%。其第二、第三阻尼更接近总体预测。普通 P=8 的均值也不精确等于总体值；有限随机相位方向、短单周期谱线密度和边界参数估计仍然存在。','','数值求根核查：直接沿用原代码的 1e−6 前向差分时，三个总体问题均在 300 次上限前达到约 1e−7–1e−6 的步长噪声平台，无法满足 1e−10 的诊断阈值；这些原始输出仍保留。随后仅在总体诊断中用矩阵指数 Fréchet 导数计算精确 ZOH 参数导数，三组均 7 次迭代达到 1e−10。总体阻尼结果与前向差分的终点相差不足 5×10⁻⁵ 个百分点。仿真 M1 的算法、模型、失败判据没有因此更改。','','## 0d. 原因与是否需修复','','1. **设置/实现：** 未发现原 P=1 配置漏加噪声、错误截频、错误丢周期或瞬态表达错误；无需因这种错误重跑并替换原 P=1 的 M1/M3。M3 的少周期稳健性另按任务 E 处理。',
    '2. **理论：** O(1) 描述渐近数量级，并不要求每个有限周期数下的偏差相同，更不要求 N=30,000 的估计均值等于总体伪真参数。P=1 时没有周期平均，随着单周期长度增长，输入信号与反馈噪声的谱功率比保持同阶，直接法的非零总体偏差仍然存在。此次延长记录显示了这一点。',
    '3. **有限样本：** 原短 P=1 第一阻尼“平均接近零”不是稳健的无偏证据；其 SD 约 28.26%，且长单周期给出显著正偏。短记录的非线性参数估计、自由留数和瞬态 nuisance 参数拟合、有限频率分辨率及有限随机相位方向会共同影响有限样本均值。这里的证据支持有限记录解释，但没有把这些成分唯一分解，不能声称已证明某一个成分是全部原因。',
    '4. **优化初值核查：** 在同一批 40 条短 P=1 数据上，将初始极点改为总体直接法伪真值，重新做线性初始化，并将步长阈值收紧到 1e−6；这仅是附加诊断，不替换主实验。结果如下。','','| 量 | 值 |','|---|---|',f"| 诊断失败数 | {init['failures']}/40 |",f"| 三个阻尼的诊断均值误差 (%) | {', '.join(f'{x:+.6f}' for x in init['mean_error_percent'][3:])} |",f"| 相比原初值，三个阻尼最大变化（百分点） | {', '.join(f'{x:.6g}' for x in init['max_abs_change_percentage_points'][3:])} |",'',
    '如果改初值和收紧停止阈值后仍复现短记录近零均值，就不能简单把该异常归咎于原初值落入一个明显错误的局部解。这个诊断也不构成全局唯一性的证明。','','## 文件与复现','','- `run_task0.py audit`：逐项配置和缓存 DFT 量纲核查。','- `run_task0.py long --workers 1 --count 40`：40 次延长单周期的 M1；逐次 JSON 和 NPZ 落盘。','- `run_task0.py population`：保留原前向差分式总体方程诊断。','- `run_task0.py population_analytic`：同一总体方程，精确 ZOH Fréchet 导数。','- `check_initialization.py`：同一短 P=1 数据的附加初值/容差诊断。','- `build_task0_report.py`：从保存的原始结果重建本报告。','- `long_data/`、`long_results/`、`population/`、`initialization_results/`：全部输入缓存和单次输出；失败没有过滤。','','所有新文件均位于 `followup_DE/task0/`；原始代码、理论、实验数据、原实验 A 输出均只读。']
    (HERE/'task0_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'summary':summaries,'initialization':{k:v for k,v in init.items() if k!='rows'}},ensure_ascii=False),flush=True)

if __name__=='__main__': main()
