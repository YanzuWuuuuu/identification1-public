"""Rebuild the follow-up report exclusively from saved trial records."""
import sys,json
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from weighted_core import statistics
HERE=Path(__file__).resolve().parent
PARAMS=['f₁','f₂','f₃','ζ₁','ζ₂','ζ₃']
PRIMARY=['M2','M3','M3w_n16_o1','M3o','M5','M5n_n16']
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def dump(p,o):p.write_text(json.dumps(o,indent=2,ensure_ascii=False,allow_nan=True),encoding='utf-8')
def fmt(x,d=4):return f'{x:.{d}f}' if x is not None and np.isfinite(x) else 'NA'
def triple(x):return '/'.join(fmt(v,3) for v in x)
def stat_table(groups,subset):
    lines=['| P | 方法/策略 | 参数 | 均值±SE (%) | SD (%) | RMSE (%) | 纳入数/40 | 失败/40 | 平均迭代或 nfev |',
           '|---:|---|---|---:|---:|---:|---:|---:|---:|']
    for (p,m),s in sorted(groups.items()):
        a=s[subset]
        for j,name in enumerate(PARAMS):
            if not a['n']:vals=['NA']*3
            else:vals=[f'{a["mean"][j]:+.4f} ± {fmt(a["se"][j])}',fmt(a['sd'][j]),fmt(a['rmse'][j])]
            lines.append(f'| {p} | {m} | {name} | '+ ' | '.join(vals)+f' | {a["n"]}/40 | {s["failures"]}/40 | {fmt(s["mean_iterations"],2)} |')
    return '\n'.join(lines)

def main():
    simrows=[read(p) for p in sorted((HERE/'sim_d'/'results').glob('*.json'))]
    erows=[read(p) for p in sorted((HERE/'short_e').glob('P*.json'))]
    rawgroups={};egroups={}
    for r in simrows:rawgroups.setdefault((r['P'],r['method']),[]).append(r)
    for r in erows:egroups.setdefault((r['P'],r['method']+'_'+r['strategy']),[]).append(r)
    assert len(simrows)==1600 and len(rawgroups)==40 and all(len(v)==40 for v in rawgroups.values()),'Incomplete D records'
    assert len(erows)==360 and len(egroups)==9 and all(len(v)==40 for v in egroups.values()),'Incomplete E records'
    gs={key:statistics(v) for key,v in rawgroups.items()};es={key:statistics(v) for key,v in egroups.items()}
    dump(HERE/'short_e'/'summary.json',[dict(P=k[0],method_strategy=k[1],**s) for k,s in sorted(es.items())])
    rules={}
    for subset in ['converged_only','all_finite']:
        details=[]
        for p in [2,4,8,16]:
            w=gs[p,'M3w_n16_o1'];m=gs[p,'M5n_n16'];ratio=np.array(w[subset]['rmse'])/np.array(m[subset]['rmse'])
            details.append(dict(P=p,ratios=ratio.tolist(),m3w_failures=w['failures'],m5n_failures=m['failures'],passes= ((ratio<=1.1)&(w['failures']<=m['failures'])).tolist()))
        A=all(all(x['passes']) for x in details if x['P'] in (2,4))
        B=all(all(np.array(x['ratios'])>1.1) for x in details)
        rules[subset]=dict(rule_A=A,rule_B=B,verdict='A' if A else 'B' if B else 'other',details=details)
    dump(HERE/'decision_rules.json',rules)
    task0=(HERE/'task0'/'task0_report.md').read_text(encoding='utf-8')
    task0='\n'.join(task0.splitlines()[2:]).replace('\n## ','\n### ')
    sections=['# 后续实验 D/E：闭环频域 SRIVC、参考 IV 与 SML\n',
              '## 1. 任务 0：P=1 异常排查（先行完成）\n\n'+task0]
    sections.append('''## 2. 仿真 D：预先固定的比较与完整统计

本轮仅研究三个支承模态，未增加柔性模态、零空间算法或间接 FRF。参数、控制器、噪声、周期长度和种子沿用实验 A；P=2/4/8/16 均为 40 次，每条数据供全部方法共用。原 A 的 P=16 只有 20 条，新增的 r=20…39 按同一公式生成并单独保存。新生成数据未覆盖旧缓存。

主比较提前固定为 **M3w_n16_o1 对 M5n_n16**：16 个最近非激励点、一次权更新。8 点和第二次更新不用于事后挑选获胜算法。M3 本轮只用激励谱线，与原 A 的 ref_band_W1（全部频点）有区别。

M2：直接法谱线；M3：均匀权参考 IV；M3o：真噪声满矩阵权参考 IV；M3w_nL_oI：L 个邻点、I 次权更新；M5：原周期 SML（P<5 时 α=.1 联合收缩，P≥5 时原始样本协方差）；M5n_nL：非激励点联合协方差 SML。所有 IV 保持相同瞬态参数和精确 ZOH。SML 沿用原有含协方差导数的 profile 目标求解器，没有改成冻结协方差的伪似然迭代。

表内误差统一为 **100[exp(log估计−log真值)−1]**。SE 是均值标准误，SD 是重复估计标准差，RMSE 是对真值的均方根误差。主表列收敛样本，但失败数量同时保留；全部有限终点的同格式统计见第 9 节。这样不能把筛掉失败后的低 RMSE 单独当作算法优势。

IV 的失败判据沿用 A：异常、非有限或达到第 80 次；停止阈值为极点对数步 <1e−5。SML 沿用原优化器成功标志和上限；其次数为函数求值 nfev，**不能与 IV 步数直接解释为效率比较**。满参数归一化 score 单独记录，未用来事后重定义失败。
''')
    quick=['### 主方法的阻尼 RMSE（收敛样本）\n','| P | 方法 | ζ₁/ζ₂/ζ₃ RMSE (%) | 失败/40 |','|---:|---|---|---:|']
    for p in [2,4,8,16]:
        for m in PRIMARY:
            s=gs[p,m];quick.append(f'| {p} | {m} | {triple(s["converged_only"]["rmse"][3:])} | {s["failures"]}/40 |')
    sections.append('\n'.join(quick)+'\n\n### 全方法、全参数统计（收敛样本）\n\n'+stat_table(gs,'converged_only'))
    gate=read(HERE/'sim_d'/'diagnostics'/'scale_parity_gate.json')
    sections.append(f'''## 3. SML 尺度校验与局部权敏感性

### 3.1 尺度校验先于批量比较

未归一化全长 DFT 的激励线 X_full=P·X_period_mean，因此 Cov(X_full)=P²·Cov(X_period_mean)，而单周期 DFT 与其 P 次均值的协方差相差 P。三种量不能混用。P=16 的独立校验：全长谱与 P 倍周期均值相对差 {gate['full_DFT_equals_P_periodmean_relative']:.3e}，正确换算后的 SML 残差相对差 {gate['residual_scale_relative']:.3e}，雅可比相对差 {gate['jacobian_scale_relative']:.3e}，重新求解的模态 log 参数最大差 {gate['estimate_scale_modal_log_difference']:.3e}。瞬态系数也随输出 DFT 单位变为 P 倍；结构参数不变。数值证明了尺度一致，**并不强制不同协方差估计的 M5n 和 M5 给出相同估计**。

### 3.2 M5n 与 M5 在 P=16 的配对比较
''')
    pairlines=['| 参数 | M5n−M5 的平均差±SE（误差百分点） | 配对差 SD | M5 RMSE (%) | M5n RMSE (%) |','|---|---:|---:|---:|---:|']
    rr5={r['replicate']:r for r in rawgroups[16,'M5']};rrn={r['replicate']:r for r in rawgroups[16,'M5n_n16']}
    pairs=[np.array(rrn[r]['error_percent'])-rr5[r]['error_percent'] for r in rr5 if rr5[r]['success'] and rrn[r]['success']]
    dif=np.array(pairs)
    for j,p in enumerate(PARAMS):pairlines.append(f'| {p} | {dif[:,j].mean():+.4f} ± {dif[:,j].std(ddof=1)/np.sqrt(len(dif)):.4f} | {dif[:,j].std(ddof=1):.4f} | {gs[16,"M5"]["converged_only"]["rmse"][j]:.4f} | {gs[16,"M5n_n16"]["converged_only"]["rmse"][j]:.4f} |')
    sections.append('\n'.join(pairlines)+f'\n\n配对均收敛数 {len(pairs)}/40；全部失败仍在总表内。')
    sections.append('''### 3.3 为什么同源噪声信息不等于同一数值权重

M3w 先在各邻点 j 用 G(f_j) 形成残差再平均；M5n 先平均 [Y−T,U] 联合协方差，再在中心 k 用 G(f_k) 投影。当邻域内频响变化快，二者即使没有收缩也不同。此外，M3w 收缩残差 4×4 协方差，M5n 收缩联合 8×8 协方差，收缩与投影不交换。因此本比较是按要求使用相同来源频点，不能宣称噪声权或正则化完全相同。

独立审计 r=0 的 n_l=16 原始协方差差异（全频 Frobenius 相对量）在 P=2/4/8/16 分别约 47.51%/27.29%/6.26%/2.61%；包含两种收缩映射后约 49.16%/29.32%/10.08%/5.87%。P=2 的邻域半宽中位数约 1.25 Hz，已跨越较窄的支承共振。这是少周期下的局部平稳近似限制，不能误归因为 P 或 P² 的尺度错误。固定中心 G 后的联合投影恒等式已核对到 2e−16 量级。

非激励谱线未混入任何激励谱线，但权重还依赖初步 β̂，有限窗口不同 DFT 频点也未必严格独立；固定 8/16 个邻点并不构成一致谱估计的渐近条件。因此这里只报告可重复的有限样本表现，不把这些计算改写成估计权重下的严格无偏证明。

### 3.4 8/16 邻点及一次/两次更新
''')
    sens=['| P | 邻点数 | 更新次数 | f₁/f₂/f₃ RMSE (%) | ζ₁/ζ₂/ζ₃ RMSE (%) | 失败/40 |','|---:|---:|---:|---|---|---:|']
    for p in [2,4,8,16]:
        for n in [8,16]:
            for o in [1,2]:
                s=gs[p,f'M3w_n{n}_o{o}'];a=s['converged_only'];sens.append(f'| {p} | {n} | {o} | {triple(a["rmse"][:3])} | {triple(a["rmse"][3:])} | {s["failures"]}/40 |')
    sections.append('\n'.join(sens)+'\n\n更新二次不是自动更好；不据此替换已固定的一次更新主方法。各参数偏差、SD 和 SE 均已列在第 2 节。')
    theory=read(HERE/'sim_d'/'oracle_theory_summary.json')
    sections.append('''## 4. M3o：理论标准差与 MC 标准差

用各次实际多正弦参考信号构造真参数处的 Z，保留全部 126 个参数（包括 24 个瞬态 nuisance 参数），计算 A=Re∑ZᴴWZ，Q=(1/2)Re∑ZᴴWΦ_vWZ，C=A⁻¹QA⁻ᵀ。Φ_v=N(σ_v²I+σ_w²G₀G₀ᴴ)，1/2 来自实参数的正频复谱方差。先完成全参数消元，再取模态子块；没有把瞬态参数当作已知。

主列把 log 参数渐近高斯方差按 exp−1 变换到百分比 SD：对方差 v，条件方差为 10000·exp(v)[exp(v)−1]，并计入不同参考设计的条件均值差异。100√v 的一阶近似也保存在 JSON。该预测仍是局部渐近近似，不含参数依赖 IV 求解的完整有限样本效应或收敛筛选；它与 MC 不吻合时不声称协方差公式已被精确验证。
''')
    tt=['| P | 参数 | 理论 SD (%) | MC SD，收敛样本 (%) | MC/理论 | M3o失败/40 |','|---:|---|---:|---:|---:|---:|']
    for row in theory:
        th=row['theory_percent_sd'];mc=row['converged_only']['sd']
        for j,p in enumerate(PARAMS):tt.append(f'| {row["P"]} | {p} | {th[j]:.4f} | {mc[j]:.4f} | {mc[j]/th[j]:.3f} | {row["failures"]}/40 |')
    sections.append('\n'.join(tt)+'\n\n独立数值审计中，oracle 的 Q=(1/2)A 恒等式相对误差约4e−13，全矩阵逆与消去120个 nuisance 坐标后的模态 Schur 协方差相对差约3e−8。P16/r0 的缩放信息矩阵条件数仍约3.3e9。理论与MC散布接近支持一阶协方差计算，但不构成有限样本精确覆盖率验证。')
    real=(HERE/'real_d'/'实测D_结果与半段差异.md').read_text(encoding='utf-8')
    real='\n'.join(real.splitlines()[2:]).replace('\n## ','\n### ')
    sections.append('## 5. 实测：完整窗口与前后两半\n\n'+real)
    sections.append('''## 6. 任务 E：少周期初值与阻尼求根

E1：原模态扰动及线性初始化；E2：相同数据的 M2 估计作为初值；E3：E2 加 LM 型 IV 矩方程求根。E3 优化的是参考 IV 的矩残差，未改成输出误差最小二乘，也没有参数正则化罚项。其矩雅可比包括参考输入随 β 改变的导数和输出模型二阶导数收缩；外层模态差分步为 1e−4，方向导数核查约2e−6，独立方向复查在2e−5以内。LM 拒绝不能降低当前尺度下实际矩残差的步，通过增大阻尼将模态对数信赖半径限制为0.25（原步长截断0.05），未截断已接受步。因此E3是包含完整矩导数与信赖域尺度变化的求根策略，不能把所有变化只归因于“加了阻尼”。

M3w 的三种策略按完整流程执行：先用该策略求均匀权 M3，再从其端点估权并用同一策略求解加权方程。因此这是初始化/求解流程稳健性比较，权重也可能随初步估计变化，不能解释为保持完全相同权矩阵后的纯求解器比较。主规格仍是 n_l=16、一次更新。初始 M3 失败的有限端点允许继续构造权重，并保留上游失败记录。

原 A 80 步和极点步阈值保持不变；E3 若不存在可接受的下降步也记为求根失败。所有结果另列满参数 score，避免把受阻尼压小的极点步误当作严格矩方程精度保证。初始简化雅可比的一个试运行未进入正式统计，保留在 short_e/pilot_quasi_jacobian/；正式 E3 已核对完整矩导数。
''')
    overview=['| P | 方法/策略 | 失败/40 | ζ₁/ζ₂/ζ₃ RMSE，收敛样本 (%) | 平均迭代 |','|---:|---|---:|---|---:|']
    for (p,m),s in sorted(es.items()):overview.append(f'| {p} | {m} | {s["failures"]}/40 | {triple(s["converged_only"]["rmse"][3:])} | {s["mean_iterations"]:.2f} |')
    sections.append('\n'.join(overview)+'''\n\n**关键负面结果：** P=2 的 M3w_E3 有39/40次失败，只有1个收敛端点；该行收敛子集的RMSE没有总体比较意义，SD/SE也无法由一个样本估计。不能据此声称LM改善了少周期稳健性。本次E2通常减少迭代，但并未消除P=1的大散布。E3的结果仅针对这里明确记录的矩残差缩放、LM实现和80步预算，不能推论所有LM/信赖域算法都无效。

E1在P=1的失败数由原A的6变成7，差异集中在r=29的80步边界：原端点score约2.46e−5，新端点约4.19e−6，新的严格score反而更小。两者初值和数据相同，满矩阵通用实现改变了浮点累加路径；这体现原极点步停止判据在弱方向上的敏感性，不是另一个统计偏差机制。未修改判据以消除这一差异。

独立检查r=0加权E3终点：缩放矩雅可比条件数约3.16e9，未阻尼正规矩阵约5.44e18；但在实际使用的阻尼下，以增广最小二乘替代正规方程求步时，模态步差仅约1e−11，实际score几乎相同。因此没有证据把39/40次上限简单归因于正规方程舍入误差；弱方向、局部非线性、矩残差尺度和保守阻尼在80步预算内的推进速度仍是限制。本轮未事后更换E3并删去这一失败批次。

### E 的全参数统计（收敛样本）

'''+stat_table(es,'converged_only'))
    ea=read(HERE/'audit_e.json');et=['### E 的停止判据与完整矩残差\n','| P | 方法/策略 | 失败/40 | 原规则成功但score>1e−5 | 成功端点最大score |','|---:|---|---:|---:|---:|']
    for r in ea['groups']:et.append(f'| {r["P"]} | {r["method"]}_{r["strategy"]} | {r["failures"]}/40 | {r["converged_score_gt_1e5"]} | {r["max_converged_score"]:.3e} |')
    sections.append('\n'.join(et)+'\n\nP2加权E3唯一的原规则成功端点，score也高于1e−5，故本轮未获得其稳健求根改善证据。这里只补充诊断，不事后修改预定失败定义。')
    sections.append('## 7. 按预先固定规则判定\n\n规则逐项作用于六个参数，不把频率和阻尼合并成一个可事后调权的综合分数。整体 A 要求 P=2 和 4 的全部六项满足且失败数不多；整体 B 要求所有 P 的全部六项均超过 1.1。否则只给分项结论。')
    rl=['| P | f₁ RMSE比 | f₂比 | f₃比 | ζ₁比 | ζ₂比 | ζ₃比 | M3w/M5n失败数 |','|---:|---:|---:|---:|---:|---:|---:|---|']
    for r in rules['converged_only']['details']:rl.append(f'| {r["P"]} | '+' | '.join(f'{x:.3f}' for x in r['ratios'])+f' | {r["m3w_failures"]}/{r["m5n_failures"]} |')
    sections.append('\n'.join(rl))
    verdict=rules['converged_only']['verdict']
    if verdict=='A':text='少周期下，带非参数权的参考 IV 与同等信息的 SML 精度相当，且不需要周期样本协方差。'
    elif verdict=='B':text='参考 IV 的作用仅限于消除直接法偏差，精度不如 SML。'
    else:text='预设整体规则 A、B 均未满足。不能写“少周期精度相当”的整体结论，也不能将全部周期和参数概括为统一劣于 SML；按上表逐项报告。'
    sections.append('**判定：'+text+'**\n\n具体为：P=2的三个阻尼RMSE比分别1.139/1.249/1.176，均超出1.1；P=4与16的六个参数均满足1.1范围；P=8只有f₁略超出（约1.110）。全部有限端点的敏感性判定为 '+rules['all_finite']['verdict']+'，详见 decision_rules.json。M3w 相对 M3 的显著改进，说明噪声矩阵加权有实际价值；这不等于已经证明它相对 SML 具有精度或效率优势。单周期无非激励点，M3w/M5n 不适用；参考 IV 虽能建立方程，少周期求根失败和大散布仍需正面报告。实测只评价一致性与收敛。')
    sections.append('''## 8. 与任务字面设定的差异及边界

1. 主方法的 n_l、权更新次数以及六参数规则聚合方式在运行前写入 PROTOCOL.md；这些原任务未指定的选择没有根据结果改变。
2. D 的 M3 按此次要求只用激励谱线；旧 A 中名为 ref_band_W1 的 M3 使用全部频点。E1/新矩阵求解与旧代码代数等价但浮点求和路径不同，弱方向上可影响边界收敛计数；完整旧结果仍保留。不能把原 A 的所有 M3 数字直接当作本轮谱线 M3。
3. P=16 从旧 A 的 20 次补到 40 次；新增数据依照相同种子公式生成于新目录。Task0 长单周期仍 settle=2，所以丢弃时间随周期长度从 12 s 增为 96 s；有效数据长度与普通 P8 一致，激励谱线数量/分辨率不同。
4. M5 沿用已有 derivative-aware profile SML，而非再次改成冻结 C_e 的步进算法；这是复用既有实现的要求。P<5 的收缩不是原始无正则 SML，已分别注明。M5n 对联合协方差采用 α=.1 收缩，M3w 对残差协方差收缩；局部频率变化和收缩映射不同的影响未隐瞒。
5. 采用全长 DFT 单位运行 M5n，以避免混用周期均值协方差；这与正确缩放后的周期均值求解在代数上等价。M5 使用周期均值，瞬态系数单位随之变化。
6. 理论标准差采用允许的等价满矩阵 sandwich 公式，按每条实际参考设计计算并保留瞬态 nuisance 参数，未只复制旧均匀权的结构参数近似。
7. 实测 IV 停止阈值统一为 A 的 1e−5/步幅0.05；前次 C 是 1e−6/0.1。SML 实测保留 160 次函数求值上限（仿真120），没有混写成相同 IV 步数。
8. P2 的两半各只有一个周期，M3w/M5/M5n 无法从该半段构造规定协方差，明确为不适用，不借用完整窗口权重，也不计作数值失败。重复的完全相同实测窗口复用一次求解，但所有180个报告上下文均展开。
9. E3 包含完整矩导数后做阻尼根求解；不是给原输出误差目标加 LM。不能下降的终点保留并记为失败。正式运行前的简化雅可比试验单独归档，未计入40次主统计。
10. 非激励谱线与激励谱线不同并不能自动证明估计权与残差严格独立；β̂ 的反馈依赖、有限窗口相关性和固定邻点数的谱估计误差仍存在。此次未改写原理论或宣称新增严格一致性定理。
11. Task0 的总体谱矩计算补齐了原 exact_limit.py 未使用的周期线掩码，并以矩阵指数 Fréchet 导数消除仅在该诊断中的有限差分收敛平台。仿真拟合模型/算法未因此改变；原诊断端点也保留。

SML 与参考 IV 的定位：原始周期 SML 需要足够重复周期估计协方差；M5n 与 M3w 都把这一需求换成非激励频点上的局部谱估计，并不意味着免去噪声统计估计。当前 M3/M3w/M3o 都使用给定控制器生成参考输入，本轮没有重新验证不用控制器的 M4。参考 IV 的消偏能力与其有限样本最优精度、求解稳定性是不同结论，应分别用上述数据判断。没有新增间接 FRF 或多实验流程的效率比较。
''')
    sections.append('## 9. 全部有限终点统计（包含数值失败，禁止选择性删样）\n\n### D\n\n'+stat_table(gs,'all_finite')+'\n\n### E\n\n'+stat_table(es,'all_finite'))
    sections.append('## 10. 可追溯文件与完整性\n\n- 原始逐次结果：sim_d/results/（1600 条）、short_e/（360 条）、task0/long_results/（40 条）。\n- 实测：real_d/raw/（140 个唯一记录，含24个不适用），summary.json 展开180个上下文。\n- 满矩阵/尺度校验：sim_d/diagnostics/、validation_e.json、real_d/diagnostics/；独立数学审计见 audit_math.md 与 audit_e.md。\n- 理论协方差逐次 NPZ：sim_d/diagnostics/*oracle_theory.npz。\n- 非有限结果和失败状态保留于 JSON；无参考真值的实测不标注准确度。\n- 原文件哈希：original_inventory.json / original_integrity_check.json；本报告脚本和所有新输出位于 followup_DE。\n- 复现命令见 README.md。')
    (HERE/'后续实验D_E_报告.md').write_text('\n\n'.join(sections)+'\n',encoding='utf-8')
    dump(HERE/'package_validation.json',dict(simulation_records=len(simrows),simulation_groups=len(gs),e_records=len(erows),e_groups=len(es),real_context_rows=read(HERE/'real_d'/'summary.json')['context_rows'],rules=rules,report_bytes=(HERE/'后续实验D_E_报告.md').stat().st_size))
    print(json.dumps({'D':len(simrows),'E':len(erows),'rules':rules['converged_only']['verdict']},ensure_ascii=False))
if __name__=='__main__':main()
