"""Read-only consolidation and deterministic Jacobian implementation self-check."""
import check_jacobian as c
import json, hashlib
from pathlib import Path
import numpy as np
import scipy.linalg as la

def main():
    res=c.fixed_residues();f=np.array([3.,10.6,25.,100.,286.573,344.945,490.39,550.])
    active=list(range(6));m,b,ctx,_=c.local_model(f,active,res)
    # Cross-check independently assembled modal/background/pole columns with
    # the original full gmodel Jacobian (keeping only C1=0 model coordinates).
    _,J=m.gmodel(ctx,b,True)
    inds=[]
    for j in active:
        for o in range(4):
            for a in range(4):inds.append(o*m.L+10*j+a)
    for off in (0,4):
        for o in range(4):
            for a in range(4):inds.append(o*m.L+60+off+a)
    inds+=list(range(m.iw,m.nb))
    X,_,_=c.matrices(f,active,res,'exact_tail_free_bg','FRF',None)
    ref=J[:,:,:,inds].reshape(-1,len(inds));ref=np.r_[ref.real,ref.imag]
    colerr=la.norm(X-ref,axis=0)/np.maximum(la.norm(ref,axis=0),1e-300)
    before=json.loads((c.OUT/'source_hashes_before.json').read_text(encoding='utf8'))
    checks=[dict(path=r['Path'],unchanged=hashlib.sha256(Path(r['Path']).read_bytes()).hexdigest().lower()==r['SHA256'].lower()) for r in before]
    result=json.loads((c.OUT/'jacobian_results.json').read_text(encoding='utf8'))
    sel=json.loads((c.OUT/'selected_design_results.json').read_text(encoding='utf8'))['results']
    bound=json.loads((c.OUT/'boundary_design_results.json').read_text(encoding='utf8'))['results']
    derivative=max(r['delay_derivative_relative_error'] for r in result['cases'])
    direct=max(r['sequential_vs_direct_max_abs_error'] for r in sel+bound)
    verify=dict(no_fit=True,no_simulation=True,no_random_generation=True,
        max_nuisance_column_relative_error=float(max(colerr)),
        max_analytic_delay_central_difference_relative_error=derivative,
        max_sequential_vs_direct_projector_absolute_error=direct,
        original_dependencies=checks)
    c.save(c.OUT/'self_check.json',verify)
    assert max(colerr)<1e-12 and derivative<1e-7 and direct<1e-10 and all(x['unchanged'] for x in checks)
    def mainrow(r):return next(t for t in r['thresholds'] if t['cut']==1e-12)
    def arr(a):return ', '.join(f'{v:.6f}' for v in a)
    lines=['# 六模态时延 Jacobian 预检查','',
        '本次只做确定性频响和无噪声 Jacobian 的单次设计检查。未生成时域数据，未调用随机数，未做任何拟合、求根或 Monte Carlo。旧文件只读，五个依赖文件检查前后的 SHA256 全部相同。',
        '', '## 结论', '',
        '在禁止 C1、保留六模态、并使用精确 ZOH 模态基元后，七个结构化时延不再出现旧自由 C0/C1 模型中的几乎完全切空间吸收。按用户预定的“投影其余所有参数（包括其他六个时延）后，任一列能量比低于 1e−3 即不通过”规则，本次全部设计的七列均通过。',
        '',
        '推荐给后续计划的候选为 FC=15、FLP=120、FSC=3，同一六模态对象；L/H 仅估计窗口内三个模态的 C0 和极点，带外模态作为已知精确机械项；时延作用于完整六模态机械响应；D0/D2 保留；六模态边界瞬态均保留，再加入每输出 1、z^{-1}、z^{-2} 共12个延迟寄存器边界候选列。L/H 为153参数，W为207参数。此结论只通过 Jacobian 设计门槛，不等于已验证有限样本精度、拟合收敛或完整时域实现。',
        '', '## 计算定义', '',
        r'对象为 $G_6(z,\tau)=E(\tau)\circ\sum_{i=1}^6 C_{0i}\beta_{0i}^{\mathrm{ZOH}}(z)$，$C_1\equiv0$，真背景为零。',
        '',
        '六个频率为 10.6、11.7、16.1、286.573、344.945、490.39 Hz；阻尼比为 0.03、0.04、0.05、0.0119、0.0132、0.0186。逐通道延迟为 (NO_j+NI_m)/5000 秒，NI=[0,1,0,1]、NO=[1,0,1,0]，公共均值200 μs。七个独立坐标为公共时延、三个正交零和输出 Helmert 对比、三个输入对比；每坐标单位100 μs。',
        '',
        'D 仅从已保存的 P20_r000.npz 读取，取激励谱线并除以 N=600000；其随机相位是既有常量，本轮不生成新相位。由同一六模态真对象代数计算 Uhat=(I+KG6)^{-1}D，未读取旧噪声 Y/U，也未积分状态方程。主检查是输出模型对参数的 M3 参考工具 Jacobian，计算导数时保持 Uhat 固定；不将其误称为带闭环灵敏度左乘的 Tyd 均值导数。FRF 检查则直接堆叠16个传函；FRF不含有限记录瞬态列，因此和输出谱检查分开报告。',
        '',
        r'解析时延列为 $t_j=\operatorname{realstack}\{-s\,[B_j\circ G]\,\widehat U\}$（乘坐标单位）。令 $X$ 为所有非时延 nuisance 列；先计算 $(I-P_X)t_j$，再对剩余六个时延的投影残差作正交投影。报告的主比值为 $r_j=\|(I-P_{[X,T_{-j}]})t_j\|^2/\|t_j\|^2$。',
        '',
        r'还报告七列联合最弱比值 $\min_v\|(I-P_X)Tv\|^2/\|Tv\|^2$，由 $T=Q_TR_T$ 后 $(I-P_X)Q_T$ 的最小奇异值平方得到。它是时延子空间对 nuisance 的最弱角度；其含义与单列条件比值不同，不应互换。',
        '',
        '所有列先做二范数归一化，SVD相对阈值主值1e−12；另查1e−8、1e−10、1e−14。全部主结果没有截掉非零参数方向；直接将“其余全部列”一次SVD投影的复核也一致。',
        '', '## 原提案与模型正确性', '',
        '| 频带 | 原局部参数数 | 原 FC60 min r | 原 FC15 min r | FC60 未表示带外响应相对范数 |',
        '|---|---:|---:|---:|---:|']
    for band in ('L','H','W'):
        old=next(r for r in result['cases'] if r['band']==band and r['case']=='local_free_bg' and r['observation']=='reference_output')
        new=next(r for r in sel if r['controller_fc']==15 and r['band']==band and r['case']=='local_free_bg')
        lines.append(f"| {band} | {old['total_count']} | {min(mainrow(old)['all_other_parameters_delay_energy']):.6f} | {min(mainrow(new)['all_other_parameters_delay_energy']):.6f} | {100*old['missing_model_relative_norm']:.3f}% |")
    lines += ['',
        'L=3–100 Hz，H=100–550 Hz，W=3–550 Hz。上表的遗漏范数在 D0=D2=0 的给定评估点计算，并不是“最优背景拟合之后的残差”。然而有限阶多项式背景不可能精确等于被省略的有理模态尾项，因此局部模型即使 Jacobian 满秩，也不能直接套用正确模型的一致性定理。这里不拟合背景来掩盖这项失配。',
        '',
        '已知带外精确项只修复稳态机械均值。随机有限记录仍有带外模态边界项，所以额外核查了全部六模态瞬态；之后又加入三个延迟边界多项式。未假设“稳态频响正确”就自动意味着完整有限记录观测模型正确。后续仍须对实际时域边界表达做无噪声恒等式验收。',
        '', '## 推荐完整边界候选：FC15', '',
        '| 频带 | 全列秩/参数数 | 七个 r：公共、输出对比1–3、输入对比1–3 | min r | 联合最弱比值 |',
        '|---|---:|---|---:|---:|']
    for r in bound:
        if r['controller_fc']!=15:continue
        t=mainrow(r)
        lines.append(f"| {r['band']} | {r['joint_rank']}/{r['total_count']} | {arr(t['all_other_parameters_delay_energy'])} | {min(t['all_other_parameters_delay_energy']):.6f} | {t['joint_weakest_energy_ratio']:.6f} |")
    lines += ['', '## 两个控制器和边界扩展的对照', '',
        '| FC | 频带 | 参数数 | 满秩 | min r | 联合最弱比值 | 归一化完整 Jacobian 条件数 |',
        '|---:|---|---:|---:|---:|---:|---:|']
    for r in bound:
        t=mainrow(r)
        lines.append(f"| {r['controller_fc']:g} | {r['band']} | {r['total_count']} | {r['joint_rank']} | {min(t['all_other_parameters_delay_energy']):.6f} | {t['joint_weakest_energy_ratio']:.6f} | {r['joint_column_normalized_condition']:.4g} |")
    lines += ['',
        '延迟列足够独立，不表示所有 nuisance 坐标都良态。L/H 加入带外瞬态和低阶边界多项式后，完整列归一化条件数约3.1×10^8/8.8×10^7；它们在窄窗上近似相关，但主阈值下没有精确冗余列。正规方程会平方这一条件数，后续应使用保留完整列空间的 tall QR/SVD 坐标化，不能隐式删除边界参数后再宣称153维原问题已收敛。W约3.4×10^2。',
        '', '## 截断敏感性', '',
        '| FC | 频带 | cutoff | nuisance秩 | min r | 联合最弱比值 |',
        '|---:|---|---:|---:|---:|---:|']
    for r in bound:
        for t in r['thresholds']:
            lines.append(f"| {r['controller_fc']:g} | {r['band']} | {t['cut']:.0e} | {t['nuisance_rank']} | {min(t['all_other_parameters_delay_energy']):.6f} | {t['joint_weakest_energy_ratio']:.6f} |")
    lines += ['',
        '1e−8 时 L 的四个弱 nuisance 方向会被截掉，表中明确保留这一现象。1e−10、1e−12、1e−14均保留完整列空间，主判定一致。因此通过结果不是依赖截断弱 nuisance 得到。',
        '', '## FRF 和实际参考激励', '',
        '| FC60、已知精确带外项+D0/D2 | FRF min r（无瞬态） | 参考输出 min r（窗口内瞬态） | 参考输出 min r（全部6瞬态） |',
        '|---|---:|---:|---:|']
    for band in ('L','H','W'):
        rf=next(r for r in result['cases'] if r['band']==band and r['case']=='exact_tail_free_bg' and r['observation']=='FRF')
        ro=next(r for r in result['cases'] if r['band']==band and r['case']=='exact_tail_free_bg' and r['observation']=='reference_output')
        ra=next(r for r in sel if r['band']==band and r['controller_fc']==60 and r['case']=='exact_tail_free_bg')
        lines.append(f"| {band} | {min(mainrow(rf)['all_other_parameters_delay_energy']):.6f} | {min(mainrow(ro)['all_other_parameters_delay_energy']):.6f} | {min(mainrow(ra)['all_other_parameters_delay_energy']):.6f} |")
    lines += ['',
        '无额外 D0/D2 的干净替代也计算并保存在 jacobian_results.json，但原背景模型已经通过，因此无需为了门槛删除背景。',
        '', '## 自查与理论边界', '',
        f"- 模态/背景/极点列逐列对照原 ModalModel.gmodel，最大相对差 {max(colerr):.3e}。",
        f"- 解析时延方向与确定性中心差分最大相对差 {derivative:.3e}。",
        f"- 分步投影和全其余列一次SVD投影最大绝对差 {direct:.3e}。",
        '- 五个只读依赖文件的 SHA256 全部未改变；没有调用随机数、时域仿真、求根或拟合。',
        '- 本地理论文档第4.1节命题7给出 C1=0 条件下不能由 C1 吸收时延的解释，但其独立性前提仍需本次检查；它不直接保证 r≥1e−3。',
        '- CT 自由 C1 模型中的切空间恒等式不能直接当成 ZOH 模型严格不可辨识定理；本轮通过的是新的 C1=0 六模态设计，不改写旧模型的失败结果。',
        '- 本轮的归一化能量比不包含实际噪声幅值、估计方差或求解稳定性证明。控制器的高频小环路和闭环稳定性另由控制器预检查给出，不能仅凭这里的比值判断。',
        '- 三个边界多项式列是针对当前真值总通道延迟0/1/2个采样周期的候选实现；本轮未证明任意分数时延的边界也由它们精确表示。后续必须先通过实际实现的无噪声模型恒等式验收，再启动拟合或统计实验。',
        '', '## 文件与复现', '',
        '- check_jacobian.py：原 FC60、三种背景/带外结构、FRF与参考输出。',
        '- check_selected_design.py：FC60/FC15、全部6模态瞬态及直接全列投影核对。',
        '- check_boundary_columns.py：再补12个延迟寄存器边界候选列。',
        '- write_report.py：只读汇总与独立 Jacobian 列核查。',
        '- jacobian_results.json / selected_design_results.json / boundary_design_results.json：完整结果、每个阈值和最弱方向（包括对应4×4时延扰动矩阵）。',
        '- self_check.json / source_hashes_before.json / source_hash_verification.json：验证与只读哈希。',
        '',
        '使用现有 Python 3.11，依次运行上述四个脚本即可复现。所有写入仅发生在当前 jacobian_check 目录。']
    (c.OUT/'六模态时延Jacobian预检查.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    manifest=[]
    for p in sorted(c.OUT.glob('*.py')):manifest.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    c.save(c.OUT/'new_script_hashes.json',manifest)
    print(json.dumps(verify,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
