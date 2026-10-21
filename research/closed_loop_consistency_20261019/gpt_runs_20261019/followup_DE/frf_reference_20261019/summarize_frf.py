"""Source-driven summary of the frozen LPM and training-selected LRM diagnostics."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[n]='1'
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,csv,hashlib
import numpy as np
from frf_core import metrics

HERE=Path(__file__).resolve().parent
ZS=('-0.5','0','0.5','1')
METHODS=('M1','M2','M3','M3w','M5','M5n')
norm=np.linalg.norm

def transfer_error(g,r):
    e=r['TY']-g@r['TU']
    return float(100*norm(e)/norm(r['TY']))

def save_csv(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w',encoding='utf-8-sig',newline='') as h:
        wr=csv.DictWriter(h,fieldnames=keys);wr.writeheader();wr.writerows(rows)

def table(header,rows):
    return '\n'.join(['| '+' | '.join(header)+' |','| '+' | '.join(['---']*len(header))+' |']+
                     ['| '+' | '.join(str(x) for x in r)+' |' for r in rows])

def main():
    allrows=[];repeat=[];fitrows=[];sources=set();missing=[]
    for z in ZS:
        old=json.loads((HERE/'results'/f'z{z}_summary.json').read_text(encoding='utf-8'))
        sources.add(old['source'])
        meta={(r['split'],r['method']):r for r in old['fits']}
        for r in old['fits']:
            if not r['available']:
                missing.append(dict(z=z,**r));continue
            sources.add(r['source'])
            obj=json.loads(Path(r['source']).read_text(encoding='utf-8'))
            for i in range(3):
                fitrows.append(dict(z_mm=z,split=r['split'],method=r['method'],mode=i+1,
                    f_hz=obj['frequency_hz'][i],zeta_percent=obj['damping_percent'][i],
                    converged=r['converged'],iterations=r['iterations'],score_normalized=r['score'],
                    iteration_unit=obj.get('iteration_unit','IV_iterations'),
                    initial_M3_failure=obj.get('initial_M3_failure'),source=r['source']))
        for label,folder,suffix in (('LPM',HERE/'results','q2_w21'),('LRM',HERE/'lrm_results','lrm')):
            refs={s:np.load(folder/f'z{z}_{s}_{suffix}.npz') for s in ('full','front','back')}
            a,b,c=(refs[s] for s in ('front','back','full'))
            repeat.append(dict(z_mm=z,reference=label,
                G_halves_percent=float(100*norm(a['G']-b['G'])/norm(c['G'])),
                TU_halves_percent=float(100*norm(a['TU']-b['TU'])/norm(c['TU'])),
                TY_halves_percent=float(100*norm(a['TY']-b['TY'])/norm(c['TY'])),
                G_random_se_percent=float(100*norm(c['se_complex'])/norm(c['G'])),
                median_se_db=float(np.nanmedian(c['se_db'])),
                median_se_phase_deg=float(np.nanmedian(c['se_phase_deg'])),
                max_TU_condition=float(np.nanmax(c['tud_condition'])),
                nonfinite_G_count=int(np.sum(~np.isfinite(c['G']).all(axis=(1,2))))))
            for fit in ('full','front','back'):
                targets=[fit] if fit=='full' else [fit,'back' if fit=='front' else 'front']
                for method in METHODS:
                    mm=meta[fit,method]
                    if not mm['available']:continue
                    g=np.load(HERE/'results'/f'z{z}_{fit}_{method}_model.npz')['G']
                    for target in targets:
                        ref=refs[target];v=metrics(g,ref['G'])
                        allrows.append(dict(z_mm=z,reference=label,fit_split=fit,reference_split=target,
                            method=method,converged=mm['converged'],iterations=mm['iterations'],
                            G_difference_percent=float(v['complex_relative_percent']),
                            inverse_free_difference_percent=transfer_error(g,ref),
                            magnitude_difference_rmse_db=v['magnitude_rmse_db'],
                            circular_phase_difference_rmse_deg=v['phase_rmse_deg']))
    save_csv('comparison_metrics.csv',allrows);save_csv('reference_repeatability.csv',repeat)
    save_csv('model_parameters_and_status.csv',fitrows)
    (HERE/'analysis_summary.json').write_text(json.dumps(dict(repeatability=repeat,comparisons=allrows,
        unavailable_models=missing),ensure_ascii=False,indent=2),encoding='utf-8')
    inventory=[]
    for s in sorted(sources):
        p=Path(s);inventory.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
    (HERE/'source_manifest.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2),encoding='utf-8')

    out=['# 旧实测数据的非参数 FRF 参照检查','',
    '日期：2026-10-19。本次只读取旧实测数据和既有参数拟合，补非参数参照、幅相图、分段检查与失败记录；没有修改理论、旧数据或已有拟合，没有新增蒙特卡洛或硬件实验。',
    '', '## 1. 结论先行','',
    '**这批旧数据能生成有参考价值的 4×4 FRF 图，但当前非参数参照还不能承担“参考 IV 比其他方法更准”的证明。** '
    '局部有理参照降低了多数位置的参照两半差，但z=0由LPM的38.4%略增至LRM的41.1%；四个位置的 LRM 前后半段整体复 FRF 差异仍为约 41%–58%。'
    '这不是硬件模型误差真值，也不能用幅值误差带的中位数较小掩盖整体不确定度。',
    '', 'z=0 时，对全8周期 LRM 参照，M1、M3、M5n 的整体复差分别约 48.5%、28.3%、36.7%；'
    '对应不做矩阵求逆的闭环关系残差约 13.8%、10.9%、12.1%。这些是描述性相容性指标；同一参照前后半差约 41.1%，因此不足以给出准确度排名。',
    '', '论文安排按本轮讨论聚焦：固定周期数、加长单周期的确定性权重一致性作为主线。'
    '本次实测图只能作为相容性和局限性证据，不证明渐近一致性，不补充新的理论首创主张。',
    '', '## 2. 数据入口与能辨识什么','',
    '- 数据为旧 `cl_data_z{-0.5,0,0.5,1}.mat`，不是后来拧紧片簧螺丝后的 sweep 新批次；z=-1 已知坏记录未使用。',
    '- 采样率5 kHz，每周期30000点（6 s）；去掉首周期，使用第2–9周期，6–54 s。前半为第2–5周期，6–30 s；后半为第6–9周期，30–54 s。频带3–25 Hz，133条激励谱线，间隔1/6 Hz。',
    '- 参考量由记录边界恢复：D=U−VQᵀ，四路D逐周期相同到数值精度。边界恒等式和周期性不等价于已证明硬件参考与噪声独立。',
    '- 同一谱线只有一个四维激励方向。重复周期降低随机噪声，却不提供四个独立方向，不能逐频点直接反演完整4×4对象。这里借助相邻频线的局部模型与激励变化恢复矩阵，因此结果依赖局部逼近。',
    '- 参数模型使用之前的固定通道时延。前半模型对后半参照避免直接重复使用同一段随机输出，但全数据先验时延未重新划分，不能称全流程盲测。单位保持记录单位。',
    '', '## 3. 非参数参照与实现核查','',
    '先以D估计闭环矩阵 T_ud 和 T_yd，再计算 G_np=T_yd T_ud⁻¹。非参数估计不输入 SRIVC 的极点、留数、时延或控制器。',
    '', '预设主分析为二阶局部多项式（LPM），21条谱线，跨度3.33 Hz，含局部瞬态项；同时保留一阶15点、二阶31点敏感性结果。'
    '这组结果不理想，完整保留而没有用另一算法替换它。无噪声一般多项式传递矩阵及瞬态的代数核验见 `method_audit.md`。',
    '', '随后增加诊断性局部有理模型（LRM）：对[U−D,Y]分别拟合二阶分子、二阶分母及二阶瞬态，恢复 T_ud 时加已知单位直通。'
    '第4路无反馈的精确结构保留。仅用四个位置的前四周期做留一激励频线预测，移除该频线在所有周期的数据；'
    '比较21/25/31点窗口和线性LRM/三次SK重加权。所有候选使用同一输出标准化权，未按参数模型吻合程度筛选。'
    '全局选中线性LRM、31点（5 Hz跨度），之后冻结。后半周期不参与配置选择。',
    '', '线性LRM含输出构造的回归列，不能宣称自动无偏；5 Hz局部近似、矩阵反演与有限周期噪声仍会影响结果。'
    '方法依据参见 Voorhoeve 等的[闭环局部参数频响估计论文](https://pure.tue.nl/ws/portalfiles/portal/93507539/Paper_MFD_LRM_paper_v501.pdf)，'
    '其局部有理模型用于降低尖锐动态的多项式插值误差，不把它当作实验真值。',
    '', '独立合成的一般局部有理MIMO+瞬态数值核验中，FRF最大绝对误差约4.1×10⁻¹⁷，留出频线预测误差约2.6×10⁻¹⁶。'
    '这排查了回归截距、复数右除及已知直通的实现错误，不证明实测局部模型无偏。全部12个全段/半段LRM和对应周期jackknife无回归失败；大误差与宽散布点全部保留。',
    '', '## 4. 参照自身的重复性','',
    '定义两半差为 100‖G_front−G_back‖F/‖G_full‖F，合并全部频率及16通道。T_ud、T_yd同理。'
    '随机SE比率为100‖SE(G)‖F/‖G_full‖F。LPM为一阶影响量传播，LRM为成对留一周期jackknife；它们不含局部插值偏差。','']
    out.append(table(['z/mm','参照','G两半差/%','T_ud两半差/%','T_yd两半差/%','G随机SE比率/%','最大cond(T_ud)'],
        [[r['z_mm'],r['reference']]+[f"{r[k]:.2f}" for k in ('G_halves_percent','TU_halves_percent','TY_halves_percent','G_random_se_percent','max_TU_condition')] for r in repeat]))
    out += ['', 'LRM全段幅值随机SE中位约0.29–0.39 dB，相位约1.67–2.73°，但个别幅值SE可达30.7 dB。'
        '中位数不能替代最弱通道/频率和整体误差检查。LPM结果尤其受闭环局部插值误差影响，不能把所有误差都归咎于最后求逆。',
        '', 'LRM分解检查进一步发现：闭环T_ud两半差约4.8%–7.6%，T_yd约3.2%–3.9%，但重构机械G的两半差为41%–58%。'
        '99.7%–99.95%的G差异能量集中在8–18 Hz支承模态频带；四个位置均以第二电机输入列贡献最大（52%–69%），第一列其次。'
        '这支持闭环重构在支承峰及弱输入方向放大差异的解释，不能由此区分局部估计误差、真实漂移与参数模型失配。'
        '频段和列分解仅用于解释，全部频点继续保留。详见 `lrm_results/diagnostic_bands.md`。',
        '', '## 5. 全记录参数模型与参照','',
        'M1=全部频点直接IV；M2=谱线直接IV；M3=谱线参考IV；M3w=邻非激励频点加权参考IV；M5=周期样本协方差SML；M5n=邻点协方差SML。'
        'M1复用实验C，其他方法复用followup D，均为同一ZOH模型类。原迭代容差、步长和终止定义不完全相同，未重跑成统一求解器。',
        '', '表内“G差”为相对参照的复Frobenius差；“闭环关系残差”为100‖T_yd−G_model T_ud‖F/‖T_yd‖F，避免对估计T_ud求逆。'
        '它仍受T_ud/T_yd估计误差影响，不是真值损失。所有失败端点纳入表格并标记。','']
    for label in ('LPM','LRM'):
        out += [f'### {label} 参照','',table(['z/mm','方法','原收敛状态','G差/%','闭环关系残差/%','幅差RMSE/dB','圆周相位差RMSE/°'],
            [[r['z_mm'],r['method'],'是' if r['converged'] else '**失败**',f"{r['G_difference_percent']:.2f}",f"{r['inverse_free_difference_percent']:.2f}",
              f"{r['magnitude_difference_rmse_db']:.2f}",f"{r['circular_phase_difference_rmse_deg']:.2f}"]
             for r in allrows if r['reference']==label and r['fit_split']==r['reference_split']=='full']), '']
    out += ['M5在z=0和1的模型与LRM参照、以及求逆前的闭环关系都有较大差异，不能因为旧优化器返回收敛就当作可靠FRF；'
            'z=0.5的M5原本失败，z=-0.5的M1原本失败，均未删除。完整频率/阻尼、原score和终止信息见 `model_parameters_and_status.csv`。'
            '这里没有重拟合SML或事后调整其参数。',
            '', '独立核对排除了M5绘图和FFT尺度错误：SML内部提取的G与本次resp求值相对差不超过3.3×10⁻¹⁶；'
            '周期DFT换到全长DFT只需同比放大瞬态系数，物理FRF完全不变。z=0的未加权输出相对残差从初值6.54%增至271.78%，'
            '同时SML目标由219500降至153646；z=1对应32.32%→332.21%，目标728745→209811。'
            '残差协方差条件数中位数约10⁵–9×10⁵。可见旧SML目标和当前数值端点没有保证良好的输出/FRF验证。'
            '保留原收敛标记，同时记为“数值终止但验证异常”；不能把这两例推广成SML理论无效。核查代码和结果见 `verify_lpm.py`、`verify_lpm_result.json`。',
            '', '## 6. 分段检查','',
            '下面只展示前半参数模型对后半LRM参照，非同段拟合分数。其他方法、同段与反向检查完整保存在 `comparison_metrics.csv`。'
            'M1后半模型原结果缺失，因此反向M1标注不可用；本次没有为补表额外重跑参数辨识。','']
    out.append(table(['z/mm','前半模型','原收敛状态','对后半G差/%','闭环关系残差/%'],
        [[r['z_mm'],r['method'],'是' if r['converged'] else '**失败**',f"{r['G_difference_percent']:.2f}",f"{r['inverse_free_difference_percent']:.2f}"]
         for r in allrows if r['reference']=='LRM' and r['fit_split']=='front' and r['reference_split']=='back' and r['method'] in ('M1','M3','M5n')]))
    out += ['', '不能把此分段差异当作独立重复实验的置信区间；只有一条连续记录，两段可能共享漂移与固定标定误差。'
        '全模型16通道均保留，弱通道没有按噪声或拟合优劣剔除。',
        '', '## 7. 论文可以怎么用，不能怎么用','',
        '- 可以：展示参数模型与纯数据局部FRF的幅相对照，明确3–25 Hz支承模态频带、参照算法和窗口；用分段散布解释哪些通道证据弱。',
        '- 可以：把本次结果作为旧数据诊断，指出“参数估计收敛”与“完整传递矩阵得到可靠验证”是两件事。',
        '- 不可以：据此宣称参考IV实测更准确、SML理论不一致，或把两者与参照的差直接当作辨识偏差。实测没有真值，参照本身散布较大。',
        '- 不可以：只画看起来较好的LRM图，隐去预设LPM结果、SML失败或极弱通道。',
        '', '更强的硬件参照需在每个频率提供至少4个线性独立的输入方向，例如4次独立相位组合的四路激励，或逐执行器激励，并同步记录参考D、总U和原始Y。'
        '在同一硬件状态下切换控制器并重复，将直接检验方法是否受反馈改变影响。增加相同周期的重复数只能改善平均，不会增加每频率输入秩。'
        '这些是后续实验设计说明，本次没有实施新实验。',
        '', '## 8. 图、复现与范围变更','',
        '- LPM图：`figures/z0_full_vs_full_main_magnitude.png`、`..._phase.png`，各位置相同命名；完整方法与分段图均保留。',
        '- LRM诊断图：文件名前缀`lrm_`，图注明“diagnostic”。独立unwrap并仅作整周显示平移，误差指标始终从原始复频响计算。',
        '- 新增LRM和求逆前残差，是因预设LPM参照重复性不足而追加的诊断；没有用参数模型选择窗口，没有替代预设结果。',
        '- 本次没有补算缺失M1反向半段模型，明确留空；有可用前半M1完成主要分段对照。',
        '- 源结果和审计：`results/`、`lrm_results/`、`data_audit.md`、`method_audit.md`；输入哈希在`source_manifest.json`。',
        '- 运行命令见`README.md`。所有代码、派生数据、图和日志限定在本目录。', '']
    (HERE/'旧数据非参数FRF参照报告.md').write_text('\n'.join(out),encoding='utf-8')
    print(json.dumps(dict(comparisons=len(allrows),reference_rows=len(repeat),parameter_rows=len(fitrows),
                         missing_models=len(missing)),ensure_ascii=False))

if __name__=='__main__':main()
