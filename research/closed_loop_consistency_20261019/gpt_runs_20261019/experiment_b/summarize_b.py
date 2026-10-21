import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
METHODS=('dir_band_W1','dir_lines_W1','ref_band_W1')
LABELS=('M1 全频点直接法','M2 谱线直接法','M3 参考 IV')
FM=np.array([10.6,11.7,16.1,286.573,344.945,490.39,718.958,755.501])
ZM=np.array([.03,.04,.05,.0119,.0132,.0186,.0322,.0346])

def main():
    methods={}; development=[]; full_data={}
    for label,method in zip(LABELS,METHODS):
        rows=[]; metas=[]
        for rep in range(20):
            base=HERE/'raw'/f'rep{rep:03d}_{method}.npz'
            retry=HERE/'raw'/f'rep{rep:03d}_{method}_retry1.npz'
            if retry.exists():
                with np.load(base) as old:
                    development.append(json.loads(str(old['metadata_json'])))
                path=retry
            else: path=base
            if not path.exists(): continue
            with np.load(path) as a:
                rows.append(a['error_percent']); metas.append(json.loads(str(a['metadata_json'])))
        data=np.array(rows); good=np.all(np.isfinite(data),axis=1) if len(data) else np.zeros(0,bool)
        full_data[method]=data
        status=np.array([m['status']=='pole_step_converged' for m in metas])
        vals=data[good]; n=len(vals)
        stats={name:value.tolist() for name,value in dict(mean=np.mean(vals,axis=0),
            std=np.std(vals,axis=0,ddof=1),se=np.std(vals,axis=0,ddof=1)/np.sqrt(n),
            rmse=np.sqrt(np.mean(vals**2,axis=0)),min=np.min(vals,axis=0),max=np.max(vals,axis=0)).items()} if n else {}
        vc=data[good & status]; nc=len(vc)
        convstats={name:value.tolist() for name,value in dict(mean=np.mean(vc,axis=0),
            std=np.std(vc,axis=0,ddof=1),se=np.std(vc,axis=0,ddof=1)/np.sqrt(nc),
            rmse=np.sqrt(np.mean(vc**2,axis=0)),min=np.min(vc,axis=0),max=np.max(vc,axis=0)).items()} if nc else {}
        methods[method]=dict(label=label,n_attempts=len(rows),n_finite=n,failures=int(np.sum(~status)),
            iteration_min=min((m['iterations'] for m in metas),default=None),
            iteration_max=max((m['iterations'] for m in metas),default=None),
            seconds_mean=float(np.mean([m['duration_seconds'] for m in metas])),
            seconds_total=float(np.sum([m['duration_seconds'] for m in metas])),
            score_max=float(max([m['score_relative'] for m in metas],default=float('nan'))),
            stats_percent=stats,n_converged=nc,stats_converged_percent=convstats,records=metas)
    paired={}
    ref=full_data['ref_band_W1']
    if len(ref)==20:
        for method in METHODS[:2]:
            a=full_data[method]
            if len(a)==20:
                ok=np.all(np.isfinite(a),axis=1)&np.all(np.isfinite(ref),axis=1)
                dif=a[ok]-ref[ok]
                assert np.allclose(dif.mean(0),a[ok].mean(0)-ref[ok].mean(0),rtol=1e-12,atol=1e-12), 'Paired means must equal differences of means on identical repetitions.'
                paired[method+' minus ref_band_W1']=dict(n=int(ok.sum()),mean_percent=dif.mean(0).tolist(),se_percent=(dif.std(0,ddof=1)/np.sqrt(ok.sum())).tolist())
    out=dict(n_mc_requested=20,frequencies_hz=FM.tolist(),damping_ratio=ZM.tolist(),paired_differences=paired,
             methods=methods,development_error_count=len(development),development_errors=development)
    (HERE/'summary.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    s=['# 实验 B：完整八模态对象的闭环对照', '',
       '本实验扩展原 `cl2.py`，不修改对象、控制器、噪声和估计方程。三个方法在每次重复中使用同一份数据及同一组模态初始扰动。', '',
       '- 对象：3 个支承模态 + 5 个柔性模态，ZOH 精确求值；离散采样 5 kHz。',
       '- 数据：每周期 30000 点、10 周期、N=300000，额外 2 周期作稳态预热；20 次蒙特卡洛。',
       '- 频带：3–800 Hz；M1/M3 为 47821 个 DFT 频点，M2 为 4783 条激励谱线。',
       '- 控制器：FC=60 Hz，FLP=120 Hz；AL=30、AH=100、SW=0.03、sd=0.2、sv=3e-7。',
       '- 参数：256 个自由留数参数 + 16 个对数频率/阻尼参数 + 64 个逐模态瞬态参数，共 336 个。',
       '- 初值：频率扰动 `log(1+0.01 ξ)`，阻尼扰动 `0.4 ξ`，线性留数和瞬态由同一数据的线性最小二乘初始化。',
       '- 随机种子：52061019+r，r=0,...,19；三种方法以百分比 `100[exp(log估计−log真值)−1]` 统计。',
       '', '## 计算正确性与运行记录', '',
       '为减少原稠密雅可比中的零元素运算，按输出通道分块累积原 IV–GN 正规方程，再对 16 个公共非线性参数作 Schur 消元。此处只是线性代数重排，不是修改估计器。与原 `cl2.out_jac` 的一次更新对照：输出、雅可比完全一致；缩放正规矩阵误差 5.78e-16，梯度误差 1.47e-15，模态步长相对误差 1.32e-12。闭环谱半径 0.9964659132 < 1。',
       '', '| 方法 | 完成/请求 | 失败 | 迭代范围 | 单次平均耗时/s | 最大归一化矩残差 |',
       '|---|---:|---:|---:|---:|---:|']
    for method in METHODS:
        m=methods[method]
        s.append(f"| {m['label']} | {m['n_attempts']}/20 | {m['failures']} | {m['iteration_min']}–{m['iteration_max']} | {m['seconds_mean']:.2f} | {m['score_max']:.3g} |")
    s+=['','失败指异常、非有限估计或达到 80 次上限。沿用原代码的收敛条件为模态对数参数最大步长 <1e-5；另外在返回解处重新计算全参数矩残差。表中残差为各列矩与工具列范数、总残差范数归一化后的最大值，仅作数值诊断，不能单独证明统计一致性。所有有限终值都纳入统计，不剔除异常值。', '']
    if development:
        s += [f'开发阶段有 {len(development)} 个参考方法调用因 NumPy 2 的批量 `solve` 右端形状要求而在第 0 次迭代报错；修正为显式单列右端后，以完全相同种子、数据及初值重跑。原始报错文件原样保留，成功重跑文件加 `_retry1`。这是运行接口错误，单独列出，不将其解释为算法收敛失败，也不隐去原记录。','']
    for section,ix in [('支承模态',range(3)),('柔性模态',range(3,8))]:
        s += [f'## {section}参数误差','',
              '单位均为百分比。均值 ± 标准误；标准差与 RMSE 另列。','',
              '| 方法 | 真频率/Hz | 频率误差均值 ± SE | 频率 SD | 频率 RMSE | 阻尼误差均值 ± SE | 阻尼 SD | 阻尼 RMSE |',
              '|---|---:|---:|---:|---:|---:|---:|---:|']
        for method in METHODS:
            m=methods[method]; st=m['stats_percent']
            for j in ix:
                s.append(f"| {m['label']} | {FM[j]:.3f} | {st['mean'][j]:.5f} ± {st['se'][j]:.5f} | {st['std'][j]:.5f} | {st['rmse'][j]:.5f} | {st['mean'][8+j]:.4f} ± {st['se'][8+j]:.4f} | {st['std'][8+j]:.4f} | {st['rmse'][8+j]:.4f} |")
        s += ['']
    s+=['## 对预期的判定','',
        '**支持“主要偏差集中在控制带宽内”，但不支持“三种方法的柔性参数都严格无偏／偏差只存在于带宽内”的强表述。**', '',
        '1. 支承模态：M1 阻尼均值偏差为 +51.13%、+22.54%、+4.86%；M2 降为 +1.39%、+1.05%、+0.815%；M3 为 −1.47%、+0.523%、+0.454%，均处于约 2 个标准误以内。参考 IV 显著消除了这组数据中的大偏差，但有限 20 次结果不是一致性证明。',
        '2. 柔性模态：M3 的 10 个频率/阻尼均值均在 1.1 个标准误内；与近似无偏相容。M2 总体也很接近真值，但 490.39 Hz 阻尼均值 +0.02823% ±0.01388%（约 2.03 SE），不能强行标注全部无偏。',
        '3. M1 的柔性误差工程量级很小，频率最大均值误差 0.00546%，阻尼最大均值误差 0.1577%；但部分已能统计分辨。例如 490.39 Hz 的频率 −0.001465% ±0.000179%（8.18 SE）、阻尼 +0.08576% ±0.01316%（6.52 SE），718.958 Hz 阻尼 −0.15767% ±0.03268%（4.82 SE）。这些不能归入“完全无偏”。',
        '4. 同一数据的配对比较进一步显示：M1−M3 在 490.39 Hz 的阻尼相对误差之差（百分比点）为 +0.07309 ±0.00929，在 718.958 Hz 为 −0.18570 ±0.01920；M2−M3 的相应差减小为 +0.01556 ±0.00488 和 −0.02934 ±0.01372。这里比较的是 100×(估计阻尼/真阻尼−1) 的差，不是阻尼比本身的绝对百分比点。谱线直接法减弱而未严格消除有限周期偏差。',
        '5. 原因解释应限于有依据的层次：闭环灵敏度在控制带宽外并不严格为零，输入扰动的输出等效噪声也随频率变化；加性模型各分块的估计方程在全频带联合求解。因此“带宽外偏差较小”有依据，“带宽外必然零偏差”没有一般保证。本次没有通过额外因素消融将以上小偏差唯一归因于某一项。',
        '', '## 收敛阈值敏感性补查', '',
        '在固定的重复 0 和 19 上，从主实验返回解继续迭代，将极点步长阈值由 1e-5 收紧为 1e-8，最多再做 20 步。6 次补查均达到补查上限，未通过更严格阈值，记录原样保留；不能将其称为严格阈值下已收敛。以参数相对误差之差（百分比点）表示继续迭代引起的最大改变：支承参数 0.001256，柔性频率 6.76e-8，柔性阻尼 8.59e-6；上述阻尼数值均不是阻尼比本身的绝对百分比点。因而在这两个预先选定重复上，小幅柔性均值差远大于继续迭代引起的参数漂移；但这不是所有 20 次运行的严格数值误差证书。原始前向差分步长仍为 1e-6，没有为获得漂亮结果修改求解器。详见 `tolerance_check.json` 和 `tolerance_*.npz`。',
        '', '## 文件与复现', '',
        '- `run_b.py`：运行及原代码等价性核查。原始 cl2 只读导入；原件 SHA256 见 `validation.json`。',
        '- `raw/repNNN_METHOD.npz`：全部参数 beta、初始值、真值、误差、每步历史和 JSON 元数据。对应 `.json` 可直接阅读。',
        '- `summary.json`：每个参数的均值、SE、SD、RMSE、最小/最大值、耗时及全部终止状态。',
        '- `run.log`：首轮日志；`reference_retry.log`：修正批量 solve 接口后的复算日志（若存在）。',
        '', '```powershell', 'python run_b.py --validate', 'python run_b.py --reps 20 --workers 2', 'python summarize_b.py', '```','']
    (HERE/'report_b.md').write_text('\n'.join(s),encoding='utf-8')
    print(json.dumps({m:{k:v for k,v in methods[m].items() if k not in ('records','stats_percent','stats_converged_percent')} for m in METHODS},ensure_ascii=False,indent=2))

if __name__=='__main__': main()
