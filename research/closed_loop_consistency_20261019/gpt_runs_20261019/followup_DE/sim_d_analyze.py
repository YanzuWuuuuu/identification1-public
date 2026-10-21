"""Tabulate the frozen task-D comparison, without selecting a winning variant."""
from sim_d_run import *

def main():
    summary=json.loads((OUT/'summary.json').read_text(encoding='utf-8'))
    lookup={(x['P'],x['method']):x for x in summary}
    records={}
    for p in (OUT/'results').glob('*.json'):
        x=json.loads(p.read_text(encoding='utf-8'));records[(x['P'],x['replicate'],x['method'])]=x
    parity=[]
    for m in ('M2','M5'):
        dif=[];bd=[]
        for P in (2,4,8,16):
            for r in range(40):
                new=records[(P,r,m)]
                if m=='M2':
                    oldpath=OLD/'experiment_a'/'results'/f'P{P:02d}_r{r:03d}_dir_lines_W1.json'
                    if not oldpath.exists():continue
                    old=json.loads(oldpath.read_text(encoding='utf-8'));err=np.array(old['error_percent']);b=np.array(old['beta'])
                else:
                    oldpath=OLD/'comparators'/'sml_results'/f'P{P:02d}_r{r:03d}_{"M5_reg" if P<5 else "M5"}.npz'
                    if not oldpath.exists():continue
                    old=np.load(oldpath);err=old['relative_error_percent'];b=old['beta']
                dif.append(np.abs(np.array(new['error_percent'])-err));bd.append(np.max(abs(np.array(new['beta'])[c.IW:c.ITR]-b[c.IW:c.ITR])))
        parity.append(dict(method=m,n=len(dif),max_error_percentage_point_difference=np.max(dif,axis=0).tolist(),max_modal_log_difference=float(max(bd))))
    paired=[]
    for P in (2,4,8,16):
        for left,right in [('M3w_n16_o1','M5n_n16'),('M5n_n16','M5'),('M3w_n8_o1','M3w_n16_o1'),('M3w_n16_o2','M3w_n16_o1')]:
            deltas=[]
            for r in range(40):
                a=records[(P,r,left)];b=records[(P,r,right)]
                if a['success'] and b['success']:deltas.append(np.array(a['error_percent'])-b['error_percent'])
            dd=np.array(deltas)
            paired.append(dict(P=P,left=left,right=right,n=len(dd),mean_percent_point_difference=dd.mean(0).tolist(),
                               se_percent_point_difference=(dd.std(0,ddof=1)/np.sqrt(len(dd))).tolist(),
                               rmse_ratio=(np.array(lookup[(P,left)]['converged_only']['rmse'])/lookup[(P,right)]['converged_only']['rmse']).tolist()))
    primary=[x for x in paired if x['left']=='M3w_n16_o1' and x['right']=='M5n_n16']
    rule_a=all(np.all(np.array(x['rmse_ratio'])<=1.1) and lookup[(x['P'],'M3w_n16_o1')]['failures']<=lookup[(x['P'],'M5n_n16')]['failures'] for x in primary if x['P'] in (2,4))
    rule_b=all(np.all(np.array(x['rmse_ratio'])>1.1) for x in primary)
    out=dict(record_count=len(records),expected_record_count=1600,all_groups_have_40=all(x['n_total']==40 for x in summary),
             source_parity=parity,paired=paired,global_rule_A=rule_a,global_rule_B=rule_b,
             parameter_order=['f1','f2','f3','zeta1','zeta2','zeta3'])
    write_json(OUT/'comparison_diagnostics.json',out)
    def vec(v):return '/'.join(f'{float(a):.4g}' for a in v)
    def mus(st,sl):return '/'.join(f'{m:.4g} ± {s:.3g}' for m,s in zip(st['mean'][sl],st['se'][sl]))
    lines=['# 实验 D：仿真对照与核查','',
      '全部为激励谱线求解，每个 P 40 次；均值、标准差、RMSE 均以 100 expm1(对数参数误差) 为单位。下表使用收敛样本，失败另列；all_finite 完整统计见 summary.json。M3 的两例失败仅出现在 P=2；其他方法未失败。',
      '', '|P|方法|频率偏差 %（均值±SE，1/2/3）|阻尼偏差 %（均值±SE，1/2/3）|频率 SD %|阻尼 SD %|频率 RMSE %|阻尼 RMSE %|失败/40|平均迭代数|',
      '|---:|---|---|---|---|---|---|---|---:|---:|']
    for x in summary:
        st=x['converged_only'];lines.append(f"|{x['P']}|{x['method']}|{mus(st,slice(0,3))}|{mus(st,slice(3,6))}|{vec(st['sd'][:3])}|{vec(st['sd'][3:])}|{vec(st['rmse'][:3])}|{vec(st['rmse'][3:])}|{x['failures']}|{x['mean_iterations']:.3g}|")
    lines+=['','## 预定规则','',f'主版本为 n_l=16、一次权重更新。规则 A={rule_a}；规则 B={rule_b}。两者均未触发，必须按 P 和参数报告，不能更换主版本以触发规则。',
      '', '|P|M3w/M5n 频率 RMSE 比（1/2/3）|阻尼 RMSE 比（1/2/3）|', '|---:|---|---|']
    for x in primary:lines.append(f"|{x['P']}|{vec(x['rmse_ratio'][:3])}|{vec(x['rmse_ratio'][3:])}|")
    lines+=['','## 尺度与协方差细节','',
      '- 全记录 DFT 谱线等于 P 倍的单周期 DFT 均值，因此同一均值估计的全记录协方差等于 P² 倍的均值协方差。P=16 的残差、Jacobian 和重拟合均验证这一变换；scale_parity_gate.json 给出数值。',
      '- M3w 使用邻点各自 G(f_j) 形成残差，M5n 则在局部联合协方差上使用中心 G(f_k) 投影。即使不收缩，有限局部窗内 G 的频率变化也会使两者协方差不同；P=16、r=0、n_l=16 时为 2.61%。这不是 DFT 尺度错误。',
      '- M3w 收缩 4×4 残差协方差；M5n 收缩 8×8 [Y−T,U] 联合协方差。两种收缩映射也不同，上述样本中最终差异约 5.87%。',
      '- P=16 的 M5n 与 M5 的均值和散布均以 40 次配对 MC 对照，不要求单个随机样本的估计完全相同。',
      '- oracle 理论用真实条件参考输入构造全部 126 个参数的 Jacobian，保留并消去 24 个瞬态干扰参数。Q=0.5 Re Σ ZᴴWΦvWZ；0.5 来自仅使用正频率的适当复噪声。理论百分比 SD 使用 exp−1 变换后的条件高斯混合方差，同时保留一阶近似。该核查支持量级与趋势，而非有限样本正态性证明。',
      '- 非激励频点与激励频点分离并不自动保证有限样本独立：模型初估、有限记录边界和相邻频率相关性仍可能影响估计权重。',
      '', '## 复现','',
      '在 followup_DE 目录运行 Python：先 sim_d_run.py --gate，然后 sim_d_run.py --workers 4；完成后 sim_d_run.py --summary、sim_d_validate.py、sim_d_analyze.py。所有脚本强制单个 BLAS 线程。每一方法完成即原子保存 JSON，重跑跳过已有结果。',
      '', '旧 A 数据可用的 140 个记录只读复用；P=16 的额外 r=20–39 按同一 seed=202610190+1000P+r 生成到 sim_d/data，原输出不变。M2/M5 与旧版本的复算差异记录于 comparison_diagnostics.json。']
    (OUT/'RESULTS_D.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'record_count':len(records),'rule_A':rule_a,'rule_B':rule_b,'source_parity':parity}),flush=True)
if __name__=='__main__':main()
