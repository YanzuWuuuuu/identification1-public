"""Render the already completed deterministic controller audit."""
from pathlib import Path
import json
import hashlib

ROOT = Path(__file__).resolve().parent
d = json.loads((ROOT/'controller_check.json').read_text(encoding='utf8'))
base = d['original_and_rolloff'][0]
rows = ['| 控制器 FC / 新低通截止 (Hz) | 高频最大奇异值 / 峰位 (Hz) | 高频最大谱半径 | 闭环最大极点模 | 稳定 |',
        '|---|---:|---:|---:|---|']
selected = [r for r in d['original_and_rolloff'] if r['added_rolloff_fc_hz'] in (None, 150., 100., 60., 40., 35.)]
selected += [d['original_and_rolloff'][-1]]
selected += [r for r in d['lower_controller_fc_alternatives'] if
    (r['controller_fc_hz']==15. and r['added_rolloff_fc_hz'] is None) or
    (r['controller_fc_hz']==20. and r['added_rolloff_fc_hz']==100.)]
for r in selected:
    cut = r['added_rolloff_fc_hz']
    label = '无新增低通' if cut is None else f'{cut:.5g}'
    rows.append(f"| {r['controller_fc_hz']:g} / {label} | {r['high_max_sigma']['value']:.6f} / {r['high_max_sigma']['frequency_hz']:.4f} | {r['high_max_rho']['value']:.6f} | {r['closed_loop']['radius']:.9f} | {'是' if r['closed_loop']['stable'] else '否'} |")
details = ['| 配置 | 10 Hz: 幅值 dB / 相位° | 30 Hz | 60 Hz | 100 Hz |', '|---|---:|---:|---:|---:|']
for label, r in [('FC60 + 150 Hz', next(r for r in selected if r['added_rolloff_fc_hz']==150.)),
                  ('FC60 + 38.93550 Hz', d['original_and_rolloff'][-1]),
                  ('FC15、无新增低通', selected[-2])]:
    pp = {x['frequency_hz']:x for x in r['points']}
    cells = [f"{pp[f]['total_controller_change_from_fc60_db']:.3f} / {pp[f]['total_controller_change_from_fc60_deg']:.3f}" for f in (10.,30.,60.,100.)]
    details.append('| '+label+' | '+' | '.join(cells)+' |')

text = '''# 控制器滚降的无噪声预检查

**结论：原 FC=60 控制器后直接串联约 150 Hz 的二阶低通不可行。它既不满足高频小回路条件，又使闭环失稳。** 在本次检查的二阶 Butterworth 家族中，没有找到同时满足稳定与 100–550 Hz 最大奇异值小于 0.1 的截止频率；更不能保持原 3–100 Hz 回路不变。

本文件只包含确定性的频率求值、状态矩阵特征值和设计边界检查。没有进行 Monte Carlo、参数拟合或时域仿真，也没有启动正式批次。代码及依赖哈希均保存在本目录。原始代码及已有结果只读。

## 1. 模型与数值定义

- 真对象保留 10.6、11.7、16.1、286.573、344.945、490.39 Hz 六个模态，删除 718.958、755.501 Hz；其余模态留数与阻尼完全沿用 cl2。
- 每个模态用精确 ZOH 离散化，采样率 5000 Hz；输入整数延迟 NI=[0,1,0,1]，输出 NO=[1,0,1,0]，保持原七维结构化时延真值。
- 原控制器 FC=60，已有 FLP=120、阻尼 0.5 的低通保持不变。“新增二阶低通”另取数字 Butterworth，按物理截止频率和 fs=5000 构造，四个输入通道共用标量响应。
- 主检验是 max σmax(GK)，另外列出 max ρ(GK)，两者不混用。用 0.02 Hz 网格覆盖 3–550 Hz，并对高频峰值作确定性的标量极值细化。
- 闭环稳定性由含全部 ZOH 模态、输入/输出整数延迟、控制器和新增滤波器状态的闭环矩阵直接求特征值，条件为全部极点模严格小于 1。
- FC 是原代码中的控制器设计参数，不等同于经过本模型精确求出的闭环带宽。

## 2. 主要结果

'''+ '\n'.join(rows)+'''

表中 38.93550 Hz 是高频最大奇异值恰好等于 0.1 的数值边界，不作为满足严格不等式的可执行配置。即使将截止进一步降低到 35 Hz 以留出幅值裕量，闭环仍失稳。

对新增截止 5–500 Hz、步长 1 Hz 的 496 个确定性配置检查，闭环最大极点模均大于 1，最小值为 1.004146656。这个有限网格不冒充对所有连续截止的解析稳定性证明，但候选 150 Hz、临界 38.93550 Hz 及其低截止邻域已有直接数值否定。没有继续自动调参寻找其他控制器结构。

## 3. 为什么“压低高频而保持 3–100 Hz 不变”在此不相容

六模态原对象在 100 Hz 处 σmax(GK)=0.668638667。新增标量低通 L 满足

\[
\sigma_{\max}(G K L)=|L|\,\sigma_{\max}(G K).
\]

所以要让包含 100 Hz 端点的高频区间满足 <0.1，必须有 |L(100)|<0.14955766，即至少降低 16.50383 dB。100 Hz 同时是要保留的低频区间端点；此处不可能同时“基本不变”和下降至少 16.5 dB。该矛盾不限于 Butterworth 截止选择。

下表给出相对原 FC=60 控制器的幅相变化；未自行给“明显改变”设定判据。

'''+ '\n'.join(details)+'''

150 Hz 新低通在 60 Hz 只衰减 0.109 dB，但新增约 33.87° 相位滞后，已经足以使当前低稳定裕量的原闭环失稳。因此只看低频幅值近似不变会误判。

## 4. 可供讨论的最小替代

若允许真实改变控制器低频回路，则采用原控制器公式、改 FC=15，保留原 FLP=120、不增加新低通，可得到：高频 max σ=0.076986965，闭环最大极点模 0.995699952。FC=20 配合新增 100 Hz 低通也可行：max σ=0.091101910，闭环最大极点模 0.995790126。

建议优先讨论 FC=15 的简单备选，因为它不引入新的控制器结构。但它显著改变 3–100 Hz 回路，不能称为“原回路保持不变”的实现；若采用，必须作为待用户确认的实验设计修正。这里没有据此生成数据或启动 MC。

原 FC=60 可以保留为高回路配置，用于定量偏差预测；FC=15 是满足高频小回路前提的候选。两者应分别使用各自真实回路和信息矩阵计算预测，不能沿用原控制器的预测值。

## 5. 独立实现核对与局限

'''
v=d['validation']
text += f"- 六模态直接状态空间频响与模态叠加频响相对差：{v['plant_state_vs_modal_relative_error']:.3e}。\n"
text += f"- 新控制器级联状态空间与标量频响相乘的最大相对差：{max(r['relative_error'] for r in v['cascade_state_vs_scalar_product']):.3e}。\n"
text += f"- 恢复旧八模态配置得到闭环极点最大模 {v['eight_mode_original_radius']:.12f}、高频 max σ={v['eight_mode_original_sigma']:.12f}，与旧审计一致。\n"
text += '- 本检查没有把带外模态近似为 D0+D2；它检查的是指定六模态真对象的精确闭环。因此后续若额外改变真实背景、模态或时延，需重新检查，不能直接沿用这里的稳定性结论。\n'
text += '- 哈希前后完全一致：'+str(d['source_hashes']==d['source_hashes_after'])+'。\n'
(ROOT/'控制器无噪声预检查.md').write_text(text, encoding='utf8')
(ROOT/'artifact_hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest()
    for p in ROOT.iterdir() if p.is_file() and p.name!='artifact_hashes.json'}, indent=2), encoding='utf8')
