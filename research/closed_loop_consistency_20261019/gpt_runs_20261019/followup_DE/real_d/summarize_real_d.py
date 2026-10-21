"""Report all real-data endpoints and independently estimated halves, without truth claims."""
import os
for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[n]='1'
import sys
sys.dont_write_bytecode=True
import json
import hashlib
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
METHODS=('M2','M3','M3w','M5','M5n')
ZPOSITIONS=('-0.5','0','0.5','1')


def key(z,p,start,method): return f'z{z}_P{p:02d}_s{start:02d}_{method}'
def read(z,p,start,method):
    path=HERE/'raw'/(key(z,p,start,method)+'.json')
    if not path.exists(): return {'status':'missing','failure':True,'missing':True}
    return json.loads(path.read_text(encoding='utf-8'))
def trio(x,digits=4):
    if x is None:return '—'
    return '/'.join(f'{a:.{digits}f}' if np.isfinite(a) else 'NaN' for a in x)
def score(x):
    if x is None:return '—'
    return f'{x:.2e}'


def main():
    requests=json.loads((HERE/'window_requests.json').read_text(encoding='utf-8'))
    raw={p.stem:json.loads(p.read_text(encoding='utf-8')) for p in (HERE/'raw').glob('*.json')}
    rows=[]
    for q in requests:
        for meth in METHODS:
            result=read(q['z'],q['nper'],q['start'],meth)
            rows.append(dict(**q,method=meth,raw_key=key(q['z'],q['nper'],q['start'],meth),result=result))
    halves=[]
    for z in ZPOSITIONS:
        for p in (2,4,8):
            for meth in METHODS:
                a,b=read(z,p//2,1,meth),read(z,p//2,1+p//2,meth)
                h=dict(z=z,P=p,method=meth,front_status=a.get('status'),back_status=b.get('status'),
                       front_failure=a.get('failure'),back_failure=b.get('failure'))
                if all('frequency_hz' in x for x in (a,b)):
                    df=np.asarray(b['frequency_hz'])-np.asarray(a['frequency_hz'])
                    dz=np.asarray(b['damping_percent'])-np.asarray(a['damping_percent'])
                    h.update(frequency_difference_hz=df.tolist(),damping_difference_percentage_points=dz.tolist(),
                             frequency_half_range_hz=(np.abs(df)/2).tolist(),
                             damping_half_range_percentage_points=(np.abs(dz)/2).tolist(),
                             both_converged=not a.get('failure') and not b.get('failure'))
                else:
                    h['unavailable']=True
                halves.append(h)
    counts={}
    for meth in METHODS:
        rr=[r for r in raw.values() if r.get('method')==meth]
        applicable=[r for r in rr if not r.get('unavailable')]
        full=[read(z,p,1,meth) for z in ZPOSITIONS for p in (2,4,8)]
        counts[meth]=dict(unique_records=len(rr),applicable=len(applicable),unavailable=len(rr)-len(applicable),
                          failures=sum(bool(r.get('failure')) for r in applicable),
                          full_window_failures=sum(bool(r.get('failure')) for r in full),
                          score_over_1e5=sum(np.isfinite(r.get('score_normalized',np.nan)) and r['score_normalized']>1e-5
                                            for r in applicable if r.get('score_normalized') is not None))
    integrity=[]
    for src in json.loads((HERE/'source_inventory.json').read_text(encoding='utf-8')):
        path=Path(src['path'])
        now=hashlib.sha256(path.read_bytes()).hexdigest()
        integrity.append(dict(path=str(path),unchanged=now==src['sha256'],sha256_before=src['sha256'],sha256_after=now))
    parity=[]
    old=HERE.parents[1]/'experiment_c'/'raw'
    for z in ZPOSITIONS:
        for p in (2,4,8):
            prior=json.loads((old/f"z{z}_P{p}_s1_zoh_{'M5_reg' if p<5 else 'M5'}.json").read_text(encoding='utf-8'))
            new=read(z,p,1,'M5')
            parity.append(dict(z=z,P=p,max_beta_difference=float(np.max(np.abs(np.asarray(prior['beta'])-np.asarray(new['beta'])))),
                               same_nfev=prior['nfev']==new['iterations']))
    agreement=[]
    for p in (2,4,8):
        dd_f=[];dd_z=[];included=[]
        for z in ZPOSITIONS:
            a,b=read(z,p,1,'M3w'),read(z,p,1,'M5n')
            if not a['failure'] and not b['failure']:
                dd_f.append(np.abs(np.asarray(a['frequency_hz'])-np.asarray(b['frequency_hz'])))
                dd_z.append(np.abs(np.asarray(a['damping_percent'])-np.asarray(b['damping_percent'])))
                included.append(z)
        agreement.append(dict(P=p,converged_positions=included,
                              max_frequency_difference_hz=np.max(dd_f,axis=0).tolist(),
                              max_damping_difference_percentage_points=np.max(dd_z,axis=0).tolist()))
    summary=dict(expected_context_rows=180,context_rows=len(rows),expected_unique_raw=140,unique_raw=len(raw),
                 counts=counts,full_and_half_contexts=rows,halves=halves,source_integrity=integrity,
                 M5_exact_parity_with_previous_C=parity,M3w_M5n_converged_agreement=agreement)
    (HERE/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# 实测 D：谱线参考 IV、非参数权与 SML', '',
           '四个位置均使用原始记录第 2 个周期起的数据；频带 3–25 Hz。模型保留前次 C 的 ZOH 三个支承模态、D0、D2、固定通道时延（统一减 100 µs）和逐模态瞬态，没有换模型。', '',
           '所有方法只使用激励谱线。M3w/M5n 的非激励频点协方差取最近 16 个总频点、收缩 α=0.1、只更新一次。M3w 收缩残差协方差，M5n 收缩联合协方差；这两种正则化映射不同。M3 的失败有限端点也会用于构造权重，并单独标记，不隐去失败。', '',
           'IV 采用 A 的 80 次上限、模态对数步阈值 1e−5、步长上限 0.05。它与前次实测 C 的 1e−6/0.1 略有不同，这是按本次统一失败规则做的调整。SML 保持前次实测最多 160 次函数求值；表内 SML 的“迭代”列为 nfev，不等同于 IV 步数。严格归一化 score 另外报告，不重新定义失败。', '',
           'M5 使用周期 DFT 均值及均值协方差（P(P−1) 分母），其初值在均值单位重新求线性系数。M5n 使用全记录 DFT 及同单位非激励频点协方差，初值也在全记录单位求解；不存在把全长协方差直接配给周期均值的 P² 尺度错误。', '',
           '三元数依次对应三个支承模态。阻尼单位为 %。没有真值，以下仅比较方法一致性、收敛和半段重复性，不判断谁更准。', '',
           '结果中应特别保留两项不利证据：M3w 在 P=2、z=−0.5/+0.5 两个完整窗口仍失败；原始 M5 在 P=8 的前两个模态与其他方法明显不一致，z=+0.5 还达到求解上限。M5 的 12 个完整窗口与前次 C 的全部 beta 逐位相同，因此该不一致不是本轮 DFT 尺度或模型改动造成的。它不能据此被判为“不准确”，也不能因优化器终止成功就视作可靠物理解。', '',
           'P=8 下 M3w/M5n 的四个位置都收敛且彼此较接近，但半段差异仍可见，例如 z=0 的 M3w 第一阶阻尼前后半段相差 −0.864 个百分点（全窗口估计约 4.754%）。不能把较小的频率差当作阻尼也已稳定的证据。', '',
           '## 完整窗口：4 个位置 × 3 个 P × 5 个方法', '',
           '| z/mm | P | 方法 | f/Hz（三模态） | ζ/%（三模态） | 状态 | 迭代或 nfev | score |',
           '|---:|---:|---|---|---|---|---:|---:|']
    for z in ZPOSITIONS:
        for p in (2,4,8):
            for meth in METHODS:
                r=read(z,p,1,meth)
                lines.append(f"| {z} | {p} | {meth} | {trio(r.get('frequency_hz'),5)} | {trio(r.get('damping_percent'),4)} | {r.get('status')} | {r.get('iterations','—')} | {score(r.get('score_normalized'))} |")
    lines += ['', '## 前后两半：分别重估', '',
              '每个 P 周期窗口拆成两个 P/2 周期窗口，各自重新拟合控制器、模型和噪声权。差值为后半减前半；|差值|/2 另存为半段散布指标，不称为标准差或置信区间。若某一半不收敛，差值仍保留并标明；P=2 的两半各只有一个周期，M3w/M5/M5n 不适用。', '',
              '| z/mm | 总 P | 方法 | Δf/Hz（三模态） | Δζ/百分点（三模态） | 前半状态 | 后半状态 |',
              '|---:|---:|---|---|---|---|---|']
    for h in halves:
        lines.append(f"| {h['z']} | {h['P']} | {h['method']} | {trio(h.get('frequency_difference_hz'),5)} | {trio(h.get('damping_difference_percentage_points'),4)} | {h['front_status']} | {h['back_status']} |")
    lines += ['', '## 失败与不适用计数（去重后的独立窗口）', '',
              '| 方法 | 完整窗口失败/12 | 所有适用窗口失败/适用数 | 不适用 | score > 1e−5 |',
              '|---|---:|---:|---:|---:|']
    for meth,c in counts.items():
        lines.append(f"| {meth} | {c['full_window_failures']}/12 | {c['failures']}/{c['applicable']} | {c['unavailable']} | {c['score_over_1e5']} |")
    lines += ['', '## 可复查性', '',
              f'36 个报告窗口上下文对应 28 个独立窗口；重复的前半与完整窗口共用完全相同的结果，避免把它们误计为独立实验。保存 {len(raw)} 个原始方法记录，报告展开为 {len(rows)} 行。P=1 半窗口的 24 个协方差方法记录标为不适用；其余 116 个方法记录全部保留。', '',
              '原始数据和前次模型文件的 SHA-256 校验：'+('全部未改动。' if all(r['unchanged'] for r in integrity) else '发现不一致，查看 summary.json。'), '',
              '逐窗口 beta、初值、原始求解器信息、频率、阻尼、score 和失败端点见 raw；FFT、周期 FFT、控制器和初值见 data；非参数协方差和邻域索引见 raw/*_cov.npz。完整 180 行上下文和 60 条半段差值见 summary.json。']
    (HERE/'实测D_结果与半段差异.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(unique_raw=len(raw),context_rows=len(rows),counts=counts,unchanged=all(r['unchanged'] for r in integrity)),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
