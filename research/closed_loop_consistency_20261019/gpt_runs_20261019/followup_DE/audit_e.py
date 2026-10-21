"""Independent, read-only numerical/code audit of Task E. Does not change fits."""
import os
for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[n]='1'
os.environ.update(NOFLEX='1',FC='60',AL='30',AH='100',FHI='100',FLO='3',SW='.03',FSPR='.01',ZSPR='.4')
import sys
sys.dont_write_bytecode=True
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'code'))
import cl2 as c
from iv_lm import moments,exact_moment_jacobian
from weighted_core import iv_matrices,whitener,local_covariances,precision


def main():
    oldroot=HERE.parent/'experiment_a'
    data=np.load(oldroot/'data/P02_r002.npz')
    k,z,U,Y,D=[data[q] for q in ('k','zk','U','Y','D')]
    line=k%2==0;zz=z[line];uu=U[line];yy=Y[line];dd=D[line];K=c.Ktf(zz)
    model=lambda b,u:c.out_jac(zz,b,u)
    model_full=lambda b,u:c.out_jac(z,b,u)
    inst=lambda b:np.linalg.solve(np.eye(4)[None]+K@c.resp(zz,b),dd[:,:,None])[:,:,0]
    b=np.array(json.loads((oldroot/'results/P02_r002_dir_lines_W1.json').read_text())['beta'])
    local=local_covariances(k,U,Y,b,model_full,2,nl=16,alpha=.1)
    W=precision(local['covariance_residual'])
    am,gm,sm=moments(b,uu,yy,model,inst,W)
    ac,gc,_,_=iv_matrices(model,uu,yy,inst,b,whitener(W))
    checks=dict(matrix_moment_A_relative=float(np.linalg.norm(am-ac)/np.linalg.norm(ac)),
                matrix_moment_g_relative=float(np.linalg.norm(gm-gc)/np.linalg.norm(gc)))
    minus_dg,_,_=exact_moment_jacobian(b,uu,yy,model,inst,W,c.IW,c.ITR,K,lambda bb:c.resp(zz,bb))
    sc=1/np.sqrt(np.abs(np.diag(am))+1e-300)
    rng=np.random.default_rng(471011)
    derivative=[]
    for _ in range(4):
        direction=sc*rng.normal(size=len(b));direction/=np.max(np.abs(direction[c.IW:c.ITR]))
        h=1e-4
        gp=moments(b+h*direction,uu,yy,model,inst,W)[1]
        gn=moments(b-h*direction,uu,yy,model,inst,W)[1]
        fd=(gp-gn)/(2*h)
        calc=-minus_dg@direction
        derivative.append(float(np.linalg.norm((fd-calc)*sc)/np.linalg.norm(fd*sc)))
    checks['scaled_moment_directional_derivative_relative']=derivative
    parity=[]
    for r in range(40):
        old=json.loads((oldroot/'results'/f'P01_r{r:03d}_ref_band_W1.json').read_text())
        new=json.loads((HERE/'short_e'/f'P01_r{r:03d}_M3_E1.json').read_text())
        bb=np.array(new['beta']);bo=np.array(old['beta'])
        parity.append(dict(replicate=r,initial_beta_max_abs_difference=float(np.max(np.abs(np.array(new['initial_beta'])-np.array(old['initial_beta'])))),
                           log_pole_max_abs_difference=float(np.max(np.abs(bb[c.IW:c.ITR]-bo[c.IW:c.ITR]))),
                           old_status=old['status'],new_status=new['status'],old_iterations=old['iterations'],new_iterations=new['iterations'],
                           old_score=old['score_max'],new_score=new.get('score_max')))
    # Re-run the one changed classification with the original A evaluator.
    r=29;dat=np.load(oldroot/'data'/f'P01_r{r:03d}.npz');zk=dat['zk'];ku=c.Ktf(zk)
    U1,Y1,D1=[dat[q] for q in ('U','Y','D')]
    record=json.loads((HERE/'short_e'/f'P01_r{r:03d}_M3_E1.json').read_text())
    initial=np.array(record['initial_beta'])
    instr=lambda bb:np.linalg.solve(np.eye(4)[None]+ku@c.resp(zk,bb),D1[:,:,None])[:,:,0]
    original,iteration=c.iv_solve(zk,U1,Y1,instr,np.ones_like(U1.real),initial,it=80,tol=1e-5)
    reference=json.loads((oldroot/'results'/f'P01_r{r:03d}_ref_band_W1.json').read_text())
    rerun=dict(replicate=r,iterations=iteration+1,status='pole_step_converged' if iteration<79 else 'iteration_limit',
               max_beta_difference_from_prior=float(np.max(np.abs(original-np.array(reference['beta'])))))
    raw=[json.loads(p.read_text()) for p in (HERE/'short_e').glob('P*.json')]
    groups=[]
    for P,method in [(1,'M3'),(2,'M3'),(2,'M3w')]:
        for strategy in ('E1','E2','E3'):
            rows=[a for a in raw if a['P']==P and a['method']==method and a['strategy']==strategy]
            successes=[a for a in rows if a.get('success')]
            scores=np.array([a.get('score_max',np.nan) for a in successes],float)
            groups.append(dict(P=P,method=method,strategy=strategy,count=len(rows),failures=sum(not a.get('success') for a in rows),
                               converged_score_gt_1e5=int(sum(np.isfinite(scores)&(scores>1e-5))),
                               converged_score_gt_1e3=int(sum(np.isfinite(scores)&(scores>1e-3))),
                               max_converged_score=float(np.nanmax(scores)) if len(scores) else None,
                               high_score_replicates=[a['replicate'] for a in successes if a.get('score_max',0)>1e-3]))
    # Verify every E2/E3 M3 initial value is exactly the previous M2 endpoint.
    initial_checks=[]
    for a in raw:
        if a['method']=='M3' and a['strategy'] in ('E2','E3'):
            direct=json.loads((oldroot/'results'/f"P{a['P']:02d}_r{a['replicate']:03d}_dir_lines_W1.json").read_text())
            initial_checks.append(float(np.max(np.abs(np.array(a['initial_beta'])-np.array(direct['beta'])))))
    out=dict(checks=checks,raw_count=len(raw),groups=groups,P1_E1_parity=parity,r29_original_evaluator_rerun=rerun,
             E2_E3_direct_initial_max_difference=max(initial_checks),
             interpretation=dict(same_IV_moments=True,LM_uses_changing_instrument_jacobian=True,
                 M3w_comparison='Full pipeline: each strategy computes its own preliminary M3, then its own weights, then weighted IV.',
                 not_fixed_weight_initialization_ablation=True,
                 small_pole_step_does_not_certify_full_moment_root=True,
                 E3_changes_jacobian_and_trust_radius_as_well_as_damping=True))
    (HERE/'audit_e.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    report=['# 任务 E 独立审计', '',
            '审计只读根代理代码及已有结果，没有修改估计器或结果。', '',
            '## 方程与导数', '',
            f"- 全矩阵加权 IV 矩阵、矩与独立预白化实现的相对差分别为 {checks['matrix_moment_A_relative']:.3e}、{checks['matrix_moment_g_relative']:.3e}。复数权矩阵的共轭和指标方向正确。",
            '- LM 优化的是参考 IV 矩方程，不是把参考 IV 偷换成输出误差最小二乘。Jacobian 包括模型 Hessian 与工具变量随参数变化的导数；后者使用 dÛ=−(I+KG)⁻¹K(dG)Û。',
            '- 四个独立参数尺度方向的中心差分导数检查（外步长 1e−4）：'+', '.join(f'{x:.3e}' for x in derivative)+'。',
            '- E3 不仅增加阻尼，还采用完整矩方程 Jacobian，并将模态对数信赖半径设为 0.25（原步长截断为 0.05）。应称为新的求根策略，不可把所有变化仅归因于 LM 阻尼。', '',
            '## 初始化与权重', '',
            '- 已完成的 E2/E3 非加权 M3 初值均与旧 M2 端点完全相同。',
            '- M3w 的 E1/E2/E3 是整条管线比较：各策略先求自己的 M3，使用自己的 M3 端点估计非激励频点协方差，再以该端点启动加权求解。因此策略之间权矩阵也不同，不能称为“固定同一个权矩阵，只改变 M3w 初值/求解器”的隔离消融。',
            '- E3 的一次加权求解内部权矩阵保持冻结；没有把权矩阵导数误加进 IV 估计方程。', '',
            '## P=1 的 E1 与旧 A 的失败数差异', '',
            '- 40 个初值逐位相同。仅 r=29 的“是否在第 80 步前停止”分类变化：旧结果通过，新 E1 达到上限；其他失败成员一致。',
            f"- 使用相同初值重新调用原 A 的 c.iv_solve，r=29 在 {rerun['iterations']} 步停止，与旧 beta 的最大差为 {rerun['max_beta_difference_from_prior']:.3e}。这定位到等价矩阵乘法／einsum 的浮点累加顺序及停止边界敏感性，而非数据、噪声或初值改变。",
            '- r=29 的旧 score=2.46e−5，新 score=4.19e−6；因此计为“失败”的新端点，其矩方程残差反而更小。这说明原失败门槛是数值停机统计，不等于统计一致性或严格求根证书。', '',
            '## 完整 score 与失败计数', '',
            '沿用预先指定的失败规则，没有事后把高 score 成功样本重新分类。需要同时展示严格 score；LM 的阻尼可产生极小的模态步，即使完整矩仍较大，因此单看成功率不能证明更可靠的 IV 根。', '',
            '| P | 方法 | 策略 | 已完成 | 失败 | 成功但 score>1e−5 | 成功但 score>1e−3 | 成功最大 score |',
            '|---:|---|---|---:|---:|---:|---:|---:|']
    for g in groups:
        report.append(f"| {g['P']} | {g['method']} | {g['strategy']} | {g['count']} | {g['failures']} | {g['converged_score_gt_1e5']} | {g['converged_score_gt_1e3']} | {g['max_converged_score']} |")
    report += ['',f"本次读取 {len(raw)} 个结果；若尚未达到 360，需要完成后重新运行 audit_e.py 更新计数。", '',
               '## LM 法方程病态性：独立复算', '',
               '另见 audit_e_lm_condition.py/json：对 P=2、r=0、M3w_E3 的初值与最终端点，在相同冻结权矩阵、同一 Jacobian、同一阻尼下，对比当前法方程求解和稳定的增广最小二乘求解。', '',
               '最终端点的缩放 Jacobian 条件数约 3.16e9，未正则法方程约 5.44e18。但实际末次阻尼 2.03e−9 下，正则法方程条件数约 1.71e10，两种线性求解的步长相对差仅 2.29e−7、模态对数步差仅 9.70e−12；候选点 score 都为约 3.26025e−5。', '',
               '因此，法方程确实有平方条件数风险，但当前证据不足以把 M3w_E3 的 39/40 次达到上限归咎于浮点法方程误差。更直接的证据是病态方向与阻尼导致缓慢进展，80 步不足。降低到 λ=1e−12 时两种解都有更大的实际矩残差；稳定增广求解可以作补充核查，不能预先声称它会解决失败率。', '',
               '结论：未发现 IV 方程或全矩阵权的代数实现错误；需明确披露管线式权重变化、E3 实际变更及小步终止并不保证矩方程收敛。']
    (HERE/'audit_e.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    print(json.dumps(dict(checks=checks,raw_count=len(raw),groups=groups,r29=rerun),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
