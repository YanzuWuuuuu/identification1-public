"""Independent finite-difference and state-space validation; no parameter fits."""
from sim_core import *
import hashlib

def relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))

def main():
    dest=ROOT/"audit"
    dest.mkdir(exist_ok=True)
    m=SimModel();truth=m.truth()
    f=np.unique(np.r_[np.geomspace(3.,850.,91),(cl.FM*np.array([.997,1.003])[:,None]).ravel()])
    c=m.full.context(f);z=np.exp(2j*np.pi*f/FS)
    G=m.resp(c,truth)
    original=cl.resp(z,cl.truth_beta())
    explicit=original*np.exp(-2j*np.pi*f[:,None,None]*TRUTH_TAU[None])
    result=dict(
        seed=712405,number_of_modes=cl.NM,NOFLEX=os.environ.get("NOFLEX"),
        frequency_hz=cl.FM,damping_ratio=cl.ZM,parameters=m.nb,
        input_delay_samples=NI,output_delay_samples=NO,truth_delay_us=TRUTH_TAU*1e6,
        reconstructed_delay_us=m.delay(truth)*1e6,
        truth_delay_max_absolute_error_s=float(abs(m.delay(truth)-TRUTH_TAU).max()),
        response_relative_error=relative(G,explicit),
        response_max_absolute_error=float(abs(G-explicit).max()),
        source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                      (Path(__file__),Path(__file__).with_name("sim_core.py"),RESEARCH/"code/cl2.py")},
        fitted_parameters=False)
    assert cl.NM==8 and os.environ["NOFLEX"]=="0"
    assert result["truth_delay_max_absolute_error_s"]<1e-15
    assert result["response_relative_error"]<1e-10

    # Independent delayed closed-loop state recursion versus explicit transfer algebra.
    AA,BB,CC,DD=plant_with_delays()
    states=np.linalg.solve(z[:,None,None]*np.eye(len(AA))-AA[None],
                           np.broadcast_to(BB,(len(z),)+BB.shape))
    transfer=CC[None]@states+DD[None]
    K=ktf(f)
    Tu=np.linalg.solve(np.eye(4)[None]+K@G,np.broadcast_to(np.eye(4),(len(f),4,4)))
    Ty=G@Tu
    result.update(closed_loop_radius=float(max(abs(la.eigvals(AA)))),
                  state_y_from_d_relative_error=relative(transfer[:,:4,:4],Ty),
                  state_u_from_d_relative_error=relative(transfer[:,4:,:4],Tu))
    assert result["state_y_from_d_relative_error"]<1e-10
    assert result["state_u_from_d_relative_error"]<1e-10

    # Populate every free linear block, including C1, D0/D2 and both transient terms,
    # so their pole/delay couplings are exercised rather than differentiated at zero.
    rng=np.random.default_rng(712405);b=truth.copy();groups={k:[] for k in ("C0","C1","D0D2","transient","poles","delays")}
    scales=np.ones(m.nb)
    for o in range(4):
        for j in range(8):
            st=o*m.full.L+10*j
            amp=max(float(np.max(abs(cl.RES[j]))),1e-12)
            groups["C0"].extend(range(st,st+4));groups["C1"].extend(range(st+4,st+8))
            groups["transient"].extend((st+8,st+9))
            b[st+4:st+8]=.15*amp*rng.normal(size=4)
            b[st+8:st+10]=1e-10*rng.normal(size=2)
            scales[st:st+8]=amp
            scales[st+8:st+10]=1e-10
        st=o*m.full.L+80
        groups["D0D2"].extend(range(st,st+8))
        b[st:st+8]=1e-8*rng.normal(size=8);scales[st:st+8]=1e-8
    groups["poles"]=list(range(m.nl,m.nl+16))
    groups["delays"]=list(range(m.nl+16,m.nb))
    U=rng.normal(size=(len(f),4))+1j*rng.normal(size=(len(f),4))
    y,J=m.output(c,b,U)
    hs=(1e-4,3e-5,1e-5,3e-6,1e-6)
    checks=[]
    for name,indices in [*groups.items(),("all",list(range(m.nb)))]:
        d=np.zeros(m.nb);d[indices]=rng.normal(size=len(indices))*scales[indices]
        dy=np.einsum("kop,p->ko",J,d)
        errors=[]
        for h in hs:
            num=(m.output(c,b+h*d,U,False)[0]-m.output(c,b-h*d,U,False)[0])/(2*h)
            errors.append(relative(num,dy))
        # Store the full step-size sweep; require the preselected h=1e-5 to pass.
        checks.append(dict(group=name,directional_relative_errors=errors,step_sizes=hs,
                           main_step=1e-5,main_error=errors[2],passed=bool(errors[2]<2e-5),
                           analytical_derivative_norm=float(np.linalg.norm(dy))))
    # Individual structural delay coordinates, not only their combined direction.
    for p in groups["delays"]:
        d=np.zeros(m.nb);d[p]=1
        h=1e-5
        num=(m.output(c,b+h*d,U,False)[0]-m.output(c,b-h*d,U,False)[0])/(2*h)
        err=relative(num,J[:,:,p])
        checks.append(dict(group=f"delay_coordinate_{p-(m.nl+16)}",main_step=h,main_error=err,passed=bool(err<2e-5)))
    result["jacobian_checks"]=checks
    result["jacobian_all_passed"]=all(v["passed"] for v in checks)
    assert result["jacobian_all_passed"]

    # Exact noiseless recursion, identical two-period settling to campaign generation.
    # No refit of poles, delays, residues, background or transient coefficients.
    a=generate(2,712405,noise=False)
    cy=m.full.context(a["f"]);y0=m.output(cy,truth,a["U"],False)[0]
    err=a["Y"]-y0
    isline=a["k"]%a["P"]==0
    result["noiseless"]=dict(
        periods=2,settle_periods=2,seed=712405,
        total_relative_residual=float(np.linalg.norm(err)/np.linalg.norm(a["Y"])),
        excited_line_relative_residual=float(np.linalg.norm(err[isline])/np.linalg.norm(a["Y"][isline])),
        nonexcited_error_relative_to_total_output=float(np.linalg.norm(err[~isline])/np.linalg.norm(a["Y"])),
        maximum_absolute_spectral_residual=float(abs(err).max()),
        transient_coefficients="truth coefficients are zero; no nuisance fit",
        delay_register_boundary_terms_in_modal_model=False,
        boundary_qualification="Finite-record delay registers can add boundary polynomials outside the modal transient span. "
        "This test uses two settling periods and repeated noiseless excitation, so boundary mismatch is numerically negligible; "
        "it does not prove omission is exact for arbitrary starting states or noisy record boundaries.")
    result["noiseless"]["passed"]=result["noiseless"]["total_relative_residual"]<1e-10
    result["all_passed"]=bool(result["noiseless"]["passed"] and result["jacobian_all_passed"])
    np.savez_compressed(dest/"simulation_validation_arrays.npz",
        frequency=f,truth=truth,test_parameters=b,U_jacobian=U,G_model=G,G_explicit=explicit,
        validation_f=a["f"],validation_k=a["k"],noiseless_U=a["U"],noiseless_Y=a["Y"],
        noiseless_model_output=y0,noiseless_residual=err,truth_delay_seconds=TRUTH_TAU)
    js(dest/"simulation_validation.json",result)
    lines=["# 完整八模态时延仿真的数值核查",
           "",
           "本文件只验证已有 `sim_core.py` 的求值、导数和数据生成，不进行参数拟合。五步 pilot 不计作正式辨识结果。",
           "",
           f"- NOFLEX={os.environ['NOFLEX']}，共 {cl.NM} 模态，模型 {m.nb} 参数。",
           f"- 模态频率：{', '.join(f'{x:.6g}' for x in cl.FM)} Hz。",
           "- 输入延迟拍数 [0,1,0,1]，输出延迟拍数 [1,0,1,0]；逐通道时延为两者之和，公共均值 200 μs。",
           f"- 与原 cl2 精确 ZOH 传函乘解析时延的相对差：{result['response_relative_error']:.3e}。",
           f"- 独立状态递推与闭环传函代数的差：Y/D 为 {result['state_y_from_d_relative_error']:.3e}，U/D 为 {result['state_u_from_d_relative_error']:.3e}。",
           f"- 闭环谱半径：{result['closed_loop_radius']:.12g}。",
           "",
           "| 导数方向 | 固定步长 1e-5 的相对误差 | 判定 |",
           "|---|---:|---|",
           *[f"| {v['group']} | {v['main_error']:.3e} | {v['passed']} |" for v in checks],
           "",
           "随机方向检查使用非零 C1、D0/D2 和两个瞬态系数，并分别检查全部七个时延坐标。"
           "JSON 保存五档差分步长的完整误差，固定主步长为 1e-5、阈值为 2e-5，未从各步长中选最小误差代替主值。",
           "",
           f"无噪声、两周期预运行后记录两周期，直接用真参数（不重新估计瞬态）得到全频点相对残差 "
           f"{result['noiseless']['total_relative_residual']:.3e}，激励线上为 "
           f"{result['noiseless']['excited_line_relative_residual']:.3e}。",
           "",
           "**边界说明：**整数时延寄存器的有限记录边界项不一定属于当前逐模态瞬态项的张成空间。"
           "本无噪声周期稳态试验中，两周期预运行使该影响在数值上可忽略；这不能证明任意初始状态或含噪记录都严格没有额外边界项。",
           "",
           f"总体核查通过：{result['all_passed']}。原文件哈希和原始数值数组分别保存在 JSON 与 NPZ。",
           ""]
    (dest/"仿真求值与导数核查.md").write_text("\n".join(lines),encoding="utf8")
    print(json.dumps({k:result[k] for k in ("number_of_modes","response_relative_error",
          "state_y_from_d_relative_error","state_u_from_d_relative_error","jacobian_all_passed","noiseless","all_passed")},
          ensure_ascii=False,default=lambda x:x.item() if isinstance(x,np.generic) else x),flush=True)

if __name__=="__main__":
    import json
    main()
