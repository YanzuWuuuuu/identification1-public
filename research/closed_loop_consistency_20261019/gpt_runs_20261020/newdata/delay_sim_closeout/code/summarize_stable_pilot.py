"""Summarize isolated solver pilots; never merge with the formal MC."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"stable_profile_pilot"


def main():
    rows = []
    for p in sorted(OUT.glob("*/result.json")):
        r = json.loads(p.read_text(encoding="utf8"))
        s = r.get("endpoint", {}).get("solver", {})
        rows.append(dict(name=p.parent.name, result_file=str(p),
                         score=s.get("score"), success=s.get("success", False),
                         status=s.get("status", r.get("status")),
                         linear_score=s.get("original_linear_score"),
                         near_bound=s.get("near_bound"),
                         calls=s.get("profile_calls"),
                         all_calls=s.get("total_profile_calls_including_diagnostics",
                                         s.get("profile_calls")),
                         seconds=s.get("total_seconds_including_diagnostics",
                                       s.get("seconds")),
                         frequency_hz=r.get("endpoint", {}).get("frequency_hz"),
                         damping_percent=r.get("endpoint", {}).get("damping_percent"),
                         error_percent=r.get("endpoint", {}).get("error_percent"),
                         delay=r.get("endpoint", {}).get("delay")))
    (OUT/"pilot_summary.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf8")
    lines = [
        "# 外层数值修复 pilot 结果", "",
        "全部端点单独保存，不计入正式 Monte Carlo。"
        "最终门槛保持原375维score<1e-6、内层线性score<1e-9、"
        "有限参数且未贴原边界；没有删去窗外模态。", "",
        "| pilot | 完整score | 线性score | 原门槛通过 | 贴边 | 求解/全部profile次数 | 秒 | 状态 |",
        "|---|---:|---:|:---:|:---:|:---:|---:|---|"]
    for r in rows:
        score = "—" if r["score"] is None else f"{r['score']:.6g}"
        linear = "—" if r["linear_score"] is None else f"{r['linear_score']:.4g}"
        seconds = "—" if r["seconds"] is None else f"{r['seconds']:.2f}"
        lines.append(f"| {r['name']} | {score} | {linear} | "
                     f"{r['success']} | {r['near_bound']} | "
                     f"{r['calls']}/{r['all_calls']} | {seconds} | {r['status']} |")
    lines += [
        "", "## 可以说明什么", "",
        "H/v2由初始score 0.37347降至6.8077e-6，"
        "W/v1由0.36473降至0.0013743。"
        "这证明此前全部停留在约0.3–0.4不能只归因于模型弱辨识；"
        "外层正规截断和共同缩步是实质性的数值障碍。", "",
        "但参数可辨识性、数值找到原矩的有效根和统计一致性仍是不同判断。"
        "未达到原375维门槛的端点继续计作失败；"
        "窗内三个模态误差较小不能代替完整验收。", "",
        "投影GN与完整中央差分收尾明确区分。"
        "central的前12次方向检查也计入实际调用与时间，"
        "不能只报告scipy的nfev。每次试探点均有eval文件，"
        "最终结果使用最后接受点，不择优取最小score的试探点。", "",
        "W/v1的360秒时间预算实际在410秒退出，"
        "原因是一次在途profile跨越时限后又执行了投影。"
        "v1源码保留；后续版本补了投影前后时间检查。"
        "此修改不改变模型、矩方程和统计判据。", "",
        "复现命令及保存规则见同目录README.md；"
        "投影等价性见../audit/stable_schur_equivalence.json。"]
    (OUT/"数值修复pilot说明.md").write_text("\n".join(lines)+"\n", encoding="utf8")


if __name__ == "__main__":
    main()
