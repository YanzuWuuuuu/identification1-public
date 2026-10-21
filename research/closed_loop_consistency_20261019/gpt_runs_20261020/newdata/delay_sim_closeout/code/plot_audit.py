"""Two deterministic diagnostic figures from saved audit data; no model fitting.

Figure contracts:
1. The specified high band does not satisfy the 0.1 loop-gain prerequisite.
2. The numerical conditioning changes with band, while projected delay-column
   energy remains small. A numerical projection does not establish an exact
   ZOH non-identifiability identity.
Backend: existing Python/matplotlib workflow. Quantitative diagnostic panels,
183 mm wide, editable PDF text and 300 dpi PNG; no smoothing or dropped points.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
OUT = ROOT / "analysis"
COLORS = ["#245A81", "#5B8065", "#80659E"]
LABELS = ["3–100 Hz", "100–550 Hz", "3–550 Hz"]
BAND_KEYS = ["low_3_100", "high_100_550", "wide_3_550"]


def csv_write(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    FIG.mkdir(exist_ok=True)
    OUT.mkdir(exist_ok=True)
    sources = [ROOT / "audit/loop_gain_grid.npz", ROOT / "audit/stability_identifiability_audit.json"]
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    audit = json.loads(sources[1].read_text(encoding="utf-8"))
    with np.load(sources[0]) as data:
        frequency = data["frequency_hz"].copy()
        gain = data["sigma_max"].copy()
    use = (frequency >= 100) & (frequency <= 550)
    fx, gy = frequency[use], gain[use]
    if not (np.isfinite(gy).all() and (gy > 0).all()):
        raise ValueError("A saved loop-gain point is nonfinite or nonpositive.")
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 9, "axes.titlesize": 9,
        "axes.linewidth": .7, "axes.spines.top": False, "axes.spines.right": False,
        "legend.fontsize": 7.5, "legend.frameon": False, "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    })
    files = []

    def save(fig, name):
        for ext in ("png", "pdf"):
            path = FIG / f"{name}.{ext}"
            fig.savefig(path, dpi=300)
            files.append(path)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 3.5))
    fig.subplots_adjust(left=.095, right=.98, bottom=.22, top=.86)
    ax.semilogy(fx, gy, color=COLORS[0], lw=1.25, label=r"$\sigma_{\max}(G_0K_0)$")
    ax.axhline(.1, color="#333333", linestyle="--", linewidth=1, label="Required upper limit: 0.1")
    for modal in (286.573, 344.945, 490.390):
        ax.axvline(modal, color="#AAAAAA", linewidth=.7, linestyle=":")
        ax.text(modal + 3.0, .96, f"{modal:.1f} Hz", transform=ax.get_xaxis_transform(),
                ha="left", va="top", fontsize=7.5, color="#555555")
    peak = audit["loop"]["bands"]["high_100_550"]["sigma_max"]
    ax.plot(peak["frequency_hz"], peak["value"], "o", color=COLORS[0], ms=3.7)
    ax.annotate(f"Maximum = {peak['value']:.3f}", (peak["frequency_hz"], peak["value"]),
                xytext=(205, .32), textcoords="data", fontsize=8,
                arrowprops={"arrowstyle": "-", "color": COLORS[0], "lw": .7},
                color=COLORS[0], ha="center")
    ax.set(xlim=(100, 550), ylim=(max(gy.min()*.65, 1e-4), max(gy.max()*2.2, 1.3)),
           xlabel="Frequency (Hz)", ylabel="Largest singular value of loop transfer")
    ax.grid(axis="y", which="major", color="#E2E2E2", linewidth=.6)
    ax.legend(loc="lower left")
    fig.suptitle("The 100–550 Hz band does not satisfy the small-loop prerequisite", y=.975, fontsize=10)
    fig.text(.095, .045, "Deterministic true-model audit; no fit and no Monte Carlo averaging.\n"
             "Modal markers: 286.6, 344.9 and 490.4 Hz. All saved in-band grid points are shown.",
             fontsize=7.4, va="bottom")
    save(fig, "audit_high_band_loop_gain")
    csv_write(OUT / "figure_high_band_loop_gain.csv",
              [dict(frequency_hz=float(x), sigma_max=float(y), reference_limit=.1) for x, y in zip(fx, gy)])

    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.9), gridspec_kw={"width_ratios": [1.06, 1]})
    fig.subplots_adjust(left=.095, right=.985, bottom=.25, top=.86, wspace=.32)
    all_values, energies = [], []
    for key, label, color, style in zip(BAND_KEYS, LABELS, COLORS, ["-", "--", "-."]):
        band = audit["information"]["bands"][key]
        singular = np.array(band["full_normalized_J_singular_values"])
        fraction = np.array(band["delay_unexplained_energy_fractions"])
        if not ((singular > 0).all() and (fraction > 0).all()
                and np.isfinite(singular).all() and np.isfinite(fraction).all()):
            raise ValueError("Saved information diagnostics cannot be plotted logarithmically without alteration.")
        normalized = singular/singular[0]
        axs[0].semilogy(np.arange(1, len(singular)+1), normalized, style, color=color, lw=1.0, label=label)
        axs[1].semilogy(np.arange(7), fraction, style, color=color, marker="o", ms=3.2, lw=.9)
        all_values.extend(dict(band=key, singular_index=j+1, column_normalized_J_singular_value=float(x),
                              relative_to_largest=float(y)) for j, (x, y) in enumerate(zip(singular, normalized)))
        energies.extend(dict(band=key, delay_coordinate=j,
                             projected_unexplained_energy_fraction=float(x)) for j, x in enumerate(fraction))
    axs[0].axhline(1e-12, color="#999999", ls=":", lw=.7)
    axs[0].text(10, 2e-12, "1e−12 rank diagnostic", fontsize=7, color="#666666")
    axs[0].set(xlim=(1, 375), ylim=(1e-17, 1.5), xlabel="Singular-value index",
               ylabel=r"$\sigma_i(J_n)/\sigma_1(J_n)$", title="a  Normalized Jacobian spectrum")
    axs[0].legend(loc="lower left", bbox_to_anchor=(.015, .035))
    axs[1].set(xlim=(-.3, 6.3), ylim=(1e-29, 1e-7), xlabel="Delay coordinate",
               ylabel="Unexplained column-energy fraction",
               xticks=np.arange(7),
               xticklabels=[r"$\tau_c$", r"$a^{(1)}$", r"$a^{(2)}$", r"$a^{(3)}$",
                            r"$b^{(1)}$", r"$b^{(2)}$", r"$b^{(3)}$"],
               title="b  Delay-column projection")
    axs[1].text(3, 1e-25, "Low band: roundoff-scale residual", ha="center", fontsize=7, color=COLORS[0])
    for ax in axs:
        ax.grid(axis="y", which="major", color="#E2E2E2", linewidth=.6)
    fig.suptitle("Window choice changes numerical conditioning and residual delay information",
                 y=.97, fontsize=9.5)
    fig.text(.095, .045,
             "375 fitted directions; each Jacobian column is normalized before projection.\n"
             r"$a^{(1:3)}, b^{(1:3)}$: zero-sum output/input delay coordinates, not individual channels." "\n"
             "Numerical projections do not prove an exact ZOH non-identifiability identity.",
             fontsize=7.0, va="bottom")
    save(fig, "audit_band_information")
    csv_write(OUT / "figure_jacobian_singular_values.csv", all_values)
    csv_write(OUT / "figure_delay_projection_energy.csv", energies)
    check = {
        "backend": "Python/matplotlib", "new_fits": False, "smoothing": False,
        "source_sha256": before, "source_unchanged": all(
            hashlib.sha256(p.read_bytes()).hexdigest() == before[str(p)] for p in sources),
        "loop_points_shown": len(fx), "high_band_audit_maximum": peak,
        "pdf_text_fonttype": 42, "artifacts": [],
        "interpretation": "Numerical deterministic diagnostics; no exact ZOH non-identifiability proof.",
    }
    for path in files:
        item = dict(path=str(path.relative_to(ROOT)), bytes=path.stat().st_size)
        if path.suffix == ".png":
            with Image.open(path) as im:
                item.update(width_pixels=im.width, height_pixels=im.height, dpi=im.info.get("dpi"))
        check["artifacts"].append(item)
    (OUT / "diagnostic_figure_checks.json").write_text(json.dumps(check, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(check, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
