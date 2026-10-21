"""Plot the held-out dynamic FIR comparison without altering actual phase.

Input: selected_results.npz beside this file (or --input).
All P arrays are frequency x output x input, ordered [Rx, Ry, Z, N].
Port directions have fixed unit Euclidean norm. Cross-response norms below
are comparison metrics in those fixed coordinates, not physical energy.
The optional delay-aligned NN error is only a labelled diagnostic: curves
always show their actual frequency response. The explicit branch delay L
gives the NN target reference a delay of 2L; it is not a uniform delay of
the whole transformed MIMO plant.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
import numpy as np


COLORS = {"static": "#808080", "nonparam": "#195EA8", "srivc": "#7947A2"}
PORTS = ["Rx", "Ry", "Z", "N"]
FLOOR = np.finfo(float).tiny


def db(x: np.ndarray) -> np.ndarray:
    return 20.0 * np.log10(np.maximum(np.abs(x), FLOOR))


def phase_deg(x: np.ndarray) -> np.ndarray:
    return np.rad2deg(np.unwrap(np.angle(x), axis=0))


def format_x(ax: plt.Axes, f: np.ndarray) -> None:
    ticks = np.asarray([3, 5, 10, 25, 50, 100, 200, 400, 800, 1000, 2000])
    ticks = ticks[(ticks >= f[0]) & (ticks <= f[-1])]
    # Preserve readable spacing for the anticipated 25--800 Hz band.
    if f[0] >= 20:
        ticks = np.asarray([25, 50, 100, 200, 400, 800, 1000, 2000])
        ticks = ticks[(ticks >= f[0]) & (ticks <= f[-1])]
    ax.set_xscale("log")
    ax.set_xlim(f[0], f[-1])
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.grid(True, which="major", color="#d4d4d4", linewidth=0.6)
    ax.grid(True, which="minor", color="#e9e9e9", linewidth=0.45, linestyle=":")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))


def save_pair(fig: plt.Figure, out: Path, name: str) -> None:
    fig.savefig(out / f"{name}.png", dpi=190, facecolor="white")
    fig.savefig(out / f"{name}.pdf", facecolor="white")
    plt.close(fig)


def label_delay(delay: int, fs: float) -> str:
    return f"reference-branch L = {delay} samples ({1000.0 * delay / fs:.3f} ms)"


def load_route_status(path: Path) -> dict[str, bool]:
    """Read only the recorded training-only decision; never change responses."""
    if not path.is_file():
        return {"nonparam": False, "srivc": False}
    record = json.loads(path.read_text(encoding="utf-8-sig"))
    result = {}
    for route in ["nonparam", "srivc"]:
        item = record.get(route, record.get("choices", {}).get(route, {}))
        if isinstance(item, dict):
            choice = item.get("choice", item)
            name = choice.get("name", "") if isinstance(choice, dict) else choice
        else:
            name = item
        result[route] = name == "static_fallback"
    return result


def route_label(d: dict, route: str, short: bool = False) -> str:
    if d.get("static_fallback", {}).get(route, False):
        return ("SRIVC route: static fallback" if route == "srivc"
                else "Nonparametric route: static fallback")
    if route == "srivc":
        return "SRIVC FIR design" if short else "FIR designed from SRIVC model"
    return "Nonparametric FIR design" if short else "FIR designed from nonparametric FRF"


def delay_summary(d: dict) -> str:
    parts = []
    for route, title in [("nonparam", "Nonparametric"), ("srivc", "SRIVC")]:
        if d.get("static_fallback", {}).get(route, False):
            parts.append(f"{title}: static fallback (no accepted dynamic correction)")
        else:
            parts.append(f"{title}: " + label_delay(d[f"delay_{route}"], d["fs"]))
    return "; ".join(parts)


def displayed_series(d: dict, short: bool = False) -> list:
    """Do not cover the baseline with a duplicate static-fallback curve."""
    result = [(d["Pstatic"], COLORS["static"], "Strong static baseline")]
    for route in ["nonparam", "srivc"]:
        if not d.get("static_fallback", {}).get(route, False):
            result.append((d[f"P_{route}"], COLORS[route], route_label(d, route, short)))
    return result


def add_fallback_legend(ax: plt.Axes, d: dict) -> None:
    for route in ["nonparam", "srivc"]:
        if d.get("static_fallback", {}).get(route, False):
            ax.plot([], [], color=COLORS[route], linestyle="--", linewidth=1.2,
                    label=route_label(d, route) + " (= gray baseline)")


def shade_band(ax: plt.Axes, f: np.ndarray, lo: float, hi: float) -> None:
    lo, hi = max(lo, f[0]), min(hi, f[-1])
    if hi > lo:
        ax.axvspan(lo, hi, color="#d8d8d8", alpha=0.35, zorder=0, linewidth=0)


def validate(data: np.lib.npyio.NpzFile) -> dict:
    required = ["f", "Pnom", "Pstatic", "P_nonparam", "P_srivc", "fs",
                "Fu_nonparam", "Fy_nonparam", "Fu_srivc", "Fy_srivc",
                "delay_nonparam", "delay_srivc"]
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Missing NPZ fields: {missing}")
    f = np.asarray(data["f"], dtype=float).reshape(-1)
    if f.size < 2 or np.any(~np.isfinite(f)) or np.any(f <= 0):
        raise ValueError("f must contain at least two finite positive frequencies")
    order = np.argsort(f)
    f = f[order]
    if np.any(np.diff(f) <= 0):
        raise ValueError("Frequency grid must not have duplicates")
    result = {"f": f}
    for key in ["Pnom", "Pstatic", "P_nonparam", "P_srivc"]:
        value = np.asarray(data[key])
        if value.shape != (f.size, 4, 4) or not np.iscomplexobj(value):
            raise ValueError(f"{key} must be a complex (N,4,4) array, got {value.shape}")
        if not np.all(np.isfinite(value)):
            raise ValueError(f"Nonfinite response values in {key}")
        result[key] = value[order]
    for key in ["Fu_nonparam", "Fy_nonparam", "Fu_srivc", "Fy_srivc"]:
        value = np.asarray(data[key])
        if value.shape != (f.size, 3) or not np.all(np.isfinite(value)):
            raise ValueError(f"{key} must be a finite (N,3) array, got {value.shape}")
        result[key] = value[order]
    fs = float(np.asarray(data["fs"]).squeeze())
    if not np.isfinite(fs) or fs <= 0:
        raise ValueError("fs must be a finite positive scalar")
    result["fs"] = fs
    for key in ["delay_nonparam", "delay_srivc"]:
        value = float(np.asarray(data[key]).squeeze())
        if not np.isfinite(value) or value < 0 or value != round(value):
            raise ValueError(f"{key} must be a nonnegative integer number of samples")
        result[key] = int(value)
    return result


def make_full_grid(d: dict, out: Path, phase: bool = False) -> None:
    f = d["f"]
    series = displayed_series(d)
    fig, axes = plt.subplots(4, 4, figsize=(17.8, 13.0), sharex=True)
    for row in range(4):
        for col in range(4):
            ax = axes[row, col]
            for values, color, label in series:
                y = phase_deg(values[:, row, col]) if phase else db(values[:, row, col])
                ax.plot(f, y, color=color, linewidth=1.15, label=label)
            ax.set_title(f"{PORTS[col]} input $\\rightarrow$ {PORTS[row]} output", fontsize=10.5)
            format_x(ax, f)
            if row == 3:
                ax.set_xlabel("Frequency (Hz)")
            if col == 0:
                ax.set_ylabel("Unwrapped phase (deg)" if phase else "Magnitude (dB)")
            ax.tick_params(labelsize=8.5)
    kind = "Continuous phase" if phase else "Magnitude"
    fig.suptitle(f"{kind}: offline held-out FRF, training-selected decoupling routes", y=0.992, fontsize=15)
    fig.text(0.5, 0.964, delay_summary(d),
             ha="center", fontsize=10)
    add_fallback_legend(axes[0, 0], d)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.952),
               frameon=False, fontsize=10)
    footer = ("Fixed unit-norm static port directions; recorded units. "
              "Actual phase is retained. Reference-branch L is not a uniform delay of the full MIMO plant."
              if phase else
              "Fixed unit-norm static port directions; recorded units. "
              "All curves are evaluated on the same held-out FRF; this is an offline check, not a hardware deployment.")
    fig.text(0.5, 0.008, footer, ha="center", fontsize=9)
    fig.subplots_adjust(left=0.065, right=0.99, top=0.895, bottom=0.055, hspace=0.38, wspace=0.23)
    save_pair(fig, out, "full_phase" if phase else "full_magnitude")


def make_cross_and_target(d: dict, out: Path) -> None:
    f = d["f"]
    series = displayed_series(d, short=True)
    fig, axes = plt.subplots(2, 2, figsize=(14.4, 9.8), sharex=True)
    for ax in axes[0]:
        shade_band(ax, f, 25, 80)
    for ax in axes[1]:
        for lo, hi in [(260, 310), (320, 380), (450, 800)]:
            shade_band(ax, f, lo, hi)
    for values, color, label in series:
        axes[0, 0].plot(f, db(np.linalg.norm(values[:, :3, 3], axis=1)), color=color,
                        label=label, linewidth=1.4)
        axes[0, 1].plot(f, db(np.linalg.norm(values[:, 3, :3], axis=1)), color=color,
                        label=label, linewidth=1.4)
        axes[1, 0].plot(f, db(values[:, 3, 3]), color=color, label=label, linewidth=1.4)
        axes[1, 1].plot(f, phase_deg(values[:, 3, 3]), color=color, label=label, linewidth=1.4)
    axes[0, 0].set_title("N input $\\rightarrow$ three logical outputs")
    axes[0, 1].set_title("Three logical inputs $\\rightarrow$ N output")
    axes[1, 0].set_title("Target channel: N $\\rightarrow$ N magnitude")
    axes[1, 1].set_title("Target channel: N $\\rightarrow$ N actual phase")
    for ax in axes[0]:
        ax.set_ylabel("Port-normalized response norm (dB)")
    axes[1, 0].set_ylabel("Magnitude (dB)")
    axes[1, 1].set_ylabel("Unwrapped phase (deg)")
    for ax in axes.flat:
        format_x(ax, f)
    for ax in axes[1]:
        ax.set_xlabel("Frequency (Hz)")
    fig.suptitle("Cross-talk and target preservation: offline held-out FRF", fontsize=15, y=0.989)
    fig.text(0.5, 0.950, delay_summary(d),
             ha="center", fontsize=10)
    add_fallback_legend(axes[0, 0], d)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.937),
               frameon=False, fontsize=10)
    errors = []
    base = d["Pstatic"][:, 3, 3]
    for name in ["nonparam", "srivc"]:
        # This comparison aligns the *reference target* to the intentional
        # two-sided delay; it never changes the plotted response phase.
        target = base * np.exp(-2j * np.pi * f * 2 * d[f"delay_{name}"] / d["fs"])
        error = np.linalg.norm(d[f"P_{name}"][:, 3, 3] - target) / max(np.linalg.norm(target), FLOOR)
        errors.append((name, 100 * error))
    fig.text(0.5, 0.093,
             f"NN complex error relative to the static target delayed by 2L: "
             f"nonparametric {errors[0][1]:.2f}%; SRIVC {errors[1][1]:.2f}% "
             f"(unweighted grid norm, {f[0]:g}--{f[-1]:g} Hz).",
             ha="center", fontsize=9)
    fig.text(0.5, 0.066,
             "Gray bands: primary cross-talk band 25--80 Hz (top); "
             "target-preservation bands 260--310, 320--380 and 450--800 Hz (bottom).",
             ha="center", fontsize=9)
    fig.text(0.5, 0.039,
             "Port norms use fixed static scales and are comparison metrics, not physical energy. "
             "The three logical directions need not be orthogonal.",
             ha="center", fontsize=9)
    fig.text(0.5, 0.012,
             "Actual response phase is retained; 2L describes the NN target reference, not a uniform MIMO delay. "
             "This is an offline check, not hardware performance.",
             ha="center", fontsize=9)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.860, bottom=0.160, hspace=0.38, wspace=0.23)
    save_pair(fig, out, "cross_and_target")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).with_name("selected_results.npz"))
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    out = args.output_dir or args.input.parent
    if not args.input.is_file():
        raise FileNotFoundError(f"Awaiting results: {args.input}")
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.spines.top": True, "axes.spines.right": True})
    with np.load(args.input, allow_pickle=False) as archive:
        d = validate(archive)
    d["static_fallback"] = load_route_status(args.input.with_name("selected_training_only.json"))
    make_full_grid(d, out, phase=False)
    make_full_grid(d, out, phase=True)
    make_cross_and_target(d, out)
    print(f"Saved six figure files in {out.resolve()}")


if __name__ == "__main__":
    main()
