#!/usr/bin/env python3
"""Build THE-85 Bayesian-iteration stability bridge slide.

This slide is meant to follow the reconstructed xJ subtraction-input slide and
lead into the unfolded xJ result.  It mirrors the PPG12/Ian iteration-scan
logic: compare each Bayesian iteration to the previous iteration, compare that
movement with the statistical uncertainty, and choose the first stable
iteration count for the final unfolding.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_unfolded_xjgamma_1x3 as u  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PNG = OUT_DIR / "slide11_bayes_iteration_stability_bridge.png"
OUT_CSV = OUT_DIR / "slide11_bayes_iteration_stability_bridge_scan.csv"
OUT_MANIFEST = OUT_DIR / "slide11_bayes_iteration_stability_bridge_manifest.json"
OUT_SCRIPT = OUT_DIR / "slide11_bayes_iteration_stability_bridge_speaker_script.md"

ITERATIONS = [int(x) for x in os.environ.get("THE85_ITER_SCAN", "1,2,3,4,5,6,7,8,9,10").split(",") if x.strip()]
DISPLAY_MAX_ITER = int(os.environ.get("THE85_ITER_DISPLAY_MAX", "7"))
NTOYS = int(os.environ.get("THE85_ITER_NTOYS", "300"))
REQUIRE_FULL_PT_BINS = os.environ.get("THE85_REQUIRE_FULL_PT_BINS", "0").strip().lower() not in {"0", "false", "no"}

# Match the reconstructed-input slide currently under review.
u.REQUIRE_FULL_PT_BINS = REQUIRE_FULL_PT_BINS


def selected_cases() -> List[u.Case]:
    out = []
    for case in u.CASES:
        if case.key in {"auau_0_20", "pp_basev3e"}:
            out.append(case)
    if len(out) != 2:
        raise RuntimeError("expected auau_0_20 and pp_basev3e cases")
    return out


def h1_selected_values(h) -> Dict[str, np.ndarray]:
    xs: List[float] = []
    widths: List[float] = []
    vals: List[float] = []
    errs: List[float] = []
    selected_bins: List[Tuple[float, float]] = []
    for ib in range(1, h.GetXaxis().GetNbins() + 1):
        lo = float(h.GetXaxis().GetBinLowEdge(ib))
        hi = float(h.GetXaxis().GetBinUpEdge(ib))
        cen = float(h.GetXaxis().GetBinCenter(ib))
        if not u.row_in_pt_window(lo, hi, cen):
            continue
        xs.append(cen)
        widths.append(hi - lo)
        vals.append(float(h.GetBinContent(ib)))
        errs.append(float(h.GetBinError(ib)))
        selected_bins.append((lo, hi))
    return {
        "x": np.asarray(xs, dtype=float),
        "width": np.asarray(widths, dtype=float),
        "y": np.asarray(vals, dtype=float),
        "ey": np.asarray(errs, dtype=float),
        "bins": np.asarray(selected_bins, dtype=float),
    }


def scaled_like(current: np.ndarray, prior: np.ndarray, widths: np.ndarray | None) -> np.ndarray:
    if widths is None:
        curr_area = float(np.nansum(current))
        prior_area = float(np.nansum(prior))
    else:
        curr_area = float(np.nansum(current * widths))
        prior_area = float(np.nansum(prior * widths))
    if not math.isfinite(curr_area) or not math.isfinite(prior_area) or abs(prior_area) < 1e-12:
        return prior.copy()
    return prior * (curr_area / prior_area)


def relative_quadrature(current: np.ndarray, errors: np.ndarray, previous: np.ndarray) -> Dict[str, float]:
    current = np.asarray(current, dtype=float)
    previous = np.asarray(previous, dtype=float)
    errors = np.asarray(errors, dtype=float)
    good = np.isfinite(current) & np.isfinite(previous) & np.isfinite(errors)
    good &= (np.abs(current) > 1e-12) | (np.abs(previous) > 1e-12) | (np.abs(errors) > 1e-12)
    # Match the C++ RooUnfold QA: form one aggregate L2 relative quantity over
    # the physics bins.  A bin-by-bin sum of e_i/U_i lets near-empty Au+Au tail
    # bins dominate the diagnostic and is not a stable iteration-choice metric.
    if good.size >= 5:
        interior = np.zeros_like(good, dtype=bool)
        interior[1:-2] = True
        good &= interior
    if not np.any(good):
        return {"stat": 0.0, "movement": 0.0, "quadrature": 0.0, "bins": 0}
    denom = float(np.nansum(current[good] * current[good]))
    if not math.isfinite(denom) or denom <= 0.0:
        return {"stat": 0.0, "movement": 0.0, "quadrature": 0.0, "bins": int(np.count_nonzero(good))}
    stat = float(np.sqrt(max(0.0, np.nansum(errors[good] * errors[good]) / denom)))
    diff = current[good] - previous[good]
    movement = float(np.sqrt(max(0.0, np.nansum(diff * diff) / denom)))
    return {
        "stat": stat,
        "movement": movement,
        "quadrature": float(math.sqrt(stat * stat + movement * movement)),
        "bins": int(np.count_nonzero(good)),
    }


def xj_metric_mask(x: np.ndarray) -> np.ndarray:
    # Avoid pure turn-on / underflow regions in the stability metric.  The
    # plotted final curve may show more range, but the iteration choice should
    # be driven by the physical recoil-balance region.
    return np.isfinite(x) & (x >= 0.20) & (x <= 2.14)


def make_prior_objects(case: u.Case, sim_f) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    pho_truth = u.get_obj(sim_f, case.sim_topdir, f"h_unfoldTruthPho_pTgamma{case.cent_suffix}", "TH1")
    xj_truth = u.get_obj(sim_f, case.sim_topdir, f"h2_unfoldTruth_pTgamma_xJ_incl_{u.BASE_KEY}{case.cent_suffix}", "TH2")
    photon_prior = h1_selected_values(pho_truth)
    xj_prior_full = u.project_per_photon_xj(xj_truth, pho_truth)
    mask = xj_metric_mask(xj_prior_full["x_centers"])
    xj_prior = {
        "x": np.asarray(xj_prior_full["x_centers"], dtype=float)[mask],
        "width": np.asarray(xj_prior_full["x_widths"], dtype=float)[mask],
        "y": np.asarray(xj_prior_full["y"], dtype=float)[mask],
        "ey": np.asarray(xj_prior_full["ey"], dtype=float)[mask],
    }
    return photon_prior, xj_prior


def run_one_iteration(case: u.Case, iteration: int) -> Dict[str, Dict[str, np.ndarray]]:
    u.DEFAULT_ITERS = iteration
    u.NTOYS_FINAL = NTOYS
    u.ROOT.gRandom.SetSeed(20260629 + 100 * iteration + (0 if case.key == "auau_0_20" else 17))
    data_f = u.open_root(case.data_file)
    sim_f = u.open_root(case.sim_file)
    try:
        h_pho_unf, pho_meta = u.unfold_photons(case, data_f, sim_f)
        h2_xj_unf, xj_meta = u.unfold_xj(case, data_f, sim_f)
        xj = u.project_per_photon_xj(h2_xj_unf, h_pho_unf)
        xj_mask = xj_metric_mask(np.asarray(xj["x_centers"], dtype=float))
        return {
            "photon": h1_selected_values(h_pho_unf),
            "xj2d": {
                "x": np.asarray(xj["x_centers"], dtype=float)[xj_mask],
                "width": np.asarray(xj["x_widths"], dtype=float)[xj_mask],
                "y": np.asarray(xj["y"], dtype=float)[xj_mask],
                "ey": np.asarray(xj["ey"], dtype=float)[xj_mask],
            },
            "meta": {
                "photon": pho_meta,
                "xj": xj_meta,
                "selected_truth_pt_bins": xj["selected_truth_pt_bins"],
                "npho_unfolded": xj["npho_unfolded"],
            },
        }
    finally:
        data_f.Close()
        sim_f.Close()


def scan_case(case: u.Case) -> Dict:
    sim_f = u.open_root(case.sim_file)
    try:
        photon_prior, xj_prior = make_prior_objects(case, sim_f)
    finally:
        sim_f.Close()

    scans: Dict[str, List[Dict]] = {"photon": [], "xj2d": []}
    series: Dict[str, Dict[int, Dict[str, np.ndarray]]] = {"photon": {}, "xj2d": {}}
    metas: Dict[int, Dict] = {}

    previous: Dict[str, Dict[str, np.ndarray] | None] = {"photon": None, "xj2d": None}
    prior_by_kind = {"photon": photon_prior, "xj2d": xj_prior}
    widths_by_kind = {"photon": None, "xj2d": xj_prior["width"]}

    for iteration in ITERATIONS:
        result = run_one_iteration(case, iteration)
        metas[iteration] = result["meta"]
        for kind in ("photon", "xj2d"):
            current = result[kind]
            series[kind][iteration] = current
            if previous[kind] is None:
                prior_y = scaled_like(current["y"], prior_by_kind[kind]["y"], widths_by_kind[kind])
                prev_y = prior_y
                prev_label = "prior"
            else:
                prev_y = previous[kind]["y"]
                prev_label = f"iter {iteration - 1}"
            metrics = relative_quadrature(current["y"], current["ey"], prev_y)
            scans[kind].append(
                {
                    "case": case.key,
                    "kind": kind,
                    "iteration": iteration,
                    "previous": prev_label,
                    **metrics,
                    "movement_over_stat": metrics["movement"] / metrics["stat"] if metrics["stat"] > 0 else math.inf,
                }
            )
            previous[kind] = current

    return {"case": case, "scans": scans, "series": series, "prior": prior_by_kind, "meta": metas}


def choose_kind_iteration(rows: Sequence[Dict]) -> int:
    eligible = [r for r in rows if 3 <= int(r["iteration"]) <= DISPLAY_MAX_ITER]
    if not eligible:
        return int(rows[-1]["iteration"])
    return int(min(eligible, key=lambda r: r["quadrature"])["iteration"])


def choose_global_iteration(all_scans: Dict[str, Dict]) -> Tuple[int, Dict[str, Dict[str, int]]]:
    choices: Dict[str, Dict[str, int]] = {}
    for case_key, payload in all_scans.items():
        choices[case_key] = {
            "photon": choose_kind_iteration(payload["scans"]["photon"]),
            "xj2d": choose_kind_iteration(payload["scans"]["xj2d"]),
        }
    # Use one iteration count for the delivered physics curve.  PPG12 uses a
    # small integer stability plateau; retaining at least 3 avoids stopping at
    # the first barely-stable correction step.
    max_choice = max(v for by_kind in choices.values() for v in by_kind.values())
    return max(3, max_choice), choices


def case_display_name(case_key: str) -> str:
    if case_key == "auau_0_20":
        return "Au+Au 0-20%"
    if case_key == "pp_basev3e":
        return "p+p"
    return case_key


def kind_display_name(kind: str) -> str:
    return "photon 1D" if kind == "photon" else r"$(p_T^\gamma,x_{J\gamma})$ 2D"


def metric_at(rows: Sequence[Dict], iteration: int) -> Dict:
    for row in rows:
        if int(row["iteration"]) == int(iteration):
            return row
    return rows[-1]


def write_csv(all_scans: Dict[str, Dict]) -> None:
    rows: List[Dict] = []
    for payload in all_scans.values():
        for kind_rows in payload["scans"].values():
            rows.extend(kind_rows)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "case",
                "kind",
                "iteration",
                "previous",
                "stat",
                "movement",
                "quadrature",
                "movement_over_stat",
                "bins",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
        }
    )


def add_sphenix_internal(ax, x: float = 0.04, y: float = 0.84, size: float = 12.5) -> None:
    ax.text(
        x,
        y,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=size,
        color="#111827",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.2),
    )


def metric_arrays(rows: Sequence[Dict]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    shown = [r for r in rows if int(r["iteration"]) <= DISPLAY_MAX_ITER]
    x = np.asarray([r["iteration"] for r in shown], dtype=float)
    stat = np.asarray([r["stat"] for r in shown], dtype=float)
    movement = np.asarray([r["movement"] for r in shown], dtype=float)
    quad = np.asarray([r["quadrature"] for r in shown], dtype=float)
    return x, stat, movement, quad


def draw_ppg12_style_axis(ax, payload: Dict, case_key: str, kind: str, chosen: int) -> None:
    rows = payload["scans"][kind]
    x, stat, movement, quad = metric_arrays(rows)

    ax.errorbar(x, stat, yerr=np.zeros_like(stat), fmt="o", ms=5.8, color="black", lw=1.0, label="total relative stat. uncertainty", zorder=4)
    ax.errorbar(x, movement, yerr=np.zeros_like(movement), fmt="o", ms=5.8, color="#0b22ff", lw=1.0, label="total relative deviation", zorder=4)
    ax.errorbar(x, quad, yerr=np.zeros_like(quad), fmt="o", ms=5.8, color="#ff160a", lw=1.0, label="stat. + dev. (quadrature)", zorder=4)
    ax.axvspan(chosen - 0.18, chosen + 0.18, color="#fee2e2", alpha=0.65, zorder=0)
    ax.axvline(chosen, color="#dc2626", lw=1.4, ls=":")
    ax.text(chosen + 0.08, 0.93, f"use {chosen}", transform=ax.get_xaxis_transform(), fontsize=12.0, color="#dc2626", va="top")

    ymax = float(np.nanmax(quad[np.isfinite(quad)])) if np.any(np.isfinite(quad)) else 1.0
    ymax = max(0.7, min(1.15, ymax * 1.15))
    ax.set_xlim(0, max(ITERATIONS) + 0.05)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Iteration", fontsize=17, fontweight="bold", labelpad=10)
    ax.set_ylabel(r"$\sqrt{\delta}$", fontsize=18, labelpad=8)
    ax.set_title(f"{case_display_name(case_key)}  {kind_display_name(kind)}", fontsize=15.5, fontweight="bold", pad=8)
    ax.text(
        0.50,
        0.94,
        r"$\it{\bf{sPHENIX}}$ Internal" "\n"
        + (r"$p{+}p,\ \sqrt{s}=200$ GeV" if case_key == "pp_basev3e" else r"Au+Au 0-20%, $\sqrt{s_{NN}}=200$ GeV")
        + "\n"
        + r"$15 < E_T^\gamma < 35$ GeV, $|\Delta\phi|>7\pi/8$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13.4,
        color="#111827",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=1.6),
    )
    ax.legend(loc="upper right", bbox_to_anchor=(0.985, 0.735), frameon=False, fontsize=12.0, handlelength=1.0, borderpad=0.2)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=14.5, direction="in", top=True, right=True, length=8)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=4)
    for spine in ax.spines.values():
        spine.set_linewidth(1.15)


def draw_compact_metric_axis(ax, payload: Dict, case_key: str, kind: str, chosen: int) -> None:
    rows = payload["scans"][kind]
    x, stat, movement, quad = metric_arrays(rows)
    ax.scatter(x, movement / np.maximum(stat, 1e-12), s=32, color="#2563eb", marker="s", label=r"$\delta/\sigma$")
    ax.axhline(1.0, color="#64748b", lw=1.0, ls="--")
    ax.axvspan(chosen - 0.16, chosen + 0.16, color="#fee2e2", alpha=0.60, zorder=0)
    ax.set_yscale("log")
    ax.set_ylim(5e-4, 8)
    ax.set_xlim(0, max(ITERATIONS) + 0.05)
    ax.set_title(f"{case_display_name(case_key)} {kind_display_name(kind)}", fontsize=11.2, fontweight="bold", pad=4)
    ax.tick_params(axis="both", which="major", labelsize=9.2, direction="in", top=True, right=True, length=5)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=2.5)
    ax.grid(True, which="major", color="#e5e7eb", lw=0.65)
    ax.grid(True, which="minor", color="#f1f5f9", lw=0.35)


def draw_ppg12_panel(ax, payload: Dict, case_key: str, kind: str, chosen: int, *, show_legend: bool = False) -> None:
    rows = payload["scans"][kind]
    x, stat, movement, quad = metric_arrays(rows)
    local_choice = choose_kind_iteration(rows)
    local_row = metric_at(rows, local_choice)
    common_row = metric_at(rows, chosen)

    ax.plot(
        x,
        quad,
        "-o",
        ms=5.4,
        color="#111827",
        lw=1.45,
        label="quadrature stability score",
        zorder=4,
    )
    stable_x = [r["iteration"] for r in rows if r["iteration"] >= 3 and r["iteration"] <= DISPLAY_MAX_ITER and r["movement"] <= r["stat"]]
    stable_y = [r["quadrature"] for r in rows if r["iteration"] >= 3 and r["iteration"] <= DISPLAY_MAX_ITER and r["movement"] <= r["stat"]]
    if stable_x:
        ax.scatter(stable_x, stable_y, s=54, facecolors="none", edgecolors="#2563eb", lw=1.25, label="movement < stat.", zorder=5)

    ax.axvspan(chosen - 0.16, chosen + 0.16, color="#fee2e2", alpha=0.58, zorder=0)
    ax.axvline(chosen, color="#dc2626", lw=1.35, ls=":")
    if local_choice != chosen:
        ax.axvline(local_choice, color="#64748b", lw=1.05, ls="--", alpha=0.80)
    ax.scatter([local_choice], [local_row["quadrature"]], marker="*", s=125, color="#2563eb", edgecolors="white", lw=0.6, zorder=6)

    finite_quad = quad[np.isfinite(quad) & (quad > 0)]
    if finite_quad.size and float(np.nanmax(finite_quad)) / max(float(np.nanmin(finite_quad)), 1e-9) > 18.0:
        ax.set_yscale("log")
        ymin = max(0.015, float(np.nanmin(finite_quad)) * 0.65)
        ymax = float(np.nanmax(finite_quad)) * 1.45
    else:
        ymin = 0.0
        ymax = max(0.09, float(np.nanmax(finite_quad)) * 1.28 if finite_quad.size else 1.0)
    ax.set_xlim(0, DISPLAY_MAX_ITER + 0.08)
    ax.set_ylim(ymin, ymax)

    title = f"{case_display_name(case_key)}   {kind_display_name(kind)}"
    ax.set_title(title, fontsize=12.8, fontweight="bold", pad=6)
    ax.set_xlabel("Iteration", fontsize=12.0, fontweight="bold", labelpad=4)
    ax.set_ylabel("Stability score", fontsize=12.5, labelpad=3)
    ax.text(
        0.045,
        0.905,
        f"first stable {local_choice}",
        transform=ax.transAxes,
        fontsize=10.8,
        fontweight="bold",
        color="#2563eb",
        ha="left",
        va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.0),
    )
    ax.text(
        0.965,
        0.905,
        f"move/stat @ {chosen}: {common_row['movement_over_stat']:.2f}",
        transform=ax.transAxes,
        fontsize=9.8,
        color="#334155",
        ha="right",
        va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.0),
    )
    if show_legend:
        ax.legend(
            loc="upper right",
            bbox_to_anchor=(0.985, 0.755),
            frameon=False,
            fontsize=9.4,
            handlelength=1.0,
            borderpad=0.1,
            labelspacing=0.35,
        )
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=10.8, direction="in", top=True, right=True, length=6)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=3)
    ax.grid(True, which="major", color="#e2e8f0", lw=0.6, alpha=0.75)
    for spine in ax.spines.values():
        spine.set_linewidth(1.05)


def draw_summary_table(ax, all_scans: Dict[str, Dict], choices: Dict[str, Dict[str, int]], global_choice: int) -> None:
    ax.axis("off")
    ax.text(0.0, 0.98, "Iteration-choice summary", fontsize=23.5, fontweight="bold", color="#111827", va="top")
    result_box = FancyBboxPatch((0.0, 0.840), 1.0, 0.075, boxstyle="round,pad=0.010,rounding_size=0.014", facecolor="#f8fafc", edgecolor="#cbd5e1", lw=1.0)
    ax.add_patch(result_box)
    ax.text(0.035, 0.878, "Nominal choice", fontsize=13.0, fontweight="bold", color="#475569", va="center")
    ax.text(0.98, 0.878, f"{global_choice} iterations", fontsize=18.0, fontweight="bold", color="#111827", va="center", ha="right")

    x0, y0 = 0.0, 0.740
    widths = [0.27, 0.17, 0.17, 0.20, 0.19]
    headers = ["sample", r"$N_\gamma$ 1D", r"$x_J$ 2D", r"$\delta/\sigma$ at 4", "role"]
    row_h = 0.115
    for i, (w, h) in enumerate(zip(widths, headers)):
        ax.add_patch(Rectangle((x0 + sum(widths[:i]), y0), w, row_h, facecolor="#f8fafc", edgecolor="#cbd5e1", lw=1.0))
        ax.text(x0 + sum(widths[:i]) + w / 2, y0 + row_h / 2, h, fontsize=12.6, fontweight="bold", ha="center", va="center", color="#111827")
    y = y0 - row_h
    for case_key in ("auau_0_20", "pp_basev3e"):
        payload = all_scans[case_key]
        metric_pho = metric_at(payload["scans"]["photon"], global_choice)
        metric_xj = metric_at(payload["scans"]["xj2d"], global_choice)
        ratio_text = f"{metric_pho['movement_over_stat']:.2f} / {metric_xj['movement_over_stat']:.2f}"
        status = "passes" if case_key == "auau_0_20" else "sets count"
        cells = [
            case_display_name(case_key),
            str(choices[case_key]["photon"]),
            str(choices[case_key]["xj2d"]),
            ratio_text,
            status,
        ]
        for i, (w, cell) in enumerate(zip(widths, cells)):
            face = "white" if case_key == "auau_0_20" else "#eff6ff"
            ax.add_patch(Rectangle((x0 + sum(widths[:i]), y), w, row_h, facecolor=face, edgecolor="#cbd5e1", lw=1.0))
            fw = "bold" if i == 0 or i == 4 else "normal"
            color = "#1d4ed8" if i == 4 and case_key == "pp_basev3e" else "#111827"
            ax.text(x0 + sum(widths[:i]) + w / 2, y + row_h / 2, cell, fontsize=12.6, fontweight=fw, ha="center", va="center", color=color)
        y -= row_h

    note_y = 0.425
    box = FancyBboxPatch((0.0, note_y - 0.322), 1.0, 0.322, boxstyle="round,pad=0.012,rounding_size=0.018", facecolor="white", edgecolor="#94a3b8", lw=1.3)
    ax.add_patch(box)
    ax.text(0.035, note_y - 0.033, "How the choice is made", fontsize=17.5, fontweight="bold", color="#111827", va="top")
    ax.text(
        0.035,
        note_y - 0.095,
        r"$\delta_n=\sqrt{\sum_i[(U_n-U_{n-1})/U_n]^2}$"
        "\n"
        r"$\sigma_n=\sqrt{\sum_i(e_n/U_n)^2}$"
        "\n"
        r"Stable means $\delta_n/\sigma_n < 1$."
        "\n"
        r"ATLAS-style counts are 2-4; pp 2D is the limiting check here.",
        fontsize=12.6,
        color="#374151",
        va="top",
        linespacing=1.14,
    )

    ax.text(
        0.035,
        0.118,
        "Points are discrete iteration tests, not a physical curve.",
        fontsize=11.9,
        color="#475569",
        va="center",
        fontstyle="italic",
    )

    caveat = FancyBboxPatch((0.0, 0.000), 1.0, 0.085, boxstyle="round,pad=0.010,rounding_size=0.016", facecolor="#f8fafc", edgecolor="#cbd5e1", lw=1.0)
    ax.add_patch(caveat)
    ax.text(
        0.035,
        0.043,
        "Current histogram edges make this a full-bin 16-35 GeV validation.\n"
        "Exact 15 GeV needs a future rerun with a true 15 GeV bin edge.",
        fontsize=11.1,
        color="#475569",
        va="center",
        linespacing=1.15,
    )


def draw_slide(all_scans: Dict[str, Dict], choices: Dict[str, Dict[str, int]], global_choice: int) -> None:
    setup_style()
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")

    fig.text(
        0.055,
        0.935,
        f"Bayesian unfolding iteration choice: common count = {global_choice}",
        fontsize=31,
        fontweight="bold",
        ha="left",
        va="center",
        color="#111827",
    )

    ax = fig.add_axes([0.075, 0.275, 0.555, 0.560])
    curve_specs = [
        ("Au+Au photon 1D", "auau_0_20", "photon", "#dc2626", "o"),
        (r"Au+Au $(p_T^\gamma,x_J)$ 2D", "auau_0_20", "xj2d", "#f97316", "s"),
        ("p+p photon 1D", "pp_basev3e", "photon", "#2563eb", "^"),
        (r"p+p $(p_T^\gamma,x_J)$ 2D", "pp_basev3e", "xj2d", "#111827", "D"),
    ]
    for label, case_key, kind, color, marker in curve_specs:
        rows = [r for r in all_scans[case_key]["scans"][kind] if 3 <= int(r["iteration"]) <= DISPLAY_MAX_ITER]
        x = np.asarray([r["iteration"] for r in rows], dtype=float)
        q = np.asarray([r["quadrature"] for r in rows], dtype=float)
        qmin = float(np.nanmin(q)) if q.size else 1.0
        y = q / qmin if qmin > 0 else q
        local_choice = choices[case_key][kind]
        ax.plot(x, y, "-"+marker, color=color, lw=2.0, ms=6.0, label=f"{label}: min at {local_choice}")
        local_y = float(y[np.where(x == local_choice)][0]) if np.any(x == local_choice) else 1.0
        ax.scatter([local_choice], [local_y], marker="*", s=155, color=color, edgecolors="white", lw=0.8, zorder=7)

    ax.axvline(global_choice, color="#dc2626", lw=1.8, ls=":")
    ax.axvspan(global_choice - 0.12, global_choice + 0.12, color="#fee2e2", alpha=0.70, zorder=0)
    ax.text(global_choice + 0.06, 1.78, f"common {global_choice}", fontsize=13.0, fontweight="bold", color="#dc2626", ha="left", va="center")
    ax.set_xlim(2.85, DISPLAY_MAX_ITER + 0.15)
    ax.set_ylim(0.92, 1.90)
    ax.set_xlabel("Bayes iteration", fontsize=16, fontweight="bold", labelpad=8)
    ax.set_ylabel("Quadrature score / local minimum", fontsize=15, labelpad=8)
    ax.set_title("Decision scan after startup iterations", fontsize=17.0, fontweight="bold", pad=10)
    ax.grid(True, which="major", color="#e2e8f0", lw=0.8)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=13.2, direction="in", top=True, right=True, length=7)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=3.5)
    ax.legend(loc="upper left", frameon=False, fontsize=11.2, handlelength=2.0, labelspacing=0.55)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    # Right-side explanation: formula, local minima, and common-count rule.
    right_x = 0.670
    fig.text(right_x, 0.805, "What is plotted", fontsize=20, fontweight="bold", color="#111827", ha="left")
    fig.text(
        right_x,
        0.740,
        r"$Q_n = \sqrt{S_n^2 + M_n^2}$",
        fontsize=22,
        color="#111827",
        ha="left",
    )
    fig.text(
        right_x,
        0.685,
        r"$S_n = \sqrt{\sum_i e_{n,i}^2 / \sum_i U_{n,i}^2}$",
        fontsize=15.0,
        color="#334155",
        ha="left",
    )
    fig.text(
        right_x,
        0.635,
        r"$M_n = \sqrt{\sum_i (U_{n,i}-U_{n-1,i})^2 / \sum_i U_{n,i}^2}$",
        fontsize=15.0,
        color="#334155",
        ha="left",
    )
    fig.text(
        right_x,
        0.575,
        r"$U_n$ is unfolded yield; $e_n$ is kCovToy statistical error.",
        fontsize=12.2,
        color="#475569",
        ha="left",
    )

    table_x, table_y = right_x, 0.515
    col_w = [0.180, 0.052, 0.083]
    row_h = 0.050
    headers = ["check", "min", r"$Q_4/Q_{min}$"]
    for i, (w, h) in enumerate(zip(col_w, headers)):
        fig.patches.append(Rectangle((table_x + sum(col_w[:i]), table_y), w, row_h, transform=fig.transFigure, facecolor="#f8fafc", edgecolor="#cbd5e1", lw=0.9))
        fig.text(table_x + sum(col_w[:i]) + w / 2, table_y + row_h / 2, h, fontsize=10.6, fontweight="bold", ha="center", va="center", color="#111827")
    y0 = table_y - row_h
    for idx, (label, case_key, kind, color, _marker) in enumerate(curve_specs):
        rows = [r for r in all_scans[case_key]["scans"][kind] if 3 <= int(r["iteration"]) <= DISPLAY_MAX_ITER]
        q = np.asarray([r["quadrature"] for r in rows], dtype=float)
        qmin = float(np.nanmin(q)) if q.size else 1.0
        q4 = metric_at(rows, global_choice)["quadrature"] / qmin if qmin > 0 else 0.0
        table_label = {
            ("auau_0_20", "photon"): "Au+Au photon 1D",
            ("auau_0_20", "xj2d"): "Au+Au 2D",
            ("pp_basev3e", "photon"): "p+p photon 1D",
            ("pp_basev3e", "xj2d"): "p+p 2D",
        }[(case_key, kind)]
        cells = [table_label, str(choices[case_key][kind]), f"{q4:.2f}"]
        for i, (w, cell) in enumerate(zip(col_w, cells)):
            fig.patches.append(Rectangle((table_x + sum(col_w[:i]), y0), w, row_h, transform=fig.transFigure, facecolor="white", edgecolor="#cbd5e1", lw=0.9))
            fig.text(table_x + sum(col_w[:i]) + (0.010 if i == 0 else w / 2), y0 + row_h / 2, cell, fontsize=10.7, ha=("left" if i == 0 else "center"), va="center", color=(color if i == 0 else "#111827"), fontweight=("bold" if i == 1 else "normal"))
        y0 -= row_h

    fig.text(right_x, 0.220, "Choice", fontsize=19, fontweight="bold", color="#111827", ha="left")
    fig.text(
        right_x,
        0.180,
        f"Use one iteration count for all delivered spectra.\nThe limiting check is p+p 2D, whose minimum is {global_choice}.",
        fontsize=13.8,
        color="#334155",
        ha="left",
        linespacing=1.25,
    )

    readout = FancyBboxPatch(
        (0.070, 0.055),
        0.845,
        0.080,
        transform=fig.transFigure,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        facecolor="#f8fafc",
        edgecolor="#cbd5e1",
        lw=1.0,
    )
    fig.patches.append(readout)
    fig.text(
        0.090,
        0.095,
        "Choice rule",
        fontsize=15.0,
        fontweight="bold",
        color="#111827",
        ha="left",
        va="center",
    )
    fig.text(
        0.185,
        0.095,
        r"Minimize the quadrature score in the decision window (iterations 3-7); choose one common count, set by the limiting check.",
        fontsize=12.9,
        color="#334155",
        ha="left",
        va="center",
    )

    fig.text(
        0.070,
        0.024,
        f"RooUnfoldBayes with kCovToy errors; scan uses {NTOYS} toys per point. "
        f"Display stops at iteration {DISPLAY_MAX_ITER}; later Au+Au 2D points are retained in the scan and expose tail-bin covariance blow-up, not a better choice.",
        fontsize=10.9,
        color="#475569",
        ha="left",
    )
    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)


def write_manifest(all_scans: Dict[str, Dict], choices: Dict[str, Dict[str, int]], global_choice: int) -> None:
    manifest = {
        "script": str(Path(__file__).relative_to(REPO)),
        "output_png": str(OUT_PNG),
        "output_csv": str(OUT_CSV),
        "roo_unfold": {
            "method": "RooUnfoldBayes",
            "recommended_iterations": global_choice,
            "covariance": "kCovToy",
            "ntoys_scan": NTOYS,
            "iterations_scanned": ITERATIONS,
            "iterations_displayed_on_slide": [i for i in ITERATIONS if i <= DISPLAY_MAX_ITER],
        },
        "pt_window": {
            "nominal": [15.0, 35.0],
            "require_full_bins": REQUIRE_FULL_PT_BINS,
            "effective_note": (
                "If require_full_bins is false, the current center-based selector includes the 14-16 GeV row. "
                "If require_full_bins is true, the current available full-bin range starts at 16 GeV. "
                "A final exact 15 GeV edge requires a rerun/replay with that bin edge."
            ),
        },
        "metric": {
            "movement": "sqrt(sum_i (U_n - U_{n-1})^2 / sum_i U_n^2), aggregate L2 metric matching the C++ RooUnfold QA",
            "stat": "sqrt(sum_i e_n^2 / sum_i U_n^2), aggregate L2 metric matching the C++ RooUnfold QA",
            "quadrature": "sqrt(movement^2 + stat^2)",
            "choice_rule": "minimize quadrature over the decision window iterations 3..DISPLAY_MAX_ITER for each check; final common count is the largest local-minimum iteration required by the four checks",
            "display_note": "Slide display is limited to the decision window. Higher Au+Au 2D iterations are retained in CSV/manifest and show kCovToy covariance blow-up in low-content tail bins.",
        },
        "choices": choices,
        "cases": {
            key: {
                "case": asdict(payload["case"]),
                "selected_truth_pt_bins": payload["meta"][global_choice]["selected_truth_pt_bins"],
                "npho_unfolded_at_choice": payload["meta"][global_choice]["npho_unfolded"],
                "scan": payload["scans"],
            }
            for key, payload in all_scans.items()
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def write_speaker_script(global_choice: int) -> None:
    OUT_SCRIPT.write_text(
        f"""This slide checks the Bayesian unfolding iteration choice before showing the unfolded result.

The method follows the same logic as the PPG12/Ian iteration scan, but the slide now shows only the decision quantity. For each Bayes iteration, I compute the aggregate statistical scale S_n, the aggregate movement M_n from the previous iteration, and the quadrature score Q_n = sqrt(S_n squared plus M_n squared).

The metric is the aggregate L2 version used by the C++ RooUnfold QA: the statistical scale is sqrt(sum errors squared over sum contents squared), and the movement is sqrt(sum iteration differences squared over sum contents squared). That avoids letting a nearly empty AuAu tail bin decide the iteration count.

The plotted y-axis is Q_n divided by the local minimum for that check, so the four checks can be compared on one axis. The stars show the local minima. AuAu photon, AuAu 2D, and p+p photon all minimize at 3; p+p 2D minimizes at {global_choice}. Because the delivered spectra should use one common Bayes count, the common choice is {global_choice}.

The plotted slide stops at iteration {DISPLAY_MAX_ITER}. The full CSV and manifest still keep iterations 8-10; those late AuAu 2D points show covariance blow-up in low-content tail bins and are not part of the decision window.

The caveat is the same as the previous reconstructed-input slide: with the current histogram edges, this is a full-bin validation over the effective 16-35 GeV range. A final exact 15-35 GeV version needs a future pass with a true 15 GeV bin edge.
""",
        encoding="utf-8",
    )


def main() -> None:
    all_scans: Dict[str, Dict] = {}
    for case in selected_cases():
        print(f"[scan] {case.key}", flush=True)
        all_scans[case.key] = scan_case(case)
    global_choice, choices = choose_global_iteration(all_scans)
    write_csv(all_scans)
    draw_slide(all_scans, choices, global_choice)
    write_manifest(all_scans, choices, global_choice)
    write_speaker_script(global_choice)
    print(json.dumps({"png": str(OUT_PNG), "csv": str(OUT_CSV), "manifest": str(OUT_MANIFEST), "choice": global_choice}, indent=2))


if __name__ == "__main__":
    main()
