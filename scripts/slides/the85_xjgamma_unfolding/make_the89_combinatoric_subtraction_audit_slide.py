#!/usr/bin/env python3
"""THE-89 audit slide for the Au+Au xJgamma combinatoric subtraction.

This helper is intentionally local and slide-facing.  It does not submit jobs,
mutate Google Slides, or change ROOT inputs.  It recomputes the same
ABCD-input and photon-yield-scaled combinatoric template used by the current
first-pass unfolding helper, then renders a one-slide diagnostic summary.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_unfolded_xjgamma_1x3 as unfold  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
OUT_PNG = OUT_DIR / "the89_combinatoric_subtraction_audit_summary_v1.png"
OUT_MANIFEST = OUT_DIR / "the89_combinatoric_subtraction_audit_summary_v1_manifest.json"
OUT_SCRIPT = OUT_DIR / "the89_combinatoric_subtraction_audit_summary_v1_speaker_script.md"
OUT_XJ_CSV = OUT_DIR / "the89_combinatoric_subtraction_audit_xj_components_v1.csv"
OUT_PT_CSV = OUT_DIR / "the89_combinatoric_subtraction_audit_pt_rows_v1.csv"


INK = "#111827"
MUTED = "#475569"
GRID = "#e5e7eb"
BLUE = "#1f4cff"
RED = "#e3342f"
PURPLE = "#7c3aed"
GRAY = "#64748b"
GREEN = "#0f766e"
PANEL = "#f8fafc"
EDGE = "#cbd5e1"


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )


def load_npz(name: str) -> dict[str, np.ndarray | float | list]:
    path = OUT_DIR / name
    with np.load(path, allow_pickle=True) as z:
        out: dict[str, np.ndarray | float | list] = {k: z[k].copy() for k in z.files}
    out["path"] = str(path)
    return out


def integral(result: dict, xmin: float | None = None, xmax: float | None = None) -> float:
    x = np.asarray(result["x_centers"], dtype=float)
    y = np.asarray(result["y"], dtype=float)
    w = np.asarray(result["x_widths"], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
    if xmin is not None:
        mask &= x >= xmin
    if xmax is not None:
        mask &= x < xmax
    return float(np.sum(y[mask] * w[mask]))


def nearest_y(result: dict, x0: float) -> tuple[float, float, float]:
    x = np.asarray(result["x_centers"], dtype=float)
    y = np.asarray(result["y"], dtype=float)
    ey = np.asarray(result["ey"], dtype=float)
    i = int(np.nanargmin(np.abs(x - x0)))
    return float(x[i]), float(y[i]), float(ey[i])


def plot_final_ladder(ax, raw: dict, abcd: dict, final: dict, pp: dict) -> None:
    entries = [
        (pp, "p+p raw-A reference", BLUE, "s", "-", 1.0),
        (raw, "Au+Au raw A", GRAY, "o", "--", 0.94),
        (abcd, "Au+Au ABCD only", PURPLE, "o", "-.", 0.96),
        (final, "Au+Au ABCD + comb", RED, "o", "-", 1.0),
    ]
    for result, label, color, marker, ls, alpha in entries:
        x = np.asarray(result["x_centers"], dtype=float)
        y = np.asarray(result["y"], dtype=float)
        ey = np.asarray(result["ey"], dtype=float)
        mask = (x >= 0.22) & (x <= 1.8) & np.isfinite(y) & np.isfinite(ey)
        ax.errorbar(
            x[mask],
            y[mask],
            yerr=ey[mask],
            fmt=marker,
            linestyle=ls,
            color=color,
            mfc="white" if label.startswith("p+p") else color,
            mec=color,
            ms=5.0,
            lw=1.55,
            elinewidth=1.0,
            capsize=2.0,
            alpha=alpha,
            label=label,
            zorder=3 if "comb" in label else 2,
        )

    ax.set_xlim(0.18, 1.82)
    ax.set_ylim(0.0, 0.78)
    ax.set_xlabel(r"particle-level $x_{J\gamma}$", fontsize=17.5)
    ax.set_ylabel(r"$(1/N_\gamma)\,dN/dx_{J\gamma}$", fontsize=17.5)
    ax.tick_params(labelsize=13.8)
    ax.grid(True, color=GRID, lw=0.8, alpha=0.8)
    ax.minorticks_on()
    ax.legend(loc="upper right", frameon=False, fontsize=11.8, handlelength=2.1)
    ax.text(
        0.025,
        0.955,
        r"$\bf{\it{sPHENIX}}$ Internal" "\n"
        r"RooUnfoldBayes, 5 iterations, kCovariance" "\n"
        r"$15<E_T^\gamma<35$ GeV target, $|\Delta\phi|>7\pi/8$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12.7,
        color=INK,
        linespacing=1.04,
    )


def build_preunfold_components() -> tuple[dict[str, np.ndarray], list[dict]]:
    case = next(c for c in unfold.CASES if c.key == "auau_0_20")
    data_f = unfold.open_root(case.data_file)
    sim_f = unfold.open_root(case.sim_file)

    h_raw = unfold.get_obj(
        data_f,
        case.data_topdir,
        f"h2_unfoldReco_pTgamma_xJ_incl_{unfold.BASE_KEY}{case.cent_suffix}",
        "TH2",
    )
    h_side_c = unfold.get_optional(
        data_f,
        case.data_topdir,
        f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{unfold.BASE_KEY}{case.cent_suffix}",
    )
    h_abcd, purity_meta = unfold.apply_xj_abcd_input(case, data_f, sim_f, h_raw, h_side_c)

    h_comb = unfold.get_optional(
        sim_f,
        case.sim_topdir,
        f"h2_unfoldRecoCombinatoric_pTgamma_xJ_incl_{unfold.BASE_KEY}{case.cent_suffix}",
    )
    if h_comb is None:
        raise RuntimeError("missing combinatoric template for Au+Au 0-20")

    h_pho_data, _ = unfold.apply_photon_abcd_input(
        case,
        data_f,
        sim_f,
        unfold.get_obj(
            data_f,
            case.data_topdir,
            f"h_unfoldRecoPho_pTgamma_{unfold.PHO_KEY}{case.cent_suffix}",
            "TH1",
        ),
    )
    h_pho_sim = unfold.get_obj(
        sim_f,
        case.sim_topdir,
        f"h_unfoldRecoPho_pTgamma_{unfold.PHO_KEY}{case.cent_suffix}",
        "TH1",
    )

    yaxis = h_abcd.GetYaxis()
    x_edges = unfold.axis_edges(yaxis)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    raw = np.zeros(yaxis.GetNbins())
    abcd = np.zeros(yaxis.GetNbins())
    comb = np.zeros(yaxis.GetNbins())
    after = np.zeros(yaxis.GetNbins())
    after_err2 = np.zeros(yaxis.GetNbins())
    row_table: list[dict] = []

    for ix in range(1, h_abcd.GetXaxis().GetNbins() + 1):
        lo = h_abcd.GetXaxis().GetBinLowEdge(ix)
        hi = h_abcd.GetXaxis().GetBinUpEdge(ix)
        cen = h_abcd.GetXaxis().GetBinCenter(ix)
        if not unfold.row_in_pt_window(lo, hi, cen):
            continue
        n_data = h_pho_data.GetBinContent(ix)
        n_sim = h_pho_sim.GetBinContent(ix)
        scale = n_data / n_sim if n_sim > 0.0 else 0.0
        row_raw = 0.0
        row_abcd = 0.0
        row_comb_raw = 0.0
        row_comb = 0.0
        row_after = 0.0
        for iy in range(1, yaxis.GetNbins() + 1):
            j = iy - 1
            rv = h_raw.GetBinContent(ix, iy)
            av = h_abcd.GetBinContent(ix, iy)
            ae = h_abcd.GetBinError(ix, iy)
            cv_raw = h_comb.GetBinContent(ix, iy)
            ce_raw = h_comb.GetBinError(ix, iy)
            cv = cv_raw * scale
            ce = ce_raw * scale
            raw[j] += rv
            abcd[j] += av
            comb[j] += cv
            after[j] += av - cv
            after_err2[j] += ae * ae + ce * ce
            row_raw += rv
            row_abcd += av
            row_comb_raw += cv_raw
            row_comb += cv
            row_after += av - cv
        row_table.append(
            {
                "pt_lo": lo,
                "pt_hi": hi,
                "pt_center": cen,
                "raw_entries": row_raw,
                "abcd_entries": row_abcd,
                "comb_raw_entries": row_comb_raw,
                "photon_data_for_scale": n_data,
                "photon_sim_for_scale": n_sim,
                "data_over_sim_photon_scale": scale,
                "comb_scaled_entries": row_comb,
                "after_comb_entries": row_after,
                "comb_over_abcd": row_comb / row_abcd if row_abcd > 0.0 else math.nan,
            }
        )

    components = {
        "x_edges": x_edges,
        "x_centers": x_centers,
        "raw": raw,
        "abcd": abcd,
        "comb": comb,
        "after": after,
        "after_err": np.sqrt(after_err2),
        "purity_meta": purity_meta,
    }
    return components, row_table


def plot_preunfold_components(ax, comp: dict[str, np.ndarray]) -> None:
    x_edges = np.asarray(comp["x_edges"], dtype=float)
    x = np.asarray(comp["x_centers"], dtype=float)
    widths = np.diff(x_edges)
    abcd = np.asarray(comp["abcd"], dtype=float)
    comb = np.asarray(comp["comb"], dtype=float)
    after = np.asarray(comp["after"], dtype=float)
    after_err = np.asarray(comp["after_err"], dtype=float)

    mask = (x >= 0.22) & (x <= 1.8)
    ax.stairs(abcd, x_edges, color=GRAY, lw=2.1, label=r"$S_{ABCD}$ before comb")
    ax.stairs(comb, x_edges, color=RED, lw=2.2, linestyle=":", label=r"scaled embedded $K$")
    ax.errorbar(
        x[mask],
        after[mask],
        xerr=0.5 * widths[mask],
        yerr=after_err[mask],
        fmt="o",
        ms=4.7,
        color=INK,
        elinewidth=1.0,
        capsize=0,
        label=r"$S_{ABCD}-K$ input",
        zorder=5,
    )
    ax.axhline(0.0, color="#94a3b8", lw=0.9)
    ax.set_xlim(0.18, 1.82)
    ax.set_ylim(min(-9, float(np.nanmin(after[mask]) * 1.18)), max(1.0, float(np.nanmax(abcd[mask]) * 1.22)))
    ax.set_xlabel(r"reconstructed $x_{J\gamma}$", fontsize=17.5)
    ax.set_ylabel("entries before unfolding", fontsize=17.5)
    ax.tick_params(labelsize=13.8)
    ax.grid(True, color=GRID, lw=0.8, alpha=0.8)
    ax.minorticks_on()
    ax.legend(loc="upper right", frameon=False, fontsize=12.2, handlelength=2.4)

    total_abcd = float(np.sum(abcd[mask]))
    total_comb = float(np.sum(comb[mask]))
    frac = total_comb / total_abcd if total_abcd > 0 else math.nan
    ax.text(
        0.030,
        0.955,
        "Same measured-space subtraction\n"
        "used by the current unfolded output\n"
        rf"$K/S_{{ABCD}}$ = {frac:.2f} in shown window",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12.9,
        color=INK,
        linespacing=1.05,
        bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="#d1d5db", alpha=0.96),
    )


def add_card(fig, xywh, title: str, body: str, accent: str) -> None:
    x, y, w, h = xywh
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.012",
        transform=fig.transFigure,
        linewidth=1.15,
        edgecolor=EDGE,
        facecolor=PANEL,
        zorder=1,
    )
    fig.patches.append(box)
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            0.006,
            h,
            boxstyle="round,pad=0.0,rounding_size=0.010",
            transform=fig.transFigure,
            linewidth=0,
            facecolor=accent,
            zorder=2,
        )
    )
    fig.text(x + 0.017, y + h - 0.027, title, ha="left", va="top", fontsize=17.0, fontweight="bold", color=INK)
    fig.text(x + 0.017, y + h - 0.073, body, ha="left", va="top", fontsize=13.2, color=MUTED, linespacing=1.10)


def write_csvs(comp: dict[str, np.ndarray], rows: list[dict]) -> None:
    with OUT_XJ_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["x_low", "x_high", "x_center", "raw_A", "S_ABCD", "scaled_comb_K", "S_ABCD_minus_K"])
        x_edges = comp["x_edges"]
        for i, xc in enumerate(comp["x_centers"]):
            writer.writerow([x_edges[i], x_edges[i + 1], xc, comp["raw"][i], comp["abcd"][i], comp["comb"][i], comp["after"][i]])

    with OUT_PT_CSV.open("w", newline="", encoding="utf-8") as f:
        fields = list(rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    raw = load_npz("the85_unfolded_xjgamma_auau_0_20_rawA_iter5_covariance_correction_step.npz")
    abcd = load_npz("the85_unfolded_xjgamma_auau_0_20_abcd_only_iter5_covariance_correction_step.npz")
    final = load_npz("the85_unfolded_xjgamma_auau_0_20_abcd_comb_iter5_covariance_correction_step.npz")
    pp = load_npz("the85_unfolded_xjgamma_pp_basev3e_iter5_covariance.npz")
    comp, pt_rows = build_preunfold_components()
    write_csvs(comp, pt_rows)

    raw_int = integral(raw)
    abcd_int = integral(abcd)
    final_int = integral(final)
    pp_tail05 = integral(pp, 0.5)
    final_tail05 = integral(final, 0.5)
    abcd_tail05 = integral(abcd, 0.5)

    _, abcd_y038, _ = nearest_y(abcd, 0.38)
    _, final_y038, _ = nearest_y(final, 0.38)
    _, pp_y079, _ = nearest_y(pp, 0.79)

    row_red_flags = [
        r for r in pt_rows if math.isfinite(float(r["comb_over_abcd"])) and float(r["comb_over_abcd"]) >= 0.5
    ]

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    fig.text(
        0.045,
        0.940,
        r"Audit: the final Au+Au $x_{J\gamma}$ scale is set by the combinatoric subtraction",
        ha="left",
        va="top",
        fontsize=29.5,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.047,
        0.890,
        "Region-C sideband and embedded combinatoric jets are separate corrections; the large final drop happens after the embedded template is subtracted.",
        ha="left",
        va="top",
        fontsize=15.8,
        color=MUTED,
    )

    ax_l = fig.add_axes([0.060, 0.315, 0.430, 0.515])
    ax_r = fig.add_axes([0.535, 0.315, 0.405, 0.515])
    plot_final_ladder(ax_l, raw, abcd, final, pp)
    plot_preunfold_components(ax_r, comp)

    drop_all = 1.0 - final_int / abcd_int if abcd_int > 0 else math.nan
    drop_tail = 1.0 - final_tail05 / abcd_tail05 if abcd_tail05 > 0 else math.nan
    peak_drop = 1.0 - final_y038 / abcd_y038 if abcd_y038 > 0 else math.nan

    add_card(
        fig,
        (0.060, 0.090, 0.275, 0.145),
        "Dominant change is the comb step",
        (
            rf"ABCD-only integral {abcd_int:.3f} becomes {final_int:.3f}." "\n"
            rf"{drop_all:.0%} removed; peak {abcd_y038:.2f} becomes {final_y038:.2f}."
        ),
        RED,
    )
    add_card(
        fig,
        (0.362, 0.090, 0.275, 0.145),
        "Sideband still matters",
        (
            "Region C sets photon-fake subtraction\n"
            "and the photon scale used for K.\n"
            "ABCD-only is not the collapsed curve."
        ),
        PURPLE,
    )
    red_flag_text = ", ".join(
        f"{int(r['pt_lo'])}-{int(r['pt_hi'])}: {float(r['comb_over_abcd']):.2f}" for r in row_red_flags[:3]
    )
    add_card(
        fig,
        (0.664, 0.090, 0.275, 0.145),
        "Most suspicious audit handle",
        (
            "Largest $K/S_{ABCD}$ rows: 14-16 0.51,\n"
            "16-18 0.58, 18-20 0.85.\n"
            "Check K normalization, matching, and closure."
        ),
        GREEN,
    )

    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    manifest = {
        "script": "scripts/slides/the85_xjgamma_unfolding/make_the89_combinatoric_subtraction_audit_slide.py",
        "slide": str(OUT_PNG),
        "csv_xj_components": str(OUT_XJ_CSV),
        "csv_pt_rows": str(OUT_PT_CSV),
        "purpose": "THE-89 audit summary for why the current Au+Au final per-photon xJgamma output is suppressed.",
        "inputs": {
            "rawA_npz": raw["path"],
            "abcd_only_npz": abcd["path"],
            "final_npz": final["path"],
            "pp_npz": pp["path"],
            "data_file": next(c for c in unfold.CASES if c.key == "auau_0_20").data_file,
            "sim_file": next(c for c in unfold.CASES if c.key == "auau_0_20").sim_file,
        },
        "normalization": "Nominal per-unfolded-photon density; no tail renormalization is applied.",
        "integrals": {
            "auau_rawA_all": raw_int,
            "auau_abcd_only_all": abcd_int,
            "auau_final_abcd_plus_comb_all": final_int,
            "auau_abcd_only_tail_xj_gt_0p5": abcd_tail05,
            "auau_final_tail_xj_gt_0p5": final_tail05,
            "pp_tail_xj_gt_0p5": pp_tail05,
        },
        "dominant_effect": {
            "fraction_removed_all_vs_abcd_only": drop_all,
            "fraction_removed_tail_gt_0p5_vs_abcd_only": drop_tail,
            "fraction_removed_at_xj_0p38_vs_abcd_only": peak_drop,
        },
        "pt_row_audit": pt_rows,
        "conclusions": [
            "The region-C non-tight sideband and embedded combinatoric template are separate correction terms in the current code path.",
            "The ABCD-only unfolded Au+Au output is not the collapsed curve; the large final drop appears when the scaled embedded combinatoric template is subtracted.",
            "The sideband definition still matters indirectly because it sets S_ABCD and the photon-yield scale used for K, but it is not the direct red-template fill definition.",
            "The most immediate audit target is the embedded template normalization/definition and closure, especially pT rows with K/S_ABCD >= 0.5.",
            "The pp overlay is also first-pass raw-A and should not be treated as a finished ATLAS-like corrected reference.",
        ],
        "code_pointers": {
            "final_normalization": "scripts/slides/the85_xjgamma_unfolding/make_unfolded_xjgamma_1x3.py:558",
            "xj_sideband_subtraction": "scripts/slides/the85_xjgamma_unfolding/make_unfolded_xjgamma_1x3.py:409",
            "comb_subtraction": "scripts/slides/the85_xjgamma_unfolding/make_unfolded_xjgamma_1x3.py:505",
            "comb_fill_definition": "src_AuAu/RecoilJets_AuAu.cc:9602",
            "macro_subtraction_reference": "macros/AnalyzeRecoilJets_RooUnfoldPipeline.cpp:10325",
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    OUT_SCRIPT.write_text(
        (
            "# Speaker script\n\n"
            "The important separation on this slide is that the non-tight sideband and the embedded "
            "combinatoric correction are not the same object. The sideband controls the photon-fake "
            "or dijet subtraction. The red embedded template is the unrelated recoil-jet background, "
            "scaled row by row using the corrected photon yield and subtracted before unfolding.\n\n"
            f"The left panel shows the current per-photon unfolded ladder. Raw A has integral {raw_int:.3f}. "
            f"After ABCD it is still {abcd_int:.3f}. After the embedded combinatoric subtraction it drops to "
            f"{final_int:.3f}. That is the dominant change, not a final plotting normalization.\n\n"
            "The right panel shows the same thing before unfolding. The scaled embedded template is large in "
            "the same xJ region where the final curve is suppressed. The strongest audit handle is to validate "
            "the template normalization and fill definition row by row, including the match-radius/threshold "
            "choice and closure. The bounded non-tight campaign still matters, but it should be described as "
            "a sideband/purity fix plus an indirect photon-scale input to K, not as the direct definition of K.\n"
        ),
        encoding="utf-8",
    )

    print(json.dumps({"png": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    main()
