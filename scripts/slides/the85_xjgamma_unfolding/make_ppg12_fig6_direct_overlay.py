#!/usr/bin/env python3
"""Overlay PPG12 Figure 6 source efficiencies with THE-85 pp output.

This is intentionally a local slide-facing diagnostic.  It reads the direct
TEfficiency readback from the PPG12 Figure 6 source ROOT and compares it to the
current THE-85 nominal-only RecoilJets pp efficiency-stage output.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ROOT


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/efficiency_stage"
PPG12_CSV = OUT_DIR / "ppg12_bdt_nom_efficiency_tefficiency_readback_full_remote_clean.csv"
THE85_ROOT = (
    REPO
    / "InputFiles/the85_pp_fig6_nominalonly_20260630/"
    / "RecoilJets_photonjet5plus10plus20_nominalonly_MERGED.root"
)
OUT_PNG = OUT_DIR / "ppg12_fig6_direct_vs_the85_nominalonly_overlay_10to35.png"
OUT_CSV = OUT_DIR / "ppg12_fig6_direct_vs_the85_nominalonly_overlay_10to35.csv"
OUT_MANIFEST = OUT_DIR / "ppg12_fig6_direct_vs_the85_nominalonly_overlay_10to35_manifest.json"

PLOT_XMIN = 10.0
PLOT_XMAX = 35.0


STAGES = [
    {
        "key": "reco",
        "label": r"$\varepsilon_{\mathrm{reco}}$",
        "color": "#111111",
        "ylim": (0.84, 1.00),
    },
    {
        "key": "reco_id",
        "label": r"$\varepsilon_{\mathrm{reco}}\times\varepsilon_{\mathrm{ID}}$",
        "color": "#d62d91",
        "ylim": (0.46, 0.78),
    },
    {
        "key": "reco_id_iso",
        "label": r"$\varepsilon_{\mathrm{reco}}\times\varepsilon_{\mathrm{ID}}\times\varepsilon_{\mathrm{iso}}$",
        "color": "#2f9638",
        "ylim": (0.22, 0.58),
    },
]


def _read_ppg12_rows() -> dict[str, dict[float, dict[str, float]]]:
    by_stage: dict[str, dict[float, dict[str, float]]] = {}
    with PPG12_CSV.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            stage = row["stage"]
            mid = float(row["pt_mid"])
            lo = float(row["pt_low"])
            hi = float(row["pt_high"])
            if hi <= PLOT_XMIN or lo >= PLOT_XMAX:
                continue
            by_stage.setdefault(stage, {})[mid] = {
                "pt_low": lo,
                "pt_high": hi,
                "pt_mid": mid,
                "eff": float(row["efficiency"]),
                "err_low": float(row["err_low"]),
                "err_high": float(row["err_high"]),
                "source_root": row["source_root"],
                "source_object": row["source_object"],
            }
    return by_stage


def _product(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    eff = a["eff"] * b["eff"]
    # Conservative uncorrelated propagation for plotting only.  PPG12 Figure 6
    # states statistical uncertainties are smaller than marker size.
    rel_low = math.hypot(a["err_low"] / a["eff"], b["err_low"] / b["eff"])
    rel_high = math.hypot(a["err_high"] / a["eff"], b["err_high"] / b["eff"])
    return {
        "pt_low": a["pt_low"],
        "pt_high": a["pt_high"],
        "pt_mid": a["pt_mid"],
        "eff": eff,
        "err_low": eff * rel_low,
        "err_high": eff * rel_high,
    }


def ratio_error(num: float, den: float, enum: float, eden: float) -> float:
    if den <= 0.0 or num < 0.0:
        return math.nan
    return math.sqrt((enum / den) ** 2 + ((num * eden) / (den * den)) ** 2)


def load_ppg12_figure6() -> dict[str, list[dict[str, float]]]:
    raw = _read_ppg12_rows()
    out: dict[str, list[dict[str, float]]] = {"reco": [], "reco_id": [], "reco_id_iso": []}
    for mid in sorted(raw["reco"]):
        out["reco"].append(raw["reco"][mid])
        out["reco_id"].append(_product(raw["reco"][mid], raw["id"][mid]))
        # This matches ppg12codeGit/plotting/plot_efficiency.C, where the green
        # Figure 6 curve is read from eff_all_eta_0 rather than recomputed.
        out["reco_id_iso"].append(raw["all"][mid])
    return out


def load_the85_nominalonly() -> dict[str, list[dict[str, float]]]:
    out: dict[str, list[dict[str, float]]] = {"reco": [], "reco_id": [], "reco_id_iso": []}
    handle = ROOT.TFile.Open(str(THE85_ROOT), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open THE-85 ROOT: {THE85_ROOT}")
    sim = handle.Get("SIM")
    if not sim:
        raise RuntimeError(f"Missing SIM directory in {THE85_ROOT}")

    h_den = sim.Get("h_photonEffTruthDen_pTgamma_0")
    h_reco = sim.Get("h_photonEffPpg12Fig6Reco_pTgamma_0")
    h_iso = sim.Get("h_photonEffPpg12Fig6RecoIso_pTgamma_0")
    h_tight_iso = sim.Get("h_photonEffPpg12Fig6RecoTightIso_pTgamma_0")
    missing = [
        name
        for name, hist in [
            ("h_photonEffTruthDen_pTgamma_0", h_den),
            ("h_photonEffPpg12Fig6Reco_pTgamma_0", h_reco),
            ("h_photonEffPpg12Fig6RecoIso_pTgamma_0", h_iso),
            ("h_photonEffPpg12Fig6RecoTightIso_pTgamma_0", h_tight_iso),
        ]
        if not hist
    ]
    if missing:
        raise RuntimeError(f"Missing THE-85 histograms: {missing}")

    ppg12_rows = load_ppg12_figure6()
    for p in ppg12_rows["reco"]:
        mid = p["pt_mid"]
        ib_den = h_den.GetXaxis().FindBin(mid)
        den = float(h_den.GetBinContent(ib_den))
        den_err = float(h_den.GetBinError(ib_den))
        if den <= 0.0:
            continue

        ib_reco = h_reco.GetXaxis().FindBin(mid)
        ib_iso = h_iso.GetXaxis().FindBin(mid)
        ib_all = h_tight_iso.GetXaxis().FindBin(mid)
        reco_num = float(h_reco.GetBinContent(ib_reco))
        reco_err = float(h_reco.GetBinError(ib_reco))
        iso_num = float(h_iso.GetBinContent(ib_iso))
        iso_err = float(h_iso.GetBinError(ib_iso))
        all_num = float(h_tight_iso.GetBinContent(ib_all))
        all_err = float(h_tight_iso.GetBinError(ib_all))

        reco_eff = reco_num / den
        reco_eff_err = ratio_error(reco_num, den, reco_err, den_err)
        id_cond_eff = all_num / iso_num if iso_num > 0.0 else math.nan
        id_cond_err = ratio_error(all_num, iso_num, all_err, iso_err) if iso_num > 0.0 else math.nan
        reco_id_eff = reco_eff * id_cond_eff if math.isfinite(id_cond_eff) else math.nan
        reco_id_err = (
            math.sqrt((id_cond_eff * reco_eff_err) ** 2 + (reco_eff * id_cond_err) ** 2)
            if math.isfinite(id_cond_eff) and math.isfinite(id_cond_err)
            else math.nan
        )

        out["reco"].append(
            {
                "pt_low": p["pt_low"],
                "pt_high": p["pt_high"],
                "pt_mid": mid,
                "eff": reco_eff,
                "err_low": reco_eff_err,
                "err_high": reco_eff_err,
                "num": reco_num,
                "den": den,
            }
        )
        out["reco_id"].append(
            {
                "pt_low": p["pt_low"],
                "pt_high": p["pt_high"],
                "pt_mid": mid,
                "eff": reco_id_eff,
                "err_low": reco_id_err,
                "err_high": reco_id_err,
                "num": all_num,
                "den": iso_num,
            }
        )
        out["reco_id_iso"].append(
            {
                "pt_low": p["pt_low"],
                "pt_high": p["pt_high"],
                "pt_mid": mid,
                "eff": all_num / den,
                "err_low": ratio_error(all_num, den, all_err, den_err),
                "err_high": ratio_error(all_num, den, all_err, den_err),
                "num": all_num,
                "den": den,
            }
        )
    handle.Close()
    for rows in out.values():
        rows.sort(key=lambda r: r["pt_mid"])
    return out


def as_arrays(rows: list[dict[str, float]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.array([r["pt_mid"] for r in rows], dtype=float)
    y = np.array([r["eff"] for r in rows], dtype=float)
    ylo = np.array([r["err_low"] for r in rows], dtype=float)
    yhi = np.array([r["err_high"] for r in rows], dtype=float)
    return x, y, ylo, yhi


def write_comparison_csv(ppg12: dict[str, list[dict[str, float]]], the85: dict[str, list[dict[str, float]]]) -> None:
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "stage",
                "pt_low",
                "pt_high",
                "pt_mid",
                "ppg12_eff",
                "the85_nominalonly_eff",
                "the85_over_ppg12",
                "the85_num",
                "the85_den",
            ],
        )
        writer.writeheader()
        for stage in ["reco", "reco_id", "reco_id_iso"]:
            by_mid = {r["pt_mid"]: r for r in the85[stage]}
            for p in ppg12[stage]:
                t = by_mid.get(p["pt_mid"])
                if not t:
                    continue
                writer.writerow(
                    {
                        "stage": stage,
                        "pt_low": p["pt_low"],
                        "pt_high": p["pt_high"],
                        "pt_mid": p["pt_mid"],
                        "ppg12_eff": p["eff"],
                        "the85_nominalonly_eff": t["eff"],
                        "the85_over_ppg12": t["eff"] / p["eff"] if p["eff"] else "",
                        "the85_num": t.get("num", ""),
                        "the85_den": t.get("den", ""),
                    }
                )


def render_overlay(ppg12: dict[str, list[dict[str, float]]], the85: dict[str, list[dict[str, float]]]) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.05,
            "xtick.major.width": 1.05,
            "ytick.major.width": 1.05,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=200)
    gs = fig.add_gridspec(2, 3, height_ratios=[3.4, 1.0], hspace=0.08, wspace=0.18)

    main_axes = []
    ratio_axes = []
    for i, stage in enumerate(STAGES):
        ax = fig.add_subplot(gs[0, i])
        rax = fig.add_subplot(gs[1, i], sharex=ax)
        main_axes.append(ax)
        ratio_axes.append(rax)

        key = stage["key"]
        color = stage["color"]
        xp, yp, ypl, yph = as_arrays(ppg12[key])
        xt, yt, ytl, yth = as_arrays(the85[key])

        ax.errorbar(
            xp,
            yp,
            yerr=[ypl, yph],
            xerr=[xp - np.array([r["pt_low"] for r in ppg12[key]]), np.array([r["pt_high"] for r in ppg12[key]]) - xp],
            fmt="o",
            ms=7.5,
            mfc="white",
            mec=color,
            mew=2.0,
            ecolor=color,
            elinewidth=1.3,
            capsize=0,
            label="PPG12 Fig. 6 source",
            zorder=4,
        )
        ax.errorbar(
            xt,
            yt,
            yerr=[ytl, yth],
            xerr=[xt - np.array([r["pt_low"] for r in the85[key]]), np.array([r["pt_high"] for r in the85[key]]) - xt],
            fmt="s",
            ms=6.8,
            mfc=color,
            mec=color,
            mew=1.2,
            ecolor=color,
            alpha=0.82,
            elinewidth=1.2,
            capsize=0,
            label="THE-85 nominal-only",
            zorder=3,
        )
        ax.set_title(stage["label"], fontsize=20, pad=10)
        ax.set_ylim(*stage["ylim"])
        ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
        ax.grid(True, which="major", color="#d8dde7", linewidth=0.9, alpha=0.65)
        ax.tick_params(axis="both", which="major", labelsize=13, direction="in", top=True, right=True, length=6)
        ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=3)
        ax.minorticks_on()
        ax.tick_params(labelbottom=False)

        ratio = yt / yp
        ratio_err = ratio * np.hypot(yth / yt, yph / yp)
        rax.axhline(1.0, color="#5b6472", lw=1.0, ls="--", zorder=1)
        rax.errorbar(
            xt,
            ratio,
            yerr=ratio_err,
            xerr=[xt - np.array([r["pt_low"] for r in the85[key]]), np.array([r["pt_high"] for r in the85[key]]) - xt],
            fmt="s",
            ms=5.5,
            mfc=color,
            mec=color,
            ecolor=color,
            alpha=0.85,
            capsize=0,
            zorder=3,
        )
        rax.set_ylim(0.55, 1.45)
        rax.set_xlabel(r"$E_{\mathrm{T}}^{\gamma,\mathrm{truth}}$ [GeV]", fontsize=15)
        rax.tick_params(axis="both", which="major", labelsize=12, direction="in", top=True, right=True, length=5)
        rax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=2.5)
        rax.minorticks_on()
        rax.grid(True, which="major", color="#d8dde7", linewidth=0.8, alpha=0.55)
        if i == 0:
            ax.set_ylabel("Efficiency", fontsize=17)
            rax.set_ylabel("THE-85 /\nPPG12", fontsize=13)
        else:
            ax.set_yticklabels([])
            rax.set_yticklabels([])

    handles, labels = main_axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.888),
        ncol=2,
        frameon=False,
        fontsize=15,
        handlelength=2.2,
        columnspacing=2.4,
    )
    fig.suptitle(
        r"Direct overlay with PPG12 Figure 6 source data, $10<E_{\mathrm{T}}^{\gamma}<35$ GeV",
        fontsize=27,
        fontweight="bold",
        y=0.965,
    )
    fig.text(
        0.5,
        0.915,
        "Open markers: PPG12 TEfficiency readback from MC_efficiency_bdt_nom.root. "
        "Filled markers: current THE-85 RecoilJets nominal-only pp output.",
        ha="center",
        va="center",
        fontsize=15.5,
        color="#344054",
    )
    fig.text(
        0.5,
        0.045,
        "Caveat: THE-85 four-component pp merge is still blocked by failed double-MbdDigitization efficiency-hist audit; this overlay uses only passing nominal components.",
        ha="center",
        va="center",
        fontsize=12.8,
        color="#667085",
    )
    fig.subplots_adjust(left=0.07, right=0.985, top=0.80, bottom=0.12)
    fig.savefig(OUT_PNG, bbox_inches="tight")
    plt.close(fig)


def write_manifest() -> None:
    manifest = {
        "artifact": str(OUT_PNG),
        "purpose": "Direct overlay of PPG12 Figure 6 source efficiency data with THE-85 nominal-only pp RecoilJets output.",
        "ppg12": {
            "csv": str(PPG12_CSV),
            "source_root_from_readback": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root",
            "macro": "ppg12codeGit/plotting/plot_efficiency.C",
            "figure_canvas": "c6 / figures/eff_photon_bdt_nom.pdf",
            "objects": {
                "reco": "eff_reco_eta_0",
                "id_conditional": "eff_id_eta_0",
                "reco_id": "eff_reco_eta_0 * eff_id_eta_0",
                "reco_id_iso": "eff_all_eta_0",
            },
            "note": "Figure 6 green curve follows eff_all_eta_0 in the PPG12 macro; magenta is reco times conditional ID.",
        },
        "the85": {
            "root": str(THE85_ROOT),
            "source": "THE-85 RecoilJets nominal-only merged ROOT, not final four-component merge.",
            "blocked_final_merge_reason": "double MbdDigitization components drained but failed required truth-binned efficiency histogram audit.",
            "histograms": {
                "denominator": "SIM/h_photonEffTruthDen_pTgamma_0",
                "reco": "SIM/h_photonEffPpg12Fig6Reco_pTgamma_0",
                "reco_iso": "SIM/h_photonEffPpg12Fig6RecoIso_pTgamma_0",
                "reco_id_iso": "SIM/h_photonEffPpg12Fig6RecoTightIso_pTgamma_0",
                "reco_id": "(Fig6 reco / truth) * (Fig6 tight+iso / Fig6 iso)",
            },
        },
        "outputs": {
            "comparison_csv": str(OUT_CSV),
            "manifest": str(OUT_MANIFEST),
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    ppg12 = load_ppg12_figure6()
    the85 = load_the85_nominalonly()
    write_comparison_csv(ppg12, the85)
    render_overlay(ppg12, the85)
    write_manifest()
    print(OUT_PNG)
    print(OUT_CSV)
    print(OUT_MANIFEST)


if __name__ == "__main__":
    main()
