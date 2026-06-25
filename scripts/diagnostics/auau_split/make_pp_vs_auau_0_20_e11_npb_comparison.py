#!/usr/bin/env python3
"""Compare pp and AuAu 0-20 E11/E33 before/after NPB preselection.

The pp row uses the repaired PPG12/baseV3E table-QA products and inclusive
sample cache.  The AuAu row uses b009 photon-trigger data and older embedded
MC that passed the reconstruction-contract E11/E33 sanity check.  This is a
diagnostic plot only; it intentionally stops before tight-BDT interpretation.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


SCRIPT_PATH = Path(__file__).resolve()
REPO = next(p for p in SCRIPT_PATH.parents if (p / "AGENTS.md").exists())
if str(REPO / "scripts") not in sys.path:
    sys.path.append(str(REPO / "scripts"))

PP_CAMPAIGN = REPO / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611"
PP_ROOT_DIR = PP_CAMPAIGN / "merged_roots"
PP_INCLUSIVE_CACHE = PP_CAMPAIGN / "inclusive_sample_hist_cache/the42_ppg12_tableqa_v1_inclusive_sample_projectx_hists.json"
PLOTTER_PATH = REPO / "scripts/plotting/pp_currentian/make_the42_ppg12_tableqa_v1_tables.py"

PP_DATA_ROOT = PP_ROOT_DIR / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
PP_SIGNAL_ROOT = PP_ROOT_DIR / "RecoilJets_photonjet5plus10plus20_MERGED.root"
PP_TOPDIR = "Photon_4_GeV_plus_MBD_NS_geq_1"

AA_DATA_ROOT = REPO / (
    "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604/"
    "b009_merged_data_roots_20260609/"
    "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
    "nonTightAuAuBDTComplement_baseVariant.root"
)
AA_SIGNAL_ROOT = REPO / (
    "dataOutput/auau_widthstudy_pt1530_wp080/combinedSimOnlyEMBEDDED/"
    "preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant/"
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
AA_INCLUSIVE_ROOT = REPO / (
    "dataOutput/auau_widthstudy_pt1530_wp080/combinedSimOnlyEMBEDDED/"
    "preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant/"
    "embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root"
)
AA_DATA_DIR = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"

OUT_DIR = REPO / "dataOutput/auauTableQA/the42_b009_data_old_good_mc_preselection_check"
OUT_PNG = OUT_DIR / "pp_vs_auau_0_20_e11e33_npb_preselection_comparison.png"
OUT_JSON = OUT_DIR / "pp_vs_auau_0_20_e11e33_npb_preselection_comparison_manifest.json"

PP_PT_TOKEN = "1535"
AA_DATA_PT_TOKENS = ("15_18", "18_20", "20_22", "22_24", "24_26", "26_28", "28_30")
AA_MC_PT_TOKENS = ("15_16", "16_18", "18_20", "20_22", "22_24", "24_26", "26_30")
AA_DATA_CENTS = ("0_10", "10_20")
AA_MC_CENT = "0_20"


def load_plotter():
    spec = importlib.util.spec_from_file_location("tableqa_plotter", PLOTTER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {PLOTTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hist_metrics(y: np.ndarray, edges: np.ndarray) -> dict[str, float]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    total = float(np.sum(y))
    sumw2 = float(np.sum(y * y))
    return {
        "sumw": total,
        "low_edge_fraction_x_lt_0p05": float(np.sum(y[centers < 0.05]) / total) if total else 0.0,
        "max_bin_fraction": float(np.max(y) / total) if total else 0.0,
        "effective_entries": float(total * total / sumw2) if sumw2 > 0 else 0.0,
    }


def norm(y: np.ndarray, edges: np.ndarray, xlim=(0.0, 1.0)) -> tuple[np.ndarray, np.ndarray]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = (centers >= xlim[0]) & (centers <= xlim[1])
    out = y.astype(float).copy()
    total = float(np.sum(out[mask]))
    if total > 0:
        out /= total
    return out, mask


def pp_arrays(plotter, files: dict, cache: dict, stage: str) -> dict:
    xlim, rebin = plotter.ppg12_axis_settings("e11_to_e33")
    data = plotter.norm_arrays(plotter.get_hist(files["data"], PP_TOPDIR, "e11_to_e33", PP_PT_TOKEN, stage), rebin, xlim)
    sig = plotter.norm_arrays(plotter.get_hist(files["signal"], "SIM", "e11_to_e33", PP_PT_TOKEN, stage), rebin, xlim)
    inc_payload = plotter.cache_hist(cache, "current_ian_jet8to40", "e11_to_e33", PP_PT_TOKEN, stage)
    inc = plotter.norm_payload(inc_payload, rebin, xlim)
    return {
        "data": data,
        "signal": sig,
        "inclusive": inc,
        "stats": {
            "inclusive": plotter.payload_stats_after_transform(inc_payload, rebin, xlim),
        },
        "label": "pp, 15-35 GeV",
    }


def sum_aa_data(stage: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
    y = None
    edges = None
    used: list[str] = []
    with uproot.open(AA_DATA_ROOT) as f:
        d = f[AA_DATA_DIR]
        for pt in AA_DATA_PT_TOKENS:
            for cent in AA_DATA_CENTS:
                name = f"h_ss_e11e33_{stage}_pT_{pt}_cent_{cent}"
                if name not in d:
                    continue
                vals, e = d[name].to_numpy(flow=False)
                vals = np.asarray(vals, dtype=float)
                e = np.asarray(e, dtype=float)
                if y is None:
                    y, edges = vals.copy(), e.copy()
                else:
                    if y.shape != vals.shape or not np.allclose(edges, e):
                        raise ValueError(f"AuAu data bin mismatch {name}")
                    y += vals
                used.append(f"{AA_DATA_DIR}/{name}")
    if y is None or edges is None:
        raise KeyError(f"No AuAu data hists for stage {stage}")
    return y, edges, used


def sum_aa_mc(path: Path, stage: str, *, hist_stage: str | None = None) -> tuple[np.ndarray, np.ndarray, list[str]]:
    y = None
    edges = None
    used: list[str] = []
    with uproot.open(path) as f:
        d = f["SIM"]
        stage_key = hist_stage or stage
        for pt in AA_MC_PT_TOKENS:
            name = f"h_ss_e11e33_{stage_key}_pT_{pt}_cent_{AA_MC_CENT}"
            if name not in d:
                continue
            vals, e = d[name].to_numpy(flow=False)
            vals = np.asarray(vals, dtype=float)
            e = np.asarray(e, dtype=float)
            if y is None:
                y, edges = vals.copy(), e.copy()
            else:
                if y.shape != vals.shape or not np.allclose(edges, e):
                    raise ValueError(f"AuAu MC bin mismatch {name}")
                y += vals
            used.append(f"SIM/{name}")
    if y is None or edges is None:
        raise KeyError(f"No AuAu MC hists for stage {stage} in {path}")
    return y, edges, used


def aa_arrays(stage: str) -> dict:
    data_y, edges, data_used = sum_aa_data(stage)
    sig_stage = "pre_sig" if stage == "pre" else stage
    inc_stage = "pre_bkg" if stage == "pre" else stage
    sig_y, sig_edges, sig_used = sum_aa_mc(AA_SIGNAL_ROOT, stage, hist_stage=sig_stage)
    inc_y, inc_edges, inc_used = sum_aa_mc(AA_INCLUSIVE_ROOT, stage, hist_stage=inc_stage)
    if not (np.allclose(edges, sig_edges) and np.allclose(edges, inc_edges)):
        raise ValueError("AuAu data/MC binning mismatch")
    data_norm, mask = norm(data_y, edges)
    sig_norm, _ = norm(sig_y, edges)
    inc_norm, _ = norm(inc_y, edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return {
        "data": (centers[mask], data_norm[mask], np.sqrt(np.maximum(data_y[mask], 0)) / max(float(np.sum(data_y[mask])), 1.0)),
        "signal": (centers[mask], sig_norm[mask], None),
        "inclusive": (centers[mask], inc_norm[mask], None),
        "raw": {
            "data": hist_metrics(data_y, edges),
            "signal": hist_metrics(sig_y, edges),
            "inclusive": hist_metrics(inc_y, edges),
        },
        "used": {
            "data": data_used,
            "signal": sig_used,
            "inclusive": inc_used,
        },
        "label": "AuAu 0-20%, photon10, 15-30 GeV",
    }


def draw_pp_panel(ax, arr: dict, title: str) -> None:
    for key, color, label in [("signal", "#d62728", "Signal MC"), ("inclusive", "#1f77b4", "Inclusive MC")]:
        x, y, _ = arr[key]
        ax.step(x, y, where="mid", color=color, linewidth=2.0, label=label)
    x, y, e = arr["data"]
    ax.errorbar(x, y, yerr=e, fmt="o", ms=4.6, mfc="white", mec="black", mew=1.1, color="black", lw=0, label="Data")
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=8)
    ax.set_xlim(0, 1)
    ax.grid(True, axis="y", color="#d8dde5", alpha=0.8)


def draw_aa_panel(ax, arr: dict, title: str) -> None:
    for key, color, label in [("signal", "#d62728", "older signal MC"), ("inclusive", "#1f77b4", "older inclusive MC")]:
        x, y, _ = arr[key]
        ax.step(x, y, where="mid", color=color, linewidth=2.0, label=label)
    x, y, e = arr["data"]
    ax.errorbar(x, y, yerr=e, fmt="o", ms=4.6, mfc="white", mec="black", mew=1.1, color="black", lw=0, label="b009 data")
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=8)
    ax.set_xlim(0, 1)
    ax.grid(True, axis="y", color="#d8dde5", alpha=0.8)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })
    plotter = load_plotter()
    pp_files = {
        "data": plotter.open_root(PP_DATA_ROOT),
        "signal": plotter.open_root(PP_SIGNAL_ROOT),
    }
    pp_cache = plotter.load_inclusive_cache(PP_INCLUSIVE_CACHE, use_stitched_inclusive=False)

    panels = {
        ("pp", "cut0"): pp_arrays(plotter, pp_files, pp_cache, "cut0"),
        ("pp", "cut1"): pp_arrays(plotter, pp_files, pp_cache, "cut1"),
        ("aa", "inclusive"): aa_arrays("inclusive"),
        ("aa", "pre"): aa_arrays("pre"),
    }

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 7.8), sharex=True, constrained_layout=False)
    draw_pp_panel(axes[0, 0], panels[("pp", "cut0")], "pp before preselection")
    draw_pp_panel(axes[0, 1], panels[("pp", "cut1")], "pp after NPB preselection")
    draw_aa_panel(axes[1, 0], panels[("aa", "inclusive")], "AuAu 0-20% before preselection")
    draw_aa_panel(axes[1, 1], panels[("aa", "pre")], "AuAu 0-20% after full preselection")

    for ax in axes.ravel():
        ax.set_ylim(bottom=0)
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.10)
    for ax in axes[:, 0]:
        ax.set_ylabel("normalized counts", fontsize=12.5)
    for ax in axes[1, :]:
        ax.set_xlabel(r"$E_{11}/E_{33}$", fontsize=13)
    axes[0, 1].legend(loc="upper right", frameon=False, fontsize=10.5)
    axes[1, 1].legend(loc="upper right", frameon=False, fontsize=10.5)

    def lowfrac_text(payload: dict, is_pp: bool) -> str:
        if is_pp:
            # Normalized arrays after visible-range normalization; use first five visible bins as a shape proxy.
            parts = []
            for key, label in [("data", "data"), ("signal", "sig"), ("inclusive", "inc")]:
                x, y, _ = payload[key]
                parts.append(f"{label} x<0.05 {float(np.sum(y[x < 0.05])):.3f}")
            return "\n".join(parts)
        return "\n".join(
            f"{label} x<0.05 {payload['raw'][key]['low_edge_fraction_x_lt_0p05']:.3f}"
            for key, label in [("data", "data"), ("signal", "sig"), ("inclusive", "inc")]
        )

    for ax, key, is_pp in [
        (axes[0, 0], ("pp", "cut0"), True),
        (axes[0, 1], ("pp", "cut1"), True),
        (axes[1, 0], ("aa", "inclusive"), False),
        (axes[1, 1], ("aa", "pre"), False),
    ]:
        ax.text(
            0.035,
            0.955,
            lowfrac_text(panels[key], is_pp),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9.1,
            bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#c7c7c7", "alpha": 0.94},
        )

    fig.suptitle("E11/E33 response to NPB preselection: pp versus AuAu 0-20%", fontsize=18, fontweight="bold", y=0.986)
    fig.text(
        0.5,
        0.943,
        "pp uses repaired PPG12/baseV3E table-QA inclusive cache; AuAu uses photon10 data and older embedded MC with validated shower-shape contract.",
        ha="center",
        va="top",
        fontsize=11.5,
    )
    fig.text(
        0.5,
        0.018,
        "Diagnostic only. pp row uses 15-35 GeV; AuAu row uses the older-good-MC 15-30 GeV overlap and stops before tight BDT.",
        ha="center",
        va="bottom",
        fontsize=10.3,
        color="#4d4d4d",
    )
    fig.subplots_adjust(left=0.07, right=0.985, top=0.89, bottom=0.09, hspace=0.28, wspace=0.14)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)

    manifest = {
        "purpose": "diagnostic comparison of full-preselection effect in pp and AuAu 0-20 for E11/E33",
        "png": str(OUT_PNG),
        "pp": {
            "data_root": str(PP_DATA_ROOT),
            "signal_root": str(PP_SIGNAL_ROOT),
            "inclusive_cache": str(PP_INCLUSIVE_CACHE),
            "topdir": PP_TOPDIR,
            "pt_token": PP_PT_TOKEN,
        },
        "auau": {
            "data_root": str(AA_DATA_ROOT),
            "data_dir": AA_DATA_DIR,
            "signal_root": str(AA_SIGNAL_ROOT),
            "inclusive_root": str(AA_INCLUSIVE_ROOT),
            "data_pt_tokens": list(AA_DATA_PT_TOKENS),
            "mc_pt_tokens": list(AA_MC_PT_TOKENS),
            "centrality": "0-20%",
        },
        "metrics": {
            "pp_cut0_low_shape_proxy": lowfrac_text(panels[("pp", "cut0")], True),
            "pp_cut1_low_shape_proxy": lowfrac_text(panels[("pp", "cut1")], True),
            "auau_inclusive": panels[("aa", "inclusive")]["raw"],
            "auau_pre": panels[("aa", "pre")]["raw"],
        },
    }
    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n")
    print(OUT_PNG)
    print(OUT_JSON)


if __name__ == "__main__":
    main()
