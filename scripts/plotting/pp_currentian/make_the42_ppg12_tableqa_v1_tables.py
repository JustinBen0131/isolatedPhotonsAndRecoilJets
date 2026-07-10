#!/usr/bin/env python3
"""Render PPG12_TABLE_QA_V1 pp checkpoint tables from RecoilJets ROOTs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)


DEFAULT_CAMPAIGN_DIR = Path("dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611")
DEFAULT_ROOT_DIR = DEFAULT_CAMPAIGN_DIR / "merged_roots"
DEFAULT_INCLUSIVE_CACHE = (
    DEFAULT_CAMPAIGN_DIR
    / "inclusive_sample_hist_cache"
    / "the42_ppg12_tableqa_v1_inclusive_sample_projectx_hists.json"
)


VARS = [
    ("weta_cogx", "weta_cogx", r"$w_{\eta}^{\mathrm{cogX}}$"),
    ("wphi_cogx", "wphi_cogx", r"$w_{\phi}^{\mathrm{cogX}}$"),
    ("et1", "et1", r"$E_T^1/E_T^{\mathrm{cluster}}$"),
    ("et2", "et2", r"$E_T^2/E_T^{\mathrm{cluster}}$"),
    ("et3", "et3", r"$E_T^3/E_T^{\mathrm{cluster}}$"),
    ("et4", "et4", r"$E_T^4/E_T^{\mathrm{cluster}}$"),
    ("e11_to_e33", "e11_to_e33", r"$E_{11}/E_{33}$"),
    ("e17_to_e77", "e17_to_e77", r"$E_{17}/E_{77}$"),
    ("e32_to_e35", "e32_to_e35", r"$E_{32}/E_{35}$"),
    ("bdt", "bdt", "bdt"),
    ("npb_score", "npb_score", "npb_score"),
]


TABLES = [
    {
        "slug": "fig13_style_pt22_28_cut0",
        "title": "Fig. 13 style: 22 < E_T < 28 GeV, no preselection",
        "pt_token": "3",
        "pt_label": r"$22<E_T<28$ GeV",
        "cut": "cut0",
        "cut_label": "no preselection",
        "include_npb_template": True,
    },
    {
        "slug": "fig19_style_pt18_22_cut1",
        "title": "Fig. 19 style: 18 < E_T < 22 GeV, preselection + NPB",
        "pt_token": "2",
        "pt_label": r"$18<E_T<22$ GeV",
        "cut": "cut1",
        "cut_label": "preselection + NPB",
        "include_npb_template": False,
    },
    {
        "slug": "fig20_style_pt10_14_cut2",
        "title": "Fig. 20 style: 10 < E_T < 14 GeV, tight BDT ID",
        "pt_token": "0",
        "pt_label": r"$10<E_T<14$ GeV",
        "cut": "cut2",
        "cut_label": "tight BDT ID",
        "include_npb_template": False,
    },
    {
        "slug": "diagnostic_pt15_35_cut0_cut1_cut2_bdt",
        "title": "15 < E_T < 35 GeV diagnostic: BDT stage flow",
        "pt_token": "1535",
        "pt_label": r"$15<E_T<35$ GeV",
        "cut": "cut0",
        "cut_label": "stage comparison",
        "include_npb_template": False,
        "stage_flow": True,
    },
]


def open_root(path: Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {path}")
    return f


def hist_to_arrays(hist: ROOT.TH1) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nb = hist.GetNbinsX()
    x = np.array([hist.GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    y = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    e = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    return x, y, e


def ppg12_axis_settings(var: str) -> tuple[tuple[float, float], int]:
    """Match ppg12codeGit/plotting/plot_showershapes_variations.C."""
    xmin, xmax, nrebin = 0.0, 1.0, 4
    if var.startswith("w"):
        xmax = 2.0
    if var == "et4":
        xmax, nrebin = 0.3, 1
    if var == "e32_to_e35":
        xmin, xmax, nrebin = 0.4, 1.0, 1
    if var == "et1":
        xmin, xmax, nrebin = 0.3, 1.0, 1
    if var in {"bdt", "npb_score"}:
        xmin, xmax, nrebin = 0.0, 1.0, 2
    return (xmin, xmax), nrebin


def rebin_arrays(
    x: np.ndarray, y: np.ndarray, e: np.ndarray, factor: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Group adjacent display bins without changing the source ROOT histograms."""
    if factor <= 1 or len(x) < factor:
        return x, y, e
    n = len(x) // factor
    if n <= 0:
        return x, y, e
    trim = n * factor
    xr = x[:trim].reshape(n, factor)
    yr = y[:trim].reshape(n, factor)
    er = e[:trim].reshape(n, factor)
    return xr.mean(axis=1), yr.sum(axis=1), np.sqrt(np.sum(er * er, axis=1))


def range_arrays(
    x: np.ndarray, y: np.ndarray, e: np.ndarray, xlim: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = (x >= xlim[0]) & (x <= xlim[1])
    return x[mask], y[mask], e[mask]


def get_hist(
    f: ROOT.TFile, topdir: str, var: str, pt_token: str, cut: str, *, prefix: str = "h2d"
) -> ROOT.TH1 | None:
    name = f"{topdir}/{prefix}_{var}_eta0_pt{pt_token}_{cut}"
    obj = f.Get(name)
    if not obj or not obj.InheritsFrom("TH1"):
        return None
    return obj


def raw_project_arrays(
    hist: ROOT.TH1 | None, display_rebin: int, xlim: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    if hist is None:
        return None
    if hist.InheritsFrom("TH2"):
        clone = hist.Clone(f"{hist.GetName()}_clone_for_plot")
        clone.SetDirectory(0)
        if display_rebin > 1:
            clone.RebinX(display_rebin)
        clone.GetXaxis().SetRangeUser(xlim[0], xlim[1])
        proj = clone.ProjectionX(f"{hist.GetName()}_px_for_plot")
        if not proj:
            return None
        proj.SetDirectory(0)
        x, y, e = hist_to_arrays(proj)
    else:
        x, y, e = hist_to_arrays(hist)
        x, y, e = rebin_arrays(x, y, e, display_rebin)
    return range_arrays(x, y, e, xlim)


def norm_arrays(
    hist: ROOT.TH1 | None, display_rebin: int, xlim: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    arrays = raw_project_arrays(hist, display_rebin, xlim)
    if arrays is None:
        return None
    x, y, e = arrays
    total = float(np.sum(y))
    if not np.isfinite(total) or total <= 0:
        return None
    return x, y / total, e / total


def norm_payload(
    payload: dict | None, display_rebin: int, xlim: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    if not payload:
        return None
    x = np.asarray(payload.get("x", []), dtype=float)
    y = np.asarray(payload.get("y", []), dtype=float)
    e = np.asarray(payload.get("e", []), dtype=float)
    if x.size == 0 or y.size != x.size or e.size != x.size:
        return None
    x, y, e = rebin_arrays(x, y, e, display_rebin)
    x, y, e = range_arrays(x, y, e, xlim)
    total = float(np.sum(y))
    if not np.isfinite(total) or total <= 0:
        return None
    return x, y / total, e / total


def payload_stats_after_transform(
    payload: dict | None, display_rebin: int, xlim: tuple[float, float]
) -> dict | None:
    if payload is None:
        return None
    x = np.asarray(payload.get("x", []), dtype=float)
    y = np.asarray(payload.get("y", []), dtype=float)
    e = np.asarray(payload.get("e", []), dtype=float)
    if x.size == 0 or y.size != x.size or e.size != x.size:
        return None
    x, y, e = rebin_arrays(x, y, e, display_rebin)
    x, y, e = range_arrays(x, y, e, xlim)
    sumw = float(np.sum(y))
    sumw2 = float(np.sum(e * e))
    return {
        "sumw_visible": sumw,
        "sumw2_visible": sumw2,
        "neff_visible": float((sumw * sumw) / sumw2) if sumw2 > 0 else 0.0,
        "max_bin_fraction_visible": float(np.max(y) / sumw) if sumw > 0 and y.size else 0.0,
        "source_samples": payload.get("source_samples"),
    }


def load_inclusive_cache(path: Path | None, *, use_stitched_inclusive: bool) -> dict | None:
    if use_stitched_inclusive:
        return None
    if path is None:
        raise RuntimeError("inclusive cache path is required unless --use-stitched-inclusive is set")
    if not path.exists():
        raise FileNotFoundError(
            f"Inclusive sample histogram cache not found: {path}. "
            "Create it with scripts/diagnostics/pp_currentian/"
            "extract_the42_tableqa_inclusive_hists.py on SDCC, then rerun this plotter."
        )
    data = json.loads(path.read_text())
    if data.get("schema") not in {
        "THE42_PPG12_TABLE_QA_INCLUSIVE_SAMPLE_HIST_CACHE_V1",
        "THE42_PPG12_TABLE_QA_INCLUSIVE_SAMPLE_HIST_CACHE_V2",
    }:
        raise RuntimeError(f"Unexpected inclusive cache schema in {path}")
    return data


def cache_hist(cache: dict | None, sample_set: str, var: str, pt_token: str, cut: str) -> dict | None:
    if cache is None:
        return None
    prefix = "h2d" if cache.get("source") == "h2d-projectx" else "h1d"
    key = f"{prefix}_{var}_eta0_pt{pt_token}_{cut}"
    return cache.get("combined", {}).get(sample_set, {}).get(key)


def draw_shape(
    ax,
    arrays,
    *,
    label,
    color,
    marker=None,
    linestyle="-",
    linewidth=1.35,
    markersize=2.9,
    markerfacecolor=None,
):
    if arrays is None:
        return False
    x, y, e = arrays
    if marker:
        ax.errorbar(
            x,
            y,
            yerr=e,
            fmt=marker,
            color=color,
            markersize=markersize,
            markerfacecolor=color if markerfacecolor is None else markerfacecolor,
            markeredgewidth=1.0,
            linewidth=0.0,
            elinewidth=0.7,
            capsize=0,
            label=label,
            zorder=4,
        )
    else:
        ax.step(x, y, where="mid", color=color, linewidth=linewidth, linestyle=linestyle, label=label)
        ax.errorbar(x, y, yerr=e, fmt="none", ecolor=color, elinewidth=0.45, capsize=0, alpha=0.65, zorder=3)
    return True


def draw_residual(ax, data_arrays, inc_arrays, xlim: tuple[float, float]) -> None:
    ax.axhline(0.0, color="#777777", linewidth=0.55, linestyle=":")
    if data_arrays is not None and inc_arrays is not None and len(data_arrays[0]) == len(inc_arrays[0]):
        x, y_data, e_data = data_arrays
        _, y_inc, _ = inc_arrays
        diff = y_data - y_inc
        ax.errorbar(
            x,
            diff,
            yerr=e_data,
            fmt="o",
            color="black",
            markersize=1.9,
            linewidth=0.0,
            elinewidth=0.55,
            capsize=0,
        )
        finite = np.isfinite(diff + e_data)
        if np.any(finite):
            ymax = float(np.max(np.abs(diff[finite]) + e_data[finite]))
            if ymax > 0:
                ax.set_ylim(-1.3 * ymax, 1.3 * ymax)
    ax.set_xlim(*xlim)
    ax.tick_params(axis="both", labelsize=6.4, direction="in", top=True, right=True)
    ax.set_ylabel("Data - Incl. MC", fontsize=6.8)
    ax.yaxis.set_label_coords(-0.16, 0.5)


def scale_arrays(arrays, scale: float):
    if arrays is None:
        return None
    x, y, e = arrays
    return x, y * scale, e * scale


def chi2_summary(data_arrays, inc_arrays) -> dict | None:
    if data_arrays is None or inc_arrays is None or len(data_arrays[0]) != len(inc_arrays[0]):
        return None
    _, d, ed = data_arrays
    _, m, em = inc_arrays
    mask = (ed * ed + em * em > 0) & (m > 0)
    ndf = int(np.count_nonzero(mask))
    if ndf <= 1:
        return None
    chi2 = float(np.sum((d[mask] - m[mask]) ** 2 / (ed[mask] ** 2 + em[mask] ** 2)))
    ndf -= 1
    pvalue = float(ROOT.TMath.Prob(chi2, ndf))
    return {"chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf, "pvalue": pvalue}


def npb_tail_scale(files: dict[str, ROOT.TFile], table: dict) -> float:
    xlim = (0.0, 2.0)
    data = raw_project_arrays(
        get_hist(files["data"], "PPG12_scaledtrigger30", "weta_cogx", table["pt_token"], "cut0"),
        4,
        xlim,
    )
    npb = raw_project_arrays(
        get_hist(files["data"], "PPG12_scaledtrigger30", "weta_cogx", table["pt_token"], "cut4"),
        4,
        xlim,
    )
    if data is None or npb is None:
        return 1.0
    x_data, y_data, _ = data
    x_npb, y_npb, _ = npb
    if len(x_data) != len(x_npb):
        return 1.0
    mask = x_data >= 1.400001
    all_tot = float(np.sum(y_data))
    npb_tot = float(np.sum(y_npb))
    all_tail = float(np.sum(y_data[mask]))
    npb_tail = float(np.sum(y_npb[mask]))
    if all_tot <= 0 or npb_tot <= 0 or npb_tail <= 0:
        return 1.0
    return float((all_tail / all_tot) / (npb_tail / npb_tot))


def add_panel_annotation(
    ax,
    table: dict,
    chi2: dict | None,
    *,
    fontsize: float = 6.6,
    linespacing: float = 1.05,
    include_chi2: bool = True,
) -> None:
    cut_text = {
        "cut0": "w/o nbkg cut",
        "cut1": "w/ nbkg cut",
        "cut2": "w/ tight cut",
    }.get(table["cut"], table["cut_label"])
    lines = [
        r"$\bf{\it{sPHENIX}}$ Internal",
        r"$p$+$p$ $\sqrt{s}=200$ GeV",
        r"$|\eta^\gamma|<0.7$",
        f"{table['pt_label']},{cut_text}",
    ]
    if include_chi2 and chi2 is not None:
        lines.extend([
            rf"$\chi^2$/ndf = {chi2['chi2']:.1f}/{chi2['ndf']} = {chi2['chi2_ndf']:.2f}",
            f"p-value = {chi2['pvalue']:.4f}",
        ])
    ax.text(
        0.04,
        0.96,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fontsize,
        linespacing=linespacing,
    )


def draw_table(
    root_dir: Path,
    out_dir: Path,
    table: dict,
    display_rebin: int,
    *,
    inclusive_cache: dict | None,
    inclusive_cache_path: Path | None,
    inclusive_sample_set: str,
    use_stitched_inclusive: bool,
) -> dict:
    data_path = root_dir / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
    sig_path = root_dir / "RecoilJets_photonjet5plus10plus20_MERGED.root"
    inc_path = root_dir / "RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root"
    files = {
        "data": open_root(data_path),
        "signal_mc": open_root(sig_path),
    }
    if use_stitched_inclusive:
        files["inclusive_mc"] = open_root(inc_path)

    ncols = 3
    nrows = int(np.ceil(len(VARS) / ncols))
    fig = plt.figure(figsize=(10.2, 4.15 * nrows), constrained_layout=False)
    fig.patch.set_facecolor("white")
    outer = fig.add_gridspec(
        nrows, ncols, left=0.075, right=0.985, bottom=0.045, top=0.965, wspace=0.30, hspace=0.28
    )
    manifest = {
        "table": table,
        "display_rebin": "ppg12_per_variable" if display_rebin <= 0 else display_rebin,
        "projection_contract": (
            "Match ppg12codeGit/plotting/plot_showershapes_variations.C: "
            "use h2d_* ProjectionX, RebinX per variable, SetRangeUser visible window, then unit normalize."
        ),
        "npb_tail_match": None,
        "inputs": {k: str(v) for k, v in {
            "data": data_path,
            "signal_mc": sig_path,
        }.items()},
        "inclusive_mc": {
            "plot_source": (
                "stitched_final_root_debug_only"
                if use_stitched_inclusive
                else "sample_level_hist_cache_no_finalStitch_rescale"
            ),
            "sample_set": inclusive_sample_set if not use_stitched_inclusive else "jet5plus8plus12plus20plus30plus40_stitched",
            "cache_path": str(inclusive_cache_path) if inclusive_cache_path else None,
            "stitched_root_not_used": None if use_stitched_inclusive else str(inc_path),
            "note": (
                "Using final stitched inclusive ROOT is diagnostic only for PPG12_TABLE_QA_V1 because "
                "those hists are preweighted in RecoilJets and finalStitch applies an additional slice weight."
            ),
        },
        "panels": [],
    }
    tail_scale = 1.0
    if table.get("include_npb_template"):
        tail_scale = npb_tail_scale(files, table)
        manifest["npb_tail_match"] = {
            "scale_factor": tail_scale,
            "source_variable": "weta_cogx",
            "tail_condition": "weta_cogx >= 1.400001 after RebinX(4), matching PPG12 macro",
            "note": "PPG12-style display scaling for NPB-tagged data template; ROOT histograms are unchanged.",
        }

    for idx, (var, axis_name, label) in enumerate(VARS):
        xlim, ppg12_rebin = ppg12_axis_settings(var)
        use_rebin = display_rebin if display_rebin > 0 else ppg12_rebin
        sub = outer[idx // 3, idx % 3].subgridspec(2, 1, height_ratios=[4.0, 1.15], hspace=0.03)
        ax = fig.add_subplot(sub[0, 0])
        rax = fig.add_subplot(sub[1, 0], sharex=ax)
        ax.set_facecolor("white")
        rax.set_facecolor("white")

        cut = table["cut"]
        data = norm_arrays(
            get_hist(files["data"], "PPG12_scaledtrigger30", var, table["pt_token"], cut),
            use_rebin,
            xlim,
        )
        sig = norm_arrays(get_hist(files["signal_mc"], "SIM", var, table["pt_token"], cut), use_rebin, xlim)
        if use_stitched_inclusive:
            inc = norm_arrays(get_hist(files["inclusive_mc"], "SIM", var, table["pt_token"], cut), use_rebin, xlim)
            inc_raw_stats = None
        else:
            inc_payload = cache_hist(inclusive_cache, inclusive_sample_set, var, table["pt_token"], cut)
            inc = norm_payload(inc_payload, use_rebin, xlim)
            inc_raw_stats = payload_stats_after_transform(inc_payload, use_rebin, xlim)
        draw_shape(ax, data, label="Data", color="black", marker="o", markersize=2.7)
        draw_shape(ax, sig, label="Signal MC", color="red")
        draw_shape(ax, inc, label="Inclusive MC", color="blue")
        if table.get("include_npb_template"):
            npb = norm_arrays(
                get_hist(files["data"], "PPG12_scaledtrigger30", var, table["pt_token"], "cut4"),
                use_rebin,
                xlim,
            )
            npb = scale_arrays(npb, tail_scale)
            draw_shape(ax, npb, label="NPB-tagged data", color="#238b1e")

        chi2 = chi2_summary(data, inc)
        add_panel_annotation(ax, table, chi2)
        if idx == 0:
            ax.legend(loc="upper right", frameon=False, fontsize=7.1, handlelength=1.5)

        ax.set_xlim(*xlim)
        ax.set_ylim(bottom=0)
        ax.tick_params(axis="both", labelsize=7.4, direction="in", top=True, right=True, labelbottom=False)
        ax.set_ylabel("normalized counts", fontsize=7.8)
        draw_residual(rax, data, inc, xlim)
        rax.set_xlabel(axis_name, fontsize=7.4)
        rax.tick_params(axis="x", labelsize=7.2)
        manifest["panels"].append({
            "variable": var,
            "axis_name": axis_name,
            "label": label,
            "xlim": xlim,
            "ppg12_rebin": ppg12_rebin,
            "applied_rebin": use_rebin,
            "chi2": chi2,
            "inclusive_mc_raw_stats": inc_raw_stats,
        })

    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"the42_ppg12_tableqa_v1_clean_{table['slug']}.png"
    fig.savefig(png, dpi=220)
    plt.close(fig)
    manifest["outputs"] = {"png": str(png)}
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-dir", type=Path, default=DEFAULT_ROOT_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_CAMPAIGN_DIR / "tables_repaired_inclusive")
    parser.add_argument("--inclusive-cache", type=Path, default=DEFAULT_INCLUSIVE_CACHE)
    parser.add_argument("--inclusive-sample-set", default="current_ian_jet8to40")
    parser.add_argument(
        "--use-stitched-inclusive",
        action="store_true",
        help="Debug-only fallback to the final stitched inclusive ROOT. Do not use for table-QA physics comparisons.",
    )
    parser.add_argument(
        "--display-rebin",
        type=int,
        default=0,
        help="0 uses PPG12 per-variable rebinning; positive values force one adjacent-bin grouping.",
    )
    args = parser.parse_args()

    table_specs = [table for table in TABLES if not table.get("stage_flow")]
    inclusive_cache = load_inclusive_cache(args.inclusive_cache, use_stitched_inclusive=args.use_stitched_inclusive)
    manifests = [
        draw_table(
            args.root_dir,
            args.out_dir,
            table,
            args.display_rebin,
            inclusive_cache=inclusive_cache,
            inclusive_cache_path=args.inclusive_cache,
            inclusive_sample_set=args.inclusive_sample_set,
            use_stitched_inclusive=args.use_stitched_inclusive,
        )
        for table in table_specs
    ]
    manifest_path = args.out_dir / "the42_ppg12_tableqa_v1_tables_manifest.json"
    manifest_path.write_text(json.dumps({"tables": manifests}, indent=2) + "\n")
    print(f"WROTE {manifest_path}")
    for m in manifests:
        print(m["outputs"]["png"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
