#!/usr/bin/env python3
"""Rebuild the old PPG12 isolation-stack diagnostic from current RecoilJets ROOTs."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]

DEFAULT_PP_DATA = REPO / "InputFiles/pp24/ppg12_photon_yield_v1_data_20260620/pp/RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
DEFAULT_PP_SIGNAL = REPO / "InputFiles/pp24/ppg12_photon_yield_v1_signal_sim_mciso_20260621/sim/RecoilJets_photonjet5plus10plus20_MERGED.root"
DEFAULT_AUAU_DATA = REPO / "InputFiles/the69_default_auau_physicsqa/RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root"
DEFAULT_AUAU_SIGNAL = REPO / "InputFiles/the69_default_auau_physicsqa/RecoilJets_embeddedPhoton12plus20_MERGED.root"


@dataclass(frozen=True)
class SampleSpec:
    tag: str
    collision_label: str
    pt_label: str
    pt_bins: tuple[str, ...]
    data_file: Path
    signal_file: Path
    data_dir: str
    signal_dir: str
    tight_pattern: str
    nontight_pattern: str
    signal_pattern: str
    output_name: str
    sideband_note: str
    include_in_summary: bool = True


def read_hist(fh: uproot.ReadOnlyDirectory, key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if key not in fh:
        raise KeyError(f"missing histogram: {key}")
    hist = fh[key]
    values, edges = hist.to_numpy(flow=False)
    variances = hist.variances(flow=False)
    if variances is None:
        variances = np.clip(values, 0.0, None)
    return values.astype(float), variances.astype(float), edges.astype(float)


def add_hists(path: Path, directory: str, pattern: str, pt_bins: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    keys = [f"{directory}/{pattern.format(pt=pt)}" for pt in pt_bins]
    with uproot.open(path) as fh:
        values_sum = None
        variances_sum = None
        edges_ref = None
        for key in keys:
            values, variances, edges = read_hist(fh, key)
            if values_sum is None:
                values_sum = np.zeros_like(values, dtype=float)
                variances_sum = np.zeros_like(variances, dtype=float)
                edges_ref = edges
            elif len(edges) != len(edges_ref) or np.max(np.abs(edges - edges_ref)) > 1e-9:
                raise ValueError(f"incompatible binning for {key}")
            values_sum += values
            variances_sum += variances
    assert values_sum is not None and variances_sum is not None and edges_ref is not None
    return values_sum, variances_sum, edges_ref, keys


def ppg12_variable_rebin(
    values: np.ndarray,
    variances: np.ndarray,
    edges: np.ndarray,
    *,
    low_region_group: int,
    high_region_group: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    new_edges = [float(edges[0])]
    rebinned_values: list[float] = []
    rebinned_variances: list[float] = []
    i = 0
    while i < len(values):
        low = edges[i]
        group = low_region_group if low < 2.5 else high_region_group
        j = min(i + group, len(values))
        rebinned_values.append(float(np.sum(values[i:j])))
        rebinned_variances.append(float(np.sum(variances[i:j])))
        new_edges.append(float(edges[j]))
        i = j
    return np.array(rebinned_values), np.array(rebinned_variances), np.array(new_edges)


def bin_width_scale(values: np.ndarray, variances: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    widths = np.diff(edges)
    return values / widths, variances / np.square(widths)


def tail_integral(values: np.ndarray, edges: np.ndarray, low: float) -> float:
    first = int(np.searchsorted(edges, low, side="right") - 1)
    first = max(0, min(first, len(values) - 1))
    return float(np.sum(values[first:]))


def old_ppg12_stack(spec: SampleSpec, tail_low: float, low_region_group: int, high_region_group: int) -> dict[str, object]:
    tight, tight_var, edges, tight_keys = add_hists(spec.data_file, spec.data_dir, spec.tight_pattern, spec.pt_bins)
    bkg, bkg_var, bkg_edges, bkg_keys = add_hists(spec.data_file, spec.data_dir, spec.nontight_pattern, spec.pt_bins)
    sig, sig_var, sig_edges, sig_keys = add_hists(spec.signal_file, spec.signal_dir, spec.signal_pattern, spec.pt_bins)
    if np.max(np.abs(edges - bkg_edges)) > 1e-9 or np.max(np.abs(edges - sig_edges)) > 1e-9:
        raise ValueError(f"{spec.tag}: data/background/signal isolation binnings do not match")

    tight, tight_var, rb_edges = ppg12_variable_rebin(
        tight, tight_var, edges, low_region_group=low_region_group, high_region_group=high_region_group
    )
    bkg, bkg_var, _ = ppg12_variable_rebin(
        bkg, bkg_var, edges, low_region_group=low_region_group, high_region_group=high_region_group
    )
    sig, sig_var, _ = ppg12_variable_rebin(
        sig, sig_var, edges, low_region_group=low_region_group, high_region_group=high_region_group
    )

    tight_density, tight_density_var = bin_width_scale(tight, tight_var, rb_edges)
    bkg_density, bkg_density_var = bin_width_scale(bkg, bkg_var, rb_edges)
    sig_density, sig_density_var = bin_width_scale(sig, sig_var, rb_edges)

    tight_tail = tail_integral(tight_density, rb_edges, tail_low)
    bkg_tail = tail_integral(bkg_density, rb_edges, tail_low)
    bkg_scale = tight_tail / bkg_tail if bkg_tail > 0 else math.nan
    bkg_scaled = bkg_density * bkg_scale
    bkg_scaled_var = bkg_density_var * bkg_scale * bkg_scale

    data_signal_density_integral = float(np.sum(tight_density) - np.sum(bkg_scaled))
    signal_density_integral = float(np.sum(sig_density))
    sig_scale = data_signal_density_integral / signal_density_integral if signal_density_integral > 0 else math.nan
    sig_scaled = sig_density * sig_scale
    sig_scaled_var = sig_density_var * sig_scale * sig_scale

    return {
        "edges": rb_edges,
        "tight": tight_density,
        "tight_err": np.sqrt(np.clip(tight_density_var, 0.0, None)),
        "background": bkg_scaled,
        "background_err": np.sqrt(np.clip(bkg_scaled_var, 0.0, None)),
        "signal_mc": sig_scaled,
        "signal_mc_err": np.sqrt(np.clip(sig_scaled_var, 0.0, None)),
        "metadata": {
            "tag": spec.tag,
            "data_file": str(spec.data_file),
            "signal_file": str(spec.signal_file),
            "tight_keys": tight_keys,
            "nontight_keys": bkg_keys,
            "signal_keys": sig_keys,
            "tail_low_GeV": tail_low,
            "tail_tight_sum_after_width_scale": tight_tail,
            "tail_background_sum_after_width_scale": bkg_tail,
            "background_tail_scale": bkg_scale,
            "data_signal_density_integral": data_signal_density_integral,
            "signal_mc_density_integral_before_scale": signal_density_integral,
            "signal_mc_scale": sig_scale,
            "sideband_note": spec.sideband_note,
            "low_region_rebin_group_below_2p5_GeV": low_region_group,
            "high_region_rebin_group_above_2p5_GeV": high_region_group,
            "old_reference": "ppg12codeGit/plotting/CONF_plots.C: rebin, Scale(width), normalize non-tight to tight tail above 6 GeV, normalize signal MC to tight-minus-background integral.",
        },
    }


def draw_step_band(
    ax: plt.Axes,
    edges: np.ndarray,
    low: np.ndarray,
    high: np.ndarray,
    facecolor: str,
    edgecolor: str,
    label: str,
    *,
    alpha: float,
    linewidth: float,
) -> None:
    """Draw a ROOT-like filled step band without per-bin rectangle borders."""

    x = np.ravel(np.column_stack([edges[:-1], edges[1:]]))
    y_low = np.repeat(low, 2)
    y_high = np.repeat(high, 2)
    ax.fill_between(x, y_low, y_high, color=facecolor, alpha=alpha, step=None, label=label, linewidth=0)
    ax.stairs(high, edges, baseline=low, fill=False, color=edgecolor, linewidth=linewidth, alpha=0.85)


def draw_stack(result: dict[str, object], spec: SampleSpec, out_path: Path) -> None:
    edges = result["edges"]
    tight = result["tight"]
    tight_err = result["tight_err"]
    bkg = result["background"]
    sig = result["signal_mc"]
    widths = np.diff(edges)
    centers = edges[:-1] + 0.5 * widths

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 1.6,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "xtick.minor.width": 1.1,
        "ytick.minor.width": 1.1,
    })

    fig, ax = plt.subplots(figsize=(7.2, 6.2), dpi=220)
    draw_step_band(ax, edges, np.zeros_like(bkg), bkg, "#f5a3a0", "#6e2c2c", "Data (Background)", alpha=0.72, linewidth=0.9)
    draw_step_band(ax, edges, bkg, bkg + sig, "#8f91f0", "#2c2a6b", "Signal MC", alpha=0.68, linewidth=0.9)
    ax.errorbar(centers, tight, yerr=tight_err, fmt="o", color="black", ecolor="black", elinewidth=1.3, capsize=0, markersize=5.8, label="Data (Signal)")

    ymax = max(float(np.nanmax(tight + tight_err)), float(np.nanmax(bkg + sig)))
    ax.set_xlim(-1.0, 15.0)
    ax.set_ylim(0, ymax * 1.18 if ymax > 0 else 1.0)
    ax.set_xlabel(r"$E_T^{\mathrm{iso,reco}}\ \mathrm{[GeV]}$", fontsize=20, loc="right")
    ax.set_ylabel("Counts / Bin Width", fontsize=20)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=17, length=7)
    ax.tick_params(axis="both", which="minor", length=4)
    ax.minorticks_on()

    ax.text(0.95, 0.93, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=19)
    ax.text(0.95, 0.84, spec.collision_label, transform=ax.transAxes, ha="right", va="top", fontsize=18)
    ax.text(0.95, 0.76, spec.pt_label, transform=ax.transAxes, ha="right", va="top", fontsize=18)
    ax.text(0.95, 0.68, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, ha="right", va="top", fontsize=18)

    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index("Data (Signal)"), labels.index("Data (Background)"), labels.index("Signal MC")]
    ax.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper right", bbox_to_anchor=(0.95, 0.58), frameon=False, fontsize=17, handlelength=1.4, handletextpad=0.6)
    fig.tight_layout(pad=1.0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def draw_summary(results: list[tuple[SampleSpec, dict[str, object], Path]], out_path: Path) -> None:
    results = [item for item in results if item[0].include_in_summary]
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 10.2), dpi=220)
    for ax, (spec, result, _) in zip(axes.flat, results):
        edges = result["edges"]
        tight = result["tight"]
        tight_err = result["tight_err"]
        bkg = result["background"]
        sig = result["signal_mc"]
        widths = np.diff(edges)
        centers = edges[:-1] + 0.5 * widths
        draw_step_band(ax, edges, np.zeros_like(bkg), bkg, "#f5a3a0", "#6e2c2c", "Data (Background)", alpha=0.72, linewidth=0.7)
        draw_step_band(ax, edges, bkg, bkg + sig, "#8f91f0", "#2c2a6b", "Signal MC", alpha=0.68, linewidth=0.7)
        ax.errorbar(centers, tight, yerr=tight_err, fmt="o", color="black", ecolor="black", elinewidth=1.0, capsize=0, markersize=3.9, label="Data (Signal)")
        ymax = max(float(np.nanmax(tight + tight_err)), float(np.nanmax(bkg + sig)))
        ax.set_xlim(-1.0, 15.0)
        ax.set_ylim(0, ymax * 1.22 if ymax > 0 else 1.0)
        ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=11, length=5)
        ax.minorticks_on()
        ax.text(0.96, 0.92, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=12.5)
        ax.text(0.96, 0.82, spec.collision_label, transform=ax.transAxes, ha="right", va="top", fontsize=12)
        ax.text(0.96, 0.73, spec.pt_label, transform=ax.transAxes, ha="right", va="top", fontsize=12)
        ax.set_xlabel(r"$E_T^{\mathrm{iso,reco}}\ \mathrm{[GeV]}$", fontsize=12)
        ax.set_ylabel("Counts / Bin Width", fontsize=12)
    axes.flat[0].legend(loc="upper right", bbox_to_anchor=(0.98, 0.63), frameon=False, fontsize=11)
    fig.suptitle("PPG12-style isolation-stack diagnostic, current local outputs", y=0.985, fontsize=18, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def build_specs(args: argparse.Namespace) -> list[SampleSpec]:
    pp_bins = ("16_18", "18_20", "20_22", "22_24", "24_26", "26_28", "28_32", "32_36")
    auau_bins = ("16_18", "18_20", "20_22", "22_24", "24_26", "26_35")
    specs: list[SampleSpec] = [
        SampleSpec(
            tag="pp_current_baseline_16to22_style_check",
            collision_label=r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV}$",
            pt_label=r"$16 < E_T^\gamma < 22\ \mathrm{GeV}$",
            pt_bins=("16_18", "18_20", "20_22"),
            data_file=args.pp_data,
            signal_file=args.pp_signal,
            data_dir="PPG12_scaledtrigger30",
            signal_dir="SIM",
            tight_pattern="h_Eiso_tight_pT_{pt}",
            nontight_pattern="h_Eiso_nonTight_pT_{pt}",
            signal_pattern="h_Eiso_tight_pT_{pt}",
            output_name="pp_current_baseline_isolation_stack_16to22_oldstyle_check.png",
            sideband_note="Direct current-output check of the old screenshot pT window.",
            include_in_summary=False,
        ),
        SampleSpec(
            tag="pp_current_baseline",
            collision_label=r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV}$",
            pt_label=r"$16 < E_T^\gamma < 36\ \mathrm{GeV}$",
            pt_bins=pp_bins,
            data_file=args.pp_data,
            signal_file=args.pp_signal,
            data_dir="PPG12_scaledtrigger30",
            signal_dir="SIM",
            tight_pattern="h_Eiso_tight_pT_{pt}",
            nontight_pattern="h_Eiso_nonTight_pT_{pt}",
            signal_pattern="h_Eiso_tight_pT_{pt}",
            output_name="pp_current_baseline_isolation_stack_16to36.png",
            sideband_note="Current pp baseline has 32-36 as the final pT bin, so this is 16-36 rather than exact 16-35.",
        )
    ]
    for cent in ("0_20", "20_50", "50_80"):
        cent_label = cent.replace("_", "-")
        specs.append(
            SampleSpec(
                tag=f"auau_default_bdt_cent_{cent}",
                collision_label=rf"$\mathrm{{Au+Au}}\ \sqrt{{s_{{NN}}}}=200\ \mathrm{{GeV}},\ {cent_label}\%$",
                pt_label=r"$16 < E_T^\gamma < 35\ \mathrm{GeV}$",
                pt_bins=auau_bins,
                data_file=args.auau_data,
                signal_file=args.auau_signal,
                data_dir="MBD_NS_geq_2_vtx_lt_150",
                signal_dir="SIM",
                tight_pattern=f"h_Eiso_tight_isoR40_pT_{{pt}}_cent_{cent}",
                nontight_pattern=f"h_Eiso_nonTight_isoR40_pT_{{pt}}_cent_{cent}",
                signal_pattern=f"h_Eiso_tight_isoR40_pT_{{pt}}_cent_{cent}",
                output_name=f"auau_default_bdt_isolation_stack_16to35_cent{cent}.png",
                sideband_note="Current AuAu non-tight histogram is nonTightAuAuBDTComplement, not a bounded sideband.",
            )
        )
    return specs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pp-data", type=Path, default=DEFAULT_PP_DATA)
    parser.add_argument("--pp-signal", type=Path, default=DEFAULT_PP_SIGNAL)
    parser.add_argument("--auau-data", type=Path, default=DEFAULT_AUAU_DATA)
    parser.add_argument("--auau-signal", type=Path, default=DEFAULT_AUAU_SIGNAL)
    parser.add_argument("--outdir", type=Path, default=REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/isolation_stack_ppg12_style_20260628")
    parser.add_argument("--tail-low", type=float, default=6.0)
    parser.add_argument("--low-region-rebin", type=int, default=2, help="Group this many source bins below 2.5 GeV. Old CONF_plots.C used 1; 2 is cleaner for current low-stat outputs.")
    parser.add_argument("--high-region-rebin", type=int, default=5, help="Group this many source bins above 2.5 GeV, matching old CONF_plots.C default.")
    args = parser.parse_args()

    manifest: dict[str, object] = {
        "description": "Current-output recreation of ppg12codeGit/plotting/CONF_plots.C isolation stack normalization.",
        "tail_normalization": "non-tight data scaled so Eiso >= tail_low matches tight data after variable rebin and bin-width scaling; signal MC scaled to tight minus scaled non-tight density integral.",
        "tail_low_GeV": args.tail_low,
        "low_region_rebin_group_below_2p5_GeV": args.low_region_rebin,
        "high_region_rebin_group_above_2p5_GeV": args.high_region_rebin,
        "outputs": [],
    }
    rendered: list[tuple[SampleSpec, dict[str, object], Path]] = []
    for spec in build_specs(args):
        result = old_ppg12_stack(spec, args.tail_low, args.low_region_rebin, args.high_region_rebin)
        out_path = args.outdir / spec.output_name
        draw_stack(result, spec, out_path)
        rendered.append((spec, result, out_path))
        metadata = dict(result["metadata"])
        metadata["output_png"] = str(out_path)
        manifest["outputs"].append(metadata)

    summary_path = args.outdir / "slide_isolation_stack_pp_auau_16to35_current_outputs.png"
    draw_summary(rendered, summary_path)
    manifest["summary_png"] = str(summary_path)
    manifest_path = args.outdir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary_path)
    print(manifest_path)
    for _, _, out_path in rendered:
        print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
