#!/usr/bin/env python3
"""Build the THE-42 WP80 AuAu shower-shape overlay slide.

This consumes merged RecoilJets ROOT outputs, sums the requested pT and
centrality slices, and renders the PPG12-style no-preselection/preselection/
tight-selection e11/e33 overlay as a full 16:9 slide PNG.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
SCRIPTS_DIR = next((p for p in THIS_FILE.parents if p.name == "scripts"), THIS_FILE.parent)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


DEFAULT_OUTDIR = Path(
    "dataOutput/auauTightBDTValidation/"
    "THE42_wp80_centlinear_ss_20260604"
)
DEFAULT_PT_BINS = "22:24,24:26,26:28"
DEFAULT_CENT_BINS = "0:10,10:20,20:30,30:40,40:50,50:60,60:80"
STAGE_SPECS = (
    {
        "key": "inclusive",
        "title": "No preselection",
        "data_tag": "inclusive",
        "signal_tag": "inclusive_sig",
        "inclusive_tag": "inclusive_bkg",
    },
    {
        "key": "pre",
        "title": "NPB preselection",
        "data_tag": "pre",
        "signal_tag": "pre_sig",
        "inclusive_tag": "pre_bkg",
    },
    {
        "key": "tight",
        "title": "Tight WP80 BDT selection",
        "data_tag": "tight",
        "signal_tag": "tight_sig",
        "inclusive_tag": "tight_bkg",
    },
)
SHOWER_SHAPE_VARS = (
    "weta",
    "wphi",
    "weta33",
    "wphi33",
    "weta35",
    "wphi53",
    "et1",
    "e11e33",
    "e32e35",
    "npbScore",
)


@dataclass
class RootMatch:
    root: Path
    object_path: str


@dataclass
class Curve:
    name: str
    tag: str
    edges: np.ndarray
    centers: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    raw_counts: np.ndarray
    raw_errors: np.ndarray
    raw_integral: float
    displayed_integral: float
    matches: list[RootMatch]


def parse_bins(spec: str) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" not in chunk:
            raise ValueError(f"Bin '{chunk}' must be formatted as low:high")
        lo_s, hi_s = chunk.split(":", 1)
        lo, hi = int(float(lo_s)), int(float(hi_s))
        if hi <= lo:
            raise ValueError(f"Invalid bin '{chunk}': high must be > low")
        out.append((lo, hi))
    if not out:
        raise ValueError(f"No bins parsed from '{spec}'")
    return out


def hist_name(var: str, tag: str, pt_bin: tuple[int, int], cent_bin: tuple[int, int] | None) -> str:
    lo, hi = pt_bin
    name = f"h_ss_{var}_{tag}_pT_{lo}_{hi}"
    if cent_bin is not None:
        clo, chi = cent_bin
        name += f"_cent_{clo}_{chi}"
    return name


def compile_trigger_regex(pattern: str | None) -> re.Pattern[str] | None:
    if not pattern:
        return None
    return re.compile(pattern)


def open_root(path: Path):
    import ROOT  # imported lazily so --help works without ROOT startup

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return f


def walk_root_dir(directory, prefix: str = ""):
    for key in directory.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        path = f"{prefix}/{name}" if prefix else name
        yield path, obj
        if obj.InheritsFrom("TDirectory"):
            yield from walk_root_dir(obj, path)


def find_hists(root_path: Path, target_name: str, trigger_regex: re.Pattern[str] | None):
    import ROOT  # noqa: F401

    matches = []
    f = open_root(root_path)
    try:
        for object_path, obj in walk_root_dir(f):
            if not obj.InheritsFrom("TH1"):
                continue
            if Path(object_path).name != target_name:
                continue
            if trigger_regex and not trigger_regex.search(object_path):
                continue
            h = obj.Clone(f"{target_name}_{len(matches)}")
            h.SetDirectory(0)
            matches.append((h, RootMatch(root=root_path, object_path=object_path)))
    finally:
        f.Close()
    return matches


def add_hist(acc, hist):
    if acc is None:
        out = hist.Clone(f"{hist.GetName()}_sum")
        out.SetDirectory(0)
        return out
    acc.Add(hist)
    return acc


def sum_stage_hist(
    roots: list[Path],
    var: str,
    tag: str,
    pt_bins: list[tuple[int, int]],
    cent_bins: list[tuple[int, int]],
    trigger_regex: re.Pattern[str] | None,
):
    acc = None
    matches: list[RootMatch] = []
    attempted: list[str] = []

    for pt_bin in pt_bins:
        found_cent_specific = False
        for cent_bin in cent_bins:
            target = hist_name(var, tag, pt_bin, cent_bin)
            attempted.append(target)
            for root in roots:
                for hist, match in find_hists(root, target, trigger_regex):
                    acc = add_hist(acc, hist)
                    matches.append(match)
                    found_cent_specific = True
        if found_cent_specific:
            continue

        target = hist_name(var, tag, pt_bin, None)
        attempted.append(target)
        for root in roots:
            for hist, match in find_hists(root, target, trigger_regex):
                acc = add_hist(acc, hist)
                matches.append(match)

    if acc is None:
        raise RuntimeError(
            f"No histograms found for var={var} tag={tag}. "
            f"Attempted names include: {attempted[:8]}"
        )
    return acc, matches, attempted


def rebin_if_needed(hist, factor: int):
    if factor <= 1:
        return hist
    nb = hist.GetNbinsX()
    if nb % factor != 0:
        raise ValueError(f"Cannot rebin {hist.GetName()} with nbins={nb} by factor={factor}")
    out = hist.Rebin(factor, f"{hist.GetName()}_rebin{factor}")
    out.SetDirectory(0)
    return out


def curve_from_hist(
    hist,
    matches: list[RootMatch],
    name: str,
    tag: str,
    x_min: float,
    x_max: float,
    rebin: int,
) -> Curve:
    hist = rebin_if_needed(hist, rebin)
    nb = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    errs = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)

    mask = (centers >= x_min) & (centers < x_max)
    edges_view = np.concatenate(([edges[np.flatnonzero(mask)[0]]], edges[np.flatnonzero(mask) + 1])) if np.any(mask) else edges
    centers_view = centers[mask]
    counts_view = counts[mask]
    errs_view = errs[mask]
    raw_integral = float(np.sum(counts))
    displayed_integral = float(np.sum(counts_view))
    values = np.zeros_like(counts_view)
    errors = np.zeros_like(errs_view)
    if displayed_integral > 0:
        values = counts_view / displayed_integral
        errors = errs_view / displayed_integral
    return Curve(
        name=name,
        tag=tag,
        edges=edges_view,
        centers=centers_view,
        values=values,
        errors=errors,
        raw_counts=counts_view,
        raw_errors=errs_view,
        raw_integral=raw_integral,
        displayed_integral=displayed_integral,
        matches=matches,
    )


def load_curve(
    roots: list[Path],
    var: str,
    tag: str,
    label: str,
    pt_bins: list[tuple[int, int]],
    cent_bins: list[tuple[int, int]],
    trigger_regex: re.Pattern[str] | None,
    x_min: float,
    x_max: float,
    rebin: int,
) -> Curve:
    hist, matches, _attempted = sum_stage_hist(roots, var, tag, pt_bins, cent_bins, trigger_regex)
    return curve_from_hist(hist, matches, label, tag, x_min, x_max, rebin)


def validate_same_binning(curves: Iterable[Curve]) -> None:
    curves = list(curves)
    ref = curves[0]
    for curve in curves[1:]:
        if len(curve.edges) != len(ref.edges) or not np.allclose(curve.edges, ref.edges):
            raise RuntimeError(f"Binning mismatch between {ref.name} and {curve.name}")


def draw_sphenix_label(ax, pt_label: str, cent_label: str, collision_label: str = "Au+Au") -> None:
    ax.text(0.045, 0.930, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes,
            ha="left", va="top", fontsize=11.5)
    ax.text(0.045, 0.860, f"{collision_label}, {pt_label}", transform=ax.transAxes,
            ha="left", va="top", fontsize=10.5)
    ax.text(0.045, 0.795, cent_label, transform=ax.transAxes,
            ha="left", va="top", fontsize=10.5)


def plot_stage(
    ax,
    rax,
    title: str,
    data: Curve | None,
    signal: Curve,
    inclusive: Curve,
    x_label: str,
    sim_only: bool,
) -> dict[str, float]:
    curves = [signal, inclusive] if sim_only else [data, signal, inclusive]
    validate_same_binning(curves)
    ax.stairs(signal.values, signal.edges, color="#D73B32", linewidth=2.0, label="Signal MC")
    ax.stairs(inclusive.values, inclusive.edges, color="#2364C8", linewidth=2.0, label="Inclusive MC")
    ymax_terms = [
        float(np.max(signal.values)) if len(signal.values) else 0.0,
        float(np.max(inclusive.values)) if len(inclusive.values) else 0.0,
    ]
    if data is not None:
        ax.errorbar(data.centers, data.values, yerr=data.errors, fmt="o", color="#111111",
                    markersize=3.6, elinewidth=1.0, capsize=1.8, label="Data", zorder=5)
        ymax_terms.append(float(np.max(data.values + data.errors)) if len(data.values) else 0.0)
    ymax = max(ymax_terms)
    ax.set_ylim(0.0, max(0.035, ymax * 1.24))
    ax.set_title(title, fontsize=16.5, fontweight="bold", pad=8)
    ax.grid(True, color="#DDE3EA", linewidth=0.8, alpha=0.85)
    ax.tick_params(axis="both", labelsize=10.5, direction="in", top=True, right=True)
    ax.set_xlim(signal.edges[0], signal.edges[-1])

    if sim_only:
        residual = signal.values - inclusive.values
        residual_err = np.sqrt(signal.errors ** 2 + inclusive.errors ** 2)
        r_color = "#374151"
    else:
        assert data is not None
        residual = data.values - inclusive.values
        residual_err = np.sqrt(data.errors ** 2 + inclusive.errors ** 2)
        r_color = "#111111"
    rax.axhline(0.0, color="#6B7280", linewidth=1.0, linestyle=":")
    rax.errorbar(signal.centers, residual, yerr=residual_err, fmt="o", color=r_color,
                 markersize=3.0, elinewidth=0.9, capsize=1.4)
    rmax = float(np.max(np.abs(residual) + residual_err)) if len(residual) else 0.01
    rax.set_ylim(-max(0.010, rmax * 1.25), max(0.010, rmax * 1.25))
    rax.grid(True, color="#E5EAF0", linewidth=0.7, alpha=0.9)
    rax.tick_params(axis="both", labelsize=9.5, direction="in", top=True, right=True)
    rax.set_xlim(signal.edges[0], signal.edges[-1])
    rax.set_xlabel(x_label, fontsize=12.5)
    metrics = {
        "signal_displayed_entries": signal.displayed_integral,
        "inclusive_displayed_entries": inclusive.displayed_integral,
        "max_abs_residual": float(np.max(np.abs(residual))) if len(residual) else math.nan,
    }
    if data is not None:
        metrics["data_displayed_entries"] = data.displayed_integral
    return metrics


def add_text_box(fig, xywh: tuple[float, float, float, float], title: str, body: str, face: str, edge: str) -> None:
    ax = fig.add_axes(xywh)
    ax.set_axis_off()
    rect = plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor=face,
                         edgecolor=edge, linewidth=1.2)
    ax.add_patch(rect)
    multiline = "\n" in body
    title_y = 0.76 if multiline else 0.70
    body_y = 0.25 if multiline else 0.34
    body_size = 11.5 if multiline else 12.2
    ax.text(0.030, title_y, title, ha="left", va="center", fontsize=13.5,
            fontweight="bold", color="#111827")
    ax.text(0.030, body_y, body, ha="left", va="center", fontsize=body_size,
            color="#263241", linespacing=1.06)


def build_slide(args: argparse.Namespace) -> dict[str, object]:
    import ROOT  # noqa: F401

    data_roots = [] if args.sim_only else [Path(p).resolve() for p in args.data_root]
    signal_roots = [Path(p).resolve() for p in args.signal_root]
    inclusive_roots = [Path(p).resolve() for p in args.inclusive_root]
    for path in data_roots + signal_roots + inclusive_roots:
        if not path.exists():
            raise FileNotFoundError(path)

    pt_bins = parse_bins(args.pt_bins)
    cent_bins = parse_bins(args.cent_bins)
    trigger_regex = compile_trigger_regex(args.trigger_regex)
    pt_label = rf"{pt_bins[0][0]} <= $p_T^\gamma$ < {pt_bins[-1][1]} GeV"
    cent_label = f"{cent_bins[0][0]}-{cent_bins[-1][1]}% centrality"
    x_label = r"$e_{11}/e_{33}$"

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.edgecolor": "#111827",
        "axes.linewidth": 1.1,
        "xtick.color": "#111827",
        "ytick.color": "#111827",
    })

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    title = (
        "WP80 BDT selection changes the simulated shower-shape comparison"
        if args.sim_only
        else "WP80 BDT selection changes the Au+Au shower-shape comparison"
    )
    fig.text(0.046, 0.958, title,
             ha="left", va="top", fontsize=24.8, fontweight="bold", color="#111827")
    fig.text(0.046, 0.908,
             rf"{pt_label}, {cent_label}; tight BDT uses {args.t80_label}",
             ha="left", va="top", fontsize=14.8, color="#263241")

    add_text_box(
        fig,
        (0.046, 0.772, 0.525, 0.065),
        "Comparison flow",
        "Shown before cuts, after NPB preselection, and after the WP80 tight-BDT split.",
        "#F4F8FC",
        "#C8D8EA",
    )
    add_text_box(
        fig,
        (0.618, 0.772, 0.336, 0.065),
        "Sample contract",
        ("Signal MC: Photon+Jet embedded 12+20\nInclusive MC: Jet embedded 12+20+30+40"
         if args.sim_only else
         "Data: Au+Au | Signal: photon-jet emb. | Inclusive: jet emb."),
        "#F7FAF8",
        "#B9DCC8",
    )

    gs = fig.add_gridspec(
        nrows=2,
        ncols=3,
        left=0.052,
        right=0.958,
        bottom=0.165,
        top=0.722,
        height_ratios=[3.1, 1.0],
        wspace=0.155,
        hspace=0.025,
    )

    stage_manifest: dict[str, object] = {}
    stage_metrics: dict[str, dict[str, float]] = {}
    legend_handles = None
    for idx, spec in enumerate(STAGE_SPECS):
        ax = fig.add_subplot(gs[0, idx])
        rax = fig.add_subplot(gs[1, idx], sharex=ax)
        data = None
        if not args.sim_only:
            data = load_curve(data_roots, args.var, spec["data_tag"], "Data", pt_bins, cent_bins,
                              trigger_regex, args.x_min, args.x_max, args.rebin)
        signal = load_curve(signal_roots, args.var, spec["signal_tag"], "Signal MC", pt_bins, cent_bins,
                            trigger_regex, args.x_min, args.x_max, args.rebin)
        inclusive = load_curve(inclusive_roots, args.var, spec["inclusive_tag"], "Inclusive MC", pt_bins,
                               cent_bins, trigger_regex, args.x_min, args.x_max, args.rebin)
        metrics = plot_stage(ax, rax, spec["title"], data, signal, inclusive, x_label, args.sim_only)
        stage_metrics[spec["key"]] = metrics
        if idx == 0:
            ax.set_ylabel("normalized counts", fontsize=12.5)
            rax.set_ylabel("Sig. - Incl. MC" if args.sim_only else "Data - Incl. MC", fontsize=10.8)
        else:
            ax.set_yticklabels([])
            rax.set_yticklabels([])
        draw_sphenix_label(ax, pt_label, cent_label, "Au+Au embedded" if args.sim_only else "Au+Au")
        if idx == 1:
            handles, labels = ax.get_legend_handles_labels()
            legend_handles = (handles, labels)
        stage_manifest[spec["key"]] = {
            "title": spec["title"],
            "tags": {
                "data": spec["data_tag"],
                "signal": spec["signal_tag"],
                "inclusive": spec["inclusive_tag"],
            },
            "metrics": metrics,
            "matches": {
                "data": [] if data is None else [m.__dict__ | {"root": str(m.root)} for m in data.matches],
                "signal": [m.__dict__ | {"root": str(m.root)} for m in signal.matches],
                "inclusive": [m.__dict__ | {"root": str(m.root)} for m in inclusive.matches],
            },
        }

    if legend_handles:
        handles, labels = legend_handles
        order = (
            [labels.index("Signal MC"), labels.index("Inclusive MC")]
            if args.sim_only else
            [labels.index("Data"), labels.index("Signal MC"), labels.index("Inclusive MC")]
        )
        handles = [handles[i] for i in order]
        labels = [labels[i] for i in order]
        leg = fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.955, 0.895),
                         ncol=3, frameon=True, fontsize=13.2, handlelength=2.2,
                         columnspacing=1.3)
        leg.get_frame().set_facecolor("white")
        leg.get_frame().set_edgecolor("#CAD4E0")
        leg.get_frame().set_linewidth(1.0)

    error_note = (
        "Simulation-only interim view: data points are intentionally omitted; lower panels show Signal MC - Inclusive MC with ROOT bin errors propagated in quadrature."
        if args.sim_only else
        "Error bars use ROOT histogram bin errors scaled by the same unit-area normalization; residual-panel errors combine Data and Inclusive MC in quadrature."
    )
    fig.text(0.052, 0.095,
             error_note,
             ha="left", va="center", fontsize=11.8, color="#4B5563")
    fig.text(0.052, 0.061,
             "Follow-up-ready histograms in the same outputs include other shower-shape variables, tight/non-tight sidebands, isolation, ABCD, and preselection-failure audits.",
             ha="left", va="center", fontsize=11.8, color="#4B5563")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    png_path = args.output_dir / args.output_name
    fig.savefig(png_path, dpi=SLIDE_DPI)
    plt.close(fig)

    manifest = {
        "campaign": "THE-42 WP80 AuAu SS overlay",
        "output_png": str(png_path.resolve()),
        "script": str(THIS_FILE),
        "plot_mode": "simulation_only" if args.sim_only else "data_signal_inclusive",
        "variable": args.var,
        "x_label": "e11/e33",
        "x_range_displayed": [args.x_min, args.x_max],
        "pt_bins": pt_bins,
        "cent_bins": cent_bins,
        "selection_definition": {
            "tight": "score > T80(c)",
            "non_tight": "score <= T80(c)",
            "T80": args.t80_label,
        },
        "roots": {
            "data": [str(p) for p in data_roots],
            "signal": [str(p) for p in signal_roots],
            "inclusive": [str(p) for p in inclusive_roots],
        },
        "normalization": "unit area over displayed x-range after summing requested pT and centrality bins",
        "stat_errors": (
            "ROOT bin errors scaled by the same normalization; sim-only residual errors are "
            "sqrt(signal_mc_err^2 + inclusive_mc_err^2)"
            if args.sim_only else
            "ROOT bin errors scaled by the same normalization; residual errors are "
            "sqrt(data_err^2 + inclusive_mc_err^2)"
        ),
        "stages": stage_manifest,
        "other_histogram_families_available": {
            "shower_shape_variables": list(SHOWER_SHAPE_VARS),
            "selection_stages": ["inclusive", "pre", "tight", "nonTight"],
            "diagnostic_families": [
                "h_Eiso_*",
                "h_isIsolated_isTight",
                "h_notIsolated_isTight",
                "h_isIsolated_notTight",
                "h_notIsolated_notTight",
                "h_preFail_*",
                "h_preRefFail_*",
                "h_tightFail_*",
                "h_tightFailMask_*",
            ],
        },
    }
    manifest_path = png_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    script_path = png_path.with_suffix(".speaker_script.md")
    script_path.write_text(build_speaker_script(png_path, manifest_path, stage_metrics, args.sim_only))
    return {"png": png_path, "manifest": manifest_path, "speaker_script": script_path}


def build_speaker_script(
    png_path: Path,
    manifest_path: Path,
    stage_metrics: dict[str, dict[str, float]],
    sim_only: bool,
) -> str:
    tight = stage_metrics.get("tight", {})
    tight_resid = tight.get("max_abs_residual", math.nan)
    resid_text = "not yet available" if not math.isfinite(tight_resid) else f"{tight_resid:.4f}"
    if sim_only:
        return (
            "# THE-42 WP80 AuAu Simulation-Only Shower-Shape Overlay Script\n\n"
            "This is the simulation-only interim version of the THE-42 shower-shape slide. I am showing the "
            "same e11 over e33 comparison through no preselection, the NPB-style preselection, and the tight "
            "WP80 BDT selection, but with the Au+Au data points intentionally omitted while the data lane is still pending.\n\n"
            "The two curves are the current THE-42 embedded samples: Photon+Jet 12+20 for the signal template, "
            "and inclusive Jet 12+20+30+40 for the background-dominated template. Each panel is unit-normalized "
            "over the displayed range after summing 22 to 28 GeV and 0 to 80 percent centrality.\n\n"
            "The important definition is unchanged from the full planned slide. Tight WP80 means the BDT score is "
            "above the centrality-dependent threshold T80 of centrality equals 0.5593 plus 0.00166 times centrality percentile.\n\n"
            f"In the tight panel, the maximum absolute displayed Signal-minus-Inclusive difference is {resid_text} "
            "in normalized-count units. The exact ROOT inputs, histogram names, normalization choice, and matched "
            f"object paths are recorded in `{manifest_path.name}` next to the PNG.\n\n"
            "Once plot-usable Au+Au data histograms are available, this same slide builder can be rerun without "
            "the simulation-only option to restore the black data points and the data-minus-inclusive residual panels.\n\n"
            f"PNG: `{png_path}`\n"
            f"Manifest: `{manifest_path}`\n"
        )
    return (
        "# THE-42 WP80 AuAu Shower-Shape Overlay Script\n\n"
        "This slide is the first direct visual check of the centrality-dependent WP80 BDT working point in the "
        "RecoilJets histogram output. I am showing the same shower-shape variable, e11 over e33, through three "
        "successive selections: no preselection, the NPB-style preselection, and then the tight WP80 BDT selection.\n\n"
        "The important definition is in the subtitle. The tight category is not a single global BDT score cut; it "
        "uses the centrality-dependent line derived from the previous working-point slide, T80 of centrality equals "
        "0.5593 plus 0.00166 times centrality percentile.\n\n"
        "The red curve is the Photon+Jet embedded signal template, the blue curve is the inclusive-jet embedded "
        "template, and the black points are the Au+Au data sample. The lower panels show data minus inclusive MC, "
        "so the audience can see whether each selection stage is moving the data away from the generic inclusive-jet "
        "shape and toward the signal-like template.\n\n"
        f"For the tight panel, the maximum absolute displayed residual against inclusive MC is {resid_text} in "
        "normalized-count units. The exact ROOT inputs, histogram names, normalization choice, and matched object "
        f"paths are recorded in `{manifest_path.name}` next to the PNG.\n\n"
        "The next step after this slide is not to stop at e11 over e33. The same output files retain the other "
        "shower-shape, isolation, ABCD, tight/non-tight, and preselection-failure histograms, so any suspicious "
        "feature of this comparison can be chased without rerunning the whole campaign.\n\n"
        f"PNG: `{png_path}`\n"
        f"Manifest: `{manifest_path}`\n"
    )


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", nargs="+", help="Merged Au+Au data ROOT file(s).")
    ap.add_argument("--signal-root", nargs="+", required=True, help="Merged Photon+Jet embedded ROOT file(s).")
    ap.add_argument("--inclusive-root", nargs="+", required=True, help="Merged inclusive-jet embedded ROOT file(s).")
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--output-name", default="the42_wp80_auau_ss_overlay_slide.png")
    ap.add_argument("--var", default="e11e33", choices=SHOWER_SHAPE_VARS)
    ap.add_argument("--pt-bins", default=DEFAULT_PT_BINS)
    ap.add_argument("--cent-bins", default=DEFAULT_CENT_BINS)
    ap.add_argument("--x-min", type=float, default=0.0)
    ap.add_argument("--x-max", type=float, default=1.0)
    ap.add_argument("--rebin", type=int, default=4, help="Rebin factor applied before plotting.")
    ap.add_argument("--trigger-regex", default=None, help="Optional regex to restrict matched ROOT object paths.")
    ap.add_argument("--sim-only", action="store_true", help="Omit Au+Au data and plot Signal MC vs Inclusive MC only.")
    ap.add_argument("--t80-label", default=r"$T_{80}(c)=0.5593+0.00166c$",
                    help="Displayed and manifest-recorded tight-BDT working-point label.")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if not args.sim_only and not args.data_root:
        raise SystemExit("THE-42 slide build failed: --data-root is required unless --sim-only is set")
    try:
        result = build_slide(args)
    except Exception as exc:
        raise SystemExit(f"THE-42 slide build failed: {exc}") from exc
    print(f"wrote_png={result['png']}")
    print(f"wrote_manifest={result['manifest']}")
    print(f"wrote_speaker_script={result['speaker_script']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
