#!/usr/bin/env python3
"""Build a slide-ready pp/AuAu shower-shape preselection comparison.

Columns are the pp reference plus central and peripheral AuAu. Rows are the
three persisted energy-sharing observables. Each cell is read vertically:
after the complete common preselection but before tight ID above, and after
tight ID below. All curves are integrated over 15 < E_T^gamma < 35 GeV.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


THIS_FILE = Path(__file__).resolve()
REPO = THIS_FILE.parents[4]
SOURCE_RENDERER = (
    REPO / "scripts/plotting/auau_bdt/render_the88_final_shower_shape_matrices.py"
)
DEFAULT_OUTDIR = (
    REPO / "dataOutput/slides/the45_jstg_20260720/the88_pp_auau_shower_shapes"
)

VARIABLES = (
    ("weta", r"$w_{\eta}^{\mathrm{COGX}}$"),
    ("e11e33", r"$E_{1\times1}/E_{3\times3}$"),
    ("e32e35", r"$E_{3\times2}/E_{3\times5}$"),
)
SYSTEMS = (
    ("pp", "p+p reference", None, "#EAF2FC", "#91B8E6"),
    ("0_20", "Au+Au 0--20%", "0_20", "#FFF3F1", "#E5A49A"),
    ("50_80", "Au+Au 50--80%", "50_80", "#F1F8F4", "#98C8B0"),
)
STAGES = (
    ("before_tight", "Before tight ID", "1", "pre"),
    ("after_tight", "After tight ID", "2", "tight"),
)
COLORS = {"Data": "#111827", "Photon simulation": "#CF352D", "Inclusive-jet simulation": "#2565D7"}
INK = "#142239"
MUTED = "#53657C"
GRID = "#DCE4EE"


def load_source_module():
    spec = importlib.util.spec_from_file_location("the88_shower_renderer", SOURCE_RENDERER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {SOURCE_RENDERER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rebin_curve(curve, factor: int):
    """Coarsen adjacent displayed bins without changing the full normalization.

    ``Curve.y`` and ``Curve.e`` are already normalized to the original
    full-candidate denominator.  Summing the contents and adding their errors
    in quadrature therefore preserves both the candidate-fraction meaning and
    the corrected unique-candidate uncertainty.
    """
    if factor <= 1:
        return curve
    count = len(curve.x)
    groups = [slice(start, min(start + factor, count)) for start in range(0, count, factor)]
    return curve.__class__(
        x=np.asarray([float(np.mean(curve.x[group])) for group in groups]),
        y=np.asarray([float(np.sum(curve.y[group])) for group in groups]),
        e=np.asarray([float(np.sqrt(np.sum(np.square(curve.e[group])))) for group in groups]),
        raw_entries=curve.raw_entries,
        effective_entries=curve.effective_entries,
        fill_multiplicity=curve.fill_multiplicity,
        total_weight=curve.total_weight,
        retained_weight=curve.retained_weight,
        retained_fraction=curve.retained_fraction,
        underflow_weight=curve.underflow_weight,
        first_bin_weight=curve.first_bin_weight,
        overflow_weight=curve.overflow_weight,
        first_bin_low=curve.first_bin_low,
        first_bin_high=curve.first_bin_high,
        object_names=curve.object_names,
    )


def add_curve(ax, curve, label: str, *, data: bool = False, data_stride: int = 1) -> None:
    color = COLORS[label]
    if data:
        stride = max(1, data_stride)
        ax.errorbar(
            curve.x[::stride],
            curve.y[::stride],
            yerr=curve.e[::stride],
            fmt="o",
            color=color,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.45,
            markersize=3.2,
            elinewidth=0.65,
            capsize=0,
            zorder=5,
        )
    else:
        ax.step(curve.x, curve.y, where="mid", color=color, linewidth=1.8, zorder=3)


def rounded_box(fig, bounds, face, edge="#CAD6E4"):
    ax = fig.add_axes(bounds)
    ax.axis("off")
    ax.add_patch(
        FancyBboxPatch(
            (0, 0),
            1,
            1,
            boxstyle="round,pad=0.005,rounding_size=0.012",
            transform=ax.transAxes,
            facecolor=face,
            edgecolor=edge,
            linewidth=1.0,
        )
    )
    return ax


def resolve_inputs(source):
    data_root, data_meta = source.root_from_pointer(source.AUAU_DATA_POINTER)
    signal_root, signal_meta = source.root_from_pointer(source.AUAU_SIGNAL_POINTER)
    inclusive_root, inclusive_meta = source.root_from_pointer(source.AUAU_INCLUSIVE_POINTER)
    pp_photon_root, pp_photon_meta = source.root_from_pointer(source.PP_PHOTONJET_POINTER)
    pp_data_pointer = source.CURRENT_ROOT / "pp_data_merged/current.json"
    pp_data_root, pp_data_meta = source.root_from_pointer(pp_data_pointer)
    pp_inclusive_pointer, pp_inclusive_meta = source.root_from_pointer(source.PP_INCLUSIVE_POINTER)
    required = source.pp_names("et1", "1")[0]
    # The July-19 current inclusive ROOT is an audited PPG12-parity artifact,
    # but this local ROOT build terminates while probing it for the archived
    # shower-shape family.  The matched July-16 compatibility ROOT is the
    # registered source for this historical shape display.  Select it first so
    # a local display-only probe cannot kill the renderer; record the explicit
    # provenance exception in the manifest instead of silently substituting it.
    if source.PP_INCLUSIVE_COMPAT_ROOT.is_file() and source.root_has_object(
        source.PP_INCLUSIVE_COMPAT_ROOT, required
    ):
        pp_inclusive_root = source.PP_INCLUSIVE_COMPAT_ROOT
        pp_inclusive_source = "registered compatible shower-shape source (current ROOT local-probe exception)"
    else:
        raise FileNotFoundError(
            "The registered compatibility pp inclusive ROOT is unavailable for "
            f"{required}; current pointer was not probed because it terminates this local ROOT build."
        )
    return {
        "data": (data_root, data_meta),
        "signal": (signal_root, signal_meta),
        "inclusive": (inclusive_root, inclusive_meta),
        "pp_photon": (pp_photon_root, pp_photon_meta),
        "pp_data": (pp_data_root, pp_data_meta),
        "pp_inclusive": (pp_inclusive_root, pp_inclusive_meta),
        "pp_inclusive_source": pp_inclusive_source,
    }


def load_pp_data_curve(source, root_path: Path, variable: str, stage: str, label: str):
    """Load the registered full-stat pp data histogram for one ID stage."""
    name = (
        f"PPG12_scaledtrigger30/h2d_{source.PP_H2D_VARIABLES[variable]}"
        f"_eta0_pt1535_cut{stage}"
    )
    handle = source.open_root(root_path)
    histogram = handle.Get(name)
    if not histogram or not histogram.InheritsFrom("TH2"):
        raise KeyError(f"Missing expected p+p data histogram {name}")
    projection = histogram.ProjectionX(f"{label}_{variable}_{stage}_projection")
    projection.SetDirectory(0)
    return source.curve_from_histogram(projection, [name])


def build_curves(source, inputs, *, rebin_factor: int = 1):
    curves = {}
    records = []
    for variable, _ in VARIABLES:
        for stage_key, _, pp_stage, auau_stage in STAGES:
            pp_data = load_pp_data_curve(
                source,
                inputs["pp_data"][0],
                variable,
                pp_stage,
                f"slide_pp_data_{variable}_{stage_key}",
            )
            pp_signal = source.load_pp_curve(
                inputs["pp_photon"][0], variable, pp_stage, f"slide_pp_photon_{variable}_{stage_key}"
            )
            pp_inclusive = source.load_pp_curve(
                inputs["pp_inclusive"][0], variable, pp_stage, f"slide_pp_inclusive_{variable}_{stage_key}"
            )
            curves[("pp", variable, stage_key)] = {
                "Data": rebin_curve(pp_data, rebin_factor),
                "Photon simulation": rebin_curve(pp_signal, rebin_factor),
                "Inclusive-jet simulation": rebin_curve(pp_inclusive, rebin_factor),
            }
            records.extend(
                [
                    source.curve_record(
                        pp_data, system="pp", variable=variable, stage=f"cut{pp_stage}", lane="data"
                    ),
                    source.curve_record(
                        pp_signal, system="pp", variable=variable, stage=f"cut{pp_stage}", lane="photon"
                    ),
                    source.curve_record(
                        pp_inclusive, system="pp", variable=variable, stage=f"cut{pp_stage}", lane="inclusive"
                    ),
                ]
            )
            for cent, _, sim_cent, _, _ in SYSTEMS[1:]:
                fill_multiplicity = source.AUAU_STAGE_FILL_MULTIPLICITY[auau_stage]
                data = source.load_curve(
                    inputs["data"][0],
                    source.data_names(variable, auau_stage, cent),
                    f"slide_data_{variable}_{cent}_{stage_key}",
                    fill_multiplicity=fill_multiplicity,
                )
                signal = source.load_curve(
                    inputs["signal"][0],
                    source.sim_names(variable, auau_stage, "sig", sim_cent),
                    f"slide_signal_{variable}_{cent}_{stage_key}",
                    fill_multiplicity=fill_multiplicity,
                )
                inclusive = source.load_curve(
                    inputs["inclusive"][0],
                    source.sim_names(variable, auau_stage, "bkg", sim_cent),
                    f"slide_inclusive_{variable}_{cent}_{stage_key}",
                    fill_multiplicity=fill_multiplicity,
                )
                curves[(cent, variable, stage_key)] = {
                    "Data": rebin_curve(data, rebin_factor),
                    "Photon simulation": rebin_curve(signal, rebin_factor),
                    "Inclusive-jet simulation": rebin_curve(inclusive, rebin_factor),
                }
                records.extend(
                    [
                        source.curve_record(
                            data,
                            system="AuAu",
                            centrality=cent,
                            variable=variable,
                            stage=auau_stage,
                            lane="data",
                        ),
                        source.curve_record(
                            signal,
                            system="AuAu",
                            centrality=cent,
                            variable=variable,
                            stage=auau_stage,
                            lane="signal",
                        ),
                        source.curve_record(
                            inclusive,
                            system="AuAu",
                            centrality=cent,
                            variable=variable,
                            stage=auau_stage,
                            lane="inclusive",
                        ),
                    ]
                )
    return curves, records


def render(curves, outdir: Path, *, output_stem: str):
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.text(
        0.050,
        0.952,
        "Tight photon ID reshapes shower inputs across p+p and Au+Au",
        ha="left",
        va="top",
        fontsize=28.0,
        fontweight="bold",
        color=INK,
    )
    legend = rounded_box(fig, [0.055, 0.835, 0.890, 0.070], "#F8FAFC")
    legend.text(
        0.025,
        0.67,
        "Read each cell vertically",
        ha="left",
        va="center",
        fontsize=13.5,
        fontweight="bold",
        color=INK,
        transform=legend.transAxes,
    )
    legend.text(
        0.025,
        0.28,
        r"both pass common preselection   |   top: before tight ID   $\downarrow$   bottom: after tight ID   |   $15<E_T^{\gamma}<35$ GeV",
        ha="left",
        va="center",
        fontsize=11.8,
        color=MUTED,
        transform=legend.transAxes,
    )
    legend.legend(
        handles=[
            Line2D([0], [0], color=COLORS["Data"], marker="o", lw=0, markersize=5.5, label="Data"),
            Line2D([0], [0], color=COLORS["Photon simulation"], lw=2.8, label="Photon simulation"),
            Line2D([0], [0], color=COLORS["Inclusive-jet simulation"], lw=2.8, label="Inclusive-jet simulation"),
        ],
        loc="center right",
        bbox_to_anchor=(0.985, 0.49),
        ncol=3,
        frameon=False,
        fontsize=11.8,
        handlelength=2.3,
        columnspacing=1.7,
    )

    left = 0.180
    right = 0.955
    top = 0.755
    bottom = 0.105
    hgap = 0.022
    vgap = 0.023
    panel_w = (right - left - 2 * hgap) / 3
    panel_h = (top - bottom - 2 * vgap) / 3
    xs = [left + col * (panel_w + hgap) for col in range(3)]
    ys = [top - (row + 1) * panel_h - row * vgap for row in range(3)]

    variable_ymax = {}
    for variable, _ in VARIABLES:
        ymax = 0.0
        for system_key, _, _, _, _ in SYSTEMS:
            for stage_key, _, _, _ in STAGES:
                for curve in curves[(system_key, variable, stage_key)].values():
                    ymax = max(ymax, float(np.max(curve.y + curve.e)))
        variable_ymax[variable] = max(0.03, 1.13 * ymax)

    for col, (system_key, system_label, _, system_face, system_edge) in enumerate(SYSTEMS):
        fig.text(
            xs[col] + panel_w / 2,
            0.785,
            system_label,
            ha="center",
            va="center",
            fontsize=16.8,
            fontweight="bold",
            color="#1D4163" if system_key == "pp" else INK,
        )
        for row, (variable, variable_label) in enumerate(VARIABLES):
            cell = rounded_box(
                fig,
                [xs[col], ys[row], panel_w, panel_h],
                system_face,
                edge=system_edge,
            )
            cell.set_zorder(0)
            if col == 0:
                fig.text(
                    0.158,
                    ys[row] + panel_h / 2,
                    variable_label,
                    ha="right",
                    va="center",
                    fontsize=14.0,
                    fontweight="bold",
                    color=INK,
                )
            inner_x = xs[col] + 0.013
            inner_w = panel_w - 0.026
            lane_h = panel_h * 0.382
            lane_specs = (
                ("before_tight", ys[row] + panel_h * 0.535, "#FCFDFE", INK),
                ("after_tight", ys[row] + panel_h * 0.080, "#F0FAF4", "#24734D"),
            )
            for stage_key, lane_y, lane_face, stage_color in lane_specs:
                stage_label = next(item[1] for item in STAGES if item[0] == stage_key)
                ax = fig.add_axes([inner_x, lane_y, inner_w, lane_h])
                ax.set_facecolor(lane_face)
                for label, curve in curves[(system_key, variable, stage_key)].items():
                    add_curve(ax, curve, label, data=(label == "Data"), data_stride=1)
                ax.set_ylim(0.0, variable_ymax[variable])
                ax.grid(axis="y", color=GRID, linewidth=0.48)
                ax.tick_params(direction="in", top=True, right=True, labelsize=7.4, pad=1.5)
                for spine in ax.spines.values():
                    spine.set_color("#4A5A6D")
                    spine.set_linewidth(0.72)
                ax.text(
                    0.025,
                    0.88,
                    stage_label,
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=8.0,
                    fontweight="bold",
                    color=stage_color,
                )
                if col != 0:
                    ax.tick_params(labelleft=False)
                if stage_key == "before_tight":
                    ax.tick_params(labelbottom=False)
                elif row == len(VARIABLES) - 1:
                    ax.set_xlabel(variable_label, fontsize=8.2, labelpad=1)
                else:
                    ax.tick_params(labelbottom=False)
                if row == 0 and col == 0 and stage_key == "before_tight":
                    ax.text(
                        0.985,
                        0.88,
                        r"$\it{\bf{sPHENIX}}$ Internal",
                        transform=ax.transAxes,
                        ha="right",
                        va="top",
                        fontsize=7.7,
                    )

    readout = rounded_box(fig, [0.055, 0.027, 0.890, 0.048], "#F2F8F5", edge="#9BCAB3")
    readout.text(
        0.025,
        0.50,
        "Central and peripheral Au+Au bracket the occupancy change; each vertical pair isolates the additional tight-ID selection.",
        ha="left",
        va="center",
        fontsize=13.0,
        color="#245E45",
        transform=readout.transAxes,
    )
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / f"{output_stem}.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)
    nodes = outdir / f"{output_stem}_layout_nodes.json"
    nodes.write_text(
        json.dumps(
            {
                "canvas": {"width": 2560, "height": 1440},
                "nodes": [
                    {"id": "title", "type": "title", "x": 128, "y": 38, "w": 2300, "h": 82},
                    {"id": "legend", "type": "card", "x": 141, "y": 137, "w": 2278, "h": 101},
                    {"id": "matrix", "type": "plot_group", "x": 180, "y": 330, "w": 2250, "h": 990},
                    {"id": "readout", "type": "card", "x": 141, "y": 1330, "w": 2278, "h": 69},
                ],
                "title_axis_x": 128,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return png, nodes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--rebin-factor",
        type=int,
        default=1,
        help="Adjacent displayed-bin aggregation factor; preserves candidate fractions and errors.",
    )
    args = parser.parse_args()
    source = load_source_module()
    inputs = resolve_inputs(source)
    if args.rebin_factor < 1:
        raise ValueError("--rebin-factor must be positive")
    suffix = "tight_id_split_3x3" if args.rebin_factor == 1 else f"tight_id_split_3x3_rebin{args.rebin_factor}"
    output_stem = f"the88_pp_auau_shower_shape_{suffix}_slide"
    curves, records = build_curves(source, inputs, rebin_factor=args.rebin_factor)
    png, nodes = render(curves, args.output_dir, output_stem=output_stem)
    notes = args.output_dir / f"{output_stem}_speaker_notes.md"
    notes.write_text(
        "# Speaker notes\n\n"
        "Read each cell vertically. Both lanes pass the complete common preselection. The upper "
        "lane is before tight photon ID and the lower lane is after tight ID. Columns compare "
        "the registered p+p reference with central "
        "0--20% and peripheral 50--80% Au+Au; rows show three shower-shape variables. "
        "The p+p column contains full-stat data, photon+jet simulation, and inclusive-jet "
        "simulation. The Au+Au columns likewise compare completed data with matched "
        "Photon12+20 and Jet12+20+30+40 embedding. Every curve is "
        "integrated over 15--35 GeV. "
        "The slide is a completed historical-production comparison and does not claim THE-111 "
        "data scoring or final data/simulation closure. "
        f"Displayed adjacent histogram bins are aggregated by a factor of {args.rebin_factor}; "
        "candidate fractions and statistical uncertainties are preserved.\n",
        encoding="utf-8",
    )
    manifest = args.output_dir / f"{output_stem}_manifest.json"
    source_paths = {}
    for key in ("data", "signal", "inclusive", "pp_data", "pp_photon", "pp_inclusive"):
        path, meta = inputs[key]
        source_paths[key] = {
            "path": str(path),
            "sha256": sha256(path),
            "current_entry_id": meta.get("current_entry_id") if isinstance(meta, dict) else None,
        }
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE88_PP_AUAU_SHOWER_SHAPE_TIGHT_ID_SPLIT_3X3_SLIDE_V1",
                "created_at": datetime.now().isoformat(timespec="seconds"),
                "generator": str(THIS_FILE),
                "source_renderer": str(SOURCE_RENDERER),
                "selection": [
                    "after complete common preselection and before tight ID",
                    "after tight ID",
                ],
                "pt_range": "15-35 GeV",
                "columns": [system[1] for system in SYSTEMS],
                "stages": [stage[1] for stage in STAGES],
                "variables": [item[0] for item in VARIABLES],
                "display_binning": {
                    "original_bin_width": 0.01,
                    "adjacent_bin_aggregation_factor": args.rebin_factor,
                    "displayed_bin_width": 0.01 * args.rebin_factor,
                    "data_marker_stride_after_rebin": 1,
                    "normalization": "full finite-candidate denominator retained before bin aggregation",
                },
                "sources": source_paths,
                "pp_inclusive_source": inputs["pp_inclusive_source"],
                "curve_records": records,
                "outputs": {"png": str(png), "layout_nodes": str(nodes), "speaker_notes": str(notes)},
                "png_sha256": sha256(png),
                "limitations": [
                    "Completed historical production; not THE-111 scored data.",
                    "Candidate-fraction shape QA, not a yield or closure test.",
                    "Exact-zero boundary and under/overflow remain in normalization and provenance but are not drawn as continuous shapes.",
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(png)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
