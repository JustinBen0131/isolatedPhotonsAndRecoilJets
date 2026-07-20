#!/usr/bin/env python3
"""Build a JSTG-ready THE-100 Au+Au shower-shape data/MC slide.

The script reads only the required histograms from the completed THE-100
data, photon-embedding, and inclusive-jet-embedding ROOT files on SDCC.  The
large ROOT files are never copied locally; a compact JSON extraction is
returned over the authenticated SSH connection and preserved beside the PNG.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import textwrap

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/slides/the45_jstg_20260720/the100_data_mc_shower_shapes"
)

ROOTS = {
    "data": (
        "/sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/"
        "the100_auau_dualview_20260714/auau/"
        "RecoilJets_auau_ALL_preselectionNewPPG12_"
        "tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root"
    ),
    "signal": (
        "/sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/"
        "the100_auau_dualview_20260714/signal/simembedded/"
        "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
        "nonTightAuAuBDTComplement_baseVariant/photonJet12and20merged_SIM/"
        "RecoilJets_embeddedPhoton12plus20_MERGED.root"
    ),
    "inclusive": (
        "/sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/"
        "the100_auau_dualview_20260714/inclusive/simembeddedinclusive/"
        "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
        "nonTightAuAuBDTComplement_baseVariant/embeddedJet12and20and30and40merged_SIM/"
        "RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
    ),
}

PT_GROUPS = (
    (r"Low $p_T^{\gamma}$: 15--21 GeV", ((15, 17), (17, 19), (19, 21))),
    (r"Mid $p_T^{\gamma}$: 21--26 GeV", ((21, 23), (23, 26))),
    (r"High $p_T^{\gamma}$: 26--35 GeV", ((26, 35),)),
)

VARIABLES = (
    {
        "key": "et1",
        "label": r"$E_1/E_{\mathrm{cluster}}$",
        "feature": "cluster_et1",
        "xlim": (0.25, 1.02),
        "rebin": 2,
    },
    {
        "key": "e11e33",
        "label": r"$E_{1\times1}/E_{3\times3}$",
        "feature": "e11_over_e33",
        "xlim": (0.00, 0.96),
        "rebin": 2,
    },
    {
        "key": "e32e35",
        "label": r"$E_{3\times2}/E_{3\times5}$",
        "feature": "e32_over_e35",
        "xlim": (0.45, 1.02),
        "rebin": 2,
    },
)

STAGES = (
    ("Before tight ID", "pre"),
    ("After tight ID", "tight"),
)
INTEGRATED_PT_LABEL = r"$15<E_T^{\gamma}<35~\mathrm{GeV}$"

WP80 = {"intercept": 0.5387310379, "slope": 0.0011102647}

INK = "#132238"
MUTED = "#506279"
GRID = "#DCE5EF"
RED = "#D53A32"
BLUE = "#2468D8"
BLACK = "#111827"
GREEN = "#159B6A"
TOP_BG = "#F6F8FB"
TIGHT_BG = "#EEF9F4"
PT_STYLES = (
    {"face": "#F1F7FF", "edge": "#92BDEC", "accent": "#2563EB"},
    {"face": "#FFF8E9", "edge": "#E3C67C", "accent": "#AD7118"},
    {"face": "#F7F1FF", "edge": "#C7A9ED", "accent": "#7C3AED"},
)


@dataclass(frozen=True)
class Curve:
    role: str
    variable: str
    pt_group: str
    stage: str
    edges: list[float]
    fraction: list[float]
    error: list[float]
    visible_integral: float
    full_integral: float
    displayed_fraction: float
    fill_multiplicity: int
    effective_entries: float
    source_keys: list[str]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.edgecolor": INK,
            "axes.linewidth": 0.9,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def rounded_box(fig, xywh, face, edge="#CAD6E4", radius=0.014, lw=1.0):
    ax = fig.add_axes(xywh)
    ax.axis("off")
    patch = FancyBboxPatch(
        (0, 0),
        1,
        1,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(patch)
    return ax


def _remote_extraction_source(roles: tuple[str, ...] = ("data", "signal", "inclusive")) -> str:
    contract = {
        "roots": {role: ROOTS[role] for role in roles},
        "pt_groups": [
            {"label": label, "bins": [list(x) for x in bins]}
            for label, bins in PT_GROUPS
        ],
        "variables": [
            {
                "key": item["key"],
                "label": item["label"],
                "xlim": list(item["xlim"]),
                "rebin": item["rebin"],
            }
            for item in VARIABLES
        ],
        "stages": [list(x) for x in STAGES],
    }
    payload = json.dumps(contract)
    return textwrap.dedent(
        f"""
        import json
        import os
        import ROOT

        ROOT.gROOT.SetBatch(True)
        ROOT.TH1.AddDirectory(False)
        contract = json.loads({payload!r})

        def hist_name(role, var, stage, lo, hi):
            if role == "data":
                return (
                    "photon_12_plus_MBD_NS_geq_2_vtx_lt_150/"
                    f"h_ss_{{var}}_{{stage}}_pT_{{lo}}_{{hi}}_cent_0_20"
                )
            suffix = "sig" if role == "signal" else "bkg"
            return f"SIM/h_ss_{{var}}_{{stage}}_{{suffix}}_pT_{{lo}}_{{hi}}_cent_0_20"

        def add_group(root_file, role, variable, stage, bins):
            total = None
            keys = []
            for lo, hi in bins:
                name = hist_name(role, variable["key"], stage, lo, hi)
                hist = root_file.Get(name)
                if not hist:
                    raise RuntimeError(f"missing ROOT object: {{name}}")
                keys.append(name)
                if total is None:
                    total = hist.Clone(f"sum_{{role}}_{{variable['key']}}_{{stage}}_{{lo}}_{{hi}}")
                    total.SetDirectory(0)
                else:
                    total.Add(hist)
            return total, keys

        def extract(hist, xlim, rebin, fill_multiplicity):
            if rebin > 1:
                hist.Rebin(rebin)
            nb = hist.GetNbinsX()
            full = float(hist.Integral(0, nb + 1))
            selected = []
            # Bin 1 is the exact/near-zero unresolved boundary component in
            # the persisted Au+Au shower-shape families.  Preserve it in the
            # denominator, but do not draw it as part of the continuous shape.
            for i in range(2, nb + 1):
                center = float(hist.GetXaxis().GetBinCenter(i))
                if xlim[0] <= center <= xlim[1]:
                    selected.append(i)
            if not selected:
                raise RuntimeError("display range contains no bins")
            first, last = selected[0], selected[-1]
            visible = float(hist.Integral(first, last))
            edges = [float(hist.GetXaxis().GetBinLowEdge(i)) for i in range(first, last + 1)]
            edges.append(float(hist.GetXaxis().GetBinUpEdge(last)))
            counts = [float(hist.GetBinContent(i)) for i in range(first, last + 1)]
            errors = [float(hist.GetBinError(i)) for i in range(first, last + 1)]
            if full > 0:
                fraction = [x / full for x in counts]
                error = [(x * fill_multiplicity**0.5) / full for x in errors]
            else:
                fraction = [0.0 for x in counts]
                error = [0.0 for x in errors]
            return {{
                "edges": edges,
                "fraction": fraction,
                "error": error,
                "visible_integral": visible,
                "full_integral": full,
                "displayed_fraction": visible / full if full > 0 else 0.0,
                "fill_multiplicity": fill_multiplicity,
                "effective_entries": float(hist.GetEntries()) / fill_multiplicity,
            }}

        output = {{
            "schema": "THE100_SHOWER_SHAPE_DATA_MC_EXTRACTION_V1",
            "campaign": "the100_auau_dualview_20260714",
            "created_at_remote": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
            "files": {{}},
            "curves": [],
        }}

        for role, path in contract["roots"].items():
            root_file = ROOT.TFile.Open(path, "READ")
            if not root_file or root_file.IsZombie():
                raise RuntimeError(f"cannot open ROOT file: {{path}}")
            if root_file.TestBit(ROOT.TFile.kRecovered):
                raise RuntimeError(f"recovered ROOT file is not accepted: {{path}}")
            output["files"][role] = {{
                "path": path,
                "bytes": os.path.getsize(path),
                "mtime": os.path.getmtime(path),
                "zombie": False,
                "recovered": False,
            }}
            try:
                for variable in contract["variables"]:
                    for pt_group in contract["pt_groups"]:
                        for stage_label, stage_key in contract["stages"]:
                            hist, keys = add_group(
                                root_file,
                                role,
                                variable,
                                stage_key,
                                pt_group["bins"],
                            )
                            fill_multiplicity = 4 if stage_key in ("pre", "tight") else 1
                            curve = extract(
                                hist,
                                variable["xlim"],
                                int(variable["rebin"]),
                                fill_multiplicity,
                            )
                            curve.update({{
                                "role": role,
                                "variable": variable["key"],
                                "pt_group": pt_group["label"],
                                "stage": stage_label,
                                "source_keys": keys,
                            }})
                            output["curves"].append(curve)
            finally:
                root_file.Close()

        print("THE100_JSON_BEGIN")
        print(json.dumps(output, separators=(",", ":")))
        print("THE100_JSON_END")
        """
    )


def extract_from_sdcc(
    timeout_seconds: int = 300,
    roles: tuple[str, ...] = ("data", "signal", "inclusive"),
) -> dict[str, object]:
    sock = subprocess.check_output(
        ["launchctl", "getenv", "SSH_AUTH_SOCK"], text=True
    ).strip()
    if not sock:
        raise RuntimeError("macOS SSH agent socket is unavailable")
    env = os.environ.copy()
    env["SSH_AUTH_SOCK"] = sock
    command = [
        "ssh",
        "-o",
        "BatchMode=yes",
        "patsfan753@ssh.sdcc.bnl.gov",
        (
            "ssh -o BatchMode=yes -o StrictHostKeyChecking=no "
            "-o UserKnownHostsFile=/dev/null "
            "sphnxuser05.sdcc.bnl.gov python3 -"
        ),
    ]
    completed = subprocess.run(
        command,
        input=_remote_extraction_source(roles),
        text=True,
        capture_output=True,
        env=env,
        timeout=timeout_seconds,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "SDCC extraction failed with code "
            f"{completed.returncode}:\n{completed.stderr[-4000:]}"
        )
    begin = completed.stdout.find("THE100_JSON_BEGIN")
    end = completed.stdout.find("THE100_JSON_END")
    if begin < 0 or end < 0 or end <= begin:
        raise RuntimeError(
            "SDCC extraction returned no marked JSON payload; tail follows:\n"
            + completed.stdout[-4000:]
        )
    text = completed.stdout[begin + len("THE100_JSON_BEGIN") : end].strip()
    payload = json.loads(text)
    expected_curves = 18 * len(roles)
    if len(payload.get("curves", [])) != expected_curves:
        raise RuntimeError(
            f"Expected {expected_curves} curves, found {len(payload.get('curves', []))}"
        )
    return payload


def curves_from_payload(payload: dict[str, object]) -> list[Curve]:
    return [Curve(**row) for row in payload["curves"]]


def lookup(curves: list[Curve]) -> dict[tuple[str, str, str, str], Curve]:
    return {(c.role, c.variable, c.pt_group, c.stage): c for c in curves}


def combine_pt_groups(curves: list[Curve]) -> list[Curve]:
    """Combine the disjoint stored pT groups before normalization."""
    grouped = lookup(curves)
    combined: list[Curve] = []
    source_labels = [label for label, _ in PT_GROUPS]
    for role in ("data", "signal", "inclusive"):
        for variable in VARIABLES:
            for stage, _ in STAGES:
                parts = [
                    grouped[(role, variable["key"], label, stage)]
                    for label in source_labels
                ]
                reference_edges = np.asarray(parts[0].edges, dtype=float)
                if any(
                    not np.allclose(reference_edges, np.asarray(part.edges, dtype=float))
                    for part in parts[1:]
                ):
                    raise RuntimeError(
                        f"incompatible edges for {role}/{variable['key']}/{stage}"
                    )
                full = float(sum(part.full_integral for part in parts))
                visible = float(sum(part.visible_integral for part in parts))
                counts = np.sum(
                    [np.asarray(part.fraction) * part.full_integral for part in parts],
                    axis=0,
                )
                absolute_errors = np.sqrt(
                    np.sum(
                        [
                            (np.asarray(part.error) * part.full_integral) ** 2
                            for part in parts
                        ],
                        axis=0,
                    )
                )
                multiplicities = {part.fill_multiplicity for part in parts}
                if len(multiplicities) != 1:
                    raise RuntimeError(
                        f"mixed fill multiplicity for {role}/{variable['key']}/{stage}"
                    )
                combined.append(
                    Curve(
                        role=role,
                        variable=variable["key"],
                        pt_group=INTEGRATED_PT_LABEL,
                        stage=stage,
                        edges=reference_edges.tolist(),
                        fraction=(counts / full).tolist() if full > 0 else [0.0] * len(counts),
                        error=(absolute_errors / full).tolist()
                        if full > 0
                        else [0.0] * len(counts),
                        visible_integral=visible,
                        full_integral=full,
                        displayed_fraction=visible / full if full > 0 else 0.0,
                        fill_multiplicity=multiplicities.pop(),
                        effective_entries=float(sum(part.effective_entries for part in parts)),
                        source_keys=[key for part in parts for key in part.source_keys],
                    )
                )
    return combined


def draw_panel(ax, lk, variable, pt_label, style) -> None:
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color(style["edge"])
        spine.set_linewidth(1.05)
    ax.axhspan(1.00, 1.44, color=TOP_BG, zorder=0)
    ax.axhspan(0.10, 0.54, color=TIGHT_BG, zorder=0)
    ax.hlines([1.00, 0.10], variable["xlim"][0], variable["xlim"][1], color=GRID, lw=0.8)

    raw_max = 0.0
    for role in ("data", "signal", "inclusive"):
        for stage, _ in STAGES:
            curve = lk[(role, variable["key"], pt_label, stage)]
            if curve.fraction:
                raw_max = max(raw_max, max(curve.fraction))
    scale = 0.34 / raw_max if raw_max > 0 else 1.0

    lanes = (("Before tight ID", 1.00), ("After tight ID", 0.10))
    for stage, base in lanes:
        for role, color, linewidth in (
            ("inclusive", BLUE, 2.1),
            ("signal", RED, 2.1),
        ):
            curve = lk[(role, variable["key"], pt_label, stage)]
            edges = np.asarray(curve.edges)
            values = base + scale * np.asarray(curve.fraction)
            ax.stairs(values, edges, color=color, linewidth=linewidth, zorder=3)

        data = lk[("data", variable["key"], pt_label, stage)]
        edges = np.asarray(data.edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        y = base + scale * np.asarray(data.fraction)
        yerr = scale * np.asarray(data.error)
        stride = max(1, len(centers) // 22)
        ax.errorbar(
            centers[::stride],
            y[::stride],
            yerr=yerr[::stride],
            fmt="o",
            color=BLACK,
            markerfacecolor=BLACK,
            markeredgecolor="white",
            markeredgewidth=0.45,
            markersize=3.1,
            elinewidth=0.65,
            capsize=0,
            zorder=5,
        )

    x0, x1 = variable["xlim"]
    label_x = x0 + 0.018 * (x1 - x0)
    ax.text(label_x, 1.42, "Before tight ID", ha="left", va="top", fontsize=12.8,
            fontweight="bold", color="#475569")
    ax.text(label_x, 0.52, "After tight ID", ha="left", va="top", fontsize=12.8,
            fontweight="bold", color=GREEN)
    ax.set_xlim(x0, x1)
    ax.set_ylim(0.05, 1.47)
    ax.grid(axis="x", color=GRID, linewidth=0.55, alpha=0.85)
    ax.tick_params(labelsize=9.0, direction="in", top=True, right=True)
    ax.set_yticks([])


def render(payload: dict[str, object], outdir: Path) -> dict[str, Path]:
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)
    curves = curves_from_payload(payload)
    lk = lookup(curves)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.text(
        0.050,
        0.952,
        r"Tight photon ID reshapes data and simulation toward compact shower cores",
        ha="left",
        va="top",
        fontsize=29.0,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.051,
        0.901,
        r"Au+Au $\sqrt{s_{NN}}=200$ GeV; preselected candidates in complete data and matched embedding",
        ha="left",
        va="top",
        fontsize=17.8,
        color=MUTED,
    )

    legend_ax = rounded_box(fig, [0.055, 0.805, 0.890, 0.055], "#F8FAFC")
    legend_ax.legend(
        handles=[
            Line2D([0], [0], color=BLACK, marker="o", markersize=5.5, lw=0, label="Au+Au data"),
            Line2D([0], [0], color=RED, lw=2.8, label="Photon embedding"),
            Line2D([0], [0], color=BLUE, lw=2.8, label="Inclusive-jet embedding"),
        ],
        loc="center left",
        bbox_to_anchor=(0.03, 0.50),
        ncol=3,
        frameon=False,
        fontsize=13.2,
        handlelength=2.6,
        columnspacing=2.2,
    )
    legend_ax.text(
        0.965,
        0.50,
        r"$T_{80}(c)=0.5387+0.001110\,c$",
        ha="right",
        va="center",
        fontsize=13.6,
        color=GREEN,
        fontweight="bold",
        transform=legend_ax.transAxes,
    )

    left = 0.165
    panel_w = 0.245
    panel_h = 0.178
    hgap = 0.030
    vgap = 0.035
    top = 0.745
    xs = [left + i * (panel_w + hgap) for i in range(3)]
    ys = [top - panel_h - i * (panel_h + vgap) for i in range(3)]

    for col, (pt_label, _) in enumerate(PT_GROUPS):
        style = PT_STYLES[col]
        rail = rounded_box(
            fig,
            [xs[col] - 0.010, ys[-1] - 0.015, panel_w + 0.020, top - ys[-1] + 0.025],
            style["face"],
            edge=style["edge"],
            radius=0.012,
            lw=0.9,
        )
        rail.patches[0].set_alpha(0.26)
        header = rounded_box(
            fig,
            [xs[col] + 0.018, top + 0.005, panel_w - 0.036, 0.035],
            style["face"],
            edge=style["edge"],
            radius=0.010,
            lw=1.0,
        )
        header.text(
            0.50,
            0.52,
            pt_label,
            ha="center",
            va="center",
            fontsize=14.0,
            fontweight="bold",
            color=style["accent"],
            transform=header.transAxes,
        )

    for row, variable in enumerate(VARIABLES):
        label = rounded_box(
            fig,
            [0.047, ys[row] + 0.022, 0.095, panel_h - 0.044],
            "white",
            edge="#CAD6E4",
            radius=0.012,
        )
        label.text(
            0.50,
            0.64,
            variable["label"],
            ha="center",
            va="center",
            fontsize=17.0,
            fontweight="bold",
            color=INK,
            transform=label.transAxes,
        )
        label.text(
            0.50,
            0.32,
            variable["feature"],
            ha="center",
            va="center",
            fontsize=10.3,
            color=MUTED,
            transform=label.transAxes,
        )
        for col, (pt_label, _) in enumerate(PT_GROUPS):
            ax = fig.add_axes([xs[col], ys[row], panel_w, panel_h])
            draw_panel(ax, lk, variable, pt_label, PT_STYLES[col])
            if row == len(VARIABLES) - 1:
                ax.set_xlabel(variable["label"], fontsize=11.5, labelpad=2)
            else:
                ax.tick_params(labelbottom=False)

    takeaway = rounded_box(fig, [0.055, 0.029, 0.890, 0.052], "#FFF7D6", edge="#E6B933")
    takeaway.text(
        0.025,
        0.50,
        "After the common preselection, tight ID further concentrates the retained candidates in "
        "the compact photon-like region.",
        ha="left",
        va="center",
        fontsize=13.8,
        color="#76510A",
        transform=takeaway.transAxes,
    )

    png = outdir / "the100_energy_sum_feature_data_mc_slide.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)

    notes = outdir / "the100_energy_sum_feature_data_mc_speaker_notes.md"
    notes.write_text(
        "# Speaker notes\n\n"
        "This is the completed THE-100 production baseline, not the new THE-111 model. "
        "Black points are Au+Au data; red is embedded Photon12+20; blue is embedded "
        "Jet12+20+30+40. The upper lane is after the complete NCB/common preselection and "
        "before the tight BDT requirement; the lower lane is after tight ID. Visible bins are "
        "normalized to the full finite-candidate population, so exact-zero, underflow, and "
        "overflow components remain in the denominator without being promoted into the "
        "continuous curve. The three columns use exact stored "
        "histogram bins: 15--21, 21--26, and 26--35 GeV. The data are a signal/background "
        "mixture and are not expected to equal either truth lane. The audience-facing point is "
        "that the tight selection moves the retained data and inclusive-jet population toward "
        "the compact photon-like region. Follow this slide with the THE-111 simulation-validation "
        "slide; do not imply that THE-111 has already been run on data.\n",
        encoding="utf-8",
    )

    nodes = outdir / "the100_energy_sum_feature_data_mc_layout_nodes.json"
    nodes.write_text(
        json.dumps(
            {
                "canvas": {"width": 2560, "height": 1440},
                "nodes": [
                    {"id": "title", "type": "title", "x": 128, "y": 38, "w": 2300, "h": 82},
                    {"id": "legend", "type": "card", "x": 141, "y": 1159, "w": 2278, "h": 79},
                    {"id": "plot_group", "type": "plot_group", "x": 120, "y": 220, "w": 2295, "h": 1115},
                    {"id": "takeaway", "type": "card", "x": 141, "y": 1324, "w": 2278, "h": 75},
                ],
                "title_axis_x": 128,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return {"png": png, "notes": notes, "layout_nodes": nodes}


def render_integrated(payload: dict[str, object], outdir: Path) -> dict[str, Path]:
    """Render one readable 15--35 GeV comparison instead of three pT columns."""
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)
    curves = combine_pt_groups(curves_from_payload(payload))
    lk = lookup(curves)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.text(
        0.050,
        0.952,
        r"Tight photon ID selects compact electromagnetic shower cores",
        ha="left",
        va="top",
        fontsize=29.0,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.051,
        0.901,
        r"Au+Au $\sqrt{s_{NN}}=200$ GeV; $15<E_T^{\gamma}<35$ GeV; data and matched embedding",
        ha="left",
        va="top",
        fontsize=17.8,
        color=MUTED,
    )

    legend_ax = rounded_box(fig, [0.055, 0.805, 0.890, 0.055], "#F8FAFC")
    legend_ax.legend(
        handles=[
            Line2D([0], [0], color=BLACK, marker="o", markersize=5.5, lw=0, label="Au+Au data"),
            Line2D([0], [0], color=RED, lw=2.8, label="Photon embedding"),
            Line2D([0], [0], color=BLUE, lw=2.8, label="Inclusive-jet embedding"),
        ],
        loc="center left",
        bbox_to_anchor=(0.03, 0.50),
        ncol=3,
        frameon=False,
        fontsize=13.2,
        handlelength=2.6,
        columnspacing=2.2,
    )
    legend_ax.text(
        0.965,
        0.50,
        r"$T_{80}(c)=0.5387+0.001110\,c$",
        ha="right",
        va="center",
        fontsize=13.6,
        color=GREEN,
        fontweight="bold",
        transform=legend_ax.transAxes,
    )

    panel_x = 0.210
    panel_w = 0.735
    panel_h = 0.190
    panel_ys = [0.585, 0.355, 0.125]
    label_x = 0.055
    label_w = 0.130
    for row, variable in enumerate(VARIABLES):
        label = rounded_box(
            fig,
            [label_x, panel_ys[row] + 0.016, label_w, panel_h - 0.032],
            "#FFFFFF",
            edge="#CAD6E4",
            radius=0.010,
            lw=1.0,
        )
        label.text(
            0.50,
            0.62,
            variable["label"],
            ha="center",
            va="center",
            fontsize=17.5,
            color=INK,
            transform=label.transAxes,
        )
        label.text(
            0.50,
            0.30,
            variable["feature"],
            ha="center",
            va="center",
            fontsize=10.8,
            color=MUTED,
            transform=label.transAxes,
        )
        ax = fig.add_axes([panel_x, panel_ys[row], panel_w, panel_h])
        draw_panel(ax, lk, variable, INTEGRATED_PT_LABEL, PT_STYLES[row])
        ax.tick_params(labelsize=11.0)
        ax.set_xlabel(variable["label"], fontsize=13.0, labelpad=2)

    takeaway = rounded_box(fig, [0.055, 0.028, 0.890, 0.052], "#FFF7D6", edge="#E6B933")
    takeaway.text(
        0.025,
        0.50,
        "Integrated over 15--35 GeV, the data retain a smooth shower-shape population and "
        "tight ID selects the compact photon-like region.",
        ha="left",
        va="center",
        fontsize=13.8,
        color="#76510A",
        transform=takeaway.transAxes,
    )

    png = outdir / "the100_energy_sum_feature_data_mc_15to35_slide.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)

    notes = outdir / "the100_energy_sum_feature_data_mc_15to35_speaker_notes.md"
    notes.write_text(
        "# Speaker notes\n\n"
        "This integrated view sums the exact stored 15--21, 21--26, and 26--35 GeV "
        "histograms before normalization. It therefore represents one 15--35 GeV candidate "
        "population rather than three independently normalized slices. Black points are "
        "completed THE-100 Au+Au data, red is embedded Photon12+20, and blue is embedded "
        "Jet12+20+30+40. The upper lane is after the complete NCB/common preselection and "
        "before tight BDT selection; the lower lane is after tight ID. Exact-zero boundary, "
        "underflow, and overflow components remain in the normalization denominator and in "
        "the JSON provenance but are not drawn as a continuous shower-shape curve. This is a "
        "historical completed-production baseline and is not THE-111 scored data.\n",
        encoding="utf-8",
    )

    nodes = outdir / "the100_energy_sum_feature_data_mc_15to35_layout_nodes.json"
    nodes.write_text(
        json.dumps(
            {
                "canvas": {"width": 2560, "height": 1440},
                "nodes": [
                    {"id": "title", "type": "title", "x": 128, "y": 38, "w": 2300, "h": 82},
                    {"id": "legend", "type": "card", "x": 141, "y": 202, "w": 2278, "h": 79},
                    {"id": "plot_group", "type": "plot_group", "x": 140, "y": 320, "w": 2278, "h": 1035},
                    {"id": "takeaway", "type": "card", "x": 141, "y": 1324, "w": 2278, "h": 75},
                ],
                "title_axis_x": 128,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return {"png": png, "notes": notes, "layout_nodes": nodes}


def write_integrated_manifest(
    payload: dict[str, object], outputs: dict[str, Path], outdir: Path
) -> Path:
    manifest = outdir / "the100_energy_sum_feature_data_mc_15to35_manifest.json"
    combined = combine_pt_groups(curves_from_payload(payload))
    combined_json = outdir / "the100_energy_sum_feature_data_mc_15to35_curves.json"
    combined_json.write_text(
        json.dumps([curve.__dict__ for curve in combined], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE100_JSTG_SHOWER_SHAPE_15TO35_SLIDE_MANIFEST_V1",
                "created_at": datetime.now().isoformat(timespec="seconds"),
                "campaign": "the100_auau_dualview_20260714",
                "status": "completed historical production baseline; not THE-111 data scoring",
                "generator": str(THIS_FILE),
                "sources": payload["files"],
                "selection": "complete NCB/common preselection compared with tight ID",
                "data_trigger_namespace": "photon_12_plus_MBD_NS_geq_2_vtx_lt_150",
                "pt_contract": "stored 15-21, 21-26, and 26-35 GeV groups summed before normalization",
                "normalization": "visible continuous bins divided by the full 15-35 GeV candidate population",
                "curves_json": str(combined_json.resolve()),
                "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
                "png_sha256": hashlib.sha256(outputs["png"].read_bytes()).hexdigest(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


def write_manifest(
    payload: dict[str, object], outputs: dict[str, Path], outdir: Path
) -> Path:
    manifest = outdir / "the100_energy_sum_feature_data_mc_manifest.json"
    curves_json = outdir / "the100_energy_sum_feature_data_mc_curves.json"
    curves_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    png_sha = hashlib.sha256(outputs["png"].read_bytes()).hexdigest()
    manifest_payload = {
        "schema": "THE100_JSTG_SHOWER_SHAPE_SLIDE_MANIFEST_V1",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "campaign": "the100_auau_dualview_20260714",
        "status": "completed production baseline; historical model and reconstruction contract",
        "generator": str(THIS_FILE),
        "sources": payload["files"],
        "sample_contract": {
            "data": "completed Au+Au data production",
            "data_trigger_namespace": "photon_12_plus_MBD_NS_geq_2_vtx_lt_150",
            "signal": "embedded Photon12+20",
            "inclusive": "embedded Jet12+20+30+40",
            "centrality": "0-20%",
            "normalization": "visible continuous bins divided by the full finite-candidate population; exact-zero boundary and under/overflow remain in the denominator",
        },
        "working_point": {
            "expression": "T80(c) = 0.5387310379 + 0.0011102647 c",
            **WP80,
        },
        "pt_groups": [
            {"label": label, "fine_bins": [list(x) for x in bins]}
            for label, bins in PT_GROUPS
        ],
        "limitations": [
            "THE-100 is the completed data-inclusive production baseline; it is not THE-111.",
            "Data are a signal/background mixture, not a truth-labeled class.",
            "The plotted comparison is candidate-fraction shape QA and is not a yield or closure test.",
            "Zero, underflow, overflow, and excluded-range fractions remain recorded in the curves JSON.",
        ],
        "outputs": {key: str(path.resolve()) for key, path in outputs.items()},
        "curves_json": str(curves_json.resolve()),
        "png_sha256": png_sha,
    }
    manifest.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--curves-json", type=Path, default=None)
    parser.add_argument("--timeout-seconds", type=int, default=300)
    parser.add_argument(
        "--refresh-data-from-sdcc",
        action="store_true",
        help="Replace only the cached data curves with a fresh SDCC extraction.",
    )
    parser.add_argument(
        "--integrated-15to35",
        action="store_true",
        help="Render one integrated 15--35 GeV comparison from the stored groups.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.curves_json:
        payload = json.loads(args.curves_json.read_text(encoding="utf-8"))
        if args.refresh_data_from_sdcc:
            refreshed = extract_from_sdcc(args.timeout_seconds, roles=("data",))
            payload["curves"] = [
                row for row in payload["curves"] if row["role"] != "data"
            ] + refreshed["curves"]
            payload["files"]["data"] = refreshed["files"]["data"]
    else:
        if args.refresh_data_from_sdcc:
            raise RuntimeError("--refresh-data-from-sdcc requires --curves-json")
        payload = extract_from_sdcc(args.timeout_seconds)
    if args.integrated_15to35:
        outputs = render_integrated(payload, args.output_dir)
        manifest = write_integrated_manifest(payload, outputs, args.output_dir)
    else:
        outputs = render(payload, args.output_dir)
        manifest = write_manifest(payload, outputs, args.output_dir)
    print(outputs["png"])
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
