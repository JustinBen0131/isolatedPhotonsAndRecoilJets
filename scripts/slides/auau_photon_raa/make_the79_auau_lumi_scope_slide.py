#!/usr/bin/env python3
"""Build a slide-facing AuAu GRL/lumi scope summary for THE-79."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "dataOutput" / "auauLumiAudit" / "THE79_auau_grl_lumi_scope_20260629"
PNG_OUT = OUT_DIR / "the79_auau_grl_lumi_scope_slide.png"
MANIFEST_OUT = OUT_DIR / "the79_auau_grl_lumi_scope_slide_manifest.json"
SCRIPT_OUT = OUT_DIR / "the79_auau_grl_lumi_scope_speaker_script.md"

W, H, DPI = 2560, 1440, 200

INK = "#111827"
MUTED = "#475569"
LIGHT = "#E2E8F0"
BLUE = "#2563EB"
RED = "#B91C1C"
GREEN = "#15803D"
ORANGE = "#D97706"
PALE_BLUE = "#EFF6FF"
PALE_RED = "#FEF2F2"
PALE_GREEN = "#F0FDF4"

AUDIT = {
    "grl_runs": 884,
    "unique_dst_roots": 177720,
    "counted_files": 177678,
    "open_failures": 42,
    "zombie_files": 0,
    "no_tree_files": 0,
    "counted_entries": 17724400785,
    "trigger14_scaled_exposure_nb": 1.064369,
}

CURRENT = {
    "paired_runs": 882,
    "matched_segments": 176577,
    "missing_runs": [68491, 72588],
    "cluster": "5552123",
    "raw_roots": 25599,
    "non_tiny_roots": 25481,
    "final_root_size_mb": 136.496,
}


def pct(num: float, den: float) -> str:
    return f"{100.0 * num / den:.2f}%"


def fmt_int(x: int) -> str:
    return f"{x:,}"


def add_round_rect(ax, x, y, w, h, fc="white", ec=LIGHT, lw=2.0, radius=0.018):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.010,rounding_size={radius}",
        transform=ax.transAxes,
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    return patch


def add_card(ax, x, y, w, h, title, value, detail, color, face):
    add_round_rect(ax, x, y, w, h, fc=face, ec=color, lw=2.2)
    ax.text(x + 0.020, y + h - 0.033, title, transform=ax.transAxes,
            ha="left", va="top", fontsize=15.0, fontweight="bold", color=color)
    ax.text(x + 0.020, y + 0.050, value, transform=ax.transAxes,
            ha="left", va="center", fontsize=25.5, fontweight="bold", color=INK)


def make_slide() -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
    })

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.patch.set_facecolor("white")
    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.axis("off")

    canvas.text(
        0.050, 0.948,
        "Au+Au data scope for the following plots",
        transform=canvas.transAxes,
        ha="left", va="top", fontsize=27, fontweight="bold", color=INK,
    )

    canvas.text(
        0.052, 0.892,
        "The following Au+Au plots use the analyzed 882-run paired-output subset.\n"
        "The previous slide audited the full 884-run production GRL; this subset is near-complete and from the same physics sample.",
        transform=canvas.transAxes,
        ha="left", va="top", fontsize=14.6, color=INK, linespacing=1.14,
    )

    add_card(
        canvas, 0.055, 0.645, 0.275, 0.145,
        "Plotted subset",
        "882 / 884 runs",
        "paired streams\n176,577 segment pairs",
        BLUE, PALE_BLUE,
    )
    add_card(
        canvas, 0.363, 0.645, 0.275, 0.145,
        "Full GRL audit",
        "177,678 files",
        "countable DST files\n17.724B entries",
        GREEN, PALE_GREEN,
    )
    add_card(
        canvas, 0.671, 0.645, 0.275, 0.145,
        "Trigger-14 scale",
        "1.064 nb$^{-1}$",
        "full-GRL DAQ scale\nnot official for slides yet",
        RED, PALE_RED,
    )

    # Coverage bars.
    ax = fig.add_axes([0.155, 0.215, 0.410, 0.355])
    labels = ["runs", "segments", "GRL files counted"]
    vals = [
        CURRENT["paired_runs"] / AUDIT["grl_runs"],
        CURRENT["matched_segments"] / AUDIT["unique_dst_roots"],
        AUDIT["counted_files"] / AUDIT["unique_dst_roots"],
    ]
    colors = [BLUE, ORANGE, GREEN]
    ypos = list(range(len(labels)))
    ax.barh(ypos, [100 * v for v in vals], color=colors, height=0.52)
    ax.set_xlim(98.8, 100.05)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=17)
    ax.set_xlabel("coverage relative to the audited 884-run GRL [%]", fontsize=15.0)
    ax.tick_params(axis="x", labelsize=13.0)
    ax.tick_params(axis="y", labelsize=14.8)
    ax.grid(axis="x", color="#CBD5E1", linewidth=1.1, alpha=0.75)
    ax.set_axisbelow(True)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#94A3B8")
    ax.spines["bottom"].set_color("#94A3B8")
    for y, v in zip(ypos, vals):
        ax.text(100 * v - 0.015, y, f"{100*v:.2f}%", ha="right", va="center",
                fontsize=15.5, fontweight="bold", color="white")
    ax.invert_yaxis()
    ax.set_title("Plotted subset and audit completeness vs full GRL", fontsize=18.0, fontweight="bold", pad=11)

    # Accounting box.
    add_round_rect(canvas, 0.615, 0.215, 0.330, 0.355, fc="white", ec=LIGHT, lw=2.0)
    canvas.text(0.635, 0.540, "Accounting notes", transform=canvas.transAxes,
                ha="left", va="top", fontsize=20.5, fontweight="bold", color=INK)
    rows = [
        f"Analyzed output: {CURRENT['raw_roots']:,} ROOT files;\n  {CURRENT['non_tiny_roots']:,} usable non-tiny files.",
        f"Missing paired GRL runs: {', '.join(str(r) for r in CURRENT['missing_runs'])}.",
        f"Full-GRL audit refs: {AUDIT['unique_dst_roots']:,}; open failures: {AUDIT['open_failures']}.",
    ]
    y = 0.490
    for row in rows:
        canvas.text(0.640, y, u"\u2022 " + row, transform=canvas.transAxes, ha="left", va="top",
                    fontsize=12.8, color=INK, linespacing=1.10)
        y -= 0.070
    canvas.text(
        0.640, 0.225,
        f"Following plots correspond to\n"
        f"approx. {AUDIT['trigger14_scaled_exposure_nb'] * CURRENT['matched_segments'] / AUDIT['unique_dst_roots']:.3f} nb$^{{-1}}$.",
        transform=canvas.transAxes,
        ha="left", va="bottom", fontsize=16.0, color=RED, linespacing=1.10,
    )

    # Thin top separator.
    canvas.add_patch(Rectangle((0.052, 0.825), 0.895, 0.004, transform=canvas.transAxes,
                               facecolor=BLUE, edgecolor="none", alpha=0.9))

    manifest = {
        "png": str(PNG_OUT),
        "current_slide_sample": CURRENT,
        "full_grl_audit": AUDIT,
        "derived": {
            "run_fraction": CURRENT["paired_runs"] / AUDIT["grl_runs"],
            "segment_fraction": CURRENT["matched_segments"] / AUDIT["unique_dst_roots"],
            "countable_file_fraction": AUDIT["counted_files"] / AUDIT["unique_dst_roots"],
            "missing_runs": AUDIT["grl_runs"] - CURRENT["paired_runs"],
            "missing_segment_or_file_refs": AUDIT["unique_dst_roots"] - CURRENT["matched_segments"],
            "segment_scaled_exposure_nb_not_official": AUDIT["trigger14_scaled_exposure_nb"] * CURRENT["matched_segments"] / AUDIT["unique_dst_roots"],
        },
        "caveats": [
            "Full-GRL Trigger-14 exposure comes from the audit thirdPass DAQ full-GRL scaled exposure.",
            "The current-list projection mapped zero current-list runs despite healthy file-level count evidence; therefore current-output luminosity is not finalized.",
            "The GRL files counted bar is the file-counting completeness of the full-GRL audit, not a current-output coverage statement.",
            "The plotted AuAu subset is the current paired-output sample, not the full 884-run GRL projection.",
        ],
    }
    fig.savefig(PNG_OUT, dpi=DPI)
    plt.close(fig)
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    SCRIPT_OUT.write_text(
        "\n".join([
            "# Speaker script",
            "",
            "This slide tells the audience exactly what data subset is used in the following Au+Au plots.",
            "Those plots use the current paired-output sample: 882 of the 884 GRL runs, with 176,577 matched CALOFITTING and ZDC segment pairs.",
            "The new audit looked at the full 884-run GRL and counted 177,678 readable DST_CALOFITTING ROOT files, corresponding to 17.724 billion counted entries.",
            "So this is not a different data period. It is essentially the same production family, with the audit covering the two missing GRL runs and about 1,143 additional file references.",
            "The full-GRL DAQ Trigger-14 scaled exposure from the audit is 1.064 inverse nanobarns.",
            "For the plotted 882-run subset, the segment-scaled exposure estimate is about 1.058 inverse nanobarns.",
            "The action item is to fix that projection, then freeze the exposure table before using it in an R_AA normalization.",
            "",
        ]),
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    print(json.dumps(make_slide(), indent=2))


if __name__ == "__main__":
    main()
