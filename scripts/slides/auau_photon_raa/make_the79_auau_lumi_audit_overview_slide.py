#!/usr/bin/env python3
"""Build a slide-facing overview of the THE-79 AuAu GRL/lumi audit."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "dataOutput" / "auauLumiAudit" / "THE79_auau_grl_lumi_scope_20260629"
PNG_OUT = OUT_DIR / "the79_auau_lumi_audit_overview_slide.png"
MANIFEST_OUT = OUT_DIR / "the79_auau_lumi_audit_overview_slide_manifest.json"
SCRIPT_OUT = OUT_DIR / "the79_auau_lumi_audit_overview_speaker_script.md"

W, H, DPI = 2560, 1440, 200

INK = "#111827"
MUTED = "#475569"
LIGHT = "#E2E8F0"
BLUE = "#2563EB"
GREEN = "#15803D"
RED = "#B91C1C"
ORANGE = "#D97706"
PALE_BLUE = "#EFF6FF"
PALE_GREEN = "#F0FDF4"
PALE_RED = "#FEF2F2"
PALE_GRAY = "#F8FAFC"

AUDIT = {
    "production_tag": "run3auau_pro001_pcdb001_v001",
    "dst_stream": "DST_CALOFITTING",
    "grl_source": "GRLs_tanner/run3auau_pro001_pcdb001_v001_dst_calofitting_grl.list",
    "audit_date": "2026-06-29",
    "production_start_date": None,
    "grl_runs": 884,
    "dst_refs": 177720,
    "counted_files": 177678,
    "open_failures": 42,
    "zombie_files": 0,
    "no_tree_files": 0,
    "counted_entries": 17724400785,
    "trigger14_raw": 25947486603,
    "trigger14_live": 24911616043,
    "trigger14_scaled": 6731068374,
    "live_raw_pct": 96.01,
    "scaled_live_pct": 27.02,
    "scaled_raw_pct": 25.94,
    "full_grl_exposure_nb": 1.064369,
    "sigma_mbd_barns": 6.324,
}

CURRENT = {
    "paired_runs": 882,
    "matched_segments": 176577,
    "missing_runs": [68491, 72588],
}


def fmt_int(x: int) -> str:
    return f"{x:,}"


def fmt_billion(x: int, ndigits: int = 3) -> str:
    return f"{x / 1e9:.{ndigits}f}B"


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


def add_metric(ax, x, y, w, h, label, value, detail, color, face):
    add_round_rect(ax, x, y, w, h, fc=face, ec=color, lw=2.3)
    ax.text(x + 0.018, y + h - 0.034, label, transform=ax.transAxes,
            ha="left", va="top", fontsize=15.5, fontweight="bold", color=color)
    ax.text(x + 0.018, y + 0.038, value, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=25.5, fontweight="bold", color=INK)


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
        "Au+Au GRL/lumi audit: available data inventory",
        transform=canvas.transAxes,
        ha="left", va="top", fontsize=31, fontweight="bold", color=INK,
    )
    canvas.add_patch(Rectangle((0.052, 0.875), 0.895, 0.004, transform=canvas.transAxes,
                               facecolor=BLUE, edgecolor="none", alpha=0.9))
    canvas.text(
        0.055, 0.835,
        f"Production tag: {AUDIT['production_tag']} {AUDIT['dst_stream']} GRL | audit {AUDIT['audit_date']}",
        transform=canvas.transAxes,
        ha="left", va="top", fontsize=15.4, fontweight="bold", color=INK,
    )

    add_metric(canvas, 0.055, 0.635, 0.205, 0.155,
               "GRL run scope", f"{AUDIT['grl_runs']} runs",
               "current Run-3 Au+Au GRL", BLUE, PALE_BLUE)
    add_metric(canvas, 0.285, 0.635, 0.205, 0.155,
               "DST inventory", fmt_int(AUDIT["dst_refs"]),
               "unique DST_CALOFITTING refs", GREEN, PALE_GREEN)
    add_metric(canvas, 0.515, 0.635, 0.205, 0.155,
               "Readable/countable", fmt_int(AUDIT["counted_files"]),
               "99.98% of full-GRL refs", GREEN, PALE_GREEN)
    add_metric(canvas, 0.745, 0.635, 0.200, 0.155,
               "Event entries", fmt_billion(AUDIT["counted_entries"], 3),
               "counted ROOT entries", RED, PALE_RED)

    # Trigger/lumi accounting table.
    add_round_rect(canvas, 0.055, 0.255, 0.455, 0.350, fc=PALE_GRAY, ec=LIGHT, lw=2.0)
    canvas.text(0.080, 0.570, "Trigger-14 DAQ scaling", transform=canvas.transAxes,
                ha="left", va="top", fontsize=21.5, fontweight="bold", color=INK)
    canvas.text(0.080, 0.528, "Reference: MBD N&S >= 2, |z| < 150 cm", transform=canvas.transAxes,
                ha="left", va="top", fontsize=13.4, color=MUTED)
    canvas.text(0.080, 0.503, f"Normalization: $\\sigma_{{\\mathrm{{MBD}}}}$ = {AUDIT['sigma_mbd_barns']:.3f} barns",
                transform=canvas.transAxes, ha="left", va="top", fontsize=13.4, color=MUTED)

    rows = [
        ("raw counts", fmt_billion(AUDIT["trigger14_raw"])),
        ("live counts", fmt_billion(AUDIT["trigger14_live"])),
        ("scaled counts", fmt_billion(AUDIT["trigger14_scaled"])),
        ("live / raw", f"{AUDIT['live_raw_pct']:.2f}%"),
        ("scaled / raw", f"{AUDIT['scaled_raw_pct']:.2f}%"),
        ("full-GRL exposure", f"{AUDIT['full_grl_exposure_nb']:.3f} nb$^{{-1}}$"),
    ]
    y0 = 0.458
    for i, (label, value) in enumerate(rows):
        y = y0 - i * 0.036
        canvas.text(0.085, y, label, transform=canvas.transAxes,
                    ha="left", va="center", fontsize=13.8, color=MUTED)
        canvas.text(0.465, y, value, transform=canvas.transAxes,
                    ha="right", va="center", fontsize=14.8, fontweight="bold", color=INK)
        if i < len(rows) - 1:
            canvas.plot([0.080, 0.470], [y - 0.020, y - 0.020], transform=canvas.transAxes,
                        color=LIGHT, lw=1.0)

    # File QA panel.
    add_round_rect(canvas, 0.550, 0.245, 0.395, 0.360, fc="white", ec=LIGHT, lw=2.0)
    canvas.text(0.575, 0.565, "File availability QA", transform=canvas.transAxes,
                ha="left", va="top", fontsize=20, fontweight="bold", color=INK)
    canvas.text(0.575, 0.532, "Counting pass status over the full-GRL file inventory", transform=canvas.transAxes,
                ha="left", va="top", fontsize=12.4, color=MUTED)

    ax = fig.add_axes([0.590, 0.405, 0.315, 0.070])
    ax.barh([0], [AUDIT["counted_files"]], color=GREEN, height=0.45)
    ax.barh([0], [AUDIT["open_failures"]], left=[AUDIT["counted_files"]], color=ORANGE, height=0.45)
    ax.set_xlim(0, AUDIT["dst_refs"])
    ax.set_yticks([])
    ax.set_xticks([0, AUDIT["dst_refs"]])
    ax.set_xticklabels(["0", fmt_int(AUDIT["dst_refs"])], fontsize=11)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.text(0.02, 0.50, "counted 99.98%", transform=ax.transAxes,
            ha="left", va="center", fontsize=12.6, fontweight="bold", color="white")

    qa_rows = [
        ("countable ROOT files", fmt_int(AUDIT["counted_files"]), GREEN),
        ("open failures", f"{AUDIT['open_failures']}", ORANGE),
        ("zombie files", f"{AUDIT['zombie_files']}", GREEN),
        ("files with no TTree", f"{AUDIT['no_tree_files']}", GREEN),
    ]
    y = 0.350
    for label, value, color in qa_rows:
        canvas.text(0.580, y, label, transform=canvas.transAxes,
                    ha="left", va="center", fontsize=12.3, color=MUTED)
        canvas.text(0.905, y, value, transform=canvas.transAxes,
                    ha="right", va="center", fontsize=13.5, fontweight="bold", color=color)
        y -= 0.034

    # Bottom definition strip.
    add_round_rect(canvas, 0.055, 0.095, 0.890, 0.110, fc="#FFF7ED", ec=ORANGE, lw=2.0)
    canvas.text(0.080, 0.175, "Open failures: definition", transform=canvas.transAxes,
                ha="left", va="top", fontsize=17.5, fontweight="bold", color=ORANGE)
    canvas.text(
        0.080, 0.140,
        "A listed DST returned OPENFAIL in the counting wrapper: the file could not be opened by ROOT/FROG/Fun4All, so no entry or Trigger-14 count was read.\n"
        "This is separate from file corruption diagnostics: zombie files = 0 and files with no TTree = 0.",
        transform=canvas.transAxes,
        ha="left", va="top", fontsize=12.2, color=INK, linespacing=1.08,
    )

    manifest = {
        "png": str(PNG_OUT),
        "audit": AUDIT,
        "current_slide_sample_for_next_slide": CURRENT,
        "interpretation": [
            "The audit establishes the currently available full-GRL file inventory and DAQ Trigger-14 scaling.",
            "The current PPG slides use the existing THE-69 paired-output subset, handled on the follow-up scope slide.",
            "No PHENIX-like R_AA production was launched by this audit.",
        ],
        "open_failure_definition": "OPENFAIL means the listed DST could not be opened by the ROOT/FROG/Fun4All counting wrapper, so no entry or Trigger-14 count was read. It is distinct from zombie and no-TTree statuses.",
    }
    fig.savefig(PNG_OUT, dpi=DPI)
    plt.close(fig)
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    SCRIPT_OUT.write_text(
        "\n".join([
            "# Speaker script",
            "",
            "This slide is the data-inventory result from the Au+Au luminosity audit.",
            f"The audit covered the current {AUDIT['grl_runs']}-run Au+Au GRL and found {fmt_int(AUDIT['dst_refs'])} unique DST_CALOFITTING references.",
            f"The file-level count pass was essentially complete: {fmt_int(AUDIT['counted_files'])} files were countable, with {AUDIT['open_failures']} open failures and no zombie or no-tree files.",
            "Here, open failure means the listed DST could not be opened by the ROOT/FROG/Fun4All counting wrapper, so it contributed no event or trigger count.",
            f"Across the countable files, the audit found {fmt_billion(AUDIT['counted_entries'], 3)} ROOT entries.",
            f"For Trigger 14, the DAQ accounting gives {fmt_billion(AUDIT['trigger14_raw'])} raw counts, {fmt_billion(AUDIT['trigger14_live'])} live counts, and {fmt_billion(AUDIT['trigger14_scaled'])} scaled counts.",
            f"Using the MBD reference normalization, that corresponds to a full-GRL scaled exposure of {AUDIT['full_grl_exposure_nb']:.3f} inverse nanobarns.",
            "The important distinction is that this slide is the full-GRL inventory. The next slide maps that inventory to the existing THE-69 output already used in today's Au+Au slides.",
            "",
        ]),
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    print(json.dumps(make_slide(), indent=2))


if __name__ == "__main__":
    main()
