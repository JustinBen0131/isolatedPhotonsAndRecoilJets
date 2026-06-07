#!/usr/bin/env python3
"""Build THE-11 BDT/isolation insight full-slide PNG candidates."""

from __future__ import annotations

import csv
import json
import textwrap
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402


REPO = Path(__file__).resolve().parents[4]
SOURCE_DIR = REPO / (
    "dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_feature_correlations/noiso_score_diagnostic"
)
CSV_PATH = SOURCE_DIR / "isolation_feature_score_correlations.csv"
SUMMARY_PATH = SOURCE_DIR / "isolation_feature_score_correlations_summary.json"
OUT_DIR = SOURCE_DIR / "the11_bdt_isolation_story_20260604"

W, H = 2560, 1440
DPI = 200

BG = "#fbfaf7"
INK = "#171717"
MUTED = "#5b6472"
TEAL = "#1f9d8a"
BLUE = "#315f9c"
RED = "#c84c4c"
GOLD = "#d79b2d"
PURPLE = "#6f4aa8"
GRAY = "#e9edf2"
SOFT_YELLOW = "#fff3c7"

SCORE_COL = "score_globalEtCent1535_bdt_noIso_ptCent7"
ISO_LABEL = {
    "reco_eiso_r30": r"$R=0.3\ E_T^{iso}$",
    "reco_eiso_r40": r"$R=0.4\ E_T^{iso}$",
}
CLASS_LABEL = {
    "all": "All candidates",
    "signal": "Truth photons",
    "background": "Jet background",
}
CLASS_COLOR = {
    "all": INK,
    "signal": RED,
    "background": BLUE,
}
FAMILY_ORDER = [
    "BDT score",
    "kinematic/context",
    "shower widths",
    "width ratios",
    "tower energy sharing",
    "E11 core ratios",
    "containment ratios",
]


@dataclass
class CorrRow:
    scope: str
    scope_label: str
    scope_type: str
    klass: str
    class_label: str
    iso_variable: str
    target_variable: str
    target_label: str
    feature_family: str
    pearson: float
    spearman: float
    n_entries: int
    spearman_n_entries: int


def setup_style() -> None:
    candidates = ["Times New Roman", "Times", "DejaVu Serif"]
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in candidates:
        if name in available:
            plt.rcParams["font.family"] = name
            break
    plt.rcParams.update(
        {
            "figure.facecolor": BG,
            "axes.facecolor": "white",
            "axes.edgecolor": "#2a2a2a",
            "axes.linewidth": 1.2,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "axes.labelcolor": INK,
            "savefig.facecolor": BG,
        }
    )


def load_rows() -> list[CorrRow]:
    rows: list[CorrRow] = []
    with CSV_PATH.open() as f:
        for row in csv.DictReader(f):
            rows.append(
                CorrRow(
                    scope=row["scope"],
                    scope_label=row["scope_label"],
                    scope_type=row["scope_type"],
                    klass=row["class"],
                    class_label=row["class_label"],
                    iso_variable=row["iso_variable"],
                    target_variable=row["target_variable"],
                    target_label=row["target_label"],
                    feature_family=row["feature_family"],
                    pearson=float(row["pearson"]),
                    spearman=float(row["spearman"]),
                    n_entries=int(float(row["n_entries"])),
                    spearman_n_entries=int(float(row["spearman_n_entries"])),
                )
            )
    return rows


def filter_rows(rows: list[CorrRow], **kwargs: str) -> list[CorrRow]:
    out = rows
    for key, value in kwargs.items():
        attr = "klass" if key == "class" else key
        out = [r for r in out if getattr(r, attr) == value]
    return out


def get_one(rows: list[CorrRow], **kwargs: str) -> CorrRow:
    out = filter_rows(rows, **kwargs)
    if len(out) != 1:
        raise RuntimeError(f"Expected one row for {kwargs}, got {len(out)}")
    return out[0]


def add_title(fig, title: str, subtitle: str) -> None:
    fig.text(0.045, 0.935, title, ha="left", va="top", fontsize=30, fontweight="bold")
    fig.text(0.045, 0.895, subtitle, ha="left", va="top", fontsize=16.5, color=MUTED)
    fig.add_artist(Rectangle((0.045, 0.866), 0.91, 0.003, transform=fig.transFigure, color="#d7dce3", lw=0))


def rounded_box(fig, xywh, facecolor="white", edgecolor="#d4dae3", lw=1.5, radius=0.018):
    patch = FancyBboxPatch(
        (xywh[0], xywh[1]),
        xywh[2],
        xywh[3],
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        transform=fig.transFigure,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        zorder=0.5,
    )
    fig.add_artist(patch)
    return patch


def add_card(fig, x, y, w, h, label, body, color, value=None):
    rounded_box(fig, (x, y, w, h), facecolor="white", edgecolor="#d8dfe8")
    fig.text(x + 0.018, y + h - 0.030, label, fontsize=16.5, fontweight="bold", color=color, va="top")
    body_wrapped = textwrap.fill(body, width=max(22, int(w * 95)))
    if value is not None:
        fig.text(x + 0.018, y + h - 0.076, value, fontsize=29, fontweight="bold", color=color, va="top")
        body_y = y + h - 0.122
    else:
        body_y = y + h - 0.064
    fig.text(x + 0.018, body_y, body_wrapped, fontsize=12.8, color=INK, va="top", linespacing=1.12)


def save_script(name: str, title: str, paragraphs: list[str]) -> Path:
    path = OUT_DIR / f"{name}_script.md"
    path.write_text("# THE-11 BDT Isolation Insight - " + title + "\n\n" + "\n\n".join(paragraphs) + "\n")
    return path


def slide_score_insight(rows: list[CorrRow]) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_title(
        fig,
        "The no-isolation BDT learns nearby-activity structure",
        "Raw cone isolation is not an input, but the score is still ordered by isolation-like shower context.",
    )

    ax = fig.add_axes([0.08, 0.28, 0.54, 0.45])
    classes = ["all", "signal", "background"]
    x = np.arange(len(classes))
    width = 0.34
    vals30 = [get_one(rows, scope="inclusive", klass=c, iso_variable="reco_eiso_r30", target_variable=SCORE_COL).pearson for c in classes]
    vals40 = [get_one(rows, scope="inclusive", klass=c, iso_variable="reco_eiso_r40", target_variable=SCORE_COL).pearson for c in classes]
    ax.axhline(0, color="#5e6670", lw=1)
    bars1 = ax.bar(x - width / 2, vals30, width=width, color=TEAL, label=r"$R=0.3\ E_T^{iso}$")
    bars2 = ax.bar(x + width / 2, vals40, width=width, color=GOLD, label=r"$R=0.4\ E_T^{iso}$")
    for bars in [bars1, bars2]:
        for b in bars:
            v = b.get_height()
            ax.text(b.get_x() + b.get_width() / 2, v - 0.025 if v < 0 else v + 0.015, f"{v:+.2f}", ha="center", va="top" if v < 0 else "bottom", fontsize=17, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([CLASS_LABEL[c] for c in classes], fontsize=16)
    ax.set_ylabel("Pearson correlation with BDT score", fontsize=17)
    ax.set_ylim(-0.52, 0.08)
    ax.tick_params(axis="y", labelsize=16)
    ax.grid(axis="y", color="#dfe4ea", lw=1, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", fontsize=14, frameon=False, ncol=2)
    ax.set_title("Score correlation with raw isolation", fontsize=20, fontweight="bold", pad=12)

    add_card(
        fig,
        0.67,
        0.635,
        0.27,
        0.145,
        "All candidates",
        "The full pool has a strong negative ordering with raw isolation.",
        INK,
        r"$\rho=-0.42$",
    )
    add_card(
        fig,
        0.67,
        0.435,
        0.27,
        0.145,
        "Truth photons",
        "For signal alone, the relation nearly vanishes.",
        RED,
        r"$\rho=-0.03$",
    )
    add_card(
        fig,
        0.67,
        0.235,
        0.27,
        0.145,
        "Jet background",
        "For background, the score still falls as nearby activity grows.",
        BLUE,
        r"$\rho=-0.28$",
    )
    rounded_box(fig, (0.08, 0.075, 0.86, 0.095), facecolor=SOFT_YELLOW, edgecolor="#eadba0")
    fig.text(
        0.105,
        0.134,
        "Interpretation:",
        fontsize=19,
        fontweight="bold",
        color=INK,
        va="center",
    )
    fig.text(
        0.245,
        0.134,
        "Not an isolation cut: high-isolation background looks less photon-like through shower/context inputs.",
        fontsize=15.5,
        color=INK,
        va="center",
    )
    out = OUT_DIR / "the11_slide01_noiso_bdt_discovers_nearby_activity.png"
    fig.savefig(out)
    plt.close(fig)
    script = save_script(
        "the11_slide01_noiso_bdt_discovers_nearby_activity",
        "Slide 1 Script - No-Isolation BDT Discovers Nearby Activity",
        [
            "First, the important point is that this BDT score is the no-isolation score. Raw cone isolation is not one of the inputs to the headline model.",
            "Even so, when I compare the score to raw isolation, the full candidate sample has a strong negative correlation. For R equals 0.3 isolation, the Pearson correlation with the BDT score is about minus 0.42.",
            "The split by truth class is the key insight. For truth photons alone the relationship is almost gone, around minus 0.03. For jet background, it is still clearly negative, around minus 0.28.",
            "So the model is not simply learning an isolation variable. It is learning shower and context structure that makes high-nearby-activity background less photon-like, while preserving a much weaker dependence for real photons.",
        ],
    )
    return out, script


def family_matrix(rows: list[CorrRow]) -> tuple[list[str], list[str], np.ndarray]:
    cols = [
        ("reco_eiso_r30", "all"),
        ("reco_eiso_r30", "signal"),
        ("reco_eiso_r30", "background"),
        ("reco_eiso_r40", "all"),
        ("reco_eiso_r40", "signal"),
        ("reco_eiso_r40", "background"),
    ]
    labels = [
        r"$R=0.3$ all",
        r"$R=0.3$ signal",
        r"$R=0.3$ background",
        r"$R=0.4$ all",
        r"$R=0.4$ signal",
        r"$R=0.4$ background",
    ]
    matrix = []
    for fam in FAMILY_ORDER:
        row = []
        for iso, klass in cols:
            subset = [
                r
                for r in rows
                if r.scope == "inclusive"
                and r.iso_variable == iso
                and r.klass == klass
                and r.feature_family == fam
            ]
            if not subset:
                row.append(np.nan)
            else:
                row.append(max(abs(r.pearson) for r in subset))
        matrix.append(row)
    return FAMILY_ORDER, labels, np.array(matrix)


def slide_family_insight(rows: list[CorrRow]) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_title(
        fig,
        "Isolation lives in shower containment and widths",
        "The strongest links are localized in specific input families; signal-only correlations remain small.",
    )
    families, labels, matrix = family_matrix(rows)

    ax = fig.add_axes([0.16, 0.29, 0.53, 0.43])
    im = ax.imshow(matrix, cmap="viridis", vmin=0, vmax=0.45, aspect="auto")
    ax.set_yticks(np.arange(len(families)))
    ax.set_yticklabels(families, fontsize=14.5)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=12.5, rotation=24, ha="right")
    ax.set_title("Max absolute Pearson correlation in each feature family", fontsize=19, fontweight="bold", pad=12)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            color = "white" if val > 0.25 else "#111827"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=13.5, color=color, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.018)
    cb.ax.tick_params(labelsize=14)
    cb.set_label(r"max $|\rho|$", fontsize=16)

    add_card(
        fig,
        0.75,
        0.60,
        0.17,
        0.15,
        "Containment",
        r"$E_{22}/E_{37}$ gives the strongest R=0.3 all-candidate link.",
        PURPLE,
        r"$|\rho|=0.36$",
    )
    add_card(
        fig,
        0.75,
        0.405,
        0.17,
        0.15,
        "Shower widths",
        r"Variables such as $w_\phi$ also track raw isolation.",
        TEAL,
        r"$|\rho|=0.28$",
    )
    add_card(
        fig,
        0.75,
        0.21,
        0.17,
        0.15,
        "Width ratios",
        r"Shape-ratio families stay nearly decoupled.",
        GOLD,
        r"$|\rho|=0.01$",
    )

    rounded_box(fig, (0.09, 0.075, 0.84, 0.09), facecolor=SOFT_YELLOW, edgecolor="#eadba0")
    fig.text(0.112, 0.133, "Model insight:", fontsize=18, fontweight="bold", va="center")
    fig.text(
        0.235,
        0.133,
        "Nearby activity appears through wider, less-contained background-like clusters.",
        fontsize=15.5,
        va="center",
    )

    out = OUT_DIR / "the11_slide02_isolation_entangled_with_shower_families.png"
    fig.savefig(out)
    plt.close(fig)
    script = save_script(
        "the11_slide02_isolation_entangled_with_shower_families",
        "Slide 2 Script - Isolation Entangled With Shower Families",
        [
            "The next question is where the isolation-like behavior is coming from, since isolation itself is not in the BDT input list.",
            "This heatmap compresses the 32 input variables into feature families. Each cell shows the strongest absolute Pearson correlation between raw isolation and any variable in that family.",
            "The largest effects are not spread evenly across the input space. They are concentrated in containment ratios and shower-width variables. For example, the containment family reaches about 0.36 in all candidates for R equals 0.3 isolation, and shower widths reach about 0.28.",
            "At the same time, the signal-only columns are mostly near zero. That is important because it says the model's isolation-like ordering is mainly a background-structure effect, not a strong distortion of the prompt-photon signal sample.",
        ],
    )
    return out, script


def slide_radius_phase_space(rows: list[CorrRow]) -> tuple[Path, Path]:
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_title(
        fig,
        "R=0.3 and R=0.4 measure one nearby-activity axis",
        "Across centrality and cluster energy, the two raw cone definitions are highly redundant.",
    )

    cent_labels = ["0-20%", "20-50%", "50-80%"]
    et_labels = ["15-17", "17-19", "19-21", "21-23", "23-25", "25-27", "27-30", "30-35"]
    classes = ["all", "signal", "background"]
    axes = []
    for k, klass in enumerate(classes):
        ax = fig.add_axes([0.07 + k * 0.295, 0.43, 0.245, 0.30])
        mat = np.full((3, 8), np.nan)
        for r in rows:
            if r.scope_type != "centrality_et" or r.klass != klass:
                continue
            if r.iso_variable != "reco_eiso_r30" or r.target_variable != "reco_eiso_r40":
                continue
            for i, cent in enumerate(cent_labels):
                for j, et in enumerate(et_labels):
                    if r.scope_label == f"{cent}, {et} GeV":
                        mat[i, j] = r.pearson
        im = ax.imshow(mat, cmap="magma", vmin=0.80, vmax=1.00, aspect="auto")
        ax.set_title(CLASS_LABEL[klass], fontsize=18, fontweight="bold", pad=8)
        ax.set_xticks(np.arange(len(et_labels)))
        ax.set_xticklabels(et_labels, fontsize=11, rotation=38, ha="right")
        ax.set_yticks(np.arange(len(cent_labels)))
        ax.set_yticklabels(cent_labels if k == 0 else ["", "", ""], fontsize=14)
        if k == 0:
            ax.set_ylabel("centrality", fontsize=15)
        ax.set_xlabel(r"cluster $E_T$ bin [GeV]", fontsize=14)
        for i in range(3):
            for j in range(8):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=9.5, color="white", fontweight="bold")
        axes.append(ax)
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.018)
    cb.ax.tick_params(labelsize=14)
    cb.set_label(r"Pearson correlation", fontsize=16)

    add_card(
        fig,
        0.08,
        0.205,
        0.25,
        0.105,
        "All candidates",
        "0.86-0.96 across cells.",
        INK,
    )
    add_card(
        fig,
        0.375,
        0.205,
        0.25,
        0.105,
        "Signal",
        "0.82-0.85 across cells.",
        RED,
    )
    add_card(
        fig,
        0.67,
        0.205,
        0.25,
        0.105,
        "Background",
        "0.88-0.94 across cells.",
        BLUE,
    )

    rounded_box(fig, (0.08, 0.065, 0.84, 0.09), facecolor=SOFT_YELLOW, edgecolor="#eadba0")
    fig.text(0.103, 0.123, "Implication:", fontsize=18, fontweight="bold", va="center")
    fig.text(
        0.215,
        0.123,
        "The cone radius is less important than the shared nearby-activity axis it measures.",
        fontsize=15.5,
        va="center",
    )

    out = OUT_DIR / "the11_slide03_isolation_radius_stable_axis.png"
    fig.savefig(out)
    plt.close(fig)
    script = save_script(
        "the11_slide03_isolation_radius_stable_axis",
        "Slide 3 Script - Isolation Radius Stable Axis",
        [
            "Finally, I checked whether R equals 0.3 and R equals 0.4 isolation are really telling us different stories.",
            "Across centrality and cluster-energy bins, the two isolation definitions are highly correlated. For all candidates the cell-by-cell Pearson correlation stays between about 0.86 and 0.96. For signal it is about 0.82 to 0.85, and for background it is about 0.88 to 0.94.",
            "So for a qualitative model explanation, I would not over-interpret the cone-radius choice. Both radii are measuring a common nearby-activity axis, and the BDT's relationship to isolation is really about that axis interacting with shower containment and width structure.",
            "That gives us a cleaner story: the current BDT is not an opaque score. It learns a photon-like shower core and penalizes background-like nearby activity indirectly through the shape variables it was allowed to use.",
        ],
    )
    return out, script


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    pngs_scripts = [
        slide_score_insight(rows),
        slide_family_insight(rows),
        slide_radius_phase_space(rows),
    ]
    summary = json.loads(SUMMARY_PATH.read_text())
    manifest = {
        "schema": "THE11_BDT_ISOLATION_INSIGHT_SLIDES_V1",
        "source_csv": str(CSV_PATH.relative_to(REPO)),
        "source_summary": str(SUMMARY_PATH.relative_to(REPO)),
        "source_rows": len(rows),
        "entries_seen": summary.get("entries_seen"),
        "files_read": summary.get("files_read"),
        "files_requested": summary.get("files_requested"),
        "missing_column_files": summary.get("missing_column_files"),
        "score_columns": summary.get("score_columns"),
        "selection_notes": summary.get("notes"),
        "validation_caveat": (
            "Targeted no-isolation rescore produced all 125/125 score caches but "
            "the formal validation merge-summary node failed; these slides use "
            "the local row-aligned ABCD-safe compact outputs."
        ),
        "google_slides_mutation": False,
        "outputs": [
            {
                "png": str(png.relative_to(REPO)),
                "script": str(script.relative_to(REPO)),
            }
            for png, script in pngs_scripts
        ],
    }
    (OUT_DIR / "the11_bdt_isolation_story_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    plan = OUT_DIR / "the11_bdt_isolation_story_plan.md"
    plan.write_text(
        "# THE-11 BDT Isolation Insight Slide Plan\n\n"
        "1. The no-isolation BDT still discovers nearby activity: compare raw isolation against the no-isolation BDT score for all, signal, and background rows.\n"
        "2. Isolation is entangled with specific shower families: show that containment ratios and shower widths carry the strongest isolation correlations, while signal-only correlations stay small.\n"
        "3. R=0.3 and R=0.4 isolation form one stable nearby-activity axis: show high cell-by-cell correlation across centrality and cluster E_T.\n\n"
        "No Google Slides mutation was performed. These are local full-slide PNG candidates for Justin review.\n"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
