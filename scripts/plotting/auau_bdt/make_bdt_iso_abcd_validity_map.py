#!/usr/bin/env python3
"""Make an ABCD-validity map for BDT/isolation interdependence."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = REPORT / "score_caches.local.list"
BDT_WP_CSV = REPO / (
    "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/"
    "the41_centdep_wp_slides_20260604/the41_centdep_bdt_wp_summary.csv"
)
ISO_FINE_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_wp_checks/"
    "auau_sliding_isolation_wp_r03_r04_eff70_80_90_fine7_no_numbers_coefficients.csv"
)
CORR_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_feature_correlations/noiso_score_diagnostic/"
    "isolation_feature_score_correlations.csv"
)
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants_abcd_validity_20260604"

SCORE = "score_ptFine_cent7"
EISO = "reco_eiso"
CENT7 = [
    (0.0, 10.0, "0-10%"),
    (10.0, 20.0, "10-20%"),
    (20.0, 30.0, "20-30%"),
    (30.0, 40.0, "30-40%"),
    (40.0, 50.0, "40-50%"),
    (50.0, 60.0, "50-60%"),
    (60.0, 80.0, "60-80%"),
]
ET8 = [
    (15.0, 17.0, "15-17"),
    (17.0, 19.0, "17-19"),
    (19.0, 21.0, "19-21"),
    (21.0, 23.0, "21-23"),
    (23.0, 25.0, "23-25"),
    (25.0, 27.0, "25-27"),
    (27.0, 30.0, "27-30"),
    (30.0, 35.0, "30-35"),
]
CENT3 = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
ET3 = [(15.0, 20.0, "15-20"), (20.0, 25.0, "20-25"), (25.0, 35.0, "25-35")]

W, H, DPI = 2560, 1440, 200
INK = "#111827"
MUTED = "#5b6472"
GREEN = "#147a4b"
RED = "#b42318"
BLUE = "#2563a6"
GOLD = "#b7791f"


def setup() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": ["Times New Roman", "DejaVu Serif"],
            "mathtext.default": "regular",
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
        }
    )


def load_arrays() -> dict[str, np.ndarray]:
    chunks: dict[str, list[np.ndarray]] = {key: [] for key in ["eiso", "score", "is_signal", "et", "cent"]}
    for line in MANIFEST.read_text().splitlines():
        if not line.strip():
            continue
        data = np.load(REPO / line.strip(), allow_pickle=True)
        chunks["eiso"].append(data[EISO].astype("float32", copy=False))
        chunks["score"].append(data[SCORE].astype("float32", copy=False))
        chunks["is_signal"].append(data["is_signal"].astype("int8", copy=False))
        chunks["et"].append(data["cluster_Et"].astype("float32", copy=False))
        chunks["cent"].append(data["centrality"].astype("float32", copy=False))
    arrays = {key: np.concatenate(vals) for key, vals in chunks.items()}
    selected = (
        np.isfinite(arrays["eiso"])
        & np.isfinite(arrays["score"])
        & (arrays["et"] >= ET3[0][0])
        & (arrays["et"] < ET3[-1][1])
        & (arrays["cent"] >= CENT3[0][0])
        & (arrays["cent"] < CENT3[-1][1])
    )
    return {key: val[selected] for key, val in arrays.items()}


def isolation_thresholds(cone: str) -> dict[tuple[float, float], float]:
    out: dict[tuple[float, float], float] = {}
    with ISO_FINE_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["cone"] == cone and abs(float(row["efficiency"]) - 0.90) < 1e-9:
                out[(float(row["cent_min"]), float(row["cent_max"]))] = float(row["threshold_gev"])
    return out


def configured_bdt_wp80() -> np.ndarray:
    out = np.full((len(CENT7), len(ET8)), np.nan)
    with BDT_WP_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["row_type"] != "et_bin_point" or row["wp_label"] != "WP80":
                continue
            cent = (float(row["centrality_min"]), float(row["centrality_max"]))
            et = (float(row["et_min"]), float(row["et_max"]))
            ci = next((i for i, (lo, hi, _) in enumerate(CENT7) if (lo, hi) == cent), None)
            ei = next((i for i, (lo, hi, _) in enumerate(ET8) if (lo, hi) == et), None)
            if ci is not None and ei is not None:
                out[ci, ei] = float(row["threshold"])
    return out


def event_cut_masks(data: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    iso = isolation_thresholds("R0.3")
    bdt = configured_bdt_wp80()
    clean = np.zeros(len(data["score"]), dtype=bool)
    bdt_pass = np.zeros(len(data["score"]), dtype=bool)
    for ci, (clo, chi, _) in enumerate(CENT7):
        cent_mask = (data["cent"] >= clo) & (data["cent"] < chi)
        clean[cent_mask] = data["eiso"][cent_mask] < iso[(clo, chi)]
        for ei, (elo, ehi, _) in enumerate(ET8):
            mask = cent_mask & (data["et"] >= elo) & (data["et"] < ehi)
            bdt_pass[mask] = data["score"][mask] >= bdt[ci, ei]
    return clean, bdt_pass


def closure_grid(data: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, dict[str, dict]]:
    clean, bdt_pass = event_cut_masks(data)
    background = data["is_signal"] == 0
    ratio = np.full((len(CENT3), len(ET3)), np.nan)
    log2_ratio = np.full_like(ratio, np.nan)
    stats: dict[str, dict] = {}
    for ci, (clo, chi, clabel) in enumerate(CENT3):
        for ei, (elo, ehi, elabel) in enumerate(ET3):
            base = background & (data["cent"] >= clo) & (data["cent"] < chi) & (data["et"] >= elo) & (data["et"] < ehi)
            # A is the signal-like ABCD corner. This is background-only MC, so
            # factorization predicts A*D/(B*C) near one.
            a = int((base & clean & bdt_pass).sum())
            b = int((base & ~clean & bdt_pass).sum())
            c = int((base & clean & ~bdt_pass).sum())
            d = int((base & ~clean & ~bdt_pass).sum())
            n = int(base.sum())
            if b > 0 and c > 0:
                raw = (a * d) / (b * c)
                ratio[ci, ei] = raw
                log2_ratio[ci, ei] = math.log(raw, 2) if raw > 0 else float("nan")
                sigma_ln = math.sqrt(sum(1.0 / max(x, 1) for x in (a, b, c, d)))
            else:
                raw = float("nan")
                sigma_ln = float("nan")
            stats[f"{clabel}_{elabel}"] = {
                "background_entries": n,
                "A_iso_pass_bdt_pass": a,
                "B_iso_fail_bdt_pass": b,
                "C_iso_pass_bdt_fail": c,
                "D_iso_fail_bdt_fail": d,
                "abcd_ratio_AD_over_BC": raw,
                "abcd_log2_ratio": float(log2_ratio[ci, ei]) if np.isfinite(log2_ratio[ci, ei]) else None,
                "sigma_ln_ratio_approx": sigma_ln if math.isfinite(sigma_ln) else None,
            }
    return ratio, log2_ratio, stats


def scope_key(cent_label: str, et_label: str) -> str:
    if et_label == "15-17":
        return f"{cent_label}, 15-17 GeV"
    return f"{cent_label}, {et_label} GeV"


def load_background_correlations() -> tuple[dict[str, dict], dict[str, dict]]:
    rows = []
    with CORR_CSV.open() as handle:
        for row in csv.DictReader(handle):
            if row["class"] != "background":
                continue
            if row["target_variable"] != "score_globalEtCent1535_bdt_noIso_ptCent7":
                continue
            if row["iso_variable"] not in {"reco_eiso_r30", "reco_eiso_r40"}:
                continue
            rows.append(row)
    by_scope: dict[tuple[str, str], dict[str, float]] = {}
    for row in rows:
        by_scope[(row["scope_label"], row["iso_variable"])] = {
            "pearson": float(row["pearson"]),
            "n_entries": int(float(row["n_entries"])),
        }

    rho30 = np.full((len(CENT3), len(ET3)), np.nan)
    rho40 = np.full_like(rho30, np.nan)
    stats: dict[str, dict] = {}
    et8_lookup = {label: (lo, hi) for lo, hi, label in ET8}
    for ci, (_, _, clabel) in enumerate(CENT3):
        for ei, (elo3, ehi3, elabel3) in enumerate(ET3):
            key = f"{clabel}_{elabel3}"
            stats[key] = {}
            for iso_var, target in [("reco_eiso_r30", rho30), ("reco_eiso_r40", rho40)]:
                vals = []
                weights = []
                for label, (elo, ehi) in et8_lookup.items():
                    if elo >= elo3 and ehi <= ehi3:
                        item = by_scope.get((scope_key(clabel, label), iso_var))
                        if item is not None:
                            vals.append(item["pearson"])
                            weights.append(item["n_entries"])
                if vals and weights:
                    avg = float(np.average(np.asarray(vals), weights=np.asarray(weights)))
                    target[ci, ei] = avg
                    stats[key][iso_var] = {
                        "weighted_pearson": avg,
                        "source_rows": len(vals),
                        "weighted_entries": int(sum(weights)),
                    }
    return {"rho30": rho30, "rho40": rho40}, stats


def style_heatmap_axis(ax, title: str, show_y: bool = True) -> None:
    ax.set_title(title, fontsize=17.0, fontweight="bold", pad=10)
    ax.set_xticks(np.arange(len(ET3)))
    ax.set_xticklabels([label for _, _, label in ET3], fontsize=12.0)
    ax.set_yticks(np.arange(len(CENT3)))
    ax.set_yticklabels([label for _, _, label in CENT3] if show_y else [], fontsize=12.0)
    ax.set_xticks(np.arange(-0.5, len(ET3), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(CENT3), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_color("#d7dee8")


def cell_text_color(val: float, norm) -> str:
    if not np.isfinite(val):
        return MUTED
    scaled = abs(norm(val) - 0.5) * 2.0
    return "white" if scaled > 0.55 else INK


def make_plot(data: dict[str, np.ndarray]) -> tuple[Path, dict]:
    ratio, log2_ratio, closure_stats = closure_grid(data)
    corrs, corr_stats = load_background_correlations()
    rho30 = corrs["rho30"]
    rho40 = corrs["rho40"]

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    gs = fig.add_gridspec(3, 4, left=0.070, right=0.950, top=0.800, bottom=0.145, wspace=0.24, hspace=0.42)
    ax_closure = fig.add_subplot(gs[:2, :2])
    ax_corr = fig.add_subplot(gs[:2, 2:])
    ax_cartoon = fig.add_subplot(gs[2, 0])
    ax_read = fig.add_subplot(gs[2, 1:])

    fig.text(0.045, 0.955, "ABCD validity check: BDT vs isolation factorization", fontsize=26.0, fontweight="bold", ha="left", va="top")
    fig.text(
        0.045,
        0.902,
        r"Background MC test. Independence means $R_{ABCD}=A D/(B C)$ stays near 1 and the BDT-isolation correlation stays near 0.",
        fontsize=13.7,
        color=MUTED,
        ha="left",
        va="top",
    )

    norm_ratio = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    ax_closure.imshow(log2_ratio, cmap="RdYlGn_r", norm=norm_ratio, aspect="auto")
    style_heatmap_axis(ax_closure, r"Closure ratio: $R_{ABCD}=A D/(B C)$")
    ax_closure.set_ylabel("centrality", fontsize=13.0)
    for i in range(len(CENT3)):
        for j in range(len(ET3)):
            raw = ratio[i, j]
            key = f"{CENT3[i][2]}_{ET3[j][2]}"
            n = closure_stats[key]["background_entries"]
            text = f"{raw:.2f}\nn={n}" if np.isfinite(raw) else f"--\nn={n}"
            ax_closure.text(j, i, text, ha="center", va="center", fontsize=15.0, fontweight="bold", color=cell_text_color(log2_ratio[i, j], norm_ratio))

    norm_rho = TwoSlopeNorm(vmin=-0.35, vcenter=0.0, vmax=0.35)
    ax_corr.imshow(rho30, cmap="RdBu_r", norm=norm_rho, aspect="auto")
    style_heatmap_axis(ax_corr, r"Residual dependence: background Pearson $\rho$")
    ax_corr.set_xlabel(r"cluster $E_T$ bin [GeV]", fontsize=13.0)
    for i in range(len(CENT3)):
        for j in range(len(ET3)):
            text = f"R0.3 {rho30[i, j]:+.2f}\nR0.4 {rho40[i, j]:+.2f}"
            ax_corr.text(j, i, text, ha="center", va="center", fontsize=14.0, fontweight="bold", color=cell_text_color(rho30[i, j], norm_rho))

    ax_cartoon.axis("off")
    ax_cartoon.set_title("ABCD corners", fontsize=14.0, fontweight="bold", pad=5)
    ax_cartoon.set_xlim(0, 1)
    ax_cartoon.set_ylim(0, 1)
    ax_cartoon.plot([0.12, 0.92], [0.50, 0.50], color=INK, lw=1.8)
    ax_cartoon.plot([0.52, 0.52], [0.10, 0.90], color=INK, lw=1.8)
    ax_cartoon.text(0.30, 0.74, "A\niso pass\nBDT pass", ha="center", va="center", fontsize=10.5, color=GREEN, fontweight="bold")
    ax_cartoon.text(0.73, 0.74, "B\niso fail\nBDT pass", ha="center", va="center", fontsize=10.5)
    ax_cartoon.text(0.30, 0.29, "C\niso pass\nBDT fail", ha="center", va="center", fontsize=10.5)
    ax_cartoon.text(0.73, 0.29, "D\niso fail\nBDT fail", ha="center", va="center", fontsize=10.5)

    ax_read.axis("off")
    ax_read.set_xlim(0, 1)
    ax_read.set_ylim(0, 1)
    readouts = [
        ("How to read it", "left panel near 1 means the ABCD product factorizes in background MC", INK),
        ("What we see", "background has a persistent negative BDT-isolation correlation; closure is not uniformly flat", RED),
        ("Use carefully", "R0.4 closure needs row-level R0.4 isolation; local row cache only has default reco Eiso", GOLD),
    ]
    y = 0.82
    for head, body, color in readouts:
        ax_read.text(0.00, y, head, fontsize=14.0, fontweight="bold", color=color, ha="left", va="top")
        ax_read.text(0.20, y, body, fontsize=13.0, color=INK if color != RED else RED, ha="left", va="top")
        y -= 0.275

    fig.text(
        0.070,
        0.063,
        r"Closure: local row-level background candidates, default/R0.3 reco $E_T^{iso}$ WP90, configured THE-41 BDT WP80.",
        fontsize=10.4,
        color=MUTED,
        ha="left",
        va="bottom",
    )
    fig.text(
        0.070,
        0.037,
        r"Correlation: full-stat no-isolation BDT aggregate CSV; broad $E_T$ cells are entry-weighted averages of finer bins.",
        fontsize=10.4,
        color=MUTED,
        ha="left",
        va="bottom",
    )

    out = OUT_DIR / "01_abcd_validity_bdt_isolation_interdependence.png"
    fig.savefig(out)
    plt.close(fig)
    stats = {
        "closure": closure_stats,
        "correlations": corr_stats,
        "closure_ratio_range": [
            float(np.nanmin(ratio)),
            float(np.nanmax(ratio)),
        ],
        "rho30_range": [
            float(np.nanmin(rho30)),
            float(np.nanmax(rho30)),
        ],
        "rho40_range": [
            float(np.nanmin(rho40)),
            float(np.nanmax(rho40)),
        ],
    }
    return out, stats


def main() -> None:
    setup()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_arrays()
    png, stats = make_plot(data)
    manifest = {
        "schema": "THE11_BDT_ISO_ABCD_VALIDITY_V1",
        "outputs": {"abcd_validity_map": str(png)},
        "sources": {
            "row_level_cache_manifest": str(MANIFEST),
            "configured_bdt_wp80_csv": str(BDT_WP_CSV),
            "isolation_wp_csv": str(ISO_FINE_CSV),
            "full_stat_noiso_correlation_csv": str(CORR_CSV),
        },
        "columns": {"score": SCORE, "isolation": EISO, "truth_label": "is_signal"},
        "definitions": {
            "A": "isolation pass and BDT pass",
            "B": "isolation fail and BDT pass",
            "C": "isolation pass and BDT fail",
            "D": "isolation fail and BDT fail",
            "closure_ratio": "A*D/(B*C), computed on background-labeled inclusive candidates",
        },
        "notes": [
            "This is an ABCD validity/interdependence diagnostic, not a candidate-yield plot.",
            "Closure uses local row-level default reco_eiso/R0.3 because reco_eiso_r40 is absent from the local row cache.",
            "BDT pass/fail uses configured THE-41 WP80 thresholds on the event's fine centrality x ET bin.",
            "Correlation values use the full-stat no-isolation BDT aggregate CSV; broad ET cells are entry-weighted averages of finer bins.",
            "Google Slides was not mutated.",
        ],
        "stats": stats,
    }
    manifest_path = OUT_DIR / "abcd_validity_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
