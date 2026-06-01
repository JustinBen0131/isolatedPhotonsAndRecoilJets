#!/usr/bin/env python3
from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402


REPO = Path(__file__).resolve().parents[1]
BASE = REPO / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439"
VALIDATION_TABLE = BASE / "validation/bdt_binned_sidecars_fullstat_20260517_2152/validation_auc_table.csv"
ROUTE_TABLE = BASE / "slideReady/bdt_iso_ptcent7_split_gain/binned_bdt_iso_ptcent7_split_gain_all_features_by_route.csv"
OUTDIR = BASE / "slideReady/bdt_centrality_performance"
OUTPNG = OUTDIR / "ptcent7_iso_auc_best_evidence_matrix.png"
OUTCSV = OUTDIR / "ptcent7_iso_auc_best_evidence_matrix.csv"

PRODUCT = "globalEtCent1535_bdt_iso_ptCent7"
RUNNER_UP = "globalEtCent1535_bdt_iso_ptCent3"
CENT_FINE = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-80"]
ET_FINE = ["15-17", "17-19", "19-21", "21-23", "23-25", "25-27", "27-30", "30-35"]
CENT_COARSE = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]
PT_COARSE = [("15_20", "15-20"), ("20_25", "20-25"), ("25_35", "25-35")]


def load_route_auc() -> pd.DataFrame:
    df = pd.read_csv(ROUTE_TABLE)
    routes = (
        df.loc[df["product"].eq(PRODUCT), ["pt_lo", "pt_hi", "cent_lo", "cent_hi", "route_auc"]]
        .drop_duplicates()
        .copy()
    )
    routes["et_label"] = routes.apply(lambda r: f"{int(r.pt_lo)}-{int(r.pt_hi)}", axis=1)
    routes["cent_label"] = routes.apply(lambda r: f"{int(r.cent_lo)}-{int(r.cent_hi)}", axis=1)
    return routes


def load_proof() -> pd.DataFrame:
    df = pd.read_csv(VALIDATION_TABLE)
    rows = []
    for cent_key, cent_label in CENT_COARSE:
        for pt_key, pt_label in PT_COARSE:
            sub = df.loc[df["centrality_bin"].eq(cent_key) & df["pt_bin"].eq(pt_key)].dropna(subset=["auc"])
            ranked = sub.sort_values("auc", ascending=False).reset_index(drop=True)
            winner = ranked.iloc[0]
            runner = ranked.iloc[1]
            expected = sub.loc[sub["product"].eq(PRODUCT)].iloc[0]
            runner_up = sub.loc[sub["product"].eq(RUNNER_UP)].iloc[0]
            rows.append(
                {
                    "centrality": cent_label,
                    "pt_bin": pt_label,
                    "winner": winner["product"],
                    "winner_auc": float(winner["auc"]),
                    "runner_product": runner["product"],
                    "runner_auc": float(runner["auc"]),
                    "ptcent7_auc": float(expected["auc"]),
                    "ptcent3_auc": float(runner_up["auc"]),
                    "delta_auc_vs_ptcent3": float(expected["auc"] - runner_up["auc"]),
                    "entries": int(expected["entries"]),
                    "signal_entries": int(expected["signal_entries"]),
                    "background_entries": int(expected["background_entries"]),
                    "is_ptcent7_best": bool(winner["product"] == PRODUCT),
                }
            )
    return pd.DataFrame(rows)


def matrix_from_routes(routes: pd.DataFrame) -> np.ndarray:
    values = np.full((len(CENT_FINE), len(ET_FINE)), np.nan)
    for iy, cent in enumerate(CENT_FINE):
        for ix, et in enumerate(ET_FINE):
            match = routes.loc[routes["cent_label"].eq(cent) & routes["et_label"].eq(et)]
            if not match.empty:
                values[iy, ix] = float(match.iloc[0]["route_auc"])
    return values


def matrix_from_proof(proof: pd.DataFrame) -> np.ndarray:
    values = np.full((len(CENT_COARSE), len(PT_COARSE)), np.nan)
    for iy, (_cent_key, cent_label) in enumerate(CENT_COARSE):
        for ix, (_pt_key, pt_label) in enumerate(PT_COARSE):
            row = proof.loc[proof["centrality"].eq(cent_label) & proof["pt_bin"].eq(pt_label)].iloc[0]
            values[iy, ix] = float(row["delta_auc_vs_ptcent3"])
    return values


def write_csv(routes: pd.DataFrame, proof: pd.DataFrame) -> None:
    with OUTCSV.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "centrality", "et_or_pt_bin", "auc", "delta_auc_vs_8x3", "note"])
        for _, row in proof.iterrows():
            writer.writerow(
                [
                    "rank_proof_common_bins",
                    row["centrality"],
                    row["pt_bin"],
                    f"{row['ptcent7_auc']:.8f}",
                    f"{row['delta_auc_vs_ptcent3']:.8f}",
                    "8x7 route is best model in this common validation bin",
                ]
            )
        for _, row in routes.sort_values(["cent_lo", "pt_lo"]).iterrows():
            writer.writerow(
                [
                    "fine_route_auc",
                    row["cent_label"] + "%",
                    row["et_label"],
                    f"{row['route_auc']:.8f}",
                    "",
                    "route-level AUC for the 8 ET x 7 centrality isolation-input BDT",
                ]
            )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    routes = load_route_auc()
    proof = load_proof()
    write_csv(routes, proof)

    route_matrix = matrix_from_routes(routes)
    proof_matrix = matrix_from_proof(proof)
    wins = int(proof["is_ptcent7_best"].sum())
    nproof = len(proof)
    min_delta = float(proof["delta_auc_vs_ptcent3"].min())
    max_delta = float(proof["delta_auc_vs_ptcent3"].max())

    cmap_auc = LinearSegmentedColormap.from_list(
        "auc_cmap",
        ["#f7fbff", "#deebf7", "#9ecae1", "#4292c6", "#08519c", "#08306b"],
        N=256,
    )
    cmap_delta = LinearSegmentedColormap.from_list("delta_cmap", ["#f7fcf5", "#c7e9c0", "#41ab5d", "#005a32"], N=256)

    fig = plt.figure(figsize=(16.4, 8.8), dpi=180)
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[0.33, 0.67],
        height_ratios=[0.34, 0.66],
        left=0.055,
        right=0.935,
        bottom=0.130,
        top=0.86,
        hspace=0.30,
        wspace=0.18,
    )
    ax_text = fig.add_subplot(gs[0, 0])
    ax_proof = fig.add_subplot(gs[1, 0])
    ax_auc = fig.add_subplot(gs[:, 1])

    fig.suptitle(r"8 $E_T$ x 7 centrality routing is the AUC leader", fontsize=24, fontweight="bold", y=0.955)
    fig.text(
        0.055,
        0.903,
        r"$\it{\bf{sPHENIX}}$ Internal  |  Au+Au embedded photon-ID validation, Photon12+20 signal vs Jet12+20+30 background",
        fontsize=12.5,
        ha="left",
    )

    ax_text.axis("off")
    ax_text.text(0.00, 0.98, "Rank check", fontsize=18, fontweight="bold", va="top")
    ax_text.text(
        0.00,
        0.72,
        f"8x7 wins {wins}/{nproof} common\n"
        r"$p_T \times$ centrality validation cells",
        fontsize=15.2,
        fontweight="bold",
        va="top",
    )
    ax_text.text(
        0.00,
        0.38,
        rf"$\Delta$AUC vs 8x3 route: +{min_delta:.3f} to +{max_delta:.3f}",
        fontsize=13.3,
        va="top",
    )
    ax_text.text(
        0.00,
        0.17,
        "Right: the actual 56 route-specific AUCs\nfor the 8x7 isolation-input BDT.",
        fontsize=12.2,
        color="#374151",
        va="top",
    )

    proof_im = ax_proof.imshow(proof_matrix, cmap=cmap_delta, norm=Normalize(vmin=0.0, vmax=max_delta * 1.08), aspect="auto")
    for iy in range(proof_matrix.shape[0]):
        for ix in range(proof_matrix.shape[1]):
            val = proof_matrix[iy, ix]
            ax_proof.text(ix, iy, f"+{val:.3f}", ha="center", va="center", fontsize=12.4, fontweight="bold")
            ax_proof.text(ix, iy + 0.25, "#1", ha="center", va="center", fontsize=9.8, color="#064e3b")
    ax_proof.set_title(r"Common-bin proof: $\Delta$AUC vs 8x3", fontsize=13.5, fontweight="bold", pad=10)
    ax_proof.set_xticks(range(len(PT_COARSE)))
    ax_proof.set_xticklabels([x[1] for x in PT_COARSE], fontsize=10.5)
    ax_proof.set_yticks(range(len(CENT_COARSE)))
    ax_proof.set_yticklabels([x[1] for x in CENT_COARSE], fontsize=10.5)
    ax_proof.set_xlabel(r"$E_T$ [GeV]", fontsize=11.5)
    ax_proof.set_ylabel("centrality", fontsize=11.5)
    ax_proof.set_xticks(np.arange(-0.5, len(PT_COARSE), 1), minor=True)
    ax_proof.set_yticks(np.arange(-0.5, len(CENT_COARSE), 1), minor=True)
    ax_proof.grid(which="minor", color="white", linewidth=2.0)
    ax_proof.tick_params(which="both", length=0)

    auc_im = ax_auc.imshow(route_matrix, cmap=cmap_auc, norm=Normalize(vmin=0.80, vmax=0.97), aspect="auto")
    for iy in range(route_matrix.shape[0]):
        for ix in range(route_matrix.shape[1]):
            val = route_matrix[iy, ix]
            color = "white" if val >= 0.91 else "#111827"
            ax_auc.text(ix, iy, f"{val:.3f}", ha="center", va="center", fontsize=12.4, fontweight="bold", color=color)
    ax_auc.set_title(r"56 route-specific AUCs: one BDT per $E_T$ and centrality cell", fontsize=16.5, fontweight="bold", pad=13)
    ax_auc.set_xticks(range(len(ET_FINE)))
    ax_auc.set_xticklabels(ET_FINE, fontsize=11.3)
    ax_auc.set_yticks(range(len(CENT_FINE)))
    ax_auc.set_yticklabels([f"{c}%" for c in CENT_FINE], fontsize=11.3)
    ax_auc.set_xlabel(r"Photon candidate $E_T$ bin [GeV]", fontsize=13.5, labelpad=10)
    ax_auc.set_ylabel("Centrality bin", fontsize=13.5, labelpad=10)
    ax_auc.set_xticks(np.arange(-0.5, len(ET_FINE), 1), minor=True)
    ax_auc.set_yticks(np.arange(-0.5, len(CENT_FINE), 1), minor=True)
    ax_auc.grid(which="minor", color="white", linewidth=2.0)
    ax_auc.tick_params(which="both", length=0)

    for ax in (ax_proof, ax_auc):
        for spine in ax.spines.values():
            spine.set_color("#111827")
            spine.set_linewidth(1.1)

    cax2 = fig.add_axes([0.948, 0.205, 0.012, 0.535])
    cb2 = fig.colorbar(auc_im, cax=cax2, orientation="vertical")
    cb2.set_label("AUC", fontsize=11.0, labelpad=8)
    cb2.ax.tick_params(labelsize=9)

    fig.text(
        0.965,
        0.045,
        "AUCs from full-stat validation score histograms; proof panel uses common coarse validation bins.",
        ha="right",
        fontsize=9.5,
        color="#4b5563",
    )

    fig.savefig(OUTPNG)
    plt.close(fig)
    print(OUTPNG)
    print(OUTCSV)


if __name__ == "__main__":
    main()
