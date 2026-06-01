#!/usr/bin/env python3
"""Correlation diagnostics between raw cone isolation and BDT inputs/scores.

This script intentionally works from validation score caches, not ROOT payloads.
It is meant for compact post-validation diagnostics: read existing `.npz` score
cache files, compute isolation-vs-feature correlations, and write CSV/PNG/JSON
artifacts suitable for slide review.
"""

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

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


PT_EDGES = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0]
CENT_EDGES = [0.0, 20.0, 50.0, 80.0]
ISO_COLUMNS = ["reco_eiso_r30", "reco_eiso_r40"]
DEFAULT_BASELINE_32 = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
    "cluster_weta35_cogx",
    "cluster_wphi53_cogx",
    "cluster_w32",
    "cluster_w52",
    "cluster_w72",
    "e11_over_e22",
    "e11_over_e13",
    "e11_over_e15",
    "e11_over_e17",
    "e11_over_e31",
    "e11_over_e51",
    "e11_over_e71",
    "e22_over_e33",
    "e22_over_e35",
    "e22_over_e37",
    "e22_over_e53",
    "cluster_weta_over_wphi",
    "cluster_weta33_over_wphi33",
    "centrality",
]

FEATURE_FAMILY = {
    "cluster_Et": "kinematic/context",
    "vertexz": "kinematic/context",
    "cluster_Eta": "kinematic/context",
    "centrality": "kinematic/context",
    "cluster_weta_cogx": "shower widths",
    "cluster_wphi_cogx": "shower widths",
    "cluster_weta33_cogx": "shower widths",
    "cluster_wphi33_cogx": "shower widths",
    "cluster_weta35_cogx": "shower widths",
    "cluster_wphi53_cogx": "shower widths",
    "cluster_w32": "shower widths",
    "cluster_w52": "shower widths",
    "cluster_w72": "shower widths",
    "cluster_weta_over_wphi": "width ratios",
    "cluster_weta33_over_wphi33": "width ratios",
    "cluster_et1": "tower energy sharing",
    "cluster_et2": "tower energy sharing",
    "cluster_et3": "tower energy sharing",
    "cluster_et4": "tower energy sharing",
    "e11_over_e33": "E11 core ratios",
    "e11_over_e22": "E11 core ratios",
    "e11_over_e13": "E11 core ratios",
    "e11_over_e15": "E11 core ratios",
    "e11_over_e17": "E11 core ratios",
    "e11_over_e31": "E11 core ratios",
    "e11_over_e51": "E11 core ratios",
    "e11_over_e71": "E11 core ratios",
    "e32_over_e35": "containment ratios",
    "e22_over_e33": "containment ratios",
    "e22_over_e35": "containment ratios",
    "e22_over_e37": "containment ratios",
    "e22_over_e53": "containment ratios",
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

FAMILY_COLORS = {
    "BDT score": "#6B7280",
    "kinematic/context": "#4E79A7",
    "shower widths": "#0072B2",
    "width ratios": "#56B4E9",
    "tower energy sharing": "#009E73",
    "E11 core ratios": "#D55E00",
    "containment ratios": "#CC79A7",
    "isolation": "#7F3C8D",
}

DISPLAY_LABELS = {
    "cluster_Et": r"cluster $E_T$",
    "cluster_Eta": r"cluster $\eta$",
    "vertexz": r"$z_{vtx}$",
    "centrality": "centrality",
    "cluster_weta_cogx": r"$w_\eta$",
    "cluster_wphi_cogx": r"$w_\phi$",
    "cluster_weta33_cogx": r"$w_\eta^{3x3}$",
    "cluster_wphi33_cogx": r"$w_\phi^{3x3}$",
    "cluster_weta35_cogx": r"$w_\eta^{3x5}$",
    "cluster_wphi53_cogx": r"$w_\phi^{5x3}$",
    "cluster_w32": r"$w_{32}$",
    "cluster_w52": r"$w_{52}$",
    "cluster_w72": r"$w_{72}$",
    "cluster_weta_over_wphi": r"$w_\eta/w_\phi$",
    "cluster_weta33_over_wphi33": r"$w_\eta^{3x3}/w_\phi^{3x3}$",
    "cluster_et1": "cluster et1",
    "cluster_et2": "cluster et2",
    "cluster_et3": "cluster et3",
    "cluster_et4": "cluster et4",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "e11_over_e22": r"$E_{11}/E_{22}$",
    "e11_over_e13": r"$E_{11}/E_{13}$",
    "e11_over_e15": r"$E_{11}/E_{15}$",
    "e11_over_e17": r"$E_{11}/E_{17}$",
    "e11_over_e31": r"$E_{11}/E_{31}$",
    "e11_over_e51": r"$E_{11}/E_{51}$",
    "e11_over_e71": r"$E_{11}/E_{71}$",
    "e32_over_e35": r"$E_{32}/E_{35}$",
    "e22_over_e33": r"$E_{22}/E_{33}$",
    "e22_over_e35": r"$E_{22}/E_{35}$",
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
    "reco_eiso_r30": r"R=0.3 $E_T^{iso}$",
    "reco_eiso_r40": r"R=0.4 $E_T^{iso}$",
}
ISO_PLAIN_LABELS = {
    "reco_eiso_r30": r"R=0.3 $E_T^{iso}$",
    "reco_eiso_r40": r"R=0.4 $E_T^{iso}$",
}


@dataclass
class Accumulator:
    n: int = 0
    sx: float = 0.0
    sy: float = 0.0
    sx2: float = 0.0
    sy2: float = 0.0
    sxy: float = 0.0

    def update(self, x: np.ndarray, y: np.ndarray) -> None:
        finite = np.isfinite(x) & np.isfinite(y)
        if not np.any(finite):
            return
        xf = x[finite].astype("float64", copy=False)
        yf = y[finite].astype("float64", copy=False)
        self.n += int(len(xf))
        self.sx += float(xf.sum())
        self.sy += float(yf.sum())
        self.sx2 += float(np.dot(xf, xf))
        self.sy2 += float(np.dot(yf, yf))
        self.sxy += float(np.dot(xf, yf))

    def corr(self) -> float:
        if self.n < 3:
            return math.nan
        cov = self.sxy - self.sx * self.sy / self.n
        vx = self.sx2 - self.sx * self.sx / self.n
        vy = self.sy2 - self.sy * self.sy / self.n
        if vx <= 0.0 or vy <= 0.0:
            return math.nan
        return cov / math.sqrt(vx * vy)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--score-cache-manifest", type=Path, required=True)
    ap.add_argument(
        "--aux-score-cache-manifest",
        type=Path,
        default=None,
        help=(
            "Optional aligned score-cache manifest used only for score columns. "
            "This supports joining raw-isolation caches with a no-isolation score rescore."
        ),
    )
    ap.add_argument("--feature-dictionary", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--score-column", action="append", default=[], help="Score column to include, with or without score_ prefix.")
    ap.add_argument("--label", default="Raw isolation vs 32-feature BDT inputs")
    ap.add_argument("--spearman-max-rows", type=int, default=750000)
    ap.add_argument("--max-cache-files", type=int, default=0, help="Debug throttle; 0 means all files.")
    return ap.parse_args()


def read_manifest(path: Path) -> list[Path]:
    paths = [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]
    if not paths:
        raise SystemExit(f"Manifest is empty: {path}")
    return paths


def read_cache_columns(path: Path, columns: list[str]) -> tuple[dict[str, np.ndarray], list[str]]:
    with np.load(path, allow_pickle=True) as data:
        missing = [col for col in columns if col not in data.files]
        frame = {col: np.asarray(data[col]) for col in columns if col in data.files}
    return frame, missing


def aligned_score_pairs(primary_paths: list[Path], aux_manifest: Path | None, max_cache_files: int = 0) -> list[tuple[Path, Path | None]]:
    if aux_manifest is None:
        return [(path, None) for path in primary_paths]
    aux_paths = read_manifest(aux_manifest)
    if max_cache_files > 0:
        aux_paths = aux_paths[:max_cache_files]
    if len(aux_paths) != len(primary_paths):
        raise SystemExit(
            f"Aligned manifests have different lengths: primary={len(primary_paths)} aux={len(aux_paths)}"
        )
    return list(zip(primary_paths, aux_paths, strict=True))


def read_baseline_features(path: Path) -> list[str]:
    if not path.is_file():
        print(f"[isoCorr][WARN] feature dictionary missing, using built-in 32-feature baseline: {path}", flush=True)
        return list(DEFAULT_BASELINE_32)
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    features = [r["feature"] for r in rows if r.get("used_in_baseline_32", "").strip().lower() == "yes"]
    if len(features) != 32:
        raise SystemExit(f"Expected 32 baseline features from {path}, got {len(features)}")
    return features


def normalize_score_columns(columns: list[str]) -> list[str]:
    out = []
    for col in columns:
        col = col.strip()
        if not col:
            continue
        out.append(col if col.startswith("score_") else f"score_{col}")
    return out


def range_key(prefix: str, lo: float, hi: float) -> str:
    return f"{prefix}_{lo:g}_{hi:g}".replace(".", "p")


def range_label(lo: float, hi: float, suffix: str = "") -> str:
    return f"{lo:g}-{hi:g}{suffix}"


def build_scopes(et: np.ndarray, cent: np.ndarray) -> list[tuple[str, str, np.ndarray, str]]:
    scopes: list[tuple[str, str, np.ndarray, str]] = [
        ("inclusive", "inclusive", np.ones(len(et), dtype=bool), "inclusive"),
    ]
    for lo, hi in zip(CENT_EDGES[:-1], CENT_EDGES[1:]):
        scopes.append((range_key("cent", lo, hi), range_label(lo, hi, "%"), (cent >= lo) & (cent < hi), "centrality"))
    for lo, hi in zip(PT_EDGES[:-1], PT_EDGES[1:]):
        scopes.append((range_key("pt", lo, hi), range_label(lo, hi, " GeV"), (et >= lo) & (et < hi), "et"))
    for clo, chi in zip(CENT_EDGES[:-1], CENT_EDGES[1:]):
        cmask = (cent >= clo) & (cent < chi)
        for plo, phi in zip(PT_EDGES[:-1], PT_EDGES[1:]):
            key = f"{range_key('cent', clo, chi)}__{range_key('pt', plo, phi)}"
            label = f"{range_label(clo, chi, '%')}, {range_label(plo, phi, ' GeV')}"
            scopes.append((key, label, cmask & (et >= plo) & (et < phi), "centrality_et"))
    return scopes


def class_masks(y: np.ndarray) -> list[tuple[str, str, np.ndarray]]:
    return [
        ("all", "all candidates", np.isin(y, [0, 1])),
        ("signal", "truth photons", y == 1),
        ("background", "jet background", y == 0),
    ]


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype="float64")
    sorted_vals = values[order]
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and sorted_vals[j] == sorted_vals[i]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0
        i = j
    return ranks


def pearson_from_arrays(x: np.ndarray, y: np.ndarray) -> float:
    finite = np.isfinite(x) & np.isfinite(y)
    if int(finite.sum()) < 3:
        return math.nan
    xf = x[finite].astype("float64", copy=False)
    yf = y[finite].astype("float64", copy=False)
    if np.std(xf) <= 0.0 or np.std(yf) <= 0.0:
        return math.nan
    return float(np.corrcoef(xf, yf)[0, 1])


def spearman_from_arrays(x: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    finite = np.isfinite(x) & np.isfinite(y)
    n = int(finite.sum())
    if n < 3:
        return math.nan, n
    rx = rankdata(x[finite].astype("float64", copy=False))
    ry = rankdata(y[finite].astype("float64", copy=False))
    return pearson_from_arrays(rx, ry), n


def correlation_row_key(scope_key: str, class_key: str, iso_col: str, target_col: str) -> tuple[str, str, str, str]:
    return scope_key, class_key, iso_col, target_col


def update_accumulators(
    acc: dict[tuple[str, str, str, str], Accumulator],
    data: dict[str, np.ndarray],
    targets: list[str],
) -> None:
    et = data["cluster_Et"].astype("float64", copy=False)
    cent = data["centrality"].astype("float64", copy=False)
    y = data["is_signal"].astype("int32", copy=False)
    base = np.isfinite(et) & np.isfinite(cent) & (et >= 15.0) & (et < 35.0) & (cent >= 0.0) & (cent < 80.0)
    for scope_key, _scope_label, scope_mask, _scope_type in build_scopes(et, cent):
        smask = base & scope_mask
        if not np.any(smask):
            continue
        for class_key, _class_label, cmask in class_masks(y):
            mask = smask & cmask
            if not np.any(mask):
                continue
            for iso_col in ISO_COLUMNS:
                x = data[iso_col][mask]
                for target_col in targets:
                    acc.setdefault(correlation_row_key(scope_key, class_key, iso_col, target_col), Accumulator()).update(
                        x, data[target_col][mask]
                    )
            acc.setdefault(correlation_row_key(scope_key, class_key, "reco_eiso_r30", "reco_eiso_r40"), Accumulator()).update(
                data["reco_eiso_r30"][mask], data["reco_eiso_r40"][mask]
            )


def validate_aux_alignment(primary: dict[str, np.ndarray], aux: dict[str, np.ndarray], path: Path) -> None:
    for col in ["is_signal", "cluster_Et", "centrality"]:
        if col not in aux:
            continue
        if len(primary[col]) != len(aux[col]):
            raise SystemExit(
                f"Aux cache row-count mismatch for {path}: {col} primary={len(primary[col])} aux={len(aux[col])}"
            )
        if col == "is_signal":
            if not np.array_equal(primary[col], aux[col]):
                raise SystemExit(f"Aux cache alignment mismatch for {path}: {col} differs")
        else:
            same = np.isclose(
                primary[col].astype("float64", copy=False),
                aux[col].astype("float64", copy=False),
                rtol=0.0,
                atol=1e-8,
            )
            if not bool(np.all(same)):
                raise SystemExit(f"Aux cache alignment mismatch for {path}: {col} differs")


def maybe_add_sample(
    samples: list[dict[str, np.ndarray]],
    data: dict[str, np.ndarray],
    columns: list[str],
    rng: np.random.Generator,
    max_rows: int,
) -> None:
    if max_rows <= 0:
        return
    n = len(data["is_signal"])
    if n == 0:
        return
    # Keep roughly a few percent initially, then downsample globally. This avoids
    # needing a full pass to know total rows while still making the Spearman
    # cross-check deterministic and memory bounded.
    take = min(n, max(1000, max_rows // 60))
    idx = rng.choice(n, size=take, replace=False) if take < n else np.arange(n)
    samples.append({col: data[col][idx].copy() for col in columns})
    total = sum(len(s["is_signal"]) for s in samples)
    if total <= max_rows * 2:
        return
    merged = {col: np.concatenate([s[col] for s in samples]) for col in columns}
    keep = rng.choice(len(merged["is_signal"]), size=max_rows, replace=False)
    samples[:] = [{col: merged[col][keep] for col in columns}]


def finalize_sample(samples: list[dict[str, np.ndarray]], columns: list[str], max_rows: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
    if not samples:
        return {col: np.array([], dtype="float32") for col in columns}
    merged = {col: np.concatenate([s[col] for s in samples]) for col in columns}
    n = len(merged["is_signal"])
    if max_rows > 0 and n > max_rows:
        keep = rng.choice(n, size=max_rows, replace=False)
        merged = {col: val[keep] for col, val in merged.items()}
    return merged


def compute_spearman_lookup(sample: dict[str, np.ndarray], targets: list[str]) -> dict[tuple[str, str, str, str], tuple[float, int]]:
    lookup: dict[tuple[str, str, str, str], tuple[float, int]] = {}
    if len(sample["is_signal"]) == 0:
        return lookup
    et = sample["cluster_Et"].astype("float64", copy=False)
    cent = sample["centrality"].astype("float64", copy=False)
    y = sample["is_signal"].astype("int32", copy=False)
    base = np.isfinite(et) & np.isfinite(cent) & (et >= 15.0) & (et < 35.0) & (cent >= 0.0) & (cent < 80.0)
    for scope_key, _scope_label, scope_mask, _scope_type in build_scopes(et, cent):
        smask = base & scope_mask
        if not np.any(smask):
            continue
        for class_key, _class_label, cmask in class_masks(y):
            mask = smask & cmask
            if not np.any(mask):
                continue
            for iso_col in ISO_COLUMNS:
                x = sample[iso_col][mask]
                for target_col in targets:
                    lookup[correlation_row_key(scope_key, class_key, iso_col, target_col)] = spearman_from_arrays(
                        x, sample[target_col][mask]
                    )
            lookup[correlation_row_key(scope_key, class_key, "reco_eiso_r30", "reco_eiso_r40")] = spearman_from_arrays(
                sample["reco_eiso_r30"][mask], sample["reco_eiso_r40"][mask]
            )
    return lookup


def scope_metadata() -> dict[str, tuple[str, str]]:
    dummy = np.zeros(1)
    return {key: (label, typ) for key, label, _mask, typ in build_scopes(dummy, dummy)}


def class_label_map() -> dict[str, str]:
    dummy = np.zeros(1, dtype="int32")
    return {key: label for key, label, _mask in class_masks(dummy)}


def target_family(target: str, score_columns: list[str]) -> str:
    if target in score_columns:
        return "BDT score"
    if target in ISO_COLUMNS:
        return "isolation"
    return FEATURE_FAMILY.get(target, "other")


def target_label(target: str, score_columns: list[str]) -> str:
    if target in score_columns:
        product = target.removeprefix("score_")
        if product == "globalEtCent1535_bdt_noIso_ptCent7":
            return "baseline BDT score"
        if "eiso" in product:
            return "isolation-input BDT score"
        return product.replace("globalEtCent1535_bdt_", "").replace("_", " ")
    return DISPLAY_LABELS.get(target, target)


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "scope",
        "scope_label",
        "scope_type",
        "class",
        "class_label",
        "iso_variable",
        "target_variable",
        "target_label",
        "feature_family",
        "pearson",
        "spearman",
        "n_entries",
        "spearman_n_entries",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def as_float(text: str) -> float:
    try:
        return float(text)
    except Exception:
        return math.nan


def plot_full_heatmap(rows: list[dict[str, str]], targets: list[str], score_columns: list[str], out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = ["all", "signal", "background"]
    ylabels = []
    mat = []
    for iso_col in ISO_COLUMNS:
        for cls in classes:
            ylabels.append(f"{ISO_PLAIN_LABELS[iso_col]}\n{cls}")
            vals = []
            for target in targets:
                match = next(
                    (
                        r
                        for r in rows
                        if r["scope"] == "inclusive"
                        and r["class"] == cls
                        and r["iso_variable"] == iso_col
                        and r["target_variable"] == target
                    ),
                    None,
                )
                vals.append(as_float(match["pearson"]) if match else math.nan)
            mat.append(vals)
    data = np.asarray(mat, dtype="float64")

    fig, ax = plt.subplots(figsize=(18.5, 7.2))
    im = ax.imshow(data, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_yticks(np.arange(len(ylabels)), ylabels, fontsize=11)
    ax.set_xticks(np.arange(len(targets)), [target_label(t, score_columns) for t in targets], rotation=55, ha="right", fontsize=10)
    fig.text(0.06, 0.965, "Raw cone isolation correlations with 32-feature BDT inputs and score", fontsize=22, fontweight="bold", ha="left")
    fig.text(0.06, 0.925, "Pearson correlation, 15 < cluster $E_T$ < 35 GeV; rows split by truth class", fontsize=13, color="#4B5563", ha="left")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if not math.isfinite(val):
                continue
            color = "white" if abs(val) > 0.55 else "#111827"
            ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=8.5, color=color)
    cbar = fig.colorbar(im, ax=ax, pad=0.012)
    cbar.set_label("Pearson correlation", fontsize=12)
    ax.tick_params(axis="both", length=0)
    fig.subplots_adjust(left=0.07, right=0.93, top=0.86, bottom=0.32)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_family_summary(rows: list[dict[str, str]], out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = ["all", "signal", "background"]
    families = [f for f in FAMILY_ORDER if any(r["feature_family"] == f for r in rows)]
    xlabels = []
    mat = []
    for iso_col in ISO_COLUMNS:
        for cls in classes:
            xlabels.append(f"{ISO_PLAIN_LABELS[iso_col]}\n{cls}")
    for fam in families:
        vals = []
        for iso_col in ISO_COLUMNS:
            for cls in classes:
                fam_vals = [
                    as_float(r["pearson"])
                    for r in rows
                    if r["scope"] == "inclusive"
                    and r["class"] == cls
                    and r["iso_variable"] == iso_col
                    and r["feature_family"] == fam
                    and r["target_variable"] != "reco_eiso_r40"
                ]
                vals.append(float(np.nanmax(np.abs(fam_vals))) if fam_vals else math.nan)
        mat.append(vals)
    data = np.asarray(mat, dtype="float64")

    fig, ax = plt.subplots(figsize=(13.5, 6.8))
    im = ax.imshow(data, cmap="viridis", vmin=0.0, vmax=max(0.75, np.nanmax(data)), aspect="auto")
    ax.set_xticks(np.arange(len(xlabels)), xlabels, fontsize=10)
    ax.set_yticks(np.arange(len(families)), families, fontsize=12)
    fig.text(0.12, 0.965, "Where raw isolation is entangled with the BDT input space", fontsize=20, fontweight="bold", ha="left")
    fig.text(0.12, 0.915, "Each cell shows the strongest absolute Pearson correlation in that feature family", fontsize=13, color="#4B5563", ha="left")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if not math.isfinite(val):
                continue
            color = "white" if val > 0.45 else "#111827"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=12, color=color, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, pad=0.014)
    cbar.set_label("max |Pearson correlation|", fontsize=12)
    ax.tick_params(axis="both", length=0)
    fig.subplots_adjust(left=0.20, right=0.90, top=0.84, bottom=0.18)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_r30_r40_cent_et(rows: list[dict[str, str]], out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    classes = ["all", "signal", "background"]
    cent_labels = [range_label(lo, hi, "%") for lo, hi in zip(CENT_EDGES[:-1], CENT_EDGES[1:])]
    pt_labels = [range_label(lo, hi, "") for lo, hi in zip(PT_EDGES[:-1], PT_EDGES[1:])]
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.2), sharey=True)
    for ax, cls in zip(axes, classes):
        mat = np.full((len(cent_labels), len(pt_labels)), np.nan)
        for i, (clo, chi) in enumerate(zip(CENT_EDGES[:-1], CENT_EDGES[1:])):
            for j, (plo, phi) in enumerate(zip(PT_EDGES[:-1], PT_EDGES[1:])):
                scope = f"{range_key('cent', clo, chi)}__{range_key('pt', plo, phi)}"
                match = next(
                    (
                        r
                        for r in rows
                        if r["scope"] == scope
                        and r["class"] == cls
                        and r["iso_variable"] == "reco_eiso_r30"
                        and r["target_variable"] == "reco_eiso_r40"
                    ),
                    None,
                )
                if match:
                    mat[i, j] = as_float(match["pearson"])
        im = ax.imshow(mat, cmap="magma", vmin=0.0, vmax=1.0, aspect="auto")
        ax.set_title(cls, fontsize=15, fontweight="bold")
        ax.set_xticks(np.arange(len(pt_labels)), pt_labels, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(np.arange(len(cent_labels)), cent_labels, fontsize=10)
        ax.set_xlabel("cluster $E_T$ bin [GeV]", fontsize=11)
        if ax is axes[0]:
            ax.set_ylabel("centrality", fontsize=11)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                val = mat[i, j]
                if math.isfinite(val):
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", color="white" if val > 0.55 else "#111827", fontsize=10)
    fig.suptitle(r"Correlation between raw R=0.3 and R=0.4 $E_T^{iso}$", fontsize=20, fontweight="bold", y=0.96)
    fig.subplots_adjust(left=0.08, right=0.89, top=0.76, bottom=0.22, wspace=0.10)
    cax = fig.add_axes([0.915, 0.22, 0.018, 0.54])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Pearson correlation", fontsize=12)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    features = read_baseline_features(args.feature_dictionary)
    score_columns = normalize_score_columns(args.score_column)
    paths = read_manifest(args.score_cache_manifest)
    if args.max_cache_files > 0:
        paths = paths[: args.max_cache_files]
    pairs = aligned_score_pairs(paths, args.aux_score_cache_manifest, args.max_cache_files)

    primary_required = ["is_signal", "cluster_Et", "centrality"] + ISO_COLUMNS + features
    aux_alignment_columns = ["is_signal", "cluster_Et", "centrality"]
    required = primary_required + score_columns
    targets = features + score_columns
    missing_columns: dict[str, list[str]] = {}
    acc: dict[tuple[str, str, str, str], Accumulator] = {}
    sample_columns = sorted(set(required))
    samples: list[dict[str, np.ndarray]] = []
    rng = np.random.default_rng(1535)
    files_read = 0
    entries_seen = 0

    for idx, (path, aux_path) in enumerate(pairs, 1):
        frame, primary_missing = read_cache_columns(path, primary_required + score_columns)
        missing = list(primary_missing)
        if score_columns and aux_path is not None:
            missing_scores = [col for col in score_columns if col not in frame]
            if missing_scores:
                aux_frame, aux_missing = read_cache_columns(aux_path, aux_alignment_columns + missing_scores)
                if aux_missing:
                    missing_columns[f"{path} + {aux_path}"] = sorted(set(missing + [f"aux:{col}" for col in aux_missing]))
                    continue
                validate_aux_alignment(frame, aux_frame, aux_path)
                for col in missing_scores:
                    frame[col] = aux_frame[col]
                missing = [col for col in missing if col not in score_columns]
        if missing:
            missing_columns[str(path)] = missing
            continue
        frame = {col: frame[col] for col in required}
        files_read += 1
        entries_seen += int(len(frame["is_signal"]))
        update_accumulators(acc, frame, targets)
        maybe_add_sample(samples, frame, sample_columns, rng, args.spearman_max_rows)
        if idx == 1 or idx % 25 == 0 or idx == len(paths):
            print(f"[isoCorr] processed {idx}/{len(paths)} files; usable={files_read}; entries_seen={entries_seen}", flush=True)

    if files_read == 0:
        preview = next(iter(missing_columns.items())) if missing_columns else ("<none>", [])
        raise SystemExit(f"No usable score caches. Example missing columns: {preview}")

    sample = finalize_sample(samples, sample_columns, args.spearman_max_rows, rng)
    spearman = compute_spearman_lookup(sample, targets)
    scopes = scope_metadata()
    labels = class_label_map()

    rows: list[dict[str, object]] = []
    for key, a in sorted(acc.items()):
        scope, cls, iso_col, target_col = key
        scope_label, scope_type = scopes.get(scope, (scope, "unknown"))
        s_val, s_n = spearman.get(key, (math.nan, 0))
        rows.append(
            {
                "scope": scope,
                "scope_label": scope_label,
                "scope_type": scope_type,
                "class": cls,
                "class_label": labels.get(cls, cls),
                "iso_variable": iso_col,
                "target_variable": target_col,
                "target_label": target_label(target_col, score_columns),
                "feature_family": target_family(target_col, score_columns),
                "pearson": f"{a.corr():.8g}" if math.isfinite(a.corr()) else "nan",
                "spearman": f"{s_val:.8g}" if math.isfinite(s_val) else "nan",
                "n_entries": a.n,
                "spearman_n_entries": s_n,
            }
        )

    csv_path = args.outdir / "isolation_feature_score_correlations.csv"
    write_csv(rows, csv_path)
    rows_as_str = load_csv_rows(csv_path)
    full_png = args.outdir / "isolation_feature_score_correlation_full_heatmap_inclusive.png"
    family_png = args.outdir / "isolation_feature_score_correlation_family_summary.png"
    cone_png = args.outdir / "r30_r40_isolation_correlation_cent_et.png"
    plot_full_heatmap(rows_as_str, targets, score_columns, full_png)
    plot_family_summary(rows_as_str, family_png)
    plot_r30_r40_cent_et(rows_as_str, cone_png)

    summary = {
        "label": args.label,
        "schema": "AUAU_ISOLATION_FEATURE_CORRELATION_V1",
        "score_cache_manifest": str(args.score_cache_manifest),
        "aux_score_cache_manifest": str(args.aux_score_cache_manifest) if args.aux_score_cache_manifest else None,
        "feature_dictionary": str(args.feature_dictionary),
        "files_requested": len(paths),
        "files_read": files_read,
        "entries_seen": entries_seen,
        "features": features,
        "score_columns": score_columns,
        "missing_column_files": len(missing_columns),
        "outputs": {
            "csv": str(csv_path),
            "full_heatmap": str(full_png),
            "family_summary": str(family_png),
            "r30_r40_cent_et": str(cone_png),
        },
        "notes": [
            "Pearson correlations use all rows available in the score-cache manifest after the 15 < cluster_Et < 35 GeV and 0 <= centrality < 80 selections.",
            f"Spearman correlations are computed on a deterministic global sample capped at {args.spearman_max_rows} rows.",
            "The main ABCD-safe BDT-score diagnostic should use a no-isolation BDT score column, not an isolation-input score column.",
        ],
    }
    (args.outdir / "isolation_feature_score_correlations_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
