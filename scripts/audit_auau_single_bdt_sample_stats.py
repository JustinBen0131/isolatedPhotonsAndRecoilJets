#!/usr/bin/env python3
"""Audit sample statistics for the single global AuAu BDT baseline.

The intended slide baseline is `centInput_pt1535`: one BDT trained over
15 < E_T < 35 GeV with E_T, centrality, base shower shapes, and 3x3 widths as
inputs.  This script counts the signal/background rows that enter that simple
configuration and reproduces the BDT train/test split at the count level.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


CENT_INPUT_PT1535_FEATURES = [
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
    "centrality",
]

DEFAULT_ET_BINS = "15,17,19,21,23,25,27,30,35"
DEFAULT_CENT_BINS = "0,10,20,30,40,50,60,80"
SAMPLE_ORDER = [
    "embeddedPhoton12",
    "embeddedPhoton20",
    "embeddedPhoton12plus20",
    "embeddedJet12",
    "embeddedJet20",
    "embeddedJet30",
    "embeddedJet12plus20",
    "embeddedJet12plus20plus30",
]


@dataclass
class SampleArrays:
    sample: str
    source_role: str
    source_id: np.ndarray
    run: np.ndarray
    evt: np.ndarray
    et: np.ndarray
    cent: np.ndarray
    label: np.ndarray
    finite_features: np.ndarray
    raw_entries: int


def parse_edges(text: str) -> list[float]:
    values = [float(tok) for tok in re.split(r"[,;:\s]+", text.strip()) if tok]
    if len(values) < 2 or any(values[i] >= values[i + 1] for i in range(len(values) - 1)):
        raise SystemExit(f"Invalid bin edges: {text!r}")
    return values


def parse_range(text: str) -> tuple[float, float]:
    values = parse_edges(text)
    if len(values) != 2:
        raise SystemExit(f"Range must contain exactly two values: {text!r}")
    return values[0], values[1]


def sample_from_path(path: str) -> tuple[str, str]:
    lower = path.lower()
    if "embeddedphoton12" in lower:
        return "embeddedPhoton12", "signal"
    if "embeddedphoton20" in lower:
        return "embeddedPhoton20", "signal"
    if "embeddedjet12" in lower:
        return "embeddedJet12", "background"
    if "embeddedjet20" in lower:
        return "embeddedJet20", "background"
    if "embeddedjet30" in lower:
        return "embeddedJet30", "background"
    return "unknown", "unknown"


def collect_root_files(source: Path, explicit_files: list[Path]) -> dict[str, list[str]]:
    files = [str(p) for p in explicit_files]
    if source:
        root = source / "extraction"
        if root.is_dir():
            files.extend(str(p) for p in root.rglob("*.root"))
        elif source.is_dir():
            files.extend(str(p) for p in source.rglob("*.root"))
    grouped: dict[str, list[str]] = {}
    for path in sorted(set(files)):
        sample, _ = sample_from_path(path)
        if sample == "unknown":
            continue
        grouped.setdefault(sample, []).append(path)
    return grouped


def vector_of_strings(paths: list[str]):
    import ROOT  # type: ignore

    vec = ROOT.std.vector("string")()
    for path in paths:
        vec.push_back(path)
    return vec


def load_root_sample(sample: str, paths: list[str], tree: str, max_files: int = 0) -> SampleArrays:
    import ROOT  # type: ignore

    if max_files > 0:
        paths = paths[:max_files]
    if not paths:
        raise SystemExit(f"No ROOT files for {sample}")
    role = sample_from_path(paths[0])[1]
    rdf = ROOT.RDataFrame(tree, vector_of_strings(paths))
    raw_entries = int(rdf.Count().GetValue())
    columns = ["run", "evt", "is_signal", "cluster_Et", "centrality"] + [
        name for name in CENT_INPUT_PT1535_FEATURES if name not in {"cluster_Et", "centrality"}
    ]
    arrays = rdf.AsNumpy(columns)
    source_code = SAMPLE_ORDER.index(sample) + 1 if sample in SAMPLE_ORDER else 0
    source_id = np.full(raw_entries, source_code, dtype=np.int16)
    run = np.asarray(arrays["run"], dtype=np.int64)
    evt = np.asarray(arrays["evt"], dtype=np.int64)
    label = np.asarray(arrays["is_signal"], dtype=np.int8)
    et = np.asarray(arrays["cluster_Et"], dtype=np.float64)
    cent = np.asarray(arrays["centrality"], dtype=np.float64)
    finite = np.isfinite(et) & np.isfinite(cent)
    for feature in CENT_INPUT_PT1535_FEATURES:
        finite &= np.isfinite(np.asarray(arrays[feature], dtype=np.float64))
    return SampleArrays(sample, role, source_id, run, evt, et, cent, label, finite, raw_entries)


def localize_cache_path(raw: str, list_path: Path) -> Path:
    path = Path(raw)
    if path.exists():
        return path
    candidate = list_path.parent / "score_caches" / path.name
    if candidate.exists():
        return candidate
    candidate = list_path.parent / path.name
    if candidate.exists():
        return candidate
    return path


def load_cache_samples(cache_list: Path) -> list[SampleArrays]:
    # Cache files do not preserve Photon12/20 vs Jet12/20 provenance, so this is
    # a fallback for local QA. Exact slide counts should use ROOT mode.
    paths = [
        localize_cache_path(line.strip(), cache_list)
        for line in cache_list.read_text().splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    et_parts: list[np.ndarray] = []
    cent_parts: list[np.ndarray] = []
    label_parts: list[np.ndarray] = []
    finite_parts: list[np.ndarray] = []
    raw_entries = 0
    source_parts: list[np.ndarray] = []
    run_parts: list[np.ndarray] = []
    evt_parts: list[np.ndarray] = []
    for path in paths:
        z = np.load(path, allow_pickle=True)
        label = np.asarray(z["is_signal"], dtype=np.int8)
        et = np.asarray(z["cluster_Et"], dtype=np.float64)
        cent = np.asarray(z["centrality"], dtype=np.float64)
        source = np.full(len(label), 0, dtype=np.int16)
        run = np.arange(raw_entries, raw_entries + len(label), dtype=np.int64)
        evt = np.zeros(len(label), dtype=np.int64)
        finite = np.isfinite(et) & np.isfinite(cent)
        for feature in CENT_INPUT_PT1535_FEATURES:
            if feature in z.files:
                finite &= np.isfinite(np.asarray(z[feature], dtype=np.float64))
        if "counts_json" in z.files:
            try:
                payload = json.loads(str(z["counts_json"].item()))
                raw_entries += int(payload.get("total_entries", 0))
            except Exception:
                raw_entries += int(len(label))
        else:
            raw_entries += int(len(label))
        et_parts.append(et)
        cent_parts.append(cent)
        label_parts.append(label)
        finite_parts.append(finite)
        source_parts.append(source)
        run_parts.append(run)
        evt_parts.append(evt)
    return [
        SampleArrays(
            "cappedValidationCache",
            "mixed",
            np.concatenate(source_parts),
            np.concatenate(run_parts),
            np.concatenate(evt_parts),
            np.concatenate(et_parts),
            np.concatenate(cent_parts),
            np.concatenate(label_parts),
            np.concatenate(finite_parts),
            raw_entries,
        )
    ]


def combine_samples(samples: list[SampleArrays]) -> list[SampleArrays]:
    out = list(samples)
    by_name = {s.sample: s for s in samples}
    for combo, parts, role in [
        ("embeddedPhoton12plus20", ["embeddedPhoton12", "embeddedPhoton20"], "signal"),
        ("embeddedJet12plus20", ["embeddedJet12", "embeddedJet20"], "background"),
        ("embeddedJet12plus20plus30", ["embeddedJet12", "embeddedJet20", "embeddedJet30"], "background"),
    ]:
        have = [by_name[name] for name in parts if name in by_name]
        if len(have) != len(parts):
            continue
        out.append(
            SampleArrays(
                combo,
                role,
                np.concatenate([s.source_id for s in have]),
                np.concatenate([s.run for s in have]),
                np.concatenate([s.evt for s in have]),
                np.concatenate([s.et for s in have]),
                np.concatenate([s.cent for s in have]),
                np.concatenate([s.label for s in have]),
                np.concatenate([s.finite_features for s in have]),
                sum(s.raw_entries for s in have),
            )
        )
    return out


def split_mask(labels: np.ndarray, test_size: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    train = np.zeros(len(labels), dtype=bool)
    test = np.zeros(len(labels), dtype=bool)
    for cls in (0, 1):
        idx = np.flatnonzero(labels == cls)
        rng.shuffle(idx)
        n_test = int(math.ceil(len(idx) * test_size))
        test[idx[:n_test]] = True
        train[idx[n_test:]] = True
    return train, test


def count_labels(labels: np.ndarray, mask: np.ndarray) -> tuple[int, int]:
    local = labels[mask]
    return int(np.sum(local == 1)), int(np.sum(local == 0))


def count_unique_events(s: SampleArrays, mask: np.ndarray) -> int:
    if len(s.run) == 0:
        return 0
    dtype = np.dtype([("source", s.source_id.dtype), ("run", s.run.dtype), ("evt", s.evt.dtype)])
    payload = np.empty(int(np.sum(mask)), dtype=dtype)
    payload["source"] = s.source_id[mask]
    payload["run"] = s.run[mask]
    payload["evt"] = s.evt[mask]
    return int(len(np.unique(payload)))


def row_for_sample(s: SampleArrays, pt_range: tuple[float, float], test_size: float, seed: int) -> dict[str, object]:
    pt = (s.et >= pt_range[0]) & (s.et < pt_range[1])
    selected = pt & s.finite_features & np.isin(s.label, [0, 1])
    valid_label = np.isin(s.label, [0, 1])
    sig_raw, bkg_raw = count_labels(s.label, np.isin(s.label, [0, 1]))
    sig_sel, bkg_sel = count_labels(s.label, selected)
    train, test = split_mask(s.label[selected], test_size, seed)
    sig_train, bkg_train = count_labels(s.label[selected], train)
    sig_test, bkg_test = count_labels(s.label[selected], test)
    return {
        "sample": s.sample,
        "role": s.source_role,
        "raw_entries": s.raw_entries,
        "raw_unique_events": count_unique_events(s, valid_label),
        "raw_signal_labels": sig_raw,
        "raw_background_labels": bkg_raw,
        "selected_15_35_entries": int(np.sum(selected)),
        "selected_15_35_unique_events": count_unique_events(s, selected),
        "selected_signal_labels": sig_sel,
        "selected_background_labels": bkg_sel,
        "train_entries": int(np.sum(train)),
        "train_signal": sig_train,
        "train_background": bkg_train,
        "test_entries": int(np.sum(test)),
        "test_signal": sig_test,
        "test_background": bkg_test,
        "finite_feature_fraction": float(np.mean(s.finite_features)) if len(s.finite_features) else math.nan,
        "selected_candidates_per_event": (
            float(np.sum(selected)) / count_unique_events(s, selected) if count_unique_events(s, selected) else math.nan
        ),
        "selected_signal_fraction": (sig_sel / (sig_sel + bkg_sel)) if (sig_sel + bkg_sel) else math.nan,
        "selected_background_fraction": (bkg_sel / (sig_sel + bkg_sel)) if (sig_sel + bkg_sel) else math.nan,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def heatmap_rows(samples: list[SampleArrays], pt_edges: list[float], cent_edges: list[float], pt_range: tuple[float, float]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for s in samples:
        selected = (
            (s.et >= pt_range[0])
            & (s.et < pt_range[1])
            & s.finite_features
            & np.isin(s.label, [0, 1])
        )
        for ip, (plo, phi) in enumerate(zip(pt_edges[:-1], pt_edges[1:])):
            for ic, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
                mask = selected & (s.et >= plo) & (s.et < phi) & (s.cent >= clo) & (s.cent < chi)
                sig, bkg = count_labels(s.label, mask)
                rows.append(
                    {
                        "sample": s.sample,
                        "role": s.source_role,
                        "et_bin": f"{plo:g}-{phi:g}",
                        "cent_bin": f"{clo:g}-{chi:g}",
                        "signal": sig,
                        "background": bkg,
                        "entries": sig + bkg,
                        "background_to_signal": (bkg / sig) if sig else math.nan,
                    }
                )
    return rows


def plot_heatmap(rows: list[dict[str, object]], sample: str, value: str, pt_edges: list[float], cent_edges: list[float], outpath: Path, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    et_labels = [f"{lo:g}-{hi:g}" for lo, hi in zip(pt_edges[:-1], pt_edges[1:])]
    cent_labels = [f"{lo:g}-{hi:g}" for lo, hi in zip(cent_edges[:-1], cent_edges[1:])]
    matrix = np.full((len(cent_labels), len(et_labels)), np.nan, dtype=float)
    for row in rows:
        if row["sample"] != sample:
            continue
        if row["et_bin"] in et_labels and row["cent_bin"] in cent_labels:
            matrix[cent_labels.index(str(row["cent_bin"])), et_labels.index(str(row["et_bin"]))] = float(row[value])
    fig, ax = plt.subplots(figsize=(9.4, 5.4), dpi=180)
    cmap = "Blues" if value == "signal" else "Reds"
    im = ax.imshow(matrix, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(et_labels)), labels=et_labels)
    ax.set_yticks(np.arange(len(cent_labels)), labels=cent_labels)
    ax.set_xlabel(r"Candidate $E_T$ bin [GeV]", fontsize=12)
    ax.set_ylabel("Centrality bin [%]", fontsize=12)
    ax.set_title(title, fontsize=15, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            if math.isfinite(val):
                ax.text(j, i, f"{val:,.0f}", ha="center", va="center", fontsize=8.4, color="black")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Candidate count", fontsize=11)
    fig.text(0.02, 0.965, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.118, 0.965, "Internal", ha="left", va="top", fontsize=13)
    fig.text(0.02, 0.92, "Au+Au embedded photon-ID training", ha="left", va="top", fontsize=10.5)
    fig.tight_layout(rect=(0.02, 0.02, 0.98, 0.91))
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath)
    plt.close(fig)


def plot_truth_composition(sample_rows: list[dict[str, object]], outpath: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [
        r
        for r in sample_rows
        if r["sample"] in {"embeddedPhoton12", "embeddedPhoton20", "embeddedJet12", "embeddedJet20", "embeddedJet30"}
    ]
    rows.sort(key=lambda r: SAMPLE_ORDER.index(str(r["sample"])))
    labels = [str(r["sample"]).replace("embedded", "") for r in rows]
    signal = np.array([float(r["selected_signal_labels"]) for r in rows])
    background = np.array([float(r["selected_background_labels"]) for r in rows])
    x = np.arange(len(rows))

    fig, ax = plt.subplots(figsize=(9.2, 5.6), dpi=180)
    ax.bar(x, signal, color="#2f7fbf", edgecolor="black", linewidth=0.7, label="Truth signal label")
    ax.bar(x, background, bottom=signal, color="#c43b4d", edgecolor="black", linewidth=0.7, label="Truth background label")
    ax.set_xticks(x, labels)
    ax.set_ylabel(r"Selected candidates, $15<E_T<35$ GeV")
    ax.set_title("Truth-label composition of the training samples", fontsize=15, fontweight="bold")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    for i, (s, b) in enumerate(zip(signal, background)):
        total = s + b
        if total <= 0:
            continue
        dominant = s if s >= b else b
        frac = dominant / total
        ax.text(i, total * 1.015, f"{frac:.1%} {'signal' if s >= b else 'background'}", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath)
    plt.close(fig)


def event_candidate_rows(sample_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for row in sample_rows:
        rows.append(
            {
                "sample": row["sample"],
                "role": row["role"],
                "raw_candidate_rows": row["raw_entries"],
                "raw_unique_events_with_candidates": row["raw_unique_events"],
                "selected_15_35_candidates": row["selected_15_35_entries"],
                "selected_15_35_unique_events": row["selected_15_35_unique_events"],
                "selected_candidates_per_event": row["selected_candidates_per_event"],
                "selected_signal_labels": row["selected_signal_labels"],
                "selected_background_labels": row["selected_background_labels"],
                "selected_signal_fraction": row["selected_signal_fraction"],
                "selected_background_fraction": row["selected_background_fraction"],
            }
        )
    return rows


def config_rows(test_size: float, seed: int) -> list[dict[str, object]]:
    return [
        {
            "field": "validated_product",
            "value": "centInput_pt1535",
        },
        {
            "field": "model_class",
            "value": "single global XGBoost/BDT classifier",
        },
        {
            "field": "training_range",
            "value": "15 < cluster_Et < 35 GeV",
        },
        {
            "field": "number_of_bdts",
            "value": 1,
        },
        {
            "field": "routing",
            "value": "none; E_T and centrality are input features",
        },
        {
            "field": "features",
            "value": ", ".join(CENT_INPUT_PT1535_FEATURES),
        },
        {
            "field": "isolation_inputs",
            "value": "none",
        },
        {
            "field": "split",
            "value": f"{100*(1-test_size):.0f}% train / {100*test_size:.0f}% test, stratified by is_signal",
        },
        {
            "field": "random_seed",
            "value": seed,
        },
    ]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=Path("/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_20260508_233049"))
    ap.add_argument("--root-files", nargs="*", type=Path, default=[])
    ap.add_argument("--score-cache-list", type=Path, default=None)
    ap.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    ap.add_argument("--outdir", type=Path, default=Path("dataOutput/auauBDTSampleAudit/single_bdt_centInput_pt1535"))
    ap.add_argument("--pt-range", default="15:35")
    ap.add_argument("--et-bins", default=DEFAULT_ET_BINS)
    ap.add_argument("--cent-bins", default=DEFAULT_CENT_BINS)
    ap.add_argument("--test-size", type=float, default=0.10)
    ap.add_argument("--random-seed", type=int, default=13)
    ap.add_argument("--max-files-per-sample", type=int, default=0)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    pt_range = parse_range(args.pt_range)
    et_edges = parse_edges(args.et_bins)
    cent_edges = parse_edges(args.cent_bins)
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.score_cache_list is not None:
        samples = load_cache_samples(args.score_cache_list)
        mode = "capped_score_cache_fallback"
    else:
        grouped = collect_root_files(args.source, args.root_files)
        if not grouped:
            raise SystemExit(f"No ROOT files found under {args.source}")
        samples = [
            load_root_sample(sample, paths, args.tree, max_files=args.max_files_per_sample)
            for sample, paths in sorted(grouped.items())
        ]
        mode = "exact_root_training_extraction"
    samples = combine_samples(samples)
    sample_rows = [row_for_sample(s, pt_range, args.test_size, args.random_seed) for s in samples]
    sample_rows.sort(key=lambda r: SAMPLE_ORDER.index(str(r["sample"])) if r["sample"] in SAMPLE_ORDER else 999)
    bins = heatmap_rows(samples, et_edges, cent_edges, pt_range)

    write_csv(args.outdir / "single_bdt_configuration.csv", config_rows(args.test_size, args.random_seed))
    write_csv(args.outdir / "sample_counts_single_bdt_centInput_pt1535.csv", sample_rows)
    write_csv(args.outdir / "event_candidate_summary_single_bdt_centInput_pt1535.csv", event_candidate_rows(sample_rows))
    write_csv(args.outdir / "counts_by_et_cent_single_bdt_centInput_pt1535.csv", bins)
    summary = {
        "schema": "RJ_AUAU_SINGLE_BDT_SAMPLE_AUDIT_V1",
        "mode": mode,
        "product": "centInput_pt1535",
        "pt_range": list(pt_range),
        "et_edges": et_edges,
        "cent_edges": cent_edges,
        "test_size": args.test_size,
        "random_seed": args.random_seed,
        "outdir": str(args.outdir),
        "samples": sample_rows,
        "caveat": "score-cache mode is capped and lacks per-sample provenance; exact slide counts require ROOT mode." if mode != "exact_root_training_extraction" else "",
    }
    (args.outdir / "audit_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (args.outdir / "README.md").write_text(
        "# Single-BDT AuAu Sample Audit\n\n"
        f"Mode: `{mode}`\n\n"
        "Baseline: `centInput_pt1535`, one global BDT over `15 < E_T < 35 GeV` "
        "with `cluster_Et` and `centrality` as features. No isolation inputs.\n"
    )

    if not args.no_plots:
        plot_truth_composition(
            sample_rows,
            args.outdir / "truth_label_composition_selected_candidates.png",
        )
        for sample, value, title in [
            ("embeddedPhoton12plus20", "signal", r"Signal counts for single global BDT, $15<E_T<35$ GeV"),
            ("embeddedJet12plus20", "background", r"Background counts for single global BDT, $15<E_T<35$ GeV"),
            ("embeddedJet12plus20plus30", "background", r"Background counts for single global BDT, $15<E_T<35$ GeV"),
        ]:
            if any(row["sample"] == sample for row in bins):
                plot_heatmap(
                    bins,
                    sample,
                    value,
                    et_edges,
                    cent_edges,
                    args.outdir / f"{sample}_{value}_counts_et_cent_heatmap.png",
                    title,
                )
    print(f"[singleBDTAudit] mode={mode} wrote {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
