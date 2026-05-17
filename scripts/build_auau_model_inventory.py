#!/usr/bin/env python3
"""Build a local AuAu photon-ID model inventory.

This script scans pulled validation/rank artifacts and writes a compact model
inventory plus a desired coverage matrix.  It is intentionally read-only with
respect to physics artifacts: it only creates summary files under
``dataOutput/auauModelInventory``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any


DEFAULT_ROOT = Path("dataOutput")
DEFAULT_OUT = Path("dataOutput/auauModelInventory")

PRIMARY_RANGE = "15-35"
SIDE_RANGES = ("5-35",)
PHYSICS_TOKEN = "auauPID"
DEFAULT_VERSION = "v001"
ROUTINGS = (
    "global_centEtInput",
    "global_etInput_noCent",
    "global_centInput_noEt",
    "etCoarse_centInput",
    "etFine_centInput",
    "cent3_etInput",
    "cent7_etInput",
    "etCoarse_cent3_routed",
    "etFine_cent3_routed",
    "etCoarse_cent7_routed",
    "etFine_cent7_routed",
)
DIRECT_ALGORITHMS = ("BDT", "MLP", "LogReg")
STACK_ALGORITHMS = ("StackLogistic", "StackGBM", "StackNN")
ALGORITHMS = (*DIRECT_ALGORITHMS, *STACK_ALGORITHMS, "OOFHeavyNN", "IsoDiagnostic")
COHORTS: dict[str, dict[str, Any]] = {
    "direct_full_no_iso": {
        "token": "directFull",
        "priority": "core",
        "algorithms": DIRECT_ALGORITHMS,
        "official_candidate": "yes",
        "diagnostic": "no",
    },
    "stack_pair_score_only": {
        "token": "stackPairScore",
        "priority": "core_control",
        "algorithms": STACK_ALGORITHMS,
        "official_candidate": "closure_pending",
        "diagnostic": "score-input closure pending",
    },
    "stack_tri_score_only": {
        "token": "stackTriScore",
        "priority": "core_control",
        "algorithms": STACK_ALGORITHMS,
        "official_candidate": "closure_pending",
        "diagnostic": "score-input closure pending",
    },
    "stack_pair_context": {
        "token": "stackPairCtx",
        "priority": "calibration",
        "algorithms": STACK_ALGORITHMS,
        "official_candidate": "closure_pending",
        "diagnostic": "score-input closure pending",
    },
    "stack_tri_context": {
        "token": "stackTriCtx",
        "priority": "calibration",
        "algorithms": STACK_ALGORITHMS,
        "official_candidate": "closure_pending",
        "diagnostic": "score-input closure pending",
    },
    "stack_pair_full_no_iso": {
        "token": "stackPairFull",
        "priority": "high_ceiling",
        "algorithms": STACK_ALGORITHMS,
        "official_candidate": "closure_pending",
        "diagnostic": "score-input closure pending",
    },
    "stack_tri_full_no_iso": {
        "token": "stackTriFull",
        "priority": "high_ceiling",
        "algorithms": STACK_ALGORITHMS,
        "official_candidate": "closure_pending",
        "diagnostic": "score-input closure pending",
    },
    "stack_bdt_score_full_no_iso": {
        "token": "stackBDTFull",
        "priority": "sanity_control",
        "algorithms": ("StackGBM", "StackNN"),
        "official_candidate": "closure_pending",
        "diagnostic": "score-input closure pending",
    },
    "stack_mlp_score_full_no_iso": {
        "token": "stackMLPFull",
        "priority": "sanity_control",
        "algorithms": ("StackGBM", "StackNN"),
        "official_candidate": "closure_pending",
        "diagnostic": "score-input closure pending",
    },
    "oof_tri_heavy_full_no_iso": {
        "token": "oofTriHeavy",
        "priority": "overfit_audit",
        "algorithms": ("OOFHeavyNN",),
        "official_candidate": "diagnostic",
        "diagnostic": "out-of-fold superlearner audit",
    },
    "iso_visible_diagnostic": {
        "token": "isoDiag",
        "priority": "diagnostic_only",
        "algorithms": ("BDT", "MLP", "StackLogistic", "StackGBM", "StackNN", "IsoDiagnostic"),
        "official_candidate": "no",
        "diagnostic": "uses reconstructed isolation-derived inputs",
    },
}
ACTIVE_COHORTS = tuple(COHORTS)
RANGE_TOKENS = {
    "15-35": "pt1535",
    "5-35": "pt0535",
    "5-40": "pt0540",
    "10-35": "pt1035",
    "15-30": "pt1530",
    "unknown": "ptUnknown",
}
ROUTING_TOKENS = {
    "global_centEtInput": "rGlobalCEt",
    "global_etInput_noCent": "rGlobalE",
    "global_centInput_noEt": "rGlobalC",
    "etCoarse_centInput": "rEC_CEt",
    "etFine_centInput": "rEF_CEt",
    "cent3_etInput": "rC3_E",
    "cent7_etInput": "rC7_E",
    "etCoarse_cent3_routed": "rECxC3",
    "etFine_cent3_routed": "rEFxC3",
    "etCoarse_cent7_routed": "rECxC7",
    "etFine_cent7_routed": "rEFxC7",
    "unknown": "rUnknown",
}
ALGORITHM_TOKENS = {
    "BDT": "bdt",
    "MLP": "mlp",
    "LogReg": "logreg",
    "StackLogistic": "stackLogistic",
    "StackGBM": "stackGBM",
    "StackNN": "stackNN",
    "Stack": "stack",
    "OOFHeavyNN": "oofHeavyNN",
    "IsoDiagnostic": "isoDiag",
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    return ap.parse_args()


def read_json(path: Path) -> Any | None:
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def fnum(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def fmt(value: Any) -> str:
    val = fnum(value)
    if val is None:
        return ""
    return f"{val:.6g}"


def metric_from_bins(bins: list[dict[str, Any]] | None, lo: float, hi: float, key: str = "auc") -> float | None:
    if not bins:
        return None
    for item in bins:
        if abs(float(item.get("lo", -999)) - lo) < 1e-6 and abs(float(item.get("hi", -999)) - hi) < 1e-6:
            return fnum(item.get(key))
    return None


def wp_fake_from_bins(bins: list[dict[str, Any]] | None, lo: float, hi: float) -> float | None:
    if not bins:
        return None
    for item in bins:
        if abs(float(item.get("lo", -999)) - lo) < 1e-6 and abs(float(item.get("hi", -999)) - hi) < 1e-6:
            wp = item.get("wp80") if isinstance(item.get("wp80"), dict) else item
            return fnum(wp.get("background_fake_rate", wp.get("wp80_fake")) if isinstance(wp, dict) else None)
    return None


def infer_range(name: str, path: Path) -> str:
    text = f"{name} {path}".lower()
    patterns = [
        ("15-35", ("pt1535", "15to35", "15-35")),
        ("15-30", ("pt1530", "15to30", "15-30")),
        ("10-35", ("pt10to35", "10to35", "10-35")),
        ("5-35", ("pt5to35", "5to35", "5-35")),
        ("5-40", ("pt5to40", "5to40", "allrange", "all-range", "expanded")),
    ]
    for label, needles in patterns:
        if any(n in text for n in needles):
            return label
    return "unknown"


def range_hint_from_products(products: dict[str, Any], path: Path) -> str:
    """Infer a report-wide pT range for rows whose product name is terse.

    Some BDT reports contain one explicit product such as ``centInput_pt1535``
    and several companion rows named only ``ptFine_cent7``.  Those companion
    rows belong to the same training/validation range, so using the report-wide
    hint prevents them from being filed under ``unknown``.
    """

    counts: dict[str, int] = {}
    path_hint = infer_range("", path)
    if path_hint != "unknown":
        counts[path_hint] = counts.get(path_hint, 0) + 1
    for product in products:
        hint = infer_range(str(product), path)
        if hint != "unknown":
            counts[hint] = counts.get(hint, 0) + 1
    if not counts:
        return "unknown"
    priority = {"15-35": 0, "5-35": 1, "5-40": 2, "10-35": 3, "15-30": 4}
    return sorted(counts, key=lambda key: (-counts[key], priority.get(key, 99), key))[0]


def infer_algorithm(name: str, path: Path, explicit: str | None = None) -> str:
    if explicit:
        low = explicit.lower()
        if low == "logistic":
            return "StackLogistic" if "stack" in str(path).lower() else "LogReg"
        if low == "gbm":
            return "StackGBM"
        if low == "nn":
            return "StackNN"
    text = f"{name} {path}".lower()
    if "oof" in text and ("super" in text or "heavy" in text):
        return "OOFHeavyNN"
    if "iso_visible" in text or "isolation_visible" in text or "isodiag" in text:
        if "stack" not in text and "bdt" not in text and "mlp" not in text:
            return "IsoDiagnostic"
    if "stack" in text:
        if "gbm" in text:
            return "StackGBM"
        if "nn" in text:
            return "StackNN"
        if "logistic" in text:
            return "StackLogistic"
        return "Stack"
    if "logreg" in text or "logistic" in text:
        return "LogReg"
    if "mlp" in text:
        return "MLP"
    return "BDT"


def infer_feature_family(name: str, path: Path, algorithm: str) -> str:
    text = f"{name} {path}".lower()
    if algorithm.startswith("Stack"):
        return "stackScores_v1"
    if "iso_kitchen" in text or "isokitchen" in text or "iso_diag" in text:
        return "fullIsoDiag_v1"
    if "widthratios" in text or "primary_ratios" in text or "primary-ratios" in text:
        return "widthRatios_v1"
    if "full" in text or "kitchen" in text or "logreg" in text:
        return "fullNoIso_v1"
    if "base3x3" in text or "3x3" in text or "ptfine" in text or "centinput" in text:
        return "base3x3_v1"
    if algorithm == "BDT":
        return "baseBDT_v1"
    if algorithm in {"MLP", "LogReg"}:
        return "fullNoIso_v1"
    return "unknownFeatures"


def split_feature_string(value: Any) -> set[str]:
    if value is None:
        return set()
    if isinstance(value, list):
        return {str(item) for item in value}
    text = str(value)
    if not text:
        return set()
    return {item.strip() for item in re.split(r"[+,; ]+", text) if item.strip()}


def infer_input_cohort(row: dict[str, Any]) -> str:
    """Classify a row by the information made available to the model.

    This is intentionally stricter than the older feature-family field.  It
    separates direct photon-ID models from score-only stackers, score+context
    calibration, full-feature stackers, one-score rescue controls, and
    isolation-visible diagnostics.
    """

    explicit = str(row.get("input_cohort") or "").strip()
    if explicit in COHORTS:
        return explicit

    algorithm = str(row.get("algorithm") or "")
    model = str(row.get("model") or "")
    source = str(row.get("source") or "")
    context = str(row.get("context") or "")
    text = f"{model} {source} {context}".lower()
    name_context = f"{model} {context}".lower()
    features = split_feature_string(row.get("features"))
    lower_features = {feature.lower() for feature in features}

    if any("reco_eiso" in feature or "eiso" in feature for feature in lower_features) or any(
        token in text for token in ("iso_visible", "isolation_visible", "iso_context", "isokitchen", "iso_kitchen")
    ):
        return "iso_visible_diagnostic"
    if algorithm == "OOFHeavyNN" or "oof_triscore" in text or "superlearner" in text:
        return "oof_tri_heavy_full_no_iso"
    if not algorithm.startswith("Stack"):
        return "direct_full_no_iso"

    has_bdt = "bdt_score" in lower_features or "bdtonly" in name_context or "bdt_mlp" in name_context
    has_mlp = "mlp_score" in lower_features or "mlponly" in name_context or "mlpscore" in name_context or "bdt_mlp" in name_context
    has_logreg = "logreg_score" in lower_features or "logreg" in name_context or "triscore" in name_context or "tri_score" in name_context
    has_full = any(
        feature in lower_features
        for feature in (
            "cluster_weta_cogx",
            "cluster_wphi_cogx",
            "cluster_weta33_cogx",
            "cluster_wphi33_cogx",
            "cluster_weta35_cogx",
            "cluster_wphi53_cogx",
            "cluster_w32",
            "e11_over_e33",
        )
    ) or "full" in text
    has_context = any(feature in lower_features for feature in ("cluster_et", "log_cluster_et", "centrality", "centrality_scaled"))
    one_score_only = ("mlpscore" in name_context or "mlp_score_full" in name_context or "mlp score" in name_context) and not has_bdt
    bdt_score_only = ("bdtscore" in name_context or "bdt_score_full" in name_context or "bdt score" in name_context) and not has_mlp

    if one_score_only and has_full:
        return "stack_mlp_score_full_no_iso"
    if bdt_score_only and has_full:
        return "stack_bdt_score_full_no_iso"
    if has_logreg and has_full:
        return "stack_tri_full_no_iso"
    if has_logreg and has_context:
        return "stack_tri_context"
    if has_logreg:
        return "stack_tri_score_only"
    if has_bdt and has_mlp and has_full:
        return "stack_pair_full_no_iso"
    if has_bdt and has_mlp and has_context:
        return "stack_pair_context"
    return "stack_pair_score_only"


def infer_train_scope(name: str, path: Path) -> str:
    text = f"{name} {path}".lower()
    if "smoke" in text or "tiny" in text:
        return "smoke"
    if "diagnostic" in text or "iso_kitchen" in text:
        return "diagnostic"
    if "capped" in text or "max_rows" in text:
        return "capped"
    if any(token in text for token in ("fullstat", "fullval", "condor", "target80", "wp80")):
        return "fullstat"
    return "unknownScope"


def infer_routing(name: str, path: Path, variant: dict[str, Any] | None = None) -> str:
    text = f"{name} {path}".lower()
    routes = None
    if variant:
        routes = variant.get("routes")

    labels = ""
    route_count = 0
    if routes:
        route_count = len(routes)
        labels = " ".join(str(r.get("label", "")) for r in routes if isinstance(r, dict)).lower()
        if "pt_" in labels and "cent_" in labels:
            return "etFine_cent7_routed" if route_count >= 40 else "etFine_cent3_routed"

    has_pt_coarse = any(token in text for token in ("ptcoarse", "etcoarse"))
    has_pt_fine = any(token in text for token in ("ptfine", "etfine", "ptbin", "ptcentdep"))
    has_pt_routing = has_pt_coarse or has_pt_fine
    has_cent7 = any(token in text for token in ("cent7", "centdepfine", "ptcentdepfine"))
    has_cent3 = any(token in text for token in ("cent3", "centdepbdts", "ptcentdep3"))

    if has_pt_coarse and has_cent7:
        return "etCoarse_cent7_routed"
    if has_pt_coarse and has_cent3:
        return "etCoarse_cent3_routed"
    if has_pt_coarse and any(token in text for token in ("centinput", "centasfeat")):
        return "etCoarse_centInput"
    if has_pt_fine and has_cent7:
        return "etFine_cent7_routed"
    if has_pt_fine and has_cent3:
        return "etFine_cent3_routed"
    if has_pt_fine and any(token in text for token in ("centinput", "centasfeat")):
        return "etFine_centInput"
    if has_cent7:
        return "cent7_etInput"
    if has_cent3:
        return "cent3_etInput"
    if "nocent" in text or "centind" in text:
        return "global_etInput_noCent"
    if "noet" in text or "etind" in text or "centonly" in text or "centralityonly" in text:
        return "global_centInput_noEt"
    if "centinput" in text or "centasfeat" in text or "centet" in text or "etcent" in text:
        return "global_centEtInput"
    return "unknown"


def row_base(
    *,
    algorithm: str,
    model: str,
    source: Path,
    split: str = "validation",
    routing: str = "unknown",
    training_range: str = "unknown",
    status: str = "local",
    stack_diagnostic: bool = False,
    notes: str = "",
) -> dict[str, Any]:
    return {
        "canonical_id": "",
        "comparable_group": "",
        "physics": PHYSICS_TOKEN,
        "pt_token": "",
        "input_cohort": "",
        "cohort_token": "",
        "cohort_priority": "",
        "official_candidate": "",
        "feature_family": "",
        "routing_token": "",
        "algorithm_token": "",
        "train_scope": "",
        "version": DEFAULT_VERSION,
        "readiness": "",
        "readiness_reason": "",
        "algorithm": algorithm,
        "model": model,
        "training_range": training_range,
        "routing": routing,
        "context": "",
        "routes": "",
        "features": "",
        "split": split,
        "status": status,
        "stack_diagnostic": "yes" if stack_diagnostic else "no",
        "diagnostic_reason": "",
        "entries": "",
        "scored_entries": "",
        "auc": "",
        "wp80_fake": "",
        "finite_fraction": "",
        "ece": "",
        "score_vs_eiso_corr": "",
        "highpt_auc_20_35": "",
        "highpt_wp80_fake_20_35": "",
        "pt_15_20_auc": "",
        "pt_20_25_auc": "",
        "pt_25_35_auc": "",
        "pt_15_18_auc": "",
        "pt_18_20_auc": "",
        "pt_20_22p5_auc": "",
        "pt_22p5_25_auc": "",
        "pt_25_30_auc": "",
        "pt_30_35_auc": "",
        "source": str(source),
        "notes": notes,
    }


def parse_bdt_validation(path: Path) -> list[dict[str, Any]]:
    data = read_json(path)
    if not isinstance(data, dict) or "products" not in data:
        return []
    if "/shards/" in str(path):
        return []
    rows: list[dict[str, Any]] = []
    products = data.get("products", {})
    range_hint = range_hint_from_products(products if isinstance(products, dict) else {}, path)
    for product, metrics in products.items():
        if not isinstance(metrics, dict):
            continue
        training_range = infer_range(product, path)
        if training_range == "unknown":
            training_range = range_hint
        row = row_base(
            algorithm="BDT",
            model=product,
            source=path,
            routing=infer_routing(product, path),
            training_range=training_range,
            status=str(data.get("status", "local")),
        )
        row["auc"] = fmt(metrics.get("auc_inclusive"))
        row["finite_fraction"] = fmt(metrics.get("finite_score_fraction"))
        row["score_vs_eiso_corr"] = fmt(metrics.get("score_eiso_pearson"))
        auc_by_pt = metrics.get("auc_by_pt") if isinstance(metrics.get("auc_by_pt"), dict) else {}
        row["pt_15_20_auc"] = fmt(auc_by_pt.get("15_20"))
        row["pt_20_25_auc"] = fmt(auc_by_pt.get("20_25"))
        row["pt_25_35_auc"] = fmt(auc_by_pt.get("25_35"))
        rows.append(row)
    return rows


def parse_mlp_validation(path: Path) -> list[dict[str, Any]]:
    data = read_json(path)
    if not isinstance(data, dict) or "products" not in data:
        return []
    if "/shards/" in str(path):
        return []
    rows: list[dict[str, Any]] = []
    products = data.get("products", {})
    range_hint = range_hint_from_products(products if isinstance(products, dict) else {}, path)
    for product, metrics in products.items():
        if not isinstance(metrics, dict) or "MLP" not in product:
            continue
        training_range = infer_range(product, path)
        if training_range == "unknown":
            training_range = range_hint
        row = row_base(
            algorithm="MLP",
            model=product,
            source=path,
            routing=infer_routing(product, path),
            training_range=training_range,
            status="READY",
        )
        row["auc"] = fmt(metrics.get("auc"))
        wp80 = (metrics.get("thresholds") or {}).get("wp080", {})
        row["wp80_fake"] = fmt(wp80.get("background_fake_rate"))
        row["finite_fraction"] = fmt(metrics.get("finite_fraction"))
        row["ece"] = fmt((metrics.get("calibration") or {}).get("ece"))
        row["score_vs_eiso_corr"] = fmt((metrics.get("correlations") or {}).get("score_vs_reco_eiso"))
        pt_bins = metrics.get("pt_bins")
        row["pt_15_20_auc"] = fmt(metric_from_bins(pt_bins, 15.0, 20.0))
        row["pt_20_25_auc"] = fmt(metric_from_bins(pt_bins, 20.0, 25.0))
        row["pt_25_35_auc"] = fmt(metric_from_bins(pt_bins, 25.0, 35.0))
        rows.append(row)
    return rows


def parse_stack_rank(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for src in reader:
            algorithm = infer_algorithm(src.get("model", ""), path, src.get("algorithm"))
            variant = src.get("variant") or src.get("model") or ""
            row = row_base(
                algorithm=algorithm,
                model=src.get("model", ""),
                source=path,
                split=src.get("split", "unknown"),
                routing=src.get("routing") or infer_routing(variant, path),
                training_range=src.get("pt_range", "").replace(":", "-") or infer_range(variant, path),
                status="READY",
                stack_diagnostic=True,
                notes="uses BDT/MLP score input",
            )
            row["input_cohort"] = src.get("input_cohort", "")
            row["context"] = src.get("context", "")
            row["routes"] = src.get("routes", "")
            row["features"] = src.get("features", "")
            row["readiness"] = src.get("readiness", "")
            row["stack_training_safety"] = src.get("stack_training_safety", "")
            row["full_stat_required"] = src.get("full_stat_required", "")
            row["matrix_wave_name"] = src.get("matrix_wave_name", "")
            row["entries"] = fmt(src.get("entries"))
            row["scored_entries"] = fmt(src.get("scored_entries"))
            for key in (
                "auc",
                "wp80_fake",
                "finite_fraction",
                "ece",
                "score_vs_eiso_corr",
                "highpt_auc_20_35",
                "highpt_wp80_fake_20_35",
                "pt_15_18_auc",
                "pt_18_20_auc",
                "pt_20_22p5_auc",
                "pt_22p5_25_auc",
                "pt_25_30_auc",
                "pt_30_35_auc",
            ):
                row[key] = fmt(src.get(key))
            rows.append(row)
    return rows


def parse_generic_rank(path: Path) -> list[dict[str, Any]]:
    if path.name == "stacked_sweep_rank_table.csv":
        return parse_stack_rank(path)
    rows: list[dict[str, Any]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        if not {"product", "auc"} & fields and not {"model", "auc"} & fields:
            return []
        for src in reader:
            model = src.get("product") or src.get("model") or src.get("candidate") or "unknown"
            algorithm = infer_algorithm(model, path)
            row = row_base(
                algorithm=algorithm,
                model=model,
                source=path,
                split=src.get("split", "validation"),
                routing=infer_routing(model, path),
                training_range=src.get("pt_range", "").replace(":", "-") or infer_range(model, path),
                status="READY",
                stack_diagnostic=algorithm.startswith("Stack"),
            )
            row["entries"] = fmt(src.get("entries"))
            row["scored_entries"] = fmt(src.get("scored_entries"))
            for dst, candidates in {
                "auc": ("auc", "inclusive_auc"),
                "wp80_fake": ("wp80_fake", "wp080_fake", "background_fake_rate"),
                "finite_fraction": ("finite_fraction",),
                "ece": ("ece",),
                "score_vs_eiso_corr": ("score_vs_eiso_corr", "score_vs_reco_eiso_corr"),
                "highpt_auc_20_35": ("highpt_auc_20_35",),
                "highpt_wp80_fake_20_35": ("highpt_wp80_fake_20_35",),
            }.items():
                for cand in candidates:
                    if cand in src:
                        row[dst] = fmt(src.get(cand))
                        break
            rows.append(row)
    return rows


def discover_rows(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(root.rglob("validation_metrics.json")):
        if "/shards/" in str(path):
            continue
        data = read_json(path)
        if not isinstance(data, dict):
            continue
        schema = str(data.get("schema", ""))
        products = data.get("products", {})
        if "MLP" in schema or any("MLP" in key for key in products):
            rows.extend(parse_mlp_validation(path))
        else:
            rows.extend(parse_bdt_validation(path))
    for path in sorted(root.rglob("*rank_table*.csv")):
        rows.extend(parse_generic_rank(path))
    return rows


def readiness_for_row(row: dict[str, Any]) -> tuple[str, str]:
    explicit = str(row.get("readiness") or "")
    if explicit in {"smoke_passed", "trained_full", "validated_full", "validated_full_oof_safe", "diagnostic_non_oof", "blocked_or_check"}:
        return explicit, "readiness supplied by matrix/validation artifact"
    reasons: list[str] = []
    finite = fnum(row.get("finite_fraction"))
    if finite is not None and finite < 0.95:
        reasons.append(f"finite_fraction={finite:.3f}<0.95")
    if row.get("stack_diagnostic") == "yes":
        safety = str(row.get("stack_training_safety") or "")
        full_required = str(row.get("full_stat_required") or "").lower() in {"true", "1", "yes"}
        if full_required and safety in {"disjoint_base_scores_heldout_test", "oof_base_scores_heldout_test"}:
            return "validated_full_oof_safe", "full-stat stack with disjoint/OOF-safe first-stage score provenance and held-out test split"
        if full_required:
            return "diagnostic_non_oof", "full-stat score-input stack but not marked OOF/disjoint-safe"
        reasons.append("score-input stack; ABCD/runtime closure pending")
    if row.get("training_range") == "unknown":
        reasons.append("training range inferred as unknown")
    if row.get("routing") == "unknown":
        reasons.append("routing inferred as unknown")
    if row.get("feature_family") == "unknownFeatures":
        reasons.append("feature family inferred as unknown")
    if reasons:
        return "check", "; ".join(reasons)
    return "validated_local", "local validation or rank evidence present"


def canonicalize_row(row: dict[str, Any]) -> dict[str, Any]:
    algorithm = str(row.get("algorithm") or "unknown")
    model = str(row.get("model") or "unknown")
    source = Path(str(row.get("source") or ""))
    training_range = str(row.get("training_range") or "unknown")
    routing = str(row.get("routing") or "unknown")
    feature_family = infer_feature_family(model, source, algorithm)
    input_cohort = infer_input_cohort(row)
    cohort = COHORTS.get(input_cohort, {})
    pt_token = RANGE_TOKENS.get(training_range, "ptUnknown")
    routing_token = ROUTING_TOKENS.get(routing, "rUnknown")
    algorithm_token = ALGORITHM_TOKENS.get(algorithm, algorithm.lower())
    cohort_token = str(cohort.get("token", input_cohort))
    train_scope = infer_train_scope(model, source)
    comparable_group = "__".join([PHYSICS_TOKEN, pt_token, cohort_token, routing_token])
    canonical_id = "__".join([comparable_group, algorithm_token, train_scope, DEFAULT_VERSION])
    row["physics"] = PHYSICS_TOKEN
    row["pt_token"] = pt_token
    row["input_cohort"] = input_cohort
    row["cohort_token"] = cohort_token
    row["cohort_priority"] = str(cohort.get("priority", "unknown"))
    row["official_candidate"] = str(cohort.get("official_candidate", "unknown"))
    row["feature_family"] = feature_family
    row["routing_token"] = routing_token
    row["algorithm_token"] = algorithm_token
    row["train_scope"] = train_scope
    row["version"] = DEFAULT_VERSION
    row["comparable_group"] = comparable_group
    row["canonical_id"] = canonical_id
    row["diagnostic_reason"] = str(cohort.get("diagnostic", ""))
    readiness, reason = readiness_for_row(row)
    row["readiness"] = readiness
    row["readiness_reason"] = reason
    return row


def canonicalize_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [canonicalize_row(row) for row in rows]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_canonical_index(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "canonical_id",
        "comparable_group",
        "readiness",
        "readiness_reason",
        "physics",
        "training_range",
        "pt_token",
        "input_cohort",
        "cohort_token",
        "cohort_priority",
        "official_candidate",
        "diagnostic_reason",
        "feature_family",
        "routing",
        "routing_token",
        "algorithm",
        "algorithm_token",
        "train_scope",
        "version",
        "model",
        "context",
        "routes",
        "split",
        "auc",
        "wp80_fake",
        "finite_fraction",
        "ece",
        "score_vs_eiso_corr",
        "highpt_auc_20_35",
        "highpt_wp80_fake_20_35",
        "source",
        "notes",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (
            str(r.get("comparable_group", "")),
            str(r.get("algorithm_token", "")),
            str(r.get("split", "")),
            str(r.get("model", "")),
        )):
            writer.writerow({field: row.get(field, "") for field in fields})


def routings_for_cohort(input_cohort: str) -> tuple[str, ...]:
    if input_cohort == "oof_tri_heavy_full_no_iso":
        return ("global_centEtInput",)
    if input_cohort == "iso_visible_diagnostic":
        return ("global_centEtInput", "etFine_cent7_routed")
    return ROUTINGS


def desired_matrix(rows: list[dict[str, Any]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    present: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    for row in rows:
        rng = row.get("training_range") or "unknown"
        routing = row.get("routing") or "unknown"
        algo = row.get("algorithm") or "unknown"
        input_cohort = row.get("input_cohort") or "unknown"
        if rng != "unknown" and routing != "unknown" and input_cohort != "unknown":
            key = (rng, input_cohort, routing, algo)
            old = present.get(key)
            if old is None or (fnum(row.get("auc")) or -1) > (fnum(old.get("auc")) or -1):
                present[key] = row
    desired_rows: list[dict[str, str]] = []
    missing: list[dict[str, str]] = []
    for rng in (PRIMARY_RANGE,) + SIDE_RANGES:
        for input_cohort in ACTIVE_COHORTS:
            cohort = COHORTS[input_cohort]
            for routing in routings_for_cohort(input_cohort):
                for algo in cohort["algorithms"]:
                    range_priority = "primary" if rng == PRIMARY_RANGE else "side"
                    present_row = present.get((rng, input_cohort, routing, algo))
                    if present_row is None:
                        status = "missing"
                    else:
                        status = str(present_row.get("readiness") or "validated_local")
                    item = {
                        "training_range": rng,
                        "input_cohort": input_cohort,
                        "cohort_priority": str(cohort["priority"]),
                        "routing": routing,
                        "algorithm": str(algo),
                        "range_priority": range_priority,
                        "status": status,
                    }
                    desired_rows.append(item)
                    if status == "missing":
                        missing.append(item)
    return desired_rows, missing


def coverage_matrix_markdown(desired: list[dict[str, str]]) -> str:
    by_key = {
        (row["training_range"], row["input_cohort"], row["routing"], row["algorithm"]): row["status"]
        for row in desired
    }
    lines: list[str] = []
    legend = "`validated_local` = local validation/rank evidence exists; `check` = artifact exists but has a validation caveat; `missing` = no local validated row."
    lines.append(legend)
    for rng in (PRIMARY_RANGE,) + SIDE_RANGES:
        lines.append(f"\n### {rng} GeV\n")
        for input_cohort in ACTIVE_COHORTS:
            cohort = COHORTS[input_cohort]
            fields = ["routing", *cohort["algorithms"]]
            rows = []
            for routing in routings_for_cohort(input_cohort):
                item = {"routing": routing}
                for algo in cohort["algorithms"]:
                    item[str(algo)] = by_key.get((rng, input_cohort, routing, str(algo)), "missing")
                rows.append(item)
            lines.append(f"\n#### {input_cohort} ({cohort['priority']})\n")
            lines.append(md_table(rows, fields))
    return "\n".join(lines)


def cohort_summary(desired: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], dict[str, int]] = {}
    for row in desired:
        key = (row["training_range"], row["input_cohort"])
        item = grouped.setdefault(key, {"total": 0, "missing": 0, "check": 0, "validated": 0})
        item["total"] += 1
        status = row["status"]
        if status == "missing":
            item["missing"] += 1
        elif status == "check":
            item["check"] += 1
        else:
            item["validated"] += 1
    rows: list[dict[str, str]] = []
    for (rng, input_cohort), counts in sorted(grouped.items()):
        cohort = COHORTS.get(input_cohort, {})
        total = counts["total"]
        ready = counts["validated"]
        rows.append(
            {
                "training_range": rng,
                "input_cohort": input_cohort,
                "cohort_priority": str(cohort.get("priority", "")),
                "validated_or_ready": str(ready),
                "check": str(counts["check"]),
                "missing": str(counts["missing"]),
                "total": str(total),
                "coverage_fraction": f"{ready / total:.3f}" if total else "",
            }
        )
    return rows


def write_training_backlog(out: Path, desired: list[dict[str, str]]) -> None:
    stages = [
        ("Stage 1: primary direct-model gaps", "15-35", {"core"}),
        ("Stage 2: primary stack controls and calibration gaps", "15-35", {"core_control", "calibration"}),
        ("Stage 3: primary high-ceiling stack gaps", "15-35", {"high_ceiling", "sanity_control", "overfit_audit"}),
        ("Stage 4: 5-35 side-study gaps", "5-35", {"core", "core_control", "calibration", "high_ceiling", "sanity_control"}),
    ]
    lines = [
        "# AuAu Model Training Backlog",
        "",
        "Generated by `scripts/build_auau_model_inventory.py`.",
        "",
        "Rows marked `missing` need training/validation. Rows marked `check` need repair,",
        "rerouting, finite-fraction fixes, or closure review before they can be used as",
        "clean comparison evidence.",
        "",
    ]
    for title, rng, cohort_priorities in stages:
        rows = [
            row for row in desired
            if row["training_range"] == rng
            and row["cohort_priority"] in cohort_priorities
            and row["status"] in {"missing", "check"}
        ]
        lines.append(f"## {title}")
        lines.append("")
        if not rows:
            lines.append("No open rows.")
            lines.append("")
            continue
        lines.append(md_table(rows, ["training_range", "input_cohort", "cohort_priority", "routing", "algorithm", "range_priority", "status"]))
        lines.append("")
    out.write_text("\n".join(lines))


def best_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keyed: dict[tuple[str, str, str, str, str], dict[str, Any]] = {}
    for row in rows:
        key = (
            row.get("training_range", ""),
            row.get("input_cohort", ""),
            row.get("routing", ""),
            row.get("algorithm", ""),
            row.get("split", ""),
        )
        old = keyed.get(key)
        if old is None or (fnum(row.get("auc")) or -1) > (fnum(old.get("auc")) or -1):
            keyed[key] = row
    return sorted(keyed.values(), key=lambda r: (
        r.get("training_range", ""),
        r.get("input_cohort", ""),
        r.get("routing", ""),
        r.get("algorithm", ""),
        r.get("split", ""),
    ))


def md_table(rows: list[dict[str, Any]], fields: list[str]) -> str:
    lines = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for row in rows:
        vals = [str(row.get(field, "")).replace("|", "\\|") for field in fields]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_markdown(out: Path, rows: list[dict[str, Any]], desired: list[dict[str, str]], missing: list[dict[str, str]]) -> None:
    active_ranges = {PRIMARY_RANGE, *SIDE_RANGES}
    summary_rows = [row for row in best_rows(rows) if row.get("training_range") in active_ranges]
    cohort_rows = cohort_summary(desired)
    text = []
    text.append("# AuAu Model Inventory\n")
    text.append("Generated by `scripts/build_auau_model_inventory.py` from local pulled artifacts.\n")
    text.append("This is an inventory, not a physics endorsement. Compare models first within matching training range, input cohort, and routing rows.\n")
    text.append("Active comparison ranges are `15-35` and `5-35`; archival rows such as older `5-40` BDT outputs remain in the CSV but are excluded from the active Markdown tables.\n")
    text.append("## Cohort Coverage Summary\n")
    text.append(md_table(cohort_rows, ["training_range", "input_cohort", "cohort_priority", "validated_or_ready", "check", "missing", "total", "coverage_fraction"]))
    text.append("\n## How To Read The Cohorts\n")
    text.append("- `direct_full_no_iso`: direct BDT/MLP/LogReg photon-ID classifiers with ABCD-safe inputs.\n")
    text.append("- `stack_pair_score_only` / `stack_tri_score_only`: only first-round classifier scores; pure complementarity/calibration checks.\n")
    text.append("- `stack_pair_context` / `stack_tri_context`: scores plus explicit `E_T`/centrality calibration context.\n")
    text.append("- `stack_pair_full_no_iso` / `stack_tri_full_no_iso`: scores plus full ABCD-safe shower-shape and energy-sharing features; highest non-isolation ceiling stackers.\n")
    text.append("- `stack_bdt_score_full_no_iso` / `stack_mlp_score_full_no_iso`: one-score rescue sanity controls.\n")
    text.append("- `iso_visible_diagnostic`: isolation-visible ceiling tests; not ABCD-safe.\n")
    text.append("## Best Local Rows By Range / Routing / Algorithm / Split\n")
    text.append(md_table(summary_rows, [
        "canonical_id",
        "readiness",
        "training_range",
        "input_cohort",
        "cohort_priority",
        "feature_family",
        "routing",
        "algorithm",
        "model",
        "split",
        "auc",
        "wp80_fake",
        "finite_fraction",
        "highpt_auc_20_35",
        "highpt_wp80_fake_20_35",
        "source",
    ]))
    text.append("\n## Desired Coverage Matrix\n")
    text.append(coverage_matrix_markdown(desired))
    text.append("\n## Missing Desired Matrix Cells\n")
    text.append(md_table(missing, ["training_range", "input_cohort", "cohort_priority", "routing", "algorithm", "range_priority", "status"]))
    text.append("\n## Decision Reminder\n")
    text.append("- Primary range is `15-35`; `5-35` is the broad side-study range. `5-40` is not part of the active comparison table.\n")
    text.append("- Primary metric for production is WP80 fake rate in the target region, with AUC as supporting evidence.\n")
    text.append("- Stackers use score inputs and require extra runtime/ABCD/closure review. Isolation-visible rows are diagnostic only.\n")
    out.write_text("\n".join(text) + "\n")


def main() -> int:
    args = parse_args()
    rows = canonicalize_rows(discover_rows(args.root))
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "model_inventory.csv", rows)
    write_canonical_index(args.outdir / "canonical_model_index.csv", rows)
    desired, missing = desired_matrix(rows)
    write_csv(args.outdir / "desired_model_matrix.csv", desired)
    write_csv(args.outdir / "missing_model_matrix.csv", missing)
    write_csv(args.outdir / "cohort_coverage_summary.csv", cohort_summary(desired))
    write_markdown(args.outdir / "model_inventory.md", rows, desired, missing)
    write_training_backlog(args.outdir / "training_backlog.md", desired)
    (args.outdir / "canonical_coverage_matrix.md").write_text(
        "# AuAu Canonical Coverage Matrix\n\n"
        + coverage_matrix_markdown(desired)
        + "\n",
    )
    print(f"[inventory] rows={len(rows)}")
    print(f"[inventory] wrote {args.outdir / 'model_inventory.csv'}")
    print(f"[inventory] wrote {args.outdir / 'canonical_model_index.csv'}")
    print(f"[inventory] wrote {args.outdir / 'cohort_coverage_summary.csv'}")
    print(f"[inventory] wrote {args.outdir / 'model_inventory.md'}")
    print(f"[inventory] wrote {args.outdir / 'training_backlog.md'}")
    print(f"[inventory] missing_cells={len(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
