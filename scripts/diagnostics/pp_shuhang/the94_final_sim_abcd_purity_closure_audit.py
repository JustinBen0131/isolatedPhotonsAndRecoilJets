#!/usr/bin/env python3
"""Forensic audit of the accepted 20-lane PPG12 inclusive-SIM contract.

The analysis keeps two deliberately separate roles:

* the timestamped PPG12 component ROOT files, used only as the historical
  lane-output reference (not as the new same-event executable oracle);
* the accepted RecoilJets lane roots and their exact additive final merge.

The script records content, Sumw2, effective population, effective weight,
and lane-substitution influence for the unsuffixed PPG12 A/B/C/D estimator
inputs.  It is diagnostic-only: it never submits, merges, transfers, promotes,
or applies a fitted normalization.  In particular, the unresolved historical
Jet8 factor is reported but never used.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import shlex
import subprocess
import sys
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[3]
SCRIPT_PATH = Path(__file__).resolve()
FOCUSED_TEST_PATH = (
    SCRIPT_PATH.parent / "tests/test_the94_final_sim_abcd_purity_closure_audit.py"
)
PURITY_EVIDENCE_PRODUCER_PATH = (
    REPO
    / "scripts/diagnostics/pp_currentian/produce_ppg12_stitched_purity_evidence.py"
)
DEFAULT_SNAPSHOT = (
    REPO
    / "dataOutput/ppg12Parity/the94_final_sim_parity_closure"
    / "accepted_20_lane_forensics_hashbound_20260720"
    / "accepted_20_lane_snapshot.json"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the94_final_sim_parity_closure"
    / "accepted_20_lane_forensics"
)
DEFAULT_CURRENT_ROOT = (
    REPO
    / "InputFiles/the94_ppg12_inclusive_archivedreco_final_20260719_1615EDT"
    / "final_merged_roots"
    / "RecoilJets_jet8plus12plus20plus30plus40_archivedRecoDI_"
    "SIplusDI_periodCombined_MERGED.root"
)
DEFAULT_MERGE_AUDIT = (
    REPO
    / "InputFiles/the94_ppg12_inclusive_archivedreco_final_20260719_1615EDT"
    / "final_merged_roots/final_root_arithmetic_audit_summary.json"
)
DEFAULT_JET8_EVIDENCE = (
    REPO
    / "dataOutput/ppg12Parity/THE-94_ppg12_pp_inclusivejet_sim_parity"
    / "jet8_reference_provenance_20260704"
)

SAMPLES = ("jet8", "jet12", "jet20", "jet30", "jet40")
PERIODS = ("0mrad", "1p5mrad")
PPG12_PERIOD_NAMES = {"0mrad": "0rad", "1p5mrad": "1p5mrad"}
INTERACTIONS = ("si", "di")
ABCD_REGIONS = ("A", "B", "C", "D")
CLASS_REGIONS = tuple(
    f"{region}_{classification}"
    for region in ABCD_REGIONS
    for classification in ("signal", "notmatch")
)
REGIONS = ABCD_REGIONS + CLASS_REGIONS
REGION_OBJECTS = {
    "A": "h_tight_iso_cluster_0",
    "B": "h_tight_noniso_cluster_0",
    "C": "h_nontight_iso_cluster_0",
    "D": "h_nontight_noniso_cluster_0",
    "A_signal": "h_tight_iso_cluster_signal_0",
    "B_signal": "h_tight_noniso_cluster_signal_0",
    "C_signal": "h_nontight_iso_cluster_signal_0",
    "D_signal": "h_nontight_noniso_cluster_signal_0",
    "A_notmatch": "h_tight_iso_cluster_notmatch_0",
    "B_notmatch": "h_tight_noniso_cluster_notmatch_0",
    "C_notmatch": "h_nontight_iso_cluster_notmatch_0",
    "D_notmatch": "h_nontight_noniso_cluster_notmatch_0",
}
EXPECTED_LANES = tuple(
    (sample, period, interaction)
    for period in PERIODS
    for sample in SAMPLES
    for interaction in INTERACTIONS
)
LANE_NAME_RE = re.compile(
    r"(?P<period>0mrad|1p5mrad)_(?P<interaction>si|di)_"
    r"(?P<sample>jet8|jet12|jet20|jet30|jet40)\.root$"
)
REPORT_ET_MIN_GEV = 10.0
REPORT_ET_MAX_GEV = 36.0


def load_module(name: str, path: Path) -> Any:
    """Load one local diagnostic implementation without duplicating it."""
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PURITY_PRODUCER = load_module(
    "ppg12_stitched_purity_evidence_for_the94_audit",
    PURITY_EVIDENCE_PRODUCER_PATH,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(payload: Any) -> str:
    """Hash a JSON-compatible payload with one deterministic encoding."""
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def implementation_binding() -> dict[str, Any]:
    """Bind the audit, its focused regression test, and reused estimator code."""
    paths = {
        "inclusive_forensics_audit": SCRIPT_PATH,
        "focused_regression_test": FOCUSED_TEST_PATH,
        "canonical_stitched_purity_producer": PURITY_EVIDENCE_PRODUCER_PATH,
        "canonical_ppg12_estimator_source": Path(
            PURITY_PRODUCER.PPG12_ESTIMATOR_SOURCE
        ).resolve(),
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"implementation binding is incomplete: {missing}")
    links = {
        role: {"path": str(path), "sha256": sha256(path)}
        for role, path in sorted(paths.items())
    }
    return {
        "schema": "ppg12-inclusive-forensics-implementation-binding/v1",
        "files": links,
        "file_set_sha256": canonical_sha256(links),
    }


def report_bin_mask(edges: np.ndarray) -> np.ndarray:
    """Return the complete 10--36 GeV PPG12 reporting range.

    The last published bin is 32--36 GeV.  Using a 35 GeV upper boundary
    silently dropped that entire bin because the source TH1 edge is 36 GeV.
    """
    edges = np.asarray(edges, dtype=float)
    if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0.0):
        raise ValueError("reported-purity bin edges are malformed")
    return (edges[:-1] >= REPORT_ET_MIN_GEV) & (
        edges[1:] <= REPORT_ET_MAX_GEV
    )


def validate_exact_stitched_purity_payload(
    path: Path, expected_edges: np.ndarray
) -> dict[str, Any]:
    """Validate one canonical seed-42/20k-toy output and expose its components.

    Numerical estimator code remains owned by
    ``produce_ppg12_stitched_purity_evidence.py``.  This audit only validates
    that producer's hash-bound output and converts the already-recorded A--D
    and cB--cD diagnostics into reader-friendly component rows.
    """
    path = path.resolve()
    payload = json.loads(path.read_text())
    observed_edges, _ = PURITY_PRODUCER._load_candidate_corrected_purity(path)
    edges = np.asarray(observed_edges, dtype=float)
    expected_edges = np.asarray(expected_edges, dtype=float)
    if not np.array_equal(edges, expected_edges):
        raise ValueError("stitched-purity evidence binning differs from inclusive snapshot")

    run_diagnostics = payload.get("run_diagnostics")
    if not isinstance(run_diagnostics, Mapping):
        raise ValueError("stitched-purity evidence lacks deterministic run diagnostics")
    first = run_diagnostics.get("first")
    repeated = run_diagnostics.get("repeated")
    if not isinstance(first, list) or not isinstance(repeated, list):
        raise ValueError("stitched-purity run diagnostics must contain two lists")
    first_sha = PURITY_PRODUCER._payload_sha256(first)
    repeated_sha = PURITY_PRODUCER._payload_sha256(repeated)
    if (
        first_sha != repeated_sha
        or first != repeated
        or run_diagnostics.get("diagnostics_sha256") != first_sha
        or run_diagnostics.get("repeated_diagnostics_sha256") != repeated_sha
    ):
        raise ValueError("seed-42/20000-toy estimator diagnostics are not deterministic")
    if len(first) != edges.size - 1:
        raise ValueError("stitched-purity diagnostics do not cover every bin")

    component_rows: list[dict[str, Any]] = []
    required = (
        "A",
        "B",
        "C",
        "D",
        "cB",
        "cC",
        "cD",
        "cB_error",
        "cC_error",
        "cD_error",
        "A_effective",
        "B_effective",
        "C_effective",
        "D_effective",
    )
    for index, row in enumerate(first, start=1):
        if not isinstance(row, Mapping) or any(
            not isinstance(row.get(name), (int, float))
            or isinstance(row.get(name), bool)
            or not math.isfinite(float(row[name]))
            for name in required
        ):
            raise ValueError(f"malformed exact estimator diagnostics in bin {index}")
        if (
            int(row.get("bin", -1)) != index
            or float(row.get("x_low", math.nan)) != float(edges[index - 1])
            or float(row.get("x_high", math.nan)) != float(edges[index])
        ):
            raise ValueError(f"exact estimator diagnostic bin identity drift in bin {index}")
        for component in ("A", "B", "C", "D"):
            effective = float(row[f"{component}_effective"])
            value = float(row[component])
            error = abs(value) / math.sqrt(effective) if effective > 0.0 else math.nan
            component_rows.append(
                {
                    "source_family": "inclusive",
                    "component": component,
                    "bin": index,
                    "x_low_gev": float(edges[index - 1]),
                    "x_high_gev": float(edges[index]),
                    "value": value,
                    "error": error,
                    "semantics": "weighted_unsuffixed_ABCD_input",
                }
            )
        for component in ("cB", "cC", "cD"):
            component_rows.append(
                {
                    "source_family": "photon",
                    "component": component,
                    "bin": index,
                    "x_low_gev": float(edges[index - 1]),
                    "x_high_gev": float(edges[index]),
                    "value": float(row[component]),
                    "error": float(row[f"{component}_error"]),
                    "semantics": "photon_signal_region_over_A_signal",
                }
            )

    algorithm = payload.get("algorithm")
    if not isinstance(algorithm, Mapping):
        raise ValueError("stitched-purity evidence lacks algorithm provenance")
    expected_algorithm_links = {
        "producer": PURITY_EVIDENCE_PRODUCER_PATH,
        "source": Path(PURITY_PRODUCER.PPG12_ESTIMATOR_SOURCE).resolve(),
    }
    algorithm_links: dict[str, dict[str, str]] = {}
    for role, expected_path in expected_algorithm_links.items():
        link = algorithm.get(role)
        if not isinstance(link, Mapping):
            raise ValueError(f"stitched-purity algorithm lacks {role} binding")
        observed_path = Path(str(link.get("path", ""))).resolve()
        observed_sha = str(link.get("sha256", ""))
        if observed_path != expected_path or observed_sha != sha256(expected_path):
            raise ValueError(f"stitched-purity {role} implementation binding is stale")
        algorithm_links[role] = {
            "path": str(observed_path),
            "sha256": observed_sha,
        }

    return {
        "status": "PASS",
        "schema": payload.get("schema"),
        "path": str(path),
        "sha256": sha256(path),
        "random_seed": payload.get("random_seed"),
        "toy_count": payload.get("toy_count"),
        "purity": payload["purity"],
        "diagnostics_sha256": first_sha,
        "component_rows": component_rows,
        "algorithm_binding": algorithm_links,
    }


def validate_and_recompute_stitched_purity_attachment(
    path: Path,
    expected_edges: np.ndarray,
    records: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Bind a full 32-lane estimator artifact to these accepted 20 lanes.

    This deliberately calls the canonical producer's lane validator,
    aggregation routine, and fixed-seed estimator.  The historical component
    audit therefore cannot accidentally grow a second implementation of the
    PPG12 quadratic/toy/Gaussian-fit semantics.
    """
    path = path.resolve()
    validated = validate_exact_stitched_purity_payload(path, expected_edges)
    payload = json.loads(path.read_text())
    inputs = payload.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ValueError("stitched-purity evidence lacks input bindings")
    contract_link = inputs.get("contract")
    if not isinstance(contract_link, Mapping):
        raise ValueError("stitched-purity evidence lacks closure-contract binding")
    contract_path = Path(str(contract_link.get("path", ""))).resolve()
    contract = PURITY_PRODUCER.assembler._read_contract(contract_path)
    contract_sha = PURITY_PRODUCER.assembler._payload_sha256(contract)
    if contract_link.get("sha256") != contract_sha:
        raise ValueError("stitched-purity closure-contract hash is stale")
    raw_links = inputs.get("lanes")
    if not isinstance(raw_links, list):
        raise ValueError("stitched-purity evidence lane links are missing")
    lanes, normalized_links, bin_edges = (
        PURITY_PRODUCER.assembler._load_and_validate_lanes(
            raw_links, contract, require_merge_input=True
        )
    )
    PURITY_PRODUCER.assembler._validate_stitched_coverage(lanes, bin_edges)
    if inputs.get("lane_set_sha256") != PURITY_PRODUCER._payload_sha256(
        normalized_links
    ):
        raise ValueError("stitched-purity lane-set hash is stale")
    if not np.array_equal(np.asarray(bin_edges, dtype=float), expected_edges):
        raise ValueError("stitched-purity lane binning differs from inclusive snapshot")

    family_counts = Counter(str(lane["family"]) for lane in lanes)
    if family_counts != Counter({"inclusive": 20, "photon": 12}):
        raise ValueError(f"stitched-purity family coverage is not 20+12: {family_counts}")
    for lane in lanes:
        if lane["family"] != "inclusive":
            continue
        record = records[
            lane_label(
                "current", lane["sample"], lane["period"], lane["interaction"]
            )
        ]
        merge_input = lane.get("merge_input")
        if not isinstance(merge_input, Mapping):
            raise ValueError(f"{lane['lane_id']} lacks exact merge-input binding")
        if merge_input.get("sha256") != record.get("sha256"):
            raise ValueError(
                f"{lane['lane_id']} does not bind the accepted inclusive ROOT"
            )

    estimator_inputs = PURITY_PRODUCER._aggregate_estimator_inputs(lanes, bin_edges)
    current_sum = sum_records(
        [
            records[lane_label("current", sample, period, interaction)]
            for sample, period, interaction in EXPECTED_LANES
        ]
    )
    for region in ("A", "B", "C", "D", "A_signal", "A_notmatch"):
        observed_sumw, observed_sumw2 = estimator_inputs["inclusive"][region]
        expected_sumw, expected_sumw2, _ = current_sum[region]
        relative = float(contract["tolerances"]["weighted_cell_relative"])
        absolute = float(contract["tolerances"]["weighted_cell_abs_floor"])
        if not np.allclose(
            np.asarray(observed_sumw), expected_sumw, rtol=relative, atol=absolute
        ):
            raise ValueError(f"32-lane estimator inclusive {region} Sumw drift")
        if not np.allclose(
            np.asarray(observed_sumw2), expected_sumw2, rtol=relative, atol=absolute
        ):
            raise ValueError(f"32-lane estimator inclusive {region} Sumw2 drift")

    recomputed = PURITY_PRODUCER._build_repeated_purity_payload(
        estimator_inputs, contract
    )
    for key in (
        "schema",
        "random_seed",
        "toy_count",
        "purity",
        "fixed_seed_repetition",
        "run_diagnostics",
    ):
        if PURITY_PRODUCER._canonical_bytes(recomputed.get(key)) != (
            PURITY_PRODUCER._canonical_bytes(payload.get(key))
        ):
            raise ValueError(f"stitched-purity evidence does not reproduce key {key}")

    validated.update(
        {
            "verification": "recomputed_twice_with_live_canonical_producer",
            "contract": {"path": str(contract_path), "sha256": contract_sha},
            "lane_count": len(lanes),
            "family_counts": dict(sorted(family_counts.items())),
            "lanes": normalized_links,
            "lane_set_sha256": PURITY_PRODUCER._payload_sha256(normalized_links),
            "inclusive_lane_root_hash_binding": "exact_20_of_20",
        }
    )
    return validated


def snapshot_record_binding(
    records: Mapping[str, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Return the exact remote-file identity frozen by a lane snapshot."""
    binding: list[dict[str, Any]] = []
    for label in sorted(records):
        record = records[label]
        file_sha256 = record.get("sha256")
        if not isinstance(file_sha256, str) or re.fullmatch(r"[0-9a-f]{64}", file_sha256) is None:
            raise ValueError(f"{label} lacks a valid remote ROOT sha256")
        size = record.get("bytes")
        mtime_ns = record.get("mtime_ns")
        if not isinstance(size, int) or isinstance(size, bool) or size <= 0:
            raise ValueError(f"{label} lacks a positive remote ROOT byte count")
        if not isinstance(mtime_ns, int) or isinstance(mtime_ns, bool) or mtime_ns <= 0:
            raise ValueError(f"{label} lacks an exact remote ROOT mtime_ns")
        path = record.get("path")
        if not isinstance(path, str) or not path:
            raise ValueError(f"{label} lacks a remote ROOT path")
        binding.append(
            {
                "label": label,
                "path": path,
                "bytes": size,
                "mtime_ns": mtime_ns,
                "sha256": file_sha256,
            }
        )
    return binding


def validate_snapshot_record_binding(
    payload: Mapping[str, Any],
    records: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Fail closed when a collected remote-file identity is absent or stale."""
    binding = snapshot_record_binding(records)
    observed = canonical_sha256(binding)
    expected = payload.get("record_set_sha256")
    if not isinstance(expected, str) or re.fullmatch(r"[0-9a-f]{64}", expected) is None:
        raise ValueError("snapshot lacks a valid record_set_sha256")
    if observed != expected:
        raise ValueError(
            f"snapshot record-set hash mismatch: expected {expected}, observed {observed}"
        )
    return {
        "status": "PASS",
        "record_count": len(binding),
        "record_set_sha256": observed,
        "records": binding,
    }


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def lane_label(prefix: str, sample: str, period: str, interaction: str) -> str:
    return f"{prefix}:{sample}:{period}:{interaction}"


def expected_lane_labels(prefix: str) -> set[str]:
    return {
        lane_label(prefix, sample, period, interaction)
        for sample, period, interaction in EXPECTED_LANES
    }


def arrays(record: Mapping[str, Any], region: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    try:
        histogram = record["histograms"][region]
    except KeyError as exc:
        raise ValueError(f"{record.get('label', '<unknown>')} has no {region} histogram") from exc
    if histogram.get("missing"):
        raise ValueError(
            f"{record.get('label', '<unknown>')} is missing {histogram.get('object', region)}"
        )
    if histogram.get("sumw2_present") is False:
        raise ValueError(
            f"{record.get('label', '<unknown>')} {region} has no stored Sumw2"
        )
    values = np.asarray(histogram["values"], dtype=float)
    variances = np.asarray(histogram["variances"], dtype=float)
    edges = np.asarray(histogram["edges"], dtype=float)
    if values.ndim != 1 or variances.shape != values.shape or edges.size != values.size + 1:
        raise ValueError(f"malformed {record.get('label')} {region} histogram")
    if np.any(~np.isfinite(values)) or np.any(~np.isfinite(variances)) or np.any(~np.isfinite(edges)):
        raise ValueError(f"non-finite {record.get('label')} {region} histogram")
    if np.any(variances < 0.0) or np.any(np.diff(edges) <= 0.0):
        raise ValueError(f"invalid variance or binning in {record.get('label')} {region}")
    return values, variances, edges


def flow_arrays(record: Mapping[str, Any], region: str) -> tuple[np.ndarray, np.ndarray]:
    """Return underflow and overflow content/Sumw2 for exact merge checks."""
    histogram = record["histograms"][region]
    if "flow_values" not in histogram or "flow_variances" not in histogram:
        raise ValueError(f"{record.get('label')} {region} lacks flow-bin evidence")
    values = np.asarray(histogram["flow_values"], dtype=float)
    variances = np.asarray(histogram["flow_variances"], dtype=float)
    if values.shape != (2,) or variances.shape != (2,):
        raise ValueError(f"malformed flow-bin evidence in {record.get('label')} {region}")
    if np.any(~np.isfinite(values)) or np.any(~np.isfinite(variances)):
        raise ValueError(f"non-finite flow-bin evidence in {record.get('label')} {region}")
    return values, variances


def effective_population(sumw: float, sumw2: float) -> float:
    if sumw == 0.0 and sumw2 == 0.0:
        return 0.0
    if sumw2 <= 0.0:
        return math.nan
    return sumw * sumw / sumw2


def effective_weight(sumw: float, sumw2: float) -> float:
    """Return the exact factor satisfying sumw = neff * effective_weight."""
    if sumw == 0.0 and sumw2 == 0.0:
        return 0.0
    if sumw == 0.0:
        return math.nan
    return sumw2 / sumw


def ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator != 0.0 else math.nan


def finite_abs(value: float) -> float:
    return abs(value) if math.isfinite(value) else math.inf


def ratio_metrics(
    actual: np.ndarray,
    reference: np.ndarray,
    expected: np.ndarray | None = None,
) -> dict[str, Any]:
    if actual.shape != reference.shape:
        raise ValueError("ratio arrays have different shapes")
    if expected is None:
        expected = np.ones(actual.shape, dtype=bool)
    expected = np.asarray(expected, dtype=bool)
    if expected.shape != actual.shape:
        raise ValueError("ratio expected-bin mask has wrong shape")
    actual = actual[expected]
    reference = reference[expected]
    finite = np.isfinite(actual) & np.isfinite(reference)
    both_zero_mask = finite & (actual == 0.0) & (reference == 0.0)
    candidate_only_zero_mask = finite & (actual == 0.0) & (reference != 0.0)
    reference_only_zero_mask = finite & (actual != 0.0) & (reference == 0.0)
    compared = finite & (actual != 0.0) & (reference != 0.0)
    ratios = actual[compared] / reference[compared]
    residuals = ratios - 1.0
    nonfinite = int(np.count_nonzero(~finite))
    candidate_only_zero = int(np.count_nonzero(candidate_only_zero_mask))
    reference_only_zero = int(np.count_nonzero(reference_only_zero_mask))
    one_sided = candidate_only_zero + reference_only_zero
    coverage_pass = bool(nonfinite == 0 and one_sided == 0)
    numeric_rms = float(np.sqrt(np.mean(residuals**2))) if ratios.size else 0.0
    numeric_max = float(np.max(np.abs(residuals))) if ratios.size else 0.0
    return {
        "expected_bins": int(np.count_nonzero(expected)),
        "compared_bins": int(ratios.size),
        "both_zero_bins": int(np.count_nonzero(both_zero_mask)),
        "candidate_only_zero_bins": candidate_only_zero,
        "reference_only_zero_bins": reference_only_zero,
        "one_sided_zero_bins": one_sided,
        "nonfinite_bins": nonfinite,
        "coverage_pass": coverage_pass,
        "ratio_mean": float(np.mean(ratios)) if ratios.size else math.nan,
        "ratio_median": float(np.median(ratios)) if ratios.size else math.nan,
        "ratio_rms_about_unity": numeric_rms if coverage_pass else math.inf,
        "ratio_max_abs_from_unity": numeric_max if coverage_pass else math.inf,
    }


def assert_same_binning(records: Iterable[Mapping[str, Any]]) -> np.ndarray:
    reference_edges: np.ndarray | None = None
    for record in records:
        for region in REGIONS:
            _, _, edges = arrays(record, region)
            if reference_edges is None:
                reference_edges = edges
            elif not np.array_equal(reference_edges, edges):
                raise ValueError(f"binning mismatch in {record['label']} {region}")
    if reference_edges is None:
        raise ValueError("no histogram records")
    return reference_edges


def records_by_label(raw_records: Any) -> dict[str, Mapping[str, Any]]:
    """Validate exact record identity before constructing a label map."""
    if not isinstance(raw_records, list):
        raise ValueError("snapshot records must be a list")
    required = (
        expected_lane_labels("ppg12")
        | expected_lane_labels("current")
        | {"ppg12:final", "current:final"}
    )
    if len(raw_records) != len(required):
        raise ValueError(
            f"snapshot must contain exactly {len(required)} records, got {len(raw_records)}"
        )
    labels = [record.get("label") if isinstance(record, Mapping) else None for record in raw_records]
    if any(not isinstance(label, str) for label in labels):
        raise ValueError("every snapshot record must have a string label")
    duplicates = sorted(label for label, count in Counter(labels).items() if count != 1)
    if duplicates:
        raise ValueError(f"duplicate record labels: {duplicates}")
    observed = set(labels)
    if observed != required:
        raise ValueError(
            f"snapshot record identity mismatch: missing={sorted(required-observed)}, "
            f"unexpected={sorted(observed-required)}"
        )
    return {str(record["label"]): record for record in raw_records}


def validate_record_set(records: Mapping[str, Mapping[str, Any]]) -> None:
    required = (
        expected_lane_labels("ppg12")
        | expected_lane_labels("current")
        | {"ppg12:final", "current:final"}
    )
    if set(records) != required:
        raise ValueError(
            f"record identity mismatch: missing={sorted(required-set(records))}, "
            f"unexpected={sorted(set(records)-required)}"
        )
    failures = [
        label
        for label in sorted(required)
        if records[label].get("error") or records[label].get("exists") is False
    ]
    if failures:
        raise ValueError(f"unreadable required records: {failures}")
    assert_same_binning(records[label] for label in sorted(required))


def sum_records(
    lane_records: Sequence[Mapping[str, Any]],
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    result: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for region in REGIONS:
        first_values, first_variances, edges = arrays(lane_records[0], region)
        storage_classes = {
            str(record["histograms"][region].get("storage_class", ""))
            for record in lane_records
        }
        storage_dtype = np.float32 if storage_classes == {"TH1F"} else np.float64
        values = np.zeros(first_values.shape, dtype=storage_dtype)
        variances = np.zeros_like(first_variances)
        for record in lane_records:
            lane_values, lane_variances, lane_edges = arrays(record, region)
            if not np.array_equal(edges, lane_edges):
                raise ValueError(f"binning mismatch while summing {record['label']} {region}")
            values += lane_values.astype(storage_dtype)
            variances += lane_variances
        result[region] = values.astype(float), variances, edges
    return result


def validate_merge_audit(path: Path, current_root: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    failures: list[str] = []
    if payload.get("status") != "PASS":
        failures.append("merge_audit_status_not_pass")
    if payload.get("inputs_fixed_order_count") != 20:
        failures.append("merge_audit_input_count_not_20")
    if payload.get("max_content_delta") != 0.0:
        failures.append("merge_audit_content_not_exact")
    if payload.get("max_sumw2_delta") != 0.0:
        failures.append("merge_audit_sumw2_not_exact")
    if payload.get("failures"):
        failures.append("merge_audit_contains_failures")
    observed_sha = sha256(current_root)
    if payload.get("output_sha256") != observed_sha:
        failures.append("merge_audit_root_sha_mismatch")
    return {
        "path": str(path),
        "root": str(current_root),
        "root_sha256": observed_sha,
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "evidence": payload,
    }


def validate_snapshot_manifest(
    payload: Mapping[str, Any],
    records: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    manifest = payload.get("merge_manifest")
    failures: list[str] = []
    if not isinstance(manifest, Mapping):
        return {"status": "FAIL", "failures": ["missing_embedded_merge_manifest"]}
    expected_paths = [
        str(records[lane_label("current", sample, period, interaction)].get("path", ""))
        for sample, period, interaction in EXPECTED_LANES
    ]
    manifest_paths = manifest.get("lane_roots_fixed_order")
    if manifest_paths != expected_paths:
        failures.append("manifest_lane_order_or_identity_mismatch")
    if len(expected_paths) != 20 or len(set(expected_paths)) != 20 or any(not p for p in expected_paths):
        failures.append("snapshot_current_lane_paths_not_exact_unique_20")
    if manifest.get("lane_count") != 20:
        failures.append("manifest_lane_count_not_20")
    if str(manifest.get("jet5", "")).lower() != "excluded":
        failures.append("manifest_jet5_not_excluded")
    merge_contract = str(manifest.get("merge_contract", "")).lower()
    unit_scale_proven = "no external scale" in merge_contract
    if not unit_scale_proven:
        failures.append("manifest_does_not_prove_no_external_scale")
    final_path = str(records["current:final"].get("path", ""))
    if str(manifest.get("final_root", "")) != final_path:
        failures.append("manifest_final_root_path_mismatch")
    final_sha256 = str(records["current:final"].get("sha256", ""))
    if str(manifest.get("final_root_sha256", "")) != final_sha256:
        failures.append("manifest_final_root_sha256_mismatch")
    lane_manifest_sha = payload.get("lane_manifest_sha256")
    if not lane_manifest_sha or lane_manifest_sha != manifest.get("lane_roots_manifest_sha256"):
        failures.append("lane_manifest_sha_mismatch")
    missing_config = [
        label
        for label in sorted(expected_lane_labels("current"))
        if not records[label].get("config", {}).get("sha256")
    ]
    if missing_config:
        failures.append("current_lane_config_digest_missing")
    config_groups: dict[str, list[str]] = {}
    for label in sorted(expected_lane_labels("current")):
        digest = records[label].get("config", {}).get("sha256")
        if digest:
            config_groups.setdefault(str(digest), []).append(label)
    return {
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "unit_scale_proven": unit_scale_proven,
        "lane_paths": expected_paths,
        "lane_manifest_sha256": lane_manifest_sha,
        "config_digest_groups": config_groups,
        "missing_config_labels": missing_config,
        "evidence": manifest,
    }


def merge_closure_rows(
    records: Mapping[str, Mapping[str, Any]],
    prefix: str,
    *,
    absolute_tolerance: float = 1.0e-9,
    relative_tolerance: float = 1.0e-12,
) -> tuple[list[dict[str, Any]], bool]:
    lane_records = [
        records[lane_label(prefix, sample, period, interaction)]
        for sample, period, interaction in EXPECTED_LANES
    ]
    summed = sum_records(lane_records)
    final = records[f"{prefix}:final"]
    rows: list[dict[str, Any]] = []
    passed = True
    for region in REGIONS:
        sum_values, sum_variances, edges = summed[region]
        final_values, final_variances, final_edges = arrays(final, region)
        if not np.array_equal(edges, final_edges):
            raise ValueError(f"{prefix} final binning mismatch in {region}")
        for index, (sumw, sumw2, target, target2) in enumerate(
            zip(sum_values, sum_variances, final_values, final_variances), start=1
        ):
            content_tolerance = max(absolute_tolerance, relative_tolerance * abs(target))
            sumw2_tolerance = max(absolute_tolerance, relative_tolerance * abs(target2))
            content_delta = float(sumw - target)
            sumw2_delta = float(sumw2 - target2)
            bin_pass = bool(
                abs(content_delta) <= content_tolerance
                and abs(sumw2_delta) <= sumw2_tolerance
            )
            passed = bool(passed and bin_pass)
            rows.append(
                {
                    "source": prefix,
                    "region": region,
                    "bin": index,
                    "x_low_gev": float(edges[index - 1]),
                    "x_high_gev": float(edges[index]),
                    "lane_sumw": float(sumw),
                    "final_sumw": float(target),
                    "content_delta": content_delta,
                    "content_tolerance": content_tolerance,
                    "lane_sumw2": float(sumw2),
                    "final_sumw2": float(target2),
                    "sumw2_delta": sumw2_delta,
                    "sumw2_tolerance": sumw2_tolerance,
                    "pass": bin_pass,
                }
            )
        storage_classes = {
            str(record["histograms"][region].get("storage_class", ""))
            for record in lane_records
        }
        flow_dtype = np.float32 if storage_classes == {"TH1F"} else np.float64
        lane_flow_values = np.zeros(2, dtype=flow_dtype)
        lane_flow_variances = np.zeros(2, dtype=float)
        for record in lane_records:
            flow_values, flow_variances = flow_arrays(record, region)
            lane_flow_values += flow_values.astype(flow_dtype)
            lane_flow_variances += flow_variances
        final_flow_values, final_flow_variances = flow_arrays(final, region)
        for flow_index, flow_name in enumerate(("underflow", "overflow")):
            target = float(final_flow_values[flow_index])
            target2 = float(final_flow_variances[flow_index])
            content_delta = float(lane_flow_values[flow_index] - target)
            sumw2_delta = float(lane_flow_variances[flow_index] - target2)
            content_tolerance = max(absolute_tolerance, relative_tolerance * abs(target))
            sumw2_tolerance = max(absolute_tolerance, relative_tolerance * abs(target2))
            bin_pass = bool(
                abs(content_delta) <= content_tolerance
                and abs(sumw2_delta) <= sumw2_tolerance
            )
            passed = bool(passed and bin_pass)
            rows.append(
                {
                    "source": prefix,
                    "region": region,
                    "bin": flow_name,
                    "x_low_gev": math.nan,
                    "x_high_gev": math.nan,
                    "lane_sumw": float(lane_flow_values[flow_index]),
                    "final_sumw": target,
                    "content_delta": content_delta,
                    "content_tolerance": content_tolerance,
                    "lane_sumw2": float(lane_flow_variances[flow_index]),
                    "final_sumw2": target2,
                    "sumw2_delta": sumw2_delta,
                    "sumw2_tolerance": sumw2_tolerance,
                    "pass": bin_pass,
                }
            )
        lane_entries = sum(
            float(record["histograms"][region].get("entries", math.nan))
            for record in lane_records
        )
        final_entries = float(final["histograms"][region].get("entries", math.nan))
        entries_tolerance = max(
            absolute_tolerance, relative_tolerance * abs(final_entries)
        )
        entries_pass = bool(
            math.isfinite(lane_entries)
            and math.isfinite(final_entries)
            and abs(lane_entries - final_entries) <= entries_tolerance
        )
        passed = bool(passed and entries_pass)
        rows.append(
            {
                "source": prefix,
                "region": region,
                "bin": "global_entries",
                "x_low_gev": math.nan,
                "x_high_gev": math.nan,
                "lane_sumw": lane_entries,
                "final_sumw": final_entries,
                "content_delta": lane_entries - final_entries,
                "content_tolerance": entries_tolerance,
                "lane_sumw2": math.nan,
                "final_sumw2": math.nan,
                "sumw2_delta": math.nan,
                "sumw2_tolerance": math.nan,
                "pass": entries_pass,
            }
        )
    return rows, passed


def lane_bin_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for sample, period, interaction in EXPECTED_LANES:
        reference = records[lane_label("ppg12", sample, period, interaction)]
        current = records[lane_label("current", sample, period, interaction)]
        ref_abcd_total = sum(arrays(reference, item)[0] for item in ABCD_REGIONS)
        cur_abcd_total = sum(arrays(current, item)[0] for item in ABCD_REGIONS)
        ref_signal_total = sum(arrays(reference, f"{item}_signal")[0] for item in ABCD_REGIONS)
        cur_signal_total = sum(arrays(current, f"{item}_signal")[0] for item in ABCD_REGIONS)
        ref_notmatch_total = sum(arrays(reference, f"{item}_notmatch")[0] for item in ABCD_REGIONS)
        cur_notmatch_total = sum(arrays(current, f"{item}_notmatch")[0] for item in ABCD_REGIONS)
        for region in REGIONS:
            reference_values, reference_variances, edges = arrays(reference, region)
            current_values, current_variances, current_edges = arrays(current, region)
            if not np.array_equal(edges, current_edges):
                raise ValueError(f"binning mismatch: {sample} {period} {interaction} {region}")
            for index, (ref, ref2, cur, cur2) in enumerate(
                zip(reference_values, reference_variances, current_values, current_variances),
                start=1,
            ):
                ref_neff = effective_population(float(ref), float(ref2))
                cur_neff = effective_population(float(cur), float(cur2))
                ref_weff = effective_weight(float(ref), float(ref2))
                cur_weff = effective_weight(float(cur), float(cur2))
                content_ratio = ratio(float(cur), float(ref))
                population_ratio = ratio(cur_neff, ref_neff)
                weight_ratio = ratio(cur_weff, ref_weff)
                reconstructed = (
                    population_ratio * weight_ratio
                    if math.isfinite(population_ratio) and math.isfinite(weight_ratio)
                    else math.nan
                )
                if region in ABCD_REGIONS:
                    ref_fraction = ratio(float(ref), float(ref_abcd_total[index - 1]))
                    cur_fraction = ratio(float(cur), float(cur_abcd_total[index - 1]))
                    ref_signal = arrays(reference, f"{region}_signal")[0][index - 1]
                    cur_signal = arrays(current, f"{region}_signal")[0][index - 1]
                    ref_notmatch = arrays(reference, f"{region}_notmatch")[0][index - 1]
                    cur_notmatch = arrays(current, f"{region}_notmatch")[0][index - 1]
                    ref_signal_fraction = ratio(float(ref_signal), float(ref))
                    cur_signal_fraction = ratio(float(cur_signal), float(cur))
                    ref_notmatch_fraction = ratio(float(ref_notmatch), float(ref))
                    cur_notmatch_fraction = ratio(float(cur_notmatch), float(cur))
                    ref_unclassified_fraction = ratio(
                        float(ref - ref_signal - ref_notmatch), float(ref)
                    )
                    cur_unclassified_fraction = ratio(
                        float(cur - cur_signal - cur_notmatch), float(cur)
                    )
                    ref_unclassified = float(ref - ref_signal - ref_notmatch)
                    cur_unclassified = float(cur - cur_signal - cur_notmatch)
                else:
                    classification = "signal" if region.endswith("_signal") else "notmatch"
                    ref_total = ref_signal_total if classification == "signal" else ref_notmatch_total
                    cur_total = cur_signal_total if classification == "signal" else cur_notmatch_total
                    ref_fraction = ratio(float(ref), float(ref_total[index - 1]))
                    cur_fraction = ratio(float(cur), float(cur_total[index - 1]))
                    ref_signal_fraction = math.nan
                    cur_signal_fraction = math.nan
                    ref_notmatch_fraction = math.nan
                    cur_notmatch_fraction = math.nan
                    ref_unclassified_fraction = math.nan
                    cur_unclassified_fraction = math.nan
                    ref_signal = math.nan
                    cur_signal = math.nan
                    ref_notmatch = math.nan
                    cur_notmatch = math.nan
                    ref_unclassified = math.nan
                    cur_unclassified = math.nan
                rows.append(
                    {
                        "sample": sample,
                        "period": period,
                        "interaction": interaction,
                        "region": region,
                        "bin": index,
                        "x_low_gev": float(edges[index - 1]),
                        "x_high_gev": float(edges[index]),
                        "ppg12_sumw": float(ref),
                        "current_sumw": float(cur),
                        "current_over_ppg12_sumw": content_ratio,
                        "ppg12_sumw2": float(ref2),
                        "current_sumw2": float(cur2),
                        "current_over_ppg12_sumw2": ratio(float(cur2), float(ref2)),
                        "ppg12_effective_population": ref_neff,
                        "current_effective_population": cur_neff,
                        "effective_population_ratio": population_ratio,
                        "ppg12_effective_weight": ref_weff,
                        "current_effective_weight": cur_weff,
                        "effective_weight_ratio": weight_ratio,
                        "population_times_weight_ratio": reconstructed,
                        "factorization_delta": (
                            reconstructed - content_ratio
                            if math.isfinite(reconstructed) and math.isfinite(content_ratio)
                            else math.nan
                        ),
                        "ppg12_histogram_entries": reference["histograms"][region].get("entries"),
                        "current_histogram_entries": current["histograms"][region].get("entries"),
                        "histogram_entries_semantics": "global_TH1_fill_population_not_per_bin",
                        "ppg12_weighted_region_fraction": ref_fraction,
                        "current_weighted_region_fraction": cur_fraction,
                        "ppg12_truth_signal_sumw": float(ref_signal),
                        "current_truth_signal_sumw": float(cur_signal),
                        "ppg12_truth_signal_fraction": ref_signal_fraction,
                        "current_truth_signal_fraction": cur_signal_fraction,
                        "ppg12_truth_notmatch_sumw": float(ref_notmatch),
                        "current_truth_notmatch_sumw": float(cur_notmatch),
                        "ppg12_truth_notmatch_fraction": ref_notmatch_fraction,
                        "current_truth_notmatch_fraction": cur_notmatch_fraction,
                        "ppg12_truth_unclassified_sumw": ref_unclassified,
                        "current_truth_unclassified_sumw": cur_unclassified,
                        "ppg12_truth_unclassified_fraction": ref_unclassified_fraction,
                        "current_truth_unclassified_fraction": cur_unclassified_fraction,
                        "per_bin_raw_fill_population": "not_stored_by_weighted_TH1",
                        "per_bin_population_proxy": "effective_population=sumw^2/sumw2",
                    }
                )
    return rows


def lane_summary_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for sample, period, interaction in EXPECTED_LANES:
        reference = records[lane_label("ppg12", sample, period, interaction)]
        current = records[lane_label("current", sample, period, interaction)]
        for region in REGIONS:
            ref_values, ref_variances, _ = arrays(reference, region)
            cur_values, cur_variances, _ = arrays(current, region)
            content = ratio_metrics(cur_values, ref_values)
            variance = ratio_metrics(cur_variances, ref_variances)
            rows.append(
                {
                    "sample": sample,
                    "period": period,
                    "interaction": interaction,
                    "region": region,
                    "ppg12_entries": reference["histograms"][region].get("entries"),
                    "current_entries": current["histograms"][region].get("entries"),
                    "entries_ratio": ratio(
                        float(current["histograms"][region].get("entries", math.nan)),
                        float(reference["histograms"][region].get("entries", math.nan)),
                    ),
                    **{f"content_{key}": value for key, value in content.items()},
                    **{f"sumw2_{key}": value for key, value in variance.items()},
                }
            )
    return rows


def central_estimators(
    histograms: Mapping[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> dict[str, np.ndarray]:
    a = histograms["A"][0]
    b = histograms["B"][0]
    c = histograms["C"][0]
    d = histograms["D"][0]
    signal = histograms["A_signal"][0]
    notmatch = histograms["A_notmatch"][0]
    with np.errstate(divide="ignore", invalid="ignore"):
        truth = np.divide(
            signal,
            a,
            out=np.full_like(signal, np.nan),
            where=a != 0.0,
        )
        classified_fraction = np.divide(
            signal + notmatch,
            a,
            out=np.full_like(signal, np.nan),
            where=a != 0.0,
        )
        unclassified_fraction = np.divide(
            a - signal - notmatch,
            a,
            out=np.full_like(signal, np.nan),
            where=a != 0.0,
        )
        raw = 1.0 - np.divide(
            b * c,
            a * d,
            out=np.full_like(a, np.nan),
            where=(a != 0.0) & (d != 0.0),
        )
    return {
        "truth": truth,
        "raw_central_algebraic": raw,
        "classified_fraction": classified_fraction,
        "unclassified_fraction": unclassified_fraction,
    }


def hybrid_sum(
    baseline: Mapping[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    current_record: Mapping[str, Any],
    reference_record: Mapping[str, Any],
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    hybrid: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for region in REGIONS:
        base_values, base_variances, edges = baseline[region]
        current_values, current_variances, current_edges = arrays(current_record, region)
        reference_values, reference_variances, reference_edges = arrays(reference_record, region)
        if not np.array_equal(edges, current_edges) or not np.array_equal(edges, reference_edges):
            raise ValueError(f"binning mismatch in hybrid {current_record['label']} {region}")
        hybrid[region] = (
            base_values - current_values + reference_values,
            base_variances - current_variances + reference_variances,
            edges,
        )
    return hybrid


def substitution_rows(
    records: Mapping[str, Mapping[str, Any]],
    *,
    by_sample: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    current_lanes = {
        (sample, period, interaction): records[lane_label("current", sample, period, interaction)]
        for sample, period, interaction in EXPECTED_LANES
    }
    ppg12_lanes = {
        (sample, period, interaction): records[lane_label("ppg12", sample, period, interaction)]
        for sample, period, interaction in EXPECTED_LANES
    }
    current_sum = sum_records(list(current_lanes.values()))
    ppg12_sum = sum_records(list(ppg12_lanes.values()))
    baseline_estimators = central_estimators(current_sum)
    oracle_estimators = central_estimators(ppg12_sum)
    edges = current_sum["A"][2]
    groups: list[tuple[str, tuple[tuple[str, str, str], ...]]] = []
    if by_sample:
        groups = [
            (
                sample,
                tuple(lane for lane in EXPECTED_LANES if lane[0] == sample),
            )
            for sample in SAMPLES
        ]
    else:
        groups = [
            (f"{sample}:{period}:{interaction}", ((sample, period, interaction),))
            for sample, period, interaction in EXPECTED_LANES
        ]

    bin_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for group_name, lanes in groups:
        hybrid = current_sum
        for lane in lanes:
            hybrid = hybrid_sum(hybrid, current_lanes[lane], ppg12_lanes[lane])
        hybrid_estimators = central_estimators(hybrid)
        group_improvements: list[float] = []
        for observable in ("truth", "raw_central_algebraic"):
            baseline = baseline_estimators[observable]
            oracle = oracle_estimators[observable]
            candidate = hybrid_estimators[observable]
            for index, (before, after, target) in enumerate(
                zip(baseline, candidate, oracle), start=1
            ):
                before_distance = finite_abs(float(before - target))
                after_distance = finite_abs(float(after - target))
                improvement = before_distance - after_distance
                if math.isfinite(improvement):
                    group_improvements.append(improvement)
                bin_rows.append(
                    {
                        "group_kind": "sample" if by_sample else "lane",
                        "group": group_name,
                        "lanes_replaced": len(lanes),
                        "observable": observable,
                        "bin": index,
                        "x_low_gev": float(edges[index - 1]),
                        "x_high_gev": float(edges[index]),
                        "baseline_current": float(before),
                        "hybrid_current_with_ppg12_group": float(after),
                        "ppg12_historical_component": float(target),
                        "baseline_abs_distance": before_distance,
                        "hybrid_abs_distance": after_distance,
                        "absolute_distance_improvement": improvement,
                        "fraction_of_abs_distance_removed": (
                            improvement / before_distance
                            if math.isfinite(before_distance) and before_distance > 0.0
                            else math.nan
                        ),
                    }
                )
        finite_improvements = np.asarray(
            [value for value in group_improvements if math.isfinite(value)], dtype=float
        )
        summary_rows.append(
            {
                "group_kind": "sample" if by_sample else "lane",
                "group": group_name,
                "lanes_replaced": len(lanes),
                "bins_compared": int(finite_improvements.size),
                "sum_abs_distance_improvement": (
                    float(np.sum(finite_improvements)) if finite_improvements.size else math.nan
                ),
                "mean_abs_distance_improvement": (
                    float(np.mean(finite_improvements)) if finite_improvements.size else math.nan
                ),
                "improved_bin_count": int(np.count_nonzero(finite_improvements > 0.0)),
                "worsened_bin_count": int(np.count_nonzero(finite_improvements < 0.0)),
            }
        )
    summary_rows.sort(
        key=lambda row: row["sum_abs_distance_improvement"], reverse=True
    )
    return bin_rows, summary_rows


def top_lane_attribution_rows(
    influence_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Select the most closure-improving one-lane substitution per estimator bin."""
    groups: dict[tuple[str, int], list[Mapping[str, Any]]] = {}
    for row in influence_rows:
        if row.get("group_kind") != "lane":
            raise ValueError("top-lane attribution requires lane-level influence rows")
        key = (str(row["observable"]), int(row["bin"]))
        groups.setdefault(key, []).append(row)
    if len(groups) == 0 or any(len(rows) != len(EXPECTED_LANES) for rows in groups.values()):
        raise ValueError("lane attribution does not contain every lane in every estimator bin")
    output: list[dict[str, Any]] = []
    for key in sorted(groups):
        rows = groups[key]
        best = max(rows, key=lambda row: float(row["absolute_distance_improvement"]))
        improvement = float(best["absolute_distance_improvement"])
        output.append(
            {
                "observable": key[0],
                "bin": key[1],
                "x_low_gev": best["x_low_gev"],
                "x_high_gev": best["x_high_gev"],
                "baseline_current": best["baseline_current"],
                "ppg12_historical_component": best["ppg12_historical_component"],
                "top_lane": best["group"],
                "hybrid_with_top_lane_replaced": best["hybrid_current_with_ppg12_group"],
                "absolute_distance_improvement": improvement,
                "fraction_of_abs_distance_removed": best["fraction_of_abs_distance_removed"],
                "attribution_status": (
                    "one_lane_substitution_improves_closure"
                    if improvement > 0.0
                    else "no_single_lane_substitution_improves_closure"
                ),
                "tested_lane_count": len(rows),
                "interpretation": "one-at-a-time substitution diagnostic; not an independent causal decomposition",
            }
        )
    return output


def component_substitution_rows(
    records: Mapping[str, Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Attribute A--D/truth-cell residuals with one-lane substitutions.

    These rows expose *which input cell* a lane changes before the nonlinear
    purity solver is evaluated.  The result is deliberately described as a
    one-at-a-time diagnostic rather than an additive causal decomposition.
    """
    components = ("A", "B", "C", "D", "A_signal", "A_notmatch")
    current_lanes = {
        lane: records[lane_label("current", *lane)] for lane in EXPECTED_LANES
    }
    reference_lanes = {
        lane: records[lane_label("ppg12", *lane)] for lane in EXPECTED_LANES
    }
    current_sum = sum_records(list(current_lanes.values()))
    reference_sum = sum_records(list(reference_lanes.values()))
    edges = current_sum["A"][2]
    in_report = report_bin_mask(edges)
    rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for lane in EXPECTED_LANES:
        label = ":".join(lane)
        lane_improvements: list[float] = []
        for component in components:
            baseline = current_sum[component][0]
            target = reference_sum[component][0]
            current_component = arrays(current_lanes[lane], component)[0]
            reference_component = arrays(reference_lanes[lane], component)[0]
            hybrid = baseline - current_component + reference_component
            for index, (before, after, oracle) in enumerate(
                zip(baseline, hybrid, target), start=1
            ):
                baseline_abs = abs(float(before - oracle))
                hybrid_abs = abs(float(after - oracle))
                if oracle != 0.0:
                    baseline_fractional = baseline_abs / abs(float(oracle))
                    hybrid_fractional = hybrid_abs / abs(float(oracle))
                    improvement = baseline_fractional - hybrid_fractional
                else:
                    baseline_fractional = math.nan
                    hybrid_fractional = math.nan
                    improvement = math.nan
                if in_report[index - 1] and math.isfinite(improvement):
                    lane_improvements.append(improvement)
                rows.append(
                    {
                        "group_kind": "lane",
                        "group": label,
                        "source_family": "inclusive",
                        "component": component,
                        "bin": index,
                        "x_low_gev": float(edges[index - 1]),
                        "x_high_gev": float(edges[index]),
                        "inside_reported_10_to_36_gev": bool(in_report[index - 1]),
                        "baseline_current": float(before),
                        "hybrid_current_with_ppg12_lane": float(after),
                        "ppg12_historical_component": float(oracle),
                        "lane_current_minus_ppg12": float(
                            current_component[index - 1]
                            - reference_component[index - 1]
                        ),
                        "baseline_abs_distance": baseline_abs,
                        "hybrid_abs_distance": hybrid_abs,
                        "baseline_fractional_distance": baseline_fractional,
                        "hybrid_fractional_distance": hybrid_fractional,
                        "fractional_distance_improvement": improvement,
                        "interpretation": (
                            "one-at-a-time historical component substitution; "
                            "not an independent causal decomposition"
                        ),
                    }
                )
        finite = np.asarray(lane_improvements, dtype=float)
        summaries.append(
            {
                "group_kind": "lane",
                "group": label,
                "source_family": "inclusive",
                "reported_component_bins_compared": int(finite.size),
                "sum_fractional_distance_improvement": (
                    float(np.sum(finite)) if finite.size else math.nan
                ),
                "mean_fractional_distance_improvement": (
                    float(np.mean(finite)) if finite.size else math.nan
                ),
                "improved_component_bin_count": int(np.count_nonzero(finite > 0.0)),
                "worsened_component_bin_count": int(np.count_nonzero(finite < 0.0)),
            }
        )
    summaries.sort(
        key=lambda row: float(row["sum_fractional_distance_improvement"]),
        reverse=True,
    )

    grouped: dict[tuple[str, int], list[dict[str, Any]]] = {}
    for row in rows:
        if row["inside_reported_10_to_36_gev"]:
            grouped.setdefault((str(row["component"]), int(row["bin"])), []).append(row)
    top_rows: list[dict[str, Any]] = []
    for (component, bin_index), candidates in sorted(grouped.items()):
        finite_candidates = [
            row
            for row in candidates
            if math.isfinite(float(row["fractional_distance_improvement"]))
        ]
        if not finite_candidates:
            continue
        best = max(
            finite_candidates,
            key=lambda row: float(row["fractional_distance_improvement"]),
        )
        top_rows.append(
            {
                "source_family": "inclusive",
                "component": component,
                "bin": bin_index,
                "x_low_gev": best["x_low_gev"],
                "x_high_gev": best["x_high_gev"],
                "top_lane": best["group"],
                "baseline_current": best["baseline_current"],
                "ppg12_historical_component": best["ppg12_historical_component"],
                "hybrid_with_top_lane_replaced": best[
                    "hybrid_current_with_ppg12_lane"
                ],
                "fractional_distance_improvement": best[
                    "fractional_distance_improvement"
                ],
                "tested_lane_count": len(candidates),
                "attribution_status": (
                    "one_lane_substitution_improves_component_closure"
                    if float(best["fractional_distance_improvement"]) > 0.0
                    else "no_single_lane_substitution_improves_component_closure"
                ),
                "interpretation": (
                    "one-at-a-time historical component substitution; not an "
                    "independent causal decomposition"
                ),
            }
        )
    return rows, summaries, top_rows


def audit_jet8_provenance(evidence_dir: Path) -> dict[str, Any]:
    manifest_path = evidence_dir / "ppg12_jet8_reference_provenance_manifest.json"
    component_path = evidence_dir / "jet8_weight_factor_comparison.csv"
    search_path = evidence_dir / "ppg12_jet8_real_ratio_search.csv"
    missing = [str(path) for path in (manifest_path, component_path, search_path) if not path.exists()]
    if missing:
        return {
            "status": "evidence_missing",
            "missing": missing,
            "external_scale_applied": None,
        }
    manifest = json.loads(manifest_path.read_text())
    component_rows = read_csv(component_path)
    component_factors = [float(row["current_over_ppg12_integral"]) for row in component_rows]
    entry_ratios = [float(row["current_over_ppg12_raw_entries"]) for row in component_rows]
    per_entry_ratios = [float(row["current_over_ppg12_integral_per_entry"]) for row in component_rows]
    search_rows = read_csv(search_path)
    recovered_documented_explanation = [
        row
        for row in search_rows
        if str(row.get("explains_2p338", "")).lower() == "true"
        and "observed" not in str(row.get("quantity", "")).lower()
        and "forced_effective" not in str(row.get("quantity", "")).lower()
    ]
    return {
        "status": "unresolved_weight_provenance",
        "observed_factor": float(manifest["observed_factor"]),
        "component_factor_min": min(component_factors),
        "component_factor_max": max(component_factors),
        "raw_entry_ratio_min": min(entry_ratios),
        "raw_entry_ratio_max": max(entry_ratios),
        "integral_per_entry_ratio_min": min(per_entry_ratios),
        "integral_per_entry_ratio_max": max(per_entry_ratios),
        "known_old_xsec_factor": float(manifest["known_old_xsec_factor"]),
        "forced_effective_xsec_pb": float(manifest["forced_effective_xsec_pb"]),
        "forced_effective_xsec_status": "numerical_only_not_found_in_ppg12_source",
        "ratio_search_row_count": len(search_rows),
        "documented_explanation_rows": recovered_documented_explanation,
        "documented_explanation_recovered": bool(recovered_documented_explanation),
        "conclusion": (
            "The approximately 2.338 factor is component-independent while raw "
            "entry counts agree at approximately 0.2 percent.  It is therefore a "
            "weight/exposure provenance difference, not a missing candidate "
            "population.  The known older Jet8 cross section is too small to "
            "explain it."
        ),
        "nominal_policy": "no_external_or_fitted_jet8_scale",
        "external_scale_applied": None,
        "evidence_sha256": {
            path.name: sha256(path)
            for path in (manifest_path, component_path, search_path)
        },
    }


def parity_summary(records: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    current_sum = sum_records(
        [
            records[lane_label("current", sample, period, interaction)]
            for sample, period, interaction in EXPECTED_LANES
        ]
    )
    ppg12_sum = sum_records(
        [
            records[lane_label("ppg12", sample, period, interaction)]
            for sample, period, interaction in EXPECTED_LANES
        ]
    )
    edges = current_sum["A"][2]
    report_mask = report_bin_mask(edges)
    if not np.any(report_mask):
        raise ValueError("no histogram bins lie inside the reported purity range")
    region_content_metrics = {
        region: ratio_metrics(current_sum[region][0], ppg12_sum[region][0], report_mask)
        for region in REGIONS
    }
    region_sumw2_metrics = {
        region: ratio_metrics(current_sum[region][1], ppg12_sum[region][1], report_mask)
        for region in REGIONS
    }

    def population_vector(values: np.ndarray, variances: np.ndarray) -> np.ndarray:
        output = np.full_like(values, np.nan)
        both_zero = (values == 0.0) & (variances == 0.0)
        output[both_zero] = 0.0
        valid = variances > 0.0
        output[valid] = values[valid] ** 2 / variances[valid]
        return output

    region_effective_population_metrics = {
        region: ratio_metrics(
            population_vector(current_sum[region][0], current_sum[region][1]),
            population_vector(ppg12_sum[region][0], ppg12_sum[region][1]),
            report_mask,
        )
        for region in REGIONS
    }
    current_estimators = central_estimators(current_sum)
    ppg12_estimators = central_estimators(ppg12_sum)
    estimator_metrics = {
        observable: ratio_metrics(
            current_estimators[observable], ppg12_estimators[observable], report_mask
        )
        for observable in ("truth", "raw_central_algebraic")
    }
    all_metrics = (
        list(region_content_metrics.values())
        + list(region_sumw2_metrics.values())
        + list(region_effective_population_metrics.values())
        + list(estimator_metrics.values())
    )
    coverage_pass = bool(all(metric["coverage_pass"] for metric in all_metrics))
    max_abs = max(
        float(metrics["ratio_max_abs_from_unity"])
        for metrics in all_metrics
    ) if all_metrics else math.inf
    rms = max(
        float(metrics["ratio_rms_about_unity"])
        for metrics in all_metrics
    ) if all_metrics else math.inf
    content_max_abs = max(
        float(metrics["ratio_max_abs_from_unity"])
        for metrics in (*region_content_metrics.values(), *estimator_metrics.values())
    )
    content_rms = max(
        float(metrics["ratio_rms_about_unity"])
        for metrics in (*region_content_metrics.values(), *estimator_metrics.values())
    )
    return {
        "reported_et_range_gev": [REPORT_ET_MIN_GEV, REPORT_ET_MAX_GEV],
        "reported_bin_count": int(np.count_nonzero(report_mask)),
        "coverage_pass": coverage_pass,
        "region_content_metrics": region_content_metrics,
        "region_sumw2_metrics": region_sumw2_metrics,
        "region_effective_population_metrics": region_effective_population_metrics,
        "central_estimator_metrics": estimator_metrics,
        "maximum_content_ratio_deviation": content_max_abs,
        "maximum_content_ratio_rms": content_rms,
        "maximum_ratio_deviation": max_abs,
        "maximum_ratio_rms": rms,
    }


def analyze(
    snapshot: Path,
    outdir: Path,
    current_root: Path,
    merge_audit: Path,
    jet8_evidence: Path,
    *,
    max_ratio_deviation: float,
    max_ratio_rms: float,
    stitched_purity_evidence: Path | None = None,
) -> dict[str, Any]:
    snapshot = snapshot.resolve()
    outdir = outdir.resolve()
    current_root = current_root.resolve()
    merge_audit = merge_audit.resolve()
    jet8_evidence = jet8_evidence.resolve()
    stitched_purity_evidence = (
        stitched_purity_evidence.resolve()
        if stitched_purity_evidence is not None
        else None
    )
    payload = json.loads(snapshot.read_text())
    records = records_by_label(payload.get("records"))
    validate_record_set(records)
    record_binding = validate_snapshot_record_binding(payload, records)
    merge_validation = validate_merge_audit(merge_audit, current_root)
    snapshot_manifest_validation = validate_snapshot_manifest(payload, records)
    current_closure, current_closure_pass = merge_closure_rows(records, "current")
    ppg12_closure, ppg12_closure_pass = merge_closure_rows(records, "ppg12")
    all_closure = current_closure + ppg12_closure
    bins = lane_bin_rows(records)
    summaries = lane_summary_rows(records)
    lane_influence, lane_influence_summary = substitution_rows(records)
    top_lane_attribution = top_lane_attribution_rows(lane_influence)
    sample_influence, sample_influence_summary = substitution_rows(records, by_sample=True)
    (
        component_influence,
        component_influence_summary,
        component_top_attribution,
    ) = component_substitution_rows(records)
    parity = parity_summary(records)
    parity_pass = bool(
        parity["coverage_pass"]
        and parity["maximum_ratio_deviation"] <= max_ratio_deviation
        and parity["maximum_ratio_rms"] <= max_ratio_rms
    )
    jet8 = audit_jet8_provenance(jet8_evidence)
    current_sum = sum_records(
        [
            records[lane_label("current", sample, period, interaction)]
            for sample, period, interaction in EXPECTED_LANES
        ]
    )
    if stitched_purity_evidence is None:
        stitched_purity: dict[str, Any] = {
            "status": "NOT_SUPPLIED",
            "role": "optional_hash_bound_32_lane_exact_estimator_attachment",
            "required_for_admission": True,
        }
        exact_component_rows: list[dict[str, Any]] = []
    else:
        stitched_purity = validate_and_recompute_stitched_purity_attachment(
            stitched_purity_evidence, current_sum["A"][2], records
        )
        exact_component_rows = list(stitched_purity.pop("component_rows"))

    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "lane_bin_contract": outdir / "lane_bin_contract.csv",
        "lane_region_summary": outdir / "lane_region_summary.csv",
        "additive_closure": outdir / "additive_closure.csv",
        "lane_substitution_influence": outdir / "lane_substitution_influence.csv",
        "lane_substitution_summary": outdir / "lane_substitution_summary.csv",
        "lane_substitution_top_attribution": outdir / "lane_substitution_top_attribution.csv",
        "lane_component_substitution_influence": outdir
        / "lane_component_substitution_influence.csv",
        "lane_component_substitution_summary": outdir
        / "lane_component_substitution_summary.csv",
        "lane_component_substitution_top_attribution": outdir
        / "lane_component_substitution_top_attribution.csv",
        "exact_estimator_components": outdir / "exact_estimator_components.csv",
        "sample_substitution_influence": outdir / "sample_substitution_influence.csv",
        "sample_substitution_summary": outdir / "sample_substitution_summary.csv",
        "jet8_provenance": outdir / "jet8_provenance_summary.json",
        "report": outdir / "forensic_report.json",
    }
    write_csv(outputs["lane_bin_contract"], bins)
    write_csv(outputs["lane_region_summary"], summaries)
    write_csv(outputs["additive_closure"], all_closure)
    write_csv(outputs["lane_substitution_influence"], lane_influence)
    write_csv(outputs["lane_substitution_summary"], lane_influence_summary)
    write_csv(outputs["lane_substitution_top_attribution"], top_lane_attribution)
    write_csv(outputs["lane_component_substitution_influence"], component_influence)
    write_csv(
        outputs["lane_component_substitution_summary"], component_influence_summary
    )
    write_csv(
        outputs["lane_component_substitution_top_attribution"],
        component_top_attribution,
    )
    write_csv(outputs["exact_estimator_components"], exact_component_rows)
    write_csv(outputs["sample_substitution_influence"], sample_influence)
    write_csv(outputs["sample_substitution_summary"], sample_influence_summary)
    outputs["jet8_provenance"].write_text(json.dumps(jet8, indent=2, sort_keys=True) + "\n")

    structural_pass = bool(
        merge_validation["status"] == "PASS"
        and snapshot_manifest_validation["status"] == "PASS"
        and current_closure_pass
        and ppg12_closure_pass
        and jet8.get("status") == "unresolved_weight_provenance"
    )
    report = {
        "schema_version": 4,
        "status": (
            "parity_closed"
            if structural_pass and parity_pass
            else "analyzed_nonclosing"
            if structural_pass
            else "structural_failure"
        ),
        "structural_pass": structural_pass,
        "parity_pass": parity_pass,
        "artifact_classification": "historical_diagnostic_only",
        "admission_eligible": False,
        "admission_ineligibility_reasons": [
            "historical_component_reference_is_not_the_same_event_executable_oracle",
            "accepted_weighted_TH1_inputs_lack_exact_per_bin_raw_fill_evidence",
        ],
        "lane_contract": {
            "lane_count": len(EXPECTED_LANES),
            "samples": SAMPLES,
            "periods": PERIODS,
            "interactions": INTERACTIONS,
            "jet5": "excluded",
            "abcd_population": "unsuffixed",
            "external_or_fitted_scale": (
                "forbidden_and_not_present"
                if snapshot_manifest_validation.get("unit_scale_proven")
                else "provenance_incomplete"
            ),
            "raw_fill_population_per_bin": (
                "unavailable_from_weighted_TH1; requires an unweighted companion "
                "or candidate-level oracle ledger"
            ),
            "reference_role": "historical_component_archive_not_same_event_executable_oracle",
        },
        "parity_thresholds": {
            "maximum_ratio_deviation": max_ratio_deviation,
            "maximum_ratio_rms": max_ratio_rms,
        },
        "parity": parity,
        "merge_validation": merge_validation,
        "snapshot_manifest_validation": snapshot_manifest_validation,
        "snapshot_sha256": sha256(snapshot),
        "snapshot_record_binding": record_binding,
        "implementation_binding": implementation_binding(),
        "exact_stitched_purity_attachment": stitched_purity,
        "current_additive_closure_pass": current_closure_pass,
        "ppg12_additive_closure_pass": ppg12_closure_pass,
        "jet8_provenance": jet8,
        "most_influential_lanes": lane_influence_summary[:5],
        "most_influential_samples": sample_influence_summary[:5],
        "most_influential_component_lanes": component_influence_summary[:5],
        "interpretation_boundary": (
            "Lane substitution compares the accepted THE-94 lanes to timestamped "
            "PPG12 component ROOTs; those files are a historical-output reference, "
            "not the same-event executable oracle. It uses truth A_signal/A and the central "
            "algebraic zero-leakage quantity 1-BC/(AD) plus separate A--D cell "
            "substitutions to attribute the inclusive residual. When supplied, the "
            "32-lane attachment binds the photon+jet cB/cC/cD inputs and recomputes "
            "the exact seed-42, 20,000-toy/Gaussian-fit PPG12 estimator through the "
            "canonical producer. Neither attachment can promote these historical "
            "TH1 decompositions to admission evidence because exact per-bin raw fill "
            "counts and a same-event executable oracle are absent."
        ),
        "outputs": {name: str(path) for name, path in outputs.items() if name != "report"},
    }
    outputs["report"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    artifact_manifest = {
        "schema": "ppg12-inclusive-lane-forensics-artifact-manifest/v2",
        "status": report["status"],
        "input_binding": {
            "snapshot": {"path": str(snapshot), "sha256": sha256(snapshot)},
            "record_set_sha256": record_binding["record_set_sha256"],
            "current_root": {"path": str(current_root), "sha256": sha256(current_root)},
            "merge_audit": {"path": str(merge_audit), "sha256": sha256(merge_audit)},
            "stitched_purity_evidence": (
                {
                    "path": str(stitched_purity_evidence),
                    "sha256": sha256(stitched_purity_evidence),
                }
                if stitched_purity_evidence is not None
                else None
            ),
        },
        "implementation_binding": report["implementation_binding"],
        "population_semantics": {
            "per_bin_effective_population": "sumw^2/sumw2",
            "global_histogram_entries": "TH1 fEntries; not a per-bin fill count",
            "per_bin_raw_fill_population": "unavailable from the accepted weighted TH1 inputs",
            "admission_use": "historical_diagnostic_only_even_when_exact_estimator_attachment_is_present",
        },
        "artifacts": {},
    }
    for name, path in sorted(outputs.items()):
        if name == "artifact_manifest":
            continue
        artifact_manifest["artifacts"][name] = {
            "path": str(path),
            "sha256": sha256(path),
        }
    artifact_manifest_path = outdir / "artifact_manifest.json"
    artifact_manifest_path.write_text(
        json.dumps(artifact_manifest, indent=2, sort_keys=True) + "\n"
    )
    report["outputs"]["artifact_manifest"] = str(artifact_manifest_path)
    return report


def remote_collect(
    lane_manifest: str,
    current_final: str,
    merge_manifest: str,
    ppg12_dir: str,
    output: Path,
) -> dict[str, Any]:
    """Collect compact read-only ROOT evidence over nested SSH."""
    payload = f'''import hashlib, json, os, re
import numpy as np
import uproot

lane_manifest = {lane_manifest!r}
current_final = {current_final!r}
merge_manifest_path = {merge_manifest!r}
ppg12_dir = {ppg12_dir!r}
samples = {json.dumps(SAMPLES)}
periods = {json.dumps(PERIODS)}
ppg12_period_names = {json.dumps(PPG12_PERIOD_NAMES)}
interactions = {json.dumps(INTERACTIONS)}
region_objects = {json.dumps(REGION_OBJECTS)}
lane_re = re.compile({LANE_NAME_RE.pattern!r})

def object_text(obj):
    if not obj:
        return ""
    try:
        return str(obj)
    except Exception:
        return ""

def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def config_digest(root_file, prefix):
    for name in (
        (prefix + "/analysis_config_yaml") if prefix else "analysis_config_yaml",
        (prefix + "/config") if prefix else "config",
        "analysis_config_yaml",
        "config",
    ):
        text = object_text(root_file[name]) if name in root_file else ""
        if text:
            return {{"object": name, "sha256": hashlib.sha256(text.encode()).hexdigest()}}
    return {{"object": None, "sha256": None}}

def read_file(label, path, prefix=""):
    row = {{"label": label, "path": path, "exists": os.path.isfile(path)}}
    if not row["exists"]:
        row["error"] = "missing"
        return row
    before = os.stat(path)
    try:
        root_file = uproot.open(path)
    except Exception as exc:
        row["error"] = "missing_or_unreadable:" + type(exc).__name__
        return row
    row["bytes"] = os.path.getsize(path)
    row["mtime"] = os.path.getmtime(path)
    row["config"] = config_digest(root_file, prefix)
    row["histograms"] = {{}}
    for region, name in region_objects.items():
        object_name = (prefix + "/" + name) if prefix else name
        if object_name not in root_file:
            row["histograms"][region] = {{"missing": True, "object": object_name}}
            continue
        histogram = root_file[object_name]
        sumw2 = np.asarray(histogram.member("fSumw2"), dtype=float)
        ncells = int(histogram.member("fNcells"))
        sumw2_present = sumw2.size == ncells
        if not sumw2_present:
            row["error"] = "missing_stored_sumw2:" + object_name
            row["histograms"][region] = {{
                "missing": False,
                "sumw2_present": False,
                "object": object_name,
            }}
            continue
        all_values = np.asarray(histogram.values(flow=True), dtype=float)
        values = all_values[1:-1].tolist()
        variances = sumw2[1:-1].tolist()
        flow_values = [float(all_values[0]), float(all_values[-1])]
        flow_variances = [float(sumw2[0]), float(sumw2[-1])]
        edges = np.asarray(histogram.axis().edges(), dtype=float).tolist()
        row["histograms"][region] = {{
            "object": object_name,
            "storage_class": histogram.classname,
            "values": values,
            "variances": variances,
            "flow_values": flow_values,
            "flow_variances": flow_variances,
            "sumw2_present": True,
            "edges": edges,
            "entries": float(histogram.member("fEntries")),
        }}
    root_file.close()
    row["sha256"] = file_sha256(path)
    after = os.stat(path)
    if before.st_size != after.st_size or before.st_mtime_ns != after.st_mtime_ns:
        raise RuntimeError("remote ROOT changed during collection: " + path)
    row["bytes"] = int(after.st_size)
    row["mtime_ns"] = int(after.st_mtime_ns)
    return row

lane_paths = [line.strip() for line in open(lane_manifest) if line.strip()]
if len(lane_paths) != 20 or len(set(lane_paths)) != 20:
    raise RuntimeError("lane manifest must contain exactly 20 unique paths")
current_paths = {{}}
for path in lane_paths:
    match = lane_re.search(path)
    if not match:
        raise RuntimeError("unparseable lane path: " + path)
    key = (match.group("sample"), match.group("period"), match.group("interaction"))
    if key in current_paths:
        raise RuntimeError("duplicate lane: " + repr(key))
    current_paths[key] = path

merge_payload = json.load(open(merge_manifest_path))
if merge_payload.get("lane_roots_fixed_order") != lane_paths:
    raise RuntimeError("merge manifest lane order/identity differs from lane manifest")
if merge_payload.get("lane_count") != 20:
    raise RuntimeError("merge manifest lane count is not 20")
if merge_payload.get("final_root") != current_final:
    raise RuntimeError("merge manifest final ROOT does not match requested final")
if str(merge_payload.get("jet5", "")).lower() != "excluded":
    raise RuntimeError("merge manifest does not exclude Jet5")
if "no external scale" not in str(merge_payload.get("merge_contract", "")).lower():
    raise RuntimeError("merge manifest does not prove no external scale")
lane_manifest_sha256 = hashlib.sha256(open(lane_manifest, "rb").read()).hexdigest()
if merge_payload.get("lane_roots_manifest_sha256") != lane_manifest_sha256:
    raise RuntimeError("merge manifest lane-manifest SHA mismatch")

records = []
for period in periods:
    ppg12_period = ppg12_period_names[period]
    for sample in samples:
        for interaction in interactions:
            component = "nom" if interaction == "si" else "double"
            records.append(read_file(
                "ppg12:%s:%s:%s" % (sample, period, interaction),
                "%s/MC_efficiency_%s_%s_bdt_nom_%s.root" % (
                    ppg12_dir, sample, component, ppg12_period
                ),
            ))
            records.append(read_file(
                "current:%s:%s:%s" % (sample, period, interaction),
                current_paths[(sample, period, interaction)],
                "SIM",
            ))
records.append(read_file("ppg12:final", ppg12_dir + "/MC_efficiency_jet_bdt_nom.root"))
records.append(read_file("current:final", current_final, "SIM"))
record_binding = [{{
    "label": row["label"],
    "path": row["path"],
    "bytes": row["bytes"],
    "mtime_ns": row["mtime_ns"],
    "sha256": row["sha256"],
}} for row in sorted(records, key=lambda item: item["label"])]
record_set_sha256 = hashlib.sha256(json.dumps(
    record_binding, sort_keys=True, separators=(",", ":"), allow_nan=False
).encode()).hexdigest()
print("JSON_RESULT_START")
print(json.dumps({{
    "records": records,
    "record_set_sha256": record_set_sha256,
    "merge_manifest": merge_payload,
    "lane_manifest": lane_manifest,
    "lane_manifest_sha256": lane_manifest_sha256,
}}, sort_keys=True))
'''
    encoded = base64.b64encode(payload.encode()).decode()
    remote_command = f"python3 -c 'import base64; exec(base64.b64decode(\"{encoded}\"))'"
    command = [
        "ssh",
        "-q",
        "-o",
        "BatchMode=yes",
        "-o",
        "ConnectTimeout=20",
        "-o",
        "ServerAliveInterval=30",
        "-o",
        "ServerAliveCountMax=4",
        "patsfan753@ssh.sdcc.bnl.gov",
        (
            "ssh -q -T -o BatchMode=yes -o ConnectTimeout=20 "
            "-o ServerAliveInterval=30 -o ServerAliveCountMax=4 -o LogLevel=ERROR "
            "-o StrictHostKeyChecking=no "
            "-o UserKnownHostsFile=/dev/null sphnxuser05.sdcc.bnl.gov "
            + shlex.quote(remote_command)
        ),
    ]
    env = os.environ.copy()
    # The desktop shell can retain a stale SSH_AUTH_SOCK after the macOS agent
    # rotates.  The SDCC policy requires the live launchd socket for the second
    # hop; prefer it whenever available and fail quickly in batch mode rather
    # than hanging at an interactive password prompt.
    launchd_socket = subprocess.check_output(
        ["launchctl", "getenv", "SSH_AUTH_SOCK"], text=True
    ).strip()
    if launchd_socket:
        env["SSH_AUTH_SOCK"] = launchd_socket
    try:
        result = subprocess.run(
            command,
            env=env,
            text=True,
            capture_output=True,
            check=False,
            timeout=900,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError("read-only remote collection exceeded 900 seconds") from exc
    if result.returncode != 0:
        raise RuntimeError(f"read-only remote collection failed: {result.stderr.strip()}")
    marker = "JSON_RESULT_START\n"
    marker_index = result.stdout.rfind(marker)
    if marker_index < 0:
        raise RuntimeError("read-only remote collection emitted no JSON marker")
    collected = json.loads(result.stdout[marker_index + len(marker) :])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(collected, indent=2, sort_keys=True) + "\n")
    return collected


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze_parser = subparsers.add_parser("analyze", help="analyze a compact lane snapshot")
    analyze_parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    analyze_parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    analyze_parser.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    analyze_parser.add_argument("--merge-audit", type=Path, default=DEFAULT_MERGE_AUDIT)
    analyze_parser.add_argument("--jet8-evidence", type=Path, default=DEFAULT_JET8_EVIDENCE)
    analyze_parser.add_argument(
        "--stitched-purity-evidence",
        type=Path,
        help=(
            "optional full 32-lane purity-repeat JSON; when supplied it is "
            "hash-validated and recomputed with the canonical seed-42/20000-toy producer"
        ),
    )
    analyze_parser.add_argument("--max-ratio-deviation", type=float, default=1.0e-4)
    analyze_parser.add_argument("--max-ratio-rms", type=float, default=2.0e-5)
    analyze_parser.add_argument(
        "--fail-on-parity",
        action="store_true",
        help="exit 2 after writing outputs when the diagnostic parity thresholds fail",
    )

    collect_parser = subparsers.add_parser(
        "collect-remote",
        help="collect compact read-only ROOT evidence from an exact lane manifest",
    )
    collect_parser.add_argument("--lane-manifest", required=True)
    collect_parser.add_argument("--current-final", required=True)
    collect_parser.add_argument("--merge-manifest", required=True)
    collect_parser.add_argument("--ppg12-dir", required=True)
    collect_parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "collect-remote":
        result = remote_collect(
            args.lane_manifest,
            args.current_final,
            args.merge_manifest,
            args.ppg12_dir,
            args.output,
        )
        print(json.dumps({"output": str(args.output), "records": len(result["records"])}, indent=2))
        return 0
    report = analyze(
        args.snapshot,
        args.outdir,
        args.current_root,
        args.merge_audit,
        args.jet8_evidence,
        max_ratio_deviation=args.max_ratio_deviation,
        max_ratio_rms=args.max_ratio_rms,
        stitched_purity_evidence=args.stitched_purity_evidence,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["structural_pass"] is False:
        return 1
    if args.fail_on_parity and report["parity_pass"] is False:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
