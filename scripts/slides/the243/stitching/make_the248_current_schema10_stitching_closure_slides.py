#!/usr/bin/env python3
"""Render THE-248 current-schema-10 pp and embedded-AuAu stitching slides.

This is a local, standalone consumer of the compact assembled stitching JSON.
It does not read ROOT files, submit work, or mutate Google Slides.  Plot images
are composed in memory; the only PNG outputs are the two requested full slides.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[4]
THE291_CONTRACT_ROLE = "canonical_source_stitch_diagnostic"
DEFAULT_OUTPUT_DIR = (
    REPO
    / "dataOutput/the243_golden_ppg_analysis_closure_deck_20260821"
    / "stitching_updates_current_20260820"
)
DEFAULT_INPUT = DEFAULT_OUTPUT_DIR / "THE248_SCHEMA10_STITCHING_ASSEMBLED_V4.json"
DEFAULT_ASSEMBLY_RECEIPT = (
    DEFAULT_OUTPUT_DIR / "THE248_SCHEMA10_STITCHING_ASSEMBLY_RECEIPT_V4.json"
)
SOURCE_PNGS = {
    "pp": DEFAULT_OUTPUT_DIR / "source_pp_stitching_slide.png",
    "auau": DEFAULT_OUTPUT_DIR / "source_auau_stitching_slide.png",
}
SOURCE_OBJECT_IDS = {
    "pp": "g3f78191bdc6_0_23",
    "auau": "g3f78191bdc6_0_28",
}
EXPECTED_SOURCE_PNG_HASHES = {
    "pp": "afa9a1b2f84526a84e1f7946ff4bf5c3cf045a1408f77b37ba20bf119a5f2b0e",
    "auau": "a39946a0b78364a5502915038224c951d8fc5a27596db7090588ccf8f04d5a1f",
}
OUTPUT_PNG_NAMES = {
    "pp": "the243_pp_current_schema10_stitch_weight_closure.png",
    "auau": "the243_auau_current_schema10_stitch_weight_closure.png",
}
SCRIPT_NAMES = {
    "pp": "the243_pp_current_schema10_stitch_weight_closure_script.md",
    "auau": "the243_auau_current_schema10_stitch_weight_closure_script.md",
}
MANIFEST_NAME = "THE248_STITCHING_SLIDE_MANIFEST.json"
ASSEMBLED_SCHEMA = "THE248Schema10StitchingAssemblyV4"
ASSEMBLED_STATUS = "PASS"
EXPECTED_ASSEMBLY_COUNTS = {
    "input_root_count": 20_006,
    "unique_input_root_count": 20_006,
    "bulk_product_count": 408,
    "foreground_product_count": 2,
    "sample_count": 14,
}
LAYOUT_NAMES = {
    "pp": "the243_pp_current_schema10_stitch_weight_closure.layout_nodes.json",
    "auau": "the243_auau_current_schema10_stitch_weight_closure.layout_nodes.json",
}

W, H = 2560, 1440
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
TIMES_BOLD_ITALIC = FONT_DIR / "Times New Roman Bold Italic.ttf"

INK = (18, 25, 39)
MUTED = (76, 87, 105)
WHITE = (255, 255, 255)
LINE = (202, 211, 225)
PP_ACCENT = (28, 111, 116)
PP_JET_ACCENT = (92, 78, 152)
AUAU_ACCENT = (26, 111, 116)
AUAU_JET_ACCENT = (92, 78, 152)
PHOTON_SOFT = (235, 248, 248)
JET_SOFT = (244, 241, 250)
NOTE_SOFT = (239, 245, 255)
OK_GREEN = (22, 101, 52)
WARN_AMBER = (146, 64, 14)

FIT_CONTRACT = {
    "name": "fixed_log_pt_cubic_diagnostic",
    "formula": "ln(y) = c0 + c1*ln(pT) + c2*ln(pT)^2 + c3*ln(pT)^3",
    "degree": 3,
    "fit_space": "unweighted least squares in ln(weighted density) versus ln(pT)",
    "point_scope": "all finite positive ownership-selected bins inside the fixed display range",
    "purpose": "diagnostic smooth reference only",
    "prohibitions": [
        "do not tune weights toward the fit",
        "do not tune ownership windows toward the fit",
        "do not use fit residuals or boundary readouts as validity gates",
    ],
}


class InputError(RuntimeError):
    """Raised when the assembled input cannot support a truthful render."""


@dataclass(frozen=True)
class Window:
    low: float
    high: float | None
    original: Any

    def contains(self, x: np.ndarray) -> np.ndarray:
        mask = x >= self.low
        if self.high is not None:
            mask &= x < self.high
        return mask

    def display(self) -> str:
        if self.high is None:
            return f"≥{self.low:g}"
        return f"{self.low:g}-{self.high:g}"


@dataclass(frozen=True)
class SampleData:
    sample_id: str
    canonical_id: str
    system: str
    family: str
    threshold: int
    display_name: str
    edges: np.ndarray
    density: np.ndarray
    sumw2: np.ndarray
    raw_counts: np.ndarray
    normalization_cross_section_pb: float
    normalization_denominator_events: float
    normalization_weight_pb_per_owned_event: float
    window: Window
    source_product_hashes: Any

    @property
    def centers(self) -> np.ndarray:
        return 0.5 * (self.edges[:-1] + self.edges[1:])

    @property
    def errors(self) -> np.ndarray:
        return np.sqrt(np.clip(self.sumw2, 0.0, None))


@dataclass(frozen=True)
class PanelConfig:
    system: str
    family: str
    title: str
    sample_label: str
    xlabel: str
    xlim: tuple[float, float]
    ratio_ylim: tuple[float, float]


PANEL_CONFIGS = {
    ("pp", "photon"): PanelConfig(
        "pp",
        "photon",
        "Photon+jet current exact-stitch contract",
        r"$p{+}p$ PYTHIA8 PhotonJet samples",
        r"max truth photon $p_T$ [GeV]",
        (5.0, 40.0),
        (0.55, 1.45),
    ),
    ("pp", "jet"): PanelConfig(
        "pp",
        "jet",
        "Inclusive-jet current exact-stitch contract",
        r"$p{+}p$ PYTHIA8 inclusive jets",
        r"max R=0.4 truth jet $p_T$ [GeV]",
        (8.0, 50.0),
        (0.55, 1.45),
    ),
    ("auau", "photon"): PanelConfig(
        "auau",
        "photon",
        "Embedded PhotonJet current ownership",
        "Embedded Au+Au Photon+Jet samples",
        r"max truth-filter photon $p_T$ [GeV]",
        (12.0, 40.0),
        (0.75, 1.25),
    ),
    ("auau", "jet"): PanelConfig(
        "auau",
        "jet",
        "Embedded inclusive-jet current ownership",
        "Embedded Au+Au inclusive jets",
        r"max R=0.4 truth jet $p_T$ [GeV]",
        (12.0, 50.0),
        (0.75, 1.25),
    ),
}

SAMPLE_STYLE = {
    "pp_photon5": ("#2ca25f", "o"),
    "pp_photon10": ("#2b8cbe", "o"),
    "pp_photon20": ("#e6550d", "o"),
    "pp_jet8": ("#d33682", "o"),
    "pp_jet12": ("#2ca02c", "o"),
    "pp_jet20": ("#1f77b4", "o"),
    "pp_jet30": ("#ff7f0e", "o"),
    "pp_jet40": ("#6f4eb2", "o"),
    "auau_photon12": ("#214cc3", "o"),
    "auau_photon20": ("#d2232a", "s"),
    "auau_jet12": ("#2b55b7", "o"),
    "auau_jet20": ("#ff7f0e", "s"),
    "auau_jet30": ("#cc2fa8", "^"),
    "auau_jet40": ("#238b45", "v"),
}

REQUIRED_CANONICAL_IDS = {
    "pp_photon5",
    "pp_photon10",
    "pp_photon20",
    "pp_jet8",
    "pp_jet12",
    "pp_jet20",
    "pp_jet30",
    "pp_jet40",
    "auau_photon12",
    "auau_photon20",
    "auau_jet12",
    "auau_jet20",
    "auau_jet30",
    "auau_jet40",
}
INTERNAL_TASK_ID_RE = re.compile(r"\bTHE(?:[-\s]?\d{1,4}[A-Z]?)\b", re.IGNORECASE)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def finite_float(value: Any, label: str, *, positive: bool = False) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise InputError(f"{label}: expected a number, got {value!r}") from exc
    if not math.isfinite(out):
        raise InputError(f"{label}: expected a finite number, got {out!r}")
    if positive and out <= 0.0:
        raise InputError(f"{label}: expected a positive number, got {out!r}")
    return out


def numeric_array(value: Any, label: str) -> np.ndarray:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise InputError(f"{label}: expected an array")
    try:
        out = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise InputError(f"{label}: array contains non-numeric values") from exc
    if out.ndim != 1 or not np.all(np.isfinite(out)):
        raise InputError(f"{label}: expected a one-dimensional finite array")
    return out


def parse_window(value: Any, label: str) -> Window:
    low: Any = None
    high: Any = None
    if isinstance(value, Mapping):
        for key in ("low_gev", "low", "lower_gev", "lower", "min_gev", "min"):
            if key in value:
                low = value[key]
                break
        for key in ("high_gev", "high", "upper_gev", "upper", "max_gev", "max"):
            if key in value:
                high = value[key]
                break
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == 2:
        low, high = value
    elif isinstance(value, str):
        nums = [float(item) for item in re.findall(r"[-+]?\d+(?:\.\d+)?", value)]
        if ">=" in value or "geq" in value.lower():
            if nums:
                low, high = nums[-1], None
        elif len(nums) >= 2:
            low, high = nums[-2], nums[-1]
    if low is None:
        raise InputError(f"{label}: cannot resolve ownership-window lower edge from {value!r}")
    low_f = finite_float(low, f"{label}.low")
    high_f: float | None
    if (
        high is None
        or (isinstance(high, str) and high.strip().lower() in {"inf", "infinity", "none", "null"})
        or (isinstance(high, (int, float)) and (not math.isfinite(float(high)) or float(high) < 0.0))
    ):
        high_f = None
    else:
        high_f = finite_float(high, f"{label}.high")
        if high_f <= low_f:
            raise InputError(f"{label}: upper edge {high_f} must exceed lower edge {low_f}")
    return Window(low=low_f, high=high_f, original=value)


def classify_sample(sample_id: str, record: Mapping[str, Any]) -> tuple[str, str, int, str, str]:
    context_fields = [
        sample_id,
        record.get("sample_id", ""),
        record.get("sample", ""),
        record.get("sample_name", ""),
        record.get("system", ""),
        record.get("collision_system", ""),
        record.get("family", ""),
        record.get("sample_family", ""),
    ]
    text = " ".join(str(item) for item in context_fields if item is not None).lower()
    compact = re.sub(r"[^a-z0-9]+", "", text)
    if any(token in compact for token in ("auau", "embedded", "embed")):
        system = "auau"
    elif re.search(r"(^|[^a-z])p\+?p([^a-z]|$)", text) or "pp" in compact:
        system = "pp"
    else:
        system = ""

    family = "photon" if "photon" in compact else ("jet" if "jet" in compact else "")
    match = re.search(r"photon(?:jet)?0*(5|10|12|20)(?!\d)", compact)
    if not match and family == "jet":
        matches = re.findall(r"jet0*(8|12|20|30|40)(?!\d)", compact)
        if matches:
            match = re.search(rf"({matches[-1]})", matches[-1])
    if match:
        threshold = int(match.group(1))
    elif "threshold_gev" in record:
        threshold = int(finite_float(record["threshold_gev"], f"samples.{sample_id}.threshold_gev"))
    else:
        threshold = -1

    if not system:
        if family == "photon" and threshold in {5, 10}:
            system = "pp"
        elif family == "jet" and threshold == 8:
            system = "pp"
    if not system or not family or threshold < 0:
        raise InputError(
            f"samples.{sample_id}: sample id/metadata must identify system, family, and threshold; "
            f"resolved system={system!r}, family={family!r}, threshold={threshold!r}"
        )
    canonical_id = f"{system}_{family}{threshold}"
    display_name = f"PhotonJet{threshold}" if family == "photon" else f"Jet{threshold}"
    return system, family, threshold, canonical_id, display_name


def load_samples(payload: Mapping[str, Any]) -> tuple[dict[str, SampleData], list[str]]:
    raw_samples = payload.get("samples")
    if not isinstance(raw_samples, Mapping):
        raise InputError("top-level 'samples' must be an object keyed by sample id")
    samples: dict[str, SampleData] = {}
    original_ids: list[str] = []
    for sample_id_raw, record_raw in raw_samples.items():
        sample_id = str(sample_id_raw)
        if not isinstance(record_raw, Mapping):
            raise InputError(f"samples.{sample_id}: expected an object")
        record = record_raw
        system, family, threshold, canonical_id, display_name = classify_sample(sample_id, record)
        if canonical_id in samples:
            raise InputError(
                f"duplicate canonical sample {canonical_id}: {samples[canonical_id].sample_id!r} and {sample_id!r}"
            )
        edges = numeric_array(record.get("bin_edges_gev"), f"samples.{sample_id}.bin_edges_gev")
        density = numeric_array(
            record.get("generator_stitching_density_pb_per_gev"),
            f"samples.{sample_id}.generator_stitching_density_pb_per_gev",
        )
        sumw2 = numeric_array(
            record.get("generator_stitching_density_sumw2_pb2_per_gev2"),
            f"samples.{sample_id}.generator_stitching_density_sumw2_pb2_per_gev2",
        )
        raw_counts = numeric_array(record.get("raw_counts"), f"samples.{sample_id}.raw_counts")
        nbins = len(edges) - 1
        if len(edges) < 2 or np.any(np.diff(edges) <= 0.0):
            raise InputError(f"samples.{sample_id}.bin_edges_gev must be strictly increasing")
        if any(len(array) != nbins for array in (density, sumw2, raw_counts)):
            raise InputError(
                f"samples.{sample_id}: bin arrays must have {nbins} entries for {len(edges)} bin edges"
            )
        if np.any(density < 0.0) or np.any(sumw2 < 0.0) or np.any(raw_counts < 0.0):
            raise InputError(f"samples.{sample_id}: density, sumw2, and raw counts must be non-negative")
        embedded_inclusive = system == "auau" and family == "jet"
        if embedded_inclusive:
            if "cross_section_weight_pb_per_event" in record:
                raise InputError(
                    f"samples.{sample_id}: legacy cross_section/generated_events weight is forbidden"
                )
            cross_section = finite_float(
                record.get("ownership_effective_cross_section_pb"),
                f"samples.{sample_id}.ownership_effective_cross_section_pb",
                positive=True,
            )
            normalization_denominator = finite_float(
                record.get("normalization_denominator_events"),
                f"samples.{sample_id}.normalization_denominator_events",
                positive=True,
            )
            weight = finite_float(
                record.get("stitching_weight_pb_per_owned_event"),
                f"samples.{sample_id}.stitching_weight_pb_per_owned_event",
                positive=True,
            )
            if (
                record.get("analysis_weight_state")
                != "CANONICAL_SOURCE_STITCH_ONLY__CENTRALITY_NOT_APPLIED"
                or record.get("nominal_downstream_analysis_ready") is not False
                or record.get("generator_stitching_channel", {}).get(
                    "centrality_reweighting_applied"
                )
                is not False
            ):
                raise InputError(
                    f"samples.{sample_id}: embedded-inclusive stitch/centrality state differs"
                )
        else:
            cross_section = finite_float(
                record.get("cross_section_pb"), f"samples.{sample_id}.cross_section_pb", positive=True
            )
            normalization_denominator = finite_float(
                record.get("generated_events"), f"samples.{sample_id}.generated_events", positive=True
            )
            weight = finite_float(
                record.get("cross_section_weight_pb_per_event"),
                f"samples.{sample_id}.cross_section_weight_pb_per_event",
                positive=True,
            )
        derived_weight = cross_section / normalization_denominator
        if not math.isclose(weight, derived_weight, rel_tol=2.0e-9, abs_tol=1.0e-18):
            raise InputError(
                f"samples.{sample_id}: accepted per-owned-event weight={weight:.12g} does not match "
                f"the bound cross section/denominator={derived_weight:.12g}"
            )
        source_hashes = record.get("source_product_hashes")
        if source_hashes in (None, {}, []):
            raise InputError(f"samples.{sample_id}.source_product_hashes is empty")
        window = parse_window(record.get("ownership_window"), f"samples.{sample_id}.ownership_window")
        samples[canonical_id] = SampleData(
            sample_id=sample_id,
            canonical_id=canonical_id,
            system=system,
            family=family,
            threshold=threshold,
            display_name=display_name,
            edges=edges,
            density=density,
            sumw2=sumw2,
            raw_counts=raw_counts,
            normalization_cross_section_pb=cross_section,
            normalization_denominator_events=normalization_denominator,
            normalization_weight_pb_per_owned_event=weight,
            window=window,
            source_product_hashes=source_hashes,
        )
        original_ids.append(sample_id)
    missing = sorted(REQUIRED_CANONICAL_IDS - set(samples))
    extra = sorted(set(samples) - REQUIRED_CANONICAL_IDS)
    if missing or extra:
        raise InputError(f"sample inventory mismatch: missing={missing}, unexpected={extra}")
    return samples, sorted(original_ids)


def validate_assembled_contract(payload: Mapping[str, Any]) -> None:
    if payload.get("schema") != ASSEMBLED_SCHEMA:
        raise InputError(
            f"assembled schema mismatch: expected={ASSEMBLED_SCHEMA!r}, actual={payload.get('schema')!r}"
        )
    if payload.get("status") != ASSEMBLED_STATUS:
        raise InputError(
            f"assembled status mismatch: expected={ASSEMBLED_STATUS!r}, actual={payload.get('status')!r}"
        )
    for field, expected in EXPECTED_ASSEMBLY_COUNTS.items():
        actual = payload.get(field)
        if isinstance(actual, bool) or not isinstance(actual, int) or actual != expected:
            raise InputError(f"assembled {field} mismatch: expected={expected}, actual={actual!r}")


def resolve_catalog(payload: Mapping[str, Any]) -> tuple[Path, str]:
    path_value = payload.get("accepted_catalog_path")
    hash_value = payload.get("accepted_catalog_sha256")
    if not path_value or not hash_value:
        raise InputError(
            "assembled input must provide top-level accepted_catalog_path and accepted_catalog_sha256"
        )
    catalog_path = Path(str(path_value)).expanduser()
    if not catalog_path.is_absolute():
        catalog_path = REPO / catalog_path
    catalog_hash = str(hash_value).lower()
    if not re.fullmatch(r"[0-9a-f]{64}", catalog_hash):
        raise InputError(f"accepted_catalog_sha256 is not a SHA-256 digest: {hash_value!r}")
    if catalog_path.exists():
        actual_hash = sha256_file(catalog_path)
        if actual_hash != catalog_hash:
            raise InputError(
                f"catalog hash mismatch for {catalog_path}: assembled={catalog_hash}, local={actual_hash}"
            )
    return catalog_path, catalog_hash


def load_input(
    path: Path,
    assembly_receipt_path: Path,
) -> tuple[dict[str, Any], dict[str, SampleData], Path, str, Path, str]:
    if not path.is_file():
        raise InputError(f"assembled input does not exist: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise InputError(f"assembled input is not valid JSON: {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise InputError("assembled input root must be a JSON object")
    validate_assembled_contract(payload)
    assembly_receipt_path = assembly_receipt_path.expanduser().resolve()
    if not assembly_receipt_path.is_file():
        raise InputError(f"assembled terminal receipt does not exist: {assembly_receipt_path}")
    try:
        terminal_receipt = json.loads(assembly_receipt_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise InputError(
            f"assembled terminal receipt is not valid JSON: {assembly_receipt_path}: {exc}"
        ) from exc
    if not isinstance(terminal_receipt, dict) or (
        terminal_receipt.get("schema") != "THE248Schema10StitchingAssemblyTerminalReceiptV4"
        or terminal_receipt.get("status") != "PASS"
        or Path(str(terminal_receipt.get("assembled_path"))).resolve() != path.resolve()
        or terminal_receipt.get("assembled_sha256") != sha256_file(path)
        or terminal_receipt.get("normalization_schema")
        != "THE248Schema10SplitStitchingNormalizationV4"
    ):
        raise InputError("assembled terminal receipt identity or hash binding differs")
    source_binding = (
        payload.get("normalization", {})
        .get("auau_embedded_inclusive", {})
        .get("source_artifact")
    )
    if not isinstance(source_binding, Mapping) or (
        terminal_receipt.get("auau_embedded_inclusive_source_stitch") != source_binding
    ):
        raise InputError("embedded-inclusive source-stitch receipt binding differs")
    try:
        from scripts.data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (
            load_source_stitch_artifact,
        )

        source = load_source_stitch_artifact(
            Path(str(source_binding.get("assembly_path"))),
            Path(str(source_binding.get("validation_receipt_path"))),
            verify_authority_file=True,
        )
    except (ImportError, OSError, ValueError) as exc:
        raise InputError(f"embedded-inclusive source-stitch chain is invalid: {exc}") from exc
    if source.binding_payload() != dict(source_binding):
        raise InputError("embedded-inclusive source-stitch chain differs after reconstruction")
    catalog_path, catalog_hash = resolve_catalog(payload)
    samples, original_ids = load_samples(payload)
    return (
        payload,
        samples,
        catalog_path,
        catalog_hash,
        assembly_receipt_path,
        sha256_file(assembly_receipt_path),
    )


def family_samples(samples: Mapping[str, SampleData], system: str, family: str) -> list[SampleData]:
    return sorted(
        (sample for sample in samples.values() if sample.system == system and sample.family == family),
        key=lambda sample: (sample.window.low, sample.threshold),
    )


def aggregate_owned_points(
    samples: Sequence[SampleData], xlim: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    by_x: dict[float, list[float]] = {}
    for sample in samples:
        x = sample.centers
        mask = sample.window.contains(x)
        mask &= (x >= xlim[0]) & (x <= xlim[1])
        mask &= np.isfinite(sample.density) & (sample.density > 0.0)
        for xx, yy, vv in zip(x[mask], sample.density[mask], sample.sumw2[mask]):
            key = round(float(xx), 10)
            acc = by_x.setdefault(key, [0.0, 0.0])
            acc[0] += float(yy)
            acc[1] += float(vv)
    if not by_x:
        raise InputError(f"no positive ownership-selected points in x range {xlim}")
    xs = np.array(sorted(by_x), dtype=float)
    ys = np.array([by_x[float(x)][0] for x in xs], dtype=float)
    variances = np.array([by_x[float(x)][1] for x in xs], dtype=float)
    return xs, ys, variances


def fit_log_pt_cubic(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, Any]:
    mask = (x > 0.0) & (y > 0.0) & np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 8:
        raise InputError("diagnostic fit requires at least eight finite positive owned bins")
    coefficients = np.polyfit(np.log(x[mask]), np.log(y[mask]), deg=int(FIT_CONTRACT["degree"]))

    def evaluate(values: np.ndarray | float) -> np.ndarray:
        arr = np.asarray(values, dtype=float)
        return np.exp(np.polyval(coefficients, np.log(arr)))

    return coefficients, evaluate


def boundary_diagnostics(samples: Sequence[SampleData], fit: Any, xlim: tuple[float, float]) -> list[dict[str, Any]]:
    diagnostics: list[dict[str, Any]] = []
    ordered = sorted(samples, key=lambda sample: (sample.window.low, sample.threshold))
    for left, right in zip(ordered[:-1], ordered[1:]):
        boundary = left.window.high
        if boundary is None or not math.isclose(boundary, right.window.low, abs_tol=1.0e-8):
            continue
        left_mask = left.window.contains(left.centers) & (left.centers >= xlim[0]) & (left.centers <= xlim[1])
        right_mask = right.window.contains(right.centers) & (right.centers >= xlim[0]) & (right.centers <= xlim[1])
        left_idx = np.flatnonzero(left_mask & (left.density > 0.0) & (left.centers < boundary))
        right_idx = np.flatnonzero(right_mask & (right.density > 0.0) & (right.centers >= boundary))
        if not len(left_idx) or not len(right_idx):
            diagnostics.append(
                {
                    "boundary_gev": boundary,
                    "left_sample": left.sample_id,
                    "right_sample": right.sample_id,
                    "status": "missing_adjacent_positive_bin",
                }
            )
            continue
        li = int(left_idx[-1])
        ri = int(right_idx[0])
        left_ref = float(fit(left.centers[li]))
        right_ref = float(fit(right.centers[ri]))
        diagnostics.append(
            {
                "boundary_gev": boundary,
                "left_sample": left.sample_id,
                "right_sample": right.sample_id,
                "left_bin_center_gev": float(left.centers[li]),
                "right_bin_center_gev": float(right.centers[ri]),
                "left_data_over_fit": float(left.density[li] / left_ref),
                "right_data_over_fit": float(right.density[ri] / right_ref),
                "status": "diagnostic_only_not_a_validity_gate",
            }
        )
    return diagnostics


def render_plot(samples: Sequence[SampleData], config: PanelConfig) -> tuple[Image.Image, dict[str, Any]]:
    x_combined, y_combined, variance_combined = aggregate_owned_points(samples, config.xlim)
    coefficients, fit = fit_log_pt_cubic(x_combined, y_combined)
    reference = fit(x_combined)
    ratios = y_combined / reference
    ratio_rms = float(np.sqrt(np.mean(np.square(ratios - 1.0))))
    boundaries = boundary_diagnostics(samples, fit, config.xlim)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(9.6, 5.9), dpi=125, facecolor="white")
    grid = fig.add_gridspec(
        2,
        1,
        height_ratios=[3.5, 1.1],
        hspace=0.055,
        left=0.15,
        right=0.965,
        top=0.90,
        bottom=0.14,
    )
    ax = fig.add_subplot(grid[0])
    ratio_ax = fig.add_subplot(grid[1], sharex=ax)

    fit_x = np.linspace(config.xlim[0], config.xlim[1], 600)
    ax.plot(
        fit_x,
        fit(fit_x),
        color="black",
        linewidth=1.8,
        linestyle="--",
        label="log-$p_T$ cubic fit (diagnostic)",
        zorder=2,
    )
    for sample in samples:
        x = sample.centers
        mask = sample.window.contains(x)
        mask &= (x >= config.xlim[0]) & (x <= config.xlim[1])
        mask &= sample.density > 0.0
        if not np.any(mask):
            continue
        color, marker = SAMPLE_STYLE[sample.canonical_id]
        ax.errorbar(
            x[mask],
            sample.density[mask],
            yerr=sample.errors[mask],
            fmt=marker,
            linestyle="none",
            markersize=4.8,
            linewidth=0.9,
            capsize=0,
            color=color,
            label=sample.display_name,
            zorder=4,
        )
        sample_ref = fit(x[mask])
        ratio_ax.errorbar(
            x[mask],
            sample.density[mask] / sample_ref,
            yerr=sample.errors[mask] / sample_ref,
            fmt=marker,
            linestyle="none",
            markersize=4.1,
            linewidth=0.8,
            capsize=0,
            color=color,
            zorder=4,
        )

    positive = y_combined[y_combined > 0.0]
    ax.set_yscale("log")
    ax.set_xlim(*config.xlim)
    ax.set_ylim(max(float(np.min(positive)) * 0.25, 1.0e-12), float(np.max(positive)) * 4.2)
    ax.set_ylabel(r"$d\sigma/dp_T$ [pb/GeV]", fontsize=20.5, fontweight="bold")
    ax.tick_params(labelsize=17.0, which="both")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.grid(axis="y", which="major", linestyle=":", linewidth=0.65, color="0.86")

    clipped_rms = min(ratio_rms, 0.44 * (config.ratio_ylim[1] - config.ratio_ylim[0]))
    ratio_ax.axhspan(1.0 - clipped_rms, 1.0 + clipped_rms, color="#dbeafe", alpha=0.72, zorder=0)
    ratio_ax.axhline(1.0, color="black", linewidth=1.0, zorder=1)
    ratio_ax.set_ylim(*config.ratio_ylim)
    ratio_ax.set_ylabel("data / fit", fontsize=18.0, fontweight="bold")
    ratio_ax.set_xlabel(config.xlabel, fontsize=21.5, fontweight="bold")
    ratio_ax.tick_params(labelsize=16.5, which="both")
    ratio_ax.grid(axis="y", which="major", linestyle=":", linewidth=0.65, color="0.86")
    ratio_ax.text(
        0.98,
        0.83,
        f"fit-shape RMS = {100.0 * ratio_rms:.1f}% (diagnostic)",
        transform=ratio_ax.transAxes,
        fontsize=15.5,
        ha="right",
        va="top",
    )

    label_box = {
        "boxstyle": "round,pad=0.28",
        "facecolor": "white",
        "edgecolor": "#cbd5e1",
        "alpha": 0.94,
        "linewidth": 0.8,
    }
    ax.text(
        0.965,
        0.955,
        r"$\it{\bf{sPHENIX}}$ Internal" + f"\n{config.sample_label}",
        transform=ax.transAxes,
        fontsize=17.5,
        ha="right",
        va="top",
        linespacing=1.05,
        bbox=label_box,
        zorder=7,
    )
    ax.text(
        0.965,
        0.735,
        config.title,
        transform=ax.transAxes,
        fontsize=18.0,
        fontweight="bold",
        ha="right",
        va="top",
        bbox=label_box,
        zorder=7,
    )
    ax.legend(
        frameon=False,
        fontsize=15.0,
        loc="lower left",
        ncol=2 if len(samples) > 2 else 1,
        handlelength=1.35,
        columnspacing=0.85,
        labelspacing=0.3,
    )

    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=125, facecolor="white")
    plt.close(fig)
    buffer.seek(0)
    image = Image.open(buffer).convert("RGB")
    metrics = {
        "fit_contract": dict(FIT_CONTRACT),
        "fit_coefficients_highest_power_first": [float(value) for value in coefficients],
        "fit_x_range_gev": list(config.xlim),
        "owned_positive_point_count": int(len(x_combined)),
        "fit_shape_rms_fraction": ratio_rms,
        "boundary_readouts": boundaries,
    }
    return image, metrics


def pil_font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    path = TIMES
    if bold and italic:
        path = TIMES_BOLD_ITALIC
    elif bold:
        path = TIMES_BOLD
    elif italic:
        path = TIMES_ITALIC
    try:
        return ImageFont.truetype(str(path), size)
    except OSError:
        return ImageFont.load_default()


FONTS = {
    "title": pil_font(70, bold=True),
    "subtitle": pil_font(38),
    "section": pil_font(42, bold=True),
    "table_header": pil_font(37, bold=True),
    "table_body": pil_font(37),
    "table_body_small": pil_font(37),
    "band_label": pil_font(39, bold=True),
    "band_body": pil_font(37),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    fill: tuple[int, int, int],
    outline: tuple[int, int, int] = LINE,
    width: int = 2,
    radius: int = 12,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def paste_contained(
    canvas: Image.Image, image: Image.Image, box: tuple[int, int, int, int]
) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    width, height = x1 - x0, y1 - y0
    scale = min(width / image.width, height / image.height)
    resized = image.resize((int(image.width * scale), int(image.height * scale)), Image.Resampling.LANCZOS)
    paste_x = x0 + (width - resized.width) // 2
    paste_y = y0 + (height - resized.height) // 2
    canvas.paste(resized, (paste_x, paste_y))
    return (paste_x, paste_y, paste_x + resized.width, paste_y + resized.height)


def wrap_pixels(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.FreeTypeFont, max_width: int) -> list[str]:
    lines: list[str] = []
    for paragraph in text.split("\n"):
        words = paragraph.split()
        line = ""
        for word in words:
            candidate = f"{line} {word}".strip()
            if line and draw.textbbox((0, 0), candidate, font=font)[2] > max_width:
                lines.append(line)
                line = word
            else:
                line = candidate
        lines.append(line)
    return lines


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font: ImageFont.FreeTypeFont,
    max_width: int,
    *,
    fill: tuple[int, int, int] = INK,
    line_gap: int = 4,
) -> int:
    x, y = xy
    for line in wrap_pixels(draw, text, font, max_width):
        draw.text((x, y), line, font=font, fill=fill)
        y += font.size + line_gap
    return y


def text_bbox(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font: ImageFont.FreeTypeFont,
) -> tuple[int, int, int, int]:
    return tuple(int(value) for value in draw.textbbox(xy, text, font=font))


def wrapped_text_bbox(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font: ImageFont.FreeTypeFont,
    max_width: int,
    *,
    line_gap: int = 4,
) -> tuple[int, int, int, int]:
    x, y = xy
    boxes: list[tuple[int, int, int, int]] = []
    for line in wrap_pixels(draw, text, font, max_width):
        boxes.append(text_bbox(draw, (x, y), line, font))
        y += font.size + line_gap
    return (
        min(box[0] for box in boxes),
        min(box[1] for box in boxes),
        max(box[2] for box in boxes),
        max(box[3] for box in boxes),
    )


def fmt_sig(value: float, digits: int = 6) -> str:
    if value == 0.0:
        return "0"
    if abs(value) >= 1.0e5 or abs(value) < 1.0e-3:
        return f"{value:.{digits - 1}e}"
    return f"{value:.{digits}g}"


def fmt_rel(value: float) -> str:
    if value >= 1.0e4:
        return f"{value:.3e}x"
    if value >= 100.0:
        return f"{value:.0f}x"
    if value >= 10.0:
        return f"{value:.1f}x"
    return f"{value:.2f}x"


def draw_weight_table(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    prefix: str,
    title: str,
    samples: Sequence[SampleData],
    accent: tuple[int, int, int],
    soft: tuple[int, int, int],
) -> list[dict[str, Any]]:
    x0, y0, x1, y1 = box
    nodes: list[dict[str, Any]] = [
        {
            "name": f"{prefix} weight table",
            "kind": "panel",
            "bbox": list(box),
            "symmetry_group": "weight tables",
            **({"title_axis_align": "left"} if prefix == "photon" else {}),
        }
    ]
    rounded(draw, box, soft, outline=accent, width=3, radius=12)
    draw.rectangle((x0 + 22, y0 + 20, x0 + 36, y0 + 70), fill=accent)
    draw.text((x0 + 57, y0 + 13), title, font=FONTS["section"], fill=INK)
    nodes.append(
        {
            "name": f"{prefix} table title",
            "kind": "text",
            "role": "audience",
            "font_px": FONTS["section"].size,
            "bbox": list(text_bbox(draw, (x0 + 57, y0 + 13), title, FONTS["section"])),
            "text": title,
        }
    )
    header_y = y0 + 80
    offsets = [26, 217, 350, 555, 735, 950]
    columns = [x0 + offset for offset in offsets]
    headers = ["sample", "window", "σ [pb]", "denom.", "w [pb/event]", "relative"]
    header_boxes: list[tuple[int, int, int, int]] = []
    for x, header in zip(columns, headers):
        draw.text((x, header_y), header, font=FONTS["table_header"], fill=MUTED)
        header_boxes.append(text_bbox(draw, (x, header_y), header, FONTS["table_header"]))
    if any(left[2] + 8 > right[0] for left, right in zip(header_boxes[:-1], header_boxes[1:])):
        raise InputError(f"{prefix} table header columns overlap")
    if header_boxes[-1][2] > x1 - 20:
        raise InputError(f"{prefix} table header exceeds panel")
    nodes.append(
        {
            "name": f"{prefix} table headers",
            "kind": "text",
            "role": "audience",
            "font_px": FONTS["table_header"].size,
            "bbox": [
                min(item[0] for item in header_boxes),
                min(item[1] for item in header_boxes),
                max(item[2] for item in header_boxes),
                max(item[3] for item in header_boxes),
            ],
            "text": " | ".join(headers),
        }
    )
    draw.line((x0 + 25, header_y + 46, x1 - 25, header_y + 46), fill=accent, width=2)
    anchor = samples[-1].normalization_weight_pb_per_owned_event
    row_start = header_y + 58
    row_step = max(39, min(43, int((y1 - row_start - 8) / len(samples))))
    row_font = FONTS["table_body"] if len(samples) <= 3 else FONTS["table_body_small"]
    for index, sample in enumerate(samples):
        y = row_start + index * row_step
        if index % 2 == 0:
            draw.rounded_rectangle((x0 + 22, y - 5, x1 - 22, y + row_step - 4), radius=6, fill=WHITE)
        color_hex = SAMPLE_STYLE[sample.canonical_id][0].lstrip("#")
        sample_color = tuple(int(color_hex[offset : offset + 2], 16) for offset in (0, 2, 4))
        values = [
            sample.display_name,
            sample.window.display(),
            fmt_sig(sample.normalization_cross_section_pb),
            f"{int(round(sample.normalization_denominator_events)):,}",
            fmt_sig(sample.normalization_weight_pb_per_owned_event),
            fmt_rel(sample.normalization_weight_pb_per_owned_event / anchor),
        ]
        value_boxes: list[tuple[int, int, int, int]] = []
        for col_index, (x, value) in enumerate(zip(columns, values)):
            draw.text((x, y), value, font=row_font, fill=sample_color if col_index == 0 else INK)
            value_boxes.append(text_bbox(draw, (x, y), value, row_font))
        if any(left[2] + 8 > right[0] for left, right in zip(value_boxes[:-1], value_boxes[1:])):
            raise InputError(f"{prefix} {sample.canonical_id} table columns overlap")
        if value_boxes[-1][2] > x1 - 20:
            raise InputError(f"{prefix} {sample.canonical_id} table row exceeds panel")
        nodes.append(
            {
                "name": f"{prefix} {sample.canonical_id} table row",
                "kind": "text",
                "role": "audience",
                "font_px": row_font.size,
                "bbox": [
                    min(item[0] for item in value_boxes),
                    min(item[1] for item in value_boxes),
                    max(item[2] for item in value_boxes),
                    max(item[3] for item in value_boxes),
                ],
                "text": " | ".join(values),
            }
        )
    return nodes


def boundary_summary(panel_name: str, diagnostics: Sequence[Mapping[str, Any]]) -> str:
    pieces: list[str] = []
    for item in diagnostics:
        boundary = item.get("boundary_gev")
        if item.get("status") != "diagnostic_only_not_a_validity_gate":
            pieces.append(f"{boundary:g} GeV: unavailable")
            continue
        pieces.append(
            f"{float(boundary):g} GeV: {float(item['left_data_over_fit']):.3f} -> "
            f"{float(item['right_data_over_fit']):.3f}"
        )
    return f"{panel_name} adjacent-bin data/fit: " + (", ".join(pieces) if pieces else "no internal boundary")


def render_slide(
    system: str,
    samples: Mapping[str, SampleData],
) -> tuple[Image.Image, dict[str, Any], dict[str, Any]]:
    photon_samples = family_samples(samples, system, "photon")
    jet_samples = family_samples(samples, system, "jet")
    photon_plot, photon_metrics = render_plot(photon_samples, PANEL_CONFIGS[(system, "photon")])
    jet_plot, jet_metrics = render_plot(jet_samples, PANEL_CONFIGS[(system, "jet")])

    canvas = Image.new("RGB", (W, H), WHITE)
    draw = ImageDraw.Draw(canvas)
    if system == "pp":
        title = "Current pp stitch-weight closure uses explicit ownership windows"
        subtitle = (
            "Per-event weight = cross section / generated events; each sample contributes only inside "
            "its declared truth-pT window."
        )
        photon_title = "PhotonJet weights; windows in GeV"
        jet_title = "Inclusive-jet weights; windows in GeV"
        photon_accent, jet_accent = PP_ACCENT, PP_JET_ACCENT
    else:
        title = "Embedded Au+Au stitch-weight closure uses explicit ownership windows"
        subtitle = (
            "Per-event weight = cross section / generated events; each sample contributes only inside "
            "its declared truth-pT window."
        )
        photon_title = "Embedded PhotonJet weights; windows in GeV"
        jet_title = "Embedded inclusive-jet weights; windows in GeV"
        photon_accent, jet_accent = AUAU_ACCENT, AUAU_JET_ACCENT

    title_xy = (72, 28)
    subtitle_xy = (74, 128)
    draw.text(title_xy, title, font=FONTS["title"], fill=INK)
    draw.text(subtitle_xy, subtitle, font=FONTS["subtitle"], fill=MUTED)
    nodes: list[dict[str, Any]] = [
        {
            "name": "slide title",
            "kind": "text",
            "role": "title",
            "title_anchor": True,
            "font_px": FONTS["title"].size,
            "bbox": list(text_bbox(draw, title_xy, title, FONTS["title"])),
            "text": title,
        },
        {
            "name": "weight definition",
            "kind": "text",
            "role": "audience",
            "font_px": FONTS["subtitle"].size,
            "bbox": list(text_bbox(draw, subtitle_xy, subtitle, FONTS["subtitle"])),
            "text": subtitle,
            "title_axis_align": "left",
        },
    ]

    plot_boxes = [(72, 180, 1246, 880), (1314, 180, 2488, 880)]
    for prefix, box, plot, accent in zip(
        ("photon", "inclusive jet"),
        plot_boxes,
        (photon_plot, jet_plot),
        (photon_accent, jet_accent),
    ):
        rounded(draw, box, WHITE, outline=accent, width=3, radius=12)
        ink_box = paste_contained(canvas, plot, (box[0] + 12, box[1] + 12, box[2] - 12, box[3] - 12))
        panel_name = f"{prefix} plot frame"
        nodes.extend(
            [
                {
                    "name": panel_name,
                    "kind": "panel",
                    "bbox": list(box),
                    "symmetry_group": "stitching plot frames",
                    **({"title_axis_align": "left"} if prefix == "photon" else {}),
                },
                {
                    "name": f"{prefix} plot ink",
                    "kind": "figure_ink",
                    "parent": panel_name,
                    "bbox": list(ink_box),
                    "outer_frame": panel_name,
                    "max_outer_frame_top_gap_px": 24,
                },
                {
                    "name": f"{prefix} plot typography",
                    "kind": "text",
                    "role": "plot_annotation",
                    "font_px": 25,
                    "bbox": list(ink_box),
                    "text": "Axes, ticks, sample legend, fit label, and diagnostic ratio readout",
                },
            ]
        )

    nodes.extend(
        draw_weight_table(
            draw,
            (72, 910, 1246, 1265),
            "photon",
            photon_title,
            photon_samples,
            photon_accent,
            PHOTON_SOFT,
        )
    )
    nodes.extend(
        draw_weight_table(
            draw,
            (1314, 910, 2488, 1265),
            "inclusive jet",
            jet_title,
            jet_samples,
            jet_accent,
            JET_SOFT,
        )
    )
    note_box = (72, 1280, 2488, 1408)
    rounded(draw, note_box, NOTE_SOFT, outline=(147, 197, 253), width=2, radius=10)
    left_label = "Owned-bin handoffs"
    left_body = "Adjacent-bin ratios compare each handoff with the fixed smooth reference."
    right_label = "Weights remain fixed"
    right_body = "The smooth fit is diagnostic only."
    left_label_xy = (94, 1288)
    left_body_xy = (94, 1334)
    right_label_xy = (1640, 1288)
    right_body_xy = (1640, 1334)
    draw.text(left_label_xy, left_label, font=FONTS["band_label"], fill=OK_GREEN)
    draw_wrapped(draw, left_body_xy, left_body, FONTS["band_body"], 1450, line_gap=2)
    draw.text(right_label_xy, right_label, font=FONTS["band_label"], fill=WARN_AMBER)
    draw_wrapped(
        draw,
        right_body_xy,
        right_body,
        FONTS["band_body"],
        820,
        line_gap=2,
    )
    nodes.extend(
        [
            {
                "name": "diagnostic readout band",
                "kind": "panel",
                "bbox": list(note_box),
                "title_axis_align": "left",
            },
            {
                "name": "owned-bin handoff label",
                "kind": "text",
                "role": "audience",
                "font_px": FONTS["band_label"].size,
                "bbox": list(text_bbox(draw, left_label_xy, left_label, FONTS["band_label"])),
                "text": left_label,
            },
            {
                "name": "owned-bin handoff explanation",
                "kind": "text",
                "role": "audience",
                "font_px": FONTS["band_body"].size,
                "bbox": list(
                    wrapped_text_bbox(
                        draw,
                        left_body_xy,
                        left_body,
                        FONTS["band_body"],
                        1450,
                        line_gap=2,
                    )
                ),
                "text": left_body,
            },
            {
                "name": "fixed-weight label",
                "kind": "text",
                "role": "audience",
                "font_px": FONTS["band_label"].size,
                "bbox": list(text_bbox(draw, right_label_xy, right_label, FONTS["band_label"])),
                "text": right_label,
            },
            {
                "name": "fixed-weight explanation",
                "kind": "text",
                "role": "audience",
                "font_px": FONTS["band_body"].size,
                "bbox": list(
                    wrapped_text_bbox(
                        draw,
                        right_body_xy,
                        right_body,
                        FONTS["band_body"],
                        820,
                        line_gap=2,
                    )
                ),
                "text": right_body,
            },
        ]
    )
    layout = {
        "schema": "slide_layout_nodes_v1",
        "slide_identity": f"{system}_current_schema10_stitch_weight_closure",
        "google_slides_object_id": SOURCE_OBJECT_IDS[system],
        "slide_size": [W, H],
        "title_axis_x": 72,
        "title_axis_tolerance_px": 6,
        "minimum_title_font_px": 67,
        "minimum_audience_font_px": 37,
        "minimum_plot_annotation_font_px": 25,
        "nodes": nodes,
    }
    return canvas, {"photon": photon_metrics, "inclusive_jet": jet_metrics}, layout


def script_text(system: str, samples: Mapping[str, SampleData]) -> str:
    photon = family_samples(samples, system, "photon")
    jets = family_samples(samples, system, "jet")
    if system == "pp":
        heading = "# pp current stitch-weight closure speaker script"
        opening = (
            "Here I am checking the current pp stitching inputs. "
            "The left panel is the maximum truth-photon pT spectrum, and the right panel is the "
            "maximum R=0.4 truth-jet pT spectrum."
        )
    else:
        heading = "# Embedded Au+Au current stitch-weight closure speaker script"
        opening = (
            "Here I am checking the current embedded Au+Au stitching inputs. "
            "The left panel uses the maximum truth-filter photon pT, and the right panel uses the "
            "maximum R=0.4 truth-jet pT."
        )

    def sample_sentence(items: Sequence[SampleData]) -> str:
        clauses = [
            f"{item.display_name} owns {item.window.display()} GeV with sigma={fmt_sig(item.normalization_cross_section_pb)} pb "
            f"and {'Npass' if item.system == 'auau' and item.family == 'jet' else 'Ngen'}="
            f"{int(round(item.normalization_denominator_events)):,}"
            for item in items
        ]
        return "; ".join(clauses) + "."

    return (
        f"{heading}\n\n"
        f"{opening}\n\n"
        "Each colored point is taken from the exact sample named in the assembled schema-10 receipt. "
        "The table makes the ownership window, cross section, accepted normalization denominator, and resulting "
        "per-event weight explicit. No fitted rescaling is applied.\n\n"
        f"For the photon side, {sample_sentence(photon)}\n\n"
        f"For the inclusive-jet side, {sample_sentence(jets)}\n\n"
        "The dashed curve is a fixed cubic polynomial in log pT fitted only as a smooth visual reference. "
        "The ratio and adjacent-bin boundary readouts are diagnostics, not acceptance tests, and they never "
        "change the sample weights or ownership windows. The takeaway is that the current inputs are shown "
        "with their provenance and stitching contract exposed for review.\n"
    )


def sample_manifest(sample: SampleData) -> dict[str, Any]:
    return {
        "assembled_sample_id": sample.sample_id,
        "canonical_id": sample.canonical_id,
        "display_name": sample.display_name,
        "system": sample.system,
        "family": sample.family,
        "threshold": sample.threshold,
        "bin_count": int(len(sample.density)),
        "bin_edges_gev": [float(value) for value in sample.edges],
        "normalization_cross_section_pb": sample.normalization_cross_section_pb,
        "normalization_denominator_events": sample.normalization_denominator_events,
        "accepted_weight_pb_per_owned_event": sample.normalization_weight_pb_per_owned_event,
        "ownership_window_original": sample.window.original,
        "ownership_window_parsed_lower_inclusive_upper_inclusive": {
            "low_gev": sample.window.low,
            "high_gev": sample.window.high,
        },
        "source_product_hashes": sample.source_product_hashes,
    }


def validate_source_pngs() -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for system, path in SOURCE_PNGS.items():
        if not path.is_file():
            raise InputError(f"source PNG is missing: {path}")
        digest = sha256_file(path)
        if digest != EXPECTED_SOURCE_PNG_HASHES[system]:
            raise InputError(
                f"source PNG hash changed for {system}: expected={EXPECTED_SOURCE_PNG_HASHES[system]}, actual={digest}"
            )
        with Image.open(path) as image:
            dimensions = [int(image.width), int(image.height)]
        out[system] = {
            "google_slides_object_id": SOURCE_OBJECT_IDS[system],
            "source_png": str(path),
            "source_png_sha256": digest,
            "source_png_dimensions_px": dimensions,
            "reuse_contract": "visual hierarchy reference only; historical plot values are not current evidence",
        }
    return out


def atomic_write_bytes(path: Path, payload: bytes) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def atomic_write_text(path: Path, payload: str) -> None:
    atomic_write_bytes(path, payload.encode("utf-8"))


def encode_png(image: Image.Image) -> bytes:
    buffer = BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    return buffer.getvalue()


def reject_internal_task_ids(text: str, label: str) -> None:
    match = INTERNAL_TASK_ID_RE.search(text)
    if match:
        raise InputError(f"{label} exposes internal task id {match.group(0)!r}")


def build(args: argparse.Namespace) -> dict[str, Any]:
    input_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    (
        payload,
        samples,
        catalog_path,
        catalog_hash,
        assembly_receipt_path,
        assembly_receipt_hash,
    ) = load_input(input_path, args.assembly_receipt)
    source_records = validate_source_pngs()
    assembled_hash = sha256_file(input_path)
    generator_path = Path(__file__).resolve()
    generator_hash = sha256_file(generator_path)

    validation_summary = {
        "ok": True,
        "input": str(input_path),
        "assembled_input_sha256": assembled_hash,
        "assembled_terminal_receipt_path": str(assembly_receipt_path),
        "assembled_terminal_receipt_sha256": assembly_receipt_hash,
        "assembled_contract": {
            "schema": payload["schema"],
            "status": payload["status"],
            **{field: payload[field] for field in EXPECTED_ASSEMBLY_COUNTS},
        },
        "accepted_catalog_path": str(catalog_path),
        "accepted_catalog_sha256": catalog_hash,
        "sample_count": len(samples),
        "sample_ids": [sample.sample_id for sample in sorted(samples.values(), key=lambda item: item.canonical_id)],
        "source_pngs": source_records,
    }
    if args.validate_only:
        return validation_summary

    output_dir.mkdir(parents=True, exist_ok=True)
    rendered: dict[str, Image.Image] = {}
    plot_metrics: dict[str, Any] = {}
    layouts: dict[str, dict[str, Any]] = {}
    for system in ("pp", "auau"):
        rendered[system], plot_metrics[system], layouts[system] = render_slide(system, samples)
        for node in layouts[system]["nodes"]:
            if node.get("kind") == "text" and node.get("role") in {"title", "audience"}:
                reject_internal_task_ids(str(node.get("text") or ""), f"{system} layout node {node.get('name')}")

    png_records: dict[str, Any] = {}
    for system in ("pp", "auau"):
        png_path = output_dir / OUTPUT_PNG_NAMES[system]
        png_payload = encode_png(rendered[system])
        atomic_write_bytes(png_path, png_payload)
        with Image.open(png_path) as check:
            if check.size != (W, H):
                raise RuntimeError(f"rendered PNG has wrong dimensions: {png_path}: {check.size}")
            check.verify()
        png_records[system] = {
            "path": str(png_path),
            "sha256": sha256_bytes(png_payload),
            "dimensions_px": [W, H],
            "google_slides_object_id": SOURCE_OBJECT_IDS[system],
        }

    script_records: dict[str, Any] = {}
    for system in ("pp", "auau"):
        script_path = output_dir / SCRIPT_NAMES[system]
        generated_script = script_text(system, samples)
        reject_internal_task_ids(generated_script, f"{system} speaker script")
        atomic_write_text(script_path, generated_script)
        script_records[system] = {
            "path": str(script_path),
            "sha256": sha256_file(script_path),
        }

    layout_records: dict[str, Any] = {}
    for system in ("pp", "auau"):
        layout_path = output_dir / LAYOUT_NAMES[system]
        atomic_write_text(layout_path, json.dumps(layouts[system], indent=2, sort_keys=True) + "\n")
        layout_records[system] = {
            "path": str(layout_path),
            "sha256": sha256_file(layout_path),
            "google_slides_object_id": SOURCE_OBJECT_IDS[system],
        }

    manifest = {
        "schema": "THE248_CURRENT_SCHEMA10_STITCHING_SLIDES_V1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "assembled_input": {
            "path": str(input_path),
            "sha256": assembled_hash,
            "declared_schema": payload.get("schema"),
        },
        "accepted_catalog": {"path": str(catalog_path), "sha256": catalog_hash},
        "generator": {"path": str(generator_path), "sha256": generator_hash},
        "source_references": source_records,
        "outputs": {
            "full_slide_pngs": png_records,
            "speaker_scripts": script_records,
            "layout_nodes": layout_records,
            "png_count_created_by_renderer": 2,
            "intermediate_pngs_created": 0,
        },
        "fit_contract": dict(FIT_CONTRACT),
        "plot_metrics": plot_metrics,
        "samples": {
            sample.canonical_id: sample_manifest(sample)
            for sample in sorted(samples.values(), key=lambda item: item.canonical_id)
        },
        "observable_contract": {
            "pp_photon": "max truth photon pT",
            "pp_inclusive_jet": "max R=0.4 truth jet pT",
            "auau_photon": "max truth-filter photon pT under the assembled current truth contract",
            "auau_inclusive_jet": "max R=0.4 truth jet pT under the assembled current truth contract",
            "normalization_channel": "generator_stitching",
            "weighted_quantity": "pp/photon raw counts use their accepted generator denominator; Au+Au embedded inclusive raw counts use ownership-effective cross section divided by same-window Npass, with Poisson count variance",
        },
        "caveats": [
            "Historical source PNGs define visual hierarchy only; their May/July values are not current-production evidence.",
            "The log-pT cubic is a fixed-form diagnostic reference; it does not tune weights or ownership windows.",
            "Adjacent-bin data/fit boundary readouts are diagnostic only and are never validity gates.",
            "Scientific acceptance still requires the upstream source-coverage, duplicate-exclusion, additivity, finite-weight, and ownership checks bound by the assembled receipt.",
            "No Google Slides content was read or mutated by this renderer.",
        ],
        "google_slides_mutation": False,
        "slide_numbers_baked_into_png": False,
    }
    manifest_path = output_dir / MANIFEST_NAME
    atomic_write_text(manifest_path, json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {
        **validation_summary,
        "outputs": {
            "pp_png": png_records["pp"]["path"],
            "auau_png": png_records["auau"]["path"],
            "pp_script": script_records["pp"]["path"],
            "auau_script": script_records["auau"]["path"],
            "pp_layout_nodes": layout_records["pp"]["path"],
            "auau_layout_nodes": layout_records["auau"]["path"],
            "manifest": str(manifest_path),
        },
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="assembled THE-248 stitching JSON")
    parser.add_argument(
        "--assembly-receipt",
        type=Path,
        default=DEFAULT_ASSEMBLY_RECEIPT,
        help="hash-bound terminal receipt for the exact V4 assembly",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="output directory")
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="validate schema, hashes, sample inventory, arrays, weights, and source references without writing outputs",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = build(args)
    except (InputError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
