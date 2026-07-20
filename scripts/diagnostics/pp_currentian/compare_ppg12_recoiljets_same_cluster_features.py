#!/usr/bin/env python3
"""Compare PPG12 slimtree clusters against RecoilJets pp photon-ID rows.

It matches signal clusters by eventnumber + truth track id, using nearest
eta/phi only to disambiguate multiple reconstructed clusters attached to the
same truth photon.  Angular displacement is recorded as reconstruction
evidence; it is not allowed to manufacture a candidate-population mismatch.
The comparator re-evaluates both preserved PPG12 BDTs and records the first
population, feature, score, route, tag, isolation, or ABCD divergence candidate
by candidate.  Classification follows the preserved executable, not an
inferred historical plotting convention.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

try:
    import ROOT  # type: ignore
except ImportError:  # Pure selection-contract tests do not need ROOT.
    ROOT = None  # type: ignore


if ROOT is not None:
    ROOT.gROOT.SetBatch(True)

RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
TRUTH_BINS = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45]
PPG12_ETA_MIN = -0.7
PPG12_ETA_MAX = 0.7
PPG12_VERTEX_ABS_MAX = 60.0
PPG12_RECO_ET_MIN = 5.0
PPG12_MC_ISO_SCALE = 1.2
PPG12_MC_ISO_SHIFT = 0.1
PPG12_RECO_ISO_MIN = -20.0
PPG12_RECO_ISO_INTERCEPT = 0.490
PPG12_RECO_ISO_SLOPE = 0.037
PPG12_NONISO_SHIFT = 0.8
PPG12_NONISO_MAX = 20.0
GEOMETRY_DIAGNOSTIC_DR = 0.02
PRESERVED_EXECUTABLE_EVIDENCE = "preserved_ppg12_executable"
PYTHON_SHADOW_EVIDENCE = "python_shadow_not_executable"
UNAVAILABLE_EXECUTABLE_EVIDENCE = "unavailable_in_ppg12_slimtree"

DEFAULT_BASE_E_MODEL = (
    "/sphenix/user/shuhangli/ppg12/FunWithxgboost/"
    "binned_models/model_base_E_split_single_tmva.root"
)
DEFAULT_BASE_V3E_MODEL = (
    "/sphenix/user/shuhangli/ppg12/FunWithxgboost/"
    "binned_models/model_base_v3E_split_single_tmva.root"
)

BASE_E_FEATURES = [
    "cluster_Et_score_input",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
]
BASE_V3E_FEATURES = [
    "cluster_Et_score_input",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
]
FEATURES = [
    ("cluster_Et_score_input", "cluster_Et"),
    ("cluster_weta_cogx", "cluster_weta_cogx"),
    ("cluster_wphi_cogx", "cluster_wphi_cogx"),
    ("vertexz", "vertexz"),
    ("cluster_Eta", "cluster_Eta"),
    ("e11_over_e33", "e11_over_e33"),
    ("cluster_et1", "cluster_et1"),
    ("cluster_et2", "cluster_et2"),
    ("cluster_et3", "cluster_et3"),
    ("cluster_et4", "cluster_et4"),
    ("e32_over_e35", "e32_over_e35"),
]


def read_executable_trace(
    trace_path: Path, response_trace_path: Path
) -> dict[tuple[int, int], dict[str, Any]]:
    """Load the instrumented executable side channel by stable tree identity."""
    rows: dict[tuple[int, int], dict[str, Any]] = {}
    with trace_path.open(newline="") as stream:
        for line_number, raw in enumerate(csv.DictReader(stream), start=2):
            try:
                key = (int(raw["tree_entry"]), int(raw["cluster_index"]))
            except (KeyError, ValueError) as exc:
                raise RuntimeError(
                    f"invalid executable trace identity at row {line_number}"
                ) from exc
            if key in rows:
                raise RuntimeError(f"duplicate executable trace candidate identity: {key}")
            rows[key] = dict(raw)
    responses: dict[tuple[int, int], dict[str, str]] = {}
    with response_trace_path.open(newline="") as stream:
        for line_number, raw in enumerate(csv.DictReader(stream), start=2):
            try:
                key = (int(raw["tree_entry"]), int(raw["cluster_index"]))
            except (KeyError, ValueError) as exc:
                raise RuntimeError(
                    f"invalid executable response identity at row {line_number}"
                ) from exc
            if key in responses:
                raise RuntimeError(f"duplicate executable response identity: {key}")
            responses[key] = dict(raw)
    for key, response in responses.items():
        if key not in rows:
            raise RuntimeError(f"response trace has no candidate row: {key}")
        for field in (
            "chain_file_index",
            "local_tree_entry",
            "runnumber",
            "eventnumber",
        ):
            if response.get(field) != rows[key].get(field):
                raise RuntimeError(f"response trace {field} mismatch for {key}")
        rows[key]["response_Et"] = response["response_Et"]
        rows[key]["response_window_pass"] = response["response_window_pass"]
    return rows


def bind_executable_trace(
    ppg12_by_event: dict[tuple[int, int], list[dict[str, Any]]],
    trace: dict[tuple[int, int], dict[str, Any]],
) -> None:
    """Replace shadow decisions only after exact trace identity/value checks."""
    seen: set[tuple[int, int]] = set()

    def require_close(row: dict[str, Any], trace_row: dict[str, Any], left: str, right: str) -> None:
        left_value = float(row[left])
        right_value = float(trace_row[right])
        tolerance = max(1.0e-6, 1.0e-5 * abs(left_value))
        if not math.isfinite(right_value) or abs(left_value - right_value) > tolerance:
            raise RuntimeError(
                f"executable trace drift for entry={row['tree_entry']} cluster={row['cluster_index']} "
                f"field={left}/{right}: shadow={left_value} executable={right_value}"
            )

    for event_rows in ppg12_by_event.values():
        for row in event_rows:
            key = (int(row["tree_entry"]), int(row["cluster_index"]))
            trace_row = trace.get(key)
            if trace_row is None:
                raise RuntimeError(f"accepted PPG12 candidate missing executable trace: {key}")
            seen.add(key)
            if int(trace_row["eventnumber"]) != int(row["eventnumber"]):
                raise RuntimeError(f"executable eventnumber mismatch for {key}")
            if int(trace_row["is_signal"]) != int(row["is_signal"]):
                raise RuntimeError(f"executable signal-status mismatch for {key}")
            for left, right in (
                ("cluster_Et", "cluster_Et"),
                ("cluster_weta_cogx", "cluster_weta_cogx"),
                ("cluster_wphi_cogx", "cluster_wphi_cogx"),
                ("vertexz", "vertexz"),
                ("cluster_Eta", "cluster_Eta"),
                ("e11_over_e33", "e11_over_e33"),
                ("cluster_et1", "cluster_et1"),
                ("cluster_et2", "cluster_et2"),
                ("cluster_et3", "cluster_et3"),
                ("cluster_et4", "cluster_et4"),
                ("e32_over_e35", "e32_over_e35"),
                ("bdt_base_E", "base_E_score"),
                ("bdt_base_v3E", "base_v3E_score"),
            ):
                require_close(row, trace_row, left, right)
            tag = 0 if int(trace_row["common_pass"]) == 0 else (
                1 if int(trace_row["tight"]) else 2 if int(trace_row["nontight"]) else 3
            )
            row.update(
                chain_file_index=int(trace_row["chain_file_index"]),
                local_tree_entry=int(trace_row["local_tree_entry"]),
                selected_model=trace_row["selected_model"],
                selected_bdt_score=float(trace_row["selected_score"]),
                ppg12_recomputed_tag=tag,
                ppg12_common_pass=int(trace_row["common_pass"]),
                raw_eiso=float(trace_row["raw_eiso"]),
                corrected_eiso=float(trace_row["corrected_eiso"]),
                iso_threshold=float(trace_row["iso_threshold"]),
                noniso_threshold=float(trace_row["noniso_threshold"]),
                is_iso=int(trace_row["is_iso"]),
                is_noniso=int(trace_row["is_noniso"]),
                logical_abcd_region=int(trace_row["logical_abcd_region"]),
                ppg12_abcd_region={0: "none", 1: "A", 2: "B", 3: "C", 4: "D"}[
                    int(trace_row["logical_abcd_region"])
                ],
                is_signal=int(trace_row["is_signal"]),
                truth_particle_index=int(trace_row["truth_particle_index"]),
                truth_class=int(trace_row["truth_class"]),
                truth_pt=float(trace_row["truth_pt"]),
                analysis_window_pass=int(trace_row["analysis_window_pass"]),
                signal_fill_A=int(trace_row["signal_fill_A"]),
                signal_fill_B=int(trace_row["signal_fill_B"]),
                signal_fill_C=int(trace_row["signal_fill_C"]),
                signal_fill_D=int(trace_row["signal_fill_D"]),
                fill_multiplicity=int(trace_row["fill_multiplicity"]),
                weight_sample=float(trace_row["sample_weight"]),
                weight_mix=float(trace_row["mix_weight"]),
                weight_lumi=float(trace_row["lumi_weight"]),
                weight_cross=float(trace_row["cross_weight"]),
                weight_vertex=float(trace_row["vertex_weight"]),
                weight_truth_vertex=float(trace_row["truth_vertex_weight"]),
                weight_trigger=float(trace_row["trigger_weight"]),
                weight_event=float(trace_row["event_weight"]),
                weight_final=float(trace_row["weight"]),
                response_Et=trace_row.get("response_Et", ""),
                response_window_pass=trace_row.get("response_window_pass", ""),
                tag_evidence_source=PRESERVED_EXECUTABLE_EVIDENCE,
                isolation_abcd_evidence_source=PRESERVED_EXECUTABLE_EVIDENCE,
                truth_response_fill_evidence_source=PRESERVED_EXECUTABLE_EVIDENCE,
                weight_evidence_source=PRESERVED_EXECUTABLE_EVIDENCE,
            )
    extra = sorted(set(trace) - seen)
    if extra:
        raise RuntimeError(
            f"executable trace contains {len(extra)} candidates outside comparator universe; "
            f"first={extra[0]}"
        )


def require_root() -> Any:
    if ROOT is None:
        raise RuntimeError(
            "PyROOT is required for ROOT-file comparison; run with the analysis environment"
        )
    return ROOT


def find_bin(x: float) -> int:
    for i in range(len(RECO_BINS) - 1):
        if RECO_BINS[i] < x < RECO_BINS[i + 1]:
            return i
    return -1


def bin_label(i: int) -> str:
    return f"{RECO_BINS[i]}-{RECO_BINS[i + 1]}"


def finite(x: Any) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def open_interval(x: float, lo: float, hi: float) -> bool:
    return finite(x) and lo < x < hi


def cpp_float_ratio(numerator: float, denominator: float) -> float:
    """Mirror IEEE floating division used by the C++ executable."""
    if math.isnan(numerator) or math.isnan(denominator):
        return float("nan")
    if denominator == 0.0:
        if numerator == 0.0:
            return float("nan")
        sign = math.copysign(1.0, numerator) * math.copysign(1.0, denominator)
        return math.copysign(float("inf"), sign)
    return numerator / denominator


def safe_div(n: float, d: float) -> float:
    return n / d if d else float("nan")


def dphi(a: float, b: float) -> float:
    x = a - b
    while x > math.pi:
        x -= 2.0 * math.pi
    while x <= -math.pi:
        x += 2.0 * math.pi
    return x


def next_ppg12_segment(
    current_segment: int,
    previous_eventnumber: int | None,
    eventnumber: int,
) -> int:
    """Recover the source-file segment from the slimtree event sequence.

    PPG12's merged slimtrees omit aborted events, so a source file does not
    contribute a fixed number of tree entries.  The EventHeader sequence is
    monotonic within a source file and resets when the next file begins.
    Dividing the merged entry index by an assumed 1,000 entries therefore
    pairs different physical events after the first omitted event.
    """
    # Multiple reconstructed candidates from one event legitimately repeat the
    # same EventHeader number.  Only a strict decrease identifies the next
    # source file; treating equality as a reset would split one event across
    # artificial segments.
    if previous_eventnumber is not None and eventnumber < previous_eventnumber:
        return current_segment + 1
    return current_segment


def recoil_source_identity(
    processed_event_ordinal: int,
    local_eventnumber: int,
    events_per_segment: int,
) -> tuple[int, int]:
    """Recover the sealed source row and local EventHeader identity.

    RecoilJets ``evt`` is the monotonic process-event ordinal across the
    five-row source slice.  ``eventnumber`` is the EventHeader sequence and
    may restart for every source file; it must never be divided to infer a
    segment.  The two values together produce a stable five-file identity.
    """
    if events_per_segment <= 0:
        raise RuntimeError("events_per_segment must be positive for exact oracle identity")
    if processed_event_ordinal <= 0:
        raise RuntimeError(
            f"invalid RecoilJets processed-event ordinal: {processed_event_ordinal}"
        )
    if local_eventnumber <= 0:
        raise RuntimeError(f"invalid local EventHeader eventnumber: {local_eventnumber}")
    segment = (processed_event_ordinal - 1) // events_per_segment
    return int(segment), int(local_eventnumber)


def ppg12_event_in_scan(
    event_key: tuple[int, int],
    wanted_events: set[tuple[int, int]] | None,
) -> bool:
    """Return whether a PPG12 event belongs to the requested scan universe.

    ``None`` is the admission-safe symmetric mode: every event in the sealed
    PPG12 source slice is scanned, including events with no RecoilJets row.
    """
    return wanted_events is None or event_key in wanted_events


def in_ppg12_analysis_window(reco_et: float, truth_pt: float) -> int:
    return int(
        finite(reco_et)
        and finite(truth_pt)
        and RECO_BINS[0] < reco_et < RECO_BINS[-1]
        and TRUTH_BINS[0] < truth_pt < TRUTH_BINS[-1]
    )


def abcd_region_code(region: str) -> int:
    return {"none": 0, "A": 1, "B": 2, "C": 3, "D": 4}.get(region, -1)


def ppg12_signal_fill_flags(region: str, analysis_window_pass: int) -> dict[str, int]:
    flags = {name: 0 for name in "ABCD"}
    if region == "A":
        flags["A"] = int(analysis_window_pass == 1)
    elif region in flags:
        flags[region] = 1
    return flags


def selected_ppg12_model(analysis_et: float) -> str:
    """Return the executable PPG12 model route for calibrated, unsmeared ET."""
    # config_bdt_nom.yaml: fallback base_E, ET bins [8, 15, 35] both base_v3E.
    if 8.0 <= analysis_et < 35.0:
        return "base_v3E"
    return "base_E"


def selected_score(scores: dict[str, float], analysis_et: float) -> tuple[str, float]:
    model = selected_ppg12_model(analysis_et)
    return model, float(scores[model])


def infer_stored_model_route(
    stored_score: float,
    base_e_score: float,
    base_v3e_score: float,
    tolerance: float = 1.0e-6,
) -> str:
    """Infer which preserved model produced a stored RecoilJets score.

    The output tree does not carry a model-name branch.  When its stored score
    agrees uniquely with one of the two independently recomputed responses,
    that equality is sufficient to recover the deployed route without
    guessing from the candidate ET.
    """
    matches = []
    if finite(stored_score) and finite(base_e_score) and abs(stored_score - base_e_score) <= tolerance:
        matches.append("base_E")
    if finite(stored_score) and finite(base_v3e_score) and abs(stored_score - base_v3e_score) <= tolerance:
        matches.append("base_v3E")
    if len(matches) == 1:
        return matches[0]
    if len(matches) == 2:
        return "ambiguous"
    return "unmatched"


def passes_ppg12_kinematics(et: float, eta: float, vertexz: float) -> bool:
    """Mirror the nominal executable's inclusive ET, eta, and vertex boundaries."""
    return (
        finite(et)
        and et >= PPG12_RECO_ET_MIN
        and open_interval(eta, PPG12_ETA_MIN, PPG12_ETA_MAX)
        and finite(vertexz)
        and abs(vertexz) <= PPG12_VERTEX_ABS_MAX
    )


def ppg12_iso_assignment(raw_topo_iso04: float, et: float) -> dict[str, float | int | str]:
    """Apply nominal PPG12 MC isolation correction and strict ABCD intervals."""
    corrected = raw_topo_iso04 * PPG12_MC_ISO_SCALE + PPG12_MC_ISO_SHIFT
    iso_max = PPG12_RECO_ISO_INTERCEPT + PPG12_RECO_ISO_SLOPE * et
    noniso_min = iso_max + PPG12_NONISO_SHIFT
    is_iso = int(open_interval(corrected, PPG12_RECO_ISO_MIN, iso_max))
    is_noniso = int(open_interval(corrected, noniso_min, PPG12_NONISO_MAX))
    return {
        "raw_eiso": raw_topo_iso04,
        "corrected_eiso": corrected,
        "iso_threshold": iso_max,
        "noniso_threshold": noniso_min,
        "is_iso": is_iso,
        "is_noniso": is_noniso,
    }


def abcd_region(tag: int, is_iso: int, is_noniso: int) -> str:
    if tag == 1 and is_iso:
        return "A"
    if tag == 1 and is_noniso:
        return "B"
    if tag == 2 and is_iso:
        return "C"
    if tag == 2 and is_noniso:
        return "D"
    return "none"


def load_rbdt(model_path: str) -> Any:
    root = require_root()
    root.TMVA.Tools.Instance()
    probe = root.TFile.Open(model_path)
    if not probe or probe.IsZombie():
        raise RuntimeError(f"Could not open PPG12 TMVA model: {model_path}")
    probe.Close()
    return root.TMVA.Experimental.RBDT("myBDT", model_path)


def eval_rbdt(model: Any, feature_names: Iterable[str], row: dict[str, float]) -> float:
    root = require_root()
    values = root.std.vector("float")()
    for name in feature_names:
        value = float(row[name])
        if not finite(value):
            raise RuntimeError(f"Non-finite BDT input {name}={value}")
        values.push_back(value)
    result = model.Compute(values)
    return float(result[0]) if len(result) else float("nan")


def ppg12_thresholds(et: float) -> tuple[float, float, float]:
    tight_min = 0.815625 - 0.0015625 * et
    nt_min = 0.7333333333333333 - 0.01333333333333333 * et
    nt_max = 0.684375 + 0.0015625 * et
    return tight_min, nt_min, nt_max


def classify(row: dict[str, float], score: float, et_for_cuts: float) -> tuple[int, dict[str, bool]]:
    """Return PPG12 tag: 0 preselection fail, 1 tight, 2 nontight, 3 neither."""
    e11e33 = row["e11_over_e33"]
    e32e35 = row["e32_over_e35"]
    wr = cpp_float_ratio(row["cluster_wphi_cogx"], row["cluster_weta_cogx"])
    common = (
        open_interval(row["cluster_prob"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 0.98)
        and wr > 0.0
        and finite(row["cluster_weta_cogx"])
        and row["cluster_weta_cogx"] < 2.0
        and row["npb_score"] > 0.5
    )
    flags: dict[str, bool] = {"common": common}
    if not common:
        return 0, flags

    tight_min, nt_min, nt_max = ppg12_thresholds(et_for_cuts)
    tight_prob = open_interval(row["cluster_prob"], 0.0, 1.0)
    tight_weta = open_interval(row["cluster_weta_cogx"], 0.0, 1.0)
    tight_wphi = open_interval(row["cluster_wphi_cogx"], 0.0, 1.0)
    tight_et1 = open_interval(row["cluster_et1"], 0.5, 1.0)
    tight_et2 = open_interval(row["cluster_et2"], 0.0, 1.0)
    tight_et3 = open_interval(row["cluster_et3"], 0.0, 1.0)
    tight_et4 = open_interval(row["cluster_et4"], 0.0, 1.0)
    tight_e11e33 = open_interval(e11e33, 0.0, 1.0)
    tight_e32e35 = open_interval(e32e35, 0.8, 1.0)
    tight_bdt = finite(score) and score > tight_min and score < 1.0

    flags.update(
        tight_prob=tight_prob,
        tight_weta=tight_weta,
        tight_wphi=tight_wphi,
        tight_et1=tight_et1,
        tight_et2=tight_et2,
        tight_et3=tight_et3,
        tight_et4=tight_et4,
        tight_e11e33=tight_e11e33,
        tight_e32e35=tight_e32e35,
        tight_bdt=tight_bdt,
        nt_bdt=finite(score) and score > nt_min and score < nt_max,
    )
    tight = all(
        flags[k]
        for k in (
            "tight_prob",
            "tight_weta",
            "tight_wphi",
            "tight_et1",
            "tight_et2",
            "tight_et3",
            "tight_et4",
            "tight_e11e33",
            "tight_e32e35",
            "tight_bdt",
        )
    )
    if tight:
        return 1, flags

    nt_shape = (
        open_interval(row["cluster_prob"], 0.0, 1.0)
        and open_interval(row["cluster_weta_cogx"], 0.0, 1.0)
        and open_interval(row["cluster_wphi_cogx"], 0.0, 1.0)
        and open_interval(row["cluster_et1"], 0.6, 1.0)
        and open_interval(row["cluster_et4"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 1.0)
        and open_interval(e32e35, 0.8, 1.0)
    )
    nfail = 0
    if not tight_weta:
        nfail += 1
    if not tight_prob:
        nfail += 1
    if not tight_bdt:
        nfail += 1
    flags["nt_shape"] = nt_shape
    flags["nfail_any"] = nfail > 0
    if nt_shape and flags["nt_bdt"] and nfail > 0:
        return 2, flags
    return 3, flags


def tower_masked(mask: Any, ieta_value: float, iphi_value: float) -> bool:
    if not mask:
        return False
    ieta = int(ieta_value)
    iphi = int(iphi_value)
    if ieta < 0 or ieta >= mask.GetNbinsX():
        return False
    if iphi < 0 or iphi >= mask.GetNbinsY():
        return False
    return mask.GetBinContent(ieta + 1, iphi + 1) > 0


@dataclass
class Agg:
    n: int = 0
    rj_common: int = 0
    ppg12_common: int = 0
    rj_tight: int = 0
    rj_nontight: int = 0
    ppg12_tight: int = 0
    ppg12_nontight: int = 0
    model_route_agree: int = 0
    tag_agree: int = 0
    common_agree: int = 0
    abcd_agree: int = 0
    stored_tag_disagree: int = 0
    stored_common_disagree: int = 0
    model_base_e: int = 0
    model_base_v3e: int = 0
    score_diff_sum: float = 0.0
    score_abs_sum: float = 0.0
    rj_score_sum: float = 0.0
    ppg12_score_sum: float = 0.0
    feature_abs: dict[str, float] = field(default_factory=lambda: defaultdict(float))
    feature_signed: dict[str, float] = field(default_factory=lambda: defaultdict(float))
    feature_n: dict[str, int] = field(default_factory=lambda: defaultdict(int))

    def add_feature(self, name: str, rj: float, ppg12: float) -> None:
        if not (finite(rj) and finite(ppg12)):
            return
        diff = float(rj) - float(ppg12)
        self.feature_abs[name] += abs(diff)
        self.feature_signed[name] += diff
        self.feature_n[name] += 1


def read_rj_rows(
    path: str,
    events_per_segment: int,
    base_e_model: Any,
    base_v3e_model: Any,
    allow_legacy_score_et_fallback: bool = False,
) -> tuple[list[dict[str, Any]], set[tuple[int, int]]]:
    root = require_root()
    f = root.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open RecoilJets ROOT: {path}")
    t = f.Get("AuAuPhotonIDTrainingTree")
    if not t:
        raise RuntimeError(f"Missing AuAuPhotonIDTrainingTree in {path}")
    branches = {branch.GetName() for branch in t.GetListOfBranches()}
    has_score_input_et = "cluster_Et_score_input" in branches
    if not has_score_input_et and not allow_legacy_score_et_fallback:
        raise RuntimeError(
            "Missing cluster_Et_score_input. Exact PPG12 score comparison refuses the "
            "historical cluster_Et proxy; pass --allow-legacy-score-et-fallback only "
            "for explicitly labeled legacy diagnostics."
        )
    required = {
        "evt",
        "eventnumber",
        "is_signal",
        "truth_track_id",
        "cluster_Et",
        "cluster_Eta",
        "cluster_Phi",
        "ppg12_kin_vertexz",
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "e11_over_e33",
        "e32_over_e35",
        "cluster_et1",
        "cluster_et2",
        "cluster_et3",
        "cluster_et4",
        "cluster_prob",
        "npb_score",
        "tight_bdt_score",
        "ppg12_common_pass",
        "ppg12_tight_tag",
        "ppg12_response_Et",
        "ppg12_logical_abcd_region",
        "ppg12_analysis_window_pass",
        "ppg12_response_window_pass",
        "ppg12_signal_fill_A",
        "ppg12_signal_fill_B",
        "ppg12_signal_fill_C",
        "ppg12_signal_fill_D",
        "ppg12_signal_fill_multiplicity",
        "ppg12_truth_class",
        "ppg12_weight_lane_code",
        "ppg12_weight_component_code",
        "ppg12_weight_slice",
        "ppg12_weight_vertex",
        "ppg12_weight_mix",
        "ppg12_weight_period",
        "ppg12_weight_final",
        "event_weight",
        "ppg12_raw_eiso",
        "ppg12_reco_eiso",
        "ppg12_iso_threshold",
        "ppg12_noniso_threshold",
        "ppg12_is_iso",
        "ppg12_is_noniso",
        "cluster_index",
        "ppg12_sample_bin",
        "ppg12_xsec_pb",
        "ppg12_xsec_weight",
        "ppg12_window_low",
        "ppg12_window_high",
        "ppg12_truth_window_pass_r04",
    }
    missing = sorted(required - branches)
    if missing:
        raise RuntimeError(f"Missing required RecoilJets oracle branches: {missing}")

    optional_float_branches = (
        "event_weight",
        "reco_eiso",
        "ppg12_raw_eiso",
        "ppg12_reco_eiso",
        "ppg12_iso_threshold",
        "ppg12_noniso_threshold",
        "ppg12_xsec_pb",
        "ppg12_xsec_weight",
        "ppg12_window_low",
        "ppg12_window_high",
        "max_truth_jet_pt_r04",
        "truth_energy_contribution",
        "ppg12_response_Et",
        "ppg12_weight_slice",
        "ppg12_weight_vertex",
        "ppg12_weight_mix",
        "ppg12_weight_period",
        "ppg12_weight_final",
    )
    optional_int_branches = (
        "run",
        "cluster_index",
        "ppg12_is_iso",
        "ppg12_is_noniso",
        "ppg12_sample_bin",
        "ppg12_truth_window_pass_r04",
        "truth_barcode",
        "ppg12_logical_abcd_region",
        "ppg12_analysis_window_pass",
        "ppg12_response_window_pass",
        "ppg12_signal_fill_A",
        "ppg12_signal_fill_B",
        "ppg12_signal_fill_C",
        "ppg12_signal_fill_D",
        "ppg12_signal_fill_multiplicity",
        "ppg12_truth_class",
        "ppg12_weight_lane_code",
        "ppg12_weight_component_code",
    )

    rows: list[dict[str, Any]] = []
    events: set[tuple[int, int]] = set()
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        evt = int(getattr(t, "evt"))
        local_eventnumber = int(getattr(t, "eventnumber"))
        segment, eventnumber = recoil_source_identity(
            evt, local_eventnumber, events_per_segment
        )
        analysis_et = float(getattr(t, "cluster_Et"))
        eta = float(getattr(t, "cluster_Eta"))
        vertexz = float(getattr(t, "ppg12_kin_vertexz"))
        if not passes_ppg12_kinematics(analysis_et, eta, vertexz):
            continue
        score_input_et = (
            float(getattr(t, "cluster_Et_score_input"))
            if has_score_input_et
            else analysis_et
        )
        row = {
            "tree_entry": i,
            "segment": segment,
            "evt": evt,
            "eventnumber": eventnumber,
            "source_event_ordinal": evt,
            "local_eventnumber": local_eventnumber,
            "truth_track_id": int(getattr(t, "truth_track_id")),
            "is_signal": int(getattr(t, "is_signal")),
            "cluster_Et": analysis_et,
            "cluster_Et_score_input": score_input_et,
            "cluster_Eta": eta,
            "cluster_Phi": float(getattr(t, "cluster_Phi")),
            "vertexz": vertexz,
            "cluster_weta_cogx": float(getattr(t, "cluster_weta_cogx")),
            "cluster_wphi_cogx": float(getattr(t, "cluster_wphi_cogx")),
            "e11_over_e33": float(getattr(t, "e11_over_e33")),
            "e32_over_e35": float(getattr(t, "e32_over_e35")),
            "cluster_et1": float(getattr(t, "cluster_et1")),
            "cluster_et2": float(getattr(t, "cluster_et2")),
            "cluster_et3": float(getattr(t, "cluster_et3")),
            "cluster_et4": float(getattr(t, "cluster_et4")),
            "cluster_prob": float(getattr(t, "cluster_prob")),
            "npb_score": float(getattr(t, "npb_score")),
            "rj_stored_bdt_score": float(getattr(t, "tight_bdt_score")),
            "rj_stored_common_pass": int(getattr(t, "ppg12_common_pass")),
            "rj_stored_tag": int(getattr(t, "ppg12_tight_tag")),
            "score_input_fallback_used": int(not has_score_input_et),
        }
        for name in optional_float_branches:
            row[name] = float(getattr(t, name)) if name in branches else float("nan")
        for name in optional_int_branches:
            row[name] = int(getattr(t, name)) if name in branches else -1

        base_e_score = eval_rbdt(base_e_model, BASE_E_FEATURES, row)
        base_v3e_score = eval_rbdt(base_v3e_model, BASE_V3E_FEATURES, row)
        model, score = selected_score(
            {"base_E": base_e_score, "base_v3E": base_v3e_score}, analysis_et
        )
        tag, flags = classify(row, score, analysis_et)
        row.update(
            rj_bdt_base_E=base_e_score,
            rj_bdt_base_v3E=base_v3e_score,
            rj_selected_model=model,
            rj_selected_bdt_score=score,
            rj_inferred_stored_model=infer_stored_model_route(
                row["rj_stored_bdt_score"], base_e_score, base_v3e_score
            ),
            rj_recomputed_tag=tag,
            rj_recomputed_common_pass=int(bool(flags.get("common"))),
        )
        row["rj_stored_abcd_region"] = abcd_region(
            row["rj_stored_tag"], row["ppg12_is_iso"], row["ppg12_is_noniso"]
        )
        row["rj_recomputed_abcd_region"] = abcd_region(
            row["rj_recomputed_tag"], row["ppg12_is_iso"], row["ppg12_is_noniso"]
        )
        row["rj_candidate_identity"] = (
            f"seg{segment}:evt{eventnumber}:trk{row['truth_track_id']}:"
            f"cluster{row['cluster_index']}"
        )
        rows.append(row)
        events.add((segment, eventnumber))
    f.Close()
    return rows, events


def read_ppg12_rows(
    path: str,
    wanted_events: set[tuple[int, int]] | None,
    max_entries: int,
    mask_path: str,
    events_per_segment: int,
    *,
    wanted_identities: set[tuple[int, int, int]] | None = None,
    population_audit: dict[tuple[int, int, int], dict[str, Any]] | None = None,
) -> dict[tuple[int, int], list[dict[str, Any]]]:
    root = require_root()
    f = root.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open PPG12 ROOT: {path}")
    t = f.Get("slimtree")
    if not t:
        raise RuntimeError(f"Missing slimtree in {path}")
    branches = {branch.GetName() for branch in t.GetListOfBranches()}
    required = {
        "eventnumber",
        "vertexz",
        "ncluster_CLUSTERINFO_CEMC",
        "cluster_Et_CLUSTERINFO_CEMC",
        "cluster_Eta_CLUSTERINFO_CEMC",
        "cluster_Phi_CLUSTERINFO_CEMC",
        "cluster_prob_CLUSTERINFO_CEMC",
        "cluster_truthtrkID_CLUSTERINFO_CEMC",
        "cluster_ietacent_CLUSTERINFO_CEMC",
        "cluster_iphicent_CLUSTERINFO_CEMC",
        "cluster_weta_cogx_CLUSTERINFO_CEMC",
        "cluster_wphi_cogx_CLUSTERINFO_CEMC",
        "cluster_et1_CLUSTERINFO_CEMC",
        "cluster_et2_CLUSTERINFO_CEMC",
        "cluster_et3_CLUSTERINFO_CEMC",
        "cluster_et4_CLUSTERINFO_CEMC",
        "cluster_e11_CLUSTERINFO_CEMC",
        "cluster_e33_CLUSTERINFO_CEMC",
        "cluster_e32_CLUSTERINFO_CEMC",
        "cluster_e35_CLUSTERINFO_CEMC",
        "cluster_npb_score_CLUSTERINFO_CEMC",
        "cluster_iso_topo_04_CLUSTERINFO_CEMC",
        "cluster_bdt_CLUSTERINFO_CEMC_base_E",
        "cluster_bdt_CLUSTERINFO_CEMC_base_v3E",
        "nparticles",
        "particle_pid",
        "particle_trkid",
        "particle_Pt",
        "particle_Eta",
        "particle_photonclass",
        "particle_truth_iso_03",
    }
    missing = sorted(required - branches)
    if missing:
        raise RuntimeError(f"Missing required executable-PPG12 oracle branches: {missing}")

    mf = root.TFile.Open(mask_path) if mask_path else None
    if mask_path and (not mf or mf.IsZombie()):
        raise RuntimeError(f"Could not open nominal PPG12 tower mask: {mask_path}")
    mask = mf.Get("mask_phisymm_tight") if mf else None
    if mask_path and not mask:
        raise RuntimeError(f"Missing mask_phisymm_tight in {mask_path}")
    out: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    stop = t.GetEntries()
    if max_entries > 0:
        stop = min(stop, max_entries)

    wanted_tracks_by_event: dict[tuple[int, int], set[int]] = defaultdict(set)
    for segment, eventnumber, truth_track_id in wanted_identities or set():
        wanted_tracks_by_event[(segment, eventnumber)].add(truth_track_id)
        if population_audit is not None:
            population_audit[(segment, eventnumber, truth_track_id)] = {
                "event_found": 0,
                "vertex_pass": 0,
                "ppg12_vertexz": float("nan"),
                "truth_signal": 0,
                "truth_cluster_count": 0,
                "tower_masked_count": 0,
                "unmasked_cluster_count": 0,
                "et_pass_count": 0,
                "eta_pass_count": 0,
                "accepted_count": 0,
            }

    segment = 0
    previous_eventnumber: int | None = None
    for ie in range(stop):
        t.GetEntry(ie)
        eventnumber = int(getattr(t, "eventnumber"))
        segment = next_ppg12_segment(segment, previous_eventnumber, eventnumber)
        previous_eventnumber = eventnumber
        event_key = (segment, eventnumber)
        # ``None`` means the symmetric oracle universe: inspect every PPG12
        # event in the sealed source slice.  A Recoil-derived filter cannot be
        # used by an admission oracle because a whole event/candidate missing
        # from RecoilJets would otherwise be invisible.
        if not ppg12_event_in_scan(event_key, wanted_events):
            continue
        vertexz = float(getattr(t, "vertexz"))
        relevant_tracks = wanted_tracks_by_event.get(event_key, set())
        if population_audit is not None:
            for truth_track_id in relevant_tracks:
                audit = population_audit[(segment, eventnumber, truth_track_id)]
                audit["event_found"] = 1
                audit["ppg12_vertexz"] = vertexz
                audit["vertex_pass"] = int(
                    finite(vertexz) and abs(vertexz) <= PPG12_VERTEX_ABS_MAX
                )

        ncluster = min(int(getattr(t, "ncluster_CLUSTERINFO_CEMC")), 20000)
        cluster_Et = getattr(t, "cluster_Et_CLUSTERINFO_CEMC")
        cluster_Eta = getattr(t, "cluster_Eta_CLUSTERINFO_CEMC")
        cluster_Phi = getattr(t, "cluster_Phi_CLUSTERINFO_CEMC")
        cluster_prob = getattr(t, "cluster_prob_CLUSTERINFO_CEMC")
        cluster_truth = getattr(t, "cluster_truthtrkID_CLUSTERINFO_CEMC")
        cluster_ietacent = getattr(t, "cluster_ietacent_CLUSTERINFO_CEMC")
        cluster_iphicent = getattr(t, "cluster_iphicent_CLUSTERINFO_CEMC")
        cluster_weta = getattr(t, "cluster_weta_cogx_CLUSTERINFO_CEMC")
        cluster_wphi = getattr(t, "cluster_wphi_cogx_CLUSTERINFO_CEMC")
        cluster_et1 = getattr(t, "cluster_et1_CLUSTERINFO_CEMC")
        cluster_et2 = getattr(t, "cluster_et2_CLUSTERINFO_CEMC")
        cluster_et3 = getattr(t, "cluster_et3_CLUSTERINFO_CEMC")
        cluster_et4 = getattr(t, "cluster_et4_CLUSTERINFO_CEMC")
        cluster_e11 = getattr(t, "cluster_e11_CLUSTERINFO_CEMC")
        cluster_e33 = getattr(t, "cluster_e33_CLUSTERINFO_CEMC")
        cluster_e32 = getattr(t, "cluster_e32_CLUSTERINFO_CEMC")
        cluster_e35 = getattr(t, "cluster_e35_CLUSTERINFO_CEMC")
        cluster_npb = getattr(t, "cluster_npb_score_CLUSTERINFO_CEMC")
        cluster_raw_iso04 = getattr(t, "cluster_iso_topo_04_CLUSTERINFO_CEMC")
        bdt_base_e = getattr(t, "cluster_bdt_CLUSTERINFO_CEMC_base_E")
        bdt_base_v3e = getattr(t, "cluster_bdt_CLUSTERINFO_CEMC_base_v3E")

        # Build the exact signal truth-track map used by the executable:
        # direct/fragmentation photons with truth isolation below 4 GeV.
        nparticles = min(int(getattr(t, "nparticles")), 20000)
        particle_pid = getattr(t, "particle_pid")
        particle_trkid = getattr(t, "particle_trkid")
        particle_pt = getattr(t, "particle_Pt")
        particle_eta = getattr(t, "particle_Eta")
        particle_class = getattr(t, "particle_photonclass")
        particle_iso03 = getattr(t, "particle_truth_iso_03")
        truth_tracks: dict[int, dict[str, float | int]] = {}
        signal_tracks: dict[int, dict[str, float | int]] = {}
        for ip in range(nparticles):
            truth_track_id = int(particle_trkid[ip])
            metadata = {
                "truth_particle_index": ip,
                "truth_pt": float(particle_pt[ip]),
                "truth_eta": float(particle_eta[ip]),
                "truth_class": int(particle_class[ip]),
                "truth_iso03": float(particle_iso03[ip]),
            }
            truth_tracks[truth_track_id] = metadata
            if (
                int(particle_pid[ip]) == 22
                and int(particle_class[ip]) < 3
                and float(particle_iso03[ip]) < 4.0
            ):
                signal_tracks[truth_track_id] = metadata

        if population_audit is not None:
            for truth_track_id in relevant_tracks:
                population_audit[(segment, eventnumber, truth_track_id)]["truth_signal"] = int(
                    truth_track_id in signal_tracks
                )

        vertex_pass = finite(vertexz) and abs(vertexz) <= PPG12_VERTEX_ABS_MAX

        for ic in range(ncluster):
            truth_track = int(cluster_truth[ic])
            audit = (
                population_audit.get((segment, eventnumber, truth_track))
                if population_audit is not None and truth_track in relevant_tracks
                else None
            )
            if audit is not None:
                audit["truth_cluster_count"] += 1
            if tower_masked(mask, float(cluster_ietacent[ic]), float(cluster_iphicent[ic])):
                if audit is not None:
                    audit["tower_masked_count"] += 1
                continue
            if audit is not None:
                audit["unmasked_cluster_count"] += 1
            raw_et = float(cluster_Et[ic])
            eta = float(cluster_Eta[ic])
            if raw_et < PPG12_RECO_ET_MIN:
                continue
            if audit is not None:
                audit["et_pass_count"] += 1
            if not open_interval(eta, PPG12_ETA_MIN, PPG12_ETA_MAX):
                continue
            if audit is not None:
                audit["eta_pass_count"] += 1
            if not vertex_pass:
                continue
            e11 = float(cluster_e11[ic])
            e33 = float(cluster_e33[ic])
            e32 = float(cluster_e32[ic])
            e35 = float(cluster_e35[ic])
            scores = {
                "base_E": float(bdt_base_e[ic]),
                "base_v3E": float(bdt_base_v3e[ic]),
            }
            model, score = selected_score(scores, raw_et)
            iso = ppg12_iso_assignment(float(cluster_raw_iso04[ic]), raw_et)
            truth = truth_tracks.get(
                truth_track,
                {
                    "truth_particle_index": -1,
                    "truth_pt": float("nan"),
                    "truth_eta": float("nan"),
                    "truth_class": -1,
                    "truth_iso03": float("nan"),
                },
            )
            row = {
                "tree_entry": ie,
                "cluster_index": ic,
                "segment": segment,
                "eventnumber": eventnumber,
                "truth_track_id": truth_track,
                "cluster_Et": raw_et,
                "cluster_Et_score_input": raw_et,
                "cluster_Eta": eta,
                "cluster_Phi": float(cluster_Phi[ic]),
                "vertexz": vertexz,
                "cluster_weta_cogx": float(cluster_weta[ic]),
                "cluster_wphi_cogx": float(cluster_wphi[ic]),
                "e11_over_e33": cpp_float_ratio(e11, e33),
                "e32_over_e35": cpp_float_ratio(e32, e35),
                "cluster_et1": float(cluster_et1[ic]),
                "cluster_et2": float(cluster_et2[ic]),
                "cluster_et3": float(cluster_et3[ic]),
                "cluster_et4": float(cluster_et4[ic]),
                "cluster_prob": float(cluster_prob[ic]),
                "npb_score": float(cluster_npb[ic]),
                "bdt_base_E": float(bdt_base_e[ic]),
                "bdt_base_v3E": float(bdt_base_v3e[ic]),
                "selected_model": model,
                "selected_bdt_score": score,
                "is_signal": int(truth_track in signal_tracks),
                **truth,
                **iso,
            }
            tag, flags = classify(row, score, raw_et)
            row["ppg12_recomputed_tag"] = tag
            row["ppg12_common_pass"] = 1 if flags.get("common") else 0
            row["ppg12_abcd_region"] = abcd_region(
                tag, int(iso["is_iso"]), int(iso["is_noniso"])
            )
            row["analysis_window_pass"] = (
                in_ppg12_analysis_window(raw_et, float(row["truth_pt"]))
                if row["is_signal"]
                else 0
            )
            row["logical_abcd_region"] = abcd_region_code(
                str(row["ppg12_abcd_region"])
            )
            fill_flags = (
                ppg12_signal_fill_flags(
                    str(row["ppg12_abcd_region"]), int(row["analysis_window_pass"])
                )
                if row["is_signal"]
                else {name: 0 for name in "ABCD"}
            )
            for name, value in fill_flags.items():
                row[f"signal_fill_{name}"] = value
            row["fill_multiplicity"] = sum(fill_flags.values())
            row["ppg12_candidate_identity"] = (
                f"seg{segment}:evt{eventnumber}:trk{truth_track}:cluster{ic}"
            )
            # These values are a literal Python transcription of the
            # preserved selection code, not outputs written by the preserved
            # RecoEffCalculator executable.  Keep the provenance explicit so
            # a downstream gate cannot promote shadow agreement to
            # executable-stage agreement.
            row["tag_evidence_source"] = PYTHON_SHADOW_EVIDENCE
            row["isolation_abcd_evidence_source"] = PYTHON_SHADOW_EVIDENCE
            row["truth_response_fill_evidence_source"] = PYTHON_SHADOW_EVIDENCE
            row["weight_evidence_source"] = UNAVAILABLE_EXECUTABLE_EVIDENCE
            row["response_Et"] = ""
            row["response_window_pass"] = ""
            row["weight_final"] = ""
            out[event_key].append(row)
            if audit is not None:
                audit["accepted_count"] += 1
    if mf:
        mf.Close()
    f.Close()
    return out


def first_population_rejection(audit: dict[str, Any] | None) -> str:
    """Return the first preserved-PPG12 population gate not passed."""
    if not audit or not int(audit.get("event_found", 0)):
        return "ppg12_event_missing"
    if not int(audit.get("vertex_pass", 0)):
        return "ppg12_vertex_rejected"
    if not int(audit.get("truth_signal", 0)):
        return "ppg12_truth_signal_rejected"
    if int(audit.get("truth_cluster_count", 0)) == 0:
        return "ppg12_truth_matched_cluster_missing"
    if int(audit.get("unmasked_cluster_count", 0)) == 0:
        return "ppg12_tower_mask_rejected"
    if int(audit.get("et_pass_count", 0)) == 0:
        return "ppg12_reco_et_rejected"
    if int(audit.get("eta_pass_count", 0)) == 0:
        return "ppg12_reco_eta_rejected"
    if int(audit.get("accepted_count", 0)) == 0:
        return "ppg12_population_unknown"
    return "recoiljets_multiplicity_surplus"


def annotate_population_mismatches(
    matches: list[tuple[dict[str, Any], dict[str, Any], float]],
    rj_only: list[dict[str, Any]],
    ppg12_only: list[dict[str, Any]],
    population_audit: dict[tuple[int, int, int], dict[str, Any]],
) -> None:
    """Attach first-gate evidence without changing candidate matching."""
    for rj in rj_only:
        identity = (
            int(rj["segment"]),
            int(rj["eventnumber"]),
            int(rj["truth_track_id"]),
        )
        audit = population_audit.get(identity)
        rj["population_mismatch_reason"] = first_population_rejection(audit)
        for key, value in (audit or {}).items():
            rj[f"ppg12_audit_{key}"] = value
    for ppg in ppg12_only:
        ppg["population_mismatch_reason"] = "recoiljets_multiplicity_deficit"


def match_rows(
    rj_rows: list[dict[str, Any]],
    ppg12_by_event: dict[tuple[int, int], list[dict[str, Any]]],
) -> tuple[
    list[tuple[dict[str, Any], dict[str, Any], float]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    matches = []
    rj_only = []
    used: set[tuple[int, int, int]] = set()
    for rj in rj_rows:
        event_key = (int(rj["segment"]), int(rj["eventnumber"]))
        candidates = [
            (idx, p)
            for idx, p in enumerate(ppg12_by_event.get(event_key, []))
            if int(p["truth_track_id"]) == int(rj["truth_track_id"])
        ]
        best = None
        best_dist = float("inf")
        for idx, p in candidates:
            key = (event_key[0], event_key[1], idx)
            if key in used:
                continue
            deta = float(rj["cluster_Eta"]) - float(p["cluster_Eta"])
            dph = dphi(float(rj["cluster_Phi"]), float(p["cluster_Phi"]))
            dist = math.hypot(deta, dph)
            if dist < best_dist:
                best_dist = dist
                best = (idx, p)
        # Candidate identity is the stable event + truth-track key.  Delta-R
        # only chooses among duplicate reconstructed clusters with that key;
        # a shifted cluster is a geometry/feature divergence, not a missing
        # candidate.  A population mismatch is retained only when one side has
        # no unused row for the same stable identity.
        if best is None:
            rj_only.append(rj)
            continue
        used.add((event_key[0], event_key[1], best[0]))
        matches.append((rj, best[1], best_dist))
    ppg12_only = []
    for (segment, eventnumber), rows in ppg12_by_event.items():
        for idx, row in enumerate(rows):
            if (int(segment), int(eventnumber), idx) not in used:
                ppg12_only.append(row)
    return matches, rj_only, ppg12_only


def write_candidate_report(
    matches: list[tuple[dict[str, Any], dict[str, Any], float]],
    rj_only: list[dict[str, Any]],
    ppg12_only: list[dict[str, Any]],
    out_csv: Path,
    *,
    lane_id: str = "",
    runtime_contract_sha256: str = "",
) -> None:
    """Write the row-level oracle evidence, including population mismatches."""
    base_fields = [
        "lane_id",
        "runtime_contract_sha256",
        "match_status",
        "candidate_identity",
        "rj_candidate_identity",
        "ppg12_candidate_identity",
        "segment",
        "eventnumber",
        "truth_track_id",
        "identity_match",
        "match_dr",
        "match_dr_exceeds_002",
        "population_mismatch_reason",
        "ppg12_audit_event_found",
        "ppg12_audit_vertex_pass",
        "ppg12_audit_ppg12_vertexz",
        "ppg12_audit_truth_signal",
        "ppg12_audit_truth_cluster_count",
        "ppg12_audit_tower_masked_count",
        "ppg12_audit_unmasked_cluster_count",
        "ppg12_audit_et_pass_count",
        "ppg12_audit_eta_pass_count",
        "ppg12_audit_accepted_count",
        "rj_tree_entry",
        "rj_source_event_ordinal",
        "rj_local_eventnumber",
        "rj_cluster_index",
        "ppg12_tree_entry",
        "ppg12_chain_file_index",
        "ppg12_local_tree_entry",
        "ppg12_cluster_index",
        "truth_particle_index",
        "rj_is_signal",
        "ppg12_is_signal",
        "signal_status_agree",
        "truth_class",
        "rj_truth_class",
        "truth_class_agree",
        "truth_pt",
        "truth_eta",
        "truth_iso03",
        "rj_cluster_Et",
        "rj_cluster_Et_score_input",
        "rj_response_Et",
        "ppg12_cluster_Et",
        "cluster_Et_delta",
        "rj_cluster_Eta",
        "ppg12_cluster_Eta",
        "cluster_Eta_delta",
        "rj_cluster_Phi",
        "ppg12_cluster_Phi",
        "cluster_Phi_delta",
        "rj_selected_model",
        "rj_inferred_stored_model",
        "ppg12_selected_model",
        "model_route_agree",
        "ppg12_tag_evidence_source",
        "ppg12_isolation_abcd_evidence_source",
        "ppg12_truth_response_fill_evidence_source",
        "ppg12_weight_evidence_source",
        "rj_bdt_base_E",
        "ppg12_bdt_base_E",
        "base_E_score_delta",
        "rj_bdt_base_v3E",
        "ppg12_bdt_base_v3E",
        "base_v3E_score_delta",
        "rj_selected_bdt_score",
        "ppg12_selected_bdt_score",
        "selected_score_delta",
        "rj_stored_bdt_score",
        "stored_minus_routed_score",
        "rj_recomputed_common_pass",
        "ppg12_common_pass",
        "common_agree",
        "stored_common_agree",
        "rj_recomputed_tag",
        "ppg12_recomputed_tag",
        "rj_tight_flag",
        "ppg12_tight_flag",
        "rj_nontight_flag",
        "ppg12_nontight_flag",
        "rj_neither_flag",
        "ppg12_neither_flag",
        "rj_preselection_fail_flag",
        "ppg12_preselection_fail_flag",
        "rj_isolation_gap_flag",
        "ppg12_isolation_gap_flag",
        "tag_agree",
        "stored_tag_agree",
        "rj_stored_common_pass",
        "rj_stored_tag",
        "rj_raw_eiso",
        "ppg12_raw_eiso",
        "rj_corrected_eiso",
        "ppg12_corrected_eiso",
        "rj_iso_threshold",
        "ppg12_iso_threshold",
        "rj_noniso_threshold",
        "ppg12_noniso_threshold",
        "rj_is_iso",
        "ppg12_is_iso",
        "rj_is_noniso",
        "ppg12_is_noniso",
        "rj_abcd_region",
        "rj_stored_abcd_region",
        "ppg12_abcd_region",
        "abcd_agree",
        "stored_abcd_agree",
        "rj_event_weight",
        "rj_xsec_pb",
        "rj_xsec_weight",
        "rj_sample_bin",
        "rj_window_low",
        "rj_window_high",
        "rj_truth_window_pass_r04",
        "rj_truth_barcode",
        "rj_truth_energy_contribution",
        "ppg12_fill_multiplicity",
        "rj_logical_abcd_region",
        "ppg12_logical_abcd_region",
        "rj_analysis_window_pass",
        "ppg12_analysis_window_pass",
        "rj_response_window_pass",
        "ppg12_response_Et",
        "ppg12_response_window_pass",
        "rj_signal_fill_A",
        "rj_signal_fill_B",
        "rj_signal_fill_C",
        "rj_signal_fill_D",
        "ppg12_signal_fill_A",
        "ppg12_signal_fill_B",
        "ppg12_signal_fill_C",
        "ppg12_signal_fill_D",
        "rj_signal_fill_multiplicity",
        "rj_weight_lane_code",
        "rj_weight_component_code",
        "rj_weight_slice",
        "rj_weight_vertex",
        "rj_weight_mix",
        "rj_weight_period",
        "rj_weight_final",
        "ppg12_weight_sample",
        "ppg12_weight_mix",
        "ppg12_weight_lumi",
        "ppg12_weight_cross",
        "ppg12_weight_vertex",
        "ppg12_weight_truth_vertex",
        "ppg12_weight_trigger",
        "ppg12_weight_event",
        "ppg12_weight_final",
        "rj_weight_factor_product",
        "rj_weight_product_delta",
        "rj_event_weight_delta",
        "rj_score_input_fallback_used",
    ]
    feature_fields: list[str] = []
    for rj_name, _ in FEATURES:
        feature_fields.extend((f"{rj_name}_rj", f"{rj_name}_ppg12", f"{rj_name}_delta"))

    def difference(left: Any, right: Any) -> float | str:
        if finite(left) and finite(right):
            return float(left) - float(right)
        return ""

    def make_record(
        status: str,
        rj: dict[str, Any] | None,
        ppg: dict[str, Any] | None,
        dist: float | None,
    ) -> dict[str, Any]:
        rj = rj or {}
        ppg = ppg or {}
        segment = rj.get("segment", ppg.get("segment", ""))
        eventnumber = rj.get("eventnumber", ppg.get("eventnumber", ""))
        truth_track_id = rj.get("truth_track_id", ppg.get("truth_track_id", ""))
        rj_cluster = rj.get("cluster_index", "missing")
        ppg_cluster = ppg.get("cluster_index", "missing")
        candidate_identity = (
            f"seg{segment}:evt{eventnumber}:trk{truth_track_id}:"
            f"rjentry{rj.get('tree_entry', 'missing')}:rjcluster{rj_cluster}:"
            f"ppgentry{ppg.get('tree_entry', 'missing')}:ppgcluster{ppg_cluster}"
        )
        weight_factors = (
            rj.get("ppg12_weight_slice"),
            rj.get("ppg12_weight_vertex"),
            rj.get("ppg12_weight_mix"),
            rj.get("ppg12_weight_period"),
        )
        weight_product: float | str = ""
        if all(finite(value) for value in weight_factors):
            weight_product = math.prod(float(value) for value in weight_factors)
        record: dict[str, Any] = {
            "lane_id": lane_id,
            "runtime_contract_sha256": runtime_contract_sha256,
            "match_status": status,
            "candidate_identity": candidate_identity,
            "rj_candidate_identity": rj.get("rj_candidate_identity", ""),
            "ppg12_candidate_identity": ppg.get("ppg12_candidate_identity", ""),
            "segment": segment,
            "eventnumber": eventnumber,
            "truth_track_id": truth_track_id,
            "identity_match": int(bool(rj and ppg)),
            "match_dr": dist if dist is not None else "",
            "match_dr_exceeds_002": int(
                dist is not None and dist > GEOMETRY_DIAGNOSTIC_DR
            ),
            "population_mismatch_reason": rj.get(
                "population_mismatch_reason", ppg.get("population_mismatch_reason", "")
            ),
            "rj_tree_entry": rj.get("tree_entry", ""),
            "rj_source_event_ordinal": rj.get("source_event_ordinal", ""),
            "rj_local_eventnumber": rj.get("local_eventnumber", ""),
            "rj_cluster_index": rj.get("cluster_index", ""),
            "ppg12_tree_entry": ppg.get("tree_entry", ""),
            "ppg12_chain_file_index": ppg.get("chain_file_index", ""),
            "ppg12_local_tree_entry": ppg.get("local_tree_entry", ""),
            "ppg12_cluster_index": ppg.get("cluster_index", ""),
            "truth_particle_index": ppg.get("truth_particle_index", ""),
            "rj_is_signal": rj.get("is_signal", ""),
            "ppg12_is_signal": ppg.get("is_signal", ""),
            "signal_status_agree": int(
                bool(rj and ppg) and rj.get("is_signal") == ppg.get("is_signal")
            ),
            "truth_class": ppg.get("truth_class", ""),
            "rj_truth_class": rj.get("ppg12_truth_class", ""),
            "truth_class_agree": int(
                bool(rj and ppg)
                and (
                    (
                        ppg.get("is_signal") == 1
                        and rj.get("ppg12_truth_class") == ppg.get("truth_class")
                    )
                    or (ppg.get("is_signal") == 0 and rj.get("is_signal") == 0)
                )
            ),
            "truth_pt": ppg.get("truth_pt", ""),
            "truth_eta": ppg.get("truth_eta", ""),
            "truth_iso03": ppg.get("truth_iso03", ""),
            "rj_cluster_Et": rj.get("cluster_Et", ""),
            "rj_cluster_Et_score_input": rj.get("cluster_Et_score_input", ""),
            "rj_response_Et": rj.get("ppg12_response_Et", ""),
            "ppg12_cluster_Et": ppg.get("cluster_Et", ""),
            "cluster_Et_delta": difference(rj.get("cluster_Et"), ppg.get("cluster_Et")),
            "rj_cluster_Eta": rj.get("cluster_Eta", ""),
            "ppg12_cluster_Eta": ppg.get("cluster_Eta", ""),
            "cluster_Eta_delta": difference(
                rj.get("cluster_Eta"), ppg.get("cluster_Eta")
            ),
            "rj_cluster_Phi": rj.get("cluster_Phi", ""),
            "ppg12_cluster_Phi": ppg.get("cluster_Phi", ""),
            "cluster_Phi_delta": difference(
                rj.get("cluster_Phi"), ppg.get("cluster_Phi")
            ),
            "rj_selected_model": rj.get("rj_selected_model", ""),
            "rj_inferred_stored_model": rj.get("rj_inferred_stored_model", ""),
            "ppg12_selected_model": ppg.get("selected_model", ""),
            "model_route_agree": int(
                bool(rj and ppg)
                and rj.get("rj_selected_model") == ppg.get("selected_model")
            ),
            "ppg12_tag_evidence_source": ppg.get(
                "tag_evidence_source", UNAVAILABLE_EXECUTABLE_EVIDENCE
            ),
            "ppg12_isolation_abcd_evidence_source": ppg.get(
                "isolation_abcd_evidence_source", UNAVAILABLE_EXECUTABLE_EVIDENCE
            ),
            "ppg12_truth_response_fill_evidence_source": ppg.get(
                "truth_response_fill_evidence_source", UNAVAILABLE_EXECUTABLE_EVIDENCE
            ),
            "ppg12_weight_evidence_source": ppg.get(
                "weight_evidence_source", UNAVAILABLE_EXECUTABLE_EVIDENCE
            ),
            "rj_bdt_base_E": rj.get("rj_bdt_base_E", ""),
            "ppg12_bdt_base_E": ppg.get("bdt_base_E", ""),
            "base_E_score_delta": difference(rj.get("rj_bdt_base_E"), ppg.get("bdt_base_E")),
            "rj_bdt_base_v3E": rj.get("rj_bdt_base_v3E", ""),
            "ppg12_bdt_base_v3E": ppg.get("bdt_base_v3E", ""),
            "base_v3E_score_delta": difference(
                rj.get("rj_bdt_base_v3E"), ppg.get("bdt_base_v3E")
            ),
            "rj_selected_bdt_score": rj.get("rj_selected_bdt_score", ""),
            "ppg12_selected_bdt_score": ppg.get("selected_bdt_score", ""),
            "selected_score_delta": difference(
                rj.get("rj_selected_bdt_score"), ppg.get("selected_bdt_score")
            ),
            "rj_stored_bdt_score": rj.get("rj_stored_bdt_score", ""),
            "stored_minus_routed_score": difference(
                rj.get("rj_stored_bdt_score"), rj.get("rj_selected_bdt_score")
            ),
            "rj_recomputed_common_pass": rj.get("rj_recomputed_common_pass", ""),
            "ppg12_common_pass": ppg.get("ppg12_common_pass", ""),
            "common_agree": int(
                bool(rj and ppg)
                and rj.get("rj_recomputed_common_pass") == ppg.get("ppg12_common_pass")
            ),
            "stored_common_agree": int(
                bool(rj and ppg)
                and rj.get("rj_stored_common_pass") == ppg.get("ppg12_common_pass")
            ),
            "rj_recomputed_tag": rj.get("rj_recomputed_tag", ""),
            "ppg12_recomputed_tag": ppg.get("ppg12_recomputed_tag", ""),
            "rj_tight_flag": int(rj.get("rj_recomputed_tag") == 1),
            "ppg12_tight_flag": int(ppg.get("ppg12_recomputed_tag") == 1),
            "rj_nontight_flag": int(rj.get("rj_recomputed_tag") == 2),
            "ppg12_nontight_flag": int(ppg.get("ppg12_recomputed_tag") == 2),
            "rj_neither_flag": int(rj.get("rj_recomputed_tag") == 3),
            "ppg12_neither_flag": int(ppg.get("ppg12_recomputed_tag") == 3),
            "rj_preselection_fail_flag": int(rj.get("rj_recomputed_tag") == 0),
            "ppg12_preselection_fail_flag": int(ppg.get("ppg12_recomputed_tag") == 0),
            "rj_isolation_gap_flag": int(
                bool(rj)
                and rj.get("ppg12_is_iso") == 0
                and rj.get("ppg12_is_noniso") == 0
            ),
            "ppg12_isolation_gap_flag": int(
                bool(ppg)
                and ppg.get("is_iso") == 0
                and ppg.get("is_noniso") == 0
            ),
            "tag_agree": int(
                bool(rj and ppg)
                and rj.get("rj_recomputed_tag") == ppg.get("ppg12_recomputed_tag")
            ),
            "stored_tag_agree": int(
                bool(rj and ppg)
                and rj.get("rj_stored_tag") == ppg.get("ppg12_recomputed_tag")
            ),
            "rj_stored_common_pass": rj.get("rj_stored_common_pass", ""),
            "rj_stored_tag": rj.get("rj_stored_tag", ""),
            "rj_raw_eiso": rj.get("ppg12_raw_eiso", ""),
            "ppg12_raw_eiso": ppg.get("raw_eiso", ""),
            "rj_corrected_eiso": rj.get("ppg12_reco_eiso", ""),
            "ppg12_corrected_eiso": ppg.get("corrected_eiso", ""),
            "rj_iso_threshold": rj.get("ppg12_iso_threshold", ""),
            "ppg12_iso_threshold": ppg.get("iso_threshold", ""),
            "rj_noniso_threshold": rj.get("ppg12_noniso_threshold", ""),
            "ppg12_noniso_threshold": ppg.get("noniso_threshold", ""),
            "rj_is_iso": rj.get("ppg12_is_iso", ""),
            "ppg12_is_iso": ppg.get("is_iso", ""),
            "rj_is_noniso": rj.get("ppg12_is_noniso", ""),
            "ppg12_is_noniso": ppg.get("is_noniso", ""),
            "rj_abcd_region": rj.get("rj_recomputed_abcd_region", ""),
            "rj_stored_abcd_region": rj.get("rj_stored_abcd_region", ""),
            "ppg12_abcd_region": ppg.get("ppg12_abcd_region", ""),
            "abcd_agree": int(
                bool(rj and ppg)
                and rj.get("rj_recomputed_abcd_region") == ppg.get("ppg12_abcd_region")
            ),
            "stored_abcd_agree": int(
                bool(rj and ppg)
                and rj.get("rj_stored_abcd_region") == ppg.get("ppg12_abcd_region")
            ),
            "rj_event_weight": rj.get("event_weight", ""),
            "rj_xsec_pb": rj.get("ppg12_xsec_pb", ""),
            "rj_xsec_weight": rj.get("ppg12_xsec_weight", ""),
            "rj_sample_bin": rj.get("ppg12_sample_bin", ""),
            "rj_window_low": rj.get("ppg12_window_low", ""),
            "rj_window_high": rj.get("ppg12_window_high", ""),
            "rj_truth_window_pass_r04": rj.get("ppg12_truth_window_pass_r04", ""),
            "rj_truth_barcode": rj.get("truth_barcode", ""),
            "rj_truth_energy_contribution": rj.get("truth_energy_contribution", ""),
            "ppg12_fill_multiplicity": ppg.get("fill_multiplicity", ""),
            "rj_logical_abcd_region": rj.get("ppg12_logical_abcd_region", ""),
            "ppg12_logical_abcd_region": ppg.get("logical_abcd_region", ""),
            "rj_analysis_window_pass": rj.get("ppg12_analysis_window_pass", ""),
            "ppg12_analysis_window_pass": ppg.get("analysis_window_pass", ""),
            "rj_response_window_pass": rj.get("ppg12_response_window_pass", ""),
            "ppg12_response_Et": ppg.get("response_Et", ""),
            "ppg12_response_window_pass": ppg.get("response_window_pass", ""),
            "rj_signal_fill_A": rj.get("ppg12_signal_fill_A", ""),
            "rj_signal_fill_B": rj.get("ppg12_signal_fill_B", ""),
            "rj_signal_fill_C": rj.get("ppg12_signal_fill_C", ""),
            "rj_signal_fill_D": rj.get("ppg12_signal_fill_D", ""),
            "ppg12_signal_fill_A": ppg.get("signal_fill_A", ""),
            "ppg12_signal_fill_B": ppg.get("signal_fill_B", ""),
            "ppg12_signal_fill_C": ppg.get("signal_fill_C", ""),
            "ppg12_signal_fill_D": ppg.get("signal_fill_D", ""),
            "rj_signal_fill_multiplicity": rj.get(
                "ppg12_signal_fill_multiplicity", ""
            ),
            "rj_weight_lane_code": rj.get("ppg12_weight_lane_code", ""),
            "rj_weight_component_code": rj.get("ppg12_weight_component_code", ""),
            "rj_weight_slice": rj.get("ppg12_weight_slice", ""),
            "rj_weight_vertex": rj.get("ppg12_weight_vertex", ""),
            "rj_weight_mix": rj.get("ppg12_weight_mix", ""),
            "rj_weight_period": rj.get("ppg12_weight_period", ""),
            "rj_weight_final": rj.get("ppg12_weight_final", ""),
            "ppg12_weight_sample": ppg.get("weight_sample", ""),
            "ppg12_weight_mix": ppg.get("weight_mix", ""),
            "ppg12_weight_lumi": ppg.get("weight_lumi", ""),
            "ppg12_weight_cross": ppg.get("weight_cross", ""),
            "ppg12_weight_vertex": ppg.get("weight_vertex", ""),
            "ppg12_weight_truth_vertex": ppg.get("weight_truth_vertex", ""),
            "ppg12_weight_trigger": ppg.get("weight_trigger", ""),
            "ppg12_weight_event": ppg.get("weight_event", ""),
            "ppg12_weight_final": ppg.get("weight_final", ""),
            "rj_weight_factor_product": weight_product,
            "rj_weight_product_delta": difference(
                rj.get("ppg12_weight_final"), weight_product
            ),
            "rj_event_weight_delta": difference(
                rj.get("event_weight"), rj.get("ppg12_weight_final")
            ),
            "rj_score_input_fallback_used": rj.get("score_input_fallback_used", ""),
        }
        for key in (
            "event_found",
            "vertex_pass",
            "ppg12_vertexz",
            "truth_signal",
            "truth_cluster_count",
            "tower_masked_count",
            "unmasked_cluster_count",
            "et_pass_count",
            "eta_pass_count",
            "accepted_count",
        ):
            record[f"ppg12_audit_{key}"] = rj.get(f"ppg12_audit_{key}", "")
        for rj_name, ppg_name in FEATURES:
            rj_value = rj.get(rj_name, "")
            ppg_value = ppg.get(ppg_name, "")
            record[f"{rj_name}_rj"] = rj_value
            record[f"{rj_name}_ppg12"] = ppg_value
            record[f"{rj_name}_delta"] = difference(rj_value, ppg_value)
        return record

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=base_fields + feature_fields)
        writer.writeheader()
        for rj, ppg, dist in matches:
            writer.writerow(make_record("matched", rj, ppg, dist))
        for rj in rj_only:
            writer.writerow(make_record("rj_only", rj, None, None))
        for ppg in ppg12_only:
            writer.writerow(make_record("ppg12_only", None, ppg, None))


def write_reports(
    matches: list[tuple[dict[str, Any], dict[str, Any], float]],
    rj_only: list[dict[str, Any]],
    ppg12_only: list[dict[str, Any]],
    out_md: Path,
    out_csv: Path,
    out_candidates_csv: Path,
    base_e_model_path: str,
    base_v3e_model_path: str,
    *,
    lane_id: str = "",
    runtime_contract_sha256: str = "",
) -> None:
    aggs = [Agg() for _ in range(len(RECO_BINS) - 1)]
    ppg12_only_by_bin = [
        {"total": 0, "common": 0, "tight": 0, "nontight": 0, "neither": 0, "preselection_fail": 0}
        for _ in range(len(RECO_BINS) - 1)
    ]
    examples = []
    for rj, ppg, dist in matches:
        ib = find_bin(float(ppg["cluster_Et"]))
        if ib < 0:
            continue
        agg = aggs[ib]
        agg.n += 1
        agg.rj_common += int(rj["rj_recomputed_common_pass"] == 1)
        agg.ppg12_common += int(ppg["ppg12_common_pass"] == 1)
        agg.rj_tight += int(rj["rj_recomputed_tag"] == 1)
        agg.rj_nontight += int(rj["rj_recomputed_tag"] == 2)
        agg.ppg12_tight += int(ppg["ppg12_recomputed_tag"] == 1)
        agg.ppg12_nontight += int(ppg["ppg12_recomputed_tag"] == 2)
        agg.model_route_agree += int(rj["rj_selected_model"] == ppg["selected_model"])
        agg.tag_agree += int(rj["rj_recomputed_tag"] == ppg["ppg12_recomputed_tag"])
        agg.common_agree += int(
            rj["rj_recomputed_common_pass"] == ppg["ppg12_common_pass"]
        )
        agg.abcd_agree += int(
            rj["rj_recomputed_abcd_region"] == ppg["ppg12_abcd_region"]
        )
        agg.stored_tag_disagree += int(rj["rj_stored_tag"] != rj["rj_recomputed_tag"])
        agg.stored_common_disagree += int(
            rj["rj_stored_common_pass"] != rj["rj_recomputed_common_pass"]
        )
        agg.model_base_e += int(ppg["selected_model"] == "base_E")
        agg.model_base_v3e += int(ppg["selected_model"] == "base_v3E")
        score_diff = float(rj["rj_selected_bdt_score"]) - float(ppg["selected_bdt_score"])
        agg.score_diff_sum += score_diff
        agg.score_abs_sum += abs(score_diff)
        agg.rj_score_sum += float(rj["rj_selected_bdt_score"])
        agg.ppg12_score_sum += float(ppg["selected_bdt_score"])
        for rj_name, ppg_name in FEATURES:
            rj_val = float(rj.get(rj_name, float("nan")))
            agg.add_feature(rj_name, rj_val, float(ppg[ppg_name]))
        if len(examples) < 30 or abs(score_diff) > min(abs(e[0]) for e in examples):
            examples.append((score_diff, dist, rj, ppg))
            examples = sorted(examples, key=lambda x: abs(x[0]), reverse=True)[:30]

    for ppg in ppg12_only:
        ib = find_bin(float(ppg["cluster_Et"]))
        if ib < 0:
            continue
        bucket = ppg12_only_by_bin[ib]
        bucket["total"] += 1
        bucket["common"] += int(ppg["ppg12_common_pass"] == 1)
        tag = int(ppg["ppg12_recomputed_tag"])
        if tag == 0:
            bucket["preselection_fail"] += 1
        elif tag == 1:
            bucket["tight"] += 1
        elif tag == 2:
            bucket["nontight"] += 1
        else:
            bucket["neither"] += 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "bin", "matches", "rj_common_frac", "ppg12_common_frac",
            "rj_nontight_common", "ppg12_nontight_common",
            "rj_tight_common", "ppg12_tight_common",
            "mean_rj_score", "mean_ppg12_score", "mean_score_diff", "mean_abs_score_diff",
            "model_route_agreement", "tag_agreement", "common_agreement", "abcd_agreement",
            "stored_tag_disagreements", "stored_common_disagreements",
            "base_v3E_rows", "base_E_rows",
        ])
        for i, agg in enumerate(aggs):
            writer.writerow([
                bin_label(i), agg.n,
                safe_div(agg.rj_common, agg.n), safe_div(agg.ppg12_common, agg.n),
                safe_div(agg.rj_nontight, agg.rj_common), safe_div(agg.ppg12_nontight, agg.ppg12_common),
                safe_div(agg.rj_tight, agg.rj_common), safe_div(agg.ppg12_tight, agg.ppg12_common),
                safe_div(agg.rj_score_sum, agg.n), safe_div(agg.ppg12_score_sum, agg.n),
                safe_div(agg.score_diff_sum, agg.n), safe_div(agg.score_abs_sum, agg.n),
                safe_div(agg.model_route_agree, agg.n), safe_div(agg.tag_agree, agg.n),
                safe_div(agg.common_agree, agg.n), safe_div(agg.abcd_agree, agg.n),
                agg.stored_tag_disagree, agg.stored_common_disagree,
                agg.model_base_v3e, agg.model_base_e,
            ])

    write_candidate_report(
        matches,
        rj_only,
        ppg12_only,
        out_candidates_csv,
        lane_id=lane_id,
        runtime_contract_sha256=runtime_contract_sha256,
    )

    mismatch_reasons = Counter(
        str(row.get("population_mismatch_reason", "unspecified"))
        for row in rj_only + ppg12_only
    )

    lines = [
        "# Same-cluster executable PPG12/RecoilJets oracle",
        "",
        f"- Matched rows: {len(matches)}",
        f"- Unmatched RecoilJets signal rows: {len(rj_only)}",
        f"- PPG12-only signal rows across the complete sealed event universe: {len(ppg12_only)}",
        f"- CSV: `{out_csv}`",
        f"- Candidate-level CSV: `{out_candidates_csv}`",
        f"- PPG12 base_E model: `{base_e_model_path}`",
        f"- PPG12 baseV3E model: `{base_v3e_model_path}`",
        "- Classification contract: calibrated unsmeared ET for model routing, "
        "thresholds, tags, ABCD assignment, and report bins.",
        "- Candidate identity contract: segment + eventnumber + truth track id. "
        "Nearest delta-R disambiguates duplicate rows only; delta-R above 0.02 is "
        "reported as a reconstruction displacement, not a population mismatch.",
        "- The executable's additive response-matrix smearing is intentionally "
        "outside this classification oracle.",
        "- RecoilJets row weights and lane metadata are exported candidate by candidate. "
        "The PPG12 stitched event weight is supplied outside the slimtree and must be "
        "checked by the 32-lane closure gate rather than inferred here.",
        "- Population mismatch reasons: "
        + (", ".join(f"{key}={value}" for key, value in sorted(mismatch_reasons.items()))
           if mismatch_reasons else "none"),
        "",
        "## Stage Summary",
        "",
        "| reco ET bin | matches | RJ common | PPG12 common | RJ NT/common | PPG12 NT/common | RJ tight/common | PPG12 tight/common | mean RJ BDT | mean PPG12 BDT | mean score diff | mean abs score diff | route agree | tag agree | ABCD agree | stored tag drift | base_v3E | base_E |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for i, agg in enumerate(aggs):
        lines.append(
            f"| {bin_label(i)} | {agg.n} | {safe_div(agg.rj_common, agg.n):.4f} | "
            f"{safe_div(agg.ppg12_common, agg.n):.4f} | "
            f"{safe_div(agg.rj_nontight, agg.rj_common):.4f} | "
            f"{safe_div(agg.ppg12_nontight, agg.ppg12_common):.4f} | "
            f"{safe_div(agg.rj_tight, agg.rj_common):.4f} | "
            f"{safe_div(agg.ppg12_tight, agg.ppg12_common):.4f} | "
            f"{safe_div(agg.rj_score_sum, agg.n):.4f} | "
            f"{safe_div(agg.ppg12_score_sum, agg.n):.4f} | "
            f"{safe_div(agg.score_diff_sum, agg.n):+.4f} | "
            f"{safe_div(agg.score_abs_sum, agg.n):.4f} | "
            f"{safe_div(agg.model_route_agree, agg.n):.4f} | "
            f"{safe_div(agg.tag_agree, agg.n):.4f} | "
            f"{safe_div(agg.abcd_agree, agg.n):.4f} | "
            f"{agg.stored_tag_disagree} | "
            f"{agg.model_base_v3e} | {agg.model_base_e} |"
        )

    lines += [
        "",
        "## PPG12-Only Signal Candidates In Same Events",
        "",
        "| reco ET bin | total | common | tight | nontight | neither | preselection fail |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for i, bucket in enumerate(ppg12_only_by_bin):
        lines.append(
            f"| {bin_label(i)} | {bucket['total']} | {bucket['common']} | "
            f"{bucket['tight']} | {bucket['nontight']} | {bucket['neither']} | "
            f"{bucket['preselection_fail']} |"
        )
    lines += ["", "## Mean Absolute Feature Differences", ""]
    for i, agg in enumerate(aggs):
        if agg.n == 0:
            continue
        lines.append(f"### {bin_label(i)}")
        lines.append("")
        lines.append("| feature | mean signed RJ-PPG12 | mean abs |")
        lines.append("| --- | ---: | ---: |")
        for name, _ in FEATURES:
            n = agg.feature_n.get(name, 0)
            lines.append(
                f"| {name} | {safe_div(agg.feature_signed.get(name, 0.0), n):+.6g} | "
                f"{safe_div(agg.feature_abs.get(name, 0.0), n):.6g} |"
            )
        lines.append("")

    lines += ["## Largest Score-Shift Examples", "", "| score diff | dR | segment | event | truth_track | RJ tag | PPG12 tag | RJ score | PPG12 score | PPG12 ET | RJ ET | RJ score-input ET | weta RJ/PPG12 | e11e33 RJ/PPG12 | et2 RJ/PPG12 | et3 RJ/PPG12 | et4 RJ/PPG12 | model |", "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- |"]
    for score_diff, dist, rj, ppg in examples:
        lines.append(
            f"| {score_diff:+.5f} | {dist:.5g} | {int(rj['segment'])} | {int(rj['eventnumber'])} | {int(rj['truth_track_id'])} | "
            f"{int(rj['rj_recomputed_tag'])} | {int(ppg['ppg12_recomputed_tag'])} | "
            f"{float(rj['rj_selected_bdt_score']):.5f} | {float(ppg['selected_bdt_score']):.5f} | "
            f"{float(ppg['cluster_Et']):.4f} | {float(rj['cluster_Et']):.4f} | {float(rj['cluster_Et_score_input']):.4f} | "
            f"{float(rj['cluster_weta_cogx']):.4f}/{float(ppg['cluster_weta_cogx']):.4f} | "
            f"{float(rj['e11_over_e33']):.4f}/{float(ppg['e11_over_e33']):.4f} | "
            f"{float(rj['cluster_et2']):.4f}/{float(ppg['cluster_et2']):.4f} | "
            f"{float(rj['cluster_et3']):.4f}/{float(ppg['cluster_et3']):.4f} | "
            f"{float(rj['cluster_et4']):.4f}/{float(ppg['cluster_et4']):.4f} | {ppg['selected_model']} |"
        )
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rj-root", required=True)
    parser.add_argument("--ppg12-root", default="/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon20/bdt_split.root")
    parser.add_argument("--mask-root", default="/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root")
    parser.add_argument("--ppg12-max-events", type=int, default=5000)
    parser.add_argument("--events-per-segment", type=int, default=1000)
    parser.add_argument("--base-e-model", default=DEFAULT_BASE_E_MODEL)
    parser.add_argument("--base-v3e-model", default=DEFAULT_BASE_V3E_MODEL)
    parser.add_argument(
        "--ppg12-executable-trace",
        type=Path,
        help="Instrumented preserved-RecoEff candidate trace CSV",
    )
    parser.add_argument(
        "--ppg12-executable-response-trace",
        type=Path,
        help="Instrumented preserved-RecoEff response-smear trace CSV",
    )
    parser.add_argument(
        "--allow-legacy-score-et-fallback",
        action="store_true",
        help=(
            "Use cluster_Et as the model-input ET only when the exact "
            "cluster_Et_score_input branch is absent; labels the output as legacy."
        ),
    )
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument(
        "--out-candidates-csv",
        help="Candidate-level evidence CSV (default: <out-csv stem>_candidates.csv)",
    )
    parser.add_argument(
        "--lane-id",
        help="Physical lane identity embedded in every executable-evidence row",
    )
    parser.add_argument(
        "--runtime-contract",
        type=Path,
        help="Exact paired-oracle contract whose SHA-256 is embedded in every row",
    )
    args = parser.parse_args()

    base_e_model = load_rbdt(args.base_e_model)
    base_v3e_model = load_rbdt(args.base_v3e_model)
    rj_rows, _rj_events = read_rj_rows(
        args.rj_root,
        args.events_per_segment,
        base_e_model,
        base_v3e_model,
        args.allow_legacy_score_et_fallback,
    )
    ppg12_by_event = read_ppg12_rows(
        args.ppg12_root,
        None,
        args.ppg12_max_events,
        args.mask_root,
        args.events_per_segment,
        wanted_identities=None,
        population_audit=(population_audit := {}),
    )
    if bool(args.ppg12_executable_trace) != bool(args.ppg12_executable_response_trace):
        raise RuntimeError(
            "candidate and response executable traces must be supplied together"
        )
    if bool(args.lane_id) != bool(args.runtime_contract):
        raise RuntimeError("--lane-id and --runtime-contract must be supplied together")
    runtime_contract_sha256 = ""
    if args.runtime_contract:
        if not args.runtime_contract.is_file() or args.runtime_contract.stat().st_size <= 0:
            raise RuntimeError(f"runtime contract is missing or empty: {args.runtime_contract}")
        runtime_contract_sha256 = hashlib.sha256(
            args.runtime_contract.read_bytes()
        ).hexdigest()
    if args.ppg12_executable_trace:
        bind_executable_trace(
            ppg12_by_event,
            read_executable_trace(
                args.ppg12_executable_trace,
                args.ppg12_executable_response_trace,
            ),
        )
    matches, rj_only, ppg12_only = match_rows(rj_rows, ppg12_by_event)
    annotate_population_mismatches(matches, rj_only, ppg12_only, population_audit)
    out_csv = Path(args.out_csv)
    out_candidates_csv = (
        Path(args.out_candidates_csv)
        if args.out_candidates_csv
        else out_csv.with_name(f"{out_csv.stem}_candidates.csv")
    )
    if out_candidates_csv.resolve() == out_csv.resolve():
        raise RuntimeError("--out-csv and --out-candidates-csv must be different files")
    write_reports(
        matches,
        rj_only,
        ppg12_only,
        Path(args.out_md),
        out_csv,
        out_candidates_csv,
        args.base_e_model,
        args.base_v3e_model,
        lane_id=args.lane_id or "",
        runtime_contract_sha256=runtime_contract_sha256,
    )
    print(
        f"matched={len(matches)} unmatched={len(rj_only)} ppg12_only={len(ppg12_only)} "
        f"out_md={args.out_md} out_csv={args.out_csv} "
        f"out_candidates_csv={out_candidates_csv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
