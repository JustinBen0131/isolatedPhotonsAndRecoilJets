from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "the94_final_sim_abcd_purity_closure_audit.py"
)
SPEC = importlib.util.spec_from_file_location("the94_forensics", SCRIPT)
assert SPEC and SPEC.loader
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def histogram(values: list[float], variances: list[float] | None = None) -> dict[str, object]:
    return {
        "values": values,
        "variances": variances or [abs(value) for value in values],
        "flow_values": [0.0, 0.0],
        "flow_variances": [0.0, 0.0],
        "sumw2_present": True,
        "edges": [10.0, 12.0, 14.0],
        "entries": float(sum(values)),
        "object": "synthetic",
    }


def lane_histograms(scale: float = 1.0) -> dict[str, dict[str, object]]:
    base = {
        "A": [100.0, 80.0],
        "B": [20.0, 16.0],
        "C": [30.0, 24.0],
        "D": [40.0, 32.0],
        "A_signal": [80.0, 64.0],
        "B_signal": [8.0, 6.4],
        "C_signal": [12.0, 9.6],
        "D_signal": [16.0, 12.8],
        "A_notmatch": [20.0, 16.0],
        "B_notmatch": [12.0, 9.6],
        "C_notmatch": [18.0, 14.4],
        "D_notmatch": [24.0, 19.2],
    }
    return {
        region: histogram(
            [scale * value for value in values],
            [scale * scale * value for value in values],
        )
        for region, values in base.items()
    }


def sum_histograms(records: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    output: dict[str, dict[str, object]] = {}
    for region in AUDIT.REGIONS:
        values = np.sum(
            [np.asarray(record["histograms"][region]["values"]) for record in records], axis=0
        )
        variances = np.sum(
            [np.asarray(record["histograms"][region]["variances"]) for record in records], axis=0
        )
        output[region] = histogram(values.tolist(), variances.tolist())
    return output


def synthetic_snapshot(*, mutate_lane: tuple[str, str, str] | None = None) -> dict[str, object]:
    records: list[dict[str, object]] = []
    current_records: list[dict[str, object]] = []
    ppg12_records: list[dict[str, object]] = []
    for lane_index, (sample, period, interaction) in enumerate(AUDIT.EXPECTED_LANES, start=1):
        scale = float(lane_index)
        ppg12_record = {
            "label": AUDIT.lane_label("ppg12", sample, period, interaction),
            "exists": True,
            "path": f"/ppg12/{period}_{interaction}_{sample}.root",
            "config": {"sha256": "ppg12-config"},
            "bytes": 1000 + lane_index,
            "mtime_ns": 1_000_000_000 + lane_index,
            "sha256": hashlib.sha256(
                f"ppg12:{sample}:{period}:{interaction}".encode()
            ).hexdigest(),
            "histograms": lane_histograms(scale),
        }
        current_histograms = lane_histograms(scale)
        if mutate_lane == (sample, period, interaction):
            current_histograms["B"]["values"][0] += 50.0
            current_histograms["B"]["variances"][0] += 50.0
            current_histograms["B"]["entries"] += 50.0
        current_record = {
            "label": AUDIT.lane_label("current", sample, period, interaction),
            "exists": True,
            "path": f"/current/{period}_{interaction}_{sample}.root",
            "config": {"sha256": f"current-config-{interaction}"},
            "bytes": 2000 + lane_index,
            "mtime_ns": 2_000_000_000 + lane_index,
            "sha256": hashlib.sha256(
                f"current:{sample}:{period}:{interaction}".encode()
            ).hexdigest(),
            "histograms": current_histograms,
        }
        ppg12_records.append(ppg12_record)
        current_records.append(current_record)
        records.extend((ppg12_record, current_record))
    records.extend(
        (
            {
                "label": "ppg12:final",
                "exists": True,
                "path": "/ppg12/final.root",
                "config": {"sha256": "ppg12-config"},
                "bytes": 3000,
                "mtime_ns": 3_000_000_000,
                "sha256": hashlib.sha256(b"ppg12:final").hexdigest(),
                "histograms": sum_histograms(ppg12_records),
            },
            {
                "label": "current:final",
                "exists": True,
                "path": "/current/final.root",
                "config": {"sha256": "current-config-si"},
                "bytes": 4000,
                "mtime_ns": 4_000_000_000,
                "sha256": hashlib.sha256(b"current:final").hexdigest(),
                "histograms": sum_histograms(current_records),
            },
        )
    )
    lane_paths = [
        f"/current/{period}_{interaction}_{sample}.root"
        for sample, period, interaction in AUDIT.EXPECTED_LANES
    ]
    payload = {
        "records": records,
        "lane_manifest_sha256": "synthetic-lane-manifest",
        "merge_manifest": {
            "lane_count": 20,
            "lane_roots_fixed_order": lane_paths,
            "lane_roots_manifest_sha256": "synthetic-lane-manifest",
            "final_root": "/current/final.root",
            "final_root_sha256": hashlib.sha256(b"current:final").hexdigest(),
            "jet5": "excluded",
            "merge_contract": "fixed-order additive; event-preweighted; no external scale",
        },
    }
    mapped = AUDIT.records_by_label(records)
    payload["record_set_sha256"] = AUDIT.canonical_sha256(
        AUDIT.snapshot_record_binding(mapped)
    )
    return payload


def merge_inputs(tmp_path: Path) -> tuple[Path, Path]:
    current_root = tmp_path / "current.root"
    current_root.write_bytes(b"synthetic-root")
    digest = hashlib.sha256(current_root.read_bytes()).hexdigest()
    merge_audit = tmp_path / "merge.json"
    merge_audit.write_text(
        json.dumps(
            {
                "status": "PASS",
                "inputs_fixed_order_count": 20,
                "max_content_delta": 0.0,
                "max_sumw2_delta": 0.0,
                "failures": [],
                "output_sha256": digest,
            }
        )
    )
    return current_root, merge_audit


def exact_estimator_payload() -> dict[str, object]:
    purity = {
        "bin_edges": [10.0, 12.0, 14.0],
        "truth": {"value": [0.8, 0.8], "error": [0.01, 0.01]},
        "raw": {"value": [0.85, 0.85], "error": [0.02, 0.02]},
        "corrected": {"value": [0.82, 0.82], "error": [0.025, 0.025]},
    }
    diagnostics = []
    for index, (low, high) in enumerate(((10.0, 12.0), (12.0, 14.0)), start=1):
        diagnostics.append(
            {
                "bin": index,
                "x_low": low,
                "x_high": high,
                "A": 100.0,
                "B": 20.0,
                "C": 30.0,
                "D": 40.0,
                "A_effective": 100.0,
                "B_effective": 20.0,
                "C_effective": 30.0,
                "D_effective": 40.0,
                "cB": 0.10,
                "cC": 0.15,
                "cD": 0.20,
                "cB_error": 0.01,
                "cC_error": 0.01,
                "cD_error": 0.01,
            }
        )
    purity_sha = AUDIT.PURITY_PRODUCER._payload_sha256(purity)
    diagnostics_sha = AUDIT.PURITY_PRODUCER._payload_sha256(diagnostics)
    return {
        "schema": "ppg12-stitched-purity-purity/v1",
        "random_seed": 42,
        "toy_count": 20000,
        "purity": purity,
        "fixed_seed_repetition": {
            "first_output_sha256": purity_sha,
            "repeated_output_sha256": purity_sha,
        },
        "run_diagnostics": {
            "first": diagnostics,
            "repeated": json.loads(json.dumps(diagnostics)),
            "diagnostics_sha256": diagnostics_sha,
            "repeated_diagnostics_sha256": diagnostics_sha,
        },
        "algorithm": {
            "producer": {
                "path": str(AUDIT.PURITY_EVIDENCE_PRODUCER_PATH),
                "sha256": AUDIT.sha256(AUDIT.PURITY_EVIDENCE_PRODUCER_PATH),
            },
            "source": {
                "path": str(AUDIT.PURITY_PRODUCER.PPG12_ESTIMATOR_SOURCE),
                "sha256": AUDIT.sha256(
                    Path(AUDIT.PURITY_PRODUCER.PPG12_ESTIMATOR_SOURCE)
                ),
            },
        },
    }


def run_analysis(tmp_path: Path, snapshot_payload: dict[str, object]) -> dict[str, object]:
    snapshot = tmp_path / "snapshot.json"
    snapshot.write_text(json.dumps(snapshot_payload))
    current_root, merge_audit = merge_inputs(tmp_path)
    jet8 = tmp_path / "jet8"
    jet8.mkdir()
    (jet8 / "ppg12_jet8_reference_provenance_manifest.json").write_text(
        json.dumps(
            {
                "observed_factor": 2.339817922469992,
                "known_old_xsec_factor": 1.1315652173913044,
                "forced_effective_xsec_pb": 4914912.348333586,
            }
        )
    )
    fields = [
        "current_over_ppg12_integral",
        "current_over_ppg12_raw_entries",
        "current_over_ppg12_integral_per_entry",
    ]
    with (jet8 / "jet8_weight_factor_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(dict(zip(fields, (2.3381, 1.0021, 2.3331))))
    (jet8 / "ppg12_jet8_real_ratio_search.csv").write_text(
        "quantity,value,explains_2p338,note\n"
        "observed_current_over_ppg12_true_period,2.3398,True,Measured failure\n"
        "forced_effective_jet8_xsec_pb,4914912,True,Numerical only\n"
    )
    return AUDIT.analyze(
        snapshot,
        tmp_path / "out",
        current_root,
        merge_audit,
        jet8,
        max_ratio_deviation=1.0e-4,
        max_ratio_rms=2.0e-5,
    )


class AuditTests(unittest.TestCase):
    def test_th1f_addition_reproduces_sequential_root_rounding(self) -> None:
        records: list[dict[str, object]] = []
        expected = np.float32(0.0)
        for index in range(20):
            value = np.float32(1.0e9 if index == 0 else 1.0)
            expected = np.float32(expected + value)
            histograms = lane_histograms()
            for region in AUDIT.REGIONS:
                histograms[region]["values"] = [float(value), 0.0]
                histograms[region]["storage_class"] = "TH1F"
            records.append({"label": f"lane-{index}", "histograms": histograms})
        summed = AUDIT.sum_records(records)
        self.assertEqual(summed["A"][0][0], float(expected))
        self.assertNotEqual(summed["A"][0][0], 1.0e9 + 19.0)

    def test_exact_20_lane_contract_closes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            report = run_analysis(Path(directory), synthetic_snapshot())
        self.assertEqual(report["status"], "parity_closed")
        self.assertIs(report["structural_pass"], True)
        self.assertIs(report["parity_pass"], True)
        self.assertEqual(report["lane_contract"]["lane_count"], 20)
        self.assertEqual(
            report["lane_contract"]["external_or_fitted_scale"],
            "forbidden_and_not_present",
        )
        self.assertEqual(report["artifact_classification"], "historical_diagnostic_only")
        self.assertIs(report["admission_eligible"], False)
        self.assertEqual(report["exact_stitched_purity_attachment"]["status"], "NOT_SUPPLIED")
        self.assertEqual(
            report["implementation_binding"]["files"]["focused_regression_test"][
                "sha256"
            ],
            AUDIT.sha256(Path(__file__).resolve()),
        )

    def test_report_mask_includes_final_32_to_36_bin(self) -> None:
        mask = AUDIT.report_bin_mask(
            np.asarray([8.0, 10.0, 12.0, 32.0, 36.0, 40.0])
        )
        np.testing.assert_array_equal(mask, [False, True, True, True, False])
        self.assertEqual(AUDIT.REPORT_ET_MAX_GEV, 36.0)

    def test_exact_estimator_attachment_exposes_photon_leakage_components(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "purity.json"
            path.write_text(json.dumps(exact_estimator_payload()))
            result = AUDIT.validate_exact_stitched_purity_payload(
                path, np.asarray([10.0, 12.0, 14.0])
            )
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["random_seed"], 42)
        self.assertEqual(result["toy_count"], 20000)
        leakage = [
            row
            for row in result["component_rows"]
            if row["source_family"] == "photon" and row["component"] == "cC"
        ]
        self.assertEqual(len(leakage), 2)
        self.assertEqual(leakage[0]["value"], 0.15)

    def test_exact_estimator_attachment_rejects_changed_repeated_diagnostics(self) -> None:
        payload = exact_estimator_payload()
        payload["run_diagnostics"]["repeated"][0]["cB"] += 1.0e-9
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "purity.json"
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "not deterministic"):
                AUDIT.validate_exact_stitched_purity_payload(
                    path, np.asarray([10.0, 12.0, 14.0])
                )

    def test_lane_substitution_identifies_mutated_lane(self) -> None:
        mutated = ("jet20", "0mrad", "di")
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            report = run_analysis(tmp_path, synthetic_snapshot(mutate_lane=mutated))
            with (tmp_path / "out/lane_substitution_summary.csv").open(newline="") as handle:
                summary = list(csv.DictReader(handle))
            with (tmp_path / "out/lane_substitution_top_attribution.csv").open(
                newline=""
            ) as handle:
                attribution = list(csv.DictReader(handle))
            with (
                tmp_path / "out/lane_component_substitution_top_attribution.csv"
            ).open(newline="") as handle:
                component_attribution = list(csv.DictReader(handle))
        self.assertEqual(report["status"], "analyzed_nonclosing")
        self.assertIs(report["parity_pass"], False)
        self.assertEqual(summary[0]["group"], "jet20:0mrad:di")
        self.assertGreater(float(summary[0]["sum_abs_distance_improvement"]), 0.0)
        raw_first_bin = next(
            row for row in attribution
            if row["observable"] == "raw_central_algebraic" and row["bin"] == "1"
        )
        self.assertEqual(raw_first_bin["top_lane"], "jet20:0mrad:di")
        self.assertEqual(raw_first_bin["tested_lane_count"], "20")
        b_first_bin = next(
            row
            for row in component_attribution
            if row["component"] == "B" and row["bin"] == "1"
        )
        self.assertEqual(b_first_bin["top_lane"], "jet20:0mrad:di")
        self.assertEqual(b_first_bin["tested_lane_count"], "20")

    def test_sumw2_only_final_merge_mismatch_fails_structure(self) -> None:
        snapshot = synthetic_snapshot()
        final = next(
            record for record in snapshot["records"] if record["label"] == "current:final"
        )
        final["histograms"]["A"]["variances"][0] += 1.0
        with tempfile.TemporaryDirectory() as directory:
            report = run_analysis(Path(directory), snapshot)
        self.assertEqual(report["status"], "structural_failure")
        self.assertIs(report["current_additive_closure_pass"], False)

    def test_sumw2_lane_parity_mismatch_fails_even_when_merge_closes(self) -> None:
        snapshot = synthetic_snapshot()
        lane = next(
            record for record in snapshot["records"]
            if record["label"] == "current:jet8:0mrad:si"
        )
        lane["histograms"]["A"]["variances"][0] *= 100.0
        current_lanes = [
            record for record in snapshot["records"]
            if str(record["label"]).startswith("current:")
            and record["label"] != "current:final"
        ]
        final = next(
            record for record in snapshot["records"] if record["label"] == "current:final"
        )
        final["histograms"] = sum_histograms(current_lanes)
        with tempfile.TemporaryDirectory() as directory:
            report = run_analysis(Path(directory), snapshot)
        self.assertEqual(report["status"], "analyzed_nonclosing")
        self.assertIs(report["parity_pass"], False)

    def test_one_sided_zero_fails_parity_coverage(self) -> None:
        direct = AUDIT.ratio_metrics(np.asarray([0.0]), np.asarray([1.0]))
        self.assertIs(direct["coverage_pass"], False)
        self.assertEqual(direct["candidate_only_zero_bins"], 1)
        self.assertEqual(direct["ratio_max_abs_from_unity"], float("inf"))
        snapshot = synthetic_snapshot()
        lane = next(
            record for record in snapshot["records"]
            if record["label"] == "current:jet8:0mrad:si"
        )
        lane["histograms"]["A"]["values"][0] = 0.0
        lane["histograms"]["A"]["variances"][0] = 0.0
        current_lanes = [
            record for record in snapshot["records"]
            if str(record["label"]).startswith("current:")
            and record["label"] != "current:final"
        ]
        final = next(
            record for record in snapshot["records"] if record["label"] == "current:final"
        )
        final["histograms"] = sum_histograms(current_lanes)
        with tempfile.TemporaryDirectory() as directory:
            report = run_analysis(Path(directory), snapshot)
        self.assertIs(report["parity_pass"], False)

    def test_duplicate_label_fails_before_map_construction(self) -> None:
        snapshot = synthetic_snapshot()
        snapshot["records"][-1]["label"] = snapshot["records"][-2]["label"]
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "duplicate record labels"):
                run_analysis(Path(directory), snapshot)

    def test_truth_purity_uses_unsuffixed_A_denominator(self) -> None:
        histograms = {
            region: (
                np.asarray(values, dtype=float),
                np.asarray(values, dtype=float),
                np.asarray([10.0, 12.0, 14.0]),
            )
            for region, values in {
                "A": [100.0, 100.0], "B": [10.0, 10.0],
                "C": [10.0, 10.0], "D": [10.0, 10.0],
                "A_signal": [60.0, 60.0], "A_notmatch": [20.0, 20.0],
                "B_signal": [0.0, 0.0], "B_notmatch": [0.0, 0.0],
                "C_signal": [0.0, 0.0], "C_notmatch": [0.0, 0.0],
                "D_signal": [0.0, 0.0], "D_notmatch": [0.0, 0.0],
            }.items()
        }
        estimators = AUDIT.central_estimators(histograms)
        np.testing.assert_allclose(estimators["truth"], [0.6, 0.6])
        np.testing.assert_allclose(estimators["classified_fraction"], [0.8, 0.8])

    def test_missing_lane_fails_closed(self) -> None:
        snapshot = synthetic_snapshot()
        snapshot["records"] = [
            record
            for record in snapshot["records"]
            if record["label"] != "current:jet40:1p5mrad:di"
        ]
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "exactly 42 records"):
                run_analysis(Path(directory), snapshot)

    def test_missing_remote_root_hash_fails_closed(self) -> None:
        snapshot = synthetic_snapshot()
        lane = next(
            record for record in snapshot["records"]
            if record["label"] == "current:jet8:0mrad:si"
        )
        del lane["sha256"]
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "lacks a valid remote ROOT sha256"):
                run_analysis(Path(directory), snapshot)

    def test_stale_record_set_hash_fails_closed(self) -> None:
        snapshot = synthetic_snapshot()
        snapshot["record_set_sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "record-set hash mismatch"):
                run_analysis(Path(directory), snapshot)

    def test_stale_final_root_hash_in_merge_manifest_fails_structure(self) -> None:
        snapshot = synthetic_snapshot()
        snapshot["merge_manifest"]["final_root_sha256"] = "f" * 64
        with tempfile.TemporaryDirectory() as directory:
            report = run_analysis(Path(directory), snapshot)
        self.assertEqual(report["status"], "structural_failure")
        self.assertIn(
            "manifest_final_root_sha256_mismatch",
            report["snapshot_manifest_validation"]["failures"],
        )

    def test_jet8_provenance_reports_weight_not_population_and_never_scales(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            evidence = Path(directory) / "jet8"
            evidence.mkdir()
            (evidence / "ppg12_jet8_reference_provenance_manifest.json").write_text(
                json.dumps(
                    {
                        "observed_factor": 2.339817922469992,
                        "known_old_xsec_factor": 1.1315652173913044,
                        "forced_effective_xsec_pb": 4914912.348333586,
                    }
                )
            )
            fields = [
                "current_over_ppg12_integral",
                "current_over_ppg12_raw_entries",
                "current_over_ppg12_integral_per_entry",
            ]
            with (evidence / "jet8_weight_factor_comparison.csv").open(
                "w", newline=""
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(
                    [
                        dict(zip(fields, (2.3381, 1.0021, 2.3331))),
                        dict(zip(fields, (2.3383, 1.0022, 2.3333))),
                    ]
                )
            (evidence / "ppg12_jet8_real_ratio_search.csv").write_text(
                "quantity,value,explains_2p338,note\nobserved,2.3398,True,test\n"
            )
            result = AUDIT.audit_jet8_provenance(evidence)
        self.assertEqual(result["status"], "unresolved_weight_provenance")
        self.assertLess(result["raw_entry_ratio_max"], 1.01)
        self.assertGreater(result["integral_per_entry_ratio_min"], 2.3)
        self.assertIsNone(result["external_scale_applied"])
        self.assertIs(result["documented_explanation_recovered"], False)
        self.assertEqual(result["nominal_policy"], "no_external_or_fitted_jet8_scale")


if __name__ == "__main__":
    unittest.main()
