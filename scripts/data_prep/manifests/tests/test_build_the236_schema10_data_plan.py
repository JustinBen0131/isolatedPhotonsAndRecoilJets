from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path
import unittest

from scripts.data_prep.manifests import build_the236_schema10_data_plan as compiler
from scripts.data_prep.manifests import seal_the236_schema10_data_sources as sealer


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Schema10SparseDataPlanTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.runtime = self._runtime_receipt()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write(self, name: str, payload: dict) -> Path:
        path = self.root / name
        path.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
        return path

    def _runtime_receipt(self) -> Path:
        bindings = {}
        for system in ("pp", "auau"):
            install = self.root / f"{system}_install"
            install.mkdir()
            binding = {
                "install_prefix": str(install),
                "release": "ana.544",
                "offline_main": "/cvmfs/sphenix/release/release_ana/ana.544",
            }
            library = install / f"libRecoilJets{'AuAu' if system == 'auau' else ''}.so"
            library.write_text(f"{system}-library\n", encoding="utf-8")
            binding["library_path"] = str(library)
            binding["library_sha256"] = digest(library)
            binding["model_sha256"] = ("5" if system == "pp" else "6") * 64
            for kind in ("wrapper", "macro", "config"):
                path = self.root / f"{system}_{kind}"
                content = f"{system}-{kind}\n"
                if kind == "config":
                    content += (
                        f"schema10_data_retention_profile: {compiler.PROFILE}\n"
                        "schema10_data_retention_contract_sha256: "
                        f"{compiler.CONTRACT_SHA256}\n"
                    )
                path.write_text(content, encoding="utf-8")
                binding[f"{kind}_path"] = str(path)
                binding[f"{kind}_sha256"] = digest(path)
            bindings[system] = binding
        return self._write(
            "runtime.json",
            {
                "schema": "THE236Schema10SparseDataRuntimeReceiptV1",
                "status": "PASS_LOCAL_BUILD_UNSUBMITTED",
                "data_retention_profile": compiler.PROFILE,
                "data_retention_contract_sha256": compiler.CONTRACT_SHA256,
                "replay_environment": {
                    "RJ_ANALYSIS_JET_NODE_SUFFIX": "_RJ",
                    "RJ_REPLAY_SCHEMA_SHA256": "1" * 64,
                    "RJ_REPLAY_SEMANTIC_SHA256": "2" * 64,
                    "RJ_REPLAY_MODEL_SHA256": "3" * 64,
                    "RJ_REPLAY_CODE_SHA256": "4" * 64,
                    "RJ_AUTO_MERGE": "0",
                    "RJ_ALLOW_NONZERO_WITH_ROOT_OUTPUT": "0",
                },
                **bindings,
            },
        )

    @staticmethod
    def _record(ordinal: int, run: int, segment: int, positive: bool = False) -> dict:
        return {
            "ordinal": ordinal,
            "source_identity": f"source-{ordinal}",
            "run": run,
            "segment": segment,
            "dst_jet": f"/sphenix/test/DST_JET/run{run}/segment{segment}.root",
            "dst_jetcalo": f"/sphenix/test/DST_JETCALO/run{run}/segment{segment}.root",
            "expected_events": 100 + ordinal,
            "photon_positive_evidence": positive,
        }

    def _source_receipts(self) -> tuple[Path, Path]:
        pp_records = [self._record(index, 50000 + index, index, index == 3) for index in range(5)]
        auau_records = [
            self._record(0, 60000, 0),
            self._record(1, 60000, 1, True),
            self._record(2, 60000, 2),
            self._record(3, 60001, 0),
            self._record(4, 60001, 1),
        ]
        pp = self._write(
            "pp.json",
            {
                "schema": "THE236PPDataSourceReceiptV1",
                "status": "PASS_FROZEN",
                "system": "pp",
                "source_authority": "THE-110",
                "historical_group_size": 20,
                "source_count": len(pp_records),
                "records": pp_records,
            },
        )
        auau = self._write(
            "auau.json",
            {
                "schema": "THE236AuAuPairedSourceReceiptV1",
                "status": "PASS_FROZEN",
                "system": "auau",
                "source_count": len(auau_records),
                "accepted_run_count": 2,
                "pairing_report_sha256": compiler.PAIRING_SHA256,
                "grouping": "RUN_BOUNDED_ORDERED_GROUPS_OF_10",
                "records": auau_records,
            },
        )
        return pp, auau

    def _compile(self, measurement: Path | None = None):
        pp, auau = self._source_receipts()
        return compiler.compile_plan(
            pp_receipt_path=pp,
            auau_receipt_path=auau,
            runtime_path=self.runtime,
            output_root="/sphenix/tg/tg01/bulk/jbennett/schema10_data_test",
            measurement_path=measurement,
            expected_counts={
                "pp_sources": 5,
                "auau_runs": 2,
                "auau_sources": 5,
                "pp_rows": 1,
                "auau_rows": 2,
                "auau_pilot": 1,
            },
            simulation_bytes=1_000,
        )

    def test_memory_rule_matches_accepted_canaries_and_stops_above_4gb(self) -> None:
        self.assertEqual(compiler.memory_request(2_219), 3_072)
        self.assertEqual(compiler.memory_request(1_827), 2_560)
        self.assertEqual(compiler.memory_request(3_276), 4_096)
        with self.assertRaisesRegex(ValueError, "above the 4096 MB stop"):
            compiler.memory_request(3_277)

    def test_compile_is_exact_run_bounded_disjoint_and_unsubmitted(self) -> None:
        payload, artifacts = self._compile()
        self.assertEqual(payload["status"], "PASS_FROZEN_UNSUBMITTED")
        self.assertFalse(payload["production_submission_authorized"])
        self.assertEqual(payload["science_schema_version"], 10)
        self.assertEqual(payload["data_retention_contract_sha256"], compiler.CONTRACT_SHA256)
        self.assertEqual(payload["source_authority"]["pp"]["row_count"], 1)
        self.assertEqual(payload["source_authority"]["auau"]["row_count"], 2)
        self.assertEqual(payload["auau_staging"]["pilot_rows"], 1)
        self.assertEqual(payload["auau_staging"]["complement_rows"], 1)
        self.assertFalse(payload["site_contract"]["getenv"])
        self.assertFalse(payload["site_contract"]["maxjobs"])
        self.assertFalse(payload["site_contract"]["automatic_retry"])
        rows = [json.loads(line) for line in artifacts["auau/rows.jsonl"].splitlines()]
        self.assertEqual([row["run"] for row in rows], ["60000", "60001"])
        self.assertEqual(sum(row["source_count"] for row in rows), 5)
        self.assertEqual(len(artifacts["pp/source_records.tsv"].splitlines()), 5)
        duplicate = json.loads(artifacts["duplicate_exclusion.json"])
        self.assertTrue(duplicate["pilot_complement_disjoint"])

    def test_measurements_recompute_memory_and_storage_admission(self) -> None:
        measurement = self._write(
            "measurement.json",
            {
                "schema": "THE236Schema10SparseDataCanaryTerminalReceiptV1",
                "status": "PASS",
                "rows": [
                    {
                        "row_id": row_id,
                        "exit_code": 0,
                        "num_job_starts": 1,
                        "holds": 0,
                        "releases": 0,
                        "retries": 0,
                        "processed_events": 1_000,
                        "scientific_output_size_bytes": 2_000,
                        "peak_memory_mb": 2_219 if row_id.startswith("pp") else 1_827,
                    }
                    for row_id in (
                        "pp.typical",
                        "pp.photon_positive",
                        "auau.typical",
                        "auau.photon_positive",
                    )
                ],
            },
        )
        payload, _ = self._compile(measurement)
        self.assertEqual(payload["status"], "READY_FOR_EXACT_CONDITIONAL_RUNDOWN_UNSUBMITTED")
        self.assertTrue(payload["storage_admission"]["admitted"])
        self.assertEqual(payload["resources"]["request_memory_mb"], {"pp": 3_072, "auau": 2_560})

    def test_duplicate_source_identity_is_rejected(self) -> None:
        pp, _ = self._source_receipts()
        payload = json.loads(pp.read_text(encoding="utf-8"))
        payload["records"][1]["source_identity"] = payload["records"][0]["source_identity"]
        pp.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate pp source identity"):
            compiler.parse_source_receipt(pp, "pp", expected_sources=5)

    def test_materialization_is_exclusive(self) -> None:
        payload, artifacts = self._compile()
        target = self.root / "packet"
        compiler.materialize(target, payload, artifacts)
        with self.assertRaises(FileExistsError):
            compiler.materialize(target, payload, artifacts)

    def test_source_sealer_writes_external_content_addressed_receipt(self) -> None:
        input_path = self.root / "normalized.jsonl"
        input_path.write_text(
            "".join(
                json.dumps(self._record(index, 61000 + index, index, index == 1), sort_keys=True)
                + "\n"
                for index in range(3)
            ),
            encoding="utf-8",
        )
        output = self.root / "sealed_pp"
        result = sealer.seal(
            system="pp",
            input_jsonl=input_path,
            output=output,
            authority_sha256="a" * 64,
            expected_sources=3,
            expected_runs=None,
        )
        receipt = json.loads(Path(result["receipt"]).read_text(encoding="utf-8"))
        self.assertEqual(receipt["source_authority"], "THE-110")
        self.assertEqual(receipt["source_count"], 3)
        records = Path(receipt["records_jsonl_path"])
        self.assertEqual(digest(records), receipt["records_jsonl_sha256"])
        parsed, _ = compiler.parse_source_receipt(
            Path(result["receipt"]), "pp", expected_sources=3
        )
        self.assertEqual([record.source_identity for record in parsed], ["source-0", "source-1", "source-2"])

    def test_worker_publication_is_candidate_validate_atomic(self) -> None:
        worker = Path(
            "/Users/patsfan753/Desktop/ThesisAnalysis/"
            "scripts/sdcc/runtime/condor/run_the236_schema10_data_row.sh"
        ).read_text(encoding="utf-8")
        self.assertIn('candidate_path="${output_path}.part"', worker)
        self.assertIn('source /opt/sphenix/core/bin/sphenix_setup.sh -n "$release"', worker)
        self.assertIn('/bin/mv -- "$candidate_path" "$output_path"', worker)
        self.assertLess(worker.index("candidate ROOT validation failed"), worker.index('/bin/mv -- "$candidate_path"'))
        # The worker spells these out in full; asserting the abbreviated form
        # made this test fail against a worker that enforces the check correctly.
        self.assertIn("processed_events!=expected_processed", worker)
        self.assertIn("processed-event count differs", worker)
        # These three assertions predate the Cling-to-PyROOT rewrite of the
        # candidate validator (the worker documents that change: the old
        # one-liner "could falsely return an empty TNamed title on valid
        # persisted metadata and turn successful science into exit 2"). The
        # contract they guard is unchanged, so assert it against the PyROOT
        # implementation instead of the retired C++ literals.
        for invariant in ("truth_photons != 0", "truth_jets != 0", "links != 0"):
            self.assertIn(invariant, worker)
        self.assertIn("int(value.GetEntries())", worker)
        self.assertIn('expected = parse_nonnegative(f"{tree}_entries")', worker)
        self.assertIn("if expected != actual:", worker)
        self.assertIn("RJ_DATA_ROOT_CONTRACT_V1 table mismatch tree=", worker)
        self.assertIn("printf '%s\\n' \"$root_contract\" >&2", worker)
        self.assertIn('event_limit="${14:-0}"', worker)
        self.assertIn('"$ClusterId" "$event_limit" "$ProcId"', worker)


if __name__ == "__main__":
    unittest.main()
