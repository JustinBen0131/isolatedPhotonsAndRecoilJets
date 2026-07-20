#!/usr/bin/env python3
"""Focused tests for the paired-oracle source-manifest builder."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "build_ppg12_paired_source_manifest.py"
REPO = Path(__file__).resolve().parents[4]
SPEC = importlib.util.spec_from_file_location("paired_source_builder", SCRIPT)
assert SPEC and SPEC.loader
BUILDER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILDER)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Fixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.runtime_dir = root / "runtime_attempt12"
        self.runtime_dir.mkdir()
        self.receipt = self._file("runtime_attempt12/build_receipt.json", "{}\n")
        self.role_paths: dict[str, Path] = {}
        for index, role in enumerate(sorted(BUILDER.REQUIRED_RUNTIME_ROLES)):
            safe = role.replace("/", "_")
            self.role_paths[role] = self._file(
                f"runtime_attempt12/assets/{index:02d}_{safe}", f"asset {role}\n"
            )
        self.runtime_manifest = self.runtime_dir / "runtime_manifest.json"
        self._write_runtime()
        self.setup = self._file("explicit/sphenix_setup.sh", "setup\n")
        self.mask = self._file("explicit/tower_mask.txt", "mask\n")
        self.config = self._file("explicit/analysis_config.yaml", "config\n")
        self.lane_tsv = root / "lane_sources.tsv"
        self._write_lanes()

    def _file(self, relative: str, text: str) -> Path:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path.resolve()

    def _write_runtime(self) -> None:
        payload = {
            "schema_version": 1,
            "runtime_profile": "new.17",
            "isolated_build": True,
            "estimator_revision": BUILDER.ESTIMATOR_REVISION,
            "build_receipt": str(self.receipt),
            "build_receipt_sha256": sha256(self.receipt),
            "files": [
                {
                    "role": role,
                    "path": str(path),
                    "sha256": sha256(path),
                }
                for role, path in reversed(list(self.role_paths.items()))
            ],
        }
        self.runtime_manifest.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def _write_lanes(self) -> None:
        rows: list[dict[str, str]] = []
        for lane_id, sample, period, interaction in reversed(BUILDER._expected_lanes()):
            # Real preserved shape: period lanes reuse one sample/interaction
            # macro and one paired G4/TRUTH_JET list contract.
            source_slug = f"{sample.lower()}_{interaction.lower()}"
            macro_body = """// Preserved aux inputs are intentionally inactive.
INPUTREADHITS::listfile[0] = \"g4.list\";
// INPUTREADHITS::listfile[1] = \"calo.list\";
/* INPUTREADHITS::listfile[2] = \"cluster.list\"; */
// INPUTREADHITS::listfile[3] = \"mbd.list\";
INPUTREADHITS::listfile[4] = \"truth.list\";
TruthJetInput *truth_input = new TruthJetInput(Jet::PARTICLE);
"""
            if interaction == "DI":
                macro_body += "truth_input->add_embedding_flag(2);\n"
            macro = self._file(
                f"lanes/{source_slug}.C", macro_body
            )
            g4 = self._file(
                f"lanes/{source_slug}_g4.list",
                "".join(
                    f"G4Hits_pythia8_{sample}_{interaction}_{index:06d}.root\n"
                    for index in range(5)
                ),
            )
            truth = self._file(
                f"lanes/{source_slug}_truth.list",
                "".join(
                    f"DST_TRUTH_JET_pythia8_{sample}_{interaction}_{index:06d}.root\n"
                    for index in range(5)
                ),
            )
            rows.append(
                {
                    "lane_id": lane_id,
                    "sample": sample,
                    "period": period,
                    "interaction": interaction,
                    "ppg_macro": str(macro),
                    "g4_full_list": str(g4),
                    "truthjet_full_list": str(truth),
                }
            )
        with self.lane_tsv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=BUILDER.LANE_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def kwargs(self) -> dict[str, object]:
        return {
            "runtime_manifest": self.runtime_manifest.resolve(),
            "expected_runtime_sha256": sha256(self.runtime_manifest),
            "setup_script": self.setup,
            "tower_mask": self.mask,
            "recoil_config": self.config,
            "lane_tsv": self.lane_tsv.resolve(),
        }


class PairedSourceManifestBuilderTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="paired-source-builder-")
        self.root = Path(self.temporary.name).resolve()
        self.fixture = Fixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_builds_exact_schema_in_canonical_order_deterministically(self) -> None:
        payload = BUILDER.build_manifest(**self.fixture.kwargs())
        expected_lane_ids = [row[0] for row in BUILDER._expected_lanes()]
        self.assertEqual(payload["schema"], "ppg12-paired-source-manifest/v1")
        self.assertEqual([row["lane_id"] for row in payload["lanes"]], expected_lane_ids)
        self.assertEqual(
            payload["common"]["base_v3e_model"],
            str(self.fixture.role_paths["ppg_apply_model_base_v3E"]),
        )
        first = self.root / "out" / "source_a.json"
        second = self.root / "out" / "source_b.json"
        self.assertEqual(BUILDER.write_manifest(first, payload), "WROTE")
        self.assertEqual(BUILDER.write_manifest(second, payload), "WROTE")
        self.assertEqual(first.read_bytes(), second.read_bytes())
        self.assertEqual(BUILDER.write_manifest(first, payload), "UNCHANGED")

    def test_output_is_accepted_by_the_existing_twelve_lane_plan(self) -> None:
        source = self.root / "source.json"
        BUILDER.write_manifest(source, BUILDER.build_manifest(**self.fixture.kwargs()))
        evidence = self.root / "plan_evidence"
        output = self.root / "paired_output"
        driver = (
            REPO
            / "scripts/sdcc/workflows/diagnostics/submit_ppg12_stitched_purity_photon_canaries.sh"
        )
        environment = {
            **os.environ,
            "RJ_PPG12_PAIRED_SOURCE_MANIFEST": str(source),
            "RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG": "paired_source_builder_test",
            "RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR": str(evidence),
            "RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT": str(output),
        }
        completed = subprocess.run(
            [str(driver), "plan"],
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr or completed.stdout)
        plan = json.loads((evidence / "photon_canary_plan.json").read_text())
        self.assertEqual(plan["schema"], "ppg12-stitched-purity-photon-canary-plan/v5")
        self.assertEqual(len(plan["lanes"]), 12)

    def test_rejects_runtime_hash_drift_and_missing_required_role(self) -> None:
        kwargs = self.fixture.kwargs()
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "explicitly authorized"):
            BUILDER.build_manifest(**{**kwargs, "expected_runtime_sha256": "0" * 64})

        payload = json.loads(self.fixture.runtime_manifest.read_text(encoding="utf-8"))
        payload["files"] = [
            row for row in payload["files"] if row["role"] != "ppg_apply_model_base_E"
        ]
        self.fixture.runtime_manifest.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        kwargs["expected_runtime_sha256"] = sha256(self.fixture.runtime_manifest)
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "lacks paired-oracle roles"):
            BUILDER.build_manifest(**kwargs)

    def test_requires_period_configs_truth_weights_and_yaml_headers(self) -> None:
        required = BUILDER.REQUIRED_RUNTIME_ROLES
        self.assertNotIn("ppg_recoeff_canonical_config", required)
        for role in (
            "ppg_recoeff_period_config_0mrad",
            "ppg_recoeff_period_config_1p5mrad",
            "ppg_recoeff_truth_vertex_reweight_0mrad",
            "ppg_recoeff_truth_vertex_reweight_1p5mrad",
            "ppg_recoeff_yaml_cpp_header_tree_receipt",
        ):
            self.assertIn(role, required)
            payload = json.loads(
                self.fixture.runtime_manifest.read_text(encoding="utf-8")
            )
            payload["files"] = [
                row for row in payload["files"] if row["role"] != role
            ]
            self.fixture.runtime_manifest.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            kwargs = self.fixture.kwargs()
            kwargs["expected_runtime_sha256"] = sha256(
                self.fixture.runtime_manifest
            )
            with self.assertRaisesRegex(
                BUILDER.ManifestBuildError, "lacks paired-oracle roles"
            ):
                BUILDER.build_manifest(**kwargs)
            self.fixture._write_runtime()

    def test_accepts_exact_cross_period_reuse_and_rejects_other_sharing(self) -> None:
        payload = BUILDER.build_manifest(**self.fixture.kwargs())
        by_id = {row["lane_id"]: row for row in payload["lanes"]}
        for photon in (5, 10, 20):
            for interaction in ("si", "di"):
                zero = by_id[f"photon:photon{photon}:0mrad:{interaction}"]
                shifted = by_id[f"photon:photon{photon}:1p5mrad:{interaction}"]
                for field in BUILDER.LANE_PATH_FIELDS:
                    self.assertEqual(zero[field], shifted[field])

        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        target = next(
            row
            for row in rows
            if row["lane_id"] == "photon:photon5:1p5mrad:si"
        )
        target["g4_full_list"] = str(
            self.fixture._file(
                "lanes/photon5_si_wrong_period_g4.list",
                "".join(
                    f"G4Hits_pythia8_Photon5_SI_wrong_{index:06d}.root\n"
                    for index in range(5)
                ),
            )
        )
        target["truthjet_full_list"] = str(
            self.fixture._file(
                "lanes/photon5_si_wrong_period_truth.list",
                "".join(
                    f"DST_TRUTH_JET_pythia8_Photon5_SI_wrong_{index:06d}.root\n"
                    for index in range(5)
                ),
            )
        )
        with self.fixture.lane_tsv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=BUILDER.LANE_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "exactly shared"):
            BUILDER.build_manifest(**self.fixture.kwargs())

    def test_rejects_nonabsolute_and_cross_contract_duplicate_lane_paths(self) -> None:
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        rows[0]["ppg_macro"] = "relative/macro.C"
        with self.fixture.lane_tsv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=BUILDER.LANE_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "must be absolute"):
            BUILDER.build_manifest(**self.fixture.kwargs())

        self.fixture._write_lanes()
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        source = next(row for row in rows if row["lane_id"] == "photon:photon5:0mrad:si")
        target = next(row for row in rows if row["lane_id"] == "photon:photon10:0mrad:si")
        target["ppg_macro"] = source["ppg_macro"]
        counterpart = next(
            row for row in rows if row["lane_id"] == "photon:photon10:1p5mrad:si"
        )
        counterpart["ppg_macro"] = source["ppg_macro"]
        with self.fixture.lane_tsv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=BUILDER.LANE_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "illegally shares"):
            BUILDER.build_manifest(**self.fixture.kwargs())

    def test_rejects_si_di_lists_with_identical_source_rows(self) -> None:
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        si = next(row for row in rows if row["lane_id"] == "photon:photon20:0mrad:si")
        di = next(row for row in rows if row["lane_id"] == "photon:photon20:0mrad:di")
        Path(di["g4_full_list"]).write_text(
            Path(si["g4_full_list"]).read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        Path(di["truthjet_full_list"]).write_text(
            Path(si["truthjet_full_list"]).read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "source rows must be distinct"):
            BUILDER.build_manifest(**self.fixture.kwargs())

    def test_relative_basenames_are_paired_by_exact_event_identity(self) -> None:
        payload = BUILDER.build_manifest(**self.fixture.kwargs())
        self.assertEqual(len(payload["lanes"]), 12)

        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        lane = next(row for row in rows if row["lane_id"] == "photon:photon10:0mrad:si")
        truth_path = Path(lane["truthjet_full_list"])
        truth_rows = truth_path.read_text(encoding="utf-8").splitlines()
        truth_rows[2] = truth_rows[2].replace("000002", "999999")
        truth_path.write_text("\n".join(truth_rows) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "event identity mismatch"):
            BUILDER.build_manifest(**self.fixture.kwargs())

    def test_rejects_unrecognized_source_prefix(self) -> None:
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        lane = next(row for row in rows if row["lane_id"] == "photon:photon5:0mrad:di")
        g4_path = Path(lane["g4_full_list"])
        g4_rows = g4_path.read_text(encoding="utf-8").splitlines()
        g4_rows[0] = g4_rows[0].replace("G4Hits_", "DST_CALO_CLUSTER_")
        g4_path.write_text("\n".join(g4_rows) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "recognized G4Hits"):
            BUILDER.build_manifest(**self.fixture.kwargs())

    def test_macro_contract_rejects_active_aux_input(self) -> None:
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        lane = next(row for row in rows if row["lane_id"] == "photon:photon5:0mrad:si")
        macro = Path(lane["ppg_macro"])
        macro.write_text(
            macro.read_text(encoding="utf-8")
            + 'INPUTREADHITS::listfile[1] = "calo.list";\n',
            encoding="utf-8",
        )
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "exactly \\{0, 4\\}"):
            BUILDER.build_manifest(**self.fixture.kwargs())

    def test_macro_contract_requires_di_embedding_and_forbids_it_in_si(self) -> None:
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        di = next(row for row in rows if row["lane_id"] == "photon:photon10:0mrad:di")
        di_macro = Path(di["ppg_macro"])
        di_macro.write_text(
            di_macro.read_text(encoding="utf-8").replace(
                "truth_input->add_embedding_flag(2);", "// removed"
            ),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "must actively call"):
            BUILDER.build_manifest(**self.fixture.kwargs())

        self.fixture._write_lanes()
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        di = next(row for row in rows if row["lane_id"] == "photon:photon10:0mrad:di")
        di_macro = Path(di["ppg_macro"])
        di_macro.write_text(
            di_macro.read_text(encoding="utf-8").replace(
                "new TruthJetInput", "existing_TruthJetInput"
            ),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "reconstruct truth jets"):
            BUILDER.build_manifest(**self.fixture.kwargs())

        self.fixture._write_lanes()
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        si = next(row for row in rows if row["lane_id"] == "photon:photon10:0mrad:si")
        si_macro = Path(si["ppg_macro"])
        si_macro.write_text(
            si_macro.read_text(encoding="utf-8")
            + "truth_input->add_embedding_flag(2);\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "must not actively call"):
            BUILDER.build_manifest(**self.fixture.kwargs())

    def test_rejects_missing_lane_and_unpaired_source_lists(self) -> None:
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        rows.pop()
        with self.fixture.lane_tsv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=BUILDER.LANE_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "exact 12 lanes"):
            BUILDER.build_manifest(**self.fixture.kwargs())

        self.fixture._write_lanes()
        with self.fixture.lane_tsv.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        truth_path = Path(rows[0]["truthjet_full_list"])
        truth_path.write_text("/truth/only/one.root\n", encoding="utf-8")
        with self.assertRaisesRegex(BUILDER.ManifestBuildError, "at least five"):
            BUILDER.build_manifest(**self.fixture.kwargs())


if __name__ == "__main__":
    unittest.main()
