#!/usr/bin/env python3
"""Focused local tests for the exact 12-lane paired-oracle Condor fanout."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[5]
DRIVER = (
    REPO
    / "scripts/sdcc/workflows/diagnostics/submit_ppg12_stitched_purity_photon_canaries.sh"
)
HELPER = REPO / "scripts/sdcc/workflows/diagnostics/ppg12_paired_oracle_condor_fanout.py"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Fixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.assets = root / "assets"
        self.assets.mkdir()
        self.common: dict[str, str] = {}
        for name in (
            "setup_script",
            "apply_bdt",
            "apply_config",
            "base_e_model",
            "base_v3e_model",
            "npb_model",
            "tower_mask",
            "recoil_runtime_manifest",
            "recoil_config",
        ):
            path = self.assets / name
            path.write_text(f"asset {name}\n", encoding="utf-8")
            self.common[name] = str(path.resolve())

        self.lane_dir = root / "lanes"
        self.lane_dir.mkdir()
        self.source_manifest = root / "paired_source.json"
        self._write_source_manifest()
        self.mock_driver = root / "mock_paired_driver.sh"
        self.mock_driver.write_text(
            """#!/usr/bin/env bash
set -Eeuo pipefail
mode=plan
token=""
lane=""
output=""
while (( $# )); do
  case "$1" in
    --run) mode=run; shift ;;
    --token) token="$2"; shift 2 ;;
    --lane-id) lane="$2"; shift 2 ;;
    --output-dir) output="$2"; shift 2 ;;
    *) shift 2 ;;
  esac
done
expected="ppg12-oracle:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
if [[ "$mode" == plan ]]; then
  printf 'PPG12_PAIRED_ORACLE_PLAN\n  lane_id: %s\n  run_token: %s\n' "$lane" "$expected"
  exit 0
fi
[[ "$token" == "$expected" && "$output" == /* && ! -e "$output" ]]
mkdir -p "$output/comparison"
printf '{}\n' > "$output/paired_oracle_contract.json"
printf 'event,candidate\n1,1\n' > "$output/comparison/paired_oracle_candidates.csv"
printf '{}\n' > "$output/comparison/executable_aggregate.json"
printf 'PASS\n' > "$output/RUN_STATE"
""",
            encoding="utf-8",
        )
        self.mock_driver.chmod(0o700)
        self.evidence = root / "evidence"
        self.output = root / "output"
        self.environment = {
            **os.environ,
            "RJ_PPG12_PAIRED_SOURCE_MANIFEST": str(self.source_manifest.resolve()),
            "RJ_PPG12_PAIRED_ORACLE_DRIVER": str(self.mock_driver.resolve()),
            "RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG": "paired_condor_unit",
            "RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR": str(self.evidence.resolve()),
            "RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT": str(self.output.resolve()),
            "RJ_CODEX_CHAT_NAME": "THE-97 | PPG12 Stitched-Purity Closure",
            "RJ_CODEX_THREAD_ID": "019f-test-thread",
        }

    def _write_source_manifest(self) -> None:
        lanes: list[dict[str, str]] = []
        for photon in (5, 10, 20):
            for interaction in ("si", "di"):
                source_slug = f"photon{photon}_{interaction}"
                macro = self.lane_dir / f"{source_slug}.C"
                g4 = self.lane_dir / f"{source_slug}_g4.list"
                truth = self.lane_dir / f"{source_slug}_truth.list"
                macro_body = """INPUTREADHITS::listfile[0] = \"g4.list\";
// INPUTREADHITS::listfile[1] = \"calo.list\";
/* INPUTREADHITS::listfile[2] = \"cluster.list\"; */
// INPUTREADHITS::listfile[3] = \"mbd.list\";
INPUTREADHITS::listfile[4] = \"truth.list\";
TruthJetInput *truth_input = new TruthJetInput(Jet::PARTICLE);
"""
                if interaction == "di":
                    macro_body += "truth_input->add_embedding_flag(2);\n"
                macro.write_text(macro_body, encoding="utf-8")
                identities = [
                    f"pythia8_Photon{photon}_{interaction.upper()}_{index:06d}.root"
                    for index in range(5)
                ]
                g4.write_text(
                    "".join(f"G4Hits_{identity}\n" for identity in identities),
                    encoding="utf-8",
                )
                truth.write_text(
                    "".join(f"DST_TRUTH_JET_{identity}\n" for identity in identities),
                    encoding="utf-8",
                )
                for period in ("0mrad", "1p5mrad"):
                    lanes.append(
                        {
                            "lane_id": f"photon:photon{photon}:{period}:{interaction}",
                            "sample": f"Photon{photon}",
                            "period": period,
                            "interaction": interaction.upper(),
                            "ppg_macro": str(macro.resolve()),
                            "g4_full_list": str(g4.resolve()),
                            "truthjet_full_list": str(truth.resolve()),
                        }
                    )
        payload = {
            "schema": "ppg12-paired-source-manifest/v1",
            "common": self.common,
            "lanes": lanes,
        }
        self.source_manifest.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    def plan(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(DRIVER), "--condor-plan"],
            env=self.environment,
            text=True,
            capture_output=True,
            check=False,
        )

    @property
    def plan_path(self) -> Path:
        return self.evidence / "photon_condor_plan.json"


class PairedOracleCondorFanoutTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="paired-condor-test-")
        self.root = Path(self.temporary.name).resolve()
        self.fixture = Fixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_plan_is_deterministic_and_exactly_twelve_jobs(self) -> None:
        first = self.fixture.plan()
        self.assertEqual(first.returncode, 0, first.stderr or first.stdout)
        paths = [
            self.fixture.plan_path,
            self.fixture.evidence / "photon_condor_queue.tsv",
            self.fixture.evidence / "run_photon_condor_lane.sh",
            self.fixture.evidence / "photon_condor.submit",
            self.fixture.evidence / "photon_condor_generation.json",
        ]
        first_hashes = {path.name: digest(path) for path in paths}
        second = self.fixture.plan()
        self.assertEqual(second.returncode, 0, second.stderr or second.stdout)
        self.assertEqual(first_hashes, {path.name: digest(path) for path in paths})

        plan = json.loads(self.fixture.plan_path.read_text(encoding="utf-8"))
        self.assertEqual(plan["schema"], "ppg12-paired-oracle-condor-plan/v1")
        self.assertEqual(len(plan["authorization"]["lanes"]), 12)
        self.assertEqual(
            len({lane["lane_id"] for lane in plan["authorization"]["lanes"]}), 12
        )
        self.assertTrue(
            all(lane["first_five_event_identity_sha256"] for lane in plan["authorization"]["lanes"])
        )
        queue_rows = (self.fixture.evidence / "photon_condor_queue.tsv").read_text().splitlines()
        self.assertEqual(len(queue_rows), 12)
        self.assertFalse(any(row.startswith("lane_id") for row in queue_rows))
        submit = (self.fixture.evidence / "photon_condor.submit").read_text()
        self.assertIn(
            'environment = "RJ_CODEX_CHAT_NAME=THE-97 | PPG12 Stitched-Purity Closure;'
            'RJ_CODEX_THREAD_ID=019f-test-thread"',
            submit,
        )
        self.assertIn(
            '+RJ_CODEX_CHAT_NAME = "THE-97 | PPG12 Stitched-Purity Closure"',
            submit,
        )
        self.assertIn("RJ_CODEX_THREAD_ID=019f-test-thread", submit)
        self.assertIn("queue lane_id,lane_key from", submit)
        self.assertNotIn("max_retries", submit.lower())
        self.assertNotIn("merge", submit.lower())

    def test_submit_provenance_rejects_condor_injection_syntax(self) -> None:
        bad_titles = (
            "THE-97;EVIL=1",
            "THE-97$(Process)",
            'THE-97\"\nqueue 999',
            "THE-97\\escape",
        )
        for index, title in enumerate(bad_titles):
            with self.subTest(title=title):
                evidence = self.root / f"bad_evidence_{index}"
                environment = {
                    **self.fixture.environment,
                    "RJ_CODEX_CHAT_NAME": title,
                    "RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG": f"bad_title_{index}",
                    "RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR": str(evidence),
                    "RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT": str(
                        self.root / f"bad_output_{index}"
                    ),
                }
                completed = subprocess.run(
                    [str(DRIVER), "--condor-plan"],
                    env=environment,
                    text=True,
                    capture_output=True,
                    check=False,
                )
                self.assertEqual(completed.returncode, 2, completed.stdout)
                self.assertIn("Condor-safe title", completed.stderr)
                self.assertFalse((evidence / "photon_condor.submit").exists())

    def test_raw_thirteenth_duplicate_lane_is_not_collapsed(self) -> None:
        source = json.loads(
            self.fixture.source_manifest.read_text(encoding="utf-8")
        )
        source["lanes"].append(dict(source["lanes"][0]))
        self.fixture.source_manifest.write_text(
            json.dumps(source, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        completed = self.fixture.plan()
        self.assertEqual(completed.returncode, 2, completed.stdout)
        self.assertIn(
            "exactly 12 raw canonical lane rows; observed=13",
            completed.stderr,
        )
        self.assertFalse(self.fixture.plan_path.exists())

    def test_generated_wrapper_runs_lanes_and_audit_is_read_only(self) -> None:
        planned = self.fixture.plan()
        self.assertEqual(planned.returncode, 0, planned.stderr or planned.stdout)
        plan = json.loads(self.fixture.plan_path.read_text(encoding="utf-8"))
        wrapper = self.fixture.evidence / "run_photon_condor_lane.sh"
        for lane in plan["authorization"]["lanes"]:
            completed = subprocess.run(
                [str(wrapper), lane["lane_id"]],
                env=self.fixture.environment,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr or completed.stdout)

        before = {
            str(path.relative_to(self.root)): digest(path)
            for path in self.root.rglob("*")
            if path.is_file()
        }
        audited = subprocess.run(
            [str(DRIVER), "--audit-condor", "--plan", str(self.fixture.plan_path)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(audited.returncode, 0, audited.stderr or audited.stdout)
        report = json.loads(audited.stdout)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["counts"], {"FAIL": 0, "MISSING": 0, "PASS": 12})
        after = {
            str(path.relative_to(self.root)): digest(path)
            for path in self.root.rglob("*")
            if path.is_file()
        }
        self.assertEqual(before, after)

    def test_wrong_token_and_provenance_fail_before_lane_output(self) -> None:
        planned = self.fixture.plan()
        self.assertEqual(planned.returncode, 0, planned.stderr or planned.stdout)
        plan = json.loads(self.fixture.plan_path.read_text(encoding="utf-8"))
        lane = plan["authorization"]["lanes"][0]
        wrong_token = subprocess.run(
            [
                str(DRIVER), "--condor-lane", "--plan", str(self.fixture.plan_path),
                "--token", "RUN_PPG12_PAIRED_CONDOR_" + "0" * 64,
                "--lane-id", lane["lane_id"],
            ],
            env=self.fixture.environment,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(wrong_token.returncode, 2)
        self.assertFalse(Path(lane["output_base"]).exists())

        wrong_environment = dict(self.fixture.environment)
        wrong_environment["RJ_CODEX_THREAD_ID"] = "different-thread"
        wrong_provenance = subprocess.run(
            [
                str(DRIVER), "--condor-lane", "--plan", str(self.fixture.plan_path),
                "--token", plan["submission_token"], "--lane-id", lane["lane_id"],
            ],
            env=wrong_environment,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(wrong_provenance.returncode, 2)
        self.assertFalse(Path(lane["output_base"]).exists())

        Path(lane["output_base"]).mkdir(parents=True)
        audited = subprocess.run(
            [str(DRIVER), "--audit-condor", "--plan", str(self.fixture.plan_path)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(audited.returncode, 2)
        report = json.loads(audited.stdout)
        self.assertEqual(report["counts"]["FAIL"], 1)
        self.assertIn("without a valid completion receipt", report["lanes"][0]["reason"])

    def test_mock_submit_is_token_bound_and_cannot_repeat(self) -> None:
        planned = self.fixture.plan()
        self.assertEqual(planned.returncode, 0, planned.stderr or planned.stdout)
        plan = json.loads(self.fixture.plan_path.read_text(encoding="utf-8"))
        invocation = self.root / "submit_calls.txt"
        mock_submit = self.root / "condor_submit"
        mock_submit.write_text(
            "#!/usr/bin/env bash\nset -euo pipefail\n"
            f"printf '%s\\n' \"$*\" >> {invocation}\n"
            "printf '12345.0 - 12345.11\\n'\n",
            encoding="utf-8",
        )
        mock_submit.chmod(0o700)
        environment = {
            **self.fixture.environment,
            "RJ_PPG12_CONDOR_SUBMIT_COMMAND": str(mock_submit.resolve()),
        }
        command = [
            str(DRIVER), "--condor-submit", "--token", plan["submission_token"]
        ]
        submitted = subprocess.run(
            command, env=environment, text=True, capture_output=True, check=False
        )
        self.assertEqual(submitted.returncode, 0, submitted.stderr or submitted.stdout)
        self.assertEqual(len(invocation.read_text().splitlines()), 1)
        receipt = json.loads(
            (self.fixture.evidence / "photon_condor_submission.json").read_text()
        )
        self.assertEqual(receipt["authorized_job_count"], 12)
        repeated = subprocess.run(
            command, env=environment, text=True, capture_output=True, check=False
        )
        self.assertEqual(repeated.returncode, 2)
        self.assertEqual(len(invocation.read_text().splitlines()), 1)

    def test_source_drift_invalidates_plan_and_audit(self) -> None:
        planned = self.fixture.plan()
        self.assertEqual(planned.returncode, 0, planned.stderr or planned.stdout)
        source = json.loads(self.fixture.source_manifest.read_text(encoding="utf-8"))
        Path(source["lanes"][0]["g4_full_list"]).write_text(
            "".join(
                f"G4Hits_pythia8_Photon5_SI_CHANGED_{index:06d}.root\n"
                for index in range(5)
            ),
            encoding="utf-8",
        )
        audited = subprocess.run(
            ["python3", str(HELPER), "audit", "--plan", str(self.fixture.plan_path)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(audited.returncode, 2)
        self.assertIn("event identity mismatch", audited.stderr)


if __name__ == "__main__":
    unittest.main()
