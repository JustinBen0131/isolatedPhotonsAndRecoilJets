#!/usr/bin/env python3
"""Fast regression checks for the simulation candidate validator."""

from pathlib import Path
import subprocess
import unittest


WORKER = Path(__file__).resolve().parents[1] / "run_the121_the122_sim_production_row.sh"


class SimulationWorkerValidatorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = WORKER.read_text(encoding="utf-8")

    def test_shell_syntax(self) -> None:
        result = subprocess.run(["bash", "-n", str(WORKER)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_site_setup_cannot_fail_on_optional_unset_variables(self) -> None:
        setup = self.source.index('source /opt/sphenix/core/bin/sphenix_setup.sh -n "$release"')
        disable = self.source.rfind("set +e; set +u", 0, setup)
        restore = self.source.index("set -u; set -e", setup)
        validator = self.source.index('python3 - "$candidate_path" <<\'PY\'', restore)
        self.assertLess(disable, setup)
        self.assertLess(setup, restore)
        self.assertLess(restore, validator)

    def test_validator_uses_pyroot_and_nested_schema10_directory(self) -> None:
        self.assertIn("import ROOT", self.source)
        self.assertNotIn("root -b -q -l -e '", self.source)
        self.assertIn('root_file.GetDirectory("ReplayFoundationV1")', self.source)
        self.assertIn('replay.Get("RJEventV1")', self.source)
        self.assertIn('replay.Get("rj_replay_complete")', self.source)
        self.assertIn('replay.Get("rj_replay_schema_version")', self.source)
        self.assertIn('events.InheritsFrom("TTree")', self.source)
        self.assertIn('complete.InheritsFrom("TNamed")', self.source)
        self.assertIn('version.InheritsFrom("TNamed")', self.source)

    def test_validation_precedes_measurement_and_publication(self) -> None:
        validator = self.source.index('python3 - "$candidate_path" <<\'PY\'')
        event_check = self.source.index("(( actual_events <= event_upper_bound ))")
        measurement = self.source.index('python3 - "$time_raw" "$measurement_candidate"')
        publication = self.source.index('/bin/mv -- "$candidate_path" "$output_path"')
        self.assertLess(validator, event_check)
        self.assertLess(event_check, measurement)
        self.assertLess(measurement, publication)


if __name__ == "__main__":
    unittest.main()
