#!/usr/bin/env python3

from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[5]
PLANNER = (
    REPO
    / "scripts/sdcc/workflows/diagnostics/plan_ppg12_stitched_purity_inclusive_oracles.py"
)
SAMPLES = ("Jet8", "Jet12", "Jet20", "Jet30", "Jet40")
PERIODS = ("0mrad", "1p5mrad")
INTERACTIONS = ("SI", "DI")


class InclusiveOraclePlannerTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_context = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmp_context.name)
        self.assets = self.tmp / "assets"
        self.assets.mkdir()
        self.shared_args = {}
        for key in (
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
            path = self.assets / key
            path.write_text(f"sealed test asset {key}\n", encoding="utf-8")
            self.shared_args[key] = str(path)
        self.spec = self.make_spec()
        self.spec_path = self.tmp / "source_spec.json"
        self.manifest_path = self.tmp / "source_manifest.json"
        self.plan_path = self.tmp / "plan.json"
        self.driver = self.tmp / "paired-driver.sh"
        self.driver.write_text(
            "#!/usr/bin/env bash\n"
            "# PPG12_PAIRED_ORACLE_CAPABILITIES: photon,inclusive\n"
            "exit 0\n",
            encoding="utf-8",
        )
        self.driver.chmod(0o755)

    def tearDown(self) -> None:
        self.tmp_context.cleanup()

    def make_spec(self) -> dict:
        lanes = []
        for period in PERIODS:
            for sample in SAMPLES:
                jet = sample.removeprefix("Jet")
                for interaction in INTERACTIONS:
                    family = (
                        "js_pp200_signal_dual"
                        if interaction == "DI"
                        else "js_pp200_signal"
                    )
                    sample_key = f"run28_jet{jet}"
                    if interaction == "DI":
                        sample_key += "_double"
                    lane_dir = self.assets / f"{sample}_{period}_{interaction}"
                    lane_dir.mkdir(exist_ok=True)
                    g4 = lane_dir / "G4Hits.matched.list"
                    truth = lane_dir / "DST_JETS.matched.list"
                    macro = lane_dir / "Fun4All_run_sim.C"
                    event_stem = (
                        f"pythia8_Jet{jet}_pythia8_Detroit"
                        if interaction == "DI"
                        else f"pythia8_Jet{jet}"
                    )
                    g4.write_text(
                        "".join(
                            f"/store/{family}/g4hits/run0028/"
                            f"G4Hits_{event_stem}-{index:05d}.root\n"
                            for index in range(7)
                        ),
                        encoding="utf-8",
                    )
                    truth.write_text(
                        "".join(
                            f"/store/{family}/nopileup/jets/run0028/jet{jet}/"
                            f"DST_TRUTH_JET_{event_stem}-{index:05d}.root\n"
                            for index in range(7)
                        ),
                        encoding="utf-8",
                    )
                    macro_body = (
                        f"// {sample} {period} {interaction} frozen macro\n"
                        "INPUTREADHITS::listfile[0] = inputFile0;\n"
                        "INPUTREADHITS::listfile[4] = inputFile4;\n"
                    )
                    if interaction == "DI":
                        macro_body += (
                            "TruthJetInput *truth_input = "
                            "new TruthJetInput(Jet::PARTICLE);\n"
                            "truth_input->add_embedding_flag(2);\n"
                        )
                    macro.write_text(
                        macro_body,
                        encoding="utf-8",
                    )
                    driver_args = dict(self.shared_args)
                    driver_args.update(
                        {
                            "ppg_macro": str(macro),
                            "g4_full_list": str(g4),
                            "truthjet_full_list": str(truth),
                        }
                    )
                    lanes.append(
                        {
                            "lane_id": (
                                f"inclusive:{sample.lower()}:{period}:{interaction.lower()}"
                            ),
                            "sample": sample,
                            "period": period,
                            "interaction": interaction,
                            "sample_key": sample_key,
                            "source_family": family,
                            "driver_args": driver_args,
                        }
                    )
        return {
            "schema": "ppg12-inclusive-paired-source-spec/v1",
            "contract": {
                "rows_per_lane": 5,
                "jet5": "excluded",
                "global_scale": 1.0,
                "jet8_scale": 1.0,
                "external_scale": 1.0,
                "ownership_gate": {
                    "variable": "max_truth_jet_pt_r04",
                    "placement": "before_downstream_inclusive_physics_fills",
                    "lower_bound": "inclusive",
                    "upper_bound": "inclusive",
                    "zero_truth_jet": "reject",
                    "windows": {
                        "Jet8": {"lower_gev": 9.0, "upper_gev": 14.0},
                        "Jet12": {"lower_gev": 14.0, "upper_gev": 21.0},
                        "Jet20": {"lower_gev": 21.0, "upper_gev": 32.0},
                        "Jet30": {"lower_gev": 32.0, "upper_gev": 42.0},
                        "Jet40": {"lower_gev": 42.0, "upper_gev": 100.0},
                    },
                },
                "candidate_et_cap": {
                    "variable": "cluster_et",
                    "placement": "before_downstream_inclusive_candidate_physics_fills",
                    "upper_bound": "inclusive",
                    "caps_gev": {
                        "Jet8": 15.0,
                        "Jet12": 23.0,
                        "Jet20": 35.0,
                        "Jet30": 45.0,
                        "Jet40": 100.0,
                    },
                },
                "execution": "foreground_only_no_scheduler_no_merge_no_promotion",
            },
            "lanes": lanes,
        }

    def run_planner(self, *args: str, expect: int = 0) -> subprocess.CompletedProcess:
        result = subprocess.run(
            [sys.executable, str(PLANNER), *args],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            expect,
            msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def write_spec(self, spec: dict | None = None) -> None:
        self.spec_path.write_text(
            json.dumps(self.spec if spec is None else spec, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )

    def freeze(self, spec: dict | None = None) -> None:
        self.write_spec(spec)
        self.run_planner(
            "freeze",
            "--spec",
            str(self.spec_path),
            "--out",
            str(self.manifest_path),
        )

    def test_freeze_verify_and_plan_are_deterministic(self) -> None:
        self.freeze()
        first = self.manifest_path.read_bytes()
        second_manifest = self.tmp / "source_manifest_second.json"
        self.run_planner(
            "freeze",
            "--spec",
            str(self.spec_path),
            "--out",
            str(second_manifest),
        )
        self.assertEqual(first, second_manifest.read_bytes())
        self.run_planner("verify", "--manifest", str(self.manifest_path))

        output_root = self.tmp / "oracle-output"
        self.run_planner(
            "plan",
            "--manifest",
            str(self.manifest_path),
            "--driver",
            str(self.driver),
            "--output-root",
            str(output_root),
            "--out",
            str(self.plan_path),
            "--require-execution-ready",
        )
        plan_first = self.plan_path.read_bytes()
        second_plan = self.tmp / "plan_second.json"
        self.run_planner(
            "plan",
            "--manifest",
            str(self.manifest_path),
            "--driver",
            str(self.driver),
            "--output-root",
            str(output_root),
            "--out",
            str(second_plan),
            "--require-execution-ready",
        )
        self.assertEqual(plan_first, second_plan.read_bytes())
        plan = json.loads(plan_first)
        self.assertTrue(plan["execution_ready"])
        self.assertEqual(plan["lane_count"], 20)
        self.assertEqual(plan["scales"], {"external": 1.0, "global": 1.0, "jet8": 1.0})
        expected_ids = [
            f"inclusive:{sample.lower()}:{period}:{interaction.lower()}"
            for period in PERIODS
            for sample in SAMPLES
            for interaction in INTERACTIONS
        ]
        self.assertEqual([lane["lane_id"] for lane in plan["lanes"]], expected_ids)
        for lane in plan["lanes"]:
            self.assertNotIn("--run", lane["argv"])
            self.assertNotIn("--token", lane["argv"])
            sample = lane["lane_id"].split(":")[1].replace("jet", "Jet")
            self.assertEqual(
                lane["ownership_window"],
                self.spec["contract"]["ownership_gate"]["windows"][sample],
            )
            self.assertEqual(
                lane["candidate_et_cap_gev"],
                self.spec["contract"]["candidate_et_cap"]["caps_gev"][sample],
            )
            self.assertRegex(lane["first_five_event_identity_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(lane["ppg_macro_contract"]["active_input_indices"], [0, 4])

    def test_endpoint_zero_truth_and_candidate_cap_contract_is_frozen(self) -> None:
        self.freeze()
        manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        ownership = manifest["contract"]["ownership_gate"]
        self.assertEqual(ownership["lower_bound"], "inclusive")
        self.assertEqual(ownership["upper_bound"], "inclusive")
        self.assertEqual(ownership["zero_truth_jet"], "reject")
        self.assertEqual(
            ownership["windows"],
            {
                "Jet8": {"lower_gev": 9.0, "upper_gev": 14.0},
                "Jet12": {"lower_gev": 14.0, "upper_gev": 21.0},
                "Jet20": {"lower_gev": 21.0, "upper_gev": 32.0},
                "Jet30": {"lower_gev": 32.0, "upper_gev": 42.0},
                "Jet40": {"lower_gev": 42.0, "upper_gev": 100.0},
            },
        )
        candidate_cap = manifest["contract"]["candidate_et_cap"]
        self.assertEqual(candidate_cap["upper_bound"], "inclusive")
        self.assertEqual(
            candidate_cap["caps_gev"],
            {
                "Jet8": 15.0,
                "Jet12": 23.0,
                "Jet20": 35.0,
                "Jet30": 45.0,
                "Jet40": 100.0,
            },
        )

        mutations = (
            ("ownership_gate", "upper_bound", "exclusive"),
            ("ownership_gate", "zero_truth_jet", "accept"),
            ("candidate_et_cap", "upper_bound", "exclusive"),
        )
        for section, key, value in mutations:
            with self.subTest(section=section, key=key):
                spec = copy.deepcopy(self.spec)
                spec["contract"][section][key] = value
                self.write_spec(spec)
                self.run_planner(
                    "freeze",
                    "--spec",
                    str(self.spec_path),
                    "--out",
                    str(self.manifest_path),
                    expect=2,
                )

    def test_rowwise_event_identity_mismatch_and_swap_fail(self) -> None:
        for mutation in ("identity", "swap"):
            with self.subTest(mutation=mutation):
                spec = copy.deepcopy(self.spec)
                lane = next(
                    row
                    for row in spec["lanes"]
                    if row["lane_id"] == "inclusive:jet12:0mrad:si"
                )
                truth = Path(lane["driver_args"]["truthjet_full_list"])
                rows = truth.read_text(encoding="utf-8").splitlines()
                if mutation == "identity":
                    rows[2] = rows[2].replace("-00002.root", "-99999.root")
                else:
                    rows[1], rows[2] = rows[2], rows[1]
                truth.write_text("\n".join(rows) + "\n", encoding="utf-8")
                self.write_spec(spec)
                result = self.run_planner(
                    "freeze",
                    "--spec",
                    str(self.spec_path),
                    "--out",
                    str(self.manifest_path),
                    expect=2,
                )
                self.assertIn("event identity mismatch", result.stderr)
                self.spec = self.make_spec()

    def test_macro_contract_rejects_active_auxiliary_input(self) -> None:
        spec = copy.deepcopy(self.spec)
        lane = next(
            row
            for row in spec["lanes"]
            if row["lane_id"] == "inclusive:jet8:0mrad:si"
        )
        macro = Path(lane["driver_args"]["ppg_macro"])
        macro.write_text(
            macro.read_text(encoding="utf-8")
            + "INPUTREADHITS::listfile[1] = inputFile1;\n",
            encoding="utf-8",
        )
        self.write_spec(spec)
        result = self.run_planner(
            "freeze",
            "--spec",
            str(self.spec_path),
            "--out",
            str(self.manifest_path),
            expect=2,
        )
        self.assertIn("exactly {0, 4}", result.stderr)

    def test_macro_contract_requires_di_truth_mode_and_forbids_si_flag2(self) -> None:
        cases = ("missing_truth_input", "missing_flag2", "si_flag2")
        for mutation in cases:
            with self.subTest(mutation=mutation):
                spec = copy.deepcopy(self.spec)
                interaction = "SI" if mutation == "si_flag2" else "DI"
                lane = next(
                    row
                    for row in spec["lanes"]
                    if row["sample"] == "Jet20"
                    and row["period"] == "0mrad"
                    and row["interaction"] == interaction
                )
                macro = Path(lane["driver_args"]["ppg_macro"])
                text = macro.read_text(encoding="utf-8")
                if mutation == "missing_truth_input":
                    text = text.replace("new TruthJetInput", "existing_TruthJetInput")
                    expected_message = "reconstruct truth jets"
                elif mutation == "missing_flag2":
                    text = text.replace("add_embedding_flag(2)", "add_embedding_flag(1)")
                    expected_message = "add_embedding_flag(2)"
                else:
                    text += "truth_input->add_embedding_flag(2);\n"
                    expected_message = "must not actively call"
                macro.write_text(text, encoding="utf-8")
                self.write_spec(spec)
                result = self.run_planner(
                    "freeze",
                    "--spec",
                    str(self.spec_path),
                    "--out",
                    str(self.manifest_path),
                    expect=2,
                )
                self.assertIn(expected_message, result.stderr)
                self.spec = self.make_spec()

    def test_missing_duplicate_and_jet5_lanes_fail(self) -> None:
        cases = {}
        missing = copy.deepcopy(self.spec)
        missing["lanes"].pop()
        cases["missing"] = missing
        duplicate = copy.deepcopy(self.spec)
        duplicate["lanes"][-1] = copy.deepcopy(duplicate["lanes"][0])
        cases["duplicate"] = duplicate
        jet5 = copy.deepcopy(self.spec)
        jet5["lanes"][0]["sample"] = "Jet5"
        jet5["lanes"][0]["sample_key"] = "run28_jet5"
        jet5["lanes"][0]["lane_id"] = "inclusive:jet5:0mrad:si"
        cases["jet5"] = jet5
        for name, spec in cases.items():
            with self.subTest(name=name):
                self.write_spec(spec)
                self.run_planner(
                    "freeze",
                    "--spec",
                    str(self.spec_path),
                    "--out",
                    str(self.manifest_path),
                    expect=2,
                )

    def test_all_nonunit_scales_fail(self) -> None:
        for key, value in (
            ("global_scale", 0.99),
            ("jet8_scale", 2.3398),
            ("external_scale", 1.000001),
        ):
            with self.subTest(key=key):
                spec = copy.deepcopy(self.spec)
                spec["contract"][key] = value
                self.write_spec(spec)
                self.run_planner(
                    "freeze",
                    "--spec",
                    str(self.spec_path),
                    "--out",
                    str(self.manifest_path),
                    expect=2,
                )

    def test_lane_scale_field_is_rejected(self) -> None:
        spec = copy.deepcopy(self.spec)
        spec["lanes"][0]["external_scale"] = 1.0
        self.write_spec(spec)
        self.run_planner(
            "freeze",
            "--spec",
            str(self.spec_path),
            "--out",
            str(self.manifest_path),
            expect=2,
        )

    def test_si_di_source_mismatch_fails(self) -> None:
        for mutation in ("declared_family", "source_rows", "sample_key"):
            with self.subTest(mutation=mutation):
                spec = copy.deepcopy(self.spec)
                di_lane = next(
                    lane for lane in spec["lanes"] if lane["interaction"] == "DI"
                )
                if mutation == "declared_family":
                    di_lane["source_family"] = "js_pp200_signal"
                elif mutation == "sample_key":
                    di_lane["sample_key"] = di_lane["sample_key"].removesuffix(
                        "_double"
                    )
                else:
                    for key in ("g4_full_list", "truthjet_full_list"):
                        path = Path(di_lane["driver_args"][key])
                        path.write_text(
                            path.read_text(encoding="utf-8").replace(
                                "/js_pp200_signal_dual/", "/js_pp200_signal/"
                            ),
                            encoding="utf-8",
                        )
                self.write_spec(spec)
                self.run_planner(
                    "freeze",
                    "--spec",
                    str(self.spec_path),
                    "--out",
                    str(self.manifest_path),
                    expect=2,
                )
                self.spec = self.make_spec()

    def test_ownership_contract_mutation_fails(self) -> None:
        for mutation in ("ownership", "jet40_upper", "candidate_cap"):
            with self.subTest(mutation=mutation):
                spec = copy.deepcopy(self.spec)
                if mutation == "ownership":
                    spec["contract"]["ownership_gate"]["windows"]["Jet30"][
                        "upper_gev"
                    ] = 47.0
                elif mutation == "jet40_upper":
                    spec["contract"]["ownership_gate"]["windows"]["Jet40"][
                        "upper_gev"
                    ] = None
                else:
                    spec["contract"]["candidate_et_cap"]["caps_gev"]["Jet20"] = 36.0
                self.write_spec(spec)
                self.run_planner(
                    "freeze",
                    "--spec",
                    str(self.spec_path),
                    "--out",
                    str(self.manifest_path),
                    expect=2,
                )

    def test_stale_source_bytes_fail_verification(self) -> None:
        self.freeze()
        lane = self.spec["lanes"][0]
        g4 = Path(lane["driver_args"]["g4_full_list"])
        g4.write_text(g4.read_text(encoding="utf-8") + "# drift\n", encoding="utf-8")
        self.run_planner(
            "verify", "--manifest", str(self.manifest_path), expect=2
        )

    def test_driver_must_advertise_inclusive_execution(self) -> None:
        self.freeze()
        self.driver.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
        output_root = self.tmp / "blocked-output"
        self.run_planner(
            "plan",
            "--manifest",
            str(self.manifest_path),
            "--driver",
            str(self.driver),
            "--output-root",
            str(output_root),
            "--out",
            str(self.plan_path),
        )
        plan = json.loads(self.plan_path.read_text(encoding="utf-8"))
        self.assertFalse(plan["execution_ready"])
        self.assertEqual(plan["status"], "plan_only_blocked")
        self.run_planner(
            "plan",
            "--manifest",
            str(self.manifest_path),
            "--driver",
            str(self.driver),
            "--output-root",
            str(output_root),
            "--require-execution-ready",
            expect=2,
        )

    def test_existing_output_root_fails(self) -> None:
        self.freeze()
        output_root = self.tmp / "already-exists"
        output_root.mkdir()
        self.run_planner(
            "plan",
            "--manifest",
            str(self.manifest_path),
            "--driver",
            str(self.driver),
            "--output-root",
            str(output_root),
            expect=2,
        )


if __name__ == "__main__":
    unittest.main()
