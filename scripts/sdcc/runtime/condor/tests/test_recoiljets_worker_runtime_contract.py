#!/usr/bin/env python3
"""Mutation checks for the frozen worker release and Fun4All status contract."""

from __future__ import annotations

from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[5]
WRAPPER = REPO / "scripts/sdcc/runtime/condor/RecoilJets_Condor.sh"
CONTROLLER = (
    REPO
    / "scripts/sdcc/workflows/diagnostics/"
    "submit_the134_shower_factorial_canaries.sh"
)
MACRO = REPO / "macros/Fun4All_recoilJets_unified_impl.C"


def contract_failures(wrapper: str, controller: str, macro: str) -> list[str]:
    failures: list[str] = []

    wrapper_requirements = {
        "declared release name": 'runtime_expected_release_name="${RJ_PINNED_RELEASE_NAME:-}"',
        "declared release prefix": 'runtime_expected_offline_main="${RJ_PINNED_OFFLINE_MAIN:-}"',
        "versioned setup": (
            'source /opt/sphenix/core/bin/sphenix_setup.sh -n '
            '"$runtime_expected_release_name"'
        ),
        "resolved release assertion": (
            '"$runtime_resolved_offline_main" == "$runtime_expected_offline_main"'
        ),
        "pre-ROOT profile": "RECOILJETS_RUNTIME_PROFILE_V1",
    }
    for label, needle in wrapper_requirements.items():
        if needle not in wrapper:
            failures.append(f"wrapper lacks {label}")

    profile_at = wrapper.find("RECOILJETS_RUNTIME_PROFILE_V1")
    root_at = wrapper.find("# ------------------------ Run ROOT macro")
    if profile_at < 0 or root_at < 0 or profile_at >= root_at:
        failures.append("runtime profile is not emitted before ROOT")

    if controller.count("RJ_PINNED_RELEASE_NAME=ana.560") != 1:
        failures.append("THE-134 p+p lane lacks one ana.560 release declaration")
    pinned_prefix = (
        "RJ_PINNED_OFFLINE_MAIN=/cvmfs/sphenix.sdcc.bnl.gov/"
        "alma9.2-gcc-14.2.0/release/release_ana/ana.560"
    )
    if controller.count(pinned_prefix) != 1:
        failures.append("THE-134 p+p lane lacks one exact ana.560 prefix")

    macro_requirements = {
        "quiet-safe status record": "std::fprintf(\n      stderr,\n      \"RECOILJETS_FUN4ALL_STATUS_V1",
        "ROOT failure propagation": "if (gSystem) gSystem->Exit(90);",
        "run return capture": "runRc = se->run(nEvents);",
        "End return capture": "const int endRc = se->End();",
    }
    for label, needle in macro_requirements.items():
        if needle not in macro:
            failures.append(f"macro lacks {label}")
    if macro.count("detail::enforce_fun4all_status(") != 2:
        failures.append("both Fun4All execution paths are not status-enforced")

    return failures


class WorkerRuntimeContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.wrapper = WRAPPER.read_text(encoding="utf-8")
        cls.controller = CONTROLLER.read_text(encoding="utf-8")
        cls.macro = MACRO.read_text(encoding="utf-8")

    def test_repository_contract(self) -> None:
        self.assertEqual(
            contract_failures(self.wrapper, self.controller, self.macro),
            [],
        )

    def assert_mutation_rejected(
        self,
        *,
        wrapper: str | None = None,
        controller: str | None = None,
        macro: str | None = None,
    ) -> None:
        self.assertTrue(
            contract_failures(
                self.wrapper if wrapper is None else wrapper,
                self.controller if controller is None else controller,
                self.macro if macro is None else macro,
            )
        )

    def test_rejects_mutable_setup_mutation(self) -> None:
        self.assert_mutation_rejected(
            wrapper=self.wrapper.replace(
                'source /opt/sphenix/core/bin/sphenix_setup.sh -n '
                '"$runtime_expected_release_name"',
                "source /opt/sphenix/core/bin/sphenix_setup.sh -n",
                1,
            )
        )

    def test_rejects_release_drift_mutation(self) -> None:
        self.assert_mutation_rejected(
            controller=self.controller.replace(
                "RJ_PINNED_RELEASE_NAME=ana.560",
                "RJ_PINNED_RELEASE_NAME=ana.561",
                1,
            )
        )

    def test_rejects_late_profile_mutation(self) -> None:
        profile_line = next(
            line
            for line in self.wrapper.splitlines()
            if "RECOILJETS_RUNTIME_PROFILE_V1" in line
        )
        mutated = self.wrapper.replace(profile_line + "\n", "", 1)
        mutated = mutated.replace(
            "# ------------------------ Run ROOT macro -------------------",
            "# ------------------------ Run ROOT macro -------------------\n"
            + profile_line,
            1,
        )
        self.assert_mutation_rejected(wrapper=mutated)

    def test_rejects_hidden_fun4all_failure_mutation(self) -> None:
        self.assert_mutation_rejected(
            macro=self.macro.replace(
                "if (gSystem) gSystem->Exit(90);",
                "if (gSystem) gSystem->Exit(0);",
                1,
            )
        )

    def test_rejects_unchecked_execution_path_mutation(self) -> None:
        self.assert_mutation_rejected(
            macro=self.macro.replace(
                'detail::enforce_fun4all_status("scaled-trigger-only", runRc, endRc);',
                "",
                1,
            )
        )


if __name__ == "__main__":
    unittest.main()
