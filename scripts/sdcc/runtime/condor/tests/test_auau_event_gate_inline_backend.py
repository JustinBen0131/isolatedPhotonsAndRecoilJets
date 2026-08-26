from __future__ import annotations

import ast
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[5]
MACRO = REPO / "macros/Fun4All_recoilJets_unified_impl.C"
VALIDATOR = REPO / "scripts/sdcc/runtime/condor/validate_auau_event_gate_companion.py"
BUILDERS = (
    REPO / "scripts/data_prep/manifests/build_the236_schema10_data_production_packet.py",
    REPO / "scripts/data_prep/manifests/build_the236_resume_packet.py",
    REPO / "scripts/data_prep/manifests/build_the236_schema10_sparse_data_canary_packet.py",
)


class InlineAuAuEventGateBackendTest(unittest.TestCase):
    def test_macro_writes_authoritative_gate_in_base_event_loop(self) -> None:
        source = MACRO.read_text(encoding="utf-8")
        self.assertIn("class Writer final : public SubsysReco", source)
        self.assertIn('TTree("AuAuEventGateV1"', source)
        self.assertIn('TTree("AuAuEventGateSummaryV1"', source)
        self.assertIn("getTriggerVector()", source)
        self.assertIn("getLiveVector()", source)
        self.assertIn("getScaledVector()", source)
        self.assertIn("std::uint64_t{1} << m_photon10_bit", source)
        self.assertIn("minimum_bias->isAuAuMinimumBias()", source)
        self.assertIn('"INLINE_SCHEMA10_AUAU_BASE_EVENT_LOOP_V1"', source)
        self.assertIn('std::string(profile) == "sparse_photon_analysis_v1"', source)
        self.assertIn("sparse AuAu DATA requires MinimumBiasClassifier", source)
        self.assertIn("sparse AuAu DATA requires absolute RJ_AUAU_EVENT_GATE_OUTPUT_CANDIDATE", source)
        registration = source.index("new rj_auau_event_gate::Writer")
        recoil_registration = source.index("se->registerSubsystem(recoilJets);", registration)
        self.assertLess(registration, recoil_registration)

    def test_validator_requires_exact_witness_join(self) -> None:
        source = VALIDATOR.read_text(encoding="utf-8")
        ast.parse(source)
        for required in (
            "AuAuEventGateV1",
            "AuAuEventGateSummaryV1",
            "gl1_missing_events",
            "minimum_bias_missing_events",
            "trigger_witness_mismatches",
            "event_sequence",
            "raw_trigger_bits",
            "live_trigger_bits",
            "intersection != len(base_rows)",
            "BASE_ROOT_PUBLISHED_LAST_AFTER_COMPANION_V1",
        ):
            self.assertIn(required, source)

    def test_every_sparse_packet_builder_carries_mandatory_validator(self) -> None:
        for builder in BUILDERS:
            source = builder.read_text(encoding="utf-8")
            with self.subTest(builder=builder.name):
                self.assertIn("AuAuInlineEventGateProductionContractV1", source)
                self.assertIn('"required_for_sparse_auau_data": True', source)
                self.assertIn('"photon10_scaled_bit": 22', source)
                self.assertIn("EVENT_GATE_VALIDATOR.read_bytes()", source)
                self.assertIn('"auau_event_gate_required": True', source)


if __name__ == "__main__":
    unittest.main()
