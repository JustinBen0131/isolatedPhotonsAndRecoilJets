from __future__ import annotations

import unittest

from photonjet.analysis.reduce import RecoilSelection
from photonjet.analysis.selection import compile_recoil_selection


class SelectionContractTest(unittest.TestCase):
    def test_boundaries_are_executed_by_the_serialized_predicates(self) -> None:
        program = compile_recoil_selection(RecoilSelection(region="inclusive"))
        self.assertTrue(program.accepts("photon", {"photon_et": 15.0, "photon_eta": 0.0}))
        self.assertFalse(program.accepts("photon", {"photon_et": 35.0, "photon_eta": 0.0}))
        self.assertFalse(program.accepts("photon", {"photon_et": 20.0, "photon_eta": 0.7}))
        self.assertFalse(
            program.accepts(
                "recoil",
                {"jet_pt": 5.0, "jet_eta": 0.0, "jet_radius": 0.4, "delta_phi": 3.0},
            )
        )
        self.assertTrue(
            program.accepts(
                "recoil",
                {"jet_pt": 5.1, "jet_eta": 0.1, "jet_radius": 0.4, "delta_phi": 3.0},
            )
        )

    def test_irrelevant_non_tight_definition_is_pruned(self) -> None:
        region_a = compile_recoil_selection(
            RecoilSelection(region="A", non_tight_definition="bounded")
        )
        region_c = compile_recoil_selection(
            RecoilSelection(region="C", non_tight_definition="complement")
        )
        a_keys = {fact.key for fact in region_a.facts}
        c_values = {fact.key: fact.value for fact in region_c.facts}
        self.assertNotIn("photon.non_tight_definition", a_keys)
        self.assertEqual(c_values["photon.non_tight_definition"], "complement")
        self.assertEqual(region_c.leader_branch, "leader_C_r04_complement_index")

    def test_every_cut_change_changes_the_program_identity(self) -> None:
        nominal = compile_recoil_selection(RecoilSelection())
        changed = compile_recoil_selection(RecoilSelection(jet_pt_min=6.0))
        self.assertNotEqual(nominal.sha256, changed.sha256)
        self.assertNotEqual(nominal.to_dict(), changed.to_dict())

    def test_abcd_predicates_execute_the_same_leader_state_that_is_serialized(self) -> None:
        program = compile_recoil_selection(RecoilSelection(region="A"))
        records = [
            {
                "photon_et": 19.0,
                "photon_encounter_ordinal": 0,
                "photon_bdt_score": 0.90,
                "photon_bdt_tight_threshold": 0.80,
                "photon_bdt_nontight_low_threshold": 0.50,
                "photon_bdt_nontight_high_threshold": 0.70,
                "photon_iso_r04": 1.0,
                "photon_iso_r04_threshold": 4.0,
                "photon_iso_r04_nonisolated_threshold": 7.0,
            },
            {
                "photon_et": 24.0,
                "photon_encounter_ordinal": 1,
                "photon_bdt_score": 0.91,
                "photon_bdt_tight_threshold": 0.80,
                "photon_bdt_nontight_low_threshold": 0.50,
                "photon_bdt_nontight_high_threshold": 0.70,
                "photon_iso_r04": 2.0,
                "photon_iso_r04_threshold": 4.0,
                "photon_iso_r04_nonisolated_threshold": 7.0,
            },
        ]
        self.assertEqual(program.choose_leader(records), 1)
        records[1]["photon_iso_r04"] = 8.0
        self.assertEqual(program.choose_leader(records), 0)


if __name__ == "__main__":
    unittest.main()
