#!/usr/bin/env python3
"""Regression checks for the deployed PPG12 pp-SIM cluster-ET contract."""

from __future__ import annotations

import math
import re
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[4]
RECOIL_CC = (REPO / "src" / "RecoilJets.cc").read_text(encoding="utf-8")
RECOIL_H = (REPO / "src" / "RecoilJets.h").read_text(encoding="utf-8")
MACRO = (REPO / "macros" / "Fun4All_recoilJets_unified_impl.C").read_text(
    encoding="utf-8"
)


def function_body(source: str, signature: str) -> str:
    start = source.index(signature)
    brace = source.index("{", start)
    depth = 0
    for pos in range(brace, len(source)):
        if source[pos] == "{":
            depth += 1
        elif source[pos] == "}":
            depth -= 1
            if depth == 0:
                return source[brace : pos + 1]
    raise AssertionError(f"unterminated function: {signature}")


def sigma_extra_gev(et: float) -> float:
    data = 0.15**2 / et + 0.05**2 / et**2 + 0.05**2
    mc = 0.185**2 / et + 0.0**2 / et**2 + 0.04**2
    return et * math.sqrt(data - mc) if data > mc else 0.0


class DeployedEtSmearContractTest(unittest.TestCase):
    def test_exact_nominal_widths(self) -> None:
        expected = {
            10.0: 0.0,
            12.0: 0.0,
            15.0: 0.170660482,
            20.0: 0.357770876,
            25.0: 0.521416340,
            30.0: 0.678785680,
            35.0: 0.833441660,
            40.0: 0.986661036,
        }
        for et, target in expected.items():
            self.assertAlmostEqual(sigma_extra_gev(et), target, places=8)

    def test_pp_sim_gate_and_kill_switch_are_explicit(self) -> None:
        self.assertIn('envFlag("RJ_PPG12_PHOTON_YIELD_ET_SMEAR", true)', RECOIL_CC)
        prepare = function_body(
            RECOIL_CC, "bool RecoilJets::preparePPG12PhotonYieldClusterEtCache()"
        )
        cuts = function_body(
            RECOIL_CC, "double RecoilJets::ppg12PhotonYieldClusterEtForCuts("
        )
        for body in (prepare, cuts):
            self.assertRegex(
                body,
                re.compile(
                    r"m_ppg12PhotonYieldEnabled\s*&&\s*m_isSim\s*&&\s*!m_isAuAu"
                ),
            )
            self.assertIn("m_ppg12PhotonYieldEtSmearEnabled", body)

    def test_one_stateful_draw_is_cached_and_response_reuses_it(self) -> None:
        prepare = function_body(
            RECOIL_CC, "bool RecoilJets::preparePPG12PhotonYieldClusterEtCache()"
        )
        response = function_body(
            RECOIL_CC, "double RecoilJets::ppg12PhotonYieldClusterEtForResponse("
        )
        self.assertIn("new TRandom3(42)", RECOIL_CC)
        self.assertEqual(prepare.count("->Gaus("), 1)
        self.assertIn("ppg12PhotonYieldTowerMasked(pho)", prepare)
        self.assertIn("m_ppg12PhotonYieldEtForCutsCache.push_back", prepare)
        self.assertNotIn("Gaus(", response)
        self.assertIn("ppg12PhotonYieldClusterEtForCuts(recoEt, candidateIndex)", response)

    def test_missing_cache_fails_closed(self) -> None:
        cuts = function_body(
            RECOIL_CC, "double RecoilJets::ppg12PhotonYieldClusterEtForCuts("
        )
        self.assertIn("m_ppg12PhotonYieldEtForCutsCacheEvent != event_count", cuts)
        self.assertIn("std::numeric_limits<double>::quiet_NaN()", cuts)
        self.assertIn("return Fun4AllReturnCodes::ABORTRUN", RECOIL_CC)

    def test_deployed_base_v3e_base_e_route(self) -> None:
        attach = function_body(
            RECOIL_CC, "void RecoilJets::attachVariantScoresToSSVars("
        )
        self.assertIn('"tight_bdt_score_base_e"', attach)
        self.assertIn("v.pt_gamma >= 8.0 && v.pt_gamma < 35.0", attach)
        self.assertIn('baseEModelFile.replace', MACRO)
        self.assertIn('"cluster_Et", "vertexz", "cluster_Eta", "e11_over_e33"', MACRO)
        self.assertIn(
            '"cluster_et1", "cluster_et2", "cluster_et3", "cluster_et4"', MACRO
        )
        self.assertIn("tight_bdt_score_base_e", RECOIL_H)


if __name__ == "__main__":
    unittest.main()

