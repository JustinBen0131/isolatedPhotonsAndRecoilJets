#!/usr/bin/env python3
"""Structural contract for detector-role-explicit PhotonClusterBuilder maps.

The ana.560-linked exhaustive channel round-trip/equivalence test belongs to
test_build_the134_ana560_calo_reco_runtime.sh.  This test instead pins every
PhotonClusterBuilder role guard and detector-sensitive call site.
"""

from __future__ import annotations

from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[5]
SOURCE = REPO_ROOT / "src" / "PhotonClusterBuilder.cc"


def validate(source: str) -> None:
    errors: list[str] = []
    normalized = re.sub(r"\s+", " ", source)

    if "->encode_key(" in source:
        errors.append("detector-dependent TowerInfoContainer::encode_key remains")
    if "->get_tower_at_key(" in source:
        errors.append("detector-dependent TowerInfoContainer::get_tower_at_key remains")

    if source.count("explicit_tower_key(") != 6:
        errors.append("expected one explicit_tower_key definition and five call sites")
    if source.count("tower_at_explicit_key(") != 7:
        errors.append("expected one tower_at_explicit_key definition and six call sites")
    if source.count("validate_tower_container_role(") != 7:
        errors.append("expected one role validator definition and six role checks")

    required_fragments = (
        "return TowerInfoDefs::encode_emcal(channel);",
        "return TowerInfoDefs::encode_hcal(channel);",
        "TowerInfoDefs::decode_emcal(tower_key)",
        "TowerInfoDefs::decode_hcal(tower_key)",
        "actual_size != expected_size",
        "actual_detector != expected_detector &&",
        "actual_detector != TowerInfoContainer::DETECTOR_INVALID",
        "actual_detector == TowerInfoContainer::DETECTOR_INVALID",
        "[tower-map][detector-explicit]",
        "TowerInfoContainer::EMCAL, 24576U",
        "TowerInfoContainer::HCAL, 1536U",
        "contradictory tower detector identity",
        "unsupported detector role for tower-key encoding",
        "unsupported detector role for tower-key decoding",
    )
    for fragment in required_fragments:
        if fragment not in source:
            errors.append(f"missing detector mapping contract: {fragment}")

    exact_role_guards = {
        "primary CEMC": (
            "Name(), m_emc_tower_node, m_emc_tower_container, "
            "TowerInfoContainer::EMCAL, 24576U"
        ),
        "primary HCALIN": (
            "Name(), m_ihcal_tower_node, m_ihcal_tower_container, "
            "TowerInfoContainer::HCAL, 1536U"
        ),
        "primary HCALOUT": (
            "Name(), m_ohcal_tower_node, m_ohcal_tower_container, "
            "TowerInfoContainer::HCAL, 1536U"
        ),
        "SUB1 CEMC retower": (
            "Name(), cemcIsoNode, m_emc_tower_container_iso, "
            "TowerInfoContainer::HCAL, 1536U"
        ),
        "SUB1 HCALIN": (
            "Name(), ihcalIsoNode, m_ihcal_tower_container_iso, "
            "TowerInfoContainer::HCAL, 1536U"
        ),
        "SUB1 HCALOUT": (
            "Name(), ohcalIsoNode, m_ohcal_tower_container_iso, "
            "TowerInfoContainer::HCAL, 1536U"
        ),
    }
    for role, fragment in exact_role_guards.items():
        count = normalized.count(fragment)
        if count != 1:
            errors.append(
                f"{role} role guard count differs: expected 1, found {count}"
            )

    exact_mapping_calls = {
        "CEMC keyed lookup": (
            "tower_at_explicit_key( m_emc_tower_container, towerinfokey, "
            "RawTowerDefs::CalorimeterId::CEMC)",
            4,
        ),
        "HCALIN keyed lookup": (
            "tower_at_explicit_key( m_ihcal_tower_container, towerinfokey, "
            "RawTowerDefs::CalorimeterId::HCALIN)",
            1,
        ),
        "HCALOUT keyed lookup": (
            "tower_at_explicit_key( m_ohcal_tower_container, towerinfokey, "
            "RawTowerDefs::CalorimeterId::HCALOUT)",
            1,
        ),
        "closest-HCAL channel map": (
            "explicit_tower_key( channel, isihcal ? "
            "RawTowerDefs::CalorimeterId::HCALIN : "
            "RawTowerDefs::CalorimeterId::HCALOUT)",
            1,
        ),
        "role-forwarded channel map": (
            "explicit_tower_key(channel, calo_id)",
            2,
        ),
        "CEMC channel map": (
            "explicit_tower_key(channel, RawTowerDefs::CalorimeterId::CEMC)",
            1,
        ),
        "CEMC QA channel map": (
            "explicit_tower_key(ch, RawTowerDefs::CalorimeterId::CEMC)",
            1,
        ),
    }
    for callsite, (fragment, expected) in exact_mapping_calls.items():
        count = normalized.count(fragment)
        if count != expected:
            errors.append(
                f"{callsite} count differs: expected {expected}, found {count}"
            )

    if errors:
        raise AssertionError("; ".join(errors))


def require_rejected(name: str, mutated: str) -> None:
    try:
        validate(mutated)
    except AssertionError:
        return
    raise AssertionError(f"mutation unexpectedly passed: {name}")


def main() -> int:
    source = SOURCE.read_text(encoding="utf-8")
    validate(source)

    mutations = {
        "generic_encode_reintroduced": source.replace(
            "explicit_tower_key(channel, RawTowerDefs::CalorimeterId::CEMC)",
            "m_emc_tower_container->encode_key(channel)",
            1,
        ),
        "detector_dependent_lookup_reintroduced": source.replace(
            "tower_at_explicit_key(\n"
            "            m_emc_tower_container, towerinfokey,\n"
            "            RawTowerDefs::CalorimeterId::CEMC)",
            "m_emc_tower_container->get_tower_at_key(towerinfokey)",
            1,
        ),
        "wrong_valid_detector_accepted": source.replace(
            "actual_detector != TowerInfoContainer::DETECTOR_INVALID",
            "false",
            1,
        ),
        "role_size_gate_removed": source.replace(
            "actual_size != expected_size",
            "false",
            1,
        ),
        "hcal_mapping_collapsed": source.replace(
            "return TowerInfoDefs::encode_hcal(channel);",
            "return 0U;",
            1,
        ),
        "invalid_detector_diagnostic_removed": source.replace(
            "[tower-map][detector-explicit]",
            "[tower-map]",
            1,
        ),
        "wrong_ohcal_sub1_role_and_size": source.replace(
            "Name(), ohcalIsoNode, m_ohcal_tower_container_iso,\n"
            "              TowerInfoContainer::HCAL, 1536U",
            "Name(), ohcalIsoNode, m_ohcal_tower_container_iso,\n"
            "              TowerInfoContainer::EMCAL, 24576U",
            1,
        ),
        "wrong_ihcal_closest_tower_role": source.replace(
            "isihcal ? RawTowerDefs::CalorimeterId::HCALIN\n"
            "                : RawTowerDefs::CalorimeterId::HCALOUT",
            "isihcal ? RawTowerDefs::CalorimeterId::CEMC\n"
            "                : RawTowerDefs::CalorimeterId::HCALOUT",
            1,
        ),
        "wrong_ihcal_shower_lookup_decoder": source.replace(
            "m_ihcal_tower_container, towerinfokey,\n"
            "                RawTowerDefs::CalorimeterId::HCALIN",
            "m_ihcal_tower_container, towerinfokey,\n"
            "                RawTowerDefs::CalorimeterId::CEMC",
            1,
        ),
    }
    for name, mutated in mutations.items():
        if mutated == source:
            raise AssertionError(f"mutation fixture did not alter source: {name}")
        require_rejected(name, mutated)

    print(
        "PASS: PhotonClusterBuilder detector-role mapping contract "
        f"source={SOURCE} mutations={len(mutations)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
