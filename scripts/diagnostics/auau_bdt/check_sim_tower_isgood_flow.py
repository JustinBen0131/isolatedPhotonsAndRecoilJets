#!/usr/bin/env python3
"""Print the local SIM TowerInfo isGood/status flow from source semantics.

This is a no-ROOT, no-SDCC diagnostic.  It verifies the relevant local source
patterns, then prints the acceptance truth table implied by:

* TowerInfov4::get_isGood()
* CaloTowerCalib::copy_tower()
* RawClusterBuilderTemplate::IsAcceptableTower(TowerInfo*)
* the embedded-MC CaloTowerStatus skip contract
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
TOWERINFO = REPO / "coresoftware_local/offline/packages/CaloBase/TowerInfov4.h"
CALIB = REPO / "coresoftware_local/offline/packages/CaloReco/CaloTowerCalib.cc"
CLUSTER = REPO / "coresoftware_local/offline/packages/CaloReco/RawClusterBuilderTemplate.cc"
CALO_CALIB_MACRO = REPO / "macros/Calo_Calib.C"
AUAU_RECOILJETS = REPO / "src_AuAu/RecoilJets_AuAu.cc"


def require_contains(path: Path, needle: str, label: str) -> None:
    text = path.read_text()
    if needle not in text:
        raise SystemExit(f"FAIL {label}: did not find expected source pattern in {path}")
    print(f"PASS {label}: {path.relative_to(REPO)}")


def require_all_contains(path: Path, needles: list[str], label: str) -> None:
    text = path.read_text()
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise SystemExit(
            f"FAIL {label}: did not find expected source pattern(s) in {path}: {missing}"
        )
    print(f"PASS {label}: {path.relative_to(REPO)}")


@dataclass(frozen=True)
class TowerCase:
    label: str
    energy: float
    is_hot: bool = False
    is_bad_chi2: bool = False
    is_no_calib: bool = False
    is_not_instr: bool = False


def status_bits(case: TowerCase) -> int:
    status = 0
    status |= int(case.is_hot) << 0
    status |= int(case.is_bad_chi2) << 2
    status |= int(case.is_not_instr) << 3
    status |= int(case.is_no_calib) << 4
    return status


def get_is_good(case: TowerCase) -> bool:
    return not (case.is_hot or case.is_bad_chi2 or case.is_no_calib or case.is_not_instr)


def builder_accepts(case: TowerCase, min_tower_e: float = 0.070, apply_tower_selection: bool = True) -> bool:
    if case.energy < min_tower_e:
        return False
    if apply_tower_selection and not get_is_good(case):
        return False
    return True


def main() -> int:
    print("SIM_TOWER_ISGOOD_FLOW_CHECK_BEGIN")
    print("source-pattern checks")
    require_contains(
        TOWERINFO,
        "bool get_isGood() const override { return !(get_isHot() || get_isBadChi2() || get_isNoCalib() || get_isNotInstr()); }",
        "TowerInfov4::get_isGood definition",
    )
    require_contains(
        CALIB,
        "_calib_towers->get_tower_at_channel(channel)->copy_tower(caloinfo_raw);",
        "CaloTowerCalib copies raw status bits to calibrated TowerInfo",
    )
    require_contains(
        CLUSTER,
        "if (!tower->get_isGood())",
        "RawClusterBuilderTemplate rejects non-good TowerInfo when tower selection is on",
    )
    require_contains(
        CALO_CALIB_MACRO,
        "[Process_Calo_Calib][isSimEmbedded] skipping CaloTowerStatus setters",
        "embedded MC skips data-style CaloTowerStatus setters",
    )
    require_all_contains(
        AUAU_RECOILJETS,
        [
            'std::getenv("RJ_EVENT_CALO_REQUIRE_ISGOOD")',
            "if (requireEventCaloGoodTowers && !tower->get_isGood()) continue;",
            "if (!std::isfinite(e)) continue;",
            "sum += e;",
        ],
        "low-calo diagnostic has controlled get_isGood switch and finite signed tower sum",
    )

    print("\ntruth table from the verified source semantics")
    print(
        f"{'case':<30} {'E':>9} {'status':>8} {'isGood':>8} "
        f"{'clusterAccept':>14}  note"
    )

    cases = [
        TowerCase("good_positive", 0.120),
        TowerCase("hot_positive", 0.120, is_hot=True),
        TowerCase("badChi2_positive", 0.120, is_bad_chi2=True),
        TowerCase("noCalib_positive", 0.120, is_no_calib=True),
        TowerCase("notInstr_positive", 0.120, is_not_instr=True),
        TowerCase("good_below_cluster_threshold", 0.050),
        TowerCase("good_negative", -0.500),
    ]

    for case in cases:
        good = get_is_good(case)
        accepted = builder_accepts(case)
        if case.energy < 0.070:
            note = "rejected by RawClusterBuilderTemplate energy threshold"
        elif not good:
            note = "rejected by get_isGood status"
        else:
            note = "accepted"
        print(
            f"{case.label:<30} {case.energy:9.3f} {status_bits(case):8d} "
            f"{str(good):>8} {str(accepted):>14}  {note}"
        )

    diagnostic_sum = sum(case.energy for case in [TowerCase("good_negative", -0.5), TowerCase("good_positive", 1.0)] if get_is_good(case))
    diagnostic_sum_reject_bad = sum(case.energy for case in [TowerCase("bad_positive", 10.0, is_hot=True), TowerCase("good_positive", 1.0)] if get_is_good(case))
    print("\nlow-calo diagnostic sum sanity")
    print(f"good finite energies [-0.5, +1.0] sum to {diagnostic_sum:+.3f}; negative good towers remain included")
    print(f"hot/bad tower +10.0 plus good +1.0 sums to {diagnostic_sum_reject_bad:+.3f}; non-good towers are rejected")
    print("SIM_TOWER_ISGOOD_FLOW_CHECK_END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
