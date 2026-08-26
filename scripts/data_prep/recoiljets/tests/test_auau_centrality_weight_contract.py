from __future__ import annotations

import json
import math
from pathlib import Path
import tempfile
import unittest

from scripts.data_prep.recoiljets.auau_centrality_weight_contract import (
    BIN_EDGES_PERCENT,
    CANONICAL_CONTRACT_STATUS,
    CANONICAL_STITCH_COMPONENT,
    CENTRALITY_COMPONENT_PREFIX,
    DIAGNOSTIC_CONTRACT_STATUS,
    AuthorityValidationError,
    ContractBuildError,
    ContractNotReadyError,
    CurrentArtifactResolutionError,
    DoubleApplicationError,
    EventWeightState,
    StaleDependencyError,
    UnsupportedCentralityError,
    UnsupportedSampleFamilyError,
    WeightApplicationOrderError,
    build_canonical_contract,
    build_current_contract,
    ensure_canonical_contract,
    load_contract_receipt,
    sha256_file,
    validate_contract_receipt,
)


SAMPLE_KEYS = (
    "auau_data_merged",
    "auau_sim_photonjet_merged",
    "auau_sim_inclusivejet_merged",
)


class CurrentContractFixture:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.current = root / "current"
        self.root_paths = {}
        self.histograms = {}
        self.five_percent = {
            "auau_data_merged": tuple(float(index + 1) for index in range(16)),
            "auau_sim_photonjet_merged": tuple(float(17 - index) for index in range(16)),
            "auau_sim_inclusivejet_merged": tuple(2.0 for _ in range(16)),
        }
        self._make()

    @staticmethod
    def _one_percent_values(five_percent_counts):
        values = []
        for count in five_percent_counts:
            values.extend([float(count) / 5.0] * 5)
        values.extend([1.0] * 20)  # 80--100% is intentionally outside the contract.
        return tuple(values)

    def _make(self) -> None:
        histogram_keys = {
            "auau_data_merged": "MBD_NS_geq_2_vtx_lt_150/h_centrality",
            "auau_sim_photonjet_merged": "SIM/h_centrality",
            "auau_sim_inclusivejet_merged": "SIM/h_centrality",
        }
        edges = tuple(float(index) for index in range(101))
        for sample_key in SAMPLE_KEYS:
            root_path = self.root / "{}.root".format(sample_key)
            root_path.write_bytes(("fixture:{}:v1".format(sample_key)).encode("ascii"))
            self.root_paths[sample_key] = root_path
            self.histograms[(root_path, histogram_keys[sample_key])] = (
                self._one_percent_values(self.five_percent[sample_key]),
                edges,
            )
            pointer_dir = self.current / sample_key
            pointer_dir.mkdir(parents=True)
            (pointer_dir / "current.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "sample_key": sample_key,
                        "current_entry_id": "fixture-{}".format(sample_key),
                        "campaign_tag": "fixture-v1",
                        "canonical_status": "canonical",
                        "root_paths": [str(root_path)],
                    },
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )

    def loader(self, path: Path, histogram_key: str):
        return self.histograms[(path, histogram_key)]

    def build(self):
        return build_current_contract(self.current, histogram_loader=self.loader)

    @staticmethod
    def authorities(calibration_sha="3" * 64):
        roles = (
            "data_target",
            "sim_stitch",
            "centrality_calibration",
            "tree_selection",
            "derivation_code",
        )
        return {
            role: {
                "authority_id": "fixture-{}-v1".format(role),
                "sha256": calibration_sha if role == "centrality_calibration" else str(index + 4) * 64,
                "version": 1,
            }
            for index, role in enumerate(roles)
        }

    def build_canonical(self, calibration_sha="3" * 64):
        return build_canonical_contract(
            authorities=self.authorities(calibration_sha),
            data_counts=self.five_percent["auau_data_merged"],
            family_counts={
                "photonjet": self.five_percent["auau_sim_photonjet_merged"],
                "inclusivejet": self.five_percent["auau_sim_inclusivejet_merged"],
            },
        )


class AuAuCentralityWeightContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.fixture = CurrentContractFixture(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_builds_family_specific_maps_from_current_pointers(self) -> None:
        contract = self.fixture.build()
        self.assertEqual(contract.status, DIAGNOSTIC_CONTRACT_STATUS)
        self.assertEqual(len(contract.data_probability), 16)
        self.assertEqual(
            tuple(item.family for item in contract.family_maps),
            ("inclusivejet", "photonjet"),
        )
        self.assertNotEqual(
            contract.family_map("photonjet").weights,
            contract.family_map("inclusivejet").weights,
        )

        # DATA count is (i+1), while inclusive SIM is flat.  Both have mean 8.5,
        # so the inclusive map is exactly DATA/mean(DATA).
        self.assertAlmostEqual(contract.weight(0.0, "inclusivejet"), 1.0 / 8.5)
        self.assertAlmostEqual(contract.weight(79.999, "inclusivejet"), 16.0 / 8.5)

        payload = contract.to_payload()
        self.assertEqual(payload["bin_edges_percent"], list(BIN_EDGES_PERCENT))
        self.assertEqual(payload["numerical_lifetime"], "current artifact dependency fingerprint")
        self.assertEqual(
            {row["sample_key"] for row in payload["dependencies"]},
            set(SAMPLE_KEYS),
        )
        for dependency in payload["dependencies"]:
            self.assertEqual(len(dependency["pointer_sha256"]), 64)
            self.assertEqual(len(dependency["roots"][0]["sha256"]), 64)

    def test_fingerprint_changes_when_a_current_root_changes(self) -> None:
        before = self.fixture.build()
        changed = self.fixture.root_paths["auau_sim_photonjet_merged"]
        changed.write_bytes(b"fixture:auau_sim_photonjet_merged:v2")
        after = self.fixture.build()
        self.assertNotEqual(before.dependency_fingerprint, after.dependency_fingerprint)
        self.assertNotEqual(before.contract_fingerprint, after.contract_fingerprint)
        self.assertEqual(
            next(
                row for row in after.dependencies if row.sample_key == "auau_sim_photonjet_merged"
            ).roots[0].sha256,
            sha256_file(changed),
        )

    def test_diagnostic_current_histogram_contract_cannot_apply(self) -> None:
        contract = self.fixture.build()
        with self.assertRaises(ContractNotReadyError):
            contract.apply_once(
                EventWeightState.after_canonical_stitch(1.0),
                10.0,
                "photonjet",
                expected_dependency_fingerprint=contract.dependency_fingerprint,
            )

    def test_calibration_authority_changes_dependency_and_stale_binding_fails(self) -> None:
        before = self.fixture.build_canonical("a" * 64)
        after = self.fixture.build_canonical("b" * 64)
        self.assertNotEqual(before.dependency_fingerprint, after.dependency_fingerprint)
        self.assertNotEqual(before.contract_fingerprint, after.contract_fingerprint)
        with self.assertRaises(StaleDependencyError):
            after.apply_once(
                EventWeightState.after_canonical_stitch(1.0),
                10.0,
                "photonjet",
                expected_dependency_fingerprint=before.dependency_fingerprint,
            )

    def test_half_open_bins_are_exact(self) -> None:
        contract = self.fixture.build()
        self.assertEqual(contract.bin_index(0.0), 0)
        self.assertEqual(contract.bin_index(4.999999), 0)
        self.assertEqual(contract.bin_index(5.0), 1)
        self.assertEqual(contract.bin_index(79.999999), 15)
        for unsupported in (-0.000001, 80.0, math.nan, math.inf, "not-centrality"):
            with self.subTest(unsupported=unsupported):
                with self.assertRaises(UnsupportedCentralityError):
                    contract.bin_index(unsupported)

    def test_apply_once_requires_stitch_then_appends_hash_bound_component(self) -> None:
        contract = self.fixture.build_canonical()
        self.assertEqual(contract.status, CANONICAL_CONTRACT_STATUS)
        state = EventWeightState.after_canonical_stitch(
            2.0,
            prior_components=("producer_event_weight",),
        )
        result = contract.apply_once(
            state,
            5.0,
            "inclusivejet",
            expected_dependency_fingerprint=contract.dependency_fingerprint,
        )
        self.assertAlmostEqual(result.value, 2.0 * (2.0 / 8.5))
        self.assertEqual(result.components[-2], CANONICAL_STITCH_COMPONENT)
        self.assertEqual(
            result.components[-1],
            "{}inclusivejet:{}".format(
                CENTRALITY_COMPONENT_PREFIX,
                contract.contract_fingerprint,
            ),
        )

    def test_wrong_order_and_double_application_fail_closed(self) -> None:
        contract = self.fixture.build_canonical()
        with self.assertRaises(WeightApplicationOrderError):
            contract.apply_once(
                EventWeightState(1.0, ("producer_event_weight",)),
                10.0,
                "photonjet",
                expected_dependency_fingerprint=contract.dependency_fingerprint,
            )
        with self.assertRaises(WeightApplicationOrderError):
            contract.apply_once(
                EventWeightState(
                    1.0,
                    ("producer_event_weight", CANONICAL_STITCH_COMPONENT, "late_component"),
                ),
                10.0,
                "photonjet",
                expected_dependency_fingerprint=contract.dependency_fingerprint,
            )

        once = contract.apply_once(
            EventWeightState.after_canonical_stitch(1.0),
            10.0,
            "photonjet",
            expected_dependency_fingerprint=contract.dependency_fingerprint,
        )
        with self.assertRaises(DoubleApplicationError):
            contract.apply_once(
                once,
                10.0,
                "photonjet",
                expected_dependency_fingerprint=contract.dependency_fingerprint,
            )

    def test_non_embedded_family_fails_closed(self) -> None:
        contract = self.fixture.build_canonical()
        with self.assertRaises(UnsupportedSampleFamilyError):
            contract.apply_once(
                EventWeightState.after_canonical_stitch(1.0),
                10.0,
                "data",
                expected_dependency_fingerprint=contract.dependency_fingerprint,
            )

    def test_zero_support_fails_contract_construction(self) -> None:
        path = self.fixture.root_paths["auau_sim_inclusivejet_merged"]
        key = "SIM/h_centrality"
        values, edges = self.fixture.histograms[(path, key)]
        bad = list(values)
        bad[25:30] = [0.0] * 5
        self.fixture.histograms[(path, key)] = (tuple(bad), edges)
        with self.assertRaisesRegex(ContractBuildError, "zero or invalid support"):
            self.fixture.build()

    def test_noncanonical_current_pointer_is_rejected(self) -> None:
        pointer = self.fixture.current / "auau_data_merged/current.json"
        payload = json.loads(pointer.read_text(encoding="utf-8"))
        payload["canonical_status"] = "not_canonical"
        pointer.write_text(json.dumps(payload) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(CurrentArtifactResolutionError, "not canonical"):
            self.fixture.build()

    def test_canonical_builder_requires_every_authority_role(self) -> None:
        authorities = self.fixture.authorities()
        del authorities["centrality_calibration"]
        with self.assertRaisesRegex(AuthorityValidationError, "centrality_calibration"):
            build_canonical_contract(
                authorities=authorities,
                data_counts=self.fixture.five_percent["auau_data_merged"],
                family_counts={
                    "photonjet": self.fixture.five_percent[
                        "auau_sim_photonjet_merged"
                    ],
                    "inclusivejet": self.fixture.five_percent[
                        "auau_sim_inclusivejet_merged"
                    ],
                },
            )

    def test_ready_receipt_load_validates_hashes_and_expected_dependency(self) -> None:
        contract = self.fixture.build_canonical()
        receipt = Path(self.temporary.name) / "canonical_contract.json"
        contract.write_receipt(receipt)
        loaded = load_contract_receipt(
            receipt,
            expected_dependency_fingerprint=contract.dependency_fingerprint,
        )
        self.assertEqual(loaded, contract)
        self.assertEqual(
            validate_contract_receipt(
                contract.to_payload(),
                expected_dependency_fingerprint=contract.dependency_fingerprint,
            ),
            contract,
        )
        with self.assertRaises(StaleDependencyError):
            load_contract_receipt(
                receipt,
                expected_dependency_fingerprint="f" * 64,
            )

    def test_ensure_creates_then_reuses_the_exact_receipt(self) -> None:
        receipt = Path(self.temporary.name) / "resolved" / "canonical_contract.json"
        arguments = {
            "authorities": self.fixture.authorities(),
            "data_counts": self.fixture.five_percent["auau_data_merged"],
            "family_counts": {
                "photonjet": self.fixture.five_percent[
                    "auau_sim_photonjet_merged"
                ],
                "inclusivejet": self.fixture.five_percent[
                    "auau_sim_inclusivejet_merged"
                ],
            },
        }

        created = ensure_canonical_contract(receipt, **arguments)
        reused = ensure_canonical_contract(receipt, **arguments)

        self.assertEqual(created.action, "CREATED")
        self.assertEqual(reused.action, "REUSED_EXACT")
        self.assertEqual(created.contract, reused.contract)
        self.assertEqual(
            load_contract_receipt(
                receipt,
                expected_dependency_fingerprint=created.contract.dependency_fingerprint,
            ),
            created.contract,
        )

    def test_ensure_rebuilds_when_data_or_calibration_changes(self) -> None:
        receipt = Path(self.temporary.name) / "canonical_contract.json"
        family_counts = {
            "photonjet": self.fixture.five_percent["auau_sim_photonjet_merged"],
            "inclusivejet": self.fixture.five_percent[
                "auau_sim_inclusivejet_merged"
            ],
        }
        initial = ensure_canonical_contract(
            receipt,
            authorities=self.fixture.authorities("a" * 64),
            data_counts=self.fixture.five_percent["auau_data_merged"],
            family_counts=family_counts,
        )

        changed_data_counts = list(self.fixture.five_percent["auau_data_merged"])
        changed_data_counts[0] += 1.0
        data_rebuild = ensure_canonical_contract(
            receipt,
            authorities=self.fixture.authorities("a" * 64),
            data_counts=changed_data_counts,
            family_counts=family_counts,
        )
        calibration_rebuild = ensure_canonical_contract(
            receipt,
            authorities=self.fixture.authorities("b" * 64),
            data_counts=changed_data_counts,
            family_counts=family_counts,
        )

        self.assertEqual(initial.action, "CREATED")
        self.assertEqual(data_rebuild.action, "REBUILT_CHANGED_DEPENDENCY")
        self.assertEqual(calibration_rebuild.action, "REBUILT_CHANGED_DEPENDENCY")
        self.assertNotEqual(
            initial.contract.dependency_fingerprint,
            data_rebuild.contract.dependency_fingerprint,
        )
        self.assertNotEqual(
            data_rebuild.contract.dependency_fingerprint,
            calibration_rebuild.contract.dependency_fingerprint,
        )
        self.assertEqual(
            load_contract_receipt(
                receipt,
                expected_dependency_fingerprint=(
                    calibration_rebuild.contract.dependency_fingerprint
                ),
            ),
            calibration_rebuild.contract,
        )


if __name__ == "__main__":
    unittest.main()
