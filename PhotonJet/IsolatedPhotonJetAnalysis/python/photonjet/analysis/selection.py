"""Executable selection contracts shared by reducers and plot annotations.

The predicate objects in this module are deliberately small. Each one owns
the field, transformation, comparison, and value used by the event loop. The
same immutable objects are serialized into the histogram receipt and later
formatted by the plotting layer. A cut therefore cannot change without also
changing the contract from which its visible label is compiled.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from typing import Any, Mapping, Sequence

from photonjet.provenance import canonical_json


@dataclass(frozen=True)
class Predicate:
    """One executable numeric requirement in a named analysis stage."""

    key: str
    stage: str
    field: str
    operator: str
    value: float | None = None
    lower: float | None = None
    upper: float | None = None
    tolerance: float | None = None
    reference: str | None = None
    lower_reference: str | None = None
    upper_reference: str | None = None
    transform: str = "identity"
    unit: str | None = None

    def __post_init__(self) -> None:
        if self.transform not in {"identity", "abs"}:
            raise ValueError(f"unsupported transform: {self.transform}")
        if self.operator not in {
            "gt", "lt", "closed_open", "near",
            "gt_field", "lt_field", "open_between_fields", "not_gt_field",
        }:
            raise ValueError(f"unsupported operator: {self.operator}")
        if self.operator in {"gt", "lt", "near"} and self.value is None:
            raise ValueError(f"{self.operator} requires value")
        if self.operator == "closed_open":
            if self.lower is None or self.upper is None or not self.lower < self.upper:
                raise ValueError("closed_open requires an increasing interval")
        if self.operator == "near" and (self.tolerance is None or self.tolerance <= 0):
            raise ValueError("near requires a positive tolerance")
        if self.operator in {"gt_field", "lt_field", "not_gt_field"} and not self.reference:
            raise ValueError(f"{self.operator} requires a reference field")
        if self.operator == "open_between_fields" and not (
            self.lower_reference and self.upper_reference
        ):
            raise ValueError("open_between_fields requires lower and upper reference fields")
        for number in (self.value, self.lower, self.upper, self.tolerance):
            if number is not None and not math.isfinite(number):
                raise ValueError(f"non-finite predicate value for {self.key}")

    def evaluate(self, record: Mapping[str, float]) -> bool:
        if self.field not in record:
            raise KeyError(f"selection field is absent: {self.field}")
        observed = float(record[self.field])
        if not math.isfinite(observed):
            return False
        if self.transform == "abs":
            observed = abs(observed)
        if self.operator == "gt":
            return observed > float(self.value)
        if self.operator == "lt":
            return observed < float(self.value)
        if self.operator == "closed_open":
            return float(self.lower) <= observed < float(self.upper)
        if self.operator == "near":
            return abs(observed - float(self.value)) < float(self.tolerance)
        if self.operator in {"gt_field", "lt_field", "not_gt_field"}:
            if self.reference not in record:
                raise KeyError(f"selection reference field is absent: {self.reference}")
            reference = float(record[self.reference])
            if not math.isfinite(reference):
                return False
            if self.operator == "gt_field":
                return observed > reference
            if self.operator == "lt_field":
                return observed < reference
            return not observed > reference
        if self.operator == "open_between_fields":
            if self.lower_reference not in record or self.upper_reference not in record:
                raise KeyError("selection interval reference field is absent")
            lower = float(record[self.lower_reference])
            upper = float(record[self.upper_reference])
            if not math.isfinite(lower) or not math.isfinite(upper) or not lower < upper:
                return False
            return lower < observed < upper
        raise AssertionError(self.operator)

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "key": self.key,
            "stage": self.stage,
            "field": self.field,
            "operator": self.operator,
            "transform": self.transform,
        }
        for name in (
            "value", "lower", "upper", "tolerance", "reference",
            "lower_reference", "upper_reference", "unit",
        ):
            value = getattr(self, name)
            if value is not None:
                result[name] = value
        return result

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "Predicate":
        allowed = {
            "key", "stage", "field", "operator", "value", "lower", "upper",
            "tolerance", "reference", "lower_reference", "upper_reference",
            "transform", "unit",
        }
        extra = sorted(set(value) - allowed)
        if extra:
            raise ValueError(f"unknown predicate fields: {extra}")
        return cls(**dict(value))


@dataclass(frozen=True)
class SelectionFact:
    """A categorical choice that affects which rows enter the payload."""

    key: str
    value: str | float

    def to_dict(self) -> dict[str, Any]:
        return {"key": self.key, "value": self.value}

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SelectionFact":
        if set(value) != {"key", "value"}:
            raise ValueError("selection fact must contain exactly key and value")
        return cls(key=str(value["key"]), value=value["value"])


@dataclass(frozen=True)
class SelectionProgram:
    """The complete, executable lineage for one recoil payload."""

    predicates: tuple[Predicate, ...]
    facts: tuple[SelectionFact, ...]
    leader_branch: str | None

    def accepts(self, stage: str, record: Mapping[str, float]) -> bool:
        return all(
            predicate.evaluate(record)
            for predicate in self.predicates
            if predicate.stage == stage
        )

    def choose_leader(self, records: Sequence[Mapping[str, float]]) -> int:
        """Execute the serialized ABCD predicates and deterministic leader rule."""

        if self.leader_branch is None:
            raise ValueError("inclusive selection has no event-leading candidate")
        eligible = [
            index for index, record in enumerate(records)
            if self.accepts("photon_class", record)
        ]
        if not eligible:
            return -1
        return min(
            eligible,
            key=lambda index: (
                -float(records[index]["photon_et"]),
                int(records[index]["photon_encounter_ordinal"]),
                index,
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "PhotonJetSelectionProgramV1",
            "predicates": [predicate.to_dict() for predicate in self.predicates],
            "facts": [fact.to_dict() for fact in self.facts],
            "leader_branch": self.leader_branch,
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(canonical_json(self.to_dict()).encode("utf-8")).hexdigest()

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SelectionProgram":
        if set(value) != {"schema", "predicates", "facts", "leader_branch"}:
            raise ValueError("selection program fields differ from PhotonJetSelectionProgramV1")
        if value["schema"] != "PhotonJetSelectionProgramV1":
            raise ValueError("unsupported selection program schema")
        predicates = value["predicates"]
        facts = value["facts"]
        if not isinstance(predicates, list) or not isinstance(facts, list):
            raise ValueError("selection predicates and facts must be lists")
        leader = value["leader_branch"]
        if leader is not None and not isinstance(leader, str):
            raise ValueError("leader_branch must be a string or null")
        return cls(
            predicates=tuple(Predicate.from_dict(item) for item in predicates),
            facts=tuple(SelectionFact.from_dict(item) for item in facts),
            leader_branch=leader,
        )


def compile_recoil_selection(selection: Any) -> SelectionProgram:
    """Compile a RecoilSelection-like object into the sole cut program."""

    region = str(selection.region)
    if region not in {"inclusive", "A", "B", "C", "D"}:
        raise ValueError(f"unsupported photon region: {region}")
    non_tight = str(selection.non_tight_definition)
    if non_tight not in {"bounded", "complement"}:
        raise ValueError(f"unsupported non-tight definition: {non_tight}")

    predicates = (
        Predicate(
            key="photon.et",
            stage="photon",
            field="photon_et",
            operator="closed_open",
            lower=float(selection.photon_et_min),
            upper=float(selection.photon_et_max),
            unit="GeV",
        ),
        Predicate(
            key="photon.abs_eta",
            stage="photon",
            field="photon_eta",
            operator="lt",
            transform="abs",
            value=float(selection.photon_abs_eta_max),
        ),
        Predicate(
            key="jet.pt",
            stage="recoil",
            field="jet_pt",
            operator="gt",
            value=float(selection.jet_pt_min),
            unit="GeV",
        ),
        Predicate(
            key="jet.abs_eta",
            stage="recoil",
            field="jet_eta",
            operator="lt",
            transform="abs",
            value=float(selection.jet_abs_eta_max),
        ),
        Predicate(
            key="jet.radius",
            stage="recoil",
            field="jet_radius",
            operator="near",
            value=float(selection.jet_radius),
            tolerance=1.0e-9,
        ),
        Predicate(
            key="recoil.delta_phi",
            stage="recoil",
            field="delta_phi",
            operator="gt",
            value=float(selection.delta_phi_min),
            unit="rad",
        ),
    )

    if region == "inclusive":
        return SelectionProgram(
            predicates=predicates,
            facts=(SelectionFact("photon.candidate_scope", "inclusive_pairs"),),
            leader_branch=None,
        )

    identity = "tight" if region in {"A", "B"} else "non_tight"
    isolation = "isolated" if region in {"A", "C"} else "nonisolated"
    facts: list[SelectionFact] = [
        SelectionFact("photon.candidate_scope", "event_leading"),
        SelectionFact("photon.abcd_region", region),
        SelectionFact("photon.id_state", identity),
        SelectionFact("photon.isolation_state", isolation),
        SelectionFact("photon.isolation_radius", 0.4),
        SelectionFact("photon.model_view", "H70"),
        SelectionFact("photon.bdt_threshold_source", "per_candidate_calibrated_threshold"),
        SelectionFact("photon.isolation_threshold_source", "per_candidate_isolation_witness"),
        SelectionFact("photon.leader_rule", "highest_et_then_encounter_ordinal"),
    ]
    if identity == "non_tight":
        facts.append(SelectionFact("photon.non_tight_definition", non_tight))
    if identity == "tight":
        classification = Predicate(
            key="photon.bdt_class",
            stage="photon_class",
            field="photon_bdt_score",
            operator="gt_field",
            reference="photon_bdt_tight_threshold",
        )
    elif non_tight == "bounded":
        classification = Predicate(
            key="photon.bdt_class",
            stage="photon_class",
            field="photon_bdt_score",
            operator="open_between_fields",
            lower_reference="photon_bdt_nontight_low_threshold",
            upper_reference="photon_bdt_nontight_high_threshold",
        )
    else:
        classification = Predicate(
            key="photon.bdt_class",
            stage="photon_class",
            field="photon_bdt_score",
            operator="not_gt_field",
            reference="photon_bdt_tight_threshold",
        )
    isolation_predicate = Predicate(
        key="photon.isolation_class",
        stage="photon_class",
        field="photon_iso_r04",
        operator="lt_field" if isolation == "isolated" else "gt_field",
        reference=(
            "photon_iso_r04_threshold"
            if isolation == "isolated"
            else "photon_iso_r04_nonisolated_threshold"
        ),
        unit="GeV",
    )
    suffix = "_complement" if region in {"C", "D"} and non_tight == "complement" else ""
    return SelectionProgram(
        predicates=(*predicates, classification, isolation_predicate),
        facts=tuple(facts),
        leader_branch=f"leader_{region}_r04{suffix}_index",
    )
