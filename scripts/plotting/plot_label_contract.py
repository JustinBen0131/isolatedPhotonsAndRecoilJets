#!/usr/bin/env python3
"""Fail-closed canonical metadata contract for collaborator-facing plots."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import re
import sys
from typing import Any, Iterable, MutableMapping


_DATASET_ROLES = {"training set", "validation set", "held-out set"}
_SAMPLE_KINDS = {"data", "simulation", "mixed"}
_DIAGNOSTIC_PLOT_KINDS = {"diagnostic", "bdt_qa"}
_DOWNSTREAM_SCIENCE_KIND_TOKENS = {
    "abcd",
    "closure",
    "covariance",
    "iteration",
    "leakage",
    "occupancy",
    "purity",
    "recoil",
    "refold",
    "response",
    "unfold",
}
FULL_ACCEPTED_SIMULATION_SCOPE = "full accepted simulation"
_PARTITION_LABEL_PATTERN = re.compile(
    r"\b(?:training|validation|held[- ]?out|holdout|buckets?)\b",
    re.IGNORECASE,
)

_SCRIPTS_ROOT = Path(__file__).resolve().parents[1]
if str(_SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_ROOT))

from data_prep.recoiljets.auau_centrality_weight_contract import (  # noqa: E402
    CANONICAL_CONTRACT_STATUS,
    load_contract_receipt,
)


def _validate_analysis_weight_receipt(
    receipt_path: str,
    expected_dependency_fingerprint: str,
) -> dict[str, object]:
    """Require an exact receipt from the current Au+Au weight provider."""
    path = Path(receipt_path)
    contract = load_contract_receipt(
        path,
        expected_dependency_fingerprint=expected_dependency_fingerprint,
        require_ready=True,
        verify_dependency_files=False,
    )
    if contract.status != CANONICAL_CONTRACT_STATUS:
        raise ValueError(
            "analysis-weight receipt is not ready for canonical Au+Au simulation application"
        )
    return {
        "path": str(path.resolve()),
        "expected_dependency_fingerprint": expected_dependency_fingerprint,
        "receipt": contract.to_payload(),
    }


@dataclass(frozen=True)
class PlotLabelContract:
    system: str
    energy_label: str
    sample_label: str
    sample_kind: str
    cut_lines: tuple[str, ...]
    input_paths: tuple[str, ...]
    centrality_label: str | None = None
    plot_kind: str = "physics"
    dataset_role: str | None = None
    simulation_scope: str | None = None
    analysis_weight_receipt_path: str | None = None
    analysis_weight_dependency_fingerprint: str | None = None
    diagnostic_weight_exception: str | None = None

    def apply_to_audit(self, audit: MutableMapping[str, Any]) -> dict[str, object]:
        """Validate this label contract and stamp its canonical audit fields.

        Nominal plot generators should use this method on the exact JSON object
        passed to ``codex_os_guard.py after-plot-render``.  The post-render
        guard independently requires these fields for every Au+Au SIM/mixed
        physics plot, so omitting this call cannot silently produce a nominal
        plot.
        """
        result = self.validate()
        audit.update(result["plot_audit_fields"])
        return result

    def validate(self) -> dict[str, object]:
        system = self.system.lower()
        if system not in {"pp", "auau"}:
            raise ValueError(f"unsupported plot system: {self.system!r}")
        if not self.energy_label.strip():
            raise ValueError("energy label is required")
        if system == "pp" and "s_{NN}" in self.energy_label:
            raise ValueError("pp plot cannot carry an s_NN energy label")
        if system == "auau":
            if "s_{NN}" not in self.energy_label:
                raise ValueError("Au+Au plot must carry an s_NN energy label")
            if not self.centrality_label or not self.centrality_label.strip():
                raise ValueError("Au+Au plot must carry an explicit centrality label")
        if not self.sample_label.strip():
            raise ValueError("sample label is required")
        sample_kind = self.sample_kind.lower()
        if sample_kind not in _SAMPLE_KINDS:
            raise ValueError(
                f"unsupported sample kind: {self.sample_kind!r}; expected data, simulation, or mixed"
            )
        if not self.cut_lines or any(not line.strip() for line in self.cut_lines):
            raise ValueError("at least one nonempty cut-definition line is required")
        if not self.input_paths:
            raise ValueError("at least one exact input path is required")
        missing = [path for path in self.input_paths if not Path(path).is_file()]
        if missing:
            raise FileNotFoundError(f"plot inputs do not exist: {missing}")
        plot_kind = self.plot_kind.strip().lower()
        dataset_role = None if self.dataset_role is None else self.dataset_role.strip().lower()
        simulation_scope = (
            None if self.simulation_scope is None else self.simulation_scope.strip().lower()
        )
        visible_text = "\n".join((self.sample_label, *self.cut_lines))
        if plot_kind == "bdt_qa":
            if dataset_role not in _DATASET_ROLES:
                raise ValueError(
                    "BDT QA plots must identify the sample as training set, validation set, or held-out set"
                )
            assert dataset_role is not None
            if dataset_role not in visible_text.lower():
                raise ValueError(
                    "BDT QA plots must visibly label the exact training set, validation set, or held-out set"
                )
            if simulation_scope is not None:
                raise ValueError(
                    "BDT QA split plots cannot claim the full accepted simulation science scope"
                )
        else:
            if dataset_role is not None:
                raise ValueError(
                    "dataset_role is allowed only for BDT-performance QA; downstream science must not inherit a training/validation/held-out split"
                )
            partition_match = _PARTITION_LABEL_PATTERN.search(visible_text)
            if partition_match:
                raise ValueError(
                    "non-BDT plot labels cannot carry training/validation/held-out/bucket partition language: {!r}".format(
                        partition_match.group(0)
                    )
                )

        # A downstream science plot cannot evade the full-universe contract by
        # appending or substituting "diagnostic" in its plot kind. The plain
        # diagnostic class remains available only for genuine sample/weight
        # diagnostics such as an unweighted centrality-shape comparison.
        is_downstream_science = (
            plot_kind not in _DIAGNOSTIC_PLOT_KINDS
            or any(token in plot_kind for token in _DOWNSTREAM_SCIENCE_KIND_TOKENS)
        )
        is_simulation_science = (
            is_downstream_science and sample_kind in {"simulation", "mixed"}
        )
        if is_simulation_science and simulation_scope != FULL_ACCEPTED_SIMULATION_SCOPE:
            raise ValueError(
                "downstream simulation science must bind simulation_scope={!r}; split-restricted simulation is BDT-QA-only".format(
                    FULL_ACCEPTED_SIMULATION_SCOPE
                )
            )
        if sample_kind == "data" and simulation_scope is not None:
            raise ValueError("data-only plots cannot carry a simulation_scope")

        has_receipt = self.analysis_weight_receipt_path is not None
        has_dependency_fingerprint = self.analysis_weight_dependency_fingerprint is not None
        has_exception = self.diagnostic_weight_exception is not None
        if has_receipt != has_dependency_fingerprint:
            raise ValueError(
                "analysis-weight receipt path and expected dependency fingerprint must be provided together"
            )
        if (has_receipt or has_dependency_fingerprint) and has_exception:
            raise ValueError(
                "analysis-weight receipt and diagnostic weight exception are mutually exclusive"
            )
        if has_exception:
            assert self.diagnostic_weight_exception is not None
            if not self.diagnostic_weight_exception.strip():
                raise ValueError("diagnostic weight exception must state a nonempty reason")
            if plot_kind not in _DIAGNOSTIC_PLOT_KINDS:
                raise ValueError(
                    "centrality-weight exceptions are allowed only for diagnostic plot kinds"
                )
            if system != "auau" or sample_kind not in {"simulation", "mixed"}:
                raise ValueError(
                    "centrality-weight diagnostic exception applies only to Au+Au simulation or mixed plots"
                )

        analysis_weight_provenance = None
        requires_auau_weight = system == "auau" and sample_kind in {"simulation", "mixed"}
        if requires_auau_weight and not has_receipt and not has_exception:
            raise ValueError(
                "Au+Au simulation or mixed plots require an exact canonical analysis-weight "
                "receipt or an explicit diagnostic-only exception"
            )
        if has_receipt:
            if not requires_auau_weight:
                raise ValueError(
                    "Au+Au centrality analysis-weight receipts apply only to simulation or mixed plots"
                )
            assert self.analysis_weight_receipt_path is not None
            assert self.analysis_weight_dependency_fingerprint is not None
            analysis_weight_provenance = _validate_analysis_weight_receipt(
                self.analysis_weight_receipt_path,
                self.analysis_weight_dependency_fingerprint,
            )
        elif has_exception:
            assert self.diagnostic_weight_exception is not None
            analysis_weight_provenance = {
                "diagnostic_exception": {
                    "reason": self.diagnostic_weight_exception.strip(),
                }
            }

        plot_audit_fields = {
            "collision_system": "Au+Au" if system == "auau" else "p+p",
            "sample_kind": sample_kind,
            "plot_kind": plot_kind,
            "sample_label": self.sample_label,
            "cut_lines": list(self.cut_lines),
            "dataset_role": dataset_role,
            "simulation_scope": simulation_scope,
            "analysis_weight_provenance": analysis_weight_provenance,
        }

        return {
            "status": "PASS",
            "contract": asdict(self),
            "analysis_weight_provenance": analysis_weight_provenance,
            "plot_audit_fields": plot_audit_fields,
        }


def canonical_cut_lines(lines: Iterable[str]) -> tuple[str, ...]:
    result = tuple(line.strip() for line in lines if line.strip())
    if not result:
        raise ValueError("empty cut-definition block")
    return result
