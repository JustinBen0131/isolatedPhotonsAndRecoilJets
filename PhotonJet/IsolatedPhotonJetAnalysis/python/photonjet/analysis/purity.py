"""Event-leading ABCD region occupancy for PhotonJetTrees_v1."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import uproot


def purity_counts(
    input_paths: Iterable[Path],
    *,
    non_tight_definition: str,
    isolation_radius: float = 0.4,
) -> dict[str, Any]:
    """Return weighted event-leading A/B/C/D counts.

    Region boundaries are strict: isolated is below the isolated threshold and
    non-isolated is above the non-isolated threshold. Candidates exactly on a
    boundary are in neither region.
    """

    if non_tight_definition not in {"bounded", "complement"}:
        raise ValueError("non_tight_definition must be bounded or complement")
    if isolation_radius not in {0.3, 0.4}:
        raise ValueError("isolation_radius must be 0.3 or 0.4")
    radius = "r03" if isolation_radius == 0.3 else "r04"
    non_tight_branch = (
        "bdt_is_nontight"
        if non_tight_definition == "bounded"
        else "bdt_is_not_tight"
    )
    branches = [
        "source_file_index",
        "event_id_hi",
        "event_id_lo",
        "photon_encounter_ordinal",
        "photon_et",
        "event_weight",
        "bdt_is_tight",
        non_tight_branch,
        f"iso_{radius}",
        f"iso_{radius}_threshold",
        f"iso_{radius}_nonisolated_threshold",
    ]
    leaders: dict[tuple[int, int, int], list[tuple[float, int, float] | None]] = {}
    for raw_path in input_paths:
        path = Path(raw_path).resolve()
        with uproot.open(path) as root:
            if "photons" not in root:
                raise ValueError(f"{path} does not contain PhotonJetTrees_v1/photons")
            tree = root["photons"]
            missing = sorted(set(branches) - set(tree.keys()))
            if missing:
                raise ValueError(f"{path} photons lacks branches: {', '.join(missing)}")
            arrays = tree.arrays(branches, library="np")
        for index in range(len(arrays["photon_et"])):
            tight = bool(arrays["bdt_is_tight"][index])
            non_tight = bool(arrays[non_tight_branch][index])
            isolation = float(arrays[f"iso_{radius}"][index])
            isolated = isolation < float(arrays[f"iso_{radius}_threshold"][index])
            nonisolated = isolation > float(
                arrays[f"iso_{radius}_nonisolated_threshold"][index]
            )
            region = -1
            if tight and isolated:
                region = 0
            if tight and nonisolated:
                region = 1
            if non_tight and isolated:
                region = 2
            if non_tight and nonisolated:
                region = 3
            if region < 0:
                continue
            key = (
                int(arrays["source_file_index"][index]),
                int(arrays["event_id_hi"][index]),
                int(arrays["event_id_lo"][index]),
            )
            current = (
                float(arrays["photon_et"][index]),
                int(arrays["photon_encounter_ordinal"][index]),
                float(arrays["event_weight"][index]),
            )
            slots = leaders.setdefault(key, [None, None, None, None])
            previous = slots[region]
            if (
                previous is None
                or current[0] > previous[0]
                or (current[0] == previous[0] and current[1] < previous[1])
            ):
                slots[region] = current
    counts = np.zeros(4, dtype=float)
    events = np.zeros(4, dtype=int)
    for slots in leaders.values():
        for region, value in enumerate(slots):
            if value is not None:
                counts[region] += value[2]
                events[region] += 1
    return {
        "schema": "PhotonJetABCDCountsV1",
        "non_tight_definition": non_tight_definition,
        "isolation_radius": isolation_radius,
        "weighted_counts": dict(zip("ABCD", counts.tolist())),
        "event_counts": dict(zip("ABCD", events.tolist())),
    }


def raw_abcd_purity(counts: dict[str, float]) -> float:
    """Compute the uncorrected ABCD estimate for diagnostics.

    Physics results must apply the declared leakage and closure treatment
    before using this quantity.
    """

    a, b, c, d = (float(counts[key]) for key in "ABCD")
    if a <= 0 or d <= 0:
        raise ValueError("raw ABCD purity requires positive A and D")
    return 1.0 - (b * c) / (a * d)
