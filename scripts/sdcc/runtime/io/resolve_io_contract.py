#!/usr/bin/env python3
"""Resolve canonical and legacy RecoilJets I/O paths for a campaign tag."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import PurePosixPath


DEFAULT_REMOTE_BASE = "/sphenix/u/{user}/scratch/thesisAnalysis"
DEFAULT_BULK_BASE = "/sphenix/tg/tg01/bulk/jbennett"


def user_name() -> str:
    return os.environ.get("USER") or "patsfan753"


def validate_token(value: str, label: str) -> None:
    bad_fragments = ("agent_context", ".codex", "codex", "THE-")
    if not value or value.strip() != value or any(ch.isspace() for ch in value):
        raise SystemExit(f"[ERROR] Unsafe {label}: {value!r}")
    if any(fragment in value for fragment in bad_fragments):
        raise SystemExit(f"[ERROR] Unsafe {label}: {value!r}")
    if value in {".", ".."} or "/" in value or "\\" in value:
        raise SystemExit(f"[ERROR] Unsafe {label}: {value!r}")


def classify_path(path: str) -> dict[str, str]:
    p = path.rstrip("/")
    name = PurePosixPath(p).name
    if "/thesisAna/simembedded_" in p or "/thesisAna/simembeddedinclusive_" in p:
        return {"classification": "legacy_flat_tg_bulk", "status": "warn"}
    if "/scratch/thesisAnalysis/output_" in p:
        return {"classification": "legacy_flat_checkout_merge", "status": "warn"}
    if name.startswith("output_") or name.startswith("tmp_") or name.endswith((".csv", ".txt")):
        return {"classification": "root_contract_clutter", "status": "fail"}
    if "/thesisAna/recoiljets/" in p or "/runs/recoiljets/" in p:
        return {"classification": "canonical_recoiljets_io", "status": "ok"}
    return {"classification": "unclassified", "status": "unknown"}


def resolve(args: argparse.Namespace) -> dict[str, object]:
    validate_token(args.campaign, "campaign")
    validate_token(args.family, "family")
    user = user_name()
    remote_base = os.environ.get("RJ_REMOTE_BASE", DEFAULT_REMOTE_BASE.format(user=user))
    bulk_base = os.environ.get("RJ_BULK_BASE", DEFAULT_BULK_BASE)
    thesis_ana = os.environ.get("RJ_THESIS_ANA_ROOT", f"{bulk_base}/thesisAna")
    recoiljets_bulk = os.environ.get("RJ_RECOILJETS_BULK_ROOT", f"{thesis_ana}/recoiljets")
    run_root = os.environ.get("RJ_RECOILJETS_RUN_ROOT", f"{remote_base}/runs/recoiljets/current")

    campaign = args.campaign
    family = args.family
    canonical = {
        "sdcc_merge_root": f"{run_root}/{campaign}",
        "tg_bulk_signal_root": f"{recoiljets_bulk}/{family}/{campaign}/signal",
        "tg_bulk_background_root": f"{recoiljets_bulk}/{family}/{campaign}/background",
        "local_signal_pull": "InputFiles/simEmbedded",
        "local_background_pull": "InputFiles/InclusiveJetSIM_EMBEDDED",
        "local_evidence_root": f"dataOutput/{family}/{campaign}",
    }
    legacy = {
        "sdcc_merge_root": f"{remote_base}/output_{campaign}",
        "tg_bulk_signal_root": f"{thesis_ana}/simembedded_{campaign}",
        "tg_bulk_background_root": f"{thesis_ana}/simembeddedinclusive_{campaign}",
    }
    return {
        "campaign": campaign,
        "family": family,
        "canonical": canonical,
        "legacy_read_compatibility": legacy,
        "compatibility_policy": "write canonical for new runs; read legacy when auditing or pulling old runs",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", help="Campaign tag, not a path.")
    parser.add_argument(
        "--family",
        default="target_wp",
        help="Campaign family such as target_wp, width_study, stacking, stitching, or manual.",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON.")
    parser.add_argument("--check-path", help="Classify a concrete path instead of resolving only.")
    args = parser.parse_args()

    resolved = resolve(args)
    if args.check_path:
        resolved["path_check"] = {"path": args.check_path, **classify_path(args.check_path)}

    if args.json:
        print(json.dumps(resolved, indent=2, sort_keys=True))
    else:
        print(f"campaign\t{resolved['campaign']}")
        print(f"family\t{resolved['family']}")
        for group in ("canonical", "legacy_read_compatibility"):
            for key, value in resolved[group].items():
                print(f"{group}.{key}\t{value}")
        if "path_check" in resolved:
            check = resolved["path_check"]
            print(f"path_check.classification\t{check['classification']}")
            print(f"path_check.status\t{check['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
