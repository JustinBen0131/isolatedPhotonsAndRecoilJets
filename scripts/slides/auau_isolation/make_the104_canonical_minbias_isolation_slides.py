#!/usr/bin/env python3
"""Generate the eight canonical-MinBias AuAu isolation slide PNGs."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


THIS_FILE = Path(__file__).resolve()
REPO = next((parent for parent in THIS_FILE.parents if (parent / "AGENTS.md").exists()), THIS_FILE.parents[4])
FIT_SCRIPT = REPO / "scripts/plotting/efficiency/make_the96_auau_isolation_efficiency_fits.py"
FIT_SLIDE_SCRIPT = REPO / "scripts/slides/auau_isolation/make_the96_r03_baseline_isolation_fit_slides.py"
DIST_SLIDE_SCRIPT = REPO / "scripts/slides/auau_isolation/make_the96_isolation_centrality_distribution_slides.py"

CANONICAL_BDT = Path(
    "/sphenix/user/patsfan753/forBlair/auauBDTdefault/"
    "auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.root"
)
CANONICAL_WP80 = "T80(c)=0.5378806890+0.0011122896*c"


def run(*arguments: str | Path) -> None:
    subprocess.run([sys.executable, *(str(argument) for argument in arguments)], check=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_file(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return resolved


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-root", type=Path, required=True)
    parser.add_argument("--phosub-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--campaign-tag", default="the104_canonical_minbias_isolation_20260715")
    args = parser.parse_args()

    base_root = require_file(args.base_root)
    phosub_root = require_file(args.phosub_root)
    output_dir = args.output_dir.expanduser().resolve()
    derivation_dir = output_dir / "derivations"
    package_dir = output_dir / "slides_25_32"
    output_dir.mkdir(parents=True, exist_ok=True)
    package_dir.mkdir(parents=True, exist_ok=True)

    for cone, radius_tag in (("R30", "r03"), ("R40", "r04")):
        fit_dir = derivation_dir / f"{radius_tag}_baseline_nominal_15to35"
        fit_slide_dir = derivation_dir / f"{radius_tag}_baseline_fit_slides"
        run(
            FIT_SCRIPT,
            "--input", base_root,
            "--outdir", fit_dir,
            "--cone", cone,
            "--algorithm", "baseline",
            "--pt-scheme", "nominal",
            "--campaign", "THE-104",
            "--campaign-tag", args.campaign_tag,
            "--output-prefix", "the104",
        )
        run(
            FIT_SLIDE_SCRIPT,
            "--cone", cone,
            "--algorithm", "baseline",
            "--pt-scheme", "nominal",
            "--source-dir", fit_dir,
            "--output-dir", fit_slide_dir,
            "--source-prefix-root", "the104",
            "--campaign-label", "THE104",
        )

    baseline_dist_dir = derivation_dir / "baseline_centrality_distributions"
    phosub_dist_dir = derivation_dir / "phosub_centrality_distributions"
    common_subtitle = (
        "Embedded Photon+Jet 12+20 | truth matched | "
        "MinimumBiasClassifier required | 0.2 GeV rebinned counts"
    )
    run(
        DIST_SLIDE_SCRIPT,
        "--input-root", base_root,
        "--output-dir", baseline_dist_dir,
        "--variant-key", "baseline",
        "--title-subject", "baseline isolation counts",
        "--subtitle", common_subtitle,
        "--campaign-tag", args.campaign_tag,
        "--algorithm", "baseline tower-cone without topocluster isolation",
        "--manifest-schema", "THE104_CANONICAL_MINBIAS_BASELINE_ISOLATION_DISTRIBUTIONS_V1",
        "--production-note", "Complete THE-104 canonical-MinBias Photon12+20 merge.",
    )
    run(
        DIST_SLIDE_SCRIPT,
        "--input-root", phosub_root,
        "--output-dir", phosub_dist_dir,
        "--variant-key", "phosub",
        "--title-subject", "PHOSUB isolation counts",
        "--subtitle", common_subtitle,
        "--campaign-tag", args.campaign_tag,
        "--algorithm", "variantB EMCal photon-cluster energy subtraction",
        "--manifest-schema", "THE104_CANONICAL_MINBIAS_PHOSUB_ISOLATION_DISTRIBUTIONS_V1",
        "--production-note", "Complete THE-104 canonical-MinBias Photon12+20 merge.",
    )

    slide_sources = {
        25: derivation_dir / "r03_baseline_fit_slides/auau_r03_baseline_isolation_pt_flat_fits_slide.png",
        26: derivation_dir / "r03_baseline_fit_slides/auau_r03_baseline_isolation_centrality_fit_test_slide.png",
        27: derivation_dir / "r04_baseline_fit_slides/auau_r04_baseline_isolation_pt_flat_fits_slide.png",
        28: derivation_dir / "r04_baseline_fit_slides/auau_r04_baseline_isolation_centrality_fit_test_slide.png",
        29: baseline_dist_dir / "auau_r03_baseline_isolation_centrality_overlay_slide.png",
        30: baseline_dist_dir / "auau_r04_baseline_isolation_centrality_overlay_slide.png",
        31: phosub_dist_dir / "auau_r03_phosub_isolation_centrality_overlay_slide.png",
        32: phosub_dist_dir / "auau_r04_phosub_isolation_centrality_overlay_slide.png",
    }
    names = {
        25: "baseline_r03_flat_fits",
        26: "baseline_r03_centrality_fits",
        27: "baseline_r04_flat_fits",
        28: "baseline_r04_centrality_fits",
        29: "baseline_r03_distributions",
        30: "baseline_r04_distributions",
        31: "phosub_r03_distributions",
        32: "phosub_r04_distributions",
    }
    packaged: list[dict[str, object]] = []
    for slide_number, source in slide_sources.items():
        source = require_file(source)
        destination = package_dir / f"slide_{slide_number}_{names[slide_number]}.png"
        shutil.copy2(source, destination)
        packaged.append(
            {
                "slide_number": slide_number,
                "png": str(destination),
                "source_png": str(source),
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
            }
        )

    manifest = output_dir / "the104_canonical_minbias_isolation_slides_25_32_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE104_CANONICAL_MINBIAS_ISOLATION_SLIDES_25_32_V1",
                "generated_at": datetime.now(timezone.utc).isoformat(),
                "campaign_tag": args.campaign_tag,
                "inputs": {
                    "base_root": str(base_root),
                    "base_root_bytes": base_root.stat().st_size,
                    "base_root_sha256": sha256(base_root),
                    "phosub_root": str(phosub_root),
                    "phosub_root_bytes": phosub_root.stat().st_size,
                    "phosub_root_sha256": sha256(phosub_root),
                },
                "physics_contract": {
                    "sample": "embedded Photon+Jet 12+20 truth-matched signal",
                    "embedded_minbias_classifier": "required before physics filling",
                    "canonical_bdt": str(CANONICAL_BDT),
                    "canonical_wp80": CANONICAL_WP80,
                    "photon_pt_edges_gev": [15, 17, 19, 21, 23, 26, 30, 35],
                    "centrality_edges_percent": list(range(0, 85, 5)),
                    "isolation_cones": [0.3, 0.4],
                },
                "slides": packaged,
                "mutation_boundary": "local PNG package only; no Google Slides mutation",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"manifest={manifest}")
    for entry in packaged:
        print(entry["png"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
