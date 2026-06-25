#!/usr/bin/env python3
"""Build the THE-58 low-calo event-quality cut JSON the trainer consumes.

Single source of truth = the THE-58 slide-5 manifest (the valley-anchored smooth
cut T(c) = a2*c^2 + a1*c + a0 and its per-5%-centrality-bin thresholds + sep
flags). This converts that cut into the envelope schema read by
``scripts/ml/training/train_auau_photon_bdt.py::load_low_calo_cut`` /
``apply_low_calo_event_quality_filter`` (required row fields: cent_lo, cent_hi,
threshold; reject an event when ``log10(E_CEMC+E_IHCal+E_OHCal+1) < threshold``
for its centrality slice ``cent_lo <= centrality < cent_hi``).

Contract preserved from the slide:
  * separable bins (0-55%): threshold = the quadratic value T(c) (these are cut).
  * merged bins (55-80%): NOT cut (band and artifact inseparable) -> a never-cut
    sentinel threshold (-1e9). log10(E_calo+1) >= 0 always, so nothing is rejected
    there, exactly matching the slide's "not cut" region.

Do NOT hand-edit the output JSON; regenerate it from the manifest if the cut
changes. This JSON is the THE-58 cut that SUPERSEDES the premature THE-32 v1 cut.
"""
from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
MANIFEST = REPO / "dataOutput/auauTightBDTValidation/THE58_lowCaloSmoothCut_20260614/the58_low_calo_smooth_cut_slide5_v3.manifest.json"
OUT = REPO / "dataOutput/auauTightBDTValidation/THE58_lowCaloSmoothCut_20260614/the58_low_calo_cut_v1.json"

NOT_CUT_SENTINEL = -1.0e9  # log10(E_calo+1) >= 0 always, so this never rejects.


def main() -> int:
    man = json.loads(MANIFEST.read_text())
    coef = [float(x) for x in man["coef"]]  # [a2, a1, a0]
    per_bin = man["per_bin"]

    envelope = []
    for b in per_bin:
        lo = float(b["lo"])
        hi = float(b["hi"])
        sep = bool(b["sep"])
        if sep:
            row = {
                "cent_lo": lo,
                "cent_hi": hi,
                "threshold": float(b["thr"]),
                "status": "ok",
                "separable": True,
                "expected_removed_frac": float(b["frac"]),
            }
        else:
            row = {
                "cent_lo": lo,
                "cent_hi": hi,
                "threshold": NOT_CUT_SENTINEL,
                "status": "not_cut",
                "separable": False,
                "expected_removed_frac": 0.0,
            }
        envelope.append(row)

    payload = {
        "schema": "THE58_LOW_CALO_CUT_V1",
        "supersedes": "THE32_LOW_CALO_CUT_V1",
        "cut_variable": "log10(E_CEMC + E_IHCal + E_OHCal + 1)",
        "centrality_source": "CentralityInfo::mbd_NS",
        "method": "valley_anchor_smooth_fit",
        "fit_degree": int(man.get("degree", 2)),
        "quadratic_coef_a2_a1_a0": coef,
        "fixed_artifact_log10E": float(man.get("artifact", 2.845)),
        "slice_width_percent": 5.0,
        "not_cut_sentinel": NOT_CUT_SENTINEL,
        "source_blind": True,
        "truth_blind": True,
        "bdt_score_blind": True,
        "boundary_convention": "cent_lo <= centrality < cent_hi; reject event when log10(E_calo+1) < threshold; equality retained; merged bins use a never-cut sentinel",
        "provenance": {
            "source_manifest": str(MANIFEST.relative_to(REPO)),
            "source_slide_gen": "scripts/slides/auau_bdt/low_calo_closure/make_the58_low_calo_smooth_cut_slide.py",
            "manifest_total": float(man.get("total", 0.0)),
            "manifest_removed": float(man.get("removed", 0.0)),
            "manifest_removed_frac": float(man.get("removed_frac", 0.0)),
        },
        "envelope": envelope,
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    n_cut = sum(1 for r in envelope if r["status"] == "ok")
    n_nocut = sum(1 for r in envelope if r["status"] == "not_cut")
    print(f"[OK] wrote {OUT}")
    print(f"     envelope rows: {len(envelope)} ({n_cut} cut, {n_nocut} not_cut); T(c) coef={coef}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
