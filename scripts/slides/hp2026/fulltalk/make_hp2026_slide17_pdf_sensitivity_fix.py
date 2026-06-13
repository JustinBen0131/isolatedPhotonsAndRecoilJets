#!/usr/bin/env python3
"""Standalone generator for HP2026 deck Slide 17: "Main result: PDF sensitivity".

Bifurcated from make_hp2026_final_result_no_fade_sequence.py (slide19_full_pdf_overview)
so the RHS can be reworked independently of Codex's file. Imports the stable helpers
(base_slide, place_full_result_plot, callout_row_card, constants) and owns the composition.

Change vs the deck slide:
  - drop the middle "Why this matters" box (vague filler: 'Beyond rate' / 'Program');
  - go from three boxes to TWO taller boxes that fill the RHS column honestly;
  - rename "Complete readout" -> "Bottom line" with concrete, non-redundant text;
  - bigger type for back-of-room readability;
  - truthful wording: the spectrum is *sensitive to* / *carries* proton-PDF information
    (the bands overlap), not claimed to discriminate/separate the PDF sets.

PNG-only candidate. Does not mutate Google Slides.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from PIL import Image

import make_hp2026_fulltalk_candidates as full
import make_hp2026_final_result_no_fade_sequence as seq


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260611-hp2026-slide17-pdf-sensitivity-fix"
PNG_PATH = OUT_DIR / "hp2026_slide17_pdf_sensitivity_candidate.png"
SCRIPT_PATH = OUT_DIR / "hp2026_slide17_pdf_sensitivity_script.md"
MANIFEST_PATH = OUT_DIR / "manifest.json"

X0, X1 = seq.RIGHT_CARD_X0, seq.RIGHT_CARD_X1


def build_slide() -> Image.Image:
    img = seq.base_slide("Main result: PDF sensitivity")
    seq.place_full_result_plot(img, mask_bottom=False)

    # Box A (blue) - what the bottom panel is, enlarged.
    seq.callout_row_card(
        img,
        (X0, 276, X1, 742),
        "Bottom panel: PDF sensitivity",
        [
            ("Method", "repeat JETPHOX with different proton PDFs"),
            ("Readout", "each PDF set traces a distinct JETPHOX/data trend"),
        ],
        seq.DATA_BLUE,
        fill=seq.SOFT_BLUE,
        title_size=49,
        label_size=42,
        body_size=45,
    )

    # Box B (yellow) - the takeaway, rewritten and enlarged to fill the rest.
    seq.callout_row_card(
        img,
        (X0, 784, X1, 1248),
        "Bottom line",
        [
            ("Sensitive", "the spectrum carries proton-PDF information, not just the overall normalization"),
            ("Baseline", "the validated p+p isolated-photon reference for γ+jet and the future Au+Au program"),
        ],
        seq.READOUT_GOLD,
        fill=seq.SOFT_YELLOW,
        title_size=49,
        label_size=42,
        body_size=45,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def write_script() -> None:
    SCRIPT_PATH.write_text(
        """# HP2026 Slide 17 Script - Main result: PDF sensitivity

This is the same measured cross section, now read all the way down to the bottom panel. There I repeat the JETPHOX calculation with several modern proton PDF sets, and each one traces a distinct trend relative to the data.

The point I want the audience to take away is that this is more than a rate measurement. Because the predictions move when I change the proton PDF, the isolated-photon spectrum carries information about the proton itself, not just the overall normalization. The bands still overlap, so I am careful to call this sensitivity rather than a clean separation of the PDF sets.

And stepping back, this is the bottom line of the talk: the first sPHENIX isolated prompt-photon cross section sets a validated p+p reference for the photon and gamma-plus-jet program, and for the future Au+Au measurements that build on it.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img = build_slide()
    img.convert("RGB").save(PNG_PATH, "PNG")
    full.write_hp2026_main_header_spec(PNG_PATH)
    write_script()
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "google_slides_mutation": False,
                "deck_slide": "HPslides_v1 Slide 17 ('Main result: PDF sensitivity')",
                "bifurcated_from": "scripts/slides/hp2026/fulltalk/make_hp2026_final_result_no_fade_sequence.py::slide19_full_pdf_overview (Codex-owned; not edited)",
                "imports": ["make_hp2026_fulltalk_candidates", "make_hp2026_final_result_no_fade_sequence (base_slide, place_full_result_plot, callout_row_card, constants)"],
                "changed": [
                    "Removed the middle 'Why this matters' box (vague filler: Beyond rate / Program).",
                    "Three boxes -> two taller boxes filling the RHS column (blue 276-742, yellow 784-1248).",
                    "Renamed 'Complete readout' -> 'Bottom line' with concrete rows: 'More than a rate' + 'Baseline'.",
                    "Body/title/label type enlarged (body 36/40 -> 45) for back-of-room readability.",
                    "Truthful wording: spectrum 'carries / is sensitive to' proton-PDF information; not claimed to discriminate/separate PDF sets (bands overlap).",
                ],
                "labels_note": "Result figure still reads 'sPHENIX Internal'; public HP version needs the released/Preliminary-labeled cross-section figure.",
                "output_png": str(PNG_PATH.relative_to(ROOT)),
                "companion_script": str(SCRIPT_PATH.relative_to(ROOT)),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG_PATH)
    print(SCRIPT_PATH)
    print(MANIFEST_PATH)


if __name__ == "__main__":
    main()
