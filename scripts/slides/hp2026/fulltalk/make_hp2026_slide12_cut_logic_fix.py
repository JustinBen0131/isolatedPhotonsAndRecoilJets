#!/usr/bin/env python3
"""Standalone generator for HP2026 deck Slide 12: "From cut logic to a BDT score".

Bifurcated from make_hp2026_bdt_yesno_sequence.py (slide09_cut_logic_to_bdt_score)
so this slide can be iterated independently while Codex works that file concurrently.
It imports only stable primitives:
  - leaf drawing helpers from the sequence module (gates, tree, score bar, text);
  - base helpers/colors from make_hp2026_fulltalk_candidates.
It owns the slide composition + header geometry here.

Fixes vs the current deck slide:
  - underline: the deck version repainted the header band and erased the shell
    divider without redrawing it. Here the shell is built with the final title up
    front, so the underline sits at the contract y=232 like slides 9 and 13.
  - card height: the three body bands are pulled into the 308..1248 region (the
    deck version ran to 1282 and crowded the footer).
  - title: stays 86 pt at (132, 76), matching the HP2026 main-header contract.

PNG-only candidate. Does not mutate Google Slides.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_fulltalk_candidates as full
import make_hp2026_bdt_yesno_sequence as seq


ROOT = full.ROOT
W, H = full.W, full.H
HDR = seq.HP2026_MAIN_HEADER
TITLE = seq.COMPRESSED_SLIDE_B_TITLE  # "From cut logic to a BDT score"

OUT_DIR = ROOT / "outputs/manual-20260611-hp2026-slide12-cut-logic-fix"
PNG_PATH = OUT_DIR / "hp2026_slide12_cut_logic_to_bdt_score_candidate.png"
SCRIPT_PATH = OUT_DIR / "hp2026_slide12_cut_logic_to_bdt_score_script.md"
MANIFEST_PATH = OUT_DIR / "manifest.json"


def build_slide() -> Image.Image:
    # Shell draws the top bars, the title (86 @ 132,76), the divider (y=232) and logo.
    # Building with the final title up front means the underline is never erased.
    img = seq.draw_slide_shell(None, title=TITLE)
    draw = ImageDraw.Draw(img, "RGBA")
    left_margin, right_margin = 132, W - 132

    # ----- Top band: fixed cuts as independent gates (contract top y=308) -----
    top = (left_margin, 308, right_margin, 500)
    full.shadow(img, top)
    draw.rounded_rectangle(top, radius=13, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((top[0], top[1], top[0] + 12, top[3]), radius=6, fill=(*full.PHOTON_DARK, 230))
    draw.text((top[0] + 42, top[1] + 26), "Fixed cuts test each handle separately", font=full.font(full.TIMES_BOLD, 43), fill=full.INK)
    gates_y = top[1] + 86
    gate_w, gate_gap, gate_x = 415, 28, top[0] + 46
    for i, (label, accent) in enumerate([
        ("compact core?", full.PHOTON_DARK),
        ("narrow shoulders?", full.SPHENIX_BLUE),
        ("not split?", full.TEAL),
    ]):
        seq.draw_simple_gate(draw, (gate_x + i * (gate_w + gate_gap), gates_y, gate_x + i * (gate_w + gate_gap) + gate_w, gates_y + 94), label, accent)
    takeaway = (top[0] + 1460, top[1] + 72, top[2] - 44, top[3] - 38)
    draw.rounded_rectangle(takeaway, radius=16, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    seq.draw_lines_centered(
        draw,
        (takeaway[0] + 28, takeaway[1] + 10, takeaway[2] - 28, takeaway[3] - 10),
        [
            ("Transparent, but rigid:", full.font(full.TIMES_BOLD, 31), full.BLUE),
            ("independent thresholds", full.font(full.TIMES_BOLD, 31), full.BLUE),
        ],
        line_gap=5,
    )

    # ----- Middle card: the decision tree (524..1008) -----
    mid = (left_margin, 524, right_margin, 1008)
    full.shadow(img, mid)
    draw.rounded_rectangle(mid, radius=13, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((mid[0], mid[1], mid[0] + 12, mid[3]), radius=6, fill=(*full.SPHENIX_BLUE, 230))
    seq.centered_text(draw, (mid[0] + 42, mid[1] + 24, mid[2] - 42, mid[1] + 82), "One learned decision tree", full.font(full.TIMES_BOLD, 48), full.INK)
    seq.draw_bridge_tree(draw, (mid[0] + 170, mid[1] + 94, mid[2] - 170, mid[3] - 58))
    tree_note = (mid[2] - 680, mid[1] + 28, mid[2] - 48, mid[1] + 116)
    draw.rounded_rectangle(tree_note, radius=14, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    seq.draw_lines_centered(
        draw,
        (tree_note[0] + 22, tree_note[1] + 8, tree_note[2] - 22, tree_note[3] - 8),
        [("Same inputs, learned conditional order.", full.font(full.TIMES_BOLD, 32), full.BLUE)],
    )

    # ----- Bottom band: boosting to one score (1032..1248, contract bottom) -----
    bot = (left_margin, 1032, right_margin, 1248)
    full.shadow(img, bot)
    draw.rounded_rectangle(bot, radius=13, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((bot[0], bot[1], bot[0] + 12, bot[3]), radius=6, fill=(*full.TEAL, 225))
    draw.text((bot[0] + 42, bot[1] + 24), "Boosting turns many trees into one score", font=full.font(full.TIMES_BOLD, 43), fill=full.INK)
    icon_y = bot[1] + 94
    icon_xs = [bot[0] + 450, bot[0] + 620, bot[0] + 790, bot[0] + 960]
    for idx, ix in enumerate(icon_xs):
        full.draw_mini_tree(draw, (ix, icon_y), scale=0.82, alpha=230)
        if idx < len(icon_xs) - 1:
            draw.text((ix + 84, icon_y + 58), "+", font=full.font(full.TIMES_BOLD, 44), fill=full.LIGHT_MUTED)
    full.draw_arrow(draw, (bot[0] + 1140, bot[1] + 142), (bot[0] + 1310, bot[1] + 142), fill=(126, 145, 164), width=7)
    score = (bot[0] + 1345, bot[1] + 94, bot[2] - 350, bot[1] + 178)
    seq.draw_score_bar(draw, score)
    score_note = (bot[2] - 318, bot[1] + 58, bot[2] - 44, bot[3] - 44)
    draw.rounded_rectangle(score_note, radius=14, fill=(238, 247, 247, 255), outline=(176, 211, 214, 255), width=2)
    seq.draw_lines_centered(
        draw,
        (score_note[0] + 20, score_note[1] + 8, score_note[2] - 20, score_note[3] - 8),
        [
            ("A stable ensemble", full.font(full.TIMES_BOLD, 29), full.BLUE),
            ("ranks each candidate", full.font(full.TIMES_BOLD, 29), full.BLUE),
            ("photon-ID score.", full.font(full.TIMES_BOLD, 29), full.BLUE),
        ],
        line_gap=3,
    )

    full.draw_hp2026_identity_footer(img)
    return img


def write_script() -> None:
    SCRIPT_PATH.write_text(
        """# HP2026 Slide 12 Script - From cut logic to a BDT score

So far I have three interpretable shower-shape handles, and the simplest thing to do is cut on each one separately: is the core compact, are the shoulders narrow, is the shower not split. That is transparent, but it is rigid, because each threshold is applied independently and the handles are actually correlated.

A single decision tree keeps the same familiar yes/no questions but lets the order be learned, so the next question can depend on the previous answer. That already exploits the correlations the fixed cuts ignore.

The BDT then boosts many shallow trees of this kind and combines them into one photon-ID score, low for background-like clusters and high for photon-like ones. That single score is the identification axis I carry into the purity measurement.
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
                "deck_slide": "HPslides_v1 Slide 12 ('From cut logic to a BDT score')",
                "bifurcated_from": "scripts/slides/hp2026/fulltalk/make_hp2026_bdt_yesno_sequence.py::slide09_cut_logic_to_bdt_score (Codex-owned; not edited)",
                "imports": ["make_hp2026_fulltalk_candidates (base helpers)", "make_hp2026_bdt_yesno_sequence (stable leaf primitives only: draw_slide_shell, draw_simple_gate, draw_bridge_tree, draw_score_bar, centered_text, draw_lines_centered)"],
                "output_png": str(PNG_PATH.relative_to(ROOT)),
                "companion_script": str(SCRIPT_PATH.relative_to(ROOT)),
                "header_contract": {"title_font_size": 86, "title_xy": [132, 76], "divider_y": 232},
                "fixes": [
                    "Title underline drawn at contract y=232 (deck version erased the shell divider via a header repaint and never redrew it).",
                    "Three body cards pulled into 308..1248 (deck version ran to 1282 and crowded the footer): top 308-500, mid 524-1008, bot 1032-1248.",
                    "Title kept at 86 @ (132,76), matching slides 9 and 13.",
                ],
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
