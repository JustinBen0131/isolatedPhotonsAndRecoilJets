#!/usr/bin/env python3
"""Symmetry audit for the HP2026 prompt-photon three-state sequence.

This is a campaign-specific declaration of geometry constraints. The reusable
logic lives in `scripts/slides/common/slide_symmetry_audit.py`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

THIS_FILE = Path(__file__).resolve()
SCRIPTS_DIR = next((p for p in THIS_FILE.parents if p.name == "scripts"), THIS_FILE.parent)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_symmetry_audit import Box, Check, SymmetryAudit, audit_layout_nodes, require_audit_passed

import make_hp2026_opening_motivation_slide as hp


ROOT = hp.ROOT
OUT_DIR = ROOT / "outputs/manual-20260609-slide5_7_prompt_photons_expand_sequence"
PNG_NAMES = [
    "slide05_state01_production_expanded.png",
    "slide06_state02_isolation_expanded.png",
    "slide07_state03_current_full_context.png",
]


def expanded_layout_checks(audit: SymmetryAudit) -> None:
    card = Box(132, 330, 2428, 1290, "slide05 production card")
    visual_band = Box(card.x0, card.y0 + 220, card.x1, card.y0 + 684, "slide05 usable production visual band")
    prompt = Box(455, 572, 2105, 992, "slide05 red prompt-photon-production box")

    inner_top = prompt.y0 + 100
    inner_bottom = prompt.y1 - 58
    direct = Box(531, inner_top, 1445, inner_bottom, "slide05 direct-photon box")
    frag = Box(1539, inner_top, 2029, inner_bottom, "slide05 fragmentation-photon box")

    audit.centered_on(prompt, visual_band, 0.5, name="slide05 red prompt box centered in usable band")
    audit.same_y_center(direct, frag, 0.5, name="slide05 direct/fragment boxes share one horizontal bisector")
    audit.centered_on_y(direct, Box(direct.x0, inner_top, direct.x1, inner_bottom, "slide05 direct usable band"), 0.5)
    audit.centered_on_y(frag, Box(frag.x0, inner_top, frag.x1, inner_bottom, "slide05 fragmentation usable band"), 0.5)
    audit.equal_padding_x(prompt, direct, frag, 1.0, name="slide05 red box has balanced inner left/right padding")
    audit.within(prompt, card, 24, name="slide05 prompt box remains inside production card")
    audit.close("slide05 production card uses footer-safe lower band", card.y1, 1290, 0)

    audit.cell_center_x(direct, 0, 2, (direct.x0 * 3 + direct.x1) / 4, 0.5, name="slide05 Compton figure on left blue-cell bisector")
    audit.cell_center_x(direct, 1, 2, (direct.x0 + direct.x1 * 3) / 4, 0.5, name="slide05 annihilation figure on right blue-cell bisector")
    audit.centered_on_x(Box(frag.cx - 1, frag.cy - 1, frag.cx + 1, frag.cy + 1, "slide05 fragmentation figure center"), frag, 0.5)


def compact_layout_checks(audit: SymmetryAudit) -> None:
    card = Box(132, 330, 1448, 1290, "slide06/07 production card")
    prompt = Box(252, 520, 1332, 920, "slide06/07 red prompt-photon-production box")
    direct = Box(296, 600, 876, 850, "slide06/07 direct-photon box")
    frag = Box(972, 600, 1280, 850, "slide06/07 fragmentation-photon box")

    audit.centered_on_x(prompt, card, 2.0, name="slide06/07 red prompt box centered in production card X")
    audit.same_y_center(direct, frag, 0.5, name="slide06/07 direct/fragment boxes share one horizontal bisector")
    audit.within(prompt, card, 24, name="slide06/07 prompt box remains inside production card")
    audit.close("slide06/07 production card uses footer-safe lower band", card.y1, 1290, 0)
    audit.cell_center_x(direct, 0, 2, (direct.x0 * 3 + direct.x1) / 4, 0.5, name="slide06/07 Compton figure on left blue-cell bisector")
    audit.cell_center_x(direct, 1, 2, (direct.x0 + direct.x1 * 3) / 4, 0.5, name="slide06/07 annihilation figure on right blue-cell bisector")
    audit.centered_on_x(Box(frag.cx - 1, frag.cy - 1, frag.cx + 1, frag.cy + 1, "slide06/07 fragmentation figure center"), frag, 0.5)


def image_checks(audit: SymmetryAudit) -> None:
    paths = [OUT_DIR / name for name in PNG_NAMES]
    for path in paths:
        audit.image_size(path, (hp.W, hp.H))

    footer_crop = (0, 1326, hp.W, hp.H)
    audit.pixel_identity(paths[0], paths[1], footer_crop, name="footer pixel identity slide05 vs slide06")
    audit.pixel_identity(paths[0], paths[2], footer_crop, name="footer pixel identity slide05 vs slide07")
    audit.pixel_identity(
        paths[1],
        paths[2],
        (120, 320, 1460, 1302),
        name="slide06 and slide07 production-card pixel identity",
    )


def emitted_node_checks(audit: SymmetryAudit) -> None:
    path = OUT_DIR / "layout_nodes.json"
    if not path.exists():
        audit.checks.append(
            Check(
                name="layout_nodes.json exists",
                kind="layout_nodes",
                ok=False,
                got=str(path),
                want="existing generator-emitted layout node manifest",
            )
        )
        return
    data = json.loads(path.read_text(encoding="utf-8"))
    nodes = data.get("nodes", [])
    if not isinstance(nodes, list):
        audit.checks.append(
            Check(
                name="layout_nodes.json nodes list",
                kind="layout_nodes",
                ok=False,
                got=type(nodes).__name__,
                want="list",
            )
        )
        return
    audit_layout_nodes(
        audit,
        nodes,
        min_audience_font_px=float(data.get("minimum_audience_font_px") or 35),
        diagram_padding_tolerance_px=8,
        caption_center_tolerance_px=5,
    )
    subtitle_limits = {
        "slide05 production card subtitle": 2428 - 70,
        "slide06 isolation expanded subtitle": 2428 - 70,
        "slide07 color-neutral subtitle": 2428 - 70,
    }
    for node in nodes:
        name = str(node.get("name") or "")
        if name not in subtitle_limits:
            continue
        box = Box(*[float(v) for v in node.get("bbox", [0, 0, 0, 0])], name)
        limit = subtitle_limits[name]
        audit.checks.append(
            Check(
                name=f"{name} right padding",
                kind="subtitle_containment",
                ok=box.x1 <= limit,
                got=round(box.x1, 3),
                want=f"<= {limit}",
                tolerance=0,
            )
        )
    focus_font_minima = {
        "slide05 slide subtitle": 56,
        "slide06 slide subtitle": 56,
        "slide07 slide subtitle": 56,
        "slide05 production card title": 76,
        "slide06 production card title": 72,
        "slide07 production card title": 72,
        "slide06 production prompt summary": 47,
        "slide07 production prompt summary": 47,
        "slide05 production background label": 40,
        "slide05 production decay line": 44,
        "slide06 production background label": 40,
        "slide06 production decay line": 42,
        "slide07 production background label": 40,
        "slide07 production decay line": 42,
        "slide06 isolated cone label": 42,
        "slide06 non-isolated cone label": 42,
        "slide07 isolated cone label": 35,
        "slide07 non-isolated cone label": 35,
        "slide06 isolation expanded equation": 64,
        "slide06 isolation expanded title": 72,
        "slide06 isolation expanded cue": 42,
        "slide06 isolation expanded takeaway line 1": 40,
        "slide06 isolation expanded takeaway line 2": 40,
        "slide07 isolation compact title": 54,
        "slide07 isolation compact equation": 38,
        "slide07 isolation compact cue": 36,
        "slide07 isolation compact takeaway line 1": 38,
        "slide07 isolation compact takeaway line 2": 38,
        "slide07 color-neutral title": 70,
        "slide07 color-neutral takeaway line 1": 42,
        "slide07 color-neutral takeaway line 2": 42,
    }
    for node in nodes:
        name = str(node.get("name") or "")
        if name not in focus_font_minima:
            continue
        audit.font_size_at_least(
            f"{name} focus readability",
            float(node.get("font_px") or 0),
            focus_font_minima[name],
            context="Prompt-photon three-state audience-facing text readability",
        )
    centered_in_boxes = {
        "slide06 isolation expanded equation": Box(1548, 524, 2370, 674, "slide06 yellow isolation equation box"),
    }
    for node in nodes:
        name = str(node.get("name") or "")
        if name not in centered_in_boxes:
            continue
        box = Box(*[float(v) for v in node.get("bbox", [0, 0, 0, 0])], name)
        parent = centered_in_boxes[name]
        audit.centered_on_x(box, parent, 3.0, name=f"{name} centered on yellow box vertical bisector")
        audit.close(
            f"{name} optically centered on yellow box horizontal bisector",
            box.cy,
            parent.cy + 8,
            3.0,
        )
    summary_containment = {
        "slide05 production prompt summary": (132 + 70, 2428 - 70),
        "slide05 production background label": (132 + 70, 2428 - 70),
        "slide05 production decay line": (132 + 70, 2428 - 70),
        "slide06 production prompt summary": (132 + 70, 1448 - 70),
        "slide06 production background label": (132 + 70, 1448 - 70),
        "slide06 production decay line": (132 + 70, 1448 - 70),
        "slide07 production prompt summary": (132 + 70, 1448 - 70),
        "slide07 production background label": (132 + 70, 1448 - 70),
        "slide07 production decay line": (132 + 70, 1448 - 70),
        "slide07 color-neutral takeaway line 1": (1496 + 52, 2428 - 34),
        "slide07 color-neutral takeaway line 2": (1496 + 52, 2428 - 34),
    }
    for node in nodes:
        name = str(node.get("name") or "")
        if name not in summary_containment:
            continue
        box = Box(*[float(v) for v in node.get("bbox", [0, 0, 0, 0])], name)
        min_x, max_x = summary_containment[name]
        audit.checks.append(
            Check(
                name=f"{name} horizontal containment",
                kind="text_containment",
                ok=box.x0 >= min_x and box.x1 <= max_x,
                got=[round(box.x0, 3), round(box.x1, 3)],
                want=[f">= {min_x}", f"<= {max_x}"],
                tolerance=0,
            )
        )


def build_audit() -> SymmetryAudit:
    audit = SymmetryAudit(
        name="HP2026 prompt photons slide 5-7 sequence symmetry",
        concept=(
            "graph-style slide QA: production cards, prompt boxes, inner physics boxes, "
            "Feynman cells, footer strips, and build-state continuity are nodes connected "
            "by centering, padding, containment, and pixel-identity constraints"
        ),
    )
    expanded_layout_checks(audit)
    compact_layout_checks(audit)
    image_checks(audit)
    emitted_node_checks(audit)
    return audit


def main() -> int:
    audit = build_audit()
    report = audit.write_json(OUT_DIR / "symmetry_audit.json")
    print(json.dumps(report, indent=2))
    require_audit_passed(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
