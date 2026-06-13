#!/usr/bin/env python3
"""Default post-render audit for serious generated slide PNG candidates.

This wraps the reusable layout-graph checks in `slide_symmetry_audit.py` with a
small CLI that slide generators can run after rendering. It intentionally keeps
deck-specific chrome, such as HP2026 headers, opt-in through separate profiles
or `codex_os_guard.py` flags.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from slide_defaults import SLIDE_DPI, SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX
from slide_symmetry_audit import Check, SymmetryAudit, audit_layout_nodes, require_audit_passed


INTERNAL_TEXT_RE = re.compile(
    r"^\s*(source|provenance|command|generated|implementation|todo|speaker note)\s*:|"
    r"\b(codex|claude)\b",
    re.IGNORECASE,
)
MARKDOWN_BULLET_RE = re.compile(r"(?m)^\s*(>|--?\s+|\*\s+)")
ASCII_TIMES_IN_FORMULA_RE = re.compile(
    r"\b\d+(?:\.\d+)?\s+x\s+\d+(?:\.\d+)?\b|"
    r"\b(median|mad|sigma|threshold|formula|t_bin|q_0\.1)\b.*\s+x\s+",
    re.IGNORECASE,
)


def pt_to_px(points: float, dpi: int = SLIDE_DPI) -> float:
    return points * dpi / 72.0


def load_layout_nodes(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, dict):
        nodes = payload.get("nodes", [])
        if not isinstance(nodes, list):
            raise ValueError(f"{path} has non-list `nodes`")
        return nodes, payload
    if isinstance(payload, list):
        return payload, {"nodes": payload}
    raise ValueError(f"{path} must be a JSON object with `nodes` or a node list")


def audit_internal_canvas_text(audit: SymmetryAudit, nodes: list[dict[str, Any]], *, allow_slide_numbers: bool) -> None:
    for node in nodes:
        if node.get("kind") != "text" or node.get("role", "audience") != "audience":
            continue
        name = str(node.get("name") or "audience text")
        text = str(node.get("text") or "")
        if not text:
            continue
        has_internal_text = bool(INTERNAL_TEXT_RE.search(text))
        audit.checks.append(
            Check(
                name=f"{name} avoids internal/provenance canvas text",
                kind="audience_text_contract",
                ok=not has_internal_text,
                got=text[:120] if has_internal_text else "audience-facing",
                want="no internal notes, provenance clutter, or agent names",
            )
        )
        is_slide_number = bool(re.fullmatch(r"\s*\d{1,3}\s*", text))
        audit.checks.append(
            Check(
                name=f"{name} avoids baked slide numbers",
                kind="slide_number_contract",
                ok=allow_slide_numbers or not is_slide_number,
                got=text.strip() if is_slide_number else "not a slide number",
                want="no baked slide number unless explicitly allowed",
            )
        )
        has_markdown_bullet = bool(MARKDOWN_BULLET_RE.search(text))
        audit.checks.append(
            Check(
                name=f"{name} avoids markdown bullet residue",
                kind="markdown_bullet_residue",
                ok=not has_markdown_bullet,
                got=text[:120] if has_markdown_bullet else "polished audience text",
                want="drawn bullets or native bullet glyphs, not markdown quote/list markers",
            )
        )
        has_ascii_times = bool(ASCII_TIMES_IN_FORMULA_RE.search(text))
        audit.checks.append(
            Check(
                name=f"{name} uses mathematical multiplication symbol when formula-like",
                kind="formula_multiplication_symbol",
                ok=not has_ascii_times,
                got=text[:120] if has_ascii_times else "no ASCII x-as-times pattern",
                want="use × or a styled math renderer, not bare x, for multiplication",
            )
        )


def audit_title_bottom_clearance(audit: SymmetryAudit, nodes: list[dict[str, Any]], *, clearance_px: float) -> None:
    title_bottoms: list[float] = []
    for node in nodes:
        if node.get("kind") != "text":
            continue
        role = str(node.get("role") or "")
        if role != "title" and not node.get("title_anchor"):
            continue
        bbox = node.get("bbox")
        if isinstance(bbox, list) and len(bbox) == 4:
            title_bottoms.append(float(bbox[3]))
    if not title_bottoms:
        return

    title_clear_y = max(title_bottoms) + clearance_px
    for node in nodes:
        if node.get("kind") != "text" or node.get("role", "audience") != "audience":
            continue
        if node.get("title_band_exception"):
            continue
        bbox = node.get("bbox")
        if not isinstance(bbox, list) or len(bbox) != 4:
            continue
        y0 = float(bbox[1])
        audit.checks.append(
            Check(
                name=f"{node.get('name') or 'audience text'} clears title bottom band",
                kind="title_bottom_clearance",
                ok=y0 >= title_clear_y,
                got=round(y0, 3),
                want=f">= {round(title_clear_y, 3)}",
                details={
                    "title_bottom_y": round(max(title_bottoms), 3),
                    "clearance_px": round(clearance_px, 3),
                    "text": str(node.get("text") or "")[:120],
                },
            )
        )


def audit_title_followup_gap(audit: SymmetryAudit, nodes: list[dict[str, Any]], *, max_gap_px: float) -> None:
    title_bottoms: list[float] = []
    for node in nodes:
        if node.get("kind") != "text":
            continue
        role = str(node.get("role") or "")
        if role != "title" and not node.get("title_anchor"):
            continue
        bbox = node.get("bbox")
        if isinstance(bbox, list) and len(bbox) == 4:
            title_bottoms.append(float(bbox[3]))
    if not title_bottoms:
        return

    title_bottom = max(title_bottoms)
    followup_candidates: list[tuple[float, dict[str, Any]]] = []
    for node in nodes:
        if node.get("kind") != "text" or node.get("role", "audience") != "audience":
            continue
        if node.get("title_band_exception"):
            continue
        bbox = node.get("bbox")
        if not isinstance(bbox, list) or len(bbox) != 4:
            continue
        y0 = float(bbox[1])
        if y0 >= title_bottom:
            followup_candidates.append((y0, node))
    if not followup_candidates:
        return

    first_y, first_node = min(followup_candidates, key=lambda item: item[0])
    gap = first_y - title_bottom
    audit.checks.append(
        Check(
            name="first audience text follows title without excessive gap",
            kind="title_followup_gap",
            ok=gap <= max_gap_px,
            got=round(gap, 3),
            want=f"<= {round(max_gap_px, 3)}",
            details={
                "title_bottom_y": round(title_bottom, 3),
                "first_text_y": round(first_y, 3),
                "first_text_name": str(first_node.get("name") or "audience text"),
                "text": str(first_node.get("text") or "")[:120],
            },
        )
    )


def audit_title_axis_alignment(
    audit: SymmetryAudit,
    nodes: list[dict[str, Any]],
    payload: dict[str, Any],
    *,
    tolerance_px: float,
) -> None:
    title_axis_raw = payload.get("title_axis_x")
    if title_axis_raw is None:
        title_lefts: list[float] = []
        for node in nodes:
            if node.get("kind") != "text":
                continue
            role = str(node.get("role") or "")
            if role != "title" and not node.get("title_anchor"):
                continue
            bbox = node.get("bbox")
            if isinstance(bbox, list) and len(bbox) == 4:
                title_lefts.append(float(bbox[0]))
        if not title_lefts:
            return
        title_axis = min(title_lefts)
    else:
        title_axis = float(title_axis_raw)

    tolerance = float(payload.get("title_axis_tolerance_px") or tolerance_px)
    for node in nodes:
        align = str(node.get("title_axis_align") or node.get("align_to_title_axis") or "").lower()
        if align not in {"left", "x0", "start"} and not node.get("title_axis_left"):
            continue
        bbox = node.get("bbox")
        if not isinstance(bbox, list) or len(bbox) != 4:
            continue
        got = float(bbox[0])
        audit.checks.append(
            Check(
                name=f"{node.get('name') or 'node'} aligns to title left axis",
                kind="title_axis_alignment",
                ok=abs(got - title_axis) <= tolerance,
                got=round(got, 3),
                want=f"{round(title_axis, 3)} +/- {round(tolerance, 3)}",
                delta=round(got - title_axis, 3),
                tolerance=tolerance,
            )
        )


def audit_text_block_vertical_fill(audit: SymmetryAudit, nodes: list[dict[str, Any]]) -> None:
    by_name = {str(node.get("name")): node for node in nodes if node.get("name")}
    for node in nodes:
        min_raw = node.get("vertical_fill_min_ratio")
        max_raw = node.get("vertical_fill_max_ratio")
        if min_raw is None and max_raw is None:
            continue
        parent_name = node.get("parent") or node.get("container")
        parent = by_name.get(str(parent_name)) if parent_name else None
        bbox = node.get("bbox")
        parent_bbox = parent.get("bbox") if isinstance(parent, dict) else None
        if not isinstance(bbox, list) or len(bbox) != 4 or not isinstance(parent_bbox, list) or len(parent_bbox) != 4:
            audit.checks.append(
                Check(
                    name=f"{node.get('name') or 'text block'} vertical fill evidence",
                    kind="text_block_vertical_fill",
                    ok=False,
                    got="missing node or parent bbox",
                    want="text block with parent bbox",
                )
            )
            continue
        child_h = float(bbox[3]) - float(bbox[1])
        parent_h = float(parent_bbox[3]) - float(parent_bbox[1])
        ratio = child_h / parent_h if parent_h > 0 else 0.0
        min_ratio = float(min_raw) if min_raw is not None else 0.0
        max_ratio = float(max_raw) if max_raw is not None else 1.0
        audit.checks.append(
            Check(
                name=f"{node.get('name') or 'text block'} uses intended vertical card space",
                kind="text_block_vertical_fill",
                ok=min_ratio <= ratio <= max_ratio,
                got=round(ratio, 3),
                want=f"{round(min_ratio, 3)} <= ratio <= {round(max_ratio, 3)}",
                details={
                    "child_height_px": round(child_h, 3),
                    "parent_height_px": round(parent_h, 3),
                    "parent": str(parent_name),
                },
            )
        )


def audit_acronym_explanation(audit: SymmetryAudit, nodes: list[dict[str, Any]]) -> None:
    for node in nodes:
        acronym = node.get("requires_acronym_expansion")
        if not acronym:
            continue
        text = str(node.get("text") or "")
        normalized = re.sub(r"\s+", " ", text).strip().lower()
        expansion_terms = node.get("acronym_expansion_terms") or []
        role_terms = node.get("functional_role_terms") or []
        expansion_ok = all(str(term).lower() in normalized for term in expansion_terms)
        role_ok = any(str(term).lower() in normalized for term in role_terms) if role_terms else True
        audit.checks.append(
            Check(
                name=f"{node.get('name') or acronym} explains acronym expansion and role",
                kind="acronym_expansion_and_role",
                ok=expansion_ok and role_ok,
                got=normalized[:160],
                want={
                    "acronym": str(acronym),
                    "expansion_terms": [str(term) for term in expansion_terms],
                    "one_functional_role_term": [str(term) for term in role_terms],
                },
            )
        )


def audit_outer_frame_top_gap(audit: SymmetryAudit, nodes: list[dict[str, Any]]) -> None:
    by_name = {str(node.get("name")): node for node in nodes if node.get("name")}
    for node in nodes:
        max_gap_raw = node.get("max_outer_frame_top_gap_px")
        outer_frame_name = node.get("outer_frame") or node.get("frame") or node.get("visual_frame")
        if max_gap_raw is None or not outer_frame_name:
            continue
        bbox = node.get("bbox")
        frame = by_name.get(str(outer_frame_name))
        frame_bbox = frame.get("bbox") if isinstance(frame, dict) else None
        if not isinstance(bbox, list) or len(bbox) != 4 or not isinstance(frame_bbox, list) or len(frame_bbox) != 4:
            audit.checks.append(
                Check(
                    name=f"{node.get('name') or 'framed content'} top gap evidence",
                    kind="outer_frame_top_gap",
                    ok=False,
                    got="missing content or frame bbox",
                    want="content bbox plus named outer_frame bbox",
                )
            )
            continue
        gap = float(bbox[1]) - float(frame_bbox[1])
        max_gap = float(max_gap_raw)
        audit.checks.append(
            Check(
                name=f"{node.get('name') or 'framed content'} starts near top of {outer_frame_name}",
                kind="outer_frame_top_gap",
                ok=0 <= gap <= max_gap,
                got=round(gap, 3),
                want=f"0 <= gap <= {round(max_gap, 3)}",
                details={
                    "frame_top_y": round(float(frame_bbox[1]), 3),
                    "content_top_y": round(float(bbox[1]), 3),
                },
            )
        )


def audit_vertical_margin_balance(audit: SymmetryAudit, nodes: list[dict[str, Any]], payload: dict[str, Any]) -> None:
    spec = payload.get("vertical_margin_balance")
    if not isinstance(spec, dict):
        return
    top_name = str(spec.get("top_node") or "")
    bottom_name = str(spec.get("bottom_node") or "")
    if not top_name or not bottom_name:
        audit.checks.append(
            Check(
                name="vertical margin balance evidence",
                kind="vertical_margin_balance",
                ok=False,
                got=spec,
                want="top_node and bottom_node names",
            )
        )
        return
    by_name = {str(node.get("name")): node for node in nodes if node.get("name")}
    top_node = by_name.get(top_name)
    bottom_node = by_name.get(bottom_name)
    top_bbox = top_node.get("bbox") if isinstance(top_node, dict) else None
    bottom_bbox = bottom_node.get("bbox") if isinstance(bottom_node, dict) else None
    slide_size = payload.get("slide_size") or [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX]
    if (
        not isinstance(top_bbox, list)
        or len(top_bbox) != 4
        or not isinstance(bottom_bbox, list)
        or len(bottom_bbox) != 4
        or not isinstance(slide_size, list)
        or len(slide_size) != 2
    ):
        audit.checks.append(
            Check(
                name="vertical margin balance evidence",
                kind="vertical_margin_balance",
                ok=False,
                got="missing top/bottom bbox or slide_size",
                want="named top and bottom node bboxes",
            )
        )
        return
    top_gap = float(top_bbox[1])
    bottom_gap = float(slide_size[1]) - float(bottom_bbox[3])
    tolerance = float(spec.get("tolerance_px") or 8.0)
    audit.checks.append(
        Check(
            name=f"{top_name} top margin balances {bottom_name} bottom margin",
            kind="vertical_margin_balance",
            ok=abs(top_gap - bottom_gap) <= tolerance,
            got={"top_gap_px": round(top_gap, 3), "bottom_gap_px": round(bottom_gap, 3)},
            want=f"abs(top-bottom) <= {round(tolerance, 3)} px",
            delta=round(top_gap - bottom_gap, 3),
            tolerance=tolerance,
        )
    )
    target_raw = spec.get("target_gap_px")
    if target_raw is not None:
        target = float(target_raw)
        top_ok = abs(top_gap - target) <= tolerance
        bottom_ok = abs(bottom_gap - target) <= tolerance
        audit.checks.append(
            Check(
                name=f"{top_name}/{bottom_name} vertical margins match target gap",
                kind="vertical_margin_target",
                ok=top_ok and bottom_ok,
                got={"top_gap_px": round(top_gap, 3), "bottom_gap_px": round(bottom_gap, 3)},
                want=f"{round(target, 3)} +/- {round(tolerance, 3)} px",
                details={"target_gap_px": round(target, 3)},
            )
        )


def build_audit(args: argparse.Namespace) -> tuple[SymmetryAudit, Path]:
    png = Path(args.png).expanduser().resolve()
    report_path = Path(args.output).expanduser().resolve() if args.output else png.with_suffix(".slide_audit.json")
    audit = SymmetryAudit(
        name=f"post-render slide audit: {png.name}",
        concept=(
            "default generated-slide audit: PNG size, audience text readability, "
            "layout-node containment, boxed-content bisectors, equal padding, "
            "colon-label styling, title anchoring, and repeated-panel symmetry"
        ),
    )
    audit.image_size(png, (args.expect_width, args.expect_height))

    if args.layout_nodes:
        layout_path = Path(args.layout_nodes).expanduser().resolve()
        try:
            nodes, payload = load_layout_nodes(layout_path)
            minimum_font_px = float(
                payload.get("minimum_audience_font_px")
                or pt_to_px(args.min_body_font_pt, args.dpi)
            )
            minimum_title_font_px = float(
                payload.get("minimum_title_font_px")
                or pt_to_px(args.min_title_font_pt, args.dpi)
            )
            minimum_plot_annotation_font_px = float(
                payload.get("minimum_plot_annotation_font_px")
                or pt_to_px(args.min_plot_annotation_font_pt, args.dpi)
            )
            audit_layout_nodes(
                audit,
                nodes,
                min_audience_font_px=minimum_font_px,
                min_plot_annotation_font_px=minimum_plot_annotation_font_px,
                min_title_font_px=minimum_title_font_px,
                title_left_max_px=args.title_left_max_px,
                title_top_max_px=args.title_top_max_px,
                diagram_padding_tolerance_px=args.diagram_tolerance_px,
                caption_center_tolerance_px=args.caption_tolerance_px,
                text_center_tolerance_px=args.text_center_tolerance_px,
                repeated_node_tolerance_px=args.repeated_node_tolerance_px,
            )
            audit_internal_canvas_text(audit, nodes, allow_slide_numbers=args.allow_slide_numbers)
            audit_title_bottom_clearance(audit, nodes, clearance_px=args.title_bottom_clearance_px)
            audit_title_followup_gap(audit, nodes, max_gap_px=args.title_followup_gap_max_px)
            audit_title_axis_alignment(audit, nodes, payload, tolerance_px=args.title_axis_tolerance_px)
            audit_text_block_vertical_fill(audit, nodes)
            audit_acronym_explanation(audit, nodes)
            audit_outer_frame_top_gap(audit, nodes)
            audit_vertical_margin_balance(audit, nodes, payload)
            audit.checks.append(
                Check(
                    name="layout node manifest loaded",
                    kind="layout_nodes",
                    ok=True,
                    got=str(layout_path),
                    want="generator-emitted layout nodes",
                    details={"node_count": len(nodes)},
                )
            )
        except Exception as exc:
            audit.checks.append(
                Check(
                    name="layout node manifest loads",
                    kind="layout_nodes",
                    ok=False,
                    got=str(exc),
                    want="valid layout node JSON",
                )
            )
    else:
        audit.checks.append(
            Check(
                name="layout node manifest not supplied",
                kind="layout_nodes_optional",
                ok=True,
                got="none",
                want="optional for legacy generators; preferred for serious slides",
            )
        )

    audit.checks.append(
        Check(
            name="deck profile remains generic unless explicitly selected",
            kind="deck_profile",
            ok=args.deck_profile in {"generic", "hp2026-main-talk"},
            got=args.deck_profile,
            want="generic or known opt-in deck profile",
        )
    )
    return audit, report_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png", required=True, help="Rendered 16:9 slide PNG")
    parser.add_argument("--layout-nodes", help="Optional generator-emitted layout_nodes.json")
    parser.add_argument("--output", help="Audit JSON path; defaults beside PNG")
    parser.add_argument("--deck-profile", default="generic", choices=["generic", "hp2026-main-talk"])
    parser.add_argument("--expect-width", type=int, default=SLIDE_WIDTH_PX)
    parser.add_argument("--expect-height", type=int, default=SLIDE_HEIGHT_PX)
    parser.add_argument("--dpi", type=int, default=SLIDE_DPI)
    # Font thresholds are policy pt-equivalents converted at render DPI, not
    # slide-specific values tuned to make a particular artifact pass.
    parser.add_argument("--min-body-font-pt", type=float, default=13.0)
    parser.add_argument("--min-plot-annotation-font-pt", type=float, default=10.0)
    parser.add_argument("--min-title-font-pt", type=float, default=24.0)
    parser.add_argument("--title-left-max-px", type=float, default=180.0)
    parser.add_argument("--title-top-max-px", type=float, default=130.0)
    parser.add_argument("--title-bottom-clearance-px", type=float, default=28.0)
    parser.add_argument("--title-followup-gap-max-px", type=float, default=220.0)
    parser.add_argument("--title-axis-tolerance-px", type=float, default=3.0)
    parser.add_argument("--diagram-tolerance-px", type=float, default=8.0)
    parser.add_argument("--caption-tolerance-px", type=float, default=5.0)
    parser.add_argument("--text-center-tolerance-px", type=float, default=10.0)
    parser.add_argument("--repeated-node-tolerance-px", type=float, default=6.0)
    parser.add_argument("--allow-slide-numbers", action="store_true")
    parser.add_argument("--require-pass", action="store_true")
    parser.add_argument("--json", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    audit, report_path = build_audit(args)
    report = audit.write_json(report_path)
    if args.require_pass:
        require_audit_passed(report)
    if args.json:
        print(json.dumps({"ok": report["ok"], "report": str(report_path)}, sort_keys=True))
    else:
        status = "OK" if report["ok"] else "FAIL"
        print(f"{status}: slide audit wrote {report_path}")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
