#!/usr/bin/env python3
"""Reusable layout-graph symmetry checks for generated slide PNGs.

This is the durable first layer of the "think like an artist" slide policy:
define usable bands, put objects on horizontal/vertical bisectors, compare
padding, and verify pixel-stable regions across build states.

It is intentionally deterministic. The future learning layer can add new
constraints from repeated Justin feedback, but slide generators should already
emit enough geometry for these checks before a candidate is shown.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

from PIL import Image, ImageChops


def _style_runs_to_bold_flags(text: str, runs: Any) -> list[bool] | None:
    if not isinstance(runs, list):
        return None
    flags: list[bool] = []
    for run in runs:
        if not isinstance(run, dict):
            return None
        content = str(run.get("text") or run.get("content") or "")
        flags.extend([bool(run.get("bold"))] * len(content))
    if len(flags) != len(text):
        return None
    return flags


def _colon_label_spans(text: str) -> list[tuple[int, int, int]]:
    spans: list[tuple[int, int, int]] = []
    offset = 0
    for line in text.splitlines(keepends=True):
        body = line.rstrip("\r\n")
        colon = body.find(":")
        if colon > 0 and colon <= 50:
            after = colon + 1
            while after < len(body) and body[after].isspace():
                after += 1
            if after < len(body):
                spans.append((offset, offset + colon, offset + after))
        offset += len(line)
    return spans


@dataclass(frozen=True)
class Box:
    """A rectangular layout node in pixel coordinates."""

    x0: float
    y0: float
    x1: float
    y1: float
    name: str = ""

    @property
    def cx(self) -> float:
        return (self.x0 + self.x1) / 2

    @property
    def cy(self) -> float:
        return (self.y0 + self.y1) / 2

    @property
    def w(self) -> float:
        return self.x1 - self.x0

    @property
    def h(self) -> float:
        return self.y1 - self.y0

    def inset(self, dx: float = 0, dy: float = 0, *, name: str | None = None) -> "Box":
        return Box(self.x0 + dx, self.y0 + dy, self.x1 - dx, self.y1 - dy, self.name if name is None else name)

    def as_int_tuple(self) -> tuple[int, int, int, int]:
        return (round(self.x0), round(self.y0), round(self.x1), round(self.y1))


@dataclass
class Check:
    name: str
    kind: str
    ok: bool
    got: Any = None
    want: Any = None
    delta: Any = None
    tolerance: float | None = None
    details: dict[str, Any] | None = None

    def to_json(self) -> dict[str, Any]:
        data = asdict(self)
        return {k: v for k, v in data.items() if v is not None}


class SymmetryAudit:
    """Collects layout-graph and pixel-continuity checks."""

    def __init__(self, *, name: str, concept: str | None = None) -> None:
        self.name = name
        self.concept = concept or (
            "layout graph: nodes=boxes/diagrams/text bands; "
            "edges=centering/equal padding/equal cells/pixel continuity"
        )
        self.checks: list[Check] = []

    def close(self, name: str, got: float, want: float, tol: float) -> None:
        delta = got - want
        self.checks.append(
            Check(
                name=name,
                kind="close",
                ok=abs(delta) <= tol,
                got=round(got, 3),
                want=round(want, 3),
                delta=round(delta, 3),
                tolerance=tol,
            )
        )

    def centered_on_x(self, child: Box, parent: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        self.close(name or f"{child.name} centered on {parent.name} X", child.cx, parent.cx, tol)

    def centered_on_y(self, child: Box, parent: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        self.close(name or f"{child.name} centered on {parent.name} Y", child.cy, parent.cy, tol)

    def centered_on(self, child: Box, parent: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        prefix = name or f"{child.name} centered on {parent.name}"
        self.centered_on_x(child, parent, tol, name=f"{prefix} X")
        self.centered_on_y(child, parent, tol, name=f"{prefix} Y")

    def same_y_center(self, a: Box, b: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        self.close(name or f"{a.name} and {b.name} share Y center", a.cy, b.cy, tol)

    def same_x_center(self, a: Box, b: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        self.close(name or f"{a.name} and {b.name} share X center", a.cx, b.cx, tol)

    def equal_padding_x(self, parent: Box, left_child: Box, right_child: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        left_pad = left_child.x0 - parent.x0
        right_pad = parent.x1 - right_child.x1
        self.close(name or f"{parent.name} left/right padding", left_pad, right_pad, tol)

    def equal_padding_y(self, parent: Box, top_child: Box, bottom_child: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        top_pad = top_child.y0 - parent.y0
        bottom_pad = parent.y1 - bottom_child.y1
        self.close(name or f"{parent.name} top/bottom padding", top_pad, bottom_pad, tol)

    def cell_center_x(self, parent: Box, index: int, count: int, got_x: float, tol: float = 1.0, *, name: str | None = None) -> None:
        want = parent.x0 + parent.w * (index + 0.5) / count
        self.close(name or f"{parent.name} cell {index + 1}/{count} X center", got_x, want, tol)

    def within(self, child: Box, parent: Box, pad: float = 0, *, name: str | None = None) -> None:
        ok = (
            child.x0 >= parent.x0 + pad
            and child.y0 >= parent.y0 + pad
            and child.x1 <= parent.x1 - pad
            and child.y1 <= parent.y1 - pad
        )
        self.checks.append(
            Check(
                name=name or f"{child.name} stays inside {parent.name}",
                kind="within",
                ok=ok,
                got=child.as_int_tuple(),
                want=parent.as_int_tuple(),
                tolerance=pad,
            )
        )

    def box_centered_y_in_band(self, child: Box, band: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        self.close(name or f"{child.name} centered in {band.name} Y band", child.cy, band.cy, tol)

    def box_centered_x_in_band(self, child: Box, band: Box, tol: float = 1.0, *, name: str | None = None) -> None:
        self.close(name or f"{child.name} centered in {band.name} X band", child.cx, band.cx, tol)

    def font_size_at_least(self, name: str, px: float, minimum_px: float, *, context: str | None = None) -> None:
        self.checks.append(
            Check(
                name=name,
                kind="font_size_minimum",
                ok=px >= minimum_px,
                got=round(px, 3),
                want=round(minimum_px, 3),
                delta=round(px - minimum_px, 3),
                details={"context": context} if context else None,
            )
        )

    def colon_label_style(self, node: dict[str, Any]) -> None:
        name = str(node.get("name") or "audience text")
        text = str(node.get("text") or "")
        spans = _colon_label_spans(text)
        if not spans or node.get("colon_style_exception"):
            return
        flags = _style_runs_to_bold_flags(text, node.get("text_runs") or node.get("runs"))
        if flags is None:
            self.checks.append(
                Check(
                    name=f"{name} colon-label style evidence",
                    kind="colon_label_style",
                    ok=False,
                    got="missing text_runs",
                    want="bold lead label plus colon, regular body after colon",
                    details={"text": text[:120]},
                )
            )
            return
        for start, colon, after in spans:
            lead_flags = [flags[idx] for idx in range(start, colon + 1) if not text[idx].isspace()]
            body_flags = [flags[idx] for idx in range(after, len(text)) if text[idx] not in "\r\n" and not text[idx].isspace()]
            ok = bool(lead_flags) and all(lead_flags) and bool(body_flags) and not any(body_flags)
            self.checks.append(
                Check(
                    name=f"{name} colon-label style",
                    kind="colon_label_style",
                    ok=ok,
                    got={
                        "lead_bold": all(lead_flags) if lead_flags else None,
                        "body_has_bold": any(body_flags) if body_flags else None,
                    },
                    want="lead label and colon bold; body regular",
                    details={"text": text[start : min(len(text), after + 80)]},
                )
            )

    def balanced_padding_y(self, parent: Box, child: Box, tol: float = 8.0, *, name: str | None = None) -> None:
        top = child.y0 - parent.y0
        bottom = parent.y1 - child.y1
        self.close(name or f"{child.name} balanced vertical padding inside {parent.name}", top, bottom, tol)

    def balanced_padding_x(self, parent: Box, child: Box, tol: float = 8.0, *, name: str | None = None) -> None:
        left = child.x0 - parent.x0
        right = parent.x1 - child.x1
        self.close(name or f"{child.name} balanced horizontal padding inside {parent.name}", left, right, tol)

    def image_size(self, path: Path, size: tuple[int, int]) -> None:
        with Image.open(path) as im:
            self.checks.append(
                Check(
                    name=f"{path.name} size",
                    kind="image_size",
                    ok=im.size == size,
                    got=list(im.size),
                    want=list(size),
                )
            )

    def pixel_identity(self, path_a: Path, path_b: Path, crop: tuple[int, int, int, int], *, name: str) -> None:
        with Image.open(path_a).convert("RGB") as a, Image.open(path_b).convert("RGB") as b:
            diff_bbox = ImageChops.difference(a.crop(crop), b.crop(crop)).getbbox()
        self.checks.append(
            Check(
                name=name,
                kind="pixel_identity",
                ok=diff_bbox is None,
                got=None if diff_bbox is None else list(diff_bbox),
                want=None,
                details={"crop": list(crop), "path_a": str(path_a), "path_b": str(path_b)},
            )
        )

    def extend(self, checks: Iterable[Check]) -> None:
        self.checks.extend(checks)

    @property
    def ok(self) -> bool:
        return all(check.ok for check in self.checks)

    def report(self) -> dict[str, Any]:
        return {
            "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
            "name": self.name,
            "ok": self.ok,
            "concept": self.concept,
            "checks": [check.to_json() for check in self.checks],
        }

    def write_json(self, path: Path) -> dict[str, Any]:
        path.parent.mkdir(parents=True, exist_ok=True)
        report = self.report()
        path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        return report


def require_audit_passed(report: dict[str, Any]) -> None:
    """Raise a compact error listing failed checks."""

    if report.get("ok"):
        return
    failed = [c for c in report.get("checks", []) if not c.get("ok")]
    lines = [f"{len(failed)} symmetry checks failed:"]
    for check in failed[:12]:
        lines.append(f"- {check.get('name')}: got {check.get('got')} want {check.get('want')}")
    if len(failed) > 12:
        lines.append(f"- ... {len(failed) - 12} more")
    raise RuntimeError("\n".join(lines))


def box_from_node(node: dict[str, Any]) -> Box:
    bbox = node.get("bbox")
    if not isinstance(bbox, (list, tuple)) or len(bbox) != 4:
        raise ValueError(f"layout node missing bbox: {node.get('name')}")
    return Box(float(bbox[0]), float(bbox[1]), float(bbox[2]), float(bbox[3]), str(node.get("name") or "node"))


def audit_layout_nodes(
    audit: SymmetryAudit,
    nodes: list[dict[str, Any]],
    *,
    slide_width_px: float = 2560,
    slide_height_px: float = 1440,
    min_audience_font_px: float = 35,
    min_plot_annotation_font_px: float = 28,
    min_title_font_px: float | None = None,
    title_left_max_px: float | None = None,
    title_top_max_px: float | None = None,
    diagram_padding_tolerance_px: float = 8,
    caption_center_tolerance_px: float = 5,
    text_center_tolerance_px: float = 10,
    repeated_node_tolerance_px: float = 6,
) -> None:
    """Apply generic slide-quality checks to generator-emitted layout nodes.

    Expected node conventions:
    - text nodes: kind=text, role=audience|title, font_px, bbox, parent=<box-name>
    - use intentional_alignment=top or vertical_alignment=top for deliberately
      top-aligned text boxes that should not be vertically centered
    - boxed non-text ink nodes: kind=diagram_ink|card_ink|figure_ink|image_ink|equation_ink,
      parent=<cell-node-name>, bbox. The parent box defines the usable band;
      the child bbox should sit on that band's bisectors unless a generator
      emits a deck-specific exception.
    - repeated panels/cards: symmetry_group=<name>, bbox
    - caption bands: name pattern matching `<prefix> process-label band`
    - process label text: name pattern matching `<prefix> process label`
    """

    by_name = {str(node.get("name")): node for node in nodes if node.get("name")}
    repeated_groups: dict[str, list[dict[str, Any]]] = {}
    for node in nodes:
        group = node.get("symmetry_group") or node.get("repeat_group")
        if group:
            repeated_groups.setdefault(str(group), []).append(node)

        role = str(node.get("role", "audience") or "audience")
        name = str(node.get("name") or "")
        is_title = role == "title" or bool(node.get("title_anchor"))
        if node.get("kind") == "text" and role in {"audience", "title", "plot_annotation", "process_marker"}:
            text_box = box_from_node(node)
            audit.checks.append(
                Check(
                    name=f"{text_box.name} stays inside slide canvas",
                    kind="canvas_containment",
                    ok=(
                        text_box.x0 >= 0
                        and text_box.y0 >= 0
                        and text_box.x1 <= slide_width_px
                        and text_box.y1 <= slide_height_px
                    ),
                    got=text_box.as_int_tuple(),
                    want=(0, 0, round(slide_width_px), round(slide_height_px)),
                    details={"text": str(node.get("text") or "")[:120]},
                )
            )
            font_px = float(node.get("font_px") or 0)
            if is_title and min_title_font_px:
                minimum_px = min_title_font_px
            elif role in {"plot_annotation", "process_marker"}:
                minimum_px = min_plot_annotation_font_px
            else:
                minimum_px = min_audience_font_px
            audit.font_size_at_least(
                str(node.get("name") or "audience text"),
                font_px,
                minimum_px,
                context=str(node.get("text") or "")[:120],
            )
            if role == "audience":
                audit.colon_label_style(node)
            if is_title:
                title_box = box_from_node(node)
                if title_left_max_px is not None:
                    audit.checks.append(
                        Check(
                            name=f"{title_box.name} left/title anchor",
                            kind="title_anchor",
                            ok=title_box.x0 <= title_left_max_px,
                            got=round(title_box.x0, 3),
                            want=f"<= {title_left_max_px}",
                        )
                    )
                if title_top_max_px is not None:
                    audit.checks.append(
                        Check(
                            name=f"{title_box.name} top/title anchor",
                            kind="title_anchor",
                            ok=title_box.y0 <= title_top_max_px,
                            got=round(title_box.y0, 3),
                            want=f"<= {title_top_max_px}",
                        )
                    )
            parent_name = node.get("parent") or node.get("container")
            parent_node = by_name.get(str(parent_name)) if parent_name else None
            if parent_node:
                parent = box_from_node(parent_node)
                child = box_from_node(node)
                audit.within(child, parent, 0, name=f"{child.name} inside {parent.name}")
                align_y = str(
                    node.get("intentional_alignment")
                    or node.get("vertical_alignment")
                    or node.get("align_y")
                    or ""
                ).lower()
                intentional_top = bool(node.get("intentional_top_aligned")) or align_y in {
                    "top",
                    "top_aligned",
                    "intentional_top",
                }
                if intentional_top:
                    audit.checks.append(
                        Check(
                            name=f"{child.name} intentionally top-aligned in {parent.name}",
                            kind="alignment_exception",
                            ok=True,
                            got=align_y or "intentional_top_aligned",
                            want="documented exception",
                        )
                    )
                else:
                    audit.box_centered_y_in_band(
                        child,
                        parent,
                        text_center_tolerance_px,
                        name=f"{child.name} vertically centered in {parent.name}",
                    )
                    audit.balanced_padding_y(
                        parent,
                        child,
                        text_center_tolerance_px,
                        name=f"{child.name} text top/bottom padding balanced in {parent.name}",
                    )

    for node in nodes:
        if node.get("kind") not in {"diagram_ink", "card_ink", "figure_ink", "image_ink", "equation_ink"}:
            continue
        parent_name = node.get("parent")
        parent_node = by_name.get(str(parent_name))
        if not parent_node:
            audit.checks.append(
                Check(
                    name=f"{node.get('name')} has parent node",
                    kind="node_link",
                    ok=False,
                    got=parent_name,
                    want="existing parent",
                )
            )
            continue
        parent = box_from_node(parent_node)
        child = box_from_node(node)
        audit.within(child, parent, 0, name=f"{child.name} inside {parent.name}")
        audit.box_centered_x_in_band(child, parent, diagram_padding_tolerance_px, name=f"{child.name} visually centered horizontally in {parent.name}")
        audit.box_centered_y_in_band(child, parent, diagram_padding_tolerance_px, name=f"{child.name} visually centered vertically in {parent.name}")
        audit.balanced_padding_x(parent, child, diagram_padding_tolerance_px, name=f"{child.name} left/right ink padding balanced in {parent.name}")
        audit.balanced_padding_y(parent, child, diagram_padding_tolerance_px, name=f"{child.name} top/bottom ink padding balanced in {parent.name}")

    for group_name, group_nodes in sorted(repeated_groups.items()):
        if len(group_nodes) < 2:
            continue
        boxes = [box_from_node(node) for node in group_nodes]
        reference = boxes[0]
        for box in boxes[1:]:
            audit.close(
                f"{group_name}: {box.name} width matches {reference.name}",
                box.w,
                reference.w,
                repeated_node_tolerance_px,
            )
            audit.close(
                f"{group_name}: {box.name} height matches {reference.name}",
                box.h,
                reference.h,
                repeated_node_tolerance_px,
            )
            audit.same_y_center(
                reference,
                box,
                repeated_node_tolerance_px,
                name=f"{group_name}: {box.name} shares horizontal bisector with {reference.name}",
            )
        for key in ("fill_color", "edge_color"):
            values = [node.get(key) for node in group_nodes if node.get(key)]
            if len(values) < 2:
                continue
            all_match = len({str(value).lower() for value in values}) == 1
            documented = all(
                bool(node.get("color_difference_intentional") or node.get("color_difference_reason"))
                for node in group_nodes
                if node.get(key)
            )
            details = None
            if documented and not all_match:
                details = {
                    "reason": [
                        node.get("color_difference_reason")
                        for node in group_nodes
                        if node.get("color_difference_reason")
                    ]
                }
            audit.checks.append(
                Check(
                    name=f"{group_name}: {key} consistency or documented semantic difference",
                    kind="symmetry_group_color",
                    ok=all_match or documented,
                    got=values,
                    want="matching colors, or documented semantic color difference",
                    details=details,
                )
            )

    for node in nodes:
        name = str(node.get("name") or "")
        if not name.endswith(" process label") or node.get("kind") != "text":
            continue
        band_name = name.removesuffix(" process label") + " process-label band"
        band_node = by_name.get(band_name)
        if not band_node:
            audit.checks.append(
                Check(
                    name=f"{name} has process-label band",
                    kind="node_link",
                    ok=False,
                    got=band_name,
                    want="existing band",
                )
            )
            continue
        label_box = box_from_node(node)
        band = box_from_node(band_node)
        audit.within(label_box, band, 0, name=f"{label_box.name} inside {band.name}")
        audit.box_centered_y_in_band(label_box, band, caption_center_tolerance_px, name=f"{label_box.name} centered in caption band")
