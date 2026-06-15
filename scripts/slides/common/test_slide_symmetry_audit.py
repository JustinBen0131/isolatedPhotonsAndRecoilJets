#!/usr/bin/env python3
"""Focused tests for the reusable slide post-render audit layer."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from PIL import Image

from slide_symmetry_audit import Box, SymmetryAudit, audit_layout_nodes


THIS_DIR = Path(__file__).resolve().parent
POST_RENDER = THIS_DIR / "post_render_slide_audit.py"


class SlideSymmetryAuditTests(unittest.TestCase):
    def test_within_records_pass_and_fail(self) -> None:
        audit = SymmetryAudit(name="within-test")
        parent = Box(0, 0, 100, 100, "parent")
        audit.within(Box(10, 10, 90, 90, "inside"), parent, 0)
        audit.within(Box(-1, 10, 90, 90, "outside"), parent, 0)
        self.assertEqual([check.kind for check in audit.checks], ["within", "within"])
        self.assertTrue(audit.checks[0].ok)
        self.assertFalse(audit.checks[1].ok)

    def test_small_audience_text_fails(self) -> None:
        audit = SymmetryAudit(name="font-test")
        nodes = [
            {"name": "card", "kind": "panel", "bbox": [0, 0, 200, 120]},
            {
                "name": "body text",
                "kind": "text",
                "role": "audience",
                "parent": "card",
                "bbox": [20, 40, 180, 80],
                "font_px": 30,
                "text": "Small body text",
            },
        ]
        audit_layout_nodes(audit, nodes, min_audience_font_px=36)
        font_checks = [check for check in audit.checks if check.kind == "font_size_minimum"]
        self.assertEqual(len(font_checks), 1)
        self.assertFalse(font_checks[0].ok)

    def test_intentional_top_alignment_is_documented_exception(self) -> None:
        audit = SymmetryAudit(name="top-align-test")
        nodes = [
            {"name": "card", "kind": "panel", "bbox": [0, 0, 200, 120]},
            {
                "name": "top aligned text",
                "kind": "text",
                "role": "audience",
                "parent": "card",
                "bbox": [20, 12, 180, 42],
                "font_px": 40,
                "text": "Intentional lead line",
                "intentional_alignment": "top",
            },
        ]
        audit_layout_nodes(audit, nodes, min_audience_font_px=36)
        self.assertTrue(any(check.kind == "alignment_exception" for check in audit.checks))
        self.assertTrue(audit.ok)

    def test_clean_wrapper_manifest_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            nodes = {
                "minimum_audience_font_px": 36,
                "minimum_title_font_px": 67,
                "nodes": [
                    {
                        "name": "slide title",
                        "kind": "text",
                        "role": "title",
                        "bbox": [100, 70, 1200, 145],
                        "font_px": 72,
                        "text": "Reference-style left title",
                    },
                    {
                        "name": "slide subtitle",
                        "kind": "text",
                        "role": "audience",
                        "bbox": [102, 178, 1700, 222],
                        "font_px": 40,
                        "text": "Readable audience subtitle",
                    },
                    {"name": "left card", "kind": "panel", "bbox": [100, 300, 900, 900], "symmetry_group": "cards"},
                    {"name": "right card", "kind": "panel", "bbox": [1000, 300, 1800, 900], "symmetry_group": "cards"},
                    {
                        "name": "left card body",
                        "kind": "text",
                        "role": "audience",
                        "parent": "left card",
                        "bbox": [180, 520, 820, 680],
                        "font_px": 40,
                        "text": "Readable audience text",
                    },
                    {
                        "name": "right card body",
                        "kind": "text",
                        "role": "audience",
                        "parent": "right card",
                        "bbox": [1080, 520, 1720, 680],
                        "font_px": 40,
                        "text": "Readable audience text",
                    },
                ],
            }
            layout = root / "layout_nodes.json"
            layout.write_text(json.dumps(nodes), encoding="utf-8")
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--layout-nodes",
                    str(layout),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
            payload = json.loads(report.read_text(encoding="utf-8"))
            self.assertTrue(payload["ok"])

    def test_reference_style_title_minimum_and_anchor(self) -> None:
        audit = SymmetryAudit(name="title-test")
        nodes = [
            {
                "name": "slide title",
                "kind": "text",
                "role": "title",
                "bbox": [220, 160, 1200, 230],
                "font_px": 60,
                "text": "Too small and too low",
            },
        ]
        audit_layout_nodes(
            audit,
            nodes,
            min_audience_font_px=36,
            min_title_font_px=67,
            title_left_max_px=180,
            title_top_max_px=130,
        )
        self.assertFalse(audit.ok)
        self.assertTrue(any(check.kind == "title_anchor" and not check.ok for check in audit.checks))
        self.assertTrue(any(check.kind == "font_size_minimum" and not check.ok for check in audit.checks))

    def test_text_overflow_past_canvas_fails(self) -> None:
        audit = SymmetryAudit(name="canvas-containment-test")
        nodes = [
            {
                "name": "overflow title",
                "kind": "text",
                "role": "title",
                "bbox": [80, 60, 2570, 140],
                "font_px": 80,
                "text": "This title extends past the right slide edge",
            },
            {
                "name": "contained body",
                "kind": "text",
                "role": "audience",
                "bbox": [120, 300, 1200, 360],
                "font_px": 42,
                "text": "Contained text",
            },
        ]
        audit_layout_nodes(audit, nodes)
        failures = [check for check in audit.checks if check.kind == "canvas_containment" and not check.ok]
        self.assertEqual(len(failures), 1)
        self.assertIn("overflow title", failures[0].name)

    def test_post_render_reports_subtitle_in_title_bottom_band(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            nodes = {
                "minimum_audience_font_px": 36,
                "minimum_title_font_px": 67,
                "nodes": [
                    {
                        "name": "slide title",
                        "kind": "text",
                        "role": "title",
                        "bbox": [80, 40, 1800, 125],
                        "font_px": 80,
                        "text": "Large reference title",
                    },
                    {
                        "name": "slide subtitle",
                        "kind": "text",
                        "role": "audience",
                        "bbox": [82, 134, 1900, 176],
                        "font_px": 40,
                        "text": "Subtitle text too close under the title",
                    },
                ],
            }
            layout = root / "layout_nodes.json"
            layout.write_text(json.dumps(nodes), encoding="utf-8")
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--layout-nodes",
                    str(layout),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            failures = [check for check in payload["checks"] if check["kind"] == "title_bottom_clearance" and not check["ok"]]
            self.assertEqual(len(failures), 1)

    def test_post_render_reports_excessive_title_followup_gap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            nodes = {
                "minimum_audience_font_px": 36,
                "minimum_title_font_px": 67,
                "nodes": [
                    {
                        "name": "slide title",
                        "kind": "text",
                        "role": "title",
                        "bbox": [80, 40, 1400, 125],
                        "font_px": 80,
                        "text": "Large reference title",
                    },
                    {
                        "name": "first audience line",
                        "kind": "text",
                        "role": "audience",
                        "bbox": [82, 390, 1700, 440],
                        "font_px": 40,
                        "text": "Audience line sits too far below title",
                    },
                ],
            }
            layout = root / "layout_nodes.json"
            layout.write_text(json.dumps(nodes), encoding="utf-8")
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--layout-nodes",
                    str(layout),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            failures = [check for check in payload["checks"] if check["kind"] == "title_followup_gap" and not check["ok"]]
            self.assertEqual(len(failures), 1)

    def test_colon_label_requires_bold_lead_and_regular_body(self) -> None:
        audit = SymmetryAudit(name="colon-style-test")
        nodes = [
            {"name": "card", "kind": "panel", "bbox": [0, 0, 300, 120]},
            {
                "name": "good label body",
                "kind": "text",
                "role": "audience",
                "parent": "card",
                "bbox": [20, 40, 280, 80],
                "font_px": 40,
                "text": "Signal: isolated photon candidates",
                "text_runs": [
                    {"text": "Signal:", "bold": True},
                    {"text": " isolated photon candidates", "bold": False},
                ],
            },
            {
                "name": "bad label body",
                "kind": "text",
                "role": "audience",
                "parent": "card",
                "bbox": [20, 40, 280, 80],
                "font_px": 40,
                "text": "Inclusive: background-rich candidates",
                "text_runs": [
                    {"text": "Inclusive:", "bold": False},
                    {"text": " background-rich candidates", "bold": True},
                ],
            },
        ]
        audit_layout_nodes(audit, nodes, min_audience_font_px=36)
        checks = [check for check in audit.checks if check.kind == "colon_label_style"]
        self.assertEqual(len(checks), 2)
        self.assertTrue(checks[0].ok)
        self.assertFalse(checks[1].ok)

    def test_plot_annotation_font_floor_and_color_difference_reason(self) -> None:
        audit = SymmetryAudit(name="plot-annotation-test")
        nodes = [
            {
                "name": "panel stat",
                "kind": "text",
                "role": "plot_annotation",
                "bbox": [10, 10, 110, 40],
                "font_px": 24,
                "text": "N=1.2M",
            },
            {
                "name": "left card",
                "kind": "panel",
                "bbox": [100, 300, 900, 900],
                "symmetry_group": "cards",
                "fill_color": "#eeeeee",
                "edge_color": "#cccccc",
                "color_difference_intentional": True,
                "color_difference_reason": "read-instructions vs interpretation",
            },
            {
                "name": "right card",
                "kind": "panel",
                "bbox": [1000, 300, 1800, 900],
                "symmetry_group": "cards",
                "fill_color": "#fff8df",
                "edge_color": "#d8c67f",
                "color_difference_intentional": True,
                "color_difference_reason": "read-instructions vs interpretation",
            },
        ]
        audit_layout_nodes(audit, nodes, min_plot_annotation_font_px=28)
        self.assertTrue(any(check.kind == "font_size_minimum" and not check.ok for check in audit.checks))
        self.assertTrue(any(check.kind == "symmetry_group_color" and check.ok for check in audit.checks))

    def test_boxed_equation_uses_bisector_checks(self) -> None:
        audit = SymmetryAudit(name="equation-bisector-test")
        nodes = [
            {"name": "equation box", "kind": "panel", "bbox": [0, 0, 300, 160]},
            {
                "name": "equation ink",
                "kind": "equation_ink",
                "parent": "equation box",
                "bbox": [30, 45, 270, 115],
            },
        ]
        audit_layout_nodes(audit, nodes)
        self.assertTrue(any("horizontally" in check.name for check in audit.checks))
        self.assertTrue(any("vertically" in check.name for check in audit.checks))
        self.assertTrue(audit.ok)

    def test_title_axis_alignment_reports_drift_below_title(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            nodes = {
                "title_axis_x": 74,
                "title_axis_tolerance_px": 2,
                "nodes": [
                    {
                        "name": "slide title",
                        "kind": "text",
                        "role": "title",
                        "bbox": [74, 30, 1800, 120],
                        "font_px": 86,
                        "text": "Aligned title",
                    },
                    {
                        "name": "good frame",
                        "kind": "panel",
                        "bbox": [74, 220, 1200, 900],
                        "title_axis_align": "left",
                    },
                    {
                        "name": "drifted frame",
                        "kind": "panel",
                        "bbox": [54, 220, 1180, 900],
                        "title_axis_align": "left",
                    },
                ],
            }
            layout = root / "layout_nodes.json"
            layout.write_text(json.dumps(nodes), encoding="utf-8")
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--layout-nodes",
                    str(layout),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            failures = [check for check in payload["checks"] if check["kind"] == "title_axis_alignment" and not check["ok"]]
            self.assertEqual(len(failures), 1)
            self.assertIn("drifted frame", failures[0]["name"])

    def test_text_block_vertical_fill_reports_underfilled_card(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            nodes = {
                "nodes": [
                    {
                        "name": "slide title",
                        "kind": "text",
                        "role": "title",
                        "bbox": [74, 30, 1800, 120],
                        "font_px": 86,
                        "text": "Aligned title",
                    },
                    {
                        "name": "takeaway card",
                        "kind": "card",
                        "bbox": [74, 1000, 2486, 1280],
                    },
                    {
                        "name": "takeaway text block",
                        "kind": "text",
                        "role": "audience",
                        "parent": "takeaway card",
                        "bbox": [150, 1090, 1900, 1180],
                        "font_px": 44,
                        "text": "Two short bullets",
                        "vertical_fill_min_ratio": 0.5,
                        "vertical_fill_max_ratio": 0.7,
                    },
                ],
            }
            layout = root / "layout_nodes.json"
            layout.write_text(json.dumps(nodes), encoding="utf-8")
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--layout-nodes",
                    str(layout),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            failures = [check for check in payload["checks"] if check["kind"] == "text_block_vertical_fill" and not check["ok"]]
            self.assertEqual(len(failures), 1)

    def test_outer_frame_top_gap_reports_loose_framed_plot(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            nodes = {
                "nodes": [
                    {
                        "name": "slide title",
                        "kind": "text",
                        "role": "title",
                        "bbox": [74, 30, 1800, 120],
                        "font_px": 86,
                        "text": "Aligned title",
                    },
                    {"name": "plot frame", "kind": "panel", "bbox": [74, 220, 1275, 1060]},
                    {
                        "name": "plot image ink",
                        "kind": "figure_ink",
                        "parent": "plot usable box",
                        "outer_frame": "plot frame",
                        "max_outer_frame_top_gap_px": 28,
                        "bbox": [94, 270, 1255, 1030],
                    },
                    {"name": "plot usable box", "kind": "figure_box", "bbox": [94, 240, 1255, 1040]},
                ],
            }
            layout = root / "layout_nodes.json"
            layout.write_text(json.dumps(nodes), encoding="utf-8")
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--layout-nodes",
                    str(layout),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            failures = [check for check in payload["checks"] if check["kind"] == "outer_frame_top_gap" and not check["ok"]]
            self.assertEqual(len(failures), 1)

    def test_vertical_margin_balance_reports_uneven_header_footer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            nodes = {
                "slide_size": [2560, 1440],
                "vertical_margin_balance": {
                    "top_node": "slide title",
                    "bottom_node": "takeaway card",
                    "tolerance_px": 4,
                },
                "nodes": [
                    {
                        "name": "slide title",
                        "kind": "text",
                        "role": "title",
                        "bbox": [74, 47, 1800, 123],
                        "font_px": 86,
                        "text": "Aligned title",
                    },
                    {
                        "name": "takeaway card",
                        "kind": "card",
                        "bbox": [74, 1098, 2486, 1370],
                    },
                ],
            }
            layout = root / "layout_nodes.json"
            layout.write_text(json.dumps(nodes), encoding="utf-8")
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--layout-nodes",
                    str(layout),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            failures = [check for check in payload["checks"] if check["kind"] == "vertical_margin_balance" and not check["ok"]]
            self.assertEqual(len(failures), 1)

    def test_vertical_margin_target_reports_mismatch_to_side_margin(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            nodes = {
                "slide_size": [2560, 1440],
                "vertical_margin_balance": {
                    "top_node": "slide title",
                    "bottom_node": "takeaway card",
                    "target_gap_px": 74,
                    "tolerance_px": 4,
                },
                "nodes": [
                    {
                        "name": "slide title",
                        "kind": "text",
                        "role": "title",
                        "bbox": [74, 47, 1800, 123],
                        "font_px": 86,
                        "text": "Aligned title",
                    },
                    {
                        "name": "takeaway card",
                        "kind": "card",
                        "bbox": [74, 1098, 2486, 1393],
                    },
                ],
            }
            layout = root / "layout_nodes.json"
            layout.write_text(json.dumps(nodes), encoding="utf-8")
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--layout-nodes",
                    str(layout),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            failures = [check for check in payload["checks"] if check["kind"] == "vertical_margin_target" and not check["ok"]]
            self.assertEqual(len(failures), 1)

    def test_white_slide_backdrop_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "white").save(png)
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr + completed.stdout)
            payload = json.loads(report.read_text(encoding="utf-8"))
            checks = [check for check in payload["checks"] if check["kind"] == "white_backdrop"]
            self.assertEqual(len(checks), 1)
            self.assertTrue(checks[0]["ok"])

    def test_tinted_slide_backdrop_fails_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            png = root / "candidate.png"
            Image.new("RGB", (2560, 1440), "#f5f7fb").save(png)
            report = root / "audit.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(POST_RENDER),
                    "--png",
                    str(png),
                    "--output",
                    str(report),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(report.read_text(encoding="utf-8"))
            failures = [check for check in payload["checks"] if check["kind"] == "white_backdrop" and not check["ok"]]
            self.assertEqual(len(failures), 1)


if __name__ == "__main__":
    unittest.main()
