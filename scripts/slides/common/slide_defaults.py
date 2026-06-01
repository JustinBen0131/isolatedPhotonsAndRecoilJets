#!/usr/bin/env python3
"""Shared defaults for local full-slide PNG generators.

These constants encode the durable ThesisAnalysis slide-generation contract:
full-slide candidates are local 16:9 PNGs, audience-facing, and provenance is
kept beside the artifact rather than baked into tiny footers.
"""
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

SLIDE_WIDTH_PX = 2560
SLIDE_HEIGHT_PX = 1440
SLIDE_DPI = 200
SLIDE_ASPECT = (16, 9)
DEFAULT_FONT_FAMILY = "Times New Roman"
DEFAULT_OUTPUT_ROOT = "dataOutput/slide_assets"

DO_NOT_BAKE_SLIDE_NUMBERS = True
DO_NOT_BAKE_TINY_PROVENANCE_FOOTERS = True


def slide_figsize(dpi: int = SLIDE_DPI) -> tuple[float, float]:
    """Return the canonical matplotlib figure size for a full-slide PNG."""
    return (SLIDE_WIDTH_PX / dpi, SLIDE_HEIGHT_PX / dpi)
