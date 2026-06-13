# Shared Slide Helpers

Use this folder for layout constants and small reusable helpers used by local
full-slide PNG generators. New slide code should prefer these defaults:

- `2560x1440` PNG output for full-slide candidates.
- Times New Roman text when composing audience-facing slide PNGs.
- No baked slide/page numbers unless Justin explicitly asks.
- No tiny provenance footers; keep evidence in chat notes, JSON/CSV sidecars,
  speaker notes, or backup material.

## Default Post-Render Audit

Use `post_render_slide_audit.py` after rendering serious slide candidates. It
checks the rendered PNG and, when the generator emits `layout_nodes.json`, the
same transferable geometry contracts that proved useful in the HP sequence:
audience text size, containment, box/text centering, repeated-card symmetry,
and internal/provenance text leakage.

Example:

```bash
python3 scripts/slides/common/post_render_slide_audit.py \
  --png dataOutput/.../candidate.png \
  --layout-nodes dataOutput/.../layout_nodes.json
python3 scripts/os/safety/codex_os_guard.py after-slide-render \
  --png dataOutput/.../candidate.png \
  --symmetry-report dataOutput/.../candidate.slide_audit.json
```

The audit is advisory-but-required for serious handoff: pass it, fix it, or
record a specific exception. Run it before any slide-worker subagent call. For
the first serious generation of a slide identity, fix checker findings first,
then call the subagent once, fix valid subagent findings, and rerun the checker
before showing Justin. For later iterations of that same slide, rerun this
checker every time and tune the checker/prompt from Justin feedback, but do
not re-call the subagent unless Justin explicitly asks or the slide identity
has materially changed. Deck-specific chrome checks, such as HP2026
header/footer identity, remain opt-in and should not become global defaults.

## Symmetry Audit

Use `slide_symmetry_audit.py` for serious slide-facing generators that draw
cards, equations, schematic objects, or progressive build states. Treat the
slide as a layout graph:

- panels, boxes, diagrams, text bands, and footer strips are nodes;
- centering, equal padding, equal visual cells, containment, and pixel identity
  across build states are edges;
- the generator or a paired campaign auditor should write a JSON report beside
  the PNG outputs.
- text nodes should include `font_px`, `bbox`, and `parent` when practical;
  use `intentional_alignment: top` only when top alignment is deliberate.
- any box that surrounds text, equations, images, figures, or diagrams should
  have both an outer node and an inner ink/text node so the checker can enforce
  the usable-band bisector method rather than relying on visual vibes.

This is the deterministic foundation for the slide-learning loop. Future
feedback can add new checks, but the baseline should already prove that objects
sit on their intended usable-band bisectors before the PNG is shown.
