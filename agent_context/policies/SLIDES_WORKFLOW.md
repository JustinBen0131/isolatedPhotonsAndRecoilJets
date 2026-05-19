# Slides Workflow

## Google Slides Safety

Prefer Google Drive/Slides connector reads, writes, and thumbnails. Browser or
Computer Use is not the default editing path.

Do not use Computer Use on Justin's browser/desktop for Google Slides unless:

1. it is genuinely useful for visual/layout checking or a narrow edit that
   cannot be done cleanly through the connector;
2. Codex first tells Justin exactly why Computer Use is needed, which deck/tab
   will be touched, and what actions will be performed;
3. Justin explicitly grants permission for that session.

Never switch to unrelated browser tabs, interact with email/calendar/chat,
submit forms, present, share, delete, or change permissions unless Justin asks.

## Backup Slide

Never remove, overwrite, or repurpose the deck's `Backup` slide. Backup is the
divider for backup material. Add new presentation slides before it and backup
material after it.

## Before Substantial Edits

- Inspect relevant past Drive presentations when available.
- Use `agent_context/SLIDE_STYLE_MAP.md` for durable style.
- Resolve deck ID, slide object ID, and relevant element IDs through connector
  reads.
- For plot replacements, first show generated PNGs in chat and wait for Justin
  to approve the image.
- Gemini / Beautify Slides may be used only as a design-critique or layout
  ideation assistant. Physics content and final layout decisions must come from
  Justin, evidence, and controlled Slides edits, never from Gemini as source of
  truth.

## Style Preferences

- Native Google Slides text uses Times New Roman: titles, body, callouts,
  captions, footers, page numbers, and editable labels.
- Prefer body text at least 14 pt when possible.
- Use line spacing around 1.2-1.8.
- Preferred update-slide style: larger 14+ pt body text, filled but not
  cramped whitespace, about 1.3 line spacing, and soft gray rounded/ovular
  callouts when helpful.
- Gray rounded callouts should avoid visible black borders; match outline to
  fill unless Justin asks otherwise.
- Colored callouts should usually be two objects: a no-text rounded background
  shape plus a transparent Times New Roman text box above it.
- RHS gray explanation boxes should use deliberate paragraph spacing.
- Bottom yellow takeaway bands should center text in a separate text box.
- For overview and conclusion slides, larger fonts and larger spacing are fine
  when they fill the space nicely without hurting clarity.

## Content Language

Match Justin's preferred slide language: concise bullets above plots, clear
hierarchy, bold lead phrases, organized sub-bullets when helpful, and a
claim/evidence/implication flow. Preserve the essence of Justin's wording while
making it easier to scan and explain aloud.

## Verification

After connector edits, fetch a fresh large thumbnail for every touched slide
and inspect it. Check no clipped text, overlaps, stale placeholders, off-slide
objects, or broken plot placement remain.
