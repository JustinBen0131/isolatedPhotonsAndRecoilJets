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

- If Justin asks about "the slides", "working point slides", "current slides",
  or a slide number without providing a deck link, first use the `WORKING
  POINT` line near the top of Today's Plan. If the live doc cannot be read,
  use `daily_cockpit.active_working_point_deck` in
  `agent_context/CODEX_WORK_REGISTER.yaml`. This working-point deck is the
  default slide context until Justin changes the line in Today's Plan or gives
  a different deck link in the prompt.
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

## Full-Slide Generation Contract

When Justin asks to `generate a slide`, interpret that as a request for a
polished full-slide replacement candidate, not merely a plot export or raw
Google Slides edit.

- Default deliverable: a complete 16:9 PNG candidate that composes the plot(s),
  title, explanatory text, and caveat/status line when needed into one
  slide-ready visual.
- For Google Slides full-slide PNG candidates in this project, render/export at
  `2560x1440` pixels unless Justin explicitly asks for a different size. Do
  not use lower-resolution `1920x1080` for slide-ready candidates.
- Regenerate or mutate the plot(s) for slide geometry: figure size, margins,
  legends, labels, panel spacing, and annotations should fit the slide cleanly
  instead of forcing a generic plot PNG into the deck.
- Use slide 22 from `WP_GammaJets_5_20_26`
  (`1-2l9IWNJMpFfunoShVrYk3ywtrIocMWMyePYnR85_zc`,
  slide `g3e2f3c3049c_1_298`) as the named exemplar for Times New Roman,
  plot-first composition, compact explanatory text, and polished internal
  analysis style.
- Visual quality comes before editability for the candidate PNG. Native Slides
  rebuilding is an optional follow-up after the candidate is approved.
- Show the candidate PNG in chat before mutating Google Slides. Include a short
  chat note with source/provenance, interpretation, caveat if applicable, and
  the next action.
- Do not bake tiny `Source:` / provenance footers into full-slide PNG
  candidates. Slides should be audience-facing by default; keep detailed
  provenance in the chat note, generated summary CSV/JSON, speaker notes, or a
  separate technical backup slide when needed.
- Do not bake slide numbers/page numbers into generated PNG candidates unless
  Justin explicitly asks. Justin manages slide numbering manually because the
  number changes with deck placement and backup-slide ordering.
- After approval, ask where to place it unless Justin already gave a target
  slide/location. Add or replace depending on the workflow Justin chooses.
- Preserve existing/backup slide state before major deck mutation. Never
  overwrite, repurpose, or remove a `Backup` divider slide.

These full-slide candidates are internal analysis artifacts as much as
presentation artifacts: they should linearize the analysis state by making clear
what is shown, why it matters, and what follows.

## Chronicle Context Rule

Chronicle is a context tool for ambiguous slide-history requests, not a source
of truth for final slide content.

- Use Chronicle when Justin refers to recent screen/deck work without enough
  detail to identify the slide, style, or plot.
- Once the target is identified, switch to Google Slides/Drive, local files,
  generated artifacts, and project ledgers for exact evidence and edits.
- Do not make Chronicle a mandatory step for every slide-generation task.

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

Before calling a full-slide candidate ready, render or open the PNG and inspect
it directly. Verify no cropped labels, cramped text, unreadable legends, awkward
white space, or stale plot/title language remain.

For dense explanatory analysis slides, also check the layout-contact points
that caused repeated slide-17 split-gain iterations:

- legends should not crowd headers or axis labels; move them into unused plot
  canvas space when that is cleaner;
- bottom callout boxes should not touch the upper plot/definition panels;
- side-by-side rounded callouts need a visible gutter and must not overlap;
- shadows should be removed when they make text or panel hierarchy worse;
- slide/page numbers should not be baked into the PNG unless Justin explicitly
  asks.

For analysis-check slides like the reco-cluster `E_T` leakage check, also apply
the slide-7 lessons recorded in `agent_context/SLIDE_STYLE_MAP.md`: use a
claim-style title, fill top callout boxes with large readable text, keep
`What is plotted` and `Main takeaway` as the primary reading path, put
sPHENIX/PYTHIA labels inside the plot canvas, use marker-only binned points
unless drawing a fit/reference, match binning across paired plots, and keep the
bottom status line short and readable rather than provenance-heavy.
