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
- Every generated full-slide PNG candidate should have a companion speaker
  script Markdown file in the same output family. The script is not an outline
  or private presenter note; it is the words Justin can practice saying out
  loud on that slide. The final chat response must include clickable links to
  both the PNG and the script.
- Do not bake tiny `Source:` / provenance footers into full-slide PNG
  candidates. Slides should be audience-facing by default; keep detailed
  provenance in the chat note, generated summary CSV/JSON, speaker notes, or a
  separate technical backup slide when needed.
- Do not put internal design notes, presenter instructions, implementation
  reminders, or text addressed to Justin/Codex on the slide canvas. Phrases
  like "read this slide as", "purpose", "speaker note", "next action", or
  internal QA/provenance language belong in the companion script, manifest, or
  chat note. On-slide text must be audience-facing: title, labels, definitions,
  claims, evidence, caveats the audience needs, and transition language only
  when it is meant to be spoken or seen by the audience.
- Do not use caret notation on a final slide canvas for text that should be
  superscripted or subscripted. Units such as inverse picobarns, inverse
  nanobarns, `p_T`, `E_T`, `x_J`, `x_{Jgamma}`, powers, indices, and isotope or
  particle labels must be rendered with true visual superscript/subscript in
  the PNG or native Slides object. If the renderer cannot format it natively,
  draw separate text runs with adjusted size and baseline. Caret notation is
  acceptable only in manifests, code, or plain-text provenance where visual
  typography is not the output.
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

## Speaker Script Contract

Use Justin's prelim narrative style as the reference voice for generated slide
scripts. The named exemplar is the Google Doc `Prelim Narrative & comments`
(`1uzYK-8p-mCeHjtqjT0hewC20lXgzgfeu77Xsb2gZAQ8`), especially the narrative for
the linked prelim `sPHENIX Subsystems` slide
(`1FlmFhtPDnp3za_I5v4yLL8gSaTeuaJxzywDF2MKmc14`, slide
`g3920ab4f8cb_4_310`).

Generated scripts should read like polished spoken narration, not like generic
speaker notes:

- Start with the exact orientation Justin would say on the slide, using natural
  first-person transitions when appropriate: for example, "At a high level",
  "To start us off", "Now", "First", "Second", "So overall", or "With that
  context in mind".
- Walk through the slide in the visual order a listener will experience it.
  Explain acronyms and definitions inline at the moment they appear.
- Keep the voice audience-facing and conversational but technically precise.
  Prefer complete spoken sentences over fragments, labels, or stage directions.
- Personalize when it matters: use "my measurement", "for this talk", or
  "what I want to emphasize" when that matches Justin's ownership of the work.
- End with a spoken takeaway or bridge that naturally carries the audience to
  the next slide. Do not use meta labels such as `Transition:` in normal
  standalone slide scripts; write the transition as a sentence Justin can say.
- Keep private reasoning, implementation notes, provenance, QA, and design
  intent out of the script body unless Justin would actually say it. Put those
  in the manifest, chat handoff, or separate technical notes.
- Format each `.md` script as `# <talk/deck> Slide <n> Script - <short title>`
  followed by short paragraphs in spoken order. Bullets are allowed only when
  they represent a natural spoken list, not an outline shortcut.

For progressive builds, keep any click or timing cues clearly separated from
the spoken script text. The default standalone-slide script should be directly
readable aloud from top to bottom.

## Progressive Build Option

Prefer separate, self-contained full-slide PNG candidates. This is the normal
and recommended presentation form, especially for PPG, JSTG,
collaboration-facing updates, and analysis-status decks, where Justin needs
simple navigation, clean screenshot/export behavior, easy slide-by-slide
review, and minimal presenter overhead.

Codex should not silently choose a progressive build sequence. A progressive
build, such as `Slide 2A` / `Slide 2B` or a short three-frame reveal, is allowed
only when Justin explicitly asks for it, or when Codex presents it as an
optional alternate beside the recommended standalone-slide version.

Use a progressive alternate only when it clearly beats separate slides for live
delivery:

- the same visual anchor must remain on screen while a single idea is revealed
  in stages;
- revealing all information at once would materially overload the audience;
- the click timing is essential to the explanation, not merely decorative;
- the sequence still behaves like one spoken slide, not a hidden extra section.

Do not use progressive builds because a slide has many objects, because
animation would look polished, because a dense slide can be decomposed, or
because it is possible. If the story is equally clear as normal consecutive
slides, generate normal consecutive slides.

When a progressive alternate is generated:

- label it as optional, not the default recommendation;
- state why it might help and why the standalone version remains safer unless
  Justin chooses the build;
- generate each frame as an individual `2560x1440` PNG candidate;
- provide a contact sheet, speaker-click script, and manifest;
- keep each frame self-contained enough to present if the sequence is later
  split into normal slides;
- do not mutate Google Slides until Justin approves the sequence and the
  insertion mode.

For approved deck insertion, treat a progressive PNG build as consecutive
duplicate-frame slides or a native Slides build depending on Justin's chosen
workflow. Do not bake slide numbers or provenance footers into any frame.

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
