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

If `daily_cockpit.context_conservation_mode.active` is true, keep the slide
handoff compressed: generate/show only the updated `2560x1440` PNG candidate
plus a terse source/caveat note unless Justin explicitly asks for more. Do not
mutate Google Slides in conservation mode unless Justin explicitly asks to
update/insert/replace the deck after reviewing the PNG. Companion speaker
scripts are optional in this mode unless Justin asks or the slide is for an
imminent talk where missing narration would create real delivery risk.

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
- Use a source-first plot policy for all slide generation, not only public
  conference talks. When an approved paper figure, collaboration-approved plot,
  public reference plot, or current vetted analysis plot already exists and is
  the correct evidence object, place that exact visual artifact on the slide
  rather than recreating it from scratch. Generated design should improve the
  presentation around the evidence object: crop, framing, shadow, layout,
  callout cards, arrows, labels, footer/header consistency, and explanatory
  scaffolding. Do not impersonate a collaboration/result plot with newly
  generated fake data, approximated curves, or decorative re-draws. If a plot
  must be regenerated from data, record the exact code, input files, tag,
  timestamps, and validation evidence in the manifest or handoff.
- Visual truthfulness is a hard requirement. Never make slide labels,
  legends, captions, callouts, titles, or spoken claims say something stronger
  or different than the actual plotted object, source file, selection, sample,
  weight, model, run tag, or approval state. If a slide uses pseudo-data,
  design scaffolding, a placeholder, an approximate recreation, an internal-only
  source, or an unvalidated output, label it that way clearly in the handoff
  and, when audience-visible, on the slide itself. If the true object cannot be
  verified, do not call the candidate slide-ready.
- Apply the same visual-quality bar to collaboration-facing decks, PPG/JSTG
  updates, internal analysis reviews, and public talks. Public talks have
  stricter release/label constraints, but internal audiences still deserve
  clean hierarchy, readable text, consistent typography, uncluttered
  composition, and digestible claim/evidence/implication structure.
- Audience readability is a hard baseline, not a final polish pass. For
  generated 16:9 `2560x1440` full-slide PNG candidates, use title typography
  equivalent to at least 24-26 pt, and larger when white space allows. Section
  labels, card titles, plot labels, and text-box lead lines should be
  equivalent to at least 16 pt. Body text inside callout cards or colored text
  boxes should be equivalent to at least 13-14 pt. If the text does not fit at
  those sizes, reduce wording, split the idea across slides, or change the
  layout; do not silently shrink audience-facing text below readability.
- Full-slide candidates should balance context, specificity, digestibility, and
  visual hierarchy without requiring Justin to catch repeated readability
  problems. Prefer fewer, sharper audience-facing statements over dense
  paragraphs. Each slide should make clear what object is shown, why it matters,
  and what conclusion the audience should carry forward.
- In explanatory boxes that define categories, samples, symbols, or comparison
  logic, structure the text as readable rows rather than one paragraph. Bold
  the lead label, for example `Signal =` or `Inclusive =`, put distinct
  definitions on separate lines, and separate the concluding question or
  takeaway into its own line. Avoid burying multiple definitions and the key
  interpretation in a single wrapped chunk.
- For panel/card shadows, use the Google Slides-style baseline unless a
  specific visual reason justifies a variation: black shadow, opacity 18%,
  angle 60 degrees, distance 4 px, blur radius 14 px. Avoid heavy drop shadows
  that compete with the plot or make internal analysis slides look decorative.
- For public HP2026 PPG12 photon slides, final PPG12 plot images must come
  from the current PPG12 paper PDF or from another explicitly public/approved
  source recorded in the manifest. IAN plots, internal ROOT/data-generated
  plots, collaborator-only figures, and screenshot crops from backup slides are
  placeholders only until released or explicitly approved for public use. If a
  candidate still shows `sPHENIX Internal` where the public version should say
  `Preliminary` or carry an approved public label, mark the slide incomplete
  rather than public-ready.
- Do not put internal design notes, presenter instructions, implementation
  reminders, or text addressed to Justin/Codex on the slide canvas. Phrases
  like "read this slide as", "purpose", "speaker note", "next action", or
  internal QA/provenance language belong in the companion script, manifest, or
  chat note. On-slide text must be audience-facing: title, labels, definitions,
  claims, evidence, caveats the audience needs, and transition language only
  when it is meant to be spoken or seen by the audience.
- For explanatory note boxes, target/result boxes, and compact interpretation
  cards, prefer a short bold lead label followed by regular-weight body text
  rather than bolding the entire sentence. Examples: `Note:` bold with the note
  body regular; `Target result:` bold with the result description regular.
  Use all-bold text only for short labels, card titles, or genuinely emphatic
  one-line takeaways.
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

### Progressive Focus Build

A `progressive focus build` is a narrow progressive-build subtype where the
same slide layout remains fixed while attention moves from one proof object to
the next. Use this only when the visual object is dense enough that showing all
panels at full emphasis would compete with the spoken sequence.

The pattern is:

- keep the exact same title, layout, footer, axes, and object placement across
  all frames;
- leave the current proof object at full contrast;
- ghost, desaturate, blur, or soften future proof objects so they remain
  visible as context but are not read yet;
- on each click, activate one additional object or annotation while preserving
  the same visual field;
- avoid adding internal labels such as `next`, `not yet`, or `click here` on
  the slide canvas.

This is useful for slides such as a BDT pipeline, ABCD sideband logic,
corrections chain, or systematics breakdown, where Justin is walking the
audience through the same diagram step by step. It is not a general animation
preference: if the content is clearer as normal consecutive standalone slides,
generate normal slides.

When generating a progressive focus build, produce each frame as an individual
`2560x1440` PNG, plus a click-aware script or clearly named per-frame scripts.
The chat handoff should state which frame is active at each step and should
include links to every PNG/script. Google Slides remains untouched until Justin
explicitly approves insertion.

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

## Iteration Learning Loop

When Justin reacts to a generated plot, full-slide candidate, or speaker script
with feedback such as "wrong", "too busy", "hard to read", "not what I meant",
"make it cleaner", "do it like yesterday", or repeated geometry/style
corrections, treat the feedback as both an artifact fix and a reusable training
signal.

For each non-trivial iteration:

1. Fix the requested artifact directly and keep the response focused on the
   requested change.
2. Identify the general failure mode in one sentence: for example, unreadable
   hierarchy, over-dense text, weak scan path, wrong plot/source object,
   bad legend placement, mismatched comparison definition, too much blue/bold
   styling, or failure to preserve a proven layout.
3. Decide whether the lesson is local or durable. Local lessons stay in the
   slide/script manifest or chat handoff; durable lessons go into
   `agent_context/SLIDE_STYLE_MAP.md`, this policy, `PLOTTING.md`, or the
   relevant generator defaults.
4. If the same correction appears twice across nearby iterations, promote it
   from a one-off preference into a default checklist item for that plot family
   or slide type.
5. Do not over-record every aesthetic opinion. Record only rules that reduce
   future repeated correction, prevent scientific miscommunication, or improve
   Justin's visual parsing.

The practical question after any failed candidate is: "What should future
Codex do differently so Justin does not need to make this same correction
again?" Apply that answer before generating the next version when it is safe
and scoped.

## Slide Candidate Self-Audit

Before showing Justin a serious slide-facing artifact, run this audit and fix
failures before handoff:

- correct source object: approved/reference plot used when one exists; no
  accidental recreation of an approved artifact;
- audience-facing canvas: no internal notes, Codex instructions, provenance
  minutiae, implementation reminders, or private next actions on the slide;
- readable typography: title, section labels, card text, axes, legend, and
  annotations are readable at meeting scale;
- true superscript/subscript: no caret notation for visual units or physics
  labels that should be typeset;
- clean hierarchy and whitespace: one main scan path, no box collisions,
  legend-axis crowding, tiny captions, or over-dense callout text;
- slide numbers/provenance: no baked slide numbers and no tiny visible
  `Source:` footers unless Justin explicitly asks;
- speaker script: companion script exists and reads as Justin's spoken
  narration, not as a generic outline or private notes;
- mutation boundary: Google Slides is untouched unless Justin approved deck
  mutation for the current target deck/slide.

If an audit item fails, regenerate or revise before presenting the artifact as
ready. If the failure repeats across iterations, promote it through the
feedback-ingestion loop.

## Feedback Ingestion

After Justin correction, record the smallest lesson that prevents recurrence.

- Classify the correction: style/narrative, physics/science,
  provenance/evidence, workflow/safety, retrieval failure, memory salience
  failure, task/status failure, tool/runtime failure, or artifact-quality
  failure.
- Decide local versus durable. Local corrections go in the artifact manifest or
  handoff. Durable corrections go in `agent_context/SLIDE_STYLE_MAP.md`, this
  policy, `PLOTTING.md`, a generator default, a negative memory trap, or a
  validator.
- Cite evidence: user correction, slide ID, PNG path, deck link, source plot,
  command output, manifest, or visual QA note.
- Do not record every small taste change. Record only lessons that reduce
  repeated user correction, prevent source/provenance mistakes, or protect the
  science story.

## Iteration Learning Loop

For serious full-slide generation or repeated slide revision, maintain a local
iteration ledger instead of leaving the learning trail in chat.

Use:

```bash
python3 scripts/os/delegation/slide_iteration_learning.py init \
  --slide-key "<deck-or-talk>-slide-<n>-<topic>" \
  --objective "<initial slide conception>"
```

The helper writes local ignored records under
`agent_context/local/slide_iterations/`. Record:

- initial conception and intended audience;
- each generated candidate PNG/script/generator attempt;
- each Justin correction, grouped by category rather than as raw annoyance;
- each Codex self-audit or Claude worker critique;
- the accepted artifact and the iteration number where it became acceptable.

After enough evidence exists, ask Claude to distill the ledger:

```bash
python3 scripts/os/delegation/slide_iteration_learning.py pack-claude \
  --run-dir agent_context/local/slide_iterations/<run_id> \
  --invoke
```

Claude's role is to return a `Slide Iteration Learning Digest`: iteration
metrics, recurring failure modes, prevention rules, a compact Codex injection
brief, and promotion candidates. This is not model-weight training. It is
structured retrieval/policy learning. Codex must verify the digest before using
it to change `SLIDE_STYLE_MAP.md`, this policy, a generator default, a
validator, or a future slide context pack.

Before generating the next related slide, Codex should load the relevant
`iteration_summary.md` or Claude digest when available and paste only the
compact `Codex Injection Brief` into the new slide-generation context. The goal
is to reduce repeated slide iterations by preventing the same concrete failure,
not to overfit every slide to the last correction.

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
