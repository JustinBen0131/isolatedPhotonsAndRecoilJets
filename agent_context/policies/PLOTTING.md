# Plotting Policy

## Defaults

- Write PNG files by default. Do not also write PDFs unless Justin asks.
- Never draw ROOT stat boxes unless Justin asks. Use `gStyle->SetOptStat(0)`
  and `SetStats(false)` on frame/data histograms.
- Reuse existing plotting infrastructure: `macros/AnalyzeRecoilJets*`,
  `macros/sPhenixStyle.*`, existing helpers, styles, axis labels, and output
  organization. After THE-23 stage 8, use `macros/MACRO_INDEX.yaml` or
  `macros/bin/thesis-macro path <id-or-basename>` to find canonical ROOT
  macro source before editing.
- Plot labels must name the actual comparison in plain physics terms. Avoid
  vague shorthand like `PPG12-like` when the legend can say what changed.

## File Placement

Before creating a new plotting script, slide-candidate generator, ROOT macro,
or diagnostic helper, load `REPO_ORGANIZATION.md`.

Default to organized side-code placement:

- reusable Python plotting utilities: `scripts/plotting/`;
- full-slide PNG builders and deck-specific slide helpers: `scripts/slides/`;
- focused audits or one-campaign diagnostics: `scripts/diagnostics/`;
- local ML/score/validation plot helpers: `scripts/ml/` when they are more ML
  than general plotting;
- offline ROOT plotting macros:
  `macros/plotting/{auau_bdt,target_wp,width_study,stitching,pp_currentian,ssqa}/`
  unless a protected runtime anchor or documented hard alias requires a
  top-level macro path.

Do not move Fun4All, SDCC, Condor, transfer, or base pipeline files while
making a plot. If a plot requires editing a protected entrypoint, keep the
existing path and follow the transfer/pipeline policy.

## sPHENIX Style

For final or slide-candidate plots, match `macros/sPhenixStyle.C` essentials:
font 42, no ROOT title/stat boxes, white background, borderless legends,
readable margins, ticks on both axes.

Canonical new-plot label:

```text
#it{#bf{sPHENIX}} Internal
```

In matplotlib, emulate this with bold italic `sPHENIX` plus upright
`Internal`. Never render plain upright `sPHENIX Internal`.

If editing an existing legacy macro whose nearby plots consistently use
`#bf{sPHENIX} #it{Internal}`, preserve that local convention rather than mixing
styles within one plot family.

## Visual QA

Any time Codex produces analysis plot images:

1. View the generated PNGs before calling them ready.
2. Check legends, TLatex labels, annotations, data markers, error bars, axis
   ranges, and text size.
3. If text or legends block data, rerun and tune coordinates/margins.
4. For templates used across many plots, inspect representative outputs from
   each template.
5. If a plot looks physically surprising, debug definitions, normalization,
   sample weights, numerator/denominator, ROOT input provenance, and histogram
   family before presenting it.

When Justin explicitly asked for generated plots, also open the containing
folder and preview the PNGs locally through the macOS `open` command when
practical, after viewing them in chat. If GUI opening requires escalation,
request narrow approval for only the needed `open` command.

## Slide Plot Contract

When Justin asks for plots for slides:

- show the candidate PNGs directly in chat first;
- provide clickable absolute paths to PNGs and folders;
- explain what is plotted, what is learned, why it does or does not belong on a
  slide, and any physics/formatting caveat;
- do not insert or replace Google Slides plots until Justin approves the image.

Generated slide-use PNGs should contain only the scientific plot: axes, data,
legend, and in-plot analysis labels. Do not bake slide titles, bullets,
takeaways, footers, or gray callouts into PNGs; those belong as editable Slides
text/shapes.

Exception: when Justin explicitly asks to `generate a slide`, follow
`agent_context/policies/SLIDES_WORKFLOW.md` instead. In that workflow the
deliverable is a polished full-slide 16:9 PNG candidate, so slide title,
narrative text, callouts, and plot geometry should be composed together and
shown in chat before any Google Slides mutation. Do not add tiny `Source:` or
provenance footers to the PNG; keep that evidence in the chat note or adjacent
summary artifacts.

For full-slide candidates, use plot-panel shadows when they improve hierarchy.
Preferred Google Slides shadow baseline: black shadow, opacity 18%, angle 60
degrees, distance 4 px, blur radius 14 px. Approximate these settings in
matplotlib or other renderers when exact Slides shadow controls are unavailable.

## Color And Readability

Use high-contrast, colorblind-aware palettes. Avoid pale yellow or low-contrast
critical lines/bars/heatmap cells. For heatmaps, annotation text should switch
between black and white based on cell brightness.
