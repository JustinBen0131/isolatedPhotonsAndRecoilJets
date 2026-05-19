# Plotting Policy

## Defaults

- Write PNG files by default. Do not also write PDFs unless Justin asks.
- Never draw ROOT stat boxes unless Justin asks. Use `gStyle->SetOptStat(0)`
  and `SetStats(false)` on frame/data histograms.
- Reuse existing plotting infrastructure: `macros/AnalyzeRecoilJets*`,
  `macros/sPhenixStyle.*`, existing helpers, styles, axis labels, and output
  organization.
- Plot labels must name the actual comparison in plain physics terms. Avoid
  vague shorthand like `PPG12-like` when the legend can say what changed.

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

## Color And Readability

Use high-contrast, colorblind-aware palettes. Avoid pale yellow or low-contrast
critical lines/bars/heatmap cells. For heatmaps, annotation text should switch
between black and white based on cell brightness.
