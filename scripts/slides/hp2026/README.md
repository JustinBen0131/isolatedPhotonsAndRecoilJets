# HP2026 Slide Generators

Home for Hard Probes 2026 talk-specific local slide generators. Keep reference
deck evidence and provenance beside generated artifacts, not as tiny text baked
into full-slide PNGs.

- `make_photon_detector_motif_variants.py`: standalone transparent PNG motif
  generator for the HP2026 title slide's photon-through-detector graphic.
- `opening/make_hp2026_opening_motivation_slide.py`: deterministic local
  `2560x1440` PNG generator for the two-slide opening sequence before any
  Google Slides insertion: Slide 2 detector/data context, then Slide 3
  isolated-prompt-photon motivation. Earlier one-slide opener variants remain
  available through explicit `--variant` choices for comparison. The
  recommended experiment/data refinement is
  `--variant experiment-dataset-standalone`, which renders two normal
  consecutive slides: sPHENIX experiment context first, p+p dataset context
  second. Use `--variant detector-data-progressive-build` only as an optional
  alternate when Justin explicitly wants a click-through build.
- `fulltalk/make_hp2026_fulltalk_candidates.py`: deterministic local
  `2560x1440` PNG batch generator for Slides 4-11 of the 11-slide HP2026 talk.
  Main physics/result plots are cropped from the current collaboration-facing
  PPG12 paper draft; generated graphics are explanatory scaffolding only.
