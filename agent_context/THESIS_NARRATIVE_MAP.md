# Thesis Narrative Map

Last updated: 2026-05-26

This is the durable spine for the thesis, talks, and analysis execution. It is
not a task board. It answers why each workstream exists and what evidence must
survive into the final story.

## North Star

Produce a defensible, high-impact sPHENIX photon and gamma-jet thesis whose
analysis chain is traceable from raw detector objects through photon ID,
background modeling, response/unfolding, and final nuclear-modification physics.

The scientific center is PPG19 Au+Au gamma-jet / `x_{J#gamma}`. PPG18 pp
gamma-jet is the validated baseline. PPG12 supplies the prompt-photon ID,
isolation, purity, BDT/NPB, and pp-infrastructure reference layer.

## Justin's Position

Justin is both physics owner and operating-system owner for this work: PPG19
chair, machine-learning engineer, analysis producer, slide/narrative author,
and collaborator-facing integrator. The OS exists because this role combines
science, production, ML, presentation, and project management in one person.

Codex should reduce cognitive load by keeping provenance, tasks, plots, jobs,
emails, and slides synchronized. Codex should not replace Justin's physics
judgment.

## Claim Hierarchy

| Layer | Thesis claim | Required evidence | Current OS surfaces |
| --- | --- | --- | --- |
| Physics target | sPHENIX can make a serious prompt-photon / gamma-jet measurement in the RHIC environment. | Detector/object definitions, photon-ID performance, isolation, purity, trigger/centrality context, approval-ready plots. | `hp26_photon_id_talk`, PPG12 notes, working-point deck. |
| pp baseline | The pp prompt-photon and gamma-jet reference is reproduced or mapped to validated PPG12/PPG18 infrastructure. | PPG12-style photon/inclusive stitching, ET/eta training controls, Shuhang/PPG12 overlay, pp score-shape validation, accepted pp slide assets. | `pp_basev3e_bdt_shuhang_overlay`, pp stitching slide tasks, `DATASET_STATUS.md`. |
| Embedded background | Embedded inclusive backgrounds are stitched with correct ownership, cross sections, and sample support. | Jet12/20/30 baseline, Jet40 diagnostic proof, effective cross sections, final ROOT/CSV/PNG evidence, invalid-output warnings. | `jet40_embedded_inclusive_stitching`, `DATASET_STATUS.md`, `STATUS_DASHBOARD.md`. |
| ML photon ID | The final photon-ID model is chosen by closure and fake-rate behavior, not AUC alone. | WP80 fake rates, isolation/ABCD closure, score-isolation correlations, feature studies, route/bin support, train/test and seed robustness. | `ml_final_model_ablation_map`, validation reports, BDT/MLP/stack diagnostics. |
| Au+Au response | Signal/background/reco-response inputs are valid enough for final gamma-jet physics. | Clean input contracts, response matrices, centrality/trigger validity, invalid centrality extension gate, unfolding readiness. | `CURRENT_PROJECT_STATE.md`, scaled-trigger validity notes, response/unfolding workstreams. |
| Final physics | The final `x_{J#gamma}` or related gamma-jet observable is interpretable, uncertainty-aware, and presentation-ready. | Unfolded distributions, systematic checks, pp/AuAu comparison, physics interpretation, slide/storyline. | future final-analysis workstreams and approval deck. |

## Non-Negotiable Evidence Rules

- No slide-ready claim without artifact path and visual QA.
- No production-ready dataset without ROOT/key/count/provenance evidence.
- No stitched-background plot without matching ownership windows and effective
  cross sections.
- No ML promotion by inclusive AUC alone; closure and WP80 fake-rate behavior
  decide.
- No broad rerun without duplicate guard and explicit user approval.
- No Google Slides mutation before a local candidate is shown unless Justin
  explicitly asks for direct deck editing.

## Current Sharp Edges

- The pp Shuhang overlay is a diagnostic discrepancy, not a solved validation
  claim, until selection/histogram/bin/model differences are understood.
- Jet40 extended-50 is validated for the targeted diagnostic scope, but the
  default production path should not silently change to Jet12+20+30+40 without a
  dedicated promotion task.
- The scaled-trigger centrality extension remains invalid until a one-run
  inclusive regression reproduces the validated inclusive result.
- Active HP2026 slide work should privilege accepted PPG12/PPG18 material first,
  then use new Justin-made plots only when they close an actual gap.

## Symbiotic Working Model

Justin:

- sets physics priority, taste, narrative ambition, and collaborator judgment;
- spots when an argument feels wrong or a plot does not tell the real story;
- decides when evidence is strong enough to present.

Codex:

- keeps the OS state coherent;
- gathers exact evidence;
- prevents duplicate/stale/risky operations;
- builds repeatable plots/scripts/checks;
- converts repeated friction into policy or automation;
- keeps the daily cockpit short and the archive durable.

## Promotion Gate

Before any work moves from analysis artifact to thesis/slide spine, it needs:

1. source data path and generation command or script;
2. physics claim it supports;
3. known limitations;
4. visual/readback QA if it is a plot or slide;
5. destination deck/chapter/section;
6. reopen condition if later evidence invalidates it.

## Evolution Rule

If a future task cannot be placed on this map, it is either:

- exploratory and belongs in backlog;
- an OS/tooling task and belongs in a separate OS workstream;
- a distraction that should be dropped;
- evidence that this map is missing a real thesis claim and should be updated.
