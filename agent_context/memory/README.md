# Memory Architecture

This folder contains tracked, neutral memory-routing registries for waking
Codex. It is not a transcript archive and it is not factual evidence by itself.
The goal is to route a task to the smallest useful old context, keep stale or
synthetic material cold, and learn whether retrieved context helped.

## Memory Temperatures

- `hot`: should surface quickly for common, high-risk, or high-value tasks.
- `warm`: useful in the right route, but not part of the normal boot surface.
- `cold`: valid but low-frequency or too specific for default retrieval.
- `quarantine`: stale, contradicted, synthetic, or provenance-weak until a
  waking check restores it.

## Retrieval Roles

- `conscious_context`: compact facts or policy pointers that should be loaded
  or verified before acting.
- `latent_nudge`: warning, method cue, analogy, visual style cue, or path
  contract that may help but cannot support a claim by itself.
- `suppress`: context that should not enter normal retrieval.
- `quarantine_candidate`: context that may need cooling until verified.

## Evidence Rules

Allowed evidence classes are `real_observed`, `human_approved`, `derived`, and
`synthetic`. Synthetic or external-model material must never surface as proof.
It can only appear as warning, critique, source lead, or proposal context with
an exact waking check.

Tracked registry records must contain sanitized pattern names and source
pointers only. Dynamic outcomes and private feedback stay ignored under
`agent_context/local/context_resonance/`.

## Memory Type Taxonomy

Use these types when deciding where a lesson belongs:

- evidence: verified local files, commands, plots, logs, user statements, or
  validator output;
- experience: observed friction, repeated correction, workflow drag, or a
  failed attempt;
- belief: provisional interpretation that still needs waking verification;
- style: Justin-facing visual hierarchy, wording, speaker-script voice, or
  formatting preference;
- policy: durable behavior rule, hard stop, or safety condition;
- skill/procedure: repeatable method, checklist, script default, or validator.

Tracked memory registries should carry compact route cues and source pointers,
not raw reflections. Detailed reflection records should use the templates in
`agent_context/templates/` and should cite concrete evidence.

## Correction Categories

When feedback is promoted into memory, assign one primary category:

- style/narrative;
- physics/science;
- provenance/evidence;
- workflow/safety;
- retrieval failure;
- memory salience failure;
- task/status failure;
- tool/runtime failure;
- artifact-quality failure.

Reusable slide/plot/script failures should normally become a style-map entry,
negative-memory trap, schema, or validator check. One-off local preferences
should stay in the artifact manifest or chat handoff.

## Plasticity

Each record can be `stable`, `plastic`, `cooling`, or `quarantine`. Stable
records are durable route hints. Plastic records are still being tuned.
Cooling records should surface less often. Quarantine records should not appear
as ordinary context.

## Feedback Loop

After a memory materially affects a task, record whether it was `useful`,
`stale`, `harmful`, `irrelevant`, or `missed`:

```bash
python3 scripts/os/context/codex_context_resonance.py record-outcome \
  --memory-id "<id>" \
  --task "<task>" \
  --result useful \
  --note "<short evidence-backed note>"
```

The feedback ledger is local-only routing evidence, not canonical project
truth. The resolver uses it to identify memories that should be promoted,
cooled, compressed, or turned into validators.

Feedback training is intentionally narrow: only durable memory, schema, and
negative-memory records are trainable. Route policies, live workstream or
artifact hints, suppressed context, synthetic material, dream output, and
external-model material remain fixed pointers; strong promotion, cooling,
quarantine, or suppression proposals require repeated evidence-backed support.
