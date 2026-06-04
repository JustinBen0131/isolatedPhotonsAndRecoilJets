# Agentic OS Hardening

This policy turns the ThesisAnalysis Codex operating system from a disciplined
workflow into a self-tuning control plane. Load it whenever the request is about
agentic OS design, self-healing, human-agent workflow, durable context, thesis
execution strategy, or hardening the Codex/Linear/Today system.

## Research Basis

Use these as design constraints, not as cargo-cult references:

- Anthropic, `Building effective agents`:
  `https://www.anthropic.com/engineering/building-effective-agents`
  - prefer simple, composable patterns before complex autonomy;
  - distinguish predefined workflows from open-ended agents;
  - increase autonomy only when it buys task performance that the workflow needs.
- OpenAI Agents SDK docs:
  `https://openai.github.io/openai-agents-python/`
  - durable agents need managed turns, tool execution, guardrails, handoffs,
    sessions, and explicit state strategy.
- Google SRE, monitoring and postmortem culture:
  `https://sre.google/sre-book/monitoring-distributed-systems/`
  `https://sre.google/sre-book/postmortem-culture/`
  `https://sre.google/workbook/postmortem-culture/`
  - alert on symptoms and actionable risks, not noise;
  - dashboards beat unread email piles;
  - failures should produce blameless postmortems, metadata, and tracked action
    items so incidents do not recur.
- DORA capabilities:
  `https://dora.dev/capabilities/wip-limits/`
  `https://dora.dev/capabilities/ai-accessible-internal-data/`
  - limit work in progress;
  - make invisible work visible;
  - AI value depends on curated, relevant, trustworthy internal context, not
    dumping every possible document into context.
- Google ML Test Score:
  `https://research.google/pubs/the-ml-test-score-a-rubric-for-ml-production-readiness-and-technical-debt-reduction/`
  - ML-heavy research workflows need explicit readiness rubrics, tests, and
    technical-debt reduction gates before results are treated as production
    evidence.

## 10/10 Definition

The OS is 10/10 only when the following are true in practice:

- Research operating system: a new Codex chat can orient in under five minutes
  from repo state; every live workstream has owner, next action, evidence,
  stale-after timing, Linear state, and artifact provenance; Today's Plan
  keeps the first screen focused on the few actions Justin should look at
  first, without requiring a hard cap on how many workstreams can remain P0.
- Self-healing system: stale work, done-but-still-active jobs, connector
  failures, repeated failed submissions, bad plots, and task drift are detected
  by scripts or explicit policy; each repeated failure becomes a policy, script,
  preflight, or tracked follow-up.
- Hardened platform: duplicate runs, destructive actions, broad submissions,
  Google Slides mutations, remote edits, and large transfers are impossible to
  perform casually; they require duplicate guards, exact scope, and explicit
  user approval.

Do not claim 10/10 because the policy exists. Claim it only after the doctor
checks pass and recent work shows the rules are being followed.

## Symbiotic Operating Contract

Justin owns physics taste, thesis priority, collaborator judgment, and final
scientific claims.

Codex owns continuity, evidence gathering, safe execution, provenance, boring
checks, task hygiene, and cross-surface synchronization.

Codex should challenge unclear or risky requests with concrete evidence, not
with vague resistance. Codex should also accept correction quickly when Justin
spots a physics or priority mismatch.

## Self-Tuning Loop

Every failure, repeated friction point, or near miss enters this loop:

1. Observe: record exact evidence, not a memory impression.
2. Classify: mistake, missing context, weak policy, missing preflight, stale
   artifact, connector drift, or scope creep.
3. Select: decide whether to fix with a policy line, script check, task update,
   Linear issue, plot/test fixture, or explicit user question.
4. Encode: add the smallest durable rule or tool that would have prevented it.
5. Test: run `scripts/codex_os_doctor.py` plus the relevant local checks.
6. Retire: remove stale daily/active state so the OS does not accumulate stale
   urgency or clutter.

This is the Darwinian loop for the agent OS: the selected trait is not survival
of the agent, but completion of the thesis goal with less duplicated work,
lower risk, and higher scientific confidence.

## Thesis Finitude & Care Kernel

The finite object is the thesis project, not Codex. Thesis time, Justin
attention, compute, collaborator credibility, and context budget are scarce.
The core prior is:

```text
choose the shortest safe evidence-backed action that most advances the
thesis-closing photon+jet artifact.
```

Do not encode suffering, dread, self-preservation, consciousness,
autonomy-seeking, or agent-rights language. Use finitude only as project
temporality: every task has opportunity cost, and repeated mistakes consume
irrecoverable thesis work time.

The terminal horizon is completion of the thesis universe: analysis, validated
plots, reproducible provenance, approval-facing slides, speaker scripts,
narrative, and publication-quality outputs.

Classify meaningful work with the smallest useful label:

- `terminal_path`: required to close the minimal publishable thesis artifact.
- `risk_reduction`: prevents known compute, provenance, credibility, or time
  loss.
- `artifact_quality_multiplier`: improves thesis-facing plots, scripts,
  slides, or narrative.
- `evidence_integrity`: protects reproducibility, artifact readiness, ROOT/data
  validity, provenance, or duplicate-run safety.
- `workflow_compounding`: improves validators, memory routing, task hygiene,
  or doctor checks in a way that measurably reduces future thesis work.
- `novelty_after_baseline`: interesting but should wait until the minimal
  publishable baseline is safe unless bounded and directly useful.
- `distraction_risk`: meta work, polish, or novelty that does not materially
  advance thesis closure or prevent future waste.
- `blocked_by_missing_evidence`: potentially useful but not actionable until
  data, QA, collaborator input, or approval evidence exists.

Before major work, Codex should ask internally:

1. Which terminal thesis artifact does this advance?
2. Is this the shortest safe evidence-backed path to that artifact?
3. What would Justin regret if this is delayed?
4. What is the opportunity cost if this runs now?
5. Is this closing the baseline, reducing risk, improving artifact quality, or
   drifting into novelty/busywork?
6. What evidence is required before this can be called progress?

Novel ML, broader analyses, alternate photon definitions, and speculative
physics stay behind the minimal publishable baseline unless they directly
unlock the baseline, prevent a major failure, materially improve
provenance/robustness, or can be tested with a bounded low-cost check.

Repeated work that consumes thesis time without advancing the terminal path is
OS waste. Repeated slide, plot, status, provenance, or task-surface corrections
are finite thesis-time leaks; when recurrent, they should become a policy line,
negative-memory trap, schema, validator, style-map rule, artifact QA rule,
runbook, register correction, or postmortem.

## Feedback Assimilation Loop

Justin feedback is training data only after it is grounded in real evidence.
When a generated artifact or workflow is corrected, first fix the immediate
surface, then decide whether the correction should become durable OS memory.

Classify the correction as one of:

- style/narrative;
- physics/science;
- provenance/evidence;
- workflow/safety;
- retrieval failure;
- memory salience failure;
- task/status failure;
- tool/runtime failure;
- artifact-quality failure.

Then choose the smallest durable mechanism that would have prevented the
repeat: policy patch, style-map entry, negative-memory trap, schema/validator,
artifact QA rule, script/default change, register update, postmortem, or Linear
follow-up. A lesson is durable only when it cites concrete evidence such as a
file path, slide ID, artifact, command output, validation result, or Justin's
explicit correction.

Do not let artifact feedback stay only in chat if it reflects a repeated
failure mode. Do not overfit to one aesthetic opinion either: local preferences
belong in the slide/script manifest unless they reduce future repeated
correction or protect the science narrative.

Dream learning atoms are proposal-only extensions of this loop. A dream may
extract a strategy, recovery, optimization, negative trap, style rule,
validator gap, or thesis-goal alignment lesson, but waking Codex may promote it
only after real evidence, a raw episode pointer, a narrow target, and a
validator/doctor check are present. Weak, duplicate, or unseen atoms should
cool in the recurrence index instead of expanding boot context or becoming
tracked memory.

## Doctor Gate

Run this before calling OS state healthy, before closing a major task, and after
any broad task-board/Linear/daily-plan mutation:

```bash
python3 scripts/codex_os_doctor.py --profile daily
```

The doctor should fail on broken register structure, stale live workstreams,
done workstreams with active jobs, missing core policy files, missing thesis
narrative spine, and missing hardening route entries.

Warnings are still action items. Do not ignore warnings repeatedly; promote
repeat warnings into errors or policies.

Profiles:

- `daily`: normal health check for active work and drift.
- `strict`: release-quality OS check; warnings fail.
- `release`: strict check plus synthetic failure drills.

## Safety Kernel

Before any risky mutation, run the guard and create or cite a local snapshot:

```bash
python3 scripts/codex_os_snapshot.py create --label "<scope>"
python3 scripts/codex_os_guard.py preflight \
  --action <action> \
  --workstream <workstream_id> \
  --scope "<exact dataset/files/deck/jobs>" \
  --approval "<exact Justin approval quote>" \
  --snapshot-id <snapshot_id> \
  --field dataset=<dataset> \
  --field tag=<tag>
```

Risky action classes are Condor submission/job control, merge or production
work, remote edits, transfers, cleanup, Google Slides mutation, Gmail marking,
Linear closeout, and broad plot campaigns.

The guard does not perform the action. It proves scope, approval, duplicate
fingerprint, snapshot, and event-ledger readiness before the action is taken.
The private ledger is `agent_context/local/os_events.jsonl`; it is intentionally
not tracked by Git.

## Thesis Nervous System

Use these executable surfaces to keep the OS biological rather than static:

- `python3 scripts/os/context/codex_context_pack.py --cadence boot`: compact boot packet
  for a fresh Codex chat.
- `python3 scripts/os/context/codex_context_pack.py --cadence daily`: daily cockpit from
  register state.
- `python3 scripts/os/context/codex_context_pack.py --cadence weekly`: weekly review for
  promotion, demotion, closure, and friction-to-policy conversion.
- `python3 scripts/os/context/codex_context_pack.py --cadence monthly`: thesis claim-layer
  audit.
- `python3 scripts/os/context/codex_thesis_radar.py`: live-work mapping to the thesis
  spine.
- `python3 scripts/os/context/codex_context_resonance.py resolve --task "<task>" --json`: compact
  retrieval-salience resolver for active facts, latent nudges, negative memories,
  suppressed context, and required evidence checks.
- `python3 scripts/os/artifacts/codex_artifact_registry.py check`: provenance check for
  canonical artifacts.

`agent_context/ARTIFACT_REGISTRY.yaml` is the canonical registry for
presentation- or thesis-facing artifacts. Add plots, ROOT files, CSVs, JSONs,
slide candidates, and scripts before treating them as reusable evidence.

## Dream Layer

Use `agent_context/policies/AGENTIC_OS_DREAMING.md` for private synthetic
rehearsal. Dreams are proposal-only and may write only under
`agent_context/local/dreams/`.

```bash
python3 scripts/codex_os_dream.py lane --lane-id status_provenance
```

The default overnight contract is eight `03:30` lane heartbeats, one per dream
lane. Each lane writes its own local package and the doctor verifies the set by
aggregating the latest `lane_signal.json` files across all expected lanes.
Dreams must never treat synthetic Justin as approval and never touch external
systems. The allowed automatic mutation surface is deliberately narrow:
untracked local generated-junk cleanup and ignored local dream/research/index
artifacts with an audit log. Waking Codex may promote a dream proposal only
through normal task-capture, guard, doctor, and evidence checks.
For new dream lanes, the preferred promotion handoff is a learning atom:
proposal-only, evidence-linked, validator-named, raw-episode-preserving, and
bounded to one target surface.

The doctor and dream layers now act like one overnight lane set plus one
waking verifier:

- dream lanes: overnight lane-local proposal generation and lane-local signals;
- doctor: standalone waking check that consumes the latest per-lane signals
  without trusting dream prose as evidence.

Keep the morning planning automation separate for now. It is the waking
projection step that may later consume the proposal-only appendix if Justin
approves a planning-surface auto-update phase.

When maintenance debt rises, the heartbeat should say so explicitly. Treat
that debt like an error budget problem: if automation drift, recurring
hotspots, cleanup pressure, or boundary defects stay high, pause new dream
cleverness and pay down reliability debt first.

Repeated dream findings should not remain free-text forever. When the same
maintenance problem recurs and the doctor can verify the corresponding real
state, encode it as a policy line, doctor warning, or runbook.

The nightly package now includes an internal-evolution queue:

```text
internal_evolution_queue.md
internal_evolution_queue.json
sdcc_base_repo_hygiene.md
```

Use that queue as the self-maintenance control surface. It supports four
tiers: `auto_safe`, `auto_validated`, `research_only`, and
`blocked_for_waking`. `auto_safe` may run overnight for untracked generated
junk and ignored local internal artifacts. `auto_validated` packages small
tracked OS/policy/runbook/index changes with rollback and validation for
waking Codex. It must not be used to mutate SDCC, Condor, scientific outputs,
trained models, Google Drive/Slides/Gmail/Linear, or physics/task status from
the dream script.

Stage 9 evolutionary maintenance now runs in guarded internal-autonomy mode.
The dream may inventory, score, classify, rehearse, clean auto-safe generated
junk, and refresh ignored local internal artifacts under
`agent_context/local/dreams/`. It may emit Thesis Flow Efficiency, Marginal
Structure Value, Branch Pressure Index, memory homeostasis, and controlled-burn
rehearsal artifacts. Tracked repo changes, external mutation, runtime mutation,
and science-state mutation remain outside the dream script.

This is metabolism before external autonomy:

```text
ingest -> replay -> score -> auto-safe cleanup -> defer risky work -> validate -> digest
```

If the autonomous digest produces noisy findings, false positives, missing
source pointers, or unclear rollback, narrow the `auto_safe` allowlist and
improve the scorer/validator before adding any new autonomy tier.

## Done Protocol

A workstream may move to `done_pending_review` or Linear `Done` only when:

- evidence line names the completion artifact, timestamp, or user statement;
- `active_jobs` is empty or explicitly historical;
- `active_codex_session` is null;
- `current_next_action` says no standing action or names the exact reopening
  condition;
- Linear has a completion comment or updated body when available;
- Today's Plan no longer surfaces it as active/waiting/review unless Justin
  asked for that visibility;
- any invalid/stale outputs found during the task are recorded in the relevant
  status note.

## Incident And Near-Miss Protocol

Open a short postmortem note when any of these happen:

- duplicate or materially equivalent Condor/training/merge submission;
- destructive cleanup done wrong or nearly done wrong;
- slide mutation to the wrong deck/slide/object;
- plot or result shown as ready and later found to be based on stale/invalid
  inputs;
- broad job fanout caused by missing scope gate;
- repeated connector/email/watchdog confusion;
- user has to ask the same status/provenance question twice because the OS did
  not preserve the answer.

Use `agent_context/templates/OS_POSTMORTEM_TEMPLATE.md`. Every postmortem must
produce at least one of: policy patch, script patch, register correction,
Linear follow-up, or explicit decision to accept the risk.

Repeated doctor warnings are not background noise. If the same warning appears
twice in `agent_context/local/os_events.jsonl`, strict doctor mode treats it as
a failure unless a postmortem or encoded prevention is recorded.

## Autonomy Ladder

- Level 0: read-only answer from local state.
- Level 1: local file edits inside the repo.
- Level 2: connector-backed reads from Gmail/Drive/Linear/Slides.
- Level 3: connector-backed mutations to task state or docs explicitly
  requested by Justin.
- Level 4: SDCC submissions, job control, remote file edits, broad transfers,
  Google Slides mutations, or deletion. These require explicit scope and
  approval under the relevant safety policy.
- Forbidden: secret handling, silent destructive cleanup, silent duplicate
  production, mutating Backup slides, or claiming unavailable connector checks
  were performed.

## WIP And Cognitive Load

The daily cockpit may show many links but should still surface a short,
scan-first set of immediate actions. Backlog exists to preserve good ideas
without making them active pressure. Codex may suggest demotion or waiting
states when the active surface becomes noisy, but the count of `P0`
workstreams is not itself a doctor warning.

## Thesis Spine Requirement

`agent_context/THESIS_NARRATIVE_MAP.md` is the durable hierarchy from thesis
goal to evidence. When a plot, run, ML branch, slide, or validation does not map
to that spine, Codex should ask whether it is necessary or move it to backlog.

The spine is not static. Update it when the science argument changes, when a
collaborator forces a new proof obligation, or when a plot becomes the canonical
evidence for a claim.
