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
6. Retire: remove stale daily/active state so the OS does not accumulate fear
   or clutter.

This is the Darwinian loop for the agent OS: the selected trait is not survival
of the agent, but completion of the thesis goal with less duplicated work,
lower risk, and higher scientific confidence.

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

- `python3 scripts/codex_context_pack.py --cadence boot`: compact boot packet
  for a fresh Codex chat.
- `python3 scripts/codex_context_pack.py --cadence daily`: daily cockpit from
  register state.
- `python3 scripts/codex_context_pack.py --cadence weekly`: weekly review for
  promotion, demotion, closure, and friction-to-policy conversion.
- `python3 scripts/codex_context_pack.py --cadence monthly`: thesis claim-layer
  audit.
- `python3 scripts/codex_thesis_radar.py`: live-work mapping to the thesis
  spine.
- `python3 scripts/codex_artifact_registry.py check`: provenance check for
  canonical artifacts.

`agent_context/ARTIFACT_REGISTRY.yaml` is the canonical registry for
presentation- or thesis-facing artifacts. Add plots, ROOT files, CSVs, JSONs,
slide candidates, and scripts before treating them as reusable evidence.

## Dream Layer

Use `agent_context/policies/AGENTIC_OS_DREAMING.md` for private synthetic
rehearsal. Dreams are proposal-only and may write only under
`agent_context/local/dreams/`.

```bash
python3 scripts/codex_os_dream.py micro
python3 scripts/codex_os_dream.py nightly
python3 scripts/codex_os_dream.py validate --latest
```

Dreams simulate tomorrow's likely questions and failure modes. They must never
treat synthetic Justin as approval, never touch external systems, and never
mutate real repo state. Waking Codex may promote a dream proposal only through
normal task-capture, guard, doctor, and evidence checks.

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
