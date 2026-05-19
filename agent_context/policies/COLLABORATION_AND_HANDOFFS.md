# Collaboration And Handoffs

## Ownership Model

- Codex owns orchestration: SDCC status strategy, Gmail pipeline emails, local
  notes, dashboards, Google Slides, Drive presentation context, integration,
  transfer-command guidance, and final sanity checks.
- Claude or any other delegated code/reference worker owns bounded lanes only
  when Justin explicitly asks for delegation: PPG12/reference inspection, local
  code changes, local static analysis, feature implementation, and compact
  handoff.
- Split work only across disjoint lanes by default. Avoid two agents editing
  the same files or solving the same unknown.

## Handoff Shape

Handoffs should be compact and structured:

- changed files;
- evidence/tests;
- physics assumptions;
- risks/open questions;
- exact next command, including SDCC upload/check command when relevant.

Do not delegate just to create motion. Delegate only when the lane is bounded,
non-overlapping, and materially advances the current task.
