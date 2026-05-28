# Ask ChatGPT Delegation

This policy defines when Codex may use the user's authenticated ChatGPT UI in
Google Chrome as a delegated research and critique instrument for
ThesisAnalysis.

## Core Rule

ChatGPT UI is an external research worker, not a source of truth, user
approval, or project memory.

Codex remains the integrator. Any ChatGPT output must be treated as:

- source leads;
- design hypotheses;
- critique prompts;
- evaluation ideas;
- alternate formulations to verify.

It is not evidence until Codex verifies it against primary sources, local
files, logs, plots, or Justin's real instructions.

## Trigger

Justin may explicitly invoke this lane with:

- "ask ChatGPT";
- "use ChatGPT";
- "deep research in ChatGPT";
- close variants that clearly mean delegated ChatGPT UI research.

Codex may also choose this lane without asking again when a task would
materially benefit from external long-form synthesis, literature scoping,
agentic-system critique, or broad research, provided the prompt is sanitized
and no sensitive data is transmitted.

## Allowed Automatic Use

When this policy is loaded and the prompt is sanitized, Codex may use Chrome or
Computer Use to:

- open or claim ChatGPT UI;
- start a new chat or continue a clearly relevant ChatGPT research chat;
- type a sanitized research prompt;
- let a long response or deep-research task run;
- collect the response;
- ask bounded follow-up questions;
- save the prompt, transcript, and Codex synthesis under local private state.

No additional confirmation is needed for those steps if the content is
non-sensitive and project-safe.

## Hard Stops

Stop and ask Justin before:

- uploading files to ChatGPT;
- pasting secrets, credentials, tokens, passwords, or private keys;
- pasting raw SDCC usernames, host chains, queue snapshots, job IDs, paths, or
  logs that identify live infrastructure;
- pasting unpublished collaborator-sensitive science details, private emails,
  personal data, or raw meeting notes;
- clicking payment, subscription, upgrade, quota-extension, account, or
  permission-changing UI;
- treating a ChatGPT answer as permission to mutate SDCC, Gmail, Drive,
  Slides, Linear, or the repo.

If ChatGPT UI indicates usage limits, payments, credits, or an upgrade gate,
stop and report that state. Do not claim that a ChatGPT UI action is free or
creditless unless the UI itself explicitly proves it during that session.

## Sanitization Contract

Prompts should describe the problem class without leaking private state. Use:

- "a scientific analysis repository" instead of private repo internals when
  possible;
- abstracted job/risk examples instead of exact SDCC cluster IDs or paths;
- summarized workflow patterns instead of raw Gmail, Drive, or Slack content;
- public or already-approved concepts for science context.

Include the safety objective directly in the prompt:

```text
Do not assume access to private project facts. Mark uncertain claims. Prefer
primary sources, named methods, and testable design recommendations.
```

## Provenance

Every delegated ChatGPT research session should be recorded under:

```text
agent_context/local/chatgpt_research/
```

Each run should include:

- `prompt.md`;
- `chatgpt_response.md`;
- `followups.md` if used;
- `codex_synthesis.md`;
- `source_leads.md`;
- `proposed_changes.md`.

These files are private local state and must not be committed.

## Verification Loop

1. Define the exact research question and thesis-OS decision it informs.
2. Sanitize the prompt.
3. Ask ChatGPT UI.
4. Extract claims, source leads, design ideas, and warnings.
5. Ask at least one follow-up when the first response is broad, uncited, or
   missing failure modes.
6. Verify important claims before implementation.
7. Convert useful output into a small policy, script, task, or dream-system
   proposal.
8. Record provenance and limitations.

## Dream-System Use

For dream-state work, ChatGPT UI is useful for:

- neuroscience source leads on sleep, memory consolidation, replay,
  abstraction, threat simulation, and emotional salience;
- ML analogs such as experience replay, self-play, world models, curriculum
  generation, adversarial evaluation, and reflection loops;
- operations analogs such as chaos engineering, incident review, monitoring,
  postmortems, and progressive delivery;
- critique of the dream pipeline's safety boundaries.

ChatGPT may suggest dream prompts and selection pressure. It must not write
real task state, impersonate Justin, or decide real scientific priorities.

## Output Standard

Codex's final synthesis should state:

- what ChatGPT was asked in sanitized form;
- which claims were useful;
- which claims were unverified or rejected;
- what local OS changes were made or proposed;
- what evidence still needs primary-source verification.
