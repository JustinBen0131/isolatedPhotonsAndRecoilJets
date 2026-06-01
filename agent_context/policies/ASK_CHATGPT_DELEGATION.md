# Ask ChatGPT Delegation

This policy defines when Codex may use the user's authenticated ChatGPT UI in
Google Chrome as a delegated research and critique instrument for
ThesisAnalysis. ChatGPT UI work must use Computer Use in the user's visible
Chrome profile by default, matching Justin's normal browser workflow.

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

Codex should treat ChatGPT as a standing research, discussion, and
brainstorming engine whenever that would materially improve a task, and always
when Justin explicitly directs ChatGPT use.

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

## Session Architecture

Use one clean research objective per chat session.

- The first message for that objective must be one consolidated prompt.
- Do not split one logical prompt across multiple partial sends.
- Do not "warm up" ChatGPT by dripping context fragments before the real ask.
- Do not send a second fragment because the first fragment was incomplete.
- If an accidental partial send happens, treat that chat as contaminated for
  that objective. Start a fresh chat and resend one complete first prompt.

This policy exists because fragmented prompt entry degrades answer quality,
pollutes the thread state, and makes provenance harder to trust the next
morning.

Before pressing Enter on the first message:

1. finish the entire prompt locally first;
2. verify it is the intended single-message prompt;
3. verify it is sanitized;
4. only then submit it once.

## Mode Routing

Public OpenAI guidance currently recommends starting with the smallest prompt
that preserves the task contract and using more reasoning only when quality
gains justify the extra time and cost. ChatGPT UI labels such as `instant`,
`thinking`, `heavy`, and `pro` are treated here as operating modes, not as
guaranteed public API contracts.

Codex should use this escalation ladder:

1. `instant`
   Use for quick fact gathering, terminology checks, simple rewrites, short
   comparisons, narrow source hunting, and low-risk brainstorming where a fast
   first pass is enough.
2. `thinking`
   Default for most useful delegated work: design critique, policy drafting,
   medium-complexity research, synthesis across a few sources, and structured
   brainstorming.
3. `heavy`
   Use when the question is multi-constraint, cross-domain, architecture-heavy,
   failure-mode sensitive, or needs a deeper literature and systems synthesis.
4. `pro`
   Use only when the UI explicitly exposes a deeper premium mode and the task
   genuinely needs maximum depth, or when `thinking`/`heavy` passes remain
   materially incomplete after bounded iteration.

Default routing rules:

- Start with the smallest mode likely to succeed.
- Prefer `thinking` over `instant` when correctness depends on synthesis rather
  than retrieval.
- Prefer `heavy` only when the first-pass answer must reason across many
  constraints, anti-patterns, and design tradeoffs.
- Do not default to `pro`.
- Escalate only when the previous answer is shallow, misses constraints, lacks
  rigor, or still leaves the design decision underdetermined.
- After the hard question is answered, de-escalate for follow-up polishing or
  extraction work.

## Prompt Construction Standard

Prompts should look like strong human power-user prompts, not a transcript of
Codex thinking. The first message should be compact, outcome-first, and
complete.

Use this structure:

1. role framing only if it sharpens the response;
2. the exact decision or research objective;
3. concise sanitized context;
4. hard constraints and exclusions;
5. required deliverable shape;
6. uncertainty and citation instructions.

Preferred characteristics:

- one objective per prompt;
- direct language;
- explicit output shape;
- explicit failure-mode request;
- ask for uncertainty marking;
- ask for source leads when research matters;
- no unnecessary chain-of-thought scaffolding;
- no giant dump when a narrower framing would do.

Avoid:

- fragmented setup messages;
- multiple unrelated asks in one first prompt;
- vague "thoughts?" prompts for high-stakes design work;
- exposing private repo, SDCC, Gmail, Drive, or collaborator-sensitive state;
- asking ChatGPT to decide permissions or mutate external systems.

## Follow-Up Discipline

Follow-ups are allowed only after the first full response lands.

- Each follow-up should be a complete single message.
- Each follow-up should have one clear purpose: clarify, deepen, challenge,
  compare, or extract.
- Prefer adversarial follow-ups such as "what would break this design?" or
  "which assumptions are weakest?" over generic "go deeper".
- If the thread drifts into multiple objectives, end it and start a new chat
  with a new consolidated first prompt.

## Allowed Automatic Use

When this policy is loaded and the prompt is sanitized, Codex must use
Computer Use in Google Chrome automatically to:

- open or claim ChatGPT UI;
- start a new chat or continue a clearly relevant ChatGPT research chat;
- type a sanitized research prompt;
- let a long response or deep-research task run;
- collect the response;
- ask bounded follow-up questions;
- save the prompt, transcript, and Codex synthesis under local private state.

No additional confirmation is needed for those steps if the content is
non-sensitive and project-safe.

Do not try the Chrome extension/browser-client path first for ChatGPT UI work.
Use it only when Justin explicitly asks for that surface or when Computer Use is
unavailable and the fallback is explained. Do not use the in-app Browser for
authenticated ChatGPT chats unless Justin explicitly asks for it.

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

For the first prompt in a delegated research session, prefer one compact
single-message ask that fits cleanly into the ChatGPT input box without manual
chunking. Put staged follow-ups in separate later messages only after the first
response lands. Do not spray one long ask across multiple accidental partial
messages.

If the UI, browser, or desktop tool causes a partial send, do not continue the
same objective in that thread. Record the failure locally and restart with a
fresh chat and the validated single-message prompt.

## Provenance

Every delegated ChatGPT research session should be recorded under:

```text
agent_context/local/chatgpt_research/
```

Each run should include:

- `prompt.md`;
- `single_message_prompt.txt`;
- `session_strategy.md`;
- `chatgpt_response.md`;
- `followups.md` if used;
- `codex_synthesis.md`;
- `source_leads.md`;
- `proposed_changes.md`.

These files are private local state and must not be committed.

## Verification Loop

1. Define the exact research question and thesis-OS decision it informs.
2. Sanitize the prompt.
3. Choose the lightest ChatGPT mode likely to succeed and record the escalation
   ladder.
4. Prepare one compact first prompt and send it as a single message.
5. Extract claims, source leads, design ideas, and warnings.
6. Ask at least one bounded follow-up when the first response is broad, uncited, or
   missing failure modes.
7. Escalate mode only if the prior answer is materially insufficient for the
   decision at hand.
8. Verify important claims before implementation.
9. Convert useful output into a small policy, script, task, or dream-system
   proposal.
10. Record provenance and limitations.

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
- which mode was used first and why;
- whether escalation was needed;
- which claims were useful;
- which claims were unverified or rejected;
- what local OS changes were made or proposed;
- what evidence still needs primary-source verification.
