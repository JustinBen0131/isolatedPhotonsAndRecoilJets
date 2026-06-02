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

Current subscription/default model policy: Justin is paying for ChatGPT Pro as
the primary online model lane. Treat ChatGPT as the go-to delegated model for
online research and critique. Do not route work to Claude, Gemini, or other
paid model subscriptions by default. Use other models only when they are already
available without extra subscription friction and there is a clear task-specific
reason they are better than ChatGPT for that interface or deliverable. If a
future subscription upgrade would materially improve the ThesisAnalysis OS,
state the expected benefit and tradeoff explicitly before recommending it.

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

Use one clean research objective per ChatGPT context. A clean context does not
always mean a brand-new conversation: reuse an existing thread only when the
prior context is an asset instead of a liability.

- The first message for that objective must be one consolidated prompt.
- The first message must be staged locally in a `.txt` file before any UI
  action. Do not compose, retype, or stream the first prompt directly into the
  ChatGPT web input.
- The UI send path is: generate/verify `single_message_prompt.txt`, generate a
  newline-safe `clipboard_prompt.txt`, copy that file to the system clipboard,
  open a fresh ChatGPT chat, paste once, visually confirm the full prompt is in
  the box, then click/send once.
- Do not split one logical prompt across multiple partial sends.
- Do not "warm up" ChatGPT by dripping context fragments before the real ask.
- Do not send a second fragment because the first fragment was incomplete.
- If an accidental partial send happens, treat that chat as contaminated for
  that objective. Start a fresh chat and resend one complete first prompt.
- If stale partial text remains in a new-chat input box, clear it and verify the
  input is empty before pasting from the staged prompt file.

This policy exists because fragmented prompt entry degrades answer quality,
pollutes the thread state, and makes provenance harder to trust the next
morning.

Before pressing Enter on the first message:

1. finish the entire prompt locally first;
2. save it to `single_message_prompt.txt`;
3. create `clipboard_prompt.txt` from the same content with hard line breaks
   collapsed unless the line breaks are semantically required;
4. verify both files are sanitized and under the length budget;
5. copy `clipboard_prompt.txt` to the clipboard;
6. paste once into a fresh ChatGPT chat;
7. visually confirm the prompt text was pasted as one message;
8. only then submit it once.

### Thread Reuse And Continuation

Before using ChatGPT UI, decide whether the right context is `fresh`,
`continue`, or `fork_with_recap`.

Use `continue` when all are true:

- the thread is about the same durable object: one paper, deck, slide family,
  campaign, code artifact, incident, or policy;
- the prior attachments, wording conventions, evidence frame, or critique
  history materially reduce context cost;
- the thread has not drifted into unrelated objectives;
- there is no stale partial prompt, accidental send, or contaminated setup;
- the current task can be stated as one complete continuation message.

Use `fresh` when any are true:

- the objective, audience, artifact, campaign, or risk class changed;
- the old thread mixed too many unrelated tasks;
- the old context contains private or sensitive details not needed now;
- the old answer has stale assumptions that would take more work to correct
  than to restate cleanly;
- a Pro/Deep escalation would benefit from a clean first prompt.

Use `fork_with_recap` when the old thread has valuable context but too much
drift. In that case, create a new chat and start with a compact context capsule
summarizing only what should survive.

For continued threads, the first Codex-authored message after returning must
still be staged locally and pasted once. Use this capsule:

```text
Continuation objective:
Relevant prior context to preserve:
What changed since the prior answer:
Current evidence or artifact frame:
Hard constraints and exclusions:
Do not assume:
Deliverable:
Failure modes to check:
```

Observed useful ChatGPT thread archetypes:

- deck/script threads: continue when the same deck PDF, screenshots, slide
  sequence, and speaker-note tone are the useful context;
- manuscript/paper project threads: reuse the project container, but keep each
  child chat scoped to one stage such as rejection response, surgical LaTeX
  revision, or broad paper strategy;
- strategic science architecture threads: use a fresh high-effort prompt with
  explicit stages, constraints, and requested output columns;
- surgical review/editing threads: use a role stack, exact allowed scope,
  hard exclusions, and patch-like output format;
- SDCC/Condor/debug threads: continue only for the same incident/campaign and
  include a fresh current-evidence plus hard-stop capsule before trusting prior
  commands or paths;
- Codex-OS/dream research threads: continue while the design objective is the
  same; start fresh or fork with recap when moving from concept research to
  implementation policy.

## Mode Routing

Public OpenAI guidance and practical UI behavior both favor starting with the
smallest mode that preserves the task contract and using more reasoning only
when quality gains justify extra time and compute. ChatGPT UI labels such as
`instant`, `thinking`, `pro`, and their submodes are treated here as operating
modes, not guaranteed public API contracts.

Codex should route by family and submode:

1. `instant`
   Use for quick fact gathering, terminology checks, simple rewrites, short
   comparisons, narrow source hunting, and low-risk brainstorming where a fast
   first pass is enough.
2. `thinking-light`
   Use for low-risk synthesis, short critiques, outline alternatives, small
   source-lead requests, and prompt polishing where a plain instant answer may
   be too shallow but deep reasoning would be wasteful.
3. `thinking-standard`
   Default for most delegated work: design critique, policy drafting,
   medium-complexity research, synthesis across a few sources, and structured
   brainstorming.
4. `thinking-extended`
   Use when the prompt has several constraints, needs failure modes, compares
   alternatives, or informs a repo policy/script change.
5. `thinking-heavy`
   Use when the question is multi-constraint, cross-domain, architecture-heavy,
   failure-mode sensitive, or needs a deeper literature and systems synthesis,
   but is still short enough for Codex to wait and collect.
6. `pro-standard`
   Use only when the UI exposes Pro and the task needs premium depth, or when a
   thinking-heavy answer remains materially incomplete after bounded iteration.
7. `pro-extended`
   Use for maximum-depth architecture/research synthesis, unusually long
   multi-domain prompts, source-lead-heavy literature scoping, or high-stakes
   OS-policy changes where standard Pro is likely to be too shallow.
8. `deep-research`
   Use only when Justin explicitly asks for Deep Research or the UI exposes a
   separate long-running research workflow that is clearly better than
   Pro-extended for the objective.

Routing matrix:

| Prompt shape | First mode | Escalate if |
| --- | --- | --- |
| Simple wording, title, short summary, narrow definition | `instant` | answer misses nuance -> `thinking-light` |
| Small critique, compact source hunt, one-page outline | `thinking-light` | lacks structure/source leads -> `thinking-standard` |
| Default policy/design/synthesis with bounded constraints | `thinking-standard` | misses constraints/failure modes -> `thinking-extended` |
| Repo-facing policy/script design, safety-sensitive workflow, several constraints | `thinking-extended` | still shallow or underdetermined -> `thinking-heavy` |
| Cross-domain architecture, dream/OS design, literature + SRE + ML synthesis | `thinking-heavy` | materially incomplete after one follow-up -> `pro-standard` |
| Maximum-depth architecture or long research job | `pro-standard` | needs broader source synthesis or prior Pro is shallow -> `pro-extended` |
| Explicit Deep Research/source-report request | `deep-research` | n/a; hand off to Justin paste-back |

Default routing rules:

- Start with the smallest mode likely to succeed.
- Prefer `thinking` over `instant` when correctness depends on synthesis rather
  than retrieval.
- Within `thinking`, prefer `standard` unless the prompt is obviously simple
  (`light`) or obviously multi-constraint/safety-sensitive
  (`extended`/`heavy`).
- Within `pro`, choose `standard` for deep single-objective synthesis and
  `extended` for maximum-depth, long, cross-domain, source-heavy work.
- Do not default to any Pro mode.
- Escalate only when the previous answer is shallow, misses constraints, lacks
  rigor, or still leaves the design decision underdetermined.
- After the hard question is answered, de-escalate for follow-up polishing or
  extraction work.
- Treat the user's own prompting pattern as signal: Justin often benefits from
  richer conceptual prompts with explicit analogies, hard boundaries, and
  deliverable shape, but Codex should convert that into a compact staged prompt
  rather than many UI fragments.

## Response Collection Policy

Codex should not spend its own working time or token budget polling long
premium ChatGPT research jobs.

- `instant`, `thinking-light`, `thinking-standard`, `thinking-extended`, and
  `thinking-heavy`: Codex may submit the sanitized staged prompt, wait for the
  response, collect it, ask bounded follow-ups, save the transcript pack, and
  synthesize the answer in the same turn.
- `pro-standard`, `pro-extended`, `deep research`, or equivalent long-running
  premium research modes: Codex submits one complete sanitized staged prompt,
  records the mode/thread/local research pack, then stops. Justin pastes the
  completed ChatGPT response back into Codex when it is done.

For a Pro/Deep handoff, Codex should tell Justin:

```text
I sent the sanitized prompt to ChatGPT in <mode>.
Local research pack: <path>
Paste the completed response back here when it finishes, and I will verify,
synthesize, and turn only useful pieces into local proposal artifacts.
```

Do not create a 15-minute polling heartbeat for Pro/Deep by default. Treat it
like a long external research job whose output returns through Justin.

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
- direct UI typing of the first prompt;
- multi-line Computer Use `type_text` entry for ChatGPT first prompts;
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
- let `instant`, `thinking`, or `heavy` responses run and collect them;
- ask bounded follow-up questions for Codex-managed modes;
- for `pro`, `extended pro`, or `deep research`, submit once and hand off to
  Justin for paste-back instead of polling;
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
5. If the mode is Pro/Deep, report the local pack and wait for Justin to paste
   the finished response before continuing.
6. Otherwise extract claims, source leads, design ideas, and warnings.
7. Ask at least one bounded follow-up when the first response is broad, uncited, or
   missing failure modes.
8. Escalate mode only if the prior answer is materially insufficient for the
   decision at hand.
9. Verify important claims before implementation.
10. Convert useful output into a small policy, script, task, or dream-system
   proposal.
11. Record provenance and limitations.

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
