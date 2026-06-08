# Delegation Kernel

This policy is the Codex-controlled boundary for external model workers in
ThesisAnalysis. Codex remains the kernel: it chooses the task, builds the
sealed context pack, invokes the worker, verifies the result, and integrates
only evidence-backed changes.

The first implemented worker is terminal Claude Code as a slide-policy worker.
Its purpose is narrow: inspect slide PNGs, generator scripts, speaker/narrative
targets, and slide policy excerpts, then return structured critique and
learning-atom proposals. It is not a second project owner.

## Authority Model

- Codex owns orchestration, final claims, repo policy, Linear, Google
  Slides/Drive, SDCC, Git, and task state.
- Workers produce reports or proposals only. Worker output is never accepted as
  proof until Codex verifies it against local files, rendered artifacts,
  primary sources, or Justin's explicit correction.
- No worker may submit jobs, remove jobs, edit SDCC, mutate Google Slides,
  update Linear, commit code, install dependencies, change credentials, or
  promote policy/memory by itself.
- Do not delegate if Codex can answer or act safely with the local context
  already loaded. Delegate only when the bounded lane can materially improve
  quality, speed, or critique.

## Sealed Context Contract

Every terminal worker invocation should start from a run directory under:

`agent_context/local/delegations/<worker_id>/<run_id>/`

This directory is local/ignored. It should contain:

- `manifest.json`: objective, inputs, policy files, command, and timestamps;
- `context_pack.md`: the exact bounded context Codex intended to expose;
- `prompt.txt`: the worker prompt passed to the model;
- `claude_stdout.txt` / `claude_stderr.txt`: raw invocation logs when run;
- `worker_report.md`: extracted report when the model returns one.

Use `scripts/os/delegation/claude_slide_worker.py` to build and invoke these
packs. Do not rely on a worker's old session state, global memory, or
continuation history. Prefer `--no-session-persistence`, a bounded budget, and
explicit named files.

## Claude Code Slide Worker

Native CLI path:

`/Users/patsfan753/.local/bin/claude`

Project-scoped agent file:

`.claude/agents/thesis-slide-policy-worker.md`

Default safe invocation shape:

```bash
python3 scripts/os/delegation/claude_slide_worker.py pack \
  --objective "Critique this slide candidate for visual clarity and narrative delivery" \
  --png /absolute/path/to/slide.png \
  --script scripts/slides/path/to/generator.py

python3 scripts/os/delegation/claude_slide_worker.py invoke \
  --pack-dir agent_context/local/delegations/claude_slide_worker/<run_id> \
  --max-budget-usd 0.50
```

The harness invokes Claude with a project agent, plan permission mode,
`--no-session-persistence`, a small budget, and explicit read-only tool
allowances. If Claude asks for broader access, authentication changes,
installation, web access, or a mutation, stop and bring the request back to
Codex/Justin.

Authentication is a user boundary. Codex may check:

```bash
/Users/patsfan753/.local/bin/claude auth status --json
```

If it reports `loggedIn: false`, Justin should complete:

```bash
/Users/patsfan753/.local/bin/claude auth login --claudeai
```

Codex must not ask for, store, or type Anthropic credentials or subscription
tokens. After Justin logs in, rerun `doctor` and the smoke invocation.

## Claude Desktop Bridge

`CLAUDE.md` is local-only guidance for Claude Desktop/Claude Code sessions in
this checkout. It is intentionally ignored by GitHub. Keep it aligned with this
policy so a desktop Claude agent enters the same compartment as the terminal
Claude slide worker:

- read `AGENTS.md`, this policy, and `COLLABORATION_AND_HANDOFFS.md`;
- treat Codex as the kernel and Claude as a bounded worker;
- for slide-policy tasks, use the same sealed context-pack directory under
  `agent_context/local/delegations/claude_slide_worker/`;
- ignore prior Claude Desktop memory or old chat state unless Codex explicitly
  includes it in the context pack;
- return the `Slide Worker Report` contract and never mutate SDCC, Slides,
  Linear, Git, repo policy, memory, or final physics state.

If a desktop Claude surface cannot use the terminal harness directly, Codex can
still create the pack locally and pass the `context_pack.md` content/path to
Claude Desktop. The output should be saved or pasted back as `worker_report.md`
for Codex verification.

## ChatGPT And Other Model Lanes

ChatGPT remains the default paid online research/critique lane when the task
needs web context, idea synthesis, or long-form outside critique. Follow
`ASK_CHATGPT_DELEGATION.md`; staged prompts and transcripts stay under
`agent_context/local/chatgpt_research/`.

Claude Code is currently the terminal/codebase-contained slide worker. It is
valuable because it can cheaply inspect local project files in a narrow context
and produce a second-pass slide critique without consuming Codex context on the
whole repository.

Future workers should be added only after the first worker has evidence:

- at least one real slide task where the worker caught a useful issue;
- one documented false-positive or waste case;
- measured iteration cost, latency, and Codex integration burden;
- a clear narrower role than Codex itself.

## Promotion Rule

Worker reports can become durable OS knowledge only through a waking Codex
integration step:

1. Codex verifies the claim against artifacts or Justin feedback.
2. Codex decides the smallest durable surface: policy, style map, template,
   validator, register note, or Linear follow-up.
3. Codex records evidence and why chat memory alone is insufficient.
4. Codex keeps the original worker run path as provenance.

If those steps are not satisfied, keep the worker output as local scratch
context only.
