# Tools And Dependencies

## Tool Discovery

- Prefer existing project scripts and local helpers before inventing new
  one-off commands.
- For text/file search, use `rg` or `rg --files` first.
- Prefer small, standard developer tools when useful, for example `ripgrep`
  for fast code search.
- For Google Docs/Slides/Sheets, Drive, Gmail, Calendar, browser, or desktop
  work, use the matching connector/plugin/skill workflow when available.
- For Office/PDF/spreadsheet/deck local artifact work, call the workspace
  dependency locator before assuming Python/Node package paths.

## Missing Tools

If a tool or binary is missing:

1. check whether the repo already has a wrapper or documented environment;
2. check `scripts/`, `macros/`, `agent_context/`, and `codex_notes/` for local
   precedent;
3. ask before installing software, fetching dependencies, or changing global
   environment state.

If a missing local CLI tool or Python package would materially improve speed,
accuracy, or reliability, ask concisely to install it and name the practical
reason. Do not silently install dependencies; use the approved escalation flow
when installation needs network access or writes outside the workspace.

## Runtime Bias

- ROOT-dependent work uses `scripts/root_in_analysis_env.sh`.
- AuAu ML work on SDCC uses the explicit `RJ_ML_PYTHON` path from
  `SDCC_OPERATIONS.md`.
- Browser/Computer Use is for visual/UI inspection or narrow UI actions, not
  for mutating Drive/Slides when a connector can do it safely.
