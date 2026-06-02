# Tools And Dependencies

## Tool Discovery

- Prefer existing project scripts and local helpers before inventing new
  one-off commands.
- For text/file search, use `rg` or `rg --files` first.
- Prefer small, standard developer tools when useful, for example `ripgrep`
  for fast code search.
- For Google Docs/Slides/Sheets, Drive, Gmail, Calendar, browser, or desktop
  work, use the matching connector/plugin/skill workflow when available.
- For delegated ChatGPT UI research, load
  `agent_context/policies/ASK_CHATGPT_DELEGATION.md` first and use Computer Use
  in the user's Chrome profile by default, with a sanitized prompt and local
  transcript provenance. Generate a research pack first, copy
  `clipboard_prompt.txt` to the clipboard, paste once, and send once. Before UI
  use, record whether the right context is a fresh chat, a continued thread, or
  a fork with recap. Never type or stream the first prompt directly into the
  ChatGPT UI. Codex may collect `instant` and `thinking_*` responses itself; for
  `pro_standard`, `pro_extended`, or `deep_research`, submit one sanitized
  staged prompt and ask Justin to paste the completed response back into Codex
  instead of polling.
- For Justin's explicit "set up a Zoom room" / "get me a Zoom invite link"
  variants, load `agent_context/policies/MEETING_CAPTURE_AND_ZOOM.md` and use
  Computer Use on Zoom for the narrow room-setup workflow: create/start the
  room, enable or verify available transcript/AI-summary capture, copy the
  invite link, and return it to Justin.
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

## Connector Session Auth

- For Gmail/Drive/Slides/Calendar connector auth errors in a long-running chat,
  do not assume the plugin is broken or needs to be revamped.
- First treat it as possible chat-session credential expiration. If applicable,
  retry once with exact IDs or a minimal connector call.
- If it still fails, tell Justin to fork/reopen the chat because connector auth
  may refresh cleanly in the forked chat. Ask for manual plugin reconnect only
  after the forked/new chat also fails.
- Never claim a connector-backed status check was performed when the connector
  was unavailable.

## Runtime Bias

- ROOT-dependent work uses `scripts/root_in_analysis_env.sh`.
- AuAu ML work on SDCC uses the explicit `RJ_ML_PYTHON` path from
  `SDCC_OPERATIONS.md`.
- For PDF text extraction, prefer the installed user-local command:
  `/Users/patsfan753/.local/bin/pdftotext`. It is Xpdf command-line tools
  4.06 for macOS ARM, installed under
  `/Users/patsfan753/.local/opt/xpdf-tools-mac-4.06/`. Use the full path
  because `~/.local/bin` may not be on `PATH` in Codex shells.
- Browser/Computer Use is for visual/UI inspection, narrow UI actions,
  approved delegated ChatGPT research, and explicit Zoom-room setup requests;
  it is not for mutating Drive/Slides when a connector can do it safely.
