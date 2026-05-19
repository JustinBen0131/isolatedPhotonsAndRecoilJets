# Hard Stops And Safety

## Duplicate Runs

If Justin asks for a new training campaign, validation, merge, RecoilJets
production, or production-style analysis that is identical or materially
equivalent to something already present, stop immediately.

Do not create configs, edit drivers, launch jobs, transfer payloads, or spend
compute. Report:

- matching existing tag/path/run name;
- evidence source, such as `STATUS_DASHBOARD.md`, `TASK_BOARD.md`, local
  `dataOutput/`, Gmail, SDCC path, ROOT file, or cluster ID;
- what differs, if anything;
- options: reuse, compare, or deliberately rerun with a new reason.

Detailed procedure: `DUPLICATE_RUN_GUARD.md`.

## Remote And Production Mutations

Ask for explicit current-task approval before:

- Condor submission, rerun, hold/release/remove;
- merge/production scripts;
- remote SDCC edits or permission changes;
- destructive cleanup;
- broad upload/download or overwrite-prone transfer;
- `--commit-push` for substantive uploader actions.

Read-only diagnostics are allowed when the task needs them and the relevant
SDCC policy is followed.

## Secrets

Never ask for, read, store, repeat, or type SDCC passwords, tokens, private
keys, or other secrets. If an interactive helper prompts for a password, Justin
enters it manually.

## Google Slides

Do not use Computer Use on Google Slides unless it is genuinely useful and
Justin gives explicit permission for that session. Never remove, overwrite, or
repurpose the `Backup` slide.

## Git And Local Files

- Preserve unrelated user changes.
- Do not run destructive git commands such as `git reset --hard` or
  `git checkout --` unless Justin explicitly asks.
- Before editing files, state the files and the intent.
- Use `apply_patch` for manual edits.

## Scientific Claims

Do not call an output valid, stale, finished, production-ready, or slide-ready
without concrete evidence. Use ROOT inspection, timestamps, job IDs, terminal
output, Gmail stage messages, or a clear user statement.
