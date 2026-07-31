#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd -P)"
helper="${repo_root}/scripts/sdcc/transfer/sftp_push_recoiljets.sh"

bash -n "$helper"
grep -Fq 'tar -cf "$tar_file" -- "${selected_remote[@]}"' "$helper"
grep -Fq 'test -d ${REMOTE_BASE}' "$helper"
grep -Fq 'before_mode=\$(stat -c %a ${REMOTE_BASE})' "$helper"
grep -Fq 'test \"\$before_mode\" = \"\$after_mode\"' "$helper"
if grep -Fq 'tar -cf "$tar_file" .' "$helper"; then
  printf 'ssh-tar helper still archives its private staging-root entry\n' >&2
  exit 1
fi

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/sftp-archive-scope.XXXXXX")"
trap 'chmod -R u+w "$tmpdir" 2>/dev/null || true; rm -rf "$tmpdir"' EXIT
stage="${tmpdir}/stage"
remote="${tmpdir}/remote"
payload="${tmpdir}/payload.tar"
mkdir -m 0700 "$stage"
mkdir -m 0755 "$remote"
mkdir -p "${stage}/scripts/sdcc/workflows/diagnostics"
printf 'bounded payload\n' \
  > "${stage}/scripts/sdcc/workflows/diagnostics/controller.sh"
(
  cd "$stage"
  COPYFILE_DISABLE=1 tar -cf "$payload" -- \
    scripts/sdcc/workflows/diagnostics/controller.sh
)

before_mode="$(stat -f '%Lp' "$remote" 2>/dev/null || stat -c '%a' "$remote")"
tar -xf "$payload" -C "$remote"
after_mode="$(stat -f '%Lp' "$remote" 2>/dev/null || stat -c '%a' "$remote")"
[[ "$before_mode" == 755 && "$after_mode" == "$before_mode" ]] || {
  printf 'selected-file archive changed remote-root mode: before=%s after=%s\n' \
    "$before_mode" "$after_mode" >&2
  exit 1
}
[[ -s "${remote}/scripts/sdcc/workflows/diagnostics/controller.sh" ]]

printf 'SFTP_PUSH_ARCHIVE_SCOPE_PASS remote_mode=%s\n' "$after_mode"
