#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
submitter="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
pp_wrapper="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor.sh"
auau_wrapper="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_AuAu.sh"

err() { printf 'ERROR: %s\n' "$*" >&2; }

# Exercise the implementation from the submitter rather than a test-local
# duplicate. The next function is the stable extraction boundary.
eval "$(
  sed -n \
    '/^normalize_frozen_snapshot_wrapper()/,/^create_pipeline_snapshot()/p' \
    "$submitter" | sed '$d'
)"

tmp="$(mktemp -d "${TMPDIR:-/tmp}/rj_snapshot_loader_test.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

assert_one() {
  local needle="$1"
  local path="$2"
  local count
  count="$(grep -Fxc "$needle" "$path" || true)"
  [[ "$count" -eq 1 ]] || {
    printf 'FAIL expected one occurrence of %q in %s, got %s\n' \
      "$needle" "$path" "$count" >&2
    exit 1
  }
}

# The live p+p wrapper must remain unchanged; only its frozen copy receives the
# snapshot-library prepend.
if grep -Fq '# RJ_SNAPSHOT_LIBRARY_PREPEND_V1' "$pp_wrapper"; then
  printf 'FAIL live p+p wrapper owns a frozen-snapshot prepend marker\n' >&2
  exit 1
fi
cp "$pp_wrapper" "$tmp/pp.sh"
normalize_frozen_snapshot_wrapper "$tmp/pp.sh" ":/release/lib64:/release/lib"
validate_frozen_snapshot_wrapper_contract "$tmp/pp.sh"
assert_one '# RJ_SNAPSHOT_LIBRARY_PREPEND_V1' "$tmp/pp.sh"
assert_one \
  '  export LD_LIBRARY_PATH="${snapshot_lib_dir}${snapshot_loader_suffix}:${LD_LIBRARY_PATH:-}"' \
  "$tmp/pp.sh"
assert_one 'snapshot_loader_suffix=":/release/lib64:/release/lib"' "$tmp/pp.sh"

# Reproduce the V4 false positive: snapshot_lib_dir exists for the builder
# override, but no LD_LIBRARY_PATH prepend exists. The materializer must add
# the marked prepend rather than treating the variable as proof of closure.
cat > "$tmp/old_false_positive.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
wrapper_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
snapshot_lib_dir="${RJ_SNAPSHOT_LIB_DIR:-${wrapper_dir}/lib}"
photon_builder_override="${snapshot_lib_dir}/libphoton_cluster_builder_override.so"
# Some frozen release lanes intentionally pair a diagnostic override.
# ------------------------ Dataset routing ------------------
rc=125
EOF
normalize_frozen_snapshot_wrapper "$tmp/old_false_positive.sh" ""
validate_frozen_snapshot_wrapper_contract "$tmp/old_false_positive.sh"
assert_one '# RJ_SNAPSHOT_LIBRARY_PREPEND_V1' "$tmp/old_false_positive.sh"
assert_one \
  '  export LD_LIBRARY_PATH="${snapshot_lib_dir}${snapshot_loader_suffix}:${LD_LIBRARY_PATH:-}"' \
  "$tmp/old_false_positive.sh"
marker_line="$(grep -Fn '# RJ_SNAPSHOT_LIBRARY_PREPEND_V1' "$tmp/old_false_positive.sh" | cut -d: -f1)"
loader_line="$(grep -Fn '# Some frozen release lanes' "$tmp/old_false_positive.sh" | cut -d: -f1)"
[[ "$marker_line" -lt "$loader_line" ]] || {
  printf 'FAIL snapshot prepend was not installed before loader logic\n' >&2
  exit 1
}

# The existing Au+Au wrapper has a real legacy prepend. Normalize it in place;
# never add a second export.
cp "$auau_wrapper" "$tmp/auau.sh"
normalize_frozen_snapshot_wrapper "$tmp/auau.sh" ""
validate_frozen_snapshot_wrapper_contract "$tmp/auau.sh"
assert_one '# RJ_SNAPSHOT_LIBRARY_PREPEND_V1' "$tmp/auau.sh"
assert_one \
  '  export LD_LIBRARY_PATH="${snapshot_lib_dir}${snapshot_loader_suffix}:${LD_LIBRARY_PATH:-}"' \
  "$tmp/auau.sh"

# A duplicate marker is a hard preflight failure, not something to repair
# ambiguously.
cp "$tmp/pp.sh" "$tmp/duplicate.sh"
printf '%s\n' '# RJ_SNAPSHOT_LIBRARY_PREPEND_V1' >> "$tmp/duplicate.sh"
if normalize_frozen_snapshot_wrapper \
  "$tmp/duplicate.sh" "" 2>"$tmp/duplicate.stderr"; then
  printf 'FAIL duplicate snapshot prepend marker was accepted\n' >&2
  exit 1
fi
grep -Fq 'duplicate prepend markers' "$tmp/duplicate.stderr"

# A correctly placed marker must not hide an export moved below ROOT/dataset
# setup. This is the exact ordering class that could re-admit V4 contamination.
cp "$tmp/old_false_positive.sh" "$tmp/export-after-dataset.sh"
python3 - "$tmp/export-after-dataset.sh" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
export = '  export LD_LIBRARY_PATH="${snapshot_lib_dir}${snapshot_loader_suffix}:${LD_LIBRARY_PATH:-}"'
rows = path.read_text().splitlines()
rows.remove(export)
rows.append(export)
path.write_text("\n".join(rows) + "\n")
PY
if validate_frozen_snapshot_wrapper_contract \
  "$tmp/export-after-dataset.sh" 2>"$tmp/export-after-dataset.stderr"; then
  printf 'FAIL export moved after dataset routing was accepted\n' >&2
  exit 1
fi
grep -Fq 'loader order must be' "$tmp/export-after-dataset.stderr"

# Exercise loader closure with a deterministic fake ldd. A staged SONAME must
# resolve into the snapshot, and the exact old mutable-install resolution must
# be rejected.
mkdir -p "$tmp/bin" "$tmp/snapshot"
: > "$tmp/snapshot/libRecoilJets.so"
: > "$tmp/snapshot/libcalo_reco.so"
ln -s libcalo_reco.so "$tmp/snapshot/libcalo_reco.so.0"
cat > "$tmp/bin/ldd" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
target="${1:?target required}"
snap="$(dirname "$target")"
if [[ "$(basename "$target")" == "libRecoilJets.so" ]]; then
  if [[ "${FAKE_MUTABLE_RESOLUTION:-0}" == "1" ]]; then
    printf '%s\n' \
      'libcalo_reco.so.0 => /sphenix/user/patsfan753/install/lib/libcalo_reco.so.0 (0x1)'
  elif [[ "${FAKE_STAGE_ESCAPE:-0}" == "1" ]]; then
    printf '%s\n' \
      'libcalo_reco.so.0 => /opt/alien/libcalo_reco.so.0 (0x1)'
  else
    printf 'libcalo_reco.so.0 => %s/libcalo_reco.so.0 (0x1)\n' "$snap"
  fi
else
  printf '%s\n' 'libc.so.6 => /lib64/libc.so.6 (0x2)'
fi
EOF
chmod +x "$tmp/bin/ldd"
PATH="$tmp/bin:$PATH" validate_snapshot_loader_closure \
  "$tmp/snapshot" pp 0 /release/lib64 /release/lib
if FAKE_MUTABLE_RESOLUTION=1 PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 0 /release/lib64 /release/lib \
    >"$tmp/mutable.stdout" 2>"$tmp/mutable.stderr"; then
  printf 'FAIL mutable private-install loader resolution was accepted\n' >&2
  exit 1
fi
grep -Fq 'mutable private install dependency' "$tmp/mutable.stderr"
if FAKE_STAGE_ESCAPE=1 PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 0 /release/lib64 /release/lib \
    >"$tmp/escape.stdout" 2>"$tmp/escape.stderr"; then
  printf 'FAIL staged SONAME escape was accepted\n' >&2
  exit 1
fi
grep -Fq 'staged dependency escaped snapshot' "$tmp/escape.stderr"

bash -n "$tmp/pp.sh"
bash -n "$tmp/old_false_positive.sh"
bash -n "$tmp/auau.sh"

printf 'PASS snapshot_loader_contract mutations=4\n'
