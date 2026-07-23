#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
submitter="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
pp_wrapper="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor.sh"
auau_wrapper="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_AuAu.sh"

err() { printf 'ERROR: %s\n' "$*" >&2; }
test_sha256() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

# Exercise the implementation from the submitter rather than a test-local
# duplicate. The next function is the stable extraction boundary.
eval "$(
  sed -n \
    '/^normalize_frozen_snapshot_wrapper()/,/^create_pipeline_snapshot()/p' \
    "$submitter" | sed '$d'
)"

tmp="$(mktemp -d "${TMPDIR:-/tmp}/rj_snapshot_loader_test.XXXXXX")"
trap 'chmod -R u+w "$tmp" 2>/dev/null || true; rm -rf "$tmp"' EXIT

[[ "$(condor_getenv_directive)" == True ]] || {
  printf 'FAIL default Condor getenv contract changed\n' >&2
  exit 1
}
[[ "$(RJ_CONDOR_SEALED_ENVIRONMENT=1 condor_getenv_directive)" == False ]] || {
  printf 'FAIL controller-only sealed Condor environment was not enabled\n' >&2
  exit 1
}
if RJ_CONDOR_SEALED_ENVIRONMENT=invalid \
  condor_getenv_directive >"$tmp/getenv.stdout" 2>"$tmp/getenv.stderr"; then
  printf 'FAIL invalid sealed-environment mode was accepted\n' >&2
  exit 1
fi
grep -Fq 'RJ_CONDOR_SEALED_ENVIRONMENT must be 0 or 1' "$tmp/getenv.stderr"
[[ "$(grep -Fxc 'getenv        = $(condor_getenv_directive)' "$submitter" || true)" == 1 ]] || {
  printf 'FAIL sealed Condor environment was not scoped to one descriptor path\n' >&2
  exit 1
}
[[ "$(grep -Fxc 'getenv        = True' "$submitter" || true)" == 3 ]] || {
  printf 'FAIL unrelated Condor descriptor getenv defaults changed\n' >&2
  exit 1
}

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

assert_one_regex() {
  local needle="$1"
  local path="$2"
  local count
  count="$(grep -Ec "$needle" "$path" || true)"
  [[ "$count" -eq 1 ]] || {
    printf 'FAIL expected one regex occurrence of %q in %s, got %s\n' \
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
pin_frozen_snapshot_release "$tmp/pp.sh" ana.560 /release/release_ana/ana.560
assert_one '# RJ_SNAPSHOT_LIBRARY_PREPEND_V1' "$tmp/pp.sh"
assert_one \
  '  export LD_LIBRARY_PATH="${snapshot_lib_dir}${snapshot_loader_suffix}:${LD_LIBRARY_PATH:-}"' \
  "$tmp/pp.sh"
assert_one 'snapshot_loader_suffix=":/release/lib64:/release/lib"' "$tmp/pp.sh"
assert_one_regex \
  '^[[:space:]]*# RJ_PINNED_SPHENIX_RELEASE_V1[[:space:]]*$' \
  "$tmp/pp.sh"
assert_one_regex \
  '^[[:space:]]*source /opt/sphenix/core/bin/sphenix_setup\.sh -n ana\.560[[:space:]]*$' \
  "$tmp/pp.sh"
[[ "$(grep -Ec '^[[:space:]]*source /opt/sphenix/core/bin/sphenix_setup\.sh -n[[:space:]]*$' "$tmp/pp.sh" || true)" == 0 ]] ||
  { printf 'FAIL p+p wrapper retained an unversioned release setup\n' >&2; exit 1; }

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
pin_frozen_snapshot_release "$tmp/auau.sh" ana.560 /release/release_ana/ana.560
assert_one '# RJ_SNAPSHOT_LIBRARY_PREPEND_V1' "$tmp/auau.sh"
assert_one \
  '  export LD_LIBRARY_PATH="${snapshot_lib_dir}${snapshot_loader_suffix}:${LD_LIBRARY_PATH:-}"' \
  "$tmp/auau.sh"
assert_one_regex \
  '^[[:space:]]*source /opt/sphenix/core/bin/sphenix_setup\.sh -n ana\.560[[:space:]]*$' \
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

# The actual THE-134 CaloReco build advertises libcalo_reco.so.0. Prove the
# frozen snapshot materializer creates that loader name as a symlink to the
# one regular candidate ELF and fails closed on SONAME or provider drift.
mkdir -p "$tmp/soname-bin" "$tmp/soname-pass"
cat > "$tmp/soname-bin/readelf" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${FAKE_SONAME_DRIFT:-0}" == "1" ]]; then
  printf ' 0x000000000000000e (SONAME) Library soname: [libcalo_reco.so.9]\n'
else
  printf ' 0x000000000000000e (SONAME) Library soname: [libcalo_reco.so.0]\n'
fi
EOF
chmod +x "$tmp/soname-bin/readelf"
printf 'one-candidate-elf\n' > "$tmp/soname-pass/libcalo_reco.so"
PATH="$tmp/soname-bin:$PATH" \
  stage_snapshot_soname_aliases \
    "$tmp/soname-pass" 1 libcalo_reco.so.0
[[ -f "$tmp/soname-pass/libcalo_reco.so" &&
   ! -L "$tmp/soname-pass/libcalo_reco.so" &&
   -L "$tmp/soname-pass/libcalo_reco.so.0" &&
   "$(readlink "$tmp/soname-pass/libcalo_reco.so.0")" == libcalo_reco.so &&
   "$tmp/soname-pass/libcalo_reco.so.0" -ef "$tmp/soname-pass/libcalo_reco.so" ]] ||
  { printf 'FAIL deterministic CaloReco SONAME alias was not staged\n' >&2; exit 1; }
soname_hash="$(test_sha256 "$tmp/soname-pass/libcalo_reco.so.0")"
api_hash="$(test_sha256 "$tmp/soname-pass/libcalo_reco.so")"
[[ "$soname_hash" == "$api_hash" ]] ||
  { printf 'FAIL CaloReco API/SONAME names do not hash to one provider\n' >&2; exit 1; }

mkdir -p "$tmp/soname-drift"
printf 'one-candidate-elf\n' > "$tmp/soname-drift/libcalo_reco.so"
if FAKE_SONAME_DRIFT=1 PATH="$tmp/soname-bin:$PATH" \
  stage_snapshot_soname_aliases \
    "$tmp/soname-drift" 1 libcalo_reco.so.0 \
    >"$tmp/soname-drift.stdout" 2>"$tmp/soname-drift.stderr"; then
  printf 'FAIL CaloReco SONAME drift was accepted\n' >&2
  exit 1
fi
grep -Fq 'Pinned CaloReco SONAME differs' "$tmp/soname-drift.stderr"

mkdir -p "$tmp/soname-collision"
printf 'one-candidate-elf\n' > "$tmp/soname-collision/libcalo_reco.so"
printf 'different-provider\n' > "$tmp/soname-collision/libcalo_reco.so.0"
if PATH="$tmp/soname-bin:$PATH" \
  stage_snapshot_soname_aliases \
    "$tmp/soname-collision" 1 libcalo_reco.so.0 \
    >"$tmp/soname-collision.stdout" 2>"$tmp/soname-collision.stderr"; then
  printf 'FAIL a second CaloReco SONAME provider was accepted\n' >&2
  exit 1
fi
grep -Fq 'SONAME alias collides with a different provider' \
  "$tmp/soname-collision.stderr"

# Exercise loader closure with a deterministic fake ldd. A staged SONAME must
# resolve into the snapshot, and the exact old mutable-install resolution must
# be rejected.
mkdir -p "$tmp/bin" "$tmp/snapshot" "$tmp/release/lib" "$tmp/release/lib64"
printf 'analysis\n' > "$tmp/snapshot/libRecoilJets.so"
printf 'caloreco\n' > "$tmp/snapshot/libcalo_reco.so"
ln -s libcalo_reco.so "$tmp/snapshot/libcalo_reco.so.0"
printf 'calo-io\n' > "$tmp/release/lib64/libcalo_io.so"
printf 'clusteriso\n' > "$tmp/release/lib64/libclusteriso.so"
printf 'jetbase\n' > "$tmp/release/lib64/libjetbase.so"
printf 'clusteriso-alternate\n' > "$tmp/release/lib/libclusteriso_alternate.so"
release_calo_io="$tmp/release/lib64/libcalo_io.so"
release_clusteriso="$tmp/release/lib64/libclusteriso.so"
release_jetbase="$tmp/release/lib64/libjetbase.so"
release_calo_io_sha="$(test_sha256 "$release_calo_io")"
release_clusteriso_sha="$(test_sha256 "$release_clusteriso")"
release_jetbase_sha="$(test_sha256 "$release_jetbase")"
cat > "$tmp/bin/ldd" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
target="${1:?target required}"
snap="$(dirname "$target")"
release_lib="${FAKE_RELEASE_LIB:?release lib required}"
if [[ "$(basename "$target")" == "libRecoilJets.so" ]]; then
  if [[ "${FAKE_MUTABLE_RESOLUTION:-0}" == "1" ]]; then
    printf '%s\n' \
      'libcalo_reco.so.0 => /sphenix/user/patsfan753/install/lib/libcalo_reco.so.0 (0x1)'
  elif [[ "${FAKE_STAGE_ESCAPE:-0}" == "1" ]]; then
    printf '%s\n' \
      'libcalo_reco.so.0 => /opt/alien/libcalo_reco.so.0 (0x1)'
  else
    printf 'libcalo_reco.so.0 => %s/libcalo_reco.so.0 (0x1)\n' "$snap"
    if [[ "${FAKE_RELEASE_COMPANION_ESCAPE:-0}" == "1" ]]; then
      printf '%s\n' \
        'libcalo_io.so.0 => /opt/alien/libcalo_io.so.0 (0x1)'
    else
      printf 'libcalo_io.so.0 => %s/libcalo_io.so (0x1)\n' "$release_lib"
      if [[ "${FAKE_RELEASE_PROVIDER_SWAP:-0}" == "1" ]]; then
        printf 'libclusteriso.so.0 => %s/../lib/libclusteriso_alternate.so (0x1)\n' "$release_lib"
      else
        printf 'libclusteriso.so.0 => %s/libclusteriso.so (0x1)\n' "$release_lib"
      fi
      printf 'libjetbase.so.0 => %s/libjetbase.so (0x1)\n' "$release_lib"
    fi
  fi
else
  printf '%s\n' 'libc.so.6 => /lib64/libc.so.6 (0x2)'
fi
EOF
chmod +x "$tmp/bin/ldd"
FAKE_RELEASE_LIB="$tmp/release/lib64" PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
  "$tmp/snapshot" pp 0 "$tmp/release/lib64" "$tmp/release/lib"
if FAKE_RELEASE_LIB="$tmp/release/lib64" FAKE_MUTABLE_RESOLUTION=1 PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 0 "$tmp/release/lib64" "$tmp/release/lib" \
    >"$tmp/mutable.stdout" 2>"$tmp/mutable.stderr"; then
  printf 'FAIL mutable private-install loader resolution was accepted\n' >&2
  exit 1
fi
grep -Fq 'mutable private install dependency' "$tmp/mutable.stderr"
if FAKE_RELEASE_LIB="$tmp/release/lib64" FAKE_STAGE_ESCAPE=1 PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 0 "$tmp/release/lib64" "$tmp/release/lib" \
    >"$tmp/escape.stdout" 2>"$tmp/escape.stderr"; then
  printf 'FAIL staged SONAME escape was accepted\n' >&2
  exit 1
fi
grep -Fq 'staged dependency escaped snapshot' "$tmp/escape.stderr"

# Single-provider mode requires only CaloReco to resolve from the immutable
# snapshot and every companion to resolve to its exact path and hash.
loader_receipt="$tmp/snapshot/snapshot_loader_receipt.json"
RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io" \
RJ_PINNED_RELEASE_CALO_IO_SHA256="$release_calo_io_sha" \
RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso" \
RJ_PINNED_RELEASE_CLUSTERISO_SHA256="$release_clusteriso_sha" \
RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase" \
RJ_PINNED_RELEASE_JETBASE_SHA256="$release_jetbase_sha" \
FAKE_RELEASE_LIB="$tmp/release/lib64" PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
  "$tmp/snapshot" pp 1 "$tmp/release/lib64" "$tmp/release/lib" 1 "$loader_receipt"
python3 - "$loader_receipt" "$release_clusteriso" <<'PY'
from pathlib import Path
import json
import sys

receipt = json.loads(Path(sys.argv[1]).read_text())
assert receipt["schema"] == "RJ_SNAPSHOT_LOADER_RECEIPT_V1"
assert receipt["status"] == "PASS"
assert receipt["mode"] == "pp"
assert receipt["pinned_calo_reco_release_companions"] is True
assert receipt["providers"]["libclusteriso.so"]["realpath"] == str(Path(sys.argv[2]).resolve())
assert receipt["providers"]["libclusteriso.so"]["observed_resolutions"] == [
    str(Path(sys.argv[2]).resolve())
]
PY
rm "$tmp/snapshot/libcalo_reco.so.0"
if RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io" \
  RJ_PINNED_RELEASE_CALO_IO_SHA256="$release_calo_io_sha" \
  RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso" \
  RJ_PINNED_RELEASE_CLUSTERISO_SHA256="$release_clusteriso_sha" \
  RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase" \
  RJ_PINNED_RELEASE_JETBASE_SHA256="$release_jetbase_sha" \
  FAKE_RELEASE_LIB="$tmp/release/lib64" PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 1 "$tmp/release/lib64" "$tmp/release/lib" 1 \
    >"$tmp/missing-soname.stdout" 2>"$tmp/missing-soname.stderr"; then
  printf 'FAIL missing CaloReco SONAME alias was accepted\n' >&2
  exit 1
fi
grep -Fq 'single-provider resolution mismatch' "$tmp/missing-soname.stderr"
ln -s libcalo_reco.so "$tmp/snapshot/libcalo_reco.so.0"

if RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io" \
  RJ_PINNED_RELEASE_CALO_IO_SHA256="$release_calo_io_sha" \
  RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso" \
  RJ_PINNED_RELEASE_CLUSTERISO_SHA256="$release_clusteriso_sha" \
  RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase" \
  RJ_PINNED_RELEASE_JETBASE_SHA256="$release_jetbase_sha" \
  FAKE_RELEASE_LIB="$tmp/release/lib64" FAKE_RELEASE_COMPANION_ESCAPE=1 PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 1 "$tmp/release/lib64" "$tmp/release/lib" 1 \
    >"$tmp/companion-escape.stdout" 2>"$tmp/companion-escape.stderr"; then
  printf 'FAIL pinned release companion escape was accepted\n' >&2
  exit 1
fi
grep -Eq 'single-provider resolution mismatch|staged dependency escaped snapshot' \
  "$tmp/companion-escape.stderr"
if RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io" \
  RJ_PINNED_RELEASE_CALO_IO_SHA256="$release_calo_io_sha" \
  RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso" \
  RJ_PINNED_RELEASE_CLUSTERISO_SHA256="$release_clusteriso_sha" \
  RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase" \
  RJ_PINNED_RELEASE_JETBASE_SHA256="$release_jetbase_sha" \
  FAKE_RELEASE_LIB="$tmp/release/lib64" FAKE_STAGE_ESCAPE=1 PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 1 "$tmp/release/lib64" "$tmp/release/lib" 1 \
    >"$tmp/pinned-calo-escape.stdout" 2>"$tmp/pinned-calo-escape.stderr"; then
  printf 'FAIL pinned CaloReco snapshot escape was accepted\n' >&2
  exit 1
fi
grep -Eq 'single-provider resolution mismatch|staged dependency escaped snapshot' \
  "$tmp/pinned-calo-escape.stderr"

if RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io" \
  RJ_PINNED_RELEASE_CALO_IO_SHA256="$release_calo_io_sha" \
  RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso" \
  RJ_PINNED_RELEASE_CLUSTERISO_SHA256="$release_clusteriso_sha" \
  RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase" \
  RJ_PINNED_RELEASE_JETBASE_SHA256="$release_jetbase_sha" \
  FAKE_RELEASE_LIB="$tmp/release/lib64" FAKE_RELEASE_PROVIDER_SWAP=1 PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 1 "$tmp/release/lib64" "$tmp/release/lib" 1 \
    >"$tmp/provider-swap.stdout" 2>"$tmp/provider-swap.stderr"; then
  printf 'FAIL same-release provider swap was accepted\n' >&2
  exit 1
fi
grep -Fq 'single-provider resolution mismatch' "$tmp/provider-swap.stderr"

if RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io" \
  RJ_PINNED_RELEASE_CALO_IO_SHA256="$release_calo_io_sha" \
  RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso" \
  RJ_PINNED_RELEASE_CLUSTERISO_SHA256="$(printf '0%.0s' {1..64})" \
  RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase" \
  RJ_PINNED_RELEASE_JETBASE_SHA256="$release_jetbase_sha" \
  FAKE_RELEASE_LIB="$tmp/release/lib64" PATH="$tmp/bin:$PATH" \
  validate_snapshot_loader_closure \
    "$tmp/snapshot" pp 1 "$tmp/release/lib64" "$tmp/release/lib" 1 \
    >"$tmp/provider-hash.stdout" 2>"$tmp/provider-hash.stderr"; then
  printf 'FAIL changed release-provider hash was accepted\n' >&2
  exit 1
fi
grep -Fq 'declared release companion hash drift' "$tmp/provider-hash.stderr"

mkdir -p "$tmp/seal/lib"
printf 'sealed-analysis\n' > "$tmp/seal/lib/libRecoilJets.so"
ln -s libRecoilJets.so "$tmp/seal/lib/libRecoilJets.so.0"
printf 'sealed-caloreco\n' > "$tmp/seal/lib/libcalo_reco.so"
ln -s libcalo_reco.so "$tmp/seal/lib/libcalo_reco.so.0"
write_and_seal_snapshot_manifest "$tmp/seal"
python3 - "$tmp/seal" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1]).resolve()
manifest = json.loads((root / "snapshot_manifest.json").read_text())
assert manifest["schema"] == "RJ_FROZEN_SNAPSHOT_MANIFEST_V1"
assert manifest["status"] == "PASS"
assert manifest["root"] == str(root)
assert {row["path"] for row in manifest["entries"]} == {
    "lib",
    "lib/libRecoilJets.so",
    "lib/libRecoilJets.so.0",
    "lib/libcalo_reco.so",
    "lib/libcalo_reco.so.0",
}
assert root.stat().st_mode & 0o222 == 0
linux_symlink_fixture_count = 0
for path in root.rglob("*"):
    if path.is_symlink():
        # Linux reports symlink lstat mode as 0777 even after chmod -R. Model
        # that explicitly: the link bits are irrelevant; the sealed parent,
        # relative target, in-root resolution, and sealed target are authority.
        linux_lstat_mode = path.lstat().st_mode | 0o777
        assert linux_lstat_mode & 0o222
        target = path.readlink()
        assert not target.is_absolute()
        resolved = path.resolve(strict=True)
        assert root in resolved.parents
        assert path.parent.stat().st_mode & 0o222 == 0
        assert resolved.stat().st_mode & 0o222 == 0
        linux_symlink_fixture_count += 1
    else:
        assert path.lstat().st_mode & 0o222 == 0
assert linux_symlink_fixture_count == 2
PY

mkdir -p "$tmp/seal-absolute/lib"
printf 'absolute-target\n' > "$tmp/seal-absolute/lib/libcalo_reco.so"
ln -s "$tmp/seal-absolute/lib/libcalo_reco.so" \
  "$tmp/seal-absolute/lib/libcalo_reco.so.0"
if write_and_seal_snapshot_manifest \
  "$tmp/seal-absolute" \
  >"$tmp/seal-absolute.stdout" 2>"$tmp/seal-absolute.stderr"; then
  printf 'FAIL absolute snapshot SONAME alias was accepted\n' >&2
  exit 1
fi
grep -Fq 'snapshot symlink must be relative' "$tmp/seal-absolute.stderr"

bash -n "$tmp/pp.sh"
bash -n "$tmp/old_false_positive.sh"
bash -n "$tmp/auau.sh"

printf 'PASS snapshot_loader_contract mutations=19 receipts=2 soname_aliases=1 linux_symlink_fixture=1\n'
