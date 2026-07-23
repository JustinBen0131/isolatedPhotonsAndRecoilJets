#!/usr/bin/env bash
# Build the THE-134 ana.560 single-provider CaloReco runtime.
#
# The authoritative CaloReco package is exported from one exact coresoftware
# commit.  Only the two hash-pinned PhotonClusterBuilder sources are overlaid.
# RawClusterBuilderTopo is then changed mechanically to use the detector-
# explicit TowerInfoDefs channel maps that a correctly typed container already
# delegates to.  No energy, geometry, threshold, ordering, or clustering
# behavior is changed.
#
# Default mode is read-only and emits action-specific authorization tokens.
# --stage-source exists so the archive/overlay/mapping contract can be tested
# without an sPHENIX runtime.  It is explicitly not runtime authority.

set -euo pipefail
umask 022

readonly SCHEMA_VERSION="THE134_ANA560_CALORECO_BUILD_CONTRACT_V2"
readonly SOURCE_SCHEMA="THE134_ANA560_CALORECO_SOURCE_MANIFEST_V2"
readonly RECEIPT_SCHEMA="THE134_ANA560_CALORECO_BUILD_RECEIPT_V2"
readonly MAPPING_PATCH_ID="THE134_RAWCLUSTERBUILDERTOPO_DETECTOR_EXPLICIT_CHANNEL_MAP_V1"
readonly ROLE_SIZE_GUARD_ID="THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_ROLE_SIZE_GUARD_V1"
readonly ABI_ADDITIONS_PROVENANCE_SHA256="5d4eca4abdaa274d308856e050d19b17e02bc62652a03dda6f0ea95719b9147e"
readonly ABI_INTENTIONAL_EXCLUSION="_ZN18TowerInfoContainer10encode_keyEj W"
readonly ABI_ROLE_GUARD_INLINE_ADDITION="_ZNK18TowerInfoContainer14get_detectoridEv W"
readonly ABI_ADDITIONS_ALLOWLIST_SHA256="700974da06ab8268fcbc750564a72f1db4eb3fcfdb8e351130d71ea1ce963ed8"
readonly EXPECTED_SONAME="libcalo_reco.so.0"
readonly EXPECTED_CORESOFTWARE_COMMIT="cba274033b5560e32600cdeaa7676b6ab4a6c971"
readonly EXPECTED_SETUP_SCRIPT="/opt/sphenix/core/bin/sphenix_setup.sh"
readonly EXPECTED_OFFLINE_MAIN="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560"
readonly RELEASE_NAME="ana.560"
readonly BUILDER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
readonly DEFAULT_ABI_ADDITIONS_ALLOWLIST="${BUILDER_DIR}/contracts/the134_ana560_calo_reco_abi_additions.allowlist"

die() {
  printf 'THE134_ANA560_CALORECO_BUILD_FAIL: %s\n' "$*" >&2
  exit 2
}

usage() {
  cat <<'EOF'
Usage:
  build_the134_ana560_calo_reco_runtime.sh \
    --output-root ABS --coresoftware-repo ABS --coresoftware-commit SHA \
    --photon-builder-cc ABS --photon-builder-cc-sha256 SHA256 \
    --photon-builder-h ABS --photon-builder-h-sha256 SHA256 \
    --baseline-calo-reco-sha256 SHA256 \
    --abi-additions-allowlist-sha256 SHA256 \
    --towerinfo-defs-header-sha256 SHA256 \
    --towerinfo-container-header-sha256 SHA256 \
    --towerinfo-containerv1-header-sha256 SHA256 \
    --calo-io-library-sha256 SHA256 \
    [--baseline-calo-reco ABS] [--towerinfo-defs-header ABS] \
    [--towerinfo-container-header ABS] [--towerinfo-containerv1-header ABS] \
    [--calo-io-library ABS] [--abi-additions-allowlist ABS] \
    [--setup-script ABS] [--expected-offline-main ABS] [--jobs N]

  build_the134_ana560_calo_reco_runtime.sh --stage-source --token TOKEN \
    <same arguments>

  build_the134_ana560_calo_reco_runtime.sh --build --token TOKEN \
    <same arguments>

Default mode is read-only and prints separate source-stage and build tokens.
--stage-source exports, overlays, patches, manifests, and seals source only;
it emits no runtime receipt and is never production authority.
--build performs the full ana.560 build and ABI/provider/runtime validation.

The output root must not exist and must be below a .recoiljets_tmp directory
or /tmp.  --setup-script must equal /opt/sphenix/core/bin/sphenix_setup.sh and
--expected-offline-main must equal
/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560.
The coresoftware revision must be exactly
cba274033b5560e32600cdeaa7676b6ab4a6c971 and available in the specified
repository.  No other commit is accepted.  Mutable checkout contents are never
copied.  The baseline library and TowerInfo authority paths default below the
exact --expected-offline-main prefix.  Their expected hashes remain mandatory.
EOF
}

sha256_file() {
  python3 - "$1" <<'PY'
from pathlib import Path
import hashlib
import sys

digest = hashlib.sha256()
with Path(sys.argv[1]).open("rb") as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b""):
        digest.update(block)
print(digest.hexdigest())
PY
}

require_sha256() {
  local label="$1"
  local value="$2"
  [[ "$value" =~ ^[0-9a-f]{64}$ ]] || die "${label} is not a lowercase SHA-256"
}

require_file_hash() {
  local label="$1"
  local path="$2"
  local expected="$3"
  [[ "$path" == /* && -f "$path" && -s "$path" ]] || \
    die "${label} is missing or empty: ${path}"
  local actual
  actual="$(sha256_file "$path")"
  [[ "$actual" == "$expected" ]] || \
    die "${label} hash differs: expected ${expected}, got ${actual}"
}

require_regular_file_hash() {
  local label="$1"
  local path="$2"
  local expected="$3"
  [[ "$path" == /* && -f "$path" ]] || \
    die "${label} is missing: ${path}"
  local actual
  actual="$(sha256_file "$path")"
  [[ "$actual" == "$expected" ]] || \
    die "${label} hash differs: expected ${expected}, got ${actual}"
}

mode="plan"
provided_token=""
output_root=""
coresoftware_repo=""
coresoftware_commit=""
photon_builder_cc=""
photon_builder_cc_sha256=""
photon_builder_h=""
photon_builder_h_sha256=""
baseline_calo_reco=""
baseline_calo_reco_sha256=""
abi_additions_allowlist="$DEFAULT_ABI_ADDITIONS_ALLOWLIST"
abi_additions_allowlist_sha256=""
towerinfo_defs_header=""
towerinfo_defs_header_sha256=""
towerinfo_container_header=""
towerinfo_container_header_sha256=""
towerinfo_containerv1_header=""
towerinfo_containerv1_header_sha256=""
calo_io_library=""
calo_io_library_sha256=""
setup_script="$EXPECTED_SETUP_SCRIPT"
expected_offline_main="$EXPECTED_OFFLINE_MAIN"
jobs=4

while (($#)); do
  case "$1" in
    --stage-source)
      [[ "$mode" == "plan" ]] || die "--stage-source and --build are mutually exclusive"
      mode="source"
      shift
      ;;
    --build)
      [[ "$mode" == "plan" ]] || die "--stage-source and --build are mutually exclusive"
      mode="build"
      shift
      ;;
    --token)
      [[ $# -ge 2 ]] || die "--token requires a value"
      provided_token="$2"
      shift 2
      ;;
    --output-root)
      [[ $# -ge 2 ]] || die "--output-root requires a value"
      output_root="$2"
      shift 2
      ;;
    --coresoftware-repo)
      [[ $# -ge 2 ]] || die "--coresoftware-repo requires a value"
      coresoftware_repo="$2"
      shift 2
      ;;
    --coresoftware-commit)
      [[ $# -ge 2 ]] || die "--coresoftware-commit requires a value"
      coresoftware_commit="$2"
      shift 2
      ;;
    --photon-builder-cc)
      [[ $# -ge 2 ]] || die "--photon-builder-cc requires a value"
      photon_builder_cc="$2"
      shift 2
      ;;
    --photon-builder-cc-sha256)
      [[ $# -ge 2 ]] || die "--photon-builder-cc-sha256 requires a value"
      photon_builder_cc_sha256="$2"
      shift 2
      ;;
    --photon-builder-h)
      [[ $# -ge 2 ]] || die "--photon-builder-h requires a value"
      photon_builder_h="$2"
      shift 2
      ;;
    --photon-builder-h-sha256)
      [[ $# -ge 2 ]] || die "--photon-builder-h-sha256 requires a value"
      photon_builder_h_sha256="$2"
      shift 2
      ;;
    --baseline-calo-reco)
      [[ $# -ge 2 ]] || die "--baseline-calo-reco requires a value"
      baseline_calo_reco="$2"
      shift 2
      ;;
    --baseline-calo-reco-sha256)
      [[ $# -ge 2 ]] || die "--baseline-calo-reco-sha256 requires a value"
      baseline_calo_reco_sha256="$2"
      shift 2
      ;;
    --abi-additions-allowlist)
      [[ $# -ge 2 ]] || die "--abi-additions-allowlist requires a value"
      abi_additions_allowlist="$2"
      shift 2
      ;;
    --abi-additions-allowlist-sha256)
      [[ $# -ge 2 ]] || die "--abi-additions-allowlist-sha256 requires a value"
      abi_additions_allowlist_sha256="$2"
      shift 2
      ;;
    --towerinfo-defs-header)
      [[ $# -ge 2 ]] || die "--towerinfo-defs-header requires a value"
      towerinfo_defs_header="$2"
      shift 2
      ;;
    --towerinfo-defs-header-sha256)
      [[ $# -ge 2 ]] || die "--towerinfo-defs-header-sha256 requires a value"
      towerinfo_defs_header_sha256="$2"
      shift 2
      ;;
    --towerinfo-container-header)
      [[ $# -ge 2 ]] || die "--towerinfo-container-header requires a value"
      towerinfo_container_header="$2"
      shift 2
      ;;
    --towerinfo-container-header-sha256)
      [[ $# -ge 2 ]] || die "--towerinfo-container-header-sha256 requires a value"
      towerinfo_container_header_sha256="$2"
      shift 2
      ;;
    --towerinfo-containerv1-header)
      [[ $# -ge 2 ]] || die "--towerinfo-containerv1-header requires a value"
      towerinfo_containerv1_header="$2"
      shift 2
      ;;
    --towerinfo-containerv1-header-sha256)
      [[ $# -ge 2 ]] || die "--towerinfo-containerv1-header-sha256 requires a value"
      towerinfo_containerv1_header_sha256="$2"
      shift 2
      ;;
    --calo-io-library)
      [[ $# -ge 2 ]] || die "--calo-io-library requires a value"
      calo_io_library="$2"
      shift 2
      ;;
    --calo-io-library-sha256)
      [[ $# -ge 2 ]] || die "--calo-io-library-sha256 requires a value"
      calo_io_library_sha256="$2"
      shift 2
      ;;
    --setup-script)
      [[ $# -ge 2 ]] || die "--setup-script requires a value"
      setup_script="$2"
      shift 2
      ;;
    --expected-offline-main)
      [[ $# -ge 2 ]] || die "--expected-offline-main requires a value"
      expected_offline_main="$2"
      shift 2
      ;;
    --jobs)
      [[ $# -ge 2 ]] || die "--jobs requires a value"
      jobs="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

[[ "$setup_script" == "$EXPECTED_SETUP_SCRIPT" ]] || \
  die "--setup-script must equal frozen THE-134 authority ${EXPECTED_SETUP_SCRIPT}"
[[ "$expected_offline_main" == "$EXPECTED_OFFLINE_MAIN" ]] || \
  die "--expected-offline-main must equal frozen THE-134 authority ${EXPECTED_OFFLINE_MAIN}"
[[ "$jobs" =~ ^[1-8]$ ]] || die "--jobs must be an integer from 1 through 8"
[[ "$output_root" == /* && "$output_root" != *$'\n'* && "$output_root" != *$'\r'* \
  && "$output_root" != *' '* ]] || die "--output-root must be absolute and contain no whitespace"
output_root="$(
  python3 - "$output_root" <<'PY'
from pathlib import Path
import sys

raw = Path(sys.argv[1])
if any(part in {".", ".."} for part in raw.parts):
    raise SystemExit("dot path components are forbidden")
parent = raw.parent.resolve(strict=True)
candidate = parent / raw.name

def below(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return path != root
    except ValueError:
        return False

tmp_root = Path("/tmp").resolve(strict=True)
allowed = below(candidate, tmp_root)
parts = raw.parts
for index, part in enumerate(parts):
    if part != ".recoiljets_tmp":
        continue
    lexical_root = Path(*parts[: index + 1])
    resolved_root = lexical_root.resolve(strict=True)
    if below(candidate, resolved_root):
        allowed = True
if not allowed:
    raise SystemExit("not below canonical /tmp or .recoiljets_tmp")
print(candidate)
PY
)" || die "--output-root failed canonical scratch-boundary validation"
[[ ! -e "$output_root" ]] || die "output root already exists: ${output_root}"

[[ "$coresoftware_repo" == /* && -d "$coresoftware_repo" ]] || \
  die "--coresoftware-repo must be an existing absolute directory"
coresoftware_repo_real="$(git -C "$coresoftware_repo" rev-parse --show-toplevel 2>/dev/null)" || \
  die "--coresoftware-repo is not a Git checkout"
coresoftware_repo_real="$(cd "$coresoftware_repo_real" && pwd -P)"
git_safe=(-c "safe.directory=${coresoftware_repo_real}")
[[ "$coresoftware_commit" =~ ^[0-9a-f]{40}$ ]] || \
  die "--coresoftware-commit must be one full lowercase commit SHA"
[[ "$coresoftware_commit" == "$EXPECTED_CORESOFTWARE_COMMIT" ]] || \
  die "--coresoftware-commit must equal frozen THE-134 authority ${EXPECTED_CORESOFTWARE_COMMIT}"
git "${git_safe[@]}" -C "$coresoftware_repo_real" cat-file -e \
  "${coresoftware_commit}^{commit}" 2>/dev/null || \
  die "coresoftware commit is not present: ${coresoftware_commit}"
resolved_commit="$(git "${git_safe[@]}" -C "$coresoftware_repo_real" rev-parse \
  "${coresoftware_commit}^{commit}")"
[[ "$resolved_commit" == "$coresoftware_commit" ]] || \
  die "coresoftware revision did not resolve to the exact requested commit"
coresoftware_tree="$(git "${git_safe[@]}" -C "$coresoftware_repo_real" rev-parse \
  "${coresoftware_commit}^{tree}")"
commit_epoch="$(git "${git_safe[@]}" -C "$coresoftware_repo_real" show -s \
  --format=%ct "$coresoftware_commit")"
[[ "$coresoftware_tree" =~ ^[0-9a-f]{40}$ && "$commit_epoch" =~ ^[0-9]+$ ]] || \
  die "could not resolve coresoftware tree or commit epoch"

required_archive_files=(
  configure.ac
  Makefile.am
  autogen.sh
  PhotonClusterBuilder.cc
  PhotonClusterBuilder.h
  RawClusterBuilderTopo.cc
  RawClusterBuilderTopo.h
)
for archive_file in "${required_archive_files[@]}"; do
  git "${git_safe[@]}" -C "$coresoftware_repo_real" cat-file -e \
    "${coresoftware_commit}:offline/packages/CaloReco/${archive_file}" 2>/dev/null || \
    die "commit lacks offline/packages/CaloReco/${archive_file}"
done

require_sha256 "PhotonClusterBuilder.cc SHA-256" "$photon_builder_cc_sha256"
require_sha256 "PhotonClusterBuilder.h SHA-256" "$photon_builder_h_sha256"
require_file_hash "PhotonClusterBuilder.cc overlay" \
  "$photon_builder_cc" "$photon_builder_cc_sha256"
require_file_hash "PhotonClusterBuilder.h overlay" \
  "$photon_builder_h" "$photon_builder_h_sha256"
[[ "$setup_script" == /* && -f "$setup_script" && -s "$setup_script" ]] || \
  die "setup script is missing or empty: ${setup_script}"
[[ "$expected_offline_main" == /* && "$expected_offline_main" == *"/${RELEASE_NAME}" ]] || \
  die "--expected-offline-main must be an absolute ${RELEASE_NAME} prefix"

baseline_calo_reco="${baseline_calo_reco:-${expected_offline_main}/lib/${EXPECTED_SONAME}}"
towerinfo_defs_header="${towerinfo_defs_header:-${expected_offline_main}/include/calobase/TowerInfoDefs.h}"
towerinfo_container_header="${towerinfo_container_header:-${expected_offline_main}/include/calobase/TowerInfoContainer.h}"
towerinfo_containerv1_header="${towerinfo_containerv1_header:-${expected_offline_main}/include/calobase/TowerInfoContainerv1.h}"
calo_io_library="${calo_io_library:-${expected_offline_main}/lib/libcalo_io.so}"

for required_path in \
  "$baseline_calo_reco" "$abi_additions_allowlist" \
  "$towerinfo_defs_header" "$towerinfo_container_header" \
  "$towerinfo_containerv1_header" "$calo_io_library"; do
  [[ "$required_path" == /* && "$required_path" != *$'\n'* \
    && "$required_path" != *$'\r'* && "$required_path" != *' '* ]] || \
    die "authority paths must be absolute and contain no whitespace: ${required_path:-<empty>}"
done
for hash_contract in \
  "baseline CaloReco library:${baseline_calo_reco_sha256}" \
  "ABI additions allowlist:${abi_additions_allowlist_sha256}" \
  "TowerInfoDefs header:${towerinfo_defs_header_sha256}" \
  "TowerInfoContainer header:${towerinfo_container_header_sha256}" \
  "TowerInfoContainerv1 header:${towerinfo_containerv1_header_sha256}" \
  "CaloIO library:${calo_io_library_sha256}"; do
  require_sha256 "${hash_contract%%:*} SHA-256" "${hash_contract#*:}"
done
require_regular_file_hash "ABI additions allowlist" \
  "$abi_additions_allowlist" "$abi_additions_allowlist_sha256"
[[ "$abi_additions_allowlist_sha256" == "$ABI_ADDITIONS_ALLOWLIST_SHA256" ]] || \
  die "ABI additions allowlist must match frozen contract SHA-256 ${ABI_ADDITIONS_ALLOWLIST_SHA256}"

builder_sha256="$(sha256_file "${BASH_SOURCE[0]}")"
setup_sha256="$(sha256_file "$setup_script")"

contract_sha256() {
  local action="$1"
  python3 - "$SCHEMA_VERSION" "$action" "$output_root" "$coresoftware_commit" \
    "$coresoftware_tree" "$commit_epoch" "$photon_builder_cc_sha256" \
    "$photon_builder_h_sha256" "$MAPPING_PATCH_ID" "$RELEASE_NAME" \
    "$ROLE_SIZE_GUARD_ID" \
    "$expected_offline_main" "$setup_sha256" "$builder_sha256" "$jobs" \
    "$baseline_calo_reco" "$baseline_calo_reco_sha256" \
    "$abi_additions_allowlist" "$abi_additions_allowlist_sha256" \
    "$ABI_ADDITIONS_PROVENANCE_SHA256" "$ABI_INTENTIONAL_EXCLUSION" \
    "$ABI_ROLE_GUARD_INLINE_ADDITION" \
    "$towerinfo_defs_header" "$towerinfo_defs_header_sha256" \
    "$towerinfo_container_header" "$towerinfo_container_header_sha256" \
    "$towerinfo_containerv1_header" "$towerinfo_containerv1_header_sha256" \
    "$calo_io_library" "$calo_io_library_sha256" <<'PY'
import hashlib
import sys

payload = "\0".join(sys.argv[1:]).encode("utf-8")
print(hashlib.sha256(payload).hexdigest())
PY
}

source_token="the134-ana560-caloreco-source:$(contract_sha256 source)"
build_token="the134-ana560-caloreco-build:$(contract_sha256 build)"

if [[ "$mode" == "plan" ]]; then
  cat <<EOF
THE134_ANA560_CALORECO_BUILD_PLAN
  schema: ${SCHEMA_VERSION}
  release: ${RELEASE_NAME}
  offline_main: ${expected_offline_main}
  pinned_setup_script: ${EXPECTED_SETUP_SCRIPT}
  pinned_offline_main: ${EXPECTED_OFFLINE_MAIN}
  coresoftware_commit: ${coresoftware_commit}
  pinned_coresoftware_commit: ${EXPECTED_CORESOFTWARE_COMMIT}
  coresoftware_tree: ${coresoftware_tree}
  package_export: git archive offline/packages/CaloReco
  photon_builder_cc_sha256: ${photon_builder_cc_sha256}
  photon_builder_h_sha256: ${photon_builder_h_sha256}
  mapping_patch: ${MAPPING_PATCH_ID}
  role_size_guard: ${ROLE_SIZE_GUARD_ID}
  baseline_calo_reco: ${baseline_calo_reco}
  baseline_calo_reco_sha256: ${baseline_calo_reco_sha256}
  abi_additions_allowlist: ${abi_additions_allowlist}
  abi_additions_allowlist_sha256: ${abi_additions_allowlist_sha256}
  abi_additions_provenance_sha256: ${ABI_ADDITIONS_PROVENANCE_SHA256}
  abi_intentional_exclusion: ${ABI_INTENTIONAL_EXCLUSION}
  towerinfo_defs_header_sha256: ${towerinfo_defs_header_sha256}
  towerinfo_container_header_sha256: ${towerinfo_container_header_sha256}
  towerinfo_containerv1_header_sha256: ${towerinfo_containerv1_header_sha256}
  calo_io_library_sha256: ${calo_io_library_sha256}
  library_contract: one full ${EXPECTED_SONAME} provides PhotonClusterBuilder and RawClusterBuilderTopo
  output_root: ${output_root}
  source_token: ${source_token}
  build_token: ${build_token}
PLAN_ONLY_NO_FILES_WRITTEN
EOF
  exit 0
fi

expected_token="$source_token"
if [[ "$mode" == "build" ]]; then
  expected_token="$build_token"
fi
[[ "$provided_token" == "$expected_token" ]] || \
  die "${mode} authorization token differs from the sealed contract"

for tool in git tar python3; do
  command -v "$tool" >/dev/null 2>&1 || die "required tool is unavailable: ${tool}"
done

clean_environment_args=(
  "HOME=${HOME:-/tmp}"
  "USER=${USER:-unknown}"
  "LOGNAME=${LOGNAME:-${USER:-unknown}}"
  "SHELL=/bin/bash"
  "PATH=/usr/bin:/bin:/usr/sbin:/sbin"
  "TMPDIR=${TMPDIR:-/tmp}"
  "LC_ALL=C"
  "LANG=C"
  "TZ=UTC"
)

if [[ "$mode" == "build" ]]; then
  release_probe="$(
    /usr/bin/env -i "${clean_environment_args[@]}" /bin/bash -s -- \
      "$setup_script" "$RELEASE_NAME" <<'BASH'
set -euo pipefail
setup_script="$1"
release_name="$2"
set +e
set +u
# shellcheck disable=SC1090
source "$setup_script" -n "$release_name" >/dev/null
setup_status=$?
set -u
set -e
[[ $setup_status -eq 0 ]]
printf '%s\n' "${OFFLINE_MAIN:-}"
BASH
  )"
  [[ "$release_probe" == "$expected_offline_main" ]] || \
    die "clean setup resolved OFFLINE_MAIN=${release_probe:-<empty>}, expected ${expected_offline_main}"
  expected_offline_main_real="$(
    python3 - "$expected_offline_main" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve(strict=True))
PY
  )" || die "could not resolve exact ana.560 prefix: ${expected_offline_main}"
  authority_contracts=(
    "${baseline_calo_reco}::${expected_offline_main}/lib/${EXPECTED_SONAME}"
    "${towerinfo_defs_header}::${expected_offline_main}/include/calobase/TowerInfoDefs.h"
    "${towerinfo_container_header}::${expected_offline_main}/include/calobase/TowerInfoContainer.h"
    "${towerinfo_containerv1_header}::${expected_offline_main}/include/calobase/TowerInfoContainerv1.h"
    "${calo_io_library}::${expected_offline_main}/lib/libcalo_io.so"
  )
  for authority_contract in "${authority_contracts[@]}"; do
    authority_path="${authority_contract%%::*}"
    expected_authority_path="${authority_contract#*::}"
    authority_real="$(
      python3 - "$authority_path" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve(strict=True))
PY
    )" || die "could not resolve ana.560 authority path: ${authority_path}"
    [[ "$authority_real" == "${expected_offline_main_real}/"* ]] || \
      die "authority path resolves outside exact ana.560 prefix: ${authority_path} -> ${authority_real}"
    expected_authority_real="$(
      python3 - "$expected_authority_path" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve(strict=True))
PY
    )" || die "could not resolve expected ana.560 authority path: ${expected_authority_path}"
    [[ "$authority_real" == "$expected_authority_real" ]] || \
      die "authority path is not the exact ana.560 artifact: ${authority_path}"
  done
  require_file_hash "baseline ana.560 CaloReco library" \
    "$baseline_calo_reco" "$baseline_calo_reco_sha256"
  require_file_hash "ana.560 TowerInfoDefs header" \
    "$towerinfo_defs_header" "$towerinfo_defs_header_sha256"
  require_file_hash "ana.560 TowerInfoContainer header" \
    "$towerinfo_container_header" "$towerinfo_container_header_sha256"
  require_file_hash "ana.560 TowerInfoContainerv1 header" \
    "$towerinfo_containerv1_header" "$towerinfo_containerv1_header_sha256"
  require_file_hash "ana.560 CaloIO library" \
    "$calo_io_library" "$calo_io_library_sha256"
fi

mkdir -p "$output_root/source"
source_package="${output_root}/source/offline/packages/CaloReco"

# Export only the committed package.  No mutable checkout file can enter the
# build through this source path.
git "${git_safe[@]}" -C "$coresoftware_repo_real" archive --format=tar \
  "$coresoftware_commit" offline/packages/CaloReco \
  | tar -xf - -C "${output_root}/source"
[[ -d "$source_package" ]] || die "git archive did not materialize the CaloReco package"

archived_photon_cc_sha256="$(sha256_file "${source_package}/PhotonClusterBuilder.cc")"
archived_photon_h_sha256="$(sha256_file "${source_package}/PhotonClusterBuilder.h")"
topo_before_sha256="$(sha256_file "${source_package}/RawClusterBuilderTopo.cc")"

install -m 0644 "$photon_builder_cc" "${source_package}/PhotonClusterBuilder.cc"
install -m 0644 "$photon_builder_h" "${source_package}/PhotonClusterBuilder.h"
require_file_hash "staged PhotonClusterBuilder.cc" \
  "${source_package}/PhotonClusterBuilder.cc" "$photon_builder_cc_sha256"
require_file_hash "staged PhotonClusterBuilder.h" \
  "${source_package}/PhotonClusterBuilder.h" "$photon_builder_h_sha256"

python3 - "${source_package}/RawClusterBuilderTopo.cc" \
  "${source_package}/THE134TowerInfoContractV1.h" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
contract_header = Path(sys.argv[2])
text = path.read_text(encoding="utf-8")

contract_header.write_text(
    """#ifndef THE134_TOWERINFO_CONTRACT_V1_H
#define THE134_TOWERINFO_CONTRACT_V1_H

#include <cstddef>

namespace the134::caloreco
{
enum class TowerInfoContractResult
{
  ACCEPT_TYPED,
  ACCEPT_LEGACY_INVALID,
  REJECT
};

inline TowerInfoContractResult evaluate_towerinfo_role_and_size(
    const int actual_detector,
    const std::size_t actual_size,
    const int expected_detector,
    const std::size_t expected_size,
    const int invalid_detector)
{
  if (actual_size != expected_size)
  {
    return TowerInfoContractResult::REJECT;
  }
  if (actual_detector == expected_detector)
  {
    return TowerInfoContractResult::ACCEPT_TYPED;
  }
  if (actual_detector == invalid_detector)
  {
    return TowerInfoContractResult::ACCEPT_LEGACY_INVALID;
  }
  return TowerInfoContractResult::REJECT;
}
}  // namespace the134::caloreco

#endif
""",
    encoding="utf-8",
)

local_include = '#include "THE134TowerInfoContractV1.h"\n'
if local_include in text:
    raise SystemExit("RawClusterBuilderTopo already contains the THE-134 contract include")
text = local_include + text

include_anchor = "#include <calobase/TowerInfoContainer.h>  // for TowerInfoContainer\n"
include_line = "#include <calobase/TowerInfoDefs.h>       // detector-explicit channel maps\n"
if text.count(include_anchor) != 1:
    raise SystemExit("RawClusterBuilderTopo TowerInfoContainer include anchor count differs")
if include_line in text:
    raise SystemExit("RawClusterBuilderTopo already contains the THE-134 mapping include")
text = text.replace(include_anchor, include_anchor + include_line, 1)

cstddef_anchor = "#include <algorithm>\n"
cstddef_line = "#include <cstddef>\n"
if text.count(cstddef_anchor) != 1:
    raise SystemExit("RawClusterBuilderTopo algorithm include anchor count differs")
if cstddef_line in text:
    raise SystemExit("RawClusterBuilderTopo already contains the THE-134 cstddef include")
text = text.replace(cstddef_anchor, cstddef_anchor + cstddef_line, 1)

namespace_anchor = """bool sort_by_pair_second(const std::pair<int, float> &a, const std::pair<int, float> &b)
{
  return (a.second > b.second);
}
"""
role_guard = """
bool validate_towerinfo_role_and_size(
    const TowerInfoContainer *container,
    const TowerInfoContainer::DETECTOR expected_detector,
    const std::size_t expected_size,
    const char *role,
    bool &legacy_warning_emitted)
{
  const auto actual_detector = container->get_detectorid();
  const auto actual_size = container->size();
  const auto result = the134::caloreco::evaluate_towerinfo_role_and_size(
      static_cast<int>(actual_detector),
      actual_size,
      static_cast<int>(expected_detector),
      expected_size,
      static_cast<int>(TowerInfoContainer::DETECTOR_INVALID));
  if (result == the134::caloreco::TowerInfoContractResult::ACCEPT_TYPED)
  {
    return true;
  }
  if (result == the134::caloreco::TowerInfoContractResult::ACCEPT_LEGACY_INVALID)
  {
    if (!legacy_warning_emitted)
    {
      std::cerr << "THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_LEGACY_INVALID_ACCEPTED"
                << " role=" << role
                << " exact_size=" << actual_size
                << " expected_detector=" << static_cast<int>(expected_detector)
                << " diagnostic_frequency=once_per_role_process"
                << std::endl;
      legacy_warning_emitted = true;
    }
    return true;
  }
  std::cerr << "THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_CONTRACT_FAIL"
            << " role=" << role
            << " expected_size=" << expected_size
            << " actual_size=" << actual_size
            << " expected_detector=" << static_cast<int>(expected_detector)
            << " actual_detector=" << static_cast<int>(actual_detector)
            << std::endl;
  return false;
}
"""
if text.count(namespace_anchor) != 1:
    raise SystemExit("RawClusterBuilderTopo anonymous-namespace anchor count differs")
if "validate_towerinfo_role_and_size" in text:
    raise SystemExit("RawClusterBuilderTopo already contains the THE-134 role/size guard")
text = text.replace(namespace_anchor, namespace_anchor + role_guard, 1)

role_call_anchor = """  if (!towerinfosOH)
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERINFO_CALIB_HCALOUT does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
"""
role_calls = """
  static bool cemc_legacy_warning_emitted = false;
  static bool hcalin_legacy_warning_emitted = false;
  static bool hcalout_legacy_warning_emitted = false;
  if (!validate_towerinfo_role_and_size(
          towerinfosEM, TowerInfoContainer::EMCAL, 24576U,
          "TOWERINFO_CALIB_CEMC", cemc_legacy_warning_emitted) ||
      !validate_towerinfo_role_and_size(
          towerinfosIH, TowerInfoContainer::HCAL, 1536U,
          "TOWERINFO_CALIB_HCALIN", hcalin_legacy_warning_emitted) ||
      !validate_towerinfo_role_and_size(
          towerinfosOH, TowerInfoContainer::HCAL, 1536U,
          "TOWERINFO_CALIB_HCALOUT", hcalout_legacy_warning_emitted))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
"""
if text.count(role_call_anchor) != 1:
    raise SystemExit("RawClusterBuilderTopo HCALOUT null-check anchor count differs")
text = text.replace(role_call_anchor, role_call_anchor + role_calls, 1)

replacements = (
    (
        "      unsigned int towerinfo_key = towerinfosEM->encode_key(iEM);\n"
        "      int ti_ieta = towerinfosEM->getTowerEtaBin(towerinfo_key);\n"
        "      int ti_iphi = towerinfosEM->getTowerPhiBin(towerinfo_key);\n",
        "      const unsigned int towerinfo_key = TowerInfoDefs::encode_emcal(iEM);\n"
        "      const int ti_ieta = static_cast<int>(TowerInfoDefs::getCaloTowerEtaBin(towerinfo_key));\n"
        "      const int ti_iphi = static_cast<int>(TowerInfoDefs::getCaloTowerPhiBin(towerinfo_key));\n",
    ),
    (
        "      unsigned int towerinfo_key = towerinfosIH->encode_key(iIH);\n"
        "      int ti_ieta = towerinfosIH->getTowerEtaBin(towerinfo_key);\n"
        "      int ti_iphi = towerinfosIH->getTowerPhiBin(towerinfo_key);\n",
        "      const unsigned int towerinfo_key = TowerInfoDefs::encode_hcal(iIH);\n"
        "      const int ti_ieta = static_cast<int>(TowerInfoDefs::getCaloTowerEtaBin(towerinfo_key));\n"
        "      const int ti_iphi = static_cast<int>(TowerInfoDefs::getCaloTowerPhiBin(towerinfo_key));\n",
    ),
    (
        "      unsigned int towerinfo_key = towerinfosOH->encode_key(iOH);\n"
        "      int ti_ieta = towerinfosOH->getTowerEtaBin(towerinfo_key);\n"
        "      int ti_iphi = towerinfosOH->getTowerPhiBin(towerinfo_key);\n",
        "      const unsigned int towerinfo_key = TowerInfoDefs::encode_hcal(iOH);\n"
        "      const int ti_ieta = static_cast<int>(TowerInfoDefs::getCaloTowerEtaBin(towerinfo_key));\n"
        "      const int ti_iphi = static_cast<int>(TowerInfoDefs::getCaloTowerPhiBin(towerinfo_key));\n",
    ),
)

for old, new in replacements:
    if text.count(old) != 1:
        raise SystemExit("RawClusterBuilderTopo channel-map source anchor count differs")
    text = text.replace(old, new, 1)

for old, _ in replacements:
    if old in text:
        raise SystemExit("RawClusterBuilderTopo retains a container-selected channel map")

path.write_text(text, encoding="utf-8")
PY

topo_after_sha256="$(sha256_file "${source_package}/RawClusterBuilderTopo.cc")"
[[ "$topo_after_sha256" != "$topo_before_sha256" ]] || \
  die "RawClusterBuilderTopo mapping patch did not change the source"

source_manifest="${output_root}/source_manifest.json"
python3 - "$source_package" "$source_manifest" "$SOURCE_SCHEMA" \
  "$coresoftware_commit" "$coresoftware_tree" "$commit_epoch" \
  "$archived_photon_cc_sha256" "$archived_photon_h_sha256" \
  "$photon_builder_cc_sha256" "$photon_builder_h_sha256" \
  "$MAPPING_PATCH_ID" "$ROLE_SIZE_GUARD_ID" \
  "$topo_before_sha256" "$topo_after_sha256" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

(
    package_arg,
    manifest_arg,
    schema,
    commit,
    tree,
    epoch,
    archived_cc,
    archived_h,
    overlay_cc,
    overlay_h,
    mapping_id,
    role_size_guard_id,
    topo_before,
    topo_after,
) = sys.argv[1:]
package = Path(package_arg)
manifest = Path(manifest_arg)

files = []
aggregate = hashlib.sha256()
for path in sorted(p for p in package.rglob("*") if p.is_file()):
    relative = path.relative_to(package).as_posix()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    files.append({"path": relative, "sha256": digest})
    aggregate.update(relative.encode("utf-8"))
    aggregate.update(b"\0")
    aggregate.update(digest.encode("ascii"))
    aggregate.update(b"\0")

payload = {
    "coresoftware": {
        "commit": commit,
        "commit_epoch": int(epoch),
        "package": "offline/packages/CaloReco",
        "tree": tree,
    },
    "files": files,
    "mapping_patch": {
        "changed_scientific_controls": [],
        "detector_maps": {
            "CEMC": "TowerInfoDefs::encode_emcal(channel)",
            "HCALIN": "TowerInfoDefs::encode_hcal(channel)",
            "HCALOUT": "TowerInfoDefs::encode_hcal(channel)",
        },
        "id": mapping_id,
        "replacement_counts": {"CEMC": 1, "HCALIN": 1, "HCALOUT": 1},
        "semantic_equivalence": (
            "For a correctly typed TowerInfoContainer these are the same "
            "TowerInfoDefs conversions selected by encode_key. Energy, tower "
            "quality, geometry, thresholds, ordering, and clustering are untouched."
        ),
        "sha256_after": topo_after,
        "sha256_before": topo_before,
    },
    "role_size_guard": {
        "accepted": {
            "CEMC": {
                "detector": "EMCAL",
                "legacy_detector_invalid": "accepted only at exact size",
                "size": 24576,
            },
            "HCALIN": {
                "detector": "HCAL",
                "legacy_detector_invalid": "accepted only at exact size",
                "size": 1536,
            },
            "HCALOUT": {
                "detector": "HCAL",
                "legacy_detector_invalid": "accepted only at exact size",
                "size": 1536,
            },
        },
        "diagnostics": {
            "frequency": "once_per_role_process",
            "failure_marker": (
                "THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_CONTRACT_FAIL"
            ),
            "legacy_invalid_marker": (
                "THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_LEGACY_INVALID_ACCEPTED"
            ),
        },
        "failure_return": "Fun4AllReturnCodes::ABORTRUN",
        "id": role_size_guard_id,
        "predicate_header": {
            "path": "THE134TowerInfoContractV1.h",
            "sha256": hashlib.sha256(
                (package / "THE134TowerInfoContractV1.h").read_bytes()
            ).hexdigest(),
        },
        "valid_path_effect": (
            "none; validation precedes geometry lookup and tower arithmetic"
        ),
        "wrong_valid_detector": "rejected",
    },
    "overlay": {
        "PhotonClusterBuilder.cc": {
            "archived_sha256": archived_cc,
            "staged_sha256": overlay_cc,
        },
        "PhotonClusterBuilder.h": {
            "archived_sha256": archived_h,
            "staged_sha256": overlay_h,
        },
    },
    "package_tree_sha256": aggregate.hexdigest(),
    "schema": schema,
}
manifest.write_text(
    json.dumps(payload, indent=2, sort_keys=True, separators=(",", ": ")) + "\n",
    encoding="utf-8",
)
PY
source_manifest_sha256="$(sha256_file "$source_manifest")"

if [[ "$mode" == "source" ]]; then
  printf '%s\n' \
    "THE134_ANA560_CALORECO_SOURCE_ONLY_NOT_RUNTIME_AUTHORITY" \
    >"${output_root}/SOURCE_ONLY_NOT_RUNTIME_AUTHORITY"
  chmod -R a-w "$output_root"
  cat <<EOF
THE134_ANA560_CALORECO_SOURCE_STAGE_PASS
  output_root: ${output_root}
  source_manifest: ${source_manifest}
  source_manifest_sha256: ${source_manifest_sha256}
  runtime_authority: false
EOF
  exit 0
fi

mkdir -p "${output_root}/work" "${output_root}/build" \
  "${output_root}/install" "${output_root}/logs" "${output_root}/evidence"
cp -a "$source_package" "${output_root}/work/CaloReco"
work_source="${output_root}/work/CaloReco"
build_dir="${output_root}/build"
install_root="${output_root}/install"
build_log="${output_root}/logs/build.log"

/usr/bin/env -i "${clean_environment_args[@]}" \
  "SOURCE_DATE_EPOCH=${commit_epoch}" "ZERO_AR_DATE=1" \
  /bin/bash -s -- "$setup_script" "$RELEASE_NAME" "$expected_offline_main" \
  "$work_source" "$build_dir" "$install_root" "$jobs" "$output_root" \
  >"$build_log" 2>&1 <<'BASH'
set -euo pipefail
setup_script="$1"
release_name="$2"
expected_offline_main="$3"
work_source="$4"
build_dir="$5"
install_root="$6"
jobs="$7"
output_root="$8"

set +e
set +u
# shellcheck disable=SC1090
source "$setup_script" -n "$release_name"
setup_status=$?
set -u
set -e
[[ $setup_status -eq 0 ]]
[[ "${OFFLINE_MAIN:-}" == "$expected_offline_main" ]]
for tool in make aclocal automake autoconf libtoolize root-config; do
  command -v "$tool" >/dev/null
done

export LC_ALL=C
export LANG=C
export TZ=UTC
export CFLAGS="${CFLAGS:-} -ffile-prefix-map=${output_root}=/the134-caloreco -fdebug-prefix-map=${output_root}=/the134-caloreco"
export CXXFLAGS="${CXXFLAGS:-} -ffile-prefix-map=${output_root}=/the134-caloreco -fdebug-prefix-map=${output_root}=/the134-caloreco"

cd "$build_dir"
/bin/bash "${work_source}/autogen.sh" --prefix="$install_root"
make -j "$jobs"
make install
BASH

runtime_link="${install_root}/lib/${EXPECTED_SONAME}"
runtime_api_link="${install_root}/lib/libcalo_reco.so"
[[ -L "$runtime_link" ]] || \
  die "build did not install ${runtime_link}"
[[ -L "$runtime_api_link" ]] || \
  die "build did not install ${runtime_api_link}"
runtime_library="$(python3 - "$runtime_link" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve(strict=True))
PY
)"
install_root_real="$(python3 - "$install_root" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve(strict=True))
PY
)"
[[ "$runtime_library" == "${install_root_real}/lib/"* && -f "$runtime_library" \
  && -s "$runtime_library" ]] || die "resolved CaloReco library escaped the install root"
runtime_api_library="$(python3 - "$runtime_api_link" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).resolve(strict=True))
PY
)"
[[ "$runtime_api_library" == "$runtime_library" ]] || \
  die "libcalo_reco.so and ${EXPECTED_SONAME} resolve to different providers"
require_file_hash "installed PhotonClusterBuilder header" \
  "${install_root}/include/caloreco/PhotonClusterBuilder.h" \
  "$photon_builder_h_sha256"

validation_log="${output_root}/logs/validation.log"
/usr/bin/env -i "${clean_environment_args[@]}" \
  /bin/bash -s -- "$setup_script" "$RELEASE_NAME" "$expected_offline_main" \
  "$install_root" "$runtime_library" "$baseline_calo_reco" \
  "$abi_additions_allowlist" "${output_root}/evidence" "$EXPECTED_SONAME" \
  "$work_source" "$calo_io_library" \
  >"$validation_log" 2>&1 <<'BASH'
set -euo pipefail
setup_script="$1"
release_name="$2"
expected_offline_main="$3"
install_root="$4"
runtime_library="$5"
baseline_library="$6"
abi_additions_allowlist="$7"
evidence_root="$8"
expected_soname="$9"
work_source="${10}"
calo_io_library="${11}"

set +e
set +u
# shellcheck disable=SC1090
source "$setup_script" -n "$release_name"
setup_status=$?
set -u
set -e
[[ $setup_status -eq 0 ]]
[[ "${OFFLINE_MAIN:-}" == "$expected_offline_main" ]]
export LD_LIBRARY_PATH="${install_root}/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="${install_root}/include:${ROOT_INCLUDE_PATH:-}"
export LC_ALL=C
export LANG=C
export TZ=UTC

for tool in readelf ldd nm c++filt root root-config c++ awk comm cmp sort; do
  command -v "$tool" >/dev/null
done

normalize_symbols() {
  nm -D --defined-only --format=posix "$1" \
    | awk 'NF >= 2 { print $1 " " $2 }' \
    | LC_ALL=C sort -u
}

normalize_symbols "$baseline_library" \
  >"${evidence_root}/baseline_dynamic_symbols.normalized.txt"
normalize_symbols "$runtime_library" \
  >"${evidence_root}/candidate_dynamic_symbols.normalized.txt"
awk '
  {
    sub(/[[:space:]]*#.*/, "")
    gsub(/^[[:space:]]+|[[:space:]]+$/, "")
    if (length($0) > 0) print
  }
' "$abi_additions_allowlist" | LC_ALL=C sort -u \
  >"${evidence_root}/abi_additions_allowlist.normalized.txt"

comm -23 \
  "${evidence_root}/baseline_dynamic_symbols.normalized.txt" \
  "${evidence_root}/candidate_dynamic_symbols.normalized.txt" \
  >"${evidence_root}/abi_removed_symbols.txt"
comm -13 \
  "${evidence_root}/baseline_dynamic_symbols.normalized.txt" \
  "${evidence_root}/candidate_dynamic_symbols.normalized.txt" \
  >"${evidence_root}/abi_added_symbols.txt"
[[ ! -s "${evidence_root}/abi_removed_symbols.txt" ]]
cmp -s \
  "${evidence_root}/abi_added_symbols.txt" \
  "${evidence_root}/abi_additions_allowlist.normalized.txt"

readelf -d "$baseline_library" >"${evidence_root}/baseline_readelf_dynamic.txt"
readelf -d "$runtime_library" >"${evidence_root}/candidate_readelf_dynamic.txt"
normalize_dynsym_metadata() {
  readelf --wide --dyn-syms "$1" \
    | awk '$1 ~ /^[0-9]+:$/ && $7 != "UND" {
        print $4 "|" $5 "|" $6 "|" $8
      }' \
    | LC_ALL=C sort -u
}
normalize_version_inventory() {
  readelf --wide --version-info "$1" \
    | sed -n 's/.*Name: \([^ ]*\).*/\1/p' \
    | LC_ALL=C sort -u
}
normalize_dynsym_metadata "$baseline_library" \
  >"${evidence_root}/baseline_dynsym_metadata.normalized.txt"
normalize_dynsym_metadata "$runtime_library" \
  >"${evidence_root}/candidate_dynsym_metadata.normalized.txt"
comm -23 \
  "${evidence_root}/baseline_dynsym_metadata.normalized.txt" \
  "${evidence_root}/candidate_dynsym_metadata.normalized.txt" \
  >"${evidence_root}/abi_removed_dynsym_metadata.txt"
[[ ! -s "${evidence_root}/abi_removed_dynsym_metadata.txt" ]]
normalize_version_inventory "$baseline_library" \
  >"${evidence_root}/baseline_version_inventory.txt"
normalize_version_inventory "$runtime_library" \
  >"${evidence_root}/candidate_version_inventory.txt"
cmp -s "${evidence_root}/baseline_version_inventory.txt" \
  "${evidence_root}/candidate_version_inventory.txt"
ldd "$runtime_library" >"${evidence_root}/ldd.txt"
nm -D --defined-only "$runtime_library" | c++filt \
  >"${evidence_root}/defined_symbols.txt"
c++ --version | head -n 1 >"${evidence_root}/compiler_version.txt"

sed -n 's/.*(SONAME).*\[\([^]]*\)\].*/\1/p' \
  "${evidence_root}/baseline_readelf_dynamic.txt" \
  >"${evidence_root}/baseline_soname.txt"
sed -n 's/.*(SONAME).*\[\([^]]*\)\].*/\1/p' \
  "${evidence_root}/candidate_readelf_dynamic.txt" \
  >"${evidence_root}/candidate_soname.txt"
sed -n 's/.*(NEEDED).*\[\([^]]*\)\].*/\1/p' \
  "${evidence_root}/baseline_readelf_dynamic.txt" | LC_ALL=C sort -u \
  >"${evidence_root}/baseline_needed.txt"
sed -n 's/.*(NEEDED).*\[\([^]]*\)\].*/\1/p' \
  "${evidence_root}/candidate_readelf_dynamic.txt" | LC_ALL=C sort -u \
  >"${evidence_root}/candidate_needed.txt"

[[ "$(wc -l <"${evidence_root}/baseline_soname.txt")" == "1" ]]
[[ "$(wc -l <"${evidence_root}/candidate_soname.txt")" == "1" ]]
[[ "$(cat "${evidence_root}/baseline_soname.txt")" == "$expected_soname" ]]
[[ "$(cat "${evidence_root}/candidate_soname.txt")" == "$expected_soname" ]]
cmp -s "${evidence_root}/baseline_needed.txt" \
  "${evidence_root}/candidate_needed.txt"
if grep -Eq '\((RPATH|RUNPATH)\)' \
    "${evidence_root}/candidate_readelf_dynamic.txt"; then
  exit 20
fi

find "$(dirname "$baseline_library")" -maxdepth 1 \
  \( -name 'libcalo_reco*.pcm' -o -name 'libcalo_reco*.rootmap' \) \
  -print | while IFS= read -r artifact; do basename "$artifact"; done \
  | LC_ALL=C sort -u >"${evidence_root}/baseline_dictionary_inventory.txt"
find "${install_root}/lib" -maxdepth 1 \
  \( -name 'libcalo_reco*.pcm' -o -name 'libcalo_reco*.rootmap' \) \
  -print | while IFS= read -r artifact; do basename "$artifact"; done \
  | LC_ALL=C sort -u >"${evidence_root}/candidate_dictionary_inventory.txt"
cmp -s \
  "${evidence_root}/baseline_dictionary_inventory.txt" \
  "${evidence_root}/candidate_dictionary_inventory.txt"

grep -Eq "\\(SONAME\\).*(\\[${expected_soname}\\])" \
  "${evidence_root}/candidate_readelf_dynamic.txt"
if grep -F 'not found' "${evidence_root}/ldd.txt" >/dev/null; then
  exit 21
fi
python3 - "${evidence_root}/ldd.txt" \
  "${evidence_root}/resolved_dependency_paths.txt" \
  "$install_root" "$expected_offline_main" <<'PY'
from pathlib import Path
import re
import sys

ldd_path, output_path, install_arg, offline_arg = sys.argv[1:]
install = Path(install_arg).resolve(strict=True)
offline = Path(offline_arg).resolve(strict=True)
offline_parts = offline.parts
platform_root = None
for index in range(len(offline_parts) - 1):
    if offline_parts[index : index + 2] == ("release", "release_ana"):
        platform_root = Path(*offline_parts[:index]).resolve(strict=True)
        break
if platform_root is None:
    raise SystemExit("could not derive immutable platform root from OFFLINE_MAIN")
system_roots = tuple(
    path.resolve(strict=True)
    for path in (Path("/lib"), Path("/lib64"), Path("/usr/lib"), Path("/usr/lib64"))
    if path.exists()
)

def below(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False

resolved = set()
for line in Path(ldd_path).read_text(encoding="utf-8").splitlines():
    match = re.search(r"=>\s+(/[^\s(]+)", line)
    if not match:
        match = re.match(r"\s*(/[^\s(]+)", line)
    if not match:
        continue
    path = Path(match.group(1)).resolve(strict=True)
    resolved.add(path)
    if not (
        below(path, install)
        or below(path, offline)
        or below(path, platform_root)
        or any(below(path, root) for root in system_roots)
    ):
        raise SystemExit(f"dependency outside immutable allowlist: {path}")

Path(output_path).write_text(
    "".join(f"{path}\n" for path in sorted(resolved, key=str)),
    encoding="utf-8",
)
PY

photon_count="$(
  grep -F -c 'PhotonClusterBuilder::process_event(PHCompositeNode*)' \
    "${evidence_root}/defined_symbols.txt" || true
)"
topo_count="$(
  grep -F -c 'RawClusterBuilderTopo::process_event(PHCompositeNode*)' \
    "${evidence_root}/defined_symbols.txt" || true
)"
[[ "$photon_count" == "1" && "$topo_count" == "1" ]]

find "${install_root}/lib" -maxdepth 1 \( -type f -o -type l \) \
  \( -name '*PhotonClusterBuilder*.so*' -o -name '*photon*override*.so*' \) \
  -print >"${evidence_root}/forbidden_providers.txt"
[[ ! -s "${evidence_root}/forbidden_providers.txt" ]]

cat >"${evidence_root}/single_provider_probe.cc" <<'CPP'
#include <dlfcn.h>
#include <link.h>

#include <climits>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace
{
std::string canonical(const char *path)
{
  char resolved[PATH_MAX] = {};
  if (!path || !*path || !realpath(path, resolved))
  {
    return {};
  }
  return resolved;
}

bool is_calo_reco(const char *path)
{
  if (!path || !*path)
  {
    return false;
  }
  const char *base = std::strrchr(path, '/');
  base = base ? base + 1 : path;
  constexpr const char *prefix = "libcalo_reco.so";
  const auto length = std::strlen(prefix);
  return std::strncmp(base, prefix, length) == 0 &&
         (base[length] == '\0' || base[length] == '.');
}

int collect(struct dl_phdr_info *info, std::size_t, void *payload)
{
  auto *providers = static_cast<std::vector<std::string> *>(payload);
  if (is_calo_reco(info->dlpi_name))
  {
    providers->push_back(canonical(info->dlpi_name));
  }
  return 0;
}

std::vector<std::string> providers()
{
  std::vector<std::string> result;
  dl_iterate_phdr(collect, &result);
  return result;
}
}  // namespace

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    return 2;
  }
  const std::string expected = canonical(argv[1]);
  if (expected.empty() || !providers().empty())
  {
    return 3;
  }
  void *handle = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if (!handle)
  {
    std::cerr << dlerror() << '\n';
    return 4;
  }
  const auto loaded = providers();
  if (loaded.size() != 1 || loaded.front() != expected)
  {
    return 5;
  }
  const char *symbols[] = {
      "_ZN20PhotonClusterBuilder13process_eventEP15PHCompositeNode",
      "_ZN21RawClusterBuilderTopo13process_eventEP15PHCompositeNode",
  };
  for (const char *symbol_name : symbols)
  {
    dlerror();
    void *symbol = dlsym(handle, symbol_name);
    if (!symbol || dlerror())
    {
      return 6;
    }
    Dl_info owner = {};
    if (!dladdr(symbol, &owner) || canonical(owner.dli_fname) != expected)
    {
      return 7;
    }
  }
  std::cout << "THE134_SINGLE_PROVIDER_PROBE_PASS provider=" << expected
            << " provider_count=1" << std::endl;
  return 0;
}
CPP
c++ -std=c++17 "${evidence_root}/single_provider_probe.cc" -ldl \
  -o "${evidence_root}/single_provider_probe"
"${evidence_root}/single_provider_probe" "$runtime_library" \
  >"${evidence_root}/single_provider_probe.txt"
grep -Fq "THE134_SINGLE_PROVIDER_PROBE_PASS provider=${runtime_library} provider_count=1" \
  "${evidence_root}/single_provider_probe.txt"

cat >"${evidence_root}/towerinfo_mapping_equivalence.cc" <<'CPP'
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoDefs.h>

#include <iostream>

int main()
{
  TowerInfoContainerv1 emcal(TowerInfoContainer::EMCAL);
  TowerInfoContainerv1 hcal(TowerInfoContainer::HCAL);
  if (emcal.size() != 24576U || hcal.size() != 1536U)
  {
    return 2;
  }
  for (unsigned int channel = 0; channel < 24576U; ++channel)
  {
    const auto container_key = emcal.encode_key(channel);
    const auto explicit_key = TowerInfoDefs::encode_emcal(channel);
    if (container_key != explicit_key ||
        emcal.decode_key(explicit_key) != channel ||
        emcal.getTowerEtaBin(container_key) !=
            TowerInfoDefs::getCaloTowerEtaBin(explicit_key) ||
        emcal.getTowerPhiBin(container_key) !=
            TowerInfoDefs::getCaloTowerPhiBin(explicit_key))
    {
      std::cerr << "EMCAL mismatch at channel " << channel << std::endl;
      return 3;
    }
  }
  for (unsigned int channel = 0; channel < 1536U; ++channel)
  {
    const auto container_key = hcal.encode_key(channel);
    const auto explicit_key = TowerInfoDefs::encode_hcal(channel);
    if (container_key != explicit_key ||
        hcal.decode_key(explicit_key) != channel ||
        hcal.getTowerEtaBin(container_key) !=
            TowerInfoDefs::getCaloTowerEtaBin(explicit_key) ||
        hcal.getTowerPhiBin(container_key) !=
            TowerInfoDefs::getCaloTowerPhiBin(explicit_key))
    {
      std::cerr << "HCAL mismatch at channel " << channel << std::endl;
      return 4;
    }
  }
  std::cout << "THE134_TOWERINFO_MAPPING_EQUIVALENCE_PASS"
            << " EMCAL=24576 HCAL=1536" << std::endl;
  return 0;
}
CPP
root_cflags=()
root_libs=()
read -r -a root_cflags <<<"$(root-config --cflags)" || true
read -r -a root_libs <<<"$(root-config --libs)" || true
c++ -std=c++17 ${root_cflags[@]+"${root_cflags[@]}"} \
  -I"${expected_offline_main}/include" \
  "${evidence_root}/towerinfo_mapping_equivalence.cc" \
  -L"${expected_offline_main}/lib" -Wl,-rpath,"${expected_offline_main}/lib" \
  "$calo_io_library" -lphool ${root_libs[@]+"${root_libs[@]}"} \
  -o "${evidence_root}/towerinfo_mapping_equivalence"
"${evidence_root}/towerinfo_mapping_equivalence" \
  >"${evidence_root}/towerinfo_mapping_equivalence.txt"
grep -Fq 'THE134_TOWERINFO_MAPPING_EQUIVALENCE_PASS EMCAL=24576 HCAL=1536' \
  "${evidence_root}/towerinfo_mapping_equivalence.txt"

cat >"${evidence_root}/towerinfo_role_size_matrix.cc" <<'CPP'
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>

#include "THE134TowerInfoContractV1.h"

#include <cstddef>
#include <iostream>

namespace
{
using Result = the134::caloreco::TowerInfoContractResult;

bool require(
    const char *name,
    const int actual_detector,
    const std::size_t actual_size,
    const int expected_detector,
    const std::size_t expected_size,
    const int invalid_detector,
    const Result expected)
{
  const auto actual = the134::caloreco::evaluate_towerinfo_role_and_size(
      actual_detector, actual_size, expected_detector, expected_size,
      invalid_detector);
  if (actual != expected)
  {
    std::cerr << "role/size matrix mismatch: " << name << std::endl;
    return false;
  }
  return true;
}
}  // namespace

int main()
{
  TowerInfoContainerv1 emcal(TowerInfoContainer::EMCAL);
  TowerInfoContainerv1 hcal(TowerInfoContainer::HCAL);
  const int emcal_id = static_cast<int>(emcal.get_detectorid());
  const int hcal_id = static_cast<int>(hcal.get_detectorid());
  const int invalid_id =
      static_cast<int>(TowerInfoContainer::DETECTOR_INVALID);
  if (emcal.size() != 24576U || hcal.size() != 1536U)
  {
    return 2;
  }
  const bool pass =
      require("typed_emcal", emcal_id, emcal.size(), emcal_id, 24576U,
              invalid_id, Result::ACCEPT_TYPED) &&
      require("typed_hcal", hcal_id, hcal.size(), hcal_id, 1536U,
              invalid_id, Result::ACCEPT_TYPED) &&
      require("wrong_valid_hcal_for_emcal", hcal_id, 24576U, emcal_id,
              24576U, invalid_id, Result::REJECT) &&
      require("wrong_valid_emcal_for_hcal", emcal_id, 1536U, hcal_id,
              1536U, invalid_id, Result::REJECT) &&
      require("legacy_invalid_emcal", invalid_id, 24576U, emcal_id, 24576U,
              invalid_id, Result::ACCEPT_LEGACY_INVALID) &&
      require("legacy_invalid_hcal", invalid_id, 1536U, hcal_id, 1536U,
              invalid_id, Result::ACCEPT_LEGACY_INVALID) &&
      require("wrong_size_emcal", emcal_id, 24575U, emcal_id, 24576U,
              invalid_id, Result::REJECT) &&
      require("wrong_size_legacy", invalid_id, 1535U, hcal_id, 1536U,
              invalid_id, Result::REJECT);
  if (!pass)
  {
    return 3;
  }
  std::cout << "THE134_TOWERINFO_ROLE_SIZE_MATRIX_PASS"
            << " typed=2 wrong_valid=2 legacy_exact=2 wrong_size=2"
            << std::endl;
  return 0;
}
CPP
c++ -std=c++17 ${root_cflags[@]+"${root_cflags[@]}"} \
  -I"${expected_offline_main}/include" -I"$work_source" \
  "${evidence_root}/towerinfo_role_size_matrix.cc" \
  -L"${expected_offline_main}/lib" -Wl,-rpath,"${expected_offline_main}/lib" \
  "$calo_io_library" -lphool ${root_libs[@]+"${root_libs[@]}"} \
  -o "${evidence_root}/towerinfo_role_size_matrix"
"${evidence_root}/towerinfo_role_size_matrix" \
  >"${evidence_root}/towerinfo_role_size_matrix.txt"
grep -Fq \
  'THE134_TOWERINFO_ROLE_SIZE_MATRIX_PASS typed=2 wrong_valid=2 legacy_exact=2 wrong_size=2' \
  "${evidence_root}/towerinfo_role_size_matrix.txt"

export THE134_CALORECO_LIBRARY="$runtime_library"
cat >"${evidence_root}/root_load_smoke.C" <<'ROOT'
#include <TSystem.h>

#include <dlfcn.h>
#include <link.h>
#include <limits.h>
#include <stdlib.h>

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace
{
std::string root_canonical(const char *path)
{
  char resolved[PATH_MAX] = {};
  if (!path || !*path || !realpath(path, resolved))
  {
    return {};
  }
  return resolved;
}

bool root_is_calo_reco(const char *path)
{
  if (!path || !*path)
  {
    return false;
  }
  const char *base = std::strrchr(path, '/');
  base = base ? base + 1 : path;
  constexpr const char *prefix = "libcalo_reco.so";
  const auto length = std::strlen(prefix);
  return std::strncmp(base, prefix, length) == 0 &&
         (base[length] == '\0' || base[length] == '.');
}

int root_collect_provider(struct dl_phdr_info *info, std::size_t, void *payload)
{
  auto *providers = static_cast<std::vector<std::string> *>(payload);
  if (root_is_calo_reco(info->dlpi_name))
  {
    providers->push_back(root_canonical(info->dlpi_name));
  }
  return 0;
}

std::vector<std::string> root_providers()
{
  std::vector<std::string> providers;
  dl_iterate_phdr(root_collect_provider, &providers);
  return providers;
}
}  // namespace

void root_load_smoke()
{
  const char *library = gSystem->Getenv("THE134_CALORECO_LIBRARY");
  if (!library || !*library)
  {
    gSystem->Exit(30);
  }
  if (!root_providers().empty())
  {
    gSystem->Exit(29);
  }
  int rc = gSystem->Load(library);
  if (rc != 0)
  {
    gSystem->Exit(31);
  }
  void *handle = dlopen(library, RTLD_NOW | RTLD_NOLOAD);
  void *symbol = handle
      ? dlsym(handle, "_ZN20PhotonClusterBuilder13process_eventEP15PHCompositeNode")
      : nullptr;
  Dl_info owner = {};
  char expected[PATH_MAX] = {};
  char actual[PATH_MAX] = {};
  const auto providers = root_providers();
  if (!handle || !symbol || !dladdr(symbol, &owner) ||
      !realpath(library, expected) || !realpath(owner.dli_fname, actual) ||
      std::strcmp(expected, actual) != 0 || providers.size() != 1U ||
      providers.front() != expected)
  {
    gSystem->Exit(32);
  }
  std::cout << "THE134_ROOT_LOAD_PROVIDER_PASS provider=" << actual
            << " load_rc=" << rc
            << " preload_provider_count=0 provider_count=1" << std::endl;
  gSystem->Exit(0);
}
ROOT
root -l -b -q "${evidence_root}/root_load_smoke.C" \
  >"${evidence_root}/root_load_smoke.txt" 2>&1
grep -Fq "THE134_ROOT_LOAD_PROVIDER_PASS provider=${runtime_library} load_rc=0" \
  "${evidence_root}/root_load_smoke.txt"
BASH

runtime_library_sha256="$(sha256_file "$runtime_library")"
receipt="${output_root}/build_receipt.json"
python3 - "$receipt" "$RECEIPT_SCHEMA" "$source_manifest" \
  "$source_manifest_sha256" "$builder_sha256" "$setup_sha256" \
  "build" "${build_token#*:}" "$build_token" "$jobs" "$output_root" \
  "$coresoftware_commit" "$coresoftware_tree" "$commit_epoch" \
  "$RELEASE_NAME" "$expected_offline_main" "$MAPPING_PATCH_ID" \
  "$ROLE_SIZE_GUARD_ID" "$runtime_library" "$runtime_library_sha256" \
  "$install_root" "$baseline_calo_reco" "$baseline_calo_reco_sha256" \
  "$abi_additions_allowlist" "$abi_additions_allowlist_sha256" \
  "$ABI_ADDITIONS_PROVENANCE_SHA256" "$ABI_INTENTIONAL_EXCLUSION" \
  "$ABI_ROLE_GUARD_INLINE_ADDITION" \
  "$towerinfo_defs_header" "$towerinfo_defs_header_sha256" \
  "$towerinfo_container_header" "$towerinfo_container_header_sha256" \
  "$towerinfo_containerv1_header" "$towerinfo_containerv1_header_sha256" \
  "$calo_io_library" "$calo_io_library_sha256" \
  "${output_root}/evidence/baseline_dynamic_symbols.normalized.txt" \
  "${output_root}/evidence/candidate_dynamic_symbols.normalized.txt" \
  "${output_root}/evidence/baseline_dynsym_metadata.normalized.txt" \
  "${output_root}/evidence/candidate_dynsym_metadata.normalized.txt" \
  "${output_root}/evidence/abi_removed_dynsym_metadata.txt" \
  "${output_root}/evidence/baseline_version_inventory.txt" \
  "${output_root}/evidence/candidate_version_inventory.txt" \
  "${output_root}/evidence/abi_additions_allowlist.normalized.txt" \
  "${output_root}/evidence/abi_added_symbols.txt" \
  "${output_root}/evidence/abi_removed_symbols.txt" \
  "${output_root}/evidence/baseline_needed.txt" \
  "${output_root}/evidence/candidate_needed.txt" \
  "${output_root}/evidence/baseline_dictionary_inventory.txt" \
  "${output_root}/evidence/candidate_dictionary_inventory.txt" \
  "${output_root}/evidence/candidate_readelf_dynamic.txt" \
  "${output_root}/evidence/ldd.txt" \
  "${output_root}/evidence/resolved_dependency_paths.txt" \
  "${output_root}/evidence/defined_symbols.txt" \
  "${output_root}/evidence/compiler_version.txt" \
  "${output_root}/evidence/single_provider_probe.cc" \
  "${output_root}/evidence/single_provider_probe" \
  "${output_root}/evidence/single_provider_probe.txt" \
  "${output_root}/evidence/root_load_smoke.C" \
  "${output_root}/evidence/root_load_smoke.txt" \
  "${output_root}/evidence/towerinfo_mapping_equivalence.cc" \
  "${output_root}/evidence/towerinfo_mapping_equivalence" \
  "${output_root}/evidence/towerinfo_mapping_equivalence.txt" \
  "${output_root}/evidence/towerinfo_role_size_matrix.cc" \
  "${output_root}/evidence/towerinfo_role_size_matrix" \
  "${output_root}/evidence/towerinfo_role_size_matrix.txt" \
  "$EXPECTED_SONAME" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

(
    receipt_arg,
    schema,
    source_manifest_arg,
    source_manifest_sha,
    builder_sha,
    setup_sha,
    authorized_action,
    build_contract_sha,
    authorized_token,
    jobs,
    canonical_output_root,
    commit,
    tree,
    epoch,
    release_name,
    offline_main,
    mapping_id,
    role_size_guard_id,
    library_arg,
    library_sha,
    install_arg,
    baseline_library_arg,
    baseline_library_sha,
    abi_allowlist_arg,
    abi_allowlist_sha,
    abi_additions_provenance_sha,
    abi_intentional_exclusion,
    abi_role_guard_inline_addition,
    towerinfo_defs_arg,
    towerinfo_defs_sha,
    towerinfo_container_arg,
    towerinfo_container_sha,
    towerinfo_containerv1_arg,
    towerinfo_containerv1_sha,
    calo_io_arg,
    calo_io_sha,
    baseline_symbols_arg,
    candidate_symbols_arg,
    baseline_dynsym_metadata_arg,
    candidate_dynsym_metadata_arg,
    removed_dynsym_metadata_arg,
    baseline_version_inventory_arg,
    candidate_version_inventory_arg,
    normalized_allowlist_arg,
    added_symbols_arg,
    removed_symbols_arg,
    baseline_needed_arg,
    candidate_needed_arg,
    baseline_dictionary_arg,
    candidate_dictionary_arg,
    readelf_arg,
    ldd_arg,
    resolved_dependencies_arg,
    symbols_arg,
    compiler_arg,
    provider_probe_source_arg,
    provider_probe_binary_arg,
    provider_probe_arg,
    root_load_source_arg,
    root_load_arg,
    mapping_equivalence_source_arg,
    mapping_equivalence_binary_arg,
    mapping_equivalence_arg,
    role_size_matrix_source_arg,
    role_size_matrix_binary_arg,
    role_size_matrix_arg,
    expected_soname,
) = sys.argv[1:]
receipt = Path(receipt_arg)
source_manifest = json.loads(Path(source_manifest_arg).read_text(encoding="utf-8"))
library = Path(library_arg)
install = Path(install_arg)
readelf_text = Path(readelf_arg).read_text(encoding="utf-8")
symbols_text = Path(symbols_arg).read_text(encoding="utf-8")

def file_receipt(path_arg: str, *, display_path=None) -> dict:
    path = Path(path_arg)
    content = path.read_bytes()
    return {
        "path": display_path if display_path is not None else path.name,
        "sha256": hashlib.sha256(content).hexdigest(),
    }

def lines(path_arg: str) -> list[str]:
    return Path(path_arg).read_text(encoding="utf-8").splitlines()

soname_matches = re.findall(r"\(SONAME\).*?\[(.*?)\]", readelf_text)
photon_count = symbols_text.count(
    "PhotonClusterBuilder::process_event(PHCompositeNode*)"
)
topo_count = symbols_text.count(
    "RawClusterBuilderTopo::process_event(PHCompositeNode*)"
)

other_shared_objects = []
for path in sorted((install / "lib").glob("*.so*")):
    if path.is_symlink() or path.resolve() == library.resolve():
        continue
    other_shared_objects.append(path.relative_to(install).as_posix())

symlinks = {}
for name in ("libcalo_reco.so", expected_soname):
    path = install / "lib" / name
    symlinks[f"install/lib/{name}"] = path.readlink().as_posix() if path.is_symlink() else None

payload = {
    "abi": {
        "additions_allowlist": {
            **file_receipt(
                abi_allowlist_arg,
                display_path=Path(abi_allowlist_arg).name,
            ),
            "expected_sha256": abi_allowlist_sha,
            "normalized": {
                **file_receipt(
                    normalized_allowlist_arg,
                    display_path="evidence/abi_additions_allowlist.normalized.txt",
                ),
                "count": len(lines(normalized_allowlist_arg)),
            },
            "provenance": {
                "immutable_prior_override_sha256": (
                    abi_additions_provenance_sha
                ),
                "ana560_role_guard_inline_additions": [
                    {
                        "authority": {
                            "path": towerinfo_container_arg,
                            "sha256": towerinfo_container_sha,
                        },
                        "reason": (
                            "detector-role guard ODR-uses the pinned ana.560 "
                            "TowerInfoContainer::get_detectorid() inline method"
                        ),
                        "symbol": abi_role_guard_inline_addition,
                    }
                ],
                "intentional_exclusions": [
                    {
                        "reason": (
                            "detector-explicit RawClusterBuilderTopo mapping "
                            "removes generic container-selected encode_key"
                        ),
                        "symbol": abi_intentional_exclusion,
                    }
                ],
            },
        },
        "added_symbols": {
            **file_receipt(
                added_symbols_arg,
                display_path="evidence/abi_added_symbols.txt",
            ),
            "count": len(lines(added_symbols_arg)),
        },
        "baseline": {
            "dynsym_metadata": {
                **file_receipt(
                    baseline_dynsym_metadata_arg,
                    display_path=(
                        "evidence/baseline_dynsym_metadata.normalized.txt"
                    ),
                ),
                "count": len(lines(baseline_dynsym_metadata_arg)),
            },
            "dynamic_symbols": {
                **file_receipt(
                    baseline_symbols_arg,
                    display_path="evidence/baseline_dynamic_symbols.normalized.txt",
                ),
                "count": len(lines(baseline_symbols_arg)),
            },
            "library": {
                "path": baseline_library_arg,
                "sha256": baseline_library_sha,
            },
            "needed": lines(baseline_needed_arg),
            "version_inventory": {
                **file_receipt(
                    baseline_version_inventory_arg,
                    display_path="evidence/baseline_version_inventory.txt",
                ),
                "entries": lines(baseline_version_inventory_arg),
            },
        },
        "candidate": {
            "dynsym_metadata": {
                **file_receipt(
                    candidate_dynsym_metadata_arg,
                    display_path=(
                        "evidence/candidate_dynsym_metadata.normalized.txt"
                    ),
                ),
                "count": len(lines(candidate_dynsym_metadata_arg)),
            },
            "dynamic_symbols": {
                **file_receipt(
                    candidate_symbols_arg,
                    display_path="evidence/candidate_dynamic_symbols.normalized.txt",
                ),
                "count": len(lines(candidate_symbols_arg)),
            },
            "needed": lines(candidate_needed_arg),
            "version_inventory": {
                **file_receipt(
                    candidate_version_inventory_arg,
                    display_path="evidence/candidate_version_inventory.txt",
                ),
                "entries": lines(candidate_version_inventory_arg),
            },
        },
        "dictionary_inventory": {
            "baseline": {
                **file_receipt(
                    baseline_dictionary_arg,
                    display_path="evidence/baseline_dictionary_inventory.txt",
                ),
                "entries": lines(baseline_dictionary_arg),
            },
            "candidate": {
                **file_receipt(
                    candidate_dictionary_arg,
                    display_path="evidence/candidate_dictionary_inventory.txt",
                ),
                "entries": lines(candidate_dictionary_arg),
            },
        },
        "needed_exact_match": lines(baseline_needed_arg) == lines(candidate_needed_arg),
        "removed_symbols": {
            **file_receipt(
                removed_symbols_arg,
                display_path="evidence/abi_removed_symbols.txt",
            ),
            "count": len(lines(removed_symbols_arg)),
        },
        "removed_dynsym_metadata": {
            **file_receipt(
                removed_dynsym_metadata_arg,
                display_path="evidence/abi_removed_dynsym_metadata.txt",
            ),
            "count": len(lines(removed_dynsym_metadata_arg)),
        },
        "rpath_runpath_absent": not bool(
            re.search(r"\((?:RPATH|RUNPATH)\)", readelf_text)
        ),
        "soname": soname_matches[0] if len(soname_matches) == 1 else None,
        "soname_expected": expected_soname,
        "status": "PASS",
        "version_inventory_exact_match": (
            lines(baseline_version_inventory_arg)
            == lines(candidate_version_inventory_arg)
        ),
        "declared_implementation_replacements": [
            {
                "path": "PhotonClusterBuilder.cc",
                "sha256": source_manifest["overlay"][
                    "PhotonClusterBuilder.cc"
                ]["staged_sha256"],
            },
            {
                "mapping_patch_id": mapping_id,
                "path": "RawClusterBuilderTopo.cc",
                "role_size_guard_id": role_size_guard_id,
                "sha256": source_manifest["mapping_patch"]["sha256_after"],
            },
        ],
    },
    "artifact": {
        "library": f"install/lib/{library.name}",
        "sha256": library_sha,
        "symlinks": symlinks,
    },
    "build": {
        "authorization": {
            "action": authorized_action,
            "canonical_output_root": canonical_output_root,
            "contract_sha256": build_contract_sha,
            "jobs": int(jobs),
            "token": authorized_token,
        },
        "builder_sha256": builder_sha,
        "compiler": Path(compiler_arg).read_text(encoding="utf-8").strip(),
        "coresoftware_commit": commit,
        "coresoftware_commit_authority": commit,
        "coresoftware_commit_pinned": True,
        "coresoftware_tree": tree,
        "offline_main_authority": offline_main,
        "offline_main_pinned": True,
        "source_date_epoch": int(epoch),
        "setup_script_authority": str(Path("/opt/sphenix/core/bin/sphenix_setup.sh")),
        "setup_script_pinned": True,
    },
    "mapping_patch": {
        "authority": {
            "CaloIO": {
                "path": calo_io_arg,
                "sha256": calo_io_sha,
            },
            "TowerInfoContainer.h": {
                "path": towerinfo_container_arg,
                "sha256": towerinfo_container_sha,
            },
            "TowerInfoContainerv1.h": {
                "path": towerinfo_containerv1_arg,
                "sha256": towerinfo_containerv1_sha,
            },
            "TowerInfoDefs.h": {
                "path": towerinfo_defs_arg,
                "sha256": towerinfo_defs_sha,
            },
        },
        "exhaustive_equivalence": {
            **file_receipt(
                mapping_equivalence_arg,
                display_path="evidence/towerinfo_mapping_equivalence.txt",
            ),
            "binary": file_receipt(
                mapping_equivalence_binary_arg,
                display_path="evidence/towerinfo_mapping_equivalence",
            ),
            "CEMC_channels": 24576,
            "HCAL_channels": 1536,
            "round_trip": True,
            "source": file_receipt(
                mapping_equivalence_source_arg,
                display_path="evidence/towerinfo_mapping_equivalence.cc",
            ),
            "status": "PASS",
        },
        "id": mapping_id,
        "scientific_controls_changed": source_manifest["mapping_patch"][
            "changed_scientific_controls"
        ],
        "sha256_after": source_manifest["mapping_patch"]["sha256_after"],
        "validated_detector_maps": source_manifest["mapping_patch"]["detector_maps"],
    },
    "runtime": {
        "dependency_policy": {
            "allowed_roots": [
                "candidate install root",
                "exact ana.560 OFFLINE_MAIN",
                "versioned ana.560 platform root",
                "/lib",
                "/lib64",
                "/usr/lib",
                "/usr/lib64",
            ],
            "resolved_paths": {
                **file_receipt(
                    resolved_dependencies_arg,
                    display_path="evidence/resolved_dependency_paths.txt",
                ),
                "entries": lines(resolved_dependencies_arg),
            },
            "status": "PASS",
        },
        "ldd_not_found": "not found" in Path(ldd_arg).read_text(encoding="utf-8"),
        "offline_main": offline_main,
        "mutable_user_dependency": bool(
            re.search(
                r"/sphenix/(?:u|user)/",
                Path(ldd_arg).read_text(encoding="utf-8"),
            )
        ),
        "release": release_name,
        "root_load": {
            **file_receipt(
                root_load_arg,
                display_path="evidence/root_load_smoke.txt",
            ),
            "load_return_code": 0,
            "preload_provider_count": 0,
            "provider_count": 1,
            "provider_realpath": str(library.resolve()),
            "source": file_receipt(
                root_load_source_arg,
                display_path="evidence/root_load_smoke.C",
            ),
            "status": "PASS",
        },
        "setup_script_sha256": setup_sha,
    },
    "role_size_guard": {
        **source_manifest["role_size_guard"],
        "executable_matrix": {
            **file_receipt(
                role_size_matrix_arg,
                display_path="evidence/towerinfo_role_size_matrix.txt",
            ),
            "binary": file_receipt(
                role_size_matrix_binary_arg,
                display_path="evidence/towerinfo_role_size_matrix",
            ),
            "cases": {
                "legacy_exact": 2,
                "typed": 2,
                "wrong_size": 2,
                "wrong_valid": 2,
            },
            "source": file_receipt(
                role_size_matrix_source_arg,
                display_path="evidence/towerinfo_role_size_matrix.cc",
            ),
            "status": "PASS",
        },
        "id": role_size_guard_id,
    },
    "schema": schema,
    "single_provider": {
        "forbidden_standalone_provider_count": len(other_shared_objects),
        "other_installed_shared_objects": other_shared_objects,
        "photon_cluster_builder_process_event_definitions": photon_count,
        "provider": f"install/lib/{library.name}",
        "provider_probe": {
            **file_receipt(
                provider_probe_arg,
                display_path="evidence/single_provider_probe.txt",
            ),
            "preload_provider_count": 0,
            "provider_count": 1,
            "provider_realpath": str(library.resolve()),
            "binary": file_receipt(
                provider_probe_binary_arg,
                display_path="evidence/single_provider_probe",
            ),
            "source": file_receipt(
                provider_probe_source_arg,
                display_path="evidence/single_provider_probe.cc",
            ),
            "status": "PASS",
        },
        "raw_cluster_builder_topo_process_event_definitions": topo_count,
    },
    "source_manifest": {
        "path": "source_manifest.json",
        "sha256": source_manifest_sha,
    },
    "status": "PASS",
}

if payload["abi"]["soname"] != expected_soname:
    raise SystemExit("receipt SONAME closure failed")
if not payload["abi"]["needed_exact_match"]:
    raise SystemExit("receipt NEEDED closure failed")
if not payload["abi"]["rpath_runpath_absent"]:
    raise SystemExit("receipt RPATH/RUNPATH closure failed")
if payload["abi"]["removed_symbols"]["count"] != 0:
    raise SystemExit("receipt ABI removal closure failed")
if payload["abi"]["removed_dynsym_metadata"]["count"] != 0:
    raise SystemExit("receipt dynamic-symbol metadata closure failed")
if payload["abi"]["added_symbols"]["sha256"] != payload["abi"]["additions_allowlist"]["normalized"]["sha256"]:
    raise SystemExit("receipt ABI additions allowlist closure failed")
if not payload["abi"]["version_inventory_exact_match"]:
    raise SystemExit("receipt symbol-version inventory closure failed")
if payload["abi"]["dictionary_inventory"]["baseline"]["entries"] != payload["abi"]["dictionary_inventory"]["candidate"]["entries"]:
    raise SystemExit("receipt dictionary inventory closure failed")
if photon_count != 1 or topo_count != 1 or other_shared_objects:
    raise SystemExit("receipt single-provider closure failed")
if payload["runtime"]["ldd_not_found"]:
    raise SystemExit("receipt dependency closure failed")
if payload["runtime"]["mutable_user_dependency"]:
    raise SystemExit("receipt contains a mutable private runtime dependency")
if payload["runtime"]["dependency_policy"]["status"] != "PASS":
    raise SystemExit("receipt dependency allowlist closure failed")
if payload["mapping_patch"]["scientific_controls_changed"]:
    raise SystemExit("mapping patch changed a scientific control")
if source_manifest["role_size_guard"]["id"] != role_size_guard_id:
    raise SystemExit("role/size guard identity differs")
if "THE134_SINGLE_PROVIDER_PROBE_PASS" not in Path(provider_probe_arg).read_text(encoding="utf-8"):
    raise SystemExit("provider probe evidence marker is missing")
if "THE134_ROOT_LOAD_PROVIDER_PASS" not in Path(root_load_arg).read_text(encoding="utf-8"):
    raise SystemExit("ROOT provider evidence marker is missing")
if "THE134_TOWERINFO_MAPPING_EQUIVALENCE_PASS" not in Path(mapping_equivalence_arg).read_text(encoding="utf-8"):
    raise SystemExit("mapping equivalence evidence marker is missing")
if "THE134_TOWERINFO_ROLE_SIZE_MATRIX_PASS" not in Path(role_size_matrix_arg).read_text(encoding="utf-8"):
    raise SystemExit("role/size matrix evidence marker is missing")
if authorized_action != "build" or authorized_token != f"the134-ana560-caloreco-build:{build_contract_sha}":
    raise SystemExit("receipt build authorization binding differs")
if Path(canonical_output_root).resolve(strict=True) != receipt.parent.resolve(strict=True):
    raise SystemExit("receipt output-root binding differs")

receipt.write_text(
    json.dumps(payload, indent=2, sort_keys=True, separators=(",", ": ")) + "\n",
    encoding="utf-8",
)
PY

receipt_sha256="$(sha256_file "$receipt")"
chmod -R a-w "$output_root"
cat <<EOF
THE134_ANA560_CALORECO_BUILD_PASS
  output_root: ${output_root}
  source_manifest: ${source_manifest}
  source_manifest_sha256: ${source_manifest_sha256}
  build_receipt: ${receipt}
  build_receipt_sha256: ${receipt_sha256}
  library: ${runtime_library}
  library_sha256: ${runtime_library_sha256}
  soname: ${EXPECTED_SONAME}
  provider_count: 1
EOF
