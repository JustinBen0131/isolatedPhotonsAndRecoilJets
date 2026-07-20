#!/usr/bin/env bash
# Build a sealed new.17 RecoilJets runtime for the paired PPG12 oracle.
#
# Default mode is a read-only plan.  The build is foreground-only, writes to a
# new contained output root, and requires the exact token printed by plan mode.
# It never submits Condor, installs into a user prefix, or mutates the checkout.

set -euo pipefail

die() {
  printf 'PPG12_ORACLE_NEW17_BUILD_FAIL: %s\n' "$*" >&2
  exit 2
}

usage() {
  cat <<'EOF'
Usage:
  build_ppg12_oracle_new17_runtime.sh --output-dir ABS [--jobs N] \
    [--photon-source-dir ABS] [--ppg-repo ABS] [--ppg-revision SHA]
  build_ppg12_oracle_new17_runtime.sh --build --token TOKEN \
    --output-dir ABS [--jobs N] [--photon-source-dir ABS] \
    [--ppg-repo ABS] [--ppg-revision SHA]

Default mode prints the immutable build contract and authorization token.
--build performs the foreground build only when TOKEN exactly matches that
contract.  ABS must be a new path below REPO/.recoiljets_tmp or /tmp.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../.." && pwd -P)"
expected_offline="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.17"

mode=plan
provided_token=""
output_dir=""
setup_script="/opt/sphenix/core/bin/sphenix_setup.sh"
jobs=4
photon_source_dir="${repo_root}/src"
ppg_repo="${repo_root}/ppg12codeGit"
# This is the last CaloAna24 source revision before the archived June PPG12
# production.  The working tree is deliberately ignored: git-show exports the
# committed source into the sealed build root.
ppg_revision="1c0ff86bf0ebabfba63a1abc4512cbe59fe48e31"

while (($#)); do
  case "$1" in
    --build) mode=build; shift ;;
    --token) [[ $# -ge 2 ]] || die "--token requires a value"; provided_token="$2"; shift 2 ;;
    --output-dir) [[ $# -ge 2 ]] || die "--output-dir requires a value"; output_dir="$2"; shift 2 ;;
    --setup-script) [[ $# -ge 2 ]] || die "--setup-script requires a value"; setup_script="$2"; shift 2 ;;
    --photon-source-dir) [[ $# -ge 2 ]] || die "--photon-source-dir requires a value"; photon_source_dir="$2"; shift 2 ;;
    --ppg-repo) [[ $# -ge 2 ]] || die "--ppg-repo requires a value"; ppg_repo="$2"; shift 2 ;;
    --ppg-revision) [[ $# -ge 2 ]] || die "--ppg-revision requires a value"; ppg_revision="$2"; shift 2 ;;
    --jobs) [[ $# -ge 2 ]] || die "--jobs requires a value"; jobs="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ "$jobs" =~ ^[1-8]$ ]] || die "--jobs must be an integer from 1 through 8"
[[ -n "$output_dir" && "$output_dir" == /* ]] || die "--output-dir must be absolute"
[[ "$output_dir" != *$'\n'* && "$output_dir" != *$'\r'* && "$output_dir" != *' '* ]] || \
  die "--output-dir contains whitespace or a newline"
case "$output_dir" in
  "${repo_root}/.recoiljets_tmp/"*|/tmp/*) ;;
  *) die "--output-dir must be below ${repo_root}/.recoiljets_tmp or /tmp" ;;
esac
[[ ! -e "$output_dir" ]] || die "output path already exists: $output_dir"
[[ "$setup_script" == /* && -f "$setup_script" && -s "$setup_script" ]] || \
  die "setup script is missing or empty: $setup_script"
[[ "$photon_source_dir" == /* && -d "$photon_source_dir" ]] || \
  die "--photon-source-dir must be an existing absolute directory: $photon_source_dir"
[[ "$ppg_repo" == /* && -d "$ppg_repo/.git" ]] || \
  die "--ppg-repo must be an existing absolute Git checkout: $ppg_repo"
[[ "$ppg_revision" =~ ^[0-9a-f]{40}$ ]] || die "--ppg-revision must be a full commit SHA"
ppg_repo_real="$(cd "$ppg_repo" && pwd -P)"
git_safe=(-c "safe.directory=${ppg_repo_real}")
git "${git_safe[@]}" -C "$ppg_repo_real" cat-file -e "${ppg_revision}^{commit}" 2>/dev/null || \
  die "--ppg-revision is not present in --ppg-repo: $ppg_revision"
ppg_source_names=(configure.ac Makefile.am autogen.sh CaloAna24.cc CaloAna24.h)
for source_name in "${ppg_source_names[@]}"; do
  git "${git_safe[@]}" -C "$ppg_repo_real" cat-file -e \
    "${ppg_revision}:anatreemaker/source/${source_name}" 2>/dev/null || \
    die "PPG12 source revision lacks anatreemaker/source/${source_name}"
done

canonical_photon_cc="${photon_source_dir}/PhotonClusterBuilder.cc"
canonical_photon_h="${photon_source_dir}/PhotonClusterBuilder.h"
recoil_source="${repo_root}/src"
macro_wrapper_source="${repo_root}/macros/Fun4All_recoilJets.C"
macro_impl_source="${repo_root}/macros/Fun4All_recoilJets_unified_impl.C"

contract_inputs=(
  "$canonical_photon_cc"
  "$canonical_photon_h"
  "${recoil_source}/configure.ac"
  "${recoil_source}/Makefile.am"
  "${recoil_source}/autogen.sh"
  "${recoil_source}/RecoilJets.cc"
  "${recoil_source}/RecoilJets.h"
  "${recoil_source}/PPG12SimWeight.h"
  "$macro_wrapper_source"
  "$macro_impl_source"
  "$setup_script"
)
for input in "${contract_inputs[@]}"; do
  [[ "$input" == /* && -f "$input" && -s "$input" ]] || die "missing contract input: $input"
done

sha256_file() {
  python3 - "$1" <<'PY'
from pathlib import Path
import hashlib
import sys

h = hashlib.sha256()
with Path(sys.argv[1]).open("rb") as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b""):
        h.update(block)
print(h.hexdigest())
PY
}

contract_token="$({
  printf '%s\n' \
    'schema_version=1' \
    'purpose=ppg12_paired_oracle_recoil_runtime' \
    'runtime_profile=new.17' \
    "offline_main=${expected_offline}" \
    "output_dir=${output_dir}" \
    "jobs=${jobs}"
  printf 'ppg_repo=%s\nppg_revision=%s\n' "$ppg_repo_real" "$ppg_revision"
  for source_name in "${ppg_source_names[@]}"; do
    printf 'ppg_source=%s sha256=%s\n' "$source_name" "$(
      git "${git_safe[@]}" -C "$ppg_repo_real" show \
        "${ppg_revision}:anatreemaker/source/${source_name}" \
        | python3 -c 'import hashlib,sys; print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())'
    )"
  done
  for input in "${contract_inputs[@]}"; do
    printf 'input=%s sha256=%s\n' "$input" "$(sha256_file "$input")"
  done
} | python3 -c 'import hashlib,sys; print("ppg12-new17:" + hashlib.sha256(sys.stdin.buffer.read()).hexdigest())')"

if [[ "$mode" == plan ]]; then
  cat <<EOF
PPG12_ORACLE_NEW17_BUILD_PLAN
  output_root: ${output_dir}
  runtime: new.17
  OFFLINE_MAIN: ${expected_offline}
  build: source-locked libCaloAna24.so plus libRecoilJets.so with renamed PPG12-oracle photon builder
  ppg_source_revision: ${ppg_revision}
  release copies: exact new.17 libcalo_reco.so, libclusteriso.so, libjetbase.so (cp -L)
  macro staging: PP wrapper/implementation with exact-count path rewrites
  validation: forbidden-route scan, readelf, ldd, and ROOT load/header smoke
  jobs: ${jobs}
  token: ${contract_token}
No files were created.  Re-run with --build --token '${contract_token}'.
EOF
  exit 0
fi

[[ "$provided_token" == "$contract_token" ]] || die "authorization token does not match the current build contract"
[[ "${RJ_PPG12_ORACLE_NEW17_CLEAN_ENV:-0}" != 1 ]] || clean_env=1

# The shell executing this file may carry an analysis release or private
# installation in multiple routing variables.  Re-exec with a blank
# environment before sourcing the one supported release.
if [[ "${clean_env:-0}" != 1 ]]; then
  self="${script_dir}/$(basename "${BASH_SOURCE[0]}")"
  exec /usr/bin/env -i \
    HOME="${HOME:-/tmp}" \
    USER="${USER:-unknown}" \
    LOGNAME="${LOGNAME:-${USER:-unknown}}" \
    SHELL=/bin/bash \
    PATH=/usr/bin:/bin:/usr/sbin:/sbin \
    RJ_PPG12_ORACLE_NEW17_CLEAN_ENV=1 \
    /bin/bash --noprofile --norc "$self" \
      --build --token "$provided_token" --output-dir "$output_dir" \
      --setup-script "$setup_script" --jobs "$jobs" \
      --photon-source-dir "$photon_source_dir" \
      --ppg-repo "$ppg_repo_real" --ppg-revision "$ppg_revision"
fi

umask 077

unset OFFLINE_MAIN MYINSTALL ROOT_INCLUDE_PATH LD_LIBRARY_PATH PYTHONPATH \
  CMAKE_PREFIX_PATH CPATH CPLUS_INCLUDE_PATH LIBRARY_PATH PKG_CONFIG_PATH
export PATH=/usr/bin:/bin:/usr/sbin:/sbin
set +e
set +u
# shellcheck disable=SC1090
source "$setup_script" -n new.17
setup_status=$?
set -u
set -e
[[ $setup_status -eq 0 ]] || die "sPHENIX new.17 setup failed"
[[ "${OFFLINE_MAIN:-}" == "$expected_offline" ]] || \
  die "setup resolved unexpected OFFLINE_MAIN: ${OFFLINE_MAIN:-<unset>}"

for command_name in python3 make aclocal automake autoconf libtoolize root root-config readelf ldd cmp git awk; do
  command -v "$command_name" >/dev/null 2>&1 || die "required build command is unavailable: $command_name"
done

forbidden_routes() {
  local value="${PATH:-}:${LD_LIBRARY_PATH:-}:${ROOT_INCLUDE_PATH:-}:${PYTHONPATH:-}:${CMAKE_PREFIX_PATH:-}:${CPATH:-}:${CPLUS_INCLUDE_PATH:-}:${LIBRARY_PATH:-}:${PKG_CONFIG_PATH:-}"
  [[ "$value" != *'/release/release_ana/'* ]] || return 1
  [[ "$value" != *'/sphenix/u/patsfan753/thesisAnalysis/install'* ]] || return 1
  [[ "$value" != *'/sphenix/user/patsfan753/install'* ]] || return 1
  [[ "$value" != *'/sphenix/u/patsfan753/thesisAnalysis_auau/install'* ]] || return 1
  return 0
}
forbidden_routes || die "new.17 setup retained an analysis-release or private-install route"
base_ld_library_path="${LD_LIBRARY_PATH:-}"
base_root_include_path="${ROOT_INCLUDE_PATH:-}"
root_libdir="$(root-config --libdir)"
[[ "$root_libdir" == /* && -d "$root_libdir" ]] || \
  die "ROOT library directory is missing or non-absolute: $root_libdir"
release_link_dirs=()
for link_dir in "${OFFLINE_MAIN}/lib" "${OFFLINE_MAIN}/lib64" "$root_libdir"; do
  [[ -d "$link_dir" ]] && release_link_dirs+=("$link_dir")
done
[[ ${#release_link_dirs[@]} -ge 2 ]] || \
  die "could not establish sealed new.17 and ROOT linker directories"
release_ldflags=""
for link_dir in "${release_link_dirs[@]}"; do
  release_ldflags+=" -L${link_dir}"
done
sealed_link_dirs="$(IFS=:; printf '%s' "${release_link_dirs[*]}")"

mkdir -p "$output_dir"
state_file="${output_dir}/BUILD_STATE"
printf 'RUNNING\n' > "$state_file"
completed=0
finish_state() {
  if [[ $completed -eq 0 ]]; then
    printf 'FAILED\n' > "$state_file"
  fi
}
trap finish_state EXIT

stage_root="${output_dir}/stage"
build_root="${output_dir}/build"
install_root="${output_dir}/install"
runtime_root="${output_dir}/runtime"
log_root="${output_dir}/logs"
recoil_stage="${stage_root}/RecoilJets"
ppg_stage="${stage_root}/PPG12CaloAna24"
mkdir -p "$stage_root" "$build_root/recoiljets" "$build_root/ppg12" \
  "$recoil_stage" "$install_root" "$runtime_root/lib" "$runtime_root/include/caloreco" \
  "$runtime_root/include/caloana" "$runtime_root/macros" "$log_root" "$ppg_stage"

for source_name in "${ppg_source_names[@]}"; do
  git "${git_safe[@]}" -C "$ppg_repo_real" show \
    "${ppg_revision}:anatreemaker/source/${source_name}" \
    > "${ppg_stage}/${source_name}"
done
chmod 700 "${ppg_stage}/autogen.sh"

for source_name in configure.ac Makefile.am autogen.sh RecoilJets.cc RecoilJets.h PPG12SimWeight.h; do
  cp -f "${recoil_source}/${source_name}" "${recoil_stage}/${source_name}"
done
cp -f "$canonical_photon_cc" "${recoil_stage}/PPG12OraclePhotonClusterBuilder.cc"
cp -f "$canonical_photon_h" "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h"

replace_exact() {
  local path="$1"
  local expected_count="$2"
  local old="$3"
  local new="$4"
  python3 - "$path" "$expected_count" "$old" "$new" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
expected = int(sys.argv[2])
old = sys.argv[3]
new = sys.argv[4]
text = path.read_text()
observed = text.count(old)
if observed != expected:
    raise SystemExit(
        f"exact rewrite count mismatch for {path}: expected {expected}, "
        f"observed {observed}, needle={old!r}"
    )
path.write_text(text.replace(old, new))
PY
}

replace_word_exact() {
  local path="$1"
  local expected_count="$2"
  local old="$3"
  local new="$4"
  python3 - "$path" "$expected_count" "$old" "$new" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
expected = int(sys.argv[2])
old = sys.argv[3]
new = sys.argv[4]
text = path.read_text()
pattern = re.compile(rf"\b{re.escape(old)}\b")
observed = len(pattern.findall(text))
if observed != expected:
    raise SystemExit(
        f"exact word-rewrite count mismatch for {path}: expected {expected}, "
        f"observed {observed}, token={old!r}"
    )
path.write_text(pattern.sub(new, text))
PY
}

# Rebuild the preserved PPG12 source against the same current new.17 headers
# and libraries used by the candidate.  The archived June binary is retained
# only as historical evidence: loading it against today's mutable new.17
# aborts on valid event-2 tower keys before the scientific comparison begins.
replace_exact "${ppg_stage}/configure.ac" 1 \
  'CXXFLAGS="$CXXFLAGS -Wall -Werror"' \
  'CXXFLAGS="$CXXFLAGS -Wall -Wno-error"'
replace_exact "${ppg_stage}/Makefile.am" 1 \
  $'-lcalotrigger_io \\ ' \
  $'-lcalotrigger_io \\'

(
  cd "${build_root}/ppg12"
  export LD_LIBRARY_PATH="${install_root}/lib:${base_ld_library_path}"
  export ROOT_INCLUDE_PATH="${install_root}/include:${base_root_include_path}"
  export CPPFLAGS="-I${install_root}/include"
  export LDFLAGS="-L${install_root}/lib${release_ldflags}"
  /bin/bash "${ppg_stage}/autogen.sh" --prefix="$install_root"
  make -j "$jobs"
  make install
) >"${log_root}/ppg12_build.log" 2>&1 || {
  tail -n 80 "${log_root}/ppg12_build.log" >&2 || true
  die "isolated source-locked libCaloAna24 build failed"
}

# These rewrites touch staged copies only.  The oracle builder is renamed so
# it can coexist with exact release libcalo_reco without an ODR/symbol
# collision.  It is compiled into libRecoilJets; release calibration and
# RawClusterBuilderTemplate remain untouched.
replace_exact "${recoil_stage}/configure.ac" 1 \
  'CXXFLAGS="$CXXFLAGS -Wall -Wextra -Werror -Wshadow"' \
  'CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Wno-error"'
replace_exact "${recoil_stage}/RecoilJets.cc" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloBase/PhotonClusterv1.h"' \
  '#include <calobase/PhotonClusterv1.h>'
replace_exact "${recoil_stage}/RecoilJets.cc" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"' \
  '// Oracle photon builder is registered by the staged Fun4All macro.'
replace_exact "${recoil_stage}/RecoilJets.h" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"' \
  '// Oracle photon builder is registered by the staged Fun4All macro.'
replace_word_exact "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" 4 \
  'PhotonClusterBuilder' 'PPG12OraclePhotonClusterBuilder'
replace_word_exact "${recoil_stage}/PPG12OraclePhotonClusterBuilder.cc" 31 \
  'PhotonClusterBuilder' 'PPG12OraclePhotonClusterBuilder'
replace_exact "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" 3 \
  'CALORECO_PHOTONCLUSTERBUILDER_H' \
  'CALOANA_PPG12ORACLEPHOTONCLUSTERBUILDER_H'
replace_exact "${recoil_stage}/Makefile.am" 1 \
  $'libRecoilJets_la_SOURCES = \\\n  RecoilJets.cc' \
  $'libRecoilJets_la_SOURCES = \\\n  RecoilJets.cc \\\n  PPG12OraclePhotonClusterBuilder.cc'
replace_exact "${recoil_stage}/Makefile.am" 1 \
  $'  -lg4detectors_io \\\n  -lphg4hit' \
  $'  -lg4detectors_io \\\n  -lphg4hit \\\n  -lg4eval \\\n  -lcdbobjects \\\n  -lffamodules \\\n  -lglobalvertex_io \\\n  -lTMVA \\\n  -lTMVAUtils'
replace_exact "${recoil_stage}/Makefile.am" 1 \
  'libRecoilJets_la_LDFLAGS = -no-undefined -version-info 0:0:0' \
  'libRecoilJets_la_LDFLAGS = $(AM_LDFLAGS) -L$(ROOTSYS)/lib -no-undefined -version-info 0:0:0'
replace_exact "${recoil_stage}/Makefile.am" 1 \
  $'pkginclude_HEADERS = \\\n  RecoilJets.h' \
  $'pkginclude_HEADERS = \\\n  RecoilJets.h \\\n  PPG12OraclePhotonClusterBuilder.h'

(
  cd "${build_root}/recoiljets"
  export LD_LIBRARY_PATH="${install_root}/lib:${base_ld_library_path}"
  export ROOT_INCLUDE_PATH="${install_root}/include:${base_root_include_path}"
  export CPPFLAGS="-I${install_root}/include"
  export LDFLAGS="-L${install_root}/lib${release_ldflags}"
  /bin/bash "${recoil_stage}/autogen.sh" --prefix="$install_root"
  make -j "$jobs"
  make install
) >"${log_root}/recoiljets_build.log" 2>&1 || {
  tail -n 80 "${log_root}/recoiljets_build.log" >&2 || true
  die "isolated libRecoilJets build failed"
}

resolve_installed_lib() {
  local name="$1"
  local candidate
  for candidate in "${install_root}/lib/${name}" "${install_root}/lib64/${name}"; do
    if [[ -r "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  die "isolated build did not install ${name}"
}

resolve_release_lib() {
  local name="$1"
  local candidate
  for candidate in "${OFFLINE_MAIN}/lib/${name}" "${OFFLINE_MAIN}/lib64/${name}"; do
    if [[ -r "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  die "new.17 does not provide ${name}"
}

built_recoil="$(resolve_installed_lib libRecoilJets.so)"
built_ppg="$(resolve_installed_lib libCaloAna24.so)"
release_calo="$(resolve_release_lib libcalo_reco.so)"
release_clusteriso="$(resolve_release_lib libclusteriso.so)"
release_jetbase="$(resolve_release_lib libjetbase.so)"

cp -L "$release_calo" "${runtime_root}/lib/libcalo_reco.so"
cp -L "$built_ppg" "${runtime_root}/lib/libCaloAna24.so"
cp -L "$built_recoil" "${runtime_root}/lib/libRecoilJets.so"
cp -L "$release_clusteriso" "${runtime_root}/lib/libclusteriso.so"
cp -L "$release_jetbase" "${runtime_root}/lib/libjetbase.so"
cmp -s "$release_calo" "${runtime_root}/lib/libcalo_reco.so" || \
  die "copied libcalo_reco differs from exact new.17 source"
cmp -s "$release_clusteriso" "${runtime_root}/lib/libclusteriso.so" || \
  die "copied libclusteriso differs from exact new.17 source"
cmp -s "$release_jetbase" "${runtime_root}/lib/libjetbase.so" || \
  die "copied libjetbase differs from exact new.17 source"
cp -f "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" \
  "${runtime_root}/include/caloana/PPG12OraclePhotonClusterBuilder.h"
# The paired-oracle manifest retains its established role/path.  This file's
# contents declare only the renamed oracle class, never the release class.
cp -f "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" \
  "${runtime_root}/include/caloreco/PhotonClusterBuilder.h"
cp -f "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" \
  "${runtime_root}/include/caloreco/PPG12OraclePhotonClusterBuilder.h"
cp -f "${recoil_stage}/RecoilJets.h" "${runtime_root}/include/caloana/RecoilJets.h"

# Preserve each copied ELF's SONAME inside the sealed library directory.
for runtime_lib in "${runtime_root}"/lib/lib*.so; do
  soname="$(readelf -d "$runtime_lib" 2>/dev/null | awk -F'[][]' '/SONAME/ {print $2; exit}' || true)"
  if [[ -n "$soname" && "$soname" != "$(basename "$runtime_lib")" ]]; then
    ln -sfn "$(basename "$runtime_lib")" "${runtime_root}/lib/${soname}"
  fi
done

calo_calib="${OFFLINE_MAIN}/rootmacros/Calo_Calib.C"
[[ -f "$calo_calib" && -s "$calo_calib" ]] || \
  die "stock new.17 Calo_Calib.C is missing: $calo_calib"

runtime_wrapper="${runtime_root}/macros/Fun4All_recoilJets.C"
runtime_impl="${runtime_root}/macros/Fun4All_recoilJets_unified_impl.C"
cp -f "$macro_wrapper_source" "$runtime_wrapper"
cp -f "$macro_impl_source" "$runtime_impl"

replace_word_exact "$runtime_impl" 16 \
  'PhotonClusterBuilder' 'PPG12OraclePhotonClusterBuilder'
replace_exact "$runtime_wrapper" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Fun4All_recoilJets_unified_impl.C"' \
  "#include \"${runtime_impl}\""
replace_exact "$runtime_impl" 3 \
  '/sphenix/u/patsfan753/thesisAnalysis/install/include' \
  "${runtime_root}/include"
replace_exact "$runtime_impl" 2 \
  '/sphenix/u/patsfan753/thesisAnalysis_auau/install/include' \
  "${runtime_root}/include"
replace_exact "$runtime_impl" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src_AuAu/RecoilJets_AuAu.h"' \
  "#include \"${runtime_root}/include/caloana/RecoilJets.h\""
replace_exact "$runtime_impl" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src/RecoilJets.h"' \
  "#include \"${runtime_root}/include/caloana/RecoilJets.h\""
replace_exact "$runtime_impl" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Calo_Calib.C"' \
  "#include \"${calo_calib}\""
replace_exact "$runtime_impl" 1 \
  'R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis_auau/install/lib/libRecoilJetsAuAu.so)' \
  'R__LOAD_LIBRARY(libRecoilJets.so)'
replace_exact "$runtime_impl" 1 \
  'R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libRecoilJets.so)' \
  'R__LOAD_LIBRARY(libRecoilJets.so)'
replace_exact "$runtime_impl" 1 \
  'R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)' \
  'R__LOAD_LIBRARY(libclusteriso.so)'
replace_exact "$runtime_impl" 1 \
  'R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)' \
  'R__LOAD_LIBRARY(libjetbase.so)'

for staged_runtime_file in "$runtime_wrapper" "$runtime_impl" \
  "${runtime_root}/include/caloana/RecoilJets.h" \
  "${runtime_root}/include/caloana/PPG12OraclePhotonClusterBuilder.h"; do
  if grep -E '/release/release_ana/|/sphenix/u/patsfan753/thesisAnalysis(_auau)?/install|/sphenix/user/patsfan753/install' \
      "$staged_runtime_file" >/dev/null; then
    die "staged runtime retains a forbidden route: $staged_runtime_file"
  fi
done

runtime_ld="${runtime_root}/lib:${base_ld_library_path}"
for runtime_lib in \
  "${runtime_root}/lib/libCaloAna24.so" \
  "${runtime_root}/lib/libcalo_reco.so" \
  "${runtime_root}/lib/libRecoilJets.so" \
  "${runtime_root}/lib/libclusteriso.so" \
  "${runtime_root}/lib/libjetbase.so"; do
  basename_noext="$(basename "$runtime_lib" .so)"
  readelf -d "$runtime_lib" >"${log_root}/readelf_${basename_noext}.txt" 2>&1 || \
    die "readelf failed for $runtime_lib"
  LD_LIBRARY_PATH="$runtime_ld" ldd "$runtime_lib" \
    >"${log_root}/ldd_${basename_noext}.txt" 2>&1 || die "ldd failed for $runtime_lib"
  if grep -E 'not found|/release/release_ana/|/sphenix/u/patsfan753/thesisAnalysis(_auau)?/install|/sphenix/user/patsfan753/install' \
      "${log_root}/ldd_${basename_noext}.txt" >/dev/null; then
    die "ldd exposed an unresolved or forbidden dependency for $runtime_lib"
  fi
  if grep -E '/release/release_new/new\.[0-9]+' "${log_root}/ldd_${basename_noext}.txt" \
      | grep -vF "$expected_offline" >/dev/null; then
    die "ldd exposed a mixed new-release dependency for $runtime_lib"
  fi
done

smoke_macro="${build_root}/smoke_new17_runtime.C"
cat > "$smoke_macro" <<EOF
#include <caloana/PPG12OraclePhotonClusterBuilder.h>
#include <TSystem.h>
#include <iostream>
void smoke_new17_runtime()
{
  const char *libraries[] = {
    "${runtime_root}/lib/libCaloAna24.so",
    "${runtime_root}/lib/libcalo_reco.so",
    "${runtime_root}/lib/libclusteriso.so",
    "${runtime_root}/lib/libjetbase.so",
    "${runtime_root}/lib/libRecoilJets.so"
  };
  for (const char *library : libraries)
  {
    const int status = gSystem->Load(library);
    std::cout << "PPG12_ORACLE_ROOT_LOAD path=" << library
              << " status=" << status << std::endl;
    if (status < 0) gSystem->Exit(91);
  }
  PPG12OraclePhotonClusterBuilder *builder = nullptr;
  if (builder != nullptr) gSystem->Exit(92);
}
EOF
(
  export LD_LIBRARY_PATH="$runtime_ld"
  export ROOT_INCLUDE_PATH="${runtime_root}/include:${base_root_include_path}"
  root -l -b -q "${smoke_macro}"
) >"${log_root}/root_smoke.log" 2>&1 || {
  tail -n 80 "${log_root}/root_smoke.log" >&2 || true
  die "ROOT load/header smoke failed"
}
[[ "$(grep -c '^PPG12_ORACLE_ROOT_LOAD ' "${log_root}/root_smoke.log")" -eq 5 ]] || \
  die "ROOT smoke did not load all five runtime libraries"

forbidden_routes || die "build contaminated the active shell with a forbidden route"

build_receipt="${output_dir}/build_receipt.json"
runtime_manifest="${output_dir}/runtime_manifest.json"
python3 - \
  "$build_receipt" "$output_dir" "$install_root" "$runtime_root" \
  "$expected_offline" "$calo_calib" "$jobs" "$sealed_link_dirs" \
  "$canonical_photon_cc" "$canonical_photon_h" \
  "$ppg_repo_real" "$ppg_revision" \
  "${ppg_stage}/CaloAna24.cc" "${ppg_stage}/CaloAna24.h" \
  "${ppg_stage}/configure.ac" "${ppg_stage}/Makefile.am" \
  "${recoil_source}/RecoilJets.cc" "${recoil_source}/RecoilJets.h" \
  "$macro_wrapper_source" "$macro_impl_source" \
  "$release_calo" "$release_clusteriso" "$release_jetbase" "$log_root" <<'PY'
from pathlib import Path
import hashlib
import json
import os
import platform
import sys

(
    receipt, output_root, install_root, runtime_root, offline_main, calo_calib,
    jobs, sealed_link_dirs, photon_cc, photon_h, ppg_repo, ppg_revision,
    ppg_cc, ppg_h, ppg_configure, ppg_makefile, recoil_cc, recoil_h, macro, impl,
    calo_source, clusteriso_source, jetbase_source, log_root,
) = sys.argv[1:]

def digest(path: str | Path) -> str:
    value = Path(path)
    h = hashlib.sha256()
    with value.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

source_paths = [
    photon_cc, photon_h, ppg_cc, ppg_h, ppg_configure, ppg_makefile,
    recoil_cc, recoil_h, macro, impl,
    calo_source, clusteriso_source, jetbase_source, calo_calib,
]
log_paths = sorted(str(path) for path in Path(log_root).iterdir() if path.is_file())
data = {
    "schema_version": 1,
    "purpose": "ppg12_paired_oracle_recoil_runtime",
    "runtime_profile": "new.17",
    "offline_main": offline_main,
    "isolated_build": True,
    "output_root": output_root,
    "install_root": install_root,
    "runtime_root": runtime_root,
    "jobs": int(jobs),
    "sealed_link_directories": sealed_link_dirs.split(":"),
    "platform": platform.platform(),
    "calo_calib": {"path": calo_calib, "sha256": digest(calo_calib)},
    "ppg12_source": {
        "repository": ppg_repo,
        "revision": ppg_revision,
        "subtree": "anatreemaker/source",
        "working_tree_ignored": True,
        "rebuilt_against_common_runtime": True,
    },
    "custom_photon_builder": {
        "class": "PPG12OraclePhotonClusterBuilder",
        "compiled_into": "libRecoilJets.so",
        "release_photon_builder_symbol_collision": False,
    },
    "sources": [
        {"path": path, "sha256": digest(path)} for path in source_paths
    ],
    "release_copies": {
        "libcalo_reco.so": {
            "source": calo_source,
            "source_realpath": os.path.realpath(calo_source),
            "sha256": digest(calo_source),
            "contains_release_calibration_and_cluster_template": True,
        },
        "libclusteriso.so": {
            "source": clusteriso_source,
            "source_realpath": os.path.realpath(clusteriso_source),
            "sha256": digest(clusteriso_source),
        },
        "libjetbase.so": {
            "source": jetbase_source,
            "source_realpath": os.path.realpath(jetbase_source),
            "sha256": digest(jetbase_source),
        },
    },
    "staged_rewrites": {
        "exact_count_enforced": True,
        "checkout_absolute_includes_removed": True,
        "macro_private_routes_removed": True,
        "staged_configure_werror_disabled": True,
        "target_linker_inherits_am_ldflags": True,
        "calo_truth_evaluator_linked_from_release": True,
        "calo_calib_resolved_by_sealed_release_path": True,
        "contained_recoiljets_tmp_paths_allowed": True,
        "full_calo_reco_rebuild": False,
        "custom_builder_renamed": True,
        "archived_ppg12_binary_reused": False,
        "ppg12_source_locked_rebuild": True,
    },
    "validation": {
        "forbidden_route_scan": "pass",
        "readelf": "pass",
        "ldd": "pass",
        "root_load_and_header_smoke": "pass",
    },
    "logs": [{"path": path, "sha256": digest(path)} for path in log_paths],
}
Path(receipt).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY

python3 - "$runtime_manifest" "$build_receipt" "$expected_offline" \
  "$runtime_wrapper" "$runtime_impl" \
  "${runtime_root}/lib/libRecoilJets.so" \
  "${runtime_root}/lib/libCaloAna24.so" \
  "${runtime_root}/lib/libcalo_reco.so" \
  "${runtime_root}/lib/libclusteriso.so" \
  "${runtime_root}/lib/libjetbase.so" \
  "${runtime_root}/include/caloreco/PhotonClusterBuilder.h" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

(
    manifest, receipt, offline_main, macro, impl, recoil, ppg, calo,
    clusteriso, jetbase, photon_header,
) = sys.argv[1:]

def digest(path: str) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

roles = [
    ("recoil_macro", macro),
    ("recoil_impl", impl),
    ("libRecoilJets.so", recoil),
    ("libCaloAna24.so", ppg),
    ("libcalo_reco.so", calo),
    ("libclusteriso.so", clusteriso),
    ("libjetbase.so", jetbase),
    ("PhotonClusterBuilder.h", photon_header),
]
data = {
    "schema_version": 1,
    "runtime_profile": "new.17",
    "offline_main": offline_main,
    "isolated_build": True,
    "build_receipt": receipt,
    "build_receipt_sha256": digest(receipt),
    "files": [
        {"role": role, "path": path, "sha256": digest(path)}
        for role, path in roles
    ],
}
Path(manifest).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY

printf 'PASS\n' > "$state_file"
completed=1
trap - EXIT
printf 'PPG12_ORACLE_NEW17_BUILD_PASS output=%s manifest=%s receipt=%s\n' \
  "$output_dir" "$runtime_manifest" "$build_receipt"
