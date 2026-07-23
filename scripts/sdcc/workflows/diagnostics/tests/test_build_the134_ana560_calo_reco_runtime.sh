#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../../.." && pwd -P)"
production_builder="${repo_root}/scripts/sdcc/workflows/diagnostics/build_the134_ana560_calo_reco_runtime.sh"
builder="$production_builder"
abi_allowlist="${repo_root}/scripts/sdcc/workflows/diagnostics/contracts/the134_ana560_calo_reco_abi_additions.allowlist"
readonly pinned_commit="cba274033b5560e32600cdeaa7676b6ab4a6c971"

fail() {
  printf 'TEST_FAIL: %s\n' "$*" >&2
  exit 1
}

sha256_file() {
  python3 - "$1" <<'PY'
from pathlib import Path
import hashlib
import sys
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
PY
}

[[ -f "$builder" ]] || fail "builder is missing"
[[ -f "$abi_allowlist" ]] || fail "durable ABI additions allowlist is missing"
bash -n "$builder"
production_help="$(bash "$production_builder" --help)"
grep -Fq '/opt/sphenix/core/bin/sphenix_setup.sh' <<<"$production_help" || \
  fail "production help omits the frozen setup authority"
grep -Fq '/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560' \
  <<<"$production_help" || fail "production help omits the frozen OFFLINE_MAIN authority"
grep -Fq "$pinned_commit" <<<"$production_help" || \
  fail "production help omits the frozen coresoftware authority"

tmp="$(mktemp -d /tmp/the134_ana560_calo_builder_test.XXXXXX)"
mutable_dependency="/tmp/the134_mutable_dependency_${$}.so"
trap 'chmod -R u+w "$tmp" 2>/dev/null || true; rm -rf "$tmp"; rm -f "$mutable_dependency"' EXIT

fixture_repo="${tmp}/coresoftware"
fixture_package="${fixture_repo}/offline/packages/CaloReco"
mkdir -p "$fixture_package"
git -C "$fixture_repo" init -q
git -C "$fixture_repo" config user.name "THE-134 Test"
git -C "$fixture_repo" config user.email "the134-test@example.invalid"

for file in configure.ac Makefile.am RawClusterBuilderTopo.h; do
  printf 'fixture %s\n' "$file" >"${fixture_package}/${file}"
done
cat >"${fixture_package}/autogen.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
prefix=""
for argument in "$@"; do
  case "$argument" in
    --prefix=*) prefix="${argument#--prefix=}" ;;
  esac
done
[[ -n "$prefix" ]]
printf '%s\n' "$prefix" >.fake_install_prefix
: >Makefile
EOF
printf 'archived photon implementation\n' >"${fixture_package}/PhotonClusterBuilder.cc"
printf 'archived photon header\n' >"${fixture_package}/PhotonClusterBuilder.h"
cat >"${fixture_package}/RawClusterBuilderTopo.cc" <<'EOF'
#include <calobase/TowerInfoContainer.h>  // for TowerInfoContainer
#include <algorithm>
#include <iostream>
#include <utility>

namespace
{
bool sort_by_pair_second(const std::pair<int, float> &a, const std::pair<int, float> &b)
{
  return (a.second > b.second);
}
}

int RawClusterBuilderTopo::process_event(PHCompositeNode *topNode)
{
  TowerInfoContainer *towerinfosEM = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainer *towerinfosIH = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towerinfosOH = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  if (!towerinfosEM)
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERINFO_CALIB_CEMC does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!towerinfosIH)
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERINFO_CALIB_HCALIN does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!towerinfosOH)
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERINFO_CALIB_HCALOUT does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_containers[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  unsigned int iEM = 0;
  unsigned int iIH = 0;
  unsigned int iOH = 0;
  {
      unsigned int towerinfo_key = towerinfosEM->encode_key(iEM);
      int ti_ieta = towerinfosEM->getTowerEtaBin(towerinfo_key);
      int ti_iphi = towerinfosEM->getTowerPhiBin(towerinfo_key);
  }
  {
      unsigned int towerinfo_key = towerinfosIH->encode_key(iIH);
      int ti_ieta = towerinfosIH->getTowerEtaBin(towerinfo_key);
      int ti_iphi = towerinfosIH->getTowerPhiBin(towerinfo_key);
  }
  {
      unsigned int towerinfo_key = towerinfosOH->encode_key(iOH);
      int ti_ieta = towerinfosOH->getTowerEtaBin(towerinfo_key);
      int ti_iphi = towerinfosOH->getTowerPhiBin(towerinfo_key);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
EOF
git -C "$fixture_repo" add offline/packages/CaloReco
git -C "$fixture_repo" commit -q -m "fixture CaloReco"
commit="$(git -C "$fixture_repo" rev-parse HEAD)"

overlay="${tmp}/overlay"
mkdir "$overlay"
printf 'authoritative H70 PhotonClusterBuilder implementation\n' \
  >"${overlay}/PhotonClusterBuilder.cc"
printf 'authoritative H70 PhotonClusterBuilder header\n' \
  >"${overlay}/PhotonClusterBuilder.h"
cc_sha="$(sha256_file "${overlay}/PhotonClusterBuilder.cc")"
h_sha="$(sha256_file "${overlay}/PhotonClusterBuilder.h")"

fake_release="${tmp}/release/release_ana/ana.560"
fake_bin="${tmp}/fake_bin"
mkdir -p "${fake_release}/lib" "${fake_release}/include/calobase" "$fake_bin"
printf 'fake exact ana.560 baseline CaloReco\n' \
  >"${fake_release}/lib/libcalo_reco.so.0"
printf 'fake exact ana.560 CaloIO\n' >"${fake_release}/lib/libcalo_io.so"
printf 'fake exact ana.560 dependency\n' >"${fake_release}/lib/libfake_dep.so"
printf 'fake TowerInfoDefs authority\n' \
  >"${fake_release}/include/calobase/TowerInfoDefs.h"
printf 'fake TowerInfoContainer authority\n' \
  >"${fake_release}/include/calobase/TowerInfoContainer.h"
printf 'fake TowerInfoContainerv1 authority\n' \
  >"${fake_release}/include/calobase/TowerInfoContainerv1.h"
baseline_sha="$(sha256_file "${fake_release}/lib/libcalo_reco.so.0")"
allowlist_sha="$(sha256_file "$abi_allowlist")"
defs_sha="$(sha256_file "${fake_release}/include/calobase/TowerInfoDefs.h")"
container_sha="$(sha256_file "${fake_release}/include/calobase/TowerInfoContainer.h")"
containerv1_sha="$(sha256_file "${fake_release}/include/calobase/TowerInfoContainerv1.h")"
calo_io_sha="$(sha256_file "${fake_release}/lib/libcalo_io.so")"

cat >"${tmp}/fake_setup.sh" <<EOF
#!/usr/bin/env bash
export OFFLINE_MAIN="${fake_release}"
export PATH="${fake_bin}:/usr/bin:/bin"
EOF

for tool in aclocal automake autoconf libtoolize; do
  cat >"${fake_bin}/${tool}" <<'EOF'
#!/usr/bin/env bash
exit 0
EOF
done

cat >"${fake_bin}/make" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if [[ " \$* " == *" install "* ]]; then
  prefix="\$(cat .fake_install_prefix)"
  mkdir -p "\${prefix}/lib" "\${prefix}/include/caloreco"
  printf 'fake candidate CaloReco runtime\n' >"\${prefix}/lib/libcalo_reco.so.0.0.0"
  ln -s libcalo_reco.so.0.0.0 "\${prefix}/lib/libcalo_reco.so.0"
  ln -s libcalo_reco.so.0 "\${prefix}/lib/libcalo_reco.so"
  cp "${overlay}/PhotonClusterBuilder.h" \
    "\${prefix}/include/caloreco/PhotonClusterBuilder.h"
fi
exit 0
EOF

cat >"${fake_bin}/readelf" <<EOF
#!/usr/bin/env bash
set -euo pipefail
library="\${@: -1}"
if [[ " \$* " == *" --dyn-syms "* ]]; then
  cat <<'OUT'
Symbol table '.dynsym' contains entries:
   Num:    Value          Size Type    Bind   Vis      Ndx Name
     1: 0000000000000000    16 FUNC    GLOBAL DEFAULT   12 _ZN20PhotonClusterBuilder13process_eventEP15PHCompositeNode
     2: 0000000000000000    16 FUNC    GLOBAL DEFAULT   12 _ZN21RawClusterBuilderTopo13process_eventEP15PHCompositeNode
OUT
  if [[ "\$library" != *"/release/"* ]]; then
    awk '
      {
        sub(/[[:space:]]*#.*/, "")
        gsub(/^[[:space:]]+|[[:space:]]+$/, "")
        if (length(\$0) > 0) {
          print " " NR + 2 ": 0000000000000000 16 FUNC GLOBAL DEFAULT 12 " \$1
        }
      }
    ' "${abi_allowlist}"
  fi
elif [[ " \$* " == *" --version-info "* ]]; then
  printf '%s\\n' 'Version needs section:' \
    '  0x0000: Version: 1  File: libfake_dep.so  Cnt: 1' \
    '  0x0010:   Name: THE134_FAKE_ABI_1  Flags: none  Version: 2'
else
  cat <<'OUT'
Dynamic section at offset 0:
 0x0000000000000001 (NEEDED)             Shared library: [libfake_dep.so]
 0x000000000000000e (SONAME)             Library soname: [libcalo_reco.so.0]
OUT
fi
EOF

cat >"${fake_bin}/ldd" <<EOF
#!/usr/bin/env bash
library="\${@: -1}"
if [[ "\$library" == *"dependency_failure"* ]]; then
  printf 'libfake_dep.so => ${mutable_dependency} (0x1)\\n'
else
  printf 'libfake_dep.so => ${fake_release}/lib/libfake_dep.so (0x1)\\n'
fi
EOF
printf 'mutable dependency that must be rejected\n' >"$mutable_dependency"

cat >"${fake_bin}/nm" <<EOF
#!/usr/bin/env bash
library="\${@: -1}"
printf '%s\n' '_ZN20PhotonClusterBuilder13process_eventEP15PHCompositeNode T 0 1'
if [[ "\$library" != *"abi_failure"* || "\$library" == *"/release/"* ]]; then
  printf '%s\n' '_ZN21RawClusterBuilderTopo13process_eventEP15PHCompositeNode T 1 1'
fi
if [[ "\$library" != *"/release/"* ]]; then
  awk '
    {
      sub(/[[:space:]]*#.*/, "")
      gsub(/^[[:space:]]+|[[:space:]]+$/, "")
      if (length(\$0) > 0) print
    }
  ' "${abi_allowlist}"
fi
if [[ "\$library" == *"abi_addition_failure"* ]]; then
  printf '%s\n' '_ZTHE134UnexpectedExportv T 2 1'
fi
EOF

cat >"${fake_bin}/c++filt" <<'EOF'
#!/usr/bin/env bash
sed \
  -e 's/_ZN20PhotonClusterBuilder13process_eventEP15PHCompositeNode/PhotonClusterBuilder::process_event(PHCompositeNode*)/' \
  -e 's/_ZN21RawClusterBuilderTopo13process_eventEP15PHCompositeNode/RawClusterBuilderTopo::process_event(PHCompositeNode*)/'
EOF

cat >"${fake_bin}/root-config" <<'EOF'
#!/usr/bin/env bash
case "${1:-}" in
  --cflags|--libs) printf '\n' ;;
  *) exit 1 ;;
esac
EOF

cat >"${fake_bin}/c++" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "--version" ]]; then
  printf 'fake-c++ ana.560 deterministic toolchain\n'
  exit 0
fi
output=""
while (($#)); do
  if [[ "$1" == "-o" ]]; then
    output="$2"
    shift 2
  else
    shift
  fi
done
[[ -n "$output" ]]
if [[ "$output" == *"single_provider_probe" ]]; then
  if [[ "$output" == *"provider_failure"* ]]; then
    printf '%s\n' '#!/usr/bin/env bash' 'exit 7' >"$output"
  else
    cat >"$output" <<'PROBE'
#!/usr/bin/env bash
printf 'THE134_SINGLE_PROVIDER_PROBE_PASS provider=%s provider_count=1\n' "$1"
PROBE
  fi
elif [[ "$output" == *"towerinfo_mapping_equivalence" ]]; then
  cat >"$output" <<'MAPPING'
#!/usr/bin/env bash
printf 'THE134_TOWERINFO_MAPPING_EQUIVALENCE_PASS EMCAL=24576 HCAL=1536\n'
MAPPING
elif [[ "$output" == *"towerinfo_role_size_matrix" ]]; then
  cat >"$output" <<'MATRIX'
#!/usr/bin/env bash
printf 'THE134_TOWERINFO_ROLE_SIZE_MATRIX_PASS typed=2 wrong_valid=2 legacy_exact=2 wrong_size=2\n'
MATRIX
else
  exit 9
fi
chmod 0755 "$output"
EOF

cat >"${fake_bin}/root" <<'EOF'
#!/usr/bin/env bash
printf 'THE134_ROOT_LOAD_PROVIDER_PASS provider=%s load_rc=0 preload_provider_count=0 provider_count=1\n' \
  "${THE134_CALORECO_LIBRARY}"
EOF
chmod 0755 "${tmp}/fake_setup.sh" "${fake_bin}"/*

# Exercise all build mechanics against a temporary byte-derived test copy whose
# two absolute filesystem authorities are mechanically substituted.  The real
# builder has no test mode or override and is separately checked above.
fixture_builder="${tmp}/build_the134_ana560_calo_reco_runtime.fixture.sh"
cp "$production_builder" "$fixture_builder"
python3 - "$fixture_builder" "$fake_release" "${tmp}/fake_setup.sh" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
replacements = {
    "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560":
        sys.argv[2],
    "/opt/sphenix/core/bin/sphenix_setup.sh": sys.argv[3],
}
for old, new in replacements.items():
    if old not in text:
        raise SystemExit(f"missing frozen authority in production builder: {old}")
    text = text.replace(old, new)
path.write_text(text, encoding="utf-8")
PY
chmod 0755 "$fixture_builder"
builder="$fixture_builder"

# The production builder has no commit override.  The local unit fixture uses
# a Git shim solely to expose its current synthetic HEAD under the frozen
# commit identity.  A real source-stage/build must contain the actual cba274
# object and never uses this shim.
git_shim_dir="${tmp}/git_shim"
mkdir "$git_shim_dir"
cat >"${git_shim_dir}/git" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
readonly expected="cba274033b5560e32600cdeaa7676b6ab4a6c971"
target="${THE134_TEST_GIT_TARGET:?missing fixture Git target}"
args=("$@")
if [[ " ${args[*]} " == *" rev-parse ${expected}^{commit} "* ]]; then
  printf '%s\n' "$expected"
  exit 0
fi
for index in "${!args[@]}"; do
  if [[ "${args[$index]}" == "$expected"* ]]; then
    args[$index]="${target}${args[$index]#"$expected"}"
  fi
done
exec /usr/bin/git "${args[@]}"
EOF
chmod 0755 "${git_shim_dir}/git"
export PATH="${git_shim_dir}:${PATH}"
export THE134_TEST_GIT_TARGET="$commit"

common_args=(
  --coresoftware-repo "$fixture_repo"
  --coresoftware-commit "$pinned_commit"
  --photon-builder-cc "${overlay}/PhotonClusterBuilder.cc"
  --photon-builder-cc-sha256 "$cc_sha"
  --photon-builder-h "${overlay}/PhotonClusterBuilder.h"
  --photon-builder-h-sha256 "$h_sha"
  --baseline-calo-reco "${fake_release}/lib/libcalo_reco.so.0"
  --baseline-calo-reco-sha256 "$baseline_sha"
  --abi-additions-allowlist "$abi_allowlist"
  --abi-additions-allowlist-sha256 "$allowlist_sha"
  --towerinfo-defs-header "${fake_release}/include/calobase/TowerInfoDefs.h"
  --towerinfo-defs-header-sha256 "$defs_sha"
  --towerinfo-container-header "${fake_release}/include/calobase/TowerInfoContainer.h"
  --towerinfo-container-header-sha256 "$container_sha"
  --towerinfo-containerv1-header "${fake_release}/include/calobase/TowerInfoContainerv1.h"
  --towerinfo-containerv1-header-sha256 "$containerv1_sha"
  --calo-io-library "${fake_release}/lib/libcalo_io.so"
  --calo-io-library-sha256 "$calo_io_sha"
  --setup-script "${tmp}/fake_setup.sh"
  --expected-offline-main "$fake_release"
  --jobs 2
)

# Exercise the unmodified production pins directly.  These checks stop before
# any Git or filesystem authority is needed and must not create output roots.
production_wrong_setup="${tmp}/production_wrong_setup"
if bash "$production_builder" --output-root "$production_wrong_setup" \
    "${common_args[@]}" \
    --expected-offline-main \
    /cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560 \
    >/dev/null 2>&1; then
  fail "production builder accepted a caller-selected setup script"
fi
[[ ! -e "$production_wrong_setup" ]] || \
  fail "production wrong-setup check mutated the output root"

production_wrong_prefix="${tmp}/production_wrong_prefix"
if bash "$production_builder" --output-root "$production_wrong_prefix" \
    "${common_args[@]}" \
    --setup-script /opt/sphenix/core/bin/sphenix_setup.sh \
    >/dev/null 2>&1; then
  fail "production builder accepted a noncanonical ana.560 lookalike prefix"
fi
[[ ! -e "$production_wrong_prefix" ]] || \
  fail "production wrong-prefix check mutated the output root"

output_one="${tmp}/one"
plan_one="$(bash "$builder" --output-root "$output_one" "${common_args[@]}")"
grep -q '^THE134_ANA560_CALORECO_BUILD_PLAN$' <<<"$plan_one" || \
  fail "plan marker is missing"
grep -q 'package_export: git archive offline/packages/CaloReco' <<<"$plan_one" || \
  fail "git archive contract is missing"
grep -q 'one full libcalo_reco.so.0 provides PhotonClusterBuilder and RawClusterBuilderTopo' \
  <<<"$plan_one" || fail "single-provider contract is missing"
grep -q "pinned_coresoftware_commit: ${pinned_commit}" <<<"$plan_one" || \
  fail "frozen coresoftware authority is missing from the plan"
grep -q "$pinned_commit" < <(bash "$builder" --help) || \
  fail "fixture help does not declare the frozen coresoftware authority"
source_token_one="$(sed -n 's/^  source_token: //p' <<<"$plan_one")"
build_token_one="$(sed -n 's/^  build_token: //p' <<<"$plan_one")"
[[ "$source_token_one" =~ ^the134-ana560-caloreco-source:[0-9a-f]{64}$ ]] || \
  fail "source token is malformed"
[[ "$build_token_one" =~ ^the134-ana560-caloreco-build:[0-9a-f]{64}$ ]] || \
  fail "build token is malformed"
[[ ! -e "$output_one" ]] || fail "plan mode mutated the output root"

bash "$builder" --stage-source --token "$source_token_one" \
  --output-root "$output_one" "${common_args[@]}" >/dev/null
[[ -f "${output_one}/SOURCE_ONLY_NOT_RUNTIME_AUTHORITY" ]] || \
  fail "source-only authority warning is missing"
[[ ! -e "${output_one}/build_receipt.json" ]] || \
  fail "source-only staging emitted runtime authority"

staged="${output_one}/source/offline/packages/CaloReco"
cmp -s "${overlay}/PhotonClusterBuilder.cc" "${staged}/PhotonClusterBuilder.cc" || \
  fail "PhotonClusterBuilder.cc overlay differs"
cmp -s "${overlay}/PhotonClusterBuilder.h" "${staged}/PhotonClusterBuilder.h" || \
  fail "PhotonClusterBuilder.h overlay differs"
grep -Fq '#include <calobase/TowerInfoDefs.h>       // detector-explicit channel maps' \
  "${staged}/RawClusterBuilderTopo.cc" || fail "TowerInfoDefs include was not added"
grep -Fq 'TowerInfoDefs::encode_emcal(iEM)' "${staged}/RawClusterBuilderTopo.cc" || \
  fail "CEMC detector-explicit map is missing"
[[ "$(grep -F -c 'TowerInfoDefs::encode_hcal(' "${staged}/RawClusterBuilderTopo.cc")" == "2" ]] || \
  fail "HCAL detector-explicit map count differs"
grep -Fq 'validate_towerinfo_role_and_size(' \
  "${staged}/RawClusterBuilderTopo.cc" || fail "role/size guard helper is missing"
grep -Fq 'towerinfosEM, TowerInfoContainer::EMCAL, 24576U,' \
  "${staged}/RawClusterBuilderTopo.cc" || fail "CEMC role/size guard differs"
[[ "$(grep -F -c 'TowerInfoContainer::HCAL, 1536U' \
  "${staged}/RawClusterBuilderTopo.cc")" == "2" ]] || \
  fail "HCAL role/size guard count differs"
grep -Fq 'actual_detector == invalid_detector' \
  "${staged}/THE134TowerInfoContractV1.h" || \
  fail "legacy DETECTOR_INVALID handling is missing"
[[ "$(grep -F -c 'static bool ' "${staged}/RawClusterBuilderTopo.cc")" == "3" ]] || \
  fail "legacy warning once-per-role state count differs"
grep -Fq 'diagnostic_frequency=once_per_role_process' \
  "${staged}/RawClusterBuilderTopo.cc" || \
  fail "legacy warning frequency marker is missing"
grep -Fq 'return Fun4AllReturnCodes::ABORTRUN;' \
  "${staged}/RawClusterBuilderTopo.cc" || fail "role/size failure is not fatal"
if grep -Fq -- '->encode_key(' "${staged}/RawClusterBuilderTopo.cc"; then
  fail "container-selected channel map survived"
fi
python3 - "${staged}/RawClusterBuilderTopo.cc" \
  "${staged}/THE134TowerInfoContractV1.h" <<'PY'
from pathlib import Path
import sys

text = Path(sys.argv[1]).read_text(encoding="utf-8")
predicate = Path(sys.argv[2]).read_text(encoding="utf-8")
size_check = predicate.index("if (actual_size != expected_size)")
typed_accept = predicate.index("if (actual_detector == expected_detector)")
legacy_accept = predicate.index("if (actual_detector == invalid_detector)")
final_reject = predicate.rindex("return TowerInfoContractResult::REJECT;")
wrong_valid_failure = text.index(
    'std::cerr << "THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_CONTRACT_FAIL"',
)
role_calls = text.index("  if (!validate_towerinfo_role_and_size(")
geometry_lookup = text.index("  _geom_containers[0] =")
first_mapping = text.index("TowerInfoDefs::encode_emcal(iEM)")
assert size_check < typed_accept < legacy_accept < final_reject
assert wrong_valid_failure < role_calls
assert role_calls < geometry_lookup < first_mapping
PY

python3 - "${output_one}/source_manifest.json" "$pinned_commit" "$cc_sha" "$h_sha" \
  "${staged}/THE134TowerInfoContractV1.h" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

path = Path(sys.argv[1])
payload = json.loads(path.read_text(encoding="utf-8"))
assert payload["schema"] == "THE134_ANA560_CALORECO_SOURCE_MANIFEST_V2"
assert payload["coresoftware"]["commit"] == sys.argv[2]
assert payload["overlay"]["PhotonClusterBuilder.cc"]["staged_sha256"] == sys.argv[3]
assert payload["overlay"]["PhotonClusterBuilder.h"]["staged_sha256"] == sys.argv[4]
assert payload["mapping_patch"]["id"] == (
    "THE134_RAWCLUSTERBUILDERTOPO_DETECTOR_EXPLICIT_CHANNEL_MAP_V1"
)
assert payload["mapping_patch"]["replacement_counts"] == {
    "CEMC": 1,
    "HCALIN": 1,
    "HCALOUT": 1,
}
assert payload["mapping_patch"]["changed_scientific_controls"] == []
assert payload["role_size_guard"] == {
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
        "failure_marker": "THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_CONTRACT_FAIL",
        "legacy_invalid_marker": (
            "THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_LEGACY_INVALID_ACCEPTED"
        ),
    },
    "failure_return": "Fun4AllReturnCodes::ABORTRUN",
    "id": "THE134_RAWCLUSTERBUILDERTOPO_TOWERINFO_ROLE_SIZE_GUARD_V1",
    "predicate_header": {
        "path": "THE134TowerInfoContractV1.h",
        "sha256": hashlib.sha256(Path(sys.argv[5]).read_bytes()).hexdigest(),
    },
    "valid_path_effect": (
        "none; validation precedes geometry lookup and tower arithmetic"
    ),
    "wrong_valid_detector": "rejected",
}
canonical = json.dumps(
    payload, indent=2, sort_keys=True, separators=(",", ": ")
) + "\n"
assert path.read_text(encoding="utf-8") == canonical
PY
python3 - "$output_one" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
for path in (root, *root.rglob("*")):
    assert path.stat().st_mode & 0o222 == 0, f"writable sealed path: {path}"
PY

manifest_one_sha="$(sha256_file "${output_one}/source_manifest.json")"

# A mutable checkout edit must not affect a commit-archive build.
printf 'uncommitted mutation that must never enter the source stage\n' \
  >>"${fixture_package}/RawClusterBuilderTopo.h"
output_two="${tmp}/two"
plan_two="$(bash "$builder" --output-root "$output_two" "${common_args[@]}")"
source_token_two="$(sed -n 's/^  source_token: //p' <<<"$plan_two")"
bash "$builder" --stage-source --token "$source_token_two" \
  --output-root "$output_two" "${common_args[@]}" >/dev/null
manifest_two_sha="$(sha256_file "${output_two}/source_manifest.json")"
[[ "$manifest_two_sha" == "$manifest_one_sha" ]] || \
  fail "mutable checkout contents changed the deterministic source manifest"

# A changed mapping anchor must be rejected rather than partially generalized.
python3 - "${fixture_package}/RawClusterBuilderTopo.cc" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
old = "      unsigned int towerinfo_key = towerinfosOH->encode_key(iOH);\n"
assert text.count(old) == 1
path.write_text(text.replace(old, "      unsigned int changed_mapping_anchor = iOH;\n"), encoding="utf-8")
PY
git -C "$fixture_repo" add offline/packages/CaloReco/RawClusterBuilderTopo.cc
git -C "$fixture_repo" commit -q -m "mutate one mapping anchor"
mutated_commit="$(git -C "$fixture_repo" rev-parse HEAD)"
export THE134_TEST_GIT_TARGET="$mutated_commit"
mutated_output="${tmp}/mutated_mapping"
mutated_args=("${common_args[@]}")
mutated_plan="$(bash "$builder" --output-root "$mutated_output" "${mutated_args[@]}")"
mutated_token="$(sed -n 's/^  source_token: //p' <<<"$mutated_plan")"
if bash "$builder" --stage-source --token "$mutated_token" \
    --output-root "$mutated_output" "${mutated_args[@]}" >/dev/null 2>&1; then
  fail "changed RawClusterBuilderTopo mapping anchor was accepted"
fi
[[ ! -e "${mutated_output}/source_manifest.json" ]] || \
  fail "rejected mapping mutation emitted a source manifest"
export THE134_TEST_GIT_TARGET="$commit"

if bash "$builder" --output-root "${tmp}/bad_hash" "${common_args[@]}" \
    --photon-builder-cc-sha256 \
    0000000000000000000000000000000000000000000000000000000000000000 \
    >/dev/null 2>&1; then
  fail "wrong overlay hash was accepted"
fi
[[ ! -e "${tmp}/bad_hash" ]] || fail "wrong overlay hash mutated an output root"

cp "$abi_allowlist" "${tmp}/mutated_allowlist"
printf '%s\n' '_ZTHE134UnauthorizedContractEditv T' \
  >>"${tmp}/mutated_allowlist"
mutated_allowlist_sha="$(sha256_file "${tmp}/mutated_allowlist")"
if bash "$builder" --output-root "${tmp}/mutated_allowlist_output" \
    "${common_args[@]}" \
    --abi-additions-allowlist "${tmp}/mutated_allowlist" \
    --abi-additions-allowlist-sha256 "$mutated_allowlist_sha" \
    >/dev/null 2>&1; then
  fail "self-consistent but unauthorized ABI allowlist was accepted"
fi
[[ ! -e "${tmp}/mutated_allowlist_output" ]] || \
  fail "unauthorized ABI allowlist mutated an output root"

if bash "$builder" --output-root "${tmp}/short_commit" \
    "${common_args[@]}" --coresoftware-commit "${commit:0:12}" \
    >/dev/null 2>&1; then
  fail "short coresoftware revision was accepted"
fi

wrong_commit_output="${tmp}/wrong_full_commit"
if bash "$builder" --output-root "$wrong_commit_output" \
    "${common_args[@]}" --coresoftware-commit "$commit" \
    >/dev/null 2>&1; then
  fail "non-authoritative full coresoftware commit was accepted"
fi
[[ ! -e "$wrong_commit_output" ]] || \
  fail "wrong full coresoftware commit mutated the output root"

mkdir "${tmp}/already_exists"
if bash "$builder" --output-root "${tmp}/already_exists" \
    "${common_args[@]}" >/dev/null 2>&1; then
  fail "existing output root was accepted"
fi

if bash "$builder" --output-root "${repo_root}/unsafe_the134_runtime" \
    "${common_args[@]}" >/dev/null 2>&1; then
  fail "unsafe repository output root was accepted"
fi

escape_output="/tmp/../tmp/$(basename "$tmp")_canonical_escape"
if bash "$builder" --output-root "$escape_output" \
    "${common_args[@]}" >/dev/null 2>&1; then
  fail "output root with a lexical escape component was accepted"
fi
[[ ! -e "/tmp/$(basename "$tmp")_canonical_escape" ]] || \
  fail "rejected lexical escape mutated its canonical output root"

wrong_action_output="${tmp}/wrong_action"
wrong_plan="$(bash "$builder" --output-root "$wrong_action_output" "${common_args[@]}")"
wrong_source_token="$(sed -n 's/^  source_token: //p' <<<"$wrong_plan")"
if bash "$builder" --build --token "$wrong_source_token" \
    --output-root "$wrong_action_output" "${common_args[@]}" >/dev/null 2>&1; then
  fail "source-stage token authorized a build"
fi
[[ ! -e "$wrong_action_output" ]] || fail "wrong-action token mutated the output root"

successful_build="${tmp}/successful_build"
successful_plan="$(
  bash "$builder" --output-root "$successful_build" "${common_args[@]}"
)"
successful_token="$(sed -n 's/^  build_token: //p' <<<"$successful_plan")"
successful_log="${tmp}/successful_build.stdout"
if ! bash "$builder" --build --token "$successful_token" \
    --output-root "$successful_build" "${common_args[@]}" \
    >"$successful_log" 2>"${successful_log}.stderr"; then
  cat "${successful_log}.stderr" >&2
  [[ -f "${successful_build}/logs/build.log" ]] && \
    cat "${successful_build}/logs/build.log" >&2
  [[ -f "${successful_build}/logs/validation.log" ]] && \
    cat "${successful_build}/logs/validation.log" >&2
  fail "successful fake build path failed"
fi
grep -Fq 'THE134_ANA560_CALORECO_BUILD_PASS' "$successful_log" || \
  fail "successful fake build did not reach PASS"
[[ -f "${successful_build}/build_receipt.json" ]] || \
  fail "successful fake build omitted its receipt"
grep -Fq 'THE134_SINGLE_PROVIDER_PROBE_PASS' \
  "${successful_build}/evidence/single_provider_probe.txt" || \
  fail "single-provider executable proof is missing"
grep -Fq 'THE134_ROOT_LOAD_PROVIDER_PASS' \
  "${successful_build}/evidence/root_load_smoke.txt" || \
  fail "ROOT provider proof is missing"
grep -Fq 'THE134_TOWERINFO_MAPPING_EQUIVALENCE_PASS EMCAL=24576 HCAL=1536' \
  "${successful_build}/evidence/towerinfo_mapping_equivalence.txt" || \
  fail "exhaustive mapping proof is missing"
grep -Fq \
  'THE134_TOWERINFO_ROLE_SIZE_MATRIX_PASS typed=2 wrong_valid=2 legacy_exact=2 wrong_size=2' \
  "${successful_build}/evidence/towerinfo_role_size_matrix.txt" || \
  fail "executable role/size matrix proof is missing"

python3 - "${successful_build}/build_receipt.json" \
  "${successful_build}/source_manifest.json" "$baseline_sha" \
  "$allowlist_sha" "$defs_sha" "$container_sha" "$containerv1_sha" \
  "$calo_io_sha" "$successful_token" "$successful_build" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

receipt_path = Path(sys.argv[1])
source_manifest_path = Path(sys.argv[2])
payload = json.loads(receipt_path.read_text(encoding="utf-8"))
assert payload["schema"] == "THE134_ANA560_CALORECO_BUILD_RECEIPT_V2"
assert payload["status"] == "PASS"
assert payload["build"]["coresoftware_commit"] == (
    "cba274033b5560e32600cdeaa7676b6ab4a6c971"
)
assert payload["build"]["coresoftware_commit_authority"] == (
    "cba274033b5560e32600cdeaa7676b6ab4a6c971"
)
assert payload["build"]["coresoftware_commit_pinned"] is True
assert payload["build"]["offline_main_authority"] == (
    payload["runtime"]["offline_main"]
)
assert payload["build"]["offline_main_pinned"] is True
assert Path(payload["build"]["setup_script_authority"]).name == "fake_setup.sh"
assert payload["build"]["setup_script_pinned"] is True
assert payload["source_manifest"]["sha256"] == hashlib.sha256(
    source_manifest_path.read_bytes()
).hexdigest()
assert payload["abi"]["status"] == "PASS"
assert payload["abi"]["baseline"]["library"]["sha256"] == sys.argv[3]
assert payload["abi"]["additions_allowlist"]["expected_sha256"] == sys.argv[4]
assert payload["abi"]["additions_allowlist"]["provenance"] == {
    "immutable_prior_override_sha256": (
        "5d4eca4abdaa274d308856e050d19b17e02bc62652a03dda6f0ea95719b9147e"
    ),
    "intentional_exclusions": [
        {
            "reason": (
                "detector-explicit RawClusterBuilderTopo mapping removes "
                "generic container-selected encode_key"
            ),
            "symbol": "_ZN18TowerInfoContainer10encode_keyEj W",
        }
    ],
}
assert payload["abi"]["removed_symbols"]["count"] == 0
assert payload["abi"]["removed_dynsym_metadata"]["count"] == 0
assert payload["abi"]["added_symbols"]["count"] == 28
assert payload["abi"]["needed_exact_match"] is True
assert payload["abi"]["rpath_runpath_absent"] is True
assert payload["abi"]["version_inventory_exact_match"] is True
assert payload["abi"]["baseline"]["version_inventory"]["entries"] == [
    "THE134_FAKE_ABI_1"
]
assert payload["abi"]["candidate"]["version_inventory"]["entries"] == [
    "THE134_FAKE_ABI_1"
]
assert payload["abi"]["dictionary_inventory"]["baseline"]["entries"] == []
assert payload["abi"]["dictionary_inventory"]["candidate"]["entries"] == []
assert [entry["path"] for entry in payload["abi"]["declared_implementation_replacements"]] == [
    "PhotonClusterBuilder.cc",
    "RawClusterBuilderTopo.cc",
]
assert payload["single_provider"]["provider_probe"]["status"] == "PASS"
assert payload["single_provider"]["provider_probe"]["preload_provider_count"] == 0
assert payload["single_provider"]["provider_probe"]["provider_count"] == 1
assert payload["runtime"]["root_load"]["load_return_code"] == 0
assert payload["runtime"]["root_load"]["preload_provider_count"] == 0
assert payload["runtime"]["root_load"]["provider_count"] == 1
assert payload["runtime"]["root_load"]["status"] == "PASS"
assert payload["runtime"]["dependency_policy"]["status"] == "PASS"
assert payload["runtime"]["dependency_policy"]["resolved_paths"]["entries"] == [
    str((
        Path(payload["runtime"]["offline_main"])
        / "lib"
        / "libfake_dep.so"
    ).resolve())
]
equivalence = payload["mapping_patch"]["exhaustive_equivalence"]
assert equivalence["status"] == "PASS"
assert equivalence["CEMC_channels"] == 24576
assert equivalence["HCAL_channels"] == 1536
assert equivalence["round_trip"] is True
authority = payload["mapping_patch"]["authority"]
assert authority["TowerInfoDefs.h"]["sha256"] == sys.argv[5]
assert authority["TowerInfoContainer.h"]["sha256"] == sys.argv[6]
assert authority["TowerInfoContainerv1.h"]["sha256"] == sys.argv[7]
assert authority["CaloIO"]["sha256"] == sys.argv[8]
assert payload["role_size_guard"]["wrong_valid_detector"] == "rejected"
assert payload["role_size_guard"]["accepted"]["CEMC"]["size"] == 24576
assert payload["role_size_guard"]["accepted"]["HCALIN"]["size"] == 1536
assert payload["role_size_guard"]["accepted"]["HCALOUT"]["size"] == 1536
assert payload["role_size_guard"]["diagnostics"]["frequency"] == (
    "once_per_role_process"
)
assert payload["role_size_guard"]["executable_matrix"]["status"] == "PASS"
assert payload["role_size_guard"]["executable_matrix"]["cases"] == {
    "legacy_exact": 2,
    "typed": 2,
    "wrong_size": 2,
    "wrong_valid": 2,
}
assert payload["build"]["authorization"] == {
    "action": "build",
    "canonical_output_root": str(Path(sys.argv[10]).resolve()),
    "contract_sha256": sys.argv[9].split(":", 1)[1],
    "jobs": 2,
    "token": sys.argv[9],
}
canonical = json.dumps(
    payload, indent=2, sort_keys=True, separators=(",", ": ")
) + "\n"
assert receipt_path.read_text(encoding="utf-8") == canonical
PY
python3 - "$successful_build" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
for path in (root, *root.rglob("*")):
    assert path.stat().st_mode & 0o222 == 0, f"writable sealed path: {path}"
PY

# A missing baseline ABI symbol must fail before runtime authority is emitted.
abi_failure="${tmp}/abi_failure"
abi_plan="$(bash "$builder" --output-root "$abi_failure" "${common_args[@]}")"
abi_token="$(sed -n 's/^  build_token: //p' <<<"$abi_plan")"
if bash "$builder" --build --token "$abi_token" \
    --output-root "$abi_failure" "${common_args[@]}" >/dev/null 2>&1; then
  fail "candidate with a removed baseline ABI symbol was accepted"
fi
[[ -s "${abi_failure}/evidence/abi_removed_symbols.txt" ]] || \
  fail "ABI-negative first-bad evidence is missing"
[[ ! -e "${abi_failure}/build_receipt.json" ]] || \
  fail "ABI-negative build emitted runtime authority"

# An unallowlisted candidate export must fail the exact additions contract.
abi_addition_failure="${tmp}/abi_addition_failure"
abi_addition_plan="$(
  bash "$builder" --output-root "$abi_addition_failure" "${common_args[@]}"
)"
abi_addition_token="$(sed -n 's/^  build_token: //p' <<<"$abi_addition_plan")"
if bash "$builder" --build --token "$abi_addition_token" \
    --output-root "$abi_addition_failure" "${common_args[@]}" \
    >/dev/null 2>&1; then
  fail "candidate with an unallowlisted ABI addition was accepted"
fi
grep -Fq '_ZTHE134UnexpectedExportv T' \
  "${abi_addition_failure}/evidence/abi_added_symbols.txt" || \
  fail "ABI-addition first-bad evidence is missing"
[[ ! -e "${abi_addition_failure}/build_receipt.json" ]] || \
  fail "unallowlisted ABI-addition build emitted runtime authority"

# A provider probe failure must also stop receipt materialization.
provider_failure="${tmp}/provider_failure"
provider_plan="$(
  bash "$builder" --output-root "$provider_failure" "${common_args[@]}"
)"
provider_token="$(sed -n 's/^  build_token: //p' <<<"$provider_plan")"
if bash "$builder" --build --token "$provider_token" \
    --output-root "$provider_failure" "${common_args[@]}" >/dev/null 2>&1; then
  fail "candidate with a failing provider proof was accepted"
fi
[[ -f "${provider_failure}/evidence/single_provider_probe" ]] || \
  fail "provider-negative probe executable is missing"
[[ ! -e "${provider_failure}/build_receipt.json" ]] || \
  fail "provider-negative build emitted runtime authority"

# A resolved dependency outside the candidate, versioned ana.560 platform, or
# system-library roots must fail before receipt materialization.
dependency_failure="${tmp}/dependency_failure"
dependency_plan="$(
  bash "$builder" --output-root "$dependency_failure" "${common_args[@]}"
)"
dependency_token="$(sed -n 's/^  build_token: //p' <<<"$dependency_plan")"
if bash "$builder" --build --token "$dependency_token" \
    --output-root "$dependency_failure" "${common_args[@]}" >/dev/null 2>&1; then
  fail "candidate with a mutable dependency was accepted"
fi
grep -Fq "$mutable_dependency" \
  "${dependency_failure}/evidence/ldd.txt" || \
  fail "dependency-negative first-bad evidence is missing"
[[ ! -e "${dependency_failure}/build_receipt.json" ]] || \
  fail "dependency-negative build emitted runtime authority"

# A caller-selected setup script must fail in plan mode before creating the
# output root.  The fixture builder has the fake setup path mechanically pinned.
cat >"${tmp}/wrong_setup.sh" <<'EOF'
#!/usr/bin/env bash
export OFFLINE_MAIN=/wrong/release_ana/ana.560
EOF
wrong_setup_output="${tmp}/wrong_setup"
wrong_setup_args=("${common_args[@]}" --setup-script "${tmp}/wrong_setup.sh")
if bash "$builder" --output-root "$wrong_setup_output" \
    "${wrong_setup_args[@]}" >/dev/null 2>&1; then
  fail "caller-selected setup authority was accepted"
fi
[[ ! -e "$wrong_setup_output" ]] || \
  fail "wrong setup authority mutated the output root"

# A lookalike ana.560 prefix must likewise fail before setup or output mutation.
wrong_prefix_output="${tmp}/wrong_prefix"
wrong_prefix="${tmp}/lookalike/release/release_ana/ana.560"
wrong_prefix_args=(
  "${common_args[@]}"
  --expected-offline-main "$wrong_prefix"
)
if bash "$builder" --output-root "$wrong_prefix_output" \
    "${wrong_prefix_args[@]}" >/dev/null 2>&1; then
  fail "noncanonical ana.560 lookalike prefix was accepted"
fi
[[ ! -e "$wrong_prefix_output" ]] || \
  fail "wrong OFFLINE_MAIN authority mutated the output root"

for invariant in \
  'git "${git_safe[@]}" -C "$coresoftware_repo_real" archive --format=tar' \
  'EXPECTED_CORESOFTWARE_COMMIT="cba274033b5560e32600cdeaa7676b6ab4a6c971"' \
  'EXPECTED_SETUP_SCRIPT="/opt/sphenix/core/bin/sphenix_setup.sh"' \
  'EXPECTED_OFFLINE_MAIN="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560"' \
  '--coresoftware-commit must equal frozen THE-134 authority' \
  '--setup-script must equal frozen THE-134 authority' \
  '--expected-offline-main must equal frozen THE-134 authority' \
  'SOURCE_DATE_EPOCH=${commit_epoch}' \
  'ZERO_AR_DATE=1' \
  'source "$setup_script" -n "$release_name"' \
  'ffile-prefix-map=${output_root}=/the134-caloreco' \
  'fdebug-prefix-map=${output_root}=/the134-caloreco' \
  '/bin/bash "${work_source}/autogen.sh" --prefix="$install_root"' \
  'make -j "$jobs"' \
  'make install' \
  'normalize_symbols "$baseline_library"' \
  'abi_additions_allowlist.normalized.txt' \
  'abi_removed_symbols.txt' \
  'baseline_needed.txt' \
  'candidate_needed.txt' \
  'baseline_dynsym_metadata.normalized.txt' \
  'candidate_dynsym_metadata.normalized.txt' \
  'baseline_version_inventory.txt' \
  'candidate_version_inventory.txt' \
  'baseline_dictionary_inventory.txt' \
  'candidate_dictionary_inventory.txt' \
  'readelf -d "$runtime_library"' \
  'ldd "$runtime_library"' \
  'resolved_dependency_paths.txt' \
  'nm -D --defined-only "$runtime_library" | c++filt' \
  'installed PhotonClusterBuilder header' \
  'libcalo_reco.so and ${EXPECTED_SONAME} resolve to different providers' \
  'export THE134_CALORECO_LIBRARY="$runtime_library"' \
  'dependency outside immutable allowlist' \
  'PhotonClusterBuilder::process_event(PHCompositeNode*)' \
  'RawClusterBuilderTopo::process_event(PHCompositeNode*)' \
  'THE134_SINGLE_PROVIDER_PROBE_PASS' \
  'THE134_TOWERINFO_MAPPING_EQUIVALENCE_PASS' \
  'THE134_TOWERINFO_ROLE_SIZE_MATRIX_PASS' \
  'if (rc != 0)' \
  'THE134_ROOT_LOAD_PROVIDER_PASS' \
  'THE134_ANA560_CALORECO_BUILD_RECEIPT_V2' \
  '"status": "PASS"' \
  'chmod -R a-w "$output_root"'; do
  grep -Fq -- "$invariant" "$production_builder" || \
    fail "missing production build invariant: ${invariant}"
done

printf 'THE134_ANA560_CALORECO_BUILDER_TEST_PASS\n'
