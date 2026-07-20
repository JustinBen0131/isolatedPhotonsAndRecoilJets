#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../../.." && pwd -P)"
builder="${repo_root}/scripts/sdcc/workflows/diagnostics/build_ppg12_oracle_new17_runtime.sh"

fail() {
  printf 'TEST_FAIL: %s\n' "$*" >&2
  exit 1
}

[[ -f "$builder" ]] || fail "builder is missing"
bash -n "$builder"

tmp="$(mktemp -d /tmp/ppg12_new17_builder_test.XXXXXX)"
trap 'rm -rf "$tmp"' EXIT
output="${tmp}/fresh_runtime"

plan="$(bash "$builder" --output-dir "$output" --setup-script /usr/bin/true --jobs 2)"
grep -q '^PPG12_ORACLE_NEW17_BUILD_PLAN$' <<<"$plan" || fail "plan marker missing"
grep -q 'runtime: new.17' <<<"$plan" || fail "new.17 contract missing"
grep -q 'source-locked libCaloAna24.so plus libRecoilJets.so with renamed PPG12-oracle photon builder' <<<"$plan" || fail "build contract missing"
grep -q 'ppg_binary_mode: source_locked_rebuild' <<<"$plan" || fail "default source-rebuild mode missing"
grep -q 'ppg_source_revision: 1c0ff86bf0ebabfba63a1abc4512cbe59fe48e31' <<<"$plan" || fail "PPG12 source revision missing"
grep -q 'exact new.17 libcalo_reco.so, libclusteriso.so, libjetbase.so (cp -L)' <<<"$plan" || fail "release-copy contract missing"
token="$(sed -n 's/^  token: //p' <<<"$plan")"
[[ "$token" =~ ^ppg12-new17:[0-9a-f]{64}$ ]] || fail "plan token is malformed"
[[ ! -e "$output" ]] || fail "plan mode mutated the output path"

photon_source="${tmp}/photon_source"
mkdir "$photon_source"
cp "${repo_root}/src/PhotonClusterBuilder.cc" "$photon_source/PhotonClusterBuilder.cc"
cp "${repo_root}/src/PhotonClusterBuilder.h" "$photon_source/PhotonClusterBuilder.h"
override_plan="$(bash "$builder" --output-dir "$output" --setup-script /usr/bin/true \
  --jobs 2 --photon-source-dir "$photon_source")"
override_token="$(sed -n 's/^  token: //p' <<<"$override_plan")"
[[ "$override_token" != "$token" ]] || fail "photon-source override did not change the sealed token"
[[ ! -e "$output" ]] || fail "override plan mode mutated the output path"

origin_root="${tmp}/origin_runtime"
origin_receipt="${origin_root}/build_receipt.json"
origin_library="${origin_root}/runtime/lib/libCaloAna24.so"
origin_manifest="${origin_root}/runtime_manifest.json"
mkdir -p "$(dirname "$origin_library")"
python3 - "$origin_receipt" "$origin_library" "$origin_manifest" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

receipt_path, library_path, manifest_path = map(Path, sys.argv[1:])
library_path.write_bytes(b"sealed-attempt12-libCaloAna24-test-fixture\n")
receipt = {
    "schema_version": 1,
    "runtime_profile": "new.17",
    "offline_main": (
        "/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/"
        "release/release_new/new.17"
    ),
    "ppg12_source": {
        "revision": "1c0ff86bf0ebabfba63a1abc4512cbe59fe48e31",
        "working_tree_ignored": True,
        "rebuilt_against_common_runtime": True,
    },
    "staged_rewrites": {
        "archived_ppg12_binary_reused": False,
        "ppg12_source_locked_rebuild": True,
    },
}
receipt_path.write_text(json.dumps(receipt, sort_keys=True) + "\n")
digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
manifest = {
    "schema_version": 1,
    "runtime_profile": "new.17",
    "offline_main": receipt["offline_main"],
    "isolated_build": True,
    "estimator_revision": "29f8223bd9b36dffab07961b597afa94185bbdf1",
    "build_receipt": str(receipt_path.resolve()),
    "build_receipt_sha256": digest(receipt_path),
    "files": [{
        "role": "libCaloAna24.so",
        "path": str(library_path.resolve()),
        "sha256": digest(library_path),
    }],
}
manifest_path.write_text(json.dumps(manifest, sort_keys=True) + "\n")
PY
origin_library_sha="$(python3 - "$origin_library" <<'PY'
from pathlib import Path
import hashlib
import sys
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)"
import_plan="$(bash "$builder" --output-dir "$output" --setup-script /usr/bin/true \
  --jobs 2 --ppg-source-runtime-manifest "$origin_manifest")"
grep -q 'ppg_binary_mode: source_locked_runtime_import' <<<"$import_plan" || \
  fail "source-locked binary-import mode missing"
grep -q "ppg_origin_library_sha256: ${origin_library_sha}" <<<"$import_plan" || \
  fail "origin library digest missing from import plan"
import_token="$(sed -n 's/^  token: //p' <<<"$import_plan")"
[[ "$import_token" =~ ^ppg12-new17:[0-9a-f]{64}$ ]] || fail "import token is malformed"
[[ "$import_token" != "$token" ]] || fail "binary import did not change the sealed token"
[[ ! -e "$output" ]] || fail "import plan mode mutated the output path"

printf 'tamper\n' >> "$origin_library"
if bash "$builder" --output-dir "$output" --setup-script /usr/bin/true \
    --ppg-source-runtime-manifest "$origin_manifest" >/dev/null 2>&1; then
  fail "tampered origin library was accepted"
fi
[[ ! -e "$output" ]] || fail "rejected import plan mutated the output path"

if bash "$builder" --build --token ppg12-new17:wrong --output-dir "$output" \
    --setup-script /usr/bin/true --jobs 2 >/dev/null 2>&1; then
  fail "wrong token was accepted"
fi
[[ ! -e "$output" ]] || fail "wrong-token build mutated the output path"

# A correct token reaches the clean-environment setup check, but /usr/bin/true
# cannot establish OFFLINE_MAIN.  The build must fail before creating output.
if bash "$builder" --build --token "$token" --output-dir "$output" \
    --setup-script /usr/bin/true --jobs 2 >/dev/null 2>&1; then
  fail "fake setup unexpectedly completed a build"
fi
[[ ! -e "$output" ]] || fail "failed setup mutated the output path"

mkdir "$output"
if bash "$builder" --output-dir "$output" --setup-script /usr/bin/true >/dev/null 2>&1; then
  fail "existing output path was accepted"
fi
rm -rf "$output"

if bash "$builder" --output-dir "${repo_root}/unsafe_runtime" \
    --setup-script /usr/bin/true >/dev/null 2>&1; then
  fail "unsafe repository destination was accepted"
fi

# Static invariants cover the remote-only portion that cannot be built on the
# macOS test host.  These are contract assertions, not substitutes for the
# builder's own readelf/ldd/ROOT runtime gates.
for invariant in \
  'exec /usr/bin/env -i' \
  'git "${git_safe[@]}" -C "$ppg_repo_real" show' \
  '/bin/bash "${ppg_stage}/autogen.sh"' \
  'built_ppg="$(resolve_installed_lib libCaloAna24.so)"' \
  '"${runtime_root}/lib/libCaloAna24.so"' \
  '"archived_ppg12_binary_reused": False' \
  '"ppg12_source_locked_rebuild": ppg_binary_mode == "source_locked_rebuild"' \
  '"ppg12_source_locked_binary_import": ppg_binary_mode == "source_locked_runtime_import"' \
  'reexec_args+=(--ppg-source-runtime-manifest "$ppg_source_runtime_manifest")' \
  'cmp -s "$built_ppg" "${runtime_root}/lib/libCaloAna24.so"' \
  'ppg_source_runtime_manifest' \
  'ppg_source_build_receipt' \
  'source "$setup_script" -n new.17' \
  'root_libdir="$(root-config --libdir)"' \
  'for link_dir in "${OFFLINE_MAIN}/lib" "${OFFLINE_MAIN}/lib64" "$root_libdir"' \
  'export LDFLAGS="-L${install_root}/lib${release_ldflags}"' \
  'libRecoilJets_la_LDFLAGS = $(AM_LDFLAGS) -L$(ROOTSYS)/lib -no-undefined -version-info 0:0:0' \
  '"sealed_link_directories": sealed_link_dirs.split(":")' \
  'calo_calib="${OFFLINE_MAIN}/rootmacros/Calo_Calib.C"' \
  '/sphenix/u/patsfan753/thesisAnalysis(_auau)?/install' \
  "not found|/release/release_ana/|/sphenix/u/patsfan753/thesisAnalysis(_auau)?/install|/sphenix/user/patsfan753/install" \
  'cp -L "$release_calo"' \
  'cp -L "$release_clusteriso"' \
  'PPG12OraclePhotonClusterBuilder.cc' \
  '"$recoil_stage" "$install_root"' \
  '/bin/bash "${recoil_stage}/autogen.sh"' \
  'replace_word_exact "$runtime_impl" 16' \
  'replace_exact "$runtime_impl"' \
  'readelf -d "$runtime_lib"' \
  'ldd "$runtime_lib"' \
  'root -l -b -q "${smoke_macro}"' \
  '"recoil_macro"' \
  '"recoil_impl"' \
  '"PhotonClusterBuilder.h"'; do
  grep -Fq "$invariant" "$builder" || fail "missing builder invariant: $invariant"
done
for invariant in \
  'expected_recoeff_sha256="e9b25fdb6dd8a6bfbbad029cb90aaddc9489fdf2846c630ea63c8c41ac771eee"' \
  'expected_recoeff_config_sha256="42b7be1628843d5b7607ab988ffb58c6d019d8ade01b3c4d528611498db95732"' \
  'expected_recoeff_period_config_0mrad_sha256="3995033c8867f4b0e21d5ebc025d36395185671da474fceec128b20db7218be2"' \
  'expected_recoeff_period_config_1p5mrad_sha256="6d2e4cc691e2fdd49271486ef193055704da00bcbd6b50ced76fcdd99cd050b8"' \
  'expected_truth_vertex_reweight_0mrad_sha256="4c2a50fa2dd4fe6e3f8b823454f19753367876edabc9fe48164447a0d07be6b9"' \
  'expected_truth_vertex_reweight_1p5mrad_sha256="1429442b2cce368bcc4b4fd613f706f27b9c194f1e835d800238e059c0804d4d"' \
  'expected_apply_bdt_sha256="bd6e7c5bc9858ddad9bc835552d818c00290bb7de3f5036f44d8bf4804734366"' \
  'expected_apply_config_sha256="b8d1bc359a647cc913f213777fc42958b532b30c37a63bb318680130eb6e321b"' \
  'expected_apply_npb_sha256="d6086dadac534013cda15cdfb69c1683776d3456d9e439903589653e8ac19eab"' \
  'apply_model_names=(base base_vr base_v0 base_v1 base_v2 base_v3 base_E base_v0E base_v1E base_v2E base_v3E)' \
  '(f"ppg_apply_model_{name}"' \
  'ppg_apply_npb_model' \
  'ppg_recoeff_period_config_0mrad' \
  'ppg_recoeff_period_config_1p5mrad' \
  'ppg_recoeff_truth_vertex_reweight_0mrad' \
  'ppg_recoeff_truth_vertex_reweight_1p5mrad' \
  'ppg_recoeff_yaml_cpp_header_tree_receipt' \
  'yaml_cpp_include_dir="/sphenix/u/shuhang98/install/include"' \
  'staged yaml-cpp header tree differs from source' \
  '#include <yaml-cpp/yaml.h>' \
  'YAML::Load("ppg12_oracle_smoke: 17")' \
  'yaml_header_smoke["ppg12_oracle_smoke"].as<int>() != 17' \
  '"period_configs": {' \
  'apply_BDT split model hash differs'; do
  grep -Fq "$invariant" "$builder" || fail "missing source/model hash invariant: $invariant"
done
grep -Fq -- '-lg4eval \' "$builder" || \
  fail "staged oracle link does not include release libg4eval"

if grep -Fq 'build_root}/caloreco' "$builder"; then
  fail "builder still rebuilds the full local CaloReco package"
fi
grep -Fq '"full_calo_reco_rebuild": False' "$builder" || \
  fail "receipt does not explicitly forbid a full CaloReco rebuild"
grep -Fq 'PPG12OraclePhotonClusterBuilder.h" 3' "$builder" || \
  fail "oracle header-guard rewrite count is not sealed to all three occurrences"
grep -Fq '"contains_release_calibration_and_cluster_template": True' "$builder" || \
  fail "receipt does not identify release calibration ownership"

printf 'PPG12_ORACLE_NEW17_BUILDER_TEST_PASS\n'
