#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
controller="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the119_replay_foundation_canaries.sh"
wrapper="${repo_root}/scripts/sdcc/workflows/diagnostics/submit_the134_shower_factorial_canaries.sh"

bash -n "$controller"
bash -n "$wrapper"
if grep -Eq '(^|[^[:alnum:]_])(/usr/bin/)?(say|afplay)([^[:alnum:]_]|$)|osascript|NSSound' \
  "$controller" "$wrapper"; then
  printf 'forbidden local audio command in THE-119/THE-134 canary surface\n' >&2
  exit 1
fi

die() {
  printf 'fixture failure: %s\n' "$*" >&2
  exit 2
}
eval "$(sed -n '/^validate_writer_extra_template(){/,/^}/p' "$controller")"
eval "$(sed -n '/^environment_template_value(){/,/^}/p' "$controller")"
eval "$(sed -n '/^render_writer_extra(){/,/^}/p' "$controller")"
eval "$(sed -n '/^render_arm_extra(){/,/^}/p' "$controller")"
eval "$(sed -n '/^sha_file(){/,/^}/p' "$controller")"
eval "$(sed -n '/^validate_pinned_runtime_provider(){/,/^}/p' "$controller")"

extra_common_template='SHARED=1'
extra_pp_template='SYSTEM=pp'
extra_auau_template='SYSTEM=auau'
writer_extra_common_template='A=1;B=__OUTPUT__/sidecar.root'
writer_extra_pp_template='WRITER_SYSTEM=pp'
writer_extra_auau_template='WRITER_SYSTEM=auau'

[[ "$(render_arm_extra pp writer /tmp/pp)" == \
  'SHARED=1;SYSTEM=pp;A=1;B=/tmp/pp/sidecar.root;WRITER_SYSTEM=pp' ]]
[[ "$(render_arm_extra auau writer /tmp/auau)" == \
  'SHARED=1;SYSTEM=auau;A=1;B=/tmp/auau/sidecar.root;WRITER_SYSTEM=auau' ]]
[[ "$(render_arm_extra pp direct /tmp/direct)" == \
  'SHARED=1;SYSTEM=pp' ]]
[[ "$(environment_template_value \
  'SYSTEM=pp;RJ_PPG12_PHOTON_YIELD=1;OTHER=ok' \
  RJ_PPG12_PHOTON_YIELD)" == 1 ]]
if environment_template_value \
  'SYSTEM=auau;OTHER=ok' RJ_PPG12_PHOTON_YIELD >/dev/null; then
  printf 'missing p+p photon-yield materialization flag was fabricated\n' >&2
  exit 1
fi

if ( validate_writer_extra_template invalid 'A=1;A=2' ) >/dev/null 2>&1; then
  printf 'duplicate writer environment keys were accepted\n' >&2
  exit 1
fi
if ( validate_writer_extra_template invalid 'A=__UNKNOWN__' ) >/dev/null 2>&1; then
  printf 'unknown writer environment placeholder was accepted\n' >&2
  exit 1
fi
if ( render_writer_extra pp writer '/tmp/unsafe path' ) >/dev/null 2>&1; then
  printf 'unsafe writer output path was accepted\n' >&2
  exit 1
fi
writer_extra_common_template='SHARED=writer-duplicate'
if ( render_arm_extra pp writer /tmp/pp ) >/dev/null 2>&1; then
  printf 'cross-arm duplicate environment keys were accepted\n' >&2
  exit 1
fi
writer_extra_common_template='A=1;B=__OUTPUT__/sidecar.root'

for required in \
  'RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1' \
  'RJ_THE134_MULTIVIEW_TRAINING_FILE=__OUTPUT__/RJPhotonTrainingViewV1.root' \
  "RJ_THE119_PP_EXTRA_ENV_TEMPLATE='RJ_REPLAY_PERIOD=0mrad" \
  "RJ_THE119_AUAU_EXTRA_ENV_TEMPLATE='RJ_REPLAY_PERIOD=AUAU_RUN24" \
  'for arm in direct writer' \
  '${arm}:pp_inclusive_sim:run28_jet8' \
  '${arm}:auau_inclusive_embedded:run28_embeddedJet12' \
  'RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS=1' \
  'RJ_PINNED_CALO_RECO_SONAME="$soname"' \
  'RJ_PINNED_CALO_RECO_SHA256="$calo_reco_sha"' \
  'RJ_PINNED_PHOTON_CLUSTER_BUILDER_HEADER_SHA256="$builder_header_sha"' \
  'RJ_PINNED_CALO_RECO_BUILD_RECEIPT_SHA256="$build_receipt_sha"' \
  'RJ_PINNED_CALO_RECO_SOURCE_MANIFEST_SHA256="$source_manifest_sha"' \
  'replacement_expected_calo_reco_sha=b32e89b3b43efa57dc825f7e0b81f126b8886fe5b83c2f9337524432539b755b' \
  'replacement_expected_builder_header_sha=255fb1b4b9a0fdb9b0ee4709483ac30a04ee8dd2813f99e6afc3914e5cd1e20f' \
  'replacement_expected_calo_reco_build_receipt_sha=ff397bfe281105454ac70b9f308efa960c81a29f3cba6903cf1aa650d4c41f0c' \
  'replacement_expected_calo_reco_source_manifest_sha=9dc58e9c0a6dc5ccc42d0baf86cb04ed17b4a0b6e5ba390b7745b1332cf8ade4' \
  'RJ_AUAU_LIBRARY_OVERRIDE="$auau_library"' \
  'readonly replacement_release_name=ana.560' \
  'readonly replacement_offline_main=/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560' \
  '8810cdfcdb1302a06567d3cf0744c12f8a9b53ae28355e8cc26b0e81d0621685' \
  '807a50cb16ba85d9d0232c06bf2ab9b369b38f3c3626554de9827cdd80225bd2' \
  'd992f11a1e6a1ccb1d74c6e09da79114c9e9742fd2f9cf5a3bc284c8aec9d507' \
  'STANDALONE_PP_REPLACEMENT_PREFLIGHT_CONSUMES_NAMESPACE' \
  '--terminal-gate-receipt "$RJ_THE119_EVIDENCE_ROOT/terminal_gate_receipt.json"'; do
  grep -F -- "$required" "$wrapper" >/dev/null
done
for required in \
  'submit(){' \
  '  preflight' \
  '  assert_fresh' \
  'environment_template_value "$final_extra" RJ_PPG12_PHOTON_YIELD' \
  'env RJ_PPG12_PHOTON_YIELD="$ppg12_photon_yield"' \
  'exec scripts/sdcc/workflows/diagnostics/submit_the119_replay_foundation_canaries.sh "$mode"'; do
  grep -F -- "$required" "$controller" "$wrapper" >/dev/null
done

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/the134-sidecar-selector.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT
eval "$(sed -n '/^submit_pp(){/,/^}/p' "$controller")"
eval "$(sed -n '/^require_replacement_namespaces()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^build_sidecar_only_keys()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^require_exact_sidecar_keys()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^require_sha256()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^sha256_file()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^require_file_hash()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^require_env_unset_or_exact()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^validate_pp_replacement_calo_reco_authority()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^configure_pp_replacement_runtime_provider()/,/^}/p' "$wrapper")"
eval "$(sed -n '/^validate_terminal_rows()/,/^}/p' "$wrapper")"

submit_fixture="${tmpdir}/submit-fixture"
mkdir -p "$submit_fixture"
printf '%s\n' \
  '#!/usr/bin/env bash' \
  'set -euo pipefail' \
  'case "${TEST_EXPECTED_PHOTON_YIELD}" in' \
  '  1)' \
  '    [[ "${RJ_PPG12_PHOTON_YIELD-}" == 1 ]] || exit 1' \
  '    [[ "${RJ_SUBMIT_EXTRA_ENV}" == *";RJ_PPG12_PHOTON_YIELD=1;"* ]] || exit 1' \
  '    ;;' \
  '  absent)' \
  '    [[ "${RJ_PPG12_PHOTON_YIELD-}" == 0 ]] || exit 1' \
  '    [[ "${RJ_SUBMIT_EXTRA_ENV}" != *"RJ_PPG12_PHOTON_YIELD="* ]] || exit 1' \
  '    ;;' \
  '  *) exit 2 ;;' \
  'esac' \
  > "${submit_fixture}/RecoilJets_Condor_submit.sh"
chmod +x "${submit_fixture}/RecoilJets_Condor_submit.sh"
want_row() { return 0; }
pp_extra() { printf '%s' "$fixture_extra"; }
base="${tmpdir}/output"
pp_cfg="${tmpdir}/analysis_config.yaml"
pp_lib="${tmpdir}/libRecoilJets.so"
schema_sha=0000000000000000000000000000000000000000000000000000000000000000
canary_nevents=1
(
  cd "$submit_fixture"
  fixture_extra='SYSTEM=pp;RJ_PPG12_PHOTON_YIELD=1;OTHER=ok'
  export TEST_EXPECTED_PHOTON_YIELD=1
  submit_pp pp_inclusive_sim isSimInclusive run28_jet8 direct
)
(
  cd "$submit_fixture"
  export RJ_PPG12_PHOTON_YIELD=1
  fixture_extra='SYSTEM=pp;OTHER=ok'
  export TEST_EXPECTED_PHOTON_YIELD=absent
  submit_pp pp_inclusive_sim isSimInclusive run28_jet8 direct
)
(
  cd "$submit_fixture"
  fixture_extra='SYSTEM=pp;RJ_PPG12_PHOTON_YIELD=1;OTHER=ok'
  export TEST_EXPECTED_PHOTON_YIELD=1
  submit_pp pp_inclusive_sim isSimInclusive run28_jet8 writer
)

provider_fixture="${tmpdir}/provider.bin"
printf 'provider-fixture\n' > "$provider_fixture"
provider_fixture_sha="$(sha256_file "$provider_fixture")"
require_sha256 fixture "$provider_fixture_sha"
require_file_hash fixture "$provider_fixture" "$provider_fixture_sha"
if ( require_file_hash fixture "$provider_fixture" \
  0000000000000000000000000000000000000000000000000000000000000000 ) \
  >/dev/null 2>&1; then
  printf 'pinned provider hash drift was accepted\n' >&2
  exit 1
fi
unset RJ_AUAU_LIBRARY_OVERRIDE

provider_root="${tmpdir}/release/release_ana/ana.560"
mkdir -p "${provider_root}/lib" "${provider_root}/lib64"
pp_lib="${tmpdir}/libRecoilJets.so"
auau_lib="${tmpdir}/libRecoilJetsAuAu.so"
calo_authority="${tmpdir}/calo-authority"
mkdir -p "${calo_authority}/install/lib"
pinned_calo_reco_library="${calo_authority}/install/lib/libcalo_reco.so.0.0.0"
original_proof_library="${tmpdir}/original-r5/install/lib/libcalo_reco.so.0.0.0"
pinned_builder_header="${tmpdir}/PhotonClusterBuilder.h"
pinned_calo_io="${provider_root}/lib/libcalo_io.so"
pinned_clusteriso="${provider_root}/lib/libclusteriso.so"
pinned_jetbase="${provider_root}/lib/libjetbase.so"
for provider in \
  "$pp_lib" "$auau_lib" "$pinned_calo_reco_library" \
  "$pinned_builder_header" "$pinned_calo_io" "$pinned_clusteriso" \
  "$pinned_jetbase"; do
  printf 'provider=%s\n' "$(basename "$provider")" > "$provider"
done
ln -s libcalo_reco.so.0.0.0 "${calo_authority}/install/lib/libcalo_reco.so"
ln -s libcalo_reco.so.0.0.0 "${calo_authority}/install/lib/libcalo_reco.so.0"
pp_library="$pp_lib"
auau_library="$auau_lib"
pinned_builder_header_sha="$(sha256_file "$pinned_builder_header")"
pinned_calo_reco_sha="$(sha256_file "$pinned_calo_reco_library")"
cat > "${calo_authority}/source_manifest.json" <<JSON
{
  "schema": "THE134_ANA560_CALORECO_SOURCE_MANIFEST_V2",
  "coresoftware": {
    "commit": "cba274033b5560e32600cdeaa7676b6ab4a6c971"
  },
  "mapping_patch": {
    "changed_scientific_controls": []
  },
  "overlay": {
    "PhotonClusterBuilder.h": {
      "staged_sha256": "${pinned_builder_header_sha}"
    }
  }
}
JSON
pinned_calo_reco_source_manifest="${calo_authority}/source_manifest.json"
pinned_calo_reco_source_manifest_sha="$(sha256_file "$pinned_calo_reco_source_manifest")"
cat > "${calo_authority}/build_receipt.json" <<JSON
{
  "schema": "THE134_ANA560_CALORECO_BUILD_RECEIPT_V3",
  "status": "PASS",
  "runtime": {
    "release": "ana.560",
    "offline_main": "${provider_root}",
    "ldd_not_found": false,
    "mutable_user_dependency": false,
    "root_load": {
      "status": "PASS",
      "load_return_code": 0,
      "preload_provider_count": 0,
      "provider_count": 1,
      "provider_realpath": "${original_proof_library}"
    }
  },
  "build": {
    "coresoftware_commit": "cba274033b5560e32600cdeaa7676b6ab4a6c971"
  },
  "artifact": {
    "sha256": "${pinned_calo_reco_sha}",
    "library": "install/lib/libcalo_reco.so.0.0.0",
    "symlinks": {
      "install/lib/libcalo_reco.so": "libcalo_reco.so.0.0.0",
      "install/lib/libcalo_reco.so.0": "libcalo_reco.so.0.0.0"
    }
  },
  "source_manifest": {
    "sha256": "${pinned_calo_reco_source_manifest_sha}",
    "path": "source_manifest.json"
  },
  "abi": {
    "status": "PASS",
    "soname": "libcalo_reco.so.0",
    "soname_expected": "libcalo_reco.so.0",
    "needed_exact_match": true,
    "rpath_runpath_exact_match": true,
    "removed_symbols": {
      "count": 0
    }
  },
  "single_provider": {
    "photon_cluster_builder_process_event_definitions": 1,
    "raw_cluster_builder_topo_process_event_definitions": 1,
    "forbidden_standalone_provider_count": 0,
    "other_installed_shared_objects": [],
    "provider_probe": {
      "status": "PASS",
      "preload_provider_count": 0,
      "provider_count": 1,
      "provider_realpath": "${original_proof_library}"
    }
  },
  "mapping_patch": {
    "id": "THE134_RAWCLUSTERBUILDERTOPO_DETECTOR_EXPLICIT_CHANNEL_MAP_V1",
    "scientific_controls_changed": []
  }
}
JSON
pinned_calo_reco_build_receipt="${calo_authority}/build_receipt.json"
pinned_calo_reco_build_receipt_sha="$(sha256_file "$pinned_calo_reco_build_receipt")"
unset \
  RJ_PP_LIBRARY_OVERRIDE \
  RJ_AUAU_LIBRARY_OVERRIDE \
  RJ_CALO_RECO_LIBRARY_OVERRIDE \
  RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE \
  RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE \
  RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS \
  RJ_PINNED_CALO_RECO_SONAME \
  RJ_PINNED_CALO_RECO_SHA256 \
  RJ_PINNED_PHOTON_CLUSTER_BUILDER_HEADER_SHA256 \
  RJ_PINNED_RELEASE_NAME \
  RJ_PINNED_OFFLINE_MAIN \
  RJ_PINNED_RELEASE_CALO_IO_PATH \
  RJ_PINNED_RELEASE_CALO_IO_SHA256 \
  RJ_PINNED_RELEASE_CLUSTERISO_PATH \
  RJ_PINNED_RELEASE_CLUSTERISO_SHA256 \
  RJ_PINNED_RELEASE_JETBASE_PATH \
  RJ_PINNED_RELEASE_JETBASE_SHA256 \
  RJ_RELEASE_CORE_LIB_DIR \
  RJ_RELEASE_CORE_LIB64_DIR \
  RJ_FORCE_RELEASE_CORE_LIBS \
  RJ_FORCE_RELEASE_CALO_IO \
  RJ_PINNED_CALO_RECO_BUILD_RECEIPT \
  RJ_PINNED_CALO_RECO_BUILD_RECEIPT_SHA256 \
  RJ_PINNED_CALO_RECO_SOURCE_MANIFEST \
  RJ_PINNED_CALO_RECO_SOURCE_MANIFEST_SHA256
configure_pp_replacement_runtime_provider \
  "$pinned_calo_reco_library" "$pinned_calo_reco_sha" \
  "$pinned_builder_header" "$pinned_builder_header_sha" \
  ana.560 "$provider_root" "${provider_root}/lib" "${provider_root}/lib64" \
  "$pinned_calo_io" "$(sha256_file "$pinned_calo_io")" \
  "$pinned_clusteriso" "$(sha256_file "$pinned_clusteriso")" \
  "$pinned_jetbase" "$(sha256_file "$pinned_jetbase")" \
  libcalo_reco.so.0 \
  "$pinned_calo_reco_build_receipt" "$pinned_calo_reco_build_receipt_sha" \
  "$pinned_calo_reco_source_manifest" "$pinned_calo_reco_source_manifest_sha"
child_provider_environment="$(env)"
for inherited in \
  "RJ_PP_LIBRARY_OVERRIDE=${pp_lib}" \
  "RJ_AUAU_LIBRARY_OVERRIDE=${auau_lib}" \
  "RJ_CALO_RECO_LIBRARY_OVERRIDE=${pinned_calo_reco_library}" \
  "RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE=${pinned_builder_header}" \
  'RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE=' \
  'RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS=1' \
  'RJ_PINNED_CALO_RECO_SONAME=libcalo_reco.so.0' \
  'RJ_PINNED_RELEASE_NAME=ana.560' \
  "RJ_PINNED_OFFLINE_MAIN=${provider_root}" \
  "RJ_PINNED_RELEASE_CALO_IO_PATH=${pinned_calo_io}" \
  "RJ_PINNED_RELEASE_CLUSTERISO_PATH=${pinned_clusteriso}" \
  "RJ_PINNED_RELEASE_JETBASE_PATH=${pinned_jetbase}" \
  "RJ_RELEASE_CORE_LIB_DIR=${provider_root}/lib" \
  "RJ_RELEASE_CORE_LIB64_DIR=${provider_root}/lib64" \
  'RJ_FORCE_RELEASE_CORE_LIBS=0' \
  'RJ_FORCE_RELEASE_CALO_IO=0' \
  "RJ_PINNED_CALO_RECO_BUILD_RECEIPT=${pinned_calo_reco_build_receipt}" \
  "RJ_PINNED_CALO_RECO_BUILD_RECEIPT_SHA256=${pinned_calo_reco_build_receipt_sha}" \
  "RJ_PINNED_CALO_RECO_SOURCE_MANIFEST=${pinned_calo_reco_source_manifest}" \
  "RJ_PINNED_CALO_RECO_SOURCE_MANIFEST_SHA256=${pinned_calo_reco_source_manifest_sha}"; do
  grep -Fx -- "$inherited" <<< "$child_provider_environment" >/dev/null ||
    { printf 'configured provider field was not inherited: %s\n' "$inherited" >&2; exit 1; }
done
pinned_calo_reco_mode=1
pinned_calo_reco_sha="$RJ_PINNED_CALO_RECO_SHA256"
pinned_builder_header_sha="$RJ_PINNED_PHOTON_CLUSTER_BUILDER_HEADER_SHA256"
pinned_calo_io_sha="$RJ_PINNED_RELEASE_CALO_IO_SHA256"
pinned_clusteriso_sha="$RJ_PINNED_RELEASE_CLUSTERISO_SHA256"
pinned_jetbase_sha="$RJ_PINNED_RELEASE_JETBASE_SHA256"
pinned_builder_library="$RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE"
pinned_calo_reco_soname="$RJ_PINNED_CALO_RECO_SONAME"
pinned_release_name="$RJ_PINNED_RELEASE_NAME"
pinned_offline_main="$RJ_PINNED_OFFLINE_MAIN"
pinned_release_lib="$RJ_RELEASE_CORE_LIB_DIR"
pinned_release_lib64="$RJ_RELEASE_CORE_LIB64_DIR"
pinned_calo_reco_build_receipt="$RJ_PINNED_CALO_RECO_BUILD_RECEIPT"
pinned_calo_reco_build_receipt_sha="$RJ_PINNED_CALO_RECO_BUILD_RECEIPT_SHA256"
pinned_calo_reco_source_manifest="$RJ_PINNED_CALO_RECO_SOURCE_MANIFEST"
pinned_calo_reco_source_manifest_sha="$RJ_PINNED_CALO_RECO_SOURCE_MANIFEST_SHA256"
validate_pinned_runtime_provider

valid_calo_sha="$pinned_calo_reco_sha"
pinned_calo_reco_sha=0000000000000000000000000000000000000000000000000000000000000000
if ( validate_pinned_runtime_provider ) >/dev/null 2>&1; then
  printf 'wrong pinned CaloReco hash was accepted\n' >&2
  exit 1
fi
pinned_calo_reco_sha="$valid_calo_sha"
RJ_AUAU_LIBRARY_OVERRIDE="${tmpdir}/mutable-default-auau.so"
if ( validate_pinned_runtime_provider ) >/dev/null 2>&1; then
  printf 'wrong AuAu snapshot companion was accepted by the runtime preflight\n' >&2
  exit 1
fi
RJ_AUAU_LIBRARY_OVERRIDE="$auau_lib"
wrong_clusteriso="$pinned_clusteriso"
pinned_clusteriso="${provider_root}/lib64/libclusteriso.so"
cp "$wrong_clusteriso" "$pinned_clusteriso"
pinned_clusteriso_sha="$(sha256_file "$pinned_clusteriso")"
if ( validate_pinned_runtime_provider ) >/dev/null 2>&1; then
  printf 'wrong release companion path was accepted\n' >&2
  exit 1
fi
pinned_clusteriso="$wrong_clusteriso"
pinned_clusteriso_sha="$(sha256_file "$pinned_clusteriso")"
valid_receipt_sha="$pinned_calo_reco_build_receipt_sha"
pinned_calo_reco_build_receipt_sha=0000000000000000000000000000000000000000000000000000000000000000
if ( validate_pinned_runtime_provider ) >/dev/null 2>&1; then
  printf 'wrong CaloReco build-receipt hash was accepted\n' >&2
  exit 1
fi
pinned_calo_reco_build_receipt_sha="$valid_receipt_sha"
valid_source_sha="$pinned_calo_reco_source_manifest_sha"
pinned_calo_reco_source_manifest_sha=0000000000000000000000000000000000000000000000000000000000000000
if ( validate_pinned_runtime_provider ) >/dev/null 2>&1; then
  printf 'wrong CaloReco source-manifest hash was accepted\n' >&2
  exit 1
fi
pinned_calo_reco_source_manifest_sha="$valid_source_sha"
ln -sfn libcalo_reco.so.0 "${calo_authority}/install/lib/libcalo_reco.so"
if ( validate_pp_replacement_calo_reco_authority \
  "$pinned_calo_reco_build_receipt" "$pinned_calo_reco_build_receipt_sha" \
  "$pinned_calo_reco_source_manifest" "$pinned_calo_reco_source_manifest_sha" \
  "$pinned_calo_reco_library" "$pinned_calo_reco_sha" \
  "$pinned_builder_header" "$pinned_builder_header_sha" \
  ana.560 "$provider_root" cba274033b5560e32600cdeaa7676b6ab4a6c971 ) \
  >/dev/null 2>&1; then
  printf 'receipt-inconsistent CaloReco loader alias was accepted\n' >&2
  exit 1
fi
ln -sfn libcalo_reco.so.0.0.0 "${calo_authority}/install/lib/libcalo_reco.so"
cp "$pinned_calo_reco_build_receipt" "${calo_authority}/mutable_receipt.json"
python3 - "${calo_authority}/mutable_receipt.json" <<'PY'
import json
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
payload = json.loads(path.read_text(encoding="utf-8"))
payload["runtime"]["mutable_user_dependency"] = True
path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
if ( validate_pp_replacement_calo_reco_authority \
  "${calo_authority}/mutable_receipt.json" \
  "$(sha256_file "${calo_authority}/mutable_receipt.json")" \
  "$pinned_calo_reco_source_manifest" "$pinned_calo_reco_source_manifest_sha" \
  "$pinned_calo_reco_library" "$pinned_calo_reco_sha" \
  "$pinned_builder_header" "$pinned_builder_header_sha" \
  ana.560 "$provider_root" cba274033b5560e32600cdeaa7676b6ab4a6c971 ) \
  >/dev/null 2>&1; then
  printf 'mutable CaloReco runtime authority was accepted\n' >&2
  exit 1
fi
unset RJ_AUAU_LIBRARY_OVERRIDE
require_env_unset_or_exact RJ_AUAU_LIBRARY_OVERRIDE /frozen/auau.so
RJ_AUAU_LIBRARY_OVERRIDE=/mutable/default/auau.so
if ( require_env_unset_or_exact RJ_AUAU_LIBRARY_OVERRIDE /frozen/auau.so ) \
  >/dev/null 2>&1; then
  printf 'mutable AuAu snapshot companion was accepted\n' >&2
  exit 1
fi
unset RJ_AUAU_LIBRARY_OVERRIDE

replacement_keys="$(build_sidecar_only_keys 0 | paste -sd, -)"
[[ "$replacement_keys" == \
  'direct:pp_inclusive_sim:run28_jet8,writer:pp_inclusive_sim:run28_jet8' ]]
full_keys="$(build_sidecar_only_keys 1 | paste -sd, -)"
[[ "$full_keys" == \
  'direct:pp_inclusive_sim:run28_jet8,direct:auau_inclusive_embedded:run28_embeddedJet12,writer:pp_inclusive_sim:run28_jet8,writer:auau_inclusive_embedded:run28_embeddedJet12' ]]
require_exact_sidecar_keys "$replacement_keys" "$replacement_keys" replacement
if ( require_exact_sidecar_keys \
  'direct:pp_inclusive_sim:run28_jet8' \
  "$replacement_keys" replacement ) >/dev/null 2>&1; then
  printf 'partial p+p replacement selector was accepted\n' >&2
  exit 1
fi
if ( require_exact_sidecar_keys \
  "${replacement_keys},writer:auau_inclusive_embedded:run28_embeddedJet12" \
  "$replacement_keys" replacement ) >/dev/null 2>&1; then
  printf 'widened p+p replacement selector was accepted\n' >&2
  exit 1
fi

replacement_tag=the134_pp_replacement_fixture
RJ_THE134_TAG="$replacement_tag"
RJ_THE134_OUTPUT_ROOT="${tmpdir}/output/${replacement_tag}"
RJ_THE134_EVIDENCE_ROOT="${tmpdir}/evidence/${replacement_tag}"
require_replacement_namespaces preflight
for special_tag in . ..; do
  RJ_THE134_TAG="$special_tag"
  RJ_THE134_OUTPUT_ROOT="${tmpdir}/special/output/${special_tag}"
  RJ_THE134_EVIDENCE_ROOT="${tmpdir}/special/evidence/${special_tag}"
  if ( require_replacement_namespaces preflight ) >/dev/null 2>&1; then
    printf 'special replacement tag %s was accepted\n' "$special_tag" >&2
    exit 1
  fi
done
RJ_THE134_TAG="$replacement_tag"
RJ_THE134_OUTPUT_ROOT="${tmpdir}/canonical/${replacement_tag}"
RJ_THE134_EVIDENCE_ROOT="${tmpdir}/canonical/sub/../${replacement_tag}"
if ( require_replacement_namespaces preflight ) >/dev/null 2>&1; then
  printf 'canonically identical replacement namespaces were accepted\n' >&2
  exit 1
fi
RJ_THE134_OUTPUT_ROOT="${tmpdir}/output/${replacement_tag}"
RJ_THE134_EVIDENCE_ROOT="${tmpdir}/evidence/${replacement_tag}"
mkdir -p "$RJ_THE134_EVIDENCE_ROOT"
if ( require_replacement_namespaces preflight ) >/dev/null 2>&1; then
  printf 'pre-existing replacement evidence namespace was accepted\n' >&2
  exit 1
fi
require_replacement_namespaces status
require_replacement_namespaces validate
require_replacement_namespaces aggregate
rm -rf "$RJ_THE134_EVIDENCE_ROOT"
mkdir -p "$RJ_THE134_OUTPUT_ROOT"
if ( require_replacement_namespaces submit ) >/dev/null 2>&1; then
  printf 'pre-existing replacement output namespace was accepted\n' >&2
  exit 1
fi
rm -rf "$RJ_THE134_OUTPUT_ROOT"
RJ_THE134_OUTPUT_ROOT="${tmpdir}/output/not-the-tag"
if ( require_replacement_namespaces preflight ) >/dev/null 2>&1; then
  printf 'mismatched replacement tag/output namespace was accepted\n' >&2
  exit 1
fi
if (
  unset RJ_THE134_EVIDENCE_ROOT
  require_replacement_namespaces preflight
) >/dev/null 2>&1; then
  printf 'replacement preflight without explicit evidence namespace was accepted\n' >&2
  exit 1
fi

RJ_THE119_EVIDENCE_ROOT="${tmpdir}/terminal"
RJ_THE119_OUTPUT_ROOT="${tmpdir}/terminal-output"
RJ_THE119_TAG=the134_pp_replacement_terminal_fixture
mkdir -p "$RJ_THE119_EVIDENCE_ROOT"
mkdir -p "$RJ_THE119_OUTPUT_ROOT"
printf '%s\n' 'fixture-preflight' > \
  "${RJ_THE119_EVIDENCE_ROOT}/preflight_receipt.txt"
sidecar_pp_replacement_mode=1
condor_q() { return 0; }
condor_history() { printf '%s\n' '4 0'; }
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
validate_terminal_rows
terminal_receipt="${RJ_THE119_EVIDENCE_ROOT}/terminal_gate_receipt.json"
[[ -s "$terminal_receipt" ]]
python3 - "$terminal_receipt" <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
assert payload["schema"] == "THE134_PP_REPLACEMENT_TERMINAL_GATE_V1"
assert payload["row_count"] == 2
assert {row["role"] for row in payload["rows"]} == {"direct", "writer"}
assert len({row["cluster_proc"] for row in payload["rows"]}) == 2
assert all(row["job_status"] == 4 for row in payload["rows"])
assert all(row["exit_code"] == 0 for row in payload["rows"])
PY
printf '%s\n' '1001 0 2 direct-row' > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'one-row replacement terminal receipt was accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8" \
  '1003 0 2 duplicate-row' > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'three-row replacement terminal receipt was accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'duplicate Condor identity was accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'duplicate replacement role was accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8_evil" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8_evil" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'suffixed replacement output paths were accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 --note=${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8 unrelated-destination" \
  "1002 0 2 --note=${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8 unrelated-destination" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'embedded replacement output paths were accepted\n' >&2
  exit 1
fi
printf '%s\n' \
  "1001 0 2 ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8 ${RJ_THE119_OUTPUT_ROOT}/direct/pp_inclusive_sim/run28_jet8" \
  "1002 0 2 --output ${RJ_THE119_OUTPUT_ROOT}/writer/pp_inclusive_sim/run28_jet8" \
  > "${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
if ( validate_terminal_rows ) >/dev/null 2>&1; then
  printf 'one replacement row containing both role paths was accepted\n' >&2
  exit 1
fi

printf 'THE119_WRITER_EXTRA_ENVIRONMENT_TEST_PASS\n'
