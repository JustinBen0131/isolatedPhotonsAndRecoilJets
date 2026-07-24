#!/usr/bin/env bash
set -euo pipefail

# Source-complete, writer-only THE-134 multi-view extraction smoke for one
# explicitly selected p+p period/SI component plus the frozen Au+Au sources.
#
# This is a narrow controller over the established RecoilJets Condor executor.
# Ordinary smoke mode owns exactly one input tuple and one Condor proc for each
# frozen training source.  Capacity mode owns exactly one seven-tuple proc for
# the fixed p+p Jet8 and Au+Au embedded-Jet12 capacity witnesses.  It never
# submits data, p+p Jet40, a direct arm, a merge, or a model training job.  The
# caller must select exactly 0mrad or 1p5mrad through RJ_THE134_PP_PERIOD; the
# other period requires a distinct tag/execution, and DI remains the separately
# typed archived-source path.  `inventory` is read-only and emits the exact
# frozen-source manifest to stdout; submission remains an explicit action.

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd -P)"
cd "$repo_root"

requested_mode="${1:-preflight}"
capacity_mode=0
mode="$requested_mode"
case "$requested_mode" in
  capacity-preflight|capacity-submit|capacity-resume-submit|capacity-status|capacity-validate)
    capacity_mode=1
    mode="${requested_mode#capacity-}"
    ;;
esac
if (( capacity_mode )); then
  tag="${RJ_THE134_EXTRACTION_TAG:-the134_h70_partition_capacity_canary_20260724_v1}"
else
  tag="${RJ_THE134_EXTRACTION_TAG:-the134_h70_multiview_extraction_smoke_20260722_v1}"
fi
output_root="${RJ_THE134_EXTRACTION_OUTPUT_ROOT:-/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/smoke/replay_foundation/${tag}}"
evidence_root="${RJ_THE134_EXTRACTION_EVIDENCE_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/evidence/qa/${tag}}"
submit_root="${RJ_THE134_EXTRACTION_SUBMIT_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/condor_sub/${tag}}"
sim_root="${RJ_THE134_SIM_LIST_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/simListFiles}"

submitter="${RJ_THE134_EXTRACTION_SUBMITTER:-${repo_root}/RecoilJets_Condor_submit.sh}"
pp_executor="${RJ_THE134_PP_EXECUTOR:-${repo_root}/RecoilJets_Condor.sh}"
auau_executor="${RJ_THE134_AUAU_EXECUTOR:-${repo_root}/RecoilJets_Condor_AuAu.sh}"
photon_cluster_header="${RJ_THE134_PHOTON_CLUSTER_HEADER:-${repo_root}/coresoftware_local/offline/packages/CaloBase/PhotonClusterv1.h}"
photon_cluster_builder_header="${RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER:-}"
calo_reco_library="${RJ_THE134_CALO_RECO_LIBRARY:-}"
calo_reco_build_receipt="${RJ_THE134_CALO_RECO_BUILD_RECEIPT:-}"
calo_reco_source_manifest="${RJ_THE134_CALO_RECO_SOURCE_MANIFEST:-}"
release_core_lib_dir="${RJ_THE134_RELEASE_CORE_LIB_DIR:-}"
release_core_lib64_dir="${RJ_THE134_RELEASE_CORE_LIB64_DIR:-}"
release_calo_io="${RJ_THE134_RELEASE_CALO_IO:-}"
release_clusteriso="${RJ_THE134_RELEASE_CLUSTERISO:-}"
release_jetbase="${RJ_THE134_RELEASE_JETBASE:-}"
pp_config="${RJ_THE134_PP_CONFIG:-${repo_root}/macros/analysis_config_the119_pp_replay_foundation.yaml}"
auau_config="${RJ_THE134_AUAU_CONFIG:-${repo_root}/macros/analysis_config_the112_auau_combined_bdt_triplet.yaml}"
pp_library="${RJ_THE134_PP_LIBRARY:-}"
auau_library="${RJ_THE134_AUAU_LIBRARY:-}"
pp_model="${RJ_THE134_PP_MODEL:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the116_models/the116_pp_matched_basev3e_15to35_20260720/models/bdt_ppg12_basev3e_15to35/pp_tight_bdt_ppg12_base_v3E_bdt_noIso_tmva.root}"
auau_model="${RJ_THE134_AUAU_MODEL:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the111_models/the111_combined_corrected_shower_ppg12_labels_20260719_1618/combined/auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.root}"
source_hash_manifest="${RJ_THE134_SOURCE_HASH_MANIFEST:-}"
pp_period="${RJ_THE134_PP_PERIOD:-}"
capacity_canary_id="${RJ_THE134_CAPACITY_CANARY_ID:-}"
capacity_preflight_receipt="${RJ_THE134_CAPACITY_PREFLIGHT_RECEIPT:-}"
capacity_preflight_receipt_sha="${RJ_THE134_CAPACITY_PREFLIGHT_RECEIPT_SHA256:-}"
capacity_full_plan="${RJ_THE134_CAPACITY_FULL_PLAN:-}"
capacity_full_plan_sha="${RJ_THE134_CAPACITY_FULL_PLAN_SHA256:-}"

resolved_config_root="${evidence_root}/resolved_configs"
observed_source_hashes="${evidence_root}/observed_source_hashes.tsv"
submission_manifest="${evidence_root}/submission_manifest.tsv"
duplicate_fingerprint="${evidence_root}/duplicate_fingerprint.sha256"
submission_journal="${evidence_root}/submission_journal.tsv"
submission_receipt="${evidence_root}/submission_receipt.tsv"
source_provenance_json="${evidence_root}/source_provenance.json"
pp_source_provenance_json="${evidence_root}/pp_source_provenance.json"
auau_source_provenance_json="${evidence_root}/auau_source_provenance.json"
root_health_join_certificate="${evidence_root}/root_health_identity_join_certificate.json"
capacity_resource_certificate="${evidence_root}/capacity_resource_certificate.json"
runtime_authority_manifest="${evidence_root}/runtime_authority_manifest.json"
runtime_authority_fingerprint="${runtime_authority_manifest}.sha256"
validator="${RJ_THE134_MATRIX_PREPARER:-${repo_root}/scripts/ml/training/prepare_the134_h70_matrix.py}"

readonly pinned_release_name="ana.560"
readonly pinned_offline_main="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560"
readonly pinned_coresoftware_commit="cba274033b5560e32600cdeaa7676b6ab4a6c971"
readonly pinned_calo_reco_soname="libcalo_reco.so.0"
readonly nominal_et_min="15.0"
readonly nominal_et_max="35.0"
readonly capture_et_min="5.0"
readonly legacy_tree_max_entries="0"
readonly event_limit_per_job="0"
readonly full_training_authority="0"
readonly extraction_cone_r="0.40"
readonly ordinary_group_size="1"
readonly capacity_group_size="7"
readonly capacity_pp_row="pp_background_jet8"
readonly capacity_auau_row="auau_background_jet12"
if (( capacity_mode )); then
  readonly execution_group_size="$capacity_group_size"
  readonly execution_row_count="2"
else
  readonly execution_group_size="$ordinary_group_size"
  readonly execution_row_count="13"
fi
readonly training_schema_text='RJ_PHOTON_TRAINING_VIEW_V1|tree=RJPhotonTrainingViewV1|identity=source,event,candidate,definition|source_stable_inputs=lane,dataset,sample,period,run,segment,input_uri,input_file,manifest|event_stable_inputs=lane,sample,run,segment,event_sequence|candidate_stable_inputs=event,encounter_ordinal,cluster_map_key|features=ordered+contract|labels=truth+source+npb|weights=component-ledger|domain=15to35'

say() { printf '[THE134-EXTRACT] %s\n' "$*"; }
die() { printf '[THE134-EXTRACT][ERROR] %s\n' "$*" >&2; exit 2; }
sha256_cmd() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$@"
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$@"
  else
    die "sha256sum or shasum -a 256 is required"
  fi
}
sha_file() { sha256_cmd "$1" | awk '{print $1}'; }
sha_text() { printf '%s' "$1" | sha256_cmd | awk '{print $1}'; }

require_sha() {
  local name="$1" value="$2"
  [[ "$value" =~ ^[0-9a-f]{64}$ ]] || die "${name} must be a lowercase 64-character SHA-256"
}

require_file_hash() {
  local label="$1" path="$2" expected="$3" actual
  require_sha "$label" "$expected"
  [[ -s "$path" ]] || die "missing ${label} input: ${path}"
  actual="$(sha_file "$path")"
  [[ "$actual" == "$expected" ]] || die "${label} hash drift: expected=${expected} actual=${actual} path=${path}"
}

validate_pp_sim_weight_contract() {
  local period="${1:-}"
  case "$period" in
    0mrad|1p5mrad) ;;
    *)
      die "RJ_THE134_PP_PERIOD must be exactly 0mrad or 1p5mrad; one tagged smoke cannot mix period-specific p+p SIM weights"
      return 2
      ;;
  esac
}

resolve_release_companion() {
  local label="$1" name="$2" path="$3" expected="$4" selected=""
  [[ -n "$path" ]] || die "RJ_THE134_${label} is required"
  if [[ -r "${release_core_lib64_dir}/${name}" ]]; then
    selected="$(cd "$(dirname "${release_core_lib64_dir}/${name}")" && pwd -P)/$(basename "${release_core_lib64_dir}/${name}")"
  elif [[ -r "${release_core_lib_dir}/${name}" ]]; then
    selected="$(cd "$(dirname "${release_core_lib_dir}/${name}")" && pwd -P)/$(basename "${release_core_lib_dir}/${name}")"
  else
    die "pinned ana.560 release companion is missing: ${name}"
  fi
  path="$(cd "$(dirname "$path")" && pwd -P)/$(basename "$path")"
  [[ "$path" == "$selected" ]] ||
    die "${label} path does not match loader-order selection: declared=${path} selected=${selected}"
  require_file_hash "$label" "$path" "$expected"
  printf '%s\n' "$selected"
}

validate_calo_reco_build_authority() {
  python3 - \
    "$calo_reco_build_receipt" "$RJ_THE134_CALO_RECO_BUILD_RECEIPT_SHA256" \
    "$calo_reco_source_manifest" "$RJ_THE134_CALO_RECO_SOURCE_MANIFEST_SHA256" \
    "$calo_reco_library" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256" \
    "$photon_cluster_builder_header" "$RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER_SHA256" \
    "$pinned_offline_main" "$pinned_release_name" "$pinned_coresoftware_commit" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

(
    receipt_arg,
    receipt_sha,
    source_arg,
    source_sha,
    library_arg,
    library_sha,
    header_arg,
    header_sha,
    offline_main,
    release_name,
    expected_commit,
) = sys.argv[1:]
receipt_path = Path(receipt_arg).resolve(strict=True)
source_path = Path(source_arg).resolve(strict=True)
library_path = Path(library_arg).resolve(strict=True)
header_path = Path(header_arg).resolve(strict=True)

def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

for label, path, expected in (
    ("build receipt", receipt_path, receipt_sha),
    ("source manifest", source_path, source_sha),
    ("CaloReco library", library_path, library_sha),
    ("PhotonClusterBuilder header", header_path, header_sha),
):
    if len(expected) != 64 or digest(path) != expected:
        raise SystemExit(f"{label} hash drift")

receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
source = json.loads(source_path.read_text(encoding="utf-8"))
if receipt.get("schema") != "THE134_ANA560_CALORECO_BUILD_RECEIPT_V3":
    raise SystemExit("CaloReco build receipt schema differs")
if receipt.get("status") != "PASS":
    raise SystemExit("CaloReco build receipt is not PASS")
if source.get("schema") != "THE134_ANA560_CALORECO_SOURCE_MANIFEST_V2":
    raise SystemExit("CaloReco source manifest schema differs")
if receipt.get("runtime", {}).get("release") != release_name:
    raise SystemExit("CaloReco receipt release differs")
if receipt.get("runtime", {}).get("offline_main") != offline_main:
    raise SystemExit("CaloReco receipt OFFLINE_MAIN differs")
if receipt.get("build", {}).get("coresoftware_commit") != expected_commit:
    raise SystemExit("CaloReco receipt coresoftware commit differs")
if source.get("coresoftware", {}).get("commit") != expected_commit:
    raise SystemExit("CaloReco source-manifest commit differs")
if receipt.get("artifact", {}).get("sha256") != library_sha:
    raise SystemExit("CaloReco receipt library digest differs")
artifact = receipt_path.parent / str(receipt.get("artifact", {}).get("library", ""))
if artifact.resolve(strict=True) != library_path:
    raise SystemExit("CaloReco receipt library path does not resolve to the declared provider")
abi = receipt.get("abi", {})
if (
    abi.get("status") != "PASS"
    or abi.get("soname") != "libcalo_reco.so.0"
    or abi.get("soname_expected") != "libcalo_reco.so.0"
    or abi.get("needed_exact_match") is not True
    or abi.get("rpath_runpath_exact_match") is not True
    or abi.get("removed_symbols", {}).get("count") != 0
):
    raise SystemExit("CaloReco ABI/SONAME receipt contract differs")
for loader_name in ("libcalo_reco.so", "libcalo_reco.so.0"):
    loader_path = receipt_path.parent / "install/lib" / loader_name
    if not loader_path.is_symlink() or loader_path.resolve(strict=True) != library_path:
        raise SystemExit(f"CaloReco builder loader alias differs: {loader_name}")
receipt_source = receipt.get("source_manifest", {})
if receipt_source.get("sha256") != source_sha:
    raise SystemExit("CaloReco receipt/source-manifest digest cross-link differs")
if (receipt_path.parent / str(receipt_source.get("path", ""))).resolve(strict=True) != source_path:
    raise SystemExit("CaloReco receipt/source-manifest path cross-link differs")
single = receipt.get("single_provider", {})
provider_probe = single.get("provider_probe", {})
if (
    single.get("photon_cluster_builder_process_event_definitions") != 1
    or single.get("raw_cluster_builder_topo_process_event_definitions") != 1
    or single.get("forbidden_standalone_provider_count") != 0
    or single.get("other_installed_shared_objects") != []
    or provider_probe.get("status") != "PASS"
    or provider_probe.get("preload_provider_count") != 0
    or provider_probe.get("provider_count") != 1
    or Path(str(provider_probe.get("provider_realpath", ""))).resolve(strict=True)
    != library_path
):
    raise SystemExit("CaloReco receipt does not prove one complete provider")
runtime = receipt.get("runtime", {})
root_load = runtime.get("root_load", {})
if (
    runtime.get("ldd_not_found") is not False
    or runtime.get("mutable_user_dependency") is not False
    or root_load.get("status") != "PASS"
    or root_load.get("load_return_code") != 0
    or root_load.get("preload_provider_count") != 0
    or root_load.get("provider_count") != 1
    or Path(str(root_load.get("provider_realpath", ""))).resolve(strict=True)
    != library_path
):
    raise SystemExit("CaloReco runtime provider proof differs")
mapping = receipt.get("mapping_patch", {})
if (
    mapping.get("id")
    != "THE134_RAWCLUSTERBUILDERTOPO_DETECTOR_EXPLICIT_CHANNEL_MAP_V1"
    or mapping.get("scientific_controls_changed") != []
):
    raise SystemExit("CaloReco mapping authority differs")
if source.get("mapping_patch", {}).get("changed_scientific_controls") != []:
    raise SystemExit("CaloReco source manifest changes scientific controls")
if source.get("overlay", {}).get("PhotonClusterBuilder.h", {}).get("staged_sha256") != header_sha:
    raise SystemExit("CaloReco source/header overlay cross-link differs")
PY
}

# row_id|system|lane|dataset|sample|role|minimum_bias_gate|photon_id_row_match
emit_matrix() {
  printf '%s\n' \
    'pp_signal_photon5|pp|pp_photon_sim|isSim|run28_photonjet5|signal|not_applicable|newPPG12' \
    'pp_signal_photon10|pp|pp_photon_sim|isSim|run28_photonjet10|signal|not_applicable|newPPG12' \
    'pp_signal_photon20|pp|pp_photon_sim|isSim|run28_photonjet20|signal|not_applicable|newPPG12' \
    'pp_background_jet8|pp|pp_inclusive_sim|isSimInclusive|run28_jet8|background|not_applicable|newPPG12' \
    'pp_background_jet12|pp|pp_inclusive_sim|isSimInclusive|run28_jet12|background|not_applicable|newPPG12' \
    'pp_background_jet20|pp|pp_inclusive_sim|isSimInclusive|run28_jet20|background|not_applicable|newPPG12' \
    'pp_background_jet30|pp|pp_inclusive_sim|isSimInclusive|run28_jet30|background|not_applicable|newPPG12' \
    'auau_signal_photon12|auau|auau_photon_embedded|isSimEmbedded|run28_embeddedPhoton12|signal|required_pass|auauBDTSideband' \
    'auau_signal_photon20|auau|auau_photon_embedded|isSimEmbedded|run28_embeddedPhoton20|signal|required_pass|auauBDTSideband' \
    'auau_background_jet12|auau|auau_inclusive_embedded|isSimEmbeddedInclusive|run28_embeddedJet12|background|required_pass|auauBDTSideband' \
    'auau_background_jet20|auau|auau_inclusive_embedded|isSimEmbeddedInclusive|run28_embeddedJet20|background|required_pass|auauBDTSideband' \
    'auau_background_jet30|auau|auau_inclusive_embedded|isSimEmbeddedInclusive|run28_embeddedJet30|background|required_pass|auauBDTSideband' \
    'auau_background_jet40|auau|auau_inclusive_embedded|isSimEmbeddedInclusive|run28_embeddedJet40|background|required_pass|auauBDTSideband'
}

emit_execution_matrix() {
  if (( capacity_mode )); then
    emit_matrix | awk -F'|' \
      -v pp="$capacity_pp_row" -v auau="$capacity_auau_row" \
      '$1==pp || $1==auau'
  else
    emit_matrix
  fi
}

validate_matrix() {
  local rows unique_rows unique_samples
  rows="$(emit_matrix | wc -l | tr -d ' ')"
  unique_rows="$(emit_matrix | cut -d'|' -f1 | sort -u | wc -l | tr -d ' ')"
  unique_samples="$(emit_matrix | cut -d'|' -f5 | sort -u | wc -l | tr -d ' ')"
  [[ "$rows" == 13 && "$unique_rows" == 13 && "$unique_samples" == 13 ]] ||
    die "source matrix must contain exactly 13 unique rows and samples"
  if emit_matrix | grep -F 'run28_jet40' >/dev/null; then
    die "p+p Jet40 is diagnostic-only and must not enter this training extraction"
  fi
  [[ "$(emit_matrix | awk -F'|' '$2=="pp" && $6=="signal" {n++} END{print n+0}')" == 3 ]] ||
    die "p+p signal source closure failed"
  [[ "$(emit_matrix | awk -F'|' '$2=="pp" && $6=="background" {n++} END{print n+0}')" == 4 ]] ||
    die "p+p background source closure failed"
  [[ "$(emit_matrix | awk -F'|' '$2=="auau" && $6=="signal" {n++} END{print n+0}')" == 2 ]] ||
    die "Au+Au signal source closure failed"
  [[ "$(emit_matrix | awk -F'|' '$2=="auau" && $6=="background" {n++} END{print n+0}')" == 4 ]] ||
    die "Au+Au background source closure failed"
  [[ "$(emit_matrix | awk -F'|' '$2=="auau" && $7!="required_pass" {n++} END{print n+0}')" == 0 ]] ||
    die "every Au+Au source must require the MinimumBias classifier"
  [[ "$(emit_execution_matrix | wc -l | tr -d ' ')" == "$execution_row_count" ]] ||
    die "execution matrix row count differs from the selected mode"
  if (( capacity_mode )); then
    [[ "$(emit_execution_matrix | awk -F'|' '$1=="pp_background_jet8" && $2=="pp" && $5=="run28_jet8" {n++} END{print n+0}')" == 1 ]] ||
      die "capacity matrix must contain exactly the frozen p+p Jet8 witness"
    [[ "$(emit_execution_matrix | awk -F'|' '$1=="auau_background_jet12" && $2=="auau" && $5=="run28_embeddedJet12" {n++} END{print n+0}')" == 1 ]] ||
      die "capacity matrix must contain exactly the frozen Au+Au embedded-Jet12 witness"
  fi
}

validate_capacity_preflight_contract() {
  (( capacity_mode )) || return 0
  [[ -n "$capacity_canary_id" && ${#capacity_canary_id} -le 128 &&
     "$capacity_canary_id" =~ ^[A-Za-z0-9_.:-]+$ ]] ||
    die "capacity mode requires a safe RJ_THE134_CAPACITY_CANARY_ID"
  require_file_hash "capacity preflight receipt" \
    "$capacity_preflight_receipt" "$capacity_preflight_receipt_sha"
  require_file_hash "capacity full extraction plan" \
    "$capacity_full_plan" "$capacity_full_plan_sha"
  [[ "$output_root" =~ ^/sphenix/.*/replay_foundation/[A-Za-z0-9_.-]+$ ]] ||
    die "capacity output root must be a fresh isolated replay_foundation namespace"
  [[ "$evidence_root" =~ ^/sphenix/u/.*/evidence/qa/[A-Za-z0-9_.-]+$ ]] ||
    die "capacity evidence root must be a fresh isolated QA namespace"
  [[ "$submit_root" =~ ^/sphenix/u/.*/(condor|condor_sub)/[A-Za-z0-9_.-]+$ ]] ||
    die "capacity submit root must be a fresh isolated Condor namespace"
  python3 - \
    "$capacity_preflight_receipt" "$capacity_full_plan" \
    "$capacity_preflight_receipt_sha" "$capacity_full_plan_sha" \
    "$capacity_pp_row" "$capacity_auau_row" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

receipt_path, plan_path = map(Path, sys.argv[1:3])
receipt_sha, plan_sha, pp_row, auau_row = sys.argv[3:]
hex64 = re.compile(r"^[0-9a-f]{64}$")

def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

if digest(receipt_path) != receipt_sha or digest(plan_path) != plan_sha:
    raise SystemExit("capacity preflight input hash drift")
receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
plan = json.loads(plan_path.read_text(encoding="utf-8"))
if (
    receipt.get("schema")
    != "THE134_FULL_MULTIVIEW_EXTRACTION_PREFLIGHT_RECEIPT_V1"
    or receipt.get("status") != "PASS"
    or receipt.get("submission_performed") is not False
    or receipt.get("authority_state") != "PREFLIGHT_RESOLVED_NOT_EARNED"
    or int(receipt.get("full_training_authority", -1)) != 0
):
    raise SystemExit("capacity preflight receipt does not represent frozen pre-submit authority")
if (
    plan.get("schema") != "THE134_FULL_MULTIVIEW_EXTRACTION_PLAN_V1"
    or plan.get("status") != "PREFLIGHT_PASS"
    or plan.get("submission_performed") is not False
    or plan.get("execution_state") != "PREFLIGHT_ONLY_NO_CONDOR_MUTATION"
):
    raise SystemExit("capacity full extraction plan is not a frozen preflight plan")
if receipt.get("artifacts", {}).get("plan", {}).get("sha256") != plan_sha:
    raise SystemExit("capacity preflight receipt does not bind the exact full plan")
partition = plan.get("execution_partition", {})
if (
    partition.get("schema") != "THE134_FULL_EXTRACTION_PARTITION_CONTRACT_V1"
    or int(partition.get("group_size", -1)) != 7
    or partition.get("capacity_canary_required_before_submission") is not True
    or partition.get("capacity_authority_earned") is not False
):
    raise SystemExit("capacity execution partition is not the frozen group-of-seven pre-gate")
for key in (
    "bundle_manifest_sha256",
    "materialization_receipt_sha256",
    "execution_partition_sha256",
):
    if not hex64.fullmatch(str(receipt.get(key, ""))):
        raise SystemExit(f"capacity receipt has malformed identity: {key}")
if receipt.get("bundle_manifest_sha256") != plan.get("input_manifests", {}).get("bundle", {}).get("sha256"):
    raise SystemExit("capacity bundle receipt identity differs between receipt and plan")
if receipt.get("materialization_receipt_sha256") != plan.get("input_manifests", {}).get("materialization", {}).get("sha256"):
    raise SystemExit("capacity materialization identity differs between receipt and plan")
if plan.get("input_manifests", {}).get("materialization", {}).get("readback") != "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT":
    raise SystemExit("capacity plan lacks the exact immutable-bundle readback authority")
for namespace_key in ("output_root", "submit_root"):
    namespace = Path(str(plan.get("campaign", {}).get(namespace_key, "")))
    if not namespace.is_absolute() or namespace.exists():
        raise SystemExit(
            f"capacity gate requires an unmaterialized full-production namespace: "
            f"{namespace_key}={namespace}"
        )
rows = {str(row.get("row_id")): row for row in plan.get("rows", [])}
expected = {
    pp_row: ("pp", "run28_jet8"),
    auau_row: ("auau", "run28_embeddedJet12"),
}
for row_id, (system, sample) in expected.items():
    row = rows.get(row_id)
    if (
        row is None
        or row.get("system") != system
        or row.get("sample") != sample
        or int(row.get("input_contract", {}).get("group_size", -1)) != 7
        or int(row.get("full_training_authority", -1)) != 0
    ):
        raise SystemExit(f"capacity witness row differs from frozen full plan: {row_id}")
PY
}

capacity_receipt_value() {
  local key="$1"
  python3 - "$capacity_preflight_receipt" "$key" <<'PY'
import json
import sys
payload = json.load(open(sys.argv[1], encoding="utf-8"))
value = payload.get(sys.argv[2], "")
if not isinstance(value, str):
    raise SystemExit(f"capacity receipt field is not text: {sys.argv[2]}")
print(value)
PY
}

yaml_value() {
  local path="$1" key="$2"
  awk -F: -v key="$key" '$1 ~ "^[[:space:]]*" key "[[:space:]]*$" {sub(/^[^:]*:[[:space:]]*/, ""); sub(/[[:space:]]*#.*/, ""); print; exit}' "$path"
}

yaml_set_scalar() {
  local path="$1" key="$2" value="$3" tmp
  tmp="${path}.tmp.$$"
  awk -v key="$key" -v value="$value" '
    BEGIN { replaced=0 }
    $0 ~ "^[[:space:]]*" key "[[:space:]]*:" {
      if (!replaced) print key ": " value
      replaced=1
      next
    }
    { print }
    END { if (!replaced) print key ": " value }
  ' "$path" > "$tmp"
  mv "$tmp" "$path"
}

resolved_config_path() {
  printf '%s/%s.yaml\n' "$resolved_config_root" "$1"
}

prepare_resolved_configs() {
  mkdir -p "$resolved_config_root"
  while IFS='|' read -r row_id system _lane _dataset _sample role _mb_gate _row_match; do
    local source target
    target="$(resolved_config_path "$row_id")"
    if [[ "$system" == pp ]]; then
      source="$pp_config"
      cp "$source" "$target"
      yaml_set_scalar "$target" pp_photonid_extract_only true
      yaml_set_scalar "$target" pp_photonid_training_tree true
      yaml_set_scalar "$target" pp_photonid_training_tree_max_entries "$legacy_tree_max_entries"
      yaml_set_scalar "$target" pp_photonid_source_role "$role"
      yaml_set_scalar "$target" pp_photonid_ppg12_filter true
      yaml_set_scalar "$target" pp_photonid_require_preselection false
      yaml_set_scalar "$target" coneR "[${extraction_cone_r}]"
      [[ "$(yaml_value "$target" pp_photonid_source_role)" == "$role" ]] ||
        die "resolved p+p source role drift for ${row_id}"
    else
      source="$auau_config"
      cp "$source" "$target"
      yaml_set_scalar "$target" auau_bdt_training_tree true
      yaml_set_scalar "$target" auau_bdt_training_tree_max_entries "$legacy_tree_max_entries"
      yaml_set_scalar "$target" auau_bdt_npb_data_tagging false
      yaml_set_scalar "$target" coneR "[${extraction_cone_r}]"
      [[ "$(yaml_value "$target" auau_bdt_training_tree)" == true ]] ||
        die "resolved Au+Au training-tree gate drift for ${row_id}"
    fi
    [[ "$(yaml_value "$target" coneR)" == "[${extraction_cone_r}]" ]] ||
      die "${row_id} must resolve exactly one nominal R=${extraction_cone_r} fanout cone"
  done < <(emit_matrix)
}

photon_id_match_count() {
  local path="$1" row_match
  row_match="$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]')"
  awk -v needle="$row_match" '
    BEGIN {inside=0; count=0}
    /^[[:space:]]*photon_id_sets[[:space:]]*:/ {inside=1; next}
    inside && /^[^[:space:]]/ {inside=0}
    inside && /^[[:space:]]*-/ {
      row=tolower($0)
      if (index(row, needle)>0) ++count
    }
    END {print count+0}
  ' "$path"
}

verify_single_owner_resolved_contract() {
  local row_id system lane dataset sample role mb_gate row_match config matches
  while IFS='|' read -r row_id system lane dataset sample role mb_gate row_match; do
    config="$(resolved_config_path "$row_id")"
    [[ "$(yaml_value "$config" coneR)" == "[${extraction_cone_r}]" ]] ||
      die "${row_id} resolved config does not contain exactly one extraction cone"
    matches="$(photon_id_match_count "$config" "$row_match")"
    [[ "$matches" == 1 ]] ||
      die "${row_id} row match '${row_match}' must resolve exactly one photon-ID row, observed=${matches}"
  done < <(emit_matrix)
}

source_hashes_for_sample() {
  local row_id="$1" system="$2" sample="$3" sample_root="${sim_root}/${sample}"
  local -a names=(
    DST_CALO_CLUSTER.matched.list
    G4Hits.matched.list
    DST_JETS.matched.list
    DST_GLOBAL.matched.list
    DST_MBD_EPD.matched.list
  )
  local -a paths=()
  local name path expected_count="" count executable_count="" current_executable_count
  local first_tuple manifest_sha tuple_sha tuple_count
  for name in "${names[@]}"; do
    path="${sample_root}/${name}"
    [[ -s "$path" ]] || die "missing frozen source list for ${sample}: ${path}"
    count="$(wc -l < "$path" | tr -d ' ')"
    [[ "$count" =~ ^[1-9][0-9]*$ ]] || die "empty or malformed source list for ${sample}: ${path}"
    if [[ -z "$expected_count" ]]; then expected_count="$count"; fi
    [[ "$count" == "$expected_count" ]] || die "five-list row-count mismatch for ${sample}"
    current_executable_count="$(awk 'NF && $0 !~ /^[[:space:]]*#/ {n++} END {print n+0}' "$path")"
    [[ "$current_executable_count" =~ ^[1-9][0-9]*$ ]] ||
      die "source list has no executable entries for ${sample}: ${path}"
    if [[ -z "$executable_count" ]]; then executable_count="$current_executable_count"; fi
    [[ "$current_executable_count" == "$executable_count" ]] ||
      die "five-list executable-row-count mismatch for ${sample}"
    paths+=( "$path" )
  done
  manifest_sha="$({
    for path in "${paths[@]}"; do
      printf '%s\t%s\t%s\t%s\n' \
        "$(basename "$path")" "$(sha_file "$path")" \
        "$(wc -l < "$path" | tr -d ' ')" \
        "$(awk 'NF && $0 !~ /^[[:space:]]*#/ {n++} END {print n+0}' "$path")"
    done
  } | sha256_cmd | awk '{print $1}')"
  tuple_count="$(paste "${paths[@]}" | awk -F'\t' '
    function executable(value) {
      sub(/^[[:space:]]+/, "", value)
      return value != "" && substr(value, 1, 1) != "#"
    }
    {
      active=0
      for (column=1; column<=5; ++column) active += executable($column)
      if (active == 0) next
      if (active != 5) exit 41
      ++count
    }
    END {print count+0}
  ')" || die "five-list executable tuple alignment failed for ${sample}"
  [[ "$tuple_count" == "$executable_count" ]] ||
    die "executable tuple count differs from per-list executable count for ${sample}"
  first_tuple="$(paste "${paths[@]}" | awk -F'\t' '
    function executable(value) {
      sub(/^[[:space:]]+/, "", value)
      return value != "" && substr(value, 1, 1) != "#"
    }
    {
      active=0
      for (column=1; column<=5; ++column) active += executable($column)
      if (active == 5) {print $0; exit}
    }
  ')"
  [[ -n "$first_tuple" ]] || die "could not resolve the first executable five-file tuple for ${sample}"
  [[ "$(awk -F'\t' '{print NF}' <<< "$first_tuple")" == 5 ]] || die "first executable tuple is not five columns for ${sample}"
  validate_one_five_file_tuple "$row_id" "$system" <(printf '%s\n' "$first_tuple")
  tuple_sha="$(printf '%s\n' "$first_tuple" | sha256_cmd | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "$manifest_sha" "$tuple_sha" "$executable_count"
}

emit_source_hash_manifest() {
  local row_id system sample hashes
  validate_matrix
  printf 'row_id\tsource_manifest_sha256\tfirst_input_tuple_sha256\texecutable_tuple_count\n'
  while IFS='|' read -r row_id system _lane _dataset sample _role _mb_gate _row_match; do
    hashes="$(source_hashes_for_sample "$row_id" "$system" "$sample")"
    printf '%s\t%s\n' "$row_id" "$hashes"
  done < <(emit_matrix)
}

write_and_verify_source_hashes() {
  local expected_manifest_sha="$1" tmp
  require_file_hash "source hash manifest" "$source_hash_manifest" "$expected_manifest_sha"
  tmp="${observed_source_hashes}.tmp.$$"
  emit_source_hash_manifest > "$tmp"
  mv "$tmp" "$observed_source_hashes"
  cmp -s "$source_hash_manifest" "$observed_source_hashes" || {
    diff -u "$source_hash_manifest" "$observed_source_hashes" >&2 || true
    die "frozen source hash manifest does not match the 13 observed source populations"
  }
}

observed_source_field() {
  local row_id="$1" column="$2"
  awk -F'\t' -v row="$row_id" -v column="$column" '$1==row {print $column}' "$observed_source_hashes"
}

compute_code_sha() {
  local logical path
  while IFS='|' read -r logical path; do
    [[ -n "$logical" && -n "$path" ]] || die "malformed aggregate-code input"
    [[ -s "$path" ]] || die "missing aggregate-code input: ${logical} (${path})"
    printf '%s  %s\n' "$(sha_file "$path")" "$logical"
  done <<EOF | sha256_cmd | awk '{print $1}'
src/RecoilJets.cc|src/RecoilJets.cc
src/RecoilJets.h|src/RecoilJets.h
src_AuAu/RecoilJets_AuAu.cc|src_AuAu/RecoilJets_AuAu.cc
src_AuAu/RecoilJets_AuAu.h|src_AuAu/RecoilJets_AuAu.h
src/RJReplayFoundationV1.h|src/RJReplayFoundationV1.h
src/RJReplayRuntimeV1.h|src/RJReplayRuntimeV1.h
src/RJShowerFactorialV1.h|src/RJShowerFactorialV1.h
src/RJPhotonTrainingViewV1.h|src/RJPhotonTrainingViewV1.h
macros/Fun4All_recoilJets_unified_impl.C|macros/Fun4All_recoilJets_unified_impl.C
RecoilJets_Condor_submit.sh|${submitter}
RecoilJets_Condor.sh|${pp_executor}
RecoilJets_Condor_AuAu.sh|${auau_executor}
external/PhotonClusterv1.h|${photon_cluster_header}
external/PhotonClusterBuilder.h|${photon_cluster_builder_header}
external/libcalo_reco.so|${calo_reco_library}
external/calo_reco_build_receipt.json|${calo_reco_build_receipt}
external/calo_reco_source_manifest.json|${calo_reco_source_manifest}
external/ana.560/libcalo_io.so|${release_calo_io}
external/ana.560/libclusteriso.so|${release_clusteriso}
external/ana.560/libjetbase.so|${release_jetbase}
scripts/sdcc/workflows/diagnostics/submit_the134_multiview_extraction_smoke.sh|scripts/sdcc/workflows/diagnostics/submit_the134_multiview_extraction_smoke.sh
EOF
}

compute_semantic_sha() {
  {
    sha256_cmd \
      src/RJReplayFoundationV1.h \
      src/RJReplayRuntimeV1.h \
      src/RJShowerFactorialV1.h \
      src/RJPhotonTrainingViewV1.h
    printf '%s\n' \
      'THE134|H70,H0,G70,G0,O70,O0,R70|nominal-domain=15<=ET<35|capture-min=5|writer-only=1|legacy-tree-max=0|source-complete=13'
  } | sha256_cmd | awk '{print $1}'
}

write_submission_manifest() {
  local code_sha="$1" replay_schema_sha="$2" training_schema_sha="$3" semantic_sha="$4"
  local tmp row_id system lane dataset sample role mb_gate row_match source_sha tuple_sha tuple_count
  local config config_sha library library_sha model model_sha sidecar row_output row_submit
  tmp="${submission_manifest}.tmp.$$"
  printf 'row_id\tsystem\tlane\tdataset\tsample\tsource_role\tminimum_bias_gate\tinput_files\tinput_jobs\tsource_manifest_sha256\tfirst_input_tuple_sha256\tresolved_config\tresolved_config_sha256\tlibrary\tlibrary_sha256\tmodel\tmodel_sha256\tcode_sha256\treplay_schema_sha256\ttraining_schema_sha256\tsemantic_sha256\tnominal_et_min_gev\tnominal_et_max_gev_exclusive\tloose_capture_et_min_gev\tlegacy_training_tree_max_entries\tevent_limit_per_job\tanalysis_output_namespace\tmultiview_sidecar\tsubmit_namespace\tscheduler_log_dir\tscheduler_stdout_dir\tscheduler_stderr_dir\texecutable_input_tuple_count\tfull_training_authority\tphoton_cluster_builder_header\tphoton_cluster_builder_header_sha256\tcalo_reco_library\tcalo_reco_library_sha256\trelease_core_lib_dir\trelease_core_lib64_dir\trelease_calo_io\trelease_calo_io_sha256\trelease_clusteriso\trelease_clusteriso_sha256\trelease_jetbase\trelease_jetbase_sha256\n' > "$tmp"
  while IFS='|' read -r row_id system lane dataset sample role mb_gate row_match; do
    source_sha="$(observed_source_field "$row_id" 2)"
    tuple_sha="$(observed_source_field "$row_id" 3)"
    tuple_count="$(observed_source_field "$row_id" 4)"
    require_sha "${row_id} source manifest" "$source_sha"
    require_sha "${row_id} first tuple" "$tuple_sha"
    config="$(resolved_config_path "$row_id")"
    config_sha="$(sha_file "$config")"
    if [[ "$system" == pp ]]; then
      library="$pp_library"
      library_sha="$RJ_THE134_PP_LIBRARY_SHA256"
      model="$pp_model"
      model_sha="$RJ_THE134_PP_MODEL_SHA256"
    else
      library="$auau_library"
      library_sha="$RJ_THE134_AUAU_LIBRARY_SHA256"
      model="$auau_model"
      model_sha="$RJ_THE134_AUAU_MODEL_SHA256"
    fi
    row_output="${output_root}/${row_id}"
    sidecar="${row_output}/${sample}/RJPhotonTrainingViewV1.root"
    row_submit="${submit_root}/${row_id}"
    [[ "$tuple_count" =~ ^[1-9][0-9]*$ ]] || die "${row_id} executable tuple count is invalid"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$row_id" "$system" "$lane" "$dataset" "$sample" "$role" "$mb_gate" \
      "$execution_group_size" \
      "$source_sha" "$tuple_sha" "$config" "$config_sha" "$library" "$library_sha" \
      "$model" "$model_sha" "$code_sha" "$replay_schema_sha" "$training_schema_sha" "$semantic_sha" \
      "$nominal_et_min" "$nominal_et_max" "$capture_et_min" "$legacy_tree_max_entries" "$event_limit_per_job" \
      "$row_output" "$sidecar" "$row_submit" \
      '/sphenix/u/patsfan753/scratch/thesisAnalysis/log' \
      '/sphenix/u/patsfan753/scratch/thesisAnalysis/stdout' \
      '/sphenix/u/patsfan753/scratch/thesisAnalysis/error' \
      "$tuple_count" "$full_training_authority" \
      "$photon_cluster_builder_header" "$RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER_SHA256" \
      "$calo_reco_library" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256" \
      "$release_core_lib_dir" "$release_core_lib64_dir" \
      "$release_calo_io" "$RJ_THE134_RELEASE_CALO_IO_SHA256" \
      "$release_clusteriso" "$RJ_THE134_RELEASE_CLUSTERISO_SHA256" \
      "$release_jetbase" "$RJ_THE134_RELEASE_JETBASE_SHA256" >> "$tmp"
  done < <(emit_matrix)
  mv "$tmp" "$submission_manifest"
  [[ "$(wc -l < "$submission_manifest" | tr -d ' ')" == 14 ]] || die "submission manifest row closure failed"
  awk -F'\t' 'NF != 46 {exit 1}' "$submission_manifest" ||
    die "submission manifest must contain exactly 46 tab-separated fields on every row"
  [[ "$(tail -n +2 "$submission_manifest" | cut -f1 | sort -u | wc -l | tr -d ' ')" == 13 ]] || die "submission manifest contains duplicate row identities"
  [[ "$(tail -n +2 "$submission_manifest" | cut -f34 | sort -u)" == "$full_training_authority" ]] ||
    die "smoke manifest must remain full_training_authority=0"
  [[ "$(tail -n +2 "$submission_manifest" | cut -f28 | sort -u | wc -l | tr -d ' ')" == 13 ]] ||
    die "every source row must own a unique multiview sidecar path"
  sha_file "$submission_manifest" > "$duplicate_fingerprint"
  : > "${evidence_root}/pp_multiview_sidecars.list"
  : > "${evidence_root}/auau_multiview_sidecars.list"
  while IFS='|' read -r row_id system _lane _dataset _sample _role _mb_gate _row_match; do
    if [[ "$system" == pp ]]; then
      manifest_field "$row_id" 28 >> "${evidence_root}/pp_multiview_sidecars.list"
    else
      manifest_field "$row_id" 28 >> "${evidence_root}/auau_multiview_sidecars.list"
    fi
  done < <(emit_execution_matrix)
  if (( capacity_mode )); then
    [[ "$(wc -l < "${evidence_root}/pp_multiview_sidecars.list" | tr -d ' ')" == 1 ]] ||
      die "capacity p+p sidecar manifest closure failed"
    [[ "$(wc -l < "${evidence_root}/auau_multiview_sidecars.list" | tr -d ' ')" == 1 ]] ||
      die "capacity Au+Au sidecar manifest closure failed"
  else
    [[ "$(wc -l < "${evidence_root}/pp_multiview_sidecars.list" | tr -d ' ')" == 7 ]] ||
      die "p+p sidecar manifest closure failed"
    [[ "$(wc -l < "${evidence_root}/auau_multiview_sidecars.list" | tr -d ' ')" == 6 ]] ||
      die "Au+Au sidecar manifest closure failed"
  fi
}

write_runtime_authority_manifest() {
  local action="${1:-ensure}"
  [[ "$action" == ensure || "$action" == verify ]] ||
    die "runtime-authority action must be ensure or verify"
  if [[ "$action" == ensure ]]; then
    if [[ -e "$runtime_authority_manifest" && ! -e "$runtime_authority_fingerprint" ]] ||
       [[ ! -e "$runtime_authority_manifest" && -e "$runtime_authority_fingerprint" ]]; then
      die "runtime-authority manifest/fingerprint pair is incomplete"
    fi
  fi
  if ! python3 - "$runtime_authority_manifest" "$action" \
    "$calo_reco_build_receipt" "$RJ_THE134_CALO_RECO_BUILD_RECEIPT_SHA256" \
    "$calo_reco_source_manifest" "$RJ_THE134_CALO_RECO_SOURCE_MANIFEST_SHA256" \
    "$calo_reco_library" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256" \
    "$release_core_lib_dir" "$release_core_lib64_dir" \
    "$release_calo_io" "$RJ_THE134_RELEASE_CALO_IO_SHA256" \
    "$release_clusteriso" "$RJ_THE134_RELEASE_CLUSTERISO_SHA256" \
    "$release_jetbase" "$RJ_THE134_RELEASE_JETBASE_SHA256" \
    "$pinned_release_name" "$pinned_offline_main" \
    "$pinned_calo_reco_soname" "$pp_period" <<'PY'
from pathlib import Path
import json
import os
import sys

(
    destination_arg,
    action,
    build_receipt,
    build_receipt_sha,
    source_manifest,
    source_manifest_sha,
    calo_reco,
    calo_reco_sha,
    release_lib,
    release_lib64,
    calo_io,
    calo_io_sha,
    clusteriso,
    clusteriso_sha,
    jetbase,
    jetbase_sha,
    release_name,
    offline_main,
    calo_reco_soname,
    pp_period,
) = sys.argv[1:]
destination = Path(destination_arg)
payload = {
    "calo_reco": {
        "build_receipt": build_receipt,
        "build_receipt_sha256": build_receipt_sha,
        "library": calo_reco,
        "library_sha256": calo_reco_sha,
        "soname": calo_reco_soname,
        "source_manifest": source_manifest,
        "source_manifest_sha256": source_manifest_sha,
    },
    "release": {
        "lib": release_lib,
        "lib64": release_lib64,
        "name": release_name,
        "offline_main": offline_main,
        "providers": {
            "libcalo_io.so": {"path": calo_io, "sha256": calo_io_sha},
            "libclusteriso.so": {"path": clusteriso, "sha256": clusteriso_sha},
            "libjetbase.so": {"path": jetbase, "sha256": jetbase_sha},
        },
    },
    "pp_sim_weight_contract": {
        "interaction": "SI",
        "mix_weight": "period_auto",
        "period": pp_period,
        "period_lumi_weight": True,
        "vertex_reweight": "period_auto",
    },
    "schema": "THE134_SINGLE_PROVIDER_RUNTIME_AUTHORITY_V2",
    "status": "PASS",
}
if action == "ensure":
    if destination.exists():
        observed = json.loads(destination.read_text(encoding="utf-8"))
        if observed != payload:
            raise SystemExit(
                "existing runtime authority manifest differs from current frozen authority"
            )
    else:
        temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, separators=(",", ": ")) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
elif action == "verify":
    observed = json.loads(destination.read_text(encoding="utf-8"))
    if observed != payload:
        raise SystemExit("runtime authority manifest differs from current frozen authority")
else:
    raise SystemExit(f"unsupported runtime authority action: {action}")
PY
  then
    die "runtime authority manifest does not match the current frozen provider and p+p weight contract"
    return 2
  fi
  if [[ "$action" == ensure ]]; then
    if [[ -e "$runtime_authority_fingerprint" ]]; then
      [[ -s "$runtime_authority_fingerprint" &&
         "$(sha_file "$runtime_authority_manifest")" == "$(cat "$runtime_authority_fingerprint")" ]] ||
        die "existing runtime-authority fingerprint differs"
    else
      sha_file "$runtime_authority_manifest" > "$runtime_authority_fingerprint"
    fi
  else
    [[ -s "$runtime_authority_fingerprint" ]] ||
      die "runtime-authority fingerprint is missing"
    [[ "$(sha_file "$runtime_authority_manifest")" == "$(cat "$runtime_authority_fingerprint")" ]] ||
      die "runtime-authority manifest fingerprint drift blocks resume"
  fi
}

require_inputs_and_hashes() {
  local actual_code actual_replay_schema actual_training_schema actual_semantic
  local pp_yaml_model auau_yaml_model
  validate_matrix
  validate_pp_sim_weight_contract "$pp_period"
  validate_capacity_preflight_contract
  [[ -x "$submitter" ]] || die "RecoilJets submitter is not executable: ${submitter}"
  [[ -x "$pp_executor" ]] || die "p+p RecoilJets executor is not executable: ${pp_executor}"
  [[ -x "$auau_executor" ]] || die "Au+Au RecoilJets executor is not executable: ${auau_executor}"
  [[ -n "$pp_library" ]] || die "RJ_THE134_PP_LIBRARY is required; mutable or historical default libraries are forbidden"
  [[ -n "$auau_library" ]] || die "RJ_THE134_AUAU_LIBRARY is required; mutable or historical default libraries are forbidden"
  [[ -n "${RJ_CODEX_CHAT_NAME:-}" ]] || die "RJ_CODEX_CHAT_NAME is required"
  [[ -n "${RJ_CODEX_THREAD_ID:-}" ]] || die "RJ_CODEX_THREAD_ID is required"
  [[ -n "$source_hash_manifest" ]] || die "RJ_THE134_SOURCE_HASH_MANIFEST is required"
  [[ -n "$photon_cluster_builder_header" ]] || die "RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER is required"
  [[ -n "$calo_reco_library" ]] || die "RJ_THE134_CALO_RECO_LIBRARY is required"
  [[ -n "$calo_reco_build_receipt" ]] || die "RJ_THE134_CALO_RECO_BUILD_RECEIPT is required"
  [[ -n "$calo_reco_source_manifest" ]] || die "RJ_THE134_CALO_RECO_SOURCE_MANIFEST is required"
  [[ -d "$release_core_lib_dir" && -d "$release_core_lib64_dir" ]] ||
    die "the exact ana.560 release lib/lib64 directories are required"
  release_core_lib_dir="$(cd "$release_core_lib_dir" && pwd -P)"
  release_core_lib64_dir="$(cd "$release_core_lib64_dir" && pwd -P)"
  [[ "$release_core_lib_dir" == "${pinned_offline_main}/lib" ]] ||
    die "RJ_THE134_RELEASE_CORE_LIB_DIR must be exactly ${pinned_offline_main}/lib"
  [[ "$release_core_lib64_dir" == "${pinned_offline_main}/lib64" ]] ||
    die "RJ_THE134_RELEASE_CORE_LIB64_DIR must be exactly ${pinned_offline_main}/lib64"
  [[ -z "${RJ_THE134_MULTIVIEW_TRAINING_FILE:-}" ]] ||
    die "controller environment must not pre-own RJ_THE134_MULTIVIEW_TRAINING_FILE"
  [[ -z "${RJ_ID_FANOUT_FILE:-}" && -z "${RJ_ID_FANOUT_DIRS_FILE:-}" ]] ||
    die "controller environment must not carry a pre-existing fanout owner"
  [[ -z "${RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE:-}" ]] ||
    die "a separate PhotonClusterBuilder override library is forbidden; the pinned CaloReco library is authoritative"
  [[ -z "${RJ_FORCE_RELEASE_CORE_LIBS:-}" ]] ||
    die "controller environment must not inherit RJ_FORCE_RELEASE_CORE_LIBS"
  [[ -z "${RJ_FORCE_RELEASE_CALO_IO:-}" ]] ||
    die "controller environment must not inherit RJ_FORCE_RELEASE_CALO_IO"
  [[ -z "${RJ_RELEASE_CALO_IO_PATH:-}" ]] ||
    die "controller environment must not inherit RJ_RELEASE_CALO_IO_PATH"

  : "${RJ_THE134_PP_LIBRARY_SHA256:?set frozen p+p library SHA-256}"
  : "${RJ_THE134_AUAU_LIBRARY_SHA256:?set frozen Au+Au library SHA-256}"
  : "${RJ_THE134_PP_CONFIG_SHA256:?set frozen p+p base-config SHA-256}"
  : "${RJ_THE134_AUAU_CONFIG_SHA256:?set frozen Au+Au base-config SHA-256}"
  : "${RJ_THE134_PP_MODEL_SHA256:?set frozen p+p model SHA-256}"
  : "${RJ_THE134_AUAU_MODEL_SHA256:?set frozen Au+Au model SHA-256}"
  : "${RJ_THE134_CODE_SHA256:?set frozen aggregate code SHA-256}"
  : "${RJ_THE134_REPLAY_SCHEMA_SHA256:?set frozen replay-schema SHA-256}"
  : "${RJ_THE134_TRAINING_SCHEMA_SHA256:?set frozen training-schema SHA-256}"
  : "${RJ_THE134_SEMANTIC_SHA256:?set frozen aggregate semantic SHA-256}"
  : "${RJ_THE134_SOURCE_HASH_MANIFEST_SHA256:?set frozen source-hash-manifest SHA-256}"
  : "${RJ_THE134_PHOTON_CLUSTER_HEADER_SHA256:?set frozen PhotonClusterv1 header SHA-256}"
  : "${RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER_SHA256:?set frozen PhotonClusterBuilder header SHA-256}"
  : "${RJ_THE134_CALO_RECO_LIBRARY_SHA256:?set frozen CaloReco library SHA-256}"
  : "${RJ_THE134_CALO_RECO_BUILD_RECEIPT_SHA256:?set frozen CaloReco build-receipt SHA-256}"
  : "${RJ_THE134_CALO_RECO_SOURCE_MANIFEST_SHA256:?set frozen CaloReco source-manifest SHA-256}"
  : "${RJ_THE134_RELEASE_CALO_IO_SHA256:?set frozen ana.560 libcalo_io SHA-256}"
  : "${RJ_THE134_RELEASE_CLUSTERISO_SHA256:?set frozen ana.560 libclusteriso SHA-256}"
  : "${RJ_THE134_RELEASE_JETBASE_SHA256:?set frozen ana.560 libjetbase SHA-256}"

  require_file_hash "p+p library" "$pp_library" "$RJ_THE134_PP_LIBRARY_SHA256"
  require_file_hash "Au+Au library" "$auau_library" "$RJ_THE134_AUAU_LIBRARY_SHA256"
  require_file_hash "p+p base config" "$pp_config" "$RJ_THE134_PP_CONFIG_SHA256"
  require_file_hash "Au+Au base config" "$auau_config" "$RJ_THE134_AUAU_CONFIG_SHA256"
  require_file_hash "p+p model" "$pp_model" "$RJ_THE134_PP_MODEL_SHA256"
  require_file_hash "Au+Au model" "$auau_model" "$RJ_THE134_AUAU_MODEL_SHA256"
  require_file_hash "PhotonClusterv1 build header" "$photon_cluster_header" "$RJ_THE134_PHOTON_CLUSTER_HEADER_SHA256"
  require_file_hash "PhotonClusterBuilder build header" "$photon_cluster_builder_header" "$RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER_SHA256"
  require_file_hash "CaloReco runtime library" "$calo_reco_library" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256"
  require_file_hash "CaloReco build receipt" "$calo_reco_build_receipt" "$RJ_THE134_CALO_RECO_BUILD_RECEIPT_SHA256"
  require_file_hash "CaloReco source manifest" "$calo_reco_source_manifest" "$RJ_THE134_CALO_RECO_SOURCE_MANIFEST_SHA256"
  release_calo_io="$(
    resolve_release_companion RELEASE_CALO_IO libcalo_io.so \
      "$release_calo_io" "$RJ_THE134_RELEASE_CALO_IO_SHA256"
  )"
  release_clusteriso="$(
    resolve_release_companion RELEASE_CLUSTERISO libclusteriso.so \
      "$release_clusteriso" "$RJ_THE134_RELEASE_CLUSTERISO_SHA256"
  )"
  release_jetbase="$(
    resolve_release_companion RELEASE_JETBASE libjetbase.so \
      "$release_jetbase" "$RJ_THE134_RELEASE_JETBASE_SHA256"
  )"
  validate_calo_reco_build_authority

  pp_yaml_model="$(yaml_value "$pp_config" tight_bdt_model_file)"
  auau_yaml_model="$(yaml_value "$auau_config" auau_tight_bdt_centInputBase3x3_model_file)"
  [[ "$pp_yaml_model" == "$pp_model" ]] || die "p+p config/model path mismatch: config=${pp_yaml_model} frozen=${pp_model}"
  [[ "$auau_yaml_model" == "$auau_model" ]] || die "Au+Au config/model path mismatch: config=${auau_yaml_model} frozen=${auau_model}"

  actual_code="$(compute_code_sha)"
  actual_replay_schema="$(sha_file src/RJReplayFoundationV1.h)"
  actual_training_schema="$(sha_text "$training_schema_text")"
  actual_semantic="$(compute_semantic_sha)"
  [[ "$actual_code" == "$RJ_THE134_CODE_SHA256" ]] || die "aggregate code hash drift: expected=${RJ_THE134_CODE_SHA256} actual=${actual_code}"
  [[ "$actual_replay_schema" == "$RJ_THE134_REPLAY_SCHEMA_SHA256" ]] || die "replay-schema hash drift"
  [[ "$actual_training_schema" == "$RJ_THE134_TRAINING_SCHEMA_SHA256" ]] || die "training-schema hash drift"
  [[ "$actual_semantic" == "$RJ_THE134_SEMANTIC_SHA256" ]] || die "aggregate semantic hash drift"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$actual_code" "$actual_replay_schema" "$actual_training_schema" "$actual_semantic" \
    "$release_core_lib_dir" "$release_core_lib64_dir" \
    "$release_calo_io" "$release_clusteriso" "$release_jetbase"
}

preflight() {
  local hashes code_sha replay_schema_sha training_schema_sha semantic_sha
  [[ ! -e "$submission_journal" && ! -e "$submission_receipt" ]] ||
    die "submission evidence already exists; use status or validate instead of rewriting preflight state"
  hashes="$(require_inputs_and_hashes)"
  IFS=$'\t' read -r code_sha replay_schema_sha training_schema_sha semantic_sha \
    release_core_lib_dir release_core_lib64_dir release_calo_io release_clusteriso release_jetbase <<< "$hashes"
  mkdir -p "$evidence_root"
  prepare_resolved_configs
  verify_single_owner_resolved_contract
  write_and_verify_source_hashes "$RJ_THE134_SOURCE_HASH_MANIFEST_SHA256"
  write_submission_manifest "$code_sha" "$replay_schema_sha" "$training_schema_sha" "$semantic_sha"
  write_runtime_authority_manifest ensure
  bash -n "$0"
  bash -n "$submitter"
  bash -n "$pp_executor"
  bash -n "$auau_executor"
  if (( capacity_mode )); then
    say "CAPACITY_PREFLIGHT_PASS selected_rows=${execution_row_count} group_size=${execution_group_size} full_manifest_rows=13 manifest=${submission_manifest} fingerprint=$(cat "$duplicate_fingerprint")"
  else
    say "PREFLIGHT_PASS rows=13 manifest=${submission_manifest} fingerprint=$(cat "$duplicate_fingerprint")"
  fi
}

manifest_field() {
  local row_id="$1" column="$2"
  awk -F'\t' -v row="$row_id" -v column="$column" '$1==row {print $column}' "$submission_manifest"
}

write_source_provenance() {
  command -v python3 >/dev/null 2>&1 || die "python3 is required to materialize source provenance"
  python3 - \
    "$submission_receipt" "$submission_manifest" \
    "$source_provenance_json" "$pp_source_provenance_json" "$auau_source_provenance_json" <<'PY'
import csv
import json
import os
import sys
from pathlib import Path

receipt_path, manifest_path, aggregate_path, pp_path, auau_path = map(Path, sys.argv[1:])
with manifest_path.open(newline="") as stream:
    manifest_rows = {row["row_id"]: row for row in csv.DictReader(stream, delimiter="\t")}
with receipt_path.open(newline="") as stream:
    receipt_rows = list(csv.DictReader(stream, delimiter="\t"))

records = []
seen_paths = set()
for receipt in receipt_rows:
    row = manifest_rows.get(receipt["row_id"])
    if row is None:
        raise SystemExit(f"receipt row absent from submission manifest: {receipt['row_id']}")
    sidecar = receipt["multiview_sidecar"]
    if sidecar in seen_paths:
        raise SystemExit(f"duplicate sidecar path in receipt: {sidecar}")
    seen_paths.add(sidecar)
    # Existing RecoilJets workers define both legacy input provenance branches
    # as the SHA-256 of their immutable staged five-file chunk list.  Preserve
    # that exact worker contract rather than inventing submit-time DST hashes.
    chunk_sha = receipt["staged_chunk_sha256"]
    records.append(
        {
            "path": sidecar,
            "row_id": receipt["row_id"],
            "system": row["system"],
            "source_sample": row["sample"],
            "input_uri_sha256": chunk_sha,
            "input_file_sha256": chunk_sha,
            "source_manifest_sha256": row["source_manifest_sha256"],
            "config_sha256": row["resolved_config_sha256"],
            "code_sha256": row["code_sha256"],
        }
    )

def write_payload(path: Path, selected: list[dict]) -> None:
    payload = {"schema": "THE134_SOURCE_PROVENANCE_V1", "inputs": selected}
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)

write_payload(aggregate_path, records)
write_payload(pp_path, [record for record in records if record["system"] == "pp"])
write_payload(auau_path, [record for record in records if record["system"] == "auau"])
PY
}

common_extra_env() {
  local row_id="$1" system="$2" lane="$3" dataset="$4" sample="$5" role="$6"
  local source_sha config config_sha model_sha sidecar model_score shower_definition si_di_role period
  source_sha="$(manifest_field "$row_id" 10)"
  config="$(manifest_field "$row_id" 12)"
  config_sha="$(manifest_field "$row_id" 13)"
  model_sha="$(manifest_field "$row_id" 17)"
  sidecar="$(manifest_field "$row_id" 28)"
  if [[ "$system" == pp ]]; then
    model_score=tight_bdt_score
    shower_definition=H70
    si_di_role=SI
    period="$pp_period"
  else
    model_score=auau_tight_bdt_score
    shower_definition=H0
    si_di_role=EMBEDDED
    period=AUAU_RUN24
  fi
  local extra
  extra="RJ_REPLAY_FOUNDATION_V1=1;RJ_REPLAY_FOUNDATION_CANARY=1;RJ_REPLAY_TRACE=0;RJ_REPLAY_LANE=${lane};RJ_REPLAY_DATASET=${dataset};RJ_REPLAY_SAMPLE=${sample};RJ_REPLAY_PERIOD=${period};RJ_REPLAY_SI_DI_ROLE=${si_di_role};RJ_REPLAY_OWNERSHIP_STATE=source_role_frozen;RJ_REPLAY_SOURCE_MANIFEST_SHA256=${source_sha};RJ_REPLAY_SOURCE_SHA256=${source_sha};RJ_REPLAY_MODEL_SHA256=${model_sha};RJ_REPLAY_MODEL_SCORE_NAME=${model_score};RJ_REPLAY_MODEL_SHOWER_DEFINITION=${shower_definition};RJ_REPLAY_CONFIG_SHA256=${config_sha};RJ_REPLAY_CODE_SHA256=${RJ_THE134_CODE_SHA256};RJ_REPLAY_SCHEMA_SHA256=${RJ_THE134_REPLAY_SCHEMA_SHA256};RJ_REPLAY_SEMANTIC_SHA256=${RJ_THE134_SEMANTIC_SHA256};RJ_REPLAY_PHOTON_CAPTURE_ET_MIN=${capture_et_min};RJ_REPLAY_JET_CONSTITUENT_PT_MIN=${capture_et_min};RJ_THE134_MULTIVIEW_TRAINING_V1=1;RJ_THE134_MULTIVIEW_TRAINING_FILE=${sidecar};RJ_THE134_EXPECTED_SOURCE_ROLE=${role}"
  if (( capacity_mode )); then
    extra="${extra};RJ_REPLAY_FOUNDATION_CAPACITY_CANARY=1;RJ_REPLAY_FOUNDATION_CAPACITY_CANARY_ID=${capacity_canary_id};RJ_REPLAY_FOUNDATION_CAPACITY_PREFLIGHT_RECEIPT=${capacity_preflight_receipt};RJ_REPLAY_FOUNDATION_CAPACITY_PREFLIGHT_RECEIPT_SHA256=${capacity_preflight_receipt_sha};RJ_REPLAY_FOUNDATION_EXECUTION_PARTITION_SHA256=$(capacity_receipt_value execution_partition_sha256);RJ_REPLAY_FOUNDATION_BUNDLE_RECEIPT_SHA256=$(capacity_receipt_value bundle_manifest_sha256);RJ_REPLAY_FOUNDATION_MATERIALIZATION_RECEIPT_SHA256=$(capacity_receipt_value materialization_receipt_sha256)"
  fi
  printf '%s' "$extra"
}

system_extra_env() {
  local system="$1" role="$2"
  if [[ "$system" == pp ]]; then
    printf '%s' \
      "RJ_PP_PHOTONID_EXTRACT_ONLY=1;RJ_PP_PHOTONID_TRAINING_TREE=1;RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES=${legacy_tree_max_entries};RJ_PP_PHOTONID_SOURCE_ROLE=${role};RJ_PP_PHOTONID_PPG12_FILTER=1;RJ_PP_PHOTONID_REQUIRE_PRESELECTION=0;RJ_SIM_ALLOW_NONE_LISTS=0"
  else
    printf '%s' \
      "RJ_AUAU_BDT_EXTRACT_ONLY=1;RJ_AUAU_BDT_TRAINING_TREE=1;RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES=${legacy_tree_max_entries};RJ_AUAU_BDT_NPB_DATA_TAGGING=0;RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=1;RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0;RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0;RJ_SIM_ALLOW_NONE_LISTS=1"
  fi
}

condor_field() {
  local path="$1" key="$2"
  awk -F= -v key="$key" '
    $1 ~ "^[[:space:]]*" key "[[:space:]]*$" {
      sub(/^[^=]*=[[:space:]]*/, ""); sub(/[[:space:]]*$/, ""); print; exit
    }
  ' "$path"
}

resolve_condor_template() {
  local value="$1" cluster="$2" proc="$3"
  local cluster_token='$(Cluster)' cluster_id_token='$(ClusterId)'
  local proc_token='$(Process)' proc_id_token='$(ProcId)'
  value="${value//$cluster_token/$cluster}"
  value="${value//$cluster_id_token/$cluster}"
  value="${value//$proc_token/$proc}"
  value="${value//$proc_id_token/$proc}"
  printf '%s\n' "$value"
}

descriptor_env_values() {
  local path="$1" key="$2" environment
  environment="$(condor_field "$path" environment)"
  awk -v environment="$environment" -v key="$key" '
    BEGIN {
      count=split(environment, fields, ";")
      prefix=key "="
      for (idx=1; idx<=count; ++idx) {
        field=fields[idx]
        sub(/^[[:space:]]+/, "", field)
        sub(/[[:space:]]+$/, "", field)
        if (substr(field, 1, length(prefix)) == prefix)
          print substr(field, length(prefix)+1)
      }
    }
  '
}

require_descriptor_env_exact() {
  local row_id="$1" path="$2" key="$3" expected="$4" value
  local -a values=()
  while IFS= read -r value; do values+=( "$value" ); done \
    < <(descriptor_env_values "$path" "$key")
  [[ "${#values[@]}" == 1 && "${values[0]}" == "$expected" ]] ||
    die "${row_id} descriptor must bind ${key}=${expected} exactly once"
}

require_descriptor_env_absent() {
  local row_id="$1" path="$2" key="$3" value
  local -a values=()
  while IFS= read -r value; do values+=( "$value" ); done \
    < <(descriptor_env_values "$path" "$key")
  [[ "${#values[@]}" == 0 ]] ||
    die "${row_id} descriptor must leave ${key} unset for frozen automatic authority"
}

analysis_tag_for_dataset() {
  case "$1" in
    isSimInclusive) echo isSimInclusive ;;
    isSimEmbedded) echo isSimEmbedded ;;
    isSimEmbeddedInclusive) echo isSimEmbeddedInclusive ;;
    *) echo isSim ;;
  esac
}

verify_sealed_snapshot_receipts() {
  local snapshot_dir="$1" loader_receipt="$2" snapshot_manifest="$3" expected_mode="$4"
  python3 - \
    "$snapshot_dir" "$loader_receipt" "$snapshot_manifest" "$expected_mode" \
    "$pinned_calo_reco_soname" \
    "$calo_reco_library" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256" \
    "$release_core_lib_dir" "$release_core_lib64_dir" \
    "$release_calo_io" "$RJ_THE134_RELEASE_CALO_IO_SHA256" \
    "$release_clusteriso" "$RJ_THE134_RELEASE_CLUSTERISO_SHA256" \
    "$release_jetbase" "$RJ_THE134_RELEASE_JETBASE_SHA256" <<'PY'
from pathlib import Path
import hashlib
import json
import os
import sys

(
    snapshot_arg,
    loader_arg,
    manifest_arg,
    expected_mode,
    expected_calo_soname,
    _calo_source,
    calo_sha,
    release_lib,
    release_lib64,
    calo_io,
    calo_io_sha,
    clusteriso,
    clusteriso_sha,
    jetbase,
    jetbase_sha,
) = sys.argv[1:]
snapshot = Path(snapshot_arg).resolve(strict=True)
loader_path = Path(loader_arg).resolve(strict=True)
manifest_path = Path(manifest_arg).resolve(strict=True)
if snapshot not in loader_path.parents or snapshot not in manifest_path.parents:
    raise SystemExit("snapshot receipt or manifest escaped the frozen snapshot")
if snapshot.stat().st_mode & 0o222:
    raise SystemExit(f"writable frozen snapshot root survived seal: {snapshot}")
for evidence_path in (loader_path, manifest_path):
    if evidence_path.stat().st_mode & 0o222:
        raise SystemExit(f"writable snapshot evidence survived seal: {evidence_path}")

def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

loader = json.loads(loader_path.read_text(encoding="utf-8"))
if (
    loader.get("schema") != "RJ_SNAPSHOT_LOADER_RECEIPT_V1"
    or loader.get("status") != "PASS"
    or loader.get("pinned_calo_reco_release_companions") is not True
    or loader.get("mode") != expected_mode
):
    raise SystemExit("snapshot loader receipt contract differs")
if loader.get("release_roots") != [release_lib64, release_lib]:
    raise SystemExit("snapshot loader receipt release roots differ")
providers = loader.get("providers", {})
calo_api = snapshot / "lib/libcalo_reco.so"
calo_soname = snapshot / "lib" / expected_calo_soname
if (
    not calo_api.is_file()
    or not calo_soname.is_symlink()
    or calo_soname.readlink().as_posix() != "libcalo_reco.so"
    or calo_soname.resolve(strict=True) != calo_api.resolve(strict=True)
):
    raise SystemExit("snapshot CaloReco API/SONAME alias closure differs")
expected = {
    "libcalo_reco.so": (calo_api.resolve(strict=True), calo_sha),
    "libcalo_io.so": (Path(calo_io).resolve(strict=True), calo_io_sha),
    "libclusteriso.so": (Path(clusteriso).resolve(strict=True), clusteriso_sha),
    "libjetbase.so": (Path(jetbase).resolve(strict=True), jetbase_sha),
}
for family, (path, expected_sha) in expected.items():
    provider = providers.get(family, {})
    if (
        provider.get("realpath") != str(path)
        or provider.get("sha256") != expected_sha
        or provider.get("observed_resolutions") != [str(path)]
    ):
        raise SystemExit(f"snapshot loader provider authority differs: {family}")
    if digest(path) != expected_sha:
        raise SystemExit(f"snapshot loader provider hash drift: {family}")

manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
if (
    manifest.get("schema") != "RJ_FROZEN_SNAPSHOT_MANIFEST_V1"
    or manifest.get("status") != "PASS"
    or manifest.get("root") != str(snapshot)
):
    raise SystemExit("snapshot manifest contract differs")
observed = []
for path in sorted(snapshot.rglob("*"), key=lambda item: item.relative_to(snapshot).as_posix()):
    if path == manifest_path:
        continue
    relative = path.relative_to(snapshot).as_posix()
    if path.is_symlink():
        target = os.readlink(path)
        if os.path.isabs(target):
            raise SystemExit(f"snapshot symlink must be relative: {relative}")
        resolved = path.resolve(strict=True)
        if snapshot not in resolved.parents:
            raise SystemExit(f"snapshot symlink escaped: {relative}")
        if path.parent.stat().st_mode & 0o222:
            raise SystemExit(f"snapshot symlink parent is writable: {relative}")
        if resolved.stat().st_mode & 0o222:
            raise SystemExit(f"snapshot symlink target is writable: {relative}")
        observed.append(
            {
                "path": relative,
                "sha256": digest(resolved),
                "symlink_target": target,
                "type": "symlink",
            }
        )
    elif path.is_file():
        if path.lstat().st_mode & 0o222:
            raise SystemExit(f"writable artifact survived snapshot seal: {relative}")
        observed.append(
            {
                "path": relative,
                "sha256": digest(path),
                "size": path.stat().st_size,
                "type": "file",
            }
        )
    elif path.is_dir():
        if path.lstat().st_mode & 0o222:
            raise SystemExit(f"writable directory survived snapshot seal: {relative}")
        observed.append({"path": relative, "type": "directory"})
    else:
        raise SystemExit(f"unsupported snapshot artifact: {relative}")
if manifest.get("entries") != observed:
    raise SystemExit("complete snapshot inventory/hash closure differs")
if manifest_path.lstat().st_mode & 0o222:
    raise SystemExit("snapshot manifest itself is writable")
PY
}

validate_five_file_tuple_count() {
  local row_id="$1" system="$2" path="$3" expected_rows="$4"
  [[ "$expected_rows" =~ ^[1-9][0-9]*$ ]] ||
    die "${row_id} expected staged tuple count must be a positive integer"
  awk -F '\t' -v tuple_system="$system" -v expected_rows="$expected_rows" '
    /^[[:space:]]*($|#)/ { next }
    {
      rows += 1
      if (NF != 5) bad = 1
      for (column = 1; column <= NF; ++column) {
        if ($column == "") bad = 1
        if (tuple_system == "pp" && ($column == "NONE" || substr($column, 1, 1) != "/")) bad = 1
        if (tuple_system == "auau" && column <= 4 && ($column == "NONE" || substr($column, 1, 1) != "/")) bad = 1
      }
      if (tuple_system == "auau" && $5 != "NONE") bad = 1
      if (tuple_system != "pp" && tuple_system != "auau") bad = 1
    }
    END { exit (bad || rows != expected_rows) }
  ' "$path" ||
    die "${row_id} staged chunk violates the typed ${system} five-column ${expected_rows}-tuple contract"
}

validate_one_five_file_tuple() {
  validate_five_file_tuple_count "$1" "$2" "$3" 1
}

validate_five_field_fanout_contract() {
  local row_id="$1" fanout_file="$2" materialized_config="$3" expected_cone="$4"
  local cone_count materialized_cone
  [[ -s "$fanout_file" ]] || {
    die "${row_id} materialized fanout contract is missing: ${fanout_file}"
    return 2
  }
  awk -F '|' '
    /^[[:space:]]*($|#)/ { next }
    {
      rows += 1
      if (NF != 5) bad = 1
      for (column = 1; column <= NF; ++column)
        if ($column == "") bad = 1
    }
    END { exit (bad || rows != 1) }
  ' "$fanout_file" || {
    die "${row_id} fanout contract must contain exactly one nonempty five-field ID row"
    return 2
  }
  [[ -s "$materialized_config" ]] || {
    die "${row_id} descriptor-bound materialized config is missing: ${materialized_config}"
    return 2
  }
  cone_count="$(grep -Ec '^[[:space:]]*coneR[[:space:]]*:' "$materialized_config" || true)"
  [[ "$cone_count" == 1 ]] || {
    die "${row_id} materialized config must contain exactly one coneR authority"
    return 2
  }
  materialized_cone="$(yaml_value "$materialized_config" coneR)"
  [[ "$materialized_cone" == "$expected_cone" || "$materialized_cone" == "0.4" ]] || {
    die "${row_id} materialized config does not carry the nominal R=${expected_cone} extraction cone"
    return 2
  }
}

# Verify the exact dry-materialized unit that will be submitted.  This is the
# pre-submission ownership proof: one descriptor, one argument row, one
# fanout row, one RecoilJets instance, and one multiview sidecar assignment.
# The tab-separated return value is consumed verbatim by submit_row.
verify_materialized_row_contract() {
  local row_id="$1" system="$2" dataset="$3" sample="$4" row_output="$5" sidecar="$6" row_submit="$7"
  local -a sub_files=() sidecar_values=() id_file_values=() id_dirs_values=() materialized_config_values=() allow_none_values=()
  local submit_file args_file args_line chunk_list chunk_sha fanout_file fanout_sha fanout_line
  local arg_sample arg_chunk arg_dataset arg_cluster arg_events arg_index arg_none arg_dest arg_extra
  local fan_dest fan_cfg fan_pre fan_tight fan_non materialized_config materialized_config_sha
  local chunk_tag analysis_tag analysis_root owner_count queue_count value
  local frozen_executable snapshot_dir snapshot_header snapshot_header_sha
  local snapshot_calo snapshot_calo_sha snapshot_analysis snapshot_analysis_sha expected_analysis_sha
  local snapshot_calo_soname
  local snapshot_impl snapshot_calo_macro expected_loader_suffix companion
  local snapshot_loader_receipt snapshot_loader_receipt_sha snapshot_manifest snapshot_manifest_sha
  local sealed_getenv

  while IFS= read -r value; do sub_files+=( "$value" ); done \
    < <(find "$row_submit" -maxdepth 1 -type f -name '*.sub' -print | sort)
  [[ "${#sub_files[@]}" == 1 ]] || die "${row_id} must dry-materialize exactly one submit descriptor, observed=${#sub_files[@]}"
  submit_file="${sub_files[0]}"
  sealed_getenv="$(condor_field "$submit_file" getenv)"
  [[ "$sealed_getenv" == False ]] ||
    die "${row_id} descriptor must disable submit-host environment inheritance"
  frozen_executable="$(condor_field "$submit_file" executable)"
  [[ -s "$frozen_executable" ]] ||
    die "${row_id} descriptor does not bind a readable frozen executable"
  snapshot_dir="$(cd "$(dirname "$frozen_executable")" && pwd -P)"
  snapshot_header="${snapshot_dir}/PhotonClusterBuilder.h"
  snapshot_calo="${snapshot_dir}/lib/libcalo_reco.so"
  snapshot_calo_soname="${snapshot_dir}/lib/${pinned_calo_reco_soname}"
  snapshot_impl="${snapshot_dir}/Fun4All_recoilJets_unified_impl.C"
  snapshot_calo_macro="${snapshot_dir}/Calo_Calib.C"
  snapshot_loader_receipt="${snapshot_dir}/snapshot_loader_receipt.json"
  snapshot_manifest="${snapshot_dir}/snapshot_manifest.json"
  expected_loader_suffix="snapshot_loader_suffix=\":${release_core_lib64_dir}:${release_core_lib_dir}\""
  if [[ "$system" == pp ]]; then
    snapshot_analysis="${snapshot_dir}/lib/libRecoilJets.so"
    expected_analysis_sha="$RJ_THE134_PP_LIBRARY_SHA256"
  else
    snapshot_analysis="${snapshot_dir}/lib/libRecoilJetsAuAu.so"
    expected_analysis_sha="$RJ_THE134_AUAU_LIBRARY_SHA256"
  fi
  require_file_hash "${row_id} snapshotted PhotonClusterBuilder header" \
    "$snapshot_header" "$RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER_SHA256"
  require_file_hash "${row_id} snapshotted CaloReco library" \
    "$snapshot_calo" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256"
  [[ -L "$snapshot_calo_soname" &&
     "$(readlink "$snapshot_calo_soname")" == "libcalo_reco.so" &&
     "$snapshot_calo_soname" -ef "$snapshot_calo" ]] ||
    die "${row_id} snapshotted CaloReco SONAME alias does not resolve to its one provider"
  require_file_hash "${row_id} snapshotted CaloReco SONAME alias" \
    "$snapshot_calo_soname" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256"
  require_file_hash "${row_id} snapshotted analysis library" \
    "$snapshot_analysis" "$expected_analysis_sha"
  [[ ! -e "${snapshot_dir}/lib/libphoton_cluster_builder_override.so" ]] ||
    die "${row_id} snapshot contains an undeclared PhotonClusterBuilder override library"
  for companion in libcalo_io.so libclusteriso.so libjetbase.so; do
    [[ ! -e "${snapshot_dir}/lib/${companion}" ]] ||
      die "${row_id} snapshot illegally duplicates release-owned ${companion}"
  done
  [[ "$(grep -Fxc "$expected_loader_suffix" "$frozen_executable" || true)" == 1 ]] ||
    die "${row_id} frozen executor lacks the exact ana.560 loader suffix"
  [[ "$(grep -Ec '^[[:space:]]*# RJ_PINNED_SPHENIX_RELEASE_V1[[:space:]]*$' "$frozen_executable" || true)" == 1 &&
     "$(grep -Ec "^[[:space:]]*source /opt/sphenix/core/bin/sphenix_setup\\.sh -n ${pinned_release_name}[[:space:]]*$" "$frozen_executable" || true)" == 1 &&
     "$(grep -Fc "$pinned_offline_main" "$frozen_executable" || true)" -ge 2 ]] ||
    die "${row_id} frozen executor lacks the exact ${pinned_release_name} runtime witness"
  [[ "$(grep -Fxc "R__LOAD_LIBRARY(${snapshot_calo})" "$snapshot_calo_macro" || true)" == 1 ]] ||
    die "${row_id} Calo_Calib macro does not load the same snapshotted CaloReco provider"
  for companion in libclusteriso.so libjetbase.so; do
    [[ "$(grep -Fxc "R__LOAD_LIBRARY(${companion})" "$snapshot_impl" || true)" == 1 ]] ||
      die "${row_id} unified macro does not bind ${companion} exactly once through ana.560"
  done
  snapshot_header_sha="$(sha_file "$snapshot_header")"
  snapshot_calo_sha="$(sha_file "$snapshot_calo")"
  snapshot_analysis_sha="$(sha_file "$snapshot_analysis")"
  verify_sealed_snapshot_receipts \
    "$snapshot_dir" "$snapshot_loader_receipt" "$snapshot_manifest" "$system"
  snapshot_loader_receipt_sha="$(sha_file "$snapshot_loader_receipt")"
  snapshot_manifest_sha="$(sha_file "$snapshot_manifest")"
  args_file="${submit_file%.sub}.args"
  [[ -s "$args_file" && "$(wc -l < "$args_file" | tr -d ' ')" == 1 ]] ||
    die "${row_id} must dry-materialize exactly one argument row"
  queue_count="$(grep -Fxc "queue arguments from ${args_file}" "$submit_file" || true)"
  [[ "$queue_count" == 1 ]] || die "${row_id} descriptor must queue exactly the frozen one-row args file"

  args_line="$(sed -n '1p' "$args_file")"
  read -r arg_sample arg_chunk arg_dataset arg_cluster arg_events arg_index arg_none arg_dest arg_extra <<< "$args_line"
  [[ -z "${arg_extra:-}" ]] || die "${row_id} argument row has unexpected extra fields"
  [[ "$arg_sample" == "$sample" && "$arg_dataset" == "$dataset" ]] ||
    die "${row_id} materialized sample/dataset drift"
  [[ "$arg_cluster" == '$(Cluster)' && "$arg_events" == "$event_limit_per_job" && "$arg_index" == 1 && "$arg_none" == NONE ]] ||
    die "${row_id} argument row violates the one-proc group-${execution_group_size} contract"
  [[ "$arg_dest" == "$row_output/"* ]] || die "${row_id} argument destination escapes its owned output namespace"
  chunk_list="$arg_chunk"
  [[ -s "$chunk_list" ]] || die "${row_id} staged chunk list is missing: ${chunk_list}"
  validate_five_file_tuple_count \
    "$row_id" "$system" "$chunk_list" "$execution_group_size"
  chunk_sha="$(sha_file "$chunk_list")"
  require_sha "${row_id} staged chunk" "$chunk_sha"

  while IFS= read -r value; do sidecar_values+=( "$value" ); done \
    < <(descriptor_env_values "$submit_file" RJ_THE134_MULTIVIEW_TRAINING_FILE)
  [[ "${#sidecar_values[@]}" == 1 && "${sidecar_values[0]}" == "$sidecar" ]] ||
    die "${row_id} descriptor must contain exactly one owner of its declared multiview sidecar"
  while IFS= read -r value; do id_file_values+=( "$value" ); done \
    < <(descriptor_env_values "$submit_file" RJ_ID_FANOUT_FILE)
  while IFS= read -r value; do id_dirs_values+=( "$value" ); done \
    < <(descriptor_env_values "$submit_file" RJ_ID_FANOUT_DIRS_FILE)
  [[ "${#id_file_values[@]}" == 1 && "${#id_dirs_values[@]}" == 1 && "${id_file_values[0]}" == "${id_dirs_values[0]}" ]] ||
    die "${row_id} descriptor must bind one identical fanout file/dirs contract"
  while IFS= read -r value; do materialized_config_values+=( "$value" ); done \
    < <(descriptor_env_values "$submit_file" RJ_CONFIG_YAML)
  [[ "${#materialized_config_values[@]}" == 1 ]] ||
    die "${row_id} descriptor must bind exactly one materialized RJ_CONFIG_YAML"
  materialized_config="${materialized_config_values[0]}"
  while IFS= read -r value; do allow_none_values+=( "$value" ); done \
    < <(descriptor_env_values "$submit_file" RJ_SIM_ALLOW_NONE_LISTS)
  [[ "${#allow_none_values[@]}" == 1 ]] ||
    die "${row_id} descriptor must bind exactly one typed RJ_SIM_ALLOW_NONE_LISTS value"
  if (( capacity_mode )); then
    require_descriptor_env_exact "$row_id" "$submit_file" \
      RJ_REPLAY_FOUNDATION_CAPACITY_CANARY 1
    require_descriptor_env_exact "$row_id" "$submit_file" \
      RJ_REPLAY_FOUNDATION_CAPACITY_CANARY_ID "$capacity_canary_id"
    require_descriptor_env_exact "$row_id" "$submit_file" \
      RJ_REPLAY_FOUNDATION_CAPACITY_PREFLIGHT_RECEIPT "$capacity_preflight_receipt"
    require_descriptor_env_exact "$row_id" "$submit_file" \
      RJ_REPLAY_FOUNDATION_CAPACITY_PREFLIGHT_RECEIPT_SHA256 "$capacity_preflight_receipt_sha"
    require_descriptor_env_exact "$row_id" "$submit_file" \
      RJ_REPLAY_FOUNDATION_EXECUTION_PARTITION_SHA256 \
      "$(capacity_receipt_value execution_partition_sha256)"
    require_descriptor_env_exact "$row_id" "$submit_file" \
      RJ_REPLAY_FOUNDATION_BUNDLE_RECEIPT_SHA256 \
      "$(capacity_receipt_value bundle_manifest_sha256)"
    require_descriptor_env_exact "$row_id" "$submit_file" \
      RJ_REPLAY_FOUNDATION_MATERIALIZATION_RECEIPT_SHA256 \
      "$(capacity_receipt_value materialization_receipt_sha256)"
  else
    require_descriptor_env_absent "$row_id" "$submit_file" \
      RJ_REPLAY_FOUNDATION_CAPACITY_CANARY
  fi
  if [[ "$system" == pp ]]; then
    [[ "${allow_none_values[0]}" == 0 ]] || die "${row_id} p+p descriptor must reject NONE lists"
    require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PHOTON_YIELD 1
    require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PHOTON_YIELD_DOUBLE 0
    require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD "$pp_period"
    require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_USE_LUMI_WEIGHT 1
    require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_STRICT_DI 0
    require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_ALLOW_ALL_SIM 0
    require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE 0
    require_descriptor_env_exact "$row_id" "$submit_file" RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE 0
    require_descriptor_env_absent "$row_id" "$submit_file" RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT
    require_descriptor_env_absent "$row_id" "$submit_file" RJ_PP_VERTEX_REWEIGHT_FILE
  else
    [[ "${allow_none_values[0]}" == 1 ]] || die "${row_id} Au+Au descriptor must authorize only its typed optional MBD list"
  fi
  fanout_file="${id_file_values[0]}"
  validate_five_field_fanout_contract "$row_id" "$fanout_file" "$materialized_config" "$extraction_cone_r"
  fanout_line="$(grep -Ev '^[[:space:]]*($|#)' "$fanout_file")"
  IFS='|' read -r fan_dest fan_cfg fan_pre fan_tight fan_non <<< "$fanout_line"
  [[ -n "$fan_dest" && -n "$fan_cfg" && -n "$fan_pre" && -n "$fan_tight" && -n "$fan_non" ]] ||
    die "${row_id} fanout row is incomplete"
  [[ "$fan_dest" == "$row_output/"* ]] || die "${row_id} fanout destination escapes its owned output namespace"
  materialized_config_sha="$(sha_file "$materialized_config")"
  require_sha "${row_id} materialized config" "$materialized_config_sha"
  owner_count=$(( ${#sub_files[@]} * $(wc -l < "$args_file" | tr -d ' ') * $(grep -Evc '^[[:space:]]*($|#)' "$fanout_file") * ${#sidecar_values[@]} ))
  [[ "$owner_count" == 1 ]] || die "${row_id} multiview sidecar owner multiplicity is ${owner_count}, expected 1"

  chunk_tag="$(basename "$chunk_list")"
  chunk_tag="${chunk_tag%.list}"
  analysis_tag="$(analysis_tag_for_dataset "$dataset")"
  analysis_root="${fan_dest}/${sample}/RecoilJets_${analysis_tag}_${fan_cfg}_${chunk_tag}.root"
  fanout_sha="$(sha_file "$fanout_file")"
  require_sha "${row_id} fanout contract" "$fanout_sha"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$submit_file" "$args_file" "$chunk_list" "$chunk_sha" \
    "$fanout_file" "$fanout_sha" "$analysis_root" "$owner_count" \
    "$snapshot_dir" "$snapshot_header" "$snapshot_header_sha" \
    "$snapshot_calo" "$snapshot_calo_sha" "$snapshot_analysis" "$snapshot_analysis_sha" \
    "$materialized_config" "$materialized_config_sha" \
    "$snapshot_loader_receipt" "$snapshot_loader_receipt_sha" \
    "$snapshot_manifest" "$snapshot_manifest_sha"
}

verify_materialization_attempt_seal() {
  local attempt_dir="${1:?materialization attempt directory required}"
  local contract_file="${attempt_dir}/materialization_contract.tsv"
  local fingerprint_file="${contract_file}.sha256"
  local expected actual

  [[ -d "$attempt_dir" && ! -L "$attempt_dir" ]] ||
    die "materialization attempt is not a real directory: ${attempt_dir}"
  [[ -f "$contract_file" && ! -L "$contract_file" && -s "$contract_file" ]] ||
    die "materialization attempt contract is missing or unsafe: ${contract_file}"
  [[ -f "$fingerprint_file" && ! -L "$fingerprint_file" && -s "$fingerprint_file" ]] ||
    die "materialization attempt fingerprint is missing or unsafe: ${fingerprint_file}"
  [[ "$(wc -l < "$fingerprint_file" | tr -d ' ')" == 1 ]] ||
    die "materialization attempt fingerprint must contain exactly one row"
  expected="$(cat "$fingerprint_file")"
  require_sha "materialization attempt fingerprint" "$expected"
  actual="$(sha_file "$contract_file")"
  [[ "$actual" == "$expected" ]] ||
    die "materialization attempt contract fingerprint drift: ${attempt_dir}"

  python3 - "$attempt_dir" <<'PY'
from pathlib import Path
import os
import stat
import sys

root = Path(sys.argv[1]).resolve(strict=True)

def assert_sealed(path: Path) -> None:
    mode = path.lstat().st_mode
    if stat.S_ISLNK(mode):
        target_text = os.readlink(path)
        target = Path(target_text)
        if target.is_absolute():
            raise SystemExit(f"materialization symlink must be relative: {path}")
        resolved = (path.parent / target).resolve(strict=True)
        try:
            resolved.relative_to(root)
        except ValueError as exc:
            raise SystemExit(
                f"materialization symlink escapes attempt root: {path} -> {target_text}"
            ) from exc
        parent_mode = path.parent.stat().st_mode
        target_mode = resolved.stat().st_mode
        if parent_mode & 0o222:
            raise SystemExit(f"materialization symlink parent is writable: {path.parent}")
        if target_mode & 0o222:
            raise SystemExit(f"materialization symlink target is writable: {resolved}")
        return
    if not (stat.S_ISREG(mode) or stat.S_ISDIR(mode)):
        raise SystemExit(f"materialization attempt contains a special file: {path}")
    if mode & 0o222:
        raise SystemExit(f"materialization attempt contains writable state: {path}")

assert_sealed(root)
for candidate in root.rglob("*"):
    assert_sealed(candidate)
PY
}

select_materialization_attempt() {
  local row_submit_root="${1:?row submit root required}"
  local entry name suffix contract_file fingerprint_file
  local completed="" completed_suffix="" incomplete_suffixes=""

  mkdir -p "$row_submit_root"
  [[ -d "$row_submit_root" && ! -L "$row_submit_root" ]] ||
    die "row submit root is not a real directory: ${row_submit_root}"

  while IFS= read -r -d '' entry; do
    name="$(basename "$entry")"
    case "$name" in
      attempt_01|attempt_02|attempt_03) ;;
      *) die "unexpected row-materialization artifact blocks retry: ${entry}" ;;
    esac
    [[ -d "$entry" && ! -L "$entry" ]] ||
      die "materialization attempt is not a real directory: ${entry}"
    contract_file="${entry}/materialization_contract.tsv"
    fingerprint_file="${contract_file}.sha256"
    if [[ -e "$contract_file" || -e "$fingerprint_file" ]]; then
      [[ -e "$contract_file" && -e "$fingerprint_file" ]] ||
        die "materialization attempt has an incomplete seal receipt: ${entry}"
      verify_materialization_attempt_seal "$entry"
      [[ -z "$completed" ]] ||
        die "multiple sealed materialization attempts block exact reuse: ${row_submit_root}"
      completed="$entry"
      completed_suffix="${name#attempt_}"
    else
      incomplete_suffixes="${incomplete_suffixes} ${name#attempt_}"
    fi
  done < <(find "$row_submit_root" -mindepth 1 -maxdepth 1 -print0)

  if [[ -n "$completed" ]]; then
    for suffix in $incomplete_suffixes; do
      (( 10#$suffix < 10#$completed_suffix )) ||
        die "an incomplete materialization attempt is newer than the sealed authority: ${row_submit_root}/attempt_${suffix}"
    done
    printf 'reuse\t%s\n' "$completed"
    return 0
  fi

  for suffix in 01 02 03; do
    entry="${row_submit_root}/attempt_${suffix}"
    if [[ ! -e "$entry" && ! -L "$entry" ]]; then
      mkdir "$entry"
      printf 'new\t%s\n' "$entry"
      return 0
    fi
  done
  die "materialization retry budget exhausted with three preserved incomplete attempts: ${row_submit_root}"
}

seal_materialization_attempt() {
  local attempt_dir="${1:?materialization attempt directory required}"
  local contract="${2:?materialization contract required}"
  local contract_file="${attempt_dir}/materialization_contract.tsv"
  local fingerprint_file="${contract_file}.sha256"

  [[ "$contract" != *$'\n'* ]] ||
    die "materialization contract must be one logical row"
  [[ ! -e "$contract_file" && ! -L "$contract_file" &&
     ! -e "$fingerprint_file" && ! -L "$fingerprint_file" ]] ||
    die "materialization attempt seal already exists: ${attempt_dir}"
  printf '%s\n' "$contract" > "$contract_file"
  sha_file "$contract_file" > "$fingerprint_file"
  python3 - "$attempt_dir" <<'PY'
from pathlib import Path
import os
import stat
import sys

root = Path(sys.argv[1]).resolve(strict=True)
for dirpath, dirnames, filenames in os.walk(root, topdown=False, followlinks=False):
    parent = Path(dirpath)
    for name in filenames + dirnames:
        path = parent / name
        mode = path.lstat().st_mode
        if stat.S_ISLNK(mode):
            continue
        os.chmod(path, stat.S_IMODE(mode) & ~0o222)
root_mode = root.lstat().st_mode
os.chmod(root, stat.S_IMODE(root_mode) & ~0o222)
PY
  verify_materialization_attempt_seal "$attempt_dir"
}

submit_row() {
  local row_id="$1" system="$2" lane="$3" dataset="$4" sample="$5" role="$6" _mb_gate="$7" row_match="$8"
  local config row_output sidecar row_submit_root row_submit library extra row_log materialize_log cluster proc cluster_proc
  local materialization_selection materialization_disposition stored_contract
  local contract submit_file args_file chunk_list chunk_sha fanout_file fanout_sha analysis_root owner_count
  local snapshot_dir snapshot_header snapshot_header_sha snapshot_calo snapshot_calo_sha
  local snapshot_analysis snapshot_analysis_sha materialized_config materialized_config_sha
  local snapshot_loader_receipt snapshot_loader_receipt_sha snapshot_manifest snapshot_manifest_sha
  local log_template out_template err_template log_path out_path err_path submitted_args expected_args args_sha
  local -a materialize_contract_env=(
    RJ_REPLAY_FOUNDATION_CANARY=1
    RJ_REPLAY_LANE="$lane"
    RJ_REPLAY_SCHEMA_SHA256="$RJ_THE134_REPLAY_SCHEMA_SHA256"
  )
  if (( capacity_mode )); then
    materialize_contract_env+=(
      RJ_REPLAY_FOUNDATION_CAPACITY_CANARY=1
      RJ_REPLAY_FOUNDATION_CAPACITY_CANARY_ID="$capacity_canary_id"
      RJ_REPLAY_FOUNDATION_CAPACITY_PREFLIGHT_RECEIPT="$capacity_preflight_receipt"
      RJ_REPLAY_FOUNDATION_CAPACITY_PREFLIGHT_RECEIPT_SHA256="$capacity_preflight_receipt_sha"
      RJ_REPLAY_FOUNDATION_EXECUTION_PARTITION_SHA256="$(capacity_receipt_value execution_partition_sha256)"
      RJ_REPLAY_FOUNDATION_BUNDLE_RECEIPT_SHA256="$(capacity_receipt_value bundle_manifest_sha256)"
      RJ_REPLAY_FOUNDATION_MATERIALIZATION_RECEIPT_SHA256="$(capacity_receipt_value materialization_receipt_sha256)"
    )
  fi
  config="$(manifest_field "$row_id" 12)"
  row_output="$(manifest_field "$row_id" 27)"
  sidecar="$(manifest_field "$row_id" 28)"
  row_submit_root="$(manifest_field "$row_id" 29)"
  if [[ "$system" == pp ]]; then library="$pp_library"; else library="$auau_library"; fi
  extra="$(common_extra_env "$row_id" "$system" "$lane" "$dataset" "$sample" "$role");$(system_extra_env "$system" "$role")"
  materialization_selection="$(select_materialization_attempt "$row_submit_root")"
  IFS=$'\t' read -r materialization_disposition row_submit <<< "$materialization_selection"
  [[ "$materialization_disposition" == new || "$materialization_disposition" == reuse ]] ||
    die "${row_id} materialization selection returned an invalid disposition"
  [[ -n "$row_submit" ]] || die "${row_id} materialization selection returned no attempt directory"
  mkdir -p "$(dirname "$sidecar")"
  materialize_log="${evidence_root}/materialize_${row_id}_$(basename "$row_submit").log"
  row_log="${evidence_root}/submit_${row_id}.log"

  if [[ "$materialization_disposition" == reuse ]]; then
    say "MATERIALIZE_REUSE row=${row_id} attempt=$(basename "$row_submit")"
    verify_materialization_attempt_seal "$row_submit"
    stored_contract="$(cat "${row_submit}/materialization_contract.tsv")"
    contract="$(verify_materialized_row_contract "$row_id" "$system" "$dataset" "$sample" "$row_output" "$sidecar" "$row_submit")"
    [[ "$contract" == "$stored_contract" ]] ||
      die "${row_id} sealed materialization contract differs from exact readback"
  else
    say "MATERIALIZE row=${row_id} dataset=${dataset} sample=${sample} role=${role} attempt=$(basename "$row_submit")"
    if [[ "$system" == pp ]]; then
      env -u RJ_FORCE_RELEASE_CORE_LIBS -u RJ_FORCE_RELEASE_CALO_IO -u RJ_RELEASE_CALO_IO_PATH \
        -u RJ_PPG12_CROSSING_PERIOD -u RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT \
        -u RJ_PP_VERTEX_REWEIGHT_FILE -u RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX \
        -u RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX \
        -u RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX \
        "${materialize_contract_env[@]}" \
        RJ_CONDOR_SEALED_ENVIRONMENT=1 \
        RJ_CODEX_CHAT_NAME="$RJ_CODEX_CHAT_NAME" RJ_CODEX_THREAD_ID="$RJ_CODEX_THREAD_ID" \
        RJ_SIM_ROOT_OVERRIDE="$sim_root" RJ_CONFIG_YAML="$config" RJ_PP_LIBRARY_OVERRIDE="$library" \
        RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE="$photon_cluster_builder_header" \
        RJ_CALO_RECO_LIBRARY_OVERRIDE="$calo_reco_library" \
        RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE= \
        RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS=1 \
        RJ_PINNED_CALO_RECO_SONAME="$pinned_calo_reco_soname" \
        RJ_PINNED_RELEASE_NAME="$pinned_release_name" \
        RJ_PINNED_OFFLINE_MAIN="$pinned_offline_main" \
        RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io" \
        RJ_PINNED_RELEASE_CALO_IO_SHA256="$RJ_THE134_RELEASE_CALO_IO_SHA256" \
        RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso" \
        RJ_PINNED_RELEASE_CLUSTERISO_SHA256="$RJ_THE134_RELEASE_CLUSTERISO_SHA256" \
        RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase" \
        RJ_PINNED_RELEASE_JETBASE_SHA256="$RJ_THE134_RELEASE_JETBASE_SHA256" \
        RJ_RELEASE_CORE_LIB_DIR="$release_core_lib_dir" \
        RJ_RELEASE_CORE_LIB64_DIR="$release_core_lib64_dir" \
        RJ_SIM_ALLOW_NONE_LISTS=0 \
        RJ_PPG12_PHOTON_YIELD=1 \
        RJ_PPG12_PHOTON_YIELD_DOUBLE=0 \
        RJ_PPG12_PERIOD="$pp_period" \
        RJ_PPG12_PERIOD_USE_LUMI_WEIGHT=1 \
        RJ_PPG12_PERIOD_STRICT_DI=0 \
        RJ_PPG12_PERIOD_ALLOW_ALL_SIM=0 \
        RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE=0 \
        RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE=0 \
        RJ_PP_PHOTONID_EXTRACT_ONLY=1 RJ_PP_PHOTONID_TRAINING_TREE=1 \
        RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES="$legacy_tree_max_entries" \
        RJ_PP_PHOTONID_SOURCE_ROLE="$role" RJ_PP_PHOTONID_PPG12_FILTER=1 \
        RJ_PP_PHOTONID_REQUIRE_PRESELECTION=0 \
        RJ_DAG_DRYRUN=1 \
        RJ_AUTO_MERGE=0 RJ_STAGE_EMAIL_MODE=none RJ_CLEAN_OUTPUT_BASE=0 \
        RJ_REQUEST_MEMORY=8000MB RJ_REQUIRE_NON_TINY_OUTPUT=1 RJ_MIN_OUTPUT_BYTES=50000 \
        RJ_FAIL_ON_MISSING_CALO_INPUT=1 RJ_VALIDATE_SIM_INPUT_PATHS=1 RJ_VALIDATE_SIM_INPUT_MAX_LINES=1 \
        RJ_PROFILE_JOB=1 RJ_JOB_HEARTBEAT_SECONDS=120 RJ_PROFILE_LABEL="${tag}_${row_id}" \
        RJ_SMOKE_OUTPUT_BASE="$row_output" RJ_SMOKE_SIM_NEVENTS="$event_limit_per_job" \
        RJ_SUBMISSION_NAMESPACE="$row_id" RJ_CONDOR_SUB_DIR="$row_submit" \
        RJ_PHOTON_ID_ROW_MATCH="$row_match" RJ_ID_FANOUT_MAX_ROWS=1 \
        RJ_SUBMIT_EXTRA_ENV="$extra" \
        "$submitter" "$dataset" condorDoAllSmoke groupSize "$execution_group_size" maxJobs 1 "SAMPLE=${sample}" \
        2>&1 | tee "$materialize_log"
    else
      env -u RJ_FORCE_RELEASE_CORE_LIBS -u RJ_FORCE_RELEASE_CALO_IO -u RJ_RELEASE_CALO_IO_PATH \
        "${materialize_contract_env[@]}" \
        RJ_CONDOR_SEALED_ENVIRONMENT=1 \
        RJ_CODEX_CHAT_NAME="$RJ_CODEX_CHAT_NAME" RJ_CODEX_THREAD_ID="$RJ_CODEX_THREAD_ID" \
        RJ_SIM_ROOT_OVERRIDE="$sim_root" RJ_CONFIG_YAML="$config" RJ_AUAU_LIBRARY_OVERRIDE="$library" \
        RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE="$photon_cluster_builder_header" \
        RJ_CALO_RECO_LIBRARY_OVERRIDE="$calo_reco_library" \
        RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE= \
        RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS=1 \
        RJ_PINNED_CALO_RECO_SONAME="$pinned_calo_reco_soname" \
        RJ_PINNED_RELEASE_NAME="$pinned_release_name" \
        RJ_PINNED_OFFLINE_MAIN="$pinned_offline_main" \
        RJ_PINNED_RELEASE_CALO_IO_PATH="$release_calo_io" \
        RJ_PINNED_RELEASE_CALO_IO_SHA256="$RJ_THE134_RELEASE_CALO_IO_SHA256" \
        RJ_PINNED_RELEASE_CLUSTERISO_PATH="$release_clusteriso" \
        RJ_PINNED_RELEASE_CLUSTERISO_SHA256="$RJ_THE134_RELEASE_CLUSTERISO_SHA256" \
        RJ_PINNED_RELEASE_JETBASE_PATH="$release_jetbase" \
        RJ_PINNED_RELEASE_JETBASE_SHA256="$RJ_THE134_RELEASE_JETBASE_SHA256" \
        RJ_RELEASE_CORE_LIB_DIR="$release_core_lib_dir" \
        RJ_RELEASE_CORE_LIB64_DIR="$release_core_lib64_dir" \
        RJ_SIM_ALLOW_NONE_LISTS=1 \
        RJ_DAG_DRYRUN=1 \
        RJ_AUTO_MERGE=0 RJ_STAGE_EMAIL_MODE=none RJ_CLEAN_OUTPUT_BASE=0 \
        RJ_REQUEST_MEMORY=8000MB RJ_REQUIRE_NON_TINY_OUTPUT=1 RJ_MIN_OUTPUT_BYTES=50000 \
        RJ_FAIL_ON_MISSING_CALO_INPUT=1 RJ_VALIDATE_SIM_INPUT_PATHS=1 RJ_VALIDATE_SIM_INPUT_MAX_LINES=1 \
        RJ_PROFILE_JOB=1 RJ_JOB_HEARTBEAT_SECONDS=120 RJ_PROFILE_LABEL="${tag}_${row_id}" \
        RJ_SMOKE_OUTPUT_BASE="$row_output" RJ_SMOKE_SIM_NEVENTS="$event_limit_per_job" \
        RJ_SUBMISSION_NAMESPACE="$row_id" RJ_CONDOR_SUB_DIR="$row_submit" \
        RJ_PHOTON_ID_ROW_MATCH="$row_match" RJ_ID_FANOUT_MAX_ROWS=1 \
        RJ_SUBMIT_EXTRA_ENV="$extra" \
        "$submitter" "$dataset" condorDoAllSmoke groupSize "$execution_group_size" maxJobs 1 "SAMPLE=${sample}" \
        2>&1 | tee "$materialize_log"
    fi

    grep -F 'RECOILJETS_SMOKETEST_DRYRUN_V1' "$materialize_log" >/dev/null ||
      die "${row_id} submitter did not attest dry materialization"
    ! grep -Eq '[0-9]+ job\(s\) submitted to cluster [0-9]+' "$materialize_log" ||
      die "${row_id} dry materialization unexpectedly submitted a Condor job"
    contract="$(verify_materialized_row_contract "$row_id" "$system" "$dataset" "$sample" "$row_output" "$sidecar" "$row_submit")"
    seal_materialization_attempt "$row_submit" "$contract"
  fi
  IFS=$'\t' read -r submit_file args_file chunk_list chunk_sha fanout_file fanout_sha analysis_root owner_count \
    snapshot_dir snapshot_header snapshot_header_sha snapshot_calo snapshot_calo_sha \
    snapshot_analysis snapshot_analysis_sha materialized_config materialized_config_sha \
    snapshot_loader_receipt snapshot_loader_receipt_sha snapshot_manifest snapshot_manifest_sha <<< "$contract"

  command -v condor_submit >/dev/null 2>&1 || die "condor_submit is required after successful dry materialization"
  say "SUBMIT row=${row_id} descriptor=${submit_file} analysis_root=${analysis_root}"
  env -u RJ_DAG_DRYRUN -u RJ_THE134_MULTIVIEW_TRAINING_FILE \
    -u RJ_ID_FANOUT_FILE -u RJ_ID_FANOUT_DIRS_FILE \
    condor_submit "$submit_file" 2>&1 | tee "$row_log"
  grep -Eq '1 job\(s\) submitted to cluster [0-9]+' "$row_log" ||
    die "${row_id} did not submit exactly one Condor job"
  cluster="$(sed -n 's/.*1 job(s) submitted to cluster \([0-9][0-9]*\).*/\1/p' "$row_log" | tail -1)"
  [[ "$cluster" =~ ^[0-9]+$ ]] || die "could not recover the exact Condor cluster for ${row_id}"
  proc=0
  cluster_proc="${cluster}.${proc}"
  # This journal append is intentionally the first durable action after exact
  # cluster recovery.  A later descriptor/provenance failure cannot hide or
  # duplicate a successfully submitted row.
  printf '%s\t%s\t%s\t%s\n' "$row_id" "$cluster_proc" "$row_log" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$submission_journal"

  expected_args="$(resolve_condor_template "$(sed -n '1p' "$args_file")" "$cluster" "$proc")"
  submitted_args="$(condor_q "$cluster" -af ProcId Args 2>/dev/null || true)"
  [[ "$(wc -l <<< "$submitted_args" | tr -d ' ')" == 1 && "$submitted_args" == "0 ${expected_args}" ]] ||
    die "${row_id} submitted ClassAd does not contain exactly proc 0 with the frozen argument row"
  args_sha="$(printf '%s\n' "$submitted_args" | sha256_cmd | awk '{print $1}')"
  require_sha "${row_id} submitted arguments" "$args_sha"

  # Re-run the same evidence check after submission.  A drift here is a
  # preserved hard failure; the durable journal prevents duplicate recovery.
  verify_materialization_attempt_seal "$row_submit"
  [[ "$(verify_materialized_row_contract "$row_id" "$system" "$dataset" "$sample" "$row_output" "$sidecar" "$row_submit")" == "$contract" ]] ||
    die "${row_id} descriptor/args/fanout ownership evidence drifted after submission"

  log_template="$(condor_field "$submit_file" log)"
  out_template="$(condor_field "$submit_file" output)"
  err_template="$(condor_field "$submit_file" error)"
  [[ -n "$log_template" && -n "$out_template" && -n "$err_template" ]] ||
    die "${row_id} submit descriptor lacks exact log/stdout/stderr templates"
  log_path="$(resolve_condor_template "$log_template" "$cluster" "$proc")"
  out_path="$(resolve_condor_template "$out_template" "$cluster" "$proc")"
  err_path="$(resolve_condor_template "$err_template" "$cluster" "$proc")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$row_id" "$cluster_proc" "$submit_file" "$args_file" "$chunk_list" "$chunk_sha" \
    "$fanout_file" "$fanout_sha" "$log_path" "$out_path" "$err_path" \
    "$analysis_root" "$sidecar" "$owner_count" "$args_sha" \
    "$snapshot_dir" "$snapshot_header" "$snapshot_header_sha" \
    "$snapshot_calo" "$snapshot_calo_sha" "$snapshot_analysis" "$snapshot_analysis_sha" \
    "$materialized_config" "$materialized_config_sha" \
    "$snapshot_loader_receipt" "$snapshot_loader_receipt_sha" \
    "$snapshot_manifest" "$snapshot_manifest_sha" \
    >> "$submission_receipt"
  write_source_provenance
}

assert_fresh_submission() {
  command -v condor_q >/dev/null 2>&1 || die "condor_q is required for duplicate preflight"
  [[ ! -e "$output_root" ]] || die "output namespace already exists: ${output_root}"
  [[ ! -e "$submit_root" ]] || die "submit namespace already exists: ${submit_root}"
  [[ ! -e "$submission_journal" ]] || die "submission journal already exists: ${submission_journal}"
  [[ ! -e "$submission_receipt" ]] || die "submission receipt already exists: ${submission_receipt}"
  [[ ! -e "$source_provenance_json" ]] || die "source provenance already exists: ${source_provenance_json}"
  if condor_q "${USER:-patsfan753}" -af Args 2>/dev/null | grep -F "$output_root" >/dev/null; then
    die "an active owned job already targets ${output_root}"
  fi
}

table_has_row() {
  local table="$1" row_id="$2"
  [[ -s "$table" ]] && awk -F'\t' -v row="$row_id" '$1==row {found=1} END {exit !found}' "$table"
}

verify_resume_contract() {
  local hashes code_sha replay_schema_sha training_schema_sha semantic_sha manifest_hash
  hashes="$(require_inputs_and_hashes)"
  IFS=$'\t' read -r code_sha replay_schema_sha training_schema_sha semantic_sha \
    release_core_lib_dir release_core_lib64_dir release_calo_io release_clusteriso release_jetbase <<< "$hashes"
  write_and_verify_source_hashes "$RJ_THE134_SOURCE_HASH_MANIFEST_SHA256"
  manifest_hash="$(sha_file "$submission_manifest")"
  [[ "$manifest_hash" == "$(cat "$duplicate_fingerprint")" ]] ||
    die "submission manifest fingerprint drift blocks resume"
  [[ "$(manifest_field pp_signal_photon5 18)" == "$code_sha" ]] || die "resume code hash differs from frozen manifest"
  [[ "$(manifest_field pp_signal_photon5 19)" == "$replay_schema_sha" ]] || die "resume replay-schema hash differs from frozen manifest"
  [[ "$(manifest_field pp_signal_photon5 20)" == "$training_schema_sha" ]] || die "resume training-schema hash differs from frozen manifest"
  [[ "$(manifest_field pp_signal_photon5 21)" == "$semantic_sha" ]] || die "resume semantic hash differs from frozen manifest"
  [[ "$(manifest_field pp_signal_photon5 15)" == "$RJ_THE134_PP_LIBRARY_SHA256" ]] || die "resume p+p library differs from frozen manifest"
  [[ "$(manifest_field auau_signal_photon12 15)" == "$RJ_THE134_AUAU_LIBRARY_SHA256" ]] || die "resume Au+Au library differs from frozen manifest"
  [[ "$(manifest_field pp_signal_photon5 17)" == "$RJ_THE134_PP_MODEL_SHA256" ]] || die "resume p+p model differs from frozen manifest"
  [[ "$(manifest_field auau_signal_photon12 17)" == "$RJ_THE134_AUAU_MODEL_SHA256" ]] || die "resume Au+Au model differs from frozen manifest"
  [[ "$(manifest_field pp_signal_photon5 35)" == "$photon_cluster_builder_header" &&
     "$(manifest_field pp_signal_photon5 36)" == "$RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER_SHA256" ]] ||
    die "resume PhotonClusterBuilder header authority differs from frozen manifest"
  [[ "$(manifest_field pp_signal_photon5 37)" == "$calo_reco_library" &&
     "$(manifest_field pp_signal_photon5 38)" == "$RJ_THE134_CALO_RECO_LIBRARY_SHA256" ]] ||
    die "resume CaloReco runtime authority differs from frozen manifest"
  [[ "$(manifest_field pp_signal_photon5 39)" == "$release_core_lib_dir" &&
     "$(manifest_field pp_signal_photon5 40)" == "$release_core_lib64_dir" &&
     "$(manifest_field pp_signal_photon5 41)" == "$release_calo_io" &&
     "$(manifest_field pp_signal_photon5 42)" == "$RJ_THE134_RELEASE_CALO_IO_SHA256" &&
     "$(manifest_field pp_signal_photon5 43)" == "$release_clusteriso" &&
     "$(manifest_field pp_signal_photon5 44)" == "$RJ_THE134_RELEASE_CLUSTERISO_SHA256" &&
     "$(manifest_field pp_signal_photon5 45)" == "$release_jetbase" &&
     "$(manifest_field pp_signal_photon5 46)" == "$RJ_THE134_RELEASE_JETBASE_SHA256" ]] ||
    die "resume ana.560 release-companion authority differs from frozen manifest"
  write_runtime_authority_manifest verify
}

submit_all() {
  preflight
  assert_fresh_submission
  mkdir -p "$evidence_root" "$submit_root"
  printf 'row_id\tcluster_proc\tsubmit_log\tsubmitted_at_utc\n' > "$submission_journal"
  printf 'row_id\tcluster_proc\tsubmit_file\targs_file\tstaged_chunk_list\tstaged_chunk_sha256\tfanout_contract_file\tfanout_contract_sha256\tcondor_log\tcondor_stdout\tcondor_stderr\tanalysis_output_root\tmultiview_sidecar\tsidecar_owner_count\tsubmitted_args_sha256\tsnapshot_dir\tsnapshot_builder_header\tsnapshot_builder_header_sha256\tsnapshot_calo_reco_library\tsnapshot_calo_reco_library_sha256\tsnapshot_analysis_library\tsnapshot_analysis_library_sha256\tmaterialized_config\tmaterialized_config_sha256\tsnapshot_loader_receipt\tsnapshot_loader_receipt_sha256\tsnapshot_manifest\tsnapshot_manifest_sha256\n' > "$submission_receipt"
  while IFS='|' read -r row_id system lane dataset sample role mb_gate row_match; do
    submit_row "$row_id" "$system" "$lane" "$dataset" "$sample" "$role" "$mb_gate" "$row_match"
  done < <(emit_execution_matrix)
  [[ "$(wc -l < "$submission_receipt" | tr -d ' ')" == "$((execution_row_count + 1))" ]] ||
    die "submission receipt is incomplete"
  awk -F'\t' 'NF != 28 {exit 1}' "$submission_receipt" ||
    die "submission receipt must contain exactly 28 tab-separated fields on every row"
  [[ "$(wc -l < "$submission_journal" | tr -d ' ')" == "$((execution_row_count + 1))" ]] ||
    die "submission journal is incomplete"
  python3 - "$source_provenance_json" "$execution_row_count" <<'PY'
import json
import sys
payload = json.load(open(sys.argv[1]))
assert payload["schema"] == "THE134_SOURCE_PROVENANCE_V1"
assert len(payload["inputs"]) == int(sys.argv[2])
PY
  say "SUBMISSION_PASS rows=${execution_row_count} group_size=${execution_group_size} capacity_mode=${capacity_mode} receipt=${submission_receipt}"
}

resume_submit() {
  local row_id system lane dataset sample role mb_gate row_match row_log
  [[ -s "$submission_manifest" && -s "$duplicate_fingerprint" ]] || die "frozen preflight manifest is required for resume"
  [[ -s "$submission_journal" && -s "$submission_receipt" ]] || die "durable journal and receipt are required for resume"
  verify_resume_contract
  while IFS='|' read -r row_id system lane dataset sample role mb_gate row_match; do
    if table_has_row "$submission_journal" "$row_id"; then
      table_has_row "$submission_receipt" "$row_id" ||
        die "${row_id} has a submitted cluster but an incomplete receipt; preserve evidence and recover manually without resubmission"
      say "RESUME_SKIP already_submitted=${row_id} cluster=$(awk -F'\t' -v row="$row_id" '$1==row {print $2}' "$submission_journal")"
      continue
    fi
    row_log="${evidence_root}/submit_${row_id}.log"
    if [[ -s "$row_log" ]] && grep -Eq 'job\(s\) submitted to cluster [0-9]+' "$row_log"; then
      die "${row_id} has an unjournaled successful submit log; manual cluster recovery is required before resume"
    fi
    submit_row "$row_id" "$system" "$lane" "$dataset" "$sample" "$role" "$mb_gate" "$row_match"
  done < <(emit_execution_matrix)
  [[ "$(wc -l < "$submission_receipt" | tr -d ' ')" == "$((execution_row_count + 1))" ]] ||
    die "resumed submission receipt is incomplete"
  awk -F'\t' 'NF != 28 {exit 1}' "$submission_receipt" ||
    die "resumed submission receipt must contain exactly 28 tab-separated fields on every row"
  [[ "$(wc -l < "$submission_journal" | tr -d ' ')" == "$((execution_row_count + 1))" ]] ||
    die "resumed submission journal is incomplete"
  say "RESUME_SUBMISSION_PASS rows=${execution_row_count} group_size=${execution_group_size} receipt=${submission_receipt}"
}

status() {
  [[ -s "$submission_manifest" ]] || die "preflight manifest is missing: ${submission_manifest}"
  say "manifest=${submission_manifest} fingerprint=$(cat "$duplicate_fingerprint")"
  if [[ ! -s "$submission_journal" ]]; then
    say "no submission journal; no owned Condor row is declared"
    return 0
  fi
  while IFS=$'\t' read -r row_id cluster_proc _submit_log _submitted_at; do
    [[ "$row_id" == row_id ]] && continue
    local state
    state="$(condor_q "$cluster_proc" -af ClusterId ProcId JobStatus HoldReason NumJobStarts ExitCode 2>/dev/null || true)"
    if [[ -n "$state" ]]; then
      printf '%s\tQUEUE\t%s\n' "$row_id" "$state"
    else
      state="$(condor_history "$cluster_proc" -limit 1 -af ClusterId ProcId JobStatus HoldReason NumJobStarts ExitCode 2>/dev/null || true)"
      printf '%s\tHISTORY\t%s\n' "$row_id" "${state:-NOT_FOUND}"
    fi
  done < "$submission_journal"
  find "$output_root" -type f \( -name '*.root' -o -name '*.log' -o -name '*.json' \) -printf '%s\t%p\n' 2>/dev/null | sort -k2 || true
}

validate_root_health_and_joins() {
  command -v python3 >/dev/null 2>&1 || die "python3 is required for ROOT health and identity validation"
  python3 - "$submission_receipt" "$submission_manifest" \
    "$root_health_join_certificate" "$pinned_calo_reco_soname" \
    "$execution_row_count" "$execution_group_size" "$capacity_mode" <<'PY'
import csv
import hashlib
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    import ROOT
except Exception as exc:  # pragma: no cover - remote environment contract
    raise SystemExit(f"PyROOT is required for non-zombie/non-recovered validation: {exc}")

receipt_path, manifest_path, certificate_path = map(Path, sys.argv[1:4])
expected_calo_soname = sys.argv[4]
expected_receipt_count = int(sys.argv[5])
expected_source_count = int(sys.argv[6])
capacity_mode = bool(int(sys.argv[7]))
with manifest_path.open(newline="") as stream:
    manifests = {row["row_id"]: row for row in csv.DictReader(stream, delimiter="\t")}
with receipt_path.open(newline="") as stream:
    receipts = list(csv.DictReader(stream, delimiter="\t"))

EXPECTED_REPLAY_TREES = {
    "RJSourceOccurrenceV1",
    "RJEventV1",
    "RJPhotonCandidateV1",
    "RJModelEvaluationV1",
    "RJShowerCellV1",
    "RJShowerFeatureViewV1",
    "RJIsolationConstituentV1",
    "RJIsolationWitnessV1",
    "RJJetV1",
    "RJJetConstituentV1",
    "RJPhotonJetPairV1",
    "RJTruthPhotonV1",
    "RJTruthJetV1",
    "RJRecoTruthLinkV1",
    "RJWeightComponentV1",
    "RJEventDisplaySnapshotV1",
}
EXPECTED_DEFINITIONS = {"H70", "H0", "G70", "G0", "O70", "O0", "R70"}
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def identity(row, prefix: str) -> tuple[int, int]:
    return int(getattr(row, f"{prefix}_hi")), int(getattr(row, f"{prefix}_lo"))


def named_title(directory, name: str) -> str:
    obj = directory.Get(name)
    if not obj or not obj.InheritsFrom("TNamed"):
        raise ValueError(f"missing TNamed metadata:{name}")
    return str(obj.GetTitle())


def open_healthy(path: Path):
    if not path.is_file():
        raise ValueError(f"missing ROOT:{path}")
    size = path.stat().st_size
    if size < 50_000:
        raise ValueError(f"ROOT below 50000 bytes:{path}:{size}")
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or not root_file.IsOpen() or root_file.IsZombie():
        raise ValueError(f"unreadable or zombie ROOT:{path}")
    if root_file.TestBit(ROOT.TFile.kRecovered):
        root_file.Close()
        raise ValueError(f"recovered ROOT is forbidden:{path}")
    return root_file, size


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_frozen_snapshot_receipts(
    receipt: dict[str, str], manifest_row: dict[str, str], snapshot_dir: Path
) -> tuple[str, str]:
    loader_path = Path(receipt["snapshot_loader_receipt"]).resolve()
    snapshot_manifest_path = Path(receipt["snapshot_manifest"]).resolve()
    for label, path, receipt_sha in (
        ("snapshot loader receipt", loader_path, receipt["snapshot_loader_receipt_sha256"]),
        ("snapshot manifest", snapshot_manifest_path, receipt["snapshot_manifest_sha256"]),
    ):
        if not path.is_file() or snapshot_dir not in path.parents:
            raise ValueError(f"missing or out-of-snapshot {label}:{path}")
        if path.stat().st_mode & 0o222:
            raise ValueError(f"writable {label} survived snapshot seal:{path}")
        if not HEX64.fullmatch(receipt_sha) or file_sha256(path) != receipt_sha:
            raise ValueError(f"{label} hash drift:{path}")

    loader = json.loads(loader_path.read_text(encoding="utf-8"))
    if (
        loader.get("schema") != "RJ_SNAPSHOT_LOADER_RECEIPT_V1"
        or loader.get("status") != "PASS"
        or loader.get("pinned_calo_reco_release_companions") is not True
        or loader.get("mode") != manifest_row["system"]
    ):
        raise ValueError("snapshot loader receipt contract differs")
    if loader.get("release_roots") != [
        manifest_row["release_core_lib64_dir"],
        manifest_row["release_core_lib_dir"],
    ]:
        raise ValueError("snapshot loader receipt release roots differ")
    expected_providers = {
        "libcalo_reco.so": (
            (snapshot_dir / "lib/libcalo_reco.so").resolve(strict=True),
            manifest_row["calo_reco_library_sha256"],
        ),
        "libcalo_io.so": (
            Path(manifest_row["release_calo_io"]).resolve(strict=True),
            manifest_row["release_calo_io_sha256"],
        ),
        "libclusteriso.so": (
            Path(manifest_row["release_clusteriso"]).resolve(strict=True),
            manifest_row["release_clusteriso_sha256"],
        ),
        "libjetbase.so": (
            Path(manifest_row["release_jetbase"]).resolve(strict=True),
            manifest_row["release_jetbase_sha256"],
        ),
    }
    calo_api = snapshot_dir / "lib/libcalo_reco.so"
    calo_soname = snapshot_dir / "lib" / expected_calo_soname
    if (
        not calo_api.is_file()
        or not calo_soname.is_symlink()
        or calo_soname.readlink().as_posix() != "libcalo_reco.so"
        or calo_soname.resolve(strict=True) != calo_api.resolve(strict=True)
    ):
        raise ValueError("snapshot CaloReco API/SONAME alias closure differs")
    providers = loader.get("providers", {})
    for family, (path, expected_sha) in expected_providers.items():
        provider = providers.get(family, {})
        if (
            provider.get("realpath") != str(path)
            or provider.get("sha256") != expected_sha
            or provider.get("observed_resolutions") != [str(path)]
            or file_sha256(path) != expected_sha
        ):
            raise ValueError(f"snapshot loader provider authority differs:{family}")

    snapshot_manifest = json.loads(
        snapshot_manifest_path.read_text(encoding="utf-8")
    )
    if (
        snapshot_manifest.get("schema") != "RJ_FROZEN_SNAPSHOT_MANIFEST_V1"
        or snapshot_manifest.get("status") != "PASS"
        or snapshot_manifest.get("root") != str(snapshot_dir)
    ):
        raise ValueError("snapshot manifest contract differs")
    if snapshot_dir.stat().st_mode & 0o222:
        raise ValueError(f"writable frozen snapshot root survived seal:{snapshot_dir}")
    observed: list[dict[str, object]] = []
    for path in sorted(
        snapshot_dir.rglob("*"),
        key=lambda item: item.relative_to(snapshot_dir).as_posix(),
    ):
        if path == snapshot_manifest_path:
            continue
        relative = path.relative_to(snapshot_dir).as_posix()
        if path.is_symlink():
            target = os.readlink(path)
            if os.path.isabs(target):
                raise ValueError(f"snapshot symlink must be relative:{relative}")
            resolved = path.resolve(strict=True)
            if snapshot_dir not in resolved.parents:
                raise ValueError(f"snapshot symlink escaped:{relative}")
            if path.parent.stat().st_mode & 0o222:
                raise ValueError(f"snapshot symlink parent is writable:{relative}")
            if resolved.stat().st_mode & 0o222:
                raise ValueError(f"snapshot symlink target is writable:{relative}")
            observed.append(
                {
                    "path": relative,
                    "sha256": file_sha256(resolved),
                    "symlink_target": target,
                    "type": "symlink",
                }
            )
        elif path.is_file():
            if path.lstat().st_mode & 0o222:
                raise ValueError(f"writable artifact survived snapshot seal:{relative}")
            observed.append(
                {
                    "path": relative,
                    "sha256": file_sha256(path),
                    "size": path.stat().st_size,
                    "type": "file",
                }
            )
        elif path.is_dir():
            if path.lstat().st_mode & 0o222:
                raise ValueError(f"writable directory survived snapshot seal:{relative}")
            observed.append({"path": relative, "type": "directory"})
        else:
            raise ValueError(f"unsupported snapshot artifact:{relative}")
    if snapshot_manifest.get("entries") != observed:
        raise ValueError("complete snapshot inventory/hash closure differs")
    return receipt["snapshot_loader_receipt_sha256"], receipt["snapshot_manifest_sha256"]


failures: list[str] = []
reports: list[dict[str, object]] = []
seen_analysis: set[str] = set()
seen_sidecars: set[str] = set()

for receipt in receipts:
    row_id = receipt["row_id"]
    report: dict[str, object] = {"row_id": row_id}
    try:
        manifest = manifests[row_id]
        if manifest["full_training_authority"] != "0":
            raise ValueError("smoke row asserted full training authority")
        if receipt["sidecar_owner_count"] != "1":
            raise ValueError("sidecar owner multiplicity is not exactly one")
        materialized_config = Path(receipt["materialized_config"]).resolve()
        materialized_config_sha = receipt["materialized_config_sha256"]
        if not materialized_config.is_file():
            raise ValueError(f"missing descriptor-bound materialized config:{materialized_config}")
        if not HEX64.fullmatch(materialized_config_sha):
            raise ValueError("materialized config receipt SHA-256 is malformed")
        if file_sha256(materialized_config) != materialized_config_sha:
            raise ValueError(f"materialized config hash drift:{materialized_config}")
        snapshot_dir = Path(receipt["snapshot_dir"]).resolve()
        snapshot_authorities = (
            (
                "builder header",
                Path(receipt["snapshot_builder_header"]).resolve(),
                receipt["snapshot_builder_header_sha256"],
                manifest["photon_cluster_builder_header_sha256"],
            ),
            (
                "CaloReco library",
                Path(receipt["snapshot_calo_reco_library"]).resolve(),
                receipt["snapshot_calo_reco_library_sha256"],
                manifest["calo_reco_library_sha256"],
            ),
            (
                "analysis library",
                Path(receipt["snapshot_analysis_library"]).resolve(),
                receipt["snapshot_analysis_library_sha256"],
                manifest["library_sha256"],
            ),
        )
        for label, path, receipt_sha, expected_sha in snapshot_authorities:
            if not path.is_file() or snapshot_dir not in path.parents:
                raise ValueError(f"missing or out-of-snapshot {label}:{path}")
            if receipt_sha != expected_sha or file_sha256(path) != expected_sha:
                raise ValueError(f"snapshotted {label} hash drift:{path}")
        if (snapshot_dir / "lib/libphoton_cluster_builder_override.so").exists():
            raise ValueError("undeclared PhotonClusterBuilder override library is present")
        loader_receipt_sha, snapshot_manifest_sha = verify_frozen_snapshot_receipts(
            receipt, manifest, snapshot_dir
        )
        analysis_text = receipt["analysis_output_root"]
        sidecar_text = receipt["multiview_sidecar"]
        if analysis_text in seen_analysis or sidecar_text in seen_sidecars:
            raise ValueError("duplicate analysis ROOT or sidecar ownership")
        seen_analysis.add(analysis_text)
        seen_sidecars.add(sidecar_text)
        analysis_path = Path(analysis_text)
        sidecar_path = Path(sidecar_text)

        analysis, analysis_size = open_healthy(analysis_path)
        replay = analysis.GetDirectory("ReplayFoundationV1")
        if not replay:
            raise ValueError("missing ReplayFoundationV1 directory")
        observed_trees = {
            str(key.GetName())
            for key in replay.GetListOfKeys()
            if str(key.GetClassName()) == "TTree"
        }
        if observed_trees != EXPECTED_REPLAY_TREES:
            raise ValueError(
                f"ReplayFoundation tree inventory drift missing={sorted(EXPECTED_REPLAY_TREES-observed_trees)} "
                f"extra={sorted(observed_trees-EXPECTED_REPLAY_TREES)}"
            )
        replay_expected = {
            "rj_replay_schema": "RJ_REPLAY_FOUNDATION_V1",
            "rj_replay_schema_version": "2",
            "rj_replay_complete": "1",
            "schema_sha256": manifest["replay_schema_sha256"],
            "semantic_sha256": manifest["semantic_sha256"],
            "source_sha256": manifest["source_manifest_sha256"],
            "model_sha256": manifest["model_sha256"],
            "config_sha256": manifest["resolved_config_sha256"],
            "code_sha256": manifest["code_sha256"],
        }
        for key, expected in replay_expected.items():
            observed = named_title(replay, key)
            if observed != expected:
                raise ValueError(f"ReplayFoundation metadata drift:{key}:{observed}!={expected}")

        source_tree = replay.Get("RJSourceOccurrenceV1")
        event_tree = replay.Get("RJEventV1")
        candidate_tree = replay.Get("RJPhotonCandidateV1")
        if int(source_tree.GetEntries()) != expected_source_count:
            raise ValueError(
                "replay source-row population differs from the execution group: "
                f"expected={expected_source_count} observed={int(source_tree.GetEntries())}"
            )
        source_ids = {identity(row, "source_occurrence_id") for row in source_tree}
        if len(source_ids) != expected_source_count:
            raise ValueError(
                "replay source identities differ from the execution group: "
                f"expected={expected_source_count} observed={len(source_ids)}"
            )
        event_to_source: dict[tuple[int, int], tuple[int, int]] = {}
        for row in event_tree:
            event_id = identity(row, "event_id")
            source_id = identity(row, "source_occurrence_id")
            if event_id in event_to_source or source_id not in source_ids:
                raise ValueError("duplicate event identity or orphan source foreign key")
            event_to_source[event_id] = source_id
        candidate_to_event: dict[tuple[int, int], tuple[int, int]] = {}
        for row in candidate_tree:
            candidate_id = identity(row, "candidate_id")
            event_id = identity(row, "event_id")
            if candidate_id in candidate_to_event or event_id not in event_to_source:
                raise ValueError("duplicate candidate identity or orphan event foreign key")
            candidate_to_event[candidate_id] = event_id

        sidecar, sidecar_size = open_healthy(sidecar_path)
        sidecar_trees = {
            str(key.GetName())
            for key in sidecar.GetListOfKeys()
            if str(key.GetClassName()) == "TTree"
        }
        if sidecar_trees != {"RJPhotonTrainingViewV1"}:
            raise ValueError(f"sidecar tree inventory drift:{sorted(sidecar_trees)}")
        sidecar_expected = {
            "rj_photon_training_schema": "RJ_PHOTON_TRAINING_VIEW_V1",
            "rj_photon_training_schema_version": "1",
            "rj_photon_training_complete": "1",
            "schema_sha256": manifest["training_schema_sha256"],
            "source_manifest_sha256": manifest["source_manifest_sha256"],
            "config_sha256": manifest["resolved_config_sha256"],
            "code_sha256": manifest["code_sha256"],
        }
        for key, expected in sidecar_expected.items():
            observed = named_title(sidecar, key)
            if observed != expected:
                raise ValueError(f"sidecar metadata drift:{key}:{observed}!={expected}")
        for feature_key in ("pp_feature_contract_sha256", "auau_feature_contract_sha256"):
            if not HEX64.fullmatch(named_title(sidecar, feature_key)):
                raise ValueError(f"invalid sidecar feature contract hash:{feature_key}")

        training_tree = sidecar.Get("RJPhotonTrainingViewV1")
        entry_count = int(training_tree.GetEntries())
        if named_title(sidecar, "rj_photon_training_entries") != str(entry_count):
            raise ValueError("sidecar completion entry count does not match its tree")
        candidate_definitions: dict[tuple[int, int], set[str]] = defaultdict(set)
        joined_rows = 0
        for row in training_tree:
            source_id = identity(row, "source_occurrence_id")
            event_id = identity(row, "event_id")
            candidate_id = identity(row, "candidate_id")
            definition = str(row.definition_name)
            if source_id not in source_ids:
                raise ValueError("sidecar source identity does not join ReplayFoundation")
            if event_to_source.get(event_id) != source_id:
                raise ValueError("sidecar event/source identity join failed")
            if candidate_to_event.get(candidate_id) != event_id:
                raise ValueError("sidecar candidate/event identity join failed")
            if definition in candidate_definitions[candidate_id]:
                raise ValueError("duplicate sidecar candidate/definition identity")
            candidate_definitions[candidate_id].add(definition)
            joined_rows += 1
        for candidate_id, definitions in candidate_definitions.items():
            if definitions != EXPECTED_DEFINITIONS:
                raise ValueError(
                    f"candidate {candidate_id} does not carry exactly seven frozen views:{sorted(definitions)}"
                )
        if joined_rows != 7 * len(candidate_definitions):
            raise ValueError("sidecar seven-view population closure failed")
        if capacity_mode and not candidate_definitions:
            raise ValueError("capacity witness contains no retained photon candidates")

        report.update(
            {
                "status": "PASS",
                "analysis_output_root": analysis_text,
                "analysis_bytes": analysis_size,
                "sidecar": sidecar_text,
                "sidecar_bytes": sidecar_size,
                "replay_tree_count": len(observed_trees),
                "replay_sources": len(source_ids),
                "replay_events": len(event_to_source),
                "replay_candidates": len(candidate_to_event),
                "sidecar_entries": entry_count,
                "sidecar_candidates": len(candidate_definitions),
                "joined_sidecar_rows": joined_rows,
                "snapshot_loader_receipt_sha256": loader_receipt_sha,
                "snapshot_manifest_sha256": snapshot_manifest_sha,
                "full_training_authority": 0,
            }
        )
        sidecar.Close()
        analysis.Close()
    except Exception as exc:
        failures.append(f"{row_id}:{exc}")
        report.update({"status": "FAIL", "failure": str(exc)})
    reports.append(report)

if len(receipts) != expected_receipt_count or len(manifests) != 13:
    failures.append(
        "row_closure:"
        f"receipts={len(receipts)}/{expected_receipt_count} manifests={len(manifests)}/13"
    )
if capacity_mode and {row.get("row_id") for row in receipts} != {
    "pp_background_jet8",
    "auau_background_jet12",
}:
    failures.append(
        "capacity_selected_rows:"
        f"{sorted(str(row.get('row_id')) for row in receipts)}"
    )

payload = {
    "schema": "THE134_SMOKE_ROOT_HEALTH_IDENTITY_JOIN_V1",
    "status": "PASS" if not failures else "FAIL",
    "scope": "capacity" if capacity_mode else "smoke",
    "full_training_authority": 0,
    "execution_group_size": expected_source_count,
    "capacity_authority_earned": bool(capacity_mode and not failures),
    "row_count": len(reports),
    "rows": reports,
    "failures": failures,
}
temporary = certificate_path.with_name(certificate_path.name + f".tmp.{os.getpid()}")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
os.replace(temporary, certificate_path)
print(json.dumps({"status": payload["status"], "rows": len(reports), "failures": failures}))
raise SystemExit(0 if not failures else 1)
PY
}

write_capacity_resource_certificate() {
  (( capacity_mode )) || return 0
  python3 - \
    "$submission_journal" "$root_health_join_certificate" \
    "$capacity_resource_certificate" "$capacity_preflight_receipt" \
    "$capacity_preflight_receipt_sha" "$capacity_full_plan" \
    "$capacity_full_plan_sha" "$capacity_canary_id" <<'PY'
from pathlib import Path
import csv
import hashlib
import json
import os
import subprocess
import sys

(
    journal_arg,
    root_certificate_arg,
    destination_arg,
    preflight_receipt_arg,
    preflight_receipt_sha,
    full_plan_arg,
    full_plan_sha,
    capacity_id,
) = sys.argv[1:]
journal_path = Path(journal_arg)
root_certificate_path = Path(root_certificate_arg)
destination = Path(destination_arg)
preflight_receipt_path = Path(preflight_receipt_arg)
full_plan_path = Path(full_plan_arg)

def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

if digest(preflight_receipt_path) != preflight_receipt_sha:
    raise SystemExit("capacity preflight receipt drifted before resource certification")
if digest(full_plan_path) != full_plan_sha:
    raise SystemExit("capacity full plan drifted before resource certification")
preflight = json.loads(preflight_receipt_path.read_text(encoding="utf-8"))
root_certificate = json.loads(root_certificate_path.read_text(encoding="utf-8"))
if (
    root_certificate.get("status") != "PASS"
    or root_certificate.get("scope") != "capacity"
    or root_certificate.get("capacity_authority_earned") is not True
    or int(root_certificate.get("execution_group_size", -1)) != 7
    or int(root_certificate.get("row_count", -1)) != 2
):
    raise SystemExit("ROOT/identity capacity certificate is not PASS")

with journal_path.open(newline="") as stream:
    journal_rows = list(csv.DictReader(stream, delimiter="\t"))
if len(journal_rows) != 2:
    raise SystemExit(f"capacity journal must contain exactly two rows: {len(journal_rows)}")

reports = []
for row in journal_rows:
    cluster_proc = row["cluster_proc"]
    completed = subprocess.run(
        ["condor_history", cluster_proc, "-limit", "1", "-json"],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise SystemExit(
            f"condor_history failed for {cluster_proc}: {completed.stderr.strip()}"
        )
    ads = json.loads(completed.stdout or "[]")
    if len(ads) != 1:
        raise SystemExit(
            f"capacity history must contain exactly one ClassAd for {cluster_proc}"
        )
    ad = ads[0]
    expected_cluster, expected_proc = map(int, cluster_proc.split(".", 1))
    if int(ad.get("ClusterId", -1)) != expected_cluster or int(ad.get("ProcId", -1)) != expected_proc:
        raise SystemExit(f"capacity history identity mismatch: {cluster_proc}")
    job_status = int(ad.get("JobStatus", -1))
    exit_code = int(ad.get("ExitCode", -1))
    starts = int(ad.get("NumJobStarts", -1))
    holds = int(ad.get("NumHolds", 0))
    request_memory_mb = int(ad.get("RequestMemory", -1))
    memory_usage_mb = int(ad.get("MemoryUsage", -1))
    resident_set_kb = int(ad.get("ResidentSetSize_RAW", ad.get("ResidentSetSize", -1)))
    wall_seconds = float(ad.get("RemoteWallClockTime", -1.0))
    if job_status != 4 or exit_code != 0:
        raise SystemExit(
            f"capacity terminal state failed: {cluster_proc} status={job_status} exit={exit_code}"
        )
    if starts != 1 or holds != 0:
        raise SystemExit(
            f"capacity execution was retried or held: {cluster_proc} starts={starts} holds={holds}"
        )
    if request_memory_mb != 8000:
        raise SystemExit(
            f"capacity RequestMemory drift: {cluster_proc} request={request_memory_mb}"
        )
    if memory_usage_mb <= 0 or memory_usage_mb > request_memory_mb:
        raise SystemExit(
            f"capacity memory usage exceeds its frozen request: {cluster_proc} "
            f"usage={memory_usage_mb} request={request_memory_mb}"
        )
    if resident_set_kb <= 0 or wall_seconds <= 0:
        raise SystemExit(
            f"capacity runtime metrics are incomplete: {cluster_proc} "
            f"rss_kb={resident_set_kb} wall={wall_seconds}"
        )
    reports.append(
        {
            "cluster_proc": cluster_proc,
            "exit_code": exit_code,
            "job_status": job_status,
            "memory_usage_mb": memory_usage_mb,
            "num_holds": holds,
            "num_job_starts": starts,
            "remote_wall_clock_seconds": wall_seconds,
            "request_memory_mb": request_memory_mb,
            "resident_set_size_kb": resident_set_kb,
            "row_id": row["row_id"],
        }
    )

payload = {
    "bundle_manifest_sha256": preflight["bundle_manifest_sha256"],
    "capacity_authority_earned": True,
    "capacity_canary_id": capacity_id,
    "execution_group_size": 7,
    "execution_partition_sha256": preflight["execution_partition_sha256"],
    "full_plan": str(full_plan_path.resolve()),
    "full_plan_sha256": full_plan_sha,
    "full_training_authority": 0,
    "materialization_receipt_sha256": preflight["materialization_receipt_sha256"],
    "preflight_receipt": str(preflight_receipt_path.resolve()),
    "preflight_receipt_sha256": preflight_receipt_sha,
    "root_health_identity_join_certificate": str(root_certificate_path.resolve()),
    "root_health_identity_join_certificate_sha256": digest(root_certificate_path),
    "rows": reports,
    "schema": "THE134_GROUP7_PARTITION_CAPACITY_CERTIFICATE_V1",
    "selected_rows": ["pp_background_jet8", "auau_background_jet12"],
    "status": "PASS",
    "submission_performed": True,
}
temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
os.replace(temporary, destination)
os.chmod(destination, 0o444)
print(
    json.dumps(
        {
            "status": "PASS",
            "capacity_authority_earned": True,
            "rows": len(reports),
            "certificate": str(destination),
            "certificate_sha256": digest(destination),
        },
        sort_keys=True,
    )
)
PY
}

validate_outputs() {
  local row_id cluster_proc queue_state history_state job_status exit_code
  [[ -s "$validator" ]] || die "THE-134 matrix preparer is missing: ${validator}"
  [[ -s "$submission_journal" &&
     "$(wc -l < "$submission_journal" | tr -d ' ')" == "$((execution_row_count + 1))" ]] ||
    die "complete ${execution_row_count}-row submission journal is required"
  [[ -s "$submission_receipt" &&
     "$(wc -l < "$submission_receipt" | tr -d ' ')" == "$((execution_row_count + 1))" ]] ||
    die "complete ${execution_row_count}-row submission receipt is required"
  [[ -s "$pp_source_provenance_json" && -s "$auau_source_provenance_json" ]] ||
    die "system-specific source provenance is incomplete"

  while IFS=$'\t' read -r row_id cluster_proc _submit_log _submitted_at; do
    [[ "$row_id" == row_id ]] && continue
    queue_state="$(condor_q "$cluster_proc" -af JobStatus HoldReason NumJobStarts 2>/dev/null || true)"
    [[ -z "$queue_state" ]] || die "${row_id} is not terminal: ${cluster_proc} ${queue_state}"
    history_state="$(condor_history "$cluster_proc" -limit 1 -af JobStatus ExitCode 2>/dev/null || true)"
    read -r job_status exit_code <<< "$history_state"
    [[ "$job_status" == 4 && "$exit_code" == 0 ]] ||
      die "${row_id} terminal gate failed: cluster=${cluster_proc} history='${history_state:-NOT_FOUND}'"
  done < "$submission_journal"

  validate_root_health_and_joins

  python3 "$validator" \
    --system pp --view H70 --scope smoke \
    --input "@${evidence_root}/pp_multiview_sidecars.list" \
    --source-provenance-json "$pp_source_provenance_json" \
    --matrix-out "${evidence_root}/pp_h70_source_complete_smoke_matrix.npz" \
    --audit-out "${evidence_root}/pp_h70_source_complete_smoke_audit.json"
  python3 "$validator" \
    --system auau --view H70 --scope smoke \
    --input "@${evidence_root}/auau_multiview_sidecars.list" \
    --source-provenance-json "$auau_source_provenance_json" \
    --matrix-out "${evidence_root}/auau_h70_source_complete_smoke_matrix.npz" \
    --audit-out "${evidence_root}/auau_h70_source_complete_smoke_audit.json"

  python3 - \
    "${evidence_root}/pp_h70_source_complete_smoke_audit.json" \
    "${evidence_root}/auau_h70_source_complete_smoke_audit.json" <<'PY'
import json
import sys
for path in sys.argv[1:]:
    payload = json.load(open(path))
    if payload.get("status") != "PASS":
        raise SystemExit(f"smoke matrix audit is not PASS: {path}")
    if payload.get("scope") != "smoke" or int(payload.get("full_training_authority", -1)) != 0:
        raise SystemExit(f"smoke audit attempted to assert full training authority: {path}")
PY

  write_capacity_resource_certificate
  if (( capacity_mode )); then
    say "CAPACITY_VALIDATION_PASS selected_rows=${execution_row_count} group_size=${execution_group_size} full_training_authority=0 capacity_certificate=${capacity_resource_certificate} root_join_certificate=${root_health_join_certificate}"
  else
    say "SMOKE_VALIDATION_PASS source_categories=13 full_training_authority=0 root_join_certificate=${root_health_join_certificate} pp_audit=${evidence_root}/pp_h70_source_complete_smoke_audit.json auau_audit=${evidence_root}/auau_h70_source_complete_smoke_audit.json"
  fi
}

case "$mode" in
  inventory) emit_source_hash_manifest ;;
  preflight) preflight ;;
  submit) submit_all ;;
  resume-submit) resume_submit ;;
  status) status ;;
  validate) validate_outputs ;;
  *) die "usage: $0 inventory|preflight|submit|resume-submit|status|validate|capacity-preflight|capacity-submit|capacity-resume-submit|capacity-status|capacity-validate" ;;
esac
