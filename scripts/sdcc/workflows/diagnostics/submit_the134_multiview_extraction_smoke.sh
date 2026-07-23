#!/usr/bin/env bash
set -euo pipefail

# Source-complete, writer-only THE-134 multi-view extraction smoke.
#
# This is a narrow controller over the established RecoilJets Condor executor.
# It owns exactly one input tuple and one Condor proc for each frozen training
# source.  It never submits data, p+p Jet40, a direct arm, a merge, or a model
# training job.  `inventory` is read-only and emits the exact frozen-source
# manifest to stdout; submission remains an explicit `submit` action.

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd -P)"
cd "$repo_root"

mode="${1:-preflight}"
tag="${RJ_THE134_EXTRACTION_TAG:-the134_h70_multiview_extraction_smoke_20260722_v1}"
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
pp_config="${RJ_THE134_PP_CONFIG:-${repo_root}/macros/analysis_config_the119_pp_replay_foundation.yaml}"
auau_config="${RJ_THE134_AUAU_CONFIG:-${repo_root}/macros/analysis_config_the112_auau_combined_bdt_triplet.yaml}"
pp_library="${RJ_THE134_PP_LIBRARY:-}"
auau_library="${RJ_THE134_AUAU_LIBRARY:-}"
pp_model="${RJ_THE134_PP_MODEL:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the116_models/the116_pp_matched_basev3e_15to35_20260720/models/bdt_ppg12_basev3e_15to35/pp_tight_bdt_ppg12_base_v3E_bdt_noIso_tmva.root}"
auau_model="${RJ_THE134_AUAU_MODEL:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the111_models/the111_combined_corrected_shower_ppg12_labels_20260719_1618/combined/auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.root}"
source_hash_manifest="${RJ_THE134_SOURCE_HASH_MANIFEST:-}"

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
validator="${RJ_THE134_MATRIX_PREPARER:-${repo_root}/scripts/ml/training/prepare_the134_h70_matrix.py}"

readonly nominal_et_min="15.0"
readonly nominal_et_max="35.0"
readonly capture_et_min="5.0"
readonly legacy_tree_max_entries="0"
readonly event_limit_per_job="0"
readonly full_training_authority="0"
readonly extraction_cone_r="0.40"
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
  local sample="$1" sample_root="${sim_root}/${sample}"
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
  tuple_sha="$(printf '%s\n' "$first_tuple" | sha256_cmd | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "$manifest_sha" "$tuple_sha" "$executable_count"
}

emit_source_hash_manifest() {
  local row_id sample hashes
  validate_matrix
  printf 'row_id\tsource_manifest_sha256\tfirst_input_tuple_sha256\texecutable_tuple_count\n'
  while IFS='|' read -r row_id _system _lane _dataset sample _role _mb_gate _row_match; do
    hashes="$(source_hashes_for_sample "$sample")"
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
  printf 'row_id\tsystem\tlane\tdataset\tsample\tsource_role\tminimum_bias_gate\tinput_files\tinput_jobs\tsource_manifest_sha256\tfirst_input_tuple_sha256\tresolved_config\tresolved_config_sha256\tlibrary\tlibrary_sha256\tmodel\tmodel_sha256\tcode_sha256\treplay_schema_sha256\ttraining_schema_sha256\tsemantic_sha256\tnominal_et_min_gev\tnominal_et_max_gev_exclusive\tloose_capture_et_min_gev\tlegacy_training_tree_max_entries\tevent_limit_per_job\tanalysis_output_namespace\tmultiview_sidecar\tsubmit_namespace\tscheduler_log_dir\tscheduler_stdout_dir\tscheduler_stderr_dir\texecutable_input_tuple_count\tfull_training_authority\tphoton_cluster_builder_header\tphoton_cluster_builder_header_sha256\tcalo_reco_library\tcalo_reco_library_sha256\n' > "$tmp"
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
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t1\t1\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$row_id" "$system" "$lane" "$dataset" "$sample" "$role" "$mb_gate" \
      "$source_sha" "$tuple_sha" "$config" "$config_sha" "$library" "$library_sha" \
      "$model" "$model_sha" "$code_sha" "$replay_schema_sha" "$training_schema_sha" "$semantic_sha" \
      "$nominal_et_min" "$nominal_et_max" "$capture_et_min" "$legacy_tree_max_entries" "$event_limit_per_job" \
      "$row_output" "$sidecar" "$row_submit" \
      '/sphenix/u/patsfan753/scratch/thesisAnalysis/log' \
      '/sphenix/u/patsfan753/scratch/thesisAnalysis/stdout' \
      '/sphenix/u/patsfan753/scratch/thesisAnalysis/error' \
      "$tuple_count" "$full_training_authority" \
      "$photon_cluster_builder_header" "$RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER_SHA256" \
      "$calo_reco_library" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256" >> "$tmp"
  done < <(emit_matrix)
  mv "$tmp" "$submission_manifest"
  [[ "$(wc -l < "$submission_manifest" | tr -d ' ')" == 14 ]] || die "submission manifest row closure failed"
  [[ "$(tail -n +2 "$submission_manifest" | cut -f1 | sort -u | wc -l | tr -d ' ')" == 13 ]] || die "submission manifest contains duplicate row identities"
  [[ "$(tail -n +2 "$submission_manifest" | cut -f34 | sort -u)" == "$full_training_authority" ]] ||
    die "smoke manifest must remain full_training_authority=0"
  [[ "$(tail -n +2 "$submission_manifest" | cut -f28 | sort -u | wc -l | tr -d ' ')" == 13 ]] ||
    die "every source row must own a unique multiview sidecar path"
  sha_file "$submission_manifest" > "$duplicate_fingerprint"
  awk -F'\t' 'NR>1 && $2=="pp" {print $28}' "$submission_manifest" > "${evidence_root}/pp_multiview_sidecars.list"
  awk -F'\t' 'NR>1 && $2=="auau" {print $28}' "$submission_manifest" > "${evidence_root}/auau_multiview_sidecars.list"
  [[ "$(wc -l < "${evidence_root}/pp_multiview_sidecars.list" | tr -d ' ')" == 7 ]] || die "p+p sidecar manifest closure failed"
  [[ "$(wc -l < "${evidence_root}/auau_multiview_sidecars.list" | tr -d ' ')" == 6 ]] || die "Au+Au sidecar manifest closure failed"
}

require_inputs_and_hashes() {
  local actual_code actual_replay_schema actual_training_schema actual_semantic
  local pp_yaml_model auau_yaml_model
  validate_matrix
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
  [[ -z "${RJ_THE134_MULTIVIEW_TRAINING_FILE:-}" ]] ||
    die "controller environment must not pre-own RJ_THE134_MULTIVIEW_TRAINING_FILE"
  [[ -z "${RJ_ID_FANOUT_FILE:-}" && -z "${RJ_ID_FANOUT_DIRS_FILE:-}" ]] ||
    die "controller environment must not carry a pre-existing fanout owner"
  [[ -z "${RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE:-}" ]] ||
    die "a separate PhotonClusterBuilder override library is forbidden; the pinned CaloReco library is authoritative"

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

  require_file_hash "p+p library" "$pp_library" "$RJ_THE134_PP_LIBRARY_SHA256"
  require_file_hash "Au+Au library" "$auau_library" "$RJ_THE134_AUAU_LIBRARY_SHA256"
  require_file_hash "p+p base config" "$pp_config" "$RJ_THE134_PP_CONFIG_SHA256"
  require_file_hash "Au+Au base config" "$auau_config" "$RJ_THE134_AUAU_CONFIG_SHA256"
  require_file_hash "p+p model" "$pp_model" "$RJ_THE134_PP_MODEL_SHA256"
  require_file_hash "Au+Au model" "$auau_model" "$RJ_THE134_AUAU_MODEL_SHA256"
  require_file_hash "PhotonClusterv1 build header" "$photon_cluster_header" "$RJ_THE134_PHOTON_CLUSTER_HEADER_SHA256"
  require_file_hash "PhotonClusterBuilder build header" "$photon_cluster_builder_header" "$RJ_THE134_PHOTON_CLUSTER_BUILDER_HEADER_SHA256"
  require_file_hash "CaloReco runtime library" "$calo_reco_library" "$RJ_THE134_CALO_RECO_LIBRARY_SHA256"

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

  printf '%s\t%s\t%s\t%s\n' "$actual_code" "$actual_replay_schema" "$actual_training_schema" "$actual_semantic"
}

preflight() {
  local hashes code_sha replay_schema_sha training_schema_sha semantic_sha
  [[ ! -e "$submission_journal" && ! -e "$submission_receipt" ]] ||
    die "submission evidence already exists; use status or validate instead of rewriting preflight state"
  hashes="$(require_inputs_and_hashes)"
  IFS=$'\t' read -r code_sha replay_schema_sha training_schema_sha semantic_sha <<< "$hashes"
  mkdir -p "$evidence_root"
  prepare_resolved_configs
  verify_single_owner_resolved_contract
  write_and_verify_source_hashes "$RJ_THE134_SOURCE_HASH_MANIFEST_SHA256"
  write_submission_manifest "$code_sha" "$replay_schema_sha" "$training_schema_sha" "$semantic_sha"
  bash -n "$0" "$submitter" "$pp_executor" "$auau_executor"
  say "PREFLIGHT_PASS rows=13 manifest=${submission_manifest} fingerprint=$(cat "$duplicate_fingerprint")"
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
    period=run28
  else
    model_score=auau_tight_bdt_score
    shower_definition=H0
    si_di_role=EMBEDDED
    period=AUAU_RUN24
  fi
  printf '%s' \
    "RJ_REPLAY_FOUNDATION_V1=1;RJ_REPLAY_FOUNDATION_CANARY=1;RJ_REPLAY_TRACE=0;RJ_REPLAY_LANE=${lane};RJ_REPLAY_DATASET=${dataset};RJ_REPLAY_SAMPLE=${sample};RJ_REPLAY_PERIOD=${period};RJ_REPLAY_SI_DI_ROLE=${si_di_role};RJ_REPLAY_OWNERSHIP_STATE=source_role_frozen;RJ_REPLAY_SOURCE_MANIFEST_SHA256=${source_sha};RJ_REPLAY_SOURCE_SHA256=${source_sha};RJ_REPLAY_MODEL_SHA256=${model_sha};RJ_REPLAY_MODEL_SCORE_NAME=${model_score};RJ_REPLAY_MODEL_SHOWER_DEFINITION=${shower_definition};RJ_REPLAY_CONFIG_SHA256=${config_sha};RJ_REPLAY_CODE_SHA256=${RJ_THE134_CODE_SHA256};RJ_REPLAY_SCHEMA_SHA256=${RJ_THE134_REPLAY_SCHEMA_SHA256};RJ_REPLAY_SEMANTIC_SHA256=${RJ_THE134_SEMANTIC_SHA256};RJ_REPLAY_PHOTON_CAPTURE_ET_MIN=${capture_et_min};RJ_REPLAY_JET_CONSTITUENT_PT_MIN=${capture_et_min};RJ_THE134_MULTIVIEW_TRAINING_V1=1;RJ_THE134_MULTIVIEW_TRAINING_FILE=${sidecar};RJ_THE134_EXPECTED_SOURCE_ROLE=${role}"
}

system_extra_env() {
  local system="$1" role="$2"
  if [[ "$system" == pp ]]; then
    printf '%s' \
      "RJ_PP_PHOTONID_EXTRACT_ONLY=1;RJ_PP_PHOTONID_TRAINING_TREE=1;RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES=${legacy_tree_max_entries};RJ_PP_PHOTONID_SOURCE_ROLE=${role};RJ_PP_PHOTONID_PPG12_FILTER=1;RJ_PP_PHOTONID_REQUIRE_PRESELECTION=0"
  else
    printf '%s' \
      "RJ_AUAU_BDT_EXTRACT_ONLY=1;RJ_AUAU_BDT_TRAINING_TREE=1;RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES=${legacy_tree_max_entries};RJ_AUAU_BDT_NPB_DATA_TAGGING=0;RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=1;RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0;RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0"
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

analysis_tag_for_dataset() {
  case "$1" in
    isSimInclusive) echo isSimInclusive ;;
    isSimEmbedded) echo isSimEmbedded ;;
    isSimEmbeddedInclusive) echo isSimEmbeddedInclusive ;;
    *) echo isSim ;;
  esac
}

# Verify the exact dry-materialized unit that will be submitted.  This is the
# pre-submission ownership proof: one descriptor, one argument row, one
# fanout row, one RecoilJets instance, and one multiview sidecar assignment.
# The tab-separated return value is consumed verbatim by submit_row.
verify_materialized_row_contract() {
  local row_id="$1" system="$2" dataset="$3" sample="$4" row_output="$5" sidecar="$6" row_submit="$7"
  local -a sub_files=() fanout_files=() sidecar_values=() id_file_values=() id_dirs_values=()
  local submit_file args_file args_line chunk_list chunk_sha fanout_file fanout_sha fanout_line
  local arg_sample arg_chunk arg_dataset arg_cluster arg_events arg_index arg_none arg_dest arg_extra
  local fan_dest fan_cfg fan_pre fan_tight fan_non fan_cone fan_sliding fan_fixed
  local chunk_tag analysis_tag analysis_root owner_count queue_count value
  local frozen_executable snapshot_dir snapshot_header snapshot_header_sha
  local snapshot_calo snapshot_calo_sha snapshot_analysis snapshot_analysis_sha expected_analysis_sha

  while IFS= read -r value; do sub_files+=( "$value" ); done \
    < <(find "$row_submit" -maxdepth 1 -type f -name '*.sub' -print | sort)
  [[ "${#sub_files[@]}" == 1 ]] || die "${row_id} must dry-materialize exactly one submit descriptor, observed=${#sub_files[@]}"
  submit_file="${sub_files[0]}"
  frozen_executable="$(condor_field "$submit_file" executable)"
  [[ -s "$frozen_executable" ]] ||
    die "${row_id} descriptor does not bind a readable frozen executable"
  snapshot_dir="$(cd "$(dirname "$frozen_executable")" && pwd -P)"
  snapshot_header="${snapshot_dir}/PhotonClusterBuilder.h"
  snapshot_calo="${snapshot_dir}/lib/libcalo_reco.so"
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
  require_file_hash "${row_id} snapshotted analysis library" \
    "$snapshot_analysis" "$expected_analysis_sha"
  [[ ! -e "${snapshot_dir}/lib/libphoton_cluster_builder_override.so" ]] ||
    die "${row_id} snapshot contains an undeclared PhotonClusterBuilder override library"
  snapshot_header_sha="$(sha_file "$snapshot_header")"
  snapshot_calo_sha="$(sha_file "$snapshot_calo")"
  snapshot_analysis_sha="$(sha_file "$snapshot_analysis")"
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
    die "${row_id} argument row violates the one-proc smoke contract"
  [[ "$arg_dest" == "$row_output/"* ]] || die "${row_id} argument destination escapes its owned output namespace"
  chunk_list="$arg_chunk"
  [[ -s "$chunk_list" ]] || die "${row_id} staged chunk list is missing: ${chunk_list}"
  [[ "$(grep -Evc '^[[:space:]]*($|#)' "$chunk_list")" == 5 ]] ||
    die "${row_id} staged chunk must contain one exact five-file input tuple"
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
  fanout_file="${id_file_values[0]}"
  [[ -s "$fanout_file" ]] || die "${row_id} materialized fanout contract is missing: ${fanout_file}"
  [[ "$(grep -Evc '^[[:space:]]*($|#)' "$fanout_file")" == 1 ]] ||
    die "${row_id} fanout contract must contain exactly one RecoilJets row"
  fanout_line="$(grep -Ev '^[[:space:]]*($|#)' "$fanout_file")"
  IFS='|' read -r fan_dest fan_cfg fan_pre fan_tight fan_non fan_cone fan_sliding fan_fixed <<< "$fanout_line"
  [[ -n "$fan_dest" && -n "$fan_cfg" && -n "$fan_pre" && -n "$fan_tight" && -n "$fan_non" ]] ||
    die "${row_id} fanout row is incomplete"
  [[ "$fan_dest" == "$row_output/"* ]] || die "${row_id} fanout destination escapes its owned output namespace"
  [[ "${fan_cone:-}" == "$extraction_cone_r" || "${fan_cone:-}" == "0.4" ]] ||
    die "${row_id} fanout row does not carry the nominal R=${extraction_cone_r} extraction cone"
  owner_count=$(( ${#sub_files[@]} * $(wc -l < "$args_file" | tr -d ' ') * $(grep -Evc '^[[:space:]]*($|#)' "$fanout_file") * ${#sidecar_values[@]} ))
  [[ "$owner_count" == 1 ]] || die "${row_id} multiview sidecar owner multiplicity is ${owner_count}, expected 1"

  chunk_tag="$(basename "$chunk_list")"
  chunk_tag="${chunk_tag%.list}"
  analysis_tag="$(analysis_tag_for_dataset "$dataset")"
  analysis_root="${fan_dest}/${sample}/RecoilJets_${analysis_tag}_${fan_cfg}_${chunk_tag}.root"
  fanout_sha="$(sha_file "$fanout_file")"
  require_sha "${row_id} fanout contract" "$fanout_sha"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$submit_file" "$args_file" "$chunk_list" "$chunk_sha" \
    "$fanout_file" "$fanout_sha" "$analysis_root" "$owner_count" \
    "$snapshot_dir" "$snapshot_header" "$snapshot_header_sha" \
    "$snapshot_calo" "$snapshot_calo_sha" "$snapshot_analysis" "$snapshot_analysis_sha"
}

submit_row() {
  local row_id="$1" system="$2" lane="$3" dataset="$4" sample="$5" role="$6" _mb_gate="$7" row_match="$8"
  local config row_output sidecar row_submit library extra row_log materialize_log cluster proc cluster_proc
  local contract submit_file args_file chunk_list chunk_sha fanout_file fanout_sha analysis_root owner_count
  local snapshot_dir snapshot_header snapshot_header_sha snapshot_calo snapshot_calo_sha
  local snapshot_analysis snapshot_analysis_sha
  local log_template out_template err_template log_path out_path err_path submitted_args expected_args args_sha
  local -a materialize_contract_env=(
    RJ_REPLAY_FOUNDATION_CANARY=1
    RJ_REPLAY_LANE="$lane"
    RJ_REPLAY_SCHEMA_SHA256="$RJ_THE134_REPLAY_SCHEMA_SHA256"
  )
  config="$(manifest_field "$row_id" 12)"
  row_output="$(manifest_field "$row_id" 27)"
  sidecar="$(manifest_field "$row_id" 28)"
  row_submit="$(manifest_field "$row_id" 29)"
  if [[ "$system" == pp ]]; then library="$pp_library"; else library="$auau_library"; fi
  extra="$(common_extra_env "$row_id" "$system" "$lane" "$dataset" "$sample" "$role");$(system_extra_env "$system" "$role")"
  mkdir -p "$row_submit" "$(dirname "$sidecar")"
  materialize_log="${evidence_root}/materialize_${row_id}.log"
  row_log="${evidence_root}/submit_${row_id}.log"

  say "MATERIALIZE row=${row_id} dataset=${dataset} sample=${sample} role=${role}"
  if [[ "$system" == pp ]]; then
    env "${materialize_contract_env[@]}" \
      RJ_CODEX_CHAT_NAME="$RJ_CODEX_CHAT_NAME" RJ_CODEX_THREAD_ID="$RJ_CODEX_THREAD_ID" \
      RJ_SIM_ROOT_OVERRIDE="$sim_root" RJ_CONFIG_YAML="$config" RJ_PP_LIBRARY_OVERRIDE="$library" \
      RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE="$photon_cluster_builder_header" \
      RJ_CALO_RECO_LIBRARY_OVERRIDE="$calo_reco_library" \
      RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE= \
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
      "$submitter" "$dataset" condorDoAllSmoke groupSize 1 maxJobs 1 "SAMPLE=${sample}" \
      2>&1 | tee "$materialize_log"
  else
    env "${materialize_contract_env[@]}" \
      RJ_CODEX_CHAT_NAME="$RJ_CODEX_CHAT_NAME" RJ_CODEX_THREAD_ID="$RJ_CODEX_THREAD_ID" \
      RJ_SIM_ROOT_OVERRIDE="$sim_root" RJ_CONFIG_YAML="$config" RJ_AUAU_LIBRARY_OVERRIDE="$library" \
      RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE="$photon_cluster_builder_header" \
      RJ_CALO_RECO_LIBRARY_OVERRIDE="$calo_reco_library" \
      RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE= \
      RJ_DAG_DRYRUN=1 \
      RJ_AUTO_MERGE=0 RJ_STAGE_EMAIL_MODE=none RJ_CLEAN_OUTPUT_BASE=0 \
      RJ_REQUEST_MEMORY=8000MB RJ_REQUIRE_NON_TINY_OUTPUT=1 RJ_MIN_OUTPUT_BYTES=50000 \
      RJ_FAIL_ON_MISSING_CALO_INPUT=1 RJ_VALIDATE_SIM_INPUT_PATHS=1 RJ_VALIDATE_SIM_INPUT_MAX_LINES=1 \
      RJ_PROFILE_JOB=1 RJ_JOB_HEARTBEAT_SECONDS=120 RJ_PROFILE_LABEL="${tag}_${row_id}" \
      RJ_SMOKE_OUTPUT_BASE="$row_output" RJ_SMOKE_SIM_NEVENTS="$event_limit_per_job" \
      RJ_SUBMISSION_NAMESPACE="$row_id" RJ_CONDOR_SUB_DIR="$row_submit" \
      RJ_PHOTON_ID_ROW_MATCH="$row_match" RJ_ID_FANOUT_MAX_ROWS=1 \
      RJ_SUBMIT_EXTRA_ENV="$extra" \
      "$submitter" "$dataset" condorDoAllSmoke groupSize 1 maxJobs 1 "SAMPLE=${sample}" \
      2>&1 | tee "$materialize_log"
  fi

  grep -F 'RECOILJETS_SMOKETEST_DRYRUN_V1' "$materialize_log" >/dev/null ||
    die "${row_id} submitter did not attest dry materialization"
  ! grep -Eq '[0-9]+ job\(s\) submitted to cluster [0-9]+' "$materialize_log" ||
    die "${row_id} dry materialization unexpectedly submitted a Condor job"
  contract="$(verify_materialized_row_contract "$row_id" "$system" "$dataset" "$sample" "$row_output" "$sidecar" "$row_submit")"
  IFS=$'\t' read -r submit_file args_file chunk_list chunk_sha fanout_file fanout_sha analysis_root owner_count \
    snapshot_dir snapshot_header snapshot_header_sha snapshot_calo snapshot_calo_sha \
    snapshot_analysis snapshot_analysis_sha <<< "$contract"

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
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$row_id" "$cluster_proc" "$submit_file" "$args_file" "$chunk_list" "$chunk_sha" \
    "$fanout_file" "$fanout_sha" "$log_path" "$out_path" "$err_path" \
    "$analysis_root" "$sidecar" "$owner_count" "$args_sha" \
    "$snapshot_dir" "$snapshot_header" "$snapshot_header_sha" \
    "$snapshot_calo" "$snapshot_calo_sha" "$snapshot_analysis" "$snapshot_analysis_sha" \
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
  IFS=$'\t' read -r code_sha replay_schema_sha training_schema_sha semantic_sha <<< "$hashes"
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
}

submit_all() {
  assert_fresh_submission
  preflight
  mkdir -p "$evidence_root" "$submit_root"
  printf 'row_id\tcluster_proc\tsubmit_log\tsubmitted_at_utc\n' > "$submission_journal"
  printf 'row_id\tcluster_proc\tsubmit_file\targs_file\tstaged_chunk_list\tstaged_chunk_sha256\tfanout_contract_file\tfanout_contract_sha256\tcondor_log\tcondor_stdout\tcondor_stderr\tanalysis_output_root\tmultiview_sidecar\tsidecar_owner_count\tsubmitted_args_sha256\tsnapshot_dir\tsnapshot_builder_header\tsnapshot_builder_header_sha256\tsnapshot_calo_reco_library\tsnapshot_calo_reco_library_sha256\tsnapshot_analysis_library\tsnapshot_analysis_library_sha256\n' > "$submission_receipt"
  while IFS='|' read -r row_id system lane dataset sample role mb_gate row_match; do
    submit_row "$row_id" "$system" "$lane" "$dataset" "$sample" "$role" "$mb_gate" "$row_match"
  done < <(emit_matrix)
  [[ "$(wc -l < "$submission_receipt" | tr -d ' ')" == 14 ]] || die "submission receipt is incomplete"
  [[ "$(wc -l < "$submission_journal" | tr -d ' ')" == 14 ]] || die "submission journal is incomplete"
  python3 - "$source_provenance_json" <<'PY'
import json
import sys
payload = json.load(open(sys.argv[1]))
assert payload["schema"] == "THE134_SOURCE_PROVENANCE_V1"
assert len(payload["inputs"]) == 13
PY
  say "SUBMISSION_PASS rows=13 receipt=${submission_receipt}"
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
  done < <(emit_matrix)
  [[ "$(wc -l < "$submission_receipt" | tr -d ' ')" == 14 ]] || die "resumed submission receipt is incomplete"
  [[ "$(wc -l < "$submission_journal" | tr -d ' ')" == 14 ]] || die "resumed submission journal is incomplete"
  say "RESUME_SUBMISSION_PASS rows=13 receipt=${submission_receipt}"
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
  python3 - "$submission_receipt" "$submission_manifest" "$root_health_join_certificate" <<'PY'
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

receipt_path, manifest_path, certificate_path = map(Path, sys.argv[1:])
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
        if int(source_tree.GetEntries()) != 1:
            raise ValueError(f"expected exactly one replay source row, observed={int(source_tree.GetEntries())}")
        source_ids = {identity(row, "source_occurrence_id") for row in source_tree}
        if len(source_ids) != 1:
            raise ValueError(f"expected exactly one replay source occurrence, observed={len(source_ids)}")
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
                "full_training_authority": 0,
            }
        )
        sidecar.Close()
        analysis.Close()
    except Exception as exc:
        failures.append(f"{row_id}:{exc}")
        report.update({"status": "FAIL", "failure": str(exc)})
    reports.append(report)

if len(receipts) != 13 or len(manifests) != 13:
    failures.append(f"row_closure:receipts={len(receipts)} manifests={len(manifests)}")

payload = {
    "schema": "THE134_SMOKE_ROOT_HEALTH_IDENTITY_JOIN_V1",
    "status": "PASS" if not failures else "FAIL",
    "scope": "smoke",
    "full_training_authority": 0,
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

validate_outputs() {
  local row_id cluster_proc queue_state history_state job_status exit_code
  [[ -s "$validator" ]] || die "THE-134 matrix preparer is missing: ${validator}"
  [[ -s "$submission_journal" && "$(wc -l < "$submission_journal" | tr -d ' ')" == 14 ]] ||
    die "complete 13-row submission journal is required"
  [[ -s "$submission_receipt" && "$(wc -l < "$submission_receipt" | tr -d ' ')" == 14 ]] ||
    die "complete 13-row submission receipt is required"
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

  say "SMOKE_VALIDATION_PASS source_categories=13 full_training_authority=0 root_join_certificate=${root_health_join_certificate} pp_audit=${evidence_root}/pp_h70_source_complete_smoke_audit.json auau_audit=${evidence_root}/auau_h70_source_complete_smoke_audit.json"
}

case "$mode" in
  inventory) emit_source_hash_manifest ;;
  preflight) preflight ;;
  submit) submit_all ;;
  resume-submit) resume_submit ;;
  status) status ;;
  validate) validate_outputs ;;
  *) die "usage: $0 inventory|preflight|submit|resume-submit|status|validate" ;;
esac
