#!/usr/bin/env bash
# Plan or submit the bounded 12-lane photon+jet half of the PPG12
# stitched-purity admission canary.
#
# This driver deliberately submits exactly one full five-file worker for each
# Photon5/10/20 x 0mrad/1p5mrad x SI/DI lane.  It never submits raw Condor
# jobs itself, never merges outputs, and never reuses an output namespace.

set -Eeuo pipefail
IFS=$'\n\t'
umask 077

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  # shellcheck disable=SC1091
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  # shellcheck disable=SC1091
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  printf '[PPG12-PHOTON-CANARY][ERROR] cannot locate recoiljets_io_paths.sh\n' >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

say() { printf '[PPG12-PHOTON-CANARY] %s\n' "$*"; }
die() { printf '[PPG12-PHOTON-CANARY][ERROR] %s\n' "$*" >&2; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  submit_ppg12_stitched_purity_photon_canaries.sh [plan]
  submit_ppg12_stitched_purity_photon_canaries.sh --submit \
    --token "SUBMIT_<campaign-tag>"

The default mode is plan.  Submission requires both the explicit --submit
flag and the exact campaign-bound token printed by plan mode.
EOF
}

mode="plan"
if [[ "${1:-}" == "plan" ]]; then shift; fi
submit_flag=0
provided_token=""
while (( $# > 0 )); do
  case "$1" in
    --submit)
      [[ "$submit_flag" -eq 0 ]] || die "--submit may be specified only once"
      submit_flag=1
      mode="submit"
      shift
      ;;
    --token)
      (( $# >= 2 )) || die "--token requires a value"
      provided_token="$2"
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

campaign_tag="${RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG:-ppg12_stitched_purity_photon_canary_$(date -u +%Y%m%dT%H%M%SZ)_$$}"
rj_validate_campaign_tag "$campaign_tag" || die "unsafe campaign tag: $campaign_tag"

expected_token="SUBMIT_${campaign_tag}"
output_root="${RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT:-$(rj_thesis_ana_root)/sim/${campaign_tag}}"
evidence_dir="${RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR:-${repo_root}/evidence/qa/${campaign_tag}}"
manifest_tsv="${evidence_dir}/photon_canary_lanes.tsv"
manifest_json="${evidence_dir}/photon_canary_plan.json"
receipt_tsv="${evidence_dir}/submission_receipts.tsv"
config_yaml="${RJ_PPG12_PHOTON_CANARY_CONFIG_YAML:-${repo_root}/macros/analysis_config.yaml}"
oracle_root_base="${RJ_PPG12_ORACLE_ROOT_BASE:-/sphenix/user/shuhangli/ppg12/FunWithxgboost}"
comparator_rel="scripts/diagnostics/pp_currentian/compare_ppg12_recoiljets_same_cluster_features.py"
replay_seed_sequence='2991264730,4256268992,2394322166,874466025,2240380304'
expected_pedestal_sequence=534

if [[ -n "${RJ_PPG12_PHOTON_CANARY_SUBMITTER:-}" ]]; then
  submitter="$RJ_PPG12_PHOTON_CANARY_SUBMITTER"
elif [[ -x "${repo_root}/RecoilJets_Condor_submit.sh" ]]; then
  # This is the canonical location in the SDCC control checkout.
  submitter="${repo_root}/RecoilJets_Condor_submit.sh"
else
  # Local source-tree fallback used for plan validation and unit tests.
  submitter="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
fi

[[ "$output_root" == /* ]] || die "output root must be absolute: $output_root"
[[ "$evidence_dir" == /* ]] || die "evidence directory must be absolute: $evidence_dir"
[[ "$oracle_root_base" == /* ]] || die "oracle ROOT base must be absolute: $oracle_root_base"
case "$output_root" in
  /sphenix/*|/tmp/*|/private/tmp/*) ;;
  *) die "output root must be an isolated /sphenix or temporary test path: $output_root" ;;
esac

for forbidden in \
  RANDOMSEED \
  RJ_PPG12_CLOSURE_CANARY \
  RJ_PPG12_CLOSURE_CANARY_ID \
  RJ_PPG12_PPSIM_REPLAY_SEEDS \
  RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE \
  RJ_PPG12_PPSIM_FIXED_RANDOMSEED \
  RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE \
  RJ_DISABLE_JES_CDB_AUDIT \
  RJ_PPG12_PEDESTAL_OVERRIDE \
  RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL \
  RJ_REQUIRE_SIM_GLOBAL \
  RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT \
  RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX \
  RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX \
  RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX
do
  if [[ -n "${!forbidden+x}" ]]; then
    die "inherited physics override is forbidden for the exact oracle canary: ${forbidden}"
  fi
  case ";${RJ_SUBMIT_EXTRA_ENV:-};" in
    *";${forbidden}="*)
      die "RJ_SUBMIT_EXTRA_ENV may not carry driver-owned oracle control: ${forbidden}"
      ;;
  esac
done

mkdir -p "$evidence_dir"

write_lane_table() {
  printf 'lane_id\tfamily\tsample_key\tperiod\tinteraction\trecoil_sample\toutput_base\tsubmission_namespace\tppg12_oracle_root\tppg12_tree\trecoil_tree\tgroup_size\tmax_jobs\tauto_merge\ttraining_tree_max_entries\tphoton_id_row\treplay_seed_sequence\texpected_pedestal_sequence\n' > "$manifest_tsv"
  local photon period interaction sample suffix lane_id lane_slug output_base namespace oracle_sample oracle_root
  for photon in 5 10 20; do
    for period in 0mrad 1p5mrad; do
      for interaction in si di; do
        suffix=""
        [[ "$interaction" == "di" ]] && suffix="_double"
        sample="run28_photonjet${photon}${suffix}"
        lane_id="photon:photon${photon}:${period}:${interaction}"
        lane_slug="photon${photon}_${period}_${interaction}"
        output_base="${output_root}/${lane_slug}"
        namespace="${campaign_tag}_${lane_slug}"
        oracle_sample="photon${photon}${suffix}"
        oracle_root="${oracle_root_base}/${oracle_sample}/bdt_split.root"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$lane_id" photon "photon${photon}" "$period" "$interaction" \
          "$sample" "$output_base" "$namespace" "$oracle_root" slimtree \
          AuAuPhotonIDTrainingTree 5 1 false 0 'newPPG12|newPPG12|newPPG12' \
          "$replay_seed_sequence" "$expected_pedestal_sequence" \
          >> "$manifest_tsv"
      done
    done
  done
}

write_json_manifest() {
  PPG12_CANARY_TSV="$manifest_tsv" \
  PPG12_CANARY_JSON="$manifest_json" \
  PPG12_CANARY_CAMPAIGN="$campaign_tag" \
  PPG12_CANARY_OUTPUT_ROOT="$output_root" \
  PPG12_CANARY_SUBMITTER="$submitter" \
  PPG12_CANARY_CONFIG="$config_yaml" \
  PPG12_CANARY_COMPARATOR="$comparator_rel" \
  PPG12_CANARY_TOKEN="$expected_token" \
  python3 - <<'PY'
import csv
import datetime as dt
import json
import os
from pathlib import Path

tsv = Path(os.environ["PPG12_CANARY_TSV"])
with tsv.open(newline="") as stream:
    lanes = list(csv.DictReader(stream, delimiter="\t"))

for lane in lanes:
    lane["group_size"] = int(lane["group_size"])
    lane["max_jobs"] = int(lane["max_jobs"])
    lane["auto_merge"] = lane["auto_merge"] == "true"
    lane["training_tree_max_entries"] = int(lane["training_tree_max_entries"])
    lane["replay_seed_sequence"] = [
        int(value) for value in lane["replay_seed_sequence"].split(",")
    ]
    lane["expected_pedestal_sequence"] = int(lane["expected_pedestal_sequence"])
    lane["submit_argv"] = [
        "isSim", "condorDoAllDirect", "groupSize", "5", "maxJobs", "1",
        f"SAMPLE={lane['recoil_sample']}",
    ]
    lane["source_contract"] = (
        "run28_double_none_g4_truthjet_none_none_g4_only_rebuild"
        if lane["interaction"] == "di"
        else "run28_single_interaction_none_g4_truthjet_none_none_g4_only_rebuild"
    )
    lane["five_column_input_graph"] = [
        "NONE", "G4Hits", "DST_TRUTH_JET", "NONE", "NONE"
    ]
    lane["registered_input_graph"] = ["G4Hits", "DST_TRUTH_JET"]
    lane["dst_global_registered"] = False
    lane["guarded_interaction_source_contract"] = True
    lane["jes_cdb_audit_disabled"] = True

payload = {
    "schema": "ppg12-stitched-purity-photon-canary-plan/v4",
    "status": "PLANNED",
    "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    "campaign_tag": os.environ["PPG12_CANARY_CAMPAIGN"],
    "dataset": "isSim",
    "canonical_submitter": os.environ["PPG12_CANARY_SUBMITTER"],
    "config_yaml": os.environ["PPG12_CANARY_CONFIG"],
    "output_root": os.environ["PPG12_CANARY_OUTPUT_ROOT"],
    "submission_token": os.environ["PPG12_CANARY_TOKEN"],
    "bounded_contract": {
        "lane_count": 12,
        "group_size": 5,
        "max_jobs_per_lane": 1,
        "input_path_validation_rows": 5,
        "full_worker_input": True,
        "automatic_merge": False,
        "candidate_tree": "AuAuPhotonIDTrainingTree",
        "candidate_tree_max_entries": 0,
        "photon_id_row": ["newPPG12", "newPPG12", "newPPG12"],
        "classification_energy": "unsmeared_calibrated_et",
        "legacy_multiplicative_cluster_smearing": "disabled",
        "response_smearing": "ppg12_additive_response_only",
        "rebuild_calo_from_g4": "all lanes",
        "g4_only": "all lanes",
        "rng_contract": "historical_fifo_replay_v2",
        "replay_seed_sequence": [
            2991264730, 4256268992, 2394322166, 874466025, 2240380304
        ],
        "expected_pedestal_sequence": 534,
        "estimator_random_seed": 42,
        "estimator_toy_count": 20000,
        "jes_cdb_audit_disabled": True,
        "si_five_column_input_graph": [
            "NONE", "G4Hits", "DST_TRUTH_JET", "NONE", "NONE"
        ],
        "di_five_column_input_graph": [
            "NONE", "G4Hits", "DST_TRUTH_JET", "NONE", "NONE"
        ],
        "si_registered_input_graph": [
            "G4Hits", "DST_TRUTH_JET"
        ],
        "di_registered_input_graph": ["G4Hits", "DST_TRUTH_JET"],
        "dst_global_registered": False,
        "table_qa": True,
        "fig13_parity_qa": True,
        "fig11_sideband_diagnostic": True,
        "random_seed_policy": "identical_historical_fifo_replay_for_both_executables_v2",
        "external_normalization": "forbidden",
    },
    "executable_oracle": {
        "ppg12_tree": "slimtree",
        "recoil_tree": "AuAuPhotonIDTrainingTree",
        "comparator": os.environ["PPG12_CANARY_COMPARATOR"],
        "join_key": ["event_identity", "candidate_identity"],
        "mapping_is_lane_explicit": True,
    },
    "lanes": lanes,
}

out = Path(os.environ["PPG12_CANARY_JSON"])
out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
}

write_lane_table
write_json_manifest

print_contract() {
  cat <<EOF
PPG12_STITCHED_PURITY_PHOTON_CANARY_V3
mode=${mode}
campaign_tag=${campaign_tag}
lane_count=12
samples=Photon5,Photon10,Photon20
periods=0mrad,1p5mrad
interaction_modes=SI,DI
group_size=5
max_jobs_per_lane=1
automatic_merge=disabled
training_tree=AuAuPhotonIDTrainingTree
training_tree_max_entries=0_unlimited
photon_id_row=newPPG12|newPPG12|newPPG12
replay_seed_sequence=${replay_seed_sequence}
expected_pedestal_sequence=${expected_pedestal_sequence}
estimator_random_seed=42
estimator_toy_count=20000
jes_cdb_audit_disabled=true
output_root=${output_root}
manifest=${manifest_json}
lane_table=${manifest_tsv}
submitter=${submitter}
submit_token=${expected_token}
EOF
}

assert_submission_preflight() {
  (( submit_flag == 1 )) || die "submission requires --submit"
  [[ "$provided_token" == "$expected_token" ]] || die "submit token mismatch; run plan and use its exact campaign-bound token"
  [[ -x "$submitter" ]] || die "canonical submitter is missing or not executable: $submitter"
  [[ -n "${RJ_PPG12_PHOTON_CANARY_CONFIG_YAML:-}" ]] || \
    die "submission requires RJ_PPG12_PHOTON_CANARY_CONFIG_YAML pointing to the archived 20260709 pp config"
  [[ -s "$config_yaml" ]] || die "canonical config is missing or empty: $config_yaml"
  [[ -n "${RJ_CODEX_CHAT_NAME:-}" ]] || die "RJ_CODEX_CHAT_NAME is required for submission provenance"
  [[ -n "${RJ_CODEX_THREAD_ID:-}" ]] || die "RJ_CODEX_THREAD_ID is required for submission provenance"
  [[ ! -e "$output_root" ]] || die "output root already exists; use a fresh campaign tag: $output_root"

  local output_base
  while IFS=$'\t' read -r _ _ _ _ _ _ output_base _rest; do
    [[ "$output_base" == "$output_root"/* ]] || die "lane output escapes campaign root: $output_base"
    [[ ! -e "$output_base" ]] || die "lane output already exists: $output_base"
  done < <(tail -n +2 "$manifest_tsv")
}

submit_lane() {
  local lane_id="$1" period="$2" interaction="$3" sample="$4" output_base="$5" namespace="$6"
  say "submit lane=${lane_id} sample=${sample} output=${output_base}"

  (
    # Strip every interaction-specific variable before constructing this lane.
    unset RJ_PPG12_PHOTON_YIELD_DOUBLE
    unset RJ_PPG12_PERIOD_STRICT_DI
    unset RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4
    unset RJ_PPG12_PPSIM_G4_ONLY
    unset RJ_PPG12_PPSIM_REPLAY_SEEDS
    unset RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE
    unset RJ_PPG12_PPSIM_FIXED_RANDOMSEED
    unset RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE
    unset RJ_DISABLE_JES_CDB_AUDIT
    unset RJ_PPG12_PEDESTAL_OVERRIDE
    unset RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL
    unset RJ_SIM_ALLOW_NONE_LISTS
    unset RJ_SIM_SIGNAL_SAMPLE_SET

    export RJ_CONFIG_YAML="$config_yaml"
    export RJ_DEST_BASE_OVERRIDE="$output_base"
    export RJ_SUBMISSION_NAMESPACE="$namespace"
    export RJ_PPG12_CLOSURE_CANARY=1
    export RJ_PPG12_CLOSURE_CANARY_ID="$lane_id"
    export RJ_PPG12_PPSIM_REPLAY_SEEDS="$replay_seed_sequence"
    export RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE="$expected_pedestal_sequence"
    export RJ_DISABLE_JES_CDB_AUDIT=1
    export RJ_PPG12_PHOTON_YIELD=1
    export RJ_PPG12_PERIOD="$period"
    export RJ_PPG12_CROSSING_PERIOD="$period"
    export RJ_PPG12_PERIOD_USE_LUMI_WEIGHT=1
    export RJ_PPG12_PERIOD_ALLOW_ALL_SIM=0
    export RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE=0
    export RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE=0
    export RJ_PPG12_TABLE_QA=1
    export RJ_PPG12_TABLE_QA_NPB_DATA_TAGGING=1
    export RJ_PPG12_FIG13_PARITY_QA=1
    export RJ_PPG12_FIG11_SB_DIAGNOSTIC=1
    export RJ_PP_PHOTONID_TRAINING_TREE=1
    export RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES=0
    export RJ_PP_PHOTONID_SOURCE_ROLE=signal
    export RJ_PP_PHOTONID_PPG12_FILTER=1
    export RJ_PP_PHOTONID_REQUIRE_PRESELECTION=0
    export RJ_PHOTON_ID_ROW_MATCH='newPPG12|newPPG12|newPPG12'
    export RJ_DISABLE_ID_FANOUT=1
    export RJ_ID_FANOUT_MAX_ROWS=1
    export RJ_DISABLE_JET_PT_INTERNALIZATION=1
    export RJ_DISABLE_DPHI_INTERNALIZATION=1
    export RJ_DIRECT_DST_DOALL=1
    export RJ_DIRECT_NEVENTS=0
    export RJ_AUTO_MERGE=0
    # This bounded canary submits only the first five matched rows. Validate
    # every path the worker can consume without paying the full-production
    # 10,000-row existence-scan cost. The submitter's separate source-contract
    # audit still checks the complete staged list for SI/DI stream identity.
    export RJ_VALIDATE_SIM_INPUT_MAX_LINES=5
    export RJ_INTERNAL_JET_PT_MINS='5.0,7.0,10.0,12.0'
    export RJ_INTERNAL_DPHI_PI_FRACTIONS='0.5,0.875'
    export RJ_REQUEST_MEMORY_MB=6000
    export RJ_REQUIRE_NON_TINY_OUTPUT=1
    export RJ_MIN_OUTPUT_BYTES=50000
    export RJ_FAIL_ON_MISSING_CALO_INPUT=0
    # Both preserved executable-oracle graphs register only G4Hits and the
    # truth-jet stream; calo, global, and MBD columns are explicitly NONE.
    export RJ_SIM_ALLOW_NONE_LISTS=1
    export RJ_PPG12_PPSIM_G4_ONLY=1

    if [[ "$interaction" == "di" ]]; then
      export RJ_SIM_SIGNAL_SAMPLE_SET=ppg12_double
      export RJ_PPG12_PHOTON_YIELD_DOUBLE=1
      export RJ_PPG12_PERIOD_STRICT_DI=1
      export RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
      # DI uses the dual-interaction G4Hits and truth-jet source set.
    else
      export RJ_PPG12_PHOTON_YIELD_DOUBLE=0
      # SI uses the single-interaction G4Hits and truth-jet source set.  The
      # preserved executable reconstructs MBD, towers, and clusters rather
      # than registering prebuilt calo/global/MBD streams.
      export RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
    fi

    "$submitter" isSim condorDoAllDirect groupSize 5 maxJobs 1 "SAMPLE=${sample}"
  )

  printf '%s\t%s\t%s\t%s\t%s\tSUCCESS\n' \
    "$lane_id" "$sample" "$period" "$interaction" "$output_base" >> "$receipt_tsv"
}

print_contract

if [[ "$mode" == "plan" ]]; then
  [[ -z "$provided_token" ]] || die "--token is valid only with submit mode"
  exit 0
fi

assert_submission_preflight
printf 'lane_id\trecoil_sample\tperiod\tinteraction\toutput_base\tstatus\n' > "$receipt_tsv"

while IFS=$'\t' read -r lane_id _family _sample_key period interaction sample output_base namespace _rest; do
  submit_lane "$lane_id" "$period" "$interaction" "$sample" "$output_base" "$namespace"
done < <(tail -n +2 "$manifest_tsv")

say "submitted exactly 12 bounded lanes; receipts=${receipt_tsv}"
