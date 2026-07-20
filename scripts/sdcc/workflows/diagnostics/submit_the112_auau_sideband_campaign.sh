#!/usr/bin/env bash
set -euo pipefail

# THE-112 is deliberately split into three gates:
#   1. rank finite score bands from historical and continuous canary surfaces;
#   2. confirm one frozen band with a small matched data/signal/inclusive canary;
#   3. submit the full matched production in separately guarded lanes.
#
# This launcher never merges output, never releases/removes/resubmits jobs, and
# never changes the classifier or WP80 surface.  Each submitting mode writes a
# one-shot receipt so an accidental second invocation fails closed.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${script_dir}/../../runtime/io/recoiljets_io_paths.sh" ]]; then
  source "${script_dir}/../../runtime/io/recoiljets_io_paths.sh"
elif [[ -f scripts/sdcc/runtime/io/recoiljets_io_paths.sh ]]; then
  source scripts/sdcc/runtime/io/recoiljets_io_paths.sh
else
  echo "[THE112][ERROR] Cannot locate recoiljets_io_paths.sh" >&2
  exit 2
fi

repo_root="$(rj_find_repo_root "$script_dir")"
cd "$repo_root"

say() { printf '[THE112] %s\n' "$*"; }
die() { printf '[THE112][ERROR] %s\n' "$1" >&2; exit "${2:-2}"; }

mode="${1:-print}"
campaign_tag="${RJ_THE112_CAMPAIGN_TAG:-the112_auau_combined_bdt_sideband_20260719}"
family="${RJ_THE112_FAMILY:-the112_auau_sideband}"
scan_yaml="${RJ_THE112_SCAN_CONFIG_YAML:-${repo_root}/macros/analysis_config_the112_auau_combined_bdt_triplet.yaml}"
build_root="${RJ_THE112_BUILD_ROOT:-${repo_root}/.recoiljets_tmp/the112_sideband_build_20260719}"
build_dir="${build_root}/build"
install_dir="${build_root}/install"
auau_library="${RJ_THE112_AUAU_LIBRARY:-${install_dir}/lib/libRecoilJetsAuAu.so}"
evidence_dir="${RJ_THE112_EVIDENCE_DIR:-${repo_root}/evidence/qa/the112_sideband_20260719}"
manifest="${evidence_dir}/campaign_manifest.tsv"

model_path="${RJ_THE112_MODEL_PATH:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the111_models/the111_combined_corrected_shower_ppg12_labels_20260719_1618/combined/auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.root}"
model_sha256="d50c69ec98558cb80730ab45fe6801d4accbcf2221c482af91e8899cede1c925"
wp_intercept="0.5544148693"
wp_slope="0.0015499421"
wp_min_et="15"
wp_max_et="35"

pairing_report="${RJ_THE112_PAIRING_REPORT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau/dst_pairing_report_auau.txt}"
pairing_report_sha256="f9404bddd6e9d83e44d7ebbf186711a616dbc83a96141e0a1f915604de419bd5"
pair_list_dir="${RJ_THE112_PAIR_LIST_DIR:-/sphenix/u/patsfan753/scratch/thesisAnalysis/dst_lists_auau}"
expected_data_runs=2776
expected_data_pairs=492280

ranker="${RJ_THE112_RANKER:-${repo_root}/scripts/diagnostics/auau_bdt/rank_the112_sideband_scan.py}"
canary_validator="${RJ_THE112_CANARY_VALIDATOR:-${repo_root}/scripts/diagnostics/auau_bdt/validate_the112_sideband_canary.py}"
# The mapped SDCC upload contract installs the protected hard runtime entrypoint
# at the checkout root.  Local worktrees retain the canonical source path.
default_submitter="${repo_root}/RecoilJets_Condor_submit.sh"
if [[ ! -f "$default_submitter" ]]; then
  default_submitter="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
fi
submitter="${RJ_THE112_SUBMITTER:-$default_submitter}"
analysis_python="${RJ_THE112_PYTHON:-${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}}"
expected_iso_views="isoR40_isSliding:0.40:true:0.0,isoR30_isSliding:0.30:true:0.0"
historical_data="${RJ_THE112_HISTORICAL_DATA_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/the100_auau_dualview_20260714/auau/RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband_baseVariant.root}"
historical_signal="${RJ_THE112_HISTORICAL_SIGNAL_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/the100_auau_dualview_20260714/signal/simembedded/preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband_baseVariant/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root}"
historical_inclusive="${RJ_THE112_HISTORICAL_INCLUSIVE_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/the100_auau_dualview_20260714/inclusive/simembeddedinclusive/preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband_baseVariant/embeddedJet12and20and30and40merged_SIM/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root}"

bulk_root="${RJ_THE112_BULK_ROOT:-$(rj_recoiljets_bulk_root)/${family}/${campaign_tag}}"
data_bulk="${bulk_root}/data"
signal_bulk="${bulk_root}/signal"
inclusive_bulk="${bulk_root}/inclusive"
canary_root="${RJ_THE112_CANARY_ROOT:-$(rj_recoiljets_bulk_root)/smoke/${family}/${campaign_tag}}"
smoke_canary_root="${canary_root}/transport_smoke"
scan_canary_root="${canary_root}/scan"
confirmation_canary_root="${canary_root}/confirmation"
merge_base="${RJ_THE112_MERGE_BASE:-$(rj_merge_root_for_campaign "$campaign_tag")}" # recorded only; automatic merge is disabled

data_group="${RJ_THE112_DATA_GROUP_SIZE:-7}"
sim_group="${RJ_THE112_SIM_GROUP_SIZE:-1}"
data_memory="${RJ_THE112_DATA_MEMORY:-6000MB}"
sim_memory="${RJ_THE112_SIM_MEMORY:-7000MB}"
smoke_events="${RJ_THE112_SMOKE_EVENTS:-20000}"
smoke_data_runs="${RJ_THE112_SMOKE_DATA_RUNS:-3}"
smoke_data_jobs="${RJ_THE112_SMOKE_DATA_JOBS:-3}"
discovery_events="${RJ_THE112_DISCOVERY_EVENTS:-20000}"
discovery_data_identities="${RJ_THE112_DISCOVERY_DATA_IDENTITIES:-384}"
discovery_data_strata="${RJ_THE112_DISCOVERY_DATA_STRATA:-16}"
discovery_signal_jobs="${RJ_THE112_DISCOVERY_SIGNAL_JOBS:-1008}"
discovery_inclusive_jobs="${RJ_THE112_DISCOVERY_INCLUSIVE_JOBS:-2016}"
allow_existing="${RJ_THE112_ALLOW_EXISTING:-0}"

export RJ_CODEX_CHAT_NAME="${RJ_CODEX_CHAT_NAME:-THE-112 | Canonical AuAu WP80 Triplet}"
# Deliberately no fallback thread id: a Codex-launched DAG without an exact
# task id is a provenance defect and must fail before submission.
export RJ_CODEX_THREAD_ID="${RJ_CODEX_THREAD_ID:-}"

usage() {
  cat <<'EOF'
Usage: submit_the112_auau_sideband_campaign.sh MODE [ARGS]

Read-only / local build modes:
  print
      Print the immutable THE-111/THE-112 contract.
  preflight
      Verify config/model/data-manifest hashes and input counts; no submission.
  build-isolated
      Build the AuAu library in a THE-112-owned hidden prefix; no submission.
  historical-scan
      Rank the exact three THE-100 bounded ROOTs with the common ranker.
  rank-scan-canary DATA_ROOT SIGNAL_ROOT INCLUSIVE_ROOT
      Rank matched THE-112 V2 continuous surfaces; no submission.
  bind-comparison-band MIN_OFFSET MAX_OFFSET RANKING_JSON
      Freeze one blinded bounded lane for dual-fanout comparison while retaining
      the full finite-score complement as the nominal/control definition.
  validate-canary STAMPED_YAML BOUNDED_ROOT COMPLEMENT_ROOT SAMPLE_KIND OUTPUT_JSON
      Audit the sliding-R40 canonical view, sliding-R30 robustness view, and
      ROOT namespaces.
  production-dryrun
      Print full data/signal/inclusive job counts using the frozen config.
  status
      Read-only matching queue and output-count report.

Submitting modes (all one-shot; no automatic merge):
  smoke-canary-submit
      Submit a three-run data + one-job Photon12 + one-job Jet12 transport
      smoke. This output is never eligible for sideband ranking.
  scan-canary-submit
      Submit the deterministic stratified discovery sample with V2 diagnostics
      enabled. This is the minimum canary eligible for sideband ranking.
  confirmation-canary-submit
      Submit the same matched canary with one ranked finite band. Requires
      RJ_THE112_CONFIRM_MIN_OFFSET, RJ_THE112_CONFIRM_MAX_OFFSET, and
      RJ_THE112_RANKING_JSON.
  production-submit-data
  production-submit-signal
  production-submit-inclusive-sample SAMPLE
      Submit one frozen full-production lane. These require
      RJ_THE112_FROZEN_CONFIG_YAML and RJ_THE112_SIDEBAND_FREEZE_JSON.
      The freeze JSON must bind the exact selected offsets, ranking packet,
      PASS confirmation-validation packets, config/model/pairing hashes, and
      canonical_iso_view=isoR40_isSliding.
      SAMPLE is run28_embeddedJet12, Jet20, Jet30, or Jet40.

No mode releases, removes, retries, or merges Condor jobs.
EOF
}

require_file() {
  local path="$1" label="$2"
  [[ -s "$path" ]] || die "Missing ${label}: ${path}"
}

require_sha256() {
  local path="$1" expected="$2" label="$3" observed
  require_file "$path" "$label"
  observed="$(sha256sum "$path" | awk '{print $1}')"
  [[ "$observed" == "$expected" ]] || \
    die "${label} SHA256 mismatch: expected ${expected}, observed ${observed} (${path})"
}

require_submission_provenance() {
  [[ -n "$RJ_CODEX_CHAT_NAME" ]] || die "RJ_CODEX_CHAT_NAME is required for submission provenance."
  [[ -n "$RJ_CODEX_THREAD_ID" ]] || die "RJ_CODEX_THREAD_ID is required for submission provenance."
}

validate_canary_dimensions() {
  python3 - "$smoke_data_runs" "$smoke_data_jobs" "$smoke_events" "$discovery_data_identities" \
      "$discovery_data_strata" "$discovery_signal_jobs" \
      "$discovery_inclusive_jobs" "$discovery_events" <<'PY'
import sys

names = (
    "smoke_data_runs",
    "smoke_data_jobs",
    "smoke_events",
    "discovery_data_identities",
    "discovery_data_strata",
    "discovery_signal_jobs",
    "discovery_inclusive_jobs",
    "discovery_events",
)
try:
    values = dict(zip(names, map(int, sys.argv[1:])))
except ValueError as exc:
    raise SystemExit(f"THE-112 canary dimensions must be integers: {exc}") from exc
if values["smoke_data_runs"] != 3:
    raise SystemExit("THE-112 transport smoke is fixed to exactly three data runs")
if values["smoke_data_jobs"] != 3:
    raise SystemExit("THE-112 transport smoke is fixed to exactly three data jobs")
if values["smoke_events"] <= 0 or values["discovery_events"] <= 0:
    raise SystemExit("THE-112 canary event caps must be positive")
if values["discovery_data_identities"] < 384:
    raise SystemExit("ranking discovery requires at least 384 paired data identities")
if values["discovery_data_strata"] < 16:
    raise SystemExit("ranking discovery requires at least 16 data run-quantile strata")
if values["discovery_data_identities"] % values["discovery_data_strata"]:
    raise SystemExit("discovery data identities must divide evenly across strata")
if values["discovery_signal_jobs"] < 1008 or values["discovery_inclusive_jobs"] < 2016:
    raise SystemExit("ranking discovery SIM lanes are below the frozen minimum sizes")
PY
}

validate_config_contract() {
  local config="$1"
  require_file "$config" "THE-112 config"
  python3 - "$config" "$model_path" "$wp_intercept" "$wp_slope" "$wp_min_et" "$wp_max_et" <<'PY'
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
text = path.read_text()
model, intercept, slope, et_min, et_max = sys.argv[2:]
checks = {
    "exact THE-111 model": model in text,
    "exact WP80 entry": (
        f'auauCentInputBase3x3BDT|centlinear|{intercept}|{slope}|{et_min}|{et_max}|1'
        in text
    ),
    "exact two-row fanout": len(re.findall(r'^\s*- \[newPPG12, auauCentInputBase3x3BDT, ', text, re.M)) == 2,
    "bounded row": "[newPPG12, auauCentInputBase3x3BDT, auauBDTSideband]" in text,
    "complement row": "[newPPG12, auauCentInputBase3x3BDT, auauBDTComplement]" in text,
    "sliding-only reconstructed isolation": bool(
        re.search(r'^isSlidingIso:\s*true\s*$', text, re.M)
        and re.search(r'^isSlidingAndFixed:\s*false\s*$', text, re.M)
        and re.search(r'^fixedGeV:\s*0(?:\.0+)?(?:\s*#.*)?$', text, re.M)
        and re.search(r'^coneR:\s*\[0\.40,\s*0\.30\]\s*$', text, re.M)
        and "fixedIso" not in text
    ),
    "canonical R40 centrality-dependent reco isolation": bool(
        re.search(r'auau_cent_iso_wp:\s*\n\s*-\s*\{aGeV:\s*7\.57,\s*bPerGeV:\s*-0\.0658', text)
        and re.search(r'auau_cent_iso_wp_r40:\s*\n\s*-\s*\{aGeV:\s*7\.57,\s*bPerGeV:\s*-0\.0658', text)
    ),
    "R30 robustness reco isolation": bool(
        re.search(r'auau_cent_iso_wp_r30:\s*\n\s*-\s*\{aGeV:\s*5\.97,\s*bPerGeV:\s*-0\.0507', text)
    ),
    "AuAu zero tower floor": bool(re.search(r'towerMin:\s*0\.0', text)),
}
failed = [name for name, passed in checks.items() if not passed]
if failed:
    raise SystemExit("THE112_CONFIG_CONTRACT_FAIL: " + ", ".join(failed))
print("THE112_CONFIG_CONTRACT_PASS")
PY
}

verify_pairing_inventory() {
  require_sha256 "$pairing_report" "$pairing_report_sha256" "AuAu pairing report"
  local run_count pair_count
  run_count="$(find "$pair_list_dir" -maxdepth 1 -type f -name 'dst_auau_jet_pair-*.list' | wc -l | tr -d ' ')"
  # xargs may split the 2,776 lists across multiple awk invocations when the
  # expanded argument vector exceeds ARG_MAX.  Sum every batch result before
  # comparing it with the frozen paired-row count.
  pair_count="$(find "$pair_list_dir" -maxdepth 1 -type f -name 'dst_auau_jet_pair-*.list' -print0 | xargs -0 awk 'NF >= 2 {n++} END {print n+0}' | awk '{total += $1} END {print total+0}')"
  [[ "$run_count" == "$expected_data_runs" ]] || \
    die "AuAu paired run-list count mismatch: expected ${expected_data_runs}, observed ${run_count}"
  [[ "$pair_count" == "$expected_data_pairs" ]] || \
    die "AuAu paired segment count mismatch: expected ${expected_data_pairs}, observed ${pair_count}"
  say "Pairing inventory PASS: ${run_count} runs / ${pair_count} matched rows"
}

write_manifest() {
  mkdir -p "$evidence_dir"
  {
    printf 'field\tvalue\n'
    printf 'campaign_tag\t%s\n' "$campaign_tag"
    printf 'model_path\t%s\n' "$model_path"
    printf 'model_sha256\t%s\n' "$model_sha256"
    printf 'wp80\tT80(c)=%s+%s*c;15<=ET<35\n' "$wp_intercept" "$wp_slope"
    printf 'pairing_report\t%s\n' "$pairing_report"
    printf 'pairing_report_sha256\t%s\n' "$pairing_report_sha256"
    printf 'data_contract\t%s runs;%s matched DST_JET+DST_JETCALO rows\n' "$expected_data_runs" "$expected_data_pairs"
    printf 'data_output\t%s\n' "$data_bulk"
    printf 'signal_output\t%s\n' "$signal_bulk"
    printf 'inclusive_output\t%s\n' "$inclusive_bulk"
    printf 'automatic_merge\tdisabled\n'
  } > "$manifest"
}

print_contract() {
  cat <<EOF
RECOILJETS_THE112_AUAU_SIDEBAND_CAMPAIGN_V1
campaign_tag=${campaign_tag}
scan_config=${scan_yaml}
isolated_library=${auau_library}
analysis_python=${analysis_python}
model=${model_path}
model_sha256=${model_sha256}
tight_wp80=T80(c)=${wp_intercept}+${wp_slope}*c, ${wp_min_et}<=ET<${wp_max_et} GeV
data_input=${expected_data_runs} runs, ${expected_data_pairs} paired DST_JET+DST_JETCALO rows
pairing_report=${pairing_report}
pairing_report_sha256=${pairing_report_sha256}
signal_samples=run28_embeddedPhoton12,run28_embeddedPhoton20
inclusive_samples=run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40
scan_grid=lower_widths{0.16,0.20,0.24,0.28,0.32,0.40};gaps{0.02,0.03,0.04,0.05,0.06,0.08}
transport_smoke=data(${smoke_data_runs} runs capped at ${smoke_data_jobs} jobs),Photon12(1 job),Jet12(1 job);not ranking eligible
discovery_canary=data(${discovery_data_identities} paired identities across ${discovery_data_strata} run strata),Photon12(${discovery_signal_jobs} jobs),Jet12(${discovery_inclusive_jobs} jobs)
canonical_reco_isolation=centrality-dependent sliding R=0.4, EisoCut=7.57-0.0658*c
robustness_reco_isolation=centrality-dependent sliding R=0.3, EisoCut=5.97-0.0507*c
truth_isolation_label_GeV=4.0
canary_v2_diagnostics=enabled
full_production_v2_diagnostics=disabled
automatic_merge=disabled
bulk_root=${bulk_root}
canary_root=${canary_root}
merge_base_record_only=${merge_base}
RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}
RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID:-UNSET_SUBMISSION_BLOCKED}
EOF
}

preflight() {
  rj_validate_campaign_tag "$campaign_tag" || exit 2
  validate_canary_dimensions
  validate_config_contract "$scan_yaml"
  require_sha256 "$model_path" "$model_sha256" "THE-111 TMVA model"
  verify_pairing_inventory
  [[ -x "$submitter" ]] || die "Missing executable RecoilJets Condor submitter: ${submitter}"
  [[ -x "$analysis_python" ]] || die "Missing THE-112 analysis Python: ${analysis_python}"
  require_file "$ranker" "THE-112 ranker"
  require_file "$canary_validator" "THE-112 canary validator"
  "$analysis_python" -c 'import numpy, uproot' || die "THE-112 analysis Python is missing numpy/uproot."
  grep -F 'RJ_INTERNAL_ISO_VIEWS_OVERRIDE' "$submitter" >/dev/null || \
    die "Submitter does not support the required ordered RJ_INTERNAL_ISO_VIEWS_OVERRIDE; the canonical R40-first contract cannot be proven."
  write_manifest
  say "Preflight PASS (read-only): exact model, WP, pairing manifest, and launcher inputs"
}

build_isolated() {
  [[ -s "${repo_root}/src_AuAu/autogen.sh" ]] || die "Missing src_AuAu/autogen.sh"
  case "$build_root" in
    "${repo_root}/.recoiljets_tmp/the112_sideband_build_"*) ;;
    *) die "Refusing build cleanup outside THE-112-owned hidden path: ${build_root}" ;;
  esac
  rm -rf "$build_dir" "$install_dir"
  mkdir -p "$build_dir" "$install_dir" "$evidence_dir"
  cp -a "${repo_root}/src_AuAu/." "$build_dir/"
  rm -rf "${build_dir}/autom4te.cache" "${build_dir}/.deps"
  rm -f "${build_dir}/config.status" "${build_dir}/config.log" \
        "${build_dir}/Makefile" "${build_dir}/libtool" "${build_dir}/stamp-h1"
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  source /opt/sphenix/core/bin/setup_local.sh "${RJ_THE112_DEP_INSTALL:-/sphenix/u/patsfan753/thesisAnalysis/install}"
  set -u
  (
    cd "$build_dir"
    bash ./autogen.sh --prefix="$install_dir"
    make -j"${RJ_THE112_BUILD_JOBS:-4}"
    make install
  ) 2>&1 | tee "${evidence_dir}/isolated_build.log"
  require_file "$auau_library" "isolated AuAu library"
  sha256sum "$auau_library" | tee "${evidence_dir}/isolated_library.sha256"
}

require_runtime_inputs() {
  local config="$1"
  rj_validate_campaign_tag "$campaign_tag" || exit 2
  validate_config_contract "$config"
  require_sha256 "$model_path" "$model_sha256" "THE-111 TMVA model"
  require_file "$auau_library" "isolated AuAu library (run build-isolated first)"
  [[ -x "$submitter" ]] || die "Missing executable RecoilJets Condor submitter: ${submitter}"
  grep -F 'RJ_INTERNAL_ISO_VIEWS_OVERRIDE' "$submitter" >/dev/null || \
    die "Submitter lacks ordered internal-isolation override support; refusing a run with an ambiguous canonical view."
}

common_env() {
  local config="$1" diagnostics="$2" role="$3"
  [[ "$role" =~ ^[A-Za-z0-9_.-]+$ ]] || die "Unsafe THE-112 submission role: ${role}"
  printf '%s\n' \
    "RJ_CONFIG_YAML=${config}" \
    "RJ_AUAU_LIBRARY_OVERRIDE=${auau_library}" \
    "RJ_INTERNAL_ISO_VIEWS_OVERRIDE=${expected_iso_views}" \
    "RJ_SUBMIT_EXTRA_ENV=RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=1;RJ_AUAU_SHOWER_SHAPE_DIAGNOSTIC_VARIANT=canonical;RJ_AUAU_DUALVIEW_DIAGNOSTICS=${diagnostics};RJ_AUAU_FIG25_CORRELATION_DIAGNOSTICS=0;RJ_THE112_CAMPAIGN_TAG=${campaign_tag};RJ_THE112_SUBMISSION_ROLE=${role}" \
    "RJ_ID_FANOUT_MAX_ROWS=2" \
    "RJ_AUTO_MERGE=0" \
    "RJ_REQUIRE_NON_TINY_OUTPUT=1" \
    "RJ_MIN_OUTPUT_BYTES=50000" \
    "RJ_PROFILE_JOB=1" \
    "RJ_PROFILE_LABEL=${campaign_tag}_${role}" \
    "RJ_JOB_HEARTBEAT_SECONDS=120" \
    "RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0" \
    "RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0" \
    "RJ_CODEX_CHAT_NAME=${RJ_CODEX_CHAT_NAME}" \
    "RJ_CODEX_THREAD_ID=${RJ_CODEX_THREAD_ID}"
}

assert_no_active_duplicate() {
  local role="$1" output_root="$2" sample="${3:-}" queue_json rc=0
  require_no_receipt "$role"
  queue_json="$(mktemp "${TMPDIR:-/tmp}/the112_condor_q.XXXXXX.json")"
  if ! condor_q "${USER:-patsfan753}" -json \
      -attributes ClusterId,ProcId,JobStatus,Environment,Args,RJCampaignTag,RJSubmissionRole \
      > "$queue_json" 2>/dev/null; then
    rm -f "$queue_json"
    die "Cannot audit the live Condor queue; duplicate protection fails closed." 4
  fi
  python3 - "$queue_json" "$campaign_tag" "$role" "$output_root" "$sample" <<'PY' || rc=$?
import json
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
campaign, role, output_root, sample = sys.argv[2:]
payload = json.loads(path.read_text() or "[]")
matches = []
for ad in payload:
    environment = str(ad.get("Environment", ""))
    arguments = str(ad.get("Args", ""))
    campaign_ad = str(ad.get("RJCampaignTag", ""))
    role_ad = str(ad.get("RJSubmissionRole", ""))
    joined = "\n".join((environment, arguments, campaign_ad, role_ad))
    campaign_match = (
        campaign_ad == campaign
        or f"RJ_THE112_CAMPAIGN_TAG={campaign}" in environment
        or campaign in arguments
    )
    role_match = (
        role_ad == role
        or f"RJ_THE112_SUBMISSION_ROLE={role}" in environment
        or role in arguments
    )
    output_match = bool(output_root and output_root in joined)
    sample_match = bool(sample and sample in joined)
    output_scope_match = output_match and (not sample or sample_match)
    if output_scope_match or (campaign_match and role_match) or (campaign_match and sample_match):
        matches.append({
            "ClusterId": ad.get("ClusterId"),
            "ProcId": ad.get("ProcId"),
            "JobStatus": ad.get("JobStatus"),
            "campaign_match": campaign_match,
            "role_match": role_match,
            "output_match": output_scope_match,
            "sample_match": sample_match,
        })
if matches:
    print(json.dumps(matches[:20], indent=2, sort_keys=True), file=sys.stderr)
    raise SystemExit(10)
PY
  rm -f "$queue_json"
  if [[ "$rc" -eq 10 ]]; then
    die "Active Condor jobs match role=${role}, output=${output_root}, sample=${sample:-none}; refusing duplicate submission." 4
  elif [[ "$rc" -ne 0 ]]; then
    die "Condor duplicate-audit parser failed for role=${role}; refusing submission." 4
  fi
}

assert_fresh_path() {
  local path="$1" label="$2"
  if [[ -e "$path" && "$allow_existing" != "1" ]]; then
    die "${label} already exists: ${path}. Use a fresh campaign tag; intentional continuation requires RJ_THE112_ALLOW_EXISTING=1 plus a distinct one-shot receipt."
  fi
}

receipt_path() {
  printf '%s/submitted_%s.receipt\n' "$evidence_dir" "$1"
}

pending_receipt_path() {
  printf '%s/submitted_%s.pending.receipt\n' "$evidence_dir" "$1"
}

require_no_receipt() {
  local role="$1" receipt pending
  receipt="$(receipt_path "$role")"
  pending="$(pending_receipt_path "$role")"
  [[ ! -e "$receipt" ]] || die "Submission receipt already exists for ${role}: ${receipt}. No automatic resubmit is permitted."
  [[ ! -e "$pending" ]] || die "Pending submission receipt already exists for ${role}: ${pending}. Inspect it before any explicit recovery; automatic resubmit is forbidden."
}

reserve_submission_receipt() {
  local role="$1" config="$2" target_root="$3"
  shift 3
  local pending binding freeze_json
  pending="$(pending_receipt_path "$role")"
  binding="${config%.yaml}.binding.json"
  freeze_json="${RJ_THE112_SIDEBAND_FREEZE_JSON:-}"
  mkdir -p "$evidence_dir"
  python3 - "$pending" "$role" "$campaign_tag" "$config" "$model_sha256" \
      "$pairing_report_sha256" "$target_root" "$RJ_CODEX_CHAT_NAME" \
      "$RJ_CODEX_THREAD_ID" "$binding" "$freeze_json" "$@" <<'PY'
import datetime
import hashlib
import os
import pathlib
import sys

(
    pending_raw,
    role,
    campaign,
    config_raw,
    model_sha,
    pairing_sha,
    target_root,
    chat_name,
    thread_id,
    binding_raw,
    freeze_raw,
    *lanes,
) = sys.argv[1:]
for label, value in (("role", role), ("campaign", campaign), ("target_root", target_root), *[("lane", item) for item in lanes]):
    if "\n" in value or "\r" in value:
        raise SystemExit(f"unsafe newline in receipt {label}")
config = pathlib.Path(config_raw)
binding = pathlib.Path(binding_raw)
freeze = pathlib.Path(freeze_raw) if freeze_raw else None
lines = [
    "schema=THE112_SUBMISSION_RECEIPT_V2",
    "state=pending",
    f"role={role}",
    f"reserved_utc={datetime.datetime.now(datetime.timezone.utc).isoformat()}",
    f"campaign_tag={campaign}",
    f"target_root={target_root}",
    f"config={config}",
    f"config_sha256={hashlib.sha256(config.read_bytes()).hexdigest()}",
    f"model_sha256={model_sha}",
    f"pairing_report_sha256={pairing_sha}",
    f"RJ_CODEX_CHAT_NAME={chat_name}",
    f"RJ_CODEX_THREAD_ID={thread_id}",
    "automatic_merge=0",
]
if binding.is_file():
    lines.extend((f"confirmation_binding={binding}", f"confirmation_binding_sha256={hashlib.sha256(binding.read_bytes()).hexdigest()}"))
if freeze and freeze.is_file():
    lines.extend((f"freeze_json={freeze}", f"freeze_json_sha256={hashlib.sha256(freeze.read_bytes()).hexdigest()}"))
lines.extend(f"planned_lane={lane}" for lane in lanes)
fd = os.open(pending_raw, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
with os.fdopen(fd, "w") as stream:
    stream.write("\n".join(lines) + "\n")
    stream.flush()
    os.fsync(stream.fileno())
PY
}

record_submission_lane() {
  local role="$1" lane="$2" target="$3" pending
  pending="$(pending_receipt_path "$role")"
  python3 - "$pending" "$lane" "$target" <<'PY'
import datetime
import os
import sys

path, lane, target = sys.argv[1:]
if any("\n" in value or "\r" in value for value in (lane, target)):
    raise SystemExit("unsafe newline in receipt lane evidence")
line = (
    f"submitted_lane={lane}\t"
    f"utc={datetime.datetime.now(datetime.timezone.utc).isoformat()}\t"
    f"target={target}\n"
).encode()
fd = os.open(path, os.O_WRONLY | os.O_APPEND)
try:
    os.write(fd, line)
    os.fsync(fd)
finally:
    os.close(fd)
PY
}

finalize_submission_receipt() {
  local role="$1" pending final
  pending="$(pending_receipt_path "$role")"
  final="$(receipt_path "$role")"
  python3 - "$pending" "$final" <<'PY'
import datetime
import os
import pathlib
import sys

pending = pathlib.Path(sys.argv[1])
final = pathlib.Path(sys.argv[2])
fd = os.open(pending, os.O_WRONLY | os.O_APPEND)
try:
    os.write(fd, f"state=completed\ncompleted_utc={datetime.datetime.now(datetime.timezone.utc).isoformat()}\n".encode())
    os.fsync(fd)
finally:
    os.close(fd)
os.link(pending, final)
pending.unlink()
PY
}

run_ranker() {
  local schema="$1" data_root="$2" signal_root="$3" inclusive_root="$4" stem="$5"
  local -a schema_args gap_args
  schema_args=()
  gap_args=()
  if [[ "$schema" == "legacy" ]]; then
    # THE-100 has 0.02-wide score bins. Only exact bin edges are legal;
    # -0.03 remains an external benchmark and is not approximated.
    schema_args+=(--legacy-iso-cut 4.0)
    gap_args+=(--gaps 0.02,0.04,0.06,0.08)
  else
    gap_args+=(--gaps 0.02,0.03,0.04,0.05,0.06,0.08)
  fi
  require_file "$ranker" "THE-112 ranker"
  [[ -x "$analysis_python" ]] || die "Missing THE-112 analysis Python: ${analysis_python}"
  require_file "$data_root" "${schema} data ROOT"
  require_file "$signal_root" "${schema} signal ROOT"
  require_file "$inclusive_root" "${schema} inclusive ROOT"
  mkdir -p "$evidence_dir"
  "$analysis_python" "$ranker" \
    --schema "$schema" \
    --data "$data_root" \
    --signal "$signal_root" \
    --inclusive "$inclusive_root" \
    --data-trigger MBD_NS_geq_2_vtx_lt_150 \
    --simulation-trigger SIM \
    --output-json "${evidence_dir}/${stem}.json" \
    --output-csv "${evidence_dir}/${stem}.csv" \
    "${schema_args[@]}" \
    "${gap_args[@]}" \
    --lower-widths 0.16,0.20,0.24,0.28,0.32,0.40 \
    --isolation-gap 0.0
}

historical_scan() {
  run_ranker legacy "$historical_data" "$historical_signal" "$historical_inclusive" \
    "the112_historical_sideband_scan"
}

rank_scan_canary() {
  [[ $# -eq 3 ]] || die "rank-scan-canary requires DATA_ROOT SIGNAL_ROOT INCLUSIVE_ROOT"
  run_ranker v2 "$1" "$2" "$3" "the112_v2_sideband_scan"
}

audit_stamped_yaml() {
  local stamped_yaml="$1" output_json="$2"
  require_file "$stamped_yaml" "stamped THE-112 YAML"
  python3 - "$stamped_yaml" "$output_json" "$expected_iso_views" <<'PY'
import json
import pathlib
import re
import sys

yaml_path = pathlib.Path(sys.argv[1])
output_path = pathlib.Path(sys.argv[2])
expected = sys.argv[3]
text = yaml_path.read_text()
match = re.search(r'(?m)^internal_iso_cone_views:\s*(.+?)\s*$', text)
observed = match.group(1).strip().strip('"\'') if match else None
failures = []
if observed != expected:
    failures.append(
        f"internal_iso_cone_views mismatch: expected {expected!r}, observed {observed!r}"
    )
if observed and "fixedIso" in observed:
    failures.append("fixed reconstructed-isolation namespace is forbidden in THE-112")
if observed and not observed.startswith("isoR40_isSliding:"):
    failures.append("first/canonical isolation view is not sliding R=0.4")
payload = {
    "schema": "THE112_STAMPED_YAML_CANONICAL_VIEW_AUDIT_V1",
    "status": "PASS" if not failures else "FAIL",
    "stamped_yaml": str(yaml_path),
    "expected_internal_iso_cone_views": expected,
    "observed_internal_iso_cone_views": observed,
    "failures": failures,
}
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(json.dumps(payload, indent=2, sort_keys=True))
if failures:
    raise SystemExit(1)
PY
}

validate_canary() {
  [[ $# -eq 5 ]] || die "validate-canary requires STAMPED_YAML BOUNDED_ROOT COMPLEMENT_ROOT SAMPLE_KIND OUTPUT_JSON"
  local stamped_yaml="$1" bounded="$2" complement="$3" sample_kind="$4" output_json="$5"
  case "$sample_kind" in data|signal|inclusive) ;; *) die "SAMPLE_KIND must be data, signal, or inclusive" ;; esac
  require_file "$canary_validator" "THE-112 canary validator"
  [[ -x "$analysis_python" ]] || die "Missing THE-112 analysis Python: ${analysis_python}"
  local stamped_audit="${output_json%.json}.stamped_yaml.json"
  audit_stamped_yaml "$stamped_yaml" "$stamped_audit"
  "$analysis_python" "$canary_validator" \
    --bounded "$bounded" \
    --complement "$complement" \
    --sample-kind "$sample_kind" \
    --json "$output_json"
  python3 - "$output_json" "$stamped_yaml" "$stamped_audit" <<'PY'
import hashlib
import json
import pathlib
import re
import sys

report_path, yaml_path, audit_path = map(pathlib.Path, sys.argv[1:])
report = json.loads(report_path.read_text())
yaml_text = yaml_path.read_text()
def scalar(key: str) -> str:
    match = re.search(rf"(?m)^{re.escape(key)}:\s*([^\s#]+)", yaml_text)
    if not match:
        raise SystemExit(f"stamped confirmation YAML is missing {key}")
    return match.group(1)
ranking_match = re.search(r"(?m)^# ranking_sha256:\s*([0-9a-f]{64})\s*$", yaml_text)
candidate_match = re.search(r"(?m)^# ranking_candidate_key:\s*(\S+)\s*$", yaml_text)
if not ranking_match or not candidate_match:
    raise SystemExit("stamped confirmation YAML is missing its ranking/candidate binding comments")
report["launcher_binding"] = {
    "schema": "THE112_CONFIRMATION_VALIDATION_BINDING_V1",
    "stamped_yaml": str(yaml_path),
    "stamped_yaml_sha256": hashlib.sha256(yaml_path.read_bytes()).hexdigest(),
    "stamped_yaml_audit": str(audit_path),
    "stamped_yaml_audit_sha256": hashlib.sha256(audit_path.read_bytes()).hexdigest(),
    "ranking_sha256": ranking_match.group(1),
    "ranking_candidate_key": candidate_match.group(1),
    "selected_offsets": {
        "min_offset": scalar("auau_nontight_bdt_relative_min_offset"),
        "max_offset": scalar("auau_nontight_bdt_relative_max_offset"),
    },
}
tmp = report_path.with_suffix(report_path.suffix + ".tmp")
tmp.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
tmp.replace(report_path)
PY
}

submit_smoke_lane() {
  local config="$1" output_base="$2" dataset="$3" sample="$4" memory="$5" role="$6"
  local -a env_args extra_args
  mapfile -t env_args < <(common_env "$config" 1 "$role")
  extra_args=()
  [[ "$dataset" != "isSimEmbeddedInclusive" ]] || extra_args+=("RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1")
  env "${env_args[@]}" "${extra_args[@]}" \
    "RJ_REQUEST_MEMORY=${memory}" \
    "RJ_SMOKE_OUTPUT_BASE=${output_base}" \
    "RJ_SMOKE_SIM_NEVENTS=${smoke_events}" \
    "$submitter" "$dataset" condorDoAllSmoke groupSize 1 maxJobs 1 "SAMPLE=${sample}"
}

submit_transport_smoke() {
  local role="$1" config="$2" output_root="$3"
  require_submission_provenance
  validate_canary_dimensions
  require_runtime_inputs "$config"
  assert_no_active_duplicate "$role" "$output_root"
  assert_fresh_path "$output_root" "${role} output root"
  verify_pairing_inventory
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/${role}_submission.log") 2>&1
  local -a env_args
  mapfile -t env_args < <(common_env "$config" 1 "$role")
  reserve_submission_receipt "$role" "$config" "$output_root" \
    "data:${output_root}/data" "signal:Photon12:${output_root}/signal" \
    "inclusive:Jet12:${output_root}/inclusive"

  say "Submitting non-ranking transport smoke: ${smoke_data_runs} data runs capped at ${smoke_data_jobs} jobs, one Photon12 job, one Jet12 job"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${data_memory}" \
    "RJ_SMOKE_OUTPUT_BASE=${output_root}/data" \
    "RJ_SMOKE_DATA_RUNS=${smoke_data_runs}" \
    "RJ_SMOKE_DATA_MAX_JOBS=${smoke_data_jobs}" \
    "RJ_SMOKE_DATA_NEVENTS=${smoke_events}" \
    "$submitter" isAuAu condor smokeTest groupSize 1
  record_submission_lane "$role" data "${output_root}/data"
  submit_smoke_lane "$config" "${output_root}/signal" isSimEmbedded run28_embeddedPhoton12 "$sim_memory" "$role"
  record_submission_lane "$role" signal_Photon12 "${output_root}/signal"
  submit_smoke_lane "$config" "${output_root}/inclusive" isSimEmbeddedInclusive run28_embeddedJet12 "$sim_memory" "$role"
  record_submission_lane "$role" inclusive_Jet12 "${output_root}/inclusive"
  finalize_submission_receipt "$role"
}

prepare_discovery_inputs() {
  local selection_root="${evidence_dir}/discovery_inputs"
  python3 - "$pair_list_dir" "${repo_root}/simListFiles" "$selection_root" \
      "$discovery_data_identities" "$discovery_data_strata" \
      "$discovery_signal_jobs" "$discovery_inclusive_jobs" <<'PY'
import hashlib
import json
import math
import os
import pathlib
import re
import shutil
import sys
import tempfile

pair_root = pathlib.Path(sys.argv[1])
sim_root = pathlib.Path(sys.argv[2])
target = pathlib.Path(sys.argv[3])
data_count, strata, photon_count, jet_count = map(int, sys.argv[4:])
if min(data_count, strata, photon_count, jet_count) <= 0 or data_count % strata:
    raise SystemExit("invalid deterministic discovery dimensions")

def sha(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def spaced_indices(total: int, count: int) -> list[int]:
    if count > total:
        raise SystemExit(f"cannot select {count} deterministic rows from {total}")
    result = [min(total - 1, math.floor((index + 0.5) * total / count)) for index in range(count)]
    if len(set(result)) != count:
        raise SystemExit("deterministic index construction produced duplicates")
    return result

parent = target.parent
parent.mkdir(parents=True, exist_ok=True)
tmp = pathlib.Path(tempfile.mkdtemp(prefix=".the112_discovery_inputs.", dir=parent))
try:
    data_lists = tmp / "data_pair_lists"
    data_lists.mkdir()
    pair_files = sorted(
        pair_root.glob("dst_auau_jet_pair-*.list"),
        key=lambda item: int(re.search(r"(\d+)\.list$", item.name).group(1)),
    )
    if not pair_files:
        raise SystemExit(f"no canonical AuAu pair lists under {pair_root}")
    per_stratum = data_count // strata
    selected_runs = []
    data_manifest = []
    for stratum in range(strata):
        begin = math.floor(stratum * len(pair_files) / strata)
        end = math.floor((stratum + 1) * len(pair_files) / strata)
        block = pair_files[begin:end]
        for local_index in spaced_indices(len(block), per_stratum):
            source = block[local_index]
            valid_rows = [
                (line_number, line)
                for line_number, line in enumerate(source.read_text().splitlines(), start=1)
                if len(line.split()) >= 2
            ]
            if not valid_rows:
                raise SystemExit(f"no paired rows in {source}")
            row_number, row = valid_rows[len(valid_rows) // 2]
            output = data_lists / source.name
            output.write_text(row + "\n")
            run = re.search(r"(\d+)\.list$", source.name).group(1)
            selected_runs.append(run)
            data_manifest.append({
                "stratum": stratum,
                "run": run,
                "source": str(source),
                "source_sha256": sha(source),
                "source_row": row_number,
                "selected_row_sha256": hashlib.sha256((row + "\n").encode()).hexdigest(),
            })
    (tmp / "data_runs.list").write_text("\n".join(selected_runs) + "\n")

    sim_destination = tmp / "simListFiles"
    sim_manifest = {}
    required_lists = (
        "DST_CALO_CLUSTER.matched.list",
        "G4Hits.matched.list",
        "DST_JETS.matched.list",
        "DST_GLOBAL.matched.list",
        "DST_MBD_EPD.matched.list",
    )
    for sample, count in (
        ("run28_embeddedPhoton12", photon_count),
        ("run28_embeddedJet12", jet_count),
    ):
        source_dir = sim_root / sample
        destination_dir = sim_destination / sample
        destination_dir.mkdir(parents=True)
        lines_by_name = {}
        for name in required_lists:
            source = source_dir / name
            if not source.is_file():
                raise SystemExit(f"missing discovery SIM list: {source}")
            lines_by_name[name] = source.read_text().splitlines()
        lengths = {len(lines) for lines in lines_by_name.values()}
        if len(lengths) != 1:
            raise SystemExit(f"unaligned matched SIM lists for {sample}: {sorted(lengths)}")
        total = lengths.pop()
        indices = spaced_indices(total, count)
        outputs = {}
        for name, lines in lines_by_name.items():
            output = destination_dir / name
            output.write_text("\n".join(lines[index] for index in indices) + "\n")
            outputs[name] = {"source": str(source_dir / name), "source_sha256": sha(source_dir / name), "selected_sha256": sha(output)}
        sim_manifest[sample] = {"source_rows": total, "selected_rows": count, "indices_sha256": hashlib.sha256(json.dumps(indices).encode()).hexdigest(), "lists": outputs}

    contract = {
        "schema": "THE112_DETERMINISTIC_STRATIFIED_DISCOVERY_INPUTS_V1",
        "status": "READY",
        "selection": {
            "data_paired_identities": data_count,
            "data_run_quantile_strata": strata,
            "signal_Photon12_rows": photon_count,
            "inclusive_Jet12_rows": jet_count,
            "rule": "deterministic midpoint quantiles; one midpoint pair row per selected data run",
        },
        "data": data_manifest,
        "simulation": sim_manifest,
    }
    (tmp / "contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")

    def tree_hashes(root: pathlib.Path) -> dict[str, str]:
        return {str(path.relative_to(root)): sha(path) for path in sorted(root.rglob("*")) if path.is_file()}

    if target.exists():
        if tree_hashes(target) != tree_hashes(tmp):
            raise SystemExit(f"existing deterministic discovery input packet differs: {target}")
        shutil.rmtree(tmp)
    else:
        os.replace(tmp, target)
finally:
    if tmp.exists():
        shutil.rmtree(tmp)
PY
  printf '%s\n' "$selection_root"
}

submit_discovery_sim_lane() {
  local config="$1" selection_root="$2" output_base="$3" dataset="$4" sample="$5" memory="$6" role="$7" jobs="$8"
  local -a env_args extra_args
  mapfile -t env_args < <(common_env "$config" 1 "$role")
  extra_args=()
  [[ "$dataset" != "isSimEmbeddedInclusive" ]] || extra_args+=("RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1")
  env "${env_args[@]}" "${extra_args[@]}" \
    "RJ_REQUEST_MEMORY=${memory}" \
    "RJ_SIM_ROOT_OVERRIDE=${selection_root}/simListFiles" \
    "RJ_DIRECT_NEVENTS=${discovery_events}" \
    "RJ_SIMEMBED_DEST_BASE=${output_base}" \
    "RJ_SIMEMBEDINCLUSIVE_DEST_BASE=${output_base}" \
    "$submitter" "$dataset" condorDoAll groupSize 1 maxJobs "$jobs" "SAMPLE=${sample}"
}

submit_discovery_canary() {
  local role="$1" config="$2" output_root="$3" selection_root
  require_submission_provenance
  validate_canary_dimensions
  require_runtime_inputs "$config"
  assert_no_active_duplicate "$role" "$output_root"
  assert_fresh_path "$output_root" "${role} output root"
  verify_pairing_inventory
  selection_root="$(prepare_discovery_inputs)"
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/${role}_submission.log") 2>&1
  local -a env_args
  mapfile -t env_args < <(common_env "$config" 1 "$role")
  reserve_submission_receipt "$role" "$config" "$output_root" \
    "data:${discovery_data_identities}:${output_root}/data" \
    "signal:Photon12:${discovery_signal_jobs}:${output_root}/signal" \
    "inclusive:Jet12:${discovery_inclusive_jobs}:${output_root}/inclusive"

  say "Submitting deterministic discovery canary: ${discovery_data_identities} paired data identities / ${discovery_data_strata} run strata, ${discovery_signal_jobs} Photon12 jobs, ${discovery_inclusive_jobs} Jet12 jobs"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${data_memory}" \
    "RJ_GOLDEN_OVERRIDE=${selection_root}/data_runs.list" \
    "RJ_LIST_DIR_OVERRIDE=${selection_root}/data_pair_lists" \
    "RJ_LIST_PREFIX_OVERRIDE=dst_auau_jet_pair" \
    "RJ_DIRECT_NEVENTS=${discovery_events}" \
    "RJ_DIRECT_MAX_JOBS=${discovery_data_identities}" \
    "RJ_DEST_BASE_OVERRIDE=${output_root}/data" \
    "$submitter" isAuAu condor all groupSize 1
  record_submission_lane "$role" data "${output_root}/data"
  submit_discovery_sim_lane "$config" "$selection_root" "${output_root}/signal" isSimEmbedded run28_embeddedPhoton12 "$sim_memory" "$role" "$discovery_signal_jobs"
  record_submission_lane "$role" signal_Photon12 "${output_root}/signal"
  submit_discovery_sim_lane "$config" "$selection_root" "${output_root}/inclusive" isSimEmbeddedInclusive run28_embeddedJet12 "$sim_memory" "$role" "$discovery_inclusive_jobs"
  record_submission_lane "$role" inclusive_Jet12 "${output_root}/inclusive"
  finalize_submission_receipt "$role"
}

generate_confirmation_config() {
  local min_offset="${RJ_THE112_CONFIRM_MIN_OFFSET:-}"
  local max_offset="${RJ_THE112_CONFIRM_MAX_OFFSET:-}"
  local ranking_json="${RJ_THE112_RANKING_JSON:-}"
  [[ -n "$min_offset" && -n "$max_offset" ]] || \
    die "Confirmation requires RJ_THE112_CONFIRM_MIN_OFFSET and RJ_THE112_CONFIRM_MAX_OFFSET."
  require_file "$ranking_json" "ranked sideband JSON"
  local slug generated binding
  slug="$(printf '%s_%s' "$min_offset" "$max_offset" | tr -- '-.' 'mp')"
  generated="${evidence_dir}/generated/analysis_config_the112_confirmation_${slug}.yaml"
  binding="${generated%.yaml}.binding.json"
  mkdir -p "$(dirname "$generated")"
  python3 - "$scan_yaml" "$generated" "$binding" "$min_offset" "$max_offset" "$ranking_json" <<'PY'
import hashlib
import json
import pathlib
import re
import sys
from decimal import Decimal, InvalidOperation

source = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])
binding_path = pathlib.Path(sys.argv[3])
lo_raw = sys.argv[4]
hi_raw = sys.argv[5]
ranking = pathlib.Path(sys.argv[6])
try:
    lo = Decimal(lo_raw)
    hi = Decimal(hi_raw)
except InvalidOperation as exc:
    raise SystemExit(f"invalid confirmation offsets: {lo_raw}, {hi_raw}") from exc
if not (lo < hi < Decimal("0")):
    raise SystemExit(
        f"confirmation offsets must satisfy min < max < 0; received {lo}, {hi}"
    )
ranking_payload = json.loads(ranking.read_text())
if ranking_payload.get("schema") != "THE112_AUAU_SIDEBAND_RANKING_V1":
    raise SystemExit("confirmation requires the machine-readable THE-112 ranking schema")
if ranking_payload.get("input_schema") != "v2":
    raise SystemExit("confirmation requires a corrected V2-surface ranking packet")
if ranking_payload.get("status") not in {"ranked", "PASS", "READY"}:
    raise SystemExit("confirmation ranking packet is not in a ranked/pass/ready state")
eligible = []
for candidate in ranking_payload.get("candidates", []):
    window = candidate.get("window", {})
    try:
        candidate_lo = Decimal(str(window.get("lower_offset")))
        candidate_hi = Decimal(str(window.get("upper_offset")))
    except InvalidOperation:
        continue
    if candidate_lo == lo and candidate_hi == hi:
        eligible.append(candidate)
if len(eligible) != 1:
    raise SystemExit(
        f"confirmation offsets [{lo}, {hi}] matched {len(eligible)} ranking candidates; expected exactly one"
    )
candidate = eligible[0]
if candidate.get("eligible_for_nominal_selection") is not True:
    raise SystemExit("requested confirmation candidate is not eligible_for_nominal_selection")
if candidate.get("passes_gates") is not True or candidate.get("gate_failures"):
    raise SystemExit("requested confirmation candidate did not pass every ranking gate")
text = source.read_text()
text, n_lo = re.subn(
    r'(?m)^auau_nontight_bdt_relative_min_offset:\s*[^\n]+$',
    f'auau_nontight_bdt_relative_min_offset: {lo_raw}',
    text,
)
text, n_hi = re.subn(
    r'(?m)^auau_nontight_bdt_relative_max_offset:\s*[^\n]+$',
    f'auau_nontight_bdt_relative_max_offset: {hi_raw}',
    text,
)
if (n_lo, n_hi) != (1, 1):
    raise SystemExit(f"failed to replace exactly one offset pair: min={n_lo} max={n_hi}")
ranking_sha = hashlib.sha256(ranking.read_bytes()).hexdigest()
header = (
    "# GENERATED THE-112 CONFIRMATION CONFIG; do not hand edit.\n"
    f"# ranking_json: {ranking}\n"
    f"# ranking_sha256: {ranking_sha}\n"
    f"# ranking_candidate_rank: {candidate.get('rank')}\n"
    f"# ranking_candidate_key: {candidate.get('window', {}).get('key')}\n"
    f"# frozen_candidate_offsets: [{lo_raw}, {hi_raw}]\n"
)
config_bytes = (header + text).encode()
tmp_output = output.with_suffix(output.suffix + ".tmp")
tmp_output.write_bytes(config_bytes)
tmp_output.replace(output)
binding_payload = {
    "schema": "THE112_CONFIRMATION_CANDIDATE_BINDING_V1",
    "status": "BOUND",
    "config": str(output),
    "config_sha256": hashlib.sha256(config_bytes).hexdigest(),
    "ranking_json": str(ranking),
    "ranking_sha256": ranking_sha,
    "ranking_schema": ranking_payload.get("schema"),
    "ranking_input_schema": ranking_payload.get("input_schema"),
    "ranking_status": ranking_payload.get("status"),
    "candidate": {
        "rank": candidate.get("rank"),
        "key": candidate.get("window", {}).get("key"),
        "lower_offset": str(lo),
        "upper_offset": str(hi),
        "passes_gates": True,
        "eligible_for_nominal_selection": True,
    },
}

tmp_binding = binding_path.with_suffix(binding_path.suffix + ".tmp")
tmp_binding.write_text(json.dumps(binding_payload, indent=2, sort_keys=True) + "\n")
tmp_binding.replace(binding_path)
PY
  validate_config_contract "$generated" >&2
  require_file "$binding" "confirmation candidate binding" >&2
  printf '%s\n' "$generated"
}

bind_comparison_band() {
  [[ $# -eq 3 ]] || die "bind-comparison-band requires MIN_OFFSET MAX_OFFSET RANKING_JSON"
  local min_offset="$1" max_offset="$2" ranking_json="$3"
  require_file "$ranking_json" "ranked sideband JSON"
  local slug generated freeze_json
  slug="$(printf '%s_%s' "$min_offset" "$max_offset" | tr -- '-.' 'mp')"
  generated="${evidence_dir}/generated/analysis_config_the112_comparison_${slug}.yaml"
  freeze_json="${generated%.yaml}.freeze.json"
  mkdir -p "$(dirname "$generated")"
  python3 - "$scan_yaml" "$generated" "$freeze_json" "$min_offset" "$max_offset" \
    "$ranking_json" "$model_sha256" "$pairing_report_sha256" <<'PY'
import hashlib
import json
import pathlib
import re
import sys
from decimal import Decimal, InvalidOperation

source = pathlib.Path(sys.argv[1])
output = pathlib.Path(sys.argv[2])
freeze_path = pathlib.Path(sys.argv[3])
lo_raw, hi_raw = sys.argv[4:6]
ranking = pathlib.Path(sys.argv[6])
model_sha, pairing_sha = sys.argv[7:9]
try:
    lo = Decimal(lo_raw)
    hi = Decimal(hi_raw)
except InvalidOperation as exc:
    raise SystemExit(f"invalid comparison offsets: {lo_raw}, {hi_raw}") from exc
if not (lo < hi < Decimal("0")):
    raise SystemExit(f"comparison offsets must satisfy min < max < 0; received {lo}, {hi}")
ranking_payload = json.loads(ranking.read_text())
if ranking_payload.get("schema") != "THE112_AUAU_SIDEBAND_RANKING_V1":
    raise SystemExit("comparison binding requires the THE-112 ranking schema")
if ranking_payload.get("input_schema") != "v2" or ranking_payload.get("status") != "ranked":
    raise SystemExit("comparison binding requires a completed corrected-V2 ranking packet")
blinding = ranking_payload.get("blinding_contract", {})
if any(
    blinding.get(key) is not False
    for key in (
        "real_data_recoil_histograms_read",
        "real_data_xjgamma_histograms_read",
        "unfolded_results_read",
    )
):
    raise SystemExit("comparison band was not selected under the frozen blinding contract")
matches = []
for candidate in ranking_payload.get("candidates", []):
    window = candidate.get("window", {})
    try:
        offsets = (
            Decimal(str(window.get("lower_offset"))),
            Decimal(str(window.get("upper_offset"))),
        )
    except InvalidOperation:
        continue
    if offsets == (lo, hi):
        matches.append(candidate)
if len(matches) != 1:
    raise SystemExit(f"comparison offsets [{lo}, {hi}] matched {len(matches)} candidates; expected one")
candidate = matches[0]
text = source.read_text()
text, n_lo = re.subn(
    r'(?m)^auau_nontight_bdt_relative_min_offset:\s*[^\n]+$',
    f'auau_nontight_bdt_relative_min_offset: {lo_raw}',
    text,
)
text, n_hi = re.subn(
    r'(?m)^auau_nontight_bdt_relative_max_offset:\s*[^\n]+$',
    f'auau_nontight_bdt_relative_max_offset: {hi_raw}',
    text,
)
if (n_lo, n_hi) != (1, 1):
    raise SystemExit(f"failed to replace exactly one offset pair: min={n_lo} max={n_hi}")
ranking_sha = hashlib.sha256(ranking.read_bytes()).hexdigest()
header = (
    "# GENERATED THE-112 BLINDED DUAL-FANOUT COMPARISON CONFIG; do not hand edit.\n"
    "# The full finite-score complement remains the nominal/control lane.\n"
    f"# ranking_json: {ranking}\n"
    f"# ranking_sha256: {ranking_sha}\n"
    f"# comparison_candidate_rank: {candidate.get('rank')}\n"
    f"# comparison_candidate_key: {candidate.get('window', {}).get('key')}\n"
    f"# frozen_comparison_offsets: [{lo_raw}, {hi_raw}]\n"
)
config_bytes = (header + text).encode()
tmp_output = output.with_suffix(output.suffix + ".tmp")
tmp_output.write_bytes(config_bytes)
tmp_output.replace(output)
payload = {
    "schema": "THE112_AUAU_SIDEBAND_PRODUCTION_FREEZE_V1",
    "status": "FROZEN_DUAL_FANOUT_COMPARISON",
    "production_permission": "dual_fanout_comparison_only",
    "full_complement_retained": True,
    "nominal_control_lane": "finite_score_complement",
    "bounded_lane_role": "blinded_comparison_pending_full_stat_validation",
    "selection_frozen_before_recoil_result_inspection": True,
    "canonical_iso_view": "isoR40_isSliding",
    "robustness_iso_view": "isoR30_isSliding",
    "config": str(output),
    "config_sha256": hashlib.sha256(config_bytes).hexdigest(),
    "model_sha256": model_sha,
    "pairing_report_sha256": pairing_sha,
    "selected_offsets": {"min_offset": str(lo), "max_offset": str(hi)},
    "selected_candidate_key": candidate.get("window", {}).get("key"),
    "ranking_packet": {"path": str(ranking), "sha256": ranking_sha},
    "ranking_candidate": {
        "rank": candidate.get("rank"),
        "key": candidate.get("window", {}).get("key"),
        "passes_all_nominal_gates": candidate.get("passes_gates"),
        "gate_failures": candidate.get("gate_failures", []),
        "summary": candidate.get("summary", {}),
    },
    "discovery_limitations": [
        "discovery simulation contains Photon12 and Jet12 only",
        "upper photon-pT bins lack inclusive-background support",
        "bounded lane is not promoted as canonical by this packet",
        "full-stat validation must occur before any bounded-lane promotion",
    ],
}
tmp_freeze = freeze_path.with_suffix(freeze_path.suffix + ".tmp")
tmp_freeze.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
tmp_freeze.replace(freeze_path)
print(json.dumps({"config": str(output), "freeze_json": str(freeze_path)}, indent=2))
PY
  validate_config_contract "$generated" >&2
  require_file "$freeze_json" "dual-fanout comparison freeze JSON" >&2
}

require_frozen_production_contract() {
  local config="${RJ_THE112_FROZEN_CONFIG_YAML:-}"
  local freeze_json="${RJ_THE112_SIDEBAND_FREEZE_JSON:-}"
  [[ -n "$config" ]] || die "Production requires RJ_THE112_FROZEN_CONFIG_YAML."
  require_file "$freeze_json" "sideband freeze JSON"
  validate_config_contract "$config" >&2
  python3 - "$freeze_json" "$config" "$model_sha256" "$pairing_report_sha256" <<'PY'
import hashlib
import json
import pathlib
import re
import sys
from decimal import Decimal, InvalidOperation

freeze_path = pathlib.Path(sys.argv[1])
config_path = pathlib.Path(sys.argv[2])
model_sha = sys.argv[3]
pairing_sha = sys.argv[4]
payload = json.loads(freeze_path.read_text())
if payload.get("schema") != "THE112_AUAU_SIDEBAND_PRODUCTION_FREEZE_V1":
    raise SystemExit("sideband freeze JSON has the wrong or missing schema")
expected = {
    "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
    "model_sha256": model_sha,
    "pairing_report_sha256": pairing_sha,
    "canonical_iso_view": "isoR40_isSliding",
}
status = payload.get("status")
if status not in {"FROZEN", "FROZEN_DUAL_FANOUT_COMPARISON"}:
    raise SystemExit("sideband freeze JSON has an unsupported status")
for key, value in expected.items():
    if payload.get(key) != value:
        raise SystemExit(
            f"sideband freeze JSON mismatch for {key}: expected {value!r}, "
            f"observed {payload.get(key)!r}"
        )

config_text = config_path.read_text()
def config_decimal(key: str) -> Decimal:
    match = re.search(rf"(?m)^{re.escape(key)}:\s*([^\s#]+)", config_text)
    if not match:
        raise SystemExit(f"frozen config is missing {key}")
    try:
        return Decimal(match.group(1))
    except InvalidOperation as exc:
        raise SystemExit(f"frozen config has invalid {key}: {match.group(1)!r}") from exc

config_lo = config_decimal("auau_nontight_bdt_relative_min_offset")
config_hi = config_decimal("auau_nontight_bdt_relative_max_offset")
selected = payload.get("selected_offsets")
if not isinstance(selected, dict):
    raise SystemExit("sideband freeze JSON must contain selected_offsets")
try:
    selected_lo = Decimal(str(selected["min_offset"]))
    selected_hi = Decimal(str(selected["max_offset"]))
except (KeyError, InvalidOperation) as exc:
    raise SystemExit("sideband freeze selected_offsets are missing or invalid") from exc
if (selected_lo, selected_hi) != (config_lo, config_hi):
    raise SystemExit(
        "frozen config offsets do not match selected_offsets: "
        f"config=({config_lo},{config_hi}) freeze=({selected_lo},{selected_hi})"
    )
if not (selected_lo < selected_hi < Decimal("0")):
    raise SystemExit("frozen selected offsets must satisfy min < max < 0")

def load_bound_packet(spec: object, label: str) -> tuple[pathlib.Path, dict]:
    if not isinstance(spec, dict) or set(("path", "sha256")) - set(spec):
        raise SystemExit(f"{label} must bind path and sha256")
    packet_path = pathlib.Path(str(spec["path"]))
    if not packet_path.is_file():
        raise SystemExit(f"{label} path is missing: {packet_path}")
    observed_sha = hashlib.sha256(packet_path.read_bytes()).hexdigest()
    if observed_sha != spec["sha256"]:
        raise SystemExit(
            f"{label} SHA256 mismatch: expected {spec['sha256']}, observed {observed_sha}"
        )
    return packet_path, json.loads(packet_path.read_text())

_, ranking = load_bound_packet(payload.get("ranking_packet"), "ranking_packet")
ranking_sha = payload["ranking_packet"]["sha256"]
if ranking.get("schema") != "THE112_AUAU_SIDEBAND_RANKING_V1":
    raise SystemExit("ranking_packet has the wrong schema")
if ranking.get("input_schema") != "v2" or ranking.get("status") not in {"ranked", "PASS", "READY"}:
    raise SystemExit("ranking_packet is not an eligible corrected-V2 ranking")
matching = []
for candidate in ranking.get("candidates", []):
    window = candidate.get("window", {})
    try:
        offsets = (
            Decimal(str(window.get("lower_offset"))),
            Decimal(str(window.get("upper_offset"))),
        )
    except InvalidOperation:
        continue
    if offsets == (selected_lo, selected_hi):
        matching.append(candidate)
if len(matching) != 1:
    raise SystemExit(f"ranking_packet contains {len(matching)} exact selected candidates; expected one")
candidate = matching[0]
candidate_key = candidate.get("window", {}).get("key")
if status == "FROZEN":
    if (
        candidate.get("eligible_for_nominal_selection") is not True
        or candidate.get("passes_gates") is not True
        or candidate.get("gate_failures")
    ):
        raise SystemExit("frozen candidate is not an eligible gate-passing ranking candidate")
    validation_specs = payload.get("confirmation_validation_packets")
    if not isinstance(validation_specs, dict):
        raise SystemExit("freeze JSON must bind confirmation_validation_packets")
    for sample_kind in ("data", "signal", "inclusive"):
        _, report = load_bound_packet(
            validation_specs.get(sample_kind),
            f"confirmation_validation_packets.{sample_kind}",
        )
        if report.get("schema") != "THE112_AUAU_SIDEBAND_CANARY_VALIDATION_V1":
            raise SystemExit(f"{sample_kind} confirmation packet has the wrong schema")
        if report.get("status") != "PASS" or report.get("sample_kind") != sample_kind:
            raise SystemExit(f"{sample_kind} confirmation packet is not an exact PASS for its lane")
        if report.get("failures"):
            raise SystemExit(f"{sample_kind} confirmation packet contains failures")
        launcher_binding = report.get("launcher_binding")
        if not isinstance(launcher_binding, dict) or launcher_binding.get("schema") != "THE112_CONFIRMATION_VALIDATION_BINDING_V1":
            raise SystemExit(f"{sample_kind} confirmation packet lacks its launcher binding")
        if launcher_binding.get("ranking_sha256") != ranking_sha:
            raise SystemExit(f"{sample_kind} confirmation packet was not validated against the frozen ranking packet")
        if launcher_binding.get("ranking_candidate_key") != candidate_key:
            raise SystemExit(f"{sample_kind} confirmation packet candidate key does not match the frozen selection")
        bound_offsets = launcher_binding.get("selected_offsets", {})
        try:
            bound_pair = (
                Decimal(str(bound_offsets["min_offset"])),
                Decimal(str(bound_offsets["max_offset"])),
            )
        except (KeyError, InvalidOperation) as exc:
            raise SystemExit(f"{sample_kind} confirmation packet has invalid bound offsets") from exc
        if bound_pair != (selected_lo, selected_hi):
            raise SystemExit(f"{sample_kind} confirmation packet offsets do not match the frozen selection")
        audit_path = pathlib.Path(str(launcher_binding.get("stamped_yaml_audit", "")))
        if not audit_path.is_file():
            raise SystemExit(f"{sample_kind} stamped-YAML audit packet is missing")
        if hashlib.sha256(audit_path.read_bytes()).hexdigest() != launcher_binding.get("stamped_yaml_audit_sha256"):
            raise SystemExit(f"{sample_kind} stamped-YAML audit SHA256 mismatch")
        audit = json.loads(audit_path.read_text())
        if audit.get("status") != "PASS":
            raise SystemExit(f"{sample_kind} stamped-YAML audit is not PASS")
else:
    required_comparison = {
        "production_permission": "dual_fanout_comparison_only",
        "full_complement_retained": True,
        "nominal_control_lane": "finite_score_complement",
        "bounded_lane_role": "blinded_comparison_pending_full_stat_validation",
        "selection_frozen_before_recoil_result_inspection": True,
    }
    for key, value in required_comparison.items():
        if payload.get(key) != value:
            raise SystemExit(
                f"dual-fanout comparison freeze mismatch for {key}: "
                f"expected {value!r}, observed {payload.get(key)!r}"
            )
    if payload.get("selected_candidate_key") != candidate_key:
        raise SystemExit("comparison freeze candidate key does not match ranking packet")
print("THE112_SIDEBAND_FREEZE_CONTRACT_PASS", file=sys.stderr)
PY
  printf '%s\n' "$config"
}

production_counts() {
  local config
  config="$(require_frozen_production_contract)"
  require_runtime_inputs "$config" >&2
  verify_pairing_inventory
  local -a env_args
  mapfile -t env_args < <(common_env "$config" 0 "production_dryrun")
  say "Full data count (expected source contract: ${expected_data_runs} runs / ${expected_data_pairs} pairs)"
  env "${env_args[@]}" "RJ_REQUEST_MEMORY=${data_memory}" \
    "$submitter" isAuAu CHECKJOBS groupSize "$data_group"
  say "Full Photon12+20 count"
  env "${env_args[@]}" "RJ_REQUEST_MEMORY=${sim_memory}" \
    "$submitter" isSimEmbedded CHECKJOBS groupSize "$sim_group"
  say "Full Jet12+20+30+40 count"
  env "${env_args[@]}" "RJ_REQUEST_MEMORY=${sim_memory}" \
    "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1" \
    "$submitter" isSimEmbeddedInclusive CHECKJOBS groupSize "$sim_group"
}

prepare_production_submit() {
  local role="$1" target_root="${2:-}" sample="${3:-}" config
  require_submission_provenance
  config="$(require_frozen_production_contract)"
  require_runtime_inputs "$config" >&2
  verify_pairing_inventory >&2
  require_no_receipt "$role"
  [[ -z "$target_root" ]] || assert_no_active_duplicate "$role" "$target_root" "$sample"
  write_manifest
  printf '%s\n' "$config"
}

submit_production_data() {
  local config role="production_data"
  config="$(prepare_production_submit "$role" "$data_bulk")"
  assert_fresh_path "$data_bulk" "production data output root"
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/${role}_submission.log") 2>&1
  local -a env_args
  mapfile -t env_args < <(common_env "$config" 0 "$role")
  reserve_submission_receipt "$role" "$config" "$data_bulk" "data:${data_bulk}"
  say "Submitting full ${expected_data_runs}-run / ${expected_data_pairs}-pair AuAu data lane"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${data_memory}" \
    "RJ_DEST_BASE_OVERRIDE=${data_bulk}" \
    "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}/data" \
    "$submitter" isAuAu condor all groupSize "$data_group"
  record_submission_lane "$role" data "$data_bulk"
  finalize_submission_receipt "$role"
}

submit_production_signal() {
  local config role="production_signal"
  config="$(prepare_production_submit "$role" "$signal_bulk")"
  assert_fresh_path "$signal_bulk" "production signal output root"
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/${role}_submission.log") 2>&1
  local -a env_args
  mapfile -t env_args < <(common_env "$config" 0 "$role")
  reserve_submission_receipt "$role" "$config" "$signal_bulk" "signal:Photon12+20:${signal_bulk}"
  say "Submitting full embedded Photon12+20 signal lane"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${sim_memory}" \
    "RJ_SIMEMBED_DEST_BASE=${signal_bulk}" \
    "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}/signal" \
    "$submitter" isSimEmbedded condorDoAll groupSize "$sim_group"
  record_submission_lane "$role" signal_Photon12+20 "$signal_bulk"
  finalize_submission_receipt "$role"
}

submit_production_inclusive_sample() {
  local sample="${1:-}"
  case "$sample" in
    run28_embeddedJet12|run28_embeddedJet20|run28_embeddedJet30|run28_embeddedJet40) ;;
    *) die "Inclusive submission requires run28_embeddedJet12, Jet20, Jet30, or Jet40." ;;
  esac
  local role="production_inclusive_${sample}" config
  config="$(prepare_production_submit "$role" "$inclusive_bulk" "$sample")"
  if find "$inclusive_bulk" -type f -path "*${sample}*" -print -quit 2>/dev/null | grep -q .; then
    die "Output already exists under ${inclusive_bulk} for ${sample}; refusing implicit resume or overwrite."
  fi
  mkdir -p "$evidence_dir"
  exec > >(tee -a "${evidence_dir}/${role}_submission.log") 2>&1
  local -a env_args
  mapfile -t env_args < <(common_env "$config" 0 "$role")
  reserve_submission_receipt "$role" "$config" "$inclusive_bulk" "inclusive:${sample}:${inclusive_bulk}"
  say "Submitting full embedded inclusive lane for ${sample} (no automatic merge)"
  env "${env_args[@]}" \
    "RJ_REQUEST_MEMORY=${sim_memory}" \
    "RJ_SIMEMBEDINCLUSIVE_DEST_BASE=${inclusive_bulk}" \
    "RJ_MERGE_OUT_BASE_OVERRIDE=${merge_base}/inclusive" \
    "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1" \
    "$submitter" isSimEmbeddedInclusive condorDoAll groupSize "$sim_group" "SAMPLE=${sample}"
  record_submission_lane "$role" "inclusive_${sample}" "$inclusive_bulk"
  finalize_submission_receipt "$role"
}

status_report() {
  print_contract
  printf '\nQUEUE (READ ONLY)\n'
  condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus HoldReason Args 2>/dev/null | \
    grep -F "$campaign_tag" || true
  printf '\nNON-TINY ROOT COUNTS\n'
  local path count
  for path in "$smoke_canary_root" "$scan_canary_root" "$confirmation_canary_root" "$data_bulk" "$signal_bulk" "$inclusive_bulk"; do
    count="$(find "$path" -type f -name '*.root' -size +50k 2>/dev/null | wc -l | tr -d ' ')"
    printf '%s\t%s\n' "$count" "$path"
  done
  printf '\nSUBMISSION RECEIPTS\n'
  find "$evidence_dir" -maxdepth 1 -type f \( -name 'submitted_*.receipt' -o -name 'submitted_*.pending.receipt' \) -print 2>/dev/null | sort || true
}

case "$mode" in
  print) print_contract ;;
  preflight) preflight ;;
  build-isolated) build_isolated ;;
  historical-scan) historical_scan ;;
  rank-scan-canary) shift; rank_scan_canary "$@" ;;
  bind-comparison-band) shift; bind_comparison_band "$@" ;;
  validate-canary) shift; validate_canary "$@" ;;
  smoke-canary-submit) submit_transport_smoke "transport_smoke" "$scan_yaml" "$smoke_canary_root" ;;
  scan-canary-submit) submit_discovery_canary "scan_discovery" "$scan_yaml" "$scan_canary_root" ;;
  confirmation-canary-submit)
    confirmation_config="$(generate_confirmation_config)"
    slug="$(sha256sum "$confirmation_config" | awk '{print substr($1,1,12)}')"
    submit_discovery_canary "confirmation_canary_${slug}" "$confirmation_config" "${confirmation_canary_root}/${slug}"
    ;;
  production-dryrun) production_counts ;;
  production-submit-data) submit_production_data ;;
  production-submit-signal) submit_production_signal ;;
  production-submit-inclusive-sample) submit_production_inclusive_sample "${2:-}" ;;
  status) status_report ;;
  help|-h|--help) usage ;;
  *) usage >&2; die "Unknown mode: ${mode}" ;;
esac
