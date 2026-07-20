#!/usr/bin/env bash
# Plan or execute the exact 12-lane source-locked paired executable canary.
# Despite the historical filename, this script does not submit RecoilJets-only
# Condor jobs. Every admitted lane runs the preserved PPG12 executable and
# RecoilJets together through run_ppg12_recoiljets_paired_oracle.sh.

set -Eeuo pipefail
IFS=$'\n\t'
umask 077

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../.." && pwd -P)"
paired_driver="${RJ_PPG12_PAIRED_ORACLE_DRIVER:-${script_dir}/run_ppg12_recoiljets_paired_oracle.sh}"
campaign_driver="${script_dir}/$(basename "${BASH_SOURCE[0]}")"
fanout_helper="${script_dir}/ppg12_paired_oracle_condor_fanout.py"

say() { printf '[PPG12-PAIRED-12] %s\n' "$*"; }
die() { printf '[PPG12-PAIRED-12][ERROR] %s\n' "$*" >&2; exit 2; }

usage() {
  cat <<'EOF'
Usage:
  submit_ppg12_stitched_purity_photon_canaries.sh [plan]
  submit_ppg12_stitched_purity_photon_canaries.sh --submit --token TOKEN
  submit_ppg12_stitched_purity_photon_canaries.sh --condor-plan
  submit_ppg12_stitched_purity_photon_canaries.sh --condor-submit --token TOKEN
  submit_ppg12_stitched_purity_photon_canaries.sh --audit-condor [--plan PATH]

Required environment:
  RJ_PPG12_PAIRED_SOURCE_MANIFEST  absolute JSON source/asset manifest

The default is a non-mutating plan. --submit is a compatibility spelling for
foreground execution of 12 paired executable lanes; it does not submit Condor.

--condor-plan deterministically generates one job per canonical lane and does
not submit. --condor-submit requires the exact generated token and invokes only
condor_submit; it has no merge, promotion, retry, release, or removal behavior.
--audit-condor is read-only and reports lane receipt coverage as JSON.
EOF
}

require_condor_environment() {
  [[ -s "$fanout_helper" ]] || die "Condor fanout helper is missing: $fanout_helper"
  [[ "$paired_driver" == /* && -x "$paired_driver" ]] || \
    die "paired driver is missing: $paired_driver"
  [[ -x "$campaign_driver" ]] || die "campaign driver is not executable: $campaign_driver"
  [[ "${RJ_PPG12_PAIRED_SOURCE_MANIFEST:-}" == /* && \
     -s "${RJ_PPG12_PAIRED_SOURCE_MANIFEST:-}" ]] || \
    die "RJ_PPG12_PAIRED_SOURCE_MANIFEST must name a non-empty absolute JSON file"
  [[ "${RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG:-}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,159}$ ]] || \
    die "Condor mode requires a safe explicit RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG"
  [[ "${RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT:-}" == /* ]] || \
    die "Condor mode requires absolute RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT"
  [[ "${RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR:-}" == /* ]] || \
    die "Condor mode requires absolute RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR"
  [[ -n "${RJ_CODEX_CHAT_NAME:-}" && -n "${RJ_CODEX_THREAD_ID:-}" ]] || \
    die "Condor mode requires explicit RJ_CODEX_CHAT_NAME and RJ_CODEX_THREAD_ID"
}

# The historical foreground interface below remains unchanged.  Condor modes
# are explicit, isolated entry points so planning cannot accidentally submit.
case "${1:-}" in
  --condor-plan)
    shift
    (( $# == 0 )) || die "--condor-plan accepts no additional arguments"
    require_condor_environment
    exec python3 "$fanout_helper" plan \
      --source-manifest "$RJ_PPG12_PAIRED_SOURCE_MANIFEST" \
      --paired-driver "$paired_driver" \
      --campaign-driver "$campaign_driver" \
      --campaign-tag "$RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG" \
      --output-root "$RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT" \
      --evidence-dir "$RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR" \
      --codex-chat-name "$RJ_CODEX_CHAT_NAME" \
      --codex-thread-id "$RJ_CODEX_THREAD_ID"
    ;;
  --condor-submit)
    shift
    provided_condor_token=""
    while (( $# )); do
      case "$1" in
        --token)
          (( $# >= 2 )) || die "--token requires a value"
          provided_condor_token="$2"
          shift 2
          ;;
        *) die "unknown --condor-submit argument: $1" ;;
      esac
    done
    [[ -n "$provided_condor_token" ]] || die "--condor-submit requires --token TOKEN"
    require_condor_environment
    condor_submit_command="${RJ_PPG12_CONDOR_SUBMIT_COMMAND:-$(command -v condor_submit || true)}"
    [[ "$condor_submit_command" == /* && -x "$condor_submit_command" ]] || \
      die "condor_submit is unavailable; set RJ_PPG12_CONDOR_SUBMIT_COMMAND"
    exec python3 "$fanout_helper" submit \
      --plan "$RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR/photon_condor_plan.json" \
      --token "$provided_condor_token" \
      --condor-submit "$condor_submit_command"
    ;;
  --audit-condor)
    shift
    audit_plan=""
    if [[ "${1:-}" == --plan ]]; then
      (( $# >= 2 )) || die "--plan requires a value"
      audit_plan="$2"
      shift 2
    else
      [[ "${RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR:-}" == /* ]] || \
        die "--audit-condor requires --plan PATH or absolute RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR"
      audit_plan="$RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR/photon_condor_plan.json"
    fi
    (( $# == 0 )) || die "--audit-condor accepts only optional --plan PATH"
    [[ "$audit_plan" == /* ]] || die "--audit-condor requires one absolute plan path"
    [[ -s "$fanout_helper" ]] || die "Condor fanout helper is missing: $fanout_helper"
    exec python3 "$fanout_helper" audit --plan "$audit_plan"
    ;;
  --condor-lane)
    shift
    [[ -s "$fanout_helper" ]] || die "Condor fanout helper is missing: $fanout_helper"
    exec python3 "$fanout_helper" lane "$@"
    ;;
esac

mode=plan
provided_token=""
if [[ "${1:-}" == plan ]]; then shift; fi
while (( $# )); do
  case "$1" in
    --submit|--run) mode=run; shift ;;
    --token) (( $# >= 2 )) || die "--token requires a value"; provided_token="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

campaign_tag="${RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG:-ppg12_paired_photon_canary_$(date -u +%Y%m%dT%H%M%SZ)_$$}"
[[ "$campaign_tag" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,159}$ ]] || die "unsafe campaign tag"
source_manifest="${RJ_PPG12_PAIRED_SOURCE_MANIFEST:-}"
[[ "$source_manifest" == /* && -s "$source_manifest" ]] || \
  die "RJ_PPG12_PAIRED_SOURCE_MANIFEST must name a non-empty absolute JSON file"
[[ -s "$fanout_helper" ]] || die "canonical paired-source validator is missing: $fanout_helper"
[[ "$paired_driver" == /* && -x "$paired_driver" ]] || die "paired driver is missing: $paired_driver"
output_root="${RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT:-/tmp/${campaign_tag}}"
evidence_dir="${RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR:-${repo_root}/evidence/qa/${campaign_tag}}"
[[ "$output_root" == /* && "$evidence_dir" == /* ]] || die "output/evidence paths must be absolute"
plan_json="${evidence_dir}/photon_canary_plan.json"
lane_tsv="${evidence_dir}/photon_canary_lanes.tsv"
receipt_tsv="${evidence_dir}/paired_execution_receipts.tsv"
mkdir -p "$evidence_dir"

python3 - "$source_manifest" "$paired_driver" "$output_root" "$campaign_tag" \
  "$plan_json" "$lane_tsv" "$fanout_helper" <<'PY'
from __future__ import annotations
import hashlib, importlib.util, json, sys
from datetime import datetime, timezone
from pathlib import Path

source_path = Path(sys.argv[1])
driver_path = Path(sys.argv[2])
output_root = Path(sys.argv[3])
campaign = sys.argv[4]
plan_path = Path(sys.argv[5])
lane_path = Path(sys.argv[6])
fanout_path = Path(sys.argv[7])

def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def asset(value: object, label: str) -> dict[str, str]:
    path = Path(str(value))
    if not path.is_absolute() or not path.is_file() or path.stat().st_size <= 0:
        raise SystemExit(f"{label} must be an existing non-empty absolute file: {path}")
    return {"path": str(path), "sha256": digest(path)}

spec = importlib.util.spec_from_file_location("ppg12_paired_oracle_fanout", fanout_path)
if spec is None or spec.loader is None:
    raise SystemExit(f"cannot load canonical source validator: {fanout_path}")
fanout = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fanout)
try:
    canonical = fanout.source_contract(source_path, str(output_root))
except fanout.FanoutError as exc:
    raise SystemExit(f"paired source manifest rejected: {exc}") from exc

common = canonical["common"]
driver = asset(driver_path, "paired driver")
lanes = []
for row in canonical["lanes"]:
    lane_id = row["lane_id"]
    sample = row["sample"]
    period = row["period"]
    interaction = row["interaction"]
    output_base = row["output_base"]
    lanes.append({
        "lane_id": lane_id,
        "sample": sample,
        "period": period,
        "interaction": interaction,
        "rows": 5,
        "execution": "source_locked_paired_executable",
        "output_base": output_base,
        "sources": row["sources"],
        "first_five_event_identity_sha256": row[
            "first_five_event_identity_sha256"
        ],
        "expected_evidence": {
            "runtime_contract": output_base + "/paired_oracle_contract.json",
            "candidate_csv": output_base + "/comparison/paired_oracle_candidates.csv",
            "executable_aggregate": output_base + "/comparison/executable_aggregate.json",
            "run_state": output_base + "/RUN_STATE",
        },
    })
auth_payload = {
    "schema": "ppg12-stitched-purity-photon-canary-plan/v5",
    "campaign_tag": campaign,
    "output_root": str(output_root),
    "source_manifest": canonical["source_manifest"],
    "source_validator": asset(fanout_path, "canonical source validator"),
    "paired_driver": driver,
    "common": common,
    "lanes": lanes,
}
encoded = json.dumps(auth_payload, sort_keys=True, separators=(",", ":")).encode()
token = "RUN_PPG12_PAIRED_" + hashlib.sha256(encoded).hexdigest()
plan = dict(auth_payload)
plan.update({
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "status": "PLANNED",
    "submission_token": token,
    "bounded_contract": {
        "lane_count": 12,
        "rows_per_lane": 5,
        "automatic_merge": False,
        "execution": "foreground_source_locked_paired_executable",
        "python_shadow_admissible": False,
        "archived_ppg12_root_admissible": False,
        "external_normalization": "forbidden",
    },
})
Path(plan_path).write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
with Path(lane_path).open("w") as stream:
    stream.write("lane_id\tsample\tperiod\tinteraction\toutput_base\tppg_macro\tg4_full_list\ttruthjet_full_list\n")
    for row in lanes:
        stream.write("\t".join((
            row["lane_id"], row["sample"], row["period"], row["interaction"],
            row["output_base"], row["sources"]["ppg_macro"]["path"],
            row["sources"]["g4_full_list"]["path"],
            row["sources"]["truthjet_full_list"]["path"],
        )) + "\n")
print(token)
PY

expected_token="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["submission_token"])' "$plan_json")"
cat <<EOF
PPG12_STITCHED_PURITY_PAIRED_PHOTON_CANARY_V5
mode=${mode}
campaign_tag=${campaign_tag}
lane_count=12
execution=source_locked_paired_executable
automatic_merge=disabled
output_root=${output_root}
manifest=${plan_json}
lane_table=${lane_tsv}
run_token=${expected_token}
EOF

if [[ "$mode" == plan ]]; then
  [[ -z "$provided_token" ]] || die "--token is valid only with --submit/--run"
  exit 0
fi
[[ "$provided_token" == "$expected_token" ]] || die "run token mismatch"
[[ ! -e "$output_root" ]] || die "output root already exists: $output_root"
mkdir -p "$output_root"
printf 'lane_id\toutput_base\tstatus\truntime_contract\truntime_contract_sha256\tcandidate_csv\tcandidate_csv_sha256\texecutable_aggregate\texecutable_aggregate_sha256\n' > "$receipt_tsv"

common_args="$({ python3 - "$plan_json" <<'PY'
import json, sys
data=json.load(open(sys.argv[1]))
for name, flag in (
 ("setup_script","--setup-script"),("apply_bdt","--apply-bdt"),
 ("apply_config","--apply-config"),("base_e_model","--base-e-model"),
 ("base_v3e_model","--base-v3e-model"),("npb_model","--npb-model"),
 ("tower_mask","--tower-mask"),("recoil_runtime_manifest","--recoil-runtime-manifest"),
 ("recoil_config","--recoil-config"),
): print(flag + "\t" + data["common"][name]["path"])
PY
} )"

while IFS=$'\t' read -r lane_id sample period interaction output_base ppg_macro g4_list truth_list; do
  [[ "$lane_id" != lane_id ]] || continue
  args=(
    --lane-id "$lane_id" --sample "$sample" --period "$period"
    --interaction "$interaction" --output-dir "$output_base"
    --ppg-macro "$ppg_macro" --g4-full-list "$g4_list"
    --truthjet-full-list "$truth_list"
  )
  while IFS=$'\t' read -r flag value; do args+=("$flag" "$value"); done <<< "$common_args"
  lane_plan="$($paired_driver "${args[@]}")"
  lane_token="$(sed -n 's/^  run_token: //p' <<< "$lane_plan")"
  [[ "$lane_token" =~ ^ppg12-oracle:[0-9a-f]{64}$ ]] || die "$lane_id: paired driver emitted no valid token"
  say "run ${lane_id}"
  "$paired_driver" --run --token "$lane_token" "${args[@]}"
  contract="${output_base}/paired_oracle_contract.json"
  candidate="${output_base}/comparison/paired_oracle_candidates.csv"
  aggregate="${output_base}/comparison/executable_aggregate.json"
  [[ "$(tr -d '[:space:]' < "${output_base}/RUN_STATE")" == PASS ]] || die "$lane_id: paired run did not pass"
  for path in "$contract" "$candidate" "$aggregate"; do [[ -s "$path" ]] || die "$lane_id: missing evidence $path"; done
  python3 - "$receipt_tsv" "$lane_id" "$output_base" "$contract" "$candidate" "$aggregate" <<'PY'
import hashlib, sys
from pathlib import Path
receipt, lane, output, contract, candidate, aggregate = sys.argv[1:]
def h(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
with open(receipt, "a") as stream:
    stream.write("\t".join((lane, output, "PASS", contract, h(contract), candidate,
        h(candidate), aggregate, h(aggregate))) + "\n")
PY
done < "$lane_tsv"

say "completed exactly 12 source-locked paired executable lanes; receipts=${receipt_tsv}"
