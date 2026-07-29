#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../.." && pwd -P)"
cd "$repo_root"

say() { printf '[THE126] %s\n' "$*"; }
die() { printf '[THE126][ERROR] %s\n' "$1" >&2; exit "${2:-2}"; }

mode="${1:-print}"
campaign_tag="${RJ_THE126_CAMPAIGN_TAG:-the126_auau_lumi10_apurva_repro_20260729_v2}"
work_root="${RJ_THE126_WORK_ROOT:-${repo_root}/.recoiljets_tmp/${campaign_tag}}"
expected_dir="${RJ_THE126_EXPECTED_DIR:-${repo_root}/condor_snapshots/audit_auau_20260709_125304/expected_lists}"
macro="${RJ_THE126_MACRO:-${repo_root}/macros/diagnostics/luminosity/Fun4All_AuAuLumi10.C}"
preserved_aggregate="${RJ_THE126_PRESERVED_AGGREGATE:-${repo_root}/dataOutput/auauLumiAudit/THE126_auau_lumi_10run_reproduction_20260728/aggregate_by_run.tsv}"
revised_csv="${RJ_THE126_REVISED_CSV:-${repo_root}/dataOutput/auauLumiAudit/THE126B_auau_lumi_or_10run_reproduction_20260729/apurva_revised_offline.csv}"
chunk_size="${RJ_THE126_CHUNK_SIZE:-10}"
request_memory="${RJ_THE126_MEMORY:-4000MB}"
dbtag="${RJ_THE126_CDB_GLOBALTAG:-newcdbtag}"
thread_id="${RJ_CODEX_THREAD_ID:-manual}"
chat_name="${RJ_CODEX_CHAT_NAME:-THE-126 | AuAu 10-Run Lumi OR Repro}"
runs=(67597 67598 67599 67600 67601 67602 67633 67634 67645 67660)
canary_runs=(67597 67633)

chunks_dir="${work_root}/chunks"
results_dir="${work_root}/results"
logs_dir="${work_root}/logs"
meta_dir="${work_root}/meta"
wrapper="${work_root}/run_chunk.sh"
all_items="${meta_dir}/all_items.tsv"
canary_items="${meta_dir}/canary_items.tsv"
remaining_items="${meta_dir}/remaining_items.tsv"
canary_sub="${meta_dir}/canary.sub"
remaining_sub="${meta_dir}/remaining.sub"
aggregate_tsv="${work_root}/aggregate_by_run.tsv"
target_tsv="${work_root}/apurva_targets.tsv"

usage() {
  cat <<'EOF'
Usage: submit_the126_auau_lumi10.sh MODE

Modes:
  print             Print the immutable campaign contract.
  prepare           Create deterministic ten-file chunks and Condor files.
  refresh-runtime   Refresh wrapper/submit files without changing chunk lists.
  smoke             Run 100 events from the first run/chunk interactively.
  canary-submit     Submit all chunks for runs 67597 and 67633.
  canary-validate   Aggregate and compare the two complete canary runs.
  remaining-submit Submit the other eight runs after canary validation.
  status            Report queue and output coverage.
  aggregate         Aggregate available chunk TSVs by run.
  regress-preserved Recompute the revised luminosities from preserved counters.

No cleanup, removal, automatic resubmission, or automatic merge is implemented.
EOF
}

require_contract() {
  [[ "$campaign_tag" == "the126_auau_lumi10_apurva_repro_20260729_v2" ]] ||
    die "Refusing non-canonical campaign tag: ${campaign_tag}"
  case "$work_root" in
    "${repo_root}/.recoiljets_tmp/${campaign_tag}") ;;
    *) die "Refusing work root outside isolated campaign path: ${work_root}" ;;
  esac
  [[ -s "$macro" ]] || die "Missing macro: ${macro}"
  [[ -d "$expected_dir" ]] || die "Missing THE-91 expected-list directory: ${expected_dir}"
  [[ "$chunk_size" =~ ^[1-9][0-9]*$ ]] || die "Invalid chunk size: ${chunk_size}"
}

print_contract() {
  cat <<EOF
THE126_AUAU_LUMI10_APURVA_REPRO_V2
campaign_tag=${campaign_tag}
work_root=${work_root}
expected_dir=${expected_dir}
macro=${macro}
macro_sha256=$(sha256sum "$macro" | awk '{print $1}')
runs=${runs[*]}
canary_runs=${canary_runs[*]}
dataset=run3auau/pro001_pcdb001_v001/DST_CALOFITTING
chunk_size=${chunk_size}
dbtag=${dbtag}
sigma_mbd_b=6.324
official_lumi_numerator=z10_trig12_or_trig14
selection_order=FlagHandler,diagnostic_input_counter,CaloStatusSkimmer,MbdReco,GlobalVertexReco,TriggerRunInfoReco,selection_counter
reference_macro_blob=d132e704544300897568017a57e0756b8a411ad8
reference_triggerqa_blob=4f5a77f28f61c6c7e42faf6c280a06e686d57789
RJ_CODEX_CHAT_NAME=${chat_name}
RJ_CODEX_THREAD_ID=${thread_id}
EOF
}

write_wrapper() {
  cat > "$wrapper" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
list="$1"
output_root="$2"
output_tsv="$3"
run="$4"
chunk="$5"
macro="$6"
dbtag="$7"
chunk_decimal="$((10#${chunk}))"

set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.561
set -u

mkdir -p "$(dirname "$output_root")" "$(dirname "$output_tsv")"
root -l -b -q "${macro}(\"${list}\",\"${output_root}\",\"${output_tsv}\",${run},${chunk_decimal},0,\"${dbtag}\")"
[[ -s "$output_root" ]]
[[ -s "$output_tsv" ]]
EOF
  chmod 0755 "$wrapper"
}

write_submit_file() {
  local items="$1" submit_file="$2"
  cat > "$submit_file" <<EOF
universe = vanilla
executable = ${wrapper}
arguments = \$(list) \$(output_root) \$(output_tsv) \$(run) \$(chunk) ${macro} ${dbtag}
output = ${logs_dir}/\$(run)_\$(chunk).out
error = ${logs_dir}/\$(run)_\$(chunk).err
log = ${logs_dir}/\$(run)_\$(chunk).log
request_cpus = 1
request_memory = ${request_memory}
should_transfer_files = NO
getenv = True
+JobBatchName = "${campaign_tag}"
+RJCodexChatName = "${chat_name}"
+RJCodexThreadId = "${thread_id}"
queue run,chunk,list,output_root,output_tsv from ${items}
EOF
}

prepare() {
  require_contract
  [[ ! -e "$work_root" ]] || die "Campaign work root already exists: ${work_root}"
  mkdir -p "$chunks_dir" "$results_dir" "$logs_dir" "$meta_dir"

  cat > "$target_tsv" <<'EOF'
run	apurva_lumi_ub_inv
67597	0.3216
67598	0.5784
67599	1.8757
67600	1.5360
67601	1.3258
67602	0.9596
67633	0.0955
67634	0.4128
67645	1.0479
67660	1.2438
EOF

  : > "$all_items"
  : > "$canary_items"
  : > "$remaining_items"
  local run src run_dir chunk_file chunk_id output_root output_tsv is_canary
  for run in "${runs[@]}"; do
    src="${expected_dir}/dst_calofitting-$(printf '%08d' "$run").list"
    [[ -s "$src" ]] || die "Missing expected list for run ${run}: ${src}"
    run_dir="${chunks_dir}/${run}"
    mkdir -p "$run_dir"
    split -d -a 4 -l "$chunk_size" --additional-suffix=.list \
      "$src" "${run_dir}/chunk_"
    is_canary=0
    [[ "$run" == "67597" || "$run" == "67633" ]] && is_canary=1
    for chunk_file in "${run_dir}"/chunk_*.list; do
      chunk_id="$(basename "$chunk_file" .list | sed 's/^chunk_//')"
      output_root="${results_dir}/${run}_${chunk_id}.root"
      output_tsv="${results_dir}/${run}_${chunk_id}.tsv"
      printf '%s\t%s\t%s\t%s\t%s\n' \
        "$run" "$chunk_id" "$chunk_file" "$output_root" "$output_tsv" |
        tee -a "$all_items" >/dev/null
      if [[ "$is_canary" == "1" ]]; then
        tail -n 1 "$all_items" >> "$canary_items"
      else
        tail -n 1 "$all_items" >> "$remaining_items"
      fi
    done
  done

  write_wrapper
  write_submit_file "$canary_items" "$canary_sub"
  write_submit_file "$remaining_items" "$remaining_sub"
  print_contract > "${meta_dir}/contract.txt"
  sha256sum "$macro" "$all_items" "$canary_items" "$remaining_items" \
    "$wrapper" "$canary_sub" "$remaining_sub" > "${meta_dir}/sha256.txt"

  local files chunks
  files="$(awk '{n+=1} END{print n+0}' "$all_items")"
  chunks="$(wc -l < "$all_items")"
  say "Prepared ${chunks} chunks covering ${files} chunk rows."
  say "Input files: $(awk '{while((getline line < $3)>0)n++; close($3)} END{print n+0}' "$all_items")"
  say "Canary jobs: $(wc -l < "$canary_items"); remaining jobs: $(wc -l < "$remaining_items")"
}

refresh_runtime() {
  require_contract
  [[ -s "$all_items" && -s "$canary_items" && -s "$remaining_items" ]] ||
    die "Prepared item files are missing; run prepare first."
  write_wrapper
  write_submit_file "$canary_items" "$canary_sub"
  write_submit_file "$remaining_items" "$remaining_sub"
  print_contract > "${meta_dir}/contract.txt"
  sha256sum "$macro" "$all_items" "$canary_items" "$remaining_items" \
    "$wrapper" "$canary_sub" "$remaining_sub" > "${meta_dir}/sha256.txt"
  say "Refreshed runtime wrapper and submit files without changing chunks."
}

smoke() {
  require_contract
  [[ -s "$all_items" && -x "$wrapper" ]] || die "Run prepare first."
  local run chunk list output_root output_tsv
  IFS=$'\t' read -r run chunk list output_root output_tsv < "$all_items"
  output_root="${results_dir}/smoke_${run}_${chunk}.root"
  output_tsv="${results_dir}/smoke_${run}_${chunk}.tsv"
  local smoke_list="${meta_dir}/smoke.list"
  head -n 1 "$list" > "$smoke_list"
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.561
  set -u
  export RJ_CODEX_CHAT_NAME="$chat_name"
  export RJ_CODEX_THREAD_ID="$thread_id"
  root -l -b -q "${macro}(\"${smoke_list}\",\"${output_root}\",\"${output_tsv}\",${run},${chunk},100,\"${dbtag}\")"
  [[ -s "$output_root" && -s "$output_tsv" ]] || die "Smoke output missing."
  cat "$output_tsv"
}

assert_no_active_duplicate() {
  if condor_q "${USER:-patsfan753}" -af JobBatchName 2>/dev/null |
      grep -Fx "$campaign_tag" >/dev/null; then
    die "Active jobs already use campaign tag ${campaign_tag}" 4
  fi
}

submit_items() {
  local submit_file="$1" label="$2"
  require_contract
  [[ -s "$submit_file" ]] || die "Missing submit file: ${submit_file}"
  assert_no_active_duplicate
  say "Submitting ${label}: $(wc -l < "$(awk '/^queue /{print $NF}' "$submit_file")") jobs"
  export RJ_CODEX_CHAT_NAME="$chat_name"
  export RJ_CODEX_THREAD_ID="$thread_id"
  condor_submit "$submit_file"
}

aggregate() {
  require_contract
  [[ -d "$results_dir" ]] || die "Missing results directory."
  python3 - "$results_dir" "$target_tsv" "$aggregate_tsv" <<'PY'
import csv
import pathlib
import sys

results_dir = pathlib.Path(sys.argv[1])
target_path = pathlib.Path(sys.argv[2])
output_path = pathlib.Path(sys.argv[3])
targets = {}
with target_path.open() as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        targets[int(row["run"])] = float(row["apurva_lumi_ub_inv"])

count_fields = [
    "input_files", "input_events", "pass_calo_status", "has_z", "z10",
    "trig12", "trig14", "z10_trig12_or_trig14", "z10_trig12",
    "z10_trig14",
]
totals = {}
chunks = {}
for path in sorted(results_dir.glob("[0-9]*_[0-9]*.tsv")):
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise SystemExit(f"invalid chunk TSV {path}: rows={len(rows)}")
    row = rows[0]
    run = int(row["run"])
    totals.setdefault(run, {field: 0 for field in count_fields})
    chunks.setdefault(run, 0)
    chunks[run] += 1
    for field in count_fields:
        totals[run][field] += int(row[field])

fields = [
    "run", "chunks", *count_fields, "sigma_mbd_b", "lumi_ub_inv",
    "lumi_trig12_ub_inv", "lumi_trig14_ub_inv",
    "apurva_lumi_ub_inv", "delta_ub_inv",
]
with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for run in sorted(totals):
        lumi = totals[run]["z10_trig12_or_trig14"] / 6.324e6
        lumi_trig12 = totals[run]["z10_trig12"] / 6.324e6
        lumi_trig14 = totals[run]["z10_trig14"] / 6.324e6
        target = targets[run]
        writer.writerow({
            "run": run,
            "chunks": chunks[run],
            **totals[run],
            "sigma_mbd_b": "6.324000",
            "lumi_ub_inv": f"{lumi:.10f}",
            "lumi_trig12_ub_inv": f"{lumi_trig12:.10f}",
            "lumi_trig14_ub_inv": f"{lumi_trig14:.10f}",
            "apurva_lumi_ub_inv": f"{target:.4f}",
            "delta_ub_inv": f"{lumi-target:+.10f}",
        })
print(output_path)
PY
  cat "$aggregate_tsv"
}

regress_preserved() {
  [[ -s "$preserved_aggregate" ]] ||
    die "Missing preserved aggregate: ${preserved_aggregate}"
  [[ -s "$revised_csv" ]] ||
    die "Missing revised Apurva CSV: ${revised_csv}"
  python3 - "$preserved_aggregate" "$revised_csv" <<'PY'
import csv
import math
import pathlib
import sys

aggregate_path = pathlib.Path(sys.argv[1])
revised_path = pathlib.Path(sys.argv[2])
run_order = [67597, 67598, 67599, 67600, 67601,
             67602, 67633, 67634, 67645, 67660]
expected_numerator = 59_428_026
expected_lumi_ub_inv = 9.397221062619

with aggregate_path.open() as handle:
    rows = {
        int(row["run"]): row
        for row in csv.DictReader(handle, delimiter="\t")
    }
with revised_path.open() as handle:
    revised = {
        int(row["Run"]): float(row["Lumi_ub_inv"])
        for row in csv.DictReader(handle)
    }

missing = [run for run in run_order if run not in rows or run not in revised]
if missing:
    raise SystemExit(f"missing regression runs: {missing}")

matches = 0
numerator = 0
for run in run_order:
    count = int(rows[run]["z10_trig12_or_trig14"])
    lumi = count / 6.324e6
    target = revised[run]
    match = f"{lumi:.4f}" == f"{target:.4f}"
    matches += int(match)
    numerator += count
    print(
        f"{run}\t{count}\t{lumi:.10f}\t{target:.4f}\t"
        f"{'MATCH' if match else 'MISMATCH'}"
    )

lumi_sum = numerator / 6.324e6
if matches != len(run_order):
    raise SystemExit(f"per-run regression failed: {matches}/{len(run_order)}")
if numerator != expected_numerator:
    raise SystemExit(
        f"numerator regression failed: {numerator} != {expected_numerator}"
    )
if not math.isclose(
    lumi_sum, expected_lumi_ub_inv, rel_tol=0.0, abs_tol=5e-13
):
    raise SystemExit(
        f"luminosity regression failed: {lumi_sum:.12f} "
        f"!= {expected_lumi_ub_inv:.12f}"
    )
print(
    f"PASS runs={len(run_order)} matches={matches} "
    f"numerator={numerator} lumi_ub_inv={lumi_sum:.12f}"
)
PY
}

validate_canaries() {
  aggregate
  local expected actual
  for run in "${canary_runs[@]}"; do
    expected="$(find "${chunks_dir}/${run}" -type f -name 'chunk_*.list' | wc -l)"
    actual="$(find "$results_dir" -maxdepth 1 -type f -name "${run}_*.tsv" | wc -l)"
    [[ "$actual" -eq "$expected" ]] ||
      die "Canary ${run} incomplete: ${actual}/${expected} chunk TSVs"
  done
  say "Both complete-run canaries have full chunk coverage."
}

status() {
  require_contract
  print_contract
  if [[ -d "$work_root" ]]; then
    printf 'prepared_chunks=%s\n' "$(wc -l < "$all_items" 2>/dev/null || echo 0)"
    printf 'result_tsv=%s\n' "$(find "$results_dir" -maxdepth 1 -type f -name '[0-9]*_[0-9]*.tsv' 2>/dev/null | wc -l)"
    printf 'result_root=%s\n' "$(find "$results_dir" -maxdepth 1 -type f -name '[0-9]*_[0-9]*.root' -size +1k 2>/dev/null | wc -l)"
  fi
  condor_q "${USER:-patsfan753}" -af ClusterId ProcId JobStatus JobBatchName Args 2>/dev/null |
    grep -F "$campaign_tag" || true
}

case "$mode" in
  print) require_contract; print_contract ;;
  prepare) prepare ;;
  refresh-runtime) refresh_runtime ;;
  smoke) smoke ;;
  canary-submit) submit_items "$canary_sub" canaries ;;
  canary-validate) validate_canaries ;;
  remaining-submit) submit_items "$remaining_sub" remaining-runs ;;
  status) status ;;
  aggregate) aggregate ;;
  regress-preserved) regress_preserved ;;
  -h|--help|help) usage ;;
  *) usage; die "Unknown mode: ${mode}" ;;
esac
