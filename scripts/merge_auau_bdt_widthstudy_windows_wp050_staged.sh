#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

stage="${1:-}"
case "$stage" in
  firstRound|secondRound|finalStitch) ;;
  *)
    cat >&2 <<'EOF'
Usage:
  ./scripts/merge_auau_bdt_widthstudy_windows_wp050_staged.sh firstRound
  ./scripts/merge_auau_bdt_widthstudy_windows_wp050_staged.sh secondRound
  ./scripts/merge_auau_bdt_widthstudy_windows_wp050_staged.sh finalStitch

This helper intentionally runs exactly one merge stage. Wait for the Condor
queue/READY email from one stage before launching the next stage.
EOF
    exit 2
    ;;
esac

campaign="${RJ_WIDTHSTUDY_CAMPAIGN_TAG:-widthstudy_windows_wp050_fixed_20260511_180500}"
windows_csv="${RJ_WIDTHSTUDY_WINDOWS:-5:35,10:35,15:35}"
merge_memory="${RJ_SIM_FIRSTROUND_REQUEST_MEMORY:-8000MB}"
merge_group="${RJ_SIM_MERGE_GROUP_SIZE:-75}"
allow_live_jobs="${RJ_MERGE_ALLOW_LIVE_JOBS:-0}"

pt_tag() {
  local lo="$1"
  local hi="$2"
  printf 'pt%sto%s' "${lo%.*}" "${hi%.*}"
}

normalize_window() {
  local item="${1//[[:space:]]/}"
  if [[ "$item" =~ ^pt([0-9]+)to([0-9]+)$ ]]; then
    printf '%s:%s\n' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
  elif [[ "$item" =~ ^([0-9]+):([0-9]+)$ ]]; then
    printf '%s\n' "$item"
  else
    echo "[ERROR] Bad window '${1}'. Use 5:35 or pt5to35." >&2
    return 2
  fi
}

live_jobs=0
if command -v condor_q >/dev/null 2>&1; then
  live_jobs="$(condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null | awk '/Total for query:/ {print $4; found=1} END {if (!found) print 0}')"
  live_jobs="${live_jobs:-0}"
fi

if [[ "$allow_live_jobs" != "1" && "$live_jobs" != "0" ]]; then
  echo "[ERROR] Refusing to run ${stage} while ${live_jobs} Condor jobs are still visible for ${USER:-patsfan753}." >&2
  echo "        Wait for the prior stage to drain, or set RJ_MERGE_ALLOW_LIVE_JOBS=1 if you know these jobs are unrelated." >&2
  exit 3
fi

log_dir="${repo_root}/condor_generated_configs/${campaign}"
mkdir -p "$log_dir"
log="${log_dir}/merge_${stage}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$log") 2>&1

echo "RECOILJETS_AUAU_BDT_WIDTHSTUDY_WINDOWS_STAGED_MERGE_V1"
echo "host=$(hostname -f 2>/dev/null || hostname)"
echo "repo_root=${repo_root}"
echo "campaign=${campaign}"
echo "stage=${stage}"
echo "windows=${windows_csv}"
echo "merge_memory=${merge_memory}"
echo "merge_group_size=${merge_group}"
echo "live_jobs_before=${live_jobs}"
echo "log=${log}"
echo

run_merge() {
  local dataset="$1"
  local bulk="$2"
  local merge_base="$3"
  local yaml="$4"
  local input_base="$bulk"
  local -a stage_args=( "$stage" )

  if [[ "$stage" != "firstRound" ]]; then
    case "$dataset" in
      isSimEmbedded) input_base="${merge_base}/simembedded" ;;
      isSimEmbeddedInclusive) input_base="${merge_base}/simembeddedinclusive" ;;
    esac
  fi

  if [[ "$stage" == "secondRound" || "$stage" == "finalStitch" ]]; then
    stage_args+=( "condor" )
  fi

  echo
  echo "====================================================================="
  echo "MERGE stage=${stage} dataset=${dataset}"
  echo "  yaml       : ${yaml}"
  echo "  input base : ${input_base}"
  echo "  output base: ${merge_base}"
  echo "====================================================================="

  [[ -f "$yaml" ]] || { echo "[ERROR] Missing YAML: ${yaml}" >&2; exit 4; }
  [[ -d "$input_base" ]] || { echo "[ERROR] Missing merge input base for stage=${stage}: ${input_base}" >&2; exit 5; }

  env \
    MERGE_CONFIG_YAML="${yaml}" \
    MERGE_SIM_INPUT_BASE_OVERRIDE="${input_base}" \
    MERGE_OUT_BASE_OVERRIDE="${merge_base}" \
    RJ_SIM_FIRSTROUND_REQUEST_MEMORY="${merge_memory}" \
    RJ_AUTO_MEMORY_RETRY_CAP_MB="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-16000}" \
    ./scripts/mergeRecoilJets.sh "${dataset}" "${stage_args[@]}" groupSize "${merge_group}"
}

IFS=',' read -r -a window_items <<< "$windows_csv"
for raw in "${window_items[@]}"; do
  window="$(normalize_window "$raw")"
  lo="${window%%:*}"
  hi="${window##*:}"
  tag="$(pt_tag "$lo" "$hi")"
  bulk_tag="${campaign}_${tag}"
  yaml="${repo_root}/condor_generated_configs/${campaign}/analysis_config_auau_bdt_widthstudy_${tag}_wp050.yaml"
  merge_base="/sphenix/u/patsfan753/scratch/thesisAnalysis/output_${bulk_tag}"
  signal_bulk="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_${bulk_tag}"
  background_bulk="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive_${bulk_tag}"

  run_merge isSimEmbedded "$signal_bulk" "$merge_base" "$yaml"
  run_merge isSimEmbeddedInclusive "$background_bulk" "$merge_base" "$yaml"
done

echo
echo "== post-stage queue snapshot =="
condor_q -nobatch "${USER:-patsfan753}" 2>/dev/null || true
echo
echo "STAGED_MERGE_LOG=${log}"
echo "NEXT_STEP="
case "$stage" in
  firstRound) echo "  Wait for firstRound merge jobs to finish/READY, then run this helper with secondRound." ;;
  secondRound) echo "  Wait for secondRound to finish, then run this helper with finalStitch." ;;
  finalStitch) echo "  Pull the final SIM outputs locally after READY/output inspection." ;;
esac
