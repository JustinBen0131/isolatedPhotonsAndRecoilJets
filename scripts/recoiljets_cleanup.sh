#!/usr/bin/env bash
set -euo pipefail

BASE="${RJ_BASE:-/sphenix/u/patsfan753/scratch/thesisAnalysis}"
BULK_BASE="${RJ_BULK_BASE:-/sphenix/tg/tg01/bulk/jbennett}"
THESIS_ANA="${BULK_BASE}/thesisAna"
THESIS_SMOKE="${BULK_BASE}/thesisAnaSmoke"
THESIS_POOLS="${BULK_BASE}/thesisAnaPools"
THESIS_POOLS_SMOKE="${BULK_BASE}/thesisAnaPoolsSmoke"

usage() {
  cat <<'EOF'
Usage:
  ./scripts/recoiljets_cleanup.sh dryrun all
  ./scripts/recoiljets_cleanup.sh apply all
  ./scripts/recoiljets_cleanup.sh dryrun local
  ./scripts/recoiljets_cleanup.sh apply local
  ./scripts/recoiljets_cleanup.sh dryrun smoke
  ./scripts/recoiljets_cleanup.sh apply smoke
  ./scripts/recoiljets_cleanup.sh dryrun dataset <dataset>
  ./scripts/recoiljets_cleanup.sh apply dataset <dataset>

Datasets:
  isPP isPPrun25 isAuAu isOO isSim isSimInclusive isSimJet5 isSimMB
  isSimEmbedded isSimEmbeddedInclusive

Notes:
  dryrun prints exact deletion candidates and does not delete.
  apply refuses bulk cleanup when the user's Condor queue is nonempty unless
  RJ_CLEAN_IGNORE_CONDOR=1 is set.
EOF
}

say() { printf '➜ %s\n' "$*"; }
warn() { printf 'WARN: %s\n' "$*" >&2; }
err() { printf 'ERROR: %s\n' "$*" >&2; }

[[ $# -ge 2 ]] || { usage; exit 2; }
MODE="$1"
TARGET="$2"
DATASET="${3:-}"

case "$MODE" in
  dryrun|apply) ;;
  *) err "First argument must be dryrun or apply, got: $MODE"; usage; exit 2 ;;
esac

case "$TARGET" in
  all|local|smoke|dataset) ;;
  *) err "Second argument must be all, local, smoke, or dataset, got: $TARGET"; usage; exit 2 ;;
esac

if [[ "$TARGET" == "dataset" && -z "$DATASET" ]]; then
  err "dataset target requires a dataset name"
  usage
  exit 2
fi

is_apply() { [[ "$MODE" == "apply" ]]; }

dryrun_prefix() {
  if is_apply; then
    printf '[delete] '
  else
    printf '[dryrun] '
  fi
}

path_allowed_for_bulk_delete() {
  local p="$1"
  case "$p" in
    "${THESIS_ANA}/pp"|"${THESIS_ANA}/pp/"*) return 0 ;;
    "${THESIS_ANA}/pp25"|"${THESIS_ANA}/pp25/"*) return 0 ;;
    "${THESIS_ANA}/auau"|"${THESIS_ANA}/auau/"*) return 0 ;;
    "${THESIS_ANA}/oo"|"${THESIS_ANA}/oo/"*) return 0 ;;
    "${THESIS_ANA}/sim"|"${THESIS_ANA}/sim/"*) return 0 ;;
    "${THESIS_ANA}/siminclusive"|"${THESIS_ANA}/siminclusive/"*) return 0 ;;
    "${THESIS_ANA}/simjet5"|"${THESIS_ANA}/simjet5/"*) return 0 ;;
    "${THESIS_ANA}/simmb"|"${THESIS_ANA}/simmb/"*) return 0 ;;
    "${THESIS_ANA}/simembedded"|"${THESIS_ANA}/simembedded/"*) return 0 ;;
    "${THESIS_ANA}/simembeddedinclusive"|"${THESIS_ANA}/simembeddedinclusive/"*) return 0 ;;
    "${THESIS_SMOKE}"|"${THESIS_SMOKE}/"*) return 0 ;;
    "${THESIS_POOLS_SMOKE}"|"${THESIS_POOLS_SMOKE}/"*) return 0 ;;
    "${THESIS_POOLS}/"*poolSmoke*) return 0 ;;
  esac
  return 1
}

remove_path() {
  local p="$1"
  [[ -n "$p" && "$p" != "/" ]] || { err "Refusing unsafe path: '$p'"; exit 60; }
  if [[ ! -e "$p" ]]; then
    return 0
  fi
  if ! path_allowed_for_bulk_delete "$p"; then
    err "Refusing to delete path outside approved cleanup roots: $p"
    exit 61
  fi
  dryrun_prefix
  printf '%s\n' "$p"
  if is_apply; then
    rm -rf -- "$p"
  fi
}

remove_file_local() {
  local p="$1"
  [[ -e "$p" ]] || return 0
  dryrun_prefix
  printf '%s\n' "$p"
  if is_apply; then
    rm -f -- "$p"
  fi
}

remove_dir_local() {
  local p="$1"
  [[ -d "$p" ]] || return 0
  case "$p" in
    "${BASE}/tmp_recoil_merge_"*|"${BASE}/.recoiljets_tmp"*|"${BASE}/condor_sub/"*smoke*|"${BASE}/condor_sub/pool_workflow_"*) ;;
    *) err "Refusing local directory outside approved patterns: $p"; exit 62 ;;
  esac
  dryrun_prefix
  printf '%s\n' "$p"
  if is_apply; then
    rm -rf -- "$p"
  fi
}

require_no_condor_jobs_for_bulk_apply() {
  if ! is_apply; then
    return 0
  fi
  [[ "${RJ_CLEAN_IGNORE_CONDOR:-0}" == "1" ]] && return 0
  command -v condor_q >/dev/null 2>&1 || {
    warn "condor_q not found; cannot verify active jobs. Set RJ_CLEAN_IGNORE_CONDOR=1 to override."
    exit 63
  }
  if condor_q -af ClusterId 2>/dev/null | grep -q .; then
    err "Refusing bulk cleanup because your Condor queue is nonempty."
    condor_q -nobatch || true
    exit 64
  fi
}

dataset_output_root() {
  case "$1" in
    isPP|pp) echo "${THESIS_ANA}/pp" ;;
    isPPrun25|pp25) echo "${THESIS_ANA}/pp25" ;;
    isAuAu|auau) echo "${THESIS_ANA}/auau" ;;
    isOO|oo) echo "${THESIS_ANA}/oo" ;;
    isSim|sim) echo "${THESIS_ANA}/sim" ;;
    isSimInclusive|isSimJet5|siminclusive|simjet5) echo "${THESIS_ANA}/siminclusive" ;;
    isSimMB|simmb) echo "${THESIS_ANA}/simmb" ;;
    isSimEmbedded|simembedded) echo "${THESIS_ANA}/simembedded" ;;
    isSimEmbeddedInclusive|simembeddedinclusive) echo "${THESIS_ANA}/simembeddedinclusive" ;;
    *) err "Unknown dataset: $1"; exit 2 ;;
  esac
}

clean_local_artifacts() {
  say "Scanning disposable local test/staging artifacts under ${BASE}"
  [[ -d "$BASE" ]] || { err "BASE does not exist: $BASE"; exit 65; }

  # Merge-stage temp directories may belong to active DAGMan/hadd workflows on
  # another submit host. Local smoke cleanup must not remove them.
  while IFS= read -r d; do remove_dir_local "$d"; done < <(
    find "$BASE" -maxdepth 1 -type d -name '.recoiljets_tmp*' 2>/dev/null | sort
  )
  while IFS= read -r d; do remove_dir_local "$d"; done < <(
    find "${BASE}/condor_sub" -maxdepth 1 -type d \( -name '*smoke*' -o -name 'pool_workflow_*' \) 2>/dev/null | sort
  )
  while IFS= read -r f; do remove_file_local "$f"; done < <(
    find "${BASE}/condor_yaml_overrides" -type f \( -name '*LOCAL*.yaml' -o -name '*LOCALSTITCHTEST*.yaml' \) 2>/dev/null | sort
  )
  while IFS= read -r f; do remove_file_local "$f"; done < <(
    find "${BASE}/condor_lists" -type f -name '*LOCAL*' 2>/dev/null | sort
  )
  say "Scanning stdout/error/log files under ${BASE}"
  for _log_dir in stdout error log; do
    [[ -d "${BASE}/${_log_dir}" ]] || continue
    while IFS= read -r f; do remove_file_local "$f"; done < <(
      find "${BASE}/${_log_dir}" -type f 2>/dev/null | sort
    )
  done
}

clean_smoke_bulk() {
  require_no_condor_jobs_for_bulk_apply
  say "Scanning smoke/pool bulk artifacts"
  remove_path "$THESIS_SMOKE"
  remove_path "$THESIS_POOLS_SMOKE"
  while IFS= read -r d; do remove_path "$d"; done < <(
    find "$THESIS_POOLS" -mindepth 1 -maxdepth 1 -type d -name '*poolSmoke*' 2>/dev/null | sort
  )
}

clean_dataset_bulk() {
  local dataset="$1"
  local root
  root="$(dataset_output_root "$dataset")"
  require_no_condor_jobs_for_bulk_apply
  say "Scanning dataset output root: ${root}"
  if [[ ! -d "$root" ]]; then
    say "No such output root yet: ${root}"
    return 0
  fi
  while IFS= read -r d; do remove_path "$d"; done < <(
    find "$root" -mindepth 1 -maxdepth 1 -print 2>/dev/null | sort
  )
}

say "Mode=${MODE} Target=${TARGET}${DATASET:+ Dataset=${DATASET}}"
case "$TARGET" in
  local)
    clean_local_artifacts
    ;;
  smoke)
    clean_smoke_bulk
    ;;
  dataset)
    clean_dataset_bulk "$DATASET"
    ;;
  all)
    clean_local_artifacts
    clean_smoke_bulk
    for ds in isPP isAuAu isSim isSimInclusive isSimEmbedded isSimEmbeddedInclusive; do
      clean_dataset_bulk "$ds"
    done
    ;;
esac

say "Cleanup ${MODE} complete."
