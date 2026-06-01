#!/usr/bin/env bash
# Shared I/O path defaults for RecoilJets SDCC workflows.

if [[ "${RJ_IO_PATHS_LOADED:-0}" == "1" ]]; then
  return 0 2>/dev/null || exit 0
fi
RJ_IO_PATHS_LOADED=1

rj_find_repo_root() {
  local start="${1:-${PWD}}"
  local dir
  dir="$(cd "$start" && pwd -P)"
  while [[ "$dir" != "/" ]]; do
    if [[ -d "$dir/macros" && -d "$dir/scripts" ]]; then
      printf '%s\n' "$dir"
      return 0
    fi
    dir="$(dirname "$dir")"
  done
  return 1
}

rj_validate_path_token() {
  local value="${1:-}"
  local label="${2:-path token}"
  case "$value" in
    ""|*[[:space:]]*|*agent_context*|*.codex*|*codex*|*THE-*|*"/../"*|*"//"*|*/..|../*)
      printf '[ERROR] Unsafe %s: %s\n' "$label" "$value" >&2
      return 2
      ;;
  esac
}

rj_validate_campaign_tag() {
  rj_validate_path_token "$1" "campaign tag"
}

rj_checkout_base() {
  printf '%s\n' "${RJ_REMOTE_BASE:-/sphenix/u/${USER:-patsfan753}/scratch/thesisAnalysis}"
}

rj_bulk_base() {
  printf '%s\n' "${RJ_BULK_BASE:-/sphenix/tg/tg01/bulk/jbennett}"
}

rj_thesis_ana_root() {
  printf '%s\n' "${RJ_THESIS_ANA_ROOT:-$(rj_bulk_base)/thesisAna}"
}

rj_recoiljets_bulk_root() {
  printf '%s\n' "${RJ_RECOILJETS_BULK_ROOT:-$(rj_thesis_ana_root)/recoiljets}"
}

rj_merge_root_for_campaign() {
  local campaign="${1:?campaign required}"
  rj_validate_campaign_tag "$campaign" || return 2
  printf '%s/%s\n' "${RJ_RECOILJETS_RUN_ROOT:-$(rj_checkout_base)/runs/recoiljets/current}" "$campaign"
}

rj_bulk_signal_root_for_campaign() {
  local family="${1:?family required}"
  local campaign="${2:?campaign required}"
  rj_validate_campaign_tag "$family" || return 2
  rj_validate_campaign_tag "$campaign" || return 2
  printf '%s/%s/%s/signal\n' "$(rj_recoiljets_bulk_root)" "$family" "$campaign"
}

rj_bulk_background_root_for_campaign() {
  local family="${1:?family required}"
  local campaign="${2:?campaign required}"
  rj_validate_campaign_tag "$family" || return 2
  rj_validate_campaign_tag "$campaign" || return 2
  printf '%s/%s/%s/background\n' "$(rj_recoiljets_bulk_root)" "$family" "$campaign"
}
