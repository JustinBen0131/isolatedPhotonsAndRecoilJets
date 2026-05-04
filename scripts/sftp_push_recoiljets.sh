#!/usr/bin/env bash
set -euo pipefail

LOCAL_BASE="/Users/patsfan753/Desktop/ThesisAnalysis"
REMOTE_BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
REMOTE_HOST="patsfan753@sftp.sdcc.bnl.gov"

LOCAL_FILES=(
  "scripts/audit_auau_grl_projection.sh"
  "scripts/estimateEmbeddedPhotonXsec.sh"
  "scripts/make_dstListsData.sh"
  "scripts/makeThesisSimLists.sh"
  "scripts/mergeRecoilJets.sh"
  "scripts/root_in_analysis_env.sh"
  "scripts/RecoilJets_Condor_AuAu.sh"
  "scripts/RecoilJets_Condor_submit.sh"
  "scripts/RecoilJets_Condor.sh"
  "macros/analysis_config.yaml"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/PrintPPStitchDiagnostics.C"
  "src/RecoilJets.cc"
  "src/RecoilJets.h"
  "src_AuAu/RecoilJets_AuAu.cc"
  "src_AuAu/RecoilJets_AuAu.h"
)

REMOTE_FILES=(
  "scripts/audit_auau_grl_projection.sh"
  "scripts/estimateEmbeddedPhotonXsec.sh"
  "scripts/make_dstListsData.sh"
  "scripts/makeThesisSimLists.sh"
  "scripts/mergeRecoilJets.sh"
  "scripts/root_in_analysis_env.sh"
  "RecoilJets_Condor_AuAu.sh"
  "RecoilJets_Condor_submit.sh"
  "RecoilJets_Condor.sh"
  "macros/analysis_config.yaml"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/PrintPPStitchDiagnostics.C"
  "src/RecoilJets.cc"
  "src/RecoilJets.h"
  "src_AuAu/RecoilJets_AuAu.cc"
  "src_AuAu/RecoilJets_AuAu.h"
)

GROUP_CONDOR=(
  "scripts/RecoilJets_Condor.sh"
  "scripts/RecoilJets_Condor_AuAu.sh"
  "scripts/RecoilJets_Condor_submit.sh"
)

GROUP_MACROS=(
  "macros/analysis_config.yaml"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/PrintPPStitchDiagnostics.C"
)

GROUP_SCRIPTS=(
  "scripts/audit_auau_grl_projection.sh"
  "scripts/estimateEmbeddedPhotonXsec.sh"
  "scripts/make_dstListsData.sh"
  "scripts/makeThesisSimLists.sh"
  "scripts/mergeRecoilJets.sh"
  "scripts/root_in_analysis_env.sh"
)

usage() {
  cat <<'EOF'
Usage:
  ./scripts/sftp_push_recoiljets.sh <group-or-file> [group-or-file ...] [--commit-push -m "message"]

Uploads selected known files from the local Mac checkout to the SDCC analysis
checkout via interactive sftp. No password is stored; sftp prompts normally.

Groups:
  condor    RecoilJets_Condor*.sh submit/wrapper files
  macros    known pipeline macros/configs
  scripts   known pipeline helper scripts
  pipeline  all known transferable pipeline files except local-only sftp helpers

Examples:
  ./scripts/sftp_push_recoiljets.sh RecoilJets.cc RecoilJets.h
  ./scripts/sftp_push_recoiljets.sh condor
  ./scripts/sftp_push_recoiljets.sh mergeRecoilJets.sh analysis_config.yaml RecoilJets_AuAu.cc
  ./scripts/sftp_push_recoiljets.sh scripts/mergeRecoilJets.sh
  ./scripts/sftp_push_recoiljets.sh RecoilJets.cc RecoilJets.h RecoilJets_Condor_submit.sh --commit-push -m "Update PP recoil jet stitching diagnostics"

Options:
  --commit-push
      After a successful SFTP upload, run git add/commit/push for exactly the
      selected transferable local files. Never uses git add .
  -m, --message
      Commit message required with --commit-push.

Known files:
EOF
  local i
  for (( i=0; i<${#LOCAL_FILES[@]}; i++ )); do
    printf '  %-45s -> %s\n' "${LOCAL_FILES[$i]}" "${REMOTE_FILES[$i]}"
  done
}

die_with_usage() {
  echo "[ERROR] $*" >&2
  echo >&2
  usage >&2
  exit 2
}

normalize_arg() {
  local x="$1"
  x="${x#./}"
  echo "$x"
}

run_commit_push() {
  echo
  echo "Git commit/push requested."
  echo "Repository : ${LOCAL_BASE}"
  echo "Message    : ${commit_message}"
  echo
  echo "Current git status:"
  (cd "$LOCAL_BASE" && git status --short)
  echo
  echo "Staging exactly the uploaded local files:"
  local f
  for f in "${selected_local[@]}"; do
    echo "  $f"
  done

  (cd "$LOCAL_BASE" && git add -- "${selected_local[@]}")

  if (cd "$LOCAL_BASE" && git diff --cached --quiet -- "${selected_local[@]}"); then
    echo
    echo "[OK] No staged changes in selected uploaded files; skipping git commit/push."
    return 0
  fi

  echo
  echo "Staged diff summary:"
  (cd "$LOCAL_BASE" && git diff --cached --stat -- "${selected_local[@]}")
  echo
  (cd "$LOCAL_BASE" && git commit -m "$commit_message" -- "${selected_local[@]}")
  (cd "$LOCAL_BASE" && git push)
  echo
  echo "[OK] Git commit and push complete."
}

find_index_for_local() {
  local rel="$1"
  local i
  for (( i=0; i<${#LOCAL_FILES[@]}; i++ )); do
    if [[ "${LOCAL_FILES[$i]}" == "$rel" ]]; then
      echo "$i"
      return 0
    fi
  done
  return 1
}

selected_local=()
selected_remote=()
selection_args=()
commit_push=0
commit_message=""

add_index() {
  local idx="$1"
  local rel="${LOCAL_FILES[$idx]}"
  local remote="${REMOTE_FILES[$idx]}"
  local existing
  for existing in "${selected_local[@]:-}"; do
    if [[ "$existing" == "$rel" ]]; then
      return 0
    fi
  done
  selected_local+=( "$rel" )
  selected_remote+=( "$remote" )
}

add_local_rel() {
  local rel="$1"
  local idx
  idx="$(find_index_for_local "$rel")" || die_with_usage "Unknown file in group definition: $rel"
  add_index "$idx"
}

add_group() {
  local group="$1"
  local f
  case "$group" in
    condor)
      for f in "${GROUP_CONDOR[@]}"; do add_local_rel "$f"; done
      ;;
    macros)
      for f in "${GROUP_MACROS[@]}"; do add_local_rel "$f"; done
      ;;
    scripts)
      for f in "${GROUP_SCRIPTS[@]}"; do add_local_rel "$f"; done
      ;;
    pipeline)
      add_local_rel "src/RecoilJets.cc"
      add_local_rel "src/RecoilJets.h"
      add_local_rel "src_AuAu/RecoilJets_AuAu.cc"
      add_local_rel "src_AuAu/RecoilJets_AuAu.h"
      add_group condor
      add_group macros
      add_group scripts
      ;;
    *)
      return 1
      ;;
  esac
}

resolve_file_arg() {
  local arg
  arg="$(normalize_arg "$1")"

  local idx
  if idx="$(find_index_for_local "$arg")"; then
    add_index "$idx"
    return 0
  fi

  if [[ "$arg" == */* ]]; then
    die_with_usage "Unknown relative path: $arg"
  fi

  local matches=()
  local i
  for (( i=0; i<${#LOCAL_FILES[@]}; i++ )); do
    if [[ "$(basename "${LOCAL_FILES[$i]}")" == "$arg" ]]; then
      matches+=( "$i" )
    fi
  done

  if (( ${#matches[@]} == 0 )); then
    die_with_usage "Unknown group, relative path, or basename: $arg"
  fi
  if (( ${#matches[@]} > 1 )); then
    echo "[ERROR] Ambiguous basename: $arg" >&2
    echo "Use one of these relative paths:" >&2
    for idx in "${matches[@]}"; do
      echo "  ${LOCAL_FILES[$idx]}" >&2
    done
    exit 2
  fi

  add_index "${matches[0]}"
}

if (( $# == 0 )); then
  usage
  exit 2
fi

case "${1:-}" in
  -h|--help|help)
    usage
    exit 0
    ;;
esac

while (( $# > 0 )); do
  case "$1" in
    --commit-push)
      commit_push=1
      shift
      ;;
    -m|--message)
      opt="$1"
      shift
      if (( $# == 0 )); then
        die_with_usage "$opt requires a commit message"
      fi
      commit_message="$1"
      shift
      ;;
    --)
      shift
      while (( $# > 0 )); do
        selection_args+=( "$1" )
        shift
      done
      ;;
    -*)
      die_with_usage "Unknown option: $1"
      ;;
    *)
      selection_args+=( "$1" )
      shift
      ;;
  esac
done

if (( commit_push )) && [[ -z "$commit_message" ]]; then
  die_with_usage "--commit-push requires -m \"commit message\""
fi

if (( ${#selection_args[@]} == 0 )); then
  die_with_usage "No groups or files selected."
fi

for arg in "${selection_args[@]}"; do
  arg="$(normalize_arg "$arg")"
  if add_group "$arg"; then
    continue
  fi
  resolve_file_arg "$arg"
done

if (( ${#selected_local[@]} == 0 )); then
  die_with_usage "No files selected."
fi

if [[ ! -d "$LOCAL_BASE" ]]; then
  echo "[ERROR] Local base does not exist: $LOCAL_BASE" >&2
  exit 3
fi

for f in "${selected_local[@]}"; do
  if [[ ! -f "${LOCAL_BASE}/${f}" ]]; then
    echo "[ERROR] Missing local file: ${LOCAL_BASE}/${f}" >&2
    exit 4
  fi
done

echo
echo "Remote host : ${REMOTE_HOST}"
echo "Local base  : ${LOCAL_BASE}"
echo "Remote base : ${REMOTE_BASE}"
echo
echo "Files to upload:"
for (( i=0; i<${#selected_local[@]}; i++ )); do
  echo "  ${LOCAL_BASE}/${selected_local[$i]}"
  echo "    -> ${REMOTE_BASE}/${selected_remote[$i]}"
done
echo
echo "This will overwrite the remote files listed above."
read -r -p "Continue? [y/N]: " confirm
case "$confirm" in
  y|Y|yes|YES|Yes)
    ;;
  *)
    echo "Aborted."
    exit 0
    ;;
esac

batch_file="$(mktemp "${TMPDIR:-/tmp}/sftp_push_recoiljets.XXXXXX")"
cleanup() {
  rm -f "$batch_file"
}
trap cleanup EXIT

{
  printf 'lcd %s\n' "$LOCAL_BASE"
  printf 'cd %s\n' "$REMOTE_BASE"
  for (( i=0; i<${#selected_local[@]}; i++ )); do
    printf 'put %s %s\n' "${selected_local[$i]}" "${selected_remote[$i]}"
  done
} > "$batch_file"

echo
echo "sftp batch commands:"
sed 's/^/  /' "$batch_file"
echo
echo "Opening interactive sftp. Enter your SDCC password when prompted."
if sftp \
    -oBatchMode=no \
    -oPreferredAuthentications=password,keyboard-interactive,publickey \
    -b "$batch_file" \
    "$REMOTE_HOST"; then
  echo
  echo "[OK] Upload complete."
else
  status=$?
  echo
  echo "[ERROR] sftp upload failed with exit code ${status}." >&2
  echo "[ERROR] No success confirmation was received from sftp." >&2
  exit "$status"
fi

echo
echo "Uploaded ${#selected_local[@]} file(s)."

if (( commit_push )); then
  run_commit_push
fi
