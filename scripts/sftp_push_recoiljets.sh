#!/usr/bin/env bash
set -euo pipefail

LOCAL_BASE="/Users/patsfan753/Desktop/ThesisAnalysis"
REMOTE_BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
REMOTE_HOST="patsfan753@sftp.sdcc.bnl.gov"

LOCAL_FILES=(
  "scripts/audit_auau_grl_projection.sh"
  "scripts/audit_auau_ml_training_smoke.py"
  "scripts/auau_tight_bdt_pipeline.sh"
  "scripts/estimateEmbeddedPhotonXsec.sh"
  "scripts/make_dstListsData.sh"
  "scripts/makeThesisSimLists.sh"
  "scripts/mergeRecoilJets.sh"
  "scripts/recoiljets_cleanup.sh"
  "scripts/root_in_analysis_env.sh"
  "scripts/RecoilJets_Condor_AuAu.sh"
  "scripts/RecoilJets_Condor_submit.sh"
  "scripts/RecoilJets_Condor.sh"
  "scripts/train_auau_jet_residual_bdt.py"
  "scripts/train_auau_photon_bdt.py"
  "scripts/validate_auau_tight_bdt_on_sim.py"
  "macros/analysis_config.yaml"
  "macros/Calo_Calib.C"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_auauTightBDTTraining.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/PrintPPStitchDiagnostics.C"
  "src/PhotonClusterBuilder.cc"
  "src/PhotonClusterBuilder.h"
  "src/RecoilJets.cc"
  "src/RecoilJets.h"
  "src_AuAu/RecoilJets_AuAu.cc"
  "src_AuAu/RecoilJets_AuAu.h"
)

REMOTE_FILES=(
  "scripts/audit_auau_grl_projection.sh"
  "scripts/audit_auau_ml_training_smoke.py"
  "scripts/auau_tight_bdt_pipeline.sh"
  "scripts/estimateEmbeddedPhotonXsec.sh"
  "scripts/make_dstListsData.sh"
  "scripts/makeThesisSimLists.sh"
  "scripts/mergeRecoilJets.sh"
  "scripts/recoiljets_cleanup.sh"
  "scripts/root_in_analysis_env.sh"
  "RecoilJets_Condor_AuAu.sh"
  "RecoilJets_Condor_submit.sh"
  "RecoilJets_Condor.sh"
  "scripts/train_auau_jet_residual_bdt.py"
  "scripts/train_auau_photon_bdt.py"
  "scripts/validate_auau_tight_bdt_on_sim.py"
  "macros/analysis_config.yaml"
  "macros/Calo_Calib.C"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_auauTightBDTTraining.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/PrintPPStitchDiagnostics.C"
  "coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.cc"
  "coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"
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
  "macros/Calo_Calib.C"
  "macros/Fun4All_recoilJets.C"
  "macros/Fun4All_recoilJets_AuAu.C"
  "macros/Fun4All_auauTightBDTTraining.C"
  "macros/Fun4All_recoilJets_unified_impl.C"
  "macros/PrintPPStitchDiagnostics.C"
)

GROUP_SCRIPTS=(
  "scripts/audit_auau_grl_projection.sh"
  "scripts/audit_auau_ml_training_smoke.py"
  "scripts/auau_tight_bdt_pipeline.sh"
  "scripts/estimateEmbeddedPhotonXsec.sh"
  "scripts/make_dstListsData.sh"
  "scripts/makeThesisSimLists.sh"
  "scripts/mergeRecoilJets.sh"
  "scripts/recoiljets_cleanup.sh"
  "scripts/root_in_analysis_env.sh"
  "scripts/train_auau_jet_residual_bdt.py"
  "scripts/train_auau_photon_bdt.py"
  "scripts/validate_auau_tight_bdt_on_sim.py"
)

usage() {
  cat <<'EOF'
Usage:
  ./scripts/sftp_push_recoiljets.sh <group-or-file> [group-or-file ...] [--commit-push -m "message"]
  ./scripts/sftp_push_recoiljets.sh status [group-or-file ...]
  ./scripts/sftp_push_recoiljets.sh diff [group-or-file ...]

Uploads selected known files from the local Mac checkout to the SDCC analysis
checkout via interactive sftp. No password is stored; sftp prompts normally.
Upload mode prints the overwrite preview and then starts sftp directly; it does
not ask for an extra y/N confirmation.

Read-only modes:
  status
      Fetch selected mapped remote files into a temp directory and print whether
      each SDCC file matches the local file. With no selection, checks pipeline.
  diff
      Fetch selected mapped remote files and print unified diffs for files that
      differ. With no selection, checks pipeline.

Groups:
  condor    RecoilJets_Condor*.sh submit/wrapper files
  macros    known pipeline macros/configs
  scripts   known pipeline helper scripts
  pipeline  all known transferable pipeline files except local-only sftp helpers
  all       alias for pipeline
  changed   all changed known transferable files in this checkout

Examples:
  ./scripts/sftp_push_recoiljets.sh status
  ./scripts/sftp_push_recoiljets.sh status changed
  ./scripts/sftp_push_recoiljets.sh status pipeline
  ./scripts/sftp_push_recoiljets.sh diff RecoilJets_Condor_submit.sh
  ./scripts/sftp_push_recoiljets.sh RecoilJets.cc RecoilJets.h
  ./scripts/sftp_push_recoiljets.sh condor
  ./scripts/sftp_push_recoiljets.sh changed
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

add_changed_known_files() {
  local line path idx any=0
  while IFS= read -r line; do
    [[ -n "$line" ]] || continue
    path="${line:3}"
    if [[ "$path" == *" -> "* ]]; then
      path="${path##* -> }"
    fi
    path="${path#\"}"
    path="${path%\"}"
    path="$(normalize_arg "$path")"
    if idx="$(find_index_for_local "$path")"; then
      add_index "$idx"
      any=1
    fi
  done < <(cd "$LOCAL_BASE" && git status --porcelain --untracked-files=all)

  (( any )) || die_with_usage "No changed known transferable files found in ${LOCAL_BASE}."
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
    pipeline|all)
      add_local_rel "src/PhotonClusterBuilder.cc"
      add_local_rel "src/PhotonClusterBuilder.h"
      add_local_rel "src/RecoilJets.cc"
      add_local_rel "src/RecoilJets.h"
      add_local_rel "src_AuAu/RecoilJets_AuAu.cc"
      add_local_rel "src_AuAu/RecoilJets_AuAu.h"
      add_group condor
      add_group macros
      add_group scripts
      ;;
    changed|changed-pipeline|changedPipeline|changed-known)
      add_changed_known_files
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

hash_file() {
  local file="$1"
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$file" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$file" | awk '{print $1}'
  else
    echo "[ERROR] Need shasum or sha256sum for status/diff mode." >&2
    exit 5
  fi
}

run_remote_compare() {
  local compare_mode="$1"
  local tmp_dir
  tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/sftp_push_recoiljets_check.XXXXXX")"
  local batch
  batch="${tmp_dir}/fetch.sftp"

  {
    printf 'cd %s\n' "$REMOTE_BASE"
    local i
    for (( i=0; i<${#selected_remote[@]}; i++ )); do
      printf -- '-get %s %s\n' "${selected_remote[$i]}" "${tmp_dir}/remote_${i}"
    done
  } > "$batch"

  echo
  echo "Remote host : ${REMOTE_HOST}"
  echo "Local base  : ${LOCAL_BASE}"
  echo "Remote base : ${REMOTE_BASE}"
  echo
  echo "Read-only ${compare_mode}: fetching ${#selected_local[@]} mapped file(s) into a temp directory."
  echo "No remote files will be modified."
  echo "Opening interactive sftp. Enter your SDCC password when prompted."

  if sftp \
      -oBatchMode=no \
      -oPreferredAuthentications=password,keyboard-interactive,publickey \
      -b "$batch" \
      "$REMOTE_HOST"; then
    :
  else
    local status=$?
    rm -rf "$tmp_dir"
    echo
    echo "[ERROR] sftp fetch failed with exit code ${status}." >&2
    exit "$status"
  fi

  echo
  printf '%-8s  %-45s  %s\n' "STATUS" "LOCAL" "REMOTE"
  printf '%-8s  %-45s  %s\n' "--------" "---------------------------------------------" "---------------------------------------------"

  local matches=0 differs=0 missing=0 checked=0
  local i local_file remote_copy local_hash remote_hash result
  for (( i=0; i<${#selected_local[@]}; i++ )); do
    local_file="${LOCAL_BASE}/${selected_local[$i]}"
    remote_copy="${tmp_dir}/remote_${i}"
    checked=$(( checked + 1 ))
    if [[ ! -f "$remote_copy" ]]; then
      result="MISSING"
      missing=$(( missing + 1 ))
    else
      local_hash="$(hash_file "$local_file")"
      remote_hash="$(hash_file "$remote_copy")"
      if [[ "$local_hash" == "$remote_hash" ]]; then
        result="MATCH"
        matches=$(( matches + 1 ))
      else
        result="DIFFER"
        differs=$(( differs + 1 ))
      fi
    fi
    printf '%-8s  %-45s  %s\n' "$result" "${selected_local[$i]}" "${selected_remote[$i]}"

    if [[ "$compare_mode" == "diff" && "$result" == "DIFFER" ]]; then
      echo
      echo "Diff for ${selected_local[$i]}:"
      diff -u "$local_file" "$remote_copy" | sed \
        -e "1s|.*|--- local:${selected_local[$i]}|" \
        -e "2s|.*|+++ sdcc:${selected_remote[$i]}|" || true
      echo
    elif [[ "$compare_mode" == "diff" && "$result" == "MISSING" ]]; then
      echo
      echo "Remote file missing for ${selected_local[$i]} -> ${selected_remote[$i]}"
      echo
    fi
  done

  rm -rf "$tmp_dir"
  echo
  echo "Summary: checked=${checked} match=${matches} differ=${differs} missing=${missing}"
  if (( differs > 0 || missing > 0 )); then
    exit 1
  fi
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

mode="upload"
case "${1:-}" in
  status|check)
    mode="status"
    shift
    ;;
  diff)
    mode="diff"
    shift
    ;;
  upload|push)
    mode="upload"
    shift
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

if [[ "$mode" != "upload" && "$commit_push" -eq 1 ]]; then
  die_with_usage "--commit-push is only valid for upload mode"
fi

if (( commit_push )) && [[ -z "$commit_message" ]]; then
  die_with_usage "--commit-push requires -m \"commit message\""
fi

if (( ${#selection_args[@]} == 0 )) && [[ "$mode" != "upload" ]]; then
  selection_args+=( "pipeline" )
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

if [[ "$mode" == "status" || "$mode" == "diff" ]]; then
  run_remote_compare "$mode"
  exit 0
fi

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
echo "Proceeding without an extra y/N prompt."

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
