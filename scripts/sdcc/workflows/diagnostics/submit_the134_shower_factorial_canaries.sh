#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd -P)"
cd "$repo_root"

requested_mode="${1:-preflight}"
sidecar_only_mode=0
sidecar_pp_replacement_mode=0
mode="$requested_mode"
case "$requested_mode" in
  sidecar-preflight|sidecar-submit|sidecar-status|sidecar-validate)
    sidecar_only_mode=1
    mode="${requested_mode#sidecar-}"
    ;;
  sidecar-pp-replacement-preflight|sidecar-pp-replacement-submit|sidecar-pp-replacement-status|sidecar-pp-replacement-validate|sidecar-pp-replacement-aggregate)
    sidecar_only_mode=1
    sidecar_pp_replacement_mode=1
    mode="${requested_mode#sidecar-pp-replacement-}"
    ;;
esac

require_replacement_namespaces() {
  local operation="$1"
  : "${RJ_THE134_TAG:?set a fresh explicit RJ_THE134_TAG for the p+p replacement pair}"
  : "${RJ_THE134_OUTPUT_ROOT:?set a fresh explicit RJ_THE134_OUTPUT_ROOT for the p+p replacement pair}"
  : "${RJ_THE134_EVIDENCE_ROOT:?set a fresh explicit RJ_THE134_EVIDENCE_ROOT for the p+p replacement pair}"
  [[ "$RJ_THE134_TAG" =~ ^[A-Za-z0-9._-]+$ ]] ||
    { echo "[THE134][ERROR] replacement tag must be one path-safe component" >&2; return 2; }
  [[ "$RJ_THE134_TAG" != "." && "$RJ_THE134_TAG" != ".." ]] ||
    { echo "[THE134][ERROR] replacement tag cannot be dot or dot-dot" >&2; return 2; }
  [[ "$RJ_THE134_OUTPUT_ROOT" == /* && "$RJ_THE134_EVIDENCE_ROOT" == /* ]] ||
    { echo "[THE134][ERROR] replacement output and evidence roots must be absolute" >&2; return 2; }
  [[ "$(basename -- "$RJ_THE134_OUTPUT_ROOT")" == "$RJ_THE134_TAG" ]] ||
    { echo "[THE134][ERROR] replacement output namespace must end in the exact tag" >&2; return 2; }
  [[ "$(basename -- "$RJ_THE134_EVIDENCE_ROOT")" == "$RJ_THE134_TAG" ]] ||
    { echo "[THE134][ERROR] replacement evidence namespace must end in the exact tag" >&2; return 2; }
  [[ "$RJ_THE134_OUTPUT_ROOT" != "$RJ_THE134_EVIDENCE_ROOT" ]] ||
    { echo "[THE134][ERROR] replacement output and evidence namespaces must differ" >&2; return 2; }
  local canonical_output canonical_evidence
  canonical_output="$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$RJ_THE134_OUTPUT_ROOT")" ||
    { echo "[THE134][ERROR] could not canonicalize replacement output namespace" >&2; return 2; }
  canonical_evidence="$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$RJ_THE134_EVIDENCE_ROOT")" ||
    { echo "[THE134][ERROR] could not canonicalize replacement evidence namespace" >&2; return 2; }
  [[ "$canonical_output" != "$canonical_evidence" ]] ||
    { echo "[THE134][ERROR] replacement namespaces canonicalize to one path" >&2; return 2; }
  [[ "$(basename -- "$canonical_output")" == "$RJ_THE134_TAG" &&
     "$(basename -- "$canonical_evidence")" == "$RJ_THE134_TAG" ]] ||
    { echo "[THE134][ERROR] canonical replacement namespaces must end in the exact tag" >&2; return 2; }
  case "$operation" in
    preflight|submit)
      [[ ! -e "$RJ_THE134_OUTPUT_ROOT" ]] ||
        { echo "[THE134][ERROR] replacement output namespace already exists" >&2; return 2; }
      [[ ! -e "$RJ_THE134_EVIDENCE_ROOT" ]] ||
        { echo "[THE134][ERROR] replacement evidence namespace already exists" >&2; return 2; }
      ;;
  esac
}

build_sidecar_only_keys() {
  local include_auau="$1"
  local arm
  for arm in direct writer; do
    printf '%s\n' "${arm}:pp_inclusive_sim:run28_jet8"
    if [[ "$include_auau" == 1 ]]; then
      printf '%s\n' "${arm}:auau_inclusive_embedded:run28_embeddedJet12"
    fi
  done
}

require_exact_sidecar_keys() {
  local observed="$1"
  local expected="$2"
  local label="$3"
  [[ "$observed" == "$expected" ]] ||
    { echo "[THE134][ERROR] ${label} owns exactly its frozen row set" >&2; return 2; }
}

if (( sidecar_pp_replacement_mode )); then
  require_replacement_namespaces "$mode"
  tag="$RJ_THE134_TAG"
elif (( sidecar_only_mode )); then
  tag="${RJ_THE134_TAG:-the134_h70_sidecar_differential_canary_20260728_v1}"
else
  tag="${RJ_THE134_TAG:-the134_h70_shower_factorial_replay_canary_20260722_v1}"
fi
build_root="${RJ_THE134_BUILD_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/.recoiljets_tmp/the134_h70_schema_v2_20260722_r3}"
pp_library="${RJ_THE134_PP_LIBRARY:-${build_root}/install/pp/lib/libRecoilJets.so}"
auau_library="${RJ_THE134_AUAU_LIBRARY:-${build_root}/install/auau/lib/libRecoilJetsAuAu.so}"
: "${RJ_THE134_PP_LIBRARY_SHA256:?set the frozen p+p library SHA-256}"
: "${RJ_THE134_AUAU_LIBRARY_SHA256:?set the frozen Au+Au library SHA-256}"

factorial_semantic_sha="$({
  sha256sum \
    src/RJReplayFoundationV1.h \
    src/RJReplayRuntimeV1.h \
    src/RJShowerFactorialV1.h
  printf '%s\n' 'THE134|H70,H0,G70,G0,O70,O0,R70|nominal-domain=15<=ET<35|schema=2'
} | sha256sum | awk '{print $1}')"
if [[ ! "$factorial_semantic_sha" =~ ^[0-9a-f]{64}$ ]]; then
  echo "[THE134][ERROR] could not derive the factorial semantic identity" >&2
  exit 2
fi

# One source witness per detector lane, each with a direct and writer arm.  The
# completed THE-119 matrix remains the base replay-infrastructure certificate;
# this bounded matrix certifies only the new multi-view shower payload.  A
# replacement mode owns only the invalidated p+p pair and cannot widen its
# selector to duplicate the preserved Au+Au rows.
if (( sidecar_pp_replacement_mode )); then
  default_only_keys="$(build_sidecar_only_keys 0 | paste -sd, -)"
elif (( sidecar_only_mode )); then
  default_only_keys="$(build_sidecar_only_keys 1 | paste -sd, -)"
else
  default_only_keys="$({
    for arm in direct writer; do
      printf '%s\n' \
        "${arm}:pp_data:pp_data" \
        "${arm}:pp_photon_sim:run28_photonjet20" \
        "${arm}:pp_inclusive_sim:run28_jet8" \
        "${arm}:auau_data:auau_data" \
        "${arm}:auau_photon_embedded:run28_embeddedPhoton12" \
        "${arm}:auau_inclusive_embedded:run28_embeddedJet12"
    done
  } | paste -sd, -)"
fi
only_keys="${RJ_THE134_ONLY_KEYS:-$default_only_keys}"
if (( sidecar_pp_replacement_mode )); then
  require_exact_sidecar_keys \
    "$only_keys" "$default_only_keys" \
    "p+p replacement direct+writer pp_inclusive_sim:run28_jet8"
fi
if (( sidecar_only_mode && ! sidecar_pp_replacement_mode )); then
  require_exact_sidecar_keys \
    "$only_keys" "$default_only_keys" \
    "sidecar differential mode"
fi

export RJ_CODEX_CHAT_NAME="THE-114+THE-134 | pp/AuAu H70 Production Gate"
export RJ_CODEX_THREAD_ID="019f80b5-dc56-7330-9ee7-56ef417547dc"
export RJ_THE119_TAG="$tag"
export RJ_THE119_OUTPUT_ROOT="${RJ_THE134_OUTPUT_ROOT:-/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/smoke/replay_foundation/${tag}}"
export RJ_THE119_EVIDENCE_ROOT="${RJ_THE134_EVIDENCE_ROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/evidence/qa/${tag}}"
export RJ_THE119_PP_LIBRARY="$pp_library"
export RJ_THE119_AUAU_LIBRARY="$auau_library"
export RJ_THE119_PP_LIBRARY_SHA256="$RJ_THE134_PP_LIBRARY_SHA256"
export RJ_THE119_AUAU_LIBRARY_SHA256="$RJ_THE134_AUAU_LIBRARY_SHA256"
export RJ_THE119_PP_MODEL_SHOWER_DEFINITION=H70
export RJ_THE119_PP_REFERENCE_MODEL_SHOWER_DEFINITION=H70
export RJ_THE119_AUAU_MODEL_SHOWER_DEFINITION=H0
export RJ_THE119_SEMANTIC_SHA256="$factorial_semantic_sha"
export RJ_THE119_ONLY_KEYS="$only_keys"
export RJ_THE119_NEVENTS="${RJ_THE134_NEVENTS:-3000}"
export RJ_THE119_PHOTON_CAPTURE_ET_MIN=5.0
export RJ_THE119_JET_CONSTITUENT_PT_MIN=5.0
export RJ_THE119_PP_CAPTURE_WITNESS_QA=1
export RJ_THE119_CODE_COMMIT="${RJ_THE134_CODE_COMMIT:-$(git rev-parse HEAD)}"
if (( sidecar_only_mode )); then
  export RJ_THE119_EXTRA_ENV_TEMPLATE=''
  export RJ_THE119_PP_WRITER_EXTRA_ENV_TEMPLATE=''
  export RJ_THE119_AUAU_WRITER_EXTRA_ENV_TEMPLATE=''
  export RJ_THE119_WRITER_EXTRA_ENV_TEMPLATE='RJ_THE134_MULTIVIEW_TRAINING_V1=1;RJ_THE134_MULTIVIEW_TRAINING_FILE=__OUTPUT__/RJPhotonTrainingViewV1.root;RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1;RJ_THE134_EXPECTED_SOURCE_ROLE=background'
  export RJ_THE119_PP_EXTRA_ENV_TEMPLATE='RJ_REPLAY_PERIOD=0mrad;RJ_REPLAY_SI_DI_ROLE=SI;RJ_REPLAY_OWNERSHIP_STATE=source_role_frozen;RJ_SIM_ALLOW_NONE_LISTS=0;RJ_PPG12_PHOTON_YIELD=1;RJ_PPG12_PHOTON_YIELD_DOUBLE=0;RJ_PPG12_PERIOD=0mrad;RJ_PPG12_PERIOD_USE_LUMI_WEIGHT=1;RJ_PPG12_PERIOD_STRICT_DI=0;RJ_PPG12_PERIOD_ALLOW_ALL_SIM=0;RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE=0;RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE=0;RJ_PP_PHOTONID_EXTRACT_ONLY=1;RJ_PP_PHOTONID_TRAINING_TREE=1;RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES=0;RJ_PP_PHOTONID_SOURCE_ROLE=background;RJ_PP_PHOTONID_PPG12_FILTER=1;RJ_PP_PHOTONID_REQUIRE_PRESELECTION=0;RJ_PINNED_RELEASE_NAME=ana.560;RJ_PINNED_OFFLINE_MAIN=/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.560'
  export RJ_THE119_AUAU_EXTRA_ENV_TEMPLATE='RJ_REPLAY_PERIOD=AUAU_RUN24;RJ_REPLAY_SI_DI_ROLE=EMBEDDED;RJ_REPLAY_OWNERSHIP_STATE=source_role_frozen;RJ_SIM_ALLOW_NONE_LISTS=1;RJ_AUAU_BDT_EXTRACT_ONLY=1;RJ_AUAU_BDT_TRAINING_TREE=1;RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES=0;RJ_AUAU_BDT_NPB_DATA_TAGGING=0;RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=1;RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION=0;RJ_AUAU_USE_TOPOCLUSTER_ISOLATION=0'
fi

validate_terminal_rows() {
  local initial_queue="${RJ_THE119_EVIDENCE_ROOT}/initial_queue.tsv"
  local preflight="${RJ_THE119_EVIDENCE_ROOT}/preflight_receipt.txt"
  local terminal_receipt="${RJ_THE119_EVIDENCE_ROOT}/terminal_gate_receipt.json"
  local terminal_rows="${RJ_THE119_EVIDENCE_ROOT}/.terminal_gate_rows.$$.tsv"
  local terminal_receipt_tmp="${terminal_receipt}.tmp.$$"
  local expected_rows=4
  local cluster proc cluster_proc queue_state history_state job_status exit_code
  local role seen_identities='|' seen_roles='|' initial_queue_sha preflight_sha
  (( sidecar_pp_replacement_mode )) && expected_rows=2
  [[ -s "$initial_queue" ]] ||
    { echo "[THE134][ERROR] missing exact initial queue receipt: $initial_queue" >&2; return 2; }
  [[ "$(wc -l < "$initial_queue" | tr -d ' ')" == "$expected_rows" ]] ||
    { echo "[THE134][ERROR] sidecar canary must own exactly ${expected_rows} rows" >&2; return 2; }
  (( sidecar_pp_replacement_mode )) && : > "$terminal_rows"
  while read -r cluster proc _status _args; do
    [[ "$cluster" =~ ^[0-9]+$ && "$proc" =~ ^[0-9]+$ ]] ||
      { echo "[THE134][ERROR] malformed exact Condor identity in $initial_queue" >&2; return 2; }
    cluster_proc="${cluster}.${proc}"
    [[ "$seen_identities" != *"|${cluster_proc}|"* ]] ||
      { echo "[THE134][ERROR] duplicate exact Condor identity: $cluster_proc" >&2; return 2; }
    seen_identities="${seen_identities}${cluster_proc}|"
    queue_state="$(condor_q "$cluster_proc" -af JobStatus HoldReason NumJobStarts 2>/dev/null || true)"
    [[ -z "$queue_state" ]] ||
      { echo "[THE134][ERROR] row is not terminal: $cluster_proc $queue_state" >&2; return 2; }
    history_state="$(condor_history "$cluster_proc" -limit 1 -af JobStatus ExitCode 2>/dev/null || true)"
    read -r job_status exit_code <<< "$history_state"
    [[ "$job_status" == 4 && "$exit_code" == 0 ]] ||
      { echo "[THE134][ERROR] terminal gate failed: $cluster_proc ${history_state:-NOT_FOUND}" >&2; return 2; }
    if (( sidecar_pp_replacement_mode )); then
      role="$(
        python3 - "$RJ_THE119_OUTPUT_ROOT" "$_args" <<'PY'
import shlex
import sys

base = sys.argv[1].rstrip("/")
tokens = shlex.split(sys.argv[2])
candidates = [
    ("direct", f"{base}/direct/pp_inclusive_sim/run28_jet8"),
    ("writer", f"{base}/writer/pp_inclusive_sim/run28_jet8"),
]
matches = [
    (index, role)
    for index, token in enumerate(tokens)
    for role, expected in candidates
    if token == expected
]
if len(matches) != 1 or matches[0][0] != len(tokens) - 1:
    raise SystemExit(2)
print(matches[0][1])
PY
      )" || {
        echo "[THE134][ERROR] replacement row lacks one exact final role/path token: $cluster_proc" >&2
        return 2
      }
      [[ "$seen_roles" != *"|${role}|"* ]] ||
        { echo "[THE134][ERROR] duplicate replacement role: $role" >&2; return 2; }
      seen_roles="${seen_roles}${role}|"
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$cluster" "$proc" "$cluster_proc" "$role" "$job_status" "$exit_code" \
        >> "$terminal_rows"
    fi
  done < "$initial_queue"
  if (( sidecar_pp_replacement_mode )); then
    [[ "$seen_roles" == *"|direct|"* && "$seen_roles" == *"|writer|"* ]] ||
      { echo "[THE134][ERROR] replacement terminal gate lacks exact direct/writer roles" >&2; return 2; }
    [[ -s "$preflight" ]] ||
      { echo "[THE134][ERROR] missing p+p replacement preflight receipt: $preflight" >&2; return 2; }
    if command -v sha256sum >/dev/null 2>&1; then
      initial_queue_sha="$(sha256sum "$initial_queue" | awk '{print $1}')"
      preflight_sha="$(sha256sum "$preflight" | awk '{print $1}')"
    else
      initial_queue_sha="$(shasum -a 256 "$initial_queue" | awk '{print $1}')"
      preflight_sha="$(shasum -a 256 "$preflight" | awk '{print $1}')"
    fi
    python3 - \
      "$terminal_receipt_tmp" "$initial_queue" "$initial_queue_sha" \
      "$preflight" "$preflight_sha" "$RJ_THE119_TAG" \
      "$RJ_THE119_OUTPUT_ROOT" "$terminal_rows" <<'PY'
import json
import pathlib
import sys

(
    receipt_tmp,
    initial_queue,
    initial_queue_sha,
    preflight,
    preflight_sha,
    tag,
    output_root,
    rows_path,
) = sys.argv[1:]
rows = []
for raw in pathlib.Path(rows_path).read_text(encoding="utf-8").splitlines():
    cluster, proc, cluster_proc, role, job_status, exit_code = raw.split("\t")
    rows.append(
        dict(
            cluster_id=int(cluster),
            proc_id=int(proc),
            cluster_proc=cluster_proc,
            role=role,
            lane="pp_inclusive_sim",
            sample="run28_jet8",
            job_status=int(job_status),
            exit_code=int(exit_code),
        )
    )
payload = dict(
    schema="THE134_PP_REPLACEMENT_TERMINAL_GATE_V1",
    tag=tag,
    output_root=str(pathlib.Path(output_root).resolve()),
    preflight_receipt=str(pathlib.Path(preflight).resolve()),
    preflight_receipt_sha256=preflight_sha,
    initial_queue_tsv=str(pathlib.Path(initial_queue).resolve()),
    initial_queue_tsv_sha256=initial_queue_sha,
    row_count=len(rows),
    rows=rows,
)
pathlib.Path(receipt_tmp).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY
    mv -f -- "$terminal_receipt_tmp" "$terminal_receipt"
    rm -f -- "$terminal_rows"
  fi
}

case "$mode" in
  preflight|submit|status)
    exec scripts/sdcc/workflows/diagnostics/submit_the119_replay_foundation_canaries.sh "$mode"
    ;;
  validate)
    (( sidecar_only_mode )) ||
      { echo "[THE134][ERROR] validate is reserved for sidecar-* modes" >&2; exit 2; }
    validate_terminal_rows
    if (( sidecar_pp_replacement_mode )); then
      exec python3 scripts/diagnostics/replay_foundation/validate_the134_sidecar_differential.py \
        --component-system pp \
        --output-root "$RJ_THE119_OUTPUT_ROOT" \
        --preflight-receipt "$RJ_THE119_EVIDENCE_ROOT/preflight_receipt.txt" \
        --terminal-gate-receipt "$RJ_THE119_EVIDENCE_ROOT/terminal_gate_receipt.json" \
        --output-json "$RJ_THE119_EVIDENCE_ROOT/sidecar_pp_component_certificate.json"
    fi
    exec python3 scripts/diagnostics/replay_foundation/validate_the134_sidecar_differential.py \
      --output-root "$RJ_THE119_OUTPUT_ROOT" \
      --preflight-receipt "$RJ_THE119_EVIDENCE_ROOT/preflight_receipt.txt" \
      --output-json "$RJ_THE119_EVIDENCE_ROOT/sidecar_differential_certificate.json"
    ;;
  aggregate)
    (( sidecar_pp_replacement_mode )) ||
      { echo "[THE134][ERROR] aggregate is reserved for sidecar-pp-replacement-* mode" >&2; exit 2; }
    : "${RJ_THE134_PRESERVED_AUAU_OUTPUT_ROOT:?set the preserved four-row output root}"
    : "${RJ_THE134_PRESERVED_AUAU_PREFLIGHT_RECEIPT:?set the preserved four-row preflight receipt}"
    exec python3 scripts/diagnostics/replay_foundation/validate_the134_sidecar_differential.py \
      --pp-component-certificate "${RJ_THE134_PP_COMPONENT_CERTIFICATE:-$RJ_THE119_EVIDENCE_ROOT/sidecar_pp_component_certificate.json}" \
      --auau-output-root "$RJ_THE134_PRESERVED_AUAU_OUTPUT_ROOT" \
      --auau-preflight-receipt "$RJ_THE134_PRESERVED_AUAU_PREFLIGHT_RECEIPT" \
      --output-json "$RJ_THE119_EVIDENCE_ROOT/sidecar_differential_aggregate_certificate.json"
    ;;
  *)
    echo "usage: $0 preflight|submit|status|sidecar-preflight|sidecar-submit|sidecar-status|sidecar-validate|sidecar-pp-replacement-preflight|sidecar-pp-replacement-submit|sidecar-pp-replacement-status|sidecar-pp-replacement-validate|sidecar-pp-replacement-aggregate" >&2
    exit 2
    ;;
esac
