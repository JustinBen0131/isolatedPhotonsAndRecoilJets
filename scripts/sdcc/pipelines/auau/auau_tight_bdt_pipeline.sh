#!/usr/bin/env bash
set -euo pipefail

RJ_REPO_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
readonly RJ_REPO_BASE
SIM_ROOT="${RJ_SIM_ROOT:-${RJ_REPO_BASE}/simListFiles}"
TRAIN_BASE="${RJ_AUAU_TIGHT_BDT_TRAIN_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining}"
LOCAL_BASE="${RJ_AUAU_TIGHT_BDT_LOCAL_BASE:-${RJ_REPO_BASE}/local_bdt_training_outputs}"
MODEL_BASE="${RJ_AUAU_BDT_MODEL_BASE:-${RJ_REPO_BASE}/bdt_models}"
MASTER_YAML="${RJ_AUAU_TIGHT_BDT_CONFIG_SRC:-${RJ_REPO_BASE}/macros/analysis_config.yaml}"
TRAIN_MACRO="${RJ_AUAU_TIGHT_BDT_MACRO:-${RJ_REPO_BASE}/macros/Fun4All_auauTightBDTTraining.C}"
TRAIN_SCRIPT="${RJ_AUAU_TIGHT_BDT_TRAIN_SCRIPT:-${RJ_REPO_BASE}/scripts/ml/training/train_auau_photon_bdt.py}"
VALIDATE_SCRIPT="${RJ_AUAU_TIGHT_BDT_VALIDATE_SCRIPT:-${RJ_REPO_BASE}/scripts/ml/validation/validate_auau_tight_bdt_on_sim.py}"
ML_PYTHON="${RJ_ML_PYTHON:-${ML_PYTHON:-python3}}"
NOTIFY_EMAILS="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"

# Canonical AuAu training uses only embedded events accepted by the offline
# MinimumBiasClassifier. Controls may opt out only with an explicit second
# switch so an unrestricted extraction cannot silently become the baseline.
require_embedded_minbias="${RJ_AUAU_TIGHT_BDT_REQUIRE_EMBEDDED_MINBIAS:-1}"
case "${require_embedded_minbias,,}" in
  1|true|yes|on) require_embedded_minbias=1 ;;
  0|false|no|off) require_embedded_minbias=0 ;;
  *) printf '[auauTightBDT][ERR] invalid RJ_AUAU_TIGHT_BDT_REQUIRE_EMBEDDED_MINBIAS=%q\n' "$require_embedded_minbias" >&2; exit 2 ;;
esac
if [[ "$require_embedded_minbias" == "0" &&
      "${RJ_AUAU_TIGHT_BDT_ALLOW_UNRESTRICTED_EMBEDDED_CONTROL:-0}" != "1" ]]; then
  printf '[auauTightBDT][ERR] canonical extraction requires MinimumBiasClassifier pass; set RJ_AUAU_TIGHT_BDT_ALLOW_UNRESTRICTED_EMBEDDED_CONTROL=1 only for a labeled control\n' >&2
  exit 2
fi
export RJ_AUAU_TIGHT_BDT_REQUIRE_EMBEDDED_MINBIAS="$require_embedded_minbias"
export RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER="$require_embedded_minbias"

DEFAULT_SIGNAL_SAMPLES=(run28_embeddedPhoton12 run28_embeddedPhoton20)
DEFAULT_BACKGROUND_SAMPLES=(run28_embeddedJet12 run28_embeddedJet20 run28_embeddedJet30 run28_embeddedJet40)
SIGNAL_SAMPLES=("${DEFAULT_SIGNAL_SAMPLES[@]}")
BACKGROUND_SAMPLES=("${DEFAULT_BACKGROUND_SAMPLES[@]}")
if [[ -n "${RJ_AUAU_TIGHT_BDT_SIGNAL_SAMPLES+x}" ]]; then
  SIGNAL_SAMPLES=()
  if [[ -n "${RJ_AUAU_TIGHT_BDT_SIGNAL_SAMPLES//[[:space:]]/}" ]]; then
    read -r -a SIGNAL_SAMPLES <<< "${RJ_AUAU_TIGHT_BDT_SIGNAL_SAMPLES}"
  fi
fi
if [[ -n "${RJ_AUAU_TIGHT_BDT_BACKGROUND_SAMPLES+x}" ]]; then
  BACKGROUND_SAMPLES=()
  if [[ -n "${RJ_AUAU_TIGHT_BDT_BACKGROUND_SAMPLES//[[:space:]]/}" ]]; then
    read -r -a BACKGROUND_SAMPLES <<< "${RJ_AUAU_TIGHT_BDT_BACKGROUND_SAMPLES}"
  fi
fi
ALL_SAMPLES=("${SIGNAL_SAMPLES[@]}" "${BACKGROUND_SAMPLES[@]}")

ts() { date +%Y%m%d_%H%M%S; }
say() { printf '\033[1;36m[auauTightBDT]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauTightBDT][WARN]\033[0m %s\n' "$*" >&2; }
err() { printf '\033[1;31m[auauTightBDT][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

require_codex_submission_provenance() {
  [[ -n "${RJ_CODEX_CHAT_NAME:-}" ]] || die "Condor submission requires RJ_CODEX_CHAT_NAME"
  [[ -n "${RJ_CODEX_THREAD_ID:-}" ]] || die "Condor submission requires RJ_CODEX_THREAD_ID"
}

join_by_comma() {
  local IFS=,
  printf '%s' "$*"
}

enforce_default_embedded_samples() {
  local expected_signal expected_background actual_signal actual_background
  expected_signal="$(join_by_comma "${DEFAULT_SIGNAL_SAMPLES[@]}")"
  expected_background="$(join_by_comma "${DEFAULT_BACKGROUND_SAMPLES[@]}")"
  actual_signal="$(join_by_comma "${SIGNAL_SAMPLES[@]}")"
  actual_background="$(join_by_comma "${BACKGROUND_SAMPLES[@]}")"
  if [[ "$actual_signal" == "$expected_signal" && "$actual_background" == "$expected_background" ]]; then
    return 0
  fi
  if [[ "${RJ_AUAU_TIGHT_BDT_ALLOW_SAMPLE_COMPOSITION_CONTROL:-0}" == "1" ]]; then
    warn "Using non-default AuAu BDT embedded sample composition: signal=${actual_signal:-<none>} background=${actual_background:-<none>}"
    return 0
  fi
  die "AuAu BDT embedded training defaults must be Photon12+20 and Jet12+20+30+40. Got signal=${actual_signal:-<none>} background=${actual_background:-<none>}. Set RJ_AUAU_TIGHT_BDT_ALLOW_SAMPLE_COMPOSITION_CONTROL=1 only for deliberate sample-composition tests."
}

enforce_default_embedded_samples

guard_generated_path() {
  local label="$1"
  local path="$2"
  [[ -n "$path" ]] || die "${label} resolved to an empty path"
  case "$path" in
    /cvmfs/*|/opt/sphenix/*|/usr/*|/bin/*|/lib/*|/lib64/*|/etc/*)
      die "${label} resolved to a protected/read-only system path: ${path}. Repo root is ${RJ_REPO_BASE}; check for environment variable collisions before rerunning."
      ;;
  esac
}

log_path_plan() {
  local label="$1"
  shift
  say "${label} path plan:"
  say "  repo root : ${RJ_REPO_BASE}"
  say "  cwd       : $(pwd -P)"
  if [[ -n "${BASE:-}" && "${BASE:-}" != "${RJ_REPO_BASE}" ]]; then
    warn "environment variable BASE='${BASE}' differs from repo root; using RJ_REPO_BASE='${RJ_REPO_BASE}'"
  fi
  local item
  for item in "$@"; do
    say "  ${item}"
  done
}

setup_sphenix_stack_env() {
  export USER="${USER:-$(id -u -n)}"
  export LOGNAME="${LOGNAME:-$USER}"
  export HOME="/sphenix/u/${LOGNAME}"
  local myinstall="/sphenix/u/${USER}/thesisAnalysis/install"
  local myinstall_auau="/sphenix/u/${USER}/thesisAnalysis_auau/install"
  set +u
  # shellcheck disable=SC1091
  source /opt/sphenix/core/bin/sphenix_setup.sh -n
  if [[ -d "$myinstall" ]]; then
    # shellcheck disable=SC1091
    source /opt/sphenix/core/bin/setup_local.sh "$myinstall" || true
  fi
  if [[ -d "$myinstall_auau" ]]; then
    # shellcheck disable=SC1091
    source /opt/sphenix/core/bin/setup_local.sh "$myinstall_auau" || true
  fi
  set -u
  return 0
}

setup_ml_python_env() {
  setup_sphenix_stack_env
  if [[ -n "${BASE:-}" && "${BASE:-}" != "${RJ_REPO_BASE}" ]]; then
    warn "sPHENIX setup exported BASE='${BASE}'; sidecar repo root remains RJ_REPO_BASE='${RJ_REPO_BASE}'"
  fi
  local ml_python_prefix=""
  local ml_python_real=""
  local ml_python_real_prefix=""
  if [[ "$ML_PYTHON" == */* ]]; then
    ml_python_prefix="$(cd "$(dirname "$ML_PYTHON")/.." && pwd -P 2>/dev/null || true)"
    ml_python_real="$(readlink -f "$ML_PYTHON" 2>/dev/null || true)"
    if [[ -n "$ml_python_real" ]]; then
      ml_python_real_prefix="$(cd "$(dirname "$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
    fi
  fi
  local -a ld_candidates=()
  [[ -n "$ml_python_prefix" ]] && ld_candidates+=("${ml_python_prefix}/lib" "${ml_python_prefix}/lib64")
  [[ -n "$ml_python_real_prefix" && "$ml_python_real_prefix" != "$ml_python_prefix" ]] && ld_candidates+=("${ml_python_real_prefix}/lib" "${ml_python_real_prefix}/lib64")
  if ((${#ld_candidates[@]})); then
    local ld_joined=""
    local d
    for d in "${ld_candidates[@]}"; do
      [[ -d "$d" ]] || continue
      ld_joined="${ld_joined:+${ld_joined}:}${d}"
    done
    [[ -n "$ld_joined" ]] && export LD_LIBRARY_PATH="${ld_joined}:${LD_LIBRARY_PATH:-}"
  fi
  unset PYTHONHOME
}

send_summary_email() {
  local subject="$1"
  local summary="$2"
  if command -v mail >/dev/null 2>&1; then
    mail -s "$subject" "$NOTIFY_EMAILS" < "$summary" || true
  fi
}

usage() {
  cat <<'EOF'
Usage:
  ./scripts/auau_tight_bdt_pipeline.sh localTest [Nevents] [NFILES=N]
  ./scripts/auau_tight_bdt_pipeline.sh smokeTest [groupSize N]
  ./scripts/auau_tight_bdt_pipeline.sh condorExtract [groupSize N]
  ./scripts/auau_tight_bdt_pipeline.sh trainFromExtraction SOURCE=/path
  ./scripts/auau_tight_bdt_pipeline.sh trainCentInput3x3FromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_bdt_pipeline.sh trainWidthStudyPt1530FromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_bdt_pipeline.sh trainWidthStudyWindowsFromExtraction SOURCE=/path [MODEL_DIR=/path] [PT_WINDOWS=5:35,10:35,15:35]
  ./scripts/auau_tight_bdt_pipeline.sh trainEtFineCentStudyFromExtraction SOURCE=/path [MODEL_DIR=/path]
  ./scripts/auau_tight_bdt_pipeline.sh trainExpandedFromExtraction SOURCE=/path [PLAN_ONLY=1]
  ./scripts/auau_tight_bdt_pipeline.sh trainExpandedFromExtractionCondor SOURCE=/path [groupSize N]
  ./scripts/auau_tight_bdt_pipeline.sh finalizeExpandedTraining SOURCE=/path MODEL_DIR=/path SUB_ROOT=/path
  ./scripts/auau_tight_bdt_pipeline.sh applyCheck MODEL_DIR=/path
  ./scripts/auau_tight_bdt_pipeline.sh validateOnSim SOURCE=/path MODEL_DIR=/path
  ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor SOURCE=/path MODEL_DIR=/path [groupSize N]
  ./scripts/auau_tight_bdt_pipeline.sh deriveWorkingPointsFromValidation VALIDATION=/path [TARGET=0.80]
  ./scripts/auau_tight_bdt_pipeline.sh generateWorkingPointConfig TEMPLATE=/path WORKING_POINTS=/path PRODUCT_MAP=a=b,c=d OUT=/path [MODEL_DIR=/path]

Sidecar AuAu tight-BDT workflow:
  extraction reads embeddedPhoton12/20 and embeddedJet12/20/30/40 samples by default,
  with optional explicit sample overrides for deliberate sample-composition controls,
  writes AuAuPhotonIDTrainingTree ROOT files, and avoids normal cfg-tag
  histogram production. validateOnSim scores those same trees with the
  final TMVA ROOT models and writes quick ROC/AUC simulation diagnostics.
  deriveWorkingPointsFromValidation converts a completed validation report into
  target-signal-efficiency BDT threshold manifests for downstream MC campaigns.
  generateWorkingPointConfig injects product-mapped working points into a frozen
  RecoilJets YAML without touching the template.
  validateOnSimCondor shards the ROOT scoring over Condor, then merges the
  score caches into the same final validation report/plots.
  trainExpandedFromExtraction builds the registry-driven expanded model study
  matrix for model selection; it does not change the normal three-model path.
  trainCentInput3x3FromExtraction trains one focused centrality-input model
  that swaps the standard shower widths for 3x3 tower-window widths.
  trainWidthStudyPt1530FromExtraction trains the focused 15-30 GeV
  centrality-input comparison: base widths, 3x3 widths, and base+3x3 widths.
  trainWidthStudyWindowsFromExtraction trains the same centrality-input
  width comparison for explicit pT windows, default 5-35, 10-35, and 15-35 GeV.
  trainEtFineCentStudyFromExtraction trains the 15-35 GeV fine-E_T
  centrality study: one centrality-input model, fine-E_T centrality-input
  models, and fine-E_T x centrality-binned models using base+3x3 widths.
  finalizeExpandedTraining re-merges an already completed expanded Condor run
  without retraining, useful after bookkeeping-only fixes.
EOF
}

wrapper_path() {
  if [[ -f "${RJ_REPO_BASE}/RecoilJets_Condor_AuAu.sh" ]]; then
    printf '%s\n' "${RJ_REPO_BASE}/RecoilJets_Condor_AuAu.sh"
  elif [[ -f "${RJ_REPO_BASE}/scripts/RecoilJets_Condor_AuAu.sh" ]]; then
    printf '%s\n' "${RJ_REPO_BASE}/scripts/RecoilJets_Condor_AuAu.sh"
  else
    die "Cannot find RecoilJets_Condor_AuAu.sh under ${RJ_REPO_BASE}"
  fi
}

sample_dataset() {
  case "$1" in
    run28_embeddedPhoton12|run28_embeddedPhoton20) printf '%s\n' "isSimEmbedded" ;;
    run28_embeddedJet12|run28_embeddedJet20|run28_embeddedJet30|run28_embeddedJet40) printf '%s\n' "isSimEmbeddedInclusive" ;;
    *) die "Unknown AuAu BDT training sample: $1" ;;
  esac
}

sample_class() {
  case "$1" in
    run28_embeddedPhoton12|run28_embeddedPhoton20) printf '%s\n' "signal" ;;
    run28_embeddedJet12|run28_embeddedJet20|run28_embeddedJet30|run28_embeddedJet40) printf '%s\n' "background" ;;
    *) die "Unknown AuAu BDT training sample: $1" ;;
  esac
}

sample_dir() {
  local sample="$1"
  local short="${sample#run28_}"
  local dir="${SIM_ROOT}/${sample}"
  [[ -d "$dir" ]] || dir="${SIM_ROOT}/${short}"
  [[ -d "$dir" ]] || die "Missing ${sample} lists under ${SIM_ROOT}. Run ./scripts/makeThesisSimLists.sh isSimEmbedded and isSimEmbeddedInclusive first."
  printf '%s\n' "$dir"
}

build_sample_master() {
  local sample="$1"
  local out="$2"
  local dir; dir="$(sample_dir "$sample")"
  local calo="${dir}/DST_CALO_CLUSTER.matched.list"
  local g4="${dir}/G4Hits.matched.list"
  local jets="${dir}/DST_JETS.matched.list"
  local glob="${dir}/DST_GLOBAL.matched.list"
  local mbd="${dir}/DST_MBD_EPD.matched.list"
  [[ -s "$calo" ]] || die "Missing $calo"
  [[ -s "$g4" ]] || die "Missing $g4"
  [[ -s "$jets" ]] || die "Missing $jets"
  [[ -s "$glob" ]] || die "Missing $glob"
  [[ -s "$mbd" ]] || die "Missing $mbd"
  paste "$calo" "$g4" "$jets" "$glob" "$mbd" | grep -E -v '^[[:space:]]*($|#)' > "$out"
  [[ -s "$out" ]] || die "Built empty master list for $sample"
}

make_training_yaml() {
  local out="$1"
  [[ -s "$MASTER_YAML" ]] || die "Missing master YAML: $MASTER_YAML"
  mkdir -p "$(dirname "$out")"
  cp -f "$MASTER_YAML" "$out"
  {
    echo
    echo "# auau_tight_bdt_pipeline.sh extraction overrides"
    echo "# Canonical embedded training requires MinimumBiasClassifier pass."
    echo "setMinBiasClassifer: true"
    echo "preselection: reference"
    echo "tight: reference"
    echo "nonTight: reference"
    echo "clusterUEpipeline: baseVariant"
    echo "vz_cut_cm: 10"
    echo "jet_pt_min: 5"
    echo "back_to_back_phi_cut: 0.875"
    echo "coneR: 0.4"
    echo "isSlidingIso: false"
    echo "fixedGeV: 4.0"
    echo "auau_bdt_training_tree: true"
    echo "auau_bdt_training_tree_max_entries: ${RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES:-0}"
    echo "auau_bdt_npb_data_tagging: false"
  } >> "$out"
}

make_root_manifest() {
  local root="$1"
  local out="$2"
  mkdir -p "$(dirname "$out")"
  if [[ -n "${RJ_AUAU_BDT_ROOT_MANIFEST:-}" ]]; then
    [[ -s "$RJ_AUAU_BDT_ROOT_MANIFEST" ]] || die "RJ_AUAU_BDT_ROOT_MANIFEST is set but missing/empty: $RJ_AUAU_BDT_ROOT_MANIFEST"
    cp -f "$RJ_AUAU_BDT_ROOT_MANIFEST" "$out"
    [[ -s "$out" ]] || die "Copied empty ROOT manifest from $RJ_AUAU_BDT_ROOT_MANIFEST"
    return 0
  fi
  if [[ "${RJ_AUAU_BDT_REUSE_EXISTING_MANIFEST:-0}" == "1" && -s "$out" ]]; then
    return 0
  fi
  find -L "$root" -type f -name '*.root' | sort -V > "$out" || true
  [[ -s "$out" ]] || die "No ROOT files found under $root"
}

validate_training_tree() {
  local manifest="$1"
  local report="$2"
  if [[ "${RJ_AUAU_BDT_SKIP_TRAINING_TREE_VALIDATE:-0}" == "1" ]]; then
    if [[ "$require_embedded_minbias" == "1" ]]; then
      die "Canonical MinimumBiasClassifier training may not skip training-tree validation. Use a fresh validated extraction, or label an unrestricted control with RJ_AUAU_TIGHT_BDT_REQUIRE_EMBEDDED_MINBIAS=0 and RJ_AUAU_TIGHT_BDT_ALLOW_UNRESTRICTED_EMBEDDED_CONTROL=1."
    fi
    mkdir -p "$(dirname "$report")"
    setup_ml_python_env
    "$ML_PYTHON" - "$manifest" "$report" <<'PY'
import json
import sys
from pathlib import Path
manifest = Path(sys.argv[1])
report = Path(sys.argv[2])
paths = [line.strip() for line in manifest.read_text().splitlines() if line.strip()]
report.write_text(json.dumps({
    "status": "SKIPPED_BY_OPERATOR",
    "reason": "RJ_AUAU_BDT_SKIP_TRAINING_TREE_VALIDATE=1; manifest came from a prior READY registry.",
    "file_count": len(paths),
    "manifest": str(manifest),
}, indent=2, sort_keys=True) + "\n")
print(f"[OK] training tree validation skipped by env; manifest files={len(paths)}")
PY
    return 0
  fi
  setup_ml_python_env
  "$ML_PYTHON" - "$manifest" "$report" <<'PY'
import json
import os
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
report = Path(sys.argv[2])
paths = [Path(x.strip()) for x in manifest.read_text().splitlines() if x.strip()]
rows = []
total = 0
signal = 0
background = 0
minbias_pass = 0
minbias_nonpass = 0
missing_minbias_branch = 0
try:
    import uproot
except Exception as exc:
    report.write_text(json.dumps({"status": "CHECK", "reason": f"uproot import failed: {exc}", "files": [str(p) for p in paths]}, indent=2) + "\n")
    print(f"[WARN] uproot validation skipped: {exc}")
    sys.exit(0)
progress_every = int(os.environ.get("RJ_AUAU_TIGHT_BDT_VALIDATE_PROGRESS_EVERY", "500") or "0")
for idx, path in enumerate(paths, 1):
    if progress_every > 0 and (idx == 1 or idx % progress_every == 0 or idx == len(paths)):
        print(f"[auauTightBDT] validating training tree {idx}/{len(paths)}: {path}", flush=True)
    with uproot.open(path) as f:
        try:
            tree = f["AuAuPhotonIDTrainingTree"]
        except Exception:
            rows.append({"file": str(path), "entries": 0, "missing_tree": True})
            continue
        n = int(tree.num_entries)
        total += n
        try:
            labels = tree["is_signal"].array(library="np")
        except Exception:
            rows.append({"file": str(path), "entries": n, "missing_is_signal": True})
            continue
        n_signal = int((labels == 1).sum())
        n_background = int((labels == 0).sum())
        signal += n_signal
        background += n_background
        require_minbias = os.environ.get("RJ_AUAU_TIGHT_BDT_REQUIRE_EMBEDDED_MINBIAS", "1") == "1"
        if require_minbias:
            if "minimum_bias_classifier_decision" not in tree.keys():
                missing_minbias_branch += 1
                decisions = None
            else:
                decisions = tree["minimum_bias_classifier_decision"].array(library="np")
                minbias_pass += int((decisions == 2).sum())
                minbias_nonpass += int((decisions != 2).sum())
        else:
            decisions = None
        sample_class = "signal" if "/signal/" in str(path) else ("background" if "/background/" in str(path) else "unknown")
        rows.append({
            "file": str(path),
            "entries": n,
            "class": sample_class,
            "signal_entries": n_signal,
            "background_entries": n_background,
            "minimum_bias_pass_entries": int((decisions == 2).sum()) if decisions is not None else None,
            "minimum_bias_nonpass_entries": int((decisions != 2).sum()) if decisions is not None else None,
            "missing_minimum_bias_branch": require_minbias and decisions is None
        })
require_minbias = os.environ.get("RJ_AUAU_TIGHT_BDT_REQUIRE_EMBEDDED_MINBIAS", "1") == "1"
base_ok = total > 0 and signal > 0 and background > 0 and all(not r.get("missing_tree") and not r.get("missing_is_signal") for r in rows)
minbias_ok = (not require_minbias) or (missing_minbias_branch == 0 and minbias_nonpass == 0 and minbias_pass == total)
status = "PASS" if base_ok and minbias_ok else "FAIL"
report.write_text(json.dumps({
    "status": status,
    "require_embedded_minimum_bias_classifier": require_minbias,
    "total_entries": total,
    "signal_entries": signal,
    "background_entries": background,
    "minimum_bias_pass_entries": minbias_pass,
    "minimum_bias_nonpass_entries": minbias_nonpass,
    "missing_minimum_bias_branch_files": missing_minbias_branch,
    "files": rows
}, indent=2, sort_keys=True) + "\n")
print(f"[OK] training tree validation: status={status} entries={total} signal={signal} background={background} minbias_pass={minbias_pass} minbias_nonpass={minbias_nonpass} missing_minbias_branch={missing_minbias_branch} files={len(paths)}")
sys.exit(0 if status == "PASS" else 3)
PY
}

write_config_snippet() {
  local model_dir="$1"
  cat > "${model_dir}/analysis_config_snippet.yaml" <<EOF
# Paste these paths into analysis_config.yaml after validation.
auau_tight_bdt_centINDcontrol_model_file: ${model_dir}/auau_tight_bdt_centINDcontrol_allCent_tmva.root
auau_tight_bdt_centAsFeat_model_file: ${model_dir}/auau_tight_bdt_centAsFeat_allCent_tmva.root
auau_tight_bdt_centDep_model_files: [${model_dir}/auau_tight_bdt_centDepBDTs_cent_000_020_tmva.root, ${model_dir}/auau_tight_bdt_centDepBDTs_cent_020_050_tmva.root, ${model_dir}/auau_tight_bdt_centDepBDTs_cent_050_080_tmva.root]
EOF
}

train_from_manifest() {
  local manifest="$1"
  local stamp="${RJ_AUAU_TIGHT_BDT_TRAIN_STAMP:-$(ts)}"
  local model_dir="${RJ_AUAU_TIGHT_BDT_MODEL_DIR:-${MODEL_BASE}/tight_${stamp}}"
  AUAU_TIGHT_BDT_LAST_MODEL_DIR="$model_dir"
  mkdir -p "$model_dir"
  setup_ml_python_env
  say "Training three tight-BDT products from @${manifest}"
  "$ML_PYTHON" "$TRAIN_SCRIPT" --task tight --tight-mode centINDcontrol \
    --input "@${manifest}" --outdir "$model_dir" \
    --prefix auau_tight_bdt_centINDcontrol
  "$ML_PYTHON" "$TRAIN_SCRIPT" --task tight --tight-mode centAsFeat \
    --input "@${manifest}" --outdir "$model_dir" \
    --prefix auau_tight_bdt_centAsFeat
  "$ML_PYTHON" "$TRAIN_SCRIPT" --task tight --tight-mode centDepBDTs \
    --input "@${manifest}" --outdir "$model_dir" \
    --prefix auau_tight_bdt_centDepBDTs \
    --cent-bins 0:20,20:50,50:80
  write_config_snippet "$model_dir"
  say "Model directory: $model_dir"
}

run_local_test() {
  local nevents="${RJ_AUAU_TIGHT_BDT_LOCAL_NEVENTS:-1000}"
  local nfiles="${RJ_AUAU_TIGHT_BDT_LOCAL_NFILES:-1}"
  if [[ "${1:-}" =~ ^[0-9]+$|^-1$ ]]; then
    nevents="$1"
    shift
  fi
  for tok in "$@"; do
    case "$tok" in
      NFILES=*) nfiles="${tok#NFILES=}" ;;
      VERBOSE=*) export RJ_VERBOSITY="${tok#VERBOSE=}" ;;
    esac
  done
  local stamp="${RJ_AUAU_TIGHT_BDT_STAMP:-$(ts)}"
  local run_root="${LOCAL_BASE}/tight_${stamp}"
  local manifest_dir="${run_root}/manifests"
  local extraction_root="${run_root}/extraction"
  local report_dir="${run_root}/reports"
  local yaml="${manifest_dir}/analysis_config_auau_tight_bdt_training.yaml"
  guard_generated_path "localTest run root" "$run_root"
  log_path_plan "localTest" \
    "run root  : ${run_root}" \
    "manifest  : ${manifest_dir}" \
    "extraction: ${extraction_root}" \
    "report    : ${report_dir}"
  mkdir -p "$manifest_dir" "$extraction_root" "$report_dir"
  make_training_yaml "$yaml"
  local wrapper; wrapper="$(wrapper_path)"
  say "localTest output root: $run_root"
  for sample in "${ALL_SAMPLES[@]}"; do
    local master="${manifest_dir}/${sample}_5col.list"
    local chunk="${manifest_dir}/${sample}_local_${nfiles}.list"
    build_sample_master "$sample" "$master"
    head -n "$nfiles" "$master" > "$chunk"
    local dataset; dataset="$(sample_dataset "$sample")"
    local klass; klass="$(sample_class "$sample")"
    local dest="${extraction_root}/${klass}/${sample}"
    mkdir -p "$dest"
    say "extract local sample=${sample} dataset=${dataset} nfiles=${nfiles} nevents=${nevents}"
    RJ_CONFIG_YAML="$yaml" \
    RJ_MACRO_PATH="$TRAIN_MACRO" \
    RJ_AUAU_BDT_EXTRACT_ONLY=1 \
    RJ_AUAU_BDT_TRAINING_TREE=1 \
    RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER="$require_embedded_minbias" \
    RJ_DISABLE_ID_FANOUT=1 \
    RJ_DISABLE_ISO_CONE_INTERNALIZATION=1 \
    RJ_DISABLE_JET_PT_INTERNALIZATION=1 \
    RJ_DISABLE_DPHI_INTERNALIZATION=1 \
    bash "$wrapper" "$sample" "$chunk" "$dataset" LOCAL "$nevents" 0 NONE "$dest"
  done
  local root_manifest="${manifest_dir}/training_roots.list"
  make_root_manifest "$extraction_root" "$root_manifest"
  validate_training_tree "$root_manifest" "${report_dir}/training_tree_validation.json"
  train_from_manifest "$root_manifest"
  say "localTest complete: $run_root"
}

run_condor_extract() {
  local mode="$1"
  shift || true
  local group_size="${RJ_AUAU_TIGHT_BDT_GROUP_SIZE:-3}"
  local max_jobs_per_sample=0
  if [[ "$mode" == "smokeTest" ]]; then
    max_jobs_per_sample="${RJ_AUAU_TIGHT_BDT_SMOKE_MAX_JOBS_PER_SAMPLE:-2}"
  fi
  while (($#)); do
    case "$1" in
      groupSize) group_size="${2:?missing value after groupSize}"; shift 2 ;;
      maxJobs) max_jobs_per_sample="${2:?missing value after maxJobs}"; shift 2 ;;
      *) shift ;;
    esac
  done
  local stamp="${RJ_AUAU_TIGHT_BDT_STAMP:-$(ts)}"
  local run_root="${TRAIN_BASE}/auauTightBDT_${stamp}"
  local sub_root="${RJ_REPO_BASE}/condor_sub/auauTightBDT_${stamp}"
  local extraction_root="${run_root}/extraction"
  local manifest_dir="${run_root}/manifests"
  local report_dir="${run_root}/reports"
  local yaml="${manifest_dir}/analysis_config_auau_tight_bdt_training.yaml"
  guard_generated_path "condorExtract run root" "$run_root"
  guard_generated_path "condorExtract submit root" "$sub_root"
  log_path_plan "${mode}" \
    "run root  : ${run_root}" \
    "submit    : ${sub_root}" \
    "manifest  : ${manifest_dir}" \
    "extraction: ${extraction_root}" \
    "report    : ${report_dir}" \
    "groupSize : ${group_size}" \
    "maxJobs   : ${max_jobs_per_sample}"
  mkdir -p "$sub_root" "$manifest_dir" "$report_dir" "$extraction_root"
  make_training_yaml "$yaml"
  local args_file="${sub_root}/extract_args.txt"
  : > "$args_file"
  local wrapper; wrapper="$(wrapper_path)"
  local chunk_idx=0
  local nevents_for_jobs="${RJ_AUAU_TIGHT_BDT_NEVENTS:-0}"
  if [[ "$mode" == "smokeTest" && -z "${RJ_AUAU_TIGHT_BDT_NEVENTS:-}" ]]; then
    nevents_for_jobs=1000
  fi
  for sample in "${ALL_SAMPLES[@]}"; do
    local master="${manifest_dir}/${sample}_5col.list"
    build_sample_master "$sample" "$master"
    local split_prefix="${sub_root}/${sample}_grp_"
    find "$sub_root" -maxdepth 1 -type f -name "$(basename "$split_prefix")*" -delete
    if (( max_jobs_per_sample > 0 )); then
      head -n "$((group_size * max_jobs_per_sample))" "$master" |
        split -l "$group_size" -d -a 5 - "$split_prefix"
    else
      split -l "$group_size" -d -a 5 "$master" "$split_prefix"
    fi
    local queued=0
    for raw in "${split_prefix}"*; do
      [[ -s "$raw" ]] || { rm -f "$raw"; continue; }
      queued=$((queued + 1))
      chunk_idx=$((chunk_idx + 1))
      local chunk="${raw}.list"
      mv "$raw" "$chunk"
      local dataset; dataset="$(sample_dataset "$sample")"
      local klass; klass="$(sample_class "$sample")"
      local dest="${extraction_root}/${klass}/${sample}"
      printf '%s %s %s %s %s %s\n' "$sample" "$chunk" "$dataset" "$nevents_for_jobs" "$chunk_idx" "$dest" >> "$args_file"
    done
    say "planned sample=${sample} chunks=${queued} cap=${max_jobs_per_sample:-0}"
  done
  [[ -s "$args_file" ]] || die "No extraction jobs planned"

  local sub="${sub_root}/auau_tight_bdt_extract.sub"
  local reqmem="${RJ_AUAU_TIGHT_BDT_REQUEST_MEMORY:-3000MB}"
  local allow_nonzero_with_output=0
  if [[ "${RJ_AUAU_TIGHT_BDT_TOLERATE_ROOT_ABORT_WITH_OUTPUT:-0}" == "1" ]]; then
    allow_nonzero_with_output=1
  fi
  local force_status="${RJ_FORCE_CALO_TOWER_STATUS_FOR_EMBEDDED:-0}"
  local event_calo_require_isgood="${RJ_EVENT_CALO_REQUIRE_ISGOOD:-1}"
  local status_audit="${RJ_CALO_STATUS_AUDIT:-0}"
  local status_input_prefix="${RJ_CALO_TOWER_STATUS_INPUT_PREFIX:-}"
  cat > "$sub" <<EOF
universe = vanilla
executable = /usr/bin/bash
arguments = ${wrapper} \$(sample) \$(chunk) \$(dataset) \$(Cluster) \$(nevents) \$(chunkidx) NONE \$(dest)
output = ${sub_root}/extract_\$(Cluster)_\$(Process).out
error = ${sub_root}/extract_\$(Cluster)_\$(Process).err
log = ${sub_root}/extract_\$(Cluster).log
request_memory = ${reqmem}
notification = Never
environment = "RJ_CONFIG_YAML=${yaml} RJ_MACRO_PATH=${TRAIN_MACRO} RJ_AUAU_BDT_EXTRACT_ONLY=1 RJ_AUAU_BDT_TRAINING_TREE=1 RJ_REQUIRE_EMBEDDED_MINBIAS_CLASSIFIER=${require_embedded_minbias} RJ_AUAU_TIGHT_BDT_REQUIRE_EMBEDDED_MINBIAS=${require_embedded_minbias} RJ_DISABLE_ID_FANOUT=1 RJ_DISABLE_ISO_CONE_INTERNALIZATION=1 RJ_DISABLE_JET_PT_INTERNALIZATION=1 RJ_DISABLE_DPHI_INTERNALIZATION=1 RJ_PROFILE_JOB=1 RJ_ALLOW_NONZERO_WITH_ROOT_OUTPUT=${allow_nonzero_with_output} RJ_FORCE_CALO_TOWER_STATUS_FOR_EMBEDDED=${force_status} RJ_EVENT_CALO_REQUIRE_ISGOOD=${event_calo_require_isgood} RJ_CALO_STATUS_AUDIT=${status_audit} RJ_CALO_TOWER_STATUS_INPUT_PREFIX=${status_input_prefix}"
queue sample,chunk,dataset,nevents,chunkidx,dest from ${args_file}
EOF

  local notify="${sub_root}/notify.sh"
  cat > "$notify" <<EOF
#!/usr/bin/env bash
set -euo pipefail
echo "[auauTightBDT] notify env setup start" >&2
export USER="\${USER:-\$(id -u -n)}"
export LOGNAME="\${LOGNAME:-\$USER}"
export HOME="/sphenix/u/\${LOGNAME}"
MYINSTALL="/sphenix/u/\${USER}/thesisAnalysis/install"
MYINSTALL_AUAU="/sphenix/u/\${USER}/thesisAnalysis_auau/install"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
echo "[auauTightBDT] notify sphenix_setup rc=\$?" >&2
if [[ -d "\$MYINSTALL" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL" || true
  echo "[auauTightBDT] notify setup_local rc=\$?" >&2
fi
if [[ -d "\$MYINSTALL_AUAU" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL_AUAU" || true
  echo "[auauTightBDT] notify setup_local_auau rc=\$?" >&2
fi
set -u
echo "[auauTightBDT] notify env setup done" >&2
ml_python="${ML_PYTHON}"
ml_python_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="\$(readlink -f "\$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
if [[ -n "\$ml_python_real" ]]; then
  ml_python_real_prefix="\$(cd "\$(dirname "\$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
fi
ld_joined=""
for d in "\$ml_python_prefix/lib" "\$ml_python_prefix/lib64" "\$ml_python_real_prefix/lib" "\$ml_python_real_prefix/lib64"; do
  [[ -n "\$d" && -d "\$d" ]] || continue
  case ":\$ld_joined:" in *":\$d:"*) ;; *) ld_joined="\${ld_joined:+\$ld_joined:}\$d" ;; esac
done
[[ -n "\$ld_joined" ]] && export LD_LIBRARY_PATH="\$ld_joined:\${LD_LIBRARY_PATH:-}"
unset PYTHONHOME
root_manifest="${manifest_dir}/training_roots.list"
expected_roots="$(wc -l < "${args_file}" | tr -d ' ')"
find -L "${extraction_root}" -type f -name '*.root' | sort -V > "\$root_manifest" || true
nroots=\$(wc -l < "\$root_manifest" | tr -d ' ')
tree_entries=0
validation_note=""
if [[ "\$nroots" != "0" ]]; then
  validation_note=\$("\$ml_python" - "\$root_manifest" 2>&1 <<'PY' || true
import sys
from pathlib import Path
try:
    import uproot
except Exception as exc:
    print(f"UPROOT_IMPORT_FAILED {exc}")
    raise SystemExit(0)
manifest = Path(sys.argv[1])
total = 0
signal = 0
background = 0
missing = 0
minbias_pass = 0
minbias_nonpass = 0
missing_minbias_branch = 0
for raw in manifest.read_text().splitlines():
    path = raw.strip()
    if not path:
        continue
    try:
        with uproot.open(path) as f:
            try:
                tree = f["AuAuPhotonIDTrainingTree"]
            except Exception:
                missing += 1
                continue
            n = int(tree.num_entries)
            total += n
            try:
                labels = tree["is_signal"].array(library="np")
            except Exception:
                missing += 1
                continue
            signal += int((labels == 1).sum())
            background += int((labels == 0).sum())
            if "minimum_bias_classifier_decision" not in tree.keys():
                missing_minbias_branch += 1
            else:
                decisions = tree["minimum_bias_classifier_decision"].array(library="np")
                minbias_pass += int((decisions == 2).sum())
                minbias_nonpass += int((decisions != 2).sum())
    except Exception as exc:
        print(f"ROOT_OPEN_FAILED {path} {exc}")
print(f"TREE_ENTRIES {total} SIGNAL_ENTRIES {signal} BACKGROUND_ENTRIES {background} MISSING_TREE_FILES {missing} MINBIAS_PASS_ENTRIES {minbias_pass} MINBIAS_NONPASS_ENTRIES {minbias_nonpass} MISSING_MINBIAS_BRANCH_FILES {missing_minbias_branch}")
PY
)
  tree_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {print \$2; exit}')
  tree_entries="\${tree_entries:-0}"
  signal_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {for (i=1; i<=NF; ++i) if (\$i=="SIGNAL_ENTRIES") {print \$(i+1); exit}}')
  background_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {for (i=1; i<=NF; ++i) if (\$i=="BACKGROUND_ENTRIES") {print \$(i+1); exit}}')
  signal_entries="\${signal_entries:-0}"
  background_entries="\${background_entries:-0}"
  minbias_pass_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {for (i=1; i<=NF; ++i) if (\$i=="MINBIAS_PASS_ENTRIES") {print \$(i+1); exit}}')
  minbias_nonpass_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {for (i=1; i<=NF; ++i) if (\$i=="MINBIAS_NONPASS_ENTRIES") {print \$(i+1); exit}}')
  missing_minbias_branch_files=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {for (i=1; i<=NF; ++i) if (\$i=="MISSING_MINBIAS_BRANCH_FILES") {print \$(i+1); exit}}')
  minbias_pass_entries="\${minbias_pass_entries:-0}"
  minbias_nonpass_entries="\${minbias_nonpass_entries:-0}"
  missing_minbias_branch_files="\${missing_minbias_branch_files:-0}"
fi
status=READY
expected_signal_samples=${#SIGNAL_SAMPLES[@]}
expected_background_samples=${#BACKGROUND_SAMPLES[@]}
if [[ "\$nroots" == "0" || "\$tree_entries" == "0" || "\$nroots" != "\$expected_roots" ]]; then status=CHECK; fi
if [[ "\$expected_signal_samples" != "0" && "\${signal_entries:-0}" == "0" ]]; then status=CHECK; fi
if [[ "\$expected_background_samples" != "0" && "\${background_entries:-0}" == "0" ]]; then status=CHECK; fi
if [[ "${require_embedded_minbias}" == "1" && ("\${minbias_pass_entries:-0}" != "\$tree_entries" || "\${minbias_nonpass_entries:-0}" != "0" || "\${missing_minbias_branch_files:-0}" != "0") ]]; then status=CHECK; fi
summary="${report_dir}/final_summary.txt"
{
  echo "RECOILJETS_STAGE_EMAIL_V1"
  echo "dataset=isSimEmbeddedAndInclusive"
  echo "stage=auauTightBDT_${mode}"
  echo "status=\${status}"
  echo "training_root=${run_root}"
  echo "root_manifest=\${root_manifest}"
  echo "root_count=\${nroots}"
  echo "expected_root_count=\${expected_roots}"
  echo "tree_entries=\${tree_entries}"
  echo "signal_entries=\${signal_entries:-0}"
  echo "background_entries=\${background_entries:-0}"
  echo "require_embedded_minimum_bias_classifier=${require_embedded_minbias}"
  echo "minimum_bias_pass_entries=\${minbias_pass_entries:-0}"
  echo "minimum_bias_nonpass_entries=\${minbias_nonpass_entries:-0}"
  echo "missing_minimum_bias_branch_files=\${missing_minbias_branch_files:-0}"
  echo "expected_signal_samples=\${expected_signal_samples}"
  echo "expected_background_samples=\${expected_background_samples}"
  if [[ -n "\${validation_note:-}" ]]; then
    echo "validation_note=\${validation_note}"
  fi
  echo "next_action=RJ_ML_PYTHON=${ML_PYTHON} ./scripts/auau_tight_bdt_pipeline.sh trainFromExtraction SOURCE=${run_root}"
} > "\$summary"
if command -v mail >/dev/null 2>&1; then
  mail -s "[RecoilJets][auauTightBDT_${mode}][\${status}]" "${NOTIFY_EMAILS}" < "\$summary" || true
fi
cat "\$summary"
EOF
  chmod +x "$notify"
  local notify_sub="${sub_root}/notify.sub"
  cat > "$notify_sub" <<EOF
universe = scheduler
executable = ${notify}
output = ${sub_root}/notify.out
error = ${sub_root}/notify.err
log = ${sub_root}/notify.log
notification = Never
queue
EOF
  local dag="${sub_root}/auau_tight_bdt_${mode}.dag"
  cat > "$dag" <<EOF
JOB EXTRACT ${sub}
FINAL NOTIFY ${notify_sub}
EOF
  say "DAG: $dag"
  say "run root: $run_root"
  say "jobs: $(wc -l < "$args_file" | tr -d ' ')  groupSize=${group_size}  nevents=${nevents_for_jobs}  request_memory=${reqmem}"
  if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
    echo "RECOILJETS_AUAU_TIGHT_BDT_DRYRUN_V1"
    echo "mode=${mode}"
    echo "run_root=${run_root}"
    echo "dag=${dag}"
    echo "jobs=$(wc -l < "$args_file" | tr -d ' ')"
    return 0
  fi
  condor_submit_dag "$dag"
}

train_from_extraction() {
  local source=""
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
    esac
  done
  [[ -n "$source" ]] || die "trainFromExtraction requires SOURCE=/path"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"
  local manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  local report_dir="${source}/reports"
  local validation_report="${report_dir}/training_tree_validation.json"
  local training_summary="${report_dir}/training_summary.txt"
  mkdir -p "$report_dir"
  make_root_manifest "$search_root" "$manifest"
  local validation_rc=0
  if validate_training_tree "$manifest" "$validation_report"; then
    validation_rc=0
  else
    validation_rc=$?
  fi

  local validation_status="UNKNOWN"
  local total_entries="0"
  local signal_entries="0"
  local background_entries="0"
  setup_ml_python_env
  eval "$("$ML_PYTHON" - "$validation_report" <<'PY'
import json
import shlex
import sys
from pathlib import Path

path = Path(sys.argv[1])
data = {}
if path.is_file():
    data = json.loads(path.read_text())
for key, default in (
    ("validation_status", data.get("status", "UNKNOWN")),
    ("total_entries", data.get("total_entries", 0)),
    ("signal_entries", data.get("signal_entries", 0)),
    ("background_entries", data.get("background_entries", 0)),
):
    print(f"{key}={shlex.quote(str(default))}")
PY
)"

  local train_rc=0
  local apply_status="SKIPPED"
  local model_dir=""
  if (( validation_rc == 0 )); then
    if train_from_manifest "$manifest"; then
      train_rc=0
    else
      train_rc=$?
    fi
    model_dir="${AUAU_TIGHT_BDT_LAST_MODEL_DIR:-}"
    if (( train_rc == 0 )) && [[ -n "$model_dir" ]]; then
      if apply_check "MODEL_DIR=${model_dir}"; then
        apply_status="PASS"
      else
        apply_status="FAIL"
      fi
    fi
  else
    train_rc=3
  fi

  local status="CHECK"
  if [[ "$validation_status" == "PASS" && "$apply_status" == "PASS" && "$train_rc" == "0" ]]; then
    status="READY"
  fi
  {
    echo "RECOILJETS_AUAU_TIGHT_BDT_TRAINING_V1"
    echo "status=${status}"
    echo "source=${source}"
    echo "root_manifest=${manifest}"
    echo "validation_report=${validation_report}"
    echo "model_dir=${model_dir:-unset}"
    echo "total_entries=${total_entries}"
    echo "signal_entries=${signal_entries}"
    echo "background_entries=${background_entries}"
    echo "validation_status=${validation_status}"
    echo "train_exit_code=${train_rc}"
    echo "applyCheck=${apply_status}"
    if [[ -n "$model_dir" ]]; then
      echo "analysis_config_snippet=${model_dir}/analysis_config_snippet.yaml"
      echo "next_action=Inspect ${model_dir}/analysis_config_snippet.yaml, add the final model paths to macros/analysis_config.yaml, then run the constrained data+MC BDT-variant validation."
    fi
  } > "$training_summary"
  cat "$training_summary"
  send_summary_email "[RecoilJets][auauTightBDT_trainFromExtraction][${status}]" "$training_summary"
  [[ "$status" == "READY" ]] || return 3
}

apply_check_registry() {
  local model_dir="$1"
  local registry="${2:-${model_dir}/model_registry.json}"
  [[ -s "$registry" ]] || die "Missing expanded model registry: $registry"
  setup_ml_python_env
  "$ML_PYTHON" - "$registry" <<'PY'
import json
import sys
from pathlib import Path
try:
    import ROOT
except Exception as exc:
    raise SystemExit(f"PyROOT import failed: {exc}")
registry = Path(sys.argv[1])
data = json.loads(registry.read_text())
bad = []
checked = 0
for spec in data.get("models", []):
    report = spec.get("report")
    if not report:
        bad.append(f"{spec.get('model_id')}: missing report")
        continue
    if report.get("status") == "skipped":
        continue
    path = Path(spec.get("output_tmva", ""))
    if not path.is_file() or path.stat().st_size <= 0:
        bad.append(f"{spec.get('model_id')}: missing {path}")
        continue
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        bad.append(f"{spec.get('model_id')}: unreadable {path}")
    else:
        checked += 1
    if f:
        f.Close()
if bad:
    raise SystemExit("Expanded applyCheck failed:\n  " + "\n  ".join(bad[:80]))
print(f"[OK] expanded applyCheck opened {checked} TMVA ROOT files from {registry}")
PY
}

validation_model_dir_from_yaml() {
  local yaml="${RJ_AUAU_BDT_VALIDATION_CONFIG:-${RJ_REPO_BASE}/macros/analysis_config_auau_bdt_validation.yaml}"
  [[ -f "$yaml" ]] || return 0
  awk -F: '
    /^[[:space:]]*auau_tight_bdt_expanded_model_dir[[:space:]]*:/ {
      sub(/^[[:space:]]+/, "", $2)
      sub(/[[:space:]]+$/, "", $2)
      print $2
      exit
    }
  ' "$yaml"
}

train_cent_input_3x3_from_extraction() {
  local source=""
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*) export RJ_AUAU_TIGHT_BDT_MODEL_DIR="${tok#MODEL_DIR=}" ;;
    esac
  done
  [[ -n "$source" ]] || die "trainCentInput3x3FromExtraction requires SOURCE=/path"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"

  local manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  local report_dir="${source}/reports"
  guard_generated_path "3x3 training report dir" "$report_dir"
  mkdir -p "$report_dir"
  make_root_manifest "$search_root" "$manifest"
  validate_training_tree "$manifest" "${report_dir}/training_tree_validation_3x3.json"

  local yaml_model_dir=""
  yaml_model_dir="$(validation_model_dir_from_yaml || true)"
  local stamp="${RJ_AUAU_TIGHT_BDT_TRAIN_STAMP:-$(ts)}"
  local model_dir="${RJ_AUAU_TIGHT_BDT_MODEL_DIR:-${yaml_model_dir:-${MODEL_BASE}/tight_3x3_${stamp}}}"
  guard_generated_path "3x3 model dir" "$model_dir"
  log_path_plan "trainCentInput3x3FromExtraction" \
    "source    : ${source}" \
    "model dir : ${model_dir}" \
    "manifest  : ${manifest}" \
    "report    : ${report_dir}" \
    "model id  : centAsFeat3x3_pt5to40"
  mkdir -p "$model_dir"
  setup_ml_python_env
  local registry="${model_dir}/model_registry_3x3.json"
  "$ML_PYTHON" "$TRAIN_SCRIPT" --task tight --campaign expanded-tight \
    --input "@${manifest}" \
    --outdir "$model_dir" \
    --cache-file "${model_dir}/training_matrix_3x3.npz" \
    --registry-output "$registry" \
    --campaign-spec-ids centAsFeat3x3_pt5to40 \
    --parallel-workers 1 \
    --n-jobs "${RJ_AUAU_BDT_XGB_N_JOBS:-4}" \
    --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}"
  apply_check_registry "$model_dir" "$registry"
  local summary="${report_dir}/cent_input_3x3_training_summary.txt"
  {
    echo "RECOILJETS_AUAU_TIGHT_BDT_CENT_INPUT_3X3_TRAINING_V1"
    echo "status=READY"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "model_file=${model_dir}/auau_tight_bdt_centAsFeat3x3_pt5to40_tmva.root"
    echo "registry=${registry}"
    echo "next_action=Rebuild src_AuAu if not already rebuilt, then run isSimEmbedded/isSimEmbeddedInclusive with macros/analysis_config_auau_bdt_validation.yaml."
  } > "$summary"
  cat "$summary"
  send_summary_email "[RecoilJets][auauTightBDT_trainCentInput3x3FromExtraction][READY]" "$summary"
}

train_width_study_pt1530_from_extraction() {
  local source=""
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*) export RJ_AUAU_TIGHT_BDT_MODEL_DIR="${tok#MODEL_DIR=}" ;;
    esac
  done
  [[ -n "$source" ]] || die "trainWidthStudyPt1530FromExtraction requires SOURCE=/path"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"

  local manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  local report_dir="${source}/reports"
  guard_generated_path "width-study training report dir" "$report_dir"
  mkdir -p "$report_dir"
  make_root_manifest "$search_root" "$manifest"
  validate_training_tree "$manifest" "${report_dir}/training_tree_validation_widthstudy_pt1530.json"

  local stamp="${RJ_AUAU_TIGHT_BDT_TRAIN_STAMP:-$(ts)}"
  local default_model_dir="/gpfs/mnt/gpfs02/sphenix/user/${USER:-patsfan753}/thesisAnalysis/bdt_models/tight_centinput_widthstudy_pt1530_current"
  local model_dir="${RJ_AUAU_TIGHT_BDT_MODEL_DIR:-$default_model_dir}"
  local cache_file="${RJ_AUAU_BDT_CACHE_FILE:-${model_dir}/training_matrix_widthstudy_pt1530.npz}"
  local registry="${model_dir}/model_registry.json"
  local spec_ids="centAsFeat_pt15to30,centAsFeat3x3_pt15to30,centAsFeatBase3x3_pt15to30"
  guard_generated_path "width-study model dir" "$model_dir"
  log_path_plan "trainWidthStudyPt1530FromExtraction" \
    "source    : ${source}" \
    "model dir : ${model_dir}" \
    "manifest  : ${manifest}" \
    "cache     : ${cache_file}" \
    "registry  : ${registry}" \
    "spec ids  : ${spec_ids}" \
    "pt window : 15 <= cluster_Et < 30 GeV"
  mkdir -p "$model_dir"
  setup_ml_python_env
  "$ML_PYTHON" "$TRAIN_SCRIPT" --task tight --campaign expanded-tight \
    --input "@${manifest}" \
    --outdir "$model_dir" \
    --cache-file "$cache_file" \
    --registry-output "$registry" \
    --extra-cent-as-feat-pt-ranges 15:30 \
    --extra-cent-as-feat-3x3-pt-ranges 15:30 \
    --extra-cent-as-feat-base3x3-pt-ranges 15:30 \
    --campaign-spec-ids "$spec_ids" \
    --parallel-workers "${RJ_AUAU_BDT_TRAIN_PARALLEL:-3}" \
    --n-jobs "${RJ_AUAU_BDT_XGB_N_JOBS:-1}" \
    --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}"
  apply_check_registry "$model_dir" "$registry"
  local summary="${report_dir}/widthstudy_pt1530_training_summary.txt"
  {
    echo "RECOILJETS_AUAU_TIGHT_BDT_WIDTHSTUDY_PT1530_TRAINING_V1"
    echo "status=READY"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "registry=${registry}"
    echo "pt_window=15:30"
    echo "model_base_widths=${model_dir}/auau_tight_bdt_centAsFeat_pt15to30_tmva.root"
    echo "model_3x3_widths=${model_dir}/auau_tight_bdt_centAsFeat3x3_pt15to30_tmva.root"
    echo "model_base_plus_3x3_widths=${model_dir}/auau_tight_bdt_centAsFeatBase3x3_pt15to30_tmva.root"
    echo "next_action=Run validateOnSimCondor with this MODEL_DIR, then run the strict 15-30 GeV WP0.80 MC validation YAML."
  } > "$summary"
  cat "$summary"
  send_summary_email "[RecoilJets][auauTightBDT_trainWidthStudyPt1530FromExtraction][READY]" "$summary"
}

train_width_study_windows_from_extraction() {
  local source=""
  local pt_windows="${RJ_AUAU_BDT_WIDTH_WINDOWS:-5:35,10:35,15:35}"
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*) export RJ_AUAU_TIGHT_BDT_MODEL_DIR="${tok#MODEL_DIR=}" ;;
      PT_WINDOWS=*) pt_windows="${tok#PT_WINDOWS=}" ;;
    esac
  done
  [[ -n "$source" ]] || die "trainWidthStudyWindowsFromExtraction requires SOURCE=/path"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"

  local manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  local report_dir="${source}/reports"
  guard_generated_path "width-window training report dir" "$report_dir"
  mkdir -p "$report_dir"
  make_root_manifest "$search_root" "$manifest"
  validate_training_tree "$manifest" "${report_dir}/training_tree_validation_widthstudy_windows.json"

  local default_model_dir="/gpfs/mnt/gpfs02/sphenix/user/${USER:-patsfan753}/thesisAnalysis/bdt_models/tight_centinput_widthstudy_windows_current"
  local model_dir="${RJ_AUAU_TIGHT_BDT_MODEL_DIR:-$default_model_dir}"
  local cache_file="${RJ_AUAU_BDT_CACHE_FILE:-${model_dir}/training_matrix_widthstudy_windows.npz}"
  local registry="${model_dir}/model_registry.json"
  guard_generated_path "width-window model dir" "$model_dir"

  local extra_ranges=""
  local spec_ids=""
  local summary_models=()
  IFS=',' read -r -a window_items <<< "$pt_windows"
  for raw in "${window_items[@]}"; do
    local item="${raw//[[:space:]]/}"
    [[ -n "$item" ]] || continue
    [[ "$item" == *:* ]] || die "PT_WINDOWS item must be lo:hi, got: $item"
    local lo="${item%%:*}"
    local hi="${item##*:}"
    local lo_tag="${lo%.*}"
    local hi_tag="${hi%.*}"
    [[ -n "$lo_tag" && -n "$hi_tag" ]] || die "Bad PT_WINDOWS item: $item"
    local tag="pt${lo_tag}to${hi_tag}"
    extra_ranges+="${extra_ranges:+,}${lo}:${hi}"
    spec_ids+="${spec_ids:+,}centAsFeat_${tag},centAsFeat3x3_${tag},centAsFeatBase3x3_${tag}"
    summary_models+=("${tag}: ${model_dir}/auau_tight_bdt_centAsFeat_${tag}_tmva.root")
    summary_models+=("${tag}: ${model_dir}/auau_tight_bdt_centAsFeat3x3_${tag}_tmva.root")
    summary_models+=("${tag}: ${model_dir}/auau_tight_bdt_centAsFeatBase3x3_${tag}_tmva.root")
  done
  [[ -n "$extra_ranges" ]] || die "No PT_WINDOWS were provided"

  log_path_plan "trainWidthStudyWindowsFromExtraction" \
    "source     : ${source}" \
    "model dir  : ${model_dir}" \
    "manifest   : ${manifest}" \
    "cache      : ${cache_file}" \
    "registry   : ${registry}" \
    "pt windows : ${extra_ranges}" \
    "spec ids   : ${spec_ids}"
  mkdir -p "$model_dir"
  setup_ml_python_env
  "$ML_PYTHON" "$TRAIN_SCRIPT" --task tight --campaign expanded-tight \
    --input "@${manifest}" \
    --outdir "$model_dir" \
    --cache-file "$cache_file" \
    --registry-output "$registry" \
    --extra-cent-as-feat-pt-ranges "$extra_ranges" \
    --extra-cent-as-feat-3x3-pt-ranges "$extra_ranges" \
    --extra-cent-as-feat-base3x3-pt-ranges "$extra_ranges" \
    --campaign-spec-ids "$spec_ids" \
    --parallel-workers "${RJ_AUAU_BDT_TRAIN_PARALLEL:-3}" \
    --n-jobs "${RJ_AUAU_BDT_XGB_N_JOBS:-1}" \
    --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}"
  apply_check_registry "$model_dir" "$registry"
  local summary="${report_dir}/widthstudy_windows_training_summary.txt"
  {
    echo "RECOILJETS_AUAU_TIGHT_BDT_WIDTHSTUDY_WINDOWS_TRAINING_V1"
    echo "status=READY"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "registry=${registry}"
    echo "pt_windows=${extra_ranges}"
    echo "trained_specs=${spec_ids}"
    for model_line in "${summary_models[@]}"; do
      echo "model=${model_line}"
    done
    echo "next_action=Run validateOnSimCondor with this MODEL_DIR, then run scripts/submit_auau_bdt_widthstudy_windows_wp050.sh for paired MC validation."
  } > "$summary"
  cat "$summary"
  send_summary_email "[RecoilJets][auauTightBDT_trainWidthStudyWindowsFromExtraction][READY]" "$summary"
}

train_etfine_centstudy_from_extraction() {
  local source=""
  local pt_bins="${RJ_AUAU_BDT_ETFINE_PT_BINS:-15,17,19,21,23,25,27,30,35}"
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*) export RJ_AUAU_TIGHT_BDT_MODEL_DIR="${tok#MODEL_DIR=}" ;;
      PT_BINS=*) pt_bins="${tok#PT_BINS=}" ;;
    esac
  done
  [[ -n "$source" ]] || die "trainEtFineCentStudyFromExtraction requires SOURCE=/path"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"

  local manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  local report_dir="${source}/reports"
  guard_generated_path "fine-E_T training report dir" "$report_dir"
  mkdir -p "$report_dir"
  make_root_manifest "$search_root" "$manifest"
  validate_training_tree "$manifest" "${report_dir}/training_tree_validation_etfine_centstudy.json"

  local default_model_dir="/gpfs/mnt/gpfs02/sphenix/user/${USER:-patsfan753}/thesisAnalysis/bdt_models/tight_etfine_centstudy_current"
  local model_dir="${RJ_AUAU_TIGHT_BDT_MODEL_DIR:-$default_model_dir}"
  local cache_file="${RJ_AUAU_BDT_CACHE_FILE:-${model_dir}/training_matrix_etfine_centstudy.npz}"
  local registry="${model_dir}/model_registry.json"
  guard_generated_path "fine-E_T model dir" "$model_dir"

  log_path_plan "trainEtFineCentStudyFromExtraction" \
    "source     : ${source}" \
    "model dir  : ${model_dir}" \
    "manifest   : ${manifest}" \
    "cache      : ${cache_file}" \
    "registry   : ${registry}" \
    "pt bins    : ${pt_bins}" \
    "products   : centInput_pt1535, ptFine_centInput, ptFine_cent3, ptFine_cent7"
  mkdir -p "$model_dir"
  setup_ml_python_env
  "$ML_PYTHON" "$TRAIN_SCRIPT" --task tight --campaign etfine-centstudy \
    --input "@${manifest}" \
    --outdir "$model_dir" \
    --cache-file "$cache_file" \
    --registry-output "$registry" \
    --pt-bins "$pt_bins" \
    --coarse-cent-bins "${RJ_AUAU_BDT_ETFINE_COARSE_CENT_BINS:-0:20,20:50,50:80}" \
    --fine-cent-bins "${RJ_AUAU_BDT_ETFINE_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}" \
    --parallel-workers "${RJ_AUAU_BDT_TRAIN_PARALLEL:-1}" \
    --n-jobs "${RJ_AUAU_BDT_XGB_N_JOBS:-1}" \
    --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}"
  apply_check_registry "$model_dir" "$registry"
  local summary="${report_dir}/etfine_centstudy_training_summary.txt"
  {
    echo "RECOILJETS_AUAU_TIGHT_BDT_ETFINE_CENTSTUDY_TRAINING_V1"
    echo "status=READY"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "registry=${registry}"
    echo "pt_bins=${pt_bins}"
    echo "expected_models=89"
    echo "model_centInput_pt1535=${model_dir}/auau_tight_bdt_centInput_pt1535_tmva.root"
    echo "next_action=Run validateOnSimCondor with this MODEL_DIR, then run the strict 15-35 GeV WP0.50 MC validation YAML."
  } > "$summary"
  cat "$summary"
  send_summary_email "[RecoilJets][auauTightBDT_trainEtFineCentStudyFromExtraction][READY]" "$summary"
}

train_expanded_from_extraction() {
  local source=""
  local plan_only=0
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      PLAN_ONLY=1|planOnly=1|DRYRUN=1) plan_only=1 ;;
    esac
  done
  [[ -n "$source" ]] || die "trainExpandedFromExtraction requires SOURCE=/path"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"
  local manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  local report_dir="${source}/reports"
  guard_generated_path "expanded training report dir" "$report_dir"
  mkdir -p "$report_dir"
  make_root_manifest "$search_root" "$manifest"
  validate_training_tree "$manifest" "${report_dir}/training_tree_validation.json"

  local stamp="${RJ_AUAU_TIGHT_BDT_TRAIN_STAMP:-$(ts)}"
  local model_dir="${RJ_AUAU_TIGHT_BDT_MODEL_DIR:-${MODEL_BASE}/tight_expanded_${stamp}}"
  local cache_file="${RJ_AUAU_BDT_CACHE_FILE:-${model_dir}/training_matrix.npz}"
  local spec_ids="${RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS:-}"
  local campaign="${RJ_AUAU_BDT_CAMPAIGN:-expanded-tight}"
  local weight_mode="${RJ_AUAU_BDT_WEIGHT_MODE:-legacy}"
  local label_contract="${RJ_AUAU_BDT_LABEL_CONTRACT:-extracted-is-signal}"
  local etcent_pt_bins="${RJ_AUAU_BDT_ETCENT_PT_BINS:-15,17,19,21,23,25,27,30,35}"
  local etcent_coarse_cent_bins="${RJ_AUAU_BDT_ETCENT_COARSE_CENT_BINS:-0:20,20:50,50:80}"
  local etcent_fine_cent_bins="${RJ_AUAU_BDT_ETCENT_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}"
  local etfine_pt_bins="${RJ_AUAU_BDT_ETFINE_PT_BINS:-15,17,19,21,23,25,27,30,35}"
  local etfine_coarse_cent_bins="${RJ_AUAU_BDT_ETFINE_COARSE_CENT_BINS:-0:20,20:50,50:80}"
  local etfine_fine_cent_bins="${RJ_AUAU_BDT_ETFINE_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}"
  # Single-pt-window base+3x3-width spec (the centAsFeatBase3x3 baseline) is only built when
  # this pt-range flag is forwarded to the trainer. Exposing it lets the DEFAULT Au+Au baseline
  # train under the canonical name centAsFeatBase3x3_pt15to35 (the SAME 14 features as
  # etfine-centstudy's centInput_pt1535). Defaults empty => no change for existing campaigns.
  local extra_base3x3_pt_ranges="${RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES:-}"
  local ppg12_expected_samples="${RJ_AUAU_BDT_PPG12_EXACT_EXPECTED_SAMPLES:-run28_embeddedPhoton12,run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40}"
  local ppg12_closure_dir="${RJ_AUAU_BDT_PPG12_EXACT_CLOSURE_DIR:-${model_dir}/slideReady/ppg12_exact_reweight_bdt}"
  local ppg12_closure_artifacts="${RJ_AUAU_BDT_PPG12_EXACT_CLOSURE_ARTIFACTS:-full}"
  local bdt_test_size="${RJ_AUAU_BDT_TEST_SIZE:-0.10}"
  local bdt_split_mode="${RJ_AUAU_BDT_SPLIT_MODE:-row}"
  local bdt_max_depth="${RJ_AUAU_BDT_MAX_DEPTH:-}"
  local bdt_random_seed="${RJ_AUAU_BDT_RANDOM_SEED:-}"
  local bdt_n_estimators="${RJ_AUAU_BDT_N_ESTIMATORS:-}"
  local bdt_learning_rate="${RJ_AUAU_BDT_LEARNING_RATE:-}"
  local bdt_subsample="${RJ_AUAU_BDT_SUBSAMPLE:-}"
  local bdt_colsample_bytree="${RJ_AUAU_BDT_COLSAMPLE_BYTREE:-}"
  local bdt_reg_alpha="${RJ_AUAU_BDT_REG_ALPHA:-}"
  local bdt_reg_lambda="${RJ_AUAU_BDT_REG_LAMBDA:-}"
  local bdt_max_bin="${RJ_AUAU_BDT_MAX_BIN:-}"
  local event_quality_cut_json="${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON:-}"
  local event_quality_audit_output="${RJ_AUAU_BDT_EVENT_QUALITY_AUDIT_OUTPUT:-${model_dir}/event_quality_filter_audit.json}"
  local event_quality_audit_only="${RJ_AUAU_BDT_EVENT_QUALITY_AUDIT_ONLY:-0}"
  local raw_eiso_pt_bins="${RJ_AUAU_BDT_EISO_CONE_PT_BINS:-15,17,19,21,23,25,27,30,35}"
  local raw_eiso_coarse_cent_bins="${RJ_AUAU_BDT_EISO_CONE_COARSE_CENT_BINS:-0:20,20:50,50:80}"
  local raw_eiso_fine_cent_bins="${RJ_AUAU_BDT_EISO_CONE_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}"
  [[ "$label_contract" == "extracted-is-signal" || "$label_contract" == "nominal-isolated-prompt" || "$label_contract" == "ppg12-source-role" ]] || die "RJ_AUAU_BDT_LABEL_CONTRACT must be extracted-is-signal, nominal-isolated-prompt, or ppg12-source-role"
  guard_generated_path "expanded training model dir" "$model_dir"
  log_path_plan "trainExpandedFromExtraction" \
    "source    : ${source}" \
    "model dir : ${model_dir}" \
    "manifest  : ${manifest}" \
    "cache     : ${cache_file}" \
    "campaign  : ${campaign}" \
    "weight    : ${weight_mode}" \
    "labels    : ${label_contract}" \
    "max depth : ${bdt_max_depth:-<trainer default>}" \
    "seed      : ${bdt_random_seed:-<trainer default>}" \
    "estimators: ${bdt_n_estimators:-<trainer default>}" \
    "learn rate: ${bdt_learning_rate:-<trainer default>}" \
    "spec ids  : ${spec_ids:-<all>}" \
    "event cut : ${event_quality_cut_json:-<disabled>}" \
    "report    : ${report_dir}" \
    "plan only : ${plan_only}"
  mkdir -p "$model_dir"
  setup_ml_python_env
  local -a args=(
    "$TRAIN_SCRIPT"
    --task tight
    --campaign "$campaign"
    --input "@${manifest}"
    --outdir "$model_dir"
    --cache-file "${cache_file}"
    --registry-output "${model_dir}/model_registry.json"
    --weight-mode "${weight_mode}"
    --label-contract "${label_contract}"
    --test-size "${bdt_test_size}"
    --split-mode "${bdt_split_mode}"
    --parallel-workers "${RJ_AUAU_BDT_TRAIN_PARALLEL:-4}"
    --n-jobs "${RJ_AUAU_BDT_XGB_N_JOBS:-1}"
    --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}"
    --minopt-majority-cap-ratio "${RJ_AUAU_BDT_MINOPT_MAJORITY_CAP_RATIO:-2.0}"
  )
  if [[ "$weight_mode" == "ppg12-exact" ]]; then
    args+=(
      --no-event-weight
      --ppg12-exact-expected-samples "$ppg12_expected_samples"
      --ppg12-exact-closure-dir "$ppg12_closure_dir"
      --ppg12-exact-closure-artifacts "$ppg12_closure_artifacts"
    )
  fi
  if [[ "$campaign" == "etcent-binned-sixpack" || "$campaign" == "etcent-binned-sixpack-noiso-ptcent7" || "$campaign" == "global-and-etcent-binned-sixpack-noiso" ]]; then
    args+=(
      --pt-bins "$etcent_pt_bins"
      --coarse-cent-bins "$etcent_coarse_cent_bins"
      --fine-cent-bins "$etcent_fine_cent_bins"
    )
  fi
  if [[ "$campaign" == "etfine-centstudy" || "$campaign" == "corrected-baseline-binned14" ]]; then
    args+=(
      --pt-bins "$etfine_pt_bins"
      --coarse-cent-bins "$etfine_coarse_cent_bins"
      --fine-cent-bins "$etfine_fine_cent_bins"
    )
  fi
  if [[ "$campaign" == "etcent-binned-eiso-cone-ablation" ]]; then
    args+=(
      --pt-bins "$raw_eiso_pt_bins"
      --coarse-cent-bins "$raw_eiso_coarse_cent_bins"
      --fine-cent-bins "$raw_eiso_fine_cent_bins"
    )
  fi
  if [[ -n "$extra_base3x3_pt_ranges" ]]; then
    args+=( --extra-cent-as-feat-base3x3-pt-ranges "$extra_base3x3_pt_ranges" )
  fi
  if [[ -n "$spec_ids" ]]; then
    args+=( --campaign-spec-ids "$spec_ids" )
  fi
  if [[ -n "$bdt_max_depth" ]]; then
    args+=( --max-depth "$bdt_max_depth" )
  fi
  if [[ -n "$bdt_random_seed" ]]; then args+=( --random-seed "$bdt_random_seed" ); fi
  if [[ -n "$bdt_n_estimators" ]]; then args+=( --n-estimators "$bdt_n_estimators" ); fi
  if [[ -n "$bdt_learning_rate" ]]; then args+=( --learning-rate "$bdt_learning_rate" ); fi
  if [[ -n "$bdt_subsample" ]]; then args+=( --subsample "$bdt_subsample" ); fi
  if [[ -n "$bdt_colsample_bytree" ]]; then args+=( --colsample-bytree "$bdt_colsample_bytree" ); fi
  if [[ -n "$bdt_reg_alpha" ]]; then args+=( --reg-alpha "$bdt_reg_alpha" ); fi
  if [[ -n "$bdt_reg_lambda" ]]; then args+=( --reg-lambda "$bdt_reg_lambda" ); fi
  if [[ -n "$bdt_max_bin" ]]; then args+=( --max-bin "$bdt_max_bin" ); fi
  if [[ -n "$event_quality_cut_json" ]]; then
    args+=( --event-quality-cut-json "$event_quality_cut_json" )
    args+=( --event-quality-audit-output "$event_quality_audit_output" )
    args+=( --require-global-event-key )
    if [[ "$event_quality_audit_only" == "1" ]]; then
      args+=( --event-quality-audit-only )
    fi
  fi
  if [[ "$plan_only" == "1" ]]; then
    args+=( --plan-only --registry-output "${model_dir}/model_registry.planned.json" )
  fi
  say "Expanded training model dir: $model_dir"
  "$ML_PYTHON" "${args[@]}"
  if [[ "$plan_only" == "1" ]]; then
    say "Expanded campaign plan only: ${model_dir}/model_registry.planned.json"
    return 0
  fi
  apply_check_registry "$model_dir"
  local summary="${report_dir}/expanded_training_summary.txt"
  {
    echo "RECOILJETS_AUAU_TIGHT_BDT_EXPANDED_TRAINING_V1"
    echo "status=READY"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "registry=${model_dir}/model_registry.json"
    echo "next_action=Run validateOnSimCondor with MODEL_DIR=${model_dir}, then inspect registry-ranked model QA."
  } > "$summary"
  cat "$summary"
  send_summary_email "[RecoilJets][auauTightBDT_trainExpandedFromExtraction][READY]" "$summary"
}

train_expanded_from_extraction_condor() {
  local source=""
  local group_size="${RJ_AUAU_BDT_EXPANDED_GROUP_SIZE:-4}"
  local reqmem="${RJ_AUAU_BDT_EXPANDED_REQUEST_MEMORY:-6000MB}"
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      groupSize=*) group_size="${tok#groupSize=}" ;;
    esac
  done
  while (($#)); do
    case "$1" in
      groupSize) group_size="${2:?missing value after groupSize}"; shift 2 ;;
      *) shift ;;
    esac
  done
  [[ -n "$source" ]] || die "trainExpandedFromExtractionCondor requires SOURCE=/path"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"
  [[ "$group_size" =~ ^[0-9]+$ && "$group_size" -gt 0 ]] || die "groupSize must be a positive integer"

  local manifest_override="${RJ_AUAU_BDT_ROOT_MANIFEST:-}"
  local manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  local report_dir="${RJ_AUAU_BDT_REPORT_DIR:-${source}/reports}"
  mkdir -p "$report_dir"
  if [[ -n "$manifest_override" ]]; then
    [[ -s "$manifest_override" ]] || die "RJ_AUAU_BDT_ROOT_MANIFEST is empty or missing: $manifest_override"
    manifest="$manifest_override"
  else
    make_root_manifest "$search_root" "$manifest"
  fi
  validate_training_tree "$manifest" "${report_dir}/training_tree_validation.json"

  local stamp="${RJ_AUAU_TIGHT_BDT_TRAIN_STAMP:-$(ts)}"
  local model_dir="${RJ_AUAU_TIGHT_BDT_MODEL_DIR:-${MODEL_BASE}/tight_expanded_${stamp}}"
  local sub_root="${RJ_REPO_BASE}/condor_sub/auauTightBDTExpanded_${stamp}"
  local shard_dir="${sub_root}/spec_shards"
  local registry_dir="${sub_root}/registries"
  local cache_file="${RJ_AUAU_BDT_CACHE_FILE:-${model_dir}/training_matrix.npz}"
  local spec_ids="${RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS:-}"
  local campaign="${RJ_AUAU_BDT_CAMPAIGN:-expanded-tight}"
  local weight_mode="${RJ_AUAU_BDT_WEIGHT_MODE:-legacy}"
  local label_contract="${RJ_AUAU_BDT_LABEL_CONTRACT:-extracted-is-signal}"
  local etcent_pt_bins="${RJ_AUAU_BDT_ETCENT_PT_BINS:-15,17,19,21,23,25,27,30,35}"
  local etcent_coarse_cent_bins="${RJ_AUAU_BDT_ETCENT_COARSE_CENT_BINS:-0:20,20:50,50:80}"
  local etcent_fine_cent_bins="${RJ_AUAU_BDT_ETCENT_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}"
  local etfine_pt_bins="${RJ_AUAU_BDT_ETFINE_PT_BINS:-15,17,19,21,23,25,27,30,35}"
  local etfine_coarse_cent_bins="${RJ_AUAU_BDT_ETFINE_COARSE_CENT_BINS:-0:20,20:50,50:80}"
  local etfine_fine_cent_bins="${RJ_AUAU_BDT_ETFINE_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}"
  # Single-pt-window base+3x3-width spec (the centAsFeatBase3x3 baseline) is only built when
  # this pt-range flag is forwarded to the trainer. Exposing it lets the DEFAULT Au+Au baseline
  # train under the canonical name centAsFeatBase3x3_pt15to35 (the SAME 14 features as
  # etfine-centstudy's centInput_pt1535). Defaults empty => no change for existing campaigns.
  local extra_base3x3_pt_ranges="${RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES:-}"
  local ppg12_expected_samples="${RJ_AUAU_BDT_PPG12_EXACT_EXPECTED_SAMPLES:-run28_embeddedPhoton12,run28_embeddedPhoton20,run28_embeddedJet12,run28_embeddedJet20,run28_embeddedJet30,run28_embeddedJet40}"
  local ppg12_closure_dir="${RJ_AUAU_BDT_PPG12_EXACT_CLOSURE_DIR:-${model_dir}/slideReady/ppg12_exact_reweight_bdt}"
  local ppg12_closure_artifacts="${RJ_AUAU_BDT_PPG12_EXACT_CLOSURE_ARTIFACTS:-full}"
  local bdt_test_size="${RJ_AUAU_BDT_TEST_SIZE:-0.10}"
  local bdt_split_mode="${RJ_AUAU_BDT_SPLIT_MODE:-row}"
  local bdt_max_depth="${RJ_AUAU_BDT_MAX_DEPTH:-}"
  local bdt_random_seed="${RJ_AUAU_BDT_RANDOM_SEED:-}"
  local bdt_n_estimators="${RJ_AUAU_BDT_N_ESTIMATORS:-}"
  local bdt_learning_rate="${RJ_AUAU_BDT_LEARNING_RATE:-}"
  local bdt_subsample="${RJ_AUAU_BDT_SUBSAMPLE:-}"
  local bdt_colsample_bytree="${RJ_AUAU_BDT_COLSAMPLE_BYTREE:-}"
  local bdt_reg_alpha="${RJ_AUAU_BDT_REG_ALPHA:-}"
  local bdt_reg_lambda="${RJ_AUAU_BDT_REG_LAMBDA:-}"
  local bdt_max_bin="${RJ_AUAU_BDT_MAX_BIN:-}"
  local event_quality_cut_json="${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON:-}"
  local event_quality_audit_output="${RJ_AUAU_BDT_EVENT_QUALITY_AUDIT_OUTPUT:-${model_dir}/event_quality_filter_audit.json}"
  local event_quality_assume_filtered="${RJ_AUAU_BDT_EVENT_QUALITY_ASSUME_FILTERED:-0}"
  local raw_eiso_pt_bins="${RJ_AUAU_BDT_EISO_CONE_PT_BINS:-15,17,19,21,23,25,27,30,35}"
  local raw_eiso_coarse_cent_bins="${RJ_AUAU_BDT_EISO_CONE_COARSE_CENT_BINS:-0:20,20:50,50:80}"
  local raw_eiso_fine_cent_bins="${RJ_AUAU_BDT_EISO_CONE_FINE_CENT_BINS:-0:10,10:20,20:30,30:40,40:50,50:60,60:80}"
  if [[ ( "$campaign" == "etcent-binned-eiso-cone-ablation" || "$weight_mode" == "ppg12-exact" ) && -z "${RJ_AUAU_BDT_EXPANDED_REQUEST_MEMORY:-}" ]]; then
    reqmem="16000MB"
  fi
  local cache_reqmem="${RJ_AUAU_BDT_CACHE_REQUEST_MEMORY:-8000MB}"
  if [[ ( "$campaign" == "etcent-binned-eiso-cone-ablation" || "$weight_mode" == "ppg12-exact" ) && -z "${RJ_AUAU_BDT_CACHE_REQUEST_MEMORY:-}" ]]; then
    cache_reqmem="16000MB"
  fi
  local staged_cache="${RJ_AUAU_BDT_STAGED_CACHE:-auto}"
  if [[ "$staged_cache" == "auto" ]]; then
    if [[ "$weight_mode" == "ppg12-exact" ]]; then
      staged_cache="1"
    else
      staged_cache="0"
    fi
  fi
  local cache_shards="${RJ_AUAU_BDT_CACHE_SHARDS:-12}"
  local cache_part_reqmem="${RJ_AUAU_BDT_CACHE_PART_REQUEST_MEMORY:-6000MB}"
  local cache_reduce_reqmem="${RJ_AUAU_BDT_CACHE_REDUCE_REQUEST_MEMORY:-12000MB}"
  local cache_part_maxjobs="${RJ_AUAU_BDT_CACHE_PART_MAXJOBS:-4}"
  local reuse_existing_cache="${RJ_AUAU_BDT_REUSE_EXISTING_CACHE:-0}"
  local expected_cache_sha256="${RJ_AUAU_BDT_EXPECT_CACHE_SHA256:-}"
  [[ "$staged_cache" == "0" || "$staged_cache" == "1" ]] || die "RJ_AUAU_BDT_STAGED_CACHE must be 0, 1, or auto"
  [[ "$reuse_existing_cache" == "0" || "$reuse_existing_cache" == "1" ]] || die "RJ_AUAU_BDT_REUSE_EXISTING_CACHE must be 0 or 1"
  [[ "$label_contract" == "extracted-is-signal" || "$label_contract" == "nominal-isolated-prompt" || "$label_contract" == "ppg12-source-role" ]] || die "RJ_AUAU_BDT_LABEL_CONTRACT must be extracted-is-signal, nominal-isolated-prompt, or ppg12-source-role"
  if [[ "$label_contract" == "ppg12-source-role" && "$staged_cache" == "1" ]]; then
    die "The PPG12 source-role label contract requires RJ_AUAU_BDT_STAGED_CACHE=0 with a copied nominal cache; the staged reducer intentionally preserves the extracted label population."
  fi
  if [[ "$reuse_existing_cache" == "1" ]]; then
    [[ "$staged_cache" == "0" ]] || die "RJ_AUAU_BDT_REUSE_EXISTING_CACHE=1 requires RJ_AUAU_BDT_STAGED_CACHE=0"
    [[ -s "$cache_file" ]] || die "Requested existing training cache is missing/empty: $cache_file"
    [[ "$expected_cache_sha256" =~ ^[0-9a-fA-F]{64}$ ]] || die "RJ_AUAU_BDT_EXPECT_CACHE_SHA256 must be the expected 64-character SHA-256 when reusing a cache"
    local observed_cache_sha256
    observed_cache_sha256="$(sha256sum "$cache_file" | awk '{print $1}')"
    [[ "${observed_cache_sha256,,}" == "${expected_cache_sha256,,}" ]] || die "Existing training-cache SHA-256 mismatch: expected=${expected_cache_sha256} observed=${observed_cache_sha256} path=${cache_file}"
  fi
  [[ "$cache_shards" =~ ^[0-9]+$ && "$cache_shards" -gt 0 ]] || die "RJ_AUAU_BDT_CACHE_SHARDS must be a positive integer"
  [[ "$cache_part_maxjobs" =~ ^[0-9]+$ && "$cache_part_maxjobs" -gt 0 ]] || die "RJ_AUAU_BDT_CACHE_PART_MAXJOBS must be a positive integer"
  guard_generated_path "expanded training model dir" "$model_dir"
  guard_generated_path "expanded training submit root" "$sub_root"
  log_path_plan "trainExpandedFromExtractionCondor" \
    "source    : ${source}" \
    "manifest  : ${manifest}" \
    "model dir : ${model_dir}" \
    "submit    : ${sub_root}" \
    "shards    : ${shard_dir}" \
    "registries: ${registry_dir}" \
    "cache     : ${cache_file}" \
    "campaign  : ${campaign}" \
    "weight    : ${weight_mode}" \
    "labels    : ${label_contract}" \
    "test size : ${bdt_test_size}" \
    "split mode: ${bdt_split_mode}" \
    "max depth : ${bdt_max_depth:-<trainer default>}" \
    "seed      : ${bdt_random_seed:-<trainer default>}" \
    "estimators: ${bdt_n_estimators:-<trainer default>}" \
    "learn rate: ${bdt_learning_rate:-<trainer default>}" \
    "spec ids  : ${spec_ids:-<all>}" \
    "event cut : ${event_quality_cut_json:-<disabled>}" \
    "groupSize : ${group_size}" \
    "requestMem: ${reqmem}" \
    "cacheMem : ${cache_reqmem}" \
    "cacheMode: $([[ "$reuse_existing_cache" == "1" ]] && echo "reuse-existing sha256=${expected_cache_sha256}" || ([[ "$staged_cache" == "1" ]] && echo "staged shards=${cache_shards} partMem=${cache_part_reqmem} reduceMem=${cache_reduce_reqmem} maxJobs=${cache_part_maxjobs}" || echo "single"))"
  mkdir -p "$model_dir" "$sub_root" "$shard_dir" "$registry_dir"

  setup_ml_python_env
  if [[ "$reuse_existing_cache" == "1" ]]; then
    "$ML_PYTHON" - "$cache_file" "$expected_cache_sha256" "$model_dir/reused_cache_provenance.json" "$label_contract" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

cache = Path(sys.argv[1])
expected = sys.argv[2].lower()
output = Path(sys.argv[3])
label_contract = sys.argv[4]
digest = hashlib.sha256()
with cache.open("rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(chunk)
observed = digest.hexdigest()
if observed != expected:
    raise SystemExit(
        f"cache changed between shell and Python preflights: expected={expected} observed={observed}"
    )
output.write_text(
    json.dumps(
        {
            "schema": "AUAU_BDT_REUSED_TRAINING_CACHE_PROVENANCE_V1",
            "status": "VERIFIED_BEFORE_TRAINING",
            "path": str(cache),
            "size_bytes": cache.stat().st_size,
            "sha256": observed,
            "label_contract_to_apply_once": label_contract,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
PY
  fi
  local planned="${model_dir}/model_registry.planned.json"
  local -a plan_args=(
    "$TRAIN_SCRIPT" --task tight --campaign "$campaign"
    --input "@${manifest}" --outdir "$model_dir"
    --plan-only --registry-output "$planned"
    --weight-mode "$weight_mode"
    --label-contract "$label_contract"
    --test-size "$bdt_test_size"
    --split-mode "$bdt_split_mode"
    --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}"
    --minopt-majority-cap-ratio "${RJ_AUAU_BDT_MINOPT_MAJORITY_CAP_RATIO:-2.0}"
  )
  if [[ "$weight_mode" == "ppg12-exact" ]]; then
    plan_args+=(
      --no-event-weight
      --ppg12-exact-expected-samples "$ppg12_expected_samples"
      --ppg12-exact-closure-dir "$ppg12_closure_dir"
      --ppg12-exact-closure-artifacts "$ppg12_closure_artifacts"
    )
  fi
  if [[ "$campaign" == "etcent-binned-sixpack" || "$campaign" == "etcent-binned-sixpack-noiso-ptcent7" || "$campaign" == "global-and-etcent-binned-sixpack-noiso" ]]; then
    plan_args+=(
      --pt-bins "$etcent_pt_bins"
      --coarse-cent-bins "$etcent_coarse_cent_bins"
      --fine-cent-bins "$etcent_fine_cent_bins"
    )
  fi
  if [[ "$campaign" == "etfine-centstudy" || "$campaign" == "corrected-baseline-binned14" ]]; then
    plan_args+=(
      --pt-bins "$etfine_pt_bins"
      --coarse-cent-bins "$etfine_coarse_cent_bins"
      --fine-cent-bins "$etfine_fine_cent_bins"
    )
  fi
  if [[ "$campaign" == "etcent-binned-eiso-cone-ablation" ]]; then
    plan_args+=(
      --pt-bins "$raw_eiso_pt_bins"
      --coarse-cent-bins "$raw_eiso_coarse_cent_bins"
      --fine-cent-bins "$raw_eiso_fine_cent_bins"
    )
  fi
  if [[ -n "$extra_base3x3_pt_ranges" ]]; then
    plan_args+=( --extra-cent-as-feat-base3x3-pt-ranges "$extra_base3x3_pt_ranges" )
  fi
  if [[ -n "$spec_ids" ]]; then
    plan_args+=( --campaign-spec-ids "$spec_ids" )
  fi
  if [[ -n "$bdt_max_depth" ]]; then
    plan_args+=( --max-depth "$bdt_max_depth" )
  fi
  if [[ -n "$bdt_random_seed" ]]; then plan_args+=( --random-seed "$bdt_random_seed" ); fi
  if [[ -n "$bdt_n_estimators" ]]; then plan_args+=( --n-estimators "$bdt_n_estimators" ); fi
  if [[ -n "$bdt_learning_rate" ]]; then plan_args+=( --learning-rate "$bdt_learning_rate" ); fi
  if [[ -n "$bdt_subsample" ]]; then plan_args+=( --subsample "$bdt_subsample" ); fi
  if [[ -n "$bdt_colsample_bytree" ]]; then plan_args+=( --colsample-bytree "$bdt_colsample_bytree" ); fi
  if [[ -n "$bdt_reg_alpha" ]]; then plan_args+=( --reg-alpha "$bdt_reg_alpha" ); fi
  if [[ -n "$bdt_reg_lambda" ]]; then plan_args+=( --reg-lambda "$bdt_reg_lambda" ); fi
  if [[ -n "$bdt_max_bin" ]]; then plan_args+=( --max-bin "$bdt_max_bin" ); fi
  if [[ -n "$event_quality_cut_json" ]]; then
    plan_args+=( --event-quality-cut-json "$event_quality_cut_json" )
    plan_args+=( --event-quality-audit-output "$event_quality_audit_output" )
    plan_args+=( --require-global-event-key )
  fi
  "$ML_PYTHON" "${plan_args[@]}"

  if [[ "$campaign" == "etcent-binned-eiso-cone-ablation" ]]; then
    "$ML_PYTHON" - "$planned" <<'PY'
import json
import sys
from collections import Counter
from pathlib import Path

planned = json.loads(Path(sys.argv[1]).read_text())
models = planned.get("models", [])
expected_counts = {
    "globalEtCent1535_bdt_eisoR30_ptCent3": 24,
    "globalEtCent1535_bdt_eisoR30_ptCent7": 56,
    "globalEtCent1535_bdt_eisoR40_ptCent3": 24,
    "globalEtCent1535_bdt_eisoR40_ptCent7": 56,
    "globalEtCent1535_bdt_eisoR30R40_ptCent3": 24,
    "globalEtCent1535_bdt_eisoR30R40_ptCent7": 56,
}
products = Counter(model.get("product") for model in models)
low_pt = [
    model.get("model_id", "")
    for model in models
    if "_pt_005_" in model.get("model_id", "")
    or "_pt_008_" in model.get("model_id", "")
    or "_pt_010_" in model.get("model_id", "")
    or "_pt_012_" in model.get("model_id", "")
    or "_pt_014_" in model.get("model_id", "")
]
pt_bins = planned.get("pt_bins")
coarse = planned.get("coarse_cent_bins")
fine = planned.get("fine_cent_bins")
expected_pt = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0]
expected_coarse = [[0.0, 20.0], [20.0, 50.0], [50.0, 80.0]]
expected_fine = [[0.0, 10.0], [10.0, 20.0], [20.0, 30.0], [30.0, 40.0], [40.0, 50.0], [50.0, 60.0], [60.0, 80.0]]
errors = []
unknown_products = sorted(set(products) - set(expected_counts))
if not products:
    errors.append("no raw-eiso products planned")
if unknown_products:
    errors.append(f"unexpected products: {unknown_products}")
for product, count in sorted(products.items()):
    expected = expected_counts.get(product)
    if expected is not None and count != expected:
        errors.append(f"unexpected count for {product}: {count}, expected {expected}")
expected_selected = sum(expected_counts[product] for product in products if product in expected_counts)
if len(models) != expected_selected or planned.get("expected_model_count") != expected_selected:
    errors.append(
        "raw-eiso model count mismatch for selected products: "
        f"expected {expected_selected}, found models={len(models)} "
        f"expected_model_count={planned.get('expected_model_count')}"
    )
if low_pt:
    errors.append("low-pT model ids found: " + ", ".join(low_pt[:8]))
if pt_bins != expected_pt:
    errors.append(f"unexpected pt_bins: {pt_bins}")
if coarse != expected_coarse:
    errors.append(f"unexpected coarse_cent_bins: {coarse}")
if fine != expected_fine:
    errors.append(f"unexpected fine_cent_bins: {fine}")
if errors:
    raise SystemExit("Raw-eiso preflight failed:\n  " + "\n  ".join(errors))
print(
    "[OK] raw-eiso preflight: "
    f"{expected_selected} selected models, products={dict(products)}, "
    "15-35 GeV pT grid, 3/7 centrality grids"
)
PY
  fi

  if [[ "$weight_mode" == "ppg12-exact" ]]; then
    "$ML_PYTHON" - "$planned" <<'PY'
import json
import sys
from collections import Counter
from pathlib import Path

planned = json.loads(Path(sys.argv[1]).read_text())
models = planned.get("models", [])
products = Counter(model.get("product") for model in models)
low_pt = [
    model.get("model_id", "")
    for model in models
    if any(tag in model.get("model_id", "") for tag in ("_pt_005_", "_pt_008_", "_pt_010_", "_pt_012_", "_pt_014_"))
]
errors = []
if planned.get("defaults", {}).get("weight_mode") != "ppg12-exact":
    errors.append(f"unexpected weight_mode: {planned.get('defaults', {}).get('weight_mode')}")
expected_pt = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0]
extended_pt = [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0, 40.0]
ppg12_ext_pt = [5.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 32.0, 36.0, 40.0]
expected_coarse = [[0.0, 20.0], [20.0, 50.0], [50.0, 80.0]]
expected_fine = [[0.0, 10.0], [10.0, 20.0], [20.0, 30.0], [30.0, 40.0], [40.0, 50.0], [50.0, 60.0], [60.0, 80.0]]
planned_pt = planned.get("pt_bins") or expected_pt
allowed_pt_grids = (expected_pt, extended_pt, ppg12_ext_pt)
campaign = planned.get("campaign")
allow_low_pt = campaign == "corrected-baseline-binned14" and planned_pt == ppg12_ext_pt
if low_pt and not allow_low_pt:
    errors.append("low-pT model ids found: " + ", ".join(low_pt[:8]))
if planned.get("defaults", {}).get("event_weight_used_for_training") is not False:
    errors.append("event_weight_used_for_training is not false")
if planned.get("defaults", {}).get("no_cross_section_weights") is not True:
    errors.append("no_cross_section_weights is not true")
ok_message = ""
if campaign == "etcent-binned-sixpack-noiso-ptcent7":
    if dict(products) != {"globalEtCent1535_bdt_noIso_ptCent7": 56}:
        errors.append(f"unexpected products/counts: {dict(products)}")
    if len(models) != 56 or planned.get("expected_model_count") != 56:
        errors.append(
            f"expected exactly 56 routed models, found models={len(models)} "
            f"expected_model_count={planned.get('expected_model_count')}"
        )
    if planned.get("pt_bins") != expected_pt:
        errors.append(f"unexpected pt_bins: {planned.get('pt_bins')}")
    if planned.get("fine_cent_bins") != expected_fine:
        errors.append(f"unexpected fine_cent_bins: {planned.get('fine_cent_bins')}")
    ok_message = "[OK] PPG12-exact preflight: 56 no-isolation ptCent7 models, 15-35 GeV grid, no event/cross-section training weights"
elif campaign == "global-sixpack":
    if dict(products) != {"globalEtCent1535_bdt_noIso": 1}:
        errors.append(f"unexpected products/counts for global 32-input BDT: {dict(products)}")
    if len(models) != 1 or planned.get("expected_model_count") != 1:
        errors.append(
            f"expected exactly 1 global 32-input BDT, found models={len(models)} "
            f"expected_model_count={planned.get('expected_model_count')}"
        )
    if models:
        model = models[0]
        features = set(model.get("features") or [])
        required = {"cluster_Et", "centrality", "cluster_weta33_cogx", "cluster_wphi33_cogx"}
        missing = sorted(required - features)
        if model.get("model_id") != "globalEtCent1535_bdt_noIso":
            errors.append(f"unexpected model_id: {model.get('model_id')}")
        if len(model.get("features") or []) != 32:
            errors.append(f"global 32-input BDT feature count is not 32: {len(model.get('features') or [])}")
        if missing:
            errors.append(f"global 32-input BDT missing required feature(s): {missing}")
    ok_message = "[OK] PPG12-exact preflight: globalEtCent1535_bdt_noIso 32-input model, no event/cross-section training weights"
elif campaign in {"etcent-binned-sixpack", "global-and-etcent-binned-sixpack-noiso"}:
    expected_counts = {
        "globalEtCent1535_bdt_noIso": 1,
        "globalEtCent1535_bdt_noIso_ptCent3": 24,
        "globalEtCent1535_bdt_noIso_ptCent7": 56,
    }
    if campaign == "etcent-binned-sixpack":
        expected_counts.pop("globalEtCent1535_bdt_noIso")
    unknown = sorted(set(products) - set(expected_counts))
    if unknown:
        errors.append(f"unexpected products for requested no-isolation routed run: {unknown}")
    if not products:
        errors.append("no products planned for 32-input BDT run")
    for product, count in sorted(products.items()):
        expected = expected_counts.get(product)
        if expected is not None and count != expected:
            errors.append(f"unexpected count for {product}: {count}, expected {expected}")
    expected_selected = sum(expected_counts[p] for p in products if p in expected_counts)
    if len(models) != expected_selected or planned.get("expected_model_count") != expected_selected:
        errors.append(
            "32-input model count mismatch: "
            f"expected {expected_selected}, found models={len(models)} "
            f"expected_model_count={planned.get('expected_model_count')}"
        )
    if planned.get("pt_bins") != expected_pt:
        errors.append(f"unexpected pt_bins: {planned.get('pt_bins')}")
    if planned.get("coarse_cent_bins") != expected_coarse:
        errors.append(f"unexpected coarse_cent_bins: {planned.get('coarse_cent_bins')}")
    if planned.get("fine_cent_bins") != expected_fine:
        errors.append(f"unexpected fine_cent_bins: {planned.get('fine_cent_bins')}")
    bad_feature_counts = [
        f"{model.get('model_id')}:{len(model.get('features') or [])}"
        for model in models
        if len(model.get("features") or []) != 32
    ]
    if bad_feature_counts:
        errors.append("32-input feature count mismatch: " + ", ".join(bad_feature_counts[:8]))
    ok_message = (
        "[OK] PPG12-exact preflight: no-isolation 32-input products "
        f"{dict(products)}, 15-35 GeV grid, no event/cross-section training weights"
    )
elif campaign == "etfine-centstudy":
    product_counts = dict(products)
    allowed_counts = ({"centInput_pt1535": 1}, {"ptFine_cent7": 56})
    if product_counts not in allowed_counts:
        errors.append(
            "unexpected products/counts for etfine-centstudy PPG12-exact run: "
            f"{product_counts}; expected one of {allowed_counts}"
        )
    expected_model_count = sum(product_counts.values())
    if len(models) != expected_model_count or planned.get("expected_model_count") != expected_model_count:
        errors.append(
            f"expected {expected_model_count} selected etfine-centstudy model(s), "
            f"found models={len(models)} expected_model_count={planned.get('expected_model_count')}"
        )
    if product_counts == {"centInput_pt1535": 1} and models:
        model = models[0]
        features = set(model.get("features") or [])
        required = {"cluster_Et", "centrality", "cluster_weta33_cogx", "cluster_wphi33_cogx"}
        missing = sorted(required - features)
        if model.get("model_id") != "centInput_pt1535":
            errors.append(f"unexpected model_id: {model.get('model_id')}")
        if missing:
            errors.append(f"global BDT missing required feature(s): {missing}")
    elif product_counts == {"ptFine_cent7": 56}:
        if planned.get("pt_bins") != expected_pt:
            errors.append(f"unexpected pt_bins: {planned.get('pt_bins')}")
        if planned.get("fine_cent_bins") != expected_fine:
            errors.append(f"unexpected fine_cent_bins: {planned.get('fine_cent_bins')}")
        for model in models:
            features = set(model.get("features") or [])
            required = {"cluster_Et", "cluster_weta33_cogx", "cluster_wphi33_cogx"}
            missing = sorted(required - features)
            if model.get("product") != "ptFine_cent7":
                errors.append(f"unexpected model product: {model.get('model_id')} -> {model.get('product')}")
            if "centrality" in features:
                errors.append(f"routed model should not include centrality as a feature: {model.get('model_id')}")
            if missing:
                errors.append(f"routed model missing required feature(s): {model.get('model_id')} -> {missing}")
            if errors:
                break
    ok_message = (
        "[OK] PPG12-exact preflight: etfine-centstudy products "
        f"{product_counts}, 15-35 GeV grid, no event/cross-section training weights"
    )
elif campaign == "corrected-baseline-shower-ladder":
    if not models:
        errors.append("no corrected-baseline shower-ladder models planned")
    allowed_products = {"baseV3E11_pt1535", "baseV3E11_cent12_pt1535", "centAsFeatBase3x3_pt15to35", "globalEtCent1535_bdt_noIso", "base14_to32_cumulative_ladder"}
    for product in products:
        if not (
            product in allowed_products
            or product.startswith("base14_plus1_")
        ):
            errors.append(f"unexpected shower-ladder product: {product}")
    for model in models:
        model_id = model.get("model_id", "")
        feats = model.get("features") or []
        pt_range = model.get("pt_range")
        cent_range = model.get("cent_range")
        if pt_range != [15.0, 35.0] or cent_range is not None:
            errors.append(f"unexpected shower-ladder phase space: {model_id} pt={pt_range} cent={cent_range}")
        if model_id == "baseV3E11_pt1535" and len(feats) != 11:
            errors.append(f"baseV3E11_pt1535 feature count is not 11: {len(feats)}")
        elif model_id == "baseV3E11_cent12_pt1535" and len(feats) != 12:
            errors.append(f"baseV3E11_cent12_pt1535 feature count is not 12: {len(feats)}")
        elif model_id == "centAsFeatBase3x3_pt15to35" and len(feats) != 14:
            errors.append(f"centAsFeatBase3x3_pt15to35 feature count is not 14: {len(feats)}")
        elif model_id.startswith("base14_plus1_") and len(feats) != 15:
            errors.append(f"{model_id} plus-one feature count is not 15: {len(feats)}")
        elif model_id.startswith("base14_ladder") and not (15 <= len(feats) <= 31):
            errors.append(f"{model_id} cumulative feature count outside 15..31: {len(feats)}")
        elif model_id == "globalEtCent1535_bdt_noIso" and len(feats) != 32:
            errors.append(f"globalEtCent1535_bdt_noIso feature count is not 32: {len(feats)}")
        if errors:
            break
    ok_message = (
        "[OK] PPG12-exact preflight: corrected-baseline shower ladder "
        f"products={dict(products)}, selected={len(models)}, no event/cross-section training weights"
    )
elif campaign == "corrected-baseline-binned14":
    if planned_pt not in allowed_pt_grids:
        errors.append(f"unexpected pt_bins: {planned.get('pt_bins')}")
    n_pt_bins = max(0, len(planned_pt) - 1)
    n_coarse = len(expected_coarse)
    n_fine = len(expected_fine)
    expected_counts = {
        "base14_perEt": n_pt_bins,
        "base14_perCent3": n_coarse,
        "base14_perCent7": n_fine,
        "base14_perEtCent3": n_pt_bins * n_coarse,
        "base14_perEtCent7": n_pt_bins * n_fine,
    }
    unknown = sorted(set(products) - set(expected_counts))
    if unknown:
        errors.append(f"unexpected corrected-baseline binned14 products: {unknown}")
    if not products:
        errors.append("no corrected-baseline binned14 products planned")
    for product, count in sorted(products.items()):
        expected = expected_counts.get(product)
        if expected is not None and count > expected:
            errors.append(f"too many models for {product}: {count}, expected <= {expected}")
    if planned.get("coarse_cent_bins") != expected_coarse:
        errors.append(f"unexpected coarse_cent_bins: {planned.get('coarse_cent_bins')}")
    if planned.get("fine_cent_bins") != expected_fine:
        errors.append(f"unexpected fine_cent_bins: {planned.get('fine_cent_bins')}")
    bad_feature_counts = [
        f"{model.get('model_id')}:{len(model.get('features') or [])}"
        for model in models
        if len(model.get("features") or []) != 14
    ]
    if bad_feature_counts:
        errors.append("corrected-baseline binned14 feature count mismatch: " + ", ".join(bad_feature_counts[:8]))
    ok_message = (
        "[OK] PPG12-exact preflight: corrected-baseline 14-feature binned products "
        f"{dict(products)}, pT grid={planned_pt}, no event/cross-section training weights"
    )
elif campaign == "corrected-baseline-iso14":
    allowed = {
        "base14_eisoR30_pt1535": 15,
        "base14_eisoR40_pt1535": 15,
        "base14_eisoR30R40_pt1535": 16,
    }
    unknown = sorted(set(products) - set(allowed))
    if unknown:
        errors.append(f"unexpected corrected-baseline iso14 products: {unknown}")
    if not products:
        errors.append("no corrected-baseline iso14 products planned")
    for model in models:
        model_id = model.get("model_id", "")
        feats = model.get("features") or []
        expected = allowed.get(model_id)
        if expected is None:
            errors.append(f"unexpected iso14 model_id: {model_id}")
        elif len(feats) != expected:
            errors.append(f"{model_id} feature count is not {expected}: {len(feats)}")
        if model.get("pt_range") != [15.0, 35.0] or model.get("cent_range") is not None:
            errors.append(f"unexpected iso14 phase space: {model_id} pt={model.get('pt_range')} cent={model.get('cent_range')}")
        if not model.get("diagnostic_only"):
            errors.append(f"iso14 model not marked diagnostic_only: {model_id}")
        if errors:
            break
    ok_message = (
        "[OK] PPG12-exact preflight: corrected-baseline base14 raw-eiso diagnostics "
        f"{dict(products)}, no event/cross-section training weights"
    )
elif campaign == "expanded-tight":
    allowed_single_products = (
        {"centAsFeatBase3x3_pt15to35": 1},
        {"centAsFeatBase3x3_pt15to40": 1},
        {"centAsFeatBase3x3_pt8to40": 1},
        {"centAsFeatBase3x3_pt5to40": 1},
    )
    allowed_single_ids = {"centAsFeatBase3x3_pt15to35", "centAsFeatBase3x3_pt15to40", "centAsFeatBase3x3_pt8to40", "centAsFeatBase3x3_pt5to40"}
    if dict(products) not in allowed_single_products:
        errors.append(f"unexpected products/counts for expanded-tight baseline: {dict(products)}")
    if len(models) != 1 or planned.get("expected_model_count") != 1:
        errors.append(
            f"expected exactly 1 centAsFeatBase3x3 model, found models={len(models)} "
            f"expected_model_count={planned.get('expected_model_count')}"
        )
    if models:
        model = models[0]
        feats = model.get("features") or []
        required = {"cluster_Et", "centrality", "cluster_weta33_cogx", "cluster_wphi33_cogx"}
        missing = sorted(required - set(feats))
        if model.get("model_id") not in allowed_single_ids:
            errors.append(f"unexpected model_id: {model.get('model_id')}")
        if len(feats) != 14:
            errors.append(f"centAsFeatBase3x3 baseline feature count is not 14: {len(feats)}")
        if missing:
            errors.append(f"centAsFeatBase3x3 baseline missing required feature(s): {missing}")
    ok_message = "[OK] PPG12-exact preflight: centAsFeatBase3x3 14-input baseline/window, no event/cross-section training weights"
else:
    errors.append(f"unexpected campaign: {campaign}")
if errors:
    raise SystemExit("PPG12-exact preflight failed:\n  " + "\n  ".join(errors))
print(ok_message)
PY
  fi

  "$ML_PYTHON" - "$planned" "$shard_dir" "$group_size" <<'PY'
import json
import sys
from pathlib import Path
planned = Path(sys.argv[1])
shard_dir = Path(sys.argv[2])
group = int(sys.argv[3])
shard_dir.mkdir(parents=True, exist_ok=True)
for old in shard_dir.glob("specs_*.list"):
    old.unlink()
models = json.loads(planned.read_text())["models"]
ids = [m["model_id"] for m in models]
for i in range(0, len(ids), group):
    out = shard_dir / f"specs_{i // group:05d}.list"
    out.write_text("\n".join(ids[i:i + group]) + "\n")
print(len(ids))
PY
  local spec_count expected_spec_count shard_count
  spec_count="$(find "$shard_dir" -maxdepth 1 -type f -name 'specs_*.list' -exec cat {} + | wc -l | tr -d ' ')"
  expected_spec_count="$("$ML_PYTHON" - "$planned" <<'PY'
import json
import sys
print(json.load(open(sys.argv[1])).get("expected_model_count", 0))
PY
)"
  shard_count="$(find "$shard_dir" -maxdepth 1 -type f -name 'specs_*.list' | wc -l | tr -d ' ')"
  [[ "$spec_count" == "$expected_spec_count" ]] || die "Expanded campaign expected ${expected_spec_count} specs, planned ${spec_count}"

  local root_shard_dir="${sub_root}/cache_root_shards"
  local cache_partial_dir="${model_dir}/cache_shards"
  local cache_partials_manifest="${sub_root}/cache_partials.list"
  local cache_part_count="0"
  if [[ "$staged_cache" == "1" ]]; then
    mkdir -p "$root_shard_dir" "$cache_partial_dir"
    "$ML_PYTHON" - "$manifest" "$root_shard_dir" "$cache_partial_dir" "$cache_partials_manifest" "$cache_shards" <<'PY'
import math
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
root_shard_dir = Path(sys.argv[2])
cache_partial_dir = Path(sys.argv[3])
partials_manifest = Path(sys.argv[4])
requested = int(sys.argv[5])

roots = [line.strip() for line in manifest.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
if not roots:
    raise SystemExit(f"empty training manifest: {manifest}")
requested = max(1, min(requested, len(roots)))
for old in root_shard_dir.glob("roots_*.list"):
    old.unlink()
chunk = int(math.ceil(len(roots) / requested))
partials = []
for idx in range(requested):
    subset = roots[idx * chunk : (idx + 1) * chunk]
    if not subset:
        continue
    root_list = root_shard_dir / f"roots_{idx:05d}.list"
    root_list.write_text("\n".join(subset) + "\n")
    partials.append(cache_partial_dir / f"training_matrix_part_{idx:05d}.npz")
partials_manifest.write_text("\n".join(str(path) for path in partials) + "\n")
print(len(partials))
PY
    cache_part_count="$(wc -l < "$cache_partials_manifest" | tr -d ' ')"
    [[ "$cache_part_count" =~ ^[0-9]+$ && "$cache_part_count" -gt 0 ]] || die "staged cache requested but no partial cache shards were created"
  fi

  local env_prelude='
export USER="${USER:-$(id -u -n)}"
export LOGNAME="${LOGNAME:-$USER}"
export HOME="/sphenix/u/${LOGNAME}"
MYINSTALL="/sphenix/u/${USER}/thesisAnalysis/install"
MYINSTALL_AUAU="/sphenix/u/${USER}/thesisAnalysis_auau/install"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
if [[ -d "$MYINSTALL" ]]; then source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL" || true; fi
if [[ -d "$MYINSTALL_AUAU" ]]; then source /opt/sphenix/core/bin/setup_local.sh "$MYINSTALL_AUAU" || true; fi
set -u
ml_python="${ML_PYTHON}"
ml_python_prefix="$(cd "$(dirname "$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="$(readlink -f "$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
if [[ -n "$ml_python_real" ]]; then ml_python_real_prefix="$(cd "$(dirname "$ml_python_real")/.." && pwd -P 2>/dev/null || true)"; fi
ld_joined=""
for d in "$ml_python_prefix/lib" "$ml_python_prefix/lib64" "$ml_python_real_prefix/lib" "$ml_python_real_prefix/lib64"; do
  [[ -n "$d" && -d "$d" ]] || continue
  case ":$ld_joined:" in *":$d:"*) ;; *) ld_joined="${ld_joined:+$ld_joined:}$d" ;; esac
done
[[ -n "$ld_joined" ]] && export LD_LIBRARY_PATH="$ld_joined:${LD_LIBRARY_PATH:-}"
unset PYTHONHOME
'

  local cache_worker="${sub_root}/expanded_cache.sh"
cat > "$cache_worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export ML_PYTHON="${ML_PYTHON}"
export RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS="${spec_ids}"
export RJ_AUAU_BDT_CAMPAIGN="${campaign}"
export RJ_AUAU_BDT_WEIGHT_MODE="${weight_mode}"
export RJ_AUAU_BDT_LABEL_CONTRACT="${label_contract}"
export RJ_AUAU_BDT_TEST_SIZE="${bdt_test_size}"
export RJ_AUAU_BDT_SPLIT_MODE="${bdt_split_mode}"
export RJ_AUAU_BDT_MAX_DEPTH="${bdt_max_depth}"
export RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON="${event_quality_cut_json}"
export RJ_AUAU_BDT_EVENT_QUALITY_AUDIT_OUTPUT="${event_quality_audit_output}"
export RJ_AUAU_BDT_EVENT_QUALITY_ASSUME_FILTERED="${event_quality_assume_filtered}"
export RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES="${extra_base3x3_pt_ranges}"
${env_prelude}
extra_args=()
extra_args+=(--weight-mode "\${RJ_AUAU_BDT_WEIGHT_MODE}")
extra_args+=(--label-contract "\${RJ_AUAU_BDT_LABEL_CONTRACT}")
extra_args+=(--test-size "\${RJ_AUAU_BDT_TEST_SIZE}")
extra_args+=(--split-mode "\${RJ_AUAU_BDT_SPLIT_MODE}")
if [[ -n "\${RJ_AUAU_BDT_MAX_DEPTH:-}" ]]; then
  extra_args+=(--max-depth "\${RJ_AUAU_BDT_MAX_DEPTH}")
fi
if [[ -n "\${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON:-}" ]]; then
  extra_args+=(--event-quality-cut-json "\${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON}")
  extra_args+=(--event-quality-audit-output "\${RJ_AUAU_BDT_EVENT_QUALITY_AUDIT_OUTPUT}")
  extra_args+=(--require-global-event-key)
  if [[ "\${RJ_AUAU_BDT_EVENT_QUALITY_ASSUME_FILTERED:-0}" == "1" ]]; then
    extra_args+=(--event-quality-assume-filtered)
  fi
fi
if [[ -n "\${RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS:-}" ]]; then
  extra_args+=(--campaign-spec-ids "\${RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS}")
fi
if [[ -n "\${RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES:-}" ]]; then
  extra_args+=(--extra-cent-as-feat-base3x3-pt-ranges "\${RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES}")
fi
if [[ "\${RJ_AUAU_BDT_WEIGHT_MODE:-legacy}" == "ppg12-exact" ]]; then
  extra_args+=(--no-event-weight)
  extra_args+=(--ppg12-exact-expected-samples "${ppg12_expected_samples}")
  extra_args+=(--ppg12-exact-closure-dir "${ppg12_closure_dir}")
  extra_args+=(--ppg12-exact-closure-artifacts "${ppg12_closure_artifacts}")
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-sixpack" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-sixpack-noiso-ptcent7" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "global-and-etcent-binned-sixpack-noiso" ]]; then
  extra_args+=(--pt-bins "${etcent_pt_bins}")
  extra_args+=(--coarse-cent-bins "${etcent_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${etcent_fine_cent_bins}")
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etfine-centstudy" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "corrected-baseline-binned14" ]]; then
  extra_args+=(--pt-bins "${etfine_pt_bins}")
  extra_args+=(--coarse-cent-bins "${etfine_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${etfine_fine_cent_bins}")
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-eiso-cone-ablation" ]]; then
  extra_args+=(--pt-bins "${raw_eiso_pt_bins}")
  extra_args+=(--coarse-cent-bins "${raw_eiso_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${raw_eiso_fine_cent_bins}")
fi
"\$ml_python" "${TRAIN_SCRIPT}" --task tight --campaign "\${RJ_AUAU_BDT_CAMPAIGN}" \\
  --input "@${manifest}" --outdir "${model_dir}" \\
  --cache-file "${cache_file}" --cache-only \\
  --registry-output "${model_dir}/model_registry.cache.json" \\
  --n-jobs "${RJ_AUAU_BDT_XGB_N_JOBS:-1}" \\
  --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}" \\
  --minopt-majority-cap-ratio "${RJ_AUAU_BDT_MINOPT_MAJORITY_CAP_RATIO:-2.0}" \\
  "\${extra_args[@]}"
EOF
  chmod +x "$cache_worker"

  local cache_part_worker="${sub_root}/expanded_cache_part.sh"
cat > "$cache_part_worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
root_manifest="\${1:?root shard manifest}"
partial_cache="\${2:?partial cache output}"
registry="\${3:?registry output}"
export ML_PYTHON="${ML_PYTHON}"
export RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS="${spec_ids}"
export RJ_AUAU_BDT_CAMPAIGN="${campaign}"
export RJ_AUAU_BDT_WEIGHT_MODE="${weight_mode}"
export RJ_AUAU_BDT_LABEL_CONTRACT="${label_contract}"
export RJ_AUAU_BDT_TEST_SIZE="${bdt_test_size}"
export RJ_AUAU_BDT_SPLIT_MODE="${bdt_split_mode}"
export RJ_AUAU_BDT_MAX_DEPTH="${bdt_max_depth}"
export RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON="${event_quality_cut_json}"
export RJ_AUAU_BDT_EVENT_QUALITY_ASSUME_FILTERED="${event_quality_assume_filtered}"
export RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES="${extra_base3x3_pt_ranges}"
export MALLOC_ARENA_MAX="${RJ_AUAU_BDT_MALLOC_ARENA_MAX:-2}"
${env_prelude}
extra_args=()
extra_args+=(--weight-mode "\${RJ_AUAU_BDT_WEIGHT_MODE}")
extra_args+=(--label-contract "\${RJ_AUAU_BDT_LABEL_CONTRACT}")
extra_args+=(--test-size "\${RJ_AUAU_BDT_TEST_SIZE}")
extra_args+=(--split-mode "\${RJ_AUAU_BDT_SPLIT_MODE}")
if [[ -n "\${RJ_AUAU_BDT_MAX_DEPTH:-}" ]]; then
  extra_args+=(--max-depth "\${RJ_AUAU_BDT_MAX_DEPTH}")
fi
if [[ -n "\${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON:-}" ]]; then
  extra_args+=(--event-quality-cut-json "\${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON}")
  extra_args+=(--event-quality-audit-output "\${partial_cache%.npz}.event_quality_filter_audit.json")
  extra_args+=(--require-global-event-key)
  if [[ "\${RJ_AUAU_BDT_EVENT_QUALITY_ASSUME_FILTERED:-0}" == "1" ]]; then
    extra_args+=(--event-quality-assume-filtered)
  fi
fi
if [[ -n "\${RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS:-}" ]]; then
  extra_args+=(--campaign-spec-ids "\${RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS}")
fi
if [[ -n "\${RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES:-}" ]]; then
  extra_args+=(--extra-cent-as-feat-base3x3-pt-ranges "\${RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES}")
fi
if [[ "\${RJ_AUAU_BDT_WEIGHT_MODE:-legacy}" == "ppg12-exact" ]]; then
  extra_args+=(--no-event-weight)
  extra_args+=(--ppg12-exact-expected-samples "${ppg12_expected_samples}")
  extra_args+=(--ppg12-exact-closure-dir "${ppg12_closure_dir}")
  extra_args+=(--ppg12-exact-closure-artifacts metadata-only)
  extra_args+=(--cache-only-skip-ppg12-exact-weights)
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-sixpack" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-sixpack-noiso-ptcent7" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "global-and-etcent-binned-sixpack-noiso" ]]; then
  extra_args+=(--pt-bins "${etcent_pt_bins}")
  extra_args+=(--coarse-cent-bins "${etcent_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${etcent_fine_cent_bins}")
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etfine-centstudy" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "corrected-baseline-binned14" ]]; then
  extra_args+=(--pt-bins "${etfine_pt_bins}")
  extra_args+=(--coarse-cent-bins "${etfine_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${etfine_fine_cent_bins}")
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-eiso-cone-ablation" ]]; then
  extra_args+=(--pt-bins "${raw_eiso_pt_bins}")
  extra_args+=(--coarse-cent-bins "${raw_eiso_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${raw_eiso_fine_cent_bins}")
fi
"\$ml_python" "${TRAIN_SCRIPT}" --task tight --campaign "\${RJ_AUAU_BDT_CAMPAIGN}" \\
  --input "@\${root_manifest}" --outdir "${model_dir}" \\
  --cache-file "\${partial_cache}" --cache-only \\
  --registry-output "\${registry}" \\
  --n-jobs "${RJ_AUAU_BDT_XGB_N_JOBS:-1}" \\
  --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}" \\
  --minopt-majority-cap-ratio "${RJ_AUAU_BDT_MINOPT_MAJORITY_CAP_RATIO:-2.0}" \\
  "\${extra_args[@]}"
EOF
  chmod +x "$cache_part_worker"

  local cache_reduce_worker="${sub_root}/expanded_cache_reduce.sh"
cat > "$cache_reduce_worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export ML_PYTHON="${ML_PYTHON}"
export MALLOC_ARENA_MAX="${RJ_AUAU_BDT_MALLOC_ARENA_MAX:-2}"
${env_prelude}
"\$ml_python" - "${TRAIN_SCRIPT}" "${cache_partials_manifest}" "${cache_file}" "${model_dir}" "${ppg12_closure_dir}" "${ppg12_expected_samples}" "${ppg12_closure_artifacts}" "${model_dir}/model_registry.cache.json" <<'PY'
import importlib.util
import csv
import json
import math
import os
import sys
from pathlib import Path
import zipfile

import numpy as np
from numpy.lib.format import write_array

trainer_path = Path(sys.argv[1])
partials_manifest = Path(sys.argv[2])
output_cache = Path(sys.argv[3])
model_dir = Path(sys.argv[4])
closure_dir = Path(sys.argv[5])
expected_samples = sys.argv[6]
closure_artifacts = sys.argv[7]
registry_output = Path(sys.argv[8])

spec = importlib.util.spec_from_file_location("train_auau_photon_bdt_runtime", trainer_path)
if spec is None or spec.loader is None:
    raise SystemExit(f"cannot load trainer module from {trainer_path}")
trainer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(trainer)

partials = [Path(line.strip()) for line in partials_manifest.read_text().splitlines() if line.strip()]
if not partials:
    raise SystemExit(f"empty staged cache partial manifest: {partials_manifest}")
missing = [str(path) for path in partials if not path.is_file() or path.stat().st_size <= 0]
if missing:
    raise SystemExit("missing staged cache partial(s):\n  " + "\n  ".join(missing[:20]))

partial_reports = []
columns_ref = None
row_counts = []
for path in partials:
    with np.load(path, allow_pickle=True) as data:
        columns = [str(item) for item in data["__columns__"].tolist()]
        if not columns:
            raise SystemExit(f"empty staged cache partial columns: {path}")
        row_count = int(len(data[columns[0]]))
        for col in columns:
            if int(len(data[col])) != row_count:
                raise SystemExit(f"staged cache partial has inconsistent column lengths: {path}:{col}")
    if columns_ref is None:
        columns_ref = sorted(columns)
    elif sorted(columns) != columns_ref:
        raise SystemExit(
            "staged cache partial schema mismatch:\n"
            f"  reference={columns_ref}\n"
            f"  path={path}\n"
            f"  columns={sorted(columns)}"
        )
    row_counts.append(row_count)
    partial_reports.append({"path": str(path), "rows": row_count, "columns": sorted(columns)})

if columns_ref is None:
    raise SystemExit(f"empty staged cache partial manifest: {partials_manifest}")

weight_col = getattr(trainer, "PPG12_EXACT_WEIGHT_COLUMN")
required_reduce_columns = ["is_signal", "cluster_Et", "cluster_Eta", "source_sample"]
missing_reduce = [col for col in required_reduce_columns if col not in columns_ref]
if missing_reduce:
    raise SystemExit("staged cache partials missing reducer columns: " + ", ".join(missing_reduce))


def concat_column(name):
    pieces = []
    for path in partials:
        with np.load(path, allow_pickle=True) as data:
            pieces.append(np.asarray(data[name]))
    if not pieces:
        return np.asarray([])
    return np.concatenate(pieces)


labels = concat_column("is_signal").astype("int32", copy=False)
cluster_et = concat_column("cluster_Et").astype("float64", copy=False)
cluster_eta = concat_column("cluster_Eta").astype("float64", copy=False)
source_sample = np.asarray([str(item) for item in concat_column("source_sample")], dtype=object)
if "centrality" in columns_ref:
    centrality = concat_column("centrality").astype("float64", copy=False)
else:
    centrality = np.full(len(labels), np.nan, dtype="float64")
n_rows = int(len(labels))


def validate_samples_numpy(labels_arr, sample_arr, expected_text):
    expected = trainer.parse_expected_samples(expected_text)
    samples = sorted(set(map(str, sample_arr.tolist())))
    missing_samples = sorted(set(expected) - set(samples))
    unexpected = sorted(set(samples) - set(expected))
    if missing_samples or unexpected:
        raise SystemExit(
            "PPG12-exact sample set mismatch: "
            f"missing={missing_samples or []} unexpected={unexpected or []} observed={samples}"
        )
    inventory = []
    mixed_label_counts = {}
    for sample in expected:
        mask = sample_arr == sample
        n_signal = int(np.sum(mask & (labels_arr == 1)))
        n_background = int(np.sum(mask & (labels_arr == 0)))
        if "Photon" in sample and n_background:
            mixed_label_counts[sample] = n_background
        elif "Jet" in sample and n_signal:
            mixed_label_counts[sample] = n_signal
        inventory.append(
            {
                "source_sample": sample,
                "n_rows": int(mask.sum()),
                "n_signal": n_signal,
                "n_background": n_background,
            }
        )
    if mixed_label_counts:
        print(
            "[WARN] PPG12-exact source_sample/truth-label mixture observed; "
            "treating source_sample as provenance and truth label as the BDT class: "
            f"{mixed_label_counts}",
            flush=True,
        )
    return {
        "expected_samples": list(expected),
        "observed_samples": samples,
        "inventory": inventory,
        "mixed_label_counts": mixed_label_counts,
        "source_sample_semantics": "provenance",
        "truth_label_semantics": "per-candidate BDT class",
        "mixed_labels_are_fatal": False,
    }, expected


def compute_ppg12_exact_weights_numpy(labels_arr, et_arr, eta_arr):
    weights = np.ones(len(labels_arr), dtype="float64")
    report = {
        "weight_mode": "ppg12-exact",
        "event_weight_used": False,
        "vertex_reweight": False,
        "centrality_event_weight": False,
        "cross_section_weight_used_for_training": False,
        "weights_computed_before_binning": True,
        "eta_range": list(trainer.PPG12_EXACT_ETA_RANGE),
        "eta_bins": int(trainer.PPG12_EXACT_N_BINS),
        "et_bins": int(trainer.PPG12_EXACT_N_BINS),
        "et_weight_cap": float(trainer.PPG12_EXACT_ET_WEIGHT_CAP),
    }
    class_counts = {}
    class_weight_factors = {}
    n_total = 0
    for cls in (0, 1):
        n_cls = int((labels_arr == cls).sum())
        class_counts[str(cls)] = n_cls
        n_total += n_cls
    if class_counts["0"] <= 0 or class_counts["1"] <= 0:
        raise SystemExit(f"PPG12-exact weights need both classes; observed counts={class_counts}")
    for cls in (0, 1):
        factor = float(n_total) / (2.0 * float(class_counts[str(cls)]))
        class_weight_factors[str(cls)] = factor
        weights[labels_arr == cls] *= factor

    eta_reports = {}
    et_reports = {}
    for cls in (0, 1):
        mask = labels_arr == cls
        eta_w, eta_report = trainer.ppg12_exact_inverse_pdf_weights(
            eta_arr[mask],
            n_bins=int(trainer.PPG12_EXACT_N_BINS),
            fixed_range=trainer.PPG12_EXACT_ETA_RANGE,
            weight_cap=None,
        )
        weights[mask] *= eta_w
        eta_reports[str(cls)] = eta_report
        et_w, et_report = trainer.ppg12_exact_inverse_pdf_weights(
            et_arr[mask],
            n_bins=int(trainer.PPG12_EXACT_N_BINS),
            fixed_range=None,
            weight_cap=float(trainer.PPG12_EXACT_ET_WEIGHT_CAP),
        )
        weights[mask] *= et_w
        et_reports[str(cls)] = et_report

    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if not finite_positive.all():
        bad = int((~finite_positive).sum())
        raise SystemExit(f"PPG12-exact weights produced {bad} non-finite/non-positive rows")
    report["class_counts"] = class_counts
    report["class_weight_factors"] = class_weight_factors
    report["eta_reweight"] = eta_reports
    report["et_reweight"] = et_reports
    report["sum_weight_class0"] = float(weights[labels_arr == 0].sum())
    report["sum_weight_class1"] = float(weights[labels_arr == 1].sum())
    report["min_weight"] = float(np.min(weights)) if len(weights) else math.nan
    report["max_weight"] = float(np.max(weights)) if len(weights) else math.nan
    report["mean_weight"] = float(np.mean(weights)) if len(weights) else math.nan
    return weights, report


def write_array_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_inventory_artifacts(labels_arr, sample_arr, et_arr, eta_arr, cent_arr, weights_arr, outdir):
    rows = []
    for sample in sorted(set(map(str, sample_arr.tolist()))):
        sample_mask = sample_arr == sample
        for cls in (0, 1):
            mask = sample_mask & (labels_arr == cls)
            rows.append(
                {
                    "source_sample": sample,
                    "class": cls,
                    "class_name": "signal" if cls == 1 else "background",
                    "n_rows": int(mask.sum()),
                    "sum_ppg12_exact_weight": float(weights_arr[mask].sum()) if mask.any() else 0.0,
                    "mean_cluster_Et": float(np.nanmean(et_arr[mask])) if mask.any() else math.nan,
                    "mean_cluster_Eta": float(np.nanmean(eta_arr[mask])) if mask.any() else math.nan,
                }
            )
    inventory_csv = outdir / "ppg12_exact_sample_inventory.csv"
    write_array_csv(
        inventory_csv,
        rows,
        ["source_sample", "class", "class_name", "n_rows", "sum_ppg12_exact_weight", "mean_cluster_Et", "mean_cluster_Eta"],
    )
    et_edges = np.asarray([15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0])
    eta_edges = np.linspace(trainer.PPG12_EXACT_ETA_RANGE[0], trainer.PPG12_EXACT_ETA_RANGE[1], int(trainer.PPG12_EXACT_N_BINS) + 1)
    cent_bins = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
    binned_rows = []
    for sample in sorted(set(map(str, sample_arr.tolist()))):
        sample_mask = sample_arr == sample
        for cls in (0, 1):
            class_mask = sample_mask & (labels_arr == cls)
            for lo, hi in zip(et_edges[:-1], et_edges[1:]):
                mask = class_mask & (et_arr >= lo) & (et_arr < hi)
                binned_rows.append({"source_sample": sample, "class": cls, "axis": "cluster_Et", "bin_low": float(lo), "bin_high": float(hi), "n_rows": int(mask.sum()), "sum_ppg12_exact_weight": float(weights_arr[mask].sum()) if mask.any() else 0.0})
            for lo, hi in zip(eta_edges[:-1], eta_edges[1:]):
                mask = class_mask & (eta_arr >= lo) & (eta_arr < hi)
                binned_rows.append({"source_sample": sample, "class": cls, "axis": "cluster_Eta", "bin_low": float(lo), "bin_high": float(hi), "n_rows": int(mask.sum()), "sum_ppg12_exact_weight": float(weights_arr[mask].sum()) if mask.any() else 0.0})
            for lo, hi in cent_bins:
                mask = class_mask & (cent_arr >= lo) & (cent_arr < hi)
                binned_rows.append({"source_sample": sample, "class": cls, "axis": "centrality", "bin_low": float(lo), "bin_high": float(hi), "n_rows": int(mask.sum()), "sum_ppg12_exact_weight": float(weights_arr[mask].sum()) if mask.any() else 0.0})
    binned_csv = outdir / "ppg12_exact_sample_inventory_binned.csv"
    write_array_csv(
        binned_csv,
        binned_rows,
        ["source_sample", "class", "axis", "bin_low", "bin_high", "n_rows", "sum_ppg12_exact_weight"],
    )
    return {"inventory_csv": str(inventory_csv), "binned_inventory_csv": str(binned_csv)}


def write_npz_streaming(path, columns, weights_arr):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + ".tmp")
    if tmp.exists():
        tmp.unlink()
    with zipfile.ZipFile(tmp, mode="w", compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
        with archive.open("__columns__.npy", mode="w", force_zip64=True) as handle:
            write_array(handle, np.asarray(columns, dtype=object), allow_pickle=True)
        for col in columns:
            if col == weight_col:
                arr = weights_arr
            elif col == "is_signal":
                arr = labels
            elif col == "cluster_Et":
                arr = cluster_et
            elif col == "cluster_Eta":
                arr = cluster_eta
            elif col == "source_sample":
                arr = source_sample
            elif col == "centrality":
                arr = centrality
            else:
                arr = concat_column(col)
            with archive.open(f"{col}.npy", mode="w", force_zip64=True) as handle:
                write_array(handle, np.asarray(arr), allow_pickle=True)
            del arr
    os.replace(tmp, path)


sample_report, _expected_samples_tuple = validate_samples_numpy(labels, source_sample, expected_samples)
weights, weight_report = compute_ppg12_exact_weights_numpy(labels, cluster_et, cluster_eta)
closure_dir_path = Path(closure_dir)
closure_dir_path.mkdir(parents=True, exist_ok=True)
inventory_paths = {}
if closure_artifacts == "full":
    inventory_paths = write_inventory_artifacts(labels, source_sample, cluster_et, cluster_eta, centrality, weights, closure_dir_path)
metadata = {
    "schema": "AUAU_BDT_PPG12_EXACT_WEIGHT_CLOSURE_V1",
    "artifact_mode": closure_artifacts,
    "sample_validation": sample_report,
    "weighting": weight_report,
    "artifacts": inventory_paths,
}
metadata_path = closure_dir_path / "ppg12_exact_reweighting_metadata.json"
metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
closure = {**metadata, "metadata_json": str(metadata_path)}
output_columns = sorted([col for col in columns_ref if col != weight_col] + [weight_col])
write_npz_streaming(output_cache, output_columns, weights)
registry_output.parent.mkdir(parents=True, exist_ok=True)
registry_output.write_text(
    json.dumps(
        {
            "schema": "AUAU_BDT_STAGED_TRAINING_CACHE_REDUCE_V1",
            "status": "CACHE_READY",
            "partial_count": len(partials),
            "rows": n_rows,
            "columns": output_columns,
            "output_cache": str(output_cache),
            "partials_manifest": str(partials_manifest),
            "partials": partial_reports,
            "ppg12_exact_closure": closure,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
print(
    json.dumps(
        {
            "schema": "AUAU_BDT_STAGED_TRAINING_CACHE_REDUCE_V1",
            "status": "CACHE_READY",
            "partial_count": len(partials),
            "rows": n_rows,
            "output_cache": str(output_cache),
        },
        sort_keys=True,
    )
)
PY
EOF
  chmod +x "$cache_reduce_worker"

  local train_worker="${sub_root}/expanded_train_worker.sh"
  cat > "$train_worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
spec_list="\${1:?spec list}"
registry="\${2:?registry output}"
export ML_PYTHON="${ML_PYTHON}"
export RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS="${spec_ids}"
export RJ_AUAU_BDT_CAMPAIGN="${campaign}"
export RJ_AUAU_BDT_WEIGHT_MODE="${weight_mode}"
export RJ_AUAU_BDT_LABEL_CONTRACT="${label_contract}"
export RJ_AUAU_BDT_TEST_SIZE="${bdt_test_size}"
export RJ_AUAU_BDT_SPLIT_MODE="${bdt_split_mode}"
export RJ_AUAU_BDT_MAX_DEPTH="${bdt_max_depth}"
export RJ_AUAU_BDT_RANDOM_SEED="${bdt_random_seed}"
export RJ_AUAU_BDT_N_ESTIMATORS="${bdt_n_estimators}"
export RJ_AUAU_BDT_LEARNING_RATE="${bdt_learning_rate}"
export RJ_AUAU_BDT_SUBSAMPLE="${bdt_subsample}"
export RJ_AUAU_BDT_COLSAMPLE_BYTREE="${bdt_colsample_bytree}"
export RJ_AUAU_BDT_REG_ALPHA="${bdt_reg_alpha}"
export RJ_AUAU_BDT_REG_LAMBDA="${bdt_reg_lambda}"
export RJ_AUAU_BDT_MAX_BIN="${bdt_max_bin}"
export RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON="${event_quality_cut_json}"
export RJ_AUAU_BDT_EVENT_QUALITY_ASSUME_FILTERED="${event_quality_assume_filtered}"
export RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES="${extra_base3x3_pt_ranges}"
${env_prelude}
extra_args=()
extra_args+=(--weight-mode "\${RJ_AUAU_BDT_WEIGHT_MODE}")
extra_args+=(--label-contract "\${RJ_AUAU_BDT_LABEL_CONTRACT}")
extra_args+=(--test-size "\${RJ_AUAU_BDT_TEST_SIZE}")
extra_args+=(--split-mode "\${RJ_AUAU_BDT_SPLIT_MODE}")
if [[ -n "\${RJ_AUAU_BDT_MAX_DEPTH:-}" ]]; then
  extra_args+=(--max-depth "\${RJ_AUAU_BDT_MAX_DEPTH}")
fi
if [[ -n "\${RJ_AUAU_BDT_RANDOM_SEED:-}" ]]; then
  extra_args+=(--random-seed "\${RJ_AUAU_BDT_RANDOM_SEED}")
fi
if [[ -n "\${RJ_AUAU_BDT_N_ESTIMATORS:-}" ]]; then
  extra_args+=(--n-estimators "\${RJ_AUAU_BDT_N_ESTIMATORS}")
fi
if [[ -n "\${RJ_AUAU_BDT_LEARNING_RATE:-}" ]]; then
  extra_args+=(--learning-rate "\${RJ_AUAU_BDT_LEARNING_RATE}")
fi
if [[ -n "\${RJ_AUAU_BDT_SUBSAMPLE:-}" ]]; then
  extra_args+=(--subsample "\${RJ_AUAU_BDT_SUBSAMPLE}")
fi
if [[ -n "\${RJ_AUAU_BDT_COLSAMPLE_BYTREE:-}" ]]; then
  extra_args+=(--colsample-bytree "\${RJ_AUAU_BDT_COLSAMPLE_BYTREE}")
fi
if [[ -n "\${RJ_AUAU_BDT_REG_ALPHA:-}" ]]; then
  extra_args+=(--reg-alpha "\${RJ_AUAU_BDT_REG_ALPHA}")
fi
if [[ -n "\${RJ_AUAU_BDT_REG_LAMBDA:-}" ]]; then
  extra_args+=(--reg-lambda "\${RJ_AUAU_BDT_REG_LAMBDA}")
fi
if [[ -n "\${RJ_AUAU_BDT_MAX_BIN:-}" ]]; then
  extra_args+=(--max-bin "\${RJ_AUAU_BDT_MAX_BIN}")
fi
if [[ -n "\${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON:-}" ]]; then
  extra_args+=(--event-quality-cut-json "\${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON}")
  extra_args+=(--require-global-event-key)
  if [[ "\${RJ_AUAU_BDT_EVENT_QUALITY_ASSUME_FILTERED:-0}" == "1" ]]; then
    extra_args+=(--event-quality-assume-filtered)
  fi
fi
if [[ -n "\${RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS:-}" ]]; then
  extra_args+=(--campaign-spec-ids "\${RJ_AUAU_BDT_CAMPAIGN_SPEC_IDS}")
fi
if [[ -n "\${RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES:-}" ]]; then
  extra_args+=(--extra-cent-as-feat-base3x3-pt-ranges "\${RJ_AUAU_BDT_EXTRA_CENT_AS_FEAT_BASE3X3_PT_RANGES}")
fi
if [[ "\${RJ_AUAU_BDT_WEIGHT_MODE:-legacy}" == "ppg12-exact" ]]; then
  extra_args+=(--no-event-weight)
  extra_args+=(--ppg12-exact-expected-samples "${ppg12_expected_samples}")
  extra_args+=(--ppg12-exact-closure-dir "${ppg12_closure_dir}")
  extra_args+=(--ppg12-exact-closure-artifacts "${ppg12_closure_artifacts}")
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-sixpack" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-sixpack-noiso-ptcent7" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "global-and-etcent-binned-sixpack-noiso" ]]; then
  extra_args+=(--pt-bins "${etcent_pt_bins}")
  extra_args+=(--coarse-cent-bins "${etcent_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${etcent_fine_cent_bins}")
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etfine-centstudy" || "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "corrected-baseline-binned14" ]]; then
  extra_args+=(--pt-bins "${etfine_pt_bins}")
  extra_args+=(--coarse-cent-bins "${etfine_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${etfine_fine_cent_bins}")
fi
if [[ "\${RJ_AUAU_BDT_CAMPAIGN:-}" == "etcent-binned-eiso-cone-ablation" ]]; then
  extra_args+=(--pt-bins "${raw_eiso_pt_bins}")
  extra_args+=(--coarse-cent-bins "${raw_eiso_coarse_cent_bins}")
  extra_args+=(--fine-cent-bins "${raw_eiso_fine_cent_bins}")
fi
"\$ml_python" "${TRAIN_SCRIPT}" --task tight --campaign "\${RJ_AUAU_BDT_CAMPAIGN}" \\
  --outdir "${model_dir}" \\
  --cache-file "${cache_file}" \\
  --campaign-spec-list "\$spec_list" \\
  --registry-output "\$registry" \\
  --parallel-workers "${RJ_AUAU_BDT_TRAIN_PARALLEL_PER_JOB:-1}" \\
  --n-jobs "${RJ_AUAU_BDT_XGB_N_JOBS:-1}" \\
  --majority-cap-ratio "${RJ_AUAU_BDT_MAJORITY_CAP_RATIO:-4.0}" \\
  --minopt-majority-cap-ratio "${RJ_AUAU_BDT_MINOPT_MAJORITY_CAP_RATIO:-2.0}" \\
  "\${extra_args[@]}"
EOF
  chmod +x "$train_worker"

  local merge_worker="${sub_root}/expanded_merge.sh"
cat > "$merge_worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export ML_PYTHON="${ML_PYTHON}"
${env_prelude}
  "\$ml_python" - "${planned}" "${registry_dir}" "${model_dir}/model_registry.json" "${source}" "${manifest}" "${model_dir}" "${report_dir}/expanded_training_summary.txt" <<'PY'
import json
import sys
from pathlib import Path
try:
    import ROOT
except Exception as exc:
    ROOT = None
planned = json.loads(Path(sys.argv[1]).read_text())
registry_dir = Path(sys.argv[2])
out_registry = Path(sys.argv[3])
source = sys.argv[4]
manifest = sys.argv[5]
model_dir = Path(sys.argv[6])
summary = Path(sys.argv[7])
reports = {}
for path in sorted(registry_dir.glob("registry_*.json")):
    data = json.loads(path.read_text())
    for spec in data.get("models", []):
        if spec.get("report"):
            reports[spec["model_id"]] = spec["report"]
output_to_spec = {Path(spec["output_tmva"]).name: spec for spec in planned.get("models", [])}
skip_hints = {}
for log in sorted(registry_dir.parent.glob("train_*.out")):
    for line in log.read_text(errors="replace").splitlines():
        if "Skipping " not in line or "class counts below minimum" not in line:
            continue
        try:
            name = line.split("Skipping ", 1)[1].split(":", 1)[0].strip()
        except IndexError:
            continue
        spec = output_to_spec.get(name)
        if not spec:
            continue
        skip_hints[spec["model_id"]] = {
            "status": "skipped",
            "skip_reason": "class counts below minimum",
            "skip_log": str(log),
            "output_tmva": spec["output_tmva"],
            "output_xgb_json": spec.get("output_xgb_json"),
            "features": spec.get("features", []),
            "model_id": spec["model_id"],
            "product": spec.get("product"),
            "role": spec.get("role"),
            "pt_range": spec.get("pt_range"),
            "cent_range": spec.get("cent_range"),
            "minority_optimized": spec.get("minority_optimized", False),
            "n_signal": 0,
            "n_background": 0,
        }
missing = []
skipped = []
trained = []
bad = []
for spec in planned.get("models", []):
    report = reports.get(spec["model_id"])
    if not report:
        report = skip_hints.get(spec["model_id"])
    spec["report"] = report
    if not report:
        missing.append(spec["model_id"])
        continue
    report_status = report.get("status", "trained")
    if report_status == "skipped":
        skipped.append(spec["model_id"])
        continue
    if report_status != "trained":
        bad.append(f"{spec['model_id']}: unknown report status {report_status}")
        continue
    trained.append(spec["model_id"])
    tmva = Path(spec["output_tmva"])
    if not tmva.is_file() or tmva.stat().st_size <= 0:
        bad.append(f"{spec['model_id']}: missing {tmva}")
        continue
    if ROOT is not None:
        f = ROOT.TFile.Open(str(tmva))
        if not f or f.IsZombie():
            bad.append(f"{spec['model_id']}: unreadable {tmva}")
        if f:
            f.Close()
status = "READY" if not missing and not bad and not skipped else ("READY_WITH_SKIPS" if not missing and not bad else "CHECK")
planned["status"] = status
planned["trained_model_count"] = len(trained)
planned["skipped_model_count"] = len(skipped)
planned["missing_model_ids"] = missing
planned["skipped_model_ids"] = skipped
planned["bad_model_files"] = bad
out_registry.write_text(json.dumps(planned, indent=2, sort_keys=True) + "\n")
summary.write_text(
    "\n".join([
        "RECOILJETS_AUAU_TIGHT_BDT_EXPANDED_TRAINING_V1",
        f"status={status}",
        f"source={source}",
        f"training_manifest={manifest}",
        f"model_dir={model_dir}",
        f"registry={out_registry}",
        f"expected_model_count={planned.get('expected_model_count', len(planned.get('models', [])))}",
        f"trained_model_count={len(trained)}",
        f"skipped_model_count={len(skipped)}",
        f"missing_model_count={len(missing)}",
        f"bad_model_count={len(bad)}",
        "skipped_model_ids=" + ",".join(skipped),
        "missing_model_ids=" + ",".join(missing),
        "next_action=Run validateOnSimCondor with this expanded MODEL_DIR and inspect registry-ranked model QA.",
        "",
    ])
)
print(summary.read_text())
sys.exit(0 if status in ("READY", "READY_WITH_SKIPS") else 3)
PY
if command -v mail >/dev/null 2>&1; then
  status=\$(awk -F= '/^status=/ {print \$2; exit}' "${report_dir}/expanded_training_summary.txt")
  mail -s "[RecoilJets][auauTightBDT_trainExpandedFromExtractionCondor][\${status:-CHECK}]" "${NOTIFY_EMAILS}" < "${report_dir}/expanded_training_summary.txt" || true
fi
EOF
  chmod +x "$merge_worker"

  local cache_sub="${sub_root}/expanded_cache.sub"
  cat > "$cache_sub" <<EOF
universe = vanilla
executable = ${cache_worker}
output = ${sub_root}/cache.out
error = ${sub_root}/cache.err
log = ${sub_root}/cache.log
request_memory = ${cache_reqmem}
notification = Never
queue
EOF

  local cache_part_nodes="${sub_root}/cache_part_nodes.txt"
  : > "$cache_part_nodes"
  if [[ "$staged_cache" == "1" ]]; then
    local part_idx=0
    local root_list partial_cache partial_registry
    for root_list in "${root_shard_dir}"/roots_*.list; do
      [[ -s "$root_list" ]] || continue
      partial_cache="$(sed -n "$((part_idx + 1))p" "$cache_partials_manifest")"
      [[ -n "$partial_cache" ]] || die "missing partial cache path for staged cache shard ${part_idx}"
      local part_label; part_label="$(printf '%05d' "$part_idx")"
      partial_registry="${cache_partial_dir}/registry_part_${part_label}.json"
      local part_sub="${sub_root}/cache_part_${part_label}.sub"
      cat > "$part_sub" <<EOF
universe = vanilla
executable = ${cache_part_worker}
arguments = ${root_list} ${partial_cache} ${partial_registry}
output = ${sub_root}/cache_part_${part_label}.out
error = ${sub_root}/cache_part_${part_label}.err
log = ${sub_root}/cache_part_${part_label}.log
request_memory = ${cache_part_reqmem}
notification = Never
queue
EOF
      printf 'CACHE_PART_%s %s\n' "$part_label" "$part_sub" >> "$cache_part_nodes"
      part_idx=$((part_idx + 1))
    done
    [[ "$part_idx" == "$cache_part_count" ]] || die "staged cache shard count mismatch: roots=${part_idx} partials=${cache_part_count}"
  fi

  local cache_reduce_sub="${sub_root}/expanded_cache_reduce.sub"
  cat > "$cache_reduce_sub" <<EOF
universe = vanilla
executable = ${cache_reduce_worker}
output = ${sub_root}/cache_reduce.out
error = ${sub_root}/cache_reduce.err
log = ${sub_root}/cache_reduce.log
request_memory = ${cache_reduce_reqmem}
notification = Never
queue
EOF

  local train_nodes="${sub_root}/train_nodes.txt"
  : > "$train_nodes"
  local idx=0
  local spec_list
  for spec_list in "${shard_dir}"/specs_*.list; do
    [[ -s "$spec_list" ]] || continue
    idx=$((idx + 1))
    local label; label="$(printf '%05d' "$idx")"
    local reg="${registry_dir}/registry_${label}.json"
    local sub="${sub_root}/train_${label}.sub"
    cat > "$sub" <<EOF
universe = vanilla
executable = ${train_worker}
arguments = ${spec_list} ${reg}
output = ${sub_root}/train_${label}.out
error = ${sub_root}/train_${label}.err
log = ${sub_root}/train_${label}.log
request_memory = ${reqmem}
notification = Never
queue
EOF
    printf 'TRAIN_%s %s\n' "$label" "$sub" >> "$train_nodes"
  done

  local merge_sub="${sub_root}/expanded_merge.sub"
  cat > "$merge_sub" <<EOF
universe = scheduler
executable = ${merge_worker}
output = ${sub_root}/merge.out
error = ${sub_root}/merge.err
log = ${sub_root}/merge.log
notification = Never
queue
EOF

  local dag="${sub_root}/auau_tight_bdt_expanded_training.dag"
  {
    if [[ "$staged_cache" == "1" ]]; then
      echo "MAXJOBS cache_part ${cache_part_maxjobs}"
      echo
      while read -r node sub; do
        [[ -n "${node:-}" ]] || continue
        echo "JOB ${node} ${sub}"
        echo "CATEGORY ${node} cache_part"
        echo "RETRY ${node} 1"
      done < "$cache_part_nodes"
      echo "JOB CACHE_REDUCE ${cache_reduce_sub}"
      while read -r node sub; do
        [[ -n "${node:-}" ]] || continue
        echo "PARENT ${node} CHILD CACHE_REDUCE"
      done < "$cache_part_nodes"
    elif [[ "$reuse_existing_cache" != "1" ]]; then
      echo "JOB CACHE ${cache_sub}"
    fi
    while read -r node sub; do
      [[ -n "${node:-}" ]] || continue
      echo "JOB ${node} ${sub}"
      if [[ "$staged_cache" == "1" ]]; then
        echo "PARENT CACHE_REDUCE CHILD ${node}"
      elif [[ "$reuse_existing_cache" != "1" ]]; then
        echo "PARENT CACHE CHILD ${node}"
      fi
      echo "RETRY ${node} 1"
    done < "$train_nodes"
    echo "FINAL MERGE ${merge_sub}"
  } > "$dag"
  say "Expanded training DAG: $dag"
  say "source: $source"
  say "model dir: $model_dir"
  say "specs=${spec_count} shards=${shard_count} groupSize=${group_size} request_memory=${reqmem}"
  if [[ "$staged_cache" == "1" ]]; then
    say "staged cache: root_shards=${cache_part_count} part_memory=${cache_part_reqmem} reduce_memory=${cache_reduce_reqmem} maxjobs=${cache_part_maxjobs}"
  elif [[ "$reuse_existing_cache" == "1" ]]; then
    say "verified existing cache: ${cache_file} sha256=${expected_cache_sha256} (no CACHE node; label contract is applied once in TRAIN)"
  fi
  if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
    echo "RECOILJETS_AUAU_TIGHT_BDT_EXPANDED_DRYRUN_V1"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "planned_registry=${planned}"
    echo "dag=${dag}"
    echo "specs=${spec_count}"
    echo "shards=${shard_count}"
    echo "staged_cache=${staged_cache}"
    echo "reuse_existing_cache=${reuse_existing_cache}"
    echo "expected_cache_sha256=${expected_cache_sha256}"
    echo "cache_part_count=${cache_part_count}"
    return 0
  fi
  condor_submit_dag "$dag"
}

finalize_expanded_training() {
  local source=""
  local model_dir=""
  local sub_root=""
  local report_dir=""
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|modelDir=*|model_dir=*) model_dir="${tok#*=}" ;;
      SUB_ROOT=*|subRoot=*|sub_root=*) sub_root="${tok#*=}" ;;
      REPORT_DIR=*|reportDir=*|report_dir=*) report_dir="${tok#*=}" ;;
    esac
  done
  [[ -n "$source" ]] || die "finalizeExpandedTraining requires SOURCE=/path"
  [[ -n "$model_dir" ]] || die "finalizeExpandedTraining requires MODEL_DIR=/path"
  [[ -n "$sub_root" ]] || die "finalizeExpandedTraining requires SUB_ROOT=/path"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"
  [[ -d "$model_dir" ]] || die "MODEL_DIR is not a directory: $model_dir"
  [[ -d "$sub_root" ]] || die "SUB_ROOT is not a directory: $sub_root"
  [[ -n "$report_dir" ]] || report_dir="${source}/reports"
  guard_generated_path "expanded finalization report dir" "$report_dir"
  guard_generated_path "expanded finalization model dir" "$model_dir"
  guard_generated_path "expanded finalization submit root" "$sub_root"
  mkdir -p "$report_dir"
  local planned="${model_dir}/model_registry.planned.json"
  local registry_dir="${sub_root}/registries"
  local out_registry="${model_dir}/model_registry.json"
  local summary="${report_dir}/expanded_training_summary.txt"
  [[ -s "$planned" ]] || die "Missing planned registry: $planned"
  [[ -d "$registry_dir" ]] || die "Missing registry shard directory: $registry_dir"
  log_path_plan "finalizeExpandedTraining" \
    "source    : ${source}" \
    "model dir : ${model_dir}" \
    "sub root  : ${sub_root}" \
    "registry  : ${out_registry}" \
    "summary   : ${summary}"
  setup_ml_python_env
  "$ML_PYTHON" - "$planned" "$registry_dir" "$out_registry" "$source" "$model_dir" "$summary" <<'PY'
import json
import sys
from pathlib import Path

try:
    import ROOT
except Exception:
    ROOT = None

planned = json.loads(Path(sys.argv[1]).read_text())
registry_dir = Path(sys.argv[2])
out_registry = Path(sys.argv[3])
source = sys.argv[4]
model_dir = Path(sys.argv[5])
summary = Path(sys.argv[6])

reports = {}
for path in sorted(registry_dir.glob("registry_*.json")):
    data = json.loads(path.read_text())
    for spec in data.get("models", []):
        if spec.get("report"):
            reports[spec["model_id"]] = spec["report"]

output_to_spec = {Path(spec["output_tmva"]).name: spec for spec in planned.get("models", [])}
skip_hints = {}
for log in sorted(registry_dir.parent.glob("train_*.out")):
    for line in log.read_text(errors="replace").splitlines():
        if "Skipping " not in line or "class counts below minimum" not in line:
            continue
        try:
            name = line.split("Skipping ", 1)[1].split(":", 1)[0].strip()
        except IndexError:
            continue
        spec = output_to_spec.get(name)
        if not spec:
            continue
        skip_hints[spec["model_id"]] = {
            "status": "skipped",
            "skip_reason": "class counts below minimum",
            "skip_log": str(log),
            "output_tmva": spec["output_tmva"],
            "output_xgb_json": spec.get("output_xgb_json"),
            "features": spec.get("features", []),
            "model_id": spec["model_id"],
            "product": spec.get("product"),
            "role": spec.get("role"),
            "pt_range": spec.get("pt_range"),
            "cent_range": spec.get("cent_range"),
            "minority_optimized": spec.get("minority_optimized", False),
            "n_signal": 0,
            "n_background": 0,
        }

missing = []
skipped = []
trained = []
bad = []
for spec in planned.get("models", []):
    report = reports.get(spec["model_id"]) or skip_hints.get(spec["model_id"])
    spec["report"] = report
    if not report:
        missing.append(spec["model_id"])
        continue
    report_status = report.get("status", "trained")
    if report_status == "skipped":
        skipped.append(spec["model_id"])
        continue
    if report_status != "trained":
        bad.append(f"{spec['model_id']}: unknown report status {report_status}")
        continue
    trained.append(spec["model_id"])
    tmva = Path(spec["output_tmva"])
    if not tmva.is_file() or tmva.stat().st_size <= 0:
        bad.append(f"{spec['model_id']}: missing {tmva}")
        continue
    if ROOT is not None:
        f = ROOT.TFile.Open(str(tmva))
        if not f or f.IsZombie():
            bad.append(f"{spec['model_id']}: unreadable {tmva}")
        if f:
            f.Close()

status = "READY" if not missing and not bad and not skipped else ("READY_WITH_SKIPS" if not missing and not bad else "CHECK")
planned["status"] = status
planned["trained_model_count"] = len(trained)
planned["skipped_model_count"] = len(skipped)
planned["missing_model_ids"] = missing
planned["skipped_model_ids"] = skipped
planned["bad_model_files"] = bad
out_registry.write_text(json.dumps(planned, indent=2, sort_keys=True) + "\n")
summary.write_text(
    "\n".join([
        "RECOILJETS_AUAU_TIGHT_BDT_EXPANDED_TRAINING_V1",
        f"status={status}",
        f"source={source}",
        f"model_dir={model_dir}",
        f"registry={out_registry}",
        f"expected_model_count={planned.get('expected_model_count', len(planned.get('models', [])))}",
        f"trained_model_count={len(trained)}",
        f"skipped_model_count={len(skipped)}",
        f"missing_model_count={len(missing)}",
        f"bad_model_count={len(bad)}",
        "skipped_model_ids=" + ",".join(skipped),
        "missing_model_ids=" + ",".join(missing),
        "next_action=Run validateOnSimCondor with this expanded MODEL_DIR and inspect registry-ranked model QA.",
        "",
    ])
)
print(summary.read_text())
sys.exit(0 if status in ("READY", "READY_WITH_SKIPS") else 3)
PY
  local status
  status=$(awk -F= '/^status=/ {print $2; exit}' "$summary" 2>/dev/null || true)
  send_summary_email "[RecoilJets][auauTightBDT_finalizeExpandedTraining][${status:-CHECK}]" "$summary"
}

apply_check() {
  local model_dir=""
  for tok in "$@"; do
    case "$tok" in
      MODEL_DIR=*) model_dir="${tok#MODEL_DIR=}" ;;
    esac
  done
  [[ -n "$model_dir" ]] || die "applyCheck requires MODEL_DIR=/path"
  local required=(
    auau_tight_bdt_centINDcontrol_allCent_tmva.root
    auau_tight_bdt_centAsFeat_allCent_tmva.root
    auau_tight_bdt_centDepBDTs_cent_000_020_tmva.root
    auau_tight_bdt_centDepBDTs_cent_020_050_tmva.root
    auau_tight_bdt_centDepBDTs_cent_050_080_tmva.root
  )
  for f in "${required[@]}"; do
    [[ -s "${model_dir}/${f}" ]] || die "Missing model product: ${model_dir}/${f}"
  done
  setup_ml_python_env
  "$ML_PYTHON" - "$model_dir" <<'PY'
import sys
from pathlib import Path
try:
    import ROOT
except Exception as exc:
    raise SystemExit(f"PyROOT import failed: {exc}")
model_dir = Path(sys.argv[1])
bad = []
for path in model_dir.glob("*_tmva.root"):
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        bad.append(str(path))
    if f:
        f.Close()
if bad:
    raise SystemExit("Unreadable TMVA ROOT files:\n  " + "\n  ".join(bad))
print(f"[OK] applyCheck opened {len(list(model_dir.glob('*_tmva.root')))} TMVA ROOT files")
PY
  say "applyCheck PASS: $model_dir"
}

validate_on_sim() {
  local source=""
  local model_dir=""
  local model_registry=""
  local outdir=""
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|modelDir=*|model_dir=*) model_dir="${tok#*=}" ;;
      MODEL_REGISTRY=*|REGISTRY=*|modelRegistry=*|model_registry=*) model_registry="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
    esac
  done
  [[ -n "$source" ]] || die "validateOnSim requires SOURCE=/path/to/auauTightBDT extraction"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"
  [[ -n "$model_dir" ]] || die "validateOnSim requires MODEL_DIR=/path/to/tight models"
  [[ -d "$model_dir" ]] || die "MODEL_DIR is not a directory: $model_dir"
  [[ -z "$model_registry" || -s "$model_registry" ]] || die "MODEL_REGISTRY is not a readable file: $model_registry"
  [[ -s "$VALIDATE_SCRIPT" ]] || die "Missing validation script: $VALIDATE_SCRIPT"

  setup_ml_python_env
  local -a args
  args=( "$VALIDATE_SCRIPT" --source "$source" --model-dir "$model_dir" )
  if [[ -n "$model_registry" ]]; then
    args+=( --model-registry "$model_registry" )
  fi
  if [[ -n "$outdir" ]]; then
    args+=( --outdir "$outdir" )
  fi
  if [[ -n "${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON:-}" ]]; then
    args+=( --event-quality-cut-json "${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON}" )
    args+=( --event-quality-audit-output "${outdir:-${source}/reports}/event_quality_filter_audit.json" )
  fi

  say "Validating tight-BDT models on embedded-sim extraction trees"
  say "  source    : $source"
  say "  model dir : $model_dir"
  [[ -n "$model_registry" ]] && say "  registry  : $model_registry"
  [[ -n "${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON:-}" ]] && say "  event cut : ${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON}"
  local rc=0
  "$ML_PYTHON" "${args[@]}" || rc=$?

  local summary=""
  if [[ -n "$outdir" && -s "${outdir}/validation_summary.txt" ]]; then
    summary="${outdir}/validation_summary.txt"
  else
    summary="$(find "${source}/reports" -path '*/model_validation_*/validation_summary.txt' -type f 2>/dev/null | sort -V | tail -n 1 || true)"
  fi
  if [[ -n "$summary" && -s "$summary" ]]; then
    local status
    status="$(awk -F= '/^status=/ {print $2; exit}' "$summary")"
    status="${status:-CHECK}"
    send_summary_email "[RecoilJets][auauTightBDT_validateOnSim][${status}]" "$summary"
    say "validation summary: $summary"
  else
    warn "No validation summary found after validateOnSim"
  fi
  return "$rc"
}

validate_on_sim_condor() {
  local source=""
  local model_dir=""
  local model_registry=""
  local outdir=""
  local group_size="${RJ_AUAU_TIGHT_BDT_VALIDATE_GROUP_SIZE:-100}"
  local total_score_max="${RJ_AUAU_TIGHT_BDT_VALIDATE_TOTAL_SCORE_MAX_ROWS:-400000}"
  local reqmem="${RJ_AUAU_TIGHT_BDT_VALIDATE_REQUEST_MEMORY:-2500MB}"
  local merge_universe="${RJ_AUAU_TIGHT_BDT_VALIDATE_MERGE_UNIVERSE:-scheduler}"
  local merge_reqmem="${RJ_AUAU_TIGHT_BDT_VALIDATE_MERGE_REQUEST_MEMORY:-}"
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|modelDir=*|model_dir=*) model_dir="${tok#*=}" ;;
      MODEL_REGISTRY=*|REGISTRY=*|modelRegistry=*|model_registry=*) model_registry="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
      groupSize=*) group_size="${tok#groupSize=}" ;;
      SCORE_MAX_ROWS=*|scoreMaxRows=*) total_score_max="${tok#*=}" ;;
    esac
  done
  while (($#)); do
    case "$1" in
      groupSize) group_size="${2:?missing value after groupSize}"; shift 2 ;;
      scoreMaxRows|totalScoreMaxRows) total_score_max="${2:?missing value after $1}"; shift 2 ;;
      *) shift ;;
    esac
  done
  [[ -n "$source" ]] || die "validateOnSimCondor requires SOURCE=/path/to/auauTightBDT extraction"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"
  [[ -n "$model_dir" ]] || die "validateOnSimCondor requires MODEL_DIR=/path/to/tight models"
  [[ -d "$model_dir" ]] || die "MODEL_DIR is not a directory: $model_dir"
  [[ -z "$model_registry" || -s "$model_registry" ]] || die "MODEL_REGISTRY is not a readable file: $model_registry"
  [[ -s "$VALIDATE_SCRIPT" ]] || die "Missing validation script: $VALIDATE_SCRIPT"
  [[ "$group_size" =~ ^[0-9]+$ && "$group_size" -gt 0 ]] || die "groupSize must be a positive integer"
  case "$merge_universe" in
    scheduler|vanilla) ;;
    *) die "RJ_AUAU_TIGHT_BDT_VALIDATE_MERGE_UNIVERSE must be scheduler or vanilla, got: ${merge_universe}" ;;
  esac

  local stamp="${RJ_AUAU_TIGHT_BDT_VALIDATE_STAMP:-$(ts)}"
  local report_root="${outdir:-${source}/reports/model_validation_condor_${stamp}}"
  local sub_root="${RJ_REPO_BASE}/condor_sub/auauTightBDTValidate_${stamp}"
  local shard_dir="${sub_root}/shards"
  local cache_dir="${report_root}/score_caches"
  local event_quality_cut_json="${RJ_AUAU_BDT_EVENT_QUALITY_CUT_JSON:-}"
  guard_generated_path "validation report root" "$report_root"
  guard_generated_path "validation submit root" "$sub_root"
  log_path_plan "validateOnSimCondor" \
    "source    : ${source}" \
    "model dir : ${model_dir}" \
    "registry  : ${model_registry:-<default model_registry.json>}" \
    "report    : ${report_root}" \
    "submit    : ${sub_root}" \
    "shards    : ${shard_dir}" \
    "caches    : ${cache_dir}" \
    "event cut : ${event_quality_cut_json:-<disabled>}" \
    "groupSize : ${group_size}" \
    "requestMem: ${reqmem}" \
    "mergeUniv : ${merge_universe}" \
    "mergeMem  : ${merge_reqmem:-<default>}" \
    "scoreMax  : ${total_score_max}"
  mkdir -p "$report_root" "$sub_root" "$shard_dir" "$cache_dir"

  local root_manifest="${source}/manifests/training_roots.list"
  local search_root="$source"
  [[ -d "${source}/extraction" ]] && search_root="${source}/extraction"
  make_root_manifest "$search_root" "$root_manifest"
  local nroots
  nroots="$(wc -l < "$root_manifest" | tr -d ' ')"
  [[ "$nroots" != "0" ]] || die "No ROOT files available for Condor validation"

  local split_prefix="${shard_dir}/roots_"
  rm -f "${split_prefix}"*
  split -l "$group_size" -d -a 5 "$root_manifest" "$split_prefix"
  local shard_count
  shard_count="$(find "$shard_dir" -maxdepth 1 -type f -name 'roots_*' | wc -l | tr -d ' ')"
  [[ "$shard_count" != "0" ]] || die "Failed to split validation manifest"
  local score_max_per_shard=0
  if [[ "$total_score_max" =~ ^[0-9]+$ && "$total_score_max" -gt 0 ]]; then
    score_max_per_shard=$(( (total_score_max + shard_count - 1) / shard_count ))
  fi

  local args_file="${sub_root}/validate_args.txt"
  : > "$args_file"
  local idx=0
  local shard
  for shard in "${split_prefix}"*; do
    [[ -s "$shard" ]] || continue
    idx=$((idx + 1))
    local shard_out="${report_root}/shards/shard_$(printf '%05d' "$idx")"
    local cache="${cache_dir}/score_cache_$(printf '%05d' "$idx").npz"
    mkdir -p "$shard_out"
    printf '%s %s %s %s\n' "$shard" "$shard_out" "$cache" "$score_max_per_shard" >> "$args_file"
  done

  local worker="${sub_root}/validate_worker.sh"
  cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
echo "[auauTightBDT] worker env setup start" >&2
export USER="\${USER:-\$(id -u -n)}"
export LOGNAME="\${LOGNAME:-\$USER}"
export HOME="/sphenix/u/\${LOGNAME}"
MYINSTALL="/sphenix/u/\${USER}/thesisAnalysis/install"
MYINSTALL_AUAU="/sphenix/u/\${USER}/thesisAnalysis_auau/install"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
echo "[auauTightBDT] worker sphenix_setup rc=\$?" >&2
if [[ -d "\$MYINSTALL" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL" || true
  echo "[auauTightBDT] worker setup_local rc=\$?" >&2
fi
if [[ -d "\$MYINSTALL_AUAU" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL_AUAU" || true
  echo "[auauTightBDT] worker setup_local_auau rc=\$?" >&2
fi
set -u
echo "[auauTightBDT] worker env setup done" >&2
manifest="\${1:?manifest}"
outdir="\${2:?outdir}"
cache="\${3:?cache}"
score_max="\${4:?score_max}"
model_registry="${model_registry}"
ml_python="${ML_PYTHON}"
ml_python_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="\$(readlink -f "\$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
if [[ -n "\$ml_python_real" ]]; then
  ml_python_real_prefix="\$(cd "\$(dirname "\$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
fi
ld_joined=""
for d in "\$ml_python_prefix/lib" "\$ml_python_prefix/lib64" "\$ml_python_real_prefix/lib" "\$ml_python_real_prefix/lib64"; do
  [[ -n "\$d" && -d "\$d" ]] || continue
  case ":\$ld_joined:" in *":\$d:"*) ;; *) ld_joined="\${ld_joined:+\$ld_joined:}\$d" ;; esac
done
[[ -n "\$ld_joined" ]] && export LD_LIBRARY_PATH="\$ld_joined:\${LD_LIBRARY_PATH:-}"
unset PYTHONHOME
model_registry_arg=()
if [[ -n "\$model_registry" ]]; then
  model_registry_arg=(--model-registry "\$model_registry")
fi
event_quality_args=()
if [[ -n "${event_quality_cut_json}" ]]; then
  event_quality_args=(--event-quality-cut-json "${event_quality_cut_json}" --event-quality-audit-output "\${outdir}/event_quality_filter_audit.json")
fi
"\$ml_python" "${VALIDATE_SCRIPT}" \\
  --source "${source}" \\
  --model-dir "${model_dir}" \\
  "\${model_registry_arg[@]}" \\
  "\${event_quality_args[@]}" \\
  --manifest "\$manifest" \\
  --outdir "\$outdir" \\
  --write-score-cache "\$cache" \\
  --score-max-rows "\$score_max" \\
  --no-plots
EOF
  chmod +x "$worker"

  local shard_sub_dir="${sub_root}/shard_subs"
  mkdir -p "$shard_sub_dir"
  rm -f "${shard_sub_dir}"/validate_shard_*.sub
  local shard_nodes_file="${sub_root}/validate_nodes.txt"
  : > "$shard_nodes_file"
  local shard_index=0
  local shard_manifest shard_out shard_cache shard_scoremax
  while read -r shard_manifest shard_out shard_cache shard_scoremax; do
    [[ -n "${shard_manifest:-}" ]] || continue
    shard_index=$((shard_index + 1))
    local shard_label
    shard_label="$(printf '%05d' "$shard_index")"
    local shard_sub="${shard_sub_dir}/validate_shard_${shard_label}.sub"
    cat > "$shard_sub" <<EOF
universe = vanilla
executable = ${worker}
arguments = ${shard_manifest} ${shard_out} ${shard_cache} ${shard_scoremax}
output = ${sub_root}/validate_${shard_label}.out
error = ${sub_root}/validate_${shard_label}.err
log = ${sub_root}/validate_${shard_label}.log
request_memory = ${reqmem}
notification = Never
queue
EOF
    printf 'VALIDATE_%s %s\n' "$shard_label" "$shard_sub" >> "$shard_nodes_file"
  done < "$args_file"
  [[ "$shard_index" -eq "$idx" ]] || die "Internal validation shard mismatch: args=${idx} submit_files=${shard_index}"

  local merge="${sub_root}/validate_merge.sh"
  cat > "$merge" <<EOF
#!/usr/bin/env bash
set -euo pipefail
echo "[auauTightBDT] merge env setup start" >&2
export USER="\${USER:-\$(id -u -n)}"
export LOGNAME="\${LOGNAME:-\$USER}"
export HOME="/sphenix/u/\${LOGNAME}"
MYINSTALL="/sphenix/u/\${USER}/thesisAnalysis/install"
MYINSTALL_AUAU="/sphenix/u/\${USER}/thesisAnalysis_auau/install"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n
echo "[auauTightBDT] merge sphenix_setup rc=\$?" >&2
if [[ -d "\$MYINSTALL" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL" || true
  echo "[auauTightBDT] merge setup_local rc=\$?" >&2
fi
if [[ -d "\$MYINSTALL_AUAU" ]]; then
  source /opt/sphenix/core/bin/setup_local.sh "\$MYINSTALL_AUAU" || true
  echo "[auauTightBDT] merge setup_local_auau rc=\$?" >&2
fi
set -u
echo "[auauTightBDT] merge env setup done" >&2
model_registry="${model_registry}"
ml_python="${ML_PYTHON}"
ml_python_prefix="\$(cd "\$(dirname "\$ml_python")/.." && pwd -P 2>/dev/null || true)"
ml_python_real="\$(readlink -f "\$ml_python" 2>/dev/null || true)"
ml_python_real_prefix=""
if [[ -n "\$ml_python_real" ]]; then
  ml_python_real_prefix="\$(cd "\$(dirname "\$ml_python_real")/.." && pwd -P 2>/dev/null || true)"
fi
ld_joined=""
for d in "\$ml_python_prefix/lib" "\$ml_python_prefix/lib64" "\$ml_python_real_prefix/lib" "\$ml_python_real_prefix/lib64"; do
  [[ -n "\$d" && -d "\$d" ]] || continue
  case ":\$ld_joined:" in *":\$d:"*) ;; *) ld_joined="\${ld_joined:+\$ld_joined:}\$d" ;; esac
done
[[ -n "\$ld_joined" ]] && export LD_LIBRARY_PATH="\$ld_joined:\${LD_LIBRARY_PATH:-}"
unset PYTHONHOME
model_registry_arg=()
if [[ -n "\$model_registry" ]]; then
  model_registry_arg=(--model-registry "\$model_registry")
fi
event_quality_args=()
if [[ -n "${event_quality_cut_json}" ]]; then
  event_quality_args=(--event-quality-cut-json "${event_quality_cut_json}" --event-quality-audit-output "${report_root}/event_quality_filter_audit.json")
fi
cache_manifest="${report_root}/score_caches.list"
find "${cache_dir}" -type f -name 'score_cache_*.npz' | sort -V > "\$cache_manifest" || true
expected=${idx}
found=\$(wc -l < "\$cache_manifest" | tr -d ' ')
summary="${report_root}/validation_summary.txt"
if [[ "\$found" != "\$expected" || "\$found" == "0" ]]; then
  {
    echo "RECOILJETS_AUAU_TIGHT_BDT_SIM_VALIDATION_V1"
    echo "status=CHECK"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    [[ -n "\$model_registry" ]] && echo "model_registry=\$model_registry"
    echo "report_dir=${report_root}"
    echo "expected_score_caches=\$expected"
    echo "found_score_caches=\$found"
    echo "notes=missing score cache shards"
    echo "next_action=Inspect ${sub_root}/validate_*.err and rerun validateOnSimCondor for failed shards."
  } > "\$summary"
else
  rc=0
  "\$ml_python" "${VALIDATE_SCRIPT}" \\
    --source "${source}" \\
    --model-dir "${model_dir}" \\
    "\${model_registry_arg[@]}" \\
    "\${event_quality_args[@]}" \\
    --merge-score-caches "\$cache_manifest" \\
    --outdir "${report_root}" || rc=\$?
fi
if [[ -s "\$summary" ]]; then
  status=\$(awk -F= '/^status=/ {print \$2; exit}' "\$summary")
  status="\${status:-CHECK}"
else
  status=CHECK
fi
if command -v mail >/dev/null 2>&1 && [[ -s "\$summary" ]]; then
  mail -s "[RecoilJets][auauTightBDT_validateOnSimCondor][\${status}]" "${NOTIFY_EMAILS}" < "\$summary" || true
fi
cat "\$summary"
[[ "\$status" == "READY" ]] || exit 3
EOF
  chmod +x "$merge"

  local merge_sub="${sub_root}/validate_merge.sub"
  cat > "$merge_sub" <<EOF
universe = ${merge_universe}
executable = ${merge}
output = ${sub_root}/validate_merge.out
error = ${sub_root}/validate_merge.err
log = ${sub_root}/validate_merge.log
notification = Never
EOF
  if [[ -n "$merge_reqmem" ]]; then
    echo "request_memory = ${merge_reqmem}" >> "$merge_sub"
  fi
  cat >> "$merge_sub" <<EOF
queue
EOF

  local dag="${sub_root}/auau_tight_bdt_validateOnSimCondor.dag"
  : > "$dag"
  while read -r node_name node_sub; do
    [[ -n "${node_name:-}" ]] || continue
    {
      echo "JOB ${node_name} ${node_sub}"
      echo "RETRY ${node_name} 1"
    } >> "$dag"
  done < "$shard_nodes_file"
  echo "FINAL MERGE ${merge_sub}" >> "$dag"
  say "Condor validation DAG: $dag"
  say "source: $source"
  say "model dir: $model_dir"
  [[ -n "$model_registry" ]] && say "model registry: $model_registry"
  say "report root: $report_root"
  say "root files=${nroots}  shards=${idx}  groupSize=${group_size}  scoreMaxPerShard=${score_max_per_shard}  request_memory=${reqmem}"
  if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
    echo "RECOILJETS_AUAU_TIGHT_BDT_VALIDATE_DRYRUN_V1"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    [[ -n "$model_registry" ]] && echo "model_registry=${model_registry}"
    echo "report_root=${report_root}"
    echo "dag=${dag}"
    echo "root_files=${nroots}"
    echo "shards=${idx}"
    return 0
  fi
  condor_submit_dag "$dag"
}

derive_working_points_from_validation() {
  local validation=""
  local target="0.80"
  local wp_mode="${RJ_AUAU_BDT_WP_MODE:-centpt}"
  local pt_bins="${RJ_AUAU_BDT_WP_PT_BINS:-}"
  local cent_bins="${RJ_AUAU_BDT_WP_CENT_BINS:-0,20,50,80}"
  local tok
  for tok in "$@"; do
    case "$tok" in
      VALIDATION=*) validation="${tok#VALIDATION=}" ;;
      TARGET=*) target="${tok#TARGET=}" ;;
      MODE=*) wp_mode="${tok#MODE=}" ;;
      PT_BINS=*) pt_bins="${tok#PT_BINS=}" ;;
      CENT_BINS=*) cent_bins="${tok#CENT_BINS=}" ;;
      *) ;;
    esac
  done
  [[ -n "$validation" ]] || die "deriveWorkingPointsFromValidation requires VALIDATION=/path/to/model_validation_report"
  [[ -d "$validation" ]] || die "Validation report directory does not exist: $validation"
  setup_ml_python_env
  say "Deriving BDT working points from validation report"
  say "  validation : $validation"
  say "  target     : $target"
  say "  mode       : $wp_mode"
  say "  pt bins    : ${pt_bins:-validator default}"
  say "  cent bins  : $cent_bins"
  "$ML_PYTHON" "$VALIDATE_SCRIPT" \
    --derive-working-points-from-report "$validation" \
    --target-signal-efficiency "$target" \
    --working-point-mode "$wp_mode" \
    --working-point-pt-bins "$pt_bins" \
    --working-point-cent-bins "$cent_bins"
}

generate_working_point_config() {
  local template=""
  local wp=""
  local product_map=""
  local out=""
  local model_dir=""
  local tok
  for tok in "$@"; do
    case "$tok" in
      TEMPLATE=*) template="${tok#TEMPLATE=}" ;;
      WORKING_POINTS=*) wp="${tok#WORKING_POINTS=}" ;;
      PRODUCT_MAP=*) product_map="${tok#PRODUCT_MAP=}" ;;
      OUT=*) out="${tok#OUT=}" ;;
      MODEL_DIR=*) model_dir="${tok#MODEL_DIR=}" ;;
      *) ;;
    esac
  done
  [[ -n "$template" ]] || die "generateWorkingPointConfig requires TEMPLATE=/path/to/template.yaml"
  [[ -n "$wp" ]] || die "generateWorkingPointConfig requires WORKING_POINTS=/path/to/bdt_working_points_target80.json"
  [[ -n "$product_map" ]] || die "generateWorkingPointConfig requires PRODUCT_MAP=variant=product[,variant=product]"
  [[ -n "$out" ]] || die "generateWorkingPointConfig requires OUT=/path/to/generated.yaml"
  local helper="${RJ_REPO_BASE}/scripts/make_auau_bdt_target_wp_config.py"
  local cmd=(python3 "$helper" --template "$template" --working-points "$wp" --out "$out" --product-map "$product_map")
  if [[ -n "$model_dir" ]]; then
    cmd+=(--model-dir "$model_dir")
  fi
  say "Generating RecoilJets YAML with validation-derived BDT working points"
  "${cmd[@]}"
}

main() {
  local mode="${1:-}"
  [[ -n "$mode" ]] || { usage; exit 2; }
  shift || true
  case "$mode" in
    -h|--help|help) usage ;;
    localTest|local) run_local_test "$@" ;;
    smokeTest) require_codex_submission_provenance; run_condor_extract "smokeTest" "$@" ;;
    condorExtract|condorDoAll) require_codex_submission_provenance; run_condor_extract "condorExtract" "$@" ;;
    trainFromExtraction) train_from_extraction "$@" ;;
    trainCentInput3x3FromExtraction) train_cent_input_3x3_from_extraction "$@" ;;
    trainWidthStudyPt1530FromExtraction) train_width_study_pt1530_from_extraction "$@" ;;
    trainWidthStudyWindowsFromExtraction) train_width_study_windows_from_extraction "$@" ;;
    trainEtFineCentStudyFromExtraction) train_etfine_centstudy_from_extraction "$@" ;;
    trainExpandedFromExtraction) train_expanded_from_extraction "$@" ;;
    trainExpandedFromExtractionCondor) require_codex_submission_provenance; train_expanded_from_extraction_condor "$@" ;;
    finalizeExpandedTraining) finalize_expanded_training "$@" ;;
    applyCheck|smokeTestApplyExisting) apply_check "$@" ;;
    validateOnSim|validateSim|simValidation) validate_on_sim "$@" ;;
    validateOnSimCondor|condorValidateOnSim|validateSimCondor|validateEtFineCentStudyOnSimCondor) require_codex_submission_provenance; validate_on_sim_condor "$@" ;;
    deriveWorkingPointsFromValidation|deriveWPFromValidation) derive_working_points_from_validation "$@" ;;
    generateWorkingPointConfig|generateTargetWPConfig) generate_working_point_config "$@" ;;
    *) usage; die "Unknown mode: $mode" ;;
  esac
}

main "$@"
