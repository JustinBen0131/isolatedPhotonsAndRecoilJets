#!/usr/bin/env bash
set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
SIM_ROOT="${RJ_SIM_ROOT:-${BASE}/simListFiles}"
TRAIN_BASE="${RJ_AUAU_TIGHT_BDT_TRAIN_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining}"
LOCAL_BASE="${RJ_AUAU_TIGHT_BDT_LOCAL_BASE:-${BASE}/local_bdt_training_outputs}"
MODEL_BASE="${RJ_AUAU_BDT_MODEL_BASE:-${BASE}/bdt_models}"
MASTER_YAML="${RJ_AUAU_TIGHT_BDT_CONFIG_SRC:-${BASE}/macros/analysis_config.yaml}"
TRAIN_MACRO="${RJ_AUAU_TIGHT_BDT_MACRO:-${BASE}/macros/Fun4All_auauTightBDTTraining.C}"
TRAIN_SCRIPT="${RJ_AUAU_TIGHT_BDT_TRAIN_SCRIPT:-${BASE}/scripts/train_auau_photon_bdt.py}"
VALIDATE_SCRIPT="${RJ_AUAU_TIGHT_BDT_VALIDATE_SCRIPT:-${BASE}/scripts/validate_auau_tight_bdt_on_sim.py}"
ML_PYTHON="${RJ_ML_PYTHON:-python3}"
NOTIFY_EMAILS="${RJ_NOTIFY_EMAILS:-just0131@gmail.com}"

SIGNAL_SAMPLES=(run28_embeddedPhoton12 run28_embeddedPhoton20)
BACKGROUND_SAMPLES=(run28_embeddedJet12 run28_embeddedJet20)
ALL_SAMPLES=("${SIGNAL_SAMPLES[@]}" "${BACKGROUND_SAMPLES[@]}")

ts() { date +%Y%m%d_%H%M%S; }
say() { printf '\033[1;36m[auauTightBDT]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[auauTightBDT][WARN]\033[0m %s\n' "$*" >&2; }
err() { printf '\033[1;31m[auauTightBDT][ERR]\033[0m %s\n' "$*" >&2; }
die() { err "$*"; exit 2; }

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
  ./scripts/auau_tight_bdt_pipeline.sh applyCheck MODEL_DIR=/path
  ./scripts/auau_tight_bdt_pipeline.sh validateOnSim SOURCE=/path MODEL_DIR=/path
  ./scripts/auau_tight_bdt_pipeline.sh validateOnSimCondor SOURCE=/path MODEL_DIR=/path [groupSize N]

Sidecar AuAu tight-BDT workflow:
  extraction reads only embeddedPhoton12/20 and embeddedJet12/20 samples,
  writes AuAuPhotonIDTrainingTree ROOT files, and avoids normal cfg-tag
  histogram production. validateOnSim scores those same trees with the
  final TMVA ROOT models and writes quick ROC/AUC simulation diagnostics.
  validateOnSimCondor shards the ROOT scoring over Condor, then merges the
  score caches into the same final validation report/plots.
EOF
}

wrapper_path() {
  if [[ -f "${BASE}/RecoilJets_Condor_AuAu.sh" ]]; then
    printf '%s\n' "${BASE}/RecoilJets_Condor_AuAu.sh"
  elif [[ -f "${BASE}/scripts/RecoilJets_Condor_AuAu.sh" ]]; then
    printf '%s\n' "${BASE}/scripts/RecoilJets_Condor_AuAu.sh"
  else
    die "Cannot find RecoilJets_Condor_AuAu.sh under ${BASE}"
  fi
}

sample_dataset() {
  case "$1" in
    run28_embeddedPhoton12|run28_embeddedPhoton20) printf '%s\n' "isSimEmbedded" ;;
    run28_embeddedJet12|run28_embeddedJet20)       printf '%s\n' "isSimEmbeddedInclusive" ;;
    *) die "Unknown AuAu BDT training sample: $1" ;;
  esac
}

sample_class() {
  case "$1" in
    run28_embeddedPhoton12|run28_embeddedPhoton20) printf '%s\n' "signal" ;;
    run28_embeddedJet12|run28_embeddedJet20)       printf '%s\n' "background" ;;
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
  find "$root" -type f -name '*.root' | sort -V > "$out" || true
  [[ -s "$out" ]] || die "No ROOT files found under $root"
}

validate_training_tree() {
  local manifest="$1"
  local report="$2"
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
        sample_class = "signal" if "/signal/" in str(path) else ("background" if "/background/" in str(path) else "unknown")
        rows.append({
            "file": str(path),
            "entries": n,
            "class": sample_class,
            "signal_entries": n_signal,
            "background_entries": n_background
        })
status = "PASS" if total > 0 and signal > 0 and background > 0 and all(not r.get("missing_tree") and not r.get("missing_is_signal") for r in rows) else "FAIL"
report.write_text(json.dumps({
    "status": status,
    "total_entries": total,
    "signal_entries": signal,
    "background_entries": background,
    "files": rows
}, indent=2, sort_keys=True) + "\n")
print(f"[OK] training tree validation: status={status} entries={total} signal={signal} background={background} files={len(paths)}")
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
  local sub_root="${BASE}/condor_sub/auauTightBDT_${stamp}"
  local extraction_root="${run_root}/extraction"
  local manifest_dir="${run_root}/manifests"
  local report_dir="${run_root}/reports"
  local yaml="${manifest_dir}/analysis_config_auau_tight_bdt_training.yaml"
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
    split -l "$group_size" -d -a 5 "$master" "$split_prefix"
    local queued=0
    for raw in "${split_prefix}"*; do
      [[ -s "$raw" ]] || { rm -f "$raw"; continue; }
      queued=$((queued + 1))
      if (( max_jobs_per_sample > 0 && queued > max_jobs_per_sample )); then
        rm -f "$raw"
        continue
      fi
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
  cat > "$sub" <<EOF
universe = vanilla
executable = /usr/bin/bash
arguments = ${wrapper} \$(sample) \$(chunk) \$(dataset) \$(Cluster) \$(nevents) \$(chunkidx) NONE \$(dest)
output = ${sub_root}/extract_\$(Cluster)_\$(Process).out
error = ${sub_root}/extract_\$(Cluster)_\$(Process).err
log = ${sub_root}/extract_\$(Cluster).log
request_memory = ${reqmem}
notification = Never
environment = "RJ_CONFIG_YAML=${yaml} RJ_MACRO_PATH=${TRAIN_MACRO} RJ_AUAU_BDT_EXTRACT_ONLY=1 RJ_AUAU_BDT_TRAINING_TREE=1 RJ_DISABLE_ID_FANOUT=1 RJ_DISABLE_ISO_CONE_INTERNALIZATION=1 RJ_DISABLE_JET_PT_INTERNALIZATION=1 RJ_DISABLE_DPHI_INTERNALIZATION=1 RJ_PROFILE_JOB=1"
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
find "${extraction_root}" -type f -name '*.root' | sort -V > "\$root_manifest" || true
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
    except Exception as exc:
        print(f"ROOT_OPEN_FAILED {path} {exc}")
print(f"TREE_ENTRIES {total} SIGNAL_ENTRIES {signal} BACKGROUND_ENTRIES {background} MISSING_TREE_FILES {missing}")
PY
)
  tree_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {print \$2; exit}')
  tree_entries="\${tree_entries:-0}"
  signal_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {for (i=1; i<=NF; ++i) if (\$i=="SIGNAL_ENTRIES") {print \$(i+1); exit}}')
  background_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {for (i=1; i<=NF; ++i) if (\$i=="BACKGROUND_ENTRIES") {print \$(i+1); exit}}')
  signal_entries="\${signal_entries:-0}"
  background_entries="\${background_entries:-0}"
fi
status=READY
if [[ "\$nroots" == "0" || "\$tree_entries" == "0" || "\${signal_entries:-0}" == "0" || "\${background_entries:-0}" == "0" ]]; then status=CHECK; fi
summary="${report_dir}/final_summary.txt"
{
  echo "RECOILJETS_STAGE_EMAIL_V1"
  echo "dataset=isSimEmbeddedAndInclusive"
  echo "stage=auauTightBDT_${mode}"
  echo "status=\${status}"
  echo "training_root=${run_root}"
  echo "root_manifest=\${root_manifest}"
  echo "root_count=\${nroots}"
  echo "tree_entries=\${tree_entries}"
  echo "signal_entries=\${signal_entries:-0}"
  echo "background_entries=\${background_entries:-0}"
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
  local outdir=""
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|modelDir=*|model_dir=*) model_dir="${tok#*=}" ;;
      OUTDIR=*|outdir=*) outdir="${tok#*=}" ;;
    esac
  done
  [[ -n "$source" ]] || die "validateOnSim requires SOURCE=/path/to/auauTightBDT extraction"
  [[ -d "$source" ]] || die "SOURCE is not a directory: $source"
  [[ -n "$model_dir" ]] || die "validateOnSim requires MODEL_DIR=/path/to/tight models"
  [[ -d "$model_dir" ]] || die "MODEL_DIR is not a directory: $model_dir"
  [[ -s "$VALIDATE_SCRIPT" ]] || die "Missing validation script: $VALIDATE_SCRIPT"

  setup_ml_python_env
  local -a args
  args=( "$VALIDATE_SCRIPT" --source "$source" --model-dir "$model_dir" )
  if [[ -n "$outdir" ]]; then
    args+=( --outdir "$outdir" )
  fi

  say "Validating tight-BDT models on embedded-sim extraction trees"
  say "  source    : $source"
  say "  model dir : $model_dir"
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
  local outdir=""
  local group_size="${RJ_AUAU_TIGHT_BDT_VALIDATE_GROUP_SIZE:-100}"
  local total_score_max="${RJ_AUAU_TIGHT_BDT_VALIDATE_TOTAL_SCORE_MAX_ROWS:-400000}"
  local reqmem="${RJ_AUAU_TIGHT_BDT_VALIDATE_REQUEST_MEMORY:-2500MB}"
  for tok in "$@"; do
    case "$tok" in
      SOURCE=*) source="${tok#SOURCE=}" ;;
      MODEL_DIR=*|MODELDIR=*|modelDir=*|model_dir=*) model_dir="${tok#*=}" ;;
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
  [[ -s "$VALIDATE_SCRIPT" ]] || die "Missing validation script: $VALIDATE_SCRIPT"
  [[ "$group_size" =~ ^[0-9]+$ && "$group_size" -gt 0 ]] || die "groupSize must be a positive integer"

  local stamp="${RJ_AUAU_TIGHT_BDT_VALIDATE_STAMP:-$(ts)}"
  local report_root="${outdir:-${source}/reports/model_validation_condor_${stamp}}"
  local sub_root="${BASE}/condor_sub/auauTightBDTValidate_${stamp}"
  local shard_dir="${sub_root}/shards"
  local cache_dir="${report_root}/score_caches"
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
"\$ml_python" "${VALIDATE_SCRIPT}" \\
  --source "${source}" \\
  --model-dir "${model_dir}" \\
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
universe = scheduler
executable = ${merge}
output = ${sub_root}/validate_merge.out
error = ${sub_root}/validate_merge.err
log = ${sub_root}/validate_merge.log
notification = Never
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
  say "report root: $report_root"
  say "root files=${nroots}  shards=${idx}  groupSize=${group_size}  scoreMaxPerShard=${score_max_per_shard}  request_memory=${reqmem}"
  if [[ "${RJ_DAG_DRYRUN:-0}" == "1" ]]; then
    echo "RECOILJETS_AUAU_TIGHT_BDT_VALIDATE_DRYRUN_V1"
    echo "source=${source}"
    echo "model_dir=${model_dir}"
    echo "report_root=${report_root}"
    echo "dag=${dag}"
    echo "root_files=${nroots}"
    echo "shards=${idx}"
    return 0
  fi
  condor_submit_dag "$dag"
}

main() {
  local mode="${1:-}"
  [[ -n "$mode" ]] || { usage; exit 2; }
  shift || true
  case "$mode" in
    -h|--help|help) usage ;;
    localTest|local) run_local_test "$@" ;;
    smokeTest) run_condor_extract "smokeTest" "$@" ;;
    condorExtract|condorDoAll) run_condor_extract "condorExtract" "$@" ;;
    trainFromExtraction) train_from_extraction "$@" ;;
    applyCheck|smokeTestApplyExisting) apply_check "$@" ;;
    validateOnSim|validateSim|simValidation) validate_on_sim "$@" ;;
    validateOnSimCondor|condorValidateOnSim|validateSimCondor) validate_on_sim_condor "$@" ;;
    *) usage; die "Unknown mode: $mode" ;;
  esac
}

main "$@"
