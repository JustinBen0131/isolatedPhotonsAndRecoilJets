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

usage() {
  cat <<'EOF'
Usage:
  ./scripts/auau_tight_bdt_pipeline.sh localTest [Nevents] [NFILES=N]
  ./scripts/auau_tight_bdt_pipeline.sh smokeTest [groupSize N]
  ./scripts/auau_tight_bdt_pipeline.sh condorExtract [groupSize N]
  ./scripts/auau_tight_bdt_pipeline.sh trainFromExtraction SOURCE=/path
  ./scripts/auau_tight_bdt_pipeline.sh applyCheck MODEL_DIR=/path

Sidecar AuAu tight-BDT workflow:
  extraction reads only embeddedPhoton12/20 and embeddedJet12/20 samples,
  writes AuAuPhotonIDTrainingTree ROOT files, and avoids normal cfg-tag
  histogram production.
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
  "$ML_PYTHON" - "$manifest" "$report" <<'PY'
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
report = Path(sys.argv[2])
paths = [Path(x.strip()) for x in manifest.read_text().splitlines() if x.strip()]
rows = []
total = 0
try:
    import uproot
except Exception as exc:
    report.write_text(json.dumps({"status": "CHECK", "reason": f"uproot import failed: {exc}", "files": [str(p) for p in paths]}, indent=2) + "\n")
    print(f"[WARN] uproot validation skipped: {exc}")
    sys.exit(0)
for path in paths:
    with uproot.open(path) as f:
        if "AuAuPhotonIDTrainingTree" not in f:
            rows.append({"file": str(path), "entries": 0, "missing_tree": True})
            continue
        n = int(f["AuAuPhotonIDTrainingTree"].num_entries)
        total += n
        rows.append({"file": str(path), "entries": n})
status = "PASS" if total > 0 and all(not r.get("missing_tree") for r in rows) else "FAIL"
report.write_text(json.dumps({"status": status, "total_entries": total, "files": rows}, indent=2, sort_keys=True) + "\n")
print(f"[OK] training tree validation: status={status} entries={total} files={len(paths)}")
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
  mkdir -p "$model_dir"
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
root_manifest="${manifest_dir}/training_roots.list"
find "${extraction_root}" -type f -name '*.root' | sort -V > "\$root_manifest" || true
nroots=\$(wc -l < "\$root_manifest" | tr -d ' ')
tree_entries=0
validation_note=""
if [[ "\$nroots" != "0" ]]; then
  validation_note=\$("${ML_PYTHON}" - "\$root_manifest" <<'PY' || true
import sys
from pathlib import Path
try:
    import uproot
except Exception as exc:
    print(f"UPROOT_IMPORT_FAILED {exc}")
    raise SystemExit(0)
manifest = Path(sys.argv[1])
total = 0
missing = 0
for raw in manifest.read_text().splitlines():
    path = raw.strip()
    if not path:
        continue
    try:
        with uproot.open(path) as f:
            if "AuAuPhotonIDTrainingTree" not in f:
                missing += 1
                continue
            total += int(f["AuAuPhotonIDTrainingTree"].num_entries)
    except Exception as exc:
        print(f"ROOT_OPEN_FAILED {path} {exc}")
print(f"TREE_ENTRIES {total} MISSING_TREE_FILES {missing}")
PY
)
  tree_entries=\$(printf '%s\n' "\$validation_note" | awk '/TREE_ENTRIES/ {print \$2; exit}')
  tree_entries="\${tree_entries:-0}"
fi
status=READY
if [[ "\$nroots" == "0" || "\$tree_entries" == "0" ]]; then status=CHECK; fi
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
  make_root_manifest "$search_root" "$manifest"
  validate_training_tree "$manifest" "${source}/reports/training_tree_validation.json"
  train_from_manifest "$manifest"
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
    *) usage; die "Unknown mode: $mode" ;;
  esac
}

main "$@"
