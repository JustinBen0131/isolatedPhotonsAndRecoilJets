#!/usr/bin/env bash
# Mirror PPG12 run-28 double-interaction SIM lists into the RecoilJets
# simListFiles layout without writing into the PPG12 checkout.

set -Eeuo pipefail

PPG12_RUN28_ROOT="${PPG12_RUN28_ROOT:-/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28}"
OUTROOT="${OUTROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/simListFiles}"
PPG12_DOUBLE_SIGNAL_DUAL_ROOT="${PPG12_DOUBLE_SIGNAL_DUAL_ROOT:-/sphenix/lustre01/sphnxpro/mdc2/js_pp200_signal_dual}"
PPG12_DOUBLE_RUN_DIR="${PPG12_DOUBLE_RUN_DIR:-run0028}"

REQUESTED_SAMPLES=()

usage() {
  cat <<'EOF'
Usage: makePPG12DoubleSimLists.sh [photon5_double|photon10_double|photon20_double ...]

Mirrors read-only PPG12 condorout/OutDir*/list files into:
  $OUTROOT/run28_photonjet5_double
  $OUTROOT/run28_photonjet10_double
  $OUTROOT/run28_photonjet20_double

Environment:
  PPG12_RUN28_ROOT  PPG12 run28 macro root
  OUTROOT           RecoilJets simListFiles root
  PPG12_DOUBLE_SIGNAL_DUAL_ROOT
                    Persistent js_pp200_signal_dual root
  PPG12_DOUBLE_RUN_DIR
                    Run directory below each persistent lane (default run0028)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    photon5_double|photon10_double|photon20_double)
      REQUESTED_SAMPLES+=( "$1" )
      shift
      ;;
    *)
      echo "[ERROR] Unsupported sample: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if (( ${#REQUESTED_SAMPLES[@]} == 0 )); then
  REQUESTED_SAMPLES=( photon5_double photon10_double photon20_double )
fi

out_sample_name() {
  case "$1" in
    photon5_double)  printf '%s\n' run28_photonjet5_double ;;
    photon10_double) printf '%s\n' run28_photonjet10_double ;;
    photon20_double) printf '%s\n' run28_photonjet20_double ;;
    *) return 1 ;;
  esac
}

sample_photonjet_dir() {
  case "$1" in
    photon5_double)  printf '%s\n' photonjet5 ;;
    photon10_double) printf '%s\n' photonjet10 ;;
    photon20_double) printf '%s\n' photonjet20 ;;
    *) return 1 ;;
  esac
}

resolve_list_entry() {
  local entry="$1"
  local base="${2:-}"

  [[ -n "$entry" ]] || return 0
  if [[ "$entry" == "NONE" || "$entry" == /* ]]; then
    printf '%s\n' "$entry"
  elif [[ -n "$base" ]]; then
    printf '%s/%s\n' "$base" "$entry"
  else
    printf '%s\n' "$entry"
  fi
}

concat_lane() {
  local src_dir="$1"
  local list_name="$2"
  local out="$3"
  local base="${4:-}"
  local line

  : > "$out"
  while IFS= read -r f; do
    while IFS= read -r line || [[ -n "$line" ]]; do
      resolve_list_entry "$line" "$base" >> "$out"
    done < "$f"
  done < <(find "${src_dir}/condorout" -mindepth 2 -maxdepth 2 -type f -name "$list_name" | sort -V)
}

validate_existing_lane() {
  local lane="$1"
  local path="$2"
  local failures=0
  local line_no=0
  local entry
  while IFS= read -r entry || [[ -n "$entry" ]]; do
    line_no=$((line_no + 1))
    [[ -n "$entry" && "$entry" != "NONE" ]] || continue
    if [[ "$entry" != /* ]]; then
      echo "[ERROR] ${lane}:${line_no} is not absolute after PPG12 double resolution: ${entry}" >&2
      failures=$((failures + 1))
    elif [[ ! -e "$entry" ]]; then
      echo "[ERROR] ${lane}:${line_no} does not exist after PPG12 double resolution: ${entry}" >&2
      failures=$((failures + 1))
    fi
    (( failures < 20 )) || break
  done < "$path"
  [[ "$failures" == "0" ]] || exit 5
}

make_none_lane() {
  local n="$1"
  local out="$2"
  awk -v n="$n" 'BEGIN{for(i=0;i<n;i++) print "NONE"}' > "$out"
}

copy_lane_to_matched() {
  local outdir="$1"
  local base="$2"
  cp "${outdir}/${base}.list" "${outdir}/${base}.matched.list"
}

for sample in "${REQUESTED_SAMPLES[@]}"; do
  src_dir="${PPG12_RUN28_ROOT}/${sample}"
  [[ -d "$src_dir" ]] || { echo "[ERROR] Missing PPG12 sample dir: $src_dir" >&2; exit 3; }
  [[ -d "${src_dir}/condorout" ]] || { echo "[ERROR] Missing PPG12 condorout dir: ${src_dir}/condorout" >&2; exit 3; }

  out_sample="$(out_sample_name "$sample")"
  photonjet_dir="$(sample_photonjet_dir "$sample")"
  g4_base="${PPG12_DOUBLE_SIGNAL_DUAL_ROOT}/g4hits/${PPG12_DOUBLE_RUN_DIR}/${photonjet_dir}"
  jets_base="${PPG12_DOUBLE_SIGNAL_DUAL_ROOT}/nopileup/jets/${PPG12_DOUBLE_RUN_DIR}/${photonjet_dir}"
  outdir="${OUTROOT}/${out_sample}"
  rm -rf "$outdir"
  mkdir -p "$outdir"

  concat_lane "$src_dir" g4hits.list "${outdir}/G4Hits.list" "$g4_base"
  concat_lane "$src_dir" dst_truth_jet.list "${outdir}/DST_JETS.list" "$jets_base"
  n_g4="$(wc -l < "${outdir}/G4Hits.list" | tr -d ' ')"
  n_jets="$(wc -l < "${outdir}/DST_JETS.list" | tr -d ' ')"
  [[ "$n_g4" -gt 0 ]] || { echo "[ERROR] Empty G4Hits lane for $sample" >&2; exit 4; }
  [[ "$n_jets" == "$n_g4" ]] || { echo "[ERROR] DST_JETS count $n_jets != G4Hits count $n_g4 for $sample" >&2; exit 4; }
  validate_existing_lane "G4Hits" "${outdir}/G4Hits.list"
  validate_existing_lane "DST_JETS" "${outdir}/DST_JETS.list"

  # PPG12 double jobs rebuild the detector-level objects from G4Hits and only
  # register the truth-jet lane. Some side list files contain stale relative
  # names, so keep the RecoilJets double input contract explicitly G4-only.
  for lane in DST_CALO_CLUSTER DST_GLOBAL DST_MBD_EPD; do
    make_none_lane "$n_g4" "${outdir}/${lane}.list"
    n_lane="$n_g4"
    [[ "$n_lane" == "$n_g4" ]] || { echo "[ERROR] ${lane} count $n_lane != G4Hits count $n_g4 for $sample" >&2; exit 4; }
  done

  for lane in DST_CALO_CLUSTER G4Hits DST_JETS DST_GLOBAL DST_MBD_EPD; do
    copy_lane_to_matched "$outdir" "$lane"
  done

  {
    echo "sample=${sample}"
    echo "out_sample=${out_sample}"
    echo "source=${src_dir}"
    echo "entries=${n_g4}"
    for lane in DST_CALO_CLUSTER G4Hits DST_JETS DST_GLOBAL DST_MBD_EPD; do
      printf '%s_count=%s\n' "$lane" "$(wc -l < "${outdir}/${lane}.matched.list" | tr -d ' ')"
      printf '%s_first=' "$lane"
      head -1 "${outdir}/${lane}.matched.list"
    done
  } > "${outdir}/summary.txt"

  echo "[OK] ${sample} -> ${outdir} entries=${n_g4}"
done
