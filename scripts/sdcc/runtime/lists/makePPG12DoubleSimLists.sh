#!/usr/bin/env bash
# Mirror PPG12 run-28 double-interaction SIM lists into the RecoilJets
# simListFiles layout without writing into the PPG12 checkout.

set -Eeuo pipefail

PPG12_RUN28_ROOT="${PPG12_RUN28_ROOT:-/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28}"
OUTROOT="${OUTROOT:-/sphenix/u/patsfan753/scratch/thesisAnalysis/simListFiles}"

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

concat_lane() {
  local src_dir="$1"
  local list_name="$2"
  local out="$3"

  : > "$out"
  while IFS= read -r f; do
    cat "$f" >> "$out"
  done < <(find "${src_dir}/condorout" -mindepth 2 -maxdepth 2 -type f -name "$list_name" | sort -V)
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
  outdir="${OUTROOT}/${out_sample}"
  rm -rf "$outdir"
  mkdir -p "$outdir"

  concat_lane "$src_dir" g4hits.list "${outdir}/G4Hits.list"
  concat_lane "$src_dir" dst_truth_jet.list "${outdir}/DST_JETS.list"
  n_g4="$(wc -l < "${outdir}/G4Hits.list" | tr -d ' ')"
  n_jets="$(wc -l < "${outdir}/DST_JETS.list" | tr -d ' ')"
  [[ "$n_g4" -gt 0 ]] || { echo "[ERROR] Empty G4Hits lane for $sample" >&2; exit 4; }
  [[ "$n_jets" == "$n_g4" ]] || { echo "[ERROR] DST_JETS count $n_jets != G4Hits count $n_g4 for $sample" >&2; exit 4; }

  concat_lane "$src_dir" dst_calo_cluster.list "${outdir}/DST_CALO_CLUSTER.list"
  concat_lane "$src_dir" dst_global.list "${outdir}/DST_GLOBAL.list"
  concat_lane "$src_dir" dst_mbd_epd.list "${outdir}/DST_MBD_EPD.list"

  for lane in DST_CALO_CLUSTER DST_GLOBAL DST_MBD_EPD; do
    n_lane="$(wc -l < "${outdir}/${lane}.list" | tr -d ' ')"
    if [[ "$n_lane" == "0" ]]; then
      make_none_lane "$n_g4" "${outdir}/${lane}.list"
      n_lane="$n_g4"
    fi
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
