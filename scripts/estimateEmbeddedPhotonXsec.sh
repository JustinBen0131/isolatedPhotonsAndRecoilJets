#!/usr/bin/env bash
# ==============================================================================
# estimateEmbeddedPhotonXsec.sh
# ==============================================================================
# PURPOSE
#   Estimate the effective generator cross sections for the embedded photon
#   samples now used by isSimEmbedded:
#
#     PhotonJet12: phpythia8_10GeV_JS_MDC2.cfg + PHPy8ParticleTrigger pT > 12
#     PhotonJet20: phpythia8_20GeV_JS_MDC2.cfg + PHPy8ParticleTrigger pT > 20
#
#   This is generator-only. It does not run detector simulation, embedding,
#   clustering, or RecoilJets. It reproduces the producer-side photon trigger
#   logic directly from PHPy8ParticleTrigger:
#
#     id == 22
#     stable-particle-only disabled
#     -1.5 < eta < 1.5
#     pT > threshold
#     abs(mother id) in {1, ..., 22}
#
# OUTPUT
#   By default, overwrites a compact output directory:
#
#     <repo>/pythia_xsec_estimate/
#       results.csv
#       summary.txt
#       run.log
#
# USAGE ON SDCC
#   From the thesisAnalysis base directory:
#
#     ./scripts/estimateEmbeddedPhotonXsec.sh --events 100000
#     ./scripts/estimateEmbeddedPhotonXsec.sh --events 1000000
#
#   Local smoke test for the Condor workflow:
#
#     ./scripts/estimateEmbeddedPhotonXsec.sh firstPass LOCAL --xsec-events 10000
#
#   Condor two-pass workflow:
#
#     ./scripts/estimateEmbeddedPhotonXsec.sh firstPass
#     ./scripts/estimateEmbeddedPhotonXsec.sh secondPass
#
#   Preferred final-estimate mode:
#
#     ./scripts/estimateEmbeddedPhotonXsec.sh --target-pass 10000 --max-raw-events 50000000
#
#   Optional:
#
#     --seed N
#     --outdir DIR
#     --sample all|PhotonJet12|PhotonJet20
#     --mode compiled      # default; faster event loop, first run compiles
#     --mode interpreted   # useful for tiny smoke tests; avoids ACLiC compile
#     --xsec-shards N      # firstPass default: 50 per sample
#     --xsec-events N      # firstPass default: 1000000 raw attempts per shard
#     --xsec-memory MEM    # firstPass default: 900MB
#     --manifest PATH      # secondPass explicit firstPass manifest
#
# NOTES
#   sigma_eff = sigma_gen * (trigger-passing generator weight sum / total
#   generator weight sum). Pythia reports sigma_gen in mb; this script also
#   prints pb.
# ==============================================================================

set -Eeuo pipefail

NEVENTS=1000000
TARGET_PASS=0
MAX_RAW_EVENTS=0
SEED=75301
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTDIR="${BASE_DIR}/pythia_xsec_estimate"
RUN_MODE="compiled"
RUN_MODE_EXPLICIT=0
STREAM_LOG="true"
SAMPLE_FILTER="all"
ACTION="run"
LOCAL_MODE=0
WORKER_ARGS=()
XSEC_SHARDS=50
XSEC_RAW_EVENTS=1000000
XSEC_MEMORY="900MB"
XSEC_BASE_SEED=75301
MANIFEST_PATH=""

usage() {
  sed -n '2,65p' "$0" | sed 's/^# \{0,1\}//g'
}

if [[ $# -gt 0 ]]; then
  case "$1" in
    run|firstPass|secondPass|xsecWorker)
      ACTION="$1"
      shift
      ;;
  esac
fi

if [[ "${ACTION}" == "xsecWorker" ]]; then
  WORKER_ARGS=("$@")
  set --
fi

if [[ "${ACTION}" == "firstPass" && "${1:-}" == "LOCAL" ]]; then
  LOCAL_MODE=1
  shift
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --events|-n)
      NEVENTS="${2:-}"
      shift 2
      ;;
    --target-pass)
      TARGET_PASS="${2:-}"
      shift 2
      ;;
    --max-raw-events)
      MAX_RAW_EVENTS="${2:-}"
      shift 2
      ;;
    --seed)
      SEED="${2:-}"
      shift 2
      ;;
    --outdir)
      OUTDIR="${2:-}"
      shift 2
      ;;
    --sample)
      SAMPLE_FILTER="${2:-}"
      shift 2
      ;;
    --mode)
      RUN_MODE="${2:-}"
      RUN_MODE_EXPLICIT=1
      shift 2
      ;;
    --quiet-log)
      STREAM_LOG="false"
      shift 1
      ;;
    --xsec-shards)
      XSEC_SHARDS="${2:-}"
      shift 2
      ;;
    --xsec-events)
      XSEC_RAW_EVENTS="${2:-}"
      shift 2
      ;;
    --xsec-memory)
      XSEC_MEMORY="${2:-}"
      shift 2
      ;;
    --xsec-base-seed)
      XSEC_BASE_SEED="${2:-}"
      shift 2
      ;;
    --manifest)
      MANIFEST_PATH="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if ! [[ "${NEVENTS}" =~ ^[0-9]+$ ]] || [[ "${NEVENTS}" -le 0 ]]; then
  echo "[ERROR] --events must be a positive integer; got '${NEVENTS}'" >&2
  exit 2
fi

if ! [[ "${TARGET_PASS}" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] --target-pass must be a non-negative integer; got '${TARGET_PASS}'" >&2
  exit 2
fi

if ! [[ "${MAX_RAW_EVENTS}" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] --max-raw-events must be a non-negative integer; got '${MAX_RAW_EVENTS}'" >&2
  exit 2
fi

if [[ "${TARGET_PASS}" -gt 0 && "${MAX_RAW_EVENTS}" -eq 0 ]]; then
  MAX_RAW_EVENTS=$(( TARGET_PASS * 100000 ))
fi

if ! [[ "${SEED}" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] --seed must be a non-negative integer; got '${SEED}'" >&2
  exit 2
fi

if [[ "${RUN_MODE}" != "compiled" && "${RUN_MODE}" != "interpreted" ]]; then
  echo "[ERROR] --mode must be 'compiled' or 'interpreted'; got '${RUN_MODE}'" >&2
  exit 2
fi

if [[ "${SAMPLE_FILTER}" != "all" && "${SAMPLE_FILTER}" != "PhotonJet12" && "${SAMPLE_FILTER}" != "PhotonJet20" ]]; then
  echo "[ERROR] --sample must be 'all', 'PhotonJet12', or 'PhotonJet20'; got '${SAMPLE_FILTER}'" >&2
  exit 2
fi

if ! [[ "${XSEC_SHARDS}" =~ ^[0-9]+$ ]] || [[ "${XSEC_SHARDS}" -le 0 ]]; then
  echo "[ERROR] --xsec-shards must be a positive integer; got '${XSEC_SHARDS}'" >&2
  exit 2
fi

if ! [[ "${XSEC_RAW_EVENTS}" =~ ^[0-9]+$ ]] || [[ "${XSEC_RAW_EVENTS}" -le 0 ]]; then
  echo "[ERROR] --xsec-events must be a positive integer; got '${XSEC_RAW_EVENTS}'" >&2
  exit 2
fi

if ! [[ "${XSEC_BASE_SEED}" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] --xsec-base-seed must be a non-negative integer; got '${XSEC_BASE_SEED}'" >&2
  exit 2
fi

setup_sphenix_env_for_worker() {
  export USER="${USER:-$(id -u -n)}"
  export LOGNAME="${LOGNAME:-${USER}}"
  export HOME="${HOME:-/sphenix/u/${LOGNAME}}"

  if [[ -r /opt/sphenix/core/bin/sphenix_setup.sh ]]; then
    set +u
    # shellcheck disable=SC1091
    source /opt/sphenix/core/bin/sphenix_setup.sh -n
    local myinstall="/sphenix/u/${USER}/thesisAnalysis/install"
    if [[ -d "${myinstall}" && -r /opt/sphenix/core/bin/setup_local.sh ]]; then
      # shellcheck disable=SC1091
      source /opt/sphenix/core/bin/setup_local.sh "${myinstall}"
    fi
    set -u
  fi
}

xsec_worker() {
  if [[ $# -ne 6 ]]; then
    echo "[ERROR] xsecWorker requires: sample shard raw_events seed workdir mode" >&2
    exit 2
  fi

  local sample="$1"
  local shard="$2"
  local raw_events="$3"
  local seed="$4"
  local workdir="$5"
  local mode="$6"

  if [[ "${sample}" != "PhotonJet12" && "${sample}" != "PhotonJet20" ]]; then
    echo "[ERROR] Invalid worker sample: ${sample}" >&2
    exit 2
  fi

  setup_sphenix_env_for_worker

  local shard_tag
  shard_tag="$(printf '%03d' "${shard}")"
  local tmp_out="${workdir}/tmp/${sample}_shard${shard_tag}"
  local result_dir="${workdir}/results"
  local log_dir="${workdir}/worker_logs"
  mkdir -p "${tmp_out}" "${result_dir}" "${log_dir}"

  echo "====================================================================="
  echo "Embedded photon xsec worker"
  echo "====================================================================="
  echo "Sample    : ${sample}"
  echo "Shard     : ${shard}"
  echo "Raw events: ${raw_events}"
  echo "Seed      : ${seed}"
  echo "Mode      : ${mode}"
  echo "Workdir   : ${workdir}"
  echo "Temp out  : ${tmp_out}"
  echo "Started   : $(date '+%Y-%m-%d %H:%M:%S')"
  echo "====================================================================="

  "${BASH:-bash}" "${SCRIPT_DIR}/estimateEmbeddedPhotonXsec.sh" run \
    --sample "${sample}" \
    --events "${raw_events}" \
    --seed "${seed}" \
    --outdir "${tmp_out}" \
    --mode "${mode}" \
    --quiet-log

  cp "${tmp_out}/results.csv" "${result_dir}/${sample}_shard${shard_tag}.csv"
  cp "${tmp_out}/summary.txt" "${log_dir}/${sample}_shard${shard_tag}.summary.txt"
  cp "${tmp_out}/run.log" "${log_dir}/${sample}_shard${shard_tag}.run.log"

  echo "[DONE] Worker result: ${result_dir}/${sample}_shard${shard_tag}.csv"
}

xsec_first_pass() {
  local timestamp
  timestamp="$(date '+%Y%m%d_%H%M%S')"
  local workdir="${BASE_DIR}/condor_snapshots/pythia_xsec_${timestamp}"
  local manifest="${BASE_DIR}/pythia_xsec_firstPass_${timestamp}.txt"
  local worker_mode="${RUN_MODE}"
  if [[ "${RUN_MODE_EXPLICIT}" -eq 0 ]]; then
    worker_mode="interpreted"
  fi

  mkdir -p "${workdir}/results" "${workdir}/worker_logs" "${workdir}/tmp"

  {
    echo "PYTHIA_XSEC_WORKDIR='${workdir}'"
    echo "PYTHIA_XSEC_SHARDS='${XSEC_SHARDS}'"
    echo "PYTHIA_XSEC_RAW_EVENTS='${XSEC_RAW_EVENTS}'"
    echo "PYTHIA_XSEC_BASE_SEED='${XSEC_BASE_SEED}'"
    echo "PYTHIA_XSEC_MODE='${worker_mode}'"
    echo "PYTHIA_XSEC_CREATED='${timestamp}'"
  } > "${manifest}"

  echo "====================================================================="
  echo "Embedded photon xsec firstPass"
  echo "====================================================================="
  echo "Mode        : $([[ "${LOCAL_MODE}" -eq 1 ]] && echo LOCAL || echo CONDOR)"
  echo "Workdir     : ${workdir}"
  echo "Manifest    : ${manifest}"
  echo "Shards      : ${XSEC_SHARDS} per sample"
  echo "Raw events  : ${XSEC_RAW_EVENTS} per shard"
  echo "Base seed   : ${XSEC_BASE_SEED}"
  echo "Worker mode : ${worker_mode}"
  echo "Memory      : ${XSEC_MEMORY}"
  echo "====================================================================="

  if [[ "${LOCAL_MODE}" -eq 1 ]]; then
    echo "[LOCAL] Running one local shard for PhotonJet12 and PhotonJet20."
    xsec_worker "PhotonJet12" 0 "${XSEC_RAW_EVENTS}" "${XSEC_BASE_SEED}" "${workdir}" "${worker_mode}"
    xsec_worker "PhotonJet20" 0 "${XSEC_RAW_EVENTS}" "$((XSEC_BASE_SEED + 1000000))" "${workdir}" "${worker_mode}"
    echo "[LOCAL DONE] Now aggregate this smoke test with:"
    echo "  ./scripts/estimateEmbeddedPhotonXsec.sh secondPass --manifest ${manifest}"
    return 0
  fi

  if ! command -v condor_submit >/dev/null 2>&1; then
    echo "[ERROR] condor_submit not found. Run this on SDCC with the Condor environment available." >&2
    exit 1
  fi

  local submit_dir="${BASE_DIR}/condor_sub"
  local log_dir="${BASE_DIR}/log"
  local stdout_dir="${BASE_DIR}/stdout"
  local error_dir="${BASE_DIR}/error"
  mkdir -p "${submit_dir}" "${log_dir}" "${stdout_dir}" "${error_dir}"

  local submit_file="${submit_dir}/pythia_xsec_${timestamp}.sub"
  {
    echo "universe = vanilla"
    echo "executable = ${SCRIPT_DIR}/estimateEmbeddedPhotonXsec.sh"
    echo "initialdir = ${BASE_DIR}"
    echo "getenv = True"
    echo "should_transfer_files = NO"
    echo "request_memory = ${XSEC_MEMORY}"
    echo "log = ${log_dir}/pythia_xsec.\$(Cluster).\$(Process).log"
    echo "output = ${stdout_dir}/pythia_xsec.\$(Cluster).\$(Process).out"
    echo "error = ${error_dir}/pythia_xsec.\$(Cluster).\$(Process).err"
    echo "stream_output = True"
    echo "stream_error = True"
    echo
    local shard
    for ((shard = 0; shard < XSEC_SHARDS; ++shard)); do
      echo "arguments = xsecWorker PhotonJet12 ${shard} ${XSEC_RAW_EVENTS} $((XSEC_BASE_SEED + shard)) ${workdir} ${worker_mode}"
      echo "queue"
    done
    for ((shard = 0; shard < XSEC_SHARDS; ++shard)); do
      echo "arguments = xsecWorker PhotonJet20 ${shard} ${XSEC_RAW_EVENTS} $((XSEC_BASE_SEED + 1000000 + shard)) ${workdir} ${worker_mode}"
      echo "queue"
    done
  } > "${submit_file}"

  echo "[INFO] Submit file: ${submit_file}"
  condor_submit "${submit_file}"
  echo "[DONE] Submitted $((2 * XSEC_SHARDS)) jobs."
  echo "[NEXT] After Condor finishes, aggregate with:"
  echo "  ./scripts/estimateEmbeddedPhotonXsec.sh secondPass --manifest ${manifest}"
}

xsec_second_pass() {
  if [[ -z "${MANIFEST_PATH}" ]]; then
    MANIFEST_PATH="$(ls -t "${BASE_DIR}"/pythia_xsec_firstPass_*.txt 2>/dev/null | head -n 1 || true)"
  fi
  if [[ -z "${MANIFEST_PATH}" || ! -r "${MANIFEST_PATH}" ]]; then
    echo "[ERROR] Could not find a readable firstPass manifest. Pass --manifest PATH." >&2
    exit 1
  fi

  # shellcheck disable=SC1090
  source "${MANIFEST_PATH}"
  local workdir="${PYTHIA_XSEC_WORKDIR:?manifest missing PYTHIA_XSEC_WORKDIR}"
  local combined="${workdir}/combined_results.csv"
  local summary_out="${workdir}/combined_summary.txt"

  if [[ ! -d "${workdir}/results" ]]; then
    echo "[ERROR] Missing result directory: ${workdir}/results" >&2
    exit 1
  fi

  local expected=$((2 * PYTHIA_XSEC_SHARDS))
  local found
  found="$(find "${workdir}/results" -maxdepth 1 -name 'PhotonJet*_shard*.csv' | wc -l | tr -d ' ')"

  echo "====================================================================="
  echo "Embedded photon xsec secondPass"
  echo "====================================================================="
  echo "Manifest : ${MANIFEST_PATH}"
  echo "Workdir  : ${workdir}"
  echo "Expected : ${expected} shard CSVs"
  echo "Found    : ${found} shard CSVs"
  echo "====================================================================="

  if [[ "${found}" -eq 0 ]]; then
    echo "[ERROR] No shard CSV files found under ${workdir}/results" >&2
    exit 1
  fi
  if [[ "${found}" -lt "${expected}" ]]; then
    echo "[WARN] Some shard CSV files are missing. Aggregating completed shards only." >&2
  fi

  {
    head -n 1 "$(find "${workdir}/results" -maxdepth 1 -name 'PhotonJet*_shard*.csv' | sort | head -n 1)"
    find "${workdir}/results" -maxdepth 1 -name 'PhotonJet*_shard*.csv' | sort | while read -r f; do
      tail -n +2 "${f}"
    done
  } > "${combined}"

  awk -F, '
    NR == 1 { next }
    {
      s = $1
      nOk[s] += $8
      nFail[s] += $9
      nPass[s] += $10
      sigmaWeighted[s] += $13 * $8
      sumW[s] += $24
      passW[s] += $25
      files[s] += 1
    }
    END {
      print "Embedded photon Pythia cross-section combined estimate"
      print "Generated: " strftime("%Y-%m-%d %H:%M:%S")
      print ""
      printf "%-12s %8s %14s %12s %14s %14s %14s %12s\n", "sample", "files", "n_ok", "n_pass", "eff_weight", "sigma_gen_mb", "sigma_eff_pb", "rel_stat"
      for (s in files) {
        sigmaGen = (nOk[s] > 0) ? sigmaWeighted[s] / nOk[s] : 0
        effW = (sumW[s] > 0) ? passW[s] / sumW[s] : 0
        sigmaEffPb = sigmaGen * effW * 1.0e9
        rel = (nPass[s] > 0) ? 1.0 / sqrt(nPass[s]) : -1
        printf "%-12s %8d %14.0f %12.0f %14.8e %14.8e %14.8e %12.6f\n", s, files[s], nOk[s], nPass[s], effW, sigmaGen, sigmaEffPb, rel
      }
      print ""
      print "Use sigma_eff_pb as the per-sample generated cross section in AnalyzeRecoilJets.h."
      print "Per-event merge weight is sigma_eff_pb divided by accepted/generated entries for that sample."
    }
  ' "${combined}" > "${summary_out}"

  cat "${summary_out}"
  echo "====================================================================="
  echo "[DONE] Combined CSV    : ${combined}"
  echo "[DONE] Combined summary: ${summary_out}"
  echo "====================================================================="
}

case "${ACTION}" in
  xsecWorker)
    xsec_worker "${WORKER_ARGS[@]}"
    exit 0
    ;;
  firstPass)
    xsec_first_pass
    exit 0
    ;;
  secondPass)
    xsec_second_pass
    exit 0
    ;;
  run)
    ;;
  *)
    echo "[ERROR] Unknown action: ${ACTION}" >&2
    usage >&2
    exit 2
    ;;
esac

if [[ -z "${CALIBRATIONROOT:-}" ]]; then
  echo "[ERROR] CALIBRATIONROOT is not set. Source the sPHENIX analysis environment first." >&2
  exit 1
fi

if ! command -v root >/dev/null 2>&1; then
  echo "[ERROR] root is not in PATH. Source the sPHENIX analysis environment first." >&2
  exit 1
fi

CFG12="${CALIBRATIONROOT}/Generators/JetStructure_TG/phpythia8_10GeV_JS_MDC2.cfg"
CFG20="${CALIBRATIONROOT}/Generators/JetStructure_TG/phpythia8_20GeV_JS_MDC2.cfg"

for cfg in "${CFG12}" "${CFG20}"; do
  if [[ ! -r "${cfg}" ]]; then
    echo "[ERROR] Cannot read Pythia config: ${cfg}" >&2
    exit 1
  fi
done

mkdir -p "${OUTDIR}"
RESULTS="${OUTDIR}/results.csv"
SUMMARY="${OUTDIR}/summary.txt"
LOG="${OUTDIR}/run.log"
MACRO_CACHE="${OUTDIR}/EstimateEmbeddedPhotonXsec.C"

ts() { date '+%Y-%m-%d %H:%M:%S'; }
say() { echo "[$(ts)] $*"; }

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/estimateEmbeddedPhotonXsec.XXXXXX")"
cleanup() {
  rm -rf "${TMPDIR}"
}
trap cleanup EXIT

MACRO="${MACRO_CACHE}"

cat > "${MACRO}" <<'EOF'
R__LOAD_LIBRARY(libpythia8)

#include <Pythia8/Pythia.h>

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  struct Sample
  {
    std::string name;
    std::string config;
    double photonPtLow = 0.0;
  };

  struct FilterBreakdown
  {
    long long eventsWithPhoton = 0;
    long long eventsWithPhotonEta = 0;
    long long eventsWithPhotonPt = 0;
    long long eventsWithPhotonParent = 0;
    long long photonCandidates = 0;
    long long photonEtaCandidates = 0;
    long long photonPtCandidates = 0;
    long long photonParentCandidates = 0;
    double maxPhotonPt = -1.0;
    double maxPhotonPtEta = -1.0;
    double maxPhotonPtEtaParent = -1.0;
  };

  bool PassProducerPhotonFilter(const Pythia8::Event& event, double ptLow, FilterBreakdown& b)
  {
    bool hasPhoton = false;
    bool hasPhotonEta = false;
    bool hasPhotonPt = false;
    bool hasPhotonParent = false;

    for (int i = 0; i < event.size(); ++i)
    {
      const auto& p = event[i];

      if (p.id() != 22)
      {
        continue;
      }

      hasPhoton = true;
      ++b.photonCandidates;
      if (p.pT() > b.maxPhotonPt)
      {
        b.maxPhotonPt = p.pT();
      }

      // Producer used SetStableParticleOnly(false), so do not require status > 0.
      if (p.eta() < -1.5 || p.eta() > 1.5)
      {
        continue;
      }
      hasPhotonEta = true;
      ++b.photonEtaCandidates;
      if (p.pT() > b.maxPhotonPtEta)
      {
        b.maxPhotonPtEta = p.pT();
      }

      if (p.pT() < ptLow)
      {
        continue;
      }
      hasPhotonPt = true;
      ++b.photonPtCandidates;

      const std::vector<int> moms = p.motherList();
      for (const int mom : moms)
      {
        if (mom < 0 || mom >= event.size())
        {
          continue;
        }
        const int absMotherId = std::abs(event[mom].id());
        if (absMotherId >= 1 && absMotherId <= 22)
        {
          hasPhotonParent = true;
          ++b.photonParentCandidates;
          if (p.pT() > b.maxPhotonPtEtaParent)
          {
            b.maxPhotonPtEtaParent = p.pT();
          }
          ++b.eventsWithPhoton;
          ++b.eventsWithPhotonEta;
          ++b.eventsWithPhotonPt;
          ++b.eventsWithPhotonParent;
          return true;
        }
      }
    }

    if (hasPhoton)
    {
      ++b.eventsWithPhoton;
    }
    if (hasPhotonEta)
    {
      ++b.eventsWithPhotonEta;
    }
    if (hasPhotonPt)
    {
      ++b.eventsWithPhotonPt;
    }
    if (hasPhotonParent)
    {
      ++b.eventsWithPhotonParent;
    }

    return false;
  }

  std::string CsvQuote(const std::string& s)
  {
    std::string out = "\"";
    for (const char c : s)
    {
      if (c == '"')
      {
        out += "\"\"";
      }
      else
      {
        out += c;
      }
    }
    out += "\"";
    return out;
  }
}

void EstimateEmbeddedPhotonXsec(long long nEvents = 1000000,
                                int seed = 75301,
                                const char* outCsv = "results.csv",
                                long long targetPass = 0,
                                long long maxRawEvents = 0,
                                const char* sampleFilterC = "all")
{
  const char* calib = std::getenv("CALIBRATIONROOT");
  if (!calib)
  {
    std::cerr << "[ERROR] CALIBRATIONROOT is not set." << std::endl;
    gSystem->Exit(1);
  }

  std::vector<Sample> samples = {
      {"PhotonJet12",
       std::string(calib) + "/Generators/JetStructure_TG/phpythia8_10GeV_JS_MDC2.cfg",
       12.0},
      {"PhotonJet20",
       std::string(calib) + "/Generators/JetStructure_TG/phpythia8_20GeV_JS_MDC2.cfg",
       20.0},
  };
  const std::string sampleFilter = sampleFilterC ? sampleFilterC : "all";

  std::ofstream csv(outCsv);
  if (!csv)
  {
    std::cerr << "[ERROR] Cannot write " << outCsv << std::endl;
    gSystem->Exit(1);
  }

  csv << "sample,config,photon_pt_low,n_requested,target_pass,max_raw_events,stop_reason,"
         "n_next_ok,n_next_failed,n_filter_pass,"
         "filter_eff_count,filter_eff_weight,sigma_gen_mb,sigma_eff_mb,sigma_eff_pb,"
         "approx_filter_rel_stat,events_with_photon,events_with_photon_eta,"
         "events_with_photon_pt,events_with_photon_parent,max_photon_pt,"
         "max_photon_pt_eta,max_photon_pt_eta_parent,total_weight_sum,"
         "trigger_pass_weight_sum\n";

  std::cout << std::setprecision(10);
  std::cout << "[INFO] Running generator-only estimate with nEvents=" << nEvents
            << " raw attempts per sample, base seed=" << seed << std::endl;
  if (targetPass > 0)
  {
    std::cout << "[INFO] Target-pass mode enabled: run each sample until nPass >= "
              << targetPass << " or raw attempts reach " << maxRawEvents << std::endl;
  }
  else
  {
    std::cout << "[INFO] Fixed-raw mode enabled: run exactly nEvents raw Pythia events per sample." << std::endl;
  }
  std::cout << "[INFO] Stages per sample: configure Pythia -> init -> event loop -> sigma_eff calculation" << std::endl;
  std::cout << "[INFO] Sample filter: " << sampleFilter << std::endl;

  for (std::size_t is = 0; is < samples.size(); ++is)
  {
    const auto& s = samples[is];
    if (sampleFilter != "all" && sampleFilter != s.name)
    {
      continue;
    }

    std::cout << "\n[INFO] Sample " << s.name << std::endl;
    std::cout << "[INFO]   config        : " << s.config << std::endl;
    std::cout << "[INFO]   photon pT cut : " << s.photonPtLow << " GeV" << std::endl;
    std::cout << "[INFO]   filter        : id==22, stable-only disabled, -1.5<eta<1.5, mother |id| in [1,22]" << std::endl;

    std::string xmlDir;
    if (const char* pythia8data = std::getenv("PYTHIA8DATA"))
    {
      xmlDir = pythia8data;
    }
    else if (const char* pythia8 = std::getenv("PYTHIA8"))
    {
      xmlDir = std::string(pythia8) + "/xmldoc";
    }

    std::cout << "[INFO]   XML dir       : " << (xmlDir.empty() ? std::string("(Pythia default)") : xmlDir) << std::endl;
    std::cout << "[INFO]   configuring Pythia..." << std::endl;

    Pythia8::Pythia pythia(xmlDir, false);
    pythia.readFile(s.config);
    pythia.readString("Random:setSeed = on");
    const int sampleSeed = seed + static_cast<int>(1000 * is);
    pythia.readString("Random:seed = " + std::to_string(sampleSeed));
    pythia.readString("Print:quiet = on");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

    std::cout << "[INFO]   seed          : " << sampleSeed << std::endl;
    std::cout << "[INFO]   calling pythia.init()..." << std::endl;
    if (!pythia.init())
    {
      std::cerr << "[ERROR] Pythia init failed for " << s.name << std::endl;
      gSystem->Exit(1);
    }
    std::cout << "[INFO]   pythia.init() complete." << std::endl;
    std::cout << "[INFO]   starting event loop." << std::endl;

    long long nOk = 0;
    long long nFailed = 0;
    long long nPass = 0;
    double sumW = 0.0;
    double sumW2 = 0.0;
    double passW = 0.0;
    double passW2 = 0.0;
    FilterBreakdown breakdown;

    const long long rawLimit = (targetPass > 0) ? maxRawEvents : nEvents;
    const long long reportEveryRaw = std::max(1LL, rawLimit / 20LL);
    const long long reportEveryPass = (targetPass > 0) ? std::max(1LL, targetPass / 20LL) : 0LL;
    long long nextPassReport = reportEveryPass;
    std::string stopReason = "raw_limit";

    for (long long iev = 0; iev < rawLimit; ++iev)
    {
      if (!pythia.next())
      {
        ++nFailed;
        continue;
      }

      ++nOk;
      const double w = pythia.info.weight();
      sumW += w;
      sumW2 += w * w;

      if (PassProducerPhotonFilter(pythia.event, s.photonPtLow, breakdown))
      {
        ++nPass;
        passW += w;
        passW2 += w * w;
      }

      if (targetPass > 0 && nPass >= nextPassReport)
      {
        std::cout << "[INFO]   " << s.name << " pass-progress:"
                  << " nPass=" << nPass << " / " << targetPass
                  << " | rawAttempts=" << (iev + 1)
                  << " nOk=" << nOk
                  << " effSoFar=" << ((nOk > 0) ? static_cast<double>(nPass) / static_cast<double>(nOk) : 0.0)
                  << std::endl;
        while (nextPassReport <= nPass)
        {
          nextPassReport += reportEveryPass;
        }
      }

      if ((iev + 1) % reportEveryRaw == 0 || iev + 1 == rawLimit)
      {
        std::cout << "[INFO]   " << s.name << " progress: " << (iev + 1)
                  << " / " << rawLimit
                  << " | nOk=" << nOk
                  << " nPass=" << nPass
                  << " effSoFar=" << ((nOk > 0) ? static_cast<double>(nPass) / static_cast<double>(nOk) : 0.0)
                  << std::endl;
      }

      if (targetPass > 0 && nPass >= targetPass)
      {
        stopReason = "target_pass";
        break;
      }
    }

    const double sigmaGenMb = pythia.info.sigmaGen();
    const double effCount = (nOk > 0) ? static_cast<double>(nPass) / static_cast<double>(nOk) : 0.0;
    const double effWeight = (sumW > 0.0) ? passW / sumW : 0.0;
    const double sigmaEffMb = sigmaGenMb * effWeight;
    const double sigmaEffPb = sigmaEffMb * 1.0e9;
    const double relStat = (nPass > 0) ? 1.0 / std::sqrt(static_cast<double>(nPass)) : -1.0;

    std::cout << "[INFO]   event loop complete; computing effective cross section." << std::endl;
    std::cout << "[INFO]   stop reason: " << stopReason
              << " | targetPass=" << targetPass
              << " | rawLimit=" << rawLimit << std::endl;
    std::cout << "[INFO]   filter breakdown:"
              << " eventsWithPhoton=" << breakdown.eventsWithPhoton
              << " eventsWithPhotonEta=" << breakdown.eventsWithPhotonEta
              << " eventsWithPhotonPt=" << breakdown.eventsWithPhotonPt
              << " eventsWithPhotonParent=" << breakdown.eventsWithPhotonParent
              << " photonCandidates=" << breakdown.photonCandidates
              << " photonEtaCandidates=" << breakdown.photonEtaCandidates
              << " photonPtCandidates=" << breakdown.photonPtCandidates
              << " photonParentCandidates=" << breakdown.photonParentCandidates
              << " maxPhotonPt=" << breakdown.maxPhotonPt
              << " maxPhotonPtEta=" << breakdown.maxPhotonPtEta
              << " maxPhotonPtEtaParent=" << breakdown.maxPhotonPtEtaParent
              << std::endl;
    std::cout << std::scientific << std::setprecision(8);
    std::cout << "[RESULT] " << s.name
              << " nOk=" << nOk
              << " nPass=" << nPass
              << " effCount=" << effCount
              << " effWeight=" << effWeight
              << " sigmaGen=" << sigmaGenMb << " mb"
              << " sigmaEff=" << sigmaEffMb << " mb"
              << " = " << sigmaEffPb << " pb"
              << " relStat~" << relStat
              << std::endl;
    std::cout << std::defaultfloat << std::setprecision(10);

    csv << s.name << ","
        << CsvQuote(s.config) << ","
        << s.photonPtLow << ","
        << nEvents << ","
        << targetPass << ","
        << rawLimit << ","
        << stopReason << ","
        << nOk << ","
        << nFailed << ","
        << nPass << ","
        << std::setprecision(12)
        << effCount << ","
        << effWeight << ","
        << sigmaGenMb << ","
        << sigmaEffMb << ","
        << sigmaEffPb << ","
        << relStat << ","
        << breakdown.eventsWithPhoton << ","
        << breakdown.eventsWithPhotonEta << ","
        << breakdown.eventsWithPhotonPt << ","
        << breakdown.eventsWithPhotonParent << ","
        << breakdown.maxPhotonPt << ","
        << breakdown.maxPhotonPtEta << ","
        << breakdown.maxPhotonPtEtaParent << ","
        << sumW << ","
        << passW << "\n";

    std::cout << "[INFO]   Pythia statistics block follows." << std::endl;
    pythia.stat();
    std::cout << "[INFO]   Done with " << s.name << "." << std::endl;
  }

  csv.close();
  std::cout << "\n[INFO] Wrote " << outCsv << std::endl;
}
EOF

echo "====================================================================="
echo "Embedded photon Pythia cross-section estimate"
echo "====================================================================="
echo "Base dir : ${BASE_DIR}"
echo "Out dir  : ${OUTDIR}"
echo "Events   : ${NEVENTS} raw attempts per sample"
echo "Target   : ${TARGET_PASS} filter-passing events per sample (0 = fixed raw attempts)"
echo "Max raw  : ${MAX_RAW_EVENTS} raw attempts per sample in target-pass mode"
echo "Seed     : ${SEED}"
echo "Mode     : ${RUN_MODE}"
echo "Sample   : ${SAMPLE_FILTER}"
echo "Config12 : ${CFG12}"
echo "Config20 : ${CFG20}"
echo "Macro    : ${MACRO}"
echo "Log      : ${LOG}"
echo "====================================================================="
echo "What will happen next:"
echo "  [1] ROOT starts."
if [[ "${RUN_MODE}" == "compiled" ]]; then
echo "  [2] ACLiC compiles ${MACRO##*/}. This can take a few minutes on first use."
echo "      Compiler warnings from Pythia headers are usually harmless."
else
echo "  [2] ROOT interprets ${MACRO##*/}. This avoids compile time but is slower for large samples."
fi
if [[ "${TARGET_PASS}" -gt 0 ]]; then
echo "  [3] Pythia initializes PhotonJet12, then runs until ${TARGET_PASS} filter passes or ${MAX_RAW_EVENTS} raw attempts."
echo "  [4] Pythia initializes PhotonJet20, then runs until ${TARGET_PASS} filter passes or ${MAX_RAW_EVENTS} raw attempts."
else
echo "  [3] Pythia initializes PhotonJet12, then runs ${NEVENTS} raw events."
echo "  [4] Pythia initializes PhotonJet20, then runs ${NEVENTS} raw events."
fi
echo "  [5] The script writes results.csv and summary.txt."
echo "====================================================================="

rm -f "${LOG}"
touch "${LOG}"

if [[ "${RUN_MODE}" == "compiled" ]]; then
  ROOT_EXPR="${MACRO}+(${NEVENTS}, ${SEED}, \"${RESULTS}\", ${TARGET_PASS}, ${MAX_RAW_EVENTS}, \"${SAMPLE_FILTER}\")"
else
  ROOT_EXPR="${MACRO}(${NEVENTS}, ${SEED}, \"${RESULTS}\", ${TARGET_PASS}, ${MAX_RAW_EVENTS}, \"${SAMPLE_FILTER}\")"
fi

say "[ROOT] Command: root -l -b -q '${ROOT_EXPR}'"
say "[ROOT] Streaming log below. A long compile is normal before the first [INFO] from the macro."
echo "---------------------------------------------------------------------"

root -l -b -q "${ROOT_EXPR}" > "${LOG}" 2>&1 &
ROOT_PID=$!

TAIL_PID=""
if [[ "${STREAM_LOG}" == "true" ]]; then
  tail -n +1 -f "${LOG}" &
  TAIL_PID=$!
fi

POLL_SECONDS=2
HEARTBEAT_SECONDS=30
SECONDS_WAITED=0
while kill -0 "${ROOT_PID}" 2>/dev/null; do
  sleep "${POLL_SECONDS}"
  SECONDS_WAITED=$((SECONDS_WAITED + POLL_SECONDS))
  if (( SECONDS_WAITED % HEARTBEAT_SECONDS == 0 )); then
    say "[heartbeat] ROOT/Pythia still running after ${SECONDS_WAITED}s; log: ${LOG}"
  fi
done

set +e
wait "${ROOT_PID}"
ROOT_STATUS=$?
set -e

if [[ -n "${TAIL_PID}" ]]; then
  sleep 1
  kill "${TAIL_PID}" 2>/dev/null || true
  wait "${TAIL_PID}" 2>/dev/null || true
fi

echo "---------------------------------------------------------------------"
say "[ROOT] Finished with exit code ${ROOT_STATUS}"

if [[ "${ROOT_STATUS}" -ne 0 ]]; then
  echo "[ERROR] ROOT/Pythia estimate failed. Full log:" >&2
  echo "  ${LOG}" >&2
  tail -80 "${LOG}" >&2 || true
  exit 1
fi

{
  echo "Embedded photon Pythia cross-section estimate"
  echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "Events per sample: ${NEVENTS}"
  echo "Target filter passes per sample: ${TARGET_PASS}"
  echo "Max raw attempts per sample: ${MAX_RAW_EVENTS}"
  echo "Seed: ${SEED}"
  echo "Mode: ${RUN_MODE}"
  echo "Sample: ${SAMPLE_FILTER}"
  echo
  echo "Producer mapping:"
  echo "  PhotonJet12 -> ${CFG12} + photon filter pT > 12 GeV, |eta| < 1.5"
  echo "  PhotonJet20 -> ${CFG20} + photon filter pT > 20 GeV, |eta| < 1.5"
  echo
  echo "Results CSV:"
  cat "${RESULTS}"
  echo
  echo "Use sigma_eff_pb / accepted_events as the per-sample merge weight."
  echo "Approximate filter-only relative stat uncertainty is approx_filter_rel_stat."
  echo
  echo "Files:"
  echo "  results: ${RESULTS}"
  echo "  log    : ${LOG}"
} > "${SUMMARY}"

cat "${SUMMARY}"

echo "====================================================================="
echo "[DONE] Summary: ${SUMMARY}"
echo "[DONE] Results: ${RESULTS}"
echo "[DONE] Log    : ${LOG}"
echo "====================================================================="
