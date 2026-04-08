#!/usr/bin/env bash
###############################################################################
# RecoilJets_Condor_submit.sh — one-stop driver for LOCAL testing, LOCAL
# IsolationAudit runs, and CONDOR submission
#
# QUICK START
#   isPP:
#     ./RecoilJets_Condor_submit.sh isPP local
#     ./RecoilJets_Condor_submit.sh isPP condor all
#
#   isPPrun25:
#     ./RecoilJets_Condor_submit.sh isPPrun25 local
#     ./RecoilJets_Condor_submit.sh isPPrun25 condor all
#
#   isAuAu:
#     ./RecoilJets_Condor_submit.sh isAuAu local
#     ./RecoilJets_Condor_submit.sh isAuAu isLocalIsoPing
#     RJ_CLUSTER_UEPIPELINE=noSub       ./RecoilJets_Condor_submit.sh isAuAu isLocalIsoPing
#     RJ_CLUSTER_UEPIPELINE=baseVariant ./RecoilJets_Condor_submit.sh isAuAu isLocalIsoPing
#     RJ_CLUSTER_UEPIPELINE=variantA    ./RecoilJets_Condor_submit.sh isAuAu isLocalIsoPing
#     RJ_CLUSTER_UEPIPELINE=variantB    ./RecoilJets_Condor_submit.sh isAuAu isLocalIsoPing
#     ./RecoilJets_Condor_submit.sh isAuAu condor all
#
#   isOO:
#     ./RecoilJets_Condor_submit.sh isOO local
#     ./RecoilJets_Condor_submit.sh isOO condor all
#
#   isSim:
#     ./RecoilJets_Condor_submit.sh isSim local
#     ./RecoilJets_Condor_submit.sh isSim condorDoAll
#
#   isSimJet5:
#     ./RecoilJets_Condor_submit.sh isSimJet5 local
#     ./RecoilJets_Condor_submit.sh isSimJet5 condorDoAll
#
#   isSimMB:
#     ./RecoilJets_Condor_submit.sh isSimMB local
#     ./RecoilJets_Condor_submit.sh isSimMB condorDoAll
#
#   isSimEmbedded:
#     ./RecoilJets_Condor_submit.sh isSimEmbedded local
#     ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll
#
# OVERVIEW
#   • This script is the top-level submission driver for:
#       isPP, isPPrun25, isAuAu, isOO, isSim, isSimJet5, isSimMB, isSimEmbedded
#   • DATA modes use per-run runlists and can do:
#       local, isLocalIsoPing, CHECKJOBS, splitGoldenRunList,
#       condor testJob, condor round, condor all
#   • SIM modes use staged matched master lists and can do:
#       local, CHECKJOBS, condorTest, condorDoAll
#   • Jobs never mix files across runs in DATA mode.
#   • SIM condorDoAll/local sweep over the YAML matrix:
#       jet_pt_min × back_to_back_dphi_min_pi_fraction × vz_cut_cm × coneR × iso mode
#       and, for AuAu-like tags, also clusterUEpipeline.
#   • DATA condor/local matrix sweeps use the same scalar-YAML override strategy.
#
# NEW: AUAU LOCAL ISOLATION AUDIT MODE
#   • Action token:
#       isLocalIsoPing
#   • Valid only for:
#       isAuAu
#   • Purpose:
#       Run the LIVE Au+Au photon-side path locally, using the existing
#       authoritative PhotonClusterBuilder -> RecoilJets isolation chain, and
#       print a detailed terminal IsolationAudit summary without introducing a
#       second isolation calculator.
#   • Selection policy:
#       - local only
#       - picks the first sorted Au+Au run that has BOTH:
#           1) a non-empty per-run DST_CALOFITTING list
#           2) active photon_10_plus_MBD_NS_geq_2_vtx_lt_150 trigger path
#              (checked via gl1_scaledown scaledown22 != -1)
#       - builds a combined local input list from grouped chunk files for that run
#       - passes env-controlled audit metadata into the live wrapper/macro path
#   • Variant control:
#       Select the photon-side UE variant with RJ_CLUSTER_UEPIPELINE:
#         noSub | baseVariant | variantA | variantB
#   • Default target:
#       RJ_ISO_AUDIT_TARGET_PER_CENT defaults to 200 inclusive photons per
#       centrality bin.
#
# DATASET TOKENS (case-insensitive)
#   DATA:
#     isPP        | pp         | PP
#     isPPrun25   | pprun25    | pp25 | PP25
#     isAuAu      | auau       | AA
#     isOO        | oo         | OO
#
#   SIM:
#     isSim         | sim         | SIM
#     isSimJet5     | simjet5     | SIMJET5
#     isSimMB       | simmb       | SIMMB
#     isSimEmbedded | simembedded | SIMEMBEDDED
#
# DATASET MATRIX
#   ┌────────────────┬─────────────────────────┬──────────────────────────────┬─────────────────────────────┐
#   │ dataset        │ wrapper                 │ macro                        │ standard bulk action        │
#   ├────────────────┼─────────────────────────┼──────────────────────────────┼─────────────────────────────┤
#   │ isPP           │ RecoilJets_Condor.sh    │ Fun4All_recoilJets.C         │ condor all                  │
#   │ isPPrun25      │ RecoilJets_Condor.sh    │ Fun4All_recoilJets.C         │ condor all                  │
#   │ isAuAu         │ RecoilJets_Condor_AuAu.sh│ Fun4All_recoilJets_AuAu.C   │ condor all                  │
#   │ isOO           │ RecoilJets_Condor_AuAu.sh│ Fun4All_recoilJets_AuAu.C   │ condor all                  │
#   │ isSim          │ RecoilJets_Condor.sh    │ Fun4All_recoilJets.C         │ condorDoAll                 │
#   │ isSimJet5      │ RecoilJets_Condor.sh    │ Fun4All_recoilJets.C         │ condorDoAll                 │
#   │ isSimMB        │ RecoilJets_Condor.sh    │ Fun4All_recoilJets.C         │ condorDoAll                 │
#   │ isSimEmbedded  │ RecoilJets_Condor_AuAu.sh│ Fun4All_recoilJets_AuAu.C   │ condorDoAll                 │
#   └────────────────┴─────────────────────────┴──────────────────────────────┴─────────────────────────────┘
#
# LOCAL / BULK COPY-PASTE MATRIX
#   DATA:
#     ./RecoilJets_Condor_submit.sh isPP local
#     ./RecoilJets_Condor_submit.sh isPP condor all
#
#     ./RecoilJets_Condor_submit.sh isPPrun25 local
#     ./RecoilJets_Condor_submit.sh isPPrun25 condor all
#
#     ./RecoilJets_Condor_submit.sh isAuAu local
#     ./RecoilJets_Condor_submit.sh isAuAu isLocalIsoPing
#     ./RecoilJets_Condor_submit.sh isAuAu condor all
#
#     ./RecoilJets_Condor_submit.sh isOO local
#     ./RecoilJets_Condor_submit.sh isOO condor all
#
#   SIM:
#     ./RecoilJets_Condor_submit.sh isSim local
#     ./RecoilJets_Condor_submit.sh isSim condorDoAll
#
#     ./RecoilJets_Condor_submit.sh isSimJet5 local
#     ./RecoilJets_Condor_submit.sh isSimJet5 condorDoAll
#
#     ./RecoilJets_Condor_submit.sh isSimMB local
#     ./RecoilJets_Condor_submit.sh isSimMB condorDoAll
#
#     ./RecoilJets_Condor_submit.sh isSimEmbedded local
#     ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll
#
# INPUT CONTRACTS
#   DATA:
#     • Golden run list + per-run list files.
#     • Per-run list file contains one ROOT file path per line.
#     • isLocalIsoPing reuses the existing Au+Au per-run list infrastructure and
#       concatenates grouped local chunks into one combined local input list.
#
#   SIM:
#     • Staged matched lists under ${BASE}/simListFiles/<sample>/.
#     • sim_init() builds a 5-column master list:
#         col1 = DST_CALO_CLUSTER
#         col2 = G4Hits
#         col3 = DST_JETS
#         col4 = DST_GLOBAL
#         col5 = DST_MBD_EPD
#     • isSimEmbedded currently still uses the same 5-column staged contract.
#
# DIRECTORY LAYOUT (fixed)
#   BASE          = /sphenix/u/patsfan753/scratch/thesisAnalysis
#   LOG_DIR       = ${BASE}/log
#   OUT_DIR       = ${BASE}/stdout
#   ERR_DIR       = ${BASE}/error
#   SUB_DIR       = ${BASE}/condor_sub
#   STAGE_ROOT    = ${BASE}/condor_lists
#   ROUND_ROOT    = ${BASE}/condor_segments
#
# DATA INPUTS
#   Golden run lists:
#     • PP      : ${BASE}/GRLs_tanner/run2pp_ana509_2024p022_v001_dst_calofitting_grl.list
#     • PP run25: ${BASE}/GRLs_tanner/run3pp_new_newcdbtag_v008_dst_calofitting_grl.list
#     • AuAu    : ${BASE}/GRLs_tanner/run3auau_new_newcdbtag_v008_dst_calofitting_grl.list
#     • OO      : ${BASE}/GRLs_tanner/run3oo_ana536_2025p010_v001_dst_calofitting_grl.list
#
#   Per-run list directories:
#     • PP      : ${BASE}/dst_lists_pp
#     • PP run25: ${BASE}/dst_lists_pp_run25
#     • AuAu    : ${BASE}/dst_lists_auau
#     • OO      : ${BASE}/dst_lists_oo
#
# OUTPUT ROOT DESTINATIONS
#   DATA:
#     • isPP        → /sphenix/tg/tg01/bulk/jbennett/thesisAna/pp
#     • isPPrun25   → /sphenix/tg/tg01/bulk/jbennett/thesisAna/pp25
#     • isAuAu      → /sphenix/tg/tg01/bulk/jbennett/thesisAna/auau
#     • isOO        → /sphenix/tg/tg01/bulk/jbennett/thesisAna/oo
#     • isLocalIsoPing local outputs:
#         ${BASE}/iso_ping_local/auau/<cfg_tag>/
#
#   SIM:
#     • isSim         → /sphenix/tg/tg01/bulk/jbennett/thesisAna/sim
#     • isSimJet5     → /sphenix/tg/tg01/bulk/jbennett/thesisAna/simjet5
#     • isSimMB       → /sphenix/tg/tg01/bulk/jbennett/thesisAna/simmb
#     • isSimEmbedded → /sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded
#
# VERBOSITY POLICY
#   • LOCAL mode:
#       RJ_VERBOSITY defaults to 10, unless overridden by VERBOSE=N.
#   • isLocalIsoPing:
#       RJ_VERBOSITY defaults to 1 unless overridden by VERBOSE=N.
#   • CONDOR mode:
#       RJ_VERBOSITY is forced to 0 in submit files.
#   • Useful debug knobs:
#       RJ_GROUP_TRACE=1
#       RJ_SUBMIT_TRACE=1
#       RJ_SUBMIT_PROGRESS_EVERY=N
#
# COMMAND SUMMARY
#   DATA:
#     $ ./RecoilJets_Condor_submit.sh <dataset> local [Nevents] [VERBOSE=N]
#     $ ./RecoilJets_Condor_submit.sh isAuAu isLocalIsoPing [VERBOSE=N]
#     $ ./RecoilJets_Condor_submit.sh <dataset> CHECKJOBS [groupSize N]
#     $ ./RecoilJets_Condor_submit.sh <dataset> splitGoldenRunList [groupSize N] [maxJobs M]
#     $ ./RecoilJets_Condor_submit.sh <dataset> condor testJob
#     $ ./RecoilJets_Condor_submit.sh <dataset> condor round <N> [groupSize N] [firstChunk]
#     $ ./RecoilJets_Condor_submit.sh <dataset> condor all [groupSize N]
#
#   SIM:
#     $ ./RecoilJets_Condor_submit.sh <simDataset> local [Nevents] [VERBOSE=N] [SAMPLE=...]
#     $ ./RecoilJets_Condor_submit.sh <simDataset> CHECKJOBS [groupSize N] [SAMPLE=...]
#     $ ./RecoilJets_Condor_submit.sh <simDataset> condorTest [SAMPLE=...]
#     $ ./RecoilJets_Condor_submit.sh <simDataset> condorDoAll [groupSize N] [SAMPLE=...]
#
# ISLOCALISOPING ENVIRONMENT / CONTROLS
#   • Required dataset:
#       isAuAu
#   • Variant selection:
#       RJ_CLUSTER_UEPIPELINE=noSub|baseVariant|variantA|variantB
#   • Audit target:
#       RJ_ISO_AUDIT_TARGET_PER_CENT=<positive integer>
#   • Internally exported metadata for the live RecoilJets_AuAu audit mode:
#       RJ_ISO_AUDIT_MODE=1
#       RJ_ISO_AUDIT_SELECTED_RUN
#       RJ_ISO_AUDIT_CHUNK_SPAN
#       RJ_ISO_AUDIT_GROUP_SIZE
#       RJ_ISO_AUDIT_COMBINED_LIST
#       RJ_ISO_AUDIT_CLUSTER_UEPIPELINE
#       RJ_ISO_AUDIT_PATH_CLASS
#       RJ_ISO_AUDIT_PHOTON_INPUT_CLUSTER_NODE
#       RJ_ISO_AUDIT_PHOTON_BUILDER_IS_AUAU
#       RJ_ISO_AUDIT_PCB_TOWER_PREFIX
#       RJ_ISO_AUDIT_PCB_EM_NODE
#       RJ_ISO_AUDIT_PCB_HI_NODE
#       RJ_ISO_AUDIT_PCB_HO_NODE
#
# CLEANING BEHAVIOR
#   • make_groups():
#       removes previous run<run8>_grp*.list for that run before regenerating.
#   • splitGoldenRunList:
#       removes previous goldenRuns_<tag>_segment*.txt before rebuilding rounds.
#   • submit_condor() in DATA mode:
#       wipes the dataset output tree under DEST_BASE before submission.
#   • condorDoAll() in SIM mode:
#       wipes previous staged group files, outputs, and tagged logs before resubmitting.
#   • isLocalIsoPing:
#       creates a combined local chunk list but does not alter condor segmentation.
#
# REQUIREMENTS
#   • condor_submit must be on PATH for Condor submissions.
#   • psql must be on PATH if using TRIGGER=<bit> or isLocalIsoPing.
#   • Required golden lists, per-run lists, and staged sim lists must exist.
###############################################################################
set -euo pipefail

# ------------------------ Pretty printing ------------------
BOLD=$'\e[1m'; DIM=$'\e[2m'; RED=$'\e[31m'; YEL=$'\e[33m'; GRN=$'\e[32m'; BLU=$'\e[34m'; RST=$'\e[0m'
say()  { printf "${BLU}➜${RST} %s\n" "$*"; }
warn() { printf "${YEL}⚠ %s${RST}\n" "$*" >&2; }
err()  { printf "${RED}✘ %s${RST}\n" "$*" >&2; }

# ------------------------ Fixed paths ----------------------
BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
# NOTE: MACRO/EXE are dataset-specific and are set in resolve_dataset().
# These defaults are placeholders only.
MACRO=""
EXE=""

# Logs (as requested)
LOG_DIR="${BASE}/log"
OUT_DIR="${BASE}/stdout"
ERR_DIR="${BASE}/error"
SUB_DIR="${BASE}/condor_sub"

# Golden lists provided by you
PP_GOLDEN="${BASE}/GRLs_tanner/run2pp_ana509_2024p022_v001_dst_calofitting_grl.list"
PP25_GOLDEN="${BASE}/GRLs_tanner/run3pp_new_newcdbtag_v008_dst_calofitting_grl.list"
AA_GOLDEN="${BASE}/GRLs_tanner/run3auau_new_newcdbtag_v008_dst_calofitting_grl.list"
OO_GOLDEN="${BASE}/GRLs_tanner/run3oo_ana536_2025p010_v001_dst_calofitting_grl.list"

# Per-run input list directories
PP_LIST_DIR="${BASE}/dst_lists_pp"
PP25_LIST_DIR="${BASE}/dst_lists_pp_run25"
AA_LIST_DIR="${BASE}/dst_lists_auau"
OO_LIST_DIR="${BASE}/dst_lists_oo"

# Where we stage grouped chunk .list files and the round files
STAGE_ROOT="${BASE}/condor_lists"
ROUND_ROOT="${BASE}/condor_segments"

# Dataset-specific output trees on TG storage
PP_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp"
PP25_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/pp25"
AA_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/auau"
OO_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/oo"

# ------------------------ Defaults -------------------------
GROUP_SIZE=7         # files per Condor job (never mixes runs)
MAX_JOBS=50000      # job budget per round file
LOCAL_EVENTS=5000   # default N for "local" if not given
TRIGGER_BIT=""      # optional: filter runs by GL1 scaledown bit (e.g., TRIGGER=26)

# ------------------------ Simulation (MC) defaults --------
# isSim expects a master list with ONE line per segment, whitespace-separated:
#   col1 = DST_CALO_CLUSTER
#   col2 = G4Hits
#   col3 = DST_JETS        (truth jets DST)
#   col4 = DST_GLOBAL      (GlobalVertexMap lives here)
#   col5 = DST_MBD_EPD     (MBD inputs; needed for reco MBD vertex)
SIM_ROOT="${BASE}/simListFiles"
SIM_SAMPLE_DEFAULT="run28_photonjet10"
SIM_SAMPLE="${SIM_SAMPLE_DEFAULT}"

# Name for the staged 5-column master list (built by sim_init from matched lists)
SIM_MASTER5_NAME="DST_SIM_MASTER_5COL.list"

# Output dir root for sim is: ${SIM_DEST_BASE}/<cfgTag>/<SIM_SAMPLE>
SIM_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/sim"
SIMJET5_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simjet5"
SIMMB_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simmb"
SIMEMBED_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded"

# Flag: set to 1 for any isSim variant (isSim, isSimJet5, isSimMB, isSimEmbedded)
IS_SIM=0

# YAML config (Fun4All default is next to the macro unless RJ_CONFIG_YAML is set)
SIM_YAML_DEFAULT="${BASE}/macros/analysis_config.yaml"
SIM_YAML_OVERRIDE_DIR="${BASE}/condor_yaml_overrides"
SIM_CFG_TAG=""

# Frozen bulk-submission runtime snapshots.
# IMPORTANT:
#   - used ONLY for DATA "condor all" and SIM "condorDoAll"
#   - each bulk submission gets a unique timestamped snapshot directory
#   - do NOT recycle one shared snapshot path, or idle old jobs can break
SNAPSHOT_ROOT="${BASE}/condor_snapshots"
BULK_FROZEN_EXE=""
BULK_FROZEN_MACRO=""

create_pipeline_snapshot() {
  local mode="$1"   # pp | auau
  local stamp="$2"

  local snap_dir="${SNAPSHOT_ROOT}/${TAG}_${stamp}"
  local snap_lib_dir="${snap_dir}/lib"
  local user_root="/sphenix/u/${USER:-$(id -u -n)}"

  local live_wrapper=""
  local live_macro=""
  local snap_wrapper=""
  local snap_macro=""

  local snap_impl="${snap_dir}/Fun4All_recoilJets_unified_impl.C"
  local snap_calo="${snap_dir}/Calo_Calib.C"
  local snap_pp_header="${snap_dir}/RecoilJets.h"
  local snap_auau_header="${snap_dir}/RecoilJets_AuAu.h"

  mkdir -p "$snap_dir" "$snap_lib_dir"

  if [[ "$mode" == "auau" ]]; then
    live_wrapper="${BASE}/RecoilJets_Condor_AuAu.sh"
    live_macro="${BASE}/macros/Fun4All_recoilJets_AuAu.C"
    snap_wrapper="${snap_dir}/RecoilJets_Condor_AuAu.sh"
    snap_macro="${snap_dir}/Fun4All_recoilJets_AuAu.C"
  else
    live_wrapper="${BASE}/RecoilJets_Condor.sh"
    live_macro="${BASE}/macros/Fun4All_recoilJets.C"
    snap_wrapper="${snap_dir}/RecoilJets_Condor.sh"
    snap_macro="${snap_dir}/Fun4All_recoilJets.C"
  fi

  cp -f "$live_wrapper" "$snap_wrapper"
  cp -f "$live_macro" "$snap_macro"
  cp -f "${BASE}/macros/Fun4All_recoilJets_unified_impl.C" "$snap_impl"
  cp -f "${BASE}/macros/Calo_Calib.C" "$snap_calo"
  cp -f "${BASE}/src/RecoilJets.h" "$snap_pp_header"
  cp -f "${BASE}/src_AuAu/RecoilJets_AuAu.h" "$snap_auau_header"

  cp -f "${user_root}/thesisAnalysis/install/lib/libcalo_reco.so" "$snap_lib_dir/"
  cp -f "${user_root}/thesisAnalysis/install/lib/libcalo_io.so" "$snap_lib_dir/"
  cp -f "${user_root}/thesisAnalysis/install/lib/libclusteriso.so" "$snap_lib_dir/"
  cp -f "${user_root}/thesisAnalysis/install/lib/libjetbase.so" "$snap_lib_dir/"
  [[ -f "${user_root}/thesisAnalysis/install/lib/libRecoilJets.so" ]] && cp -f "${user_root}/thesisAnalysis/install/lib/libRecoilJets.so" "$snap_lib_dir/"
  [[ -f "${user_root}/thesisAnalysis_auau/install/lib/libRecoilJetsAuAu.so" ]] && cp -f "${user_root}/thesisAnalysis_auau/install/lib/libRecoilJetsAuAu.so" "$snap_lib_dir/"

  # Copy companion ROOT PCM dictionaries so R__LOAD_LIBRARY doesn't spew missing-PCM errors
  cp -f "${user_root}/thesisAnalysis/install/lib/"*_rdict.pcm "$snap_lib_dir/" 2>/dev/null || true
  cp -f "${user_root}/thesisAnalysis_auau/install/lib/"*_rdict.pcm "$snap_lib_dir/" 2>/dev/null || true

  sed -i "s|#include \"/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Fun4All_recoilJets_unified_impl.C\"|#include \"${snap_impl}\"|" "$snap_macro"
  sed -i "s|#include \"/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Calo_Calib.C\"|#include \"${snap_calo}\"|" "$snap_impl"
  sed -i "s|#include \"/sphenix/u/patsfan753/scratch/thesisAnalysis/src/RecoilJets.h\"|#include \"${snap_pp_header}\"|" "$snap_impl"
  sed -i "s|#include \"/sphenix/u/patsfan753/scratch/thesisAnalysis/src_AuAu/RecoilJets_AuAu.h\"|#include \"${snap_auau_header}\"|" "$snap_impl"

  sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libcalo_reco.so)|" "$snap_impl"
  sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_io.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libcalo_io.so)|" "$snap_impl"
  sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libclusteriso.so)|" "$snap_impl"
  sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libjetbase.so)|" "$snap_impl"
  sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libRecoilJets.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libRecoilJets.so)|" "$snap_impl"
  sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis_auau/install/lib/libRecoilJetsAuAu.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libRecoilJetsAuAu.so)|" "$snap_impl"

  chmod +x "$snap_wrapper"

  BULK_FROZEN_EXE="$snap_wrapper"
  BULK_FROZEN_MACRO="$snap_macro"

  say "Pipeline snapshot created:"
  say "  snapshot dir : ${snap_dir}"
  say "  frozen exe   : ${BULK_FROZEN_EXE}"
  say "  frozen macro : ${BULK_FROZEN_MACRO}"
}

cleanup_bulk_snapshots_for_tag() {
  mkdir -p "$SNAPSHOT_ROOT"
  local n_old
  n_old=$(find "$SNAPSHOT_ROOT" -mindepth 1 -maxdepth 1 -type d -name "${TAG}_*" 2>/dev/null | wc -l | awk '{print $1}')
  if (( n_old > 0 )); then
    say "Removing ${n_old} older bulk snapshot dir(s) for tag ${TAG}"
    find "$SNAPSHOT_ROOT" -mindepth 1 -maxdepth 1 -type d -name "${TAG}_*" -exec rm -rf {} +
  fi
}

# ------------------------ Helpers --------------------------
usage() {
  cat <<USAGE
${BOLD}Usage:${RST}

${BOLD}DATA modes:${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> local [Nevents] [VERBOSE=N]${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> CHECKJOBS [groupSize N]${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> splitGoldenRunList [groupSize N] [maxJobs M]${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor testJob${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor round <N> [groupSize N] [firstChunk]${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor all [groupSize N]${RST}

${BOLD}SIM mode (photonJet10/20 merged):${RST}
  ${BOLD}$0 isSim local [Nevents] [VERBOSE=N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim CHECKJOBS [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorTest [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorDoAll [groupSize N] [SAMPLE=run28_photonjet10]${RST}

${BOLD}SIM mode (photonJet5 single-slice):${RST}
  ${BOLD}$0 isSimJet5 local [Nevents] [VERBOSE=N] [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimJet5 CHECKJOBS [groupSize N] [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimJet5 condorTest [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimJet5 condorDoAll [groupSize N] [SAMPLE=run28_jet5]${RST}

${BOLD}SIM mode (MinBias DETROIT):${RST}
  ${BOLD}$0 isSimMB local [Nevents] [VERBOSE=N] [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB CHECKJOBS [groupSize N] [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB condorTest [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB condorDoAll [groupSize N] [SAMPLE=run28_detroit]${RST}

Examples:
  $0 isPP  local 5000
  $0 isAuAu splitGoldenRunList groupSize 3 maxJobs 12000
  $0 isAuAu condor round 1
  $0 isPP  condor all groupSize 4
  $0 isPP  CHECKJOBS groupSize 4

  $0 isSim CHECKJOBS groupSize 5
  $0 isSim local 5000
  $0 isSim condorTest
  $0 isSim condorDoAll groupSize 5

  $0 isSimJet5 local 5000
  $0 isSimJet5 condorDoAll groupSize 5
  $0 isSimMB local 5000
  $0 isSimMB condorDoAll groupSize 5
USAGE
  exit 2
}

need_cmd() { command -v "$1" >/dev/null 2>&1 || { err "Missing required command: $1"; exit 3; }; }

# ------------------------ SIM YAML sweep helpers --------------------------
trim_ws() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf "%s" "$s"
}

sim_yaml_master_path() {
  if [[ -n "${RJ_CONFIG_YAML:-}" ]]; then
    local p; p="$(trim_ws "${RJ_CONFIG_YAML}")"
    [[ -n "$p" ]] && { echo "$p"; return; }
  fi
  echo "$SIM_YAML_DEFAULT"
}

yaml_get_values() {
  local key="$1" file="$2"
  local line rhs inside
  line=$(grep -E "^[[:space:]]*${key}:" "$file" | head -n 1 || true)
  [[ -n "$line" ]] || { err "Missing YAML key '${key}' in ${file}"; exit 70; }
  line="${line%%#*}"
  rhs="${line#*:}"
  rhs="$(trim_ws "$rhs")"

  if [[ "$rhs" == \[*\] ]]; then
    inside="${rhs#\[}"
    inside="${inside%\]}"
    IFS=',' read -ra parts <<< "$inside"
    for p in "${parts[@]}"; do
      p="$(trim_ws "$p")"
      [[ -n "$p" ]] && echo "$p"
    done
  else
    [[ -n "$rhs" ]] && echo "$rhs"
  fi
}

sim_is_close() { awk -v x="$1" -v y="$2" 'BEGIN{d=x-y; if(d<0)d=-d; exit !(d<1e-6)}'; }

sim_pt_tag() {
  local pt="$1"
  if [[ "$pt" =~ ^([0-9]+)\.0+$ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    local s="$pt"
    s="${s//./p}"
    s="${s//-/m}"
    echo "$s"
  fi
}

sim_b2b_tag() {
  local frac="$1"
  if sim_is_close "$frac" "0.5"; then
    echo "pi_2"
  elif sim_is_close "$frac" "0.875"; then
    echo "7pi_8"
  else
    local s="$frac"
    s="${s//./p}"
    s="${s//-/m}"
    echo "piFrac${s}"
  fi
}

sim_vz_tag() {
  local vz="$1"
  if [[ "$vz" =~ ^([0-9]+)\.0+$ ]]; then
    echo "vz${BASH_REMATCH[1]}"
  else
    local s="$vz"
    s="${s//./p}"
    s="${s//-/m}"
    echo "vz${s}"
  fi
}

sim_cone_tag() {
  local r="$1"
  local r100
  r100=$(awk -v r="$r" 'BEGIN{v=int(r*100+0.5); printf "%d", v}')
  echo "isoR${r100}"
}

sim_iso_tag() {
  local sliding="$1" fixed="$2"
  if [[ "$sliding" == "true" ]]; then
    echo "isSliding"
  else
    if [[ "$fixed" =~ ^([0-9]+)\.0+$ ]]; then
      echo "fixedIso${BASH_REMATCH[1]}GeV"
    else
      local s="$fixed"
      s="${s//./p}"
      echo "fixedIso${s}GeV"
    fi
  fi
}

# Build parallel arrays: iso_tags[], iso_sliding[], iso_fixed[]
# from YAML keys: isSlidingAndFixed, isSlidingIso, fixedGeV[]
build_iso_modes() {
  local yaml_file="$1"
  iso_tags=()
  iso_sliding=()
  iso_fixed=()

  local _both; _both="$(yaml_get_values "isSlidingAndFixed" "$yaml_file" 2>/dev/null || echo "false")"
  local _slide; _slide="$(yaml_get_values "isSlidingIso" "$yaml_file" 2>/dev/null || echo "false")"
  local -a _fixeds
  mapfile -t _fixeds < <( yaml_get_values "fixedGeV" "$yaml_file" 2>/dev/null )
  (( ${#_fixeds[@]} )) || _fixeds=( "2.0" )

  _both="$(trim_ws "$_both")"
  _slide="$(trim_ws "$_slide")"

  if [[ "$_both" == "true" ]]; then
    # sliding mode entry
    iso_tags+=( "$(sim_iso_tag "true" "0")" )
    iso_sliding+=( "true" )
    iso_fixed+=( "${_fixeds[0]}" )
    # one fixed entry per fixedGeV value
    for fv in "${_fixeds[@]}"; do
      iso_tags+=( "$(sim_iso_tag "false" "$fv")" )
      iso_sliding+=( "false" )
      iso_fixed+=( "$fv" )
    done
  elif [[ "$_slide" == "true" ]]; then
    iso_tags+=( "$(sim_iso_tag "true" "0")" )
    iso_sliding+=( "true" )
    iso_fixed+=( "${_fixeds[0]}" )
  else
    for fv in "${_fixeds[@]}"; do
      iso_tags+=( "$(sim_iso_tag "false" "$fv")" )
      iso_sliding+=( "false" )
      iso_fixed+=( "$fv" )
    done
  fi
}

# Determine UE pipeline matrix modes based on dataset type.
# For AuAu/OO: reads clusterUEpipeline array from YAML; adds to naming tag.
# For pp/SIM:  always "noSub", NOT included in naming tag.
read_uepipe_modes() {
  local yaml="$1" tag="$2"
  uepipe_modes=()
  uepipe_in_tag=0
  if [[ "$tag" == "auau" || "$tag" == "oo" || "$tag" == "simembedded" || "$tag" == "simembeddedinclusive" ]]; then
    mapfile -t uepipe_modes < <( yaml_get_values "clusterUEpipeline" "$yaml" 2>/dev/null )
    # Translate legacy bool values
    local -a cleaned=()
    for v in "${uepipe_modes[@]}"; do
      v="$(trim_ws "$v")"
      case "$v" in
        true|1)  cleaned+=( "variantA" ) ;;
        false|0) cleaned+=( "noSub" ) ;;
        noSub|baseVariant|variantA|variantB) cleaned+=( "$v" ) ;;
        *) err "Unknown clusterUEpipeline value '$v'"; exit 73 ;;
      esac
    done
    uepipe_modes=( "${cleaned[@]}" )
    (( ${#uepipe_modes[@]} )) || uepipe_modes=( "noSub" )
    uepipe_in_tag=1
  else
    uepipe_modes=( "noSub" )
    uepipe_in_tag=0
  fi
}

sim_make_yaml_override() {
  local master="$1" pt="$2" frac="$3" vz="$4" cone="$5" sliding="$6" fixed="$7" uepipe="$8" tag="$9" stamp="${10:-}"
  mkdir -p "$SIM_YAML_OVERRIDE_DIR"
  local out="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${tag}${stamp:+_${stamp}}.yaml"

  [[ -s "$master" ]] || { err "Master YAML not found or empty: $master"; exit 71; }
  grep -Eq '^[[:space:]]*jet_pt_min:' "$master" || { err "YAML missing key 'jet_pt_min' in $master"; exit 71; }
  grep -Eq '^[[:space:]]*back_to_back_dphi_min_pi_fraction:' "$master" || { err "YAML missing key 'back_to_back_dphi_min_pi_fraction' in $master"; exit 71; }
  grep -Eq '^[[:space:]]*vz_cut_cm:' "$master" || { err "YAML missing key 'vz_cut_cm' in $master"; exit 71; }
  grep -Eq '^[[:space:]]*coneR:' "$master" || { err "YAML missing key 'coneR' in $master"; exit 71; }

  sed -E \
    -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${pt}|" \
    -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${frac}|" \
    -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${vz}|" \
    -e "s|^([[:space:]]*coneR:).*|\\1 ${cone}|" \
    -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${sliding}|" \
    -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${fixed}|" \
    -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${uepipe}|" \
    "$master" > "$out"

  echo "$out"
}

# Dataset resolver (sets globals: DATASET, GOLDEN, LIST_DIR, DEST_BASE, TAG, STAGE_DIR, ROUND_DIR)
resolve_dataset() {
  case "${1:-}" in
    isPP|pp|PP)
      DATASET="isPP"
      GOLDEN="$PP_GOLDEN"
      LIST_DIR="$PP_LIST_DIR"
      LIST_PREFIX="dst_jetcalo"
      DEST_BASE="$PP_DEST_BASE"
      TAG="pp"
      MACRO="${BASE}/macros/Fun4All_recoilJets.C"
      EXE="${BASE}/RecoilJets_Condor.sh"
      ;;
    isPPrun25|pprun25|pp25|PP25)
      DATASET="isPPrun25"
      GOLDEN="$PP25_GOLDEN"
      LIST_DIR="$PP25_LIST_DIR"
      LIST_PREFIX="dst_calofitting"
      DEST_BASE="$PP25_DEST_BASE"
      TAG="pp25"
      MACRO="${BASE}/macros/Fun4All_recoilJets.C"
      EXE="${BASE}/RecoilJets_Condor.sh"
      ;;
    isAuAu|auau|AA)
      DATASET="isAuAu"
      GOLDEN="$AA_GOLDEN"
      LIST_DIR="$AA_LIST_DIR"
      LIST_PREFIX="dst_calofitting"
      DEST_BASE="$AA_DEST_BASE"
      TAG="auau"
      MACRO="${BASE}/macros/Fun4All_recoilJets_AuAu.C"
      EXE="${BASE}/RecoilJets_Condor_AuAu.sh"
      ;;
    isOO|oo|OO)
      DATASET="isOO"
      GOLDEN="$OO_GOLDEN"
      LIST_DIR="$OO_LIST_DIR"
      LIST_PREFIX="dst_calofitting"
      DEST_BASE="$OO_DEST_BASE"
      TAG="oo"
      MACRO="${BASE}/macros/Fun4All_recoilJets_AuAu.C"
      EXE="${BASE}/RecoilJets_Condor_AuAu.sh"
      ;;
    isSim|sim|SIM)
      DATASET="isSim"
      GOLDEN=""        # not used in sim mode
      LIST_DIR=""      # not used in sim mode
      LIST_PREFIX=""   # not used in sim mode
      DEST_BASE="$SIM_DEST_BASE"   # output dir will be DEST_BASE/SIM_SAMPLE
      TAG="sim"
      MACRO="${BASE}/macros/Fun4All_recoilJets.C"
      EXE="${BASE}/RecoilJets_Condor.sh"
      IS_SIM=1
      ;;
    isSimEmbedded|issimembedded|simembedded|SIMEMBEDDED)
      DATASET="isSimEmbedded"
      GOLDEN=""
      LIST_DIR=""
      LIST_PREFIX=""
      DEST_BASE="$SIMEMBED_DEST_BASE"
      TAG="simembedded"
      MACRO="${BASE}/macros/Fun4All_recoilJets_AuAu.C"
      EXE="${BASE}/RecoilJets_Condor_AuAu.sh"
      IS_SIM=1
      SIM_SAMPLE_DEFAULT="run28_embeddedPhoton20"
      SIM_SAMPLE="$SIM_SAMPLE_DEFAULT"
      ;;
    isSimEmbeddedInclusive|issimembeddedinclusive|simembeddedinclusive|SIMEMBEDDEDINCLUSIVE)
      DATASET="isSimEmbeddedInclusive"
      GOLDEN=""
      LIST_DIR=""
      LIST_PREFIX=""
      DEST_BASE="$SIMEMBED_DEST_BASE"
      TAG="simembeddedinclusive"
      MACRO="${BASE}/macros/Fun4All_recoilJets_AuAu.C"
      EXE="${BASE}/RecoilJets_Condor_AuAu.sh"
      IS_SIM=1
      SIM_SAMPLE_DEFAULT="run28_embeddedJet20"
      SIM_SAMPLE="$SIM_SAMPLE_DEFAULT"
      ;;
    isSimJet5|isSimjet5|simjet5|SIMJET5)
      DATASET="isSimJet5"
      GOLDEN=""
      LIST_DIR=""
      LIST_PREFIX=""
      DEST_BASE="$SIMJET5_DEST_BASE"
      TAG="simjet5"
      MACRO="${BASE}/macros/Fun4All_recoilJets.C"
      EXE="${BASE}/RecoilJets_Condor.sh"
      IS_SIM=1
      SIM_SAMPLE_DEFAULT="run28_jet5"
      SIM_SAMPLE="$SIM_SAMPLE_DEFAULT"
      ;;
    isSimMB|isSimMb|issimmb|simmb|SIMMB)
      DATASET="isSimMB"
      GOLDEN=""
      LIST_DIR=""
      LIST_PREFIX=""
      DEST_BASE="$SIMMB_DEST_BASE"
      TAG="simmb"
      MACRO="${BASE}/macros/Fun4All_recoilJets.C"
      EXE="${BASE}/RecoilJets_Condor.sh"
      IS_SIM=1
      SIM_SAMPLE_DEFAULT="run28_detroit"
      SIM_SAMPLE="$SIM_SAMPLE_DEFAULT"
      ;;
    *)
      usage
      ;;
  esac

  STAGE_DIR="${STAGE_ROOT}/${TAG}"
  ROUND_DIR="${ROUND_ROOT}/${TAG}"
  mkdir -p "$STAGE_DIR" "$ROUND_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR" "$SUB_DIR"
}

# Count ceil(n/d)
ceil_div() { local n="$1" d="$2"; echo $(( (n + d - 1) / d )); }

# Format run to 8 digits
run8() { printf "%08d" "$((10#$1))"; }

# Check if a GL1 trigger bit is active for a run (scaledown != -1).
# Usage: is_trigger_active <run8> <bit> ; returns 0 if active, 1 otherwise.
is_trigger_active() {
  local run8="$1" bit="$2"
  [[ -z "$bit" ]] && return 0
  local run_dec=$((10#$run8))
  local val
  val=$(psql -h sphnxdaqdbreplica -d daq -At -q -F $'\t' \
        -c "SELECT scaledown${bit} FROM gl1_scaledown WHERE runnumber=${run_dec};" \
        | tr -d '[:space:]')
  [[ -n "$val" && "$val" != "-1" ]]
}

# IsolationAudit local-only run selector:
#   - scans the existing Au+Au golden run manifest in sorted order
#   - requires BOTH:
#       1) per-run DST_CALOFITTING list exists and is non-empty
#       2) photon_10_plus_MBD_NS_geq_2_vtx_lt_150 trigger path is active
pick_first_iso_ping_run() {
  local trig_bit="$1"

  [[ -s "$GOLDEN" ]] || return 1

  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    local r8
    r8="$(run8 "$rn")"

    local src="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
    [[ -s "$src" ]] || continue

    if is_trigger_active "$r8" "$trig_bit"; then
      echo "$r8"
      return 0
    fi
  done < <(grep -E '^[0-9]+' "$GOLDEN" | sort -n)

  return 1
}

# Build per-run grouped list files; prints absolute paths to grouped lists, one per line
#   make_groups <run8> <groupSize>  → writes STAGE_DIR/run<run8>_grpXXX.list files
make_groups() {
  local r8="$1" gs="$2"
  local src="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
  [[ -s "$src" ]] || { warn "List is missing/empty for run ${r8} → $src"; return 1; }

  # Clean old groups for this run
  rm -f "${STAGE_DIR}/run${r8}_grp"*.list 2>/dev/null || true
  rm -f "${STAGE_DIR}/run${r8}_LOCAL_"*.list 2>/dev/null || true

  local nfiles; nfiles=$(wc -l < "$src" | awk '{print $1}')
  local ngroups; ngroups=$(ceil_div "$nfiles" "$gs")

  if [[ "${RJ_GROUP_TRACE:-0}" == "1" ]]; then
    say "[group] run=${r8} src=$(basename "$src") files=${nfiles} groupSize=${gs} groups=${ngroups}" >&2
  fi

  local start=1
  local g=1
  while (( g <= ngroups )); do
    local end=$(( start + gs - 1 ))
    local out="${STAGE_DIR}/run${r8}_grp$(printf "%03d" "$g").list"
    # sed is 1-indexed on lines; clamp 'end' to nfiles
    if (( end > nfiles )); then end="$nfiles"; fi
    sed -n "${start},${end}p" "$src" > "$out"
    echo "$out"
    start=$(( end + 1 ))
    g=$(( g + 1 ))
  done

  if [[ "${RJ_GROUP_TRACE:-0}" == "1" ]]; then
    say "[group] run=${r8} finished grouping under ${STAGE_DIR}" >&2
  fi
}

# ------------------------ Simulation helpers --------------------------
# Initializes paths for isSim mode and prepares a cleaned master list.
sim_init() {
  SIM_DIR="${SIM_ROOT}/${SIM_SAMPLE}"
  [[ -d "$SIM_DIR" ]] || { err "Sim sample directory not found: $SIM_DIR"; exit 20; }

  SIM_STAGE_DIR="${STAGE_DIR}/${SIM_SAMPLE}"
  mkdir -p "$SIM_STAGE_DIR"

  # Build a staged 5-column master list from the matched lists created by makeThesisSimLists.sh
  local calo="${SIM_DIR}/DST_CALO_CLUSTER.matched.list"
  local g4="${SIM_DIR}/G4Hits.matched.list"
  local jets="${SIM_DIR}/DST_JETS.matched.list"
  local glob="${SIM_DIR}/DST_GLOBAL.matched.list"
  local mbd="${SIM_DIR}/DST_MBD_EPD.matched.list"

  [[ -s "$calo" ]] || { err "Missing: $calo"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$g4"   ]] || { err "Missing: $g4";   err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$jets" ]] || { err "Missing: $jets"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$glob" ]] || { err "Missing: $glob"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$mbd"  ]] || { err "Missing: $mbd";  err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }

  SIM_MASTER_LIST="${SIM_STAGE_DIR}/sim_${SIM_SAMPLE}_${SIM_MASTER5_NAME}"
  local _ninput; _ninput=$(wc -l < "$calo" | tr -d ' ')
  if [[ "${ACTION:-}" != "CHECKJOBS" ]]; then
    say "    [sim_init] paste 5 matched lists (${_ninput} lines each) → master list…" >&2
  fi
  paste "$calo" "$g4" "$jets" "$glob" "$mbd" > "$SIM_MASTER_LIST"
  if [[ "${ACTION:-}" != "CHECKJOBS" ]]; then
    say "    [sim_init] paste done ($(wc -l < "$SIM_MASTER_LIST" | tr -d ' ') lines)" >&2
  fi

  # Clean master list: strip blank lines + comment lines (keeps ALL columns)
  SIM_CLEAN_LIST="${SIM_STAGE_DIR}/sim_${SIM_SAMPLE}_PAIR_MASTER_CLEAN.list"
  if [[ "${ACTION:-}" != "CHECKJOBS" ]]; then
    say "    [sim_init] cleaning master list (strip blanks/comments)…" >&2
  fi
  grep -E -v '^[[:space:]]*($|#)' "$SIM_MASTER_LIST" > "$SIM_CLEAN_LIST" || true
  if [[ "${ACTION:-}" != "CHECKJOBS" ]]; then
    say "    [sim_init] clean list ready: $(wc -l < "$SIM_CLEAN_LIST" | tr -d ' ') lines → $(basename "$SIM_CLEAN_LIST")" >&2
  fi
  [[ -s "$SIM_CLEAN_LIST" ]] || { err "Sim master list empty after cleaning: $SIM_MASTER_LIST"; exit 22; }

  SIM_OUT_DIR="${DEST_BASE}/${SIM_SAMPLE}"
  [[ "${ACTION:-}" != "CHECKJOBS" ]] && mkdir -p "$SIM_OUT_DIR"

  if [[ -n "${SIM_CFG_TAG:-}" ]]; then
    SIM_JOB_PREFIX="sim_${SIM_SAMPLE}_${SIM_CFG_TAG}"
  else
    SIM_JOB_PREFIX="sim_${SIM_SAMPLE}"
  fi
}

# Build grouped chunk lists for isSim (one job per chunk list)
make_sim_groups() {
  local gs="$1"
  sim_init

  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp"*.list 2>/dev/null || true

  local _nclean; _nclean=$(wc -l < "$SIM_CLEAN_LIST" | tr -d ' ')
  local _nexpect=$(( (_nclean + gs - 1) / gs ))
  say "    [make_sim_groups] splitting ${_nclean} lines into chunks of ${gs} (expect ~${_nexpect} groups)…" >&2

  # Use split(1) for O(n) grouping instead of sed-in-a-loop (critical for 200k+ file samples)
  local prefix="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp_raw_"
  split -l "$gs" -d -a 5 "$SIM_CLEAN_LIST" "$prefix"
  say "    [make_sim_groups] split done, renaming chunk files…" >&2

  # Rename split's numeric suffixes to our grpNNN.list naming convention
  local g=0
  local first_group=""
  local last_group=""
  for raw in "${prefix}"*; do
    [[ -s "$raw" ]] || { rm -f "$raw"; continue; }
    (( g+=1 ))
    local out="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp$(printf "%03d" "$g").list"
    mv "$raw" "$out"
    [[ -z "$first_group" ]] && first_group="$out"
    last_group="$out"
    echo "$out"
  done
  say "    [make_sim_groups] renamed ${g} group files" >&2
  if [[ "${RJ_GROUP_TRACE:-0}" == "1" && "$g" -gt 0 ]]; then
    say "    [make_sim_groups] first=$(basename "$first_group")  last=$(basename "$last_group")" >&2
  fi
}

# Dry-run job count for isSim
check_jobs_sim() {
  local gs="$GROUP_SIZE"
  local master_yaml
  local n_cfg=1
  local per_cfg_jobs=0
  local total_jobs=0
  local -a samples=()

  master_yaml="$(sim_yaml_master_path)"
  [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

  mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
  mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
  mapfile -t sim_vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
  mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
  (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
  (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
  (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
  (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
  build_iso_modes "$master_yaml"
  read_uepipe_modes "$master_yaml" "$TAG"

  if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
    case "$DATASET" in
      isSimEmbedded)          samples=( "run28_embeddedPhoton10" "run28_embeddedPhoton20" ) ;;
      isSimEmbeddedInclusive) samples=( "run28_embeddedJet10" "run28_embeddedJet20" ) ;;
      isSimJet5)              samples=( "run28_jet5" ) ;;
      isSimMB)                samples=( "run28_detroit" ) ;;
      *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
    esac
  else
    samples=( "${SIM_SAMPLE}" )
  fi

  n_cfg=$(( ${#sim_pts[@]} * ${#sim_fracs[@]} * ${#sim_vzs[@]} * ${#sim_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))

  say "CHECKJOBS (dataset=${DATASET}, tag=${TAG})"
  say "  YAML master         : ${master_yaml}"
  say "  groupSize           : ${gs}"
  say "  samples             : [${samples[*]}]"
  say "  master list source  : ${SIM_ROOT}"
  echo

  say "${BOLD}Matrix dimensions:${RST}"
  say "  jet_pt_min                        : [${sim_pts[*]}]  (${#sim_pts[@]} values)"
  say "  back_to_back_dphi_min_pi_fraction : [${sim_fracs[*]}]  (${#sim_fracs[@]} values)"
  say "  vz_cut_cm                         : [${sim_vzs[*]}]  (${#sim_vzs[@]} values)"
  say "  coneR                             : [${sim_cones[*]}]  (${#sim_cones[@]} values)"
  say "  iso modes                         : [${iso_tags[*]}]  (${#iso_tags[@]} values)"
  say "  clusterUEpipeline                 : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} values)"
  say "  cfg combos (product)              : ${BOLD}${n_cfg}${RST}"
  echo

  say "${BOLD}Per-sample input file counts:${RST}"
  for samp in "${samples[@]}"; do
    SIM_SAMPLE="$samp"
    sim_init
    local nfiles
    local njobs
    nfiles=$(wc -l < "$SIM_CLEAN_LIST" | awk '{print $1}')
    njobs=$(ceil_div "$nfiles" "$gs")
    per_cfg_jobs=$(( per_cfg_jobs + njobs ))
    say "  sample=${BOLD}${SIM_SAMPLE}${RST}  input_files=${nfiles}  jobs_per_cfg=${njobs}  master_list=${SIM_MASTER_LIST}"
  done
  echo

  say "${BOLD}Full cfg tag list (${n_cfg} entries × ${#samples[@]} samples = ${BOLD}$((n_cfg * ${#samples[@]}))${RST} submit blocks):${RST}"
  printf "  ${BOLD}%3s │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s${RST}\n" \
         "#" "cfg_tag" "pt" "frac" "vz" "cone" "iso" "uepipe"
  printf "  ────┼────────────────────────────────────────────────────────────────────────┼───── ────── ─────── ─────── ────────────────── ────────────\n"
  local cfg_num=0
  for _pt in "${sim_pts[@]}"; do
    for _frac in "${sim_fracs[@]}"; do
      for _vz in "${sim_vzs[@]}"; do
        for _cone in "${sim_cones[@]}"; do
          for (( _ci=0; _ci<${#iso_tags[@]}; _ci++ )); do
            for _ue in "${uepipe_modes[@]}"; do
              (( cfg_num+=1 ))
              local _tag
              _tag="jetMinPt$(sim_pt_tag "$_pt")_$(sim_b2b_tag "$_frac")_$(sim_vz_tag "$_vz")_$(sim_cone_tag "$_cone")_${iso_tags[$_ci]}"
              (( uepipe_in_tag )) && _tag="${_tag}_${_ue}"
              printf "  %3d │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s\n" \
                     "$cfg_num" "$_tag" "$_pt" "$_frac" "$_vz" "$_cone" "${iso_tags[$_ci]}" "$_ue"
            done
          done
        done
      done
    done
  done
  echo

  total_jobs=$(( n_cfg * per_cfg_jobs ))
  say "${BOLD}Job count summary:${RST}"
  say "  cfg combos             : ${n_cfg}"
  say "  jobs per combo (Σsamp) : ${per_cfg_jobs}"
  say "  ─────────────────────────────────"
  say "  ${BOLD}TOTAL CONDOR JOBS        : ${BOLD}${total_jobs}${RST}"
  echo

  say "${BOLD}groupSize sensitivity (total Condor jobs submitted):${RST}"
  echo
  printf "  ${BOLD}%-12s │ %-14s${RST}\n" "groupSize" "TOTAL JOBS"
  printf "  ─────────────┼───────────────\n"
  for _gs_try in 4 6 7 8 10; do
    local _per_cfg_jobs_try=0
    for samp in "${samples[@]}"; do
      SIM_SAMPLE="$samp"
      sim_init
      local _nft
      local _nj
      _nft=$(wc -l < "$SIM_CLEAN_LIST" | awk '{print $1}')
      _nj=$(ceil_div "$_nft" "$_gs_try")
      _per_cfg_jobs_try=$(( _per_cfg_jobs_try + _nj ))
    done
    local _tj=$(( n_cfg * _per_cfg_jobs_try ))
    local _marker=""
    [[ "$_gs_try" -eq "$gs" ]] && _marker=" ${BOLD}← current${RST}"
    printf "  %-12s │ %14s%s\n" "$_gs_try" "$_tj" "$_marker"
  done
  echo

  say "Output tree: each tag becomes a subdirectory under ${DEST_BASE}/"
  local _ex_tag
  _ex_tag="jetMinPt$(sim_pt_tag "${sim_pts[0]}")_$(sim_b2b_tag "${sim_fracs[0]}")_$(sim_vz_tag "${sim_vzs[0]}")_$(sim_cone_tag "${sim_cones[0]}")_${iso_tags[0]}"
  (( uepipe_in_tag )) && _ex_tag="${_ex_tag}_${uepipe_modes[0]}"
  say "  e.g. ${DIM}${DEST_BASE}/${_ex_tag}/run28_photonjet10/*.root${RST}"
}

# Wipe previous sim outputs + sim-tagged logs/out/err (ONLY used before condorDoAll)
wipe_sim_artifacts() {
  sim_init

  say "WIPING previous isSim artifacts for sample ${BOLD}${SIM_SAMPLE}${RST}"
  say "  Output ROOT dir : ${SIM_OUT_DIR}"
  say "  Logs prefix     : ${SIM_JOB_PREFIX}"

  rm -f "${SIM_OUT_DIR}/"*.root 2>/dev/null || true
  find "${SIM_OUT_DIR}" -maxdepth 1 -name "*_LOCAL_*.root" -delete 2>/dev/null || true
  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp"*.list 2>/dev/null || true
  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_LOCAL_"*.list 2>/dev/null || true
  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_condorTest_"*.list 2>/dev/null || true

  rm -f "${LOG_DIR}/${SIM_JOB_PREFIX}.job."*.log 2>/dev/null || true
  rm -f "${OUT_DIR}/${SIM_JOB_PREFIX}.job."*.out 2>/dev/null || true
  rm -f "${ERR_DIR}/${SIM_JOB_PREFIX}.job."*.err 2>/dev/null || true
}

# Split golden run list into “round” files so that each round’s sum of jobs
# (sum over runs of ceil(nFiles/GROUP_SIZE)) ≤ MAX_JOBS
split_golden() {
  local gs="$1" cap="$2"
  say "Splitting golden run list → groups of ${BOLD}${gs}${RST} files/job, cap ${BOLD}${cap}${RST} jobs/round"
  [[ -s "$GOLDEN" ]] || { err "Golden run list not found or empty: $GOLDEN"; exit 4; }

  # Clean previous rounds for this dataset
  rm -f "${ROUND_DIR}/goldenRuns_${TAG}_segment"*.txt 2>/dev/null || true

  local seg=1 jobs_in_seg=0 runs_in_seg=0
  local segfile="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
  : > "$segfile"

  local total_jobs=0 total_runs=0
  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    local r8; r8="$(run8 "$rn")"
    local lf="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
    if [[ ! -s "$lf" ]]; then
      warn "No per-run list for ${r8}; skipping"
      continue
    fi
    local nfiles; nfiles=$(wc -l < "$lf" | awk '{print $1}')
    local nj; nj=$(ceil_div "$nfiles" "$gs")

    # If adding this run would exceed the cap, flush segment
    if (( jobs_in_seg + nj > cap )) && (( runs_in_seg > 0 )); then
      say "Round ${BOLD}${seg}${RST}: runs=${runs_in_seg}, jobs=${jobs_in_seg} → ${segfile}"
      seg=$(( seg + 1 ))
      segfile="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
      : > "$segfile"
      jobs_in_seg=0; runs_in_seg=0
    fi

    echo "$r8" >> "$segfile"
    jobs_in_seg=$(( jobs_in_seg + nj ))
    runs_in_seg=$(( runs_in_seg + 1 ))
    total_jobs=$(( total_jobs + nj ))
    total_runs=$(( total_runs + 1 ))
  done < "$GOLDEN"

  say "Round ${BOLD}${seg}${RST}: runs=${runs_in_seg}, jobs=${jobs_in_seg} → ${segfile}"
  say "Total runs=${BOLD}${total_runs}${RST}, total condensed jobs≈${BOLD}${total_jobs}${RST}"
}

# Dry-run job count for a hypothetical "condor all" at current GROUP_SIZE
#   check_jobs_all
# Reads the golden list and per-run lists; prints totals only. No side effects.
check_jobs_all() {
  local gs="$GROUP_SIZE"
  [[ -s "$GOLDEN" ]] || { err "Golden run list not found or empty: $GOLDEN"; exit 4; }

  # Read YAML matrix dimensions (same logic as condor all / local paths)
  local data_yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  local -a ck_pts ck_fracs ck_vzs ck_cones
  mapfile -t ck_pts   < <( yaml_get_values "jet_pt_min" "$data_yaml_src" )
  mapfile -t ck_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$data_yaml_src" )
  mapfile -t ck_vzs   < <( yaml_get_values "vz_cut_cm" "$data_yaml_src" )
  mapfile -t ck_cones < <( yaml_get_values "coneR" "$data_yaml_src" )
  (( ${#ck_pts[@]} ))   || { err "No values found for jet_pt_min in $data_yaml_src"; exit 72; }
  (( ${#ck_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $data_yaml_src"; exit 72; }
  (( ${#ck_vzs[@]} ))   || { err "No values found for vz_cut_cm in $data_yaml_src"; exit 72; }
  (( ${#ck_cones[@]} )) || { err "No values found for coneR in $data_yaml_src"; exit 72; }
  build_iso_modes "$data_yaml_src"
  read_uepipe_modes "$data_yaml_src" "$TAG"

  local total_runs=0 total_files=0 base_jobs=0 missing=0

  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    local r8; r8="$(run8 "$rn")"
    local lf="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"

    if [[ ! -s "$lf" ]]; then
      warn "No per-run list for ${r8}; skipping"
      ((missing++))
      continue
    fi

    local nfiles; nfiles=$(wc -l < "$lf" | awk '{print $1}')
    local nj;     nj=$(ceil_div "$nfiles" "$gs")

    (( total_runs  += 1    ))
    (( total_files += nfiles ))
    (( base_jobs   += nj   ))
  done < "$GOLDEN"

  local n_matrix=$(( ${#ck_pts[@]} * ${#ck_fracs[@]} * ${#ck_vzs[@]} * ${#ck_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
  local total_jobs=$(( n_matrix * base_jobs ))

  say "CHECKJOBS (dataset=${DATASET}, tag=${TAG})"
  say "  YAML source          : ${data_yaml_src}"
  say "  groupSize            : ${gs}"
  echo
  say "${BOLD}Matrix dimensions:${RST}"
  say "  jet_pt_min           : [${ck_pts[*]}]  (${#ck_pts[@]} values)"
  say "  b2b_dphi_pi_fraction : [${ck_fracs[*]}]  (${#ck_fracs[@]} values)"
  say "  vz_cut_cm            : [${ck_vzs[*]}]  (${#ck_vzs[@]} values)"
  say "  coneR                : [${ck_cones[*]}]  (${#ck_cones[@]} values)"
  say "  iso modes            : [${iso_tags[*]}]  (${#iso_tags[@]} values)"
  say "  clusterUEpipeline    : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} values)"
  say "  matrix combos        : ${BOLD}${n_matrix}${RST}"
  echo
  say "${BOLD}Per-matrix-cell base:${RST}"
  say "  golden runs (w/lists): ${total_runs}"
  say "  total input files    : ${total_files}"
  say "  base jobs per combo  : ${base_jobs}"
  (( missing > 0 )) && warn "  runs skipped (missing lists): ${missing}"
  echo

  say "${BOLD}Full cfg tag list (${n_matrix} entries):${RST}"
  echo
  printf "  ${BOLD}%3s │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s${RST}\n" \
         "#" "cfg_tag" "pt" "frac" "vz" "cone" "iso" "uepipe"
  printf "  ────┼────────────────────────────────────────────────────────────────────────┼───── ────── ─────── ─────── ────────────────── ────────────\n"
  local cfg_num=0
  for _pt in "${ck_pts[@]}"; do
  for _frac in "${ck_fracs[@]}"; do
  for _vz in "${ck_vzs[@]}"; do
    for _cone in "${ck_cones[@]}"; do
    for (( _ci=0; _ci<${#iso_tags[@]}; _ci++ )); do
    for _ue in "${uepipe_modes[@]}"; do
      (( ++cfg_num ))
      local _dpt; _dpt="jetMinPt$(sim_pt_tag "$_pt")"
      local _dfrac; _dfrac="$(sim_b2b_tag "$_frac")"
      local _dvz; _dvz="$(sim_vz_tag "$_vz")"
      local _dcone; _dcone="$(sim_cone_tag "$_cone")"
      local _tag="${_dpt}_${_dfrac}_${_dvz}_${_dcone}_${iso_tags[$_ci]}"
      (( uepipe_in_tag )) && _tag="${_tag}_${_ue}"
      printf "  %3d │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s\n" \
             "$cfg_num" "$_tag" "$_pt" "$_frac" "$_vz" "$_cone" "${iso_tags[$_ci]}" "$_ue"
    done
    done
    done
  done
  done
  done
  echo

  say "${BOLD}Job count summary (groupSize=${gs}):${RST}"
  say "  matrix combos          : ${n_matrix}"
  say "  base jobs per combo    : ${base_jobs}"
  say "  ─────────────────────────────────"
  say "  ${BOLD}TOTAL CONDOR JOBS        : ${total_jobs}${RST}"
  echo

  # groupSize sensitivity table
  say "${BOLD}groupSize sensitivity (total Condor jobs submitted):${RST}"
  echo
  printf "  ${BOLD}%-12s │ %-14s${RST}\n" "groupSize" "TOTAL JOBS"
  printf "  ─────────────┼───────────────\n"
  for _gs_try in 4 6 7 8 10; do
    local _bj=0
    while IFS= read -r _rn; do
      [[ -z "$_rn" || "$_rn" =~ ^# ]] && continue
      local _r8t; _r8t="$(run8 "$_rn")"
      local _lft="${LIST_DIR}/${LIST_PREFIX}-${_r8t}.list"
      [[ -s "$_lft" ]] || continue
      local _nft; _nft=$(wc -l < "$_lft" | awk '{print $1}')
      _bj=$(( _bj + $( ceil_div "$_nft" "$_gs_try" ) ))
    done < "$GOLDEN"
    local _tj=$(( n_matrix * _bj ))
    local _marker=""
    [[ "$_gs_try" -eq "$gs" ]] && _marker=" ${BOLD}← current${RST}"
    printf "  %-12s │ %14s%s\n" "$_gs_try" "$_tj" "$_marker"
  done
  echo

  say "Output tree: ${DEST_BASE}/<cfg_tag>/<run8>/*.root"
  local _ex_tag="jetMinPt$(sim_pt_tag "${ck_pts[0]}")_$(sim_b2b_tag "${ck_fracs[0]}")_$(sim_vz_tag "${ck_vzs[0]}")_$(sim_cone_tag "${ck_cones[0]}")_${iso_tags[0]}"
  (( uepipe_in_tag )) && _ex_tag="${_ex_tag}_${uepipe_modes[0]}"
  say "  e.g. ${DIM}${DEST_BASE}/${_ex_tag}/00067599/*.root${RST}"
}

# Submit a set of runs (from a round file or the whole golden list) to Condor
#   submit_condor <runs_source> [firstChunk]
submit_condor() {
  local source="$1"
  local first_chunk="${2:-}"

  [[ -s "$source" ]] || { err "Run source not found or empty: $source"; exit 5; }

  # Clean stale .sub files for this TAG
  rm -f "${SUB_DIR}/RecoilJets_${TAG}_"*.sub 2>/dev/null || true

  # -------------------------------------------------------------------
  # DATA output wipe (pp/auau/oo only): unconditional on every submit
  # -------------------------------------------------------------------
  if [[ "$DATASET" != "isSim" ]]; then
    [[ -n "${DEST_BASE:-}" ]] || { err "DEST_BASE is empty; refusing to wipe."; exit 60; }
    [[ "$DEST_BASE" != "/" ]] || { err "DEST_BASE is '/' ; refusing to wipe."; exit 60; }

    case "$DEST_BASE" in
      */thesisAna/pp|*/thesisAna/pp25|*/thesisAna/auau|*/thesisAna/oo) ;;
      */thesisAna/pp/*|*/thesisAna/pp25/*|*/thesisAna/auau/*|*/thesisAna/oo/*) ;;
      *)
        err "Refusing to wipe DEST_BASE='$DEST_BASE' (not an expected thesisAna/{pp|pp25|auau|oo} path)"
        exit 61
        ;;
    esac

    say "WIPING OUTPUT TREE (dataset=${DATASET}) → ${DEST_BASE}"
    mkdir -p "$DEST_BASE"
    find "$DEST_BASE" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true
  fi

  local stamp; stamp="$(date +%Y%m%d_%H%M%S)"
  local sub="${SUB_DIR}/RecoilJets_${TAG}_${stamp}.sub"
  local exe_to_use="${BULK_FROZEN_EXE:-${EXE}}"
  local macro_env=""
  [[ -n "${BULK_FROZEN_MACRO:-}" ]] && macro_env=";RJ_MACRO_PATH=${BULK_FROZEN_MACRO}"

  # Snapshot YAML at submit time so idle jobs are immune to later edits
  local yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  local yaml_snap="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${stamp}.yaml"
  local source_runs
  source_runs=$(grep -cE '^[0-9]+' "$source" 2>/dev/null || true)
  mkdir -p "$SIM_YAML_OVERRIDE_DIR"
  cp -f "$yaml_src" "$yaml_snap"
  say "YAML snapshot: ${yaml_snap}"
  say "Submit context: source=${source}  runs=${source_runs:-0}  groupSize=${GROUP_SIZE}  firstChunk=${first_chunk:-none}"
  say "Submit environment: RJ_DATASET=${DATASET}  RJ_VERBOSITY=0  RJ_CONFIG_YAML=${yaml_snap}${macro_env}"

  cat > "$sub" <<SUB
universe      = vanilla
executable    = ${exe_to_use}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/job.\$(Cluster).\$(Process).err
request_memory= 2000MB
should_transfer_files = NO
stream_output = True
stream_error  = True
# Force dataset & quiet macro on Condor (YAML frozen at submit time):
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=0;RJ_CONFIG_YAML=${yaml_snap}${macro_env}
SUB

  local queued=0
  local run_counter=0
  local t0; t0="$(date +%s)"

  # Controls:
  #   RJ_SUBMIT_TRACE=1               -> print one line per run
  #   RJ_SUBMIT_PROGRESS_EVERY=25     -> print progress every N runs (default 25)
  local submit_trace="${RJ_SUBMIT_TRACE:-0}"
  local progress_every="${RJ_SUBMIT_PROGRESS_EVERY:-25}"
  [[ "$progress_every" =~ ^[0-9]+$ ]] || progress_every=25
  (( progress_every > 0 )) || progress_every=25

  say "Building Condor submit file: $(basename "$sub")"
  say "Progress: export RJ_SUBMIT_TRACE=1 for per-run lines (or set RJ_SUBMIT_PROGRESS_EVERY=N)"

  while IFS= read -r r8_raw; do
    [[ -z "$r8_raw" || "$r8_raw" =~ ^# ]] && continue
    r8="$(run8 "$r8_raw")"   # ← zero-pad to 8 digits
    (( run_counter+=1 ))

    if (( submit_trace )); then
      say "[submit] run ${r8} (run #${run_counter})"
    elif (( run_counter % progress_every == 0 )); then
      local now elapsed
      now="$(date +%s)"
      elapsed=$(( now - t0 ))
      say "[submit] processed ${run_counter} runs, queued ${queued} jobs so far (elapsed ${elapsed}s)"
    fi

    # Optional trigger filter: require GL1 scaledown bit to be active
    if [[ -n "${TRIGGER_BIT}" ]]; then
      if ! is_trigger_active "$r8" "$TRIGGER_BIT"; then
        say "Skipping run ${r8} (trigger bit ${TRIGGER_BIT} not active)"
        continue
      fi
    fi

    # Build group lists for this run
    mapfile -t groups < <( make_groups "$r8" "$GROUP_SIZE" )

    (( ${#groups[@]} )) || { warn "No groups produced for run $r8; skipping"; continue; }

    if (( submit_trace )); then
      say "    files per job=${GROUP_SIZE} -> groups produced=${#groups[@]}"
    fi

    # If "firstChunk" is requested, only submit the first group (smoke test)
    if [[ "$first_chunk" == "firstChunk" ]]; then
      groups=( "${groups[0]}" )
      (( submit_trace )) && say "    firstChunk enabled -> submitting 1 group for run ${r8}"
    fi

    local gidx=0
    for glist in "${groups[@]}"; do
      (( gidx+=1 ))
      printf 'arguments = %s %s %s $(Cluster) 0 %d NONE %s\nqueue\n\n' \
             "$r8" "$glist" "$DATASET" "$gidx" "$DEST_BASE" >> "$sub"
      (( queued+=1 ))
    done

  done < "$source"

  local submit_elapsed
  submit_elapsed=$(( $(date +%s) - t0 ))
  say "Submitting ${BOLD}${queued}${RST} jobs → $(basename "$sub")"
  say "Submit summary: runsProcessed=${run_counter}  groupSize=${GROUP_SIZE}  firstChunk=${first_chunk:-none}  elapsed=${submit_elapsed}s"
  say "Wrapper path: ${exe_to_use}"
  need_cmd condor_submit
  condor_submit "$sub"
}

# ------------------------ Parse CLI ------------------------
[[ $# -ge 1 ]] || usage
resolve_dataset "$1"

# Parse remaining tokens (order-agnostic):
ACTION=""
DRYRUN=0
SIM_SAMPLE_EXPLICIT=0
GROUP_SIZE_EXPLICIT=0
tokens=( "${@:2}" )
for (( idx=0; idx<${#tokens[@]}; idx++ )); do
  tok="${tokens[$idx]}"
  case "$tok" in
    local|isLocalIsoPing|condor|splitGoldenRunList|condorTest|condorDoAll)
      ACTION="$tok"
      ;;
    groupSize)
      next=$((idx+1))
      [[ $next -lt ${#tokens[@]} ]] || { err "groupSize requires a value"; exit 2; }
      GROUP_SIZE="${tokens[$next]}"
      GROUP_SIZE_EXPLICIT=1
      idx=$next
      ;;
    maxJobs)
      next=$((idx+1))
      [[ $next -lt ${#tokens[@]} ]] || { err "maxJobs requires a value"; exit 2; }
      MAX_JOBS="${tokens[$next]}"
      idx=$next
      ;;
    TRIGGER=*)
      TRIGGER_BIT="${tok#TRIGGER=}"
      [[ "$TRIGGER_BIT" =~ ^[0-9]+$ ]] || { err "TRIGGER must be an integer bit index (e.g., TRIGGER=26)"; exit 2; }
      ;;
    SAMPLE=*)
      SIM_SAMPLE="${tok#SAMPLE=}"
      SIM_SAMPLE_EXPLICIT=1
      ;;
    CHECKJOBS)
      if [[ -n "$ACTION" ]]; then
        DRYRUN=1
      else
        ACTION="CHECKJOBS"
      fi
      ;;
    *)
      :  # ignore unrecognized tokens
      ;;
  esac
done

# Default action if none provided
if [[ -z "$ACTION" ]]; then
  if [[ "$IS_SIM" -eq 1 ]]; then
    ACTION="CHECKJOBS"
  else
    ACTION="condor"
  fi
fi

# If a trigger filter is requested, require psql
if [[ -n "${TRIGGER_BIT}" ]]; then
  need_cmd psql
fi

# Only print the full banner when we're actually running work (not CHECKJOBS)
if [[ "$ACTION" != "CHECKJOBS" ]]; then
  yaml_src_print="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  say "Action=${BOLD}${ACTION}${RST}  Dataset=${BOLD}${DATASET}${RST}  Tag=${TAG}"
  say "Executable=${EXE}"
  say "Macro=${MACRO}"
  say "YAML source=${yaml_src_print}"
  if [[ "$IS_SIM" -eq 1 ]]; then
    say "Sim sample default=${SIM_SAMPLE_DEFAULT}  current=${SIM_SAMPLE}  explicit=${SIM_SAMPLE_EXPLICIT}"
  else
    say "Input lists dir=${LIST_DIR}"
    say "Golden runs=${GOLDEN}"
  fi
  say "Stage dir=${STAGE_DIR}"
  say "Rounds dir=${ROUND_DIR}"
  say "Dest base=${DEST_BASE}"
  if [[ -n "${TRIGGER_BIT}" ]]; then
    say "Trigger filter      : bit=${TRIGGER_BIT} (scaledown != -1 required)"
  fi
  say "Group size=${GROUP_SIZE} (explicit=${GROUP_SIZE_EXPLICIT})  Max jobs/round=${MAX_JOBS}"
  say "Trace knobs         : RJ_GROUP_TRACE=${RJ_GROUP_TRACE:-0}  RJ_SUBMIT_TRACE=${RJ_SUBMIT_TRACE:-0}  RJ_SUBMIT_PROGRESS_EVERY=${RJ_SUBMIT_PROGRESS_EVERY:-25}"
  echo
fi

# ------------------------ Actions --------------------------
case "$ACTION" in
  CHECKJOBS)
    # Dry-run only
    if [[ "$IS_SIM" -eq 1 ]]; then
      check_jobs_sim
    else
      check_jobs_all
    fi
    exit 0
    ;;

  local)
    # Parse optional [Nevents] and VERBOSE=N from tokens after 'local'
    nevt="$LOCAL_EVENTS"
    RJV="10"                     # default for local runs
    user_nevt_set=0
    rest=( "${@:3}" )
    for t in "${rest[@]}"; do
      if [[ "$t" =~ ^VERBOSE=([0-9]+)$ ]]; then
        RJV="${BASH_REMATCH[1]}"
      elif [[ "$t" =~ ^[0-9]+$ ]]; then
        nevt="$t"
        user_nevt_set=1
      fi
    done

    if [[ "$IS_SIM" -eq 1 ]]; then
      # Default SIM local smoke test should be fast unless user explicitly sets Nevents
      if (( user_nevt_set == 0 )); then
        nevt="100"
      fi

      # Mirror condorDoAll grouping behavior for a faithful smoke test:
      # groupSize defaults to 7 unless explicitly provided by the user.
      gs_local="$GROUP_SIZE"
      if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
        gs_local="7"
      fi

      master_yaml="$(sim_yaml_master_path)"
      [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

      mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
      mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
      mapfile -t sim_vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
      mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
      (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
      (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
      (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
      (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
      build_iso_modes "$master_yaml"
      read_uepipe_modes "$master_yaml" "$TAG"

      samples=()
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
        case "$DATASET" in
          isSimEmbedded)          samples=( "run28_embeddedPhoton10" "run28_embeddedPhoton20" ) ;;
          isSimEmbeddedInclusive) samples=( "run28_embeddedJet10" "run28_embeddedJet20" ) ;;
          isSimJet5)              samples=( "run28_jet5" ) ;;
          isSimMB)                samples=( "run28_detroit" ) ;;
          *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
        esac
      else
        samples=( "${SIM_SAMPLE}" )
      fi

      SIM_DEST_BASE_RESOLVED="$DEST_BASE"

      say "SIM local smoke test (mirrors condorDoAll matrix)"
      say "  YAML master : ${master_yaml}"
      say "  groupSize   : ${gs_local} (condorDoAll default)"
      say "  events      : ${nevt} (default for SIM local)"
      say "  samples     : ${samples[*]}"
      say "  dest base   : ${SIM_DEST_BASE_RESOLVED}"
      echo

      for pt in "${sim_pts[@]}"; do
        for frac in "${sim_fracs[@]}"; do
          for vz in "${sim_vzs[@]}"; do
          for cone in "${sim_cones[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
          for uepipe in "${uepipe_modes[@]}"; do
          SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$pt")_$(sim_b2b_tag "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso_tags[$iso_idx]}"
          (( uepipe_in_tag )) && SIM_CFG_TAG="${SIM_CFG_TAG}_${uepipe}"
          DEST_BASE="${SIM_DEST_BASE_RESOLVED}/${SIM_CFG_TAG}"
          yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt" "$frac" "$vz" "$cone" "${iso_sliding[$iso_idx]}" "${iso_fixed[$iso_idx]}" "${uepipe}" "$SIM_CFG_TAG" "LOCAL")"

          for samp in "${samples[@]}"; do
            SIM_SAMPLE="$samp"
            GROUP_SIZE="$gs_local"

            # For local mode we only need the first input line — skip creating all group files
            sim_init
            say "  [local] sample=${SIM_SAMPLE} ($(wc -l < "$SIM_CLEAN_LIST" | tr -d ' ') entries) — extracting first line, skipping full grouping"

            tmp="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_LOCAL_firstfile_grp001.list"
            head -n 1 "$SIM_CLEAN_LIST" > "$tmp"
            [[ -s "$tmp" ]] || { err "No sim entries (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"; exit 30; }

            in_line="$(head -n 1 "$tmp" 2>/dev/null || true)"
            chunk_base="$(basename "$tmp")"
            chunk_tag="${chunk_base%.list}"
            out_root_preview="${DEST_BASE}/${SIM_SAMPLE}/RecoilJets_${DATASET}_${chunk_tag}.root"

            say "----------------------------------------"
            say "SIM local: tag=${SIM_CFG_TAG}  sample=${SIM_SAMPLE}"
            say "  jet_pt_min=${pt}  back_to_back_pi_fraction=${frac}  vz_cut_cm=${vz}  coneR=${cone}  iso=${iso_tags[$iso_idx]}  uepipe=${uepipe}"
            say "  YAML override: ${yaml_override}"
            say "  grp001 list   : ${tmp}"
            say "  input line    : ${in_line}"
            say "  temp list     : ${tmp}"
            say "  out ROOT      : ${out_root_preview}"
            say "  wrapper env   : RJ_VERBOSITY=${RJV} RJ_CONFIG_YAML=${yaml_override}"
            say "  wrapper args  : sample=${SIM_SAMPLE} dataset=${DATASET} mode=LOCAL nevents=${nevt} chunk=1 dest=${DEST_BASE}"
            find "${DEST_BASE}/${SIM_SAMPLE}" -type f -name "*_LOCAL_*.root" -delete 2>/dev/null || true
            say "Invoking wrapper locally…"

            RJ_VERBOSITY="$RJV" RJ_CONFIG_YAML="$yaml_override" bash "$EXE" "$SIM_SAMPLE" "$tmp" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
            echo
          done
          done
          done
          done
          done
        done
      done
    else
      say "Local test on ${DATASET}  (events=${nevt}, RJ_VERBOSITY=${RJV})"

      [[ -s "$GOLDEN" ]] || { err "Golden list is empty: $GOLDEN"; exit 6; }

      # Pick first golden run; if TRIGGER_BIT set, pick the first run with that bit active
      r8=""
      while IFS= read -r rn; do
        [[ -z "$rn" || "$rn" =~ ^# ]] && continue
        cand="$(run8 "$rn")"
        if [[ -n "${TRIGGER_BIT}" ]]; then
          if is_trigger_active "$cand" "$TRIGGER_BIT"; then r8="$cand"; break; fi
        else
          r8="$cand"; break
        fi
      done < "$GOLDEN"

      [[ -n "$r8" ]] || { err "No run in $GOLDEN satisfies the requested trigger filter (TRIGGER=${TRIGGER_BIT:-none})."; exit 6; }

      src="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
      [[ -s "$src" ]] || { err "Per-run list missing or empty: $src"; exit 7; }

      # Build a 1-file list (first file only)
      tmp="${STAGE_DIR}/run${r8}_LOCAL_firstfile.list"
      head -n 1 "$src" > "$tmp"
      say "Temp list → $tmp"

      # Local debug helpers (safe defaults):
      #   RJ_CRASH_BACKTRACE=1 prints a stack trace on SIGSEGV/SIGABRT.
      #   RJ_F4A_VERBOSE=1 enables Fun4AllServer subsystem-level tracing (chatty).
      #   RJ_STEP_EVENTS=1 runs 1 event at a time with per-event banners (only for small nevt).
      RJ_CRASH_BACKTRACE_LOCAL=1
      RJ_F4A_VERBOSE_LOCAL=0
      RJ_STEP_EVENTS_LOCAL=0
      if [[ "$nevt" =~ ^[0-9]+$ ]] && (( nevt > 0 && nevt <= 50 )); then
        RJ_F4A_VERBOSE_LOCAL=1
        RJ_STEP_EVENTS_LOCAL=1
      fi

      # Read YAML matrix arrays for DATA local
      data_yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
      mapfile -t data_pts   < <( yaml_get_values "jet_pt_min" "$data_yaml_src" )
      mapfile -t data_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$data_yaml_src" )
      mapfile -t data_vzs   < <( yaml_get_values "vz_cut_cm" "$data_yaml_src" )
      mapfile -t data_cones < <( yaml_get_values "coneR" "$data_yaml_src" )
      (( ${#data_pts[@]} ))   || { err "No values found for jet_pt_min in $data_yaml_src"; exit 72; }
      (( ${#data_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $data_yaml_src"; exit 72; }
      (( ${#data_vzs[@]} ))   || { err "No values found for vz_cut_cm in $data_yaml_src"; exit 72; }
      (( ${#data_cones[@]} )) || { err "No values found for coneR in $data_yaml_src"; exit 72; }
      build_iso_modes "$data_yaml_src"
      read_uepipe_modes "$data_yaml_src" "$TAG"
      DATA_DEST_BASE_SAVED="$DEST_BASE"

      for data_pt in "${data_pts[@]}"; do
      for data_frac in "${data_fracs[@]}"; do
      for data_vz in "${data_vzs[@]}"; do
        for data_cone in "${data_cones[@]}"; do
        for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
        for uepipe in "${uepipe_modes[@]}"; do
        dpt_tag="jetMinPt$(sim_pt_tag "$data_pt")"
        dfrac_tag="$(sim_b2b_tag "$data_frac")"
        dvz_tag="$(sim_vz_tag "$data_vz")"
        dcone_tag="$(sim_cone_tag "$data_cone")"
        data_cfg_tag="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_tags[$iso_idx]}"
        (( uepipe_in_tag )) && data_cfg_tag="${data_cfg_tag}_${uepipe}"
        yaml_override="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${data_cfg_tag}_LOCAL.yaml"
        mkdir -p "$SIM_YAML_OVERRIDE_DIR"
        sed -E \
          -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${data_pt}|" \
          -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${data_frac}|" \
          -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${data_vz}|" \
          -e "s|^([[:space:]]*coneR:).*|\\1 ${data_cone}|" \
          -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${iso_sliding[$iso_idx]}|" \
          -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${iso_fixed[$iso_idx]}|" \
          -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${uepipe}|" \
          "$data_yaml_src" > "$yaml_override"
        DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"

        say "----------------------------------------"
        say "DATA local: pt=${data_pt}  frac=${data_frac}  vz=${data_vz}  coneR=${data_cone}  iso=${iso_tags[$iso_idx]}  uepipe=${uepipe}  tag=${data_cfg_tag}"
        say "  YAML override: ${yaml_override}"
        say "  DEST_BASE    : ${DEST_BASE}"
        say "  wrapper env  : RJ_DATASET=${DATASET} RJ_VERBOSITY=${RJV} RJ_CONFIG_YAML=${yaml_override} RJ_CRASH_BACKTRACE=${RJ_CRASH_BACKTRACE_LOCAL} RJ_F4A_VERBOSE=${RJ_F4A_VERBOSE_LOCAL} RJ_STEP_EVENTS=${RJ_STEP_EVENTS_LOCAL}"
        say "  wrapper args : run=${r8} dataset=${DATASET} mode=LOCAL nevents=${nevt} chunk=1 dest=${DEST_BASE}"
        find "${DEST_BASE}" -type f -name "*_LOCAL_*.root" -delete 2>/dev/null || true
        say "Invoking wrapper locally…"

        RJ_DATASET="$DATASET" RJ_VERBOSITY="$RJV" \
        RJ_CONFIG_YAML="$yaml_override" \
        RJ_CRASH_BACKTRACE="$RJ_CRASH_BACKTRACE_LOCAL" \
        RJ_F4A_VERBOSE="$RJ_F4A_VERBOSE_LOCAL" \
        RJ_STEP_EVENTS="$RJ_STEP_EVENTS_LOCAL" \
        bash "$EXE" "$r8" "$tmp" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
        echo
        done
        done
        done
      done
      done
      done
    fi
    ;;

  isLocalIsoPing)
    [[ "$IS_SIM" -eq 0 ]] || { err "isLocalIsoPing is DATA-only"; exit 2; }
    [[ "$DATASET" == "isAuAu" ]] || { err "isLocalIsoPing is implemented only for Au+Au"; exit 2; }

    need_cmd psql

    # Sensible default:
    #   200 inclusive photons per centrality bin gives stable means/medians/fractions
    #   without forcing an excessively long local run.
    ISO_PING_TARGET_PER_CENT="${RJ_ISO_AUDIT_TARGET_PER_CENT:-200}"
    [[ "$ISO_PING_TARGET_PER_CENT" =~ ^[0-9]+$ ]] || { err "RJ_ISO_AUDIT_TARGET_PER_CENT must be a positive integer"; exit 2; }
    (( ISO_PING_TARGET_PER_CENT > 0 )) || { err "RJ_ISO_AUDIT_TARGET_PER_CENT must be > 0"; exit 2; }

    ISO_PING_REQUIRED_CENT_BINS="${RJ_ISO_AUDIT_REQUIRED_CENT_BINS:-3}"
    [[ "$ISO_PING_REQUIRED_CENT_BINS" =~ ^[0-9]+$ ]] || { err "RJ_ISO_AUDIT_REQUIRED_CENT_BINS must be a positive integer"; exit 2; }
    (( ISO_PING_REQUIRED_CENT_BINS > 0 )) || { err "RJ_ISO_AUDIT_REQUIRED_CENT_BINS must be > 0"; exit 2; }

    ISO_PING_TRIGGER_BIT=22
    ISO_PING_TRIGGER_KEY="photon_10_plus_MBD_NS_geq_2_vtx_lt_150"

    RJV="1"
    rest=( "${@:3}" )
    for t in "${rest[@]}"; do
      if [[ "$t" =~ ^VERBOSE=([0-9]+)$ ]]; then
        RJV="${BASH_REMATCH[1]}"
      fi
    done

    data_yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
    [[ -s "$data_yaml_src" ]] || { err "YAML source missing: $data_yaml_src"; exit 72; }

    mapfile -t data_pts   < <( yaml_get_values "jet_pt_min" "$data_yaml_src" )
    mapfile -t data_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$data_yaml_src" )
    mapfile -t data_vzs   < <( yaml_get_values "vz_cut_cm" "$data_yaml_src" )
    mapfile -t data_cones < <( yaml_get_values "coneR" "$data_yaml_src" )
    (( ${#data_pts[@]} ))   || { err "No values found for jet_pt_min in $data_yaml_src"; exit 72; }
    (( ${#data_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $data_yaml_src"; exit 72; }
    (( ${#data_vzs[@]} ))   || { err "No values found for vz_cut_cm in $data_yaml_src"; exit 72; }
    (( ${#data_cones[@]} )) || { err "No values found for coneR in $data_yaml_src"; exit 72; }

    build_iso_modes "$data_yaml_src"
    read_uepipe_modes "$data_yaml_src" "$TAG"

    selected_uepipe="${RJ_CLUSTER_UEPIPELINE:-}"
    selected_uepipe="$(trim_ws "$selected_uepipe")"
    [[ -n "$selected_uepipe" ]] || selected_uepipe="${uepipe_modes[0]}"

    case "$selected_uepipe" in
      noSub|baseVariant|variantA|variantB) ;;
      *) err "Unsupported RJ_CLUSTER_UEPIPELINE for isLocalIsoPing: ${selected_uepipe}"; exit 2 ;;
    esac

    pt0="${data_pts[0]}"
    frac0="${data_fracs[0]}"
    vz0="${data_vzs[0]}"
    cone0="${data_cones[0]}"

    iso_idx=0
    dpt_tag="jetMinPt$(sim_pt_tag "$pt0")"
    dfrac_tag="$(sim_b2b_tag "$frac0")"
    dvz_tag="$(sim_vz_tag "$vz0")"
    dcone_tag="$(sim_cone_tag "$cone0")"
    data_cfg_tag="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_tags[$iso_idx]}_${selected_uepipe}"

    yaml_override="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${data_cfg_tag}_ISOPING.yaml"
    mkdir -p "$SIM_YAML_OVERRIDE_DIR"
    sed -E \
      -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${pt0}|" \
      -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${frac0}|" \
      -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${vz0}|" \
      -e "s|^([[:space:]]*coneR:).*|\\1 ${cone0}|" \
      -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${iso_sliding[$iso_idx]}|" \
      -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${iso_fixed[$iso_idx]}|" \
      -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${selected_uepipe}|" \
      "$data_yaml_src" > "$yaml_override"

    r8="$(pick_first_iso_ping_run "$ISO_PING_TRIGGER_BIT")"
    [[ -n "$r8" ]] || { err "No Au+Au run with a non-empty per-run list and active ${ISO_PING_TRIGGER_KEY} trigger"; exit 6; }

    mapfile -t groups < <( make_groups "$r8" "$GROUP_SIZE" )
    (( ${#groups[@]} )) || { err "No chunk lists were produced for run ${r8}"; exit 9; }

    combined="${STAGE_DIR}/run${r8}_isLocalIsoPing_allchunks.list"
    : > "$combined"
    for grp in "${groups[@]}"; do
      cat "$grp" >> "$combined"
    done

    first_group="$(basename "${groups[0]}")"
    last_group="$(basename "${groups[$(( ${#groups[@]} - 1 ))]}")"
    chunk_span="${first_group} .. ${last_group}"

    tower_prefix="${RJ_TOWERINFO_PREFIX:-TOWERINFO_CALIB}"
    audit_photon_input_node="CLUSTERINFO_CEMC"
    audit_photon_builder_is_auau="false"
    audit_pcb_em_node="${tower_prefix}_CEMC"
    audit_pcb_hi_node="${tower_prefix}_HCALIN"
    audit_pcb_ho_node="${tower_prefix}_HCALOUT"

    case "$selected_uepipe" in
      variantA|variantB)
        audit_photon_input_node="CLUSTERINFO_CEMC_PHOSUB"
        audit_pcb_em_node="${tower_prefix}_CEMC_PHOSUB"
        audit_pcb_hi_node="${tower_prefix}_HCALIN_SUB1"
        audit_pcb_ho_node="${tower_prefix}_HCALOUT_SUB1"
        ;;
      baseVariant)
        audit_photon_builder_is_auau="true"
        ;;
    esac

    local_dest="${BASE}/iso_ping_local/${TAG}/${data_cfg_tag}"
    mkdir -p "$local_dest"

    say "IsolationAudit local ping"
    say "  run selection   : first sorted Au+Au run with list + active ${ISO_PING_TRIGGER_KEY} (bit ${ISO_PING_TRIGGER_BIT})"
    say "  selected run    : ${r8}"
    say "  chunk span      : ${chunk_span}  (nChunks=${#groups[@]}, groupSize=${GROUP_SIZE})"
    say "  YAML override   : ${yaml_override}"
    say "  path class      : ${selected_uepipe}"
    say "  target / cent   : ${ISO_PING_TARGET_PER_CENT}"
    say "  required bins   : first ${ISO_PING_REQUIRED_CENT_BINS}"
    say "  combined list   : ${combined}"
    say "  output dir      : ${local_dest}"
    say "Invoking wrapper locally…"

    RJ_DATASET="$DATASET" \
    RJ_VERBOSITY="$RJV" \
    RJ_CONFIG_YAML="$yaml_override" \
    RJ_CLUSTER_UEPIPELINE="$selected_uepipe" \
    RJ_ISO_AUDIT_MODE="1" \
    RJ_ISO_AUDIT_TARGET_PER_CENT="$ISO_PING_TARGET_PER_CENT" \
    RJ_ISO_AUDIT_REQUIRED_CENT_BINS="$ISO_PING_REQUIRED_CENT_BINS" \
    RJ_ISO_AUDIT_DATASET_TOKEN="$TAG" \
    RJ_ISO_AUDIT_SELECTED_RUN="$r8" \
    RJ_ISO_AUDIT_CHUNK_SPAN="$chunk_span" \
    RJ_ISO_AUDIT_GROUP_SIZE="$GROUP_SIZE" \
    RJ_ISO_AUDIT_COMBINED_LIST="$combined" \
    RJ_ISO_AUDIT_CLUSTER_UEPIPELINE="$selected_uepipe" \
    RJ_ISO_AUDIT_PATH_CLASS="$selected_uepipe" \
    RJ_ISO_AUDIT_PHOTON_INPUT_CLUSTER_NODE="$audit_photon_input_node" \
    RJ_ISO_AUDIT_PHOTON_BUILDER_IS_AUAU="$audit_photon_builder_is_auau" \
    RJ_ISO_AUDIT_PCB_TOWER_PREFIX="$tower_prefix" \
    RJ_ISO_AUDIT_PCB_EM_NODE="$audit_pcb_em_node" \
    RJ_ISO_AUDIT_PCB_HI_NODE="$audit_pcb_hi_node" \
    RJ_ISO_AUDIT_PCB_HO_NODE="$audit_pcb_ho_node" \
    RJ_CRASH_BACKTRACE="1" \
    RJ_F4A_VERBOSE="0" \
    RJ_STEP_EVENTS="0" \
    bash "$EXE" "$r8" "$combined" "$DATASET" LOCAL 0 1 NONE "$local_dest"
    ;;

  condorTest)
    # One verbose Condor job on first sim file
    [[ "$IS_SIM" -eq 1 ]] || { err "condorTest is only valid for isSim variants"; exit 2; }
    need_cmd condor_submit

    master_yaml="$(sim_yaml_master_path)"
    [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

    mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
    mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
    mapfile -t sim_vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
    mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
    build_iso_modes "$master_yaml"

    SIM_DEST_BASE_RESOLVED="$DEST_BASE"

    pt0="${sim_pts[0]}"
    frac0="${sim_fracs[0]}"
    vz0="${sim_vzs[0]}"
    cone0="${sim_cones[0]}"
    SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$pt0")_$(sim_b2b_tag "$frac0")_$(sim_vz_tag "$vz0")_$(sim_cone_tag "$cone0")_${iso_tags[0]}"
    DEST_BASE="${SIM_DEST_BASE_RESOLVED}/${SIM_CFG_TAG}"
    stamp="$(date +%Y%m%d_%H%M%S)"
    yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt0" "$frac0" "$vz0" "$cone0" "${iso_sliding[0]}" "${iso_fixed[0]}" "noSub" "$SIM_CFG_TAG" "$stamp")"

    sim_init

    tmp="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_condorTest_firstfile.list"
    head -n 1 "$SIM_CLEAN_LIST" > "$tmp"

    # Clean stale .sub files for SIM
    rm -f "${SUB_DIR}/RecoilJets_sim_"*.sub 2>/dev/null || true
    sub="${SUB_DIR}/RecoilJets_sim_${SIM_CFG_TAG}_${SIM_SAMPLE}_${stamp}_TEST.sub"

    cat > "$sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).err
request_memory= 2000MB
should_transfer_files = NO
stream_output = True
stream_error  = True
environment   = RJ_VERBOSITY=10;RJ_CONFIG_YAML=${yaml_override}
arguments     = ${SIM_SAMPLE} ${tmp} ${DATASET} \$(Cluster) 0 1 NONE ${DEST_BASE}
queue
SUB

    say "Submitting 1 ${DATASET} condorTest job (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG}) → $(basename "$sub")"
    say "Output ROOT dir: ${DEST_BASE}/${SIM_SAMPLE}"
    say "YAML override: ${yaml_override}"
    condor_submit "$sub"
    ;;

  condorDoAll)
    # Submit all sim files in grouped chunks. MUST wipe outputs/logs first.
    [[ "$IS_SIM" -eq 1 ]] || { err "condorDoAll is only valid for isSim variants"; exit 2; }
    gs_doall="$GROUP_SIZE"
    if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
      gs_doall="7"
    fi

    master_yaml="$(sim_yaml_master_path)"
    [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

    mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
    mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
    mapfile -t sim_vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
    mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
    build_iso_modes "$master_yaml"
    read_uepipe_modes "$master_yaml" "$TAG"

    samples=()
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
      case "$DATASET" in
        isSimEmbedded)          samples=( "run28_embeddedPhoton10" "run28_embeddedPhoton20" ) ;;
        isSimEmbeddedInclusive) samples=( "run28_embeddedJet10" "run28_embeddedJet20" ) ;;
        isSimJet5)              samples=( "run28_jet5" ) ;;
        isSimMB)                samples=( "run28_detroit" ) ;;
        *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
      esac
    else
      samples=( "${SIM_SAMPLE}" )
    fi

    # If CHECKJOBS was also provided, do a dry-run count for the FULL condorDoAll matrix and exit.
    if [[ "${DRYRUN:-0}" -eq 1 ]]; then
      n_cfg=$(( ${#sim_pts[@]} * ${#sim_fracs[@]} * ${#sim_vzs[@]} * ${#sim_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
      per_cfg_jobs=0

      say "CHECKJOBS (${DATASET} condorDoAll matrix)"
      say "  YAML master         : ${master_yaml}"
      say "  groupSize (baseline): ${gs_doall}"
      echo
      say "${BOLD}Matrix dimensions:${RST}"
      say "  jet_pt_min                        : [${sim_pts[*]}]  (${#sim_pts[@]} values)"
      say "  back_to_back_dphi_min_pi_fraction : [${sim_fracs[*]}]  (${#sim_fracs[@]} values)"
      say "  vz_cut_cm                         : [${sim_vzs[*]}]  (${#sim_vzs[@]} values)"
      say "  coneR                             : [${sim_cones[*]}]  (${#sim_cones[@]} values)"
      say "  iso modes                         : [${iso_tags[*]}]  (${#iso_tags[@]} values)"
      say "  clusterUEpipeline                 : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} values)"
      say "  cfg combos (product)              : ${BOLD}${n_cfg}${RST}"
      say "  samples                           : [${samples[*]}]  (${#samples[@]} values)"
      echo

      say "${BOLD}Per-sample input file counts:${RST}"
      for samp in "${samples[@]}"; do
        SIM_SAMPLE="$samp"
        sim_init
        nfiles=$(wc -l < "$SIM_CLEAN_LIST" | awk '{print $1}')
        njobs=$(ceil_div "$nfiles" "$gs_doall")
        say "  sample=${BOLD}${SIM_SAMPLE}${RST}  input_files=${nfiles}  jobs_per_cfg=${njobs}"
        per_cfg_jobs=$(( per_cfg_jobs + njobs ))
      done
      echo

      say "${BOLD}Full cfg tag list (${n_cfg} entries × ${#samples[@]} samples = ${BOLD}$((n_cfg * ${#samples[@]}))${RST} submit blocks):${RST}"
      cfg_num=0
      for pt in "${sim_pts[@]}"; do
        for frac in "${sim_fracs[@]}"; do
          for vz in "${sim_vzs[@]}"; do
          for cone in "${sim_cones[@]}"; do
          for (( _ci=0; _ci<${#iso_tags[@]}; _ci++ )); do
          for uepipe in "${uepipe_modes[@]}"; do
            (( cfg_num+=1 ))
            _tag="jetMinPt$(sim_pt_tag "$pt")_$(sim_b2b_tag "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso_tags[$_ci]}"
            (( uepipe_in_tag )) && _tag="${_tag}_${uepipe}"
            printf "  ${DIM}%3d${RST} │ %-70s │ pt=%-5s frac=%-6s vz=%-5s cone=%-5s iso=%-16s uepipe=%s\n" \
                   "$cfg_num" "$_tag" "$pt" "$frac" "$vz" "$cone" "${iso_tags[$_ci]}" "$uepipe"
          done
          done
          done
          done
        done
      done
      echo

      total_jobs=$(( n_cfg * per_cfg_jobs ))
      say "${BOLD}Job count summary:${RST}"
      say "  cfg combos             : ${n_cfg}"
      say "  jobs per combo (Σsamp) : ${per_cfg_jobs}"
      say "  ─────────────────────────────────"
      say "  ${BOLD}TOTAL CONDOR JOBS        : ${BOLD}${total_jobs}${RST}"
      echo
      say "Output tree: each tag becomes a subdirectory under ${DEST_BASE}/"
      say "  e.g. ${DIM}${DEST_BASE}/jetMinPt5_7pi_8_vz30_isoR30_fixedIso2GeV/<sample>/*.root${RST}"
      exit 0
    fi

    SIM_DEST_BASE_RESOLVED="$DEST_BASE"

    need_cmd condor_submit
    doall_stamp="$(date +%Y%m%d_%H%M%S)"
    # Freeze pipeline for this bulk submission
    cleanup_bulk_snapshots_for_tag
    case "$DATASET" in
      isSimEmbedded|isSimEmbeddedInclusive) create_pipeline_snapshot "auau" "$doall_stamp" ;;
      *)                                    create_pipeline_snapshot "pp" "$doall_stamp" ;;
    esac
    # Clean stale .sub files and YAML overrides for SIM
    rm -f "${SUB_DIR}/RecoilJets_sim_"*.sub "${SUB_DIR}/RecoilJets_${TAG}_"*.sub 2>/dev/null || true
    rm -f "${SIM_YAML_OVERRIDE_DIR}/analysis_config_jetMinPt"*.yaml "${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_"*.yaml 2>/dev/null || true
    for pt in "${sim_pts[@]}"; do
      for frac in "${sim_fracs[@]}"; do
        for vz in "${sim_vzs[@]}"; do
        for cone in "${sim_cones[@]}"; do
        for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
        for uepipe in "${uepipe_modes[@]}"; do
        SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$pt")_$(sim_b2b_tag "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso_tags[$iso_idx]}"
        (( uepipe_in_tag )) && SIM_CFG_TAG="${SIM_CFG_TAG}_${uepipe}"
        DEST_BASE="${SIM_DEST_BASE_RESOLVED}/${SIM_CFG_TAG}"
        yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt" "$frac" "$vz" "$cone" "${iso_sliding[$iso_idx]}" "${iso_fixed[$iso_idx]}" "${uepipe}" "$SIM_CFG_TAG" "$doall_stamp")"

        for samp in "${samples[@]}"; do
          SIM_SAMPLE="$samp"
          GROUP_SIZE="$gs_doall"

          wipe_sim_artifacts

          mapfile -t groups < <( make_sim_groups "$GROUP_SIZE" )
          (( ${#groups[@]} )) || { err "No sim groups produced (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"; exit 30; }

          stamp="$(date +%Y%m%d_%H%M%S)"
          sub="${SUB_DIR}/RecoilJets_sim_${SIM_CFG_TAG}_${SIM_SAMPLE}_${stamp}.sub"
          exe_for_sub="${BULK_FROZEN_EXE:-${EXE}}"
          macro_env_for_sub=""
          [[ -n "${BULK_FROZEN_MACRO:-}" ]] && macro_env_for_sub=";RJ_MACRO_PATH=${BULK_FROZEN_MACRO}"

          cat > "$sub" <<SUB
universe      = vanilla
executable    = ${exe_for_sub}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).err
request_memory= 2000MB
should_transfer_files = NO
stream_output = True
stream_error  = True
environment   = RJ_VERBOSITY=0;RJ_CONFIG_YAML=${yaml_override}${macro_env_for_sub}
SUB

          gidx=0
          for glist in "${groups[@]}"; do
            (( gidx+=1 ))
            printf 'arguments = %s %s %s $(Cluster) 0 %d NONE %s\nqueue\n\n' \
                   "$SIM_SAMPLE" "$glist" "$DATASET" "$gidx" "$DEST_BASE" >> "$sub"
          done

          say "Submitting ${DATASET} condorDoAll (tag=${SIM_CFG_TAG}, sample=${SIM_SAMPLE}, groupSize=${GROUP_SIZE}, uepipe=${uepipe}) → jobs=${BOLD}${#groups[@]}${RST}"
          say "Output ROOT dir: ${DEST_BASE}/${SIM_SAMPLE}"
          say "YAML override: ${yaml_override}"
          condor_submit "$sub"
        done
        done
        done
        done
        done
      done
    done
    ;;

  splitGoldenRunList)
    [[ "$IS_SIM" -eq 1 ]] && { err "splitGoldenRunList is not used in isSim mode"; exit 2; }
    split_golden "$GROUP_SIZE" "$MAX_JOBS"
    ;;

  condor)
    [[ "$IS_SIM" -eq 1 ]] && { err "Use: isSim condorTest | isSim condorDoAll (not 'condor')"; exit 2; }

    # Read YAML matrix arrays for DATA
    data_yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
    mapfile -t data_pts   < <( yaml_get_values "jet_pt_min" "$data_yaml_src" )
    mapfile -t data_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$data_yaml_src" )
    mapfile -t data_vzs   < <( yaml_get_values "vz_cut_cm" "$data_yaml_src" )
    mapfile -t data_cones < <( yaml_get_values "coneR" "$data_yaml_src" )
    (( ${#data_pts[@]} ))   || { err "No values found for jet_pt_min in $data_yaml_src"; exit 72; }
    (( ${#data_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $data_yaml_src"; exit 72; }
    (( ${#data_vzs[@]} ))   || { err "No values found for vz_cut_cm in $data_yaml_src"; exit 72; }
    (( ${#data_cones[@]} )) || { err "No values found for coneR in $data_yaml_src"; exit 72; }
    build_iso_modes "$data_yaml_src"
    read_uepipe_modes "$data_yaml_src" "$TAG"
    DATA_DEST_BASE_SAVED="$DEST_BASE"

    # DATA ONLY: condor testJob | condor round <N> [firstChunk] | condor all
    submode="${3:-}"
    case "$submode" in
      testJob)
        [[ -s "$GOLDEN" ]] || { err "Golden list is empty: $GOLDEN"; exit 6; }
        r8=""
        while IFS= read -r rn; do
          [[ -z "$rn" || "$rn" =~ ^# ]] && continue
          cand="$(run8 "$rn")"
          if [[ -n "${TRIGGER_BIT}" ]]; then
            if is_trigger_active "$cand" "$TRIGGER_BIT"; then r8="$cand"; break; fi
          else
            r8="$cand"; break
          fi
        done < "$GOLDEN"
        [[ -n "$r8" ]] || { err "No run in $GOLDEN satisfies the requested trigger filter (TRIGGER=${TRIGGER_BIT:-none})."; exit 6; }

        mapfile -t groups < <( make_groups "$r8" 1 )
        (( ${#groups[@]} )) || { err "No groups produced for run $r8"; exit 9; }
        glist="${groups[0]}"

        # Use first values for testJob (mirrors SIM condorTest)
        pt0="${data_pts[0]}"
        frac0="${data_fracs[0]}"
        vz0="${data_vzs[0]}"
        cone0="${data_cones[0]}"
        dpt_tag="jetMinPt$(sim_pt_tag "$pt0")"
        dfrac_tag="$(sim_b2b_tag "$frac0")"
        dvz_tag="$(sim_vz_tag "$vz0")"
        dcone_tag="$(sim_cone_tag "$cone0")"
        data_cfg_tag="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_tags[0]}"
        (( uepipe_in_tag )) && data_cfg_tag="${data_cfg_tag}_${uepipe_modes[0]}"

        stamp="$(date +%Y%m%d_%H%M%S)"
        sub="${SUB_DIR}/RecoilJets_${TAG}_${data_cfg_tag}_${stamp}_TEST.sub"

        # Snapshot YAML at submit time, pinning all matrix axes to first values
        yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
        yaml_snap="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${data_cfg_tag}_${stamp}_TEST.yaml"
        mkdir -p "$SIM_YAML_OVERRIDE_DIR"
        sed -E \
          -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${pt0}|" \
          -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${frac0}|" \
          -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${vz0}|" \
          -e "s|^([[:space:]]*coneR:).*|\\1 ${cone0}|" \
          -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${iso_sliding[0]}|" \
          -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${iso_fixed[0]}|" \
          -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${uepipe_modes[0]}|" \
          "$yaml_src" > "$yaml_snap"
        say "YAML snapshot (pt=${pt0}, frac=${frac0}, vz=${vz0}, coneR=${cone0}, iso=${iso_tags[0]}, uepipe=${uepipe_modes[0]}): ${yaml_snap}"
        DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"

        cat > "$sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/job.\$(Cluster).\$(Process).err
request_memory= 2000MB
should_transfer_files = NO
stream_output = True
stream_error  = True
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=10;RJ_CONFIG_YAML=${yaml_snap}
arguments     = ${r8} ${glist} ${DATASET} \$(Cluster) 0 1 NONE ${DEST_BASE}
queue
SUB
        say "Submitting 1 test job on run ${BOLD}${r8}${RST} (first chunk, groupSize=1, vz=${vz0}) → $(basename "$sub")"
        condor_submit "$sub"
        ;;
      round)
        seg="${4:?round number required}"
        firstChunk="${5:-}"
        round_file="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
        [[ -s "$round_file" ]] || { err "Round file not found: $round_file. Run 'splitGoldenRunList' first."; exit 8; }

        # Clean stale YAML overrides before fresh submission
        rm -f "${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_"*.yaml 2>/dev/null || true

        for data_pt in "${data_pts[@]}"; do
        for data_frac in "${data_fracs[@]}"; do
        for data_vz in "${data_vzs[@]}"; do
          for data_cone in "${data_cones[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
          for uepipe in "${uepipe_modes[@]}"; do
          dpt_tag="jetMinPt$(sim_pt_tag "$data_pt")"
          dfrac_tag="$(sim_b2b_tag "$data_frac")"
          dvz_tag="$(sim_vz_tag "$data_vz")"
          dcone_tag="$(sim_cone_tag "$data_cone")"
          data_cfg_tag="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_tags[$iso_idx]}"
          (( uepipe_in_tag )) && data_cfg_tag="${data_cfg_tag}_${uepipe}"
          yaml_override="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${data_cfg_tag}.yaml"
          mkdir -p "$SIM_YAML_OVERRIDE_DIR"
          sed -E \
            -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${data_pt}|" \
            -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${data_frac}|" \
            -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${data_vz}|" \
            -e "s|^([[:space:]]*coneR:).*|\\1 ${data_cone}|" \
            -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${iso_sliding[$iso_idx]}|" \
            -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${iso_fixed[$iso_idx]}|" \
            -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${uepipe}|" \
            "$data_yaml_src" > "$yaml_override"
          export RJ_CONFIG_YAML="$yaml_override"
          DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"

          say "DATA condor round (pt=${data_pt}, frac=${data_frac}, vz=${data_vz}, coneR=${data_cone}, iso=${iso_tags[$iso_idx]}, uepipe=${uepipe}, tag=${data_cfg_tag})"
          say "  YAML override: ${yaml_override}"
          say "  DEST_BASE    : ${DEST_BASE}"
          submit_condor "$round_file" "$firstChunk"
          done
          done
          done
        done
        done
        done
        ;;
      all|"")
        # Clean stale YAML overrides before fresh submission
        rm -f "${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_"*.yaml 2>/dev/null || true

        # Freeze pipeline for this bulk submission
        cleanup_bulk_snapshots_for_tag
        case "$DATASET" in
          isAuAu|isOO) create_pipeline_snapshot "auau" "$(date +%Y%m%d_%H%M%S)" ;;
          *)           create_pipeline_snapshot "pp" "$(date +%Y%m%d_%H%M%S)" ;;
        esac

        n_matrix=$(( ${#data_pts[@]} * ${#data_fracs[@]} * ${#data_vzs[@]} * ${#data_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
        say "═══════════════════════════════════════════════════════════════"
        say "${BOLD}CONDOR ALL: ${n_matrix} matrix configuration(s), groupSize=${GROUP_SIZE}${RST}"
        say "  dataset       : ${DATASET}"
        say "  YAML source   : ${data_yaml_src}"
        say "  golden list   : $(basename "$GOLDEN") ($(grep -cE '^[0-9]+' "$GOLDEN") runs)"
        say "═══════════════════════════════════════════════════════════════"
        echo

        cell_num=0
        total_queued_all=0
        t0_all="$(date +%s)"

        for data_pt in "${data_pts[@]}"; do
        for data_frac in "${data_fracs[@]}"; do
        for data_vz in "${data_vzs[@]}"; do
          for data_cone in "${data_cones[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
          for uepipe in "${uepipe_modes[@]}"; do
          dpt_tag="jetMinPt$(sim_pt_tag "$data_pt")"
          dfrac_tag="$(sim_b2b_tag "$data_frac")"
          dvz_tag="$(sim_vz_tag "$data_vz")"
          dcone_tag="$(sim_cone_tag "$data_cone")"
          data_cfg_tag="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_tags[$iso_idx]}"
          (( uepipe_in_tag )) && data_cfg_tag="${data_cfg_tag}_${uepipe}"
          yaml_override="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${data_cfg_tag}.yaml"
          mkdir -p "$SIM_YAML_OVERRIDE_DIR"
          sed -E \
            -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${data_pt}|" \
            -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${data_frac}|" \
            -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${data_vz}|" \
            -e "s|^([[:space:]]*coneR:).*|\\1 ${data_cone}|" \
            -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${iso_sliding[$iso_idx]}|" \
            -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${iso_fixed[$iso_idx]}|" \
            -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${uepipe}|" \
            "$data_yaml_src" > "$yaml_override"
          export RJ_CONFIG_YAML="$yaml_override"
          DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"

          (( ++cell_num ))
          say "───────────────────────────────────────────────────────────────"
          say "${BOLD}[${cell_num}/${n_matrix}] ${data_cfg_tag}${RST}"
          say "  YAML override : ${yaml_override}"
          say "  DEST_BASE     : ${DEST_BASE}"
          tmp_src="${ROUND_DIR}/ALL_${TAG}_${data_cfg_tag}_$(date +%s).txt"
          grep -E '^[0-9]+' "$GOLDEN" > "$tmp_src"

          submit_condor "$tmp_src" ""
          done
          done
          done
        done
        done
        done

        elapsed_all=$(( $(date +%s) - t0_all ))
        say "═══════════════════════════════════════════════════════════════"
        say "${BOLD}CONDOR ALL complete: ${cell_num} configurations submitted (${elapsed_all}s)${RST}"
        say "═══════════════════════════════════════════════════════════════"
        ;;
      *)
        usage
        ;;
    esac
    ;;

  *)
    usage
    ;;
esac

say "Done."
