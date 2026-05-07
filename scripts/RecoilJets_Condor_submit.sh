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
#     ./RecoilJets_Condor_submit.sh isSim checkModels
#     ./RecoilJets_Condor_submit.sh isSim condorDoAllFromScratch
#     ./RecoilJets_Condor_submit.sh isSim condorDoAll
#     ./RecoilJets_Condor_submit.sh isSim condorHistFromPool
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
#     ./RecoilJets_Condor_submit.sh isSimEmbedded condorHistFromPool
#
#   AuAu embedded photon-ID BDT training extraction:
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainTightBDT localTest
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainTightBDT smokeTestFirstPass
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainTightBDT smokeTestSecondPass SOURCE=/path/to/extraction
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainTightBDT smokeTestApplyExisting MODEL_DIR=/path/to/tight/models
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainNPB local
#     RJ_AUAU_BDT_NPB_DATA_LOCAL=1 ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainNPB local 5000 VERBOSE=10
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainJetMLResidual local
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainMLAll local 5000 VERBOSE=10
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainMLAll local 5000 NFILES=5 VERBOSE=10
#     RJ_ML_PYTHON=/path/to/python ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainMLAll resume /path/to/mlIntegration_run
#
# OVERVIEW
#   • This script is the top-level submission driver for:
#       isPP, isPPrun25, isAuAu, isOO, isSim, isSimJet5, isSimMB, isSimEmbedded
#   • DATA modes use per-run runlists and can do:
#       local, isLocalIsoPing, CHECKJOBS, splitGoldenRunList,
#       condor testJob, condor round, condor allFromScratch, condor all
#   • SIM modes use staged matched master lists and can do:
#       local, checkModels, CHECKJOBS, condorTest,
#       condorDoAllFromScratch, condorDoAll, condorDoAllDirect
#   • Jobs never mix files across runs in DATA mode.
#   • Direct DST/local sweeps use the YAML matrix:
#       jet_pt_min × back_to_back_dphi_min_pi_fraction × vz_cut_cm × coneR × iso mode
#       and, for AuAu-like tags, also clusterUEpipeline.
#   • Pool workflows split the expensive DST pass from replay-only variations:
#       condorDoAllFromScratch / condor allFromScratch capture reusable pools,
#       while condorDoAll / condor all replay existing pools into the same
#       cfg-tagged output tree expected by the merge/get scripts.
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
#     ./RecoilJets_Condor_submit.sh isSim condorHistFromPool
#
#     ./RecoilJets_Condor_submit.sh isSimJet5 local
#     ./RecoilJets_Condor_submit.sh isSimJet5 condorDoAll
#
#     ./RecoilJets_Condor_submit.sh isSimMB local
#     ./RecoilJets_Condor_submit.sh isSimMB condorDoAll
#
#     ./RecoilJets_Condor_submit.sh isSimEmbedded local
#     ./RecoilJets_Condor_submit.sh isSimEmbedded condorDoAll
#     ./RecoilJets_Condor_submit.sh isSimEmbedded condorHistFromPool
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
#     $ ./RecoilJets_Condor_submit.sh <simDataset> condorHistFromPool [groupSize N] [SAMPLE=...]
#     $ ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainTightBDT local [Nevents] [VERBOSE=N]
#     $ ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainTightBDT condorDoAll [groupSize N]
#     $ ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainNPB local [Nevents] [VERBOSE=N]
#     $ ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainJetMLResidual local [Nevents] [VERBOSE=N]
#     $ ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainMLAll local [NeventsPerSample] [NFILES=N] [VERBOSE=N]
#     $ RJ_ML_PYTHON=/path/to/python ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainMLAll resume /path/to/mlIntegration_run
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
AA_GOLDEN="${BASE}/GRLs_tanner/run3auau_pro001_pcdb001_v001_dst_calofitting_grl.list"
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

SCALED_TRIGGER_RUNLIST="${BASE}/dst_lists_auau/scaledEffRuns_MBD_NS_geq_2_vtx_lt_150__Pho10_12.list"
SCALED_TRIGGER_CONFIG_SRC="${BASE}/dst_lists_auau/trigger_scaled_efficiency_studies_auau.txt"
SCALED_TRIGGER_CONFIG_EXPECTED="${BASE}/dst_lists_auau/scaledEffConfig_MBD_NS_geq_2_vtx_lt_150__Pho10_12.txt"
SCALED_TRIGGER_OUTPUT_SUFFIX="${RJ_SCALED_TRIGGER_OUTPUT_SUFFIX:-_scaledTriggerStudy}"

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
AUAU_BDT_DEST_BASE="${RJ_AUAU_BDT_DEST_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdt_training}"
POOL_DEST_ROOT="${RJ_POOL_DEST_ROOT:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaPools}"
AUAU_BDT_LOCAL_BASE="${RJ_AUAU_BDT_LOCAL_BASE:-${BASE}/local_bdt_training_outputs}"
AUAU_BDT_MODEL_BASE="${RJ_AUAU_BDT_MODEL_BASE:-${BASE}/bdt_models}"
AUAU_BDT_SIGNAL_SAMPLES_DEFAULT="run28_embeddedPhoton12 run28_embeddedPhoton20"
AUAU_BDT_BACKGROUND_SAMPLES_DEFAULT="run28_embeddedJet12 run28_embeddedJet20"

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
  say "Keeping existing bulk snapshot dirs for tag ${TAG}; live Condor jobs may still reference them."
  say "Manual cleanup later: remove ${SNAPSHOT_ROOT}/${TAG}_* only after all matching jobs have left the queue."
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
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor smokeTest [groupSize N]${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor poolSmoke [groupSize N] [maxJobs 12]${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor allFromScratch [groupSize N]${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor all [groupSize N]${RST}
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor allDirect [groupSize N]${RST}

${BOLD}SIM mode (photonJet10/20 merged):${RST}
  ${BOLD}$0 isSim local [Nevents] [VERBOSE=N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim checkModels [VERBOSE=N]${RST}
  ${BOLD}$0 isSim CHECKJOBS [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorTest [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorDoAllSmoke [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorDoAllFromScratch [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorDoAll [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorHistFromPool [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorDoAllDirect [groupSize N] [SAMPLE=run28_photonjet10]${RST}

${BOLD}SIM mode (embedded photon AuAu-like signal):${RST}
  ${BOLD}$0 isSimEmbedded local [Nevents] [VERBOSE=N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded CHECKJOBS [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded condorDoAllSmoke [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded condorDoAllFromScratch [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded condorDoAll [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded condorHistFromPool [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded condorDoAllDirect [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}

${BOLD}SIM mode (embedded inclusive-jet AuAu-like background):${RST}
  ${BOLD}$0 isSimEmbeddedInclusive local [Nevents] [VERBOSE=N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive CHECKJOBS [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive condorDoAllSmoke [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive condorDoAllFromScratch [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive condorDoAll [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive condorHistFromPool [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive condorDoAllDirect [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}

${BOLD}SIM mode (photonJet5 single-slice):${RST}
  ${BOLD}$0 isSimJet5 local [Nevents] [VERBOSE=N] [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimJet5 CHECKJOBS [groupSize N] [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimJet5 condorTest [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimJet5 condorDoAll [groupSize N] [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimJet5 condorHistFromPool [groupSize N] [SAMPLE=run28_jet5]${RST}

${BOLD}SIM mode (MinBias DETROIT):${RST}
  ${BOLD}$0 isSimMB local [Nevents] [VERBOSE=N] [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB CHECKJOBS [groupSize N] [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB condorTest [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB condorDoAll [groupSize N] [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB condorHistFromPool [groupSize N] [SAMPLE=run28_detroit]${RST}

${BOLD}AuAu embedded photon-ID BDT training:${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT localTest [Nevents] [VERBOSE=N]${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT smokeTestFirstPass [groupSize N]${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT smokeTestSecondPass SOURCE=/path/to/finished/extraction${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT smokeTestApplyExisting MODEL_DIR=/path/to/tight/models [Nevents] [NFILES=N] [VERBOSE=N]${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainNPB local [Nevents] [VERBOSE=N]${RST}
  ${BOLD}RJ_AUAU_BDT_NPB_DATA_LOCAL=1 $0 isSimEmbeddedAndInclusive trainNPB local [Nevents] [VERBOSE=N]${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainJetMLResidual local [Nevents] [VERBOSE=N]${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainMLAll local [NeventsPerSample] [NFILES=N] [VERBOSE=N]${RST}
  ${BOLD}RJ_ML_PYTHON=/path/to/python $0 isSimEmbeddedAndInclusive trainMLAll resume /path/to/mlIntegration_run${RST}

  Signal samples default to:     ${AUAU_BDT_SIGNAL_SAMPLES_DEFAULT}
  Background samples default to: ${AUAU_BDT_BACKGROUND_SAMPLES_DEFAULT}
  Override with RJ_AUAU_BDT_SIGNAL_SAMPLES / RJ_AUAU_BDT_BACKGROUND_SAMPLES.
  ML training uses RJ_ML_PYTHON if set; otherwise it uses python3.
  Local ML event targets assume ${RJ_SIM_LOCAL_EVENTS_PER_FILE:-1000} events/input file and group enough input lines
  into one local Fun4All invocation per sample; override with NFILES=N or RJ_SIM_LOCAL_EVENTS_PER_FILE.

Examples:
  $0 isPP  local 5000
  $0 isAuAu splitGoldenRunList groupSize 3 maxJobs 12000
  $0 isAuAu condor round 1
  $0 isPP  condor smokeTest
  $0 isPP  condor poolSmoke groupSize 7 maxJobs 12
  $0 isPP  condor all groupSize 4
  $0 isPP  CHECKJOBS groupSize 4

  $0 isSim CHECKJOBS groupSize 5
  $0 isSim local 5000
  $0 isSim checkModels
  $0 isSim condorTest
  $0 isSim condorDoAllSmoke
  $0 isSimEmbedded condorDoAllSmoke
  $0 isSimEmbeddedInclusive condorDoAllSmoke
  $0 isSim condorDoAllFromScratch groupSize 5
  $0 isSim condorDoAll groupSize 5
  $0 isSim condorHistFromPool groupSize 5

  $0 isSimJet5 local 5000
  $0 isSimJet5 condorDoAll groupSize 5
  $0 isSimJet5 condorHistFromPool groupSize 5
  $0 isSimMB local 5000
  $0 isSimMB condorDoAll groupSize 5
  $0 isSimMB condorHistFromPool groupSize 5
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

memory_request_to_mb() {
  local request="$1"
  local digits="${request//[^0-9]/}"
  local mb="${digits:-0}"
  case "$request" in
    *[Gg][Bb]|*[Gg]) mb=$(( mb * 1024 )) ;;
  esac
  printf '%s\n' "$mb"
}

truthy_env() {
  case "${1:-}" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

condor_pool_submit_extra_lines() {
  local stage="$1"
  local is_smoke="${2:-0}"
  local use_late="${RJ_CONDOR_USE_LATE_MATERIALIZATION:-1}"
  local max_idle="" max_materialize="" max_runtime=""

  if truthy_env "$use_late"; then
    case "$stage" in
      capture)
        max_idle="${RJ_CONDOR_CAPTURE_MAX_IDLE:-${RJ_CONDOR_MAX_IDLE:-2000}}"
        max_materialize="${RJ_CONDOR_CAPTURE_MAX_MATERIALIZE:-${RJ_CONDOR_MAX_MATERIALIZE:-5000}}"
        ;;
      replay)
        max_idle="${RJ_CONDOR_REPLAY_MAX_IDLE:-${RJ_CONDOR_MAX_IDLE:-1000}}"
        max_materialize="${RJ_CONDOR_REPLAY_MAX_MATERIALIZE:-${RJ_CONDOR_MAX_MATERIALIZE:-3000}}"
        ;;
      *)
        max_idle="${RJ_CONDOR_MAX_IDLE:-1000}"
        max_materialize="${RJ_CONDOR_MAX_MATERIALIZE:-3000}"
        ;;
    esac
    if [[ "$max_idle" =~ ^[0-9]+$ && "$max_idle" -gt 0 ]]; then
      printf 'max_idle     = %s\n' "$max_idle"
    fi
    if [[ "$max_materialize" =~ ^[0-9]+$ && "$max_materialize" -gt 0 ]]; then
      printf 'max_materialize = %s\n' "$max_materialize"
    fi
  fi

  if [[ "$is_smoke" == "1" ]]; then
    case "$stage" in
      capture) max_runtime="${RJ_SMOKE_CAPTURE_MAX_RUNTIME_SECONDS:-${RJ_SMOKE_MAX_RUNTIME_SECONDS:-1200}}" ;;
      replay)  max_runtime="${RJ_SMOKE_REPLAY_MAX_RUNTIME_SECONDS:-${RJ_SMOKE_MAX_RUNTIME_SECONDS:-1200}}" ;;
      *)       max_runtime="${RJ_SMOKE_MAX_RUNTIME_SECONDS:-1200}" ;;
    esac
    if [[ "$max_runtime" =~ ^[0-9]+$ && "$max_runtime" -gt 0 ]]; then
      printf 'periodic_remove = (JobStatus == 2) && ((time() - EnteredCurrentStatus) > %s)\n' "$max_runtime"
    fi
  fi
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

notify_emails_csv() {
  if [[ -n "${RJ_NOTIFY_EMAILS:-}" ]]; then
    printf "%s\n" "$RJ_NOTIFY_EMAILS"
    return 0
  fi
  local yaml
  yaml="$(sim_yaml_master_path)"
  [[ -f "$yaml" ]] || return 0
  local -a emails=()
  mapfile -t emails < <( yaml_get_values "notify_emails" "$yaml" 2>/dev/null || true )
  local out=""
  local e
  for e in "${emails[@]}"; do
    e="$(trim_ws "$e")"
    [[ -z "$e" ]] && continue
    out="${out:+${out},}${e}"
  done
  printf "%s\n" "$out"
}

submit_dag_with_notify() {
  local dag="$1"
  if [[ "${RJ_DAG_DRYRUN:-0}" == "1" || "${RJ_DAG_DRYRUN:-0}" == "true" || "${RJ_DAG_DRYRUN:-0}" == "TRUE" ]]; then
    say "RJ_DAG_DRYRUN=1 → DAG files were built but not submitted."
    say "  dag      : ${dag}"
    say "  manifest : $(dirname "$dag")/manifest.txt"
    echo "RECOILJETS_DAG_DRYRUN_V1 dag=${dag} manifest=$(dirname "$dag")/manifest.txt"
    return 0
  fi
  if [[ "${RJ_DAG_DRYRUN:-0}" != "1" && "${RJ_DAG_DRYRUN:-0}" != "true" && "${RJ_DAG_DRYRUN:-0}" != "TRUE" ]]; then
    need_cmd condor_submit_dag
  fi
  local emails
  emails="$(notify_emails_csv)"
  if [[ -n "$emails" ]]; then
    if [[ "${RJ_DAGMAN_NATIVE_EMAIL:-0}" == "1" || "${RJ_DAGMAN_NATIVE_EMAIL:-0}" == "true" || "${RJ_DAGMAN_NATIVE_EMAIL:-0}" == "TRUE" ]]; then
      say "Submitting DAG with native DAGMan error email plus parseable FINAL summary → ${emails}"
      condor_submit_dag -notification Error -append "notify_user = ${emails}" "$dag"
    else
      say "Submitting DAG with quiet DAGMan native mail; parseable FINAL summary will notify → ${emails}"
      condor_submit_dag -notification Never "$dag"
    fi
  else
    say "Submitting DAG without email notification; notify_emails/RJ_NOTIFY_EMAILS is empty."
    condor_submit_dag "$dag"
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

matrix_cfg_tag() {
  local pt="$1" frac="$2" vz="$3" cone="$4" iso_idx="$5" uepipe="$6"
  local tag
  tag="jetMinPt$(sim_pt_tag "$pt")_$(sim_b2b_tag "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso_base_tags[$iso_idx]}"
  (( ${uepipe_in_tag:-0} )) && tag="${tag}_${uepipe}"
  tag="${tag}_${iso_selection_tags[$iso_idx]}"
  echo "$tag"
}

working_point_cfg_tag() {
  local iso_idx="$1" uepipe="$2"
  local tag="wp_${iso_selection_tags[$iso_idx]}"
  (( ${uepipe_in_tag:-0} )) && tag="${uepipe}_${tag}"
  echo "$tag"
}

analysis_view_id() {
  local pt="$1" frac="$2" vz="$3" cone="$4" iso_idx="$5"
  printf 'jetPt%s_%s_%s_%s_%s\n' \
    "$(sim_pt_tag "$pt")" \
    "$(sim_b2b_tag "$frac")" \
    "$(sim_vz_tag "$vz")" \
    "${iso_base_tags[$iso_idx]}" \
    "$(sim_cone_tag "$cone")"
}

read_replay_cones() {
  local yaml="$1"
  replay_cones=()
  mapfile -t replay_cones < <( yaml_get_values "coneR" "$yaml" 2>/dev/null )
  (( ${#replay_cones[@]} )) || replay_cones=( "0.4" )
}

read_capture_cones() {
  local yaml="$1"
  capture_cones=()
  local -a stored_cones_tmp
  mapfile -t stored_cones_tmp < <( yaml_get_values "stored_isolation_cones" "$yaml" 2>/dev/null )
  if (( ${#stored_cones_tmp[@]} > 0 )); then
    capture_cones=( "${stored_cones_tmp[0]}" )
  fi
  if (( ${#capture_cones[@]} == 0 )); then
    mapfile -t capture_cones < <( yaml_get_values "captured_coneR" "$yaml" 2>/dev/null )
  fi
  if (( ${#capture_cones[@]} == 0 )); then
    local -a replay_tmp
    mapfile -t replay_tmp < <( yaml_get_values "coneR" "$yaml" 2>/dev/null )
    (( ${#replay_tmp[@]} )) && capture_cones=( "${replay_tmp[0]}" )
  fi
  (( ${#capture_cones[@]} )) || capture_cones=( "0.4" )
}

pool_schema_tag() {
  local yaml="$1"
  local line rhs
  line="$(grep -E "^[[:space:]]*(pool_schema_version|schema_version):" "$yaml" | head -n 1 || true)"
  rhs="${line#*:}"
  rhs="$(trim_ws "${rhs%%#*}")"
  [[ -n "$rhs" ]] || rhs="2"
  rhs="${rhs//./p}"
  printf "schema%s\n" "$rhs"
}

pool_capture_cfg_tag() {
  local cone="$1" uepipe="$2" yaml="$3"
  local tag
  local -a stored
  mapfile -t stored < <( yaml_get_values "stored_isolation_cones" "$yaml" 2>/dev/null )
  if (( ${#stored[@]} > 1 )); then
    tag="pool_$(pool_schema_tag "$yaml")_isoRmulti"
  else
    tag="pool_$(pool_schema_tag "$yaml")_$(sim_cone_tag "$cone")"
  fi
  (( ${uepipe_in_tag:-0} )) && tag="${tag}_${uepipe}"
  echo "$tag"
}

sim_analysis_tag_for_dataset() {
  case "$1" in
    isSimJet5) echo "isSimJet5" ;;
    isSimMB) echo "isSimMB" ;;
    isSimEmbedded) echo "isSimEmbedded" ;;
    isSimEmbeddedInclusive) echo "isSimEmbeddedInclusive" ;;
    *) echo "isSim" ;;
  esac
}

sanitize_node_name() {
  local s="$1"
  s="${s//[^A-Za-z0-9_]/_}"
  [[ "$s" =~ ^[A-Za-z] ]] || s="N_${s}"
  printf "%s\n" "$s"
}

iso_group_key() {
  local iso_idx="$1"
  printf '%s|%s|%s\n' "${iso_base_tags[$iso_idx]}" "${iso_sliding[$iso_idx]}" "${iso_fixed[$iso_idx]}"
}

iso_idx_is_group_leader() {
  local iso_idx="$1"
  local key
  key="$(iso_group_key "$iso_idx")"
  local j
  for (( j=0; j<iso_idx; j++ )); do
    [[ "$(iso_group_key "$j")" == "$key" ]] && return 1
  done
  return 0
}

iso_group_count() {
  local c=0 i
  for (( i=0; i<${#iso_tags[@]}; i++ )); do
    if iso_idx_is_group_leader "$i"; then (( c+=1 )); fi
  done
  printf '%d\n' "$c"
}

emit_id_fanout_dirs_file() {
  local out_file="$1" dest_root="$2" pt="$3" frac="$4" vz="$5" cone="$6" iso_idx="$7" uepipe="$8"
  local key cfg j
  key="$(iso_group_key "$iso_idx")"
  : > "$out_file"
  for (( j=0; j<${#iso_tags[@]}; j++ )); do
    [[ "$(iso_group_key "$j")" == "$key" ]] || continue
    cfg="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$j" "$uepipe")"
    printf '%s|%s|%s|%s|%s\n' \
      "${dest_root}/${cfg}" "$cfg" "${iso_preselection[$j]}" "${iso_tight[$j]}" "${iso_nonTight[$j]}" >> "$out_file"
  done
}

emit_pool_replay_fanout_dirs_file() {
  local out_file="$1" dest_root="$2" master="$3" stamp="$4" cones_text="$5" uepipe="$6" pts_text="$7" fracs_text="$8" vzs_text="$9"
  local pt frac vz cone iso_idx cfg wp_cfg view_id inline_spec out_cfg materialize
  local public_layout="${RJ_POOL_REPLAY_PUBLIC_LAYOUT:-legacy}"
  if [[ "${RJ_POOL_REPLAY_USE_INTERNAL_VIEWS:-0}" == "1" || "${RJ_POOL_REPLAY_USE_INTERNAL_VIEWS:-0}" == "true" || "${RJ_POOL_REPLAY_USE_INTERNAL_VIEWS:-0}" == "TRUE" ]]; then
    public_layout="views"
  fi
  local -a cones_arr pts_arr fracs_arr vzs_arr
  read -r -a cones_arr <<< "$cones_text"
  read -r -a pts_arr <<< "$pts_text"
  read -r -a fracs_arr <<< "$fracs_text"
  read -r -a vzs_arr <<< "$vzs_text"
  local total=$(( ${#cones_arr[@]} * ${#pts_arr[@]} * ${#fracs_arr[@]} * ${#vzs_arr[@]} * ${#iso_tags[@]} ))
  local count=0
  local emit_trace=0
  if [[ "${trace_pool:-0}" == "1" || "${RJ_POOL_TRACE:-0}" == "1" || "${RJ_POOL_TRACE:-0}" == "true" || "${RJ_POOL_TRACE:-0}" == "TRUE" ]]; then
    emit_trace=1
  fi
  local emit_every="${trace_every:-${RJ_POOL_TRACE_EVERY:-250}}"
  [[ "$emit_every" =~ ^[0-9]+$ && "$emit_every" -gt 0 ]] || emit_every=250
  (( emit_trace )) && say "[poolFanout] writing ${total} inline replay view rows → ${out_file}"
  : > "$out_file"
  for cone in "${cones_arr[@]}"; do
    for pt in "${pts_arr[@]}"; do
      for frac in "${fracs_arr[@]}"; do
        for vz in "${vzs_arr[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
            cfg="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$iso_idx" "$uepipe")"
            wp_cfg="$(working_point_cfg_tag "$iso_idx" "$uepipe")"
            if [[ "$public_layout" == "views" ]]; then
              out_cfg="$wp_cfg"
              view_id="$(analysis_view_id "$pt" "$frac" "$vz" "$cone" "$iso_idx")"
              materialize="physics"
            else
              out_cfg="$cfg"
              view_id="legacy"
              materialize="legacy"
            fi
            inline_spec="INLINE_VIEW_V1;jet_pt_min=${pt};back_to_back_dphi_min_pi_fraction=${frac};vz_cut_cm=${vz};coneR=${cone};isSlidingIso=${iso_sliding[$iso_idx]};fixedGeV=${iso_fixed[$iso_idx]};clusterUEpipeline=${uepipe};cfg_tag=${cfg};view_id=${view_id}"
            printf '%s|%s|%s|%s|%s|%s|%s|%s\n' \
              "${dest_root}/${out_cfg}" "$cfg" "$inline_spec" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}" "$view_id" "$materialize" >> "$out_file"
            (( count+=1 ))
            (( emit_trace && (count == 1 || count % emit_every == 0 || count == total) )) && say "[poolFanout] wrote ${count}/${total} rows (${cfg})"
          done
        done
      done
    done
  done
  return 0
}

fanout_dest_allowed() {
  local fan_dest="$1"
  case "$fan_dest" in
    */thesisAna/pp/*|*/thesisAna/pp25/*|*/thesisAna/auau/*|*/thesisAna/oo/*|*/thesisAna/sim/*|*/thesisAna/simembedded/*|*/thesisAna/simembeddedinclusive/*|*/thesisAna/simjet5/*|*/thesisAna/simmb/*)
      return 0
      ;;
    */thesisAnaSmoke/pp_smokeTest_*/*|*/thesisAnaSmoke/pp25_smokeTest_*/*|*/thesisAnaSmoke/auau_smokeTest_*/*|*/thesisAnaSmoke/oo_smokeTest_*/*|*/thesisAnaSmoke/sim_smokeTest_*/*|*/thesisAnaSmoke/simembedded_smokeTest_*/*|*/thesisAnaSmoke/simembeddedinclusive_smokeTest_*/*|*/thesisAnaSmoke/simjet5_smokeTest_*/*|*/thesisAnaSmoke/simmb_smokeTest_*/*)
      return 0
      ;;
  esac
  return 1
}

clean_fanout_output_dirs_from_file() {
  local fanout_file="$1" label="$2" trace="${3:-0}" every="${4:-25}"
  shift 4 || true
  local -a leaf_dirs=( "$@" )
  declare -A cleaned_dest=()
  local fan_line fan_dest
  local rows=0 cleaned=0
  while IFS= read -r fan_line; do
    [[ -z "${fan_line:-}" || "${fan_line:0:1}" == "#" ]] && continue
    IFS='|' read -r fan_dest _fan_cfg _fan_yaml _fan_pre _fan_tight _fan_nonTight <<< "$fan_line"
    [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
    (( rows+=1 ))
    [[ -n "${cleaned_dest[$fan_dest]:-}" ]] && continue
    cleaned_dest["$fan_dest"]=1
    fanout_dest_allowed "$fan_dest" || { err "Refusing to wipe fanout DEST_BASE='$fan_dest'"; exit 62; }
    (( cleaned+=1 ))
    (( trace && (cleaned == 1 || cleaned % every == 0) )) && say "[poolFanout] cleaning unique output ${cleaned} (${label}): ${fan_dest}"
    if (( ${#leaf_dirs[@]} > 0 )); then
      local leaf
      for leaf in "${leaf_dirs[@]}"; do
        mkdir -p "${fan_dest}/${leaf}"
        rm -f "${fan_dest}/${leaf}/"*.root 2>/dev/null || true
      done
    else
      mkdir -p "$fan_dest"
      find "$fan_dest" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true
    fi
  done < "$fanout_file"
  (( trace )) && say "[poolFanout] cleaned ${cleaned} unique output roots from ${rows} fanout rows (${label})"
  CLEAN_FANOUT_LAST_COUNT="$cleaned"
}

append_fanout_cfgs_to_manifest() {
  local fanout_file="$1" manifest="$2"
  [[ -s "$fanout_file" ]] || return 0
  awk -F'|' '
    NF >= 2 && $1 !~ /^#/ && $1 != "" {
      out=$1
      sub(/\/$/, "", out)
      n=split(out, parts, "/")
      if (parts[n] != "") print "cfg_tag=" parts[n]
      if ($2 != "") print "view_cfg_tag=" $2
      if (NF >= 7 && $7 != "") print "view_id=" $7
    }' "$fanout_file" >> "$manifest"
}

split_pool_replay_fanout_shards() {
  local fanout_file="$1" shard_dir="$2" label="$3" roots_per_shard="${4:-2}" trace="${5:-0}"
  [[ -s "$fanout_file" ]] || { err "Cannot shard empty fanout file: ${fanout_file}"; exit 95; }
  [[ "$roots_per_shard" =~ ^[0-9]+$ && "$roots_per_shard" -gt 0 ]] || roots_per_shard=2
  mkdir -p "$shard_dir"

  local safe_label shard_idx roots_in_shard rows unique_roots line dest shard_file i
  safe_label="$(sanitize_node_name "$label")"
  rm -f "${shard_dir}/${safe_label}_shard"*.txt 2>/dev/null || true

  declare -A root_to_shard=()
  shard_idx=0
  roots_in_shard=0
  rows=0
  unique_roots=0

  while IFS= read -r line; do
    [[ -z "${line:-}" || "${line:0:1}" == "#" ]] && continue
    dest="${line%%|*}"
    [[ -n "$dest" ]] || continue
    if [[ -z "${root_to_shard[$dest]:-}" ]]; then
      if (( shard_idx == 0 || roots_in_shard >= roots_per_shard )); then
        (( shard_idx+=1 ))
        roots_in_shard=0
        shard_file="${shard_dir}/${safe_label}_shard$(printf "%03d" "$shard_idx").txt"
        : > "$shard_file"
      fi
      root_to_shard["$dest"]="$shard_idx"
      (( roots_in_shard+=1 ))
      (( unique_roots+=1 ))
    fi
    shard_file="${shard_dir}/${safe_label}_shard$(printf "%03d" "${root_to_shard[$dest]}").txt"
    printf '%s\n' "$line" >> "$shard_file"
    (( rows+=1 ))
  done < "$fanout_file"

  (( shard_idx > 0 )) || { err "Fanout sharding produced no shards for ${fanout_file}"; exit 95; }
  (( trace )) && say "[poolFanout] sharded ${rows} rows / ${unique_roots} output roots into ${shard_idx} replay shard(s), rootsPerShard=${roots_per_shard} (${label})" >&2
  FANOUT_SHARD_LAST_COUNT="$shard_idx"
  FANOUT_SHARD_LAST_ROOTS="$unique_roots"
  FANOUT_SHARD_LAST_ROWS="$rows"
  for (( i=1; i<=shard_idx; i++ )); do
    shard_file="${shard_dir}/${safe_label}_shard$(printf "%03d" "$i").txt"
    [[ -s "$shard_file" ]] && printf '%s\n' "$shard_file"
  done
}

check_dag_worker_job_budget() {
  local label="$1" capture_jobs="$2" replay_jobs="$3" dryrun="${4:-0}" min_jobs="${5:-0}"
  local max_jobs="${RJ_DAG_MAX_WORKER_JOBS:-50000}"
  local warn_jobs="${RJ_DAG_WARN_WORKER_JOBS:-40000}"
  [[ "$max_jobs" =~ ^[0-9]+$ && "$max_jobs" -gt 0 ]] || max_jobs=50000
  [[ "$warn_jobs" =~ ^[0-9]+$ && "$warn_jobs" -gt 0 ]] || warn_jobs=40000
  [[ "$min_jobs" =~ ^[0-9]+$ ]] || min_jobs=0
  local total_jobs=$(( capture_jobs + replay_jobs ))

  say "${BOLD}Worker-job budget check:${RST} ${label}"
  say "  capture worker jobs : ${capture_jobs}"
  say "  replay worker jobs  : ${replay_jobs}"
  say "  total worker jobs   : ${total_jobs}"
  (( min_jobs > 0 )) && say "  target minimum jobs : ${min_jobs}"
  say "  warning / hard cap  : ${warn_jobs} / ${max_jobs}"

  if (( min_jobs > 0 && total_jobs < min_jobs )); then
    warn "${label} is below the preferred worker-job target (${total_jobs} < ${min_jobs}). This is allowed, but jobs may be too coarse to justify lower memory requests; consider lowering groupSize/RJ_POOL_REPLAY_GROUP_SIZE or RJ_REPLAY_OUTPUT_ROOTS_PER_SHARD if runtime or slot matching looks poor."
  fi
  if (( total_jobs > max_jobs )); then
    if (( dryrun )); then
      warn "Dry-run DAG exceeds RJ_DAG_MAX_WORKER_JOBS=${max_jobs}; not submitting because RJ_DAG_DRYRUN=1."
    elif [[ "${RJ_DAG_ALLOW_OVER_BUDGET:-0}" == "1" || "${RJ_DAG_ALLOW_OVER_BUDGET:-0}" == "true" || "${RJ_DAG_ALLOW_OVER_BUDGET:-0}" == "TRUE" ]]; then
      warn "DAG exceeds RJ_DAG_MAX_WORKER_JOBS=${max_jobs}, but RJ_DAG_ALLOW_OVER_BUDGET is enabled."
    else
      err "Refusing to submit ${label}: planned worker jobs (${total_jobs}) exceed RJ_DAG_MAX_WORKER_JOBS=${max_jobs}."
      err "Reduce job count by increasing groupSize / RJ_POOL_REPLAY_GROUP_SIZE, increasing RJ_REPLAY_OUTPUT_ROOTS_PER_SHARD, or splitting the production into smaller rounds."
      exit 64
    fi
  elif (( total_jobs >= warn_jobs )); then
    warn "${label} is in the high job-count zone (${total_jobs} worker jobs). This is allowed, but keep it near the 40k-50k comfort band."
  fi
}

publish_production_manifest() {
  local manifest="$1" out_base="$2"
  [[ -s "$manifest" && -n "$out_base" ]] || return 0
  mkdir -p "$out_base"
  cp "$manifest" "${out_base}/.recoiljets_latest_manifest.txt"
}

write_dag_final_summary_files() {
  local dag_dir="$1" workflow="$2" dataset="$3" pool_base="$4" out_base="$5" manifest="$6" dag_file="${7:-}"
  local emails
  emails="$(notify_emails_csv)"
  local exec_file="${dag_dir}/final_summary_${workflow}.sh"
  local sub_file="${dag_dir}/final_summary_${workflow}.sub"
  mkdir -p "$dag_dir"

  cat > "$exec_file" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
workflow="${1:?workflow required}"
dataset="${2:?dataset required}"
pool_base="${3:?pool base required}"
out_base="${4:?output base required}"
manifest="${5:?manifest required}"
emails="${6:-}"
dag_file="${7:-}"
[[ "$emails" == "NONE" ]] && emails=""
dag_dir="$(dirname "$manifest")"

pool_count=0
out_count=0
[[ -d "$pool_base" ]] && pool_count=$(find "$pool_base" -type f -name '*.root' 2>/dev/null | wc -l | awk '{print $1}')
[[ -d "$out_base" ]] && out_count=$(find "$out_base" -type f -name '*.root' 2>/dev/null | wc -l | awk '{print $1}')

status="READY"
status_note="DAG reached its FINAL summary node."
dagman_out=""
nodes_log=""
rescue_count=0
if [[ -n "$dag_file" && "$dag_file" != "NONE" ]]; then
  dagman_out="${dag_file}.dagman.out"
  nodes_log="${dag_file}.nodes.log"
  if compgen -G "${dag_file}.rescue*" >/dev/null 2>&1; then
    rescue_count=$(compgen -G "${dag_file}.rescue*" | wc -l | awk '{print $1}')
    status="FAILED"
    status_note="DAGMan rescue file(s) were produced; inspect the DAG logs before continuing."
  elif [[ -s "$dagman_out" ]] && grep -Eiq 'DAG_STATUS_NODE_FAILED|ERROR: the following Node|failed with status|Node return val:[[:space:]]*[1-9]' "$dagman_out"; then
    status="FAILED"
    status_note="At least one DAG node failed; inspect the failed worker stdout/err before continuing."
  elif [[ -s "$dagman_out" ]] && grep -Eiq 'DAG abort|aborted|Job was held|held job|ULOG_JOB_HELD|DAG_STATUS_RM' "$dagman_out"; then
    status="CHECK"
    status_note="DAGMan log contains error/hold-like text; inspect logs before treating outputs as final."
  fi
fi

body="$(mktemp "${TMPDIR:-/tmp}/recoiljets_dag_summary.XXXXXX")"
profile_lines="$(mktemp "${TMPDIR:-/tmp}/recoiljets_profile_lines.XXXXXX")"
heartbeat_lines="$(mktemp "${TMPDIR:-/tmp}/recoiljets_heartbeat_lines.XXXXXX")"
profile_summary="${dag_dir}/smoke_profile_summary_${workflow}.txt"
profile_raw="${dag_dir}/smoke_profile_rows_${workflow}.txt"
report_dir="${out_base%/}/_pipeline_reports/${workflow}"
report_summary="${report_dir}/final_summary.txt"
report_manifest="${report_dir}/manifest.txt"
report_profile="${report_dir}/smoke_profile_summary.txt"
report_profile_rows="${report_dir}/smoke_profile_rows.txt"
report_job_profile="${report_dir}/job_profile_summary.txt"
report_job_profile_rows="${report_dir}/job_profile_rows.txt"
report_job_heartbeat_rows="${report_dir}/job_heartbeat_rows.txt"
report_merge_profile="${report_dir}/merge_profile_summary.txt"
report_merge_profile_rows="${report_dir}/merge_profile_rows.txt"
report_root_summary="${report_dir}/root_output_summary.txt"
report_hist_summary="${report_dir}/histogram_sanity_summary.txt"
report_failure_summary="${report_dir}/failure_summary.txt"
report_tuning_json="${report_dir}/tuning_inputs.json"
report_smoke_run_stats="${report_dir}/smoke_selected_runs_stats.txt"
manifest_get() {
  local key="$1"
  [[ -s "$manifest" ]] || return 0
  awk -F= -v key="$key" '$1 == key { sub(/^[^=]*=/, ""); print; exit }' "$manifest" 2>/dev/null || true
}
manifest_group_size="$(manifest_get group_size)"
manifest_replay_group_size="$(manifest_get replay_group_size)"
manifest_request_memory="$(manifest_get request_memory)"
manifest_request_memory_mb="$(manifest_get request_memory_mb)"
manifest_capture_request_memory="$(manifest_get capture_request_memory)"
manifest_replay_request_memory="$(manifest_get replay_request_memory)"
manifest_capture_request_memory_mb="$(manifest_get capture_request_memory_mb)"
manifest_replay_request_memory_mb="$(manifest_get replay_request_memory_mb)"
manifest_replay_output_roots_per_shard="$(manifest_get replay_output_roots_per_shard)"
manifest_capture_nevents="$(manifest_get capture_nevents)"
manifest_smoke_cap="$(manifest_get smoke_capture_job_cap)"
manifest_smoke_run_stats="$(manifest_get smoke_selected_run_stats)"
manifest_photon_rows="$(manifest_get photon_id_sets_expanded)"
manifest_replay_cones="$(manifest_get replay_cones)"
manifest_capture_cones="$(manifest_get capture_cones)"
manifest_profile_job="$(manifest_get profile_job)"
expected_capture_jobs=0
expected_replay_jobs=0
expected_total_jobs=0
if [[ -s "$manifest" ]]; then
  expected_capture_jobs="$(awk '
    /^capture_node=/ {
      for (i = 1; i <= NF; ++i) {
        if ($i ~ /^jobs=/) {
          sub(/^jobs=/, "", $i)
          s += $i + 0
        }
      }
    }
    END { print s + 0 }
  ' "$manifest")"
  expected_replay_jobs="$(awk '
    /^replay_node=/ {
      for (i = 1; i <= NF; ++i) {
        if ($i ~ /^jobs=/) {
          sub(/^jobs=/, "", $i)
          s += $i + 0
        }
      }
    }
    END { print s + 0 }
  ' "$manifest")"
fi
expected_total_jobs=$(( expected_capture_jobs + expected_replay_jobs ))

mapfile -t profile_globs_status < <([[ -s "$manifest" ]] && awk -F= '/^profile_glob=/ {print $2}' "$manifest" | sort -u || true)
if (( ${#profile_globs_status[@]} > 0 )); then
  for profile_glob in "${profile_globs_status[@]}"; do
    while IFS= read -r profile_file; do
      [[ -s "$profile_file" ]] || continue
      grep -h 'RECOILJETS_JOB_PROFILE_V1' "$profile_file" >> "$profile_lines" || true
      grep -h 'RECOILJETS_JOB_HEARTBEAT_V1' "$profile_file" >> "$heartbeat_lines" || true
    done < <(compgen -G "$profile_glob" || true)
  done
fi
profile_row_count=0
profile_failure_count=0
if [[ -s "$profile_lines" ]]; then
  profile_row_count="$(grep -c 'RECOILJETS_JOB_PROFILE_V1' "$profile_lines" || echo 0)"
  profile_failure_count="$(awk '
    {
      code = ""
      for (i = 1; i <= NF; ++i) {
        if ($i ~ /^exit_code=/) {
          code = $i
          sub(/^exit_code=/, "", code)
        }
      }
      if (code != "" && code != "0") failed += 1
    }
    END { print failed + 0 }
  ' "$profile_lines")"
fi
if (( expected_capture_jobs > 0 && pool_count < expected_capture_jobs )); then
  status="FAILED"
  status_note="Fewer pool ROOT files were produced than planned capture jobs; inspect worker stdout/err and DAG logs before continuing."
elif (( expected_replay_jobs > 0 && out_count == 0 )); then
  status="FAILED"
  status_note="Replay jobs were planned but no output ROOT files were produced; inspect replay worker stdout/err and DAG logs."
elif [[ "$manifest_profile_job" == "1" ]] && (( expected_total_jobs > 0 && profile_row_count < expected_total_jobs )); then
  status="CHECK"
  status_note="Fewer worker profile rows were found than planned DAG worker jobs; inspect profile globs and worker stdout before tuning."
fi
if (( profile_failure_count > 0 )); then
  status="FAILED"
  status_note="One or more worker profiles reported non-zero exit codes; inspect failed_profile lines and worker logs before continuing."
fi
{
  echo "RECOILJETS_STAGE_EMAIL_V1"
  echo "status=${status}"
  echo "status_note=${status_note}"
  echo "dataset=${dataset}"
  echo "stage=${workflow}"
  echo "stage_type=pool_dag"
  echo "pool_base=${pool_base}"
  echo "out_base=${out_base}"
  echo "manifest=${manifest}"
  echo "pool_root_files=${pool_count}"
  echo "output_root_files=${out_count}"
  echo "expected_capture_jobs=${expected_capture_jobs}"
  echo "expected_replay_jobs=${expected_replay_jobs}"
  echo "profile_rows=${profile_row_count}"
  echo "profile_failures=${profile_failure_count}"
  echo "dag_file=${dag_file}"
  echo "dagman_out=${dagman_out}"
  echo "nodes_log=${nodes_log}"
  echo "rescue_file_count=${rescue_count}"
  echo "profile_summary=${profile_summary}"
  echo "published_report_dir=${report_dir}"
  echo "published_final_summary=${report_summary}"
  echo "published_manifest=${report_manifest}"
  echo "published_profile_summary=${report_profile}"
  echo "published_profile_rows=${report_profile_rows}"
  echo "published_heartbeat_rows=${report_job_heartbeat_rows}"
  [[ -n "$manifest_smoke_run_stats" ]] && echo "published_smoke_selected_run_stats=${report_smoke_run_stats}"
  echo "next_action=If status is READY, continue with the normal merge/pull/plot step for this dataset. If status is CHECK or FAILED, inspect dagman_out, nodes_log, and the job .err files first."
  echo
  if [[ -s "$manifest" ]]; then
    echo "manifest:"
    sed -n '1,160p' "$manifest"
    if (( ${#profile_globs_status[@]} > 0 )); then
      echo
      echo "job_profiles:"
      [[ -s "$profile_lines" ]] && cat "$profile_lines"
      if [[ -s "$profile_lines" ]]; then
        {
          echo "RECOILJETS_SMOKE_PROFILE_SUMMARY_V1"
          echo "workflow=${workflow}"
          echo "dataset=${dataset}"
          echo "manifest=${manifest}"
          awk \
            -v dataset="${dataset}" \
            -v group_size="${manifest_group_size}" \
            -v replay_group_size="${manifest_replay_group_size}" \
            -v request_memory="${manifest_request_memory}" \
            -v request_memory_mb="${manifest_request_memory_mb}" \
            -v capture_request_memory="${manifest_capture_request_memory}" \
            -v replay_request_memory="${manifest_replay_request_memory}" \
            -v capture_request_memory_mb="${manifest_capture_request_memory_mb}" \
            -v replay_request_memory_mb="${manifest_replay_request_memory_mb}" \
            -v replay_output_roots_per_shard="${manifest_replay_output_roots_per_shard}" \
            -v capture_nevents="${manifest_capture_nevents}" \
            -v smoke_cap="${manifest_smoke_cap}" \
            -v photon_rows="${manifest_photon_rows}" \
            -v replay_cones="${manifest_replay_cones}" \
            -v capture_cones="${manifest_capture_cones}" '
            function getv(k,    i, s) {
              for (i = 1; i <= NF; ++i) {
                s = $i
                if (s ~ "^" k "=") {
                  sub("^" k "=", "", s)
                  return s
                }
              }
              return ""
            }
            function mb_from_kb(kb) { return int((kb + 1023) / 1024) }
            function rec_mem_mb(max_mb) {
              if (max_mb <= 0) return 0
              return int(((max_mb * 3 / 2) + 512 + 499) / 500) * 500
            }
            function max(a,b) { return (a > b ? a : b) }
            function int_or_zero(x) { return (x == "" ? 0 : x + 0) }
            {
              stage = getv("stage"); if (stage == "") stage = "unknown"
              exit_code = getv("exit_code")
              elapsed = getv("elapsed_seconds") + 0
              rss = getv("max_rss_kb") + 0
              inputs = getv("input_files") + 0
              outputs = getv("output_files") + 0
              bytes = getv("output_bytes") + 0
              req_mb = getv("request_memory_mb") + 0
              user_cpu = getv("user_cpu_s") + 0
              system_cpu = getv("system_cpu_s") + 0
              cpu_pct = getv("cpu_percent") + 0
              major_faults = getv("major_page_faults") + 0
              minor_faults = getv("minor_page_faults") + 0
              fs_in = getv("fs_inputs") + 0
              fs_out = getv("fs_outputs") + 0
              n[stage] += 1
              total += 1
              elapsed_sum[stage] += elapsed
              user_cpu_sum[stage] += user_cpu
              system_cpu_sum[stage] += system_cpu
              cpu_pct_sum[stage] += cpu_pct
              major_fault_sum[stage] += major_faults
              minor_fault_sum[stage] += minor_faults
              fs_input_sum[stage] += fs_in
              fs_output_sum[stage] += fs_out
              input_sum[stage] += inputs
              output_sum[stage] += outputs
              bytes_sum[stage] += bytes
              if (!(stage in elapsed_min) || elapsed < elapsed_min[stage]) elapsed_min[stage] = elapsed
              if (elapsed > elapsed_max[stage]) elapsed_max[stage] = elapsed
              if (rss > rss_max[stage]) rss_max[stage] = rss
              if (req_mb > request_max[stage]) request_max[stage] = req_mb
              if (rss > overall_rss_max) overall_rss_max = rss
              if (req_mb > overall_request_max) overall_request_max = req_mb
              overall_elapsed += elapsed
              overall_user_cpu += user_cpu
              overall_system_cpu += system_cpu
              overall_cpu_pct += cpu_pct
              overall_major_faults += major_faults
              overall_minor_faults += minor_faults
              overall_fs_inputs += fs_in
              overall_fs_outputs += fs_out
              overall_inputs += inputs
              overall_outputs += outputs
              overall_bytes += bytes
              if (exit_code != "0") {
                failed += 1
                failed_stage[stage] += 1
                failed_lines = failed_lines "\nfailed_profile=" $0
              }
            }
            END {
              print "profile_rows=" total
              print "profile_failures=" failed + 0
              print "group_size=" group_size
              print "replay_group_size=" replay_group_size
              print "request_memory=" request_memory
              print "request_memory_mb=" request_memory_mb
              print "capture_request_memory=" capture_request_memory
              print "replay_request_memory=" replay_request_memory
              print "capture_request_memory_mb=" capture_request_memory_mb
              print "replay_request_memory_mb=" replay_request_memory_mb
              print "replay_output_roots_per_shard=" replay_output_roots_per_shard
              print "capture_nevents_per_worker=" capture_nevents
              print "smoke_capture_job_cap=" smoke_cap
              print "photon_id_rows=" photon_rows
              print "capture_cones=" capture_cones
              print "replay_cones=" replay_cones
              for (stage in n) {
                avg_elapsed = elapsed_sum[stage] / n[stage]
                max_mb = mb_from_kb(rss_max[stage])
                rec_mb = rec_mem_mb(max_mb)
                sec_per_input = (input_sum[stage] > 0 ? elapsed_sum[stage] / input_sum[stage] : 0)
                avg_cpu_pct = (n[stage] > 0 ? cpu_pct_sum[stage] / n[stage] : 0)
                stage_req = request_max[stage] + 0
                stage_usage = (stage_req > 0 && max_mb > 0 ? max_mb / stage_req : 0)
                printf("stage_summary stage=%s jobs=%d failed=%d elapsed_avg_s=%.1f elapsed_min_s=%d elapsed_max_s=%d input_files=%d sec_per_input_file=%.2f user_cpu_total_s=%.1f system_cpu_total_s=%.1f cpu_avg_pct=%.1f fs_inputs=%d fs_outputs=%d major_page_faults=%d minor_page_faults=%d max_rss_mb=%d requested_mb=%d usage_fraction=%.2f output_files=%d output_mb=%.1f headroom_request_memory_mb=%d\n",
                       stage, n[stage], failed_stage[stage] + 0, avg_elapsed,
                       elapsed_min[stage], elapsed_max[stage], input_sum[stage],
                       sec_per_input, user_cpu_sum[stage], system_cpu_sum[stage],
                       avg_cpu_pct, fs_input_sum[stage], fs_output_sum[stage],
                       major_fault_sum[stage], minor_fault_sum[stage],
                       max_mb, stage_req, stage_usage, output_sum[stage], bytes_sum[stage] / 1048576.0, rec_mb)
              }
              overall_max_mb = mb_from_kb(overall_rss_max)
              recommended_mb = rec_mem_mb(overall_max_mb)
              overall_sec_per_input = (overall_inputs > 0 ? overall_elapsed / overall_inputs : 0)
              overall_cpu_avg = (total > 0 ? overall_cpu_pct / total : 0)
              printf("overall_summary jobs=%d elapsed_total_s=%.1f input_files=%d sec_per_input_file=%.2f user_cpu_total_s=%.1f system_cpu_total_s=%.1f cpu_avg_pct=%.1f fs_inputs=%d fs_outputs=%d major_page_faults=%d minor_page_faults=%d max_rss_mb=%d output_files=%d output_mb=%.1f headroom_request_memory_mb=%d\n",
                     total, overall_elapsed, overall_inputs, overall_sec_per_input,
                     overall_user_cpu, overall_system_cpu, overall_cpu_avg,
                     overall_fs_inputs, overall_fs_outputs, overall_major_faults,
                     overall_minor_faults, overall_max_mb, overall_outputs,
                     overall_bytes / 1048576.0, recommended_mb)
              req_mb = (overall_request_max > 0 ? overall_request_max : int_or_zero(request_memory_mb))
              if (req_mb > 0 && overall_max_mb > 0) {
                usage = overall_max_mb / req_mb
                printf("memory_headroom requested_mb=%d max_rss_mb=%d usage_fraction=%.2f spare_mb=%d\n",
                       req_mb, overall_max_mb, usage, req_mb - overall_max_mb)
              } else {
                usage = 0
                print "memory_headroom requested_mb=unknown max_rss_mb=" overall_max_mb " usage_fraction=unknown spare_mb=unknown"
              }
              cap_events = int_or_zero(capture_nevents)
              smoke_jobs = int_or_zero(smoke_cap)
              if (cap_events == 0) {
                print "resource_confidence=high_for_this_exact_grouping_full_event_workers"
              } else if (cap_events < 10000) {
                print "resource_confidence=fast_pipeline_check_only"
              } else {
                print "resource_confidence=single_pass_tuning_smoke"
              }
              if (failed > 0) {
                print "measurement_status=invalid_worker_failures"
                print "measurement_note=Worker failures make timing and memory summaries incomplete. Inspect failed_profile lines and worker logs before inferring production settings."
                printf "%s\n", failed_lines
              } else {
                if (req_mb > 0 && usage > 0.90) {
                  print "measurement_status=complete_memory_tight"
                  print "measurement_note=Measurements were collected, but max_rss is close to requested memory. Use the report locally to infer whether to increase memory or reduce groupSize."
                } else if (cap_events > 0 && cap_events < 10000) {
                  print "measurement_status=complete_small_event_cap"
                  print "measurement_note=Measurements were collected, but the event cap is below the intended one-pass tuning default. Treat runtime scaling cautiously."
                } else if (total < 6) {
                  print "measurement_status=complete_low_profile_count"
                  print "measurement_note=Measurements were collected from too few worker profiles for robust groupSize inference."
                } else {
                  print "measurement_status=complete_one_pass_tuning"
                  print "measurement_note=Measurements are complete enough for local inference of production groupSize and request_memory."
                }
                print "one_pass_smoke_target=smokeTest is a capped fast pipeline check by default; set RJ_SMOKE_DATA_NEVENTS=0 or RJ_SMOKE_SIM_NEVENTS=0 for full worker-input tuning, and use legacy poolSmoke/RJ_POOL_SMOKE_NEVENTS for DATA-only capped pool checks."
              }
            }' "$profile_lines"
        } > "$profile_summary"
        echo
        echo "profile_summary_report:"
        cat "$profile_summary"
      else
        {
          echo "RECOILJETS_SMOKE_PROFILE_SUMMARY_V1"
          echo "workflow=${workflow}"
          echo "dataset=${dataset}"
          echo "profile_rows=0"
          echo "recommendation=No RECOILJETS_JOB_PROFILE_V1 lines were found. Inspect worker stdout/err paths from the submit files."
        } > "$profile_summary"
        echo
        echo "profile_summary_report:"
        cat "$profile_summary"
      fi
    fi
  else
    echo "manifest missing or empty: ${manifest}"
  fi
} > "$body"

mkdir -p "$report_dir"
[[ -s "$profile_lines" ]] && cp "$profile_lines" "$profile_raw"
cp "$body" "$report_summary"
[[ -s "$manifest" ]] && cp "$manifest" "$report_manifest"
[[ -n "$manifest_smoke_run_stats" && -s "$manifest_smoke_run_stats" ]] && cp "$manifest_smoke_run_stats" "$report_smoke_run_stats"
[[ -s "$profile_summary" ]] && cp "$profile_summary" "$report_profile"
[[ -s "$profile_raw" ]] && cp "$profile_raw" "$report_profile_rows"
[[ -s "$profile_summary" ]] && cp "$profile_summary" "$report_job_profile"
[[ -s "$profile_raw" ]] && cp "$profile_raw" "$report_job_profile_rows"
[[ -s "$heartbeat_lines" ]] && cp "$heartbeat_lines" "$report_job_heartbeat_rows"
if [[ ! -s "$report_job_heartbeat_rows" ]]; then
  {
    echo "RECOILJETS_JOB_HEARTBEAT_ROWS_V1"
    echo "status=no_heartbeat_rows_found"
    echo "note=Short jobs may finish before the first RJ_JOB_HEARTBEAT_SECONDS interval."
  } > "$report_job_heartbeat_rows"
fi
{
  echo "RECOILJETS_ROOT_OUTPUT_SUMMARY_V1"
  echo "dataset=${dataset}"
  echo "workflow=${workflow}"
  echo "pool_root_files=${pool_count}"
  echo "output_root_files=${out_count}"
  [[ -d "$out_base" ]] && find "$out_base" -type f -name '*.root' -printf '%s %p\n' 2>/dev/null | sort -nr | head -50 || true
} > "$report_root_summary"
{
  echo "RECOILJETS_HISTOGRAM_SANITY_SUMMARY_V1"
  echo "dataset=${dataset}"
  echo "workflow=${workflow}"
  echo "status=deferred_root_level_check"
  echo "note=Final DAG node confirmed ROOT file presence. Histogram-level semantic checks run after reports/ROOT outputs are pulled locally."
} > "$report_hist_summary"
{
  echo "RECOILJETS_FAILURE_SUMMARY_V1"
  echo "dataset=${dataset}"
  echo "workflow=${workflow}"
  echo "status=${status}"
  echo "status_note=${status_note}"
  echo "rescue_file_count=${rescue_count}"
  echo "dagman_out=${dagman_out}"
  echo "nodes_log=${nodes_log}"
  [[ -s "$profile_summary" ]] && grep -E 'profile_failures=|failed_profile=|measurement_status=' "$profile_summary" || true
} > "$report_failure_summary"
[[ -s "$profile_summary" ]] && cp "$profile_summary" "$report_merge_profile"
[[ -s "$profile_raw" ]] && cp "$profile_raw" "$report_merge_profile_rows"
profile_rows_json=0
failure_count_json=0
rss_max_json=0
elapsed_max_json=0
if [[ -s "$profile_summary" ]]; then
  profile_rows_json="$(awk -F= '$1=="profile_rows"{print $2; exit}' "$profile_summary")"
  failure_count_json="$(awk -F= '$1=="profile_failures"{print $2; exit}' "$profile_summary")"
  rss_max_json="$(awk -F'[ =]' '/overall_summary/ {for(i=1;i<=NF;i++) if($i=="max_rss_mb"){print $(i+1); exit}}' "$profile_summary")"
  elapsed_max_json="$(awk -F'[ =]' '/stage_summary/ {for(i=1;i<=NF;i++) if($i=="elapsed_max_s" && $(i+1)>m)m=$(i+1)} END{print m+0}' "$profile_summary")"
fi
rec_mem_json=0
if [[ "${rss_max_json:-0}" =~ ^[0-9]+$ && "${rss_max_json:-0}" -gt 0 ]]; then
  rec_mem_json=$(( (rss_max_json * 135 + 99) / 100 + 512 ))
fi
current_group_json="${manifest_group_size:-0}"
if [[ ! "$current_group_json" =~ ^[0-9]+$ || "$current_group_json" -le 0 ]]; then current_group_json=1; fi
cat > "$report_tuning_json" <<JSON
{
  "dataset": "${dataset}",
  "workflow": "${workflow}",
  "current_groupSize": ${current_group_json},
  "current_replay_groupSize": ${manifest_replay_group_size:-0},
  "current_request_memory_mb": ${manifest_request_memory_mb:-0},
  "current_capture_request_memory_mb": ${manifest_capture_request_memory_mb:-0},
  "current_replay_request_memory_mb": ${manifest_replay_request_memory_mb:-0},
  "replay_output_roots_per_shard": ${manifest_replay_output_roots_per_shard:-0},
  "profile_rows": ${profile_rows_json:-0},
  "failure_count": ${failure_count_json:-0},
  "rss_max_mb": ${rss_max_json:-0},
  "elapsed_max_s": ${elapsed_max_json:-0},
  "input_files_per_job": null,
  "output_bytes_per_job": null,
  "recommended_memory_floor_mb": ${rec_mem_json},
  "recommended_groupSize_candidates": [$(( current_group_json > 1 ? current_group_json / 2 : 1 )), ${current_group_json}, $(( (current_group_json * 3 + 1) / 2 )), $(( current_group_json * 2 ))],
  "notes": "Use this as tuning input only. Codex should infer final dataset-specific groupSize and request_memory after pulling reports."
}
JSON

cat "$body"

if [[ -n "$emails" ]]; then
  subject="[RecoilJets][${dataset}][${workflow}][${status}]"
  mail_recipients="${emails//,/ }"
  if command -v mail >/dev/null 2>&1; then
    # shellcheck disable=SC2086
    mail -s "$subject" $mail_recipients < "$body" || true
  elif command -v mailx >/dev/null 2>&1; then
    # shellcheck disable=SC2086
    mailx -s "$subject" $mail_recipients < "$body" || true
  else
    echo "[notify] mail/mailx not found; notification skipped for ${emails}" >&2
  fi
fi
rm -f "$body" "$profile_lines" "$heartbeat_lines"
EOS
  chmod +x "$exec_file"

  cat > "$sub_file" <<SUB
universe      = vanilla
executable    = ${exec_file}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/final_summary_${workflow}.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/final_summary_${workflow}.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/final_summary_${workflow}.\$(Cluster).\$(Process).err
request_memory= 512MB
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
arguments     = ${workflow} ${dataset} ${pool_base} ${out_base} ${manifest} ${emails:-NONE} ${dag_file:-NONE}
queue
SUB
  printf "%s\n" "$sub_file"
}

selection_mode_normalize() {
  selection_mode_normalize_for_key "" "$1"
}

selection_mode_normalize_for_key() {
  local key="$1"
  local mode
  mode="$(trim_ws "$2")"
  case "$key:$mode" in
    preselection:variantA|preselection:VariantA|preselection:varianta|preselection:newPPG12|preselection:NewPPG12|preselection:newppg12) echo "newPPG12"; return 0 ;;
    preselection:variantB|preselection:VariantB|preselection:variantb|preselection:noPreCriteria|preselection:NoPreCriteria|preselection:noprecriteria) echo "noPreCriteria"; return 0 ;;
    preselection:variantC|preselection:VariantC|preselection:variantc|preselection:onlyNPB|preselection:OnlyNPB|preselection:onlynpb) echo "onlyNPB"; return 0 ;;
    preselection:variantD|preselection:VariantD|preselection:variantd|preselection:refPlusNPB|preselection:RefPlusNPB|preselection:refplusnpb) echo "refPlusNPB"; return 0 ;;
    preselection:variantE|preselection:VariantE|preselection:variante|preselection:auauOnlyNPB|preselection:AuAuOnlyNPB|preselection:auauonlynpb) echo "auauOnlyNPB"; return 0 ;;
    tight:variantA|tight:VariantA|tight:varianta|tight:newPPG12|tight:NewPPG12|tight:newppg12) echo "newPPG12"; return 0 ;;
    tight:variantB|tight:VariantB|tight:variantb|tight:auauEmbeddedBDT|tight:AuAuEmbeddedBDT|tight:auauembeddedbdt) echo "auauEmbeddedBDT"; return 0 ;;
    nonTight:variantA|nonTight:VariantA|nonTight:varianta|nonTight:bdtSideband|nonTight:BDTSideband|nonTight:bdtsideband|nonTight:newPPG12|nonTight:NewPPG12|nonTight:newppg12) echo "newPPG12"; return 0 ;;
    nonTight:variantB|nonTight:VariantB|nonTight:variantb|nonTight:auauBDTSideband|nonTight:AuAuBDTSideband|nonTight:auaubdtsideband) echo "auauBDTSideband"; return 0 ;;
    nonTight:variantC|nonTight:VariantC|nonTight:variantc|nonTight:auauBDTComplement|nonTight:AuAuBDTComplement|nonTight:auaubdtcomplement) echo "auauBDTComplement"; return 0 ;;
  esac
  case "$mode" in
    ""|reference|Reference) echo "reference" ;;
    variantA|VariantA|varianta|newPPG12|NewPPG12|newppg12) echo "newPPG12" ;;
    auauEmbeddedBDT|AuAuEmbeddedBDT|auauembeddedbdt) echo "auauEmbeddedBDT" ;;
    auauBDTSideband|AuAuBDTSideband|auaubdtsideband) echo "auauBDTSideband" ;;
    auauBDTComplement|AuAuBDTComplement|auaubdtcomplement) echo "auauBDTComplement" ;;
    variantB|VariantB|variantb) echo "variantB" ;;
    variantC|VariantC|variantc) echo "variantC" ;;
    variantD|VariantD|variantd) echo "variantD" ;;
    variantE|VariantE|variante) echo "variantE" ;;
    centINDcontrol|CentINDControl|centindcontrol|centINDControl) echo "centINDcontrol" ;;
    centAsFeat|CentAsFeat|centasfeat|centAsFeature) echo "centAsFeat" ;;
    centDepBDTs|CentDepBDTs|centdepbdts|centDepBDT) echo "centDepBDTs" ;;
    *) echo "$mode" ;;
  esac
}

selection_mode_tag() {
  local key="$1"
  local mode
  mode="$(selection_mode_normalize_for_key "$key" "$2")"
  case "$mode" in
    reference) echo "${key}Reference" ;;
    newPPG12) echo "${key}NewPPG12" ;;
    noPreCriteria) echo "${key}NoPreCriteria" ;;
    onlyNPB) echo "${key}OnlyNPB" ;;
    refPlusNPB) echo "${key}RefPlusNPB" ;;
    auauOnlyNPB) echo "${key}AuAuOnlyNPB" ;;
    auauEmbeddedBDT) echo "${key}AuAuEmbeddedBDT" ;;
    auauBDTSideband) echo "${key}AuAuBDTSideband" ;;
    auauBDTComplement) echo "${key}AuAuBDTComplement" ;;
    variantB) echo "${key}VariantB" ;;
    variantC) echo "${key}VariantC" ;;
    variantD) echo "${key}VariantD" ;;
    variantE) echo "${key}VariantE" ;;
    *)
      local first="${mode:0:1}"
      local rest="${mode:1}"
      echo "${key}${first^^}${rest}"
      ;;
  esac
}

yaml_get_photon_id_sets() {
  local yaml="$1"
  awk '
    function trim(s){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }
    BEGIN{inset=0}
    {
      line=$0
      sub(/#.*/, "", line)
      if (line ~ /^[[:space:]]*photon_id_sets[[:space:]]*:/) { inset=1; next }
      if (inset && line ~ /^[[:alnum:]_][[:alnum:]_[:space:]-]*:/) { exit }
      if (!inset || line !~ /\[/) next
      sub(/^.*\[/, "", line)
      sub(/\].*$/, "", line)
      n=split(line, a, ",")
      if (n >= 3) print trim(a[1]) "|" trim(a[2]) "|" trim(a[3])
    }' "$yaml"
}

# Build parallel arrays:
#   iso_tags[]            = full iso + selection label
#   iso_base_tags[]       = iso-only tag (placed before clusterUEpipeline in cfg_tag)
#   iso_selection_tags[]  = selection-only tag (placed after clusterUEpipeline in cfg_tag)
#   iso_sliding[], iso_fixed[]
#   iso_preselection[], iso_tight[], iso_nonTight[]
# from YAML keys:
#   isSlidingAndFixed, isSlidingIso, fixedGeV[],
#   preselection[], tight[], nonTight[]
build_iso_modes() {
  local yaml_file="$1"
  iso_tags=()
  iso_base_tags=()
  iso_selection_tags=()
  iso_sliding=()
  iso_fixed=()
  iso_preselection=()
  iso_tight=()
  iso_nonTight=()

  local -a _base_tags
  local -a _base_sliding
  local -a _base_fixed
  _base_tags=()
  _base_sliding=()
  _base_fixed=()

  local _both; _both="$(yaml_get_values "isSlidingAndFixed" "$yaml_file" 2>/dev/null || echo "false")"
  local _slide; _slide="$(yaml_get_values "isSlidingIso" "$yaml_file" 2>/dev/null || echo "false")"
  local -a _fixeds
  mapfile -t _fixeds < <( yaml_get_values "fixedGeV" "$yaml_file" 2>/dev/null )
  (( ${#_fixeds[@]} )) || _fixeds=( "2.0" )

  local -a _id_rows
  mapfile -t _id_rows < <( yaml_get_photon_id_sets "$yaml_file" )
  (( ${#_id_rows[@]} > 0 )) || { err "YAML must define photon_id_sets with rows like [preselection, tight, nonTight]: $yaml_file"; exit 74; }

  _both="$(trim_ws "$_both")"
  _slide="$(trim_ws "$_slide")"

  if [[ "$_both" == "true" ]]; then
    _base_tags+=( "$(sim_iso_tag "true" "0")" )
    _base_sliding+=( "true" )
    _base_fixed+=( "${_fixeds[0]}" )
    for fv in "${_fixeds[@]}"; do
      _base_tags+=( "$(sim_iso_tag "false" "$fv")" )
      _base_sliding+=( "false" )
      _base_fixed+=( "$fv" )
    done
  elif [[ "$_slide" == "true" ]]; then
    _base_tags+=( "$(sim_iso_tag "true" "0")" )
    _base_sliding+=( "true" )
    _base_fixed+=( "${_fixeds[0]}" )
  else
    for fv in "${_fixeds[@]}"; do
      _base_tags+=( "$(sim_iso_tag "false" "$fv")" )
      _base_sliding+=( "false" )
      _base_fixed+=( "$fv" )
    done
  fi

  local _i
  local _pre
  local _tight
  local _nonTight
  for (( _i=0; _i<${#_base_tags[@]}; _i++ )); do
    local _row
    for _row in "${_id_rows[@]}"; do
      IFS='|' read -r _pre _tight _nonTight <<< "$_row"
      local _pre_norm
      local _tight_norm
      local _nonTight_norm
      local selection_tag
      _pre_norm="$(selection_mode_normalize_for_key "preselection" "$_pre")"
      _tight_norm="$(selection_mode_normalize_for_key "tight" "$_tight")"
      _nonTight_norm="$(selection_mode_normalize_for_key "nonTight" "$_nonTight")"
      selection_tag="$(selection_mode_tag "preselection" "$_pre_norm")_$(selection_mode_tag "tight" "$_tight_norm")_$(selection_mode_tag "nonTight" "$_nonTight_norm")"

      iso_tags+=( "${_base_tags[$_i]}_${selection_tag}" )
      iso_base_tags+=( "${_base_tags[$_i]}" )
      iso_selection_tags+=( "${selection_tag}" )
      iso_sliding+=( "${_base_sliding[$_i]}" )
      iso_fixed+=( "${_base_fixed[$_i]}" )
      iso_preselection+=( "${_pre_norm}" )
      iso_tight+=( "${_tight_norm}" )
      iso_nonTight+=( "${_nonTight_norm}" )
    done
  done
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

yaml_set_scalar_in_place() {
  local file="$1" key="$2" value="$3"
  local tmp="${file}.tmp.$$"
  if grep -Eq "^[[:space:]]*${key}:" "$file"; then
    sed -E "s|^([[:space:]]*${key}:).*|\\1 ${value}|" "$file" > "$tmp"
    mv "$tmp" "$file"
  else
    {
      printf '\n'
      printf '%s: %s\n' "$key" "$value"
    } >> "$file"
  fi
}

pin_photon_id_scalars_in_yaml() {
  local file="$1" preselection="$2" tight="$3" nonTight="$4"
  yaml_set_scalar_in_place "$file" "preselection" "$preselection"
  yaml_set_scalar_in_place "$file" "tight" "$tight"
  yaml_set_scalar_in_place "$file" "nonTight" "$nonTight"
}

sim_make_yaml_override() {
  local master="$1" pt="$2" frac="$3" vz="$4" cone="$5" sliding="$6" fixed="$7" uepipe="$8" preselection="$9" tight="${10}" nonTight="${11}" tag="${12}" stamp="${13:-}" force_sliding_and_fixed="${14:-}"
  mkdir -p "$SIM_YAML_OVERRIDE_DIR"
  local out="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${tag}${stamp:+_${stamp}}.yaml"

  [[ -s "$master" ]] || { err "Master YAML not found or empty: $master"; exit 71; }
  grep -Eq '^[[:space:]]*jet_pt_min:' "$master" || { err "YAML missing key 'jet_pt_min' in $master"; exit 71; }
  grep -Eq '^[[:space:]]*back_to_back_dphi_min_pi_fraction:' "$master" || { err "YAML missing key 'back_to_back_dphi_min_pi_fraction' in $master"; exit 71; }
  grep -Eq '^[[:space:]]*vz_cut_cm:' "$master" || { err "YAML missing key 'vz_cut_cm' in $master"; exit 71; }
  grep -Eq '^[[:space:]]*coneR:' "$master" || { err "YAML missing key 'coneR' in $master"; exit 71; }

  local -a sed_args
  sed_args=(
    -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${pt}|"
    -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${frac}|"
    -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${vz}|"
    -e "s|^([[:space:]]*coneR:).*|\\1 ${cone}|"
    -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${sliding}|"
    -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${fixed}|"
    -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${uepipe}|"
  )
  if [[ -n "$force_sliding_and_fixed" ]]; then
    sed_args+=( -e "s|^([[:space:]]*isSlidingAndFixed:).*|\\1 ${force_sliding_and_fixed}|" )
  fi

  sed -E "${sed_args[@]}" "$master" > "$out"
  pin_photon_id_scalars_in_yaml "$out" "$preselection" "$tight" "$nonTight"

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
      DEST_BASE="${RJ_SIMEMBED_DEST_BASE:-$SIMEMBED_DEST_BASE}"
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
      DEST_BASE="${RJ_SIMEMBED_DEST_BASE:-$SIMEMBED_DEST_BASE}"
      TAG="simembeddedinclusive"
      MACRO="${BASE}/macros/Fun4All_recoilJets_AuAu.C"
      EXE="${BASE}/RecoilJets_Condor_AuAu.sh"
      IS_SIM=1
      SIM_SAMPLE_DEFAULT="run28_embeddedJet20"
      SIM_SAMPLE="$SIM_SAMPLE_DEFAULT"
      ;;
    isSimEmbeddedAndInclusive|issimembeddedandinclusive|simembeddedandinclusive|SIMEMBEDDEDANDINCLUSIVE)
      DATASET="isSimEmbeddedAndInclusive"
      GOLDEN=""
      LIST_DIR=""
      LIST_PREFIX=""
      DEST_BASE="$AUAU_BDT_DEST_BASE"
      TAG="simembeddedandinclusive"
      MACRO="${BASE}/macros/Fun4All_recoilJets_AuAu.C"
      EXE="${BASE}/RecoilJets_Condor_AuAu.sh"
      IS_SIM=1
      SIM_SAMPLE_DEFAULT="run28_embeddedPhoton20"
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
  if [[ -n "${RJ_DEST_BASE_OVERRIDE:-}" ]]; then
    DEST_BASE="$RJ_DEST_BASE_OVERRIDE"
  fi
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

  if [[ -n "${SIM_STAGE_NAMESPACE:-}" ]]; then
    SIM_STAGE_DIR="${STAGE_DIR}/${SIM_STAGE_NAMESPACE}/${SIM_SAMPLE}"
  else
    SIM_STAGE_DIR="${STAGE_DIR}/${SIM_SAMPLE}"
  fi
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
      isSimEmbedded)          samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
      isSimEmbeddedInclusive) samples=( "run28_embeddedJet12" "run28_embeddedJet20" ) ;;
      isSimJet5)              samples=( "run28_jet5" ) ;;
      isSimMB)                samples=( "run28_detroit" ) ;;
      *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
    esac
  else
    samples=( "${SIM_SAMPLE}" )
  fi

  n_cfg=$(( ${#sim_pts[@]} * ${#sim_fracs[@]} * ${#sim_vzs[@]} * ${#sim_cones[@]} * $(iso_group_count) * ${#uepipe_modes[@]} ))

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
  say "  iso groups submitted              : $(iso_group_count) upstream mode(s)"
  say "  photon-ID fanout outputs          : ${#iso_tags[@]} cfg tag(s) across those iso groups"
  if (( ${uepipe_in_tag:-0} )); then
    say "  clusterUEpipeline                 : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} values; tagged)"
  else
    say "  clusterUEpipeline                 : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} forced value; not tagged for ${TAG})"
  fi
  say "  upstream DST jobs configs         : ${BOLD}${n_cfg}${RST}"
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

  say "${BOLD}Full upstream ID-fanout cell list (${n_cfg} entries × ${#samples[@]} samples = ${BOLD}$((n_cfg * ${#samples[@]}))${RST} submit blocks):${RST}"
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
              _tag="$(matrix_cfg_tag "$_pt" "$_frac" "$_vz" "$_cone" "$_ci" "$_ue")"
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
  say "  upstream configs       : ${n_cfg}"
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
  _ex_tag="$(matrix_cfg_tag "${sim_pts[0]}" "${sim_fracs[0]}" "${sim_vzs[0]}" "${sim_cones[0]}" 0 "${uepipe_modes[0]}")"
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

  local n_matrix=$(( ${#ck_pts[@]} * ${#ck_fracs[@]} * ${#ck_vzs[@]} * ${#ck_cones[@]} * $(iso_group_count) * ${#uepipe_modes[@]} ))
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
  say "  iso groups submitted : $(iso_group_count) upstream mode(s)"
  say "  photon-ID fanout cfgs: ${#iso_tags[@]} cfg tag(s) across those iso groups"
  if (( ${uepipe_in_tag:-0} )); then
    say "  clusterUEpipeline    : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} values; tagged)"
  else
    say "  clusterUEpipeline    : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} forced value; not tagged for ${TAG})"
  fi
  say "  upstream matrix cells: ${BOLD}${n_matrix}${RST}"
  echo
  say "${BOLD}Per-matrix-cell base:${RST}"
  say "  golden runs (w/lists): ${total_runs}"
  say "  total input files    : ${total_files}"
  say "  base jobs per combo  : ${base_jobs}"
  (( missing > 0 )) && warn "  runs skipped (missing lists): ${missing}"
  echo

  say "${BOLD}Full upstream ID-fanout cell list (${n_matrix} entries):${RST}"
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
    iso_idx_is_group_leader "$_ci" || continue
    for _ue in "${uepipe_modes[@]}"; do
      (( ++cfg_num ))
      local _dpt; _dpt="jetMinPt$(sim_pt_tag "$_pt")"
      local _dfrac; _dfrac="$(sim_b2b_tag "$_frac")"
      local _dvz; _dvz="$(sim_vz_tag "$_vz")"
      local _dcone; _dcone="$(sim_cone_tag "$_cone")"
      local _tag
      _tag="$(matrix_cfg_tag "$_pt" "$_frac" "$_vz" "$_cone" "$_ci" "$_ue")"
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
  say "  upstream matrix cells  : ${n_matrix}"
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
  local _ex_tag
  _ex_tag="$(matrix_cfg_tag "${ck_pts[0]}" "${ck_fracs[0]}" "${ck_vzs[0]}" "${ck_cones[0]}" 0 "${uepipe_modes[0]}")"
  say "  e.g. ${DIM}${DEST_BASE}/${_ex_tag}/00067599/*.root${RST}"
}

workflow_check() {
  local yaml_src
  yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  [[ -s "$yaml_src" ]] || { err "Workflow check YAML not found or empty: $yaml_src"; exit 72; }

  local -a pts fracs vzs cones
  mapfile -t pts   < <( yaml_get_values "jet_pt_min" "$yaml_src" )
  mapfile -t fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$yaml_src" )
  mapfile -t vzs   < <( yaml_get_values "vz_cut_cm" "$yaml_src" )
  mapfile -t cones < <( yaml_get_values "coneR" "$yaml_src" )
  (( ${#pts[@]} ))   || { err "No values found for jet_pt_min in $yaml_src"; exit 72; }
  (( ${#fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $yaml_src"; exit 72; }
  (( ${#vzs[@]} ))   || { err "No values found for vz_cut_cm in $yaml_src"; exit 72; }
  (( ${#cones[@]} )) || { err "No values found for coneR in $yaml_src"; exit 72; }

  build_iso_modes "$yaml_src"
  read_uepipe_modes "$yaml_src" "$TAG"

  local n_cfg=$(( ${#pts[@]} * ${#fracs[@]} * ${#vzs[@]} * ${#cones[@]} * $(iso_group_count) * ${#uepipe_modes[@]} ))
  local ex_tag
  ex_tag="$(matrix_cfg_tag "${pts[0]}" "${fracs[0]}" "${vzs[0]}" "${cones[0]}" 0 "${uepipe_modes[0]}")"
  local override
  override="$(sim_make_yaml_override "$yaml_src" "${pts[0]}" "${fracs[0]}" "${vzs[0]}" "${cones[0]}" "${iso_sliding[0]}" "${iso_fixed[0]}" "${uepipe_modes[0]}" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}" "$ex_tag" "WORKFLOWCHECK" "false")"

  say "${BOLD}Workflow check${RST}"
  say "  dataset          : ${DATASET}"
  say "  YAML             : ${yaml_src}"
  say "  photon_id_sets   : ${#iso_selection_tags[@]} ID/iso rows after expansion"
  say "  matrix cfg count : ${n_cfg}"
  say "  example cfg_tag  : ${ex_tag}"
  say "  example override : ${override}"
  say "  override scalars : preselection=$(grep -E '^[[:space:]]*preselection:' "$override" | tail -1 | awk -F: '{print $2}' | xargs)"
  say "                     tight=$(grep -E '^[[:space:]]*tight:' "$override" | tail -1 | awk -F: '{print $2}' | xargs)"
  say "                     nonTight=$(grep -E '^[[:space:]]*nonTight:' "$override" | tail -1 | awk -F: '{print $2}' | xargs)"
  say "  notify_emails    : $(grep -E '^[[:space:]]*notify_emails:' "$yaml_src" | head -1 | cut -d: -f2- | xargs || true)"
  say "  note             : pool replay is active; from-scratch workflows capture only DST-dependent axes, then fan out replay-only ID/cut/bin variations."
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
  local args_file="${SUB_DIR}/RecoilJets_${TAG}_${stamp}.args"
  local submit_stage_dir="${STAGE_DIR}/${TAG}_${stamp}"
  local exe_to_use="${BULK_FROZEN_EXE:-${EXE}}"
  local macro_env=""
  [[ -n "${BULK_FROZEN_MACRO:-}" ]] && macro_env=";RJ_MACRO_PATH=${BULK_FROZEN_MACRO}"
  local submit_extra_env="${RJ_SUBMIT_EXTRA_ENV:-}"
  [[ -n "$submit_extra_env" && "$submit_extra_env" != \;* ]] && submit_extra_env=";${submit_extra_env}"
  local request_memory="${RJ_REQUEST_MEMORY:-2000MB}"
  mkdir -p "$submit_stage_dir"
  : > "$args_file"
  say "Submit chunk-list stage: ${submit_stage_dir}"

  # Snapshot YAML at submit time so idle jobs are immune to later edits
  local yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  local yaml_snap="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${stamp}.yaml"
  local source_runs
  source_runs=$(grep -cE '^[0-9]+' "$source" 2>/dev/null || true)
  mkdir -p "$SIM_YAML_OVERRIDE_DIR"
  cp -f "$yaml_src" "$yaml_snap"
  say "YAML snapshot: ${yaml_snap}"
  say "Submit context: source=${source}  runs=${source_runs:-0}  groupSize=${GROUP_SIZE}  firstChunk=${first_chunk:-none}"
  say "Submit environment: RJ_DATASET=${DATASET}  RJ_VERBOSITY=0  RJ_CONFIG_YAML=${yaml_snap}${macro_env}${submit_extra_env}"

  cat > "$sub" <<SUB
universe      = vanilla
executable    = ${exe_to_use}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/job.\$(Cluster).\$(Process).err
request_memory= ${request_memory}
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
# Force dataset & quiet macro on Condor (YAML frozen at submit time):
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=0;RJ_CONFIG_YAML=${yaml_snap}${macro_env}${submit_extra_env}
queue arguments from ${args_file}
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

    # Build group lists for this run in this submit's frozen staging directory.
    # This prevents later submissions from overwriting chunk lists still used by queued jobs.
    mapfile -t groups < <( STAGE_DIR="$submit_stage_dir" make_groups "$r8" "$GROUP_SIZE" )

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
      printf '%s %s %s $(Cluster) 0 %d NONE %s\n' \
             "$r8" "$glist" "$DATASET" "$gidx" "$DEST_BASE" >> "$args_file"
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

data_analysis_tag_for_dataset() {
  case "$1" in
    isAuAu) echo "isAuAu" ;;
    isOO) echo "isOO" ;;
    *) echo "isPP" ;;
  esac
}

select_largest_stat_data_runs() {
  local count="${1:-10}"
  local out_file="${2:?output run list required}"
  local stats_file="${3:-}"
  [[ "$count" =~ ^[0-9]+$ && "$count" -gt 0 ]] || count=10
  mkdir -p "$(dirname "$out_file")"
  : > "$out_file"
  [[ -z "$stats_file" ]] || { mkdir -p "$(dirname "$stats_file")"; : > "$stats_file"; }
  local ranked
  ranked="$(mktemp "${TMPDIR:-/tmp}/recoiljets_smoke_runs.XXXXXX")"
  local ranked_top
  ranked_top="$(mktemp "${TMPDIR:-/tmp}/recoiljets_smoke_runs_top.XXXXXX")"
  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    local r8 src nfiles
    r8="$(run8 "$rn")"
    if [[ -n "${TRIGGER_BIT}" ]] && ! is_trigger_active "$r8" "$TRIGGER_BIT"; then
      continue
    fi
    src="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
    [[ -s "$src" ]] || continue
    nfiles="$(grep -Evc '^[[:space:]]*($|#)' "$src" 2>/dev/null || echo 0)"
    [[ "$nfiles" =~ ^[0-9]+$ && "$nfiles" -gt 0 ]] || continue
    printf '%012d %s\n' "$nfiles" "$r8" >> "$ranked"
  done < "$GOLDEN"
  sort -nr "$ranked" | head -n "$count" > "$ranked_top"
  awk '{print $2}' "$ranked_top" > "$out_file"
  if [[ -n "$stats_file" ]]; then
    awk '{printf "run=%s input_files=%d\n", $2, $1+0}' "$ranked_top" > "$stats_file"
  fi
  rm -f "$ranked" "$ranked_top"
  [[ -s "$out_file" ]]
}

submit_data_pool_workflow() {
  local from_scratch="$1"  # 1: DST->pool->replay DAG, 0: replay existing pools
  local smoke_capture_cap="${2:-0}"
  local workflow_flavor="${3:-}"
  local master_yaml="$data_yaml_src"
  local workflow_stamp
  workflow_stamp="$(date +%Y%m%d_%H%M%S)"
  local dag_dryrun=0
  if [[ "${RJ_DAG_DRYRUN:-0}" == "1" || "${RJ_DAG_DRYRUN:-0}" == "true" || "${RJ_DAG_DRYRUN:-0}" == "TRUE" ]]; then
    dag_dryrun=1
  fi
  [[ "$smoke_capture_cap" =~ ^[0-9]+$ ]] || smoke_capture_cap=0
  local is_smoke=0
  if [[ "$workflow_flavor" == "poolSmoke" || "$workflow_flavor" == "smokeTest" || "$workflow_flavor" == "smoke" || "$smoke_capture_cap" -gt 0 ]]; then
    is_smoke=1
  fi
  local trace_pool=0
  if (( is_smoke || dag_dryrun )) || [[ "${RJ_POOL_TRACE:-0}" == "1" || "${RJ_POOL_TRACE:-0}" == "true" || "${RJ_POOL_TRACE:-0}" == "TRUE" ]]; then
    trace_pool=1
  fi
  local trace_every="${RJ_POOL_TRACE_EVERY:-25}"
  [[ "$trace_every" =~ ^[0-9]+$ && "$trace_every" -gt 0 ]] || trace_every=25
  local pool_base out_base
  if (( is_smoke )); then
    if [[ "$workflow_flavor" == "smokeTest" ]]; then
      pool_base="${RJ_POOL_OUTPUT_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaPoolsSmoke/${TAG}_smokeTest_${workflow_stamp}}"
      out_base="${RJ_POOL_SMOKE_OUTPUT_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaSmoke/${TAG}_smokeTest_${workflow_stamp}}"
    else
      pool_base="${RJ_POOL_OUTPUT_BASE:-${POOL_DEST_ROOT}/${TAG}_poolSmoke_${workflow_stamp}}"
      out_base="${RJ_POOL_SMOKE_OUTPUT_BASE:-${DATA_DEST_BASE_SAVED%/}/_poolSmoke_${workflow_stamp}}"
    fi
  else
    pool_base="${RJ_POOL_INPUT_BASE:-${RJ_POOL_OUTPUT_BASE:-${POOL_DEST_ROOT}/${TAG}}}"
    out_base="${DATA_DEST_BASE_SAVED%/}"
  fi
  local dag_dir="${SUB_DIR}/pool_workflow_${TAG}_${workflow_stamp}"
  local dag="${dag_dir}/RecoilJets_pool_${DATASET}_${workflow_stamp}.dag"
  local manifest="${dag_dir}/manifest.txt"
  local pool_macro="${BASE}/macros/Fun4All_recoilJets_poolReplay.C"
  local replay_gs="${RJ_POOL_REPLAY_GROUP_SIZE:-}"
  if [[ -z "$replay_gs" ]]; then
    case "$DATASET" in
      isAuAu|isOO) replay_gs="${RJ_DATA_REPLAY_GROUPSIZE_AUAU:-7}" ;;
      *)           replay_gs="${RJ_DATA_REPLAY_GROUPSIZE_PP:-$GROUP_SIZE}" ;;
    esac
  fi
  if (( is_smoke )) && [[ -z "${RJ_POOL_REPLAY_GROUP_SIZE:-}" ]] && (( smoke_capture_cap > 0 )); then
    replay_gs="$smoke_capture_cap"
  fi
  [[ "$replay_gs" =~ ^[0-9]+$ && "$replay_gs" -gt 0 ]] || replay_gs="$GROUP_SIZE"
  local profile_job="${RJ_PROFILE_JOB:-0}"
  (( is_smoke )) && profile_job=1
  local capture_nevents=0
  if (( is_smoke )); then
    if [[ "$workflow_flavor" == "smokeTest" ]]; then
      capture_nevents="${RJ_SMOKE_DATA_NEVENTS:-3000}"
    else
      capture_nevents="${RJ_POOL_SMOKE_NEVENTS:-20000}"
    fi
  fi
  [[ "$capture_nevents" == "-1" || "$capture_nevents" =~ ^[0-9]+$ ]] || { err "Smoke nEvents must be -1 or a non-negative integer, got '${capture_nevents}'"; exit 2; }
  local request_memory="${RJ_REQUEST_MEMORY:-2000MB}"
  local request_memory_digits="${request_memory//[^0-9]/}"
  local request_memory_mb="${request_memory_digits:-0}"
  case "$request_memory" in
    *[Gg][Bb]|*[Gg]) request_memory_mb=$(( request_memory_mb * 1024 )) ;;
  esac
  local capture_request_memory="${RJ_CAPTURE_REQUEST_MEMORY:-$request_memory}"
  local replay_request_memory="${RJ_REPLAY_REQUEST_MEMORY:-}"
  if [[ -z "$replay_request_memory" ]]; then
    case "$DATASET" in
      isAuAu|isOO) replay_request_memory="${RJ_DEFAULT_REPLAY_REQUEST_MEMORY_AUAU:-2500MB}" ;;
      *)           replay_request_memory="${RJ_DEFAULT_REPLAY_REQUEST_MEMORY_DATA:-1800MB}" ;;
    esac
  fi
  local capture_request_memory_mb replay_request_memory_mb
  capture_request_memory_mb="$(memory_request_to_mb "$capture_request_memory")"
  replay_request_memory_mb="$(memory_request_to_mb "$replay_request_memory")"
  local replay_roots_per_shard="${RJ_REPLAY_OUTPUT_ROOTS_PER_SHARD:-${RJ_POOL_REPLAY_OUTPUT_ROOTS_PER_SHARD:-999999}}"
  [[ "$replay_roots_per_shard" =~ ^[0-9]+$ && "$replay_roots_per_shard" -gt 0 ]] || replay_roots_per_shard=999999
  local replay_max_open_outputs="${RJ_REPLAY_MAX_OPEN_OUTPUTS:-$replay_roots_per_shard}"
  local job_heartbeat_seconds="${RJ_JOB_HEARTBEAT_SECONDS:-0}"
  if (( is_smoke )) && [[ "$job_heartbeat_seconds" == "0" ]]; then
    job_heartbeat_seconds="${RJ_SMOKE_JOB_HEARTBEAT_SECONDS:-120}"
  fi
  local capture_submit_extra replay_submit_extra
  capture_submit_extra="$(condor_pool_submit_extra_lines capture "$is_smoke")"
  replay_submit_extra="$(condor_pool_submit_extra_lines replay "$is_smoke")"

  if [[ "${RJ_DAG_DRYRUN:-0}" != "1" && "${RJ_DAG_DRYRUN:-0}" != "true" && "${RJ_DAG_DRYRUN:-0}" != "TRUE" ]]; then
    need_cmd condor_submit_dag
  fi
  [[ -f "$pool_macro" ]] || { err "Pool replay macro not found: $pool_macro"; exit 92; }
  mkdir -p "$dag_dir" "$SIM_YAML_OVERRIDE_DIR"
  : > "$dag"

  read_replay_cones "$master_yaml"
  read_capture_cones "$master_yaml"
  local -a data_replay_cones=( "${replay_cones[@]}" )
  local -a data_capture_cones=( "${capture_cones[@]}" )
  local workflow_golden="$GOLDEN"
  local smoke_selected_runs=""
  local smoke_selected_run_stats=""
  if [[ "$workflow_flavor" == "smokeTest" ]]; then
    local smoke_run_count="${RJ_SMOKE_DATA_RUNS:-10}"
    smoke_selected_runs="${dag_dir}/smoke_selected_runs_${TAG}.list"
    smoke_selected_run_stats="${dag_dir}/smoke_selected_runs_${TAG}_stats.txt"
    select_largest_stat_data_runs "$smoke_run_count" "$smoke_selected_runs" "$smoke_selected_run_stats" || { err "smokeTest could not select any DATA runs from ${GOLDEN}"; exit 99; }
    workflow_golden="$smoke_selected_runs"
  fi
  local golden_rows=0
  [[ -s "$workflow_golden" ]] && golden_rows="$(grep -Evc '^[[:space:]]*($|#)' "$workflow_golden" 2>/dev/null || echo 0)"
  local smoke_capture_jobs_per_run=0
  if (( is_smoke )); then
    smoke_capture_jobs_per_run="${RJ_SMOKE_CAPTURE_JOBS_PER_RUN:-1}"
    [[ "$smoke_capture_jobs_per_run" =~ ^[0-9]+$ ]] || smoke_capture_jobs_per_run=1
  fi

  if [[ "$from_scratch" -eq 1 ]]; then
    cleanup_bulk_snapshots_for_tag
    case "$DATASET" in
      isAuAu|isOO) create_pipeline_snapshot "auau" "$workflow_stamp" ;;
      *)           create_pipeline_snapshot "pp" "$workflow_stamp" ;;
    esac
  fi

  {
    echo "workflow=$([[ "$from_scratch" -eq 1 ]] && echo dataAllFromScratch || echo dataAllFromPool)"
    echo "dataset=${DATASET}"
    echo "smoke=${is_smoke}"
    echo "smoke_capture_job_cap=${smoke_capture_cap}"
    (( is_smoke )) && echo "smoke_capture_jobs_per_run=${smoke_capture_jobs_per_run}"
    [[ -n "$smoke_selected_runs" ]] && echo "smoke_selected_runs_file=${smoke_selected_runs}"
    [[ -n "$smoke_selected_runs" ]] && echo "smoke_selected_runs=${golden_rows}"
    [[ -n "$smoke_selected_run_stats" ]] && echo "smoke_selected_run_stats=${smoke_selected_run_stats}"
    echo "capture_nevents=${capture_nevents}"
    echo "replay_nevents=0"
    echo "profile_job=${profile_job}"
    echo "pool_profile_counts=1"
    echo "job_heartbeat_seconds=${job_heartbeat_seconds}"
    echo "condor_use_late_materialization=${RJ_CONDOR_USE_LATE_MATERIALIZATION:-1}"
    echo "condor_capture_max_idle=${RJ_CONDOR_CAPTURE_MAX_IDLE:-${RJ_CONDOR_MAX_IDLE:-2000}}"
    echo "condor_capture_max_materialize=${RJ_CONDOR_CAPTURE_MAX_MATERIALIZE:-${RJ_CONDOR_MAX_MATERIALIZE:-5000}}"
    echo "condor_replay_max_idle=${RJ_CONDOR_REPLAY_MAX_IDLE:-${RJ_CONDOR_MAX_IDLE:-1000}}"
    echo "condor_replay_max_materialize=${RJ_CONDOR_REPLAY_MAX_MATERIALIZE:-${RJ_CONDOR_MAX_MATERIALIZE:-3000}}"
    echo "request_memory=${request_memory}"
    echo "request_memory_mb=${request_memory_mb}"
    echo "capture_request_memory=${capture_request_memory}"
    echo "replay_request_memory=${replay_request_memory}"
    echo "capture_request_memory_mb=${capture_request_memory_mb}"
    echo "replay_request_memory_mb=${replay_request_memory_mb}"
    echo "yaml=${master_yaml}"
    echo "pool_base=${pool_base}"
    echo "output_base=${out_base}"
    echo "group_size=${GROUP_SIZE}"
    echo "replay_group_size=${replay_gs}"
    echo "replay_output_roots_per_shard=${replay_roots_per_shard}"
    echo "replay_max_open_outputs=${replay_max_open_outputs}"
    echo "replay_cones=${data_replay_cones[*]}"
    echo "capture_cones=${data_capture_cones[*]}"
    echo "uepipe_modes=${uepipe_modes[*]}"
    echo "photon_id_sets_expanded=${#iso_selection_tags[@]}"
    echo "submitted_at=${workflow_stamp}"
  } > "$manifest"

  declare -A cap_node_by_tag=()
  declare -A pool_list_by_tag_run=()
  local capture_count=0
  local replay_count=0
  local capture_job_total=0
  local replay_job_total=0
  local analysis_tag
  analysis_tag="$(data_analysis_tag_for_dataset "$DATASET")"

  say "${BOLD}DATA pool workflow requested${RST}"
  say "  mode        : $([[ "$from_scratch" -eq 1 ]] && echo allFromScratch || echo all-from-existing-pools)$([[ "$is_smoke" -eq 1 ]] && echo ' smoke' || true)"
  say "  dataset     : ${DATASET}"
  say "  pool base   : ${pool_base}"
  say "  output base : ${out_base}"
  say "  memory      : capture=${capture_request_memory} replay=${replay_request_memory}"
  say "  replay split: poolFiles/job=${replay_gs} outputRoots/shard=${replay_roots_per_shard} maxOpenOutputs=${replay_max_open_outputs}"
  (( is_smoke )) && say "  smoke cap   : ${smoke_capture_cap} capture jobs total, ${smoke_capture_jobs_per_run}/selected-run; profiling enabled"
  (( is_smoke )) && say "  smoke events: capture nEvents=${capture_nevents} per worker; replay nEvents=0 over captured pools"
  say "  DAG dir     : ${dag_dir}"
  if (( trace_pool )); then
    say "  trace       : enabled (RJ_POOL_TRACE_EVERY=${trace_every}; dryrun=${dag_dryrun}; smoke=${is_smoke})"
    say "  golden runs : ${golden_rows}"
    say "  capture axes: storedIsolation=[${data_capture_cones[*]}] clusterUE=[${uepipe_modes[*]}]"
    say "  replay axes : coneR=[${data_replay_cones[*]}] jetPt=[${data_pts[*]}] dphi=[${data_fracs[*]}] vz=[${data_vzs[*]}] ID/iso rows=${#iso_selection_tags[@]}"
  fi

  local smoke_capture_queued_total=0
  if [[ "$from_scratch" -eq 1 ]]; then
    for data_cone in "${data_capture_cones[@]}"; do
      for uepipe in "${uepipe_modes[@]}"; do
        local cap_tag
        cap_tag="$(pool_capture_cfg_tag "$data_cone" "$uepipe" "$master_yaml")"
        (( trace_pool )) && say "[poolSmoke] capture planning begin: cap_tag=${cap_tag} cone=${data_cone} clusterUE=${uepipe}"
        local cap_yaml
        cap_yaml="$(sim_make_yaml_override "$master_yaml" "${data_pts[0]}" "${data_fracs[0]}" "${data_vzs[0]}" "$data_cone" "${iso_sliding[0]}" "${iso_fixed[0]}" "$uepipe" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}" "$cap_tag" "$workflow_stamp" "false")"
        local cap_stage_dir="${STAGE_DIR}/${TAG}_${workflow_stamp}_capture_${cap_tag}"
        mkdir -p "$cap_stage_dir"
        local cap_root="${pool_base}/${cap_tag}"
        case "$cap_root" in
          */thesisAnaPools/*|*/thesisAnaPoolsSmoke/*) ;;
          *) err "Refusing to wipe pool capture root outside thesisAnaPools: ${cap_root}"; exit 62 ;;
        esac
        mkdir -p "$cap_root"
        (( trace_pool )) && say "[poolSmoke] cleaning capture pool root: ${cap_root}"
        find "$cap_root" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true
        (( trace_pool )) && say "[poolSmoke] capture pool root clean: ${cap_root}"

        local cap_sub="${dag_dir}/capture_${cap_tag}.sub"
        local cap_args="${dag_dir}/capture_${cap_tag}.args"
        : > "$cap_args"
        local exe_for_sub="${BULK_FROZEN_EXE:-${EXE}}"
        local macro_env_for_sub=""
        [[ -n "${BULK_FROZEN_MACRO:-}" ]] && macro_env_for_sub=";RJ_MACRO_PATH=${BULK_FROZEN_MACRO}"
        local cap_prefix="capture_${TAG}_${workflow_stamp}_${cap_tag}"
        cat > "$cap_sub" <<SUB
universe      = vanilla
executable    = ${exe_for_sub}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${cap_prefix}.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${cap_prefix}.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${cap_prefix}.\$(Cluster).\$(Process).err
request_memory= ${capture_request_memory}
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=0;RJ_CONFIG_YAML=${cap_yaml}${macro_env_for_sub};RJ_POOL_MODE=capture;RJ_PROFILE_JOB=${profile_job};RJ_POOL_PROFILE_COUNTS=1;RJ_JOB_HEARTBEAT_SECONDS=${job_heartbeat_seconds};RJ_PROFILE_STAGE=capture;RJ_PROFILE_LABEL=${cap_tag};RJ_REQUEST_MEMORY_MB=${capture_request_memory_mb}
${capture_submit_extra}
queue arguments from ${cap_args}
SUB

        local queued=0
        local runs_seen=0 runs_with_groups=0
        while IFS= read -r rn; do
          if (( is_smoke && smoke_capture_cap > 0 && smoke_capture_queued_total >= smoke_capture_cap )); then
            (( trace_pool )) && say "[poolSmoke] capture cap reached (${smoke_capture_queued_total}/${smoke_capture_cap}); stopping capture scan for ${cap_tag}"
            break
          fi
          [[ -z "$rn" || "$rn" =~ ^# ]] && continue
          (( runs_seen+=1 ))
          local r8
          r8="$(run8 "$rn")"
          if [[ -n "${TRIGGER_BIT}" ]] && ! is_trigger_active "$r8" "$TRIGGER_BIT"; then
            (( trace_pool && (runs_seen == 1 || runs_seen % trace_every == 0) )) && say "[poolSmoke] capture scan ${cap_tag}: run ${runs_seen}/${golden_rows} ${r8} skipped by trigger bit ${TRIGGER_BIT}"
            continue
          fi
          (( trace_pool && (runs_seen == 1 || runs_seen % trace_every == 0) )) && say "[poolSmoke] capture scan ${cap_tag}: run ${runs_seen}/${golden_rows} ${r8} grouping DST lists"
          mapfile -t groups < <( STAGE_DIR="$cap_stage_dir" make_groups "$r8" "$GROUP_SIZE" )
          (( ${#groups[@]} )) || continue
          (( runs_with_groups+=1 ))
          if (( is_smoke && smoke_capture_cap > 0 )); then
            local remaining=$(( smoke_capture_cap - smoke_capture_queued_total ))
            (( remaining > 0 )) || break
            if (( smoke_capture_jobs_per_run > 0 && ${#groups[@]} > smoke_capture_jobs_per_run )); then
              say "Smoke per-run cap: trimming capture groups for ${cap_tag}/${r8}: ${#groups[@]} → ${smoke_capture_jobs_per_run}"
              groups=( "${groups[@]:0:$smoke_capture_jobs_per_run}" )
            fi
            if (( ${#groups[@]} > remaining )); then
              say "Smoke cap: trimming capture groups for ${cap_tag}/${r8}: ${#groups[@]} → ${remaining}"
              groups=( "${groups[@]:0:$remaining}" )
            fi
          fi
          local expected_pool_list="${dag_dir}/pool_expected_${cap_tag}_${r8}.list"
          : > "$expected_pool_list"
          local gidx=0
          for glist in "${groups[@]}"; do
            (( gidx+=1 ))
            local chunk_base chunk_tag
            chunk_base="$(basename "$glist")"
            chunk_tag="${chunk_base%.list}"
            printf '%s\n' "${cap_root}/${r8}/RecoilJets_${analysis_tag}_${chunk_tag}.root" >> "$expected_pool_list"
            printf '%s %s %s $(Cluster) %s %d NONE %s\n' \
                   "$r8" "$glist" "$DATASET" "$capture_nevents" "$gidx" "$cap_root" >> "$cap_args"
            (( queued+=1 ))
            (( smoke_capture_queued_total+=1 ))
          done
          pool_list_by_tag_run["${cap_tag}|${r8}"]="$expected_pool_list"
          (( trace_pool )) && say "[poolSmoke] capture queue ${cap_tag}/${r8}: groups=${#groups[@]} queuedForTag=${queued} queuedTotal=${smoke_capture_queued_total}"
        done < "$workflow_golden"

        if (( queued == 0 )); then
          warn "No DATA pool capture jobs queued for ${cap_tag}; skipping capture node"
          rm -f "$cap_sub"
          continue
        fi

        local cap_node
        cap_node="CAP_$(sanitize_node_name "${cap_tag}")"
        printf 'JOB %s %s\n' "$cap_node" "$cap_sub" >> "$dag"
        cap_node_by_tag["$cap_tag"]="$cap_node"
        (( capture_count+=1 ))
        (( capture_job_total+=queued ))
        echo "capture_node=${cap_node} tag=${cap_tag} jobs=${queued} yaml=${cap_yaml} output=${cap_root}" >> "$manifest"
        echo "profile_glob=${OUT_DIR}/${cap_prefix}.*.out" >> "$manifest"
        (( trace_pool )) && say "[poolSmoke] capture node ready: node=${cap_node} jobs=${queued} runsSeen=${runs_seen} runsWithGroups=${runs_with_groups} sub=${cap_sub}"
      done
    done
  fi

  local cell_num=0
  for data_cone in "${data_capture_cones[@]}"; do
  for uepipe in "${uepipe_modes[@]}"; do
      local cap_tag fanout_dirs replay_sub replay_node fanout_count
      cap_tag="$(pool_capture_cfg_tag "$data_cone" "$uepipe" "$master_yaml")"
      fanout_dirs="${SIM_YAML_OVERRIDE_DIR}/pool_replay_fanout_${TAG}_${cap_tag}_${workflow_stamp}.txt"
      (( trace_pool )) && say "[poolSmoke] replay fanout begin: cap_tag=${cap_tag} outputBase=${out_base}"
      emit_pool_replay_fanout_dirs_file "$fanout_dirs" "$out_base" "$master_yaml" "$workflow_stamp" "${data_replay_cones[*]}" "$uepipe" "${data_pts[*]}" "${data_fracs[*]}" "${data_vzs[*]}"
      fanout_count="$(wc -l < "$fanout_dirs" | awk '{print $1}')"
      (( trace_pool )) && say "[poolSmoke] replay fanout written: cap_tag=${cap_tag} fanoutOutputs=${fanout_count} file=${fanout_dirs}"
      append_fanout_cfgs_to_manifest "$fanout_dirs" "$manifest"

      (( trace_pool )) && say "[poolSmoke] cleaning fanout output dirs begin: cap_tag=${cap_tag} rows=${fanout_count}"
      clean_fanout_output_dirs_from_file "$fanout_dirs" "$cap_tag" "$trace_pool" "$trace_every"
      (( trace_pool )) && say "[poolSmoke] cleaning fanout output dirs done: cap_tag=${cap_tag} uniqueRoots=${CLEAN_FANOUT_LAST_COUNT:-0}"

      local -a fanout_shards=()
      mapfile -t fanout_shards < <(split_pool_replay_fanout_shards "$fanout_dirs" "${dag_dir}/fanoutShards/${cap_tag}" "$cap_tag" "$replay_roots_per_shard" "$trace_pool")
      (( ${#fanout_shards[@]} )) || { err "No replay fanout shards were produced for ${cap_tag}"; exit 95; }

      local shard_i=0
      local fanout_shard shard_rows shard_roots
      for fanout_shard in "${fanout_shards[@]}"; do
      (( shard_i+=1 ))
      shard_rows="$(wc -l < "$fanout_shard" | awk '{print $1}')"
      shard_roots="$(awk -F'|' 'NF && $1 !~ /^#/ && $1 != "" {seen[$1]=1} END{for(k in seen)c++; print c+0}' "$fanout_shard")"
      replay_sub="${dag_dir}/replay_${cap_tag}_shard$(printf "%03d" "$shard_i").sub"
      local replay_args="${dag_dir}/replay_${cap_tag}_shard$(printf "%03d" "$shard_i").args"
      : > "$replay_args"
      local replay_prefix="replay_${TAG}_${workflow_stamp}_${cap_tag}_shard$(printf "%03d" "$shard_i")"
      cat > "$replay_sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${replay_prefix}.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${replay_prefix}.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${replay_prefix}.\$(Cluster).\$(Process).err
request_memory= ${replay_request_memory}
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=0;RJ_CONFIG_YAML=${master_yaml};RJ_MACRO_PATH=${pool_macro};RJ_ID_FANOUT_DIRS_FILE=${fanout_shard};RJ_REPLAY_MAX_OPEN_OUTPUTS=${replay_max_open_outputs};RJ_REPLAY_OUTPUT_ROOTS_PER_SHARD=${replay_roots_per_shard};RJ_PROFILE_JOB=${profile_job};RJ_POOL_PROFILE_COUNTS=1;RJ_JOB_HEARTBEAT_SECONDS=${job_heartbeat_seconds};RJ_PROFILE_STAGE=replay;RJ_PROFILE_LABEL=${cap_tag}_shard$(printf "%03d" "$shard_i");RJ_REQUEST_MEMORY_MB=${replay_request_memory_mb}
${replay_submit_extra}
queue arguments from ${replay_args}
SUB

      local queued=0
      local replay_runs_seen=0 replay_runs_with_pools=0
      while IFS= read -r rn; do
        [[ -z "$rn" || "$rn" =~ ^# ]] && continue
        (( replay_runs_seen+=1 ))
        local r8
        r8="$(run8 "$rn")"
        if [[ -n "${TRIGGER_BIT}" ]] && ! is_trigger_active "$r8" "$TRIGGER_BIT"; then
          (( trace_pool && (replay_runs_seen == 1 || replay_runs_seen % trace_every == 0) )) && say "[poolSmoke] replay scan ${cap_tag}: run ${replay_runs_seen}/${golden_rows} ${r8} skipped by trigger bit ${TRIGGER_BIT}"
          continue
        fi

        local pool_all
        if [[ "$from_scratch" -eq 1 ]]; then
          pool_all="${pool_list_by_tag_run["${cap_tag}|${r8}"]:-}"
          [[ -n "$pool_all" && -s "$pool_all" ]] || continue
        else
          local pool_dir="${pool_base}/${cap_tag}/${r8}"
          [[ -d "$pool_dir" ]] || { err "Missing pool directory for ${cap_tag}/${r8}: $pool_dir"; exit 93; }
          pool_all="${dag_dir}/pool_existing_${cap_tag}_${r8}.list"
          find "$pool_dir" -maxdepth 1 -type f -name '*.root' | sort > "$pool_all"
          [[ -s "$pool_all" ]] || { err "No pool ROOT files found in $pool_dir"; exit 94; }
        fi
        (( replay_runs_with_pools+=1 ))
        (( trace_pool )) && say "[poolSmoke] replay queue source ${cap_tag}/${r8}: poolList=${pool_all} replayGroupSize=${replay_gs}"

        local pool_stage_dir="${dag_dir}/poolReplayLists/${cap_tag}/${r8}"
        mkdir -p "$pool_stage_dir"
        rm -f "${pool_stage_dir}/pool_grp"*.list "${pool_stage_dir}/pool_grp_raw_"* 2>/dev/null || true
        split -l "$replay_gs" -d -a 5 "$pool_all" "${pool_stage_dir}/pool_grp_raw_"
        local gidx=0
        for raw in "${pool_stage_dir}/pool_grp_raw_"*; do
          [[ -s "$raw" ]] || { rm -f "$raw"; continue; }
          (( gidx+=1 ))
          local out_list="${pool_stage_dir}/pool_grp$(printf "%03d" "$gidx").list"
          mv "$raw" "$out_list"
          printf '%s %s %s $(Cluster) 0 %d NONE %s\n' \
                 "$r8" "$out_list" "$DATASET" "$gidx" "${out_base}/${cap_tag}" >> "$replay_args"
          (( queued+=1 ))
        done
        (( trace_pool )) && say "[poolSmoke] replay queue ${cap_tag}/${r8}: replayChunks=${gidx} queuedForTag=${queued}"
      done < "$workflow_golden"

      if (( queued == 0 )); then
        warn "No DATA pool replay jobs queued for ${cap_tag} shard ${shard_i}; skipping replay node"
        rm -f "$replay_sub"
        continue
      fi

      replay_node="REP_$(sanitize_node_name "${cap_tag}_shard$(printf "%03d" "$shard_i")")"
      if [[ "$from_scratch" -eq 1 ]]; then
        local cap_node="${cap_node_by_tag[$cap_tag]:-}"
        if [[ -z "$cap_node" && "$is_smoke" -eq 1 ]]; then
          warn "Smoke cap did not include capture jobs for ${cap_tag}; skipping replay shard ${shard_i}"
          rm -f "$replay_sub"
          continue
        fi
        [[ -n "$cap_node" ]] || { err "Internal DAG build error: missing capture node for ${cap_tag}"; exit 97; }
        printf 'JOB %s %s\n' "$replay_node" "$replay_sub" >> "$dag"
        printf 'PARENT %s CHILD %s\n' "$cap_node" "$replay_node" >> "$dag"
      else
        printf 'JOB %s %s\n' "$replay_node" "$replay_sub" >> "$dag"
      fi
      (( replay_count+=1 ))
      (( replay_job_total+=queued ))
      (( cell_num+=shard_rows ))
      echo "replay_node=${replay_node} cap_tag=${cap_tag} shard=${shard_i}/${#fanout_shards[@]} jobs=${queued} fanout_outputs=${shard_rows} fanout_roots=${shard_roots} fanout=${fanout_shard} parent_fanout=${fanout_dirs}" >> "$manifest"
      echo "profile_glob=${OUT_DIR}/${replay_prefix}.*.out" >> "$manifest"
      (( trace_pool )) && say "[poolSmoke] replay node ready: node=${replay_node} jobs=${queued} runsSeen=${replay_runs_seen} runsWithPools=${replay_runs_with_pools} shard=${shard_i}/${#fanout_shards[@]} fanoutRows=${shard_rows} fanoutRoots=${shard_roots} sub=${replay_sub}"
      done
  done
  done

  if (( is_smoke )) && (( capture_count == 0 || replay_count == 0 )); then
    err "poolSmoke produced capture_count=${capture_count}, replay_count=${replay_count}; refusing to submit an empty smoke DAG"
    exit 98
  fi

  local final_sub workflow_name
  if (( is_smoke )); then
    workflow_name="dataPoolSmoke_${TAG}_${workflow_stamp}"
  else
    workflow_name="$([[ "$from_scratch" -eq 1 ]] && echo "dataPoolFromScratch_${TAG}_${workflow_stamp}" || echo "dataPoolReplay_${TAG}_${workflow_stamp}")"
  fi
  final_sub="$(write_dag_final_summary_files "$dag_dir" "$workflow_name" "$DATASET" "$pool_base" "$out_base" "$manifest" "$dag")"
  printf 'FINAL FINAL_SUMMARY %s\n' "$final_sub" >> "$dag"

  say "DATA pool DAG build summary:"
  say "  capture nodes : ${capture_count}"
  say "  replay nodes  : ${replay_count}"
  say "  capture jobs  : ${capture_job_total}"
  say "  replay jobs   : ${replay_job_total}"
  say "  manifest      : ${manifest}"
  say "  dag           : ${dag}"
  {
    echo "capture_worker_jobs=${capture_job_total}"
    echo "replay_worker_jobs=${replay_job_total}"
    echo "total_worker_jobs=$(( capture_job_total + replay_job_total ))"
    echo "dag_max_worker_jobs=${RJ_DAG_MAX_WORKER_JOBS:-50000}"
    echo "dag_warn_worker_jobs=${RJ_DAG_WARN_WORKER_JOBS:-40000}"
  } >> "$manifest"
  local data_min_worker_jobs=0
  if (( ! is_smoke )); then
    data_min_worker_jobs="${RJ_DAG_MIN_WORKER_JOBS:-${RJ_DAG_MIN_WORKER_JOBS_DATA:-25000}}"
  fi
  echo "dag_min_worker_jobs=${data_min_worker_jobs}" >> "$manifest"
  check_dag_worker_job_budget "DATA pool workflow ${TAG}" "$capture_job_total" "$replay_job_total" "$dag_dryrun" "$data_min_worker_jobs"
  (( trace_pool )) && say "[poolSmoke] publishing latest manifest marker under ${out_base}"
  publish_production_manifest "$manifest" "$out_base"
  (( trace_pool )) && say "[poolSmoke] DAG build complete; submit step next (dryrun=${dag_dryrun})"
  if (( dag_dryrun )) && [[ "$workflow_flavor" == "smokeTest" ]]; then
    echo "RECOILJETS_SMOKETEST_DRYRUN_V1"
    echo "dataset=${DATASET}"
    echo "mode=condor smokeTest"
    echo "runs_or_samples=${golden_rows}"
    echo "manifest=${manifest}"
    echo "dag=${dag}"
  fi
  submit_dag_with_notify "$dag"
}

# ------------------------ AuAu embedded BDT training helpers --------------------------
auau_bdt_words() {
  local text="$1"
  local word
  for word in $text; do
    [[ -n "$word" ]] && printf "%s\n" "$word"
  done
}

auau_bdt_signal_samples() {
  auau_bdt_words "${RJ_AUAU_BDT_SIGNAL_SAMPLES:-$AUAU_BDT_SIGNAL_SAMPLES_DEFAULT}"
}

auau_bdt_background_samples() {
  auau_bdt_words "${RJ_AUAU_BDT_BACKGROUND_SAMPLES:-$AUAU_BDT_BACKGROUND_SAMPLES_DEFAULT}"
}

auau_ml_default_nfiles_for_events() {
  local events="$1"
  local events_per_file="${RJ_SIM_LOCAL_EVENTS_PER_FILE:-1000}"
  [[ "$events_per_file" =~ ^[0-9]+$ && "$events_per_file" -ge 1 ]] || events_per_file=1000
  [[ "$events" =~ ^[0-9]+$ && "$events" -ge 1 ]] || events=1
  echo $(( (events + events_per_file - 1) / events_per_file ))
}

auau_bdt_parse_local_controls() {
  AUAU_BDT_LOCAL_EVENTS="${LOCAL_EVENTS}"
  AUAU_BDT_LOCAL_VERBOSITY="10"
  AUAU_BDT_LOCAL_NFILES="${RJ_AUAU_BDT_LOCAL_NFILES:-}"
  local t
  for t in "${tokens[@]}"; do
    if [[ "$t" =~ ^VERBOSE=([0-9]+)$ ]]; then
      AUAU_BDT_LOCAL_VERBOSITY="${BASH_REMATCH[1]}"
    elif [[ "$t" =~ ^NFILES=([0-9]+)$ ]]; then
      AUAU_BDT_LOCAL_NFILES="${BASH_REMATCH[1]}"
    elif [[ "$t" =~ ^FILES=([0-9]+)$ ]]; then
      AUAU_BDT_LOCAL_NFILES="${BASH_REMATCH[1]}"
    elif [[ "$t" =~ ^[0-9]+$ ]]; then
      AUAU_BDT_LOCAL_EVENTS="$t"
    fi
  done
  [[ "$AUAU_BDT_LOCAL_EVENTS" =~ ^[0-9]+$ ]] || { err "BDT local event count must be an integer"; exit 2; }
  if [[ -z "$AUAU_BDT_LOCAL_NFILES" ]]; then
    AUAU_BDT_LOCAL_NFILES="$(auau_ml_default_nfiles_for_events "$AUAU_BDT_LOCAL_EVENTS")"
  fi
  [[ "$AUAU_BDT_LOCAL_NFILES" =~ ^[0-9]+$ && "$AUAU_BDT_LOCAL_NFILES" -ge 1 ]] || { err "BDT local NFILES must be a positive integer"; exit 2; }
}

auau_bdt_make_training_yaml() {
  local task="$1"
  local stamp="$2"
  local master_yaml="$3"

  [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

  local -a pts fracs vzs cones fixeds uepipes
  mapfile -t pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
  mapfile -t fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
  mapfile -t vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
  mapfile -t cones < <( yaml_get_values "coneR" "$master_yaml" )
  mapfile -t fixeds < <( yaml_get_values "fixedGeV" "$master_yaml" 2>/dev/null )
  (( ${#pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
  (( ${#fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
  (( ${#vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
  (( ${#cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
  (( ${#fixeds[@]} )) || fixeds=( "2.0" )

  read_uepipe_modes "$master_yaml" "simembedded"
  uepipes=( "${uepipe_modes[@]}" )
  (( ${#uepipes[@]} )) || uepipes=( "noSub" )

  local out
  out="$(sim_make_yaml_override \
    "$master_yaml" \
    "${pts[0]}" "${fracs[0]}" "${vzs[0]}" "${cones[0]}" \
    "false" "${fixeds[0]}" "${uepipes[0]}" \
    "reference" "reference" "reference" \
    "auauBDTTraining_${task}" "$stamp" "false")"

  local tmp="${out}.tmp"
  local npb_data_tagging="false"
  [[ "$task" == "npbData" ]] && npb_data_tagging="true"
  sed -E \
    -e "s|^([[:space:]]*auau_bdt_training_tree:).*|\\1 true|" \
    -e "s|^([[:space:]]*auau_bdt_training_tree_max_entries:).*|\\1 ${RJ_AUAU_BDT_TRAINING_TREE_MAX_ENTRIES:-0}|" \
    -e "s|^([[:space:]]*auau_bdt_npb_data_tagging:).*|\\1 ${npb_data_tagging}|" \
    "$out" > "$tmp"
  mv -f "$tmp" "$out"
  echo "$out"
}

auau_bdt_cfg_tag_from_yaml() {
  local yaml="$1"
  local -a pts fracs vzs cones
  mapfile -t pts   < <( yaml_get_values "jet_pt_min" "$yaml" )
  mapfile -t fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$yaml" )
  mapfile -t vzs   < <( yaml_get_values "vz_cut_cm" "$yaml" )
  mapfile -t cones < <( yaml_get_values "coneR" "$yaml" )
  build_iso_modes "$yaml"
  read_uepipe_modes "$yaml" "simembedded"
  matrix_cfg_tag "${pts[0]}" "${fracs[0]}" "${vzs[0]}" "${cones[0]}" 0 "${uepipe_modes[0]}"
}

auau_bdt_collect_roots() {
  local root_dir="$1"
  local manifest="$2"
  mkdir -p "$(dirname "$manifest")"
  find "$root_dir" \
    -type f \
    -name "RecoilJets_*.root" \
    ! -path "*/_bad_root_files_*/*" \
    ! -path "*/bad_root_quarantine/*" \
    | sort > "$manifest" || true
  [[ -s "$manifest" ]] || { err "No RecoilJets ROOT files found under ${root_dir}"; exit 31; }
}

auau_ml_validate_manifest_tree() {
  local manifest="$1"
  local tree="$2"
  local label="$3"
  local report="${4:-}"
  [[ -s "$manifest" ]] || { err "Cannot validate empty manifest: $manifest"; exit 31; }
  if [[ -n "$report" ]]; then
    mkdir -p "$(dirname "$report")"
  fi
  auau_ml_run_python - "$manifest" "$tree" "$label" "$report" <<'PY'
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
tree_name = sys.argv[2]
label = sys.argv[3]
report = Path(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else None

try:
    import uproot
except Exception as exc:
    raise SystemExit(f"uproot is required to validate {label}: {exc}") from exc

paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip() and not line.startswith("#")]
missing = []
entries = []
for path in paths:
    try:
        with uproot.open(path) as root_file:
            if tree_name not in root_file:
                missing.append(str(path))
                continue
            entries.append((str(path), int(root_file[tree_name].num_entries)))
    except Exception as exc:  # noqa: BLE001
        missing.append(f"{path} ({exc})")

total_entries = sum(count for _, count in entries)
empty = [path for path, count in entries if count == 0]
lines = [
    f"label={label}",
    f"tree={tree_name}",
    f"files={len(paths)}",
    f"files_with_tree={len(entries)}",
    f"files_missing_tree={len(missing)}",
    f"files_empty_tree={len(empty)}",
    f"entries_total={total_entries}",
]
if missing:
    lines.append("missing_tree_preview:")
    lines.extend(f"  {item}" for item in missing[:20])
if empty:
    lines.append("empty_tree_preview:")
    lines.extend(f"  {item}" for item in empty[:20])
text = "\n".join(lines) + "\n"
print(text, end="")
if report:
    report.write_text(text)
if missing:
    raise SystemExit(f"{label}: {len(missing)} file(s) are missing {tree_name}; refusing to train on stale/incomplete extraction.")
PY
}

auau_ml_python_label() {
  echo "${RJ_ML_PYTHON:-python3}"
}

auau_ml_run_python() {
  local py="${RJ_ML_PYTHON:-python3}"
  local -a pycmd
  read -r -a pycmd <<< "$py"
  "${pycmd[@]}" "$@"
}

auau_ml_check_python_deps() {
  if [[ "${RJ_ML_SKIP_PYTHON_PREFLIGHT:-0}" == "1" ]]; then
    warn "Skipping ML Python dependency preflight because RJ_ML_SKIP_PYTHON_PREFLIGHT=1"
    return 0
  fi

  say "Checking ML Python environment: $(auau_ml_python_label)"
  if auau_ml_run_python - <<'PY'
import importlib
import sys

required = ["uproot", "pandas", "numpy", "sklearn", "xgboost", "ROOT"]
missing = []
for module in required:
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append((module, str(exc)))

if missing:
    print("[ERROR] Missing ML Python dependencies:", file=sys.stderr)
    for module, exc in missing:
        print(f"  {module}: {exc}", file=sys.stderr)
    sys.exit(17)

print("[OK] ML Python dependencies available: " + ", ".join(required))
PY
  then
    return 0
  fi

  err "ML Python dependency preflight failed."
  err "Set RJ_ML_PYTHON to a Python with uproot, pandas, numpy, scikit-learn, xgboost, and PyROOT."
  err "Example: RJ_ML_PYTHON=/path/to/python ./RecoilJets_Condor_submit.sh isSimEmbeddedAndInclusive trainMLAll resume /path/to/mlIntegration_run"
  return 17
}

auau_bdt_train_from_manifest() {
  local task="$1"
  local manifest="$2"
  local outdir="$3"
  mkdir -p "$outdir"
  if [[ "$task" == "npb" && -z "${RJ_AUAU_BDT_NPB_INPUTS:-}" ]]; then
    warn "NPB extraction smoke test finished, but no real NPB-labeled data inputs were supplied."
    warn "PPG12-style NPB training needs npb_label=0 rows from timing/streak-tagged data clusters."
    warn "Set RJ_AUAU_BDT_NPB_INPUTS to ROOT files with AuAuPhotonIDTrainingTree/npb_label, then rerun."
    say "Physics-cluster manifest from embedded sim: ${manifest}"
    return 0
  fi

  local -a train_args
  train_args=( "--task" "$task" "--input" "@${manifest}" "--outdir" "$outdir" )
  if [[ "$task" == "npb" && -n "${RJ_AUAU_BDT_NPB_INPUTS:-}" ]]; then
    local npb_manifest="${outdir}/npb_labeled_inputs.list"
    auau_bdt_words "${RJ_AUAU_BDT_NPB_INPUTS}" > "$npb_manifest"
    train_args=( "--task" "$task" "--input" "@${manifest}" "@${npb_manifest}" "--missing-label-value" "1" "--outdir" "$outdir" )
  fi

  say "Training AuAu ${task} BDT:"
  say "  manifest : ${manifest}"
  say "  outdir   : ${outdir}"
  if [[ "$task" == "tight" && "${RJ_AUAU_BDT_LOCAL_WIDE:-1}" == "1" ]]; then
    local cent_bins="${RJ_AUAU_TIGHT_BDT_CENT_BINS:-0:20,20:50,50:80}"
    say "  product: centINDcontrol (PPG12 base_v3E feature order, no centrality feature, all centralities)"
    auau_ml_run_python "${BASE}/scripts/train_auau_photon_bdt.py" "${train_args[@]}" \
      --tight-mode centINDcontrol \
      --prefix "auau_tight_bdt_centINDcontrol" \
      --cent-bins "$cent_bins" \
      --test-size 0.2 \
      --random-seed 42 \
      --n-estimators 750 \
      --max-depth 5 \
      --learning-rate 0.1 \
      --subsample 0.5 \
      --colsample-bytree 0.6 \
      --reg-alpha 5.0 \
      --reg-lambda 0.3 \
      --grow-policy lossguide \
      --max-bin 256 \
      --background-subsample-fraction 0.3 \
      --background-subsample-et-threshold 15 \
      --background-subsample-bins 20 \
      --background-subsample-seed 42 \
      --background-subsample-flatten

    say "  product: centAsFeat (same model family with centrality appended as a feature)"
    auau_ml_run_python "${BASE}/scripts/train_auau_photon_bdt.py" "${train_args[@]}" \
      --tight-mode centAsFeat \
      --prefix "auau_tight_bdt_centAsFeat" \
      --cent-bins "$cent_bins" \
      --test-size 0.2 \
      --random-seed 42 \
      --n-estimators 750 \
      --max-depth 5 \
      --learning-rate 0.1 \
      --subsample 0.5 \
      --colsample-bytree 0.6 \
      --reg-alpha 5.0 \
      --reg-lambda 0.3 \
      --grow-policy lossguide \
      --max-bin 256 \
      --background-subsample-fraction 0.3 \
      --background-subsample-et-threshold 15 \
      --background-subsample-bins 20 \
      --background-subsample-seed 42 \
      --background-subsample-flatten

    say "  product: centDepBDTs (one PPG12-like BDT per configured centrality bin: ${cent_bins})"
    auau_ml_run_python "${BASE}/scripts/train_auau_photon_bdt.py" "${train_args[@]}" \
      --tight-mode centDepBDTs \
      --prefix "auau_tight_bdt_centDepBDTs" \
      --cent-bins "$cent_bins" \
      --test-size 0.2 \
      --random-seed 42 \
      --n-estimators 750 \
      --max-depth 5 \
      --learning-rate 0.1 \
      --subsample 0.5 \
      --colsample-bytree 0.6 \
      --reg-alpha 5.0 \
      --reg-lambda 0.3 \
      --grow-policy lossguide \
      --max-bin 256 \
      --background-subsample-fraction 0.3 \
      --background-subsample-et-threshold 15 \
      --background-subsample-bins 20 \
      --background-subsample-seed 42 \
      --background-subsample-flatten
  else
    auau_ml_run_python "${BASE}/scripts/train_auau_photon_bdt.py" "${train_args[@]}"
  fi
}

auau_bdt_csv_to_yaml_list() {
  local csv="$1"
  local out="["
  local IFS=','
  local item first=1
  for item in $csv; do
    [[ -n "$item" ]] || continue
    if [[ "$first" -eq 0 ]]; then out+=", "; fi
    out+="$item"
    first=0
  done
  out+="]"
  printf '%s\n' "$out"
}

auau_tight_bdt_centdep_model_list() {
  local model_dir="$1"
  local cent_bins="${RJ_AUAU_TIGHT_BDT_CENT_BINS:-0:20,20:50,50:80}"
  local csv="" item lo hi tag model
  local IFS=','
  for item in $cent_bins; do
    lo="${item%%:*}"
    hi="${item##*:}"
    printf -v tag "cent_%03d_%03d" "$lo" "$hi"
    model="${model_dir}/auau_tight_bdt_centDepBDTs_${tag}_tmva.root"
    [[ -s "$model" ]] || { err "Missing centrality-dependent tight BDT model: $model"; exit 33; }
    csv+="${csv:+,}${model}"
  done
  auau_bdt_csv_to_yaml_list "$csv"
}

auau_tight_bdt_apply_product() {
  local product="$1" model_dir="$2" master_yaml="$3" stamp="$4" output_base="$5"
  shift 5
  local -a samples=( "$@" )
  local yaml model

  yaml="$(auau_ml_first_matrix_yaml "$master_yaml" "tightBDTApply_${product}_${stamp}" "$stamp" "reference" "$product" "$product")"
  ml_yaml_set_scalar "$yaml" "auau_bdt_training_tree" "false"
  ml_yaml_set_scalar "$yaml" "jet_ml_training_tree" "false"
  ml_yaml_set_scalar "$yaml" "jet_ml_correction_enabled" "false"

  case "$product" in
    centINDcontrol)
      model="${model_dir}/auau_tight_bdt_centINDcontrol_allCent_tmva.root"
      [[ -s "$model" ]] || { err "Missing centINDcontrol tight BDT model: $model"; exit 33; }
      ml_yaml_set_scalar "$yaml" "auau_tight_bdt_centINDcontrol_model_file" "$model"
      ml_yaml_set_scalar "$yaml" "auau_tight_bdt_model_file" "$model"
      ;;
    centAsFeat)
      model="${model_dir}/auau_tight_bdt_centAsFeat_allCent_tmva.root"
      [[ -s "$model" ]] || { err "Missing centAsFeat tight BDT model: $model"; exit 33; }
      ml_yaml_set_scalar "$yaml" "auau_tight_bdt_centAsFeat_model_file" "$model"
      ml_yaml_set_scalar "$yaml" "auau_tight_bdt_model_file" "$model"
      ;;
    centDepBDTs)
      ml_yaml_set_scalar "$yaml" "auau_tight_bdt_centDep_model_files" "$(auau_tight_bdt_centdep_model_list "$model_dir")"
      ;;
    *)
      err "Unknown tight BDT apply product: $product"
      exit 2
      ;;
  esac

  say "Re-applying tight BDT product ${product}"
  say "  YAML : ${yaml}"
  local samp
  for samp in "${samples[@]}"; do
    RJ_CONFIG_YAML="$yaml" \
    RJ_SIM_LOCAL_GROUPED=1 \
    RJ_SIM_LOCAL_NFILES="$AUAU_BDT_LOCAL_NFILES" \
    RJ_LOCAL_SIM_OUTPUT_BASE="${output_base}/${product}" \
    "$0" isSimEmbedded local "$AUAU_BDT_LOCAL_EVENTS" "VERBOSE=${AUAU_BDT_LOCAL_VERBOSITY}" "SAMPLE=${samp}"
  done
}

auau_bdt_run_local() {
  local task="$1"
  auau_bdt_parse_local_controls
  if [[ "$task" == "tight" || -n "${RJ_AUAU_BDT_NPB_INPUTS:-}" || "${RJ_AUAU_BDT_NPB_DATA_LOCAL:-0}" == "1" ]]; then
    auau_ml_check_python_deps
  fi

  local stamp master_yaml train_yaml cfg_tag local_root manifest outdir
  local npb_data_manifest=""
  stamp="$(date +%Y%m%d_%H%M%S)"
  master_yaml="$(sim_yaml_master_path)"
  train_yaml="$(auau_bdt_make_training_yaml "$task" "$stamp" "$master_yaml")"
  cfg_tag="$(auau_bdt_cfg_tag_from_yaml "$train_yaml")"
  local_root="${AUAU_BDT_LOCAL_BASE}/${task}_${stamp}"
  manifest="${local_root}/training_roots.list"
  outdir="${AUAU_BDT_MODEL_BASE}/${task}_${stamp}"

  local -a signal_samples background_samples
  mapfile -t signal_samples < <( auau_bdt_signal_samples )
  mapfile -t background_samples < <( auau_bdt_background_samples )

  say "AuAu embedded ${task} BDT LOCAL extraction"
  say "  YAML override      : ${train_yaml}"
  say "  cfg tag            : ${cfg_tag}"
  say "  local output root  : ${local_root}"
  say "  signal samples     : ${signal_samples[*]}"
  say "  background samples : ${background_samples[*]}"
  say "  events/sample      : ${AUAU_BDT_LOCAL_EVENTS}"
  say "  grouped files/sample: ${AUAU_BDT_LOCAL_NFILES} (assuming ${RJ_SIM_LOCAL_EVENTS_PER_FILE:-1000} events/input file unless overridden)"
  say "  verbosity          : ${AUAU_BDT_LOCAL_VERBOSITY}"
  echo

  local samp
  for samp in "${signal_samples[@]}"; do
    say "[BDT local] signal extraction sample=${samp}"
    RJ_CONFIG_YAML="$train_yaml" \
    RJ_SIM_LOCAL_GROUPED=1 \
    RJ_SIM_LOCAL_NFILES="$AUAU_BDT_LOCAL_NFILES" \
    RJ_LOCAL_SIM_OUTPUT_BASE="${local_root}/simembedded" \
    "$0" isSimEmbedded local "$AUAU_BDT_LOCAL_EVENTS" "VERBOSE=${AUAU_BDT_LOCAL_VERBOSITY}" "SAMPLE=${samp}"
  done
  for samp in "${background_samples[@]}"; do
    say "[BDT local] background extraction sample=${samp}"
    RJ_CONFIG_YAML="$train_yaml" \
    RJ_SIM_LOCAL_GROUPED=1 \
    RJ_SIM_LOCAL_NFILES="$AUAU_BDT_LOCAL_NFILES" \
    RJ_LOCAL_SIM_OUTPUT_BASE="${local_root}/simembeddedinclusive" \
    "$0" isSimEmbeddedInclusive local "$AUAU_BDT_LOCAL_EVENTS" "VERBOSE=${AUAU_BDT_LOCAL_VERBOSITY}" "SAMPLE=${samp}"
  done

  auau_bdt_collect_roots "$local_root" "$manifest"
  auau_ml_validate_manifest_tree "$manifest" "AuAuPhotonIDTrainingTree" "${task} photon rows" "${local_root}/qa_photon_tree_validation.txt"

  if [[ "$task" == "npb" && "${RJ_AUAU_BDT_NPB_DATA_LOCAL:-0}" == "1" ]]; then
    local data_yaml data_root
    data_yaml="$(auau_bdt_make_training_yaml "npbData" "$stamp" "$master_yaml")"
    data_root="${local_root}/auau_data_npb"
    say "AuAu data NPB LOCAL extraction (PPG12 timing-tagged NPB label)"
    say "  YAML override : ${data_yaml}"
    say "  output root   : ${data_root}"
    say "  label policy  : npb_label=0 for tagged data NPB rows"
    RJ_CONFIG_YAML="$data_yaml" \
    RJ_LOCAL_DATA_OUTPUT_BASE="$data_root" \
    "$0" isAuAu local "$AUAU_BDT_LOCAL_EVENTS" "VERBOSE=${AUAU_BDT_LOCAL_VERBOSITY}"
    npb_data_manifest="${local_root}/npb_data_training_roots.list"
    auau_bdt_collect_roots "$data_root" "$npb_data_manifest"
    local npb_data_report="${local_root}/qa_npb_data_tree_validation.txt"
    auau_ml_validate_manifest_tree "$npb_data_manifest" "AuAuPhotonIDTrainingTree" "AuAu data NPB rows" "$npb_data_report"
    local npb_data_entries
    npb_data_entries="$(awk -F= '/^entries_total=/{print $2}' "$npb_data_report" 2>/dev/null | tail -n 1)"
    if [[ "${npb_data_entries:-0}" == "0" ]]; then
      warn "AuAu data NPB extraction produced zero tagged rows. Skipping NPB training for this local smoke pass."
      npb_data_manifest=""
    fi
  fi

  if [[ -n "$npb_data_manifest" ]]; then
    RJ_AUAU_BDT_NPB_INPUTS="@${npb_data_manifest}" auau_bdt_train_from_manifest "$task" "$manifest" "$outdir"
  else
    auau_bdt_train_from_manifest "$task" "$manifest" "$outdir"
  fi

  if [[ "$task" == "tight" && "${RJ_AUAU_TIGHT_BDT_APPLY_LOCAL:-1}" == "1" ]]; then
    mkdir -p "${local_root}/apply"
    say "AuAu tight BDT LOCAL re-apply sanity"
    auau_tight_bdt_apply_product "centINDcontrol" "$outdir" "$master_yaml" "$stamp" "${local_root}/apply" "${signal_samples[@]}"
    auau_tight_bdt_apply_product "centAsFeat" "$outdir" "$master_yaml" "$stamp" "${local_root}/apply" "${signal_samples[@]}"
    auau_tight_bdt_apply_product "centDepBDTs" "$outdir" "$master_yaml" "$stamp" "${local_root}/apply" "${signal_samples[@]}"
    {
      echo "run_base=${local_root}"
      echo "model_dir=${outdir}"
      echo "manifest=${manifest}"
      echo "cent_bins=${RJ_AUAU_TIGHT_BDT_CENT_BINS:-0:20,20:50,50:80}"
      echo "products=centINDcontrol centAsFeat centDepBDTs"
    } > "${local_root}/tight_bdt_local_summary.env"
    say "Tight BDT local pass complete."
    say "  run base : ${local_root}"
    say "  models   : ${outdir}"
    say "  pull from Mac with:"
    say "    ./scripts/sftp_get_recoiljets_outputs.sh tightBDTSmoke ${local_root}"
  fi
}

auau_bdt_run_condor_do_all() {
  local task="$1"
  local gs_doall="$GROUP_SIZE"
  if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
    gs_doall="7"
  fi

  local stamp master_yaml train_yaml cfg_tag outdir
  stamp="$(date +%Y%m%d_%H%M%S)"
  master_yaml="$(sim_yaml_master_path)"
  train_yaml="$(auau_bdt_make_training_yaml "$task" "$stamp" "$master_yaml")"
  cfg_tag="$(auau_bdt_cfg_tag_from_yaml "$train_yaml")"
  outdir="${AUAU_BDT_MODEL_BASE}/${task}_${stamp}"

  local -a signal_samples background_samples
  mapfile -t signal_samples < <( auau_bdt_signal_samples )
  mapfile -t background_samples < <( auau_bdt_background_samples )

  say "AuAu embedded ${task} BDT CONDOR extraction"
  say "  YAML override      : ${train_yaml}"
  say "  cfg tag            : ${cfg_tag}"
  say "  output base        : ${AUAU_BDT_DEST_BASE}"
  say "  groupSize          : ${gs_doall}"
  say "  signal samples     : ${signal_samples[*]}"
  say "  background samples : ${background_samples[*]}"
  echo

  local samp
  for samp in "${signal_samples[@]}"; do
    say "[BDT condor] signal extraction sample=${samp}"
    RJ_CONFIG_YAML="$train_yaml" \
    RJ_SIMEMBED_DEST_BASE="$AUAU_BDT_DEST_BASE" \
    "$0" isSimEmbedded condorDoAll groupSize "$gs_doall" maxJobs "$MAX_JOBS" "SAMPLE=${samp}"
  done
  for samp in "${background_samples[@]}"; do
    say "[BDT condor] background extraction sample=${samp}"
    RJ_CONFIG_YAML="$train_yaml" \
    RJ_SIMEMBED_DEST_BASE="$AUAU_BDT_DEST_BASE" \
    "$0" isSimEmbeddedInclusive condorDoAll groupSize "$gs_doall" maxJobs "$MAX_JOBS" "SAMPLE=${samp}"
  done

  say "Condor extraction was submitted. After jobs finish, train from the produced ROOT files with:"
  echo
  printf '  mkdir -p %q\n' "$outdir"
  printf '  find %q -path %q -name %q | sort > %q\n' "$AUAU_BDT_DEST_BASE" "*/${cfg_tag}/*" "RecoilJets_*.root" "${outdir}/training_roots.list"
  if [[ "$task" == "npb" ]]; then
    printf '  # Append real data-tagged NPB ROOT files with npb_label=0 rows before training.\n'
    printf '  # Embedded sim files missing npb_label are treated as physics-side npb_label=1 by --missing-label-value 1.\n'
    printf '  %q --task %q --input @%q @/path/to/npb_labeled_data_roots.list --missing-label-value 1 --outdir %q\n' "${BASE}/scripts/train_auau_photon_bdt.py" "$task" "${outdir}/training_roots.list" "$outdir"
  else
    printf '  %q --task %q --input @%q --outdir %q\n' "${BASE}/scripts/train_auau_photon_bdt.py" "$task" "${outdir}/training_roots.list" "$outdir"
  fi
}

auau_tight_bdt_smoke_second_pass() {
  auau_ml_check_python_deps
  local source_base="${RJ_AUAU_TIGHT_BDT_SMOKE_SOURCE:-}"
  local t expect_source=0
  for t in "${tokens[@]}"; do
    if [[ "$expect_source" -eq 1 ]]; then
      source_base="$t"
      expect_source=0
      continue
    fi
    if [[ "$t" == "SOURCE" || "$t" == "source" ]]; then
      expect_source=1
    elif [[ "$t" =~ ^SOURCE=(.+)$ ]]; then
      source_base="${BASH_REMATCH[1]}"
    elif [[ "$t" == /* || "$t" == ./* ]]; then
      source_base="$t"
    fi
  done
  [[ -n "$source_base" ]] || { err "smokeTestSecondPass needs SOURCE=/path/to/finished/extraction or a path argument"; exit 2; }
  [[ -d "$source_base" ]] || { err "smokeTestSecondPass source directory does not exist: $source_base"; exit 2; }

  local stamp outbase manifest model_dir
  stamp="$(date +%Y%m%d_%H%M%S)"
  outbase="${AUAU_BDT_LOCAL_BASE}/tightSmokeSecond_${stamp}"
  manifest="${outbase}/training_roots.list"
  model_dir="${AUAU_BDT_MODEL_BASE}/tightSmokeSecond_${stamp}"
  mkdir -p "$outbase" "$model_dir"

  auau_bdt_collect_roots "$source_base" "$manifest"
  auau_ml_validate_manifest_tree "$manifest" "AuAuPhotonIDTrainingTree" "tight smoke second-pass photon rows" "${outbase}/qa_photon_tree_validation.txt"
  RJ_AUAU_BDT_LOCAL_WIDE=1 auau_bdt_train_from_manifest "tight" "$manifest" "$model_dir"

  {
    echo "source_base=${source_base}"
    echo "run_base=${outbase}"
    echo "model_dir=${model_dir}"
    echo "manifest=${manifest}"
    echo "cent_bins=${RJ_AUAU_TIGHT_BDT_CENT_BINS:-0:20,20:50,50:80}"
  } > "${outbase}/tight_bdt_smoke_second_pass.env"

  say "Tight BDT smoke second pass complete."
  say "  source   : ${source_base}"
  say "  run base : ${outbase}"
  say "  models   : ${model_dir}"
  say "  pull from Mac with:"
  say "    ./scripts/sftp_get_recoiljets_outputs.sh tightBDTSmoke ${outbase}"
}

auau_tight_bdt_smoke_apply_existing() {
  auau_bdt_parse_local_controls

  local model_dir="${RJ_AUAU_TIGHT_BDT_MODEL_DIR:-${RJ_AUAU_TIGHT_BDT_APPLY_MODEL_DIR:-}}"
  local t expect_model_dir=0
  for t in "${tokens[@]}"; do
    if [[ "$expect_model_dir" -eq 1 ]]; then
      model_dir="$t"
      expect_model_dir=0
      continue
    fi
    if [[ "$t" == "MODEL_DIR" || "$t" == "MODELDIR" || "$t" == "model_dir" || "$t" == "modelDir" ]]; then
      expect_model_dir=1
    elif [[ "$t" =~ ^MODEL_DIR=(.+)$ ]]; then
      model_dir="${BASH_REMATCH[1]}"
    elif [[ "$t" =~ ^MODELDIR=(.+)$ ]]; then
      model_dir="${BASH_REMATCH[1]}"
    elif [[ "$t" =~ ^model_dir=(.+)$ ]]; then
      model_dir="${BASH_REMATCH[1]}"
    elif [[ "$t" =~ ^modelDir=(.+)$ ]]; then
      model_dir="${BASH_REMATCH[1]}"
    fi
  done

  [[ -n "$model_dir" ]] || { err "smokeTestApplyExisting needs MODEL_DIR=/path/to/tight/models"; exit 2; }
  [[ -d "$model_dir" ]] || { err "Tight-BDT model directory does not exist: $model_dir"; exit 2; }

  local stamp master_yaml outbase
  stamp="$(date +%Y%m%d_%H%M%S)"
  master_yaml="$(sim_yaml_master_path)"
  outbase="${AUAU_BDT_LOCAL_BASE}/tightSmokeApply_${stamp}"
  mkdir -p "$outbase"

  local -a signal_samples
  mapfile -t signal_samples < <( auau_bdt_signal_samples )

  say "AuAu tight BDT APPLY-ONLY smoke using existing models"
  say "  model dir          : ${model_dir}"
  say "  local output root  : ${outbase}"
  say "  signal samples     : ${signal_samples[*]}"
  say "  events/sample      : ${AUAU_BDT_LOCAL_EVENTS}"
  say "  grouped files/sample: ${AUAU_BDT_LOCAL_NFILES} (assuming ${RJ_SIM_LOCAL_EVENTS_PER_FILE:-1000} events/input file unless overridden)"
  say "  verbosity          : ${AUAU_BDT_LOCAL_VERBOSITY}"
  echo

  auau_tight_bdt_apply_product "centINDcontrol" "$model_dir" "$master_yaml" "$stamp" "$outbase" "${signal_samples[@]}"
  auau_tight_bdt_apply_product "centAsFeat" "$model_dir" "$master_yaml" "$stamp" "$outbase" "${signal_samples[@]}"
  auau_tight_bdt_apply_product "centDepBDTs" "$model_dir" "$master_yaml" "$stamp" "$outbase" "${signal_samples[@]}"

  {
    echo "run_base=${outbase}"
    echo "model_dir=${model_dir}"
    echo "cent_bins=${RJ_AUAU_TIGHT_BDT_CENT_BINS:-0:20,20:50,50:80}"
    echo "products=centINDcontrol centAsFeat centDepBDTs"
  } > "${outbase}/tight_bdt_apply_existing_summary.env"

  say "Tight BDT apply-only smoke complete."
  say "  run base : ${outbase}"
  say "  models   : ${model_dir}"
  say "  pull from Mac with:"
  say "    ./scripts/sftp_get_recoiljets_outputs.sh tightBDTSmoke ${outbase}"
}

auau_jetml_make_training_yaml() {
  local stamp="$1"
  local master_yaml="$2"
  [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

  local -a pts fracs vzs cones fixeds uepipes
  mapfile -t pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
  mapfile -t fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
  mapfile -t vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
  mapfile -t cones < <( yaml_get_values "coneR" "$master_yaml" )
  mapfile -t fixeds < <( yaml_get_values "fixedGeV" "$master_yaml" 2>/dev/null )
  (( ${#pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
  (( ${#fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
  (( ${#vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
  (( ${#cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
  (( ${#fixeds[@]} )) || fixeds=( "2.0" )

  read_uepipe_modes "$master_yaml" "simembedded"
  uepipes=( "${uepipe_modes[@]}" )
  (( ${#uepipes[@]} )) || uepipes=( "noSub" )

  local out tmp
  out="$(sim_make_yaml_override \
    "$master_yaml" \
    "${pts[0]}" "${fracs[0]}" "${vzs[0]}" "${cones[0]}" \
    "false" "${fixeds[0]}" "${uepipes[0]}" \
    "reference" "reference" "reference" \
    "jetMLTraining" "$stamp" "false")"

  tmp="${out}.tmp"
  sed -E \
    -e "s|^([[:space:]]*auau_bdt_training_tree:).*|\\1 false|" \
    -e "s|^([[:space:]]*jet_ml_training_tree:).*|\\1 true|" \
    -e "s|^([[:space:]]*jet_ml_training_tree_max_entries:).*|\\1 ${RJ_JET_ML_TRAINING_TREE_MAX_ENTRIES:-0}|" \
    -e "s|^([[:space:]]*jet_ml_correction_enabled:).*|\\1 false|" \
    "$out" > "$tmp"
  mv -f "$tmp" "$out"
  echo "$out"
}

auau_jetml_run_local() {
  auau_bdt_parse_local_controls
  auau_ml_check_python_deps

  local stamp master_yaml train_yaml cfg_tag local_root manifest outdir
  stamp="$(date +%Y%m%d_%H%M%S)"
  master_yaml="$(sim_yaml_master_path)"
  train_yaml="$(auau_jetml_make_training_yaml "$stamp" "$master_yaml")"
  cfg_tag="$(auau_bdt_cfg_tag_from_yaml "$train_yaml")"
  local_root="${AUAU_BDT_LOCAL_BASE}/jetResidual_${stamp}"
  manifest="${local_root}/jetml_training_roots.list"
  outdir="${AUAU_BDT_MODEL_BASE}/jetResidual_${stamp}"

  local -a signal_samples background_samples
  mapfile -t signal_samples < <( auau_bdt_signal_samples )
  mapfile -t background_samples < <( auau_bdt_background_samples )

  say "AuAu embedded JetML residual LOCAL extraction"
  say "  YAML override      : ${train_yaml}"
  say "  cfg tag            : ${cfg_tag}"
  say "  local output root  : ${local_root}"
  say "  photon+jet samples : ${signal_samples[*]}"
  say "  inclusive samples  : ${background_samples[*]}"
  say "  events/sample      : ${AUAU_BDT_LOCAL_EVENTS}"
  say "  grouped files/sample: ${AUAU_BDT_LOCAL_NFILES} (assuming ${RJ_SIM_LOCAL_EVENTS_PER_FILE:-1000} events/input file unless overridden)"
  say "  verbosity          : ${AUAU_BDT_LOCAL_VERBOSITY}"
  echo

  local samp
  for samp in "${signal_samples[@]}"; do
    say "[JetML local] photon+jet extraction sample=${samp}"
    RJ_CONFIG_YAML="$train_yaml" \
    RJ_SIM_LOCAL_GROUPED=1 \
    RJ_SIM_LOCAL_NFILES="$AUAU_BDT_LOCAL_NFILES" \
    RJ_LOCAL_SIM_OUTPUT_BASE="${local_root}/simembedded" \
    "$0" isSimEmbedded local "$AUAU_BDT_LOCAL_EVENTS" "VERBOSE=${AUAU_BDT_LOCAL_VERBOSITY}" "SAMPLE=${samp}"
  done
  for samp in "${background_samples[@]}"; do
    say "[JetML local] inclusive extraction sample=${samp}"
    RJ_CONFIG_YAML="$train_yaml" \
    RJ_SIM_LOCAL_GROUPED=1 \
    RJ_SIM_LOCAL_NFILES="$AUAU_BDT_LOCAL_NFILES" \
    RJ_LOCAL_SIM_OUTPUT_BASE="${local_root}/simembeddedinclusive" \
    "$0" isSimEmbeddedInclusive local "$AUAU_BDT_LOCAL_EVENTS" "VERBOSE=${AUAU_BDT_LOCAL_VERBOSITY}" "SAMPLE=${samp}"
  done

  auau_bdt_collect_roots "$local_root" "$manifest"
  auau_ml_validate_manifest_tree "$manifest" "JetResidualMLTrainingTree" "JetML residual rows" "${local_root}/qa_jetml_tree_validation.txt"
  mkdir -p "$outdir"
  say "Training JetML residual model:"
  say "  manifest : ${manifest}"
  say "  outdir   : ${outdir}"
  if [[ "${RJ_JET_ML_LOCAL_WIDE:-1}" == "1" ]]; then
    say "  wide scan: core low-bias feature set"
    auau_ml_run_python "${BASE}/scripts/train_auau_jet_residual_bdt.py" \
      --input "@${manifest}" \
      --outdir "${outdir}/core" \
      --model-name "jetResidualBDT_core_allCent_tmva.root" \
      --features "reco_areaSub_pt,raw_pt,jet_area,jet_eta,centrality,R"

    say "  wide scan: minimal kinematic/background feature set"
    auau_ml_run_python "${BASE}/scripts/train_auau_jet_residual_bdt.py" \
      --input "@${manifest}" \
      --outdir "${outdir}/minimal" \
      --model-name "jetResidualBDT_minimal_allCent_tmva.root" \
      --features "reco_areaSub_pt,jet_eta,centrality,R"

    say "  wide scan: core plus vertex and recoil-angle diagnostics"
    auau_ml_run_python "${BASE}/scripts/train_auau_jet_residual_bdt.py" \
      --input "@${manifest}" \
      --outdir "${outdir}/core_dphi_vz" \
      --model-name "jetResidualBDT_coreDphiVz_allCent_tmva.root" \
      --features "reco_areaSub_pt,raw_pt,jet_area,jet_eta,centrality,R,vertexz,dphi_gamma_jet"
  else
    auau_ml_run_python "${BASE}/scripts/train_auau_jet_residual_bdt.py" --input "@${manifest}" --outdir "$outdir"
  fi
}

auau_jetml_run_condor_do_all() {
  local gs_doall="$GROUP_SIZE"
  if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
    gs_doall="7"
  fi

  local stamp master_yaml train_yaml cfg_tag outdir
  stamp="$(date +%Y%m%d_%H%M%S)"
  master_yaml="$(sim_yaml_master_path)"
  train_yaml="$(auau_jetml_make_training_yaml "$stamp" "$master_yaml")"
  cfg_tag="$(auau_bdt_cfg_tag_from_yaml "$train_yaml")"
  outdir="${AUAU_BDT_MODEL_BASE}/jetResidual_${stamp}"

  local -a signal_samples background_samples
  mapfile -t signal_samples < <( auau_bdt_signal_samples )
  mapfile -t background_samples < <( auau_bdt_background_samples )

  say "AuAu embedded JetML residual CONDOR extraction"
  say "  YAML override      : ${train_yaml}"
  say "  cfg tag            : ${cfg_tag}"
  say "  output base        : ${AUAU_BDT_DEST_BASE}"
  say "  groupSize          : ${gs_doall}"
  say "  photon+jet samples : ${signal_samples[*]}"
  say "  inclusive samples  : ${background_samples[*]}"
  echo

  local samp
  for samp in "${signal_samples[@]}"; do
    say "[JetML condor] photon+jet extraction sample=${samp}"
    RJ_CONFIG_YAML="$train_yaml" \
    RJ_SIMEMBED_DEST_BASE="$AUAU_BDT_DEST_BASE" \
    "$0" isSimEmbedded condorDoAll groupSize "$gs_doall" "SAMPLE=${samp}"
  done
  for samp in "${background_samples[@]}"; do
    say "[JetML condor] inclusive extraction sample=${samp}"
    RJ_CONFIG_YAML="$train_yaml" \
    RJ_SIMEMBED_DEST_BASE="$AUAU_BDT_DEST_BASE" \
    "$0" isSimEmbeddedInclusive condorDoAll groupSize "$gs_doall" "SAMPLE=${samp}"
  done

  say "Condor extraction was submitted. After jobs finish, train from the produced ROOT files with:"
  echo
  printf '  mkdir -p %q\n' "$outdir"
  printf '  find %q -path %q -name %q | sort > %q\n' "$AUAU_BDT_DEST_BASE" "*/${cfg_tag}/*" "RecoilJets_*.root" "${outdir}/jetml_training_roots.list"
  printf '  %q --input @%q --outdir %q\n' "${BASE}/scripts/train_auau_jet_residual_bdt.py" "${outdir}/jetml_training_roots.list" "$outdir"
}

ml_yaml_set_scalar() {
  local yaml="$1" key="$2" value="$3"
  local tmp="${yaml}.tmp"
  sed -E "s|^([[:space:]]*${key}:).*|\\1 ${value}|" "$yaml" > "$tmp"
  mv -f "$tmp" "$yaml"
}

auau_ml_first_matrix_yaml() {
  local master_yaml="$1" tag="$2" stamp="$3" preselection="$4" tight="$5" nonTight="$6"
  local -a pts fracs vzs cones fixeds uepipes
  mapfile -t pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
  mapfile -t fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
  mapfile -t vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
  mapfile -t cones < <( yaml_get_values "coneR" "$master_yaml" )
  mapfile -t fixeds < <( yaml_get_values "fixedGeV" "$master_yaml" 2>/dev/null )
  (( ${#pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
  (( ${#fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
  (( ${#vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
  (( ${#cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
  (( ${#fixeds[@]} )) || fixeds=( "2.0" )

  read_uepipe_modes "$master_yaml" "simembedded"
  uepipes=( "${uepipe_modes[@]}" )
  (( ${#uepipes[@]} )) || uepipes=( "noSub" )

  sim_make_yaml_override \
    "$master_yaml" \
    "${pts[0]}" "${fracs[0]}" "${vzs[0]}" "${cones[0]}" \
    "false" "${fixeds[0]}" "${uepipes[0]}" \
    "$preselection" "$tight" "$nonTight" \
    "$tag" "$stamp" "false"
}

auau_ml_parse_integration_controls() {
  AUAU_ML_LOCAL_EVENTS="${LOCAL_EVENTS}"
  AUAU_ML_LOCAL_NFILES="${RJ_ML_LOCAL_NFILES:-}"
  AUAU_ML_LOCAL_VERBOSITY="10"
  AUAU_ML_RESUME_BASE="${RJ_ML_RESUME_BASE:-}"
  local t
  local expect_resume_path=0
  for t in "${tokens[@]}"; do
    if [[ "$expect_resume_path" -eq 1 ]]; then
      AUAU_ML_RESUME_BASE="$t"
      expect_resume_path=0
      continue
    fi
    if [[ "$t" =~ ^VERBOSE=([0-9]+)$ ]]; then
      AUAU_ML_LOCAL_VERBOSITY="${BASH_REMATCH[1]}"
    elif [[ "$t" =~ ^NFILES=([0-9]+)$ ]]; then
      AUAU_ML_LOCAL_NFILES="${BASH_REMATCH[1]}"
    elif [[ "$t" =~ ^FILES=([0-9]+)$ ]]; then
      AUAU_ML_LOCAL_NFILES="${BASH_REMATCH[1]}"
    elif [[ "$t" == "resume" ]]; then
      expect_resume_path=1
    elif [[ "$t" =~ ^(BASE|RUN_BASE|RESUME_BASE)=(.+)$ ]]; then
      AUAU_ML_RESUME_BASE="${BASH_REMATCH[2]}"
    elif [[ "$t" =~ ^[0-9]+$ ]]; then
      AUAU_ML_LOCAL_EVENTS="$t"
    fi
  done
  [[ "$AUAU_ML_LOCAL_EVENTS" =~ ^[0-9]+$ ]] || { err "trainMLAll event count must be an integer"; exit 2; }
  if [[ -z "$AUAU_ML_LOCAL_NFILES" ]]; then
    AUAU_ML_LOCAL_NFILES="$(auau_ml_default_nfiles_for_events "$AUAU_ML_LOCAL_EVENTS")"
  fi
  [[ "$AUAU_ML_LOCAL_NFILES" =~ ^[0-9]+$ && "$AUAU_ML_LOCAL_NFILES" -ge 1 ]] || { err "NFILES must be a positive integer"; exit 2; }
}

auau_ml_run_sample_set_local() {
  local dataset="$1"
  local yaml="$2"
  local output_base="$3"
  local nevents="$4"
  local verbosity="$5"
  shift 5
  local -a samples=( "$@" )

  local samp
  for samp in "${samples[@]}"; do
    say "[ML integration] dataset=${dataset} sample=${samp}"
    RJ_CONFIG_YAML="$yaml" \
    RJ_SIM_LOCAL_GROUPED=1 \
    RJ_SIM_LOCAL_NFILES="$AUAU_ML_LOCAL_NFILES" \
    RJ_LOCAL_SIM_OUTPUT_BASE="$output_base" \
    "$0" "$dataset" local "$nevents" "VERBOSE=${verbosity}" "SAMPLE=${samp}"
  done
}

auau_ml_root_count() {
  local dir="$1"
  if [[ -d "$dir" ]]; then
    find "$dir" -type f -name "RecoilJets_*.root" | wc -l | tr -d '[:space:]'
  else
    echo 0
  fi
}

auau_ml_run_all_local() {
  auau_ml_parse_integration_controls

  local stamp run_base master_yaml row_yaml jetml_yaml tight_manifest jetml_manifest
  local tight_model_dir npb_model_dir jetml_model_dir tight_apply_yaml jetml_apply_yaml
  if [[ "${TRAIN_MODE:-local}" == "resume" ]]; then
    [[ -n "${AUAU_ML_RESUME_BASE:-}" ]] || { err "trainMLAll resume requires a run directory path"; exit 2; }
    run_base="$AUAU_ML_RESUME_BASE"
    [[ -d "$run_base" ]] || { err "Resume directory does not exist: $run_base"; exit 2; }
    stamp="$(basename "$run_base")"
    stamp="${stamp#mlIntegration_}"
  else
    stamp="$(date +%Y%m%d_%H%M%S)"
    run_base="${RJ_ML_LOCAL_BASE:-${BASE}/local_ml_pipeline_tests}/mlIntegration_${stamp}"
  fi
  master_yaml="$(sim_yaml_master_path)"
  mkdir -p "$run_base"/{extract,models,apply,manifests,qa}

  auau_ml_check_python_deps

  local -a signal_samples background_samples
  mapfile -t signal_samples < <( auau_bdt_signal_samples )
  mapfile -t background_samples < <( auau_bdt_background_samples )

  say "AuAu ML integrated LOCAL pipeline"
  say "  output base        : ${run_base}"
  if [[ "${TRAIN_MODE:-local}" == "resume" ]]; then
    say "  resume mode        : reuse existing extraction ROOT files when present"
  fi
  say "  target events/sample: ${AUAU_ML_LOCAL_EVENTS}"
  say "  grouped files/sample: ${AUAU_ML_LOCAL_NFILES} (assuming ${RJ_SIM_LOCAL_EVENTS_PER_FILE:-1000} events/input file unless overridden)"
  say "  photon+jet samples : ${signal_samples[*]}"
  say "  inclusive samples  : ${background_samples[*]}"
  say "  reference style    : PPG12 fixed feature order + TMVA/RBDT export + analysis re-apply"
  echo

  row_yaml="$(auau_bdt_make_training_yaml "allMLPhotonRows" "$stamp" "$master_yaml")"
  say "Stage A1: extract shared photon-ID rows for tight BDT and NPB sanity"
  say "  row YAML: ${row_yaml}"
  local photon_root_count
  photon_root_count="$(auau_ml_root_count "${run_base}/extract/photon_rows")"
  if [[ "${TRAIN_MODE:-local}" == "resume" && "$photon_root_count" -gt 0 ]]; then
    say "  found existing photon-row ROOT files: ${photon_root_count}; skipping photon extraction"
  else
    auau_ml_run_sample_set_local "isSimEmbedded" "$row_yaml" "${run_base}/extract/photon_rows/simembedded" "$AUAU_ML_LOCAL_EVENTS" "$AUAU_ML_LOCAL_VERBOSITY" "${signal_samples[@]}"
    auau_ml_run_sample_set_local "isSimEmbeddedInclusive" "$row_yaml" "${run_base}/extract/photon_rows/simembeddedinclusive" "$AUAU_ML_LOCAL_EVENTS" "$AUAU_ML_LOCAL_VERBOSITY" "${background_samples[@]}"
  fi

  tight_manifest="${run_base}/manifests/photon_training_roots.list"
  auau_bdt_collect_roots "${run_base}/extract/photon_rows" "$tight_manifest"
  auau_ml_validate_manifest_tree "$tight_manifest" "AuAuPhotonIDTrainingTree" "shared photon-ID rows" "${run_base}/qa/photon_tree_validation.txt"

  tight_model_dir="${run_base}/models/tight"
  npb_model_dir="${run_base}/models/npb"
  say "Stage B1: train tight BDT variants from shared photon rows"
  RJ_AUAU_BDT_LOCAL_WIDE=1 auau_bdt_train_from_manifest "tight" "$tight_manifest" "$tight_model_dir"

  say "Stage B2: run NPB training sanity or PPG12 data-tagged training"
  local npb_data_manifest=""
  if [[ "${RJ_AUAU_BDT_NPB_DATA_LOCAL:-0}" == "1" ]]; then
    local npb_data_yaml="${run_base}/manifests/npb_data_yaml.path"
    local data_yaml data_root
    data_yaml="$(auau_bdt_make_training_yaml "npbData" "$stamp" "$master_yaml")"
    printf "%s\n" "$data_yaml" > "$npb_data_yaml"
    data_root="${run_base}/extract/npb_data"
    say "  extracting AuAu timing-tagged NPB rows with PPG12-style labels"
    say "  row YAML: ${data_yaml}"
    RJ_CONFIG_YAML="$data_yaml" \
    RJ_LOCAL_DATA_OUTPUT_BASE="$data_root" \
    "$0" isAuAu local "$AUAU_ML_LOCAL_EVENTS" "VERBOSE=${AUAU_ML_LOCAL_VERBOSITY}"
    npb_data_manifest="${run_base}/manifests/npb_data_training_roots.list"
    auau_bdt_collect_roots "$data_root" "$npb_data_manifest"
    local npb_data_report="${run_base}/qa/npb_data_tree_validation.txt"
    auau_ml_validate_manifest_tree "$npb_data_manifest" "AuAuPhotonIDTrainingTree" "AuAu data NPB rows" "$npb_data_report"
    local npb_data_entries
    npb_data_entries="$(awk -F= '/^entries_total=/{print $2}' "$npb_data_report" 2>/dev/null | tail -n 1)"
    if [[ "${npb_data_entries:-0}" == "0" ]]; then
      warn "AuAu data NPB extraction produced zero tagged rows. Skipping NPB training/apply for this local smoke pass."
      npb_data_manifest=""
      auau_bdt_train_from_manifest "npb" "$tight_manifest" "$npb_model_dir"
    else
      RJ_AUAU_BDT_NPB_INPUTS="@${npb_data_manifest}" auau_bdt_train_from_manifest "npb" "$tight_manifest" "$npb_model_dir"
    fi
  else
    auau_bdt_train_from_manifest "npb" "$tight_manifest" "$npb_model_dir"
  fi

  jetml_yaml="$(auau_jetml_make_training_yaml "$stamp" "$master_yaml")"
  say "Stage A2: extract JetML residual rows"
  say "  row YAML: ${jetml_yaml}"
  local jetml_root_count
  jetml_root_count="$(auau_ml_root_count "${run_base}/extract/jetml_rows")"
  if [[ "${TRAIN_MODE:-local}" == "resume" && "$jetml_root_count" -gt 0 ]]; then
    say "  found existing JetML-row ROOT files: ${jetml_root_count}; skipping JetML extraction"
  else
    auau_ml_run_sample_set_local "isSimEmbedded" "$jetml_yaml" "${run_base}/extract/jetml_rows/simembedded" "$AUAU_ML_LOCAL_EVENTS" "$AUAU_ML_LOCAL_VERBOSITY" "${signal_samples[@]}"
    auau_ml_run_sample_set_local "isSimEmbeddedInclusive" "$jetml_yaml" "${run_base}/extract/jetml_rows/simembeddedinclusive" "$AUAU_ML_LOCAL_EVENTS" "$AUAU_ML_LOCAL_VERBOSITY" "${background_samples[@]}"
  fi

  jetml_manifest="${run_base}/manifests/jetml_training_roots.list"
  auau_bdt_collect_roots "${run_base}/extract/jetml_rows" "$jetml_manifest"
  auau_ml_validate_manifest_tree "$jetml_manifest" "JetResidualMLTrainingTree" "JetML residual rows" "${run_base}/qa/jetml_tree_validation.txt"

  jetml_model_dir="${run_base}/models/jetResidual"
  mkdir -p "$jetml_model_dir"
  say "Stage B3: train JetML residual feature-set scan"
  RJ_JET_ML_LOCAL_WIDE=1 \
  AUAU_BDT_MODEL_BASE="${run_base}/models" \
  auau_ml_run_python "${BASE}/scripts/train_auau_jet_residual_bdt.py" \
    --input "@${jetml_manifest}" \
    --outdir "${jetml_model_dir}/core" \
    --model-name "jetResidualBDT_core_allCent_tmva.root" \
    --features "reco_areaSub_pt,raw_pt,jet_area,jet_eta,centrality,R"
  auau_ml_run_python "${BASE}/scripts/train_auau_jet_residual_bdt.py" \
    --input "@${jetml_manifest}" \
    --outdir "${jetml_model_dir}/minimal" \
    --model-name "jetResidualBDT_minimal_allCent_tmva.root" \
    --features "reco_areaSub_pt,jet_eta,centrality,R"
  auau_ml_run_python "${BASE}/scripts/train_auau_jet_residual_bdt.py" \
    --input "@${jetml_manifest}" \
    --outdir "${jetml_model_dir}/core_dphi_vz" \
    --model-name "jetResidualBDT_coreDphiVz_allCent_tmva.root" \
    --features "reco_areaSub_pt,raw_pt,jet_area,jet_eta,centrality,R,vertexz,dphi_gamma_jet"

  local tight_model="${tight_model_dir}/auau_tight_bdt_nominal_allCent_tmva.root"
  local npb_model="${npb_model_dir}/auau_npb_bdt_allCent_tmva.root"
  local jetml_model="${jetml_model_dir}/core/jetResidualBDT_core_allCent_tmva.root"

  if [[ -s "$npb_model" ]]; then
    local npb_apply_yaml
    npb_apply_yaml="$(auau_ml_first_matrix_yaml "$master_yaml" "mlApply_npb_${stamp}" "$stamp" "auauOnlyNPB" "reference" "reference")"
    ml_yaml_set_scalar "$npb_apply_yaml" "auau_npb_model_file" "$npb_model"
    ml_yaml_set_scalar "$npb_apply_yaml" "auau_bdt_training_tree" "false"
    ml_yaml_set_scalar "$npb_apply_yaml" "jet_ml_training_tree" "false"
    ml_yaml_set_scalar "$npb_apply_yaml" "jet_ml_correction_enabled" "false"

    say "Stage C0: re-apply trained NPB model through preselection=auauOnlyNPB"
    say "  model: ${npb_model}"
    say "  YAML : ${npb_apply_yaml}"
    auau_ml_run_sample_set_local "isSimEmbedded" "$npb_apply_yaml" "${run_base}/apply/npb/simembedded" "$AUAU_ML_LOCAL_EVENTS" "$AUAU_ML_LOCAL_VERBOSITY" "${signal_samples[@]}"
    auau_ml_run_sample_set_local "isSimEmbeddedInclusive" "$npb_apply_yaml" "${run_base}/apply/npb/simembeddedinclusive" "$AUAU_ML_LOCAL_EVENTS" "$AUAU_ML_LOCAL_VERBOSITY" "${background_samples[@]}"
  else
    warn "Skipping NPB apply stage because model is missing: ${npb_model}"
  fi

  if [[ -s "$tight_model" ]]; then
    tight_apply_yaml="$(auau_ml_first_matrix_yaml "$master_yaml" "mlApply_tight_${stamp}" "$stamp" "reference" "auauEmbeddedBDT" "auauBDTSideband")"
    ml_yaml_set_scalar "$tight_apply_yaml" "auau_tight_bdt_model_file" "$tight_model"
    ml_yaml_set_scalar "$tight_apply_yaml" "auau_bdt_training_tree" "false"
    ml_yaml_set_scalar "$tight_apply_yaml" "jet_ml_training_tree" "false"
    ml_yaml_set_scalar "$tight_apply_yaml" "jet_ml_correction_enabled" "false"

    say "Stage C1: re-apply trained tight BDT model through production analysis path"
    say "  model: ${tight_model}"
    say "  YAML : ${tight_apply_yaml}"
    auau_ml_run_sample_set_local "isSimEmbedded" "$tight_apply_yaml" "${run_base}/apply/tight" "$AUAU_ML_LOCAL_EVENTS" "$AUAU_ML_LOCAL_VERBOSITY" "${signal_samples[@]}"
  else
    warn "Skipping tight apply stage because model is missing: ${tight_model}"
  fi

  if [[ -s "$jetml_model" ]]; then
    jetml_apply_yaml="$(auau_ml_first_matrix_yaml "$master_yaml" "mlApply_jetResidual_${stamp}" "$stamp" "reference" "reference" "reference")"
    ml_yaml_set_scalar "$jetml_apply_yaml" "auau_bdt_training_tree" "false"
    ml_yaml_set_scalar "$jetml_apply_yaml" "jet_ml_training_tree" "false"
    ml_yaml_set_scalar "$jetml_apply_yaml" "jet_ml_correction_enabled" "true"
    ml_yaml_set_scalar "$jetml_apply_yaml" "jet_ml_model_file" "$jetml_model"
    ml_yaml_set_scalar "$jetml_apply_yaml" "jet_ml_features" "[reco_areaSub_pt, raw_pt, jet_area, jet_eta, centrality, R]"

    say "Stage C2: re-apply trained JetML core residual model through production analysis path"
    say "  model: ${jetml_model}"
    say "  YAML : ${jetml_apply_yaml}"
    auau_ml_run_sample_set_local "isSimEmbedded" "$jetml_apply_yaml" "${run_base}/apply/jetResidual" "$AUAU_ML_LOCAL_EVENTS" "$AUAU_ML_LOCAL_VERBOSITY" "${signal_samples[@]}"
  else
    warn "Skipping JetML apply stage because model is missing: ${jetml_model}"
  fi

  {
    echo "run_base=${run_base}"
    echo "target_events_per_sample=${AUAU_ML_LOCAL_EVENTS}"
    echo "grouped_input_files_per_sample=${AUAU_ML_LOCAL_NFILES}"
    echo "assumed_events_per_input_file=${RJ_SIM_LOCAL_EVENTS_PER_FILE:-1000}"
    echo "photon_training_manifest=${tight_manifest}"
    echo "npb_data_training_manifest=${npb_data_manifest:-}"
    echo "jetml_training_manifest=${jetml_manifest}"
    echo "tight_model=${tight_model}"
    echo "npb_model=${npb_model}"
    echo "jetml_model=${jetml_model}"
    echo "qa_dir=${run_base}/qa"
    echo "pull_command=./scripts/sftp_get_recoiljets_outputs.sh mlIntegrationLatest"
  } > "${run_base}/README.local_ml_pipeline.txt"

  say "Integrated local ML pipeline complete."
  say "  run base : ${run_base}"
  say "  pull from Mac with:"
  say "    ./scripts/sftp_get_recoiljets_outputs.sh mlIntegrationLatest"
}

scaled_trigger_prepare_artifacts() {
  [[ "$DATASET" == "isAuAu" ]] || { err "scaledTriggerStudy is valid only for isAuAu"; exit 2; }
  [[ -s "$SCALED_TRIGGER_RUNLIST" ]] || {
    err "Missing scaled-trigger run list: ${SCALED_TRIGGER_RUNLIST}"
    err "Run first on SDCC: ./scripts/make_dstListsData.sh auau QA scaledTriggerAna"
    exit 80
  }
  [[ -s "$SCALED_TRIGGER_CONFIG_SRC" ]] || {
    err "Missing scaled-trigger CONFIG table: ${SCALED_TRIGGER_CONFIG_SRC}"
    err "Run first on SDCC: ./scripts/make_dstListsData.sh auau QA scaledTriggerAna"
    exit 80
  }
  if [[ ! -s "$SCALED_TRIGGER_CONFIG_EXPECTED" ]]; then
    cp -f "$SCALED_TRIGGER_CONFIG_SRC" "$SCALED_TRIGGER_CONFIG_EXPECTED"
    say "Created merge config alias: ${SCALED_TRIGGER_CONFIG_EXPECTED}"
  fi
}

scaled_trigger_prepare_single_yaml() {
  local master_yaml="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  local -a data_pts data_fracs data_vzs data_cones
  mapfile -t data_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
  mapfile -t data_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
  mapfile -t data_vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
  mapfile -t data_cones < <( yaml_get_values "coneR" "$master_yaml" )
  (( ${#data_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
  (( ${#data_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
  (( ${#data_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
  (( ${#data_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }

  build_iso_modes "$master_yaml"
  read_uepipe_modes "$master_yaml" "$TAG"

  local data_pt="${data_pts[0]}"
  local data_frac="${data_fracs[0]}"
  local data_vz="${data_vzs[0]}"
  local data_cone="${data_cones[0]}"
  local uepipe="${uepipe_modes[0]}"
  SCALED_TRIGGER_VZ_CUT="$data_vz"

  local dpt_tag dfrac_tag dvz_tag dcone_tag
  dpt_tag="jetMinPt$(sim_pt_tag "$data_pt")"
  dfrac_tag="$(sim_b2b_tag "$data_frac")"
  dvz_tag="$(sim_vz_tag "$data_vz")"
  dcone_tag="$(sim_cone_tag "$data_cone")"

  SCALED_TRIGGER_CFG_TAG="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_base_tags[0]}"
  (( uepipe_in_tag )) && SCALED_TRIGGER_CFG_TAG="${SCALED_TRIGGER_CFG_TAG}_${uepipe}"
  SCALED_TRIGGER_CFG_TAG="${SCALED_TRIGGER_CFG_TAG}_${iso_selection_tags[0]}"
  SCALED_TRIGGER_OUTPUT_TAG="${SCALED_TRIGGER_CFG_TAG}${SCALED_TRIGGER_OUTPUT_SUFFIX}"

  SCALED_TRIGGER_YAML="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${SCALED_TRIGGER_CFG_TAG}_scaledTriggerStudy.yaml"
  mkdir -p "$SIM_YAML_OVERRIDE_DIR"
  sed -E \
    -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${data_pt}|" \
    -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${data_frac}|" \
    -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${data_vz}|" \
    -e "s|^([[:space:]]*coneR:).*|\\1 ${data_cone}|" \
    -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${iso_sliding[0]}|" \
    -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${iso_fixed[0]}|" \
    -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${uepipe}|" \
    "$master_yaml" > "$SCALED_TRIGGER_YAML"
  pin_photon_id_scalars_in_yaml "$SCALED_TRIGGER_YAML" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}"

  say "scaledTriggerStudy single config:"
  say "  cfg_tag       : ${SCALED_TRIGGER_CFG_TAG}"
  say "  output tag    : ${SCALED_TRIGGER_OUTPUT_TAG}"
  say "  YAML override : ${SCALED_TRIGGER_YAML}"
  say "  run list      : ${SCALED_TRIGGER_RUNLIST} ($(grep -cE '^[0-9]+' "$SCALED_TRIGGER_RUNLIST") runs)"
}

# ------------------------ Parse CLI ------------------------
[[ $# -ge 1 ]] || usage
resolve_dataset "$1"

# Parse remaining tokens (order-agnostic):
ACTION=""
TRAIN_MODE=""
DRYRUN=0
SIM_SAMPLE_EXPLICIT=0
GROUP_SIZE_EXPLICIT=0
MAX_JOBS_EXPLICIT=0
tokens=( "${@:2}" )
for (( idx=0; idx<${#tokens[@]}; idx++ )); do
  tok="${tokens[$idx]}"
  case "$tok" in
    trainTightBDT|trainNPB|trainJetMLResidual|trainMLAll)
      ACTION="$tok"
      ;;
    scaledTriggerStudy)
      ACTION="$tok"
      ;;
    local|localTest|condorDoAll|condorDoAllSmoke|condorDoAllDirect|condorDoAllFromScratch|condorHistFromPool|resume|smokeTestFirstPass|smokeTestSecondPass|smokeTestApplyExisting)
      if [[ "$ACTION" == trainTightBDT || "$ACTION" == trainNPB || "$ACTION" == trainJetMLResidual || "$ACTION" == trainMLAll || "$ACTION" == scaledTriggerStudy ]]; then
        TRAIN_MODE="$tok"
      else
        ACTION="$tok"
      fi
      ;;
    checkModels|workflowCheck|isLocalIsoPing|condor|splitGoldenRunList|condorTest)
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
      MAX_JOBS_EXPLICIT=1
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
if [[ "$ACTION" == trainTightBDT || "$ACTION" == trainNPB || "$ACTION" == trainJetMLResidual || "$ACTION" == trainMLAll ]]; then
  [[ "$DATASET" == "isSimEmbeddedAndInclusive" ]] || { err "${ACTION} is valid only as: $0 isSimEmbeddedAndInclusive ${ACTION} <local|condorDoAll>"; exit 2; }
  [[ -n "$TRAIN_MODE" ]] || TRAIN_MODE="local"
fi
if [[ "$ACTION" == "scaledTriggerStudy" ]]; then
  [[ "$DATASET" == "isAuAu" ]] || { err "scaledTriggerStudy is valid only as: $0 isAuAu scaledTriggerStudy <local|condorDoAll>"; exit 2; }
  [[ -n "$TRAIN_MODE" ]] || TRAIN_MODE="local"
  [[ "$TRAIN_MODE" == "local" || "$TRAIN_MODE" == "condorDoAll" ]] || { err "scaledTriggerStudy mode must be local or condorDoAll, got '${TRAIN_MODE}'"; exit 2; }
  if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
    GROUP_SIZE=20
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
  if [[ -n "${TRAIN_MODE:-}" ]]; then
    say "Training mode=${BOLD}${TRAIN_MODE}${RST}"
  fi
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
  scaledTriggerStudy)
    scaled_trigger_prepare_artifacts
    scaled_trigger_prepare_single_yaml
    export RJ_CONFIG_YAML="$SCALED_TRIGGER_YAML"
    export RJ_SCALED_TRIGGER_RUNLIST="$SCALED_TRIGGER_RUNLIST"
    export RJ_SCALED_TRIGGER_STUDY_ONLY=1
    export RJ_REQUEST_MEMORY="${RJ_SCALED_TRIGGER_REQUEST_MEMORY:-1000MB}"
    export RJ_SUBMIT_EXTRA_ENV="RJ_SCALED_TRIGGER_STUDY_ONLY=1;RJ_SCALED_TRIGGER_RUNLIST=${SCALED_TRIGGER_RUNLIST};RJ_SCALED_TRIGGER_VZ_MAX_CM=${RJ_SCALED_TRIGGER_VZ_MAX_CM:-${SCALED_TRIGGER_VZ_CUT}}"

    DEST_BASE="${AA_DEST_BASE}/${SCALED_TRIGGER_OUTPUT_TAG}"

    case "$TRAIN_MODE" in
      local)
        RJV="10"
        nevt="$LOCAL_EVENTS"
        for tok in "${tokens[@]}"; do
          if [[ "$tok" =~ ^[0-9]+$ ]]; then nevt="$tok"; fi
          if [[ "$tok" == VERBOSE=* ]]; then RJV="${tok#VERBOSE=}"; fi
        done

        r8="$(awk 'NF && $1 !~ /^#/ { printf "%08d\n", $1; exit }' "$SCALED_TRIGGER_RUNLIST")"
        [[ -n "$r8" ]] || { err "No run found in ${SCALED_TRIGGER_RUNLIST}"; exit 81; }
        mkdir -p "$STAGE_DIR"
        mapfile -t groups < <( make_groups "$r8" "$GROUP_SIZE" )
        (( ${#groups[@]} )) || { err "No input groups produced for run ${r8}"; exit 82; }
        glist="${groups[0]}"

        say "scaledTriggerStudy local smoke test"
        say "  run          : ${r8}"
        say "  groupSize    : ${GROUP_SIZE}"
        say "  list chunk   : ${glist}"
        say "  events       : ${nevt}"
        say "  DEST_BASE    : ${DEST_BASE}"
        say "  memory target: local run (condor default would be ${RJ_REQUEST_MEMORY})"
        say "Invoking wrapper locally..."

        RJ_DATASET="$DATASET" RJ_VERBOSITY="$RJV" \
        RJ_CONFIG_YAML="$SCALED_TRIGGER_YAML" \
        RJ_SCALED_TRIGGER_STUDY_ONLY=1 \
        RJ_SCALED_TRIGGER_RUNLIST="$SCALED_TRIGGER_RUNLIST" \
        RJ_SCALED_TRIGGER_VZ_MAX_CM="${RJ_SCALED_TRIGGER_VZ_MAX_CM:-${SCALED_TRIGGER_VZ_CUT}}" \
        bash "$EXE" "$r8" "$glist" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
        ;;

      condorDoAll)
        cleanup_bulk_snapshots_for_tag
        create_pipeline_snapshot "auau" "$(date +%Y%m%d_%H%M%S)"
        say "scaledTriggerStudy condorDoAll"
        say "  runs         : $(grep -cE '^[0-9]+' "$SCALED_TRIGGER_RUNLIST")"
        say "  groupSize    : ${GROUP_SIZE}"
        say "  request mem  : ${RJ_REQUEST_MEMORY}"
        say "  DEST_BASE    : ${DEST_BASE}"
        submit_condor "$SCALED_TRIGGER_RUNLIST" ""
        ;;
    esac
    ;;

  trainTightBDT|trainNPB)
    task="tight"
    [[ "$ACTION" == "trainNPB" ]] && task="npb"
    case "${TRAIN_MODE:-local}" in
      local|localTest)
        auau_bdt_run_local "$task"
        ;;
      condorDoAll)
        auau_bdt_run_condor_do_all "$task"
        ;;
      smokeTestFirstPass)
        [[ "$task" == "tight" ]] || { err "smokeTestFirstPass is currently only implemented for trainTightBDT"; exit 2; }
        if [[ -z "${RJ_MAX_JOBS:-}" ]]; then
          MAX_JOBS="${RJ_AUAU_TIGHT_BDT_SMOKE_MAX_JOBS_PER_SAMPLE:-4}"
        fi
        if [[ -z "${RJ_AUAU_BDT_DEST_BASE:-}" ]]; then
          AUAU_BDT_DEST_BASE="${BASE}/condor_tight_bdt_smoke/tightSmokeFirst_$(date +%Y%m%d_%H%M%S)"
        fi
        say "tight BDT smokeTestFirstPass caps extraction at MAX_JOBS=${MAX_JOBS} per sample."
        say "With default groupSize=7 and 1000 events/input file this is about $(( MAX_JOBS * 7 ))k events per sample."
        say "Unique smoke output base: ${AUAU_BDT_DEST_BASE}"
        auau_bdt_run_condor_do_all "$task"
        ;;
      smokeTestSecondPass)
        [[ "$task" == "tight" ]] || { err "smokeTestSecondPass is currently only implemented for trainTightBDT"; exit 2; }
        auau_tight_bdt_smoke_second_pass
        ;;
      smokeTestApplyExisting)
        [[ "$task" == "tight" ]] || { err "smokeTestApplyExisting is currently only implemented for trainTightBDT"; exit 2; }
        auau_tight_bdt_smoke_apply_existing
        ;;
      *)
        err "${ACTION} mode must be local, localTest, condorDoAll, smokeTestFirstPass, smokeTestSecondPass, or smokeTestApplyExisting; got '${TRAIN_MODE}'"
        exit 2
        ;;
    esac
    ;;

  trainJetMLResidual)
    case "${TRAIN_MODE:-local}" in
      local)
        auau_jetml_run_local
        ;;
      condorDoAll)
        auau_jetml_run_condor_do_all
        ;;
      *)
        err "${ACTION} mode must be local or condorDoAll, got '${TRAIN_MODE}'"
        exit 2
        ;;
    esac
    ;;

  trainMLAll)
    case "${TRAIN_MODE:-local}" in
      local|resume)
        auau_ml_run_all_local
        ;;
      *)
        err "${ACTION} currently supports local or resume mode, got '${TRAIN_MODE}'"
        exit 2
        ;;
    esac
    ;;

  CHECKJOBS)
    # Dry-run only
    if [[ "$IS_SIM" -eq 1 ]]; then
      check_jobs_sim
    else
      check_jobs_all
    fi
    exit 0
    ;;

  workflowCheck)
    workflow_check
    exit 0
    ;;

  checkModels)
    [[ "$DATASET" == "isSim" ]] || { err "checkModels is only valid as: $0 isSim checkModels"; exit 2; }

    nevt="${RJ_CHECKMODELS_EVENTS:-1000}"
    RJV="${RJ_CHECKMODELS_VERBOSITY:-2}"
    rest=( "${@:3}" )
    for t in "${rest[@]}"; do
      if [[ "$t" =~ ^VERBOSE=([0-9]+)$ ]]; then
        RJV="${BASH_REMATCH[1]}"
      elif [[ "$t" =~ ^[0-9]+$ ]]; then
        nevt="$t"
      fi
    done
    [[ "$nevt" =~ ^[0-9]+$ ]] || { err "checkModels event count must be an integer"; exit 2; }

    master_yaml="$(sim_yaml_master_path)"
    [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

    mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
    mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
    mapfile -t sim_vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
    read_replay_cones "$master_yaml"
    read_capture_cones "$master_yaml"
    sim_cones=( "${replay_cones[@]}" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
    (( ${#capture_cones[@]} )) || { err "No capture cone could be derived from $master_yaml"; exit 72; }
    read_uepipe_modes "$master_yaml" "$TAG"
    check_pt="${sim_pts[0]}"
    check_frac="${sim_fracs[0]}"
    check_vz="${sim_vzs[0]}"
    check_cone="${sim_cones[0]}"
    check_fixed="2.0"
    check_sliding="false"
    check_iso_base_tag="$(sim_iso_tag "$check_sliding" "$check_fixed")"

    SIM_SAMPLE="run28_photonjet20"
    SIM_SAMPLE_EXPLICIT=1
    SIM_STAGE_NAMESPACE="model_check"
    SIM_DEST_BASE_RESOLVED="${BASE}/local_sim_outputs/model_check"
    [[ "$SIM_DEST_BASE_RESOLVED" == "${BASE}/local_sim_outputs/model_check" ]] || { err "Refusing to clean unexpected model-check path: ${SIM_DEST_BASE_RESOLVED}"; exit 74; }
    say "Cleaning previous model-check outputs: ${SIM_DEST_BASE_RESOLVED}"
    rm -rf "${SIM_DEST_BASE_RESOLVED}"
    mkdir -p "${SIM_DEST_BASE_RESOLVED}"
    model_check_manifest="${SIM_DEST_BASE_RESOLVED}/checkModels_manifest.txt"
    : > "$model_check_manifest"

    say "SIM model check (local photon20 diagnostic)"
    say "  YAML master : ${master_yaml}"
    say "  events      : ${nevt}"
    say "  sample      : ${SIM_SAMPLE}"
    say "  scalar axes : jet_pt_min=${check_pt}  back_to_back_pi_fraction=${check_frac}  vz_cut_cm=${check_vz}  coneR=${check_cone}  iso=${check_iso_base_tag}"
    say "  ID passes   : reference/reference/reference, newPPG12/newPPG12/newPPG12"
    say "  dest base   : ${SIM_DEST_BASE_RESOLVED}"
    echo

    check_preselections=( "reference" "newPPG12" )
    check_tights=( "reference" "newPPG12" )
    check_nonTights=( "reference" "newPPG12" )

    for check_idx in "${!check_preselections[@]}"; do
      for uepipe in "${uepipe_modes[@]}"; do
          preselection="${check_preselections[$check_idx]}"
          tight="${check_tights[$check_idx]}"
          nonTight="${check_nonTights[$check_idx]}"
          selection_tag="$(selection_mode_tag "preselection" "$preselection")_$(selection_mode_tag "tight" "$tight")_$(selection_mode_tag "nonTight" "$nonTight")"
          SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$check_pt")_$(sim_b2b_tag "$check_frac")_$(sim_vz_tag "$check_vz")_$(sim_cone_tag "$check_cone")_${check_iso_base_tag}"
          (( uepipe_in_tag )) && SIM_CFG_TAG="${SIM_CFG_TAG}_${uepipe}"
          SIM_CFG_TAG="${SIM_CFG_TAG}_${selection_tag}"
          DEST_BASE="${SIM_DEST_BASE_RESOLVED}/${SIM_CFG_TAG}"
          yaml_override="$(sim_make_yaml_override "$master_yaml" "$check_pt" "$check_frac" "$check_vz" "$check_cone" "$check_sliding" "$check_fixed" "${uepipe}" "$preselection" "$tight" "$nonTight" "$SIM_CFG_TAG" "CHECKMODELS" "false")"

          GROUP_SIZE="1"
          sim_init

          tmp="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_CHECKMODELS_firstfile_grp001.list"
          head -n 1 "$SIM_CLEAN_LIST" > "$tmp"
          [[ -s "$tmp" ]] || { err "No sim entries (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"; exit 30; }

          in_line="$(head -n 1 "$tmp" 2>/dev/null || true)"
          chunk_base="$(basename "$tmp")"
          chunk_tag="${chunk_base%.list}"
          out_root_preview="${DEST_BASE}/${SIM_SAMPLE}/RecoilJets_${DATASET}_${chunk_tag}.root"

          say "----------------------------------------"
          say "SIM checkModels: tag=${SIM_CFG_TAG}  sample=${SIM_SAMPLE}"
          say "  jet_pt_min=${check_pt}  back_to_back_pi_fraction=${check_frac}  vz_cut_cm=${check_vz}  coneR=${check_cone}  iso=${check_iso_base_tag}  uepipe=${uepipe}"
          say "  photon ID    : preselection=${preselection}  tight=${tight}  nonTight=${nonTight}"
          say "  YAML override: ${yaml_override}"
          say "  input line    : ${in_line}"
          say "  temp list     : ${tmp}"
          say "  out ROOT      : ${out_root_preview}"
          say "  wrapper env   : RJ_VERBOSITY=${RJV} RJ_CONFIG_YAML=${yaml_override}"
          find "${DEST_BASE}/${SIM_SAMPLE}" -type f -name "*_CHECKMODELS_*.root" -delete 2>/dev/null || true
          say "Invoking wrapper locally..."

          RJ_VERBOSITY="$RJV" RJ_CONFIG_YAML="$yaml_override" bash "$EXE" "$SIM_SAMPLE" "$tmp" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
          if [[ -s "$out_root_preview" ]]; then
            printf "%s\n" "$out_root_preview" >> "$model_check_manifest"
          else
            warn "Expected output ROOT was not produced or is empty: ${out_root_preview}"
          fi
          echo
      done
    done

    mapfile -t model_check_roots < "$model_check_manifest"
    say "Model-check ROOT files produced: ${#model_check_roots[@]}"
    for rf in "${model_check_roots[@]}"; do
      say "  ${rf}"
    done
    say "Model-check manifest: ${model_check_manifest}"

    summary_macro="${BASE}/macros/SummarizeModelCheck.C"
    if [[ -s "$summary_macro" ]]; then
      summary_log="${SIM_DEST_BASE_RESOLVED}/checkModels_summary.log"
      root_arg="${summary_macro}(\"${model_check_manifest}\",\"${summary_log}\")"
      say "Running final model-check summary..."
      say "Summary log: ${summary_log}"
      if [[ -n "${RJ_ROOT_CMD:-}" ]]; then
        ${RJ_ROOT_CMD} -l -b -q "$root_arg" 2>&1 | tee "$summary_log"
      elif command -v root >/dev/null 2>&1; then
        root -l -b -q "$root_arg" 2>&1 | tee "$summary_log"
      else
        warn "ROOT is not on PATH; skipping ${summary_macro}"
      fi
      say "SFTP summary log from SSH: ${summary_log}"
    else
      warn "Summary macro not found: ${summary_macro}"
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
      sim_local_nfiles="${RJ_SIM_LOCAL_NFILES:-1}"
      [[ "$sim_local_nfiles" =~ ^[0-9]+$ && "$sim_local_nfiles" -ge 1 ]] || { err "RJ_SIM_LOCAL_NFILES must be a positive integer"; exit 2; }
      sim_local_grouped="${RJ_SIM_LOCAL_GROUPED:-0}"
      [[ "$sim_local_grouped" == "1" || "$sim_local_grouped" == "true" ]] && sim_local_grouped=1 || sim_local_grouped=0

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
          isSimEmbedded)          samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
          isSimEmbeddedInclusive) samples=( "run28_embeddedJet12" "run28_embeddedJet20" ) ;;
          isSimJet5)              samples=( "run28_jet5" ) ;;
          isSimMB)                samples=( "run28_detroit" ) ;;
          *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
        esac
      else
        samples=( "${SIM_SAMPLE}" )
      fi

      SIM_DEST_BASE_RESOLVED="${RJ_LOCAL_SIM_OUTPUT_BASE:-${BASE}/local_sim_outputs/${TAG}}"
      mkdir -p "${SIM_DEST_BASE_RESOLVED}"

      say "SIM local smoke test (mirrors condorDoAll matrix)"
      say "  YAML master : ${master_yaml}"
      say "  groupSize   : ${gs_local} (condorDoAll default)"
      say "  events      : ${nevt} (target for each local invocation)"
      say "  files/sample: ${sim_local_nfiles}"
      if (( sim_local_grouped )); then
        say "  local mode  : grouped list (one Fun4All invocation per sample/config)"
      else
        say "  local mode  : sequential one-file smoke jobs"
      fi
      say "  samples     : ${samples[*]}"
      say "  dest base   : ${SIM_DEST_BASE_RESOLVED}"
      echo

      for pt in "${sim_pts[@]}"; do
        for frac in "${sim_fracs[@]}"; do
          for vz in "${sim_vzs[@]}"; do
          for cone in "${sim_cones[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
          for uepipe in "${uepipe_modes[@]}"; do
          SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$pt")_$(sim_b2b_tag "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso_base_tags[$iso_idx]}"
          (( uepipe_in_tag )) && SIM_CFG_TAG="${SIM_CFG_TAG}_${uepipe}"
          SIM_CFG_TAG="${SIM_CFG_TAG}_${iso_selection_tags[$iso_idx]}"
          DEST_BASE="${SIM_DEST_BASE_RESOLVED}/${SIM_CFG_TAG}"
          yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt" "$frac" "$vz" "$cone" "${iso_sliding[$iso_idx]}" "${iso_fixed[$iso_idx]}" "${uepipe}" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}" "$SIM_CFG_TAG" "LOCAL")"

          for samp in "${samples[@]}"; do
            SIM_SAMPLE="$samp"
            GROUP_SIZE="$gs_local"

            sim_init
            n_available="$(wc -l < "$SIM_CLEAN_LIST" | tr -d ' ')"
            say "  [local] sample=${SIM_SAMPLE} (${n_available} entries) — using first ${sim_local_nfiles} input line(s)"
            mkdir -p "${DEST_BASE}/${SIM_SAMPLE}"
            find "${DEST_BASE}/${SIM_SAMPLE}" -type f -name "*_LOCAL_*.root" -delete 2>/dev/null || true

            if (( sim_local_grouped )); then
              tmp="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_LOCAL_files001-$(printf '%03d' "$sim_local_nfiles")_grp001.list"
              head -n "$sim_local_nfiles" "$SIM_CLEAN_LIST" > "$tmp"
              [[ -s "$tmp" ]] || { err "No sim entries (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"; exit 30; }

              chunk_base="$(basename "$tmp")"
              chunk_tag="${chunk_base%.list}"
              out_root_preview="${DEST_BASE}/${SIM_SAMPLE}/RecoilJets_${DATASET}_${chunk_tag}.root"
              n_group_lines="$(wc -l < "$tmp" | tr -d ' ')"
              first_input_line="$(head -n 1 "$tmp" 2>/dev/null || true)"
              last_input_line="$(tail -n 1 "$tmp" 2>/dev/null || true)"

              say "----------------------------------------"
              say "SIM local: tag=${SIM_CFG_TAG}  sample=${SIM_SAMPLE} groupedFiles=${n_group_lines}/${sim_local_nfiles}"
              say "  jet_pt_min=${pt}  back_to_back_pi_fraction=${frac}  vz_cut_cm=${vz}  coneR=${cone}  iso=${iso_tags[$iso_idx]}  uepipe=${uepipe}"
              say "  YAML override: ${yaml_override}"
              say "  grouped list : ${tmp}"
              say "  first input  : ${first_input_line}"
              say "  last input   : ${last_input_line}"
              say "  out ROOT      : ${out_root_preview}"
              say "  wrapper env   : RJ_VERBOSITY=${RJV} RJ_CONFIG_YAML=${yaml_override}"
              say "  wrapper args  : sample=${SIM_SAMPLE} dataset=${DATASET} mode=LOCAL nevents=${nevt} chunk=1 dest=${DEST_BASE}"
              say "Invoking wrapper locally…"

              RJ_VERBOSITY="$RJV" RJ_CONFIG_YAML="$yaml_override" bash "$EXE" "$SIM_SAMPLE" "$tmp" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
              echo
            else
              local_file_idx=0
              while IFS= read -r in_line; do
                [[ -n "$in_line" ]] || continue
                (( local_file_idx+=1 ))

                tmp="${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_LOCAL_file$(printf '%03d' "$local_file_idx")_grp001.list"
                printf "%s\n" "$in_line" > "$tmp"
                [[ -s "$tmp" ]] || { err "No sim entries (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"; exit 30; }

                chunk_base="$(basename "$tmp")"
                chunk_tag="${chunk_base%.list}"
                out_root_preview="${DEST_BASE}/${SIM_SAMPLE}/RecoilJets_${DATASET}_${chunk_tag}.root"

                say "----------------------------------------"
                say "SIM local: tag=${SIM_CFG_TAG}  sample=${SIM_SAMPLE} file=${local_file_idx}/${sim_local_nfiles}"
                say "  jet_pt_min=${pt}  back_to_back_pi_fraction=${frac}  vz_cut_cm=${vz}  coneR=${cone}  iso=${iso_tags[$iso_idx]}  uepipe=${uepipe}"
                say "  YAML override: ${yaml_override}"
                say "  grp001 list   : ${tmp}"
                say "  input line    : ${in_line}"
                say "  temp list     : ${tmp}"
                say "  out ROOT      : ${out_root_preview}"
                say "  wrapper env   : RJ_VERBOSITY=${RJV} RJ_CONFIG_YAML=${yaml_override}"
                say "  wrapper args  : sample=${SIM_SAMPLE} dataset=${DATASET} mode=LOCAL nevents=${nevt} chunk=${local_file_idx} dest=${DEST_BASE}"
                say "Invoking wrapper locally…"

                RJ_VERBOSITY="$RJV" RJ_CONFIG_YAML="$yaml_override" bash "$EXE" "$SIM_SAMPLE" "$tmp" "$DATASET" LOCAL "$nevt" "$local_file_idx" NONE "$DEST_BASE"
                echo
              done < <(head -n "$sim_local_nfiles" "$SIM_CLEAN_LIST")
            fi
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
      DATA_DEST_BASE_SAVED="${RJ_LOCAL_DATA_OUTPUT_BASE:-$DEST_BASE}"

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
        data_cfg_tag="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_base_tags[$iso_idx]}"
        (( uepipe_in_tag )) && data_cfg_tag="${data_cfg_tag}_${uepipe}"
        data_cfg_tag="${data_cfg_tag}_${iso_selection_tags[$iso_idx]}"
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
        pin_photon_id_scalars_in_yaml "$yaml_override" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}"
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
    data_cfg_tag="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_base_tags[$iso_idx]}_${selected_uepipe}_${iso_selection_tags[$iso_idx]}"

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
    pin_photon_id_scalars_in_yaml "$yaml_override" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}"

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
    read_replay_cones "$master_yaml"
    read_capture_cones "$master_yaml"
    sim_cones=( "${replay_cones[@]}" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
    (( ${#capture_cones[@]} )) || { err "No capture cone could be derived from $master_yaml"; exit 72; }
    build_iso_modes "$master_yaml"
    read_uepipe_modes "$master_yaml" "$TAG"

    SIM_DEST_BASE_RESOLVED="$DEST_BASE"

    pt0="${sim_pts[0]}"
    frac0="${sim_fracs[0]}"
    vz0="${sim_vzs[0]}"
    cone0="${sim_cones[0]}"
    SIM_CFG_TAG="jetMinPt$(sim_pt_tag "$pt0")_$(sim_b2b_tag "$frac0")_$(sim_vz_tag "$vz0")_$(sim_cone_tag "$cone0")_${iso_base_tags[0]}"
    (( uepipe_in_tag )) && SIM_CFG_TAG="${SIM_CFG_TAG}_${uepipe_modes[0]}"
    SIM_CFG_TAG="${SIM_CFG_TAG}_${iso_selection_tags[0]}"
    DEST_BASE="${SIM_DEST_BASE_RESOLVED}/${SIM_CFG_TAG}"
    stamp="$(date +%Y%m%d_%H%M%S)"
    yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt0" "$frac0" "$vz0" "$cone0" "${iso_sliding[0]}" "${iso_fixed[0]}" "${uepipe_modes[0]}" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}" "$SIM_CFG_TAG" "$stamp")"

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

  condorDoAllFromScratch|condorDoAllSmoke)
    [[ "$IS_SIM" -eq 1 ]] || { err "${ACTION} is currently implemented for SIM-style pool capture/replay DAGs."; exit 2; }
    sim_smoke=0
    [[ "$ACTION" == "condorDoAllSmoke" ]] && sim_smoke=1
    if (( sim_smoke )); then
      set -E
      trap 'rc=$?; err "SIM smoke DAG build failed (dataset=${DATASET:-unknown}, action=${ACTION:-unknown}, line=${LINENO}, rc=${rc}, command=${BASH_COMMAND}, dag_dir=${dag_dir:-unset})"; exit "$rc"' ERR
    fi
    if (( sim_smoke )) && [[ "${MAX_JOBS_EXPLICIT:-0}" -eq 0 ]]; then
      case "$DATASET" in
        isSimEmbedded)          MAX_JOBS="${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE_EMBEDDED:-${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE:-12}}" ;;
        isSimEmbeddedInclusive) MAX_JOBS="${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE_EMBEDDED_INCLUSIVE:-${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE_EMBEDDED:-${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE:-12}}}" ;;
        *)                      MAX_JOBS="${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE:-12}" ;;
      esac
    fi
    dag_dryrun=0
    if [[ "${RJ_DAG_DRYRUN:-0}" == "1" || "${RJ_DAG_DRYRUN:-0}" == "true" || "${RJ_DAG_DRYRUN:-0}" == "TRUE" ]]; then
      dag_dryrun=1
    fi

    master_yaml="$(sim_yaml_master_path)"
    [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

    if (( sim_smoke )) && [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
      case "$DATASET" in
        isSimEmbedded)          GROUP_SIZE="${RJ_SMOKE_GROUPSIZE_EMBEDDED:-2}" ;;
        isSimEmbeddedInclusive) GROUP_SIZE="${RJ_SMOKE_GROUPSIZE_EMBEDDED_INCLUSIVE:-${RJ_SMOKE_GROUPSIZE_EMBEDDED:-2}}" ;;
        *)                      GROUP_SIZE="${RJ_SMOKE_GROUPSIZE_SIM:-2}" ;;
      esac
    fi
    request_memory_was_explicit=0
    [[ -n "${RJ_REQUEST_MEMORY:-}" ]] && request_memory_was_explicit=1
    if (( sim_smoke )) && [[ -z "${RJ_REQUEST_MEMORY:-}" ]]; then
      case "$DATASET" in
        isSimEmbedded)          RJ_REQUEST_MEMORY="${RJ_SMOKE_REQUEST_MEMORY_EMBEDDED:-2500MB}" ;;
        isSimEmbeddedInclusive) RJ_REQUEST_MEMORY="${RJ_SMOKE_REQUEST_MEMORY_EMBEDDED_INCLUSIVE:-${RJ_SMOKE_REQUEST_MEMORY_EMBEDDED:-2500MB}}" ;;
        *)                      RJ_REQUEST_MEMORY="${RJ_SMOKE_REQUEST_MEMORY_SIM:-1800MB}" ;;
      esac
      export RJ_REQUEST_MEMORY
    fi
    request_memory="${RJ_REQUEST_MEMORY:-2000MB}"
    request_memory_digits="${request_memory//[^0-9]/}"
    request_memory_mb="${request_memory_digits:-0}"
    case "$request_memory" in
      *[Gg][Bb]|*[Gg]) request_memory_mb=$(( request_memory_mb * 1024 )) ;;
    esac
    capture_request_memory="${RJ_CAPTURE_REQUEST_MEMORY:-$request_memory}"
    replay_request_memory="${RJ_REPLAY_REQUEST_MEMORY:-}"
    if [[ -z "$replay_request_memory" ]]; then
      if (( sim_smoke )) && (( request_memory_was_explicit == 0 )); then
        case "$DATASET" in
          isSimEmbedded)          replay_request_memory="${RJ_SMOKE_REPLAY_REQUEST_MEMORY_EMBEDDED:-2500MB}" ;;
          isSimEmbeddedInclusive) replay_request_memory="${RJ_SMOKE_REPLAY_REQUEST_MEMORY_EMBEDDED_INCLUSIVE:-${RJ_SMOKE_REPLAY_REQUEST_MEMORY_EMBEDDED:-2500MB}}" ;;
          *)                      replay_request_memory="${RJ_SMOKE_REPLAY_REQUEST_MEMORY_SIM:-2200MB}" ;;
        esac
      else
        replay_request_memory="$request_memory"
      fi
    fi
    replay_roots_per_shard="${RJ_REPLAY_OUTPUT_ROOTS_PER_SHARD:-${RJ_POOL_REPLAY_OUTPUT_ROOTS_PER_SHARD:-999999}}"
    [[ "$replay_roots_per_shard" =~ ^[0-9]+$ && "$replay_roots_per_shard" -gt 0 ]] || replay_roots_per_shard=999999
    replay_max_open_outputs="${RJ_REPLAY_MAX_OPEN_OUTPUTS:-$replay_roots_per_shard}"
    capture_request_memory_mb="$(memory_request_to_mb "$capture_request_memory")"
    replay_request_memory_mb="$(memory_request_to_mb "$replay_request_memory")"
    profile_job="${RJ_PROFILE_JOB:-0}"
    (( sim_smoke )) && profile_job=1
    sim_capture_nevents="${RJ_SMOKE_SIM_NEVENTS:-3000}"
    job_heartbeat_seconds="${RJ_JOB_HEARTBEAT_SECONDS:-0}"
    if (( sim_smoke )) && [[ "$job_heartbeat_seconds" == "0" ]]; then
      job_heartbeat_seconds="${RJ_SMOKE_JOB_HEARTBEAT_SECONDS:-120}"
    fi
    capture_submit_extra="$(condor_pool_submit_extra_lines capture "$sim_smoke")"
    replay_submit_extra="$(condor_pool_submit_extra_lines replay "$sim_smoke")"

    gs_capture="$GROUP_SIZE"
    gs_replay="${RJ_POOL_REPLAY_GROUP_SIZE:-20}"
    if (( sim_smoke )) && [[ -z "${RJ_POOL_REPLAY_GROUP_SIZE:-}" ]]; then
      case "$DATASET" in
        isSimEmbedded)          gs_replay="${RJ_SMOKE_REPLAY_GROUPSIZE_EMBEDDED:-5}" ;;
        isSimEmbeddedInclusive) gs_replay="${RJ_SMOKE_REPLAY_GROUPSIZE_EMBEDDED_INCLUSIVE:-${RJ_SMOKE_REPLAY_GROUPSIZE_EMBEDDED:-5}}" ;;
        *)                      gs_replay="${RJ_SMOKE_REPLAY_GROUPSIZE_SIM:-5}" ;;
      esac
    fi
    if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 1 ]]; then
      gs_replay="$GROUP_SIZE"
    fi
    [[ "$gs_replay" =~ ^[0-9]+$ && "$gs_replay" -gt 0 ]] || gs_replay=20

    mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
    mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
    mapfile -t sim_vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
    read_replay_cones "$master_yaml"
    read_capture_cones "$master_yaml"
    sim_cones=( "${replay_cones[@]}" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
    (( ${#capture_cones[@]} )) || { err "No capture cone could be derived from $master_yaml"; exit 72; }
    build_iso_modes "$master_yaml"
    read_uepipe_modes "$master_yaml" "$TAG"

    samples=()
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
      case "$DATASET" in
        isSimEmbedded)          samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
        isSimEmbeddedInclusive) samples=( "run28_embeddedJet12" "run28_embeddedJet20" ) ;;
        isSimJet5)              samples=( "run28_jet5" ) ;;
        isSimMB)                samples=( "run28_detroit" ) ;;
        *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
      esac
    else
      samples=( "${SIM_SAMPLE}" )
    fi

    (( dag_dryrun )) || need_cmd condor_submit_dag

    workflow_stamp="$(date +%Y%m%d_%H%M%S)"
    if (( sim_smoke )); then
      pool_base="${RJ_POOL_OUTPUT_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaPoolsSmoke/${TAG}_smokeTest_${workflow_stamp}}"
      out_base="${RJ_SMOKE_OUTPUT_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaSmoke/${TAG}_smokeTest_${workflow_stamp}}"
    else
      pool_base="${RJ_POOL_OUTPUT_BASE:-${POOL_DEST_ROOT}/${TAG}}"
      out_base="${DEST_BASE%/}"
    fi
    dag_dir="${SUB_DIR}/pool_workflow_${TAG}_${workflow_stamp}"
    mkdir -p "$dag_dir" "$SIM_YAML_OVERRIDE_DIR"
    dag="${dag_dir}/RecoilJets_pool_${DATASET}_${workflow_stamp}.dag"
    manifest="${dag_dir}/manifest.txt"
    : > "$dag"

    cleanup_bulk_snapshots_for_tag
    case "$DATASET" in
      isSimEmbedded|isSimEmbeddedInclusive) create_pipeline_snapshot "auau" "$workflow_stamp" ;;
      *)                                    create_pipeline_snapshot "pp" "$workflow_stamp" ;;
    esac

    pool_macro="${BASE}/macros/Fun4All_recoilJets_poolReplay.C"
    [[ -f "$pool_macro" ]] || { err "Pool replay macro not found: $pool_macro"; exit 92; }

    {
      echo "workflow=${ACTION}"
      echo "dataset=${DATASET}"
      echo "smoke=${sim_smoke}"
      echo "yaml=${master_yaml}"
      echo "pool_base=${pool_base}"
      echo "output_base=${out_base}"
      echo "capture_group_size=${gs_capture}"
      echo "group_size=${gs_capture}"
      echo "replay_group_size=${gs_replay}"
      echo "smoke_capture_job_cap=${MAX_JOBS}"
      echo "request_memory=${request_memory}"
      echo "request_memory_mb=${request_memory_mb}"
      echo "capture_request_memory=${capture_request_memory}"
      echo "replay_request_memory=${replay_request_memory}"
      echo "capture_request_memory_mb=${capture_request_memory_mb}"
      echo "replay_request_memory_mb=${replay_request_memory_mb}"
      echo "replay_output_roots_per_shard=${replay_roots_per_shard}"
      echo "replay_max_open_outputs=${replay_max_open_outputs}"
      echo "capture_nevents=${sim_capture_nevents}"
      echo "profile_job=${profile_job}"
      echo "pool_profile_counts=1"
      echo "job_heartbeat_seconds=${job_heartbeat_seconds}"
      echo "condor_use_late_materialization=${RJ_CONDOR_USE_LATE_MATERIALIZATION:-1}"
      echo "condor_capture_max_idle=${RJ_CONDOR_CAPTURE_MAX_IDLE:-${RJ_CONDOR_MAX_IDLE:-2000}}"
      echo "condor_capture_max_materialize=${RJ_CONDOR_CAPTURE_MAX_MATERIALIZE:-${RJ_CONDOR_MAX_MATERIALIZE:-5000}}"
      echo "condor_replay_max_idle=${RJ_CONDOR_REPLAY_MAX_IDLE:-${RJ_CONDOR_MAX_IDLE:-1000}}"
      echo "condor_replay_max_materialize=${RJ_CONDOR_REPLAY_MAX_MATERIALIZE:-${RJ_CONDOR_MAX_MATERIALIZE:-3000}}"
      echo "samples=${samples[*]}"
      echo "replay_cones=${sim_cones[*]}"
      echo "capture_cones=${capture_cones[*]}"
      echo "uepipe_modes=${uepipe_modes[*]}"
      echo "photon_id_sets_expanded=${#iso_selection_tags[@]}"
      echo "submitted_at=${workflow_stamp}"
    } > "$manifest"

    declare -A cap_node_by_key=()
    declare -A pool_list_by_key=()
    capture_count=0
    replay_count=0
    capture_job_total=0
    replay_job_total=0

    say "${BOLD}${ACTION} DAG requested${RST}"
    say "  dataset       : ${DATASET}"
    say "  pool base     : ${pool_base}"
    say "  output base   : ${out_base}"
    (( sim_smoke )) && say "  smoke outputs : isolated thesisAnaSmoke/thesisAnaPoolsSmoke directories"
    say "  memory        : capture=${capture_request_memory} replay=${replay_request_memory}"
    (( sim_smoke )) && say "  smoke cap     : max ${MAX_JOBS} capture jobs/sample, nevents=${sim_capture_nevents}/worker, heartbeat=${job_heartbeat_seconds}s"
    say "  replay split  : poolFiles/job=${gs_replay} outputRoots/shard=${replay_roots_per_shard} maxOpenOutputs=${replay_max_open_outputs}"
    say "  profiling     : ${profile_job}"
    say "  capture axes  : storedIsolation=[${capture_cones[*]}], clusterUE=[${uepipe_modes[*]}]"
    say "  replay axes   : pT=[${sim_pts[*]}], dphi=[${sim_fracs[*]}], vz=[${sim_vzs[*]}], ID rows=${#iso_selection_tags[@]}"
    say "  DAG dir       : ${dag_dir}"

    for cone in "${capture_cones[@]}"; do
      for uepipe in "${uepipe_modes[@]}"; do
        POOL_CAPTURE_TAG="$(pool_capture_cfg_tag "$cone" "$uepipe" "$master_yaml")"
        (( sim_smoke )) && say "[simSmoke] capture planning begin: cap_tag=${POOL_CAPTURE_TAG} cone=${cone} clusterUE=${uepipe}"
        cap_yaml="$(sim_make_yaml_override "$master_yaml" "${sim_pts[0]}" "${sim_fracs[0]}" "${sim_vzs[0]}" "$cone" "${iso_sliding[0]}" "${iso_fixed[0]}" "$uepipe" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}" "$POOL_CAPTURE_TAG" "$workflow_stamp" "false")"

        for samp in "${samples[@]}"; do
          SIM_SAMPLE="$samp"
          SIM_CFG_TAG="$POOL_CAPTURE_TAG"
          DEST_BASE="${pool_base}/${POOL_CAPTURE_TAG}"
          GROUP_SIZE="$gs_capture"
          SIM_STAGE_NAMESPACE="${TAG}_${workflow_stamp}_capture_${POOL_CAPTURE_TAG}"
          export SIM_STAGE_NAMESPACE
          mapfile -t groups < <( make_sim_groups "$GROUP_SIZE" )
          (( ${#groups[@]} )) || { err "No sim capture groups produced (sample=${SIM_SAMPLE}, tag=${POOL_CAPTURE_TAG})"; exit 30; }
          (( sim_smoke )) && say "[simSmoke] capture groups ready: cap_tag=${POOL_CAPTURE_TAG} sample=${SIM_SAMPLE} groups=${#groups[@]} groupSize=${GROUP_SIZE}"
          if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 && "${#groups[@]}" -gt "$MAX_JOBS" ]]; then
            say "Capping pool capture group list for sample=${SIM_SAMPLE}, tag=${POOL_CAPTURE_TAG}: ${#groups[@]} → ${MAX_JOBS} jobs"
            groups=( "${groups[@]:0:$MAX_JOBS}" )
            (( sim_smoke )) && say "[simSmoke] capture cap applied: cap_tag=${POOL_CAPTURE_TAG} sample=${SIM_SAMPLE} cappedGroups=${#groups[@]}"
          fi

          pool_out_dir="${pool_base}/${POOL_CAPTURE_TAG}/${SIM_SAMPLE}"
          mkdir -p "$pool_out_dir"
          rm -f "${pool_out_dir}/"*.root 2>/dev/null || true

          expected_pool_list="${dag_dir}/pool_expected_${POOL_CAPTURE_TAG}_${SIM_SAMPLE}.list"
          : > "$expected_pool_list"
          analysis_tag="$(sim_analysis_tag_for_dataset "$DATASET")"
          for glist in "${groups[@]}"; do
            chunk_base="$(basename "$glist")"
            chunk_tag="${chunk_base%.list}"
            printf '%s\n' "${pool_out_dir}/RecoilJets_${analysis_tag}_${chunk_tag}.root" >> "$expected_pool_list"
          done
          (( sim_smoke )) && say "[simSmoke] expected pool list written: cap_tag=${POOL_CAPTURE_TAG} sample=${SIM_SAMPLE} file=${expected_pool_list}"

          cap_sub="${dag_dir}/capture_${POOL_CAPTURE_TAG}_${SIM_SAMPLE}.sub"
          cap_args="${dag_dir}/capture_${POOL_CAPTURE_TAG}_${SIM_SAMPLE}.args"
          : > "$cap_args"
          exe_for_sub="${BULK_FROZEN_EXE:-${EXE}}"
          macro_env_for_sub=""
          [[ -n "${BULK_FROZEN_MACRO:-}" ]] && macro_env_for_sub=";RJ_MACRO_PATH=${BULK_FROZEN_MACRO}"
          cap_prefix="capture_${workflow_stamp}_${SIM_SAMPLE}_${POOL_CAPTURE_TAG}"
          cat > "$cap_sub" <<SUB
universe      = vanilla
executable    = ${exe_for_sub}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${cap_prefix}.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${cap_prefix}.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${cap_prefix}.\$(Cluster).\$(Process).err
request_memory= ${capture_request_memory}
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
environment   = RJ_VERBOSITY=0;RJ_CONFIG_YAML=${cap_yaml}${macro_env_for_sub};RJ_POOL_MODE=capture;RJ_PROFILE_JOB=${profile_job};RJ_POOL_PROFILE_COUNTS=1;RJ_JOB_HEARTBEAT_SECONDS=${job_heartbeat_seconds};RJ_PROFILE_STAGE=capture;RJ_PROFILE_LABEL=${POOL_CAPTURE_TAG}_${SIM_SAMPLE};RJ_REQUEST_MEMORY_MB=${capture_request_memory_mb}
${capture_submit_extra}
queue arguments from ${cap_args}
SUB
          gidx=0
          for glist in "${groups[@]}"; do
            (( gidx+=1 ))
            printf '%s %s %s $(Cluster) %s %d NONE %s\n' \
                   "$SIM_SAMPLE" "$glist" "$DATASET" "$sim_capture_nevents" "$gidx" "${pool_base}/${POOL_CAPTURE_TAG}" >> "$cap_args"
          done

          cap_node="CAP_$(sanitize_node_name "${POOL_CAPTURE_TAG}_${SIM_SAMPLE}")"
          printf 'JOB %s %s\n' "$cap_node" "$cap_sub" >> "$dag"
          cap_node_by_key["${POOL_CAPTURE_TAG}|${SIM_SAMPLE}"]="$cap_node"
          pool_list_by_key["${POOL_CAPTURE_TAG}|${SIM_SAMPLE}"]="$expected_pool_list"
          (( capture_count+=1 ))
          (( capture_job_total+=${#groups[@]} ))
          {
            echo "capture_node=${cap_node} tag=${POOL_CAPTURE_TAG} sample=${SIM_SAMPLE} jobs=${#groups[@]} yaml=${cap_yaml} output=${pool_out_dir}"
            echo "profile_glob=${OUT_DIR}/${cap_prefix}.*.out"
          } >> "$manifest"
          (( sim_smoke )) && say "[simSmoke] capture node ready: node=${cap_node} sample=${SIM_SAMPLE} jobs=${#groups[@]} sub=${cap_sub}"
        done
      done
    done

    for cone in "${capture_cones[@]}"; do
      for uepipe in "${uepipe_modes[@]}"; do
        POOL_CAPTURE_TAG="$(pool_capture_cfg_tag "$cone" "$uepipe" "$master_yaml")"
        fanout_dirs="${SIM_YAML_OVERRIDE_DIR}/pool_replay_fanout_${TAG}_${POOL_CAPTURE_TAG}_${workflow_stamp}.txt"
        (( sim_smoke )) && say "[simSmoke] replay fanout begin: cap_tag=${POOL_CAPTURE_TAG} outputBase=${out_base}"
        emit_pool_replay_fanout_dirs_file "$fanout_dirs" "$out_base" "$master_yaml" "$workflow_stamp" "${sim_cones[*]}" "$uepipe" "${sim_pts[*]}" "${sim_fracs[*]}" "${sim_vzs[*]}"
        fanout_count="$(wc -l < "$fanout_dirs" | awk '{print $1}')"
        (( fanout_count > 0 )) || { err "Pool replay fanout is empty for ${POOL_CAPTURE_TAG}: ${fanout_dirs}"; exit 95; }
        (( sim_smoke )) && say "[simSmoke] replay fanout written: cap_tag=${POOL_CAPTURE_TAG} fanoutOutputs=${fanout_count} file=${fanout_dirs}"
        append_fanout_cfgs_to_manifest "$fanout_dirs" "$manifest"

        (( sim_smoke )) && say "[simSmoke] cleaning fanout output dirs begin: cap_tag=${POOL_CAPTURE_TAG} rows=${fanout_count}"
        clean_fanout_output_dirs_from_file "$fanout_dirs" "${POOL_CAPTURE_TAG}" 1 25 "${samples[@]}"
        (( sim_smoke )) && say "[simSmoke] cleaning fanout output dirs done: cap_tag=${POOL_CAPTURE_TAG} uniqueRoots=${CLEAN_FANOUT_LAST_COUNT:-0}"

        fanout_shards=()
        mapfile -t fanout_shards < <(split_pool_replay_fanout_shards "$fanout_dirs" "${dag_dir}/fanoutShards/${POOL_CAPTURE_TAG}" "$POOL_CAPTURE_TAG" "$replay_roots_per_shard" "$sim_smoke")
        (( ${#fanout_shards[@]} )) || { err "No replay fanout shards were produced for ${POOL_CAPTURE_TAG}"; exit 95; }

        for samp in "${samples[@]}"; do
          cap_key="${POOL_CAPTURE_TAG}|${samp}"
          cap_node="${cap_node_by_key[$cap_key]:-}"
          pool_all="${pool_list_by_key[$cap_key]:-}"
          [[ -n "$cap_node" && -s "$pool_all" ]] || { err "Internal DAG build error: missing capture node/pool list for ${cap_key}"; exit 97; }

          pool_stage_dir="${dag_dir}/poolReplayLists/${POOL_CAPTURE_TAG}/${samp}"
          mkdir -p "$pool_stage_dir"
          rm -f "${pool_stage_dir}/pool_grp"*.list "${pool_stage_dir}/pool_grp_raw_"* 2>/dev/null || true
          split -l "$gs_replay" -d -a 5 "$pool_all" "${pool_stage_dir}/pool_grp_raw_"
          groups=()
          gidx=0
          for raw in "${pool_stage_dir}/pool_grp_raw_"*; do
            [[ -s "$raw" ]] || { rm -f "$raw"; continue; }
            (( gidx+=1 ))
            out_list="${pool_stage_dir}/pool_grp$(printf "%03d" "$gidx").list"
            mv "$raw" "$out_list"
            groups+=( "$out_list" )
          done
          (( ${#groups[@]} )) || { err "No replay pool groups produced for ${POOL_CAPTURE_TAG}/${samp}"; exit 95; }

          shard_i=0
          for fanout_shard in "${fanout_shards[@]}"; do
          (( shard_i+=1 ))
          shard_rows="$(wc -l < "$fanout_shard" | awk '{print $1}')"
          shard_roots="$(awk -F'|' 'NF && $1 !~ /^#/ && $1 != "" {seen[$1]=1} END{for(k in seen)c++; print c+0}' "$fanout_shard")"
          replay_sub="${dag_dir}/replay_${POOL_CAPTURE_TAG}_${samp}_shard$(printf "%03d" "$shard_i").sub"
          replay_args="${dag_dir}/replay_${POOL_CAPTURE_TAG}_${samp}_shard$(printf "%03d" "$shard_i").args"
          : > "$replay_args"
          replay_prefix="replay_${workflow_stamp}_${samp}_${POOL_CAPTURE_TAG}_shard$(printf "%03d" "$shard_i")"
          cat > "$replay_sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${replay_prefix}.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${replay_prefix}.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${replay_prefix}.\$(Cluster).\$(Process).err
request_memory= ${replay_request_memory}
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
environment   = RJ_VERBOSITY=0;RJ_CONFIG_YAML=${master_yaml};RJ_MACRO_PATH=${pool_macro};RJ_ID_FANOUT_DIRS_FILE=${fanout_shard};RJ_REPLAY_MAX_OPEN_OUTPUTS=${replay_max_open_outputs};RJ_REPLAY_OUTPUT_ROOTS_PER_SHARD=${replay_roots_per_shard};RJ_PROFILE_JOB=${profile_job};RJ_POOL_PROFILE_COUNTS=1;RJ_JOB_HEARTBEAT_SECONDS=${job_heartbeat_seconds};RJ_PROFILE_STAGE=replay;RJ_PROFILE_LABEL=${POOL_CAPTURE_TAG}_${samp}_shard$(printf "%03d" "$shard_i");RJ_REQUEST_MEMORY_MB=${replay_request_memory_mb}
${replay_submit_extra}
queue arguments from ${replay_args}
SUB
          g=0
          for glist in "${groups[@]}"; do
            (( g+=1 ))
            printf '%s %s %s $(Cluster) 0 %d NONE %s\n' \
                   "$samp" "$glist" "$DATASET" "$g" "${out_base}/${POOL_CAPTURE_TAG}" >> "$replay_args"
          done
          replay_node="REP_$(sanitize_node_name "${POOL_CAPTURE_TAG}_${samp}_shard$(printf "%03d" "$shard_i")")"
          printf 'JOB %s %s\n' "$replay_node" "$replay_sub" >> "$dag"
          printf 'PARENT %s CHILD %s\n' "$cap_node" "$replay_node" >> "$dag"
          (( replay_count+=1 ))
          (( replay_job_total+=${#groups[@]} ))
          echo "replay_node=${replay_node} cap_tag=${POOL_CAPTURE_TAG} sample=${samp} shard=${shard_i}/${#fanout_shards[@]} jobs=${#groups[@]} fanout_outputs=${shard_rows} fanout_roots=${shard_roots} parent=${cap_node} fanout=${fanout_shard} parent_fanout=${fanout_dirs}" >> "$manifest"
          echo "profile_glob=${OUT_DIR}/${replay_prefix}.*.out" >> "$manifest"
          (( sim_smoke )) && say "[simSmoke] replay node ready: node=${replay_node} sample=${samp} jobs=${#groups[@]} parent=${cap_node} shard=${shard_i}/${#fanout_shards[@]} fanoutRows=${shard_rows} fanoutRoots=${shard_roots} sub=${replay_sub}"
          done
        done
      done
    done

    if (( sim_smoke )) && (( capture_count == 0 || replay_count == 0 )); then
      err "condorDoAllSmoke produced capture_count=${capture_count}, replay_count=${replay_count}; refusing to submit an empty SIM smoke DAG"
      exit 98
    fi

    sim_workflow_name="$([[ "$sim_smoke" -eq 1 ]] && echo "simPoolSmoke_${TAG}_${workflow_stamp}" || echo "poolFromScratch_${TAG}_${workflow_stamp}")"
    final_sub="$(write_dag_final_summary_files "$dag_dir" "$sim_workflow_name" "$DATASET" "$pool_base" "$out_base" "$manifest" "$dag")"
    printf 'FINAL FINAL_SUMMARY %s\n' "$final_sub" >> "$dag"

    say "DAG build summary:"
    say "  capture nodes : ${capture_count}"
    say "  replay nodes  : ${replay_count}"
    say "  capture jobs  : ${capture_job_total}"
    say "  replay jobs   : ${replay_job_total}"
    say "  manifest      : ${manifest}"
    say "  dag           : ${dag}"
    {
      echo "capture_worker_jobs=${capture_job_total}"
      echo "replay_worker_jobs=${replay_job_total}"
      echo "total_worker_jobs=$(( capture_job_total + replay_job_total ))"
      echo "dag_max_worker_jobs=${RJ_DAG_MAX_WORKER_JOBS:-50000}"
      echo "dag_warn_worker_jobs=${RJ_DAG_WARN_WORKER_JOBS:-40000}"
    } >> "$manifest"
    sim_min_worker_jobs=0
    if (( ! sim_smoke )); then
      sim_min_worker_jobs="${RJ_DAG_MIN_WORKER_JOBS:-${RJ_DAG_MIN_WORKER_JOBS_SIM:-15000}}"
    fi
    echo "dag_min_worker_jobs=${sim_min_worker_jobs}" >> "$manifest"
    check_dag_worker_job_budget "SIM pool workflow ${TAG}" "$capture_job_total" "$replay_job_total" "$dag_dryrun" "$sim_min_worker_jobs"
    publish_production_manifest "$manifest" "$out_base"
    if (( dag_dryrun && sim_smoke )); then
      echo "RECOILJETS_SMOKETEST_DRYRUN_V1"
      echo "dataset=${DATASET}"
      echo "mode=${ACTION}"
      echo "runs_or_samples=${samples[*]}"
      echo "manifest=${manifest}"
      echo "dag=${dag}"
    fi
    submit_dag_with_notify "$dag"
    (( sim_smoke )) && trap - ERR
    ;;

  condorHistFromPool)
    [[ "$IS_SIM" -eq 1 ]] || { err "condorHistFromPool is currently defined for SIM-style staged pool replay only."; exit 2; }

    master_yaml="$(sim_yaml_master_path)"
    [[ -s "$master_yaml" ]] || { err "Master YAML not found or empty: $master_yaml"; exit 72; }

    gs_doall="$GROUP_SIZE"
    if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
      gs_doall="20"
    fi
    replay_request_memory="${RJ_REPLAY_REQUEST_MEMORY:-}"
    if [[ -z "$replay_request_memory" ]]; then
      case "$DATASET" in
        isSimEmbedded|isSimEmbeddedInclusive) replay_request_memory="${RJ_DEFAULT_REPLAY_REQUEST_MEMORY_EMBEDDED:-3000MB}" ;;
        *)                                   replay_request_memory="${RJ_DEFAULT_REPLAY_REQUEST_MEMORY_SIM:-2500MB}" ;;
      esac
    fi
    replay_request_memory_mb="$(memory_request_to_mb "$replay_request_memory")"
    replay_roots_per_shard="${RJ_REPLAY_OUTPUT_ROOTS_PER_SHARD:-${RJ_POOL_REPLAY_OUTPUT_ROOTS_PER_SHARD:-999999}}"
    [[ "$replay_roots_per_shard" =~ ^[0-9]+$ && "$replay_roots_per_shard" -gt 0 ]] || replay_roots_per_shard=999999
    replay_max_open_outputs="${RJ_REPLAY_MAX_OPEN_OUTPUTS:-$replay_roots_per_shard}"
    profile_job="${RJ_PROFILE_JOB:-0}"

    mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
    mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
    mapfile -t sim_vzs   < <( yaml_get_values "vz_cut_cm" "$master_yaml" )
    read_replay_cones "$master_yaml"
    read_capture_cones "$master_yaml"
    sim_cones=( "${replay_cones[@]}" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
    (( ${#capture_cones[@]} )) || { err "No capture cone could be derived from $master_yaml"; exit 72; }
    build_iso_modes "$master_yaml"
    read_uepipe_modes "$master_yaml" "$TAG"

    samples=()
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
      case "$DATASET" in
        isSimEmbedded)          samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
        isSimEmbeddedInclusive) samples=( "run28_embeddedJet12" "run28_embeddedJet20" ) ;;
        isSimJet5)              samples=( "run28_jet5" ) ;;
        isSimMB)                samples=( "run28_detroit" ) ;;
        *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
      esac
    else
      samples=( "${SIM_SAMPLE}" )
    fi

    pool_base="${RJ_POOL_INPUT_BASE:-${POOL_DEST_ROOT}/${TAG}}"

    say "${BOLD}condorHistFromPool requested${RST}"
    say "  dataset        : ${DATASET}"
    say "  YAML           : ${master_yaml}"
    say "  pool base      : ${pool_base}"
    say "  output base    : ${DEST_BASE}"
    say "  groupSize      : ${gs_doall}"
    say "  request mem    : replay=${replay_request_memory}"
    say "  replay split   : outputRoots/shard=${replay_roots_per_shard} maxOpenOutputs=${replay_max_open_outputs}"
    say "  photon_id_sets : ${#iso_selection_tags[@]} ID rows after expansion"
    say "  iso groups     : $(iso_group_count)"
    say "  clusterUE modes: [${uepipe_modes[*]}]"

    need_cmd condor_submit_dag
    replay_stamp="$(date +%Y%m%d_%H%M%S)"
    dag_dir="${SUB_DIR}/pool_replay_${TAG}_${replay_stamp}"
    dag="${dag_dir}/RecoilJets_poolReplay_${DATASET}_${replay_stamp}.dag"
    manifest="${dag_dir}/manifest.txt"
    mkdir -p "$dag_dir"
    : > "$dag"
    pool_macro="${BASE}/macros/Fun4All_recoilJets_poolReplay.C"
    [[ -f "$pool_macro" ]] || { err "Pool replay macro not found: $pool_macro"; exit 92; }
    {
      echo "workflow=condorHistFromPool"
      echo "dataset=${DATASET}"
      echo "yaml=${master_yaml}"
      echo "pool_base=${pool_base}"
      echo "output_base=${DEST_BASE}"
      echo "group_size=${gs_doall}"
      echo "replay_group_size=${gs_doall}"
      echo "replay_request_memory=${replay_request_memory}"
      echo "replay_request_memory_mb=${replay_request_memory_mb}"
      echo "replay_output_roots_per_shard=${replay_roots_per_shard}"
      echo "replay_max_open_outputs=${replay_max_open_outputs}"
      echo "profile_job=${profile_job}"
      echo "samples=${samples[*]}"
      echo "submitted_at=${replay_stamp}"
    } > "$manifest"
    replay_count=0
    replay_job_total=0

    DEST_BASE="${DEST_BASE%/}"
    for cone in "${capture_cones[@]}"; do
      for uepipe in "${uepipe_modes[@]}"; do
        POOL_CAPTURE_TAG="$(pool_capture_cfg_tag "$cone" "$uepipe" "$master_yaml")"
        fanout_dirs="${SIM_YAML_OVERRIDE_DIR}/pool_replay_fanout_${TAG}_${POOL_CAPTURE_TAG}_${replay_stamp}.txt"
        emit_pool_replay_fanout_dirs_file "$fanout_dirs" "$DEST_BASE" "$master_yaml" "$replay_stamp" "${sim_cones[*]}" "$uepipe" "${sim_pts[*]}" "${sim_fracs[*]}" "${sim_vzs[*]}"
        fanout_count="$(wc -l < "$fanout_dirs" | awk '{print $1}')"
        append_fanout_cfgs_to_manifest "$fanout_dirs" "$manifest"

        clean_fanout_output_dirs_from_file "$fanout_dirs" "${POOL_CAPTURE_TAG}" 1 25 "${samples[@]}"
        fanout_shards=()
        mapfile -t fanout_shards < <(split_pool_replay_fanout_shards "$fanout_dirs" "${dag_dir}/fanoutShards/${POOL_CAPTURE_TAG}" "$POOL_CAPTURE_TAG" "$replay_roots_per_shard" 1)
        (( ${#fanout_shards[@]} )) || { err "No replay fanout shards were produced for ${POOL_CAPTURE_TAG}"; exit 95; }

        for samp in "${samples[@]}"; do
          pool_dir="${pool_base}/${POOL_CAPTURE_TAG}/${samp}"
          [[ -d "$pool_dir" ]] || { err "Missing pool directory for ${POOL_CAPTURE_TAG}/${samp}: $pool_dir"; exit 93; }
          pool_stage_dir="${STAGE_DIR}/poolReplay_${replay_stamp}/${POOL_CAPTURE_TAG}/${samp}"
          mkdir -p "$pool_stage_dir"
          pool_all="${pool_stage_dir}/pool_all.list"
          find "$pool_dir" -maxdepth 1 -type f -name '*.root' | sort > "$pool_all"
          [[ -s "$pool_all" ]] || { err "No pool ROOT files found in $pool_dir"; exit 94; }

          rm -f "${pool_stage_dir}/pool_grp"*.list 2>/dev/null || true
          split -l "$gs_doall" -d -a 5 "$pool_all" "${pool_stage_dir}/pool_grp_raw_"
          groups=()
          gidx=0
          for raw in "${pool_stage_dir}/pool_grp_raw_"*; do
            [[ -s "$raw" ]] || { rm -f "$raw"; continue; }
            (( gidx+=1 ))
            out_list="${pool_stage_dir}/pool_grp$(printf "%03d" "$gidx").list"
            mv "$raw" "$out_list"
            groups+=( "$out_list" )
          done
          (( ${#groups[@]} )) || { err "No pool replay groups produced for ${POOL_CAPTURE_TAG}/${samp}"; exit 95; }
          if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 && "${#groups[@]}" -gt "$MAX_JOBS" ]]; then
            groups=( "${groups[@]:0:$MAX_JOBS}" )
          fi

          SIM_SAMPLE="$samp"
          shard_i=0
          for fanout_shard in "${fanout_shards[@]}"; do
          (( shard_i+=1 ))
          shard_rows="$(wc -l < "$fanout_shard" | awk '{print $1}')"
          shard_roots="$(awk -F'|' 'NF && $1 !~ /^#/ && $1 != "" {seen[$1]=1} END{for(k in seen)c++; print c+0}' "$fanout_shard")"
          SIM_JOB_PREFIX="poolReplay_${samp}_${POOL_CAPTURE_TAG}_shard$(printf "%03d" "$shard_i")"
          sub="${dag_dir}/RecoilJets_poolReplay_${POOL_CAPTURE_TAG}_${samp}_shard$(printf "%03d" "$shard_i").sub"
          args_file="${dag_dir}/RecoilJets_poolReplay_${POOL_CAPTURE_TAG}_${samp}_shard$(printf "%03d" "$shard_i").args"
          : > "$args_file"

          cat > "$sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).err
request_memory= ${replay_request_memory}
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
environment   = RJ_VERBOSITY=0;RJ_CONFIG_YAML=${master_yaml};RJ_MACRO_PATH=${pool_macro};RJ_ID_FANOUT_DIRS_FILE=${fanout_shard};RJ_REPLAY_MAX_OPEN_OUTPUTS=${replay_max_open_outputs};RJ_REPLAY_OUTPUT_ROOTS_PER_SHARD=${replay_roots_per_shard};RJ_PROFILE_JOB=${profile_job};RJ_PROFILE_STAGE=replay;RJ_PROFILE_LABEL=${POOL_CAPTURE_TAG}_${samp}_shard$(printf "%03d" "$shard_i");RJ_REQUEST_MEMORY_MB=${replay_request_memory_mb}
queue arguments from ${args_file}
SUB
          g=0
          for glist in "${groups[@]}"; do
            (( g+=1 ))
            printf '%s %s %s $(Cluster) 0 %d NONE %s\n' \
                   "$samp" "$glist" "$DATASET" "$g" "${DEST_BASE}/${POOL_CAPTURE_TAG}" >> "$args_file"
          done

          say "Submitting pool replay (pool tag=${POOL_CAPTURE_TAG}, sample=${samp}) → jobs=${BOLD}${#groups[@]}${RST}"
          say "  shard       : ${shard_i}/${#fanout_shards[@]} (${shard_rows} rows, ${shard_roots} output roots)"
          say "  pool dir    : ${pool_dir}"
          replay_node="REP_$(sanitize_node_name "${POOL_CAPTURE_TAG}_${samp}_shard$(printf "%03d" "$shard_i")")"
          printf 'JOB %s %s\n' "$replay_node" "$sub" >> "$dag"
          (( replay_count+=1 ))
          (( replay_job_total+=${#groups[@]} ))
          echo "replay_node=${replay_node} cap_tag=${POOL_CAPTURE_TAG} sample=${samp} shard=${shard_i}/${#fanout_shards[@]} jobs=${#groups[@]} fanout_outputs=${shard_rows} fanout_roots=${shard_roots} pool=${pool_dir} fanout=${fanout_shard} parent_fanout=${fanout_dirs}" >> "$manifest"
          echo "profile_glob=${OUT_DIR}/${SIM_JOB_PREFIX}.job.*.out" >> "$manifest"
          done
        done
      done
    done
    (( replay_count > 0 )) || { err "No pool replay DAG nodes were produced"; exit 96; }
    final_sub="$(write_dag_final_summary_files "$dag_dir" "poolReplay_${TAG}_${replay_stamp}" "$DATASET" "$pool_base" "$DEST_BASE" "$manifest" "$dag")"
    printf 'FINAL FINAL_SUMMARY %s\n' "$final_sub" >> "$dag"
    say "Pool replay DAG build summary:"
    say "  replay nodes : ${replay_count}"
    say "  replay jobs  : ${replay_job_total}"
    say "  manifest     : ${manifest}"
    say "  dag          : ${dag}"
    {
      echo "capture_worker_jobs=0"
      echo "replay_worker_jobs=${replay_job_total}"
      echo "total_worker_jobs=${replay_job_total}"
      echo "dag_max_worker_jobs=${RJ_DAG_MAX_WORKER_JOBS:-50000}"
      echo "dag_warn_worker_jobs=${RJ_DAG_WARN_WORKER_JOBS:-40000}"
    } >> "$manifest"
    sim_min_worker_jobs="${RJ_DAG_MIN_WORKER_JOBS:-${RJ_DAG_MIN_WORKER_JOBS_SIM:-15000}}"
    echo "dag_min_worker_jobs=${sim_min_worker_jobs}" >> "$manifest"
    check_dag_worker_job_budget "SIM pool replay ${TAG}" 0 "$replay_job_total" 0 "$sim_min_worker_jobs"
    publish_production_manifest "$manifest" "$DEST_BASE"
    submit_dag_with_notify "$dag"
    ;;

  condorDoAllDirect)
    [[ "$IS_SIM" -eq 1 ]] || { err "condorDoAllDirect is only valid for isSim variants"; exit 2; }
    warn "condorDoAllDirect reopens DSTs. Prefer condorDoAllFromScratch once to build pools, then condorDoAll for replay."
    _direct_args=( "$DATASET" condorDoAll groupSize "$GROUP_SIZE" )
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]]; then
      _direct_args+=( "SAMPLE=${SIM_SAMPLE}" )
    fi
    if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 ]]; then
      _direct_args+=( maxJobs "$MAX_JOBS" )
    fi
    RJ_DIRECT_DST_DOALL=1 "$0" "${_direct_args[@]}"
    ;;

  condorDoAll)
    # Submit all sim files in grouped chunks. MUST wipe outputs/logs first.
    [[ "$IS_SIM" -eq 1 ]] || { err "condorDoAll is only valid for isSim variants"; exit 2; }
    if [[ "${RJ_DIRECT_DST_DOALL:-0}" != "1" && "${RJ_POOL_MODE:-}" != "capture" && "${RJ_POOL_MODE:-}" != "captureOnly" && "${RJ_POOL_MODE:-}" != "capture_only" ]]; then
      say "${BOLD}condorDoAll now replays existing AnalysisPool ROOT files.${RST}"
      say "Use ${BOLD}${0} ${DATASET} condorDoAllFromScratch${RST} to build pools from DST first."
      _pool_replay_args=( "$DATASET" condorHistFromPool groupSize "$GROUP_SIZE" )
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]]; then
        _pool_replay_args+=( "SAMPLE=${SIM_SAMPLE}" )
      fi
      if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 ]]; then
        _pool_replay_args+=( maxJobs "$MAX_JOBS" )
      fi
      "$0" "${_pool_replay_args[@]}"
      exit $?
    fi
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
        isSimEmbedded)          samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
        isSimEmbeddedInclusive) samples=( "run28_embeddedJet12" "run28_embeddedJet20" ) ;;
        isSimJet5)              samples=( "run28_jet5" ) ;;
        isSimMB)                samples=( "run28_detroit" ) ;;
        *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
      esac
    else
      samples=( "${SIM_SAMPLE}" )
    fi

    # If CHECKJOBS was also provided, do a dry-run count for the FULL condorDoAll matrix and exit.
    if [[ "${DRYRUN:-0}" -eq 1 ]]; then
      n_cfg=$(( ${#sim_pts[@]} * ${#sim_fracs[@]} * ${#sim_vzs[@]} * ${#sim_cones[@]} * $(iso_group_count) * ${#uepipe_modes[@]} ))
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
      say "  iso groups submitted              : $(iso_group_count) upstream mode(s)"
      say "  photon-ID fanout outputs          : ${#iso_tags[@]} cfg tag(s) across those iso groups"
      if (( ${uepipe_in_tag:-0} )); then
        say "  clusterUEpipeline                 : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} values; tagged)"
      else
        say "  clusterUEpipeline                 : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} forced value; not tagged for ${TAG})"
      fi
      say "  upstream DST job configs          : ${BOLD}${n_cfg}${RST}"
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

      say "${BOLD}Full upstream ID-fanout cell list (${n_cfg} entries × ${#samples[@]} samples = ${BOLD}$((n_cfg * ${#samples[@]}))${RST} submit blocks):${RST}"
      cfg_num=0
      for pt in "${sim_pts[@]}"; do
        for frac in "${sim_fracs[@]}"; do
          for vz in "${sim_vzs[@]}"; do
          for cone in "${sim_cones[@]}"; do
          for (( _ci=0; _ci<${#iso_tags[@]}; _ci++ )); do
          iso_idx_is_group_leader "$_ci" || continue
          for uepipe in "${uepipe_modes[@]}"; do
            (( cfg_num+=1 ))
            _tag="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$_ci" "$uepipe")"
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
      say "  upstream configs       : ${n_cfg}"
      say "  jobs per combo (Σsamp) : ${per_cfg_jobs}"
      say "  ─────────────────────────────────"
      say "  ${BOLD}TOTAL CONDOR JOBS        : ${BOLD}${total_jobs}${RST}"
      echo
      _ex_tag="$(matrix_cfg_tag "${sim_pts[0]}" "${sim_fracs[0]}" "${sim_vzs[0]}" "${sim_cones[0]}" 0 "${uepipe_modes[0]}")"
      say "Output tree: each tag becomes a subdirectory under ${DEST_BASE}/"
      say "  e.g. ${DIM}${DEST_BASE}/${_ex_tag}/<sample>/*.root${RST}"
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
    # Clean stale .sub files only. Keep YAML overrides/snapshots because live Condor jobs may still reference them.
    rm -f "${SUB_DIR}/RecoilJets_sim_"*.sub "${SUB_DIR}/RecoilJets_${TAG}_"*.sub 2>/dev/null || true
    SIM_STAGE_NAMESPACE="${TAG}_${doall_stamp}"
    export SIM_STAGE_NAMESPACE
    pool_capture_submit=0
    case "${RJ_POOL_MODE:-}" in
      capture|captureOnly|capture_only) pool_capture_submit=1 ;;
    esac
    say "SIM chunk-list stage namespace: ${SIM_STAGE_NAMESPACE}"
    for pt in "${sim_pts[@]}"; do
      for frac in "${sim_fracs[@]}"; do
        for vz in "${sim_vzs[@]}"; do
        for cone in "${sim_cones[@]}"; do
        for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
        iso_idx_is_group_leader "$iso_idx" || continue
        for uepipe in "${uepipe_modes[@]}"; do
        SIM_CFG_TAG="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$iso_idx" "$uepipe")"
        DEST_BASE="${SIM_DEST_BASE_RESOLVED}/${SIM_CFG_TAG}"
        fanout_dirs="${SIM_YAML_OVERRIDE_DIR}/id_fanout_${SIM_CFG_TAG}_${doall_stamp}.txt"
        if [[ "$pool_capture_submit" -eq 0 ]]; then
          emit_id_fanout_dirs_file "$fanout_dirs" "$SIM_DEST_BASE_RESOLVED" "$pt" "$frac" "$vz" "$cone" "$iso_idx" "$uepipe"
        else
          : > "$fanout_dirs"
        fi
        yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt" "$frac" "$vz" "$cone" "${iso_sliding[$iso_idx]}" "${iso_fixed[$iso_idx]}" "${uepipe}" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}" "$SIM_CFG_TAG" "$doall_stamp")"

        for samp in "${samples[@]}"; do
          SIM_SAMPLE="$samp"
          GROUP_SIZE="$gs_doall"
          sim_init

          if [[ "$pool_capture_submit" -eq 1 ]]; then
            SIM_OUT_DIR="${DEST_BASE}/${SIM_SAMPLE}"
            mkdir -p "$SIM_OUT_DIR"
            rm -f "${SIM_OUT_DIR}/"*.root 2>/dev/null || true
            find "${SIM_OUT_DIR}" -maxdepth 1 -name "*_LOCAL_*.root" -delete 2>/dev/null || true
          else
            while IFS='|' read -r fan_dest _fan_cfg _fan_pre _fan_tight _fan_nonTight; do
              [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
              SIM_OUT_DIR="${fan_dest}/${SIM_SAMPLE}"
              mkdir -p "$SIM_OUT_DIR"
              rm -f "${SIM_OUT_DIR}/"*.root 2>/dev/null || true
              find "${SIM_OUT_DIR}" -maxdepth 1 -name "*_LOCAL_*.root" -delete 2>/dev/null || true
            done < "$fanout_dirs"
          fi
          rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp"*.list 2>/dev/null || true
          rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_LOCAL_"*.list 2>/dev/null || true
          rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_condorTest_"*.list 2>/dev/null || true

          mapfile -t groups < <( make_sim_groups "$GROUP_SIZE" )
          (( ${#groups[@]} )) || { err "No sim groups produced (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"; exit 30; }
          if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 && "${#groups[@]}" -gt "$MAX_JOBS" ]]; then
            say "Capping ${DATASET} group list for sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG}: ${#groups[@]} → ${MAX_JOBS} jobs"
            groups=( "${groups[@]:0:$MAX_JOBS}" )
          fi

          stamp="$(date +%Y%m%d_%H%M%S)"
          sub="${SUB_DIR}/RecoilJets_sim_${SIM_CFG_TAG}_${SIM_SAMPLE}_${stamp}.sub"
          args_file="${SUB_DIR}/RecoilJets_sim_${SIM_CFG_TAG}_${SIM_SAMPLE}_${stamp}.args"
          : > "$args_file"
          exe_for_sub="${BULK_FROZEN_EXE:-${EXE}}"
          macro_env_for_sub=""
          [[ -n "${BULK_FROZEN_MACRO:-}" ]] && macro_env_for_sub=";RJ_MACRO_PATH=${BULK_FROZEN_MACRO}"
          fanout_env_for_sub=""
          [[ "$pool_capture_submit" -eq 0 ]] && fanout_env_for_sub=";RJ_ID_FANOUT_DIRS_FILE=${fanout_dirs}"
          pool_env_for_sub=""
          [[ "$pool_capture_submit" -eq 1 ]] && pool_env_for_sub=";RJ_POOL_MODE=capture"

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
notification  = Never
environment   = RJ_VERBOSITY=0;RJ_CONFIG_YAML=${yaml_override}${macro_env_for_sub}${fanout_env_for_sub}${pool_env_for_sub}
queue arguments from ${args_file}
SUB

          gidx=0
          for glist in "${groups[@]}"; do
            (( gidx+=1 ))
            printf '%s %s %s $(Cluster) 0 %d NONE %s\n' \
                   "$SIM_SAMPLE" "$glist" "$DATASET" "$gidx" "$DEST_BASE" >> "$args_file"
          done

          if [[ "$pool_capture_submit" -eq 1 ]]; then
            say "Submitting ${DATASET} pool capture (tag=${SIM_CFG_TAG}, sample=${SIM_SAMPLE}, groupSize=${GROUP_SIZE}, uepipe=${uepipe}) → jobs=${BOLD}${#groups[@]}${RST}"
            say "Pool output dir: ${DEST_BASE}/${SIM_SAMPLE}"
          else
            say "Submitting ${DATASET} condorDoAll ID-fanout (base tag=${SIM_CFG_TAG}, sample=${SIM_SAMPLE}, groupSize=${GROUP_SIZE}, uepipe=${uepipe}) → jobs=${BOLD}${#groups[@]}${RST}"
            say "Fanout outputs: $(wc -l < "$fanout_dirs" | awk '{print $1}') cfg tags from one DST pass"
            say "Fanout dirs file: ${fanout_dirs}"
          fi
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
        data_cfg_tag="${dpt_tag}_${dfrac_tag}_${dvz_tag}_${dcone_tag}_${iso_base_tags[0]}"
        (( uepipe_in_tag )) && data_cfg_tag="${data_cfg_tag}_${uepipe_modes[0]}"
        data_cfg_tag="${data_cfg_tag}_${iso_selection_tags[0]}"

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
        pin_photon_id_scalars_in_yaml "$yaml_snap" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}"
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
notification  = Never
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=10;RJ_CONFIG_YAML=${yaml_snap}
arguments     = ${r8} ${glist} ${DATASET} \$(Cluster) 0 1 NONE ${DEST_BASE}
queue
SUB
        say "Submitting 1 test job on run ${BOLD}${r8}${RST} (first chunk, groupSize=1, vz=${vz0}) → $(basename "$sub")"
        condor_submit "$sub"
        ;;
      poolSmoke)
        smoke_jobs="${RJ_POOL_SMOKE_JOBS:-12}"
        if [[ "${MAX_JOBS_EXPLICIT:-0}" -eq 1 ]]; then
          smoke_jobs="$MAX_JOBS"
        fi
        [[ "$smoke_jobs" =~ ^[0-9]+$ && "$smoke_jobs" -gt 0 ]] || { err "poolSmoke job cap must be a positive integer, got '${smoke_jobs}'"; exit 2; }
        say "${BOLD}DATA poolSmoke requested${RST}"
        say "  dataset        : ${DATASET}"
        say "  groupSize      : ${GROUP_SIZE}"
        say "  capture jobs   : ${smoke_jobs} total cap"
        say "  replay grouping: ${RJ_POOL_REPLAY_GROUP_SIZE:-${smoke_jobs}} pool files/job unless overridden"
        say "  profiling      : enabled via RJ_PROFILE_JOB=1 in capture/replay workers"
        say "  outputs        : isolated _poolSmoke_<timestamp> directories, not production cfg dirs"
        submit_data_pool_workflow 1 "$smoke_jobs" "poolSmoke"
        ;;
      smokeTest)
        if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
          case "$DATASET" in
            isAuAu|isOO) GROUP_SIZE="${RJ_SMOKE_GROUPSIZE_AUAU:-1}" ;;
            *)           GROUP_SIZE="${RJ_SMOKE_GROUPSIZE_PP:-2}" ;;
          esac
        fi
        if [[ -z "${RJ_REQUEST_MEMORY:-}" ]]; then
          case "$DATASET" in
            isAuAu|isOO) RJ_REQUEST_MEMORY="${RJ_SMOKE_REQUEST_MEMORY_AUAU:-2500MB}" ;;
            *)           RJ_REQUEST_MEMORY="${RJ_SMOKE_REQUEST_MEMORY_PP:-1800MB}" ;;
          esac
          export RJ_REQUEST_MEMORY
        fi
        smoke_jobs="${RJ_SMOKE_CAPTURE_JOBS:-${RJ_SMOKE_DATA_RUNS:-10}}"
        if [[ "${MAX_JOBS_EXPLICIT:-0}" -eq 1 ]]; then
          smoke_jobs="$MAX_JOBS"
        fi
        [[ "$smoke_jobs" =~ ^[0-9]+$ && "$smoke_jobs" -gt 0 ]] || { err "DATA smokeTest capture cap must be a positive integer, got '${smoke_jobs}'"; exit 2; }
        say "${BOLD}DATA smokeTest requested${RST}"
        say "  dataset       : ${DATASET}"
        say "  selected runs : ${RJ_SMOKE_DATA_RUNS:-10} largest-statistics golden runs"
        say "  groupSize     : ${GROUP_SIZE}"
        say "  capture jobs  : ${smoke_jobs} total cap"
        say "  jobs/run cap  : ${RJ_SMOKE_CAPTURE_JOBS_PER_RUN:-1}"
        say "  nEvents       : ${RJ_SMOKE_DATA_NEVENTS:-3000} per capture worker (0 means full worker input)"
        say "  request mem   : ${RJ_REQUEST_MEMORY}"
        say "  heartbeats    : ${RJ_SMOKE_JOB_HEARTBEAT_SECONDS:-120}s while workers run"
        say "  outputs       : isolated thesisAnaSmoke/thesisAnaPoolsSmoke directories"
        submit_data_pool_workflow 1 "$smoke_jobs" "smokeTest"
        ;;
      round)
        seg="${4:?round number required}"
        firstChunk="${5:-}"
        round_file="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
        [[ -s "$round_file" ]] || { err "Round file not found: $round_file. Run 'splitGoldenRunList' first."; exit 8; }

        # Keep existing YAML overrides/snapshots; live Condor jobs may still reference them.
        mkdir -p "$SIM_YAML_OVERRIDE_DIR"

        for data_pt in "${data_pts[@]}"; do
        for data_frac in "${data_fracs[@]}"; do
        for data_vz in "${data_vzs[@]}"; do
          for data_cone in "${data_cones[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
          iso_idx_is_group_leader "$iso_idx" || continue
          for uepipe in "${uepipe_modes[@]}"; do
          data_cfg_tag="$(matrix_cfg_tag "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe")"
          yaml_override="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${data_cfg_tag}.yaml"
          fanout_dirs="${SIM_YAML_OVERRIDE_DIR}/id_fanout_${TAG}_${data_cfg_tag}.txt"
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
          pin_photon_id_scalars_in_yaml "$yaml_override" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}"
          emit_id_fanout_dirs_file "$fanout_dirs" "$DATA_DEST_BASE_SAVED" "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe"
          export RJ_CONFIG_YAML="$yaml_override"
          DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"

          say "DATA condor round (pt=${data_pt}, frac=${data_frac}, vz=${data_vz}, coneR=${data_cone}, iso=${iso_tags[$iso_idx]}, uepipe=${uepipe}, tag=${data_cfg_tag})"
          say "  YAML override: ${yaml_override}"
          say "  Fanout cfgs  : $(wc -l < "$fanout_dirs" | awk '{print $1}') output cfg tags from one DST pass"
          say "  Fanout file  : ${fanout_dirs}"
          while IFS='|' read -r fan_dest _fan_cfg _fan_pre _fan_tight _fan_nonTight; do
            [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
            case "$fan_dest" in
              */thesisAna/pp/*|*/thesisAna/pp25/*|*/thesisAna/auau/*|*/thesisAna/oo/*) ;;
              *) err "Refusing to wipe fanout DEST_BASE='$fan_dest'"; exit 62 ;;
            esac
            mkdir -p "$fan_dest"
            find "$fan_dest" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true
          done < "$fanout_dirs"
          _old_submit_extra_env="${RJ_SUBMIT_EXTRA_ENV:-}"
          export RJ_SUBMIT_EXTRA_ENV="${_old_submit_extra_env:+${_old_submit_extra_env};}RJ_ID_FANOUT_DIRS_FILE=${fanout_dirs}"
          submit_condor "$round_file" "$firstChunk"
          export RJ_SUBMIT_EXTRA_ENV="$_old_submit_extra_env"
          done
          done
          done
        done
        done
        done
        ;;
      allFromScratch)
        submit_data_pool_workflow 1
        ;;
      allDirect)
        warn "condor allDirect reopens DSTs. Prefer 'condor allFromScratch' once, then 'condor all' for pool replay."
        _direct_args=( "$DATASET" condor all groupSize "$GROUP_SIZE" )
        if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 ]]; then
          _direct_args+=( maxJobs "$MAX_JOBS" )
        fi
        RJ_DIRECT_DST_DOALL=1 "$0" "${_direct_args[@]}"
        ;;
      all|"")
        if [[ "${RJ_DIRECT_DST_DOALL:-0}" != "1" && "${RJ_POOL_MODE:-}" != "capture" && "${RJ_POOL_MODE:-}" != "captureOnly" && "${RJ_POOL_MODE:-}" != "capture_only" ]]; then
          say "${BOLD}condor all now replays existing AnalysisPool ROOT files for DATA.${RST}"
          say "Use ${BOLD}${0} ${DATASET} condor allFromScratch${RST} to build pools from DST first."
          submit_data_pool_workflow 0
        else
        # Keep existing YAML overrides/snapshots; live Condor jobs may still reference them.
        mkdir -p "$SIM_YAML_OVERRIDE_DIR"

        # Freeze pipeline for this bulk submission
        cleanup_bulk_snapshots_for_tag
        case "$DATASET" in
          isAuAu|isOO) create_pipeline_snapshot "auau" "$(date +%Y%m%d_%H%M%S)" ;;
          *)           create_pipeline_snapshot "pp" "$(date +%Y%m%d_%H%M%S)" ;;
        esac

        n_matrix=$(( ${#data_pts[@]} * ${#data_fracs[@]} * ${#data_vzs[@]} * ${#data_cones[@]} * $(iso_group_count) * ${#uepipe_modes[@]} ))
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
          iso_idx_is_group_leader "$iso_idx" || continue
          for uepipe in "${uepipe_modes[@]}"; do
          data_cfg_tag="$(matrix_cfg_tag "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe")"
          yaml_override="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${data_cfg_tag}.yaml"
          fanout_dirs="${SIM_YAML_OVERRIDE_DIR}/id_fanout_${TAG}_${data_cfg_tag}.txt"
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
          pin_photon_id_scalars_in_yaml "$yaml_override" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}"
          emit_id_fanout_dirs_file "$fanout_dirs" "$DATA_DEST_BASE_SAVED" "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe"
          export RJ_CONFIG_YAML="$yaml_override"
          DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"

          (( ++cell_num ))
          say "───────────────────────────────────────────────────────────────"
          say "${BOLD}[${cell_num}/${n_matrix}] ${data_cfg_tag}${RST}"
          say "  YAML override : ${yaml_override}"
          say "  Fanout cfgs   : $(wc -l < "$fanout_dirs" | awk '{print $1}') output cfg tags from one DST pass"
          say "  Fanout file   : ${fanout_dirs}"
          tmp_src="${ROUND_DIR}/ALL_${TAG}_${data_cfg_tag}_$(date +%s).txt"
          grep -E '^[0-9]+' "$GOLDEN" > "$tmp_src"

          while IFS='|' read -r fan_dest _fan_cfg _fan_pre _fan_tight _fan_nonTight; do
            [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
            case "$fan_dest" in
              */thesisAna/pp/*|*/thesisAna/pp25/*|*/thesisAna/auau/*|*/thesisAna/oo/*) ;;
              *) err "Refusing to wipe fanout DEST_BASE='$fan_dest'"; exit 62 ;;
            esac
            mkdir -p "$fan_dest"
            find "$fan_dest" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true
          done < "$fanout_dirs"

          _old_submit_extra_env="${RJ_SUBMIT_EXTRA_ENV:-}"
          export RJ_SUBMIT_EXTRA_ENV="${_old_submit_extra_env:+${_old_submit_extra_env};}RJ_ID_FANOUT_DIRS_FILE=${fanout_dirs}"
          submit_condor "$tmp_src" ""
          export RJ_SUBMIT_EXTRA_ENV="$_old_submit_extra_env"
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
        fi
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
