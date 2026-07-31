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
#       isPP, isPPrun25, isAuAu, isOO, isSim, isSimJet5, isSimMB,
#       isSimEmbedded, isSimEmbeddedInclusive
#   • DATA modes use per-run runlists and can do:
#       local, isLocalIsoPing, CHECKJOBS, orchestrationSelfTest, splitGoldenRunList,
#       condor testJob, condor round, condor smokeTest, condor all
#   • SIM modes use staged matched master lists and can do:
#       local, checkModels, CHECKJOBS, orchestrationSelfTest, condorTest,
#       condorDoAllSmoke, condorDoAll, condorDoAllDirect
#   • Jobs never mix files across runs in DATA mode.
#   • Production sweeps use the legacy RecoilJets histogram engine with direct
#     photon-ID fanout inside each Fun4All DST pass. Isolation cone/mode are
#     internal histogram views by default, so current one-UE productions submit
#     15 final cfg ROOT files while each file contains its configured iso/cone
#     views. Au+Au defaults to sliding R=0.4 plus sliding R=0.3 only.
#     jet_pt_min and back_to_back_dphi_min_pi_fraction are also internal recoil
#     histogram scans by default, so they do not force repeated DST passes.
#     vz_cut_cm is collapsed by dataset default unless RJ_VZ_SCAN_ALL=1
#     or RJ_VZ_CUT_OVERRIDE=<cm> is set.
#     Each fanout output is still an ordinary cfg-tagged ROOT file filled by
#     RecoilJets/RecoilJets_AuAu, so AnalyzeRecoilJets* sees the legacy output
#     contract unchanged.
#   • Pool/replay histogram production has been removed from this submitter.
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
#     isSimJet5     | isSimInclusive | simjet5 | siminclusive | SIMJET5
#     isSimMB       | simmb       | SIMMB
#     isSimEmbedded          | simembedded          | SIMEMBEDDED
#     isSimEmbeddedInclusive | simembeddedinclusive | SIMEMBEDDEDINCLUSIVE
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
#   │ isSimEmbeddedInclusive│ RecoilJets_Condor_AuAu.sh│ Fun4All_recoilJets_AuAu.C│ condorDoAll                 │
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
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive local
#     ./RecoilJets_Condor_submit.sh isSimEmbeddedInclusive condorDoAll
#
# INPUT CONTRACTS
#   DATA:
#     • Golden run list + per-run list files.
#     • PP Run-24 PPG12 parity uses paired per-run lists:
#       col1 = DST_Jet, col2 = DST_JETCALO, both from ana521_2025p007_v001.
#       This mirrors ppg12codeGit/anatreemaker/macro_maketree/data/ana521/run.sh.
#     • Other DATA per-run list files contain one ROOT file path per line,
#       except Au+Au lists that opt into paired DST_ZDC_RAW for MB classification.
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
#     • isSimEmbedded and isSimEmbeddedInclusive use the same 5-column staged contract.
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
#     • PP      : ${BASE}/GRLs_tanner/run2pp_ana521_2025p007_v001_ppg12_runlist.list
#                 or, if absent, Shuhang's PPG12 ana521 runList.txt
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
#     • isSimInclusive → /sphenix/tg/tg01/bulk/jbennett/thesisAna/siminclusive
#     • isSimJet5      → legacy alias for isSimInclusive
#     • isSimMB       → /sphenix/tg/tg01/bulk/jbennett/thesisAna/simmb
#     • isSimEmbedded → /sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded
#     • isSimEmbeddedInclusive → /sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive
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
#       production output cleanup is explicit via RJ_CLEAN_OUTPUT_BASE=1.
#   • condorDoAll() in SIM mode:
#       production output cleanup is explicit via RJ_CLEAN_OUTPUT_BASE=1.
#   • local mode:
#       removes disposable *_LOCAL* YAML/list/root artifacts before running.
#   • isLocalIsoPing:
#       creates a combined local chunk list but does not alter condor segmentation.
#
# REQUIREMENTS
#   • condor_submit must be on PATH for Condor submissions.
#   • psql must be on PATH if using TRIGGER=<bit> or isLocalIsoPing.
#   • Required golden lists, per-run lists, and staged sim lists must exist.
###############################################################################
set -euo pipefail

# Shared sPHENIX runtime artifacts must remain inspectable on Lustre/GPFS.
# Credentials and private agent state do not belong under this tree.
umask 0022

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
# Allow independent submit-file namespaces when one campaign is intentionally
# distributed across multiple schedds.  The default remains unchanged for all
# existing callers; a campaign may set RJ_CONDOR_SUB_DIR per submit host to
# prevent concurrent submitters from removing or overwriting one another's
# transient .sub/.args files.
SUB_DIR="${RJ_CONDOR_SUB_DIR:-${BASE}/condor_sub}"
CLEANUP_HELPER="${BASE}/scripts/recoiljets_cleanup.sh"

# Golden lists provided by you
PP_GOLDEN_DEFAULT="${BASE}/GRLs_tanner/run2pp_ana521_2025p007_v001_ppg12_runlist.list"
PP_GOLDEN_PPG12_SOURCE="/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/runList.txt"
PP_GOLDEN="${RJ_PP_GOLDEN:-$PP_GOLDEN_DEFAULT}"
if [[ -z "${RJ_PP_GOLDEN:-}" && ! -f "$PP_GOLDEN" && -f "$PP_GOLDEN_PPG12_SOURCE" ]]; then
  PP_GOLDEN="$PP_GOLDEN_PPG12_SOURCE"
fi
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
SCALED_TRIGGER_JET_PT_MIN="${RJ_SCALED_TRIGGER_JET_PT_MIN:-5.0}"
SCALED_TRIGGER_DPHI_FRAC="${RJ_SCALED_TRIGGER_DPHI_FRAC:-0.875}"
SCALED_TRIGGER_VZ_CUT_CM="${RJ_SCALED_TRIGGER_VZ_CUT_CM:-60}"
SCALED_TRIGGER_CONE_R="${RJ_SCALED_TRIGGER_CONE_R:-0.40}"
SCALED_TRIGGER_IS_SLIDING="${RJ_SCALED_TRIGGER_IS_SLIDING:-true}"
SCALED_TRIGGER_FIXED_GEV="${RJ_SCALED_TRIGGER_FIXED_GEV:-4.0}"
SCALED_TRIGGER_UEPIPELINE="${RJ_SCALED_TRIGGER_UEPIPELINE:-baseVariant}"
SCALED_TRIGGER_CENT_EDGES="${RJ_SCALED_TRIGGER_CENT_EDGES:-0,20,50,80}"
SCALED_TRIGGER_CENT_OUTPUT_SUFFIX="${RJ_SCALED_TRIGGER_CENT_OUTPUT_SUFFIX:-_scaledTriggerCentStudy_cent0_20_50_80}"

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
SIM_ROOT="${RJ_SIM_ROOT_OVERRIDE:-${BASE}/simListFiles}"
SIM_SAMPLE_DEFAULT="run28_photonjet10"
SIM_SAMPLE="${SIM_SAMPLE_DEFAULT}"

# Name for the staged 5-column master list (built by sim_init from matched lists)
SIM_MASTER5_NAME="DST_SIM_MASTER_5COL.list"

# Output dir root for sim is: ${SIM_DEST_BASE}/<cfgTag>/<SIM_SAMPLE>
SIM_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/sim"
SIMINCLUSIVE_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/siminclusive"
SIMJET5_DEST_BASE="$SIMINCLUSIVE_DEST_BASE"
SIMMB_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simmb"
SIMEMBED_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded"
SIMEMBEDINCLUSIVE_DEST_BASE="/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembeddedinclusive"
AUAU_BDT_DEST_BASE="${RJ_AUAU_BDT_DEST_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAna/simembedded_bdt_training}"
AUAU_BDT_LOCAL_BASE="${RJ_AUAU_BDT_LOCAL_BASE:-${BASE}/local_bdt_training_outputs}"
AUAU_BDT_MODEL_BASE="${RJ_AUAU_BDT_MODEL_BASE:-${BASE}/bdt_models}"
AUAU_BDT_SIGNAL_SAMPLES_DEFAULT="run28_embeddedPhoton12 run28_embeddedPhoton20"
AUAU_BDT_BACKGROUND_SAMPLES_DEFAULT="run28_embeddedJet12 run28_embeddedJet20"

# Keep the long-standing embedded-inclusive default as Jet12+Jet20. Set
# RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES=1 for the Jet12+Jet20+Jet30 stitching
# study, or RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=1 for Jet12+20+30+40.
simembeddedinclusive_sample_list() {
  if [[ "${RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES:-0}" == "1" || "${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40:-0}" == "1" ]]; then
    printf "%s\n" "run28_embeddedJet12" "run28_embeddedJet20" "run28_embeddedJet30" "run28_embeddedJet40"
  elif [[ "${RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES:-0}" == "1" || "${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30:-0}" == "1" ]]; then
    printf "%s\n" "run28_embeddedJet12" "run28_embeddedJet20" "run28_embeddedJet30"
  else
    printf "%s\n" "run28_embeddedJet12" "run28_embeddedJet20"
  fi
}

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
BULK_FROZEN_SNAPSHOT_LOADER_RECEIPT=""
BULK_FROZEN_SNAPSHOT_LOADER_RECEIPT_SHA256=""
BULK_FROZEN_SNAPSHOT_MANIFEST=""
BULK_FROZEN_SNAPSHOT_MANIFEST_SHA256=""
PPG12_ARCHIVED_OFFLINE_MAIN="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.541"
BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST=""
BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST_SHA256=""

ppg12_archived_di_sample() {
  case "${1:-}" in
    run28_photonjet5_double|run28_photonjet10_double|run28_photonjet20_double|\
    run28_jet8_double|run28_jet12_double|run28_jet20_double|run28_jet30_double|run28_jet40_double) return 0 ;;
  esac
  return 1
}

ppg12_archived_di_lane() {
  case "${1:-${DATASET:-}}" in
    isSim|isSimInclusive) ;;
    *) return 1 ;;
  esac
  ppg12_archived_di_sample "${2:-${SIM_SAMPLE:-}}"
}

ppg12_archived_di_campaign_requested() {
  case "${DATASET:-}" in
    isSim|isSimInclusive) ;;
    *) return 1 ;;
  esac
  [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]] || return 1
  ppg12_archived_di_sample "${SIM_SAMPLE:-}"
}

validate_ppg12_archived_canary_manifest() {
  local manifest="${1:?accepted full-group canary manifest required}"
  local expected_sha="${2:?accepted full-group canary manifest SHA-256 required}"
  [[ "$manifest" == /* && -s "$manifest" ]] || {
    err "Archived PPG12 DI production requires an absolute readable accepted-canary manifest."
    return 96
  }
  [[ "$expected_sha" =~ ^[0-9a-f]{64}$ ]] || {
    err "Archived PPG12 DI accepted-canary manifest SHA-256 is invalid."
    return 96
  }
  [[ "$(sha256sum "$manifest" | awk '{print $1}')" == "$expected_sha" ]] || {
    err "Archived PPG12 DI accepted-canary manifest hash drift."
    return 96
  }
  python3 - "$manifest" "$PPG12_ARCHIVED_OFFLINE_MAIN" <<'PY'
import json
from pathlib import Path
import sys

path = Path(sys.argv[1]).resolve()
expected_offline = sys.argv[2]
data = json.loads(path.read_text(encoding="utf-8"))
software = data.get("software", {})
contract = data.get("contract", {})
comparison = contract.get("comparison", {})
strict = contract.get("strict_post_audit", {})
if data.get("status") != "full_group_root_and_comparator_gates_passed":
    raise SystemExit("accepted canary status is not the full-group production gate")
if software.get("offline_main") != expected_offline:
    raise SystemExit("accepted canary does not bind the fixed ana.541 release")
if contract.get("mode") != "full_group" or contract.get("dataset") != "isSimInclusive":
    raise SystemExit("accepted canary is not an inclusive full-group result")
if contract.get("sample") != "run28_jet20_double" or contract.get("truth_jets_mode") != "DST":
    raise SystemExit("accepted canary sample/truth contract drift")
if contract.get("expected_processed_events") != 5000 or contract.get("full_group_rows") != 5:
    raise SystemExit("accepted canary event/source coverage drift")
if strict.get("status") != "passed":
    raise SystemExit("accepted canary strict post-audit is not passed")
status = str(comparison.get("physics_closure_status", ""))
ac_status = str(comparison.get("A_C", {}).get("physics_status", ""))
if "A_C_closed" not in status and not ac_status.startswith("closed_"):
    raise SystemExit("accepted canary lacks A/C physics closure")
PY
}

normalize_frozen_snapshot_wrapper() {
  local wrapper="${1:?frozen wrapper required}"
  local release_prefix="${2:-}"
  python3 - "$wrapper" "$release_prefix" <<'PY'
from pathlib import Path
import sys

wrapper = Path(sys.argv[1])
release_prefix = sys.argv[2]
text = wrapper.read_text()

marker = "# RJ_SNAPSHOT_LIBRARY_PREPEND_V1"
canonical_export = (
    '  export LD_LIBRARY_PATH="${snapshot_lib_dir}${snapshot_loader_suffix}:'
    '${LD_LIBRARY_PATH:-}"'
)
legacy_export = '  export LD_LIBRARY_PATH="${snapshot_lib_dir}:${LD_LIBRARY_PATH:-}"'
suffix_prefix = 'snapshot_loader_suffix='
dataset_anchor = "# ------------------------ Dataset routing"
loader_anchors = (
    "# Some frozen release lanes intentionally pair",
    "# A bounded diagnostic may replace",
    dataset_anchor,
)
root_include = (
    '  export ROOT_INCLUDE_PATH="$wrapper_dir:$wrapper_dir/include:'
    '${ROOT_INCLUDE_PATH:-}"'
)

if dataset_anchor not in text:
    raise SystemExit(f"snapshot wrapper insertion anchor not found in {wrapper}")
if text.count(marker) > 1:
    raise SystemExit(f"snapshot wrapper has duplicate prepend markers: {wrapper}")
if text.count(canonical_export) > 1 or text.count(legacy_export) > 1:
    raise SystemExit(f"snapshot wrapper has duplicate LD_LIBRARY_PATH prepends: {wrapper}")

suffix_literal = release_prefix.replace("\\", "\\\\").replace('"', '\\"')
suffix_line = f'snapshot_loader_suffix="{suffix_literal}"'

if text.count(marker) == 1 and text.count(canonical_export) == 1:
    suffix_rows = [
        row for row in text.splitlines() if row.startswith(suffix_prefix)
    ]
    if len(suffix_rows) != 1:
        raise SystemExit(
            f"snapshot wrapper marker/export lacks one loader suffix: {wrapper}"
        )
    text = text.replace(suffix_rows[0], suffix_line, 1)
elif text.count(marker) == 0 and text.count(legacy_export) == 1:
    # Normalize the already-correct legacy Au+Au prepend rather than adding a
    # second one. The stable marker lets later preflight distinguish a real
    # prepend from an unrelated snapshot_lib_dir assignment.
    text = text.replace(
        legacy_export,
        f"{suffix_line}\n{marker}\n{canonical_export}",
        1,
    )
elif (
    text.count(marker) == 0
    and text.count(canonical_export) == 0
    and text.count(legacy_export) == 0
):
    # This is the old p+p false-positive shape: snapshot_lib_dir may already
    # exist solely for the optional builder override, but there is no loader
    # prepend. Install the complete contract before any loader/dataset logic.
    insert_at = min(
        (text.index(anchor) for anchor in loader_anchors if anchor in text),
        default=-1,
    )
    if insert_at < 0:
        raise SystemExit(f"snapshot wrapper loader anchor not found in {wrapper}")
    block = f"""# Frozen Condor snapshots carry a sibling lib/ directory with copied local
# analysis libraries. Put it first so DT_NEEDED SONAME lookups and explicit
# ROOT loads resolve to the same immutable snapshot.
wrapper_dir="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd -P)"
snapshot_lib_dir="${{RJ_SNAPSHOT_LIB_DIR:-${{wrapper_dir}}/lib}}"
{suffix_line}
{marker}
if [[ -d "$snapshot_lib_dir" ]]; then
{canonical_export}
  echo "[INFO] Snapshot lib prepended: $snapshot_lib_dir"
fi

"""
    text = text[:insert_at] + block + text[insert_at:]
else:
    raise SystemExit(
        f"snapshot wrapper has an incomplete prepend marker/export contract: {wrapper}"
    )

if text.count(root_include) == 0:
    root_block = f"""if [[ -d "$wrapper_dir" ]]; then
{root_include}
fi

"""
    text = text.replace(dataset_anchor, root_block + dataset_anchor, 1)
elif text.count(root_include) != 1:
    raise SystemExit(f"snapshot wrapper has duplicate ROOT include prepends: {wrapper}")

wrapper.write_text(text)
PY
}

pin_frozen_snapshot_release() {
  local wrapper="${1:?frozen wrapper required}"
  local release_name="${2:?release name required}"
  local offline_main="${3:?offline prefix required}"
  python3 - "$wrapper" "$release_name" "$offline_main" <<'PY'
from pathlib import Path
import re
import sys

wrapper = Path(sys.argv[1])
release_name = sys.argv[2]
offline_main = sys.argv[3]
text = wrapper.read_text()
marker = "# RJ_PINNED_SPHENIX_RELEASE_V1"
setup = "source /opt/sphenix/core/bin/sphenix_setup.sh -n"
setup_re = re.compile(
    r"^(?P<indent>[ \t]*)source /opt/sphenix/core/bin/sphenix_setup\.sh -n[ \t]*$",
    re.MULTILINE,
)
matches = list(setup_re.finditer(text))
if text.count(marker) != 0:
    raise SystemExit(f"frozen wrapper already has a pinned-release marker: {wrapper}")
if len(matches) != 1:
    raise SystemExit(
        f"frozen wrapper must have exactly one unversioned setup call, got "
        f"{len(matches)}: {wrapper}"
    )
indent = matches[0].group("indent")
replacement = f"""{indent}{marker}
{indent}unset LD_PRELOAD
{indent}{setup} {release_name}
{indent}if [[ "${{OFFLINE_MAIN:-}}" != "{offline_main}" ]]; then
{indent}  echo "[FATAL] Frozen snapshot resolved OFFLINE_MAIN='${{OFFLINE_MAIN:-<unset>}}', expected '{offline_main}'."
{indent}  exit 96
{indent}fi"""
text = text[: matches[0].start()] + replacement + text[matches[0].end() :]
wrapper.write_text(text)
PY
}

condor_getenv_directive() {
  case "${RJ_CONDOR_SEALED_ENVIRONMENT:-0}" in
    0|false|FALSE|no|NO|'')
      printf 'True\n'
      ;;
    1|true|TRUE|yes|YES)
      printf 'False\n'
      ;;
    *)
      err "RJ_CONDOR_SEALED_ENVIRONMENT must be 0 or 1"
      return 2
      ;;
  esac
}

snapshot_sha256_file() {
  local path="${1:?snapshot artifact required}"
  python3 - "$path" <<'PY'
from pathlib import Path
import hashlib
import sys

value = hashlib.sha256()
with Path(sys.argv[1]).open("rb") as stream:
    for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
        value.update(block)
print(value.hexdigest())
PY
}

stage_snapshot_soname_aliases() {
  local snap_lib_dir="${1:?snapshot library directory required}"
  local pinned_calo_reco="${2:-0}"
  local expected_calo_reco_soname="${3:-}"
  local snap_so soname alias_path target_name

  if ! command -v readelf >/dev/null 2>&1; then
    if [[ "$pinned_calo_reco" == "1" ]]; then
      err "Pinned CaloReco snapshot requires readelf to prove and stage its SONAME"
      return 2
    fi
    warn "readelf not available; snapshot SONAME links were not generated"
    return 0
  fi

  for snap_so in "${snap_lib_dir}"/lib*.so*; do
    [[ -f "$snap_so" && ! -L "$snap_so" ]] || continue
    soname="$(readelf -d "$snap_so" 2>/dev/null | awk -F'[][]' '/SONAME/ {print $2; exit}' || true)"
    [[ -n "$soname" ]] || continue
    [[ "$soname" == "$(basename "$soname")" && "$soname" == lib*.so* ]] || {
      err "Snapshot ELF declares an unsafe SONAME: ${snap_so} => ${soname}"
      return 2
    }
    target_name="$(basename "$snap_so")"
    [[ "$soname" == "$target_name" ]] && continue
    alias_path="${snap_lib_dir}/${soname}"
    if [[ -e "$alias_path" || -L "$alias_path" ]]; then
      [[ -L "$alias_path" &&
         "$(readlink "$alias_path")" == "$target_name" &&
         "$alias_path" -ef "$snap_so" ]] || {
        err "Snapshot SONAME alias collides with a different provider: ${alias_path}"
        return 2
      }
    else
      ln -s "$target_name" "$alias_path"
    fi
  done

  if [[ "$pinned_calo_reco" == "1" ]]; then
    [[ "$expected_calo_reco_soname" =~ ^libcalo_reco\.so(\.[0-9]+)+$ ]] || {
      err "Pinned CaloReco requires an explicit versioned SONAME contract"
      return 2
    }
    snap_so="${snap_lib_dir}/libcalo_reco.so"
    soname="$(readelf -d "$snap_so" 2>/dev/null | awk -F'[][]' '/SONAME/ {print $2; exit}' || true)"
    [[ "$soname" == "$expected_calo_reco_soname" ]] || {
      err "Pinned CaloReco SONAME differs: observed=${soname:-<unset>} expected=${expected_calo_reco_soname}"
      return 2
    }
    alias_path="${snap_lib_dir}/${expected_calo_reco_soname}"
    [[ -L "$alias_path" &&
       "$(readlink "$alias_path")" == "libcalo_reco.so" &&
       "$alias_path" -ef "$snap_so" ]] || {
      err "Pinned CaloReco SONAME alias does not resolve to the one snapshot provider"
      return 2
    }
  fi
}

validate_frozen_snapshot_wrapper_contract() {
  local wrapper="${1:?frozen wrapper required}"
  local marker="# RJ_SNAPSHOT_LIBRARY_PREPEND_V1"
  local export_line='  export LD_LIBRARY_PATH="${snapshot_lib_dir}${snapshot_loader_suffix}:${LD_LIBRARY_PATH:-}"'
  local marker_count export_count suffix_count suffix_line marker_line export_line_number loader_line dataset_line boundary_line

  marker_count="$(grep -Fxc "$marker" "$wrapper" || true)"
  export_count="$(grep -Fxc "$export_line" "$wrapper" || true)"
  suffix_count="$(grep -Ec '^snapshot_loader_suffix="[^"]*"$' "$wrapper" || true)"
  if [[ "$marker_count" -ne 1 || "$export_count" -ne 1 || "$suffix_count" -ne 1 ]]; then
    err "Frozen wrapper must contain exactly one marked snapshot-library prepend: ${wrapper}"
    return 2
  fi

  suffix_line="$(grep -nE '^snapshot_loader_suffix="[^"]*"$' "$wrapper" | cut -d: -f1)"
  marker_line="$(grep -Fn "$marker" "$wrapper" | cut -d: -f1)"
  export_line_number="$(grep -Fn "$export_line" "$wrapper" | cut -d: -f1)"
  dataset_line="$(grep -Fn '# ------------------------ Dataset routing' "$wrapper" | cut -d: -f1)"
  if [[ -z "$dataset_line" ]]; then
    err "Frozen wrapper lacks the dataset-routing boundary: ${wrapper}"
    return 2
  fi
  loader_line="$(
    grep -nE \
      '^# (Some frozen release lanes intentionally pair|A bounded diagnostic may replace)' \
      "$wrapper" | head -n 1 | cut -d: -f1 || true
  )"
  boundary_line="${loader_line:-$dataset_line}"
  if ! [[ "$suffix_line" -lt "$marker_line" &&
          "$marker_line" -lt "$export_line_number" &&
          "$export_line_number" -lt "$boundary_line" &&
          "$export_line_number" -lt "$dataset_line" ]]; then
    err "Frozen wrapper loader order must be suffix < marker < export < ROOT loader < dataset routing: ${wrapper}"
    return 2
  fi

  if ! grep -Fq 'echo "[INFO] Snapshot lib prepended: ' "$wrapper"; then
    err "Frozen wrapper lacks the stable snapshot-library runtime witness: ${wrapper}"
    return 2
  fi
  return 0
}

write_and_seal_snapshot_manifest() {
  local snap_dir="${1:?snapshot directory required}"
  local manifest="${snap_dir}/snapshot_manifest.json"
  python3 - "$snap_dir" "$manifest" <<'PY'
from pathlib import Path
import hashlib
import json
import os
import sys

root = Path(sys.argv[1]).resolve()
manifest = Path(sys.argv[2])
if not root.is_dir() or manifest.exists():
    raise SystemExit(f"snapshot manifest precondition failed: {root}")

entries = []
for path in sorted(root.rglob("*"), key=lambda item: item.relative_to(root).as_posix()):
    relative = path.relative_to(root).as_posix()
    stat_result = path.lstat()
    if path.is_symlink():
        target = os.readlink(path)
        if os.path.isabs(target):
            raise SystemExit(f"snapshot symlink must be relative: {relative} -> {target}")
        resolved = path.resolve(strict=True)
        if root not in resolved.parents:
            raise SystemExit(f"snapshot symlink escapes its root: {relative} -> {target}")
        digest = hashlib.sha256(resolved.read_bytes()).hexdigest()
        entries.append(
            {
                "path": relative,
                "sha256": digest,
                "symlink_target": target,
                "type": "symlink",
            }
        )
    elif path.is_file():
        entries.append(
            {
                "path": relative,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "size": stat_result.st_size,
                "type": "file",
            }
        )
    elif path.is_dir():
        entries.append({"path": relative, "type": "directory"})
    else:
        raise SystemExit(f"unsupported snapshot artifact type: {relative}")

payload = {
    "entries": entries,
    "root": str(root),
    "schema": "RJ_FROZEN_SNAPSHOT_MANIFEST_V1",
    "status": "PASS",
}
manifest.write_text(
    json.dumps(payload, indent=2, sort_keys=True, separators=(",", ": ")) + "\n",
    encoding="utf-8",
)
PY
  chmod -R a-w "$snap_dir"
  [[ -s "$manifest" ]] || {
    err "Frozen snapshot manifest was not created: ${manifest}"
    return 2
  }
  BULK_FROZEN_SNAPSHOT_MANIFEST="$manifest"
  BULK_FROZEN_SNAPSHOT_MANIFEST_SHA256="$(snapshot_sha256_file "$manifest")"
}

validate_snapshot_loader_closure() {
  local snap_lib_dir="${1:?snapshot lib directory required}"
  local mode="${2:?snapshot mode required}"
  local use_release_core="${3:?release-core flag required}"
  local release_core_lib64="${4:?release lib64 directory required}"
  local release_core_lib="${5:?release lib directory required}"
  local pinned_calo_reco="${6:-0}"
  local receipt="${7:-}"
  local release_prefix=""
  local report target report_is_temporary=1
  local expected_calo_io="${RJ_PINNED_RELEASE_CALO_IO_PATH:-}"
  local expected_calo_io_sha="${RJ_PINNED_RELEASE_CALO_IO_SHA256:-}"
  local expected_clusteriso="${RJ_PINNED_RELEASE_CLUSTERISO_PATH:-}"
  local expected_clusteriso_sha="${RJ_PINNED_RELEASE_CLUSTERISO_SHA256:-}"
  local expected_jetbase="${RJ_PINNED_RELEASE_JETBASE_PATH:-}"
  local expected_jetbase_sha="${RJ_PINNED_RELEASE_JETBASE_SHA256:-}"
  local -a targets=()

  command -v ldd >/dev/null 2>&1 || {
    err "ldd is required for frozen snapshot loader preflight"
    return 2
  }
  [[ -d "$snap_lib_dir" ]] || {
    err "Frozen snapshot library directory is missing: ${snap_lib_dir}"
    return 2
  }

  if [[ "$mode" == "auau" ]]; then
    targets+=("${snap_lib_dir}/libRecoilJetsAuAu.so")
  else
    targets+=("${snap_lib_dir}/libRecoilJets.so")
  fi
  for target in \
    "${snap_lib_dir}/libcalo_reco.so" \
    "${snap_lib_dir}/libcalo_io.so" \
    "${snap_lib_dir}/libclusteriso.so" \
    "${snap_lib_dir}/libjetbase.so"; do
    [[ -f "$target" ]] && targets+=("$target")
  done
  for target in "${targets[@]}"; do
    [[ -f "$target" ]] || {
      err "Frozen snapshot loader target is missing: ${target}"
      return 2
    }
  done

  if [[ "$use_release_core" == "1" ]]; then
    release_prefix=":${release_core_lib64}:${release_core_lib}"
  fi
  if [[ "$pinned_calo_reco" == "1" ]]; then
    local expected_path expected_sha
    for expected_path in \
      "$expected_calo_io" "$expected_clusteriso" "$expected_jetbase"; do
      [[ "$expected_path" == /* && -f "$expected_path" && -s "$expected_path" ]] || {
        err "Pinned release companion is missing or not absolute: ${expected_path:-<unset>}"
        return 2
      }
    done
    for expected_sha in \
      "$expected_calo_io_sha" "$expected_clusteriso_sha" "$expected_jetbase_sha"; do
      [[ "$expected_sha" =~ ^[0-9a-f]{64}$ ]] || {
        err "Pinned release companion has an invalid SHA-256: ${expected_sha:-<unset>}"
        return 2
      }
    done
    # These companions are loaded explicitly by the frozen Fun4All macro and
    # therefore need not appear in any staged library's DT_NEEDED closure.
    # Inspect the exact declared providers directly so loader health and
    # provider identity are still proven without manufacturing a dependency.
    targets+=(
      "$expected_calo_io"
      "$expected_clusteriso"
      "$expected_jetbase"
    )
  fi
  if [[ -n "$receipt" ]]; then
    [[ "$receipt" == /* && ! -e "$receipt" ]] || {
      err "Snapshot loader receipt must be a fresh absolute path: ${receipt}"
      return 2
    }
    report="${receipt%.json}.ldd.txt"
    [[ ! -e "$report" ]] || {
      err "Snapshot loader report already exists: ${report}"
      return 2
    }
    report_is_temporary=0
  else
    report="$(mktemp "${TMPDIR:-/tmp}/rj_snapshot_ldd.XXXXXX")"
  fi
  for target in "${targets[@]}"; do
    printf '@@TARGET %s\n' "$target" >> "$report"
    if ! LD_LIBRARY_PATH="${snap_lib_dir}${release_prefix}:${LD_LIBRARY_PATH:-}" \
      ldd "$target" >> "$report" 2>&1; then
      err "ldd failed for frozen snapshot target: ${target}"
      (( report_is_temporary )) && rm -f "$report"
      return 2
    fi
  done

  if ! python3 - "$report" "$mode" "$snap_lib_dir" "$use_release_core" \
    "$release_core_lib64" "$release_core_lib" "$pinned_calo_reco" "$receipt" \
    "$expected_calo_io" "$expected_calo_io_sha" \
    "$expected_clusteriso" "$expected_clusteriso_sha" \
    "$expected_jetbase" "$expected_jetbase_sha" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

report = Path(sys.argv[1])
mode = sys.argv[2]
snapshot = Path(sys.argv[3]).resolve()
use_release = sys.argv[4] == "1"
release_roots = (Path(sys.argv[5]).resolve(), Path(sys.argv[6]).resolve())
pinned_calo_reco = sys.argv[7] == "1"
receipt = Path(sys.argv[8]) if sys.argv[8] else None
expected_args = sys.argv[9:]
expected_companions = {
    "libcalo_io.so": (Path(expected_args[0]), expected_args[1]),
    "libclusteriso.so": (Path(expected_args[2]), expected_args[3]),
    "libjetbase.so": (Path(expected_args[4]), expected_args[5]),
}
staged_names = {entry.name for entry in snapshot.iterdir()}
staged_families = {
    name.split(".so", 1)[0] + ".so" for name in staged_names if ".so" in name
}
local_core = ("libcalo_reco.so", "libcalo_io.so", "libclusteriso.so", "libjetbase.so")
targets = 0
failures = []
observed: dict[str, set[str]] = {family: set() for family in local_core}

def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

def beneath(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False

expected_real: dict[str, Path] = {
    "libcalo_reco.so": (snapshot / "libcalo_reco.so").resolve(),
}
if pinned_calo_reco:
    for family, (declared, expected_sha) in expected_companions.items():
        try:
            real = declared.resolve(strict=True)
        except (FileNotFoundError, RuntimeError):
            failures.append(f"missing declared release companion: {family} => {declared}")
            continue
        expected_real[family] = real
        if not any(beneath(real, root) for root in release_roots):
            failures.append(f"declared release companion escaped pinned release: {family} => {real}")
        elif digest(real) != expected_sha:
            failures.append(f"declared release companion hash drift: {family} => {real}")

for raw in report.read_text(errors="replace").splitlines():
    line = raw.strip()
    if line.startswith("@@TARGET "):
        targets += 1
        target_text = line.removeprefix("@@TARGET ")
        target_path = Path(target_text)
        try:
            target_real = target_path.resolve(strict=True)
        except (FileNotFoundError, RuntimeError):
            failures.append(f"missing inspected loader target: {target_text}")
            continue
        target_family = next(
            (
                candidate
                for candidate in local_core
                if target_path.name == candidate
                or target_path.name.startswith(candidate + ".")
            ),
            None,
        )
        if target_family is not None:
            observed[target_family].add(str(target_real))
            expected = expected_real.get(target_family)
            if pinned_calo_reco and expected is not None and target_real != expected:
                failures.append(
                    f"single-provider target mismatch: {target_text}; expected {expected}"
                )
        continue
    if "=> not found" in line:
        failures.append(f"unresolved dependency: {line}")
        continue
    match = re.match(r"(\S+)\s+=>\s+(\S+)\s+", line)
    if not match:
        continue
    name, resolved_text = match.groups()
    resolved = Path(resolved_text).resolve()
    family = next(
        (
            candidate
            for candidate in local_core
            if name == candidate or name.startswith(candidate + ".")
        ),
        None,
    )
    if family is not None:
        observed[family].add(str(resolved))
    if re.match(
        r"^/sphenix/(?:u|user)/[^/]+/"
        r"(?:(?:thesisAnalysis|thesisAnalysis_auau)/)?install/lib(?:64)?/",
        resolved_text,
    ):
        failures.append(f"mutable private install dependency: {name} => {resolved_text}")
    staged_family = name.split(".so", 1)[0] + ".so" if ".so" in name else name
    if (name in staged_names or staged_family in staged_families) and not beneath(resolved, snapshot):
        failures.append(f"staged dependency escaped snapshot: {name} => {resolved_text}")
    if use_release and family is not None:
        if pinned_calo_reco:
            expected = expected_real.get(family)
            if expected is not None and resolved != expected:
                failures.append(
                    f"single-provider resolution mismatch: {name} => {resolved_text}; "
                    f"expected {expected}"
                )
        elif not any(beneath(resolved, root) for root in release_roots):
            failures.append(
                f"release companion escaped pinned release: {name} => {resolved_text}"
            )

if targets == 0:
    failures.append("no frozen loader targets were inspected")
if pinned_calo_reco:
    for family in local_core:
        if not observed[family]:
            failures.append(
                f"selected provider was not observed in frozen-target dependency closure: {family}"
            )
if failures:
    raise SystemExit("\n".join(failures))

if receipt is not None:
    providers = {}
    for family in local_core:
        real = expected_real.get(family)
        if real is None:
            raise SystemExit(f"no selected provider was recorded for {family}")
        providers[family] = {
            "observed_resolutions": sorted(observed[family]),
            "realpath": str(real),
            "sha256": digest(real),
        }
        if family in expected_companions:
            providers[family]["declared_path"] = str(expected_companions[family][0])
    payload = {
        "ldd_report": {
            "path": report.name,
            "sha256": digest(report),
        },
        "mode": mode,
        "pinned_calo_reco_release_companions": pinned_calo_reco,
        "providers": providers,
        "release_roots": [str(root) for root in release_roots],
        "schema": "RJ_SNAPSHOT_LOADER_RECEIPT_V1",
        "snapshot_lib": str(snapshot),
        "status": "PASS",
        "targets_inspected": targets,
    }
    receipt.write_text(
        json.dumps(payload, indent=2, sort_keys=True, separators=(",", ": ")) + "\n",
        encoding="utf-8",
    )
PY
  then
    err "Frozen snapshot dynamic-loader closure failed; inspect ${report}"
    cat "$report" >&2
    (( report_is_temporary )) && rm -f "$report"
    return 2
  fi
  if (( report_is_temporary )); then
    rm -f "$report"
  else
    [[ -s "$receipt" ]] || {
      err "Frozen snapshot loader receipt was not created: ${receipt}"
      return 2
    }
  fi
  return 0
}

create_pipeline_snapshot() {
  local mode="$1"   # pp | auau
  local stamp="$2"

  local snap_dir="${SNAPSHOT_ROOT}/${TAG}_${stamp}"
  local snap_lib_dir="${snap_dir}/lib"
  local user_root="/sphenix/u/${USER:-$(id -u -n)}"
  local pp_library_source="${RJ_PP_LIBRARY_OVERRIDE:-${user_root}/thesisAnalysis/install/lib/libRecoilJets.so}"
  local auau_library_source="${RJ_AUAU_LIBRARY_OVERRIDE:-${user_root}/thesisAnalysis_auau/install/lib/libRecoilJetsAuAu.so}"
  local calo_reco_library_source="${RJ_CALO_RECO_LIBRARY_OVERRIDE:-${user_root}/thesisAnalysis/install/lib/libcalo_reco.so}"
  local photon_cluster_builder_header_source="${RJ_PHOTON_CLUSTER_BUILDER_HEADER_OVERRIDE:-${user_root}/thesisAnalysis/install/include/caloreco/PhotonClusterBuilder.h}"
  local photon_cluster_builder_library_source="${RJ_PHOTON_CLUSTER_BUILDER_LIBRARY_OVERRIDE:-}"

  local live_wrapper=""
  local live_macro=""
  local snap_wrapper=""
  local snap_macro=""

  local snap_impl="${snap_dir}/Fun4All_recoilJets_unified_impl.C"
  local snap_calo="${snap_dir}/Calo_Calib.C"
  local snap_pp_header="${snap_dir}/RecoilJets.h"
  local snap_auau_header="${snap_dir}/RecoilJets_AuAu.h"
  local snap_photon_cluster_builder_header="${snap_dir}/PhotonClusterBuilder.h"
  local use_release_core_libs=0
  local use_pinned_calo_reco_release_companions=0
  local use_ppg12_archived_runtime=0
  local ppg12_canary_manifest=""
  local ppg12_canary_manifest_sha256=""
  local release_core_lib_dir="${RJ_RELEASE_CORE_LIB_DIR:-/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.558/lib}"
  local release_core_lib64_dir="${RJ_RELEASE_CORE_LIB64_DIR:-/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.558/lib64}"
  local pinned_release_name="${RJ_PINNED_RELEASE_NAME:-}"
  local pinned_offline_main="${RJ_PINNED_OFFLINE_MAIN:-}"
  local pinned_calo_reco_soname="${RJ_PINNED_CALO_RECO_SONAME:-}"
  local snapshot_loader_receipt="${snap_dir}/snapshot_loader_receipt.json"

  if [[ "$mode" == "pp" ]] && ppg12_archived_di_campaign_requested; then
    use_ppg12_archived_runtime=1
    use_release_core_libs=1
    release_core_lib_dir="${PPG12_ARCHIVED_OFFLINE_MAIN}/lib"
    release_core_lib64_dir="${PPG12_ARCHIVED_OFFLINE_MAIN}/lib64"
    ppg12_canary_manifest="${RJ_PPG12_DI_ARCHIVED_CANARY_MANIFEST:-}"
    ppg12_canary_manifest_sha256="${RJ_PPG12_DI_ARCHIVED_CANARY_MANIFEST_SHA256:-}"
    validate_ppg12_archived_canary_manifest "$ppg12_canary_manifest" "$ppg12_canary_manifest_sha256" || return $?
  fi

  if env_truthy "${RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS:-0}"; then
    (( use_ppg12_archived_runtime == 0 )) || {
      err "RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS cannot be combined with the archived PPG12 runtime"
      return 2
    }
    ! env_truthy "${RJ_FORCE_RELEASE_CORE_LIBS:-0}" || {
      err "pinned CaloReco plus release companions is mutually exclusive with RJ_FORCE_RELEASE_CORE_LIBS"
      return 2
    }
    [[ ! "${RJ_FORCE_RELEASE_CALO_IO:-0}" =~ ^(1|true|TRUE|yes|YES|on|ON)$ ]] || {
      err "pinned CaloReco plus release companions supersedes RJ_FORCE_RELEASE_CALO_IO"
      return 2
    }
    use_pinned_calo_reco_release_companions=1
    use_release_core_libs=1
    [[ "$pinned_release_name" =~ ^ana\.[0-9]+$ ]] || {
      err "Pinned CaloReco mode requires RJ_PINNED_RELEASE_NAME=ana.NNN"
      return 2
    }
    [[ "$pinned_offline_main" == /*/release/release_ana/"$pinned_release_name" ]] || {
      err "Pinned CaloReco mode requires one exact RJ_PINNED_OFFLINE_MAIN release prefix"
      return 2
    }
    [[ "$pinned_calo_reco_soname" =~ ^libcalo_reco\.so(\.[0-9]+)+$ ]] || {
      err "Pinned CaloReco mode requires RJ_PINNED_CALO_RECO_SONAME"
      return 2
    }
    [[ "$(cd "$release_core_lib_dir" && pwd -P)" == "${pinned_offline_main}/lib" &&
       "$(cd "$release_core_lib64_dir" && pwd -P)" == "${pinned_offline_main}/lib64" ]] || {
      err "Pinned release companion directories do not match RJ_PINNED_OFFLINE_MAIN"
      return 2
    }
  fi

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
  if (( use_pinned_calo_reco_release_companions )); then
    pin_frozen_snapshot_release \
      "$snap_wrapper" "$pinned_release_name" "$pinned_offline_main" || return $?
  fi
  if [[ ! -r "$photon_cluster_builder_header_source" ]]; then
    err "PhotonClusterBuilder snapshot header is missing: ${photon_cluster_builder_header_source}"
    exit 2
  fi
  cp -f "$photon_cluster_builder_header_source" "$snap_photon_cluster_builder_header"
  if [[ -n "$photon_cluster_builder_library_source" ]]; then
    if [[ ! -r "$photon_cluster_builder_library_source" ]]; then
      err "PhotonClusterBuilder override library is missing: ${photon_cluster_builder_library_source}"
      exit 2
    fi
    cp -f "$photon_cluster_builder_library_source" \
      "$snap_lib_dir/libphoton_cluster_builder_override.so"
  fi

  if (( use_ppg12_archived_runtime )); then
    # Reuse only the detector-support payload proven by the accepted full-group
    # canary. Current campaign macro/library/schema/model bytes remain current.
    python3 - "$ppg12_canary_manifest" "$snap_dir" "$snap_lib_dir" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys

manifest = Path(sys.argv[1]).resolve()
snap_dir = Path(sys.argv[2]).resolve()
snap_lib = Path(sys.argv[3]).resolve()
data = json.loads(manifest.read_text(encoding="utf-8"))
artifacts = [row for row in data.get("artifacts", []) if isinstance(row, dict)]
required = {
    "/source/ppg12_runtime/lib/libCaloWaveformSim.so.0.0.0": snap_lib / "libCaloWaveformSim.so.0.0.0",
    "/source/ppg12_runtime/lib/libphg4hit.so.0.0.0": snap_lib / "libphg4hit.so.0.0.0",
    "/source/ppg12_runtime/lib/libg4testbench.so.0.0.0": snap_lib / "libg4testbench.so.0.0.0",
    "/source/ppg12_runtime/include/calowaveformsim/CaloWaveformSim.h": snap_dir / "include/calowaveformsim/CaloWaveformSim.h",
}

def digest(path: Path) -> str:
    h = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

for suffix, destination in required.items():
    rows = [row for row in artifacts if str(row.get("path", "")).endswith(suffix)]
    if len(rows) != 1:
        raise SystemExit(f"accepted canary does not uniquely record {suffix}")
    source = Path(str(rows[0]["path"])).resolve()
    expected = str(rows[0].get("sha256", ""))
    if not source.is_file() or len(expected) != 64 or digest(source) != expected:
        raise SystemExit(f"accepted support artifact is missing or hash-drifted: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)

for stem in ("libCaloWaveformSim", "libphg4hit", "libg4testbench"):
    (snap_lib / f"{stem}.so.0").symlink_to(f"{stem}.so.0.0.0")
    (snap_lib / f"{stem}.so").symlink_to(f"{stem}.so.0.0.0")
PY
  fi

  if (( use_pinned_calo_reco_release_companions )); then
    local pinned_companion_so
    local -a pinned_release_companions=(libcalo_io.so libclusteriso.so libjetbase.so)
    [[ -r "$calo_reco_library_source" ]] || {
      err "Pinned CaloReco runtime library is missing: ${calo_reco_library_source}"
      return 2
    }
    for pinned_companion_so in "${pinned_release_companions[@]}"; do
      if [[ ! -r "${release_core_lib_dir}/${pinned_companion_so}" && ! -r "${release_core_lib64_dir}/${pinned_companion_so}" ]]; then
        err "Pinned CaloReco runtime requires ${pinned_companion_so} in ${release_core_lib_dir} or ${release_core_lib64_dir}"
        return 2
      fi
    done
    cp -f "$calo_reco_library_source" "$snap_lib_dir/libcalo_reco.so"
    say "RJ_PINNED_CALO_RECO_RELEASE_COMPANIONS=1: snapshotting one custom CaloReco and pinning companion libraries to the declared release."
  elif [[ "$mode" != "auau" ]] && { (( use_ppg12_archived_runtime )) || env_truthy "${RJ_FORCE_RELEASE_CORE_LIBS:-0}"; }; then
    use_release_core_libs=1
    local release_core_so
    local -a required_release_core=(libcalo_reco.so libclusteriso.so libjetbase.so)
    (( use_ppg12_archived_runtime )) && required_release_core+=(libcalo_io.so)
    for release_core_so in "${required_release_core[@]}"; do
      if [[ ! -r "${release_core_lib_dir}/${release_core_so}" && ! -r "${release_core_lib64_dir}/${release_core_so}" ]]; then
        err "RJ_FORCE_RELEASE_CORE_LIBS requested, but ${release_core_so} is not readable in ${release_core_lib_dir} or ${release_core_lib64_dir}"
        exit 2
      fi
    done
    say "RJ_FORCE_RELEASE_CORE_LIBS=1: using release CaloReco/ClusterIso/JetBase instead of private core library snapshots."
  else
    if [[ ! -r "$calo_reco_library_source" ]]; then
      err "CaloReco snapshot library is missing: ${calo_reco_library_source}"
      exit 2
    fi
    cp -f "$calo_reco_library_source" "$snap_lib_dir/libcalo_reco.so"
    cp -f "${user_root}/thesisAnalysis/install/lib/libcalo_io.so" "$snap_lib_dir/"
    cp -f "${user_root}/thesisAnalysis/install/lib/libclusteriso.so" "$snap_lib_dir/"
    cp -f "${user_root}/thesisAnalysis/install/lib/libjetbase.so" "$snap_lib_dir/"
  fi
  if [[ -f "$pp_library_source" ]]; then
    cp -f "$pp_library_source" "$snap_lib_dir/libRecoilJets.so"
  elif [[ "$mode" != "auau" ]]; then
    err "p+p snapshot library is missing: ${pp_library_source}"
    exit 2
  fi
  if [[ -f "$auau_library_source" ]]; then
    cp -f "$auau_library_source" "$snap_lib_dir/libRecoilJetsAuAu.so"
  elif [[ "$mode" == "auau" ]]; then
    err "AuAu snapshot library is missing: ${auau_library_source}"
    exit 2
  fi

  # Copy only the ROOT PCM dictionaries defined by the archived libphg4hit.
  # The PPG12 DI support libraries were built against ana.541, so omitting
  # these PCMs prevents ROOT from loading libphg4hit.  Copying every ana.541
  # PCM is also invalid: it makes Cling discover dictionaries whose matching
  # libraries are not in this bounded snapshot and breaks macro materialization.
  if (( use_ppg12_archived_runtime )); then
    local -a archived_ppg12_di_pcm_allowlist=(
      EicEventHeader_Dict_rdict.pcm
      EicEventHeaderv1_Dict_rdict.pcm
      PHG4EventHeader_Dict_rdict.pcm
      PHG4EventHeaderv1_Dict_rdict.pcm
      PHG4HitContainer_Dict_rdict.pcm
      PHG4HitEval_Dict_rdict.pcm
      PHG4Hit_Dict_rdict.pcm
      PHG4Hitv1_Dict_rdict.pcm
      PHG4InEvent_Dict_rdict.pcm
      PHG4Particle_Dict_rdict.pcm
      PHG4Particlev1_Dict_rdict.pcm
      PHG4Particlev2_Dict_rdict.pcm
      PHG4Particlev3_Dict_rdict.pcm
      PHG4Shower_Dict_rdict.pcm
      PHG4Showerv1_Dict_rdict.pcm
      PHG4TruthInfoContainer_Dict_rdict.pcm
      PHG4VtxPoint_Dict_rdict.pcm
      PHG4VtxPointv1_Dict_rdict.pcm
      PHG4VtxPointv2_Dict_rdict.pcm
    )
    local archived_pcm_name
    local archived_pcm_source
    for archived_pcm_name in "${archived_ppg12_di_pcm_allowlist[@]}"; do
      archived_pcm_source=""
      if [[ -r "${release_core_lib_dir}/${archived_pcm_name}" ]]; then
        archived_pcm_source="${release_core_lib_dir}/${archived_pcm_name}"
      elif [[ -r "${release_core_lib64_dir}/${archived_pcm_name}" ]]; then
        archived_pcm_source="${release_core_lib64_dir}/${archived_pcm_name}"
      else
        err "Archived PPG12 DI runtime is missing required libphg4hit PCM: ${archived_pcm_name}"
        exit 2
      fi
      cp -f "$archived_pcm_source" "${snap_lib_dir}/${archived_pcm_name}"
    done
    if [[ "$(find "$snap_lib_dir" -maxdepth 1 -type f -name '*_rdict.pcm' | wc -l)" -ne "${#archived_ppg12_di_pcm_allowlist[@]}" ]]; then
      err "Archived PPG12 DI runtime staged an unexpected ROOT PCM inventory"
      exit 2
    fi
  elif (( use_pinned_calo_reco_release_companions )); then
    # In the single-provider mode, copy dictionaries only from the three
    # immutable campaign build roots.  Never admit mutable private-install
    # dictionaries alongside release-owned CaloReco companions.
    cp -f "$(dirname "$pp_library_source")/"*_rdict.pcm "$snap_lib_dir/" 2>/dev/null || true
    cp -f "$(dirname "$auau_library_source")/"*_rdict.pcm "$snap_lib_dir/" 2>/dev/null || true
    cp -f "$(dirname "$calo_reco_library_source")/"*_rdict.pcm "$snap_lib_dir/" 2>/dev/null || true
  else
    cp -f "${user_root}/thesisAnalysis/install/lib/"*_rdict.pcm "$snap_lib_dir/" 2>/dev/null || true
    cp -f "${user_root}/thesisAnalysis_auau/install/lib/"*_rdict.pcm "$snap_lib_dir/" 2>/dev/null || true
    if [[ -n "${RJ_CALO_RECO_LIBRARY_OVERRIDE:-}" ]]; then
      cp -f "$(dirname "$calo_reco_library_source")/"*_rdict.pcm "$snap_lib_dir/" 2>/dev/null || true
    fi
  fi

  # Preserve dynamic-loader identity inside the frozen snapshot. The copied
  # files are bare linker names, while dependent ELFs request their SONAMEs.
  # Pinned CaloReco fails closed unless its exact versioned loader alias is a
  # symlink to the one snapshotted inode.
  stage_snapshot_soname_aliases \
    "$snap_lib_dir" \
    "$use_pinned_calo_reco_release_companions" \
    "$pinned_calo_reco_soname" || return $?

  sed -i "s|#include \"/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Fun4All_recoilJets_unified_impl.C\"|#include \"${snap_impl}\"|" "$snap_macro"
  sed -i "s|#include \"/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Calo_Calib.C\"|#include \"${snap_calo}\"|" "$snap_impl"
  sed -i "s|#include \"/sphenix/u/patsfan753/scratch/thesisAnalysis/src/RecoilJets.h\"|#include \"${snap_pp_header}\"|" "$snap_impl"
  sed -i "s|#include \"/sphenix/u/patsfan753/scratch/thesisAnalysis/src_AuAu/RecoilJets_AuAu.h\"|#include \"${snap_auau_header}\"|" "$snap_impl"
  sed -i "s|#include \"/sphenix/u/patsfan753/thesisAnalysis/install/include/caloreco/PhotonClusterBuilder.h\"|#include \"${snap_photon_cluster_builder_header}\"|" "$snap_impl"

  if (( use_pinned_calo_reco_release_companions )); then
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libcalo_reco.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_io.so)|R__LOAD_LIBRARY(libcalo_io.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)|R__LOAD_LIBRARY(libclusteriso.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)|R__LOAD_LIBRARY(libjetbase.so)|" "$snap_impl"
  elif (( use_release_core_libs )); then
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)|R__LOAD_LIBRARY(libcalo_reco.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_io.so)|R__LOAD_LIBRARY(libcalo_io.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)|R__LOAD_LIBRARY(libclusteriso.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)|R__LOAD_LIBRARY(libjetbase.so)|" "$snap_impl"
  else
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libcalo_reco.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_io.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libcalo_io.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libclusteriso.so)|" "$snap_impl"
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libjetbase.so)|" "$snap_impl"
  fi
  sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libRecoilJets.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libRecoilJets.so)|" "$snap_impl"
  sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis_auau/install/lib/libRecoilJetsAuAu.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libRecoilJetsAuAu.so)|" "$snap_impl"
  if (( use_pinned_calo_reco_release_companions )); then
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libcalo_reco.so)|" "$snap_calo"
  elif (( use_release_core_libs )); then
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)|R__LOAD_LIBRARY(libcalo_reco.so)|" "$snap_calo"
  else
    sed -i "s|R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)|R__LOAD_LIBRARY(${snap_lib_dir}/libcalo_reco.so)|" "$snap_calo"
  fi
  if [[ "${RJ_FORCE_RELEASE_CALO_IO:-0}" =~ ^(1|true|TRUE|yes|YES|on|ON)$ ]]; then
    local release_calo_io="${RJ_RELEASE_CALO_IO_PATH:-/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/ana.558/lib/libcalo_io.so.0}"
    if [[ ! -r "$release_calo_io" ]]; then
      err "RJ_FORCE_RELEASE_CALO_IO requested, but libcalo_io is not readable: ${release_calo_io}"
      exit 2
    fi
    # Keep ROOT to a single calo_io dictionary provider. Preloading the release
    # library while the frozen macro R__LOAD_LIBRARYs the private snapshot copy
    # duplicates classes and can abort before event processing starts.
    sed -i "s|R__LOAD_LIBRARY(${snap_lib_dir}/libcalo_io.so)|R__LOAD_LIBRARY(${release_calo_io})|" "$snap_impl"
  fi

  local release_prefix=""
  if (( use_release_core_libs )); then
    release_prefix=":${release_core_lib64_dir}:${release_core_lib_dir}"
  fi
  normalize_frozen_snapshot_wrapper "$snap_wrapper" "$release_prefix"

  chmod +x "$snap_wrapper"
  if ! bash -n "$snap_wrapper"; then
    err "Frozen wrapper failed bash -n: ${snap_wrapper}"
    exit 2
  fi
  if ! grep -Eq '^[[:space:]]*rc=125([[:space:]]|$)' "$snap_wrapper"; then
    err "Frozen wrapper lacks the rc sentinel required to prevent rc-unbound holds: ${snap_wrapper}"
    exit 2
  fi
  validate_frozen_snapshot_wrapper_contract "$snap_wrapper" || return $?
  if (( use_pinned_calo_reco_release_companions )); then
    [[ "$(grep -Ec '^[[:space:]]*# RJ_PINNED_SPHENIX_RELEASE_V1[[:space:]]*$' "$snap_wrapper" || true)" == 1 &&
       "$(grep -Ec "^[[:space:]]*source /opt/sphenix/core/bin/sphenix_setup\\.sh -n ${pinned_release_name}[[:space:]]*$" "$snap_wrapper" || true)" == 1 &&
       "$(grep -Fc "$pinned_offline_main" "$snap_wrapper" || true)" -ge 2 ]] || {
      err "Frozen wrapper lacks the exact pinned ${pinned_release_name} runtime witness"
      return 2
    }
  fi
  validate_snapshot_loader_closure \
    "$snap_lib_dir" \
    "$mode" \
    "$use_release_core_libs" \
    "$release_core_lib64_dir" \
    "$release_core_lib_dir" \
    "$use_pinned_calo_reco_release_companions" \
    "$(
      if (( use_pinned_calo_reco_release_companions )); then
        printf '%s' "$snapshot_loader_receipt"
      fi
    )" || return $?
  if (( use_pinned_calo_reco_release_companions )); then
    BULK_FROZEN_SNAPSHOT_LOADER_RECEIPT="$snapshot_loader_receipt"
    BULK_FROZEN_SNAPSHOT_LOADER_RECEIPT_SHA256="$(
      snapshot_sha256_file "$snapshot_loader_receipt"
    )"
  else
    BULK_FROZEN_SNAPSHOT_LOADER_RECEIPT=""
    BULK_FROZEN_SNAPSHOT_LOADER_RECEIPT_SHA256=""
  fi

  if (( use_ppg12_archived_runtime )); then
    local accepted_manifest_copy="${snap_dir}/accepted_full_group_canary_manifest.json"
    local runtime_manifest="${snap_dir}/ppg12_archived_runtime_manifest.sha256"
    cp -f "$ppg12_canary_manifest" "$accepted_manifest_copy"
    (
      cd "$snap_dir"
      sha256sum \
        accepted_full_group_canary_manifest.json \
        RecoilJets_Condor.sh \
        Fun4All_recoilJets.C \
        Fun4All_recoilJets_unified_impl.C \
        Calo_Calib.C \
        RecoilJets.h \
        lib/libRecoilJets.so \
        lib/libCaloWaveformSim.so.0.0.0 \
        lib/libCaloWaveformSim.so.0 \
        lib/libphg4hit.so.0.0.0 \
        lib/libphg4hit.so.0 \
        lib/libg4testbench.so.0.0.0 \
        lib/libg4testbench.so.0 \
        include/calowaveformsim/CaloWaveformSim.h \
        > "$runtime_manifest"
      find lib -maxdepth 1 -type f -name '*_rdict.pcm' -print0 \
        | sort -z \
        | xargs -0 -r sha256sum \
        >> "$runtime_manifest"
    )
    [[ -s "$runtime_manifest" ]] || {
      err "Archived PPG12 DI runtime manifest was not created: ${runtime_manifest}"
      return 96
    }
    BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST="$runtime_manifest"
    BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST_SHA256="$(sha256sum "$runtime_manifest" | awk '{print $1}')"
  else
    BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST=""
    BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST_SHA256=""
  fi

  write_and_seal_snapshot_manifest "$snap_dir" || return $?
  BULK_FROZEN_EXE="$snap_wrapper"
  BULK_FROZEN_MACRO="$snap_macro"

  say "Pipeline snapshot created:"
  say "  snapshot dir : ${snap_dir}"
  say "  frozen exe   : ${BULK_FROZEN_EXE}"
  say "  frozen macro : ${BULK_FROZEN_MACRO}"
  if [[ "$mode" == "auau" ]]; then
    say "  AuAu library : ${auau_library_source}"
  fi
  return 0
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
  ${BOLD}$0 <isPP|isPPrun25|isAuAu|isOO> condor all [groupSize N]${RST}
  ${DIM}$0 <isPP|isPPrun25|isAuAu|isOO> condor allFromScratch [groupSize N]          # alias for direct fanout${RST}

${BOLD}SIM mode (photonJet10/20 merged):${RST}
  ${BOLD}$0 isSim local [Nevents] [VERBOSE=N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim checkModels [VERBOSE=N]${RST}
  ${BOLD}$0 isSim CHECKJOBS [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorTest [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorDoAllSmoke [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${BOLD}$0 isSim condorDoAll [groupSize N] [SAMPLE=run28_photonjet10]${RST}
  ${DIM}$0 isSim condorDoAllFromScratch [groupSize N] [SAMPLE=run28_photonjet10]  # alias for direct fanout${RST}

${BOLD}SIM mode (embedded photon AuAu-like signal):${RST}
  ${BOLD}$0 isSimEmbedded local [Nevents] [VERBOSE=N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded CHECKJOBS [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded condorDoAllSmoke [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${BOLD}$0 isSimEmbedded condorDoAll [groupSize N] [SAMPLE=run28_embeddedPhoton20]${RST}
  ${DIM}$0 isSimEmbedded condorDoAllFromScratch [groupSize N] [SAMPLE=run28_embeddedPhoton20]  # alias for direct fanout${RST}

${BOLD}SIM mode (embedded inclusive-jet AuAu-like background):${RST}
  ${BOLD}$0 isSimEmbeddedInclusive local [Nevents] [VERBOSE=N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive CHECKJOBS [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive condorDoAllSmoke [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}
  ${BOLD}$0 isSimEmbeddedInclusive condorDoAll [groupSize N] [SAMPLE=run28_embeddedJet20]${RST}
  ${DIM}$0 isSimEmbeddedInclusive condorDoAllFromScratch [groupSize N] [SAMPLE=run28_embeddedJet20]  # alias for direct fanout${RST}

${BOLD}SIM mode (inclusive jet 5/8/12/20/30/40):${RST}
  ${BOLD}$0 isSimInclusive local [Nevents] [VERBOSE=N] [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimInclusive CHECKJOBS [groupSize N] [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimInclusive condorTest [SAMPLE=run28_jet5]${RST}
  ${BOLD}$0 isSimInclusive condorDoAll [groupSize N] [SAMPLE=run28_jet5]${RST}
  ${DIM}$0 isSimJet5 ...  # legacy alias for isSimInclusive${RST}

${BOLD}SIM mode (MinBias DETROIT):${RST}
  ${BOLD}$0 isSimMB local [Nevents] [VERBOSE=N] [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB CHECKJOBS [groupSize N] [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB condorTest [SAMPLE=run28_detroit]${RST}
  ${BOLD}$0 isSimMB condorDoAll [groupSize N] [SAMPLE=run28_detroit]${RST}

${BOLD}AuAu embedded photon-ID BDT training:${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT localTest [Nevents] [NFILES=N]${RST}  ${DIM}# forwards to scripts/auau_tight_bdt_pipeline.sh${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT smokeTest [groupSize N]${RST}         ${DIM}# sidecar extraction DAG${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT trainFromExtraction SOURCE=/path${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT trainExpandedFromExtraction SOURCE=/path [PLAN_ONLY=1]${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT trainExpandedFromExtractionCondor SOURCE=/path [groupSize N]${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT applyCheck MODEL_DIR=/path/to/tight/models${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT validateOnSim SOURCE=/path MODEL_DIR=/path${RST}
  ${BOLD}$0 isSimEmbeddedAndInclusive trainTightBDT validateOnSimCondor SOURCE=/path MODEL_DIR=/path [groupSize N]${RST}
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
  $0 isPP  condor all groupSize 4
  $0 isPP  CHECKJOBS groupSize 4

  $0 isSim CHECKJOBS groupSize 5
  $0 isSim local 5000
  $0 isSim checkModels
  $0 isSim condorTest
  $0 isSim condorDoAllSmoke
  $0 isSimEmbedded condorDoAllSmoke
  $0 isSimEmbeddedInclusive condorDoAllSmoke
  $0 isSim condorDoAll groupSize 5

  $0 isSimJet5 local 5000
  $0 isSimJet5 condorDoAll groupSize 5
  $0 isSimMB local 5000
  $0 isSimMB condorDoAll groupSize 5

  RJ_CLEAN_OUTPUT_BASE=1 $0 isAuAu condor all groupSize 7
  RJ_CLEAN_OUTPUT_BASE=1 $0 isSimEmbedded condorDoAll groupSize 7

Production histogram jobs use the original RecoilJets modules with direct
photon-ID fanout plus internal jet-pT, dphi, and iso/cone histogram views.
Use RJ_ID_FANOUT_MAX_ROWS=5 or 10 only if a smoke test shows memory pressure.
Pool/replay histogram production was removed.
Automatic production DAGs send one final READY/CHECK/FAILED email by default;
nested merge DAGs run quiet with strict output validation. Use RJ_AUTO_MERGE=0
for analysis-only submission, or run mergeRecoilJets.sh manually for per-cfg
debug emails.
Production output cleanup is explicit via RJ_CLEAN_OUTPUT_BASE=1 or
./scripts/recoiljets_cleanup.sh dryrun|apply <target>.
USAGE
  exit 2
}

need_cmd() { command -v "$1" >/dev/null 2>&1 || { err "Missing required command: $1"; exit 3; }; }

cleanup_helper_available() {
  [[ -x "$CLEANUP_HELPER" ]]
}

cleanup_apply_or_dryrun_mode() {
  if dag_dryrun_enabled; then
    printf '%s\n' "dryrun"
  else
    printf '%s\n' "apply"
  fi
}

cleanup_dataset_outputs_if_requested() {
  [[ "${RJ_CLEAN_OUTPUT_BASE:-0}" == "1" ]] || return 0
  cleanup_helper_available || { err "RJ_CLEAN_OUTPUT_BASE=1 requested but cleanup helper is missing/executable bit is off: ${CLEANUP_HELPER}"; exit 66; }
  local clean_mode
  clean_mode="$(cleanup_apply_or_dryrun_mode)"
  say "RJ_CLEAN_OUTPUT_BASE=1 → ${clean_mode} cleanup for dataset ${DATASET}"
  "$CLEANUP_HELPER" "$clean_mode" dataset "$DATASET"
}

cleanup_local_artifacts_before_local_run() {
  cleanup_helper_available || {
    warn "Cleanup helper missing/executable bit is off; skipping automatic local artifact cleanup: ${CLEANUP_HELPER}"
    return 0
  }
  say "Cleaning disposable *_LOCAL* artifacts before local test"
  "$CLEANUP_HELPER" apply local
}

write_cleanup_manifest() {
  local manifest="$1"
  shift || true
  mkdir -p "$(dirname "$manifest")"
  {
    echo "RECOILJETS_CLEANUP_MANIFEST_V1"
    echo "created=$(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "dataset=${DATASET}"
    echo "tag=${TAG}"
    echo "action=${ACTION:-}"
    echo "dest_base=${DEST_BASE:-}"
    echo "stage_dir=${STAGE_DIR:-}"
    echo "yaml_override_dir=${SIM_YAML_OVERRIDE_DIR:-}"
    echo "sub_dir=${SUB_DIR:-}"
    local item
    for item in "$@"; do
      echo "artifact=${item}"
    done
  } > "$manifest"
}

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

memory_request_from_env_or_default() {
  local default_request="$1"
  if [[ -n "${RJ_REQUEST_MEMORY:-}" ]]; then
    printf '%s\n' "$RJ_REQUEST_MEMORY"
    return 0
  fi
  if [[ -n "${RJ_REQUEST_MEMORY_MB:-}" ]]; then
    case "$RJ_REQUEST_MEMORY_MB" in
      *[Mm][Bb]|*[Gg][Bb]|*[Gg])
        printf '%s\n' "$RJ_REQUEST_MEMORY_MB"
        ;;
      *)
        printf '%sMB\n' "$RJ_REQUEST_MEMORY_MB"
        ;;
    esac
    return 0
  fi
  printf '%s\n' "$default_request"
}

condor_auto_memory_retry_block() {
  local base_mb="${1:-0}"
  [[ "$base_mb" =~ ^[0-9]+$ ]] || base_mb=0
  local cap_mb="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-8000}"
  local retries="${RJ_AUTO_MEMORY_RETRY_MAX_RELEASES:-2}"
  [[ "$cap_mb" =~ ^[0-9]+$ ]] || cap_mb=8000
  [[ "$retries" =~ ^[0-9]+$ ]] || retries=2
  if [[ "${RJ_AUTO_MEMORY_RETRY:-1}" == "0" || "$base_mb" -le 0 || "$base_mb" -ge "$cap_mb" || "$retries" -le 0 ]]; then
    printf 'request_memory= %s\n' "$base_mb"
    return 0
  fi
  cat <<EOT
+RJBaseRequestMemoryMb = ${base_mb}
+RJMemoryRetryCapMb = ${cap_mb}
+RJMemoryRetryMaxReleases = ${retries}
request_memory = ifThenElse(MemoryUsage =!= undefined, ifThenElse(int(MemoryUsage * 1.35 + 512) > ${cap_mb}, ${cap_mb}, ifThenElse(int(MemoryUsage * 1.35 + 512) > ${base_mb}, int(MemoryUsage * 1.35 + 512), ${base_mb})), ${base_mb})
periodic_release = (JobStatus == 5) && (NumJobStarts <= ${retries}) && (RequestMemory < ${cap_mb}) && ((HoldReasonCode == 34) || ((HoldReason =!= undefined) && regexp("memory|Memory|cgroup|request_memory|request memory", HoldReason)))
EOT
}

condor_late_materialization_block() {
  local max_materialize="${RJ_CONDOR_MAX_MATERIALIZE:-0}"
  local max_idle="${RJ_CONDOR_MAX_IDLE:-0}"
  [[ "$max_materialize" =~ ^[0-9]+$ ]] || {
    err "RJ_CONDOR_MAX_MATERIALIZE must be a non-negative integer"
    return 2
  }
  [[ "$max_idle" =~ ^[0-9]+$ ]] || {
    err "RJ_CONDOR_MAX_IDLE must be a non-negative integer"
    return 2
  }
  if (( max_materialize == 0 && max_idle == 0 )); then
    return 0
  fi
  if (( max_materialize < 1 || max_materialize > 256 )); then
    err "RJ_CONDOR_MAX_MATERIALIZE must be in [1,256] when enabled"
    return 2
  fi
  if (( max_idle < 1 || max_idle > max_materialize )); then
    err "RJ_CONDOR_MAX_IDLE must be in [1,RJ_CONDOR_MAX_MATERIALIZE]"
    return 2
  fi
  printf 'max_materialize = %s\nmax_idle = %s\n' \
    "$max_materialize" "$max_idle"
}

condor_worker_failure_hold_block() {
  [[ "${RJ_HOLD_FAILED_WORKERS:-1}" == "0" ]] && return 0
  cat <<'EOT'
# Keep failed workers visible for diagnosis instead of letting them silently finish.
on_exit_hold = (ExitBySignal == True) || (ExitCode =!= 0)
on_exit_hold_reason = "RecoilJets worker exited nonzero or by signal; inspect stdout/err before release/remove"
EOT
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

env_truthy() {
  case "${1:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
  esac
  return 1
}

dataset_is_sim_like() {
  case "${1:-${DATASET:-}}" in
    isSim|isSimEmbedded|isSimEmbeddedInclusive|isSimEmbeddedAndInclusive|isSimInclusive|isSimJet5|isSimMB) return 0 ;;
  esac
  return 1
}

append_submit_extra_env_var() {
  local extra="$1"
  local key="$2"
  local val="${!key:-}"
  [[ -n "$val" ]] || { printf '%s' "$extra"; return; }
  case ";${extra};" in
    *";${key}="*) printf '%s' "$extra"; return ;;
  esac
  printf '%s' "${extra:+${extra};}${key}=${val}"
}

remove_submit_extra_env_var() {
  local extra="$1"
  local key="$2"
  local out=""
  local old_ifs="$IFS"
  local part part_key
  IFS=';'
  for part in $extra; do
    [[ -n "$part" ]] || continue
    part_key="${part%%=*}"
    [[ "$part_key" == "$key" ]] && continue
    out="${out:+${out};}${part}"
  done
  IFS="$old_ifs"
  printf '%s' "$out"
}

append_submit_extra_env_literal() {
  local extra="$1"
  local key="$2"
  local value="$3"
  printf '%s' "${extra:+${extra};}${key}=${value}"
}

finalize_ppg12_archived_di_submit_env() {
  local extra="$1"
  local dataset="${2:-${DATASET:-}}"
  local sample="${3:-${SIM_SAMPLE:-}}"
  local key

  extra="$(remove_submit_extra_env_var "$extra" RJ_PPG12_DI_ARCHIVED_RECO_CHAIN)"
  if ! ppg12_archived_di_lane "$dataset" "$sample"; then
    append_submit_extra_env_literal "$extra" RJ_PPG12_DI_ARCHIVED_RECO_CHAIN 0
    return 0
  fi

  for key in \
    RJ_TRUTH_JETS_MODE \
    RJ_PPG12_DI_ARCHIVED_RELEASE \
    RJ_PPG12_DI_RUNTIME_MANIFEST \
    RJ_PPG12_DI_RUNTIME_MANIFEST_SHA256 \
    RJ_PPG12_PERIOD_ALLOW_ALL_SIM \
    RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE \
    RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE \
    RJ_PPG12_PERIOD_USE_LUMI_WEIGHT \
    RJ_PPG12_PHOTON_YIELD \
    RJ_PPG12_PHOTON_YIELD_DOUBLE \
    RJ_PPG12_PERIOD_STRICT_DI \
    RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4 \
    RJ_PPG12_PPSIM_G4_ONLY \
    RJ_SIM_ALLOW_NONE_LISTS; do
    extra="$(remove_submit_extra_env_var "$extra" "$key")"
  done

  [[ -n "$BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST" &&
     "$BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST" == /* &&
     -s "$BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST" ]] || {
    err "Archived PPG12 DI submission has no immutable runtime manifest from create_pipeline_snapshot."
    return 96
  }
  [[ "$BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST_SHA256" =~ ^[0-9a-f]{64}$ ]] || {
    err "Archived PPG12 DI submission has no valid runtime-manifest digest."
    return 96
  }

  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_DI_ARCHIVED_RECO_CHAIN 1)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_TRUTH_JETS_MODE DST)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_DI_ARCHIVED_RELEASE "$PPG12_ARCHIVED_OFFLINE_MAIN")"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_DI_RUNTIME_MANIFEST "$BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST")"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_DI_RUNTIME_MANIFEST_SHA256 "$BULK_PPG12_ARCHIVED_RUNTIME_MANIFEST_SHA256")"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PERIOD_ALLOW_ALL_SIM 0)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE 0)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE 0)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PERIOD_USE_LUMI_WEIGHT 1)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PHOTON_YIELD 1)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PHOTON_YIELD_DOUBLE 1)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PERIOD_STRICT_DI 1)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4 1)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_PPG12_PPSIM_G4_ONLY 1)"
  extra="$(append_submit_extra_env_literal "$extra" RJ_SIM_ALLOW_NONE_LISTS 1)"
  printf '%s' "$extra"
}

submit_extra_env_var_is_truthy() {
  local extra="$1"
  local key="$2"
  local old_ifs="$IFS"
  local part part_key part_value
  IFS=';'
  for part in $extra; do
    [[ -n "$part" ]] || continue
    part_key="${part%%=*}"
    [[ "$part_key" == "$key" ]] || continue
    part_value="${part#*=}"
    IFS="$old_ifs"
    env_truthy "$part_value"
    return
  done
  IFS="$old_ifs"
  return 1
}

ppg12_period_sim_uses_auto_mix_weight() {
  dataset_is_sim_like "${DATASET:-}" || return 1
  env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}" || return 1
  [[ -n "${RJ_PPG12_PERIOD:-}" ]] || return 1
  env_truthy "${RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE:-0}" && return 1
  [[ -n "${RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT:-}" || "${RJ_SUBMIT_EXTRA_ENV:-}" == *RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT=* ]] || return 1
  return 0
}

effective_require_non_tiny_output() {
  if env_truthy "${RJ_REQUIRE_NON_TINY_OUTPUT:-0}"; then
    printf '1'
    return
  fi
  if env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}" \
    || env_truthy "${RJ_PPG12_TABLE_QA:-0}" \
    || env_truthy "${RJ_PPG12_FIG11_SB_DIAGNOSTIC:-0}" \
    || env_truthy "${RJ_PPG12_FIG13_PARITY_QA:-0}" \
    || env_truthy "${RJ_PPG12_FIG7_TRIGGER_DIAGNOSTIC:-0}" \
    || env_truthy "${RJ_PPG12_FIG13_BIT30_DIAGNOSTIC:-0}"; then
    printf '1'
    return
  fi
  printf '0'
}

build_submit_extra_env_fragment() {
  local extra="${RJ_SUBMIT_EXTRA_ENV:-}"
  local sim_dataset=0
  dataset_is_sim_like "${DATASET:-}" && sim_dataset=1
  if (( ! sim_dataset )) && [[ "${DATASET:-}" == "isPP" ]] \
    && [[ -z "${RJ_PPG12_FIG13_BIT30_DIAGNOSTIC+x}" ]] \
    && { env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}" \
      || env_truthy "${RJ_PPG12_TABLE_QA:-0}" \
      || env_truthy "${RJ_PPG12_FIG13_PARITY_QA:-0}"; }; then
    # PPG12 pp-data parity uses direct GL1 ScaledVector bit 30 as the
    # canonical Photon 4 GeV + MBD N&S trigger, not the TriggerAnalyzer alias.
    RJ_PPG12_FIG13_BIT30_DIAGNOSTIC=1
  fi
  # Load-bearing worker-side analysis modes must reach Condor workers, not
  # only the submit shell. RJ_SUBMIT_EXTRA_ENV remains the general override.
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PHOTON_YIELD)"
  if (( sim_dataset )); then
    if env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}" &&
       { env_truthy "${RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX:-0}" ||
         env_truthy "${RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX:-0}" ||
         env_truthy "${RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX:-0}" ||
         submit_extra_env_var_is_truthy "$extra" RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX ||
         submit_extra_env_var_is_truthy "$extra" RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX ||
         submit_extra_env_var_is_truthy "$extra" RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX; }; then
      printf '%s\n' \
        "ERROR: PPG12 truth-for-reconstructed-object vertex mode is forbidden; use reconstructed MBD z for reco kinematics/Eiso and truth z only for SI/DI weights." >&2
      return 97
    fi
    extra="$(remove_submit_extra_env_var "$extra" RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX)"
    extra="$(remove_submit_extra_env_var "$extra" RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX)"
    extra="$(remove_submit_extra_env_var "$extra" RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PHOTON_YIELD_DOUBLE)"
    if ppg12_period_sim_uses_auto_mix_weight; then
      extra="$(remove_submit_extra_env_var "$extra" RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT)"
    else
      extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT)"
    fi
  fi
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_CROSSING_PERIOD)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_FILTER_DATA)"
  if (( sim_dataset )); then
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_USE_LUMI_WEIGHT)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_STRICT_DI)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_ALLOW_ALL_SIM)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_ALLOW_MIX_OVERRIDE)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PERIOD_ALLOW_VERTEX_FILE_OVERRIDE)"
  fi
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_FIG8_BUILD_NOSPLIT)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_FIG8_CLUSTER_NODE)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_FIG8_FALLBACK_TO_SPLIT)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_MC_ISO_SCALE)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_MC_ISO_SHIFT)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_MBD_T0_CORRECTION_FILE)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_NPB_DATA_TAGGING)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_NPB_TIME_SAMPLE_NS)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_NPB_DELTA_T_CUT)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_NPB_WETA_MIN)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_NPB_AWAY_JET_PT_MIN)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_TABLE_QA_NPB_AWAY_JET_DPHI_MIN)"
  if (( sim_dataset )); then
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_FIG11_SB_DIAGNOSTIC)"
  fi
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_FIG13_PARITY_QA)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_FIG7_TRIGGER_DIAGNOSTIC)"
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_FIG13_BIT30_DIAGNOSTIC)"
  if (( sim_dataset )); then
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_CLOSURE_CANARY)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_CLOSURE_CANARY_ID)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PPSIM_REPLAY_SEEDS)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE)"
    extra="$(append_submit_extra_env_var "$extra" RJ_DISABLE_JES_CDB_AUDIT)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PPSIM_G4_ONLY)"
    extra="$(append_submit_extra_env_var "$extra" RJ_SIM_ALLOW_NONE_LISTS)"
    extra="$(append_submit_extra_env_var "$extra" RJ_FORCE_RELEASE_CORE_LIBS)"
    extra="$(append_submit_extra_env_var "$extra" RJ_RELEASE_CORE_LIB_DIR)"
    extra="$(append_submit_extra_env_var "$extra" RJ_RELEASE_CORE_LIB64_DIR)"
    extra="$(append_submit_extra_env_var "$extra" RJ_FORCE_RELEASE_CALO_IO)"
    extra="$(append_submit_extra_env_var "$extra" RJ_RELEASE_CALO_IO_PATH)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PP_VERTEX_REWEIGHT_FILE)"
    extra="$(append_submit_extra_env_var "$extra" RJ_PP_VERTEX_REWEIGHT_HIST)"
  fi
  extra="$(append_submit_extra_env_var "$extra" RJ_PPG12_PP_DATA_PAIRED)"
  # A fixed reconstructed AuAu-isolation study is a guarded opt-in. If the
  # submit-side contract accepts that opt-in, propagate the same flag to the
  # worker so the C++ hard stop evaluates the identical authorization state.
  extra="$(append_submit_extra_env_var "$extra" RJ_ALLOW_FIXED_RECO_ISO_VIEWS)"
  if (( sim_dataset )); then
    extra="$(finalize_ppg12_archived_di_submit_env "$extra" "${DATASET:-}" "${SIM_SAMPLE:-}")" || return $?
  fi
  [[ -n "$extra" && "$extra" != \;* ]] && extra=";${extra}"
  printf '%s' "$extra"
}

dataset_default_vz_cut() {
  case "${DATASET:-}" in
    isAuAu|isOO|isSimEmbedded|isSimEmbeddedInclusive|isSimEmbeddedAndInclusive)
      printf '%s\n' "10"
      ;;
    isPP|isPPrun25|isSim|isSimJet5|isSimInclusive|isSimMB|"")
      printf '%s\n' "60"
      ;;
    *)
      printf '%s\n' "60"
      ;;
  esac
}

vz_value_in_list() {
  local needle="$1"; shift || true
  local v
  for v in "$@"; do
    sim_is_close "$needle" "$v" && return 0
  done
  return 1
}

vz_selection_reason() {
  if env_truthy "${RJ_VZ_SCAN_ALL:-0}"; then
    printf '%s\n' "full YAML scan via RJ_VZ_SCAN_ALL=1"
  elif [[ -n "${RJ_VZ_CUT_OVERRIDE:-}" ]]; then
    printf '%s\n' "explicit override via RJ_VZ_CUT_OVERRIDE=${RJ_VZ_CUT_OVERRIDE}"
  else
    printf '%s\n' "dataset default for ${DATASET:-unknown}: $(dataset_default_vz_cut) cm"
  fi
}

vz_selection_reason_for_yaml() {
  local yaml="$1"
  local -a yaml_vzs=()
  mapfile -t yaml_vzs < <( yaml_get_values "vz_cut_cm" "$yaml" )
  if env_truthy "${RJ_VZ_SCAN_ALL:-0}"; then
    printf '%s\n' "full YAML scan via RJ_VZ_SCAN_ALL=1"
  elif [[ -n "${RJ_VZ_CUT_OVERRIDE:-}" ]]; then
    printf '%s\n' "explicit override via RJ_VZ_CUT_OVERRIDE=${RJ_VZ_CUT_OVERRIDE}"
  elif (( ${#yaml_vzs[@]} == 1 )); then
    printf '%s\n' "scalar YAML value"
  else
    printf '%s\n' "dataset default for ${DATASET:-unknown}: $(dataset_default_vz_cut) cm"
  fi
}

dataset_vz_values() {
  local yaml="$1"
  local -a yaml_vzs=()
  mapfile -t yaml_vzs < <( yaml_get_values "vz_cut_cm" "$yaml" )
  (( ${#yaml_vzs[@]} )) || { err "No values found for vz_cut_cm in $yaml"; exit 72; }

  if env_truthy "${RJ_VZ_SCAN_ALL:-0}"; then
    printf '%s\n' "${yaml_vzs[@]}"
    return 0
  fi

  local selected="${RJ_VZ_CUT_OVERRIDE:-}"
  if [[ -z "$selected" ]]; then
    if (( ${#yaml_vzs[@]} == 1 )); then
      selected="${yaml_vzs[0]}"
    else
      selected="$(dataset_default_vz_cut)"
    fi
  fi
  selected="$(trim_ws "$selected")"
  [[ -n "$selected" ]] || { err "Selected vz_cut_cm is empty for dataset=${DATASET:-unknown}"; exit 72; }
  if ! [[ "$selected" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
    err "Selected vz_cut_cm must be numeric, got '${selected}'"
    exit 72
  fi

  if ! vz_value_in_list "$selected" "${yaml_vzs[@]}"; then
    if env_truthy "${RJ_ALLOW_VZ_CUT_OVERRIDE_OUTSIDE_YAML:-0}"; then
      warn "Selected vz_cut_cm=${selected} is not listed in ${yaml}; allowing because RJ_ALLOW_VZ_CUT_OVERRIDE_OUTSIDE_YAML=1"
    else
      err "Selected vz_cut_cm=${selected} is not in YAML vz_cut_cm list [${yaml_vzs[*]}] for dataset=${DATASET:-unknown}"
      err "Use RJ_VZ_SCAN_ALL=1 to scan YAML values, RJ_VZ_CUT_OVERRIDE=<listed value>, or RJ_ALLOW_VZ_CUT_OVERRIDE_OUTSIDE_YAML=1 for an explicit non-YAML value."
      exit 72
    fi
  fi
  printf '%s\n' "$selected"
}

say_vz_selection_summary() {
  local yaml="$1"; shift
  local -a selected=( "$@" )
  local -a available=()
  mapfile -t available < <( yaml_get_values "vz_cut_cm" "$yaml" )
  say "  vz YAML available                  : [${available[*]}]"
  say "  vz selected for production         : [${selected[*]}]  ($(vz_selection_reason_for_yaml "$yaml"))"
}

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

jetpt_internal_enabled() {
  case "${RJ_DISABLE_JET_PT_INTERNALIZATION:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 1 ;;
  esac
  case "${RJ_INTERNALIZE_JET_PT:-1}" in
    0|false|FALSE|no|NO|off|OFF) return 1 ;;
  esac
  return 0
}

jetpt_scan_csv() {
  local out="" v
  for v in "$@"; do
    v="$(trim_ws "$v")"
    [[ -z "$v" ]] && continue
    out="${out:+${out},}${v}"
  done
  printf '%s\n' "$out"
}

jetpt_nominal_value() {
  (( "$#" > 0 )) && printf '%s\n' "$1"
}

jetpt_submit_values() {
  if jetpt_internal_enabled && (( "$#" > 0 )); then
    jetpt_nominal_value "$@"
  else
    printf '%s\n' "$@"
  fi
}

jetpt_tag_component() {
  local pt="$1"
  if jetpt_internal_enabled; then
    echo "jetMinPtScan"
  else
    echo "jetMinPt$(sim_pt_tag "$pt")"
  fi
}

jetpt_env_fragment() {
  if ! jetpt_internal_enabled; then
    return 0
  fi
  local csv
  csv="$(jetpt_scan_csv "$@")"
  [[ -n "$csv" ]] && printf ';RJ_INTERNAL_JET_PT_MINS=%s' "$csv"
}

embedded_inclusive_stitch_env_fragment() {
  local names=(
    RJ_SIMEMBEDDED_PHOTON12_LO
    RJ_SIMEMBEDDED_PHOTON12_HI
    RJ_SIMEMBEDDED_PHOTON20_LO
    RJ_SIMEMBEDDED_PHOTON20_HI
    RJ_SIMEMBEDDED_SIGMA_PHOTON12_PB
    RJ_SIMEMBEDDED_SIGMA_PHOTON20_PB
    RJ_SIMEMBEDDEDINCLUSIVE_JET12_LO
    RJ_SIMEMBEDDEDINCLUSIVE_JET12_HI
    RJ_SIMEMBEDDEDINCLUSIVE_JET20_LO
    RJ_SIMEMBEDDEDINCLUSIVE_JET20_HI
    RJ_SIMEMBEDDEDINCLUSIVE_JET30_LO
    RJ_SIMEMBEDDEDINCLUSIVE_JET30_HI
    RJ_SIMEMBEDDEDINCLUSIVE_JET40_LO
    RJ_SIMEMBEDDEDINCLUSIVE_SIGMA_JET12_PB
    RJ_SIMEMBEDDEDINCLUSIVE_SIGMA_JET20_PB
    RJ_SIMEMBEDDEDINCLUSIVE_SIGMA_JET30_PB
    RJ_SIMEMBEDDEDINCLUSIVE_SIGMA_JET40_PB
    RJ_RECO_CLUSTER_ET_FINE_DIAG
    RJ_RECO_CLUSTER_ET_FINE_MAX
    RJ_PP_NPB_SCORE_MIN_ET
    RJ_PP_NPB_SCORE_MAX_ET
    RJ_PP_PHOTONID_EXTRACT_ONLY
    RJ_PP_PHOTONID_TRAINING_TREE
    RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES
    RJ_PP_PHOTONID_SOURCE_ROLE
    RJ_PP_PHOTONID_PPG12_FILTER
    RJ_PP_PHOTONID_REQUIRE_PRESELECTION
    RJ_PP_INCLUSIVE_CLUSTER_ET_MAX_OVERRIDE
    RJ_PHOTON_ID_ROW_MATCH
    RJ_CURRENT_IAN_RAW_PHOTON_ID_ROW_MATCH
    RJ_THE44_PYTHIA_AUTOPSY
    RJ_THE44_PYTHIA_AUTOPSY_MAX_ENTRIES
    RJ_THE44_PYTHIA_AUTOPSY_HIGH_BDT_MIN
    RJ_THE44_PYTHIA_AUTOPSY_MID_BDT_MIN
    RJ_THE44_PYTHIA_AUTOPSY_MID_BDT_MAX
    RJ_THE44_PYTHIA_AUTOPSY_MIN_PT
    RJ_THE44_PYTHIA_AUTOPSY_MAX_PT
    RJ_THE44_PYTHIA_AUTOPSY_PARTICLE_CONE
    RJ_THE44_PYTHIA_AUTOPSY_PARTICLE_MIN_PT
    RJ_THE44_PYTHIA_AUTOPSY_MAX_PARTICLES
    RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION
    RJ_AUAU_USE_TOPOCLUSTER_ISOLATION
    RJ_AUAU_TOPOCLUSTER_EXCLUDE_CANDIDATE
    RJ_CODEX_CHAT_NAME
    RJ_CODEX_THREAD_ID
  )
  local name value
  for name in "${names[@]}"; do
    value="${!name:-}"
    [[ -n "$value" ]] && printf ';%s=%s' "$name" "$value"
  done
  return 0
}

embedded_inclusive_stitch_env_args() {
  local names=(
    RJ_SIMEMBEDDED_PHOTON12_LO
    RJ_SIMEMBEDDED_PHOTON12_HI
    RJ_SIMEMBEDDED_PHOTON20_LO
    RJ_SIMEMBEDDED_PHOTON20_HI
    RJ_SIMEMBEDDED_SIGMA_PHOTON12_PB
    RJ_SIMEMBEDDED_SIGMA_PHOTON20_PB
    RJ_SIMEMBEDDEDINCLUSIVE_JET12_LO
    RJ_SIMEMBEDDEDINCLUSIVE_JET12_HI
    RJ_SIMEMBEDDEDINCLUSIVE_JET20_LO
    RJ_SIMEMBEDDEDINCLUSIVE_JET20_HI
    RJ_SIMEMBEDDEDINCLUSIVE_JET30_LO
    RJ_SIMEMBEDDEDINCLUSIVE_JET30_HI
    RJ_SIMEMBEDDEDINCLUSIVE_JET40_LO
    RJ_SIMEMBEDDEDINCLUSIVE_SIGMA_JET12_PB
    RJ_SIMEMBEDDEDINCLUSIVE_SIGMA_JET20_PB
    RJ_SIMEMBEDDEDINCLUSIVE_SIGMA_JET30_PB
    RJ_SIMEMBEDDEDINCLUSIVE_SIGMA_JET40_PB
    RJ_RECO_CLUSTER_ET_FINE_DIAG
    RJ_RECO_CLUSTER_ET_FINE_MAX
    RJ_PP_NPB_SCORE_MIN_ET
    RJ_PP_NPB_SCORE_MAX_ET
    RJ_PP_PHOTONID_EXTRACT_ONLY
    RJ_PP_PHOTONID_TRAINING_TREE
    RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES
    RJ_PP_PHOTONID_SOURCE_ROLE
    RJ_PP_PHOTONID_PPG12_FILTER
    RJ_PP_PHOTONID_REQUIRE_PRESELECTION
    RJ_PP_INCLUSIVE_CLUSTER_ET_MAX_OVERRIDE
    RJ_PHOTON_ID_ROW_MATCH
    RJ_CURRENT_IAN_RAW_PHOTON_ID_ROW_MATCH
    RJ_THE44_PYTHIA_AUTOPSY
    RJ_THE44_PYTHIA_AUTOPSY_MAX_ENTRIES
    RJ_THE44_PYTHIA_AUTOPSY_HIGH_BDT_MIN
    RJ_THE44_PYTHIA_AUTOPSY_MID_BDT_MIN
    RJ_THE44_PYTHIA_AUTOPSY_MID_BDT_MAX
    RJ_THE44_PYTHIA_AUTOPSY_MIN_PT
    RJ_THE44_PYTHIA_AUTOPSY_MAX_PT
    RJ_THE44_PYTHIA_AUTOPSY_PARTICLE_CONE
    RJ_THE44_PYTHIA_AUTOPSY_PARTICLE_MIN_PT
    RJ_THE44_PYTHIA_AUTOPSY_MAX_PARTICLES
    RJ_AUAU_BUILD_TOPOCLUSTER_ISOLATION
    RJ_AUAU_USE_TOPOCLUSTER_ISOLATION
    RJ_AUAU_TOPOCLUSTER_EXCLUDE_CANDIDATE
    RJ_CODEX_CHAT_NAME
    RJ_CODEX_THREAD_ID
  )
  local name value
  for name in "${names[@]}"; do
    value="${!name:-}"
    [[ -n "$value" ]] && printf '%s=%s\n' "$name" "$value"
  done
  return 0
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

dphi_internal_enabled() {
  case "${RJ_DISABLE_DPHI_INTERNALIZATION:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 1 ;;
  esac
  case "${RJ_INTERNALIZE_DPHI:-1}" in
    0|false|FALSE|no|NO|off|OFF) return 1 ;;
  esac
  return 0
}

dphi_scan_csv() {
  local out="" v
  for v in "$@"; do
    v="$(trim_ws "$v")"
    [[ -z "$v" ]] && continue
    out="${out:+${out},}${v}"
  done
  printf '%s\n' "$out"
}

dphi_nominal_value() {
  local v
  for v in "$@"; do
    if sim_is_close "$v" "0.875"; then
      printf '%s\n' "$v"
      return 0
    fi
  done
  (( "$#" > 0 )) && printf '%s\n' "$1"
}

dphi_submit_values() {
  if dphi_internal_enabled && (( "$#" > 0 )); then
    dphi_nominal_value "$@"
  else
    printf '%s\n' "$@"
  fi
}

dphi_tag_component() {
  local frac="$1"
  if dphi_internal_enabled; then
    echo "dphiScan"
  else
    sim_b2b_tag "$frac"
  fi
}

dphi_env_fragment() {
  if ! dphi_internal_enabled; then
    return 0
  fi
  local csv
  csv="$(dphi_scan_csv "$@")"
  [[ -n "$csv" ]] && printf ';RJ_INTERNAL_DPHI_PI_FRACTIONS=%s' "$csv"
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

ppg12_photon_yield_enabled() {
  case "${RJ_PPG12_PHOTON_YIELD:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

iso_view_internal_enabled() {
  id_fanout_enabled || return 1
  ppg12_photon_yield_enabled && return 1
  [[ "${RJ_DIRECT_DST_DOALL:-0}" == "1" ]] && return 1
  case "${RJ_DISABLE_ISO_CONE_INTERNALIZATION:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 1 ;;
  esac
  return 0
}

iso_view_fixed_value_for_tag() {
  case "${TAG:-}" in
    auau|oo|simembedded|simembeddedinclusive) echo "${RJ_INTERNAL_FIXED_ISO_GEV_AUAU:-4.0}" ;;
    *)                                       echo "${RJ_INTERNAL_FIXED_ISO_GEV_PP:-2.0}" ;;
  esac
}

iso_view_fixed_label() {
  local fixed="$1"
  fixed="$(trim_ws "$fixed")"
  if [[ "$fixed" =~ ^([0-9]+)\.0+$ ]]; then
    fixed="${BASH_REMATCH[1]}"
  else
    fixed="${fixed//./p}"
    fixed="${fixed//-/m}"
  fi
  printf '%s\n' "$fixed"
}

iso_view_env_fragment() {
  iso_view_internal_enabled || return 0
  local override="${RJ_INTERNAL_ISO_VIEWS_OVERRIDE:-}"
  if [[ -n "$(trim_ws "$override")" ]]; then
    case "${TAG:-}" in
      auau|oo|simembedded|simembeddedinclusive)
        case "${RJ_ALLOW_FIXED_RECO_ISO_VIEWS:-0}" in
          1|true|TRUE|yes|YES|on|ON) ;;
          *)
            case "$override" in
              "isoR40_isSliding:0.40:true:0.0"|\
              "isoR40_isSliding:0.40:true:0.0,isoR30_isSliding:0.30:true:0.0") ;;
              *) die "AuAu isolation override must be canonical sliding R=0.4, optionally followed by sliding R=0.3; fixed/extra/mislabeled views require explicit RJ_ALLOW_FIXED_RECO_ISO_VIEWS=1 authorization" ;;
            esac
            ;;
        esac
        ;;
    esac
    printf ';RJ_INTERNAL_ISO_VIEWS=%s' "$override"
    return 0
  fi

  # Canonical Au+Au reconstruction always uses the centrality-dependent
  # sliding isolation.  R=0.4 is the first/canonical view and R=0.3 is a
  # robustness view.  A fixed reconstructed-isolation study must use the
  # guarded explicit override above.  This does not alter the independent
  # fixed E_T^iso,truth < 4 GeV truth-label definition.
  case "${TAG:-}" in
    auau|oo|simembedded|simembeddedinclusive)
      printf ';RJ_INTERNAL_ISO_VIEWS=isoR40_isSliding:0.40:true:0.0,isoR30_isSliding:0.30:true:0.0'
      return 0
      ;;
  esac

  local fixed fixed_label
  fixed="$(iso_view_fixed_value_for_tag)"
  fixed_label="$(iso_view_fixed_label "$fixed")"
  printf ';RJ_INTERNAL_ISO_VIEWS=isoR30_fixedIso%sGeV:0.30:false:%s,isoR40_fixedIso%sGeV:0.40:false:%s,isoR30_isSliding:0.30:true:%s,isoR40_isSliding:0.40:true:%s' \
    "$fixed_label" "$fixed" "$fixed_label" "$fixed" "$fixed" "$fixed"
}

iso_view_env_value() {
  local frag
  frag="$(iso_view_env_fragment)"
  printf '%s\n' "${frag#;RJ_INTERNAL_ISO_VIEWS=}"
}

cone_submit_values() {
  if ppg12_photon_yield_enabled; then
    printf '%s\n' "0.40"
    return 0
  fi
  if iso_view_internal_enabled; then
    local v
    for v in "$@"; do
      if sim_is_close "$v" "0.40"; then
        printf '%s\n' "$v"
        return 0
      fi
    done
    (( "$#" > 0 )) && printf '%s\n' "$1"
  else
    printf '%s\n' "$@"
  fi
}

matrix_cfg_tag() {
  local pt="$1" frac="$2" vz="$3" cone="$4" iso_idx="$5" uepipe="$6"
  local tag
  if iso_view_internal_enabled; then
    tag="${iso_selection_tags[$iso_idx]}"
    (( ${uepipe_in_tag:-0} )) && tag="${tag}_${uepipe}"
    echo "$tag"
    return 0
  fi
  tag="$(jetpt_tag_component "$pt")_$(dphi_tag_component "$frac")_$(sim_vz_tag "$vz")_$(sim_cone_tag "$cone")_${iso_base_tags[$iso_idx]}"
  (( ${uepipe_in_tag:-0} )) && tag="${tag}_${uepipe}"
  tag="${tag}_${iso_selection_tags[$iso_idx]}"
  echo "$tag"
}

direct_fanout_shard_tag() {
  local pt="$1" frac="$2" vz="$3" uepipe="$4" shard_idx="$5"
  local tag
  tag="$(jetpt_tag_component "$pt")_$(dphi_tag_component "$frac")_$(sim_vz_tag "$vz")_isoConeFanout"
  (( ${uepipe_in_tag:-0} )) && tag="${tag}_${uepipe}"
  printf '%s_shard%03d\n' "$tag" "$shard_idx"
}

sim_analysis_tag_for_dataset() {
  case "$1" in
    isSimInclusive|isSimJet5) echo "isSimInclusive" ;;
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

id_fanout_max_rows() {
  local n="${RJ_ID_FANOUT_MAX_ROWS:-15}"
  [[ "$n" =~ ^[0-9]+$ && "$n" -gt 0 ]] || { err "RJ_ID_FANOUT_MAX_ROWS must be a positive integer, got '${n}'"; exit 2; }
  printf '%d\n' "$n"
}

iso_idx_is_group_leader() {
  local iso_idx="$1"
  if iso_view_internal_enabled && id_fanout_enabled; then
    local max_rows
    max_rows="$(id_fanout_max_rows)"
    (( iso_idx % max_rows == 0 ))
    return $?
  fi
  local key
  key="$(iso_group_key "$iso_idx")"
  local j
  for (( j=0; j<iso_idx; j++ )); do
    [[ "$(iso_group_key "$j")" == "$key" ]] && return 1
  done
  return 0
}

iso_group_count() {
  if iso_view_internal_enabled && id_fanout_enabled; then
    local max_rows
    max_rows="$(id_fanout_max_rows)"
    ceil_div "${#iso_tags[@]}" "$max_rows"
    return 0
  fi
  local c=0 i
  for (( i=0; i<${#iso_tags[@]}; i++ )); do
    if iso_idx_is_group_leader "$i"; then (( c+=1 )); fi
  done
  printf '%d\n' "$c"
}

id_fanout_enabled() {
  # Fig.11 S/B background production needs one PPG12 diagnostic TH2 output.
  # The generic photon-ID fanout writes many secondary cfg ROOTs with long
  # preselection/tight/non-tight names; those are not part of the Fig.11
  # contract and can exceed filesystem filename limits.
  if env_truthy "${RJ_PPG12_FIG11_SB_DIAGNOSTIC:-0}"; then
    return 1
  fi
  case "${RJ_DISABLE_ID_FANOUT:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 1 ;;
  esac
  return 0
}

auto_merge_enabled() {
  case "${RJ_AUTO_MERGE:-1}" in
    0|false|FALSE|no|NO|off|OFF) return 1 ;;
  esac
  return 0
}

dag_dryrun_enabled() {
  case "${RJ_DAG_DRYRUN:-0}" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
  esac
  return 1
}

dag_collect_enabled() {
  [[ -n "${RJ_COLLECT_DAG_FILE:-}" ]]
}

iso_cone_fanout_enabled() {
  # Production no longer uses separate iso/cone output fanout shards. ConeR and
  # fixed/sliding isolation are internal histogram views inside each cfg ROOT.
  return 1
}

direct_fanout_shard_count() {
  local n_cones="$1"
  local n_rows="$2"
  if iso_cone_fanout_enabled; then
    local total=$(( n_cones * n_rows ))
    local max_outputs
    max_outputs="$(id_fanout_max_rows)"
    ceil_div "$total" "$max_outputs"
  else
    iso_submit_count
  fi
}

direct_worker_nevents() {
  local n="${RJ_DIRECT_NEVENTS:-0}"
  [[ "$n" == "-1" ]] && n=0
  [[ "$n" =~ ^[0-9]+$ ]] || { err "RJ_DIRECT_NEVENTS must be -1 or a non-negative integer, got '${RJ_DIRECT_NEVENTS}'"; exit 2; }
  printf '%s\n' "$n"
}

checkjobs_max_output_tags() {
  local n="${RJ_CHECKJOBS_MAX_OUTPUT_TAGS:-200}"
  [[ "$n" =~ ^[0-9]+$ && "$n" -gt 0 ]] || { err "RJ_CHECKJOBS_MAX_OUTPUT_TAGS must be a positive integer, got '${n}'"; exit 2; }
  printf '%s\n' "$n"
}

iso_submit_count() {
  if id_fanout_enabled; then
    iso_group_count
  else
    printf '%d\n' "${#iso_tags[@]}"
  fi
}

emit_id_fanout_dirs_file() {
  local out_file="$1" dest_root="$2" pt="$3" frac="$4" vz="$5" cone="$6" iso_idx="$7" uepipe="$8"
  local key cfg j
  : > "$out_file"
  if iso_view_internal_enabled; then
    local max_rows start end
    max_rows="$(id_fanout_max_rows)"
    start="$iso_idx"
    end=$(( start + max_rows ))
    (( end > ${#iso_tags[@]} )) && end="${#iso_tags[@]}"
    for (( j=start; j<end; j++ )); do
      cfg="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$j" "$uepipe")"
      printf '%s|%s|%s|%s|%s\n' \
        "${dest_root}/${cfg}" "$cfg" "${iso_preselection[$j]}" "${iso_tight[$j]}" "${iso_nonTight[$j]}" >> "$out_file"
    done
    return 0
  fi

  key="$(iso_group_key "$iso_idx")"
  for (( j=0; j<${#iso_tags[@]}; j++ )); do
    [[ "$(iso_group_key "$j")" == "$key" ]] || continue
    cfg="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$j" "$uepipe")"
    printf '%s|%s|%s|%s|%s\n' \
      "${dest_root}/${cfg}" "$cfg" "${iso_preselection[$j]}" "${iso_tight[$j]}" "${iso_nonTight[$j]}" >> "$out_file"
  done
}

iso_group_leader_ordinal() {
  local iso_idx="$1"
  local ord=0 j
  for (( j=0; j<=iso_idx; j++ )); do
    iso_idx_is_group_leader "$j" && (( ord+=1 ))
  done
  printf '%d\n' "$ord"
}

direct_fanout_cell_ordinal() {
  local cone="$1" iso_idx="$2"
  shift 2
  local -a cones=( "$@" )
  local cone_pos=-1 i
  for (( i=0; i<${#cones[@]}; i++ )); do
    if [[ "${cones[$i]}" == "$cone" ]]; then
      cone_pos="$i"
      break
    fi
  done
  (( cone_pos >= 0 )) || { err "Internal fanout error: cone ${cone} not found in cone list [${cones[*]}]"; exit 2; }
  local groups ord
  groups="$(iso_group_count)"
  ord="$(iso_group_leader_ordinal "$iso_idx")"
  printf '%d\n' $(( cone_pos * groups + ord ))
}

emit_direct_fanout_shard_file() {
  local out_file="$1" dest_root="$2" pt="$3" frac="$4" vz="$5" uepipe="$6" shard_idx="$7"
  shift 7
  local -a cones=( "$@" )
  local max_outputs
  max_outputs="$(id_fanout_max_rows)"
  : > "$out_file"
  local row=0 written=0 cone j cfg row_shard
  for cone in "${cones[@]}"; do
    for (( j=0; j<${#iso_tags[@]}; j++ )); do
      (( row+=1 ))
      row_shard=$(( (row - 1) / max_outputs + 1 ))
      (( row_shard == shard_idx )) || continue
      cfg="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$j" "$uepipe")"
      printf '%s|%s|%s|%s|%s|%s|%s|%s\n' \
        "${dest_root}/${cfg}" "$cfg" "${iso_preselection[$j]}" "${iso_tight[$j]}" "${iso_nonTight[$j]}" \
        "$cone" "${iso_sliding[$j]}" "${iso_fixed[$j]}" >> "$out_file"
      (( written+=1 ))
    done
  done
  (( written > 0 )) || { err "Internal fanout error: shard ${shard_idx} wrote no rows to ${out_file}"; exit 2; }
}

fanout_dest_allowed() {
  local fan_dest="$1"
  case "$fan_dest" in
    */thesisAna/pp/*|*/thesisAna/pp25/*|*/thesisAna/auau/*|*/thesisAna/oo/*|*/thesisAna/sim/*|*/thesisAna/simembedded/*|*/thesisAna/simembeddedinclusive/*|*/thesisAna/siminclusive/*|*/thesisAna/simjet5/*|*/thesisAna/simmb/*|*/thesisAna/recoiljets/*)
      return 0
      ;;
    */thesisAnaSmoke/pp_smokeTest_*/*|*/thesisAnaSmoke/pp25_smokeTest_*/*|*/thesisAnaSmoke/auau_smokeTest_*/*|*/thesisAnaSmoke/oo_smokeTest_*/*|*/thesisAnaSmoke/sim_smokeTest_*/*|*/thesisAnaSmoke/simembedded_smokeTest_*/*|*/thesisAnaSmoke/simembeddedinclusive_smokeTest_*/*|*/thesisAnaSmoke/siminclusive_smokeTest_*/*|*/thesisAnaSmoke/simjet5_smokeTest_*/*|*/thesisAnaSmoke/simmb_smokeTest_*/*)
      return 0
      ;;
  esac
  return 1
}

declare -a RJ_DAG_COLLECTED_NODES=()

resolve_codex_thread_id() {
  printf '%s\n' "${RJ_CODEX_THREAD_ID:-${CODEX_THREAD_ID:-}}"
}

resolve_codex_chat_name() {
  if [[ -n "${RJ_CODEX_CHAT_NAME:-}" ]]; then
    printf '%s\n' "$RJ_CODEX_CHAT_NAME"
    return 0
  fi
  if [[ -n "${CODEX_CHAT_NAME:-}" ]]; then
    printf '%s\n' "$CODEX_CHAT_NAME"
    return 0
  fi
  if [[ -n "${CODEX_THREAD_NAME:-}" ]]; then
    printf '%s\n' "$CODEX_THREAD_NAME"
    return 0
  fi
  if [[ -n "${CODEX_THREAD_TITLE:-}" ]]; then
    printf '%s\n' "$CODEX_THREAD_TITLE"
    return 0
  fi

  local thread_id
  thread_id="$(resolve_codex_thread_id)"
  if [[ -n "$thread_id" && -r "${CODEX_SESSION_INDEX:-${HOME}/.codex/session_index.jsonl}" ]]; then
    awk -v id="$thread_id" '
      index($0, "\"id\":\"" id "\"") {
        line = $0
        sub(/^.*"thread_name":"?/, "", line)
        sub(/".*$/, "", line)
        if (length(line) > 0) name = line
      }
      END { if (length(name) > 0) print name }
    ' "${CODEX_SESSION_INDEX:-${HOME}/.codex/session_index.jsonl}"
    return 0
  fi

  printf '%s\n' "unknown"
}

condor_submit_queue_args_file() {
  local sub="$1"
  awk '
    /^[[:space:]]*queue[[:space:]]+arguments[[:space:]]+from[[:space:]]+/ {
      for (i = 1; i <= NF; ++i) {
        if (tolower($i) == "from" && i < NF) {
          print $(i + 1)
          exit
        }
      }
    }
  ' "$sub"
}

stage_collected_condor_submit() {
  local sub="$1"
  [[ -s "$sub" ]] || { err "DAG collection refused missing/empty submit file: $sub"; return 66; }
  [[ -n "${RJ_COLLECT_DAG_FILE:-}" ]] || { err "DAG collection missing RJ_COLLECT_DAG_FILE"; return 66; }

  local stable_dir="${RJ_COLLECT_SUB_DIR:-$(dirname "$RJ_COLLECT_DAG_FILE")/analysis_submit_files}"
  mkdir -p "$stable_dir"

  local stable_sub="${stable_dir}/$(basename "$sub")"
  cp -f "$sub" "$stable_sub"

  local args_path
  args_path="$(condor_submit_queue_args_file "$sub" || true)"
  if [[ -n "$args_path" ]]; then
    [[ -s "$args_path" ]] || {
      err "DAG collection refused submit with missing/empty queue args file: sub=$sub args=$args_path"
      return 67
    }
    local stable_args="${stable_dir}/$(basename "$args_path")"
    cp -f "$args_path" "$stable_args"

    local tmp_sub="${stable_sub}.tmp"
    while IFS= read -r line || [[ -n "$line" ]]; do
      if [[ "$line" =~ ^[[:space:]]*queue[[:space:]]+arguments[[:space:]]+from[[:space:]]+ ]]; then
        printf 'queue arguments from %s\n' "$stable_args"
      else
        printf '%s\n' "$line"
      fi
    done < "$stable_sub" > "$tmp_sub"
    mv -f "$tmp_sub" "$stable_sub"
  fi

  [[ -s "$stable_sub" ]] || { err "DAG collection produced missing/empty stable submit file: $stable_sub"; return 68; }
  printf '%s\n' "$stable_sub"
}

validate_auto_dag_submit_files() {
  local dag="$1"
  [[ -s "$dag" ]] || { err "Automatic workflow DAG is missing/empty: $dag"; return 69; }

  local checked=0
  local missing=0
  local kind node sub args_path
  while read -r kind node sub _rest; do
    case "$kind" in
      JOB|FINAL) ;;
      *) continue ;;
    esac
    (( checked += 1 ))
    if [[ ! -s "$sub" ]]; then
      err "Automatic workflow DAG references missing/empty submit file: kind=$kind node=$node sub=$sub"
      (( missing += 1 ))
      continue
    fi
    args_path="$(condor_submit_queue_args_file "$sub" || true)"
    if [[ -n "$args_path" && ! -s "$args_path" ]]; then
      err "Automatic workflow DAG references missing/empty queue args file: kind=$kind node=$node sub=$sub args=$args_path"
      (( missing += 1 ))
    fi
  done < "$dag"

  (( checked > 0 )) || { err "Automatic workflow DAG has no JOB/FINAL submit nodes: $dag"; return 70; }
  (( missing == 0 )) || return 71
  say "Automatic workflow DAG preflight: ${checked} JOB/FINAL submit files present"
}

submit_or_collect_condor() {
  local sub="$1"
  local label="${2:-job}"
  if dag_collect_enabled; then
    local node
    node="$(sanitize_node_name "${RJ_COLLECT_NODE_PREFIX:-RJ}_${label}_${#RJ_DAG_COLLECTED_NODES[@]}")"
    local dag_sub
    dag_sub="$(stage_collected_condor_submit "$sub")" || exit $?
    printf 'JOB %s %s\n' "$node" "$dag_sub" >> "$RJ_COLLECT_DAG_FILE"
    RJ_DAG_COLLECTED_NODES+=( "$node" )
    say "Added Condor submit to orchestration DAG: node=${node} sub=${dag_sub}"
    return 0
  fi
  if dag_dryrun_enabled; then
    say "DRYRUN: validated Condor submit file without submission: ${sub}"
    return 0
  fi
  need_cmd condor_submit
  condor_submit "$sub"
}

write_auto_stage_runner() {
  local runner="$1"
  cat > "$runner" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
recoiljets_dagman_failure_re='failed with|DAG abort|aborted|Job was held|held job|POST Script failed|Node Status:[[:space:]]*STATUS_ERROR|Error:[[:space:]].*failed|return value [1-9][0-9]*|strict validation status=(CHECK|FAILED)|FINAL Node failed'
stage_key="${1:?stage key required}"
args_file="${2:?stage argument file required}"
meta_file="${3:-}"
[[ -s "$args_file" ]] || { echo "[auto-stage] missing/empty argument file: ${args_file}" >&2; exit 19; }
mapfile -t stage_cmd < "$args_file"
(( ${#stage_cmd[@]} > 0 )) || { echo "[auto-stage] no command arguments in ${args_file}" >&2; exit 19; }
emails=""
dataset="unknown"
dag_file=""
next_stage=""
if [[ -n "$meta_file" && -s "$meta_file" ]]; then
  mapfile -t meta < "$meta_file"
  emails="${meta[0]:-}"
  dataset="${meta[1]:-unknown}"
  dag_file="${meta[2]:-}"
  next_stage="${meta[3]:-}"
  codex_chat_name="${meta[4]:-unknown}"
  codex_thread_id="${meta[5]:-unknown}"
fi
poll_seconds="${RJ_ORCH_POLL_SECONDS:-120}"
[[ "$poll_seconds" =~ ^[0-9]+$ && "$poll_seconds" -gt 0 ]] || poll_seconds=120
log_file="${TMPDIR:-/tmp}/recoiljets_auto_${stage_key}_$$_submit.log"

send_stage_mail() {
  local status="$1" note="$2" clusters="${3:-}"
  [[ -n "$emails" ]] || return 0
  local subject="[RecoilJets][${stage_key}][${status}]"
  local msg
  msg="$(mktemp "${TMPDIR:-/tmp}/recoiljets_auto_stage.XXXXXX")"
  {
    echo "RECOILJETS_STAGE_EMAIL_V1"
    echo "==================== CODEX SUBMISSION ===================="
    echo "codex_chat_name=${codex_chat_name:-unknown}"
    echo "codex_thread_id=${codex_thread_id:-unknown}"
    echo "submitted_from=$(pwd)"
    echo "=========================================================="
    echo "status=${status}"
    echo "status_note=${note}"
    echo "stage=${stage_key}"
    echo "stage_type=analysis_to_merge_stage"
    echo "dataset=${dataset}"
    echo "email_policy=dataset-level stage boundary emails; nested per-cfg merge emails suppressed"
    echo "dag_file=${dag_file}"
    echo "stage_command_log=${log_file}"
    echo "clusters=${clusters}"
    echo "next_stage=${next_stage}"
    echo
    echo "message:"
    echo "  ${note}"
  } > "$msg"
  local timeout_s="${RJ_STAGE_EMAIL_TIMEOUT_SECONDS:-25}"
  [[ "$timeout_s" =~ ^[0-9]+$ && "$timeout_s" -gt 0 ]] || timeout_s=25
  if command -v mail >/dev/null 2>&1; then
    if command -v timeout >/dev/null 2>&1; then
      timeout "$timeout_s" mail -s "$subject" "$emails" < "$msg" || echo "[auto-stage] mail failed/timed out for ${emails}; status=${status}" >&2
    else
      mail -s "$subject" "$emails" < "$msg" || echo "[auto-stage] mail failed for ${emails}; status=${status}" >&2
    fi
  elif command -v mailx >/dev/null 2>&1; then
    if command -v timeout >/dev/null 2>&1; then
      timeout "$timeout_s" mailx -s "$subject" "$emails" < "$msg" || echo "[auto-stage] mailx failed/timed out for ${emails}; status=${status}" >&2
    else
      mailx -s "$subject" "$emails" < "$msg" || echo "[auto-stage] mailx failed for ${emails}; status=${status}" >&2
    fi
  else
    echo "[auto-stage] mail/mailx not found; notification skipped for ${emails}" >&2
    cat "$msg" >&2
  fi
  rm -f "$msg"
}

echo "[auto-stage] stage=${stage_key}"
printf '[auto-stage] command:'
printf ' %q' "${stage_cmd[@]}"
printf '\n'
send_stage_mail "STARTED" "Dataset=${dataset}: upstream DAG node reached ${stage_key}; launching this merge/check stage now."
"${stage_cmd[@]}" 2>&1 | tee "$log_file"
mapfile -t clusters < <(grep -Eo 'submitted to cluster [0-9]+' "$log_file" | awk '{print $4}' | sort -u)
if (( ${#clusters[@]} == 0 )); then
  echo "[auto-stage] no Condor cluster was detected in submit output for ${stage_key}" >&2
  send_stage_mail "FAILED" "Dataset=${dataset}: ${stage_key} did not report a submitted Condor cluster; inspect stage_command_log."
  exit 20
fi
echo "[auto-stage] clusters=${clusters[*]}"
while :; do
  active=0
  for c in "${clusters[@]}"; do
    if condor_q "$c" -autoformat ClusterId ProcId 2>/dev/null | grep -q .; then
      active=1
      break
    fi
  done
  (( active == 0 )) && break
  echo "[auto-stage] waiting stage=${stage_key} active_clusters=${clusters[*]} poll=${poll_seconds}s"
  sleep "$poll_seconds"
done

mapfile -t dagman_logs < <(
  grep -Eo '/[^[:space:]]+\.dag\.dagman\.out' "$log_file" | sort -u || true
)
status=0
for dlog in "${dagman_logs[@]}"; do
  dag="${dlog%.dagman.out}"
  if compgen -G "${dag}.rescue*" >/dev/null 2>&1; then
    echo "[auto-stage] rescue file found for ${stage_key}: ${dag}.rescue*" >&2
    status=30
  fi
  if [[ -s "$dlog" ]] && grep -Eiq "$recoiljets_dagman_failure_re" "$dlog"; then
    echo "[auto-stage] error/hold-like text found in ${dlog}" >&2
    status=31
  fi
done
if (( status != 0 )); then
  echo "[auto-stage] stage=${stage_key} finished with CHECK status; stopping parent DAG" >&2
  send_stage_mail "CHECK" "Dataset=${dataset}: ${stage_key} completed but strict validation found rescue/error/hold-like evidence; parent DAG will stop." "${clusters[*]}"
  exit "$status"
fi
echo "[auto-stage] stage=${stage_key} complete"
send_stage_mail "READY" "Dataset=${dataset}: ${stage_key} finished cleanly. ${next_stage}" "${clusters[*]}"
EOS
  chmod +x "$runner"
}

write_auto_final_notify() {
  local script="$1"
  cat > "$script" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
recoiljets_dagman_failure_re='failed with|DAG abort|aborted|Job was held|held job|POST Script failed|Node Status:[[:space:]]*STATUS_ERROR|Error:[[:space:]].*failed|return value [1-9][0-9]*|strict validation status=(CHECK|FAILED)|FINAL Node failed'
meta_file="${1:?metadata file required}"
[[ -s "$meta_file" ]] || { echo "[auto-notify] missing/empty metadata file: ${meta_file}" >&2; exit 19; }
mapfile -t meta < "$meta_file"
(( ${#meta[@]} >= 5 )) || { echo "[auto-notify] metadata file has too few rows: ${meta_file}" >&2; exit 19; }
emails="${meta[0]}"
stage_key="${meta[1]}"
dataset="${meta[2]}"
dag_file="${meta[3]}"
next_action="${meta[4]}"
final_output_base="${meta[5]:-}"
codex_chat_name="${meta[6]:-unknown}"
codex_thread_id="${meta[7]:-unknown}"
dagman_out="${dag_file}.dagman.out"
nodes_log="${dag_file}.nodes.log"
status="READY"
status_note="Automatic RecoilJets workflow reached the final notification node."
rescue_count=0
if compgen -G "${dag_file}.rescue*" >/dev/null 2>&1; then
  rescue_count=$(compgen -G "${dag_file}.rescue*" | wc -l | awk '{print $1}')
  status="FAILED"
  status_note="Top-level DAGMan rescue file(s) were produced; inspect the DAG logs before continuing."
elif [[ -s "$dagman_out" ]] && grep -Eiq "$recoiljets_dagman_failure_re" "$dagman_out"; then
  status="CHECK"
  status_note="Top-level DAGMan log contains error/hold-like text; inspect logs before treating outputs as ready."
fi
subject="[RecoilJets][${stage_key}][${status}]"
msg="$(mktemp "${TMPDIR:-/tmp}/recoiljets_auto_notify.XXXXXX")"
{
  echo "RECOILJETS_STAGE_EMAIL_V1"
  echo "==================== CODEX SUBMISSION ===================="
  echo "codex_chat_name=${codex_chat_name:-unknown}"
  echo "codex_thread_id=${codex_thread_id:-unknown}"
  echo "submitted_from=$(pwd)"
  echo "=========================================================="
  echo "status=${status}"
  echo "status_note=${status_note}"
  echo "stage=${stage_key}"
  echo "stage_type=analysis_to_merge_dag"
  echo "dataset=${dataset}"
  echo "email_policy=one final auto-workflow email; nested merge DAG emails suppressed with strict validation"
  echo "dag_file=${dag_file}"
  echo "dagman_out=${dagman_out}"
  echo "nodes_log=${nodes_log}"
  echo "rescue_file_count=${rescue_count}"
  if [[ -n "$final_output_base" ]]; then
    echo "final_output_base=${final_output_base}"
  fi
  echo "next_action=${next_action}"
  echo
  echo "message:"
  echo "  Automatic RecoilJets workflow finished its remote stages for dataset=${dataset}."
  if [[ -n "$final_output_base" ]]; then
    echo "  Final output base: ${final_output_base}"
  fi
  echo "  Next action: ${next_action}"
} > "$msg"
timeout_s="${RJ_STAGE_EMAIL_TIMEOUT_SECONDS:-25}"
[[ "$timeout_s" =~ ^[0-9]+$ && "$timeout_s" -gt 0 ]] || timeout_s=25
if command -v mail >/dev/null 2>&1; then
  if command -v timeout >/dev/null 2>&1; then
    timeout "$timeout_s" mail -s "$subject" "$emails" < "$msg" || echo "[auto-notify] mail command failed/timed out for ${emails}; status=${status}" >&2
  else
    mail -s "$subject" "$emails" < "$msg" || echo "[auto-notify] mail command failed for ${emails}; status=${status}" >&2
  fi
elif command -v mailx >/dev/null 2>&1; then
  if command -v timeout >/dev/null 2>&1; then
    timeout "$timeout_s" mailx -s "$subject" "$emails" < "$msg" || echo "[auto-notify] mailx command failed/timed out for ${emails}; status=${status}" >&2
  else
    mailx -s "$subject" "$emails" < "$msg" || echo "[auto-notify] mailx command failed for ${emails}; status=${status}" >&2
  fi
else
  echo "[auto-notify] mail/mailx not found; notification skipped for ${emails}" >&2
  cat "$msg" >&2
fi
rm -f "$msg"
EOS
  chmod +x "$script"
}

notify_emails_csv_from_yaml() {
  if [[ -n "${RJ_NOTIFY_EMAILS:-}" ]]; then
    printf '%s\n' "$RJ_NOTIFY_EMAILS"
    return 0
  fi
  local yaml="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  [[ -f "$yaml" ]] || return 0
  local line
  line="$(grep -E '^[[:space:]]*notify_emails:' "$yaml" | head -1 | sed 's/^[^:]*://; s/[][]//g; s/[[:space:]]//g' || true)"
  printf '%s\n' "$line"
}

add_auto_stage_node() {
  local dag="$1" node="$2" runner="$3" stage_key="$4"
  shift 4
  local sub="${dag%/*}/${node}.sub"
  local args_file="${dag%/*}/${node}.args"
  local meta_file="${dag%/*}/${node}.meta"
  local emails
  local codex_chat_name
  local codex_thread_id
  local next_stage_note
  emails="$(notify_emails_csv_from_yaml)"
  codex_chat_name="$(resolve_codex_chat_name)"
  codex_thread_id="$(resolve_codex_thread_id)"
  case "$node" in
    DATA_PERRUN)
      next_stage_note="If READY, analysis outputs merged per run and the parent DAG will start DATA_SLICERUNS next." ;;
    DATA_SLICERUNS)
      if [[ "${RJ_AUTO_FINAL_ADDCHUNKS:-1}" == "0" ]]; then
        next_stage_note="If READY, remote data sliceRuns outputs are ready; the final user-controlled addChunks/pull step is next."
      else
        next_stage_note="If READY, remote data sliceRuns outputs are ready and the parent DAG will start DATA_FINAL_ADDCHUNKS next."
      fi ;;
    DATA_FINAL_ADDCHUNKS)
      next_stage_note="If READY, final data ROOT files exist and are ready to pull or inspect." ;;
    SIM_FIRSTROUND)
      next_stage_note="If READY, remote SIM firstRound outputs are ready and the parent DAG will start SIM_SECONDROUND next." ;;
    SIM_SECONDROUND)
      next_stage_note="If READY, final SIM sample ROOT outputs exist and the parent DAG will start SIM_FINALSTITCH next." ;;
    SIM_FINALSTITCH)
      next_stage_note="If READY, final stitched SIM ROOT outputs exist and are ready to pull or inspect." ;;
    *)
      next_stage_note="If READY, the next DAG dependency is eligible to run." ;;
  esac
  printf '%s\n' "$@" > "$args_file"
  printf '%s\n%s\n%s\n%s\n%s\n%s\n' "$emails" "$DATASET" "$dag" "$next_stage_note" "$codex_chat_name" "$codex_thread_id" > "$meta_file"
  cat > "$sub" <<EOT
universe   = scheduler
executable = /bin/true
log        = $LOG_DIR/auto.${node}.\$(Cluster).\$(Process).log
notification = Never
queue
EOT
  printf 'JOB %s %s NOOP\n' "$node" "$sub" >> "$dag"
  printf 'SCRIPT POST %s %s %s %s %s\n' "$node" "$runner" "$stage_key" "$args_file" "$meta_file" >> "$dag"
}

add_auto_final_node() {
  local dag="$1" node="$2" notify_script="$3" stage_key="$4" dataset="$5" next_action="$6" final_output_base="${7:-}"
  local emails
  local codex_chat_name
  local codex_thread_id
  emails="$(notify_emails_csv_from_yaml)"
  [[ -n "$emails" ]] || return 0
  codex_chat_name="$(resolve_codex_chat_name)"
  codex_thread_id="$(resolve_codex_thread_id)"
  local sub="${dag%/*}/${node}.sub"
  local meta_file="${dag%/*}/${node}.meta"
  printf '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n' "$emails" "$stage_key" "$dataset" "$dag" "$next_action" "$final_output_base" "$codex_chat_name" "$codex_thread_id" > "$meta_file"
  cat > "$sub" <<EOT
universe   = scheduler
executable = /bin/true
log        = $LOG_DIR/auto.${node}.\$(Cluster).\$(Process).log
notification = Never
queue
EOT
  printf 'FINAL %s %s NOOP\n' "$node" "$sub" >> "$dag"
  printf 'SCRIPT POST %s %s %s\n' "$node" "$notify_script" "$meta_file" >> "$dag"
}

orchestration_self_test() {
  local failure_re='failed with|DAG abort|aborted|Job was held|held job|POST Script failed|Node Status:[[:space:]]*STATUS_ERROR|Error:[[:space:]].*failed|return value [1-9][0-9]*|strict validation status=(CHECK|FAILED)|FINAL Node failed'
  local tmp
  tmp="$(mktemp -d "${TMPDIR:-/tmp}/recoiljets_orch_selftest.XXXXXX")"
  trap 'rm -rf "$tmp"' RETURN

  cat > "${tmp}/normal.dagman.out" <<'EOF'
05/07/26 21:03:45 Daemon Log is logging: D_ALWAYS D_ERROR D_STATUS
05/07/26 21:03:45 DAGMAN_LOG_ON_NFS_IS_ERROR setting: False
05/07/26 21:03:45 DAG status: 0 (DAG_STATUS_OK)
05/07/26 21:05:25 Node MERGE job proc (123.0.0) completed successfully.
EOF

  cat > "${tmp}/failed_post.dagman.out" <<'EOF'
05/07/26 21:05:45 POST Script of node DATA_PERRUN failed with status 31
05/07/26 21:05:45 Node Status: STATUS_ERROR
EOF

  cat > "${tmp}/held_job.dagman.out" <<'EOF'
05/07/26 21:05:45 Job was held.
EOF

  cat > "${tmp}/strict_failed.dagman.out" <<'EOF'
[notify] data_perRun_auau_all: strict validation status=CHECK; failing merge stage DAG
EOF

  local failures=0
  if grep -Eiq "$failure_re" "${tmp}/normal.dagman.out"; then
    echo "RECOILJETS_ORCH_SELFTEST_FAIL normal_log_false_positive=${tmp}/normal.dagman.out"
    failures=$((failures + 1))
  else
    echo "RECOILJETS_ORCH_SELFTEST_PASS normal_dagman_noise_ignored"
  fi

  local f
  for f in failed_post held_job strict_failed; do
    if grep -Eiq "$failure_re" "${tmp}/${f}.dagman.out"; then
      echo "RECOILJETS_ORCH_SELFTEST_PASS ${f}_detected"
    else
      echo "RECOILJETS_ORCH_SELFTEST_FAIL ${f}_missed=${tmp}/${f}.dagman.out"
      failures=$((failures + 1))
    fi
  done

  if (( failures > 0 )); then
    echo "RECOILJETS_ORCH_SELFTEST_V1 status=FAIL failures=${failures}"
    return 1
  fi
  echo "RECOILJETS_ORCH_SELFTEST_V1 status=PASS"
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
    IFS='|' read -r fan_dest _fan_cfg _fan_pre _fan_tight _fan_nonTight <<< "$fan_line"
    [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
    (( rows+=1 ))
    [[ -n "${cleaned_dest[$fan_dest]:-}" ]] && continue
    cleaned_dest["$fan_dest"]=1
    fanout_dest_allowed "$fan_dest" || { err "Refusing to wipe fanout DEST_BASE='$fan_dest'"; exit 62; }
    (( cleaned+=1 ))
    (( trace && (cleaned == 1 || cleaned % every == 0) )) && say "[idFanout] cleaning unique output ${cleaned} (${label}): ${fan_dest}"
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
  (( trace )) && say "[idFanout] cleaned ${cleaned} unique output roots from ${rows} fanout rows (${label})"
  CLEAN_FANOUT_LAST_COUNT="$cleaned"
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
    tight:auauNoCentBDT|tight:auaunocentbdt) echo "auauNoCentBDT"; return 0 ;;
    tight:auauCentInputBDT|tight:auaucentinputbdt) echo "auauCentInputBDT"; return 0 ;;
    tight:auauCentInput3x3BDT|tight:auaucentinput3x3bdt) echo "auauCentInput3x3BDT"; return 0 ;;
    tight:auauCentInputBase3x3BDT|tight:auaucentinputbase3x3bdt) echo "auauCentInputBase3x3BDT"; return 0 ;;
    tight:auauCentInputMinOptBDT|tight:auaucentinputminoptbdt) echo "auauCentInputMinOptBDT"; return 0 ;;
    tight:auauCent3BDT|tight:auaucent3bdt) echo "auauCent3BDT"; return 0 ;;
    tight:auauCent7BDT|tight:auaucent7bdt) echo "auauCent7BDT"; return 0 ;;
    tight:auauPtBinCentInputBDT|tight:auauptbincentinputbdt) echo "auauPtBinCentInputBDT"; return 0 ;;
    tight:auauPtCent3BDT|tight:auauptcent3bdt) echo "auauPtCent3BDT"; return 0 ;;
    tight:auauPtCent7BDT|tight:auauptcent7bdt) echo "auauPtCent7BDT"; return 0 ;;
    tight:auauEtFineCentInputBDT|tight:auauetfinecentinputbdt) echo "auauEtFineCentInputBDT"; return 0 ;;
    tight:auauEtFineCent3BDT|tight:auauetfinecent3bdt) echo "auauEtFineCent3BDT"; return 0 ;;
    tight:auauEtFineCent7BDT|tight:auauetfinecent7bdt) echo "auauEtFineCent7BDT"; return 0 ;;
    tight:auauCentInputMLP|tight:auaucentinputmlp) echo "auauCentInputMLP"; return 0 ;;
    tight:auauNoCentBase3x3MLP|tight:auaunocentbase3x3mlp) echo "auauNoCentBase3x3MLP"; return 0 ;;
    tight:auauCentInputBase3x3MLP|tight:auaucentinputbase3x3mlp) echo "auauCentInputBase3x3MLP"; return 0 ;;
    tight:auauBDTMLPStack|tight:auaubdtmlpstack|tight:bdtmlpstack) echo "auauBDTMLPStack"; return 0 ;;
    nonTight:variantA|nonTight:VariantA|nonTight:varianta|nonTight:bdtSideband|nonTight:BDTSideband|nonTight:bdtsideband|nonTight:newPPG12|nonTight:NewPPG12|nonTight:newppg12) echo "newPPG12"; return 0 ;;
    nonTight:variantB|nonTight:VariantB|nonTight:variantb|nonTight:auauBDTSideband|nonTight:AuAuBDTSideband|nonTight:auaubdtsideband) echo "auauBDTSideband"; return 0 ;;
    nonTight:variantC|nonTight:VariantC|nonTight:variantc|nonTight:auauBDTComplement|nonTight:AuAuBDTComplement|nonTight:auaubdtcomplement) echo "auauBDTComplement"; return 0 ;;
    nonTight:auauMLPSideband|nonTight:auaumlpsideband) echo "auauMLPSideband"; return 0 ;;
    nonTight:auauMLPComplement|nonTight:auaumlpcomplement) echo "auauMLPComplement"; return 0 ;;
    nonTight:auauBDTMLPStackSideband|nonTight:auaubdtmlpstacksideband) echo "auauBDTMLPStackSideband"; return 0 ;;
    nonTight:auauBDTMLPStackComplement|nonTight:auaubdtmlpstackcomplement) echo "auauBDTMLPStackComplement"; return 0 ;;
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
    auauNoCentBDT|auaunocentbdt) echo "auauNoCentBDT" ;;
    auauCentInputBDT|auaucentinputbdt) echo "auauCentInputBDT" ;;
    auauCentInput3x3BDT|auaucentinput3x3bdt) echo "auauCentInput3x3BDT" ;;
    auauCentInputBase3x3BDT|auaucentinputbase3x3bdt) echo "auauCentInputBase3x3BDT" ;;
    auauCentInputMinOptBDT|auaucentinputminoptbdt) echo "auauCentInputMinOptBDT" ;;
    auauCent3BDT|auaucent3bdt) echo "auauCent3BDT" ;;
    auauCent7BDT|auaucent7bdt) echo "auauCent7BDT" ;;
    auauPtBinCentInputBDT|auauptbincentinputbdt) echo "auauPtBinCentInputBDT" ;;
    auauPtCent3BDT|auauptcent3bdt) echo "auauPtCent3BDT" ;;
    auauPtCent7BDT|auauptcent7bdt) echo "auauPtCent7BDT" ;;
    auauEtFineCentInputBDT|auauetfinecentinputbdt) echo "auauEtFineCentInputBDT" ;;
    auauEtFineCent3BDT|auauetfinecent3bdt) echo "auauEtFineCent3BDT" ;;
    auauEtFineCent7BDT|auauetfinecent7bdt) echo "auauEtFineCent7BDT" ;;
    auauCentInputMLP|auaucentinputmlp) echo "auauCentInputMLP" ;;
    auauNoCentBase3x3MLP|auaunocentbase3x3mlp) echo "auauNoCentBase3x3MLP" ;;
    auauCentInputBase3x3MLP|auaucentinputbase3x3mlp) echo "auauCentInputBase3x3MLP" ;;
    auauBDTMLPStack|auaubdtmlpstack|bdtmlpstack) echo "auauBDTMLPStack" ;;
    auauMLPSideband|auaumlpsideband) echo "auauMLPSideband" ;;
    auauMLPComplement|auaumlpcomplement) echo "auauMLPComplement" ;;
    auauBDTMLPStackSideband|auaubdtmlpstacksideband) echo "auauBDTMLPStackSideband" ;;
    auauBDTMLPStackComplement|auaubdtmlpstackcomplement) echo "auauBDTMLPStackComplement" ;;
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
    auauNoCentBDT) echo "${key}AuAuNoCentBDT" ;;
    auauCentInputBDT) echo "${key}AuAuCentInputBDT" ;;
    auauCentInput3x3BDT) echo "${key}AuAuCentInput3x3BDT" ;;
    auauCentInputBase3x3BDT) echo "${key}AuAuCentInputBase3x3BDT" ;;
    auauCentInputMinOptBDT) echo "${key}AuAuCentInputMinOptBDT" ;;
    auauCent3BDT) echo "${key}AuAuCent3BDT" ;;
    auauCent7BDT) echo "${key}AuAuCent7BDT" ;;
    auauPtBinCentInputBDT) echo "${key}AuAuPtBinCentInputBDT" ;;
    auauPtCent3BDT) echo "${key}AuAuPtCent3BDT" ;;
    auauPtCent7BDT) echo "${key}AuAuPtCent7BDT" ;;
    auauEtFineCentInputBDT) echo "${key}AuAuEtFineCentInputBDT" ;;
    auauEtFineCent3BDT) echo "${key}AuAuEtFineCent3BDT" ;;
    auauEtFineCent7BDT) echo "${key}AuAuEtFineCent7BDT" ;;
    auauCentInputMLP) echo "${key}AuAuCentInputMLP" ;;
    auauNoCentBase3x3MLP) echo "${key}AuAuNoCentBase3x3MLP" ;;
    auauCentInputBase3x3MLP) echo "${key}AuAuCentInputBase3x3MLP" ;;
    auauBDTMLPStack) echo "${key}AuAuBDTMLPStack" ;;
    auauMLPSideband) echo "${key}AuAuMLPSideband" ;;
    auauMLPComplement) echo "${key}AuAuMLPComplement" ;;
    auauBDTMLPStackSideband) echo "${key}AuAuBDTMLPStackSideband" ;;
    auauBDTMLPStackComplement) echo "${key}AuAuBDTMLPStackComplement" ;;
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

# Build parallel arrays. In optimized internal-iso-view production,
# iso_tags[] collapses to photon-ID selection labels only; cone/isolation
# choices are written as suffixed histograms inside each cfg ROOT file.
# In direct validation modes, the old scalar iso/cone cfg behavior is kept.
#   iso_tags[]            = full iso + selection label, or selection-only label
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

  if [[ -n "${RJ_PHOTON_ID_ROW_MATCH:-}" ]]; then
    local _match_lc="${RJ_PHOTON_ID_ROW_MATCH,,}"
    local -a _filtered_id_rows=()
    local _filter_row _filter_pre _filter_tight _filter_non _filter_tag _filter_norm
    for _filter_row in "${_id_rows[@]}"; do
      IFS='|' read -r _filter_pre _filter_tight _filter_non <<< "$_filter_row"
      _filter_pre="$(selection_mode_normalize_for_key "preselection" "$_filter_pre")"
      _filter_tight="$(selection_mode_normalize_for_key "tight" "$_filter_tight")"
      _filter_non="$(selection_mode_normalize_for_key "nonTight" "$_filter_non")"
      _filter_tag="$(selection_mode_tag "preselection" "$_filter_pre")_$(selection_mode_tag "tight" "$_filter_tight")_$(selection_mode_tag "nonTight" "$_filter_non")"
      _filter_norm="${_filter_pre}|${_filter_tight}|${_filter_non}|${_filter_tag}"
      [[ "${_filter_norm,,}" == *"${_match_lc}"* ]] && _filtered_id_rows+=( "$_filter_row" )
    done
    (( ${#_filtered_id_rows[@]} > 0 )) || { err "RJ_PHOTON_ID_ROW_MATCH='${RJ_PHOTON_ID_ROW_MATCH}' matched no photon_id_sets rows in ${yaml_file}"; exit 75; }
    say "RJ_PHOTON_ID_ROW_MATCH=${RJ_PHOTON_ID_ROW_MATCH} -> ${#_filtered_id_rows[@]}/${#_id_rows[@]} photon-ID row(s)"
    _id_rows=( "${_filtered_id_rows[@]}" )
  fi

  _both="$(trim_ws "$_both")"
  _slide="$(trim_ws "$_slide")"

  case "${TAG:-}" in
    auau|oo|simembedded|simembeddedinclusive)
      case "${RJ_ALLOW_FIXED_RECO_ISO_VIEWS:-0}" in
        1|true|TRUE|yes|YES|on|ON) ;;
        *)
          [[ "$_slide" == "true" ]] || {
            err "AuAu reconstructed isolation must be centrality-dependent sliding (isSlidingIso: true)";
            exit 76;
          }
          [[ "$_both" != "true" ]] || {
            err "Fixed reconstructed AuAu isolation views are disabled unless RJ_ALLOW_FIXED_RECO_ISO_VIEWS=1 is explicitly authorized";
            exit 76;
          }
          if (( ${#_fixeds[@]} != 1 )) || [[ ! "$(trim_ws "${_fixeds[0]}")" =~ ^[+-]?0+([.]0+)?$ ]]; then
            err "AuAu sliding-only configs must stamp fixedGeV: 0.0 as an inert compatibility sentinel";
            exit 76;
          fi
          ;;
      esac
      ;;
  esac

  if ppg12_photon_yield_enabled; then
    _both="false"
    _slide="true"
    _fixeds=( "2.0" )
  fi

  if iso_view_internal_enabled; then
    local _fixed_internal
    _fixed_internal="$(iso_view_fixed_value_for_tag)"
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

      iso_tags+=( "${selection_tag}" )
      iso_base_tags+=( "isoViewScan" )
      iso_selection_tags+=( "${selection_tag}" )
      case "${TAG:-}" in
        auau|oo|simembedded|simembeddedinclusive)
          # The stamped base row must match the canonical view even though the
          # internal loop later evaluates the R=0.4/R=0.3 sliding pair.
          # Single-view training and pre-loop code intentionally see this
          # sliding contract as well.
          iso_sliding+=( "true" )
          iso_fixed+=( "0.0" )
          ;;
        *)
          iso_sliding+=( "false" )
          iso_fixed+=( "${_fixed_internal}" )
          ;;
      esac
      iso_preselection+=( "${_pre_norm}" )
      iso_tight+=( "${_tight_norm}" )
      iso_nonTight+=( "${_nonTight_norm}" )
    done
    return 0
  fi

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

force_ppg12_photon_yield_yaml_contract() {
  local file="$1"
  ppg12_photon_yield_enabled || return 0

  yaml_set_scalar_in_place "$file" "coneR" "0.40"
  yaml_set_scalar_in_place "$file" "isSlidingIso" "true"
  yaml_set_scalar_in_place "$file" "isSlidingAndFixed" "false"
  yaml_set_scalar_in_place "$file" "fixedGeV" "2.0"
  yaml_set_scalar_in_place "$file" "isolation_wp" "{aGeV: 0.490, bPerGeV: 0.037, sideGapGeV: 0.8, truthIsoGeV: 4.0, towerMin: 0.0}"
  yaml_set_scalar_in_place "$file" "pp_iso_wp_r40" "{aGeV: 0.490, bPerGeV: 0.037, sideGapGeV: 0.8}"
}

pin_photon_id_scalars_in_yaml() {
  local file="$1" preselection="$2" tight="$3" nonTight="$4"
  yaml_set_scalar_in_place "$file" "preselection" "$preselection"
  yaml_set_scalar_in_place "$file" "tight" "$tight"
  yaml_set_scalar_in_place "$file" "nonTight" "$nonTight"
}

truthy_to_yaml_bool() {
  local value="${1:-0}"
  case "$value" in
    1|true|TRUE|yes|YES|on|ON) printf '%s\n' "true" ;;
    *)                         printf '%s\n' "false" ;;
  esac
}

propagate_pp_photonid_controls_to_yaml() {
  local file="$1"

  # RecoilJets reads the pp photon-ID extraction controls from env vars, but
  # the Fun4All steering macro also installs ProcessEnvSetter modules from the
  # generated YAML.  If these keys are absent from the per-worker YAML, those
  # setters reset the controls to their default false values inside Condor.
  # Stamp the requested controls into the YAML snapshot so local and Condor
  # extraction follow the same contract.
  if [[ -n "${RJ_PP_PHOTONID_EXTRACT_ONLY:-}" ]]; then
    yaml_set_scalar_in_place "$file" "pp_photonid_extract_only" \
      "$(truthy_to_yaml_bool "$RJ_PP_PHOTONID_EXTRACT_ONLY")"
  fi
  if [[ -n "${RJ_PP_PHOTONID_TRAINING_TREE:-}" ]]; then
    yaml_set_scalar_in_place "$file" "pp_photonid_training_tree" \
      "$(truthy_to_yaml_bool "$RJ_PP_PHOTONID_TRAINING_TREE")"
  fi
  if [[ -n "${RJ_PP_PHOTONID_PPG12_FILTER:-}" ]]; then
    yaml_set_scalar_in_place "$file" "pp_photonid_ppg12_filter" \
      "$(truthy_to_yaml_bool "$RJ_PP_PHOTONID_PPG12_FILTER")"
  fi
  if [[ -n "${RJ_PP_PHOTONID_REQUIRE_PRESELECTION:-}" ]]; then
    yaml_set_scalar_in_place "$file" "pp_photonid_require_preselection" \
      "$(truthy_to_yaml_bool "$RJ_PP_PHOTONID_REQUIRE_PRESELECTION")"
  fi
  if [[ -n "${RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES:-}" ]]; then
    yaml_set_scalar_in_place "$file" "pp_photonid_training_tree_max_entries" \
      "${RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES}"
  fi
  if [[ -n "${RJ_PP_PHOTONID_SOURCE_ROLE:-}" ]]; then
    yaml_set_scalar_in_place "$file" "pp_photonid_source_role" \
      "${RJ_PP_PHOTONID_SOURCE_ROLE}"
  fi
}

propagate_reco_cluster_et_bins_to_yaml() {
  local file="$1"
  local bins="${RJ_RECO_CLUSTER_ET_ANALYSIS_BINS:-}"
  [[ -n "$bins" ]] || return 0
  [[ "$bins" == \[*\] ]] || bins="[${bins}]"
  yaml_set_scalar_in_place "$file" "jes3_photon_pt_bins" "$bins"
  yaml_set_scalar_in_place "$file" "unfold_reco_photon_pt_bins" "$bins"
  yaml_set_scalar_in_place "$file" "unfold_truth_photon_pt_bins" "$bins"
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
  force_ppg12_photon_yield_yaml_contract "$out"
  pin_photon_id_scalars_in_yaml "$out" "$preselection" "$tight" "$nonTight"
  propagate_pp_photonid_controls_to_yaml "$out"
  propagate_reco_cluster_et_bins_to_yaml "$out"
  echo "$out"
}

# Dataset resolver (sets globals: DATASET, GOLDEN, LIST_DIR, DEST_BASE, TAG, STAGE_DIR, ROUND_DIR)
resolve_dataset() {
  case "${1:-}" in
    isPP|pp|PP)
      DATASET="isPP"
      GOLDEN="$PP_GOLDEN"
      LIST_DIR="$PP_LIST_DIR"
      LIST_PREFIX="${RJ_PP_LIST_PREFIX:-dst_ppg12_pair}"
      export RJ_PPG12_PP_DATA_PAIRED="${RJ_PPG12_PP_DATA_PAIRED:-1}"
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
      LIST_PREFIX="${RJ_AUAU_LIST_PREFIX:-dst_auau_jet_pair}"
      export RJ_AUAU_DATA_PAIRED="${RJ_AUAU_DATA_PAIRED:-1}"
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
      DEST_BASE="${RJ_SIMEMBEDINCLUSIVE_DEST_BASE:-$SIMEMBEDINCLUSIVE_DEST_BASE}"
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
    isSimInclusive|issiminclusive|siminclusive|SIMINCLUSIVE)
      DATASET="isSimInclusive"
      GOLDEN=""
      LIST_DIR=""
      LIST_PREFIX=""
      DEST_BASE="$SIMINCLUSIVE_DEST_BASE"
      TAG="siminclusive"
      MACRO="${BASE}/macros/Fun4All_recoilJets.C"
      EXE="${BASE}/RecoilJets_Condor.sh"
      IS_SIM=1
      SIM_SAMPLE_DEFAULT="run28_jet5"
      SIM_SAMPLE="$SIM_SAMPLE_DEFAULT"
      ;;
    isSimJet5|isSimjet5|simjet5|SIMJET5)
      DATASET="isSimJet5"
      GOLDEN=""
      LIST_DIR=""
      LIST_PREFIX=""
      DEST_BASE="$SIMJET5_DEST_BASE"
      TAG="siminclusive"
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
  if [[ -n "${RJ_GOLDEN_OVERRIDE:-}" ]]; then
    GOLDEN="$RJ_GOLDEN_OVERRIDE"
  fi
  if [[ -n "${RJ_LIST_DIR_OVERRIDE:-}" ]]; then
    LIST_DIR="$RJ_LIST_DIR_OVERRIDE"
  fi
  if [[ -n "${RJ_LIST_PREFIX_OVERRIDE:-}" ]]; then
    LIST_PREFIX="$RJ_LIST_PREFIX_OVERRIDE"
  fi
  mkdir -p "$STAGE_DIR" "$ROUND_DIR" "$LOG_DIR" "$OUT_DIR" "$ERR_DIR" "$SUB_DIR"
}

# Count ceil(n/d)
ceil_div() { local n="$1" d="$2"; echo $(( (n + d - 1) / d )); }

# Format run to 8 digits
run8() { printf "%08d" "$((10#$1))"; }

ppg12_pp_strict_list_coverage_enabled() {
  [[ "$DATASET" == "isPP" ]] || return 1
  case "${RJ_PPG12_PP_STRICT_LIST_COVERAGE:-1}" in
    0|false|FALSE|no|NO|off|OFF) return 1 ;;
    *) return 0 ;;
  esac
}

validate_ppg12_pp_list_coverage() {
  local source="${1:?run source required}"
  [[ -s "$source" ]] || { err "PPG12 pp data source not found or empty: $source"; return 4; }

  if [[ "$LIST_PREFIX" != "dst_ppg12_pair" ]]; then
    case "${RJ_PPG12_PP_ALLOW_LEGACY_LIST_PREFIX:-0}" in
      1|true|TRUE|yes|YES|on|ON) ;;
      *)
        err "PPG12 pp strict intake requires LIST_PREFIX=dst_ppg12_pair; got '${LIST_PREFIX}'. Set RJ_PPG12_PP_ALLOW_LEGACY_LIST_PREFIX=1 only for an explicit diagnostic."
        return 3
        ;;
    esac
  fi

  local total=0 missing=0 rn r8 src
  local missing_tmp
  missing_tmp="$(mktemp "${TMPDIR:-/tmp}/rj_ppg12_pp_missing.XXXXXX")" || return 5

  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    r8="$(run8 "$rn")"
    (( total += 1 ))
    src="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
    if [[ ! -s "$src" ]]; then
      (( missing += 1 ))
      printf '%s\t%s\n' "$r8" "$src" >> "$missing_tmp"
    fi
  done < "$source"

  if (( total == 0 )); then
    rm -f "$missing_tmp"
    err "PPG12 pp data source contains no run numbers: $source"
    return 6
  fi

  if (( missing > 0 )); then
    err "PPG12 pp data intake is incomplete: ${missing}/${total} requested runs lack submit-ready ${LIST_PREFIX}-<RUN8>.list files."
    sed -n '1,20p' "$missing_tmp" >&2 || true
    rm -f "$missing_tmp"
    return 7
  fi

  rm -f "$missing_tmp"
  say "PPG12 pp data intake coverage OK: runs=${total} list_prefix=${LIST_PREFIX} list_dir=${LIST_DIR}"
  return 0
}

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
resolve_data_input_list() {
  local src="$1" out="$2" r8="${3:-unknown}"
  [[ -s "$src" ]] || { warn "resolve input requested for missing/empty list: $src"; return 1; }

  local tmpdir macro log
  tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/rj_resolve_inputs.XXXXXX")" || return 1
  macro="${tmpdir}/rj_resolve_input_list.C"
  log="${out}.resolve.log"

  cat > "$macro" <<'ROOTMACRO'
#include <frog/FROG.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>

namespace {
bool is_protocol(const std::string& s) { return s.find("://") != std::string::npos; }
bool is_absolute(const std::string& s) { return !s.empty() && s[0] == '/'; }
bool exists_local(const std::string& s) { return gSystem && !gSystem->AccessPathName(s.c_str()); }

std::string physical_run3auau_pro001_path(const std::string& token)
{
  const std::string calo_prefix = "DST_CALOFITTING_run3auau_pro001_pcdb001_v001-";
  const std::string zdc_prefix = "DST_ZDC_RAW_run3auau_pro001_pcdb001_v001-";
  std::string stream;
  std::size_t offset = std::string::npos;
  if (token.rfind(calo_prefix, 0) == 0)
  {
    stream = "DST_CALOFITTING";
    offset = calo_prefix.size();
  }
  else if (token.rfind(zdc_prefix, 0) == 0)
  {
    stream = "DST_ZDC_RAW";
    offset = zdc_prefix.size();
  }
  else
  {
    return "";
  }

  if (token.size() < offset + 8) return "";
  const std::string run_str = token.substr(offset, 8);
  int run = 0;
  try
  {
    run = std::stoi(run_str);
  }
  catch (...)
  {
    return "";
  }
  const int lo = (run / 100) * 100;
  const int hi = lo + 100;
  char bucket[64];
  std::snprintf(bucket, sizeof(bucket), "run_%08d_%08d", lo, hi);
  const std::string path =
      "/sphenix/lustre01/sphnxpro/production/run3auau/physics/"
      "pro001_pcdb001_v001/" + stream + "/" + bucket + "/" + token;
  return exists_local(path) ? path : "";
}

std::string resolve_one(const std::string& token, int column, const std::string& run, bool& changed)
{
  if (token.empty() || token == "NONE" || token[0] == '#') return token;
  const bool required_calo = (column == 0 || token.find("DST_CALOFITTING") != std::string::npos);
  if (is_protocol(token)) return token;
  if (is_absolute(token))
  {
    if (required_calo && !exists_local(token))
    {
      std::cerr << "[resolve][FAIL] run=" << run
                << " required CALO path is absolute but inaccessible: " << token << std::endl;
      throw 10;
    }
    return token;
  }
  if (exists_local(token)) return token;
  const std::string physical = physical_run3auau_pro001_path(token);
  if (!physical.empty())
  {
    if (physical != token) changed = true;
    return physical;
  }

  FROG frog;
  const char* raw = frog.location(token.c_str());
  if (!raw || !std::string(raw).size())
  {
    if (!required_calo) return token;
    std::cerr << "[resolve][FAIL] run=" << run << " column=" << column
              << " FROG returned empty for " << token << std::endl;
    throw 11;
  }
  std::string resolved(raw);
  if (!is_protocol(resolved) && !is_absolute(resolved) && !exists_local(resolved))
  {
    if (!required_calo) return token;
    std::cerr << "[resolve][FAIL] run=" << run << " column=" << column
              << " FROG unresolved: " << token << " -> " << resolved << std::endl;
    throw 12;
  }
  if (required_calo && !is_protocol(resolved) && !exists_local(resolved))
  {
    std::cerr << "[resolve][FAIL] run=" << run
              << " required CALO path resolved but is inaccessible: "
              << token << " -> " << resolved << std::endl;
    throw 13;
  }
  if (resolved != token) changed = true;
  return resolved;
}
}  // namespace

void rj_resolve_input_list(const char* in_path, const char* out_path, const char* run)
{
  std::ifstream in(in_path);
  if (!in)
  {
    std::cerr << "[resolve][FAIL] cannot open input list: " << in_path << std::endl;
    gSystem->Exit(20);
  }
  std::ofstream out(out_path);
  if (!out)
  {
    std::cerr << "[resolve][FAIL] cannot open output list: " << out_path << std::endl;
    gSystem->Exit(21);
  }

  std::string line;
  std::size_t lines = 0;
  std::size_t kept_lines = 0;
  std::size_t skipped_lines = 0;
  std::size_t changed_lines = 0;
  const char* filter_env = std::getenv("RJ_RESOLVE_FILTER_BAD_INPUT_LINES");
  const bool filter_bad = filter_env &&
                          (std::string(filter_env) == "1" ||
                           std::string(filter_env) == "true" ||
                           std::string(filter_env) == "TRUE" ||
                           std::string(filter_env) == "yes" ||
                           std::string(filter_env) == "YES");

  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#')
    {
      out << line << '\n';
      continue;
    }
    ++lines;
    std::istringstream ss(line);
    std::vector<std::string> cols;
    std::string tok;
    while (ss >> tok) cols.push_back(tok);
    if (cols.empty())
    {
      out << line << '\n';
      continue;
    }
    bool changed = false;
    try
    {
      for (std::size_t i = 0; i < cols.size(); ++i)
      {
        cols[i] = resolve_one(cols[i], static_cast<int>(i), run, changed);
      }
    }
    catch (int rc)
    {
      if (!filter_bad) gSystem->Exit(rc);
      ++skipped_lines;
      continue;
    }
    for (std::size_t i = 0; i < cols.size(); ++i)
    {
      if (i) out << '\t';
      out << cols[i];
    }
    out << '\n';
    ++kept_lines;
    if (changed) ++changed_lines;
  }
  if (kept_lines == 0)
  {
    std::cerr << "[resolve][FAIL] run=" << run
              << " kept_lines=0 skipped_lines=" << skipped_lines
              << " input=" << in_path << std::endl;
    gSystem->Exit(30);
  }

  std::cout << "[resolve][OK] run=" << run
            << " lines=" << lines
            << " kept_lines=" << kept_lines
            << " skipped_lines=" << skipped_lines
            << " changed_lines=" << changed_lines
            << " input=" << in_path
            << " output=" << out_path << std::endl;
}
ROOTMACRO

  if ! root -l -b -q "${macro}(\"${src}\",\"${out}\",\"${r8}\")" > "$log" 2>&1; then
    warn "Input resolver failed for run ${r8}; see ${log}"
    rm -rf "$tmpdir"
    return 1
  fi
  rm -rf "$tmpdir"
  [[ -s "$out" ]] || { warn "Input resolver produced empty output for run ${r8}: ${out}"; return 1; }
  if [[ "${RJ_GROUP_TRACE:-0}" == "1" ]]; then
    tail -5 "$log" >&2 || true
  fi
}

resolve_golden_run_list() {
  local source="${1:?run source required}"
  local out_file="${2:?output run list required}"
  local stats_file="${3:-}"
  [[ -s "$source" ]] || { err "Run source not found or empty: $source"; return 5; }
  mkdir -p "$(dirname "$out_file")"
  : > "$out_file"
  if [[ -n "$stats_file" ]]; then
    mkdir -p "$(dirname "$stats_file")"
    : > "$stats_file"
  fi

  local stamp probe_stage old_stage total=0 pass=0 fail=0 r8_raw r8 src resolved
  stamp="$(date +%Y%m%d_%H%M%S)"
  probe_stage="${STAGE_DIR}/resolveProbe_${stamp}"
  mkdir -p "$probe_stage"
  old_stage="${STAGE_DIR}"

  while IFS= read -r r8_raw; do
    [[ -z "$r8_raw" || "$r8_raw" =~ ^# ]] && continue
    r8="$(run8 "$r8_raw")"
    (( total += 1 ))
    src="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
    if [[ ! -s "$src" ]]; then
      (( fail += 1 ))
      [[ -z "$stats_file" ]] || printf '%s\tFAIL\tmissing_list\t%s\n' "$r8" "$src" >> "$stats_file"
      continue
    fi
    resolved="${probe_stage}/run${r8}_resolved.list"
    if RJ_RESOLVE_FILTER_BAD_INPUT_LINES="${RJ_RESOLVE_FILTER_BAD_INPUT_LINES:-1}" resolve_data_input_list "$src" "$resolved" "$r8"; then
      (( pass += 1 ))
      printf '%s\n' "$r8" >> "$out_file"
      [[ -z "$stats_file" ]] || printf '%s\tPASS\tinput_files=%s\tresolved_files=%s\t%s\n' "$r8" "$(wc -l < "$src" | tr -d ' ')" "$(grep -vcE '^[[:space:]]*(#|$)' "$resolved" 2>/dev/null || echo 0)" "$resolved" >> "$stats_file"
    else
      (( fail += 1 ))
      [[ -z "$stats_file" ]] || printf '%s\tFAIL\tresolve_failed\t%s\n' "$r8" "${resolved}.resolve.log" >> "$stats_file"
    fi
  done < "$source"
  STAGE_DIR="$old_stage"
  say "Resolver preflight summary: source=${source} total=${total} pass=${pass} fail=${fail}"
  say "Resolved run list: ${out_file}"
  [[ -z "$stats_file" ]] || say "Resolver stats: ${stats_file}"
  (( pass > 0 ))
}

make_groups() {
  local r8="$1" gs="$2"
  local src="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
  [[ -s "$src" ]] || { warn "List is missing/empty for run ${r8} → $src"; return 1; }

  local src_for_group="$src"
  if [[ "${RJ_RESOLVE_GROUP_INPUTS:-0}" == "1" || "${RJ_RESOLVE_GROUP_INPUTS:-0}" == "true" || "${RJ_RESOLVE_GROUP_INPUTS:-0}" == "TRUE" ]]; then
    local resolved_src="${STAGE_DIR}/run${r8}_resolved.list"
    RJ_RESOLVE_FILTER_BAD_INPUT_LINES="${RJ_RESOLVE_FILTER_BAD_INPUT_LINES:-1}" resolve_data_input_list "$src" "$resolved_src" "$r8" || return 1
    src_for_group="$resolved_src"
  fi

  # Clean old groups for this run
  rm -f "${STAGE_DIR}/run${r8}_grp"*.list 2>/dev/null || true
  rm -f "${STAGE_DIR}/run${r8}_LOCAL_"*.list 2>/dev/null || true

  local nfiles; nfiles=$(wc -l < "$src_for_group" | awk '{print $1}')
  local ngroups; ngroups=$(ceil_div "$nfiles" "$gs")

  if [[ "${RJ_GROUP_TRACE:-0}" == "1" ]]; then
    say "[group] run=${r8} src=$(basename "$src_for_group") files=${nfiles} groupSize=${gs} groups=${ngroups}" >&2
  fi

  local start=1
  local g=1
  while (( g <= ngroups )); do
    local end=$(( start + gs - 1 ))
    local out="${STAGE_DIR}/run${r8}_grp$(printf "%03d" "$g").list"
    # sed is 1-indexed on lines; clamp 'end' to nfiles
    if (( end > nfiles )); then end="$nfiles"; fi
    sed -n "${start},${end}p" "$src_for_group" > "$out"
    echo "$out"
    start=$(( end + 1 ))
    g=$(( g + 1 ))
  done

  if [[ "${RJ_GROUP_TRACE:-0}" == "1" ]]; then
    say "[group] run=${r8} finished grouping under ${STAGE_DIR}" >&2
  fi
}

# ------------------------ Simulation helpers --------------------------
sim_path_validation_requested() {
  env_truthy "${RJ_VALIDATE_SIM_INPUT_PATHS:-0}" && return 0
  env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}" && return 0
  env_truthy "${RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4:-0}" && return 0
  return 1
}

sim_requires_global_lane() {
  # The deterministic stitched-purity closure canary reproduces the preserved
  # executable graphs exactly.  Neither SI nor DI registers DST_GLOBAL; the
  # vertex inputs are reconstructed inside the SI graph and absent for DI.
  env_truthy "${RJ_PPG12_CLOSURE_CANARY:-0}" && return 1
  ppg12_archived_di_campaign_requested && return 1
  env_truthy "${RJ_REQUIRE_SIM_GLOBAL:-0}" && return 0
  env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}" && return 0
  env_truthy "${RJ_PPG12_TABLE_QA:-0}" && return 0
  env_truthy "${RJ_PPG12_FIG7_TRIGGER_DIAGNOSTIC:-0}" && return 0
  env_truthy "${RJ_PPG12_FIG13_PARITY_QA:-0}" && return 0
  [[ -n "${RJ_PPG12_PERIOD:-}" ]] && return 0
  return 1
}

validate_sim_clean_list_paths() {
  local list="$1"
  local allow_none_lists="$2"
  local failures=0
  local line_no=0
  local max_lines="${RJ_VALIDATE_SIM_INPUT_MAX_LINES:-32}"
  local hard_max_lines="${RJ_LOGIN_NODE_MAX_PATH_VALIDATION_LINES:-256}"
  if [[ ! "$max_lines" =~ ^[0-9]+$ ]] || (( max_lines < 1 )); then
    err "RJ_VALIDATE_SIM_INPUT_MAX_LINES must be a positive integer"
    exit 24
  fi
  if [[ ! "$hard_max_lines" =~ ^[0-9]+$ ]] \
     || (( hard_max_lines < 1 || hard_max_lines > 256 )); then
    err "RJ_LOGIN_NODE_MAX_PATH_VALIDATION_LINES must be in [1,256]"
    exit 24
  fi
  if (( max_lines > hard_max_lines )); then
    err "SIM input validation request ${max_lines} exceeds login-node hard cap ${hard_max_lines}"
    exit 24
  fi
  local line col_idx p
  local -a cols
  while IFS= read -r line || [[ -n "$line" ]]; do
    line_no=$((line_no + 1))
    IFS=$'\t' read -r -a cols <<< "$line"
    if (( ${#cols[@]} != 5 )); then
      err "SIM input validation: line ${line_no} has ${#cols[@]} columns, expected 5 in ${list}"
      failures=$((failures + 1))
      (( failures < 20 )) || break
      continue
    fi
    for col_idx in 0 1 2 3 4; do
      p="${cols[$col_idx]}"
      if [[ "$p" == "NONE" ]]; then
        if (( col_idx == 3 )) && sim_requires_global_lane; then
          err "SIM input validation: line ${line_no} column 4 is NONE, but this PPG12 SIM contract requires a real DST_GLOBAL reco-vertex stream"
          failures=$((failures + 1))
          (( failures < 20 )) || break 2
          continue
        fi
        if (( allow_none_lists )); then
          continue
        fi
        err "SIM input validation: line ${line_no} column $((col_idx + 1)) is NONE but RJ_SIM_ALLOW_NONE_LISTS is not enabled"
        failures=$((failures + 1))
      elif [[ "$p" != /* ]]; then
        err "SIM input validation: line ${line_no} column $((col_idx + 1)) is not an absolute path: ${p}"
        failures=$((failures + 1))
      elif [[ ! -e "$p" ]]; then
        err "SIM input validation: line ${line_no} column $((col_idx + 1)) path does not exist: ${p}"
        failures=$((failures + 1))
      fi
      (( failures < 20 )) || break 2
    done
    if (( line_no >= max_lines )); then
      say "    [sim_init] validation capped at ${line_no} line(s) by RJ_VALIDATE_SIM_INPUT_MAX_LINES=${max_lines}" >&2
      break
    fi
  done < "$list"

  if (( failures > 0 )); then
    err "SIM input validation failed with ${failures} displayed failure(s); refusing to submit zero-event workers."
    err "For PPG12 double samples, regenerate lists with makePPG12DoubleSimLists.sh so relative G4/TRUTH entries resolve to js_pp200_signal_dual paths."
    exit 24
  fi
}

# Historical replay controls belong only to one explicitly typed bounded
# canary: either the stitched-purity closure oracle or the replay-foundation
# DI-neutrality witness. Ordinary production and every Au+Au mode retain their
# natural RNG/pedestal contract.
validate_ppg12_closure_canary_controls() {
  local sample="${SIM_SAMPLE:-}"
  local closure_canary=0 di_neutrality_canary=0
  env_truthy "${RJ_PPG12_CLOSURE_CANARY:-0}" && closure_canary=1
  env_truthy "${RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY:-0}" && di_neutrality_canary=1

  local inherited_extra=";${RJ_SUBMIT_EXTRA_ENV:-};"
  case "$inherited_extra" in
    *';RJ_PPG12_CLOSURE_CANARY='*|*';RJ_PPG12_CLOSURE_CANARY_ID='*|\
    *';RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY='*|*';RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY_ID='*|\
    *';RJ_PPG12_PPSIM_REPLAY_SEEDS='*|*';RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE='*|\
    *';RJ_PPG12_PPSIM_FIXED_RANDOMSEED='*|*';RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE='*)
      err "PPG12 historical canary controls must be explicit submit-shell variables, not inherited through RJ_SUBMIT_EXTRA_ENV."
      return 98
      ;;
  esac

  local replay_seeds_set=0 expected_pedestal_set=0
  [[ -n "${RJ_PPG12_PPSIM_REPLAY_SEEDS+x}" ]] && replay_seeds_set=1
  [[ -n "${RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE+x}" ]] && expected_pedestal_set=1

  if (( closure_canary && di_neutrality_canary )); then
    err "PPG12 closure and replay-foundation DI-neutrality canaries are mutually exclusive."
    return 98
  fi

  if (( ! closure_canary && ! di_neutrality_canary )); then
    if (( replay_seeds_set || expected_pedestal_set )); then
      err "PPG12 historical replay controls require an explicit closure or replay-foundation DI-neutrality canary."
      return 98
    fi
    return 0
  fi

  if (( di_neutrality_canary )); then
    [[ "$sample" =~ ^run28_(photonjet(5|10|20)|jet(8|12|20|30|40))_double$ ]] || {
      err "Replay-foundation DI-neutrality canary is restricted to frozen Run-28 double-interaction pp lanes; sample=${sample:-<unset>}."
      return 98
    }
    if ! env_truthy "${RJ_REPLAY_FOUNDATION_CANARY:-0}" ||
       ! env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}" ||
       ! env_truthy "${RJ_PPG12_PHOTON_YIELD_DOUBLE:-0}" ||
       ! env_truthy "${RJ_PPG12_PERIOD_STRICT_DI:-0}" ||
       ! env_truthy "${RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4:-0}" ||
       ! env_truthy "${RJ_PPG12_PPSIM_G4_ONLY:-0}"; then
      err "Replay-foundation DI-neutrality canary requires the pp archived-DI G4 rebuild canary path."
      return 98
    fi
    if [[ -n "${RJ_PPG12_PPSIM_FIXED_RANDOMSEED+x}" ||
          -n "${RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE+x}" ]]; then
      err "Replay-foundation DI-neutrality canary forbids synthetic fixed-seed controls."
      return 98
    fi
    if (( ! replay_seeds_set || ! expected_pedestal_set )) ||
       [[ "${RJ_PPG12_PPSIM_REPLAY_SEEDS:-}" != "2991264730,4256268992,2394322166,874466025,2240380304" ]] ||
       [[ "${RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE:-}" != "534" ]]; then
      err "Replay-foundation DI-neutrality canary requires the exact historical five-seed FIFO and pedestal sequence 534."
      return 98
    fi
    local di_canary_id="${RJ_REPLAY_FOUNDATION_DI_NEUTRALITY_CANARY_ID:-}"
    if [[ -z "$di_canary_id" || ${#di_canary_id} -gt 128 || ! "$di_canary_id" =~ ^[A-Za-z0-9_.:-]+$ ]]; then
      err "Replay-foundation DI-neutrality canary requires a safe nonempty canary ID."
      return 98
    fi
    if [[ -n "${RANDOMSEED+x}" || -n "${RJ_PPG12_PEDESTAL_OVERRIDE+x}" ||
          -n "${RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL+x}" ]]; then
      err "Replay-foundation DI-neutrality canary rejects RANDOMSEED and pedestal overrides."
      return 98
    fi
    return 0
  fi

  case "$inherited_extra" in
    *';RANDOMSEED='*|*';RJ_PPG12_PEDESTAL_OVERRIDE='*|\
    *';RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL='*|*';RJ_DISABLE_JES_CDB_AUDIT='*|\
    *';RJ_REQUIRE_SIM_GLOBAL='*)
      err "PPG12 closure canary rejects RNG, pedestal, and JES controls inherited through RJ_SUBMIT_EXTRA_ENV."
      return 98
      ;;
  esac

  [[ "$sample" =~ ^run28_(photonjet(5|10|20)|jet(8|12|20|30|40))(_double)?$ ]] || {
    err "PPG12 closure canary is restricted to the frozen Run-28 pp stitched-purity lanes; sample=${sample:-<unset>}."
    return 98
  }
  if ! env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}"; then
    err "PPG12 closure canary requires RJ_PPG12_PHOTON_YIELD=1."
    return 98
  fi
  if env_truthy "${RJ_REQUIRE_SIM_GLOBAL:-0}"; then
    err "PPG12 closure canary forbids RJ_REQUIRE_SIM_GLOBAL; the preserved executable graph has no DST_GLOBAL lane."
    return 98
  fi
  if [[ -n "${RJ_PPG12_PPSIM_FIXED_RANDOMSEED+x}" || \
        -n "${RJ_PPG12_PPSIM_FIXED_PEDESTAL_SEQUENCE+x}" ]]; then
    err "PPG12 closure canary forbids the disproven synthetic fixed-seed controls; use the exact historical FIFO replay controls."
    return 98
  fi
  if (( ! replay_seeds_set || ! expected_pedestal_set )); then
    err "PPG12 closure canary requires RJ_PPG12_PPSIM_REPLAY_SEEDS and RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE."
    return 98
  fi
  [[ "${RJ_PPG12_PPSIM_REPLAY_SEEDS:-}" == "2991264730,4256268992,2394322166,874466025,2240380304" ]] || {
    err "PPG12 closure canary requires the exact five-seed historical FIFO replay sequence."
    return 98
  }
  [[ "${RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE:-}" == "534" ]] || {
    err "PPG12 closure canary requires the replayed pedestal sequence 534; got ${RJ_PPG12_PPSIM_EXPECT_PEDESTAL_SEQUENCE:-<unset>}."
    return 98
  }
  if ! env_truthy "${RJ_DISABLE_JES_CDB_AUDIT:-0}"; then
    err "PPG12 closure canary requires RJ_DISABLE_JES_CDB_AUDIT=1 so no pre-graph JES CDB lookup is introduced."
    return 98
  fi

  local canary_id="${RJ_PPG12_CLOSURE_CANARY_ID:-}"
  if [[ -z "$canary_id" || ${#canary_id} -gt 128 || ! "$canary_id" =~ ^[A-Za-z0-9_.:-]+$ ]]; then
    err "PPG12 closure canary requires a safe nonempty RJ_PPG12_CLOSURE_CANARY_ID (letters, digits, dot, underscore, colon, or hyphen; max 128 characters)."
    return 98
  fi
  if [[ -n "${RANDOMSEED+x}" || -n "${RJ_PPG12_PEDESTAL_OVERRIDE+x}" || \
        -n "${RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL+x}" ]]; then
    err "PPG12 closure canary rejects RANDOMSEED, RJ_PPG12_PEDESTAL_OVERRIDE, and RJ_PPG12_DI_ARCHIVED_EXPECT_PEDESTAL; use only the historical replay controls."
    return 98
  fi
}

# Fail closed when a Run-28 PPG12 simulation lane's interaction-mode flags,
# sample identity, and staged source paths disagree.  DI is not defined by an
# environment flag alone: it requires a *_double sample backed by the
# js_pp200_signal_dual G4Hits, truth-jet, and global streams.  Conversely, an
# SI lane must not silently consume those dual-interaction streams.
validate_ppg12_sim_source_contract() {
  local sample="${SIM_SAMPLE:-}"
  local list="${SIM_CLEAN_LIST:-}"

  [[ "$sample" =~ ^run28_(photonjet(5|10|20)|jet(5|8|12|20|30|40))(_double)?$ ]] || return 0
  validate_ppg12_closure_canary_controls || return $?

  local closure_canary=0
  env_truthy "${RJ_PPG12_CLOSURE_CANARY:-0}" && closure_canary=1
  local archived_di=0
  ppg12_archived_di_campaign_requested && archived_di=1

  local contract_requested=0
  env_truthy "${RJ_PPG12_PHOTON_YIELD:-0}" && contract_requested=1
  env_truthy "${RJ_PPG12_PHOTON_YIELD_DOUBLE:-0}" && contract_requested=1
  env_truthy "${RJ_PPG12_PERIOD_STRICT_DI:-0}" && contract_requested=1
  env_truthy "${RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4:-0}" && contract_requested=1
  env_truthy "${RJ_PPG12_PPSIM_G4_ONLY:-0}" && contract_requested=1
  (( contract_requested )) || return 0

  [[ -s "$list" ]] || {
    err "PPG12 SIM source contract: cleaned five-column list is missing or empty: ${list:-<unset>}"
    return 98
  }

  local sample_is_double=0
  [[ "$sample" == *_double ]] && sample_is_double=1

  local flag_double=0 flag_strict_di=0 flag_rebuild=0 flag_g4_only=0
  env_truthy "${RJ_PPG12_PHOTON_YIELD_DOUBLE:-0}" && flag_double=1
  env_truthy "${RJ_PPG12_PERIOD_STRICT_DI:-0}" && flag_strict_di=1
  env_truthy "${RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4:-0}" && flag_rebuild=1
  env_truthy "${RJ_PPG12_PPSIM_G4_ONLY:-0}" && flag_g4_only=1

  if (( closure_canary )); then
    if (( sample_is_double )); then
      if (( ! flag_double || ! flag_strict_di || ! flag_rebuild || ! flag_g4_only )); then
        err "PPG12 closure DI graph requires double=1, strict_DI=1, rebuild=1, and g4_only=1 for sample=${sample}."
        return 98
      fi
      if ! awk -F '\t' '
        BEGIN { bad=0 }
        NF != 5 ||
        $1 != "NONE" ||
        $2 !~ /\/js_pp200_signal_dual\/g4hits\// ||
        $3 !~ /\/js_pp200_signal_dual\/nopileup\/jets\// ||
        $4 != "NONE" ||
        $5 != "NONE" {
          if (bad < 5) {
            printf "PPG12 closure DI graph mismatch row %d: calo=%s g4=%s truthjet=%s global=%s mbd=%s\n", NR, $1, $2, $3, $4, $5 > "/dev/stderr"
          }
          bad++
        }
        END { exit bad == 0 ? 0 : 1 }
      ' "$list"; then
        err "PPG12 closure DI graph must be exactly NONE,g4,truthjet,NONE,NONE with dual-interaction G4/truth-jet sources."
        return 98
      fi
    else
      if (( flag_double || flag_strict_di || ! flag_rebuild || ! flag_g4_only )); then
        err "PPG12 closure SI graph requires double=0, strict_DI=0, rebuild=1, and g4_only=1 for sample=${sample}."
        return 98
      fi
      if ! awk -F '\t' '
        BEGIN { bad=0 }
        NF != 5 ||
        $1 != "NONE" ||
        $2 !~ /\/js_pp200_signal\/g4hits\// ||
        $3 !~ /\/js_pp200_signal\/nopileup\/jets\// ||
        $4 != "NONE" ||
        $5 != "NONE" {
          if (bad < 5) {
            printf "PPG12 closure SI graph mismatch row %d: calo=%s g4=%s truthjet=%s global=%s mbd=%s\n", NR, $1, $2, $3, $4, $5 > "/dev/stderr"
          }
          bad++
        }
        END { exit bad == 0 ? 0 : 1 }
      ' "$list"; then
        err "PPG12 closure SI graph must be exactly NONE,g4,truthjet,NONE,NONE with single-interaction G4/truth-jet sources."
        return 98
      fi
    fi
    say "    [sim_init] exact PPG12 closure graph passed: id=${RJ_PPG12_CLOSURE_CANARY_ID} sample=${sample} mode=$([[ $sample_is_double -eq 1 ]] && printf DI || printf SI) rng=historical_fifo_replay_v2 pedestal=534" >&2
    return 0
  fi

  if (( sample_is_double )); then
    if (( ! flag_double || ! flag_strict_di || ! flag_rebuild || ! flag_g4_only )); then
      err "PPG12 DI source contract: sample=${sample} requires RJ_PPG12_PHOTON_YIELD_DOUBLE=1, RJ_PPG12_PERIOD_STRICT_DI=1, RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1, and RJ_PPG12_PPSIM_G4_ONLY=1."
      return 98
    fi
    if (( archived_di )); then
      if ! awk -F '\t' '
        BEGIN { bad=0 }
        NF != 5 ||
        $1 != "NONE" ||
        $2 !~ /\/js_pp200_signal_dual\/g4hits\// ||
        $3 !~ /\/js_pp200_signal_dual\/nopileup\/jets\// ||
        $4 != "NONE" ||
        $5 != "NONE" {
          if (bad < 5) {
            printf "PPG12 archived DI graph mismatch row %d: CALO=%s G4Hits=%s DST_JETS=%s DST_GLOBAL=%s MBD=%s\n", NR, $1, $2, $3, $4, $5 > "/dev/stderr"
          }
          bad++
        }
        END { exit bad == 0 ? 0 : 1 }
      ' "$list"; then
        err "PPG12 archived DI source contract requires exactly NONE,G4Hits,DST_JETS,NONE,NONE with dual-interaction G4/truth-jet sources."
        return 98
      fi
    elif ! awk -F '\t' '
      BEGIN { bad=0 }
      NF != 5 ||
      $2 !~ /\/js_pp200_signal_dual\/g4hits\// ||
      $3 !~ /\/js_pp200_signal_dual\/nopileup\/jets\// ||
      $4 !~ /\/js_pp200_signal_dual\/nopileup\/global\// {
        if (bad < 5) {
          printf "PPG12 DI source mismatch at row %d: G4Hits=%s DST_JETS=%s DST_GLOBAL=%s\n", NR, $2, $3, $4 > "/dev/stderr"
        }
        bad++
      }
      END { exit bad == 0 ? 0 : 1 }
    ' "$list"; then
      err "PPG12 DI source contract: refusing sample=${sample}; every row must use js_pp200_signal_dual G4Hits, truth-jet, and global streams."
      return 98
    fi
  else
    if (( flag_double || flag_strict_di || flag_rebuild || flag_g4_only )); then
      err "PPG12 SI source contract: sample=${sample} is unsuffixed but one or more DI-only flags are enabled; use the matching *_double sample or disable the DI-only flags."
      return 98
    fi
    if ! awk -F '\t' '
      BEGIN { bad=0 }
      NF != 5 ||
      $2 ~ /\/js_pp200_signal_dual\// ||
      $3 ~ /\/js_pp200_signal_dual\// ||
      $4 ~ /\/js_pp200_signal_dual\// { bad=1; exit }
      END { exit bad == 0 ? 0 : 1 }
    ' "$list"; then
      err "PPG12 SI source contract: refusing sample=${sample}; an SI sample must not consume js_pp200_signal_dual streams."
      return 98
    fi
  fi

  say "    [sim_init] PPG12 interaction/source contract passed: sample=${sample} mode=$([[ $sample_is_double -eq 1 ]] && printf DI || printf SI) rows=$(wc -l < "$list" | tr -d ' ')" >&2
}

# Canonical JSON hash of
# agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml.  This
# intentionally couples broad Run-28 pp-SIM production to the reviewed local
# closure contract.  A contract edit therefore requires a deliberate
# submitter update and a new admission; an old admission cannot drift forward.
PPG12_STITCHED_PURITY_CONTRACT_SHA256="a59bdb9d50f798885f53002f4d4945c8a94193c60d21f45d501360300ff12b14"

ppg12_sha256_file() {
  local path="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$path" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$path" | awk '{print $1}'
  else
    python3 - "$path" <<'PY'
import hashlib
import sys

h = hashlib.sha256()
with open(sys.argv[1], "rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        h.update(block)
print(h.hexdigest())
PY
  fi
}

the134_full_extraction_requested() {
  env_truthy "${RJ_THE134_MULTIVIEW_TRAINING_V1:-0}" && return 0
  env_truthy "${RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1:-0}" && return 0
  env_truthy "${RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT:-0}" && return 0
  [[ -n "${RJ_THE134_EXTRACTION_EXECUTION_MANIFEST:-}" ]] && return 0
  return 1
}

# THE-134's source-complete training extraction is not a PPG12 stitched-purity
# production.  It does, however, deliberately reuse the frozen PPG12 Run-28
# source/period/weight implementation.  Admit that one operation only when the
# independent full-extraction controller exposes its exact, hash-pinned,
# unexpired execution and authorization receipts for the current row.
#
# This is not a flag-only exemption: the receipts bind all thirteen rows, the
# group-of-seven partition, immutable submitter/config/code/source identities,
# output namespaces, duplicate guard, quota guard, and Justin's exact-scope
# authorization.  Any missing or mutated binding fails before condor_submit.
validate_the134_full_extraction_admission() {
  [[ "${ACTION:-}" == "condorDoAll" ]] || {
    err "THE-134 full extraction admission permits only condorDoAll."
    return 99
  }
  if [[ "${GROUP_SIZE_EXPLICIT:-0}" -ne 1 || "${GROUP_SIZE:-0}" -ne 7 ]]; then
    err "THE-134 full extraction admission requires explicit groupSize 7."
    return 99
  fi
  if auto_merge_enabled; then
    err "THE-134 full extraction admission requires RJ_AUTO_MERGE=0."
    return 99
  fi
  for key in \
    RJ_THE134_MULTIVIEW_TRAINING_V1 \
    RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1 \
    RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT
  do
    env_truthy "${!key:-0}" || {
      err "THE-134 full extraction admission requires ${key}=1."
      return 99
    }
  done
  if [[ "${DATASET:-}" == "isSim" || "${DATASET:-}" == "isSimInclusive" ]]; then
    for key in RJ_PP_PHOTONID_EXTRACT_ONLY RJ_PP_PHOTONID_TRAINING_TREE; do
      env_truthy "${!key:-0}" || {
        err "THE-134 p+p full extraction admission requires ${key}=1."
        return 99
      }
    done
    [[ "${RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES:-}" == "0" ]] || {
      err "THE-134 p+p full extraction admission requires an untruncated training tree."
      return 99
    }
  fi

  local execution_manifest="${RJ_THE134_EXTRACTION_EXECUTION_MANIFEST:-}"
  local execution_sha="${RJ_THE134_EXTRACTION_EXECUTION_MANIFEST_SHA256:-}"
  local authorization_receipt="${RJ_THE134_EXTRACTION_AUTHORIZATION_RECEIPT:-}"
  local authorization_sha="${RJ_THE134_EXTRACTION_AUTHORIZATION_RECEIPT_SHA256:-}"
  local row_id="${RJ_THE134_EXTRACTION_ROW_ID:-}"
  local row_fingerprint="${RJ_THE134_EXTRACTION_ROW_FINGERPRINT_SHA256:-}"
  if [[ ! -f "$execution_manifest" || -L "$execution_manifest" ||
        ! "$execution_sha" =~ ^[0-9a-f]{64}$ ||
        "$(ppg12_sha256_file "$execution_manifest")" != "$execution_sha" ]]; then
    err "THE-134 full extraction requires an exact regular execution manifest."
    return 99
  fi
  if [[ ! -f "$authorization_receipt" || -L "$authorization_receipt" ||
        ! "$authorization_sha" =~ ^[0-9a-f]{64}$ ||
        "$(ppg12_sha256_file "$authorization_receipt")" != "$authorization_sha" ]]; then
    err "THE-134 full extraction requires an exact regular authorization receipt."
    return 99
  fi
  if [[ -z "$row_id" || ! "$row_id" =~ ^[a-z0-9_]+$ ||
        ! "$row_fingerprint" =~ ^[0-9a-f]{64}$ ]]; then
    err "THE-134 full extraction row identity is missing or malformed."
    return 99
  fi
  [[ -s "${SIM_CLEAN_LIST:-}" ]] || {
    err "THE-134 full extraction cleaned source list is missing."
    return 99
  }

  local observed_jobs
  observed_jobs="$(( ($(wc -l < "$SIM_CLEAN_LIST") + 6) / 7 ))"
  if ! env \
    THE134_EXECUTION_MANIFEST="$execution_manifest" \
    THE134_EXECUTION_SHA256="$execution_sha" \
    THE134_AUTHORIZATION_RECEIPT="$authorization_receipt" \
    THE134_AUTHORIZATION_SHA256="$authorization_sha" \
    THE134_ROW_ID="$row_id" \
    THE134_ROW_FINGERPRINT_SHA256="$row_fingerprint" \
    THE134_EXPECTED_JOB_COUNT="$observed_jobs" \
    THE134_SUBMITTER_PATH="${RJ_THE134_EXTRACTION_SUBMITTER_PATH:-$0}" \
    THE134_ACTION="${ACTION:-}" \
    THE134_DATASET="${DATASET:-}" \
    THE134_SAMPLE="${SIM_SAMPLE:-}" \
    THE134_OUTPUT_NAMESPACE="${RJ_DEST_BASE_OVERRIDE:-}" \
    THE134_SUBMIT_NAMESPACE="${RJ_CONDOR_SUB_DIR:-}" \
    THE134_SUBMISSION_NAMESPACE="${RJ_SUBMISSION_NAMESPACE:-}" \
    THE134_CONFIG_PATH="${RJ_CONFIG_YAML:-}" \
    python3 - <<'PY'
import hashlib
import json
import os
import re
import time
from pathlib import Path

HEX64 = re.compile(r"^[0-9a-f]{64}$")

def fail(message):
    print(
        f"ERROR: THE-134 full extraction admission rejected: {message}",
        file=os.sys.stderr,
    )
    raise SystemExit(99)

def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def read_json(path, label):
    try:
        payload = json.loads(path.read_text())
    except Exception as exc:
        fail(f"cannot read {label}: {exc}")
    if not isinstance(payload, dict):
        fail(f"{label} is not a JSON object")
    return payload

def require_sha(value, label):
    if not isinstance(value, str) or not HEX64.fullmatch(value):
        fail(f"{label} is not a SHA-256")
    return value

execution_path = Path(os.environ["THE134_EXECUTION_MANIFEST"]).resolve()
authorization_path = Path(os.environ["THE134_AUTHORIZATION_RECEIPT"]).resolve()
execution_sha = require_sha(
    os.environ["THE134_EXECUTION_SHA256"], "execution manifest hash"
)
authorization_sha = require_sha(
    os.environ["THE134_AUTHORIZATION_SHA256"], "authorization receipt hash"
)
if file_sha256(execution_path) != execution_sha:
    fail("execution manifest hash drift")
if file_sha256(authorization_path) != authorization_sha:
    fail("authorization receipt hash drift")

execution = read_json(execution_path, "execution manifest")
authorization = read_json(authorization_path, "authorization receipt")
if (
    execution.get("schema") != "THE134_FULL_EXTRACTION_EXECUTION_MANIFEST_V1"
    or execution.get("status") != "READY_TO_SUBMIT_EXACT_GROUP7"
):
    fail("execution manifest schema/status differs")
if (
    authorization.get("schema")
    != "THE134_FULL_EXTRACTION_SUBMISSION_AUTHORIZATION_V1"
    or authorization.get("status") != "PASS_EXACT_SUBMISSION_AUTHORIZED"
):
    fail("authorization schema/status differs")
campaign = execution.get("campaign")
if not isinstance(campaign, dict) or authorization.get("campaign") != campaign:
    fail("execution and authorization campaigns differ")
bindings = execution.get("bindings")
if not isinstance(bindings, dict):
    fail("execution bindings are missing")
authorization_binding = bindings.get("authorization")
if (
    not isinstance(authorization_binding, dict)
    or Path(str(authorization_binding.get("path", ""))).resolve()
    != authorization_path
    or authorization_binding.get("sha256") != authorization_sha
):
    fail("execution authorization binding differs")
approval = authorization.get("approval")
if (
    not isinstance(approval, dict)
    or approval.get("explicit") is not True
    or approval.get("scope") != "THE134_FULL_SOURCE_COMPLETE_GROUP7_EXTRACTION"
):
    fail("exact-scope approval is missing")
expires = authorization.get("expires_at_unix")
if not isinstance(expires, int) or expires <= int(time.time()):
    fail("authorization is expired")
duplicate = authorization.get("duplicate_guard")
quota = authorization.get("quota_guard")
capacity = authorization.get("capacity_guard")
if (
    not isinstance(duplicate, dict)
    or duplicate.get("status") != "PASS_NO_ACTIVE_DUPLICATE"
    or duplicate.get("active_matching_jobs") != 0
):
    fail("duplicate guard is not PASS")
if (
    not isinstance(quota, dict)
    or quota.get("status") != "PASS"
    or quota.get("authoritative") is not True
):
    fail("quota guard is not authoritative PASS")
if capacity != {"status": "PASS", "group_size": 7, "request_memory_mb": 3000}:
    fail("capacity guard differs")
counts = execution.get("counts")
if (
    not isinstance(counts, dict)
    or counts.get("row_count") != 13
    or counts.get("group_size") != 7
    or counts.get("request_memory_mb") != 3000
):
    fail("execution limits differ")
authority = execution.get("authority")
if (
    not isinstance(authority, dict)
    or authority.get("submission_authority") is not True
    or authority.get("broad_production_authority") is not False
    or authority.get("the121_authority") is not False
    or authority.get("the122_authority") is not False
    or authority.get("canonical_promotion") is not False
):
    fail("execution authority boundary differs")

row_id = os.environ["THE134_ROW_ID"]
row_fingerprint = require_sha(
    os.environ["THE134_ROW_FINGERPRINT_SHA256"], "row fingerprint"
)
rows = execution.get("rows")
if (
    not isinstance(rows, list)
    or len(rows) != counts.get("row_count")
    or any(not isinstance(candidate, dict) for candidate in rows)
    or len({candidate.get("row_id") for candidate in rows}) != len(rows)
    or sum(
        candidate.get("expected_job_count", 0)
        for candidate in rows
        if isinstance(candidate.get("expected_job_count"), int)
    )
    != counts.get("job_count")
):
    fail("execution rows are missing")
matching = [row for row in rows if isinstance(row, dict) and row.get("row_id") == row_id]
if len(matching) != 1:
    fail("execution row is missing or duplicated")
row = matching[0]
if row.get("row_fingerprint_sha256") != row_fingerprint:
    fail("row fingerprint differs")
if (
    row.get("dataset") != os.environ["THE134_DATASET"]
    or row.get("sample") != os.environ["THE134_SAMPLE"]
    or row.get("expected_job_count") != int(os.environ["THE134_EXPECTED_JOB_COUNT"])
):
    fail("row dataset/sample/job-count differs")
if (
    row.get("analysis_output_namespace") != os.environ["THE134_OUTPUT_NAMESPACE"]
    or row.get("submit_namespace") != os.environ["THE134_SUBMIT_NAMESPACE"]
    or row.get("evidence_namespace")
    != f"{campaign.get('evidence_root')}/{row_id}"
    or os.environ["THE134_SUBMISSION_NAMESPACE"] != row_id
):
    fail("row namespace differs")
expected_sidecar = (
    f"{row['analysis_output_namespace']}/training_views/"
    "$(Cluster).$(Process).root"
)
if row.get("training_sidecar_template") != expected_sidecar:
    fail("row sidecar template differs")
submitter = row.get("submitter")
submitter_path = Path(os.environ["THE134_SUBMITTER_PATH"]).resolve()
if (
    not isinstance(submitter, dict)
    or Path(str(submitter.get("path", ""))).resolve() != submitter_path
    or file_sha256(submitter_path)
    != require_sha(submitter.get("sha256"), "row submitter hash")
):
    fail("row submitter binding differs")
if os.environ["THE134_ACTION"] != "condorDoAll":
    fail("row action differs")
config_path = Path(os.environ["THE134_CONFIG_PATH"]).resolve()
if (
    not config_path.is_file()
    or file_sha256(config_path)
    != require_sha(row.get("config_sha256"), "row config hash")
):
    fail("row configuration binding differs")
sealed = row.get("materialization_environment")
if not isinstance(sealed, dict):
    fail("row materialization environment is missing")
expected_environment = {
    "RJ_DAG_DRYRUN": "1",
    "RJ_AUTO_MERGE": "0",
    "RJ_AUTO_MEMORY_RETRY": "0",
    "RJ_AUTO_MEMORY_RETRY_MAX_RELEASES": "0",
    "RJ_REQUEST_MEMORY": "3000MB",
    "RJ_CONDOR_MAX_MATERIALIZE": "20",
    "RJ_CONDOR_MAX_IDLE": "5",
    "RJ_HOLD_FAILED_WORKERS": "1",
    "RJ_VALIDATE_SIM_INPUT_PATHS": "1",
    "RJ_VALIDATE_SIM_INPUT_MAX_LINES": "32",
    "RJ_LOGIN_NODE_MAX_PATH_VALIDATION_LINES": "256",
    "RJ_LOGIN_NODE_MAX_GROUP_FILES_PER_ROW": "2048",
    "RJ_DEST_BASE_OVERRIDE": os.environ["THE134_OUTPUT_NAMESPACE"],
    "RJ_CONDOR_SUB_DIR": os.environ["THE134_SUBMIT_NAMESPACE"],
    "RJ_SUBMISSION_NAMESPACE": row_id,
    "RJ_THE134_MULTIVIEW_TRAINING_V1": "1",
    "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1": "1",
    "RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT": "1",
}
pp_environment = {
    "RJ_PP_PHOTONID_EXTRACT_ONLY": "1",
    "RJ_PP_PHOTONID_TRAINING_TREE": "1",
    "RJ_PP_PHOTONID_TRAINING_TREE_MAX_ENTRIES": "0",
    "RJ_PPG12_PHOTON_YIELD": "1",
    "RJ_PPG12_PHOTON_YIELD_DOUBLE": "0",
}
if row.get("system") == "pp":
    expected_environment.update(pp_environment)
elif row.get("system") == "auau":
    unexpected = sorted(key for key in pp_environment if key in sealed)
    if unexpected:
        fail("Au+Au row sealed environment contains p+p-only bindings")
else:
    fail("row system is unsupported")
for key, value in expected_environment.items():
    if sealed.get(key) != value:
        fail(f"row sealed environment differs for {key}")
for row_field, environment_key in (
    ("full_source_manifest_sha256", "RJ_REPLAY_SOURCE_MANIFEST_SHA256"),
    ("config_sha256", "RJ_REPLAY_CONFIG_SHA256"),
    ("code_sha256", "RJ_REPLAY_CODE_SHA256"),
):
    if sealed.get(environment_key) != row.get(row_field):
        fail(f"row identity binding differs for {environment_key}")
profile = sealed.get("RJ_PROFILE_LABEL")
if profile != f"{campaign.get('tag')}_{row_id}":
    fail("row profile/campaign binding differs")

print(
    "THE134_FULL_EXTRACTION_ADMISSION_PASS "
    f"row={row_id} jobs={row['expected_job_count']} execution={execution_sha}"
)
PY
  then
    return 99
  fi
  # The admission already proved that this cleaned source list maps to the
  # exact row count sealed in the execution manifest.  Keep that count in the
  # current shell so group materialization cannot fall back to the legacy
  # MAX_JOBS=50000 ceiling.
  THE134_ADMITTED_EXPECTED_JOB_COUNT="$observed_jobs"
  say "    [sim_init] THE-134 source-complete sidecar extraction admitted: row=${row_id} execution_sha256=${execution_sha}" >&2
}

# Fail closed before any broad Run-28 pp photon+jet or inclusive-jet
# submission.  A deliberately bounded admission canary is the sole manifest
# exemption: exactly one five-file group per lane, isolated output namespace,
# and no automatic merge.  Production must bind the exact passing admission,
# current source list, current YAML, and every frozen physics-contract hash.
validate_ppg12_stitched_purity_admission() {
  local sample="${SIM_SAMPLE:-}"
  local list="${SIM_CLEAN_LIST:-}"
  local yaml

  case "${ACTION:-}" in
    condorDoAll|condorDoAllDirect) ;;
    *) return 0 ;;
  esac
  [[ "$sample" =~ ^run28_(photonjet(5|10|20)|jet(8|12|20|30|40))(_double)?$ ]] || return 0

  # THE-119 exercises the general replay-foundation writer on exactly one
  # isolated source group.  It is not a PPG12 stitched-purity production and
  # must not inherit that campaign's historical-RNG admission packet.  Keep
  # this exemption fail-closed and bounded so it cannot open broad running.
  if env_truthy "${RJ_REPLAY_FOUNDATION_CANARY:-0}"; then
    local replay_capacity_canary=0
    env_truthy "${RJ_REPLAY_FOUNDATION_CAPACITY_CANARY:-0}" &&
      replay_capacity_canary=1
    if (( replay_capacity_canary )); then
      if [[ "${GROUP_SIZE_EXPLICIT:-0}" -ne 1 || "${GROUP_SIZE:-0}" -ne 7 ||
            "${MAX_JOBS_EXPLICIT:-0}" -ne 1 || "${MAX_JOBS:-0}" -ne 1 ]]; then
        err "Replay-foundation capacity canary requires explicit groupSize 7 and maxJobs 1."
        return 99
      fi
      local capacity_id="${RJ_REPLAY_FOUNDATION_CAPACITY_CANARY_ID:-}"
      local capacity_receipt="${RJ_REPLAY_FOUNDATION_CAPACITY_PREFLIGHT_RECEIPT:-}"
      local capacity_receipt_sha="${RJ_REPLAY_FOUNDATION_CAPACITY_PREFLIGHT_RECEIPT_SHA256:-}"
      local execution_partition_sha="${RJ_REPLAY_FOUNDATION_EXECUTION_PARTITION_SHA256:-}"
      local bundle_receipt_sha="${RJ_REPLAY_FOUNDATION_BUNDLE_RECEIPT_SHA256:-}"
      local materialization_receipt_sha="${RJ_REPLAY_FOUNDATION_MATERIALIZATION_RECEIPT_SHA256:-}"
      if [[ -z "$capacity_id" || ${#capacity_id} -gt 128 ||
            ! "$capacity_id" =~ ^[A-Za-z0-9_.:-]+$ ]]; then
        err "Replay-foundation capacity canary requires a safe nonempty capacity ID."
        return 99
      fi
      if [[ ! -s "$capacity_receipt" ||
            ! "$capacity_receipt_sha" =~ ^[0-9a-f]{64}$ ||
            "$(ppg12_sha256_file "$capacity_receipt")" != "$capacity_receipt_sha" ]]; then
        err "Replay-foundation capacity canary requires an exact SHA-pinned preflight receipt."
        return 99
      fi
      if [[ ! "$execution_partition_sha" =~ ^[0-9a-f]{64}$ ||
            ! "$bundle_receipt_sha" =~ ^[0-9a-f]{64}$ ||
            ! "$materialization_receipt_sha" =~ ^[0-9a-f]{64}$ ]]; then
        err "Replay-foundation capacity canary requires frozen execution-partition, bundle, and materialization SHA-256 identities."
        return 99
      fi
    elif [[ "${GROUP_SIZE_EXPLICIT:-0}" -ne 1 || "${GROUP_SIZE:-0}" -ne 1 ||
            "${MAX_JOBS_EXPLICIT:-0}" -ne 1 || "${MAX_JOBS:-0}" -ne 1 ]]; then
      err "Replay-foundation canary requires explicit groupSize 1 and maxJobs 1."
      return 99
    fi
    if auto_merge_enabled; then
      err "Replay-foundation canary requires RJ_AUTO_MERGE=0."
      return 99
    fi
    if [[ -z "${RJ_REPLAY_LANE:-}" ||
          ! "${RJ_REPLAY_SCHEMA_SHA256:-}" =~ ^[0-9a-fA-F]{64}$ ||
          ! "${RJ_DEST_BASE_OVERRIDE:-}" =~ ^/sphenix/.*/replay_foundation/ ]]; then
      err "Replay-foundation canary requires a lane, schema hash, and isolated replay_foundation output path."
      return 99
    fi
    if (( replay_capacity_canary )); then
      say "    [sim_init] bounded replay-foundation capacity canary admitted: id=${capacity_id} lane=${RJ_REPLAY_LANE} sample=${sample} groupSize=7 maxJobs=1 autoMerge=off partition=${execution_partition_sha}" >&2
    else
      say "    [sim_init] bounded replay-foundation canary admitted: lane=${RJ_REPLAY_LANE} sample=${sample} groupSize=1 maxJobs=1 autoMerge=off" >&2
    fi
    return 0
  fi

  if the134_full_extraction_requested; then
    validate_the134_full_extraction_admission
    return $?
  fi

  if env_truthy "${RJ_PPG12_CLOSURE_CANARY:-0}"; then
    if [[ "${GROUP_SIZE_EXPLICIT:-0}" -ne 1 || "${GROUP_SIZE:-0}" -ne 5 || \
          "${MAX_JOBS_EXPLICIT:-0}" -ne 1 || "${MAX_JOBS:-0}" -ne 1 ]]; then
      err "PPG12 closure admission canary requires explicit groupSize 5 and maxJobs 1."
      return 99
    fi
    if auto_merge_enabled; then
      err "PPG12 closure admission canary requires RJ_AUTO_MERGE=0."
      return 99
    fi
    if [[ -z "${RJ_PPG12_CLOSURE_CANARY_ID:-}" || -z "${RJ_SUBMISSION_NAMESPACE:-}" || -z "${RJ_DEST_BASE_OVERRIDE:-}" ]]; then
      err "PPG12 closure admission canary requires RJ_PPG12_CLOSURE_CANARY_ID, RJ_SUBMISSION_NAMESPACE, and an isolated RJ_DEST_BASE_OVERRIDE."
      return 99
    fi
    if [[ ! "${RJ_SUBMISSION_NAMESPACE}" =~ ^[A-Za-z0-9_.-]+$ || ${#RJ_SUBMISSION_NAMESPACE} -gt 160 ]]; then
      err "PPG12 closure admission canary requires a safe bounded RJ_SUBMISSION_NAMESPACE."
      return 99
    fi
    case "${RJ_DEST_BASE_OVERRIDE}" in
      /sphenix/*|/tmp/*) ;;
      *)
        err "PPG12 closure admission canary output must be an absolute isolated /sphenix or /tmp path."
        return 99
        ;;
    esac
    if [[ -e "${RJ_DEST_BASE_OVERRIDE}" ]]; then
      err "PPG12 closure admission canary output already exists; use a fresh isolated path: ${RJ_DEST_BASE_OVERRIDE}"
      return 99
    fi
    say "    [sim_init] bounded PPG12 stitched-purity admission canary passed: id=${RJ_PPG12_CLOSURE_CANARY_ID} sample=${sample} groupSize=5 maxJobs=1 autoMerge=off rng=historical_fifo_replay_v2 pedestal=534 jesCdbAudit=off" >&2
    return 0
  fi

  if auto_merge_enabled; then
    err "PPG12 stitched-purity production requires RJ_AUTO_MERGE=0; merging is permitted only after the complete production audit."
    return 99
  fi

  local admission="${RJ_PPG12_CLOSURE_ADMISSION_MANIFEST:-}"
  local expected_admission_sha="${RJ_PPG12_CLOSURE_ADMISSION_SHA256:-}"
  local runtime_manifest="${RJ_PPG12_CLOSURE_RUNTIME_MANIFEST:-}"
  [[ -s "$admission" ]] || {
    err "PPG12 stitched-purity production is blocked: set RJ_PPG12_CLOSURE_ADMISSION_MANIFEST to a passing admission_manifest.json."
    return 99
  }
  [[ "$expected_admission_sha" =~ ^[0-9a-fA-F]{64}$ ]] || {
    err "PPG12 stitched-purity production is blocked: RJ_PPG12_CLOSURE_ADMISSION_SHA256 is required."
    return 99
  }
  [[ -s "$runtime_manifest" ]] || {
    err "PPG12 stitched-purity production is blocked: RJ_PPG12_CLOSURE_RUNTIME_MANIFEST must name the source-locked executed-runtime manifest."
    return 99
  }
  local actual_admission_sha
  actual_admission_sha="$(ppg12_sha256_file "$admission")"
  local actual_admission_sha_lower expected_admission_sha_lower
  actual_admission_sha_lower="$(printf '%s' "$actual_admission_sha" | tr '[:upper:]' '[:lower:]')"
  expected_admission_sha_lower="$(printf '%s' "$expected_admission_sha" | tr '[:upper:]' '[:lower:]')"
  [[ "$actual_admission_sha_lower" == "$expected_admission_sha_lower" ]] || {
    err "PPG12 stitched-purity admission file hash changed: expected=${expected_admission_sha} actual=${actual_admission_sha}."
    return 99
  }

  [[ -s "$list" ]] || { err "PPG12 stitched-purity source list is missing: ${list:-<unset>}"; return 99; }
  yaml="$(sim_yaml_master_path)"
  [[ -s "$yaml" ]] || { err "PPG12 stitched-purity YAML is missing: ${yaml:-<unset>}"; return 99; }

  local list_sha yaml_sha interaction family sample_key period lane_id
  list_sha="$(ppg12_sha256_file "$list")"
  yaml_sha="$(ppg12_sha256_file "$yaml")"
  interaction="si"
  [[ "$sample" == *_double ]] && interaction="di"
  period="${RJ_PPG12_PERIOD:-}"
  case "$period" in
    0mrad|1p5mrad) ;;
    *) err "PPG12 stitched-purity production requires RJ_PPG12_PERIOD=0mrad or 1p5mrad."; return 99 ;;
  esac
  sample_key="${sample#run28_}"
  sample_key="${sample_key%_double}"
  if [[ "$sample_key" == photonjet* ]]; then
    family="photon"
    sample_key="photon${sample_key#photonjet}"
  else
    family="inclusive"
  fi
  lane_id="${family}:${sample_key}:${period}:${interaction}"

  if ! env \
    PPG12_ADMISSION_PATH="$admission" \
    PPG12_LANE_ID="$lane_id" \
    PPG12_LIST_SHA256="$list_sha" \
    PPG12_YAML_SHA256="$yaml_sha" \
    PPG12_YAML_PATH="$yaml" \
    PPG12_RUNTIME_MANIFEST_PATH="$runtime_manifest" \
    PPG12_CONTRACT_SHA256="$PPG12_STITCHED_PURITY_CONTRACT_SHA256" \
    python3 - <<'PY'
import json
import hashlib
import os
import re
import sys
from pathlib import Path

hex64 = re.compile(r"^[0-9a-f]{64}$")

PROVENANCE_FIELDS = (
    "implementation_sha256", "source_set_sha256", "config_contract_sha256",
    "model_set_sha256", "reconstruction_contract_sha256",
    "ownership_contract_sha256", "weights_contract_sha256",
    "estimator_contract_sha256",
)
LANE_FIELDS = (
    "source_list_sha256", "config_sha256", "reconstruction_sha256",
    "model_set_sha256", "ownership_sha256", "weight_sha256",
    "estimator_sha256",
    "external_scale",
)
RUNTIME_ROLE_GROUPS = {
    "config_sha256": (
        "lane_config", "ppg_apply_bdt_config", "ppg_recoeff_period_config",
    ),
    "reconstruction_sha256": (
        "recoil_macro", "recoil_impl", "libRecoilJets.so", "libCaloAna24.so",
        "libcalo_reco.so", "libclusteriso.so", "libjetbase.so",
        "PhotonClusterBuilder.h",
    ),
    "model_set_sha256": (
        "ppg_apply_model_base_E", "ppg_apply_model_base_v3E",
        "ppg_apply_npb_model",
    ),
    "ownership_sha256": ("recoil_impl", "libRecoilJets.so"),
    "weight_sha256": (
        "recoil_impl", "ppg_recoeff_cross_section_header",
        "ppg_recoeff_truth_vertex_header", "ppg_recoeff_period_config",
        "ppg_recoeff_truth_vertex_reweight",
    ),
    "estimator_sha256": (
        "ppg_apply_bdt_macro", "ppg_recoeff_yaml_cpp_header_tree_receipt",
        "ppg_recoeff_source_macro",
        "ppg_recoeff_macro", "ppg_calculate_photon_yield",
    ),
}
RUNTIME_REQUIRED_VALUES = {
    "schema_version": 1,
    "runtime_profile": "new.17",
    "isolated_build": True,
    "estimator_revision": "29f8223bd9b36dffab07961b597afa94185bbdf1",
}

def fail(message: str) -> None:
    print(f"ERROR: PPG12 stitched-purity admission rejected: {message}", file=sys.stderr)
    raise SystemExit(99)

def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def payload_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()

def read_json(path: Path, label: str) -> dict:
    try:
        with path.open() as handle:
            payload = json.load(handle)
    except Exception as exc:
        fail(f"cannot read {label} JSON {path}: {exc}")
    if not isinstance(payload, dict):
        fail(f"{label} must contain a JSON object: {path}")
    return payload

def resolve_path(owner: Path, raw_path: object, label: str) -> Path:
    if not isinstance(raw_path, str) or not raw_path:
        fail(f"{label} path is missing")
    linked = Path(raw_path).expanduser()
    if not linked.is_absolute():
        linked = owner.parent / linked
    try:
        return linked.resolve(strict=True)
    except (FileNotFoundError, OSError) as exc:
        fail(f"{label} path cannot be resolved: {linked}: {exc}")

def resolve_file_link(owner: Path, link: object, label: str):
    if not isinstance(link, dict):
        fail(f"{label} link is missing")
    expected_sha = link.get("sha256")
    if not isinstance(expected_sha, str) or not hex64.fullmatch(expected_sha):
        fail(f"{label} link has an invalid SHA-256")
    linked = resolve_path(owner, link.get("path"), label)
    actual_sha = file_sha256(linked)
    if actual_sha != expected_sha:
        fail(
            f"{label} hash drift: path={linked} expected={expected_sha} "
            f"actual={actual_sha}"
        )
    return linked, actual_sha

def resolve_link(owner: Path, link: object, label: str):
    linked, _ = resolve_file_link(owner, link, label)
    return linked, read_json(linked, label)

def runtime_contract(runtime_path: Path, config_path: Path):
    runtime = read_json(runtime_path, "executed runtime manifest")
    for field, expected in RUNTIME_REQUIRED_VALUES.items():
        if runtime.get(field) != expected:
            fail(
                f"executed runtime manifest {field} drift: "
                f"expected={expected!r} observed={runtime.get(field)!r}"
            )
    receipt_path = resolve_path(
        runtime_path, runtime.get("build_receipt"), "runtime build receipt"
    )
    receipt_sha = runtime.get("build_receipt_sha256")
    if not isinstance(receipt_sha, str) or not hex64.fullmatch(receipt_sha):
        fail("executed runtime manifest has an invalid build-receipt SHA-256")
    if file_sha256(receipt_path) != receipt_sha:
        fail("executed runtime build receipt hash drift")
    raw_files = runtime.get("files")
    if not isinstance(raw_files, list) or not raw_files:
        fail("executed runtime manifest file list is empty")
    by_role = {}
    observed_paths = set()
    for index, row in enumerate(raw_files):
        role = row.get("role") if isinstance(row, dict) else None
        if not isinstance(role, str) or not role or role in by_role:
            fail(f"executed runtime manifest role {index} is missing or duplicated")
        linked, actual_sha = resolve_file_link(
            runtime_path, row, f"executed runtime role {role}"
        )
        if str(linked) in observed_paths:
            fail(f"executed runtime aliases one physical file under role {role}")
        observed_paths.add(str(linked))
        by_role[role] = {"role": role, "path": str(linked), "sha256": actual_sha}
    if "lane_config" in by_role:
        fail("executed runtime manifest may not override the submitted lane config")
    by_role["lane_config"] = {
        "role": "lane_config", "path": str(config_path),
        "sha256": file_sha256(config_path),
    }
    source_sets = {}
    field_hashes = {}
    for field, roles in sorted(RUNTIME_ROLE_GROUPS.items()):
        missing = sorted(set(roles) - set(by_role))
        if missing:
            fail(f"executed runtime lacks roles required for {field}: {missing}")
        files = [by_role[role] for role in sorted(roles)]
        portable = [
            {"role": row["role"], "sha256": row["sha256"]} for row in files
        ]
        set_sha = payload_sha256(portable)
        source_sets[field] = {
            "roles": sorted(roles),
            "files": files,
            "portable_role_set_sha256": set_sha,
        }
        field_hashes[field] = set_sha
    return runtime, receipt_path, receipt_sha, source_sets, field_hashes

def frozen_snapshot(manifest: dict, label: str) -> dict:
    provenance = manifest.get("provenance")
    rows = manifest.get("lanes")
    if not isinstance(provenance, dict) or not isinstance(rows, list):
        fail(f"linked {label} manifest lacks provenance or lanes")
    frozen_lanes = {}
    for index, row in enumerate(rows):
        if not isinstance(row, dict) or not isinstance(row.get("lane_id"), str):
            fail(f"linked {label} manifest has invalid lane {index}")
        lane_id = row["lane_id"]
        if lane_id in frozen_lanes:
            fail(f"linked {label} manifest duplicates lane {lane_id}")
        frozen_lanes[lane_id] = {field: row.get(field) for field in LANE_FIELDS}
    return {
        "provenance": {field: provenance.get(field) for field in PROVENANCE_FIELDS},
        "lanes": dict(sorted(frozen_lanes.items())),
        "abcd_population": manifest.get("abcd_population"),
        "external_scale": manifest.get("external_scale"),
        "random_seed": manifest.get("random_seed"),
        "toy_count": manifest.get("toy_count"),
    }

path = Path(os.environ["PPG12_ADMISSION_PATH"]).expanduser().resolve()
admission = read_json(path, "admission manifest")

if admission.get("schema") != "ppg12-stitched-purity-admission/v1" or admission.get("status") != "PASS":
    fail("manifest is not a passing v1 admission")
if admission.get("contract_sha256") != os.environ["PPG12_CONTRACT_SHA256"]:
    fail("admission was emitted under a different closure contract")

# A pass-shaped admission JSON is not sufficient.  Resolve every artifact the
# real closure gate bound, verify its bytes, and require the admission's frozen
# snapshots to be reproducible from those exact manifests.
reference_path, reference_manifest = resolve_link(
    path, admission.get("reference_manifest"), "admission reference_manifest"
)
candidate_path, candidate_manifest = resolve_link(
    path, admission.get("candidate_manifest"), "admission candidate_manifest"
)
merge_path, merge_audit = resolve_link(
    path, admission.get("merge_audit"), "admission merge_audit"
)
if reference_manifest.get("schema") != "ppg12-stitched-purity-manifest/v1" or reference_manifest.get("role") != "reference":
    fail("linked reference manifest has the wrong schema or role")
if candidate_manifest.get("schema") != "ppg12-stitched-purity-manifest/v1" or candidate_manifest.get("role") not in ("candidate", "production"):
    fail("linked candidate manifest has the wrong schema or role")
if admission.get("reference_frozen") != frozen_snapshot(reference_manifest, "reference"):
    fail("reference_frozen does not match the linked reference manifest")
if admission.get("candidate_frozen") != frozen_snapshot(candidate_manifest, "candidate"):
    fail("candidate_frozen does not match the linked candidate manifest")

# Recompute every frozen provenance source set from the current files.  Broad
# admission must never trust caller-supplied hash strings that can simply be
# copied from an old admission after implementation/model/contract drift.
assembly = candidate_manifest.get("assembly")
if not isinstance(assembly, dict):
    fail("linked candidate manifest lacks canonical assembly evidence")
provenance_path, provenance_payload = resolve_link(
    candidate_path, assembly.get("provenance"), "candidate provenance"
)
source_sets = provenance_payload.get("source_sets")
manifest_provenance = candidate_manifest.get("provenance")
if not isinstance(source_sets, dict) or not isinstance(manifest_provenance, dict):
    fail("candidate provenance lacks hash-bound source sets")
if set(source_sets) != set(PROVENANCE_FIELDS):
    fail("candidate provenance source-set coverage is not exact")
for field in PROVENANCE_FIELDS:
    source_set = source_sets.get(field)
    files = source_set.get("files") if isinstance(source_set, dict) else None
    if not isinstance(files, list) or not files:
        fail(f"candidate provenance source set is empty: {field}")
    normalized = []
    for index, link in enumerate(files):
        linked, actual_sha = resolve_file_link(
            provenance_path, link, f"candidate provenance {field} file {index}"
        )
        normalized.append({"path": str(linked), "sha256": actual_sha})
    normalized.sort(key=lambda row: row["path"])
    if len({row["path"] for row in normalized}) != len(normalized):
        fail(f"candidate provenance source set duplicates a file: {field}")
    current_set_sha = payload_sha256(normalized)
    if source_set.get("set_sha256") != current_set_sha:
        fail(f"candidate provenance source-set hash drift: {field}")
    if manifest_provenance.get(field) != current_set_sha:
        fail(f"candidate manifest provenance drift: {field}")

if merge_audit.get("schema") != "ppg12-stitched-purity-merge-audit/v1":
    fail("linked merge audit has the wrong schema")
merge_candidate_path, _ = resolve_link(
    merge_path, merge_audit.get("candidate_manifest"),
    "merge-audit candidate_manifest",
)
if merge_candidate_path != candidate_path:
    fail("linked merge audit is bound to a different candidate manifest")
audits = merge_audit.get("audits")
if not isinstance(audits, list) or not audits:
    fail("linked merge audit has no audited merge families")
for index, audit in enumerate(audits):
    if not isinstance(audit, dict) or str(audit.get("status", "")).upper() != "PASS" or audit.get("failures") not in ([], None):
        fail(f"linked merge audit family {index} is not passing")

gate_report_path = path.parent / "gate_report.json"
gate_report = read_json(gate_report_path, "admission gate report")
expected_gate_sha = admission.get("gate_report_payload_sha256")
if not isinstance(expected_gate_sha, str) or not hex64.fullmatch(expected_gate_sha):
    fail("admission gate_report_payload_sha256 is invalid")
actual_gate_sha = payload_sha256(gate_report)
if actual_gate_sha != expected_gate_sha:
    fail(
        "admission gate-report payload hash drift: "
        f"expected={expected_gate_sha} actual={actual_gate_sha}"
    )
if (
    gate_report.get("schema") != "ppg12-stitched-purity-gate-report/v1"
    or gate_report.get("mode") != "admit"
    or gate_report.get("status") != "PASS"
    or gate_report.get("failure_count") != 0
    or gate_report.get("failures") != []
):
    fail("linked admission gate report is not a zero-failure admit PASS")
report_contract = gate_report.get("contract")
if not isinstance(report_contract, dict) or report_contract.get("sha256") != os.environ["PPG12_CONTRACT_SHA256"]:
    fail("linked admission gate report has a stale closure contract")
if resolve_path(gate_report_path, gate_report.get("reference_manifest"), "gate-report reference_manifest") != reference_path:
    fail("admission gate report references a different reference manifest")
if resolve_path(gate_report_path, gate_report.get("candidate_manifest"), "gate-report candidate_manifest") != candidate_path:
    fail("admission gate report references a different candidate manifest")

expected_lanes = {
    f"inclusive:jet{sample}:{period}:{interaction}"
    for sample in (8, 12, 20, 30, 40)
    for period in ("0mrad", "1p5mrad")
    for interaction in ("si", "di")
} | {
    f"photon:photon{sample}:{period}:{interaction}"
    for sample in (5, 10, 20)
    for period in ("0mrad", "1p5mrad")
    for interaction in ("si", "di")
}

reference = admission.get("reference_frozen")
if not isinstance(reference, dict) or not isinstance(reference.get("lanes"), dict):
    fail("reference_frozen is missing")
if set(reference["lanes"]) != expected_lanes:
    fail("reference admission does not contain the exact 32-lane contract")

candidate = admission.get("candidate_frozen")
if not isinstance(candidate, dict):
    fail("candidate_frozen is missing")
if candidate.get("abcd_population") != "unsuffixed":
    fail("candidate admission does not use unsuffixed inclusive A/B/C/D")
if candidate.get("external_scale") != 1.0:
    fail("candidate admission uses a forbidden global normalization")
if candidate.get("random_seed") != 42 or candidate.get("toy_count") != 20000:
    fail("candidate admission does not use seed 42 and 20000 toys")

lanes = candidate.get("lanes")
if not isinstance(lanes, dict) or set(lanes) != expected_lanes:
    missing = sorted(expected_lanes - set(lanes or {})) if isinstance(lanes, dict) else sorted(expected_lanes)
    extra = sorted(set(lanes or {}) - expected_lanes) if isinstance(lanes, dict) else []
    fail(f"admission lane set is not the exact 32-lane contract; missing={missing} extra={extra}")

for lane_name, lane in lanes.items():
    if not isinstance(lane, dict) or lane.get("external_scale") != 1.0:
        fail(f"lane {lane_name} has a forbidden external scale")
    for field in (
        "source_list_sha256", "config_sha256", "reconstruction_sha256",
        "model_set_sha256", "ownership_sha256", "weight_sha256",
        "estimator_sha256",
    ):
        value = lane.get(field)
        if not isinstance(value, str) or not hex64.fullmatch(value):
            fail(f"lane {lane_name} has invalid {field}")

lane_id = os.environ["PPG12_LANE_ID"]
lane = lanes.get(lane_id)
if lane is None:
    fail(f"current lane is absent: {lane_id}")
if lane["source_list_sha256"] != os.environ["PPG12_LIST_SHA256"]:
    fail(f"source-list hash drift for {lane_id}")

config_path = Path(os.environ["PPG12_YAML_PATH"]).expanduser().resolve()
if file_sha256(config_path) != os.environ["PPG12_YAML_SHA256"]:
    fail(f"submitted configuration bytes changed during admission for {lane_id}")
runtime_path = Path(os.environ["PPG12_RUNTIME_MANIFEST_PATH"]).expanduser().resolve()
runtime_payload, receipt_path, receipt_sha, runtime_source_sets, runtime_hashes = (
    runtime_contract(runtime_path, config_path)
)
for field, current_hash in runtime_hashes.items():
    if lane.get(field) != current_hash:
        fail(f"executed runtime contract drift for {lane_id}: {field}")

# Re-hash the current lane's canonical extraction and require it to bind the
# same source list, submitted config, and role-labelled executed runtime.  The
# caller cannot replace these hashes with arbitrary environment strings.
raw_lane_links = assembly.get("lanes")
if not isinstance(raw_lane_links, list):
    fail("candidate manifest assembly lane links are missing")
raw_lane_link = next(
    (row for row in raw_lane_links if isinstance(row, dict) and row.get("lane_id") == lane_id),
    None,
)
raw_lane_path, raw_lane = resolve_link(
    candidate_path, raw_lane_link, f"candidate raw lane {lane_id}"
)
extraction = raw_lane.get("extraction")
evidence = extraction.get("evidence") if isinstance(extraction, dict) else None
if not isinstance(evidence, dict) or set(evidence) != {"source_list", "event_set", "config"}:
    fail(f"candidate raw lane lacks extraction evidence: {lane_id}")
source_path, source_sha = resolve_file_link(
    raw_lane_path, evidence.get("source_list"),
    f"candidate lane {lane_id} evidence source_list",
)
if source_sha != os.environ["PPG12_LIST_SHA256"] or raw_lane.get("source_list_sha256") != source_sha:
    fail(f"candidate lane source-list evidence drift for {lane_id}")
_, event_sha = resolve_file_link(
    raw_lane_path, evidence.get("event_set"),
    f"candidate lane {lane_id} evidence event_set",
)
if raw_lane.get("event_set_sha256") != event_sha:
    fail(f"candidate lane event-set evidence drift for {lane_id}")
admitted_config_path, admitted_config_sha = resolve_file_link(
    raw_lane_path, evidence.get("config"),
    f"candidate lane {lane_id} evidence config",
)
if admitted_config_path != config_path or admitted_config_sha != os.environ["PPG12_YAML_SHA256"]:
    fail(f"candidate lane submitted-config evidence drift for {lane_id}")

runtime_evidence = extraction.get("runtime_evidence")
if not isinstance(runtime_evidence, dict):
    fail(f"candidate raw lane lacks runtime evidence: {lane_id}")
admitted_runtime_path, admitted_runtime_sha = resolve_file_link(
    raw_lane_path, runtime_evidence.get("manifest"),
    f"candidate lane {lane_id} runtime manifest",
)
if admitted_runtime_path != runtime_path or admitted_runtime_sha != file_sha256(runtime_path):
    fail(f"candidate lane executed-runtime manifest drift for {lane_id}")
admitted_receipt_path, admitted_receipt_sha = resolve_file_link(
    raw_lane_path, runtime_evidence.get("build_receipt"),
    f"candidate lane {lane_id} runtime build receipt",
)
if admitted_receipt_path != receipt_path or admitted_receipt_sha != receipt_sha:
    fail(f"candidate lane runtime build-receipt drift for {lane_id}")
if runtime_evidence.get("manifest_payload_sha256") != payload_sha256(runtime_payload):
    fail(f"candidate lane runtime-manifest payload hash drift for {lane_id}")
if extraction.get("runtime_source_sets") != runtime_source_sets:
    fail(f"candidate lane role-labelled runtime source sets drift for {lane_id}")
for field, current_hash in runtime_hashes.items():
    if raw_lane.get(field) != current_hash or lane.get(field) != current_hash:
        fail(f"candidate lane runtime field drift for {lane_id}: {field}")

print(f"PPG12_STITCHED_PURITY_ADMISSION_PASS lane={lane_id}")
PY
  then
    return 99
  fi
  say "    [sim_init] PPG12 stitched-purity production admission passed: lane=${lane_id} admission_sha256=${actual_admission_sha}" >&2
}

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
  local allow_none_lists=0
  if env_truthy "${RJ_PPG12_CLOSURE_CANARY:-0}" \
    || env_truthy "${RJ_SIM_ALLOW_NONE_LISTS:-0}" \
    || env_truthy "${RJ_PPG12_PPSIM_G4_ONLY:-0}"; then
    allow_none_lists=1
  fi

  [[ -s "$g4"   ]] || { err "Missing: $g4";   err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }
  [[ -s "$jets" ]] || { err "Missing: $jets"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23; }

  SIM_MASTER_LIST="${SIM_STAGE_DIR}/sim_${SIM_SAMPLE}_${SIM_MASTER5_NAME}"
  local _ninput; _ninput=$(wc -l < "$g4" | tr -d ' ')
  local _none_dir="${SIM_STAGE_DIR}/none_placeholders"
  local _none_calo="${_none_dir}/DST_CALO_CLUSTER.NONE.list"
  local _none_glob="${_none_dir}/DST_GLOBAL.NONE.list"
  local _none_mbd="${_none_dir}/DST_MBD_EPD.NONE.list"
  make_none_sim_list() {
    local out="$1"
    mkdir -p "$(dirname "$out")"
    awk -v n="${_ninput}" 'BEGIN{for(i=0;i<n;i++) print "NONE"}' > "$out"
  }
  if env_truthy "${RJ_PPG12_CLOSURE_CANARY:-0}"; then
    # Freeze the five-column executable-oracle graph independently of which
    # optional matched lists happen to be present in the sample directory.
    make_none_sim_list "$_none_glob"
    glob="$_none_glob"
    make_none_sim_list "$_none_calo"
    make_none_sim_list "$_none_mbd"
    calo="$_none_calo"
    mbd="$_none_mbd"
  fi
  if ppg12_archived_di_campaign_requested; then
    # The hash-bound ana.541 DI runtime reconstructs detector inputs from the
    # dual-interaction G4 stream. Its accepted graph has no prebuilt CALO,
    # GLOBAL, or MBD lanes.
    make_none_sim_list "$_none_calo"
    make_none_sim_list "$_none_glob"
    make_none_sim_list "$_none_mbd"
    calo="$_none_calo"
    glob="$_none_glob"
    mbd="$_none_mbd"
  fi
  if [[ ! -s "$calo" ]]; then
    if (( allow_none_lists )); then
      make_none_sim_list "$_none_calo"
      calo="$_none_calo"
    else
      err "Missing: $calo"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23
    fi
  fi
  if [[ ! -s "$glob" ]]; then
    if sim_requires_global_lane; then
      err "Missing: $glob"
      err "This PPG12 SIM contract requires DST_GLOBAL; run makePPG12DoubleSimLists.sh or makeThesisSimLists.sh with the global lane before submitting."
      exit 23
    fi
    if (( allow_none_lists )); then
      make_none_sim_list "$_none_glob"
      glob="$_none_glob"
    else
      err "Missing: $glob"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23
    fi
  fi
  if [[ ! -s "$mbd" ]]; then
    if (( allow_none_lists )); then
      make_none_sim_list "$_none_mbd"
      mbd="$_none_mbd"
    else
      err "Missing: $mbd"; err "Run makeThesisSimLists.sh for ${SIM_SAMPLE}"; exit 23
    fi
  fi
  if [[ "${ACTION:-}" != "CHECKJOBS" ]]; then
    say "    [sim_init] paste 5 matched lists (${_ninput} G4 lines; allow_NONE=${allow_none_lists}) → master list…" >&2
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
  if sim_path_validation_requested; then
    [[ "${ACTION:-}" != "CHECKJOBS" ]] && say "    [sim_init] validating SIM input paths before submission…" >&2
    validate_sim_clean_list_paths "$SIM_CLEAN_LIST" "$allow_none_lists"
  fi
  validate_ppg12_sim_source_contract
  validate_ppg12_stitched_purity_admission

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
  local max_groups="${2:-0}"
  # This function's stdout is a typed stream containing only group-list
  # paths.  Initialization diagnostics must never become Condor arguments.
  sim_init >&2

  rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp"*.list 2>/dev/null || true

  local _nclean; _nclean=$(wc -l < "$SIM_CLEAN_LIST" | tr -d ' ')
  local _nexpect=$(( (_nclean + gs - 1) / gs ))
  local _nsource="$_nclean"
  if [[ "$max_groups" =~ ^[0-9]+$ && "$max_groups" -gt 0 && "$_nexpect" -gt "$max_groups" ]]; then
    local _ncap=$(( max_groups * gs ))
    (( _ncap > _nclean )) && _ncap="$_nclean"
    _nsource="$_ncap"
    say "    [make_sim_groups] limiting split source for maxJobs=${max_groups}: ${_nclean} → ${_nsource} lines" >&2
  fi
  say "    [make_sim_groups] splitting ${_nsource} lines into chunks of ${gs} (expect ~${_nexpect} groups before cap)…" >&2

  # One bounded Python process preserves split(1)'s exact byte/order semantics
  # without spawning one mv(1) process per logical group on the login node.
  local hard_max_groups="${RJ_LOGIN_NODE_MAX_GROUP_FILES_PER_ROW:-2048}"
  if [[ ! "$hard_max_groups" =~ ^[0-9]+$ ]] \
     || (( hard_max_groups < 1 || hard_max_groups > 4096 )); then
    err "RJ_LOGIN_NODE_MAX_GROUP_FILES_PER_ROW must be in [1,4096]"
    return 25
  fi
  python3 - \
    "$SIM_CLEAN_LIST" "$SIM_STAGE_DIR" "$SIM_JOB_PREFIX" \
    "$gs" "$max_groups" "$hard_max_groups" "${RJ_GROUP_TRACE:-0}" <<'PY'
import os
import sys
from pathlib import Path

source = Path(sys.argv[1])
stage = Path(sys.argv[2])
prefix = sys.argv[3]
group_size = int(sys.argv[4])
max_groups = int(sys.argv[5])
hard_max_groups = int(sys.argv[6])
trace = sys.argv[7] == "1"

if group_size < 1 or max_groups < 0:
    raise SystemExit("invalid group-size/max-groups contract")
if max_groups > hard_max_groups:
    raise SystemExit(
        f"requested max-groups {max_groups} exceeds hard cap {hard_max_groups}"
    )
created: list[Path] = []
stream = None
try:
    with source.open("rb") as handle:
        for line_index, line in enumerate(handle):
            group_index = line_index // group_size + 1
            if max_groups and group_index > max_groups:
                break
            if group_index > hard_max_groups:
                raise RuntimeError(
                    f"manifest requires more than {hard_max_groups} group files"
                )
            if line_index % group_size == 0:
                if stream is not None:
                    stream.close()
                destination = stage / f"{prefix}_grp{group_index:03d}.list"
                descriptor = os.open(
                    destination,
                    os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                    0o644,
                )
                os.fchmod(descriptor, 0o644)
                stream = os.fdopen(descriptor, "wb")
                created.append(destination)
            stream.write(line)
    if stream is not None:
        stream.close()
        stream = None
except Exception:
    if stream is not None:
        stream.close()
    for destination in created:
        try:
            destination.unlink()
        except FileNotFoundError:
            pass
    raise

if not created:
    raise SystemExit("no SIM groups were produced")
for destination in created:
    print(destination)
print(f"    [make_sim_groups] wrote {len(created)} group files in one bounded process", file=sys.stderr)
if trace:
    print(
        f"    [make_sim_groups] first={created[0].name} last={created[-1].name}",
        file=sys.stderr,
    )
PY
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
  mapfile -t sim_vzs   < <( dataset_vz_values "$master_yaml" )
  mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
  (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
  (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
  (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
  (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
  local -a sim_view_cones=( "${sim_cones[@]}" )
  mapfile -t sim_cones < <( cone_submit_values "${sim_view_cones[@]}" )
  local -a sim_submit_pts sim_submit_fracs
  mapfile -t sim_submit_pts < <( jetpt_submit_values "${sim_pts[@]}" )
  mapfile -t sim_submit_fracs < <( dphi_submit_values "${sim_fracs[@]}" )
  build_iso_modes "$master_yaml"
  read_uepipe_modes "$master_yaml" "$TAG"

  if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
    case "$DATASET" in
      isSimEmbedded)          samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
      isSimEmbeddedInclusive) mapfile -t samples < <(simembeddedinclusive_sample_list) ;;
      isSimInclusive|isSimJet5) samples=( "run28_jet5" "run28_jet8" "run28_jet12" "run28_jet20" "run28_jet30" "run28_jet40" ) ;;
      isSimMB)                samples=( "run28_detroit" ) ;;
      *)                      samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" ) ;;
    esac
  else
    samples=( "${SIM_SAMPLE}" )
  fi

  local iso_submit_n
  iso_submit_n="$(iso_submit_count)"
  local fanout_shards_per_axis=0
  local final_cfg_n=0
  if iso_cone_fanout_enabled; then
    fanout_shards_per_axis="$(direct_fanout_shard_count "${#sim_cones[@]}" "${#iso_tags[@]}")"
    n_cfg=$(( ${#sim_submit_pts[@]} * ${#sim_submit_fracs[@]} * ${#sim_vzs[@]} * fanout_shards_per_axis * ${#uepipe_modes[@]} ))
    final_cfg_n=$(( ${#sim_submit_pts[@]} * ${#sim_submit_fracs[@]} * ${#sim_vzs[@]} * ${#sim_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
  else
    n_cfg=$(( ${#sim_submit_pts[@]} * ${#sim_submit_fracs[@]} * ${#sim_vzs[@]} * ${#sim_cones[@]} * iso_submit_n * ${#uepipe_modes[@]} ))
    final_cfg_n=$(( ${#sim_submit_pts[@]} * ${#sim_submit_fracs[@]} * ${#sim_vzs[@]} * ${#sim_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
  fi

  say "CHECKJOBS (dataset=${DATASET}, tag=${TAG})"
  say "  YAML master         : ${master_yaml}"
  say "  groupSize           : ${gs}"
  say "  samples             : [${samples[*]}]"
  say "  master list source  : ${SIM_ROOT}"
  echo

  say "${BOLD}Matrix dimensions:${RST}"
  if jetpt_internal_enabled; then
    say "  jet_pt_min                        : [${sim_pts[*]}]  (${#sim_pts[@]} values; internal hist scan)"
  else
    say "  jet_pt_min                        : [${sim_pts[*]}]  (${#sim_pts[@]} cfg-tag values)"
  fi
  if dphi_internal_enabled; then
    say "  back_to_back_dphi_min_pi_fraction : [${sim_fracs[*]}]  (${#sim_fracs[@]} values; internal hist scan)"
  else
    say "  back_to_back_dphi_min_pi_fraction : [${sim_fracs[*]}]  (${#sim_fracs[@]} cfg-tag values)"
  fi
  say "  vz_cut_cm                         : [${sim_vzs[*]}]  (${#sim_vzs[@]} cfg-tag value(s); $(vz_selection_reason_for_yaml "$master_yaml"))"
  say_vz_selection_summary "$master_yaml" "${sim_vzs[@]}"
  if iso_view_internal_enabled; then
    say "  coneR                             : [${sim_view_cones[*]}]  (${#sim_view_cones[@]} internal iso/cone view values; submit scalar=${sim_cones[*]})"
    say "  iso/cone views                    : $(iso_view_env_value)"
  else
    say "  coneR                             : [${sim_cones[*]}]  (${#sim_cones[@]} values)"
  fi
  if id_fanout_enabled; then
    if iso_cone_fanout_enabled; then
      say "  iso/cone fanout shards            : ${fanout_shards_per_axis} upstream shard(s) per pt/dphi/vz/UE axis (cap=$(id_fanout_max_rows) outputs/pass)"
      say "  final cone×iso×ID outputs         : $(( ${#sim_cones[@]} * ${#iso_tags[@]} )) cfg tag(s) per pt/dphi/vz/UE axis"
    else
      say "  photon-ID fanout shards           : ${iso_submit_n} upstream shard(s), cap=${RJ_ID_FANOUT_MAX_ROWS:-15} cfg outputs/pass"
      say "  photon-ID cfg outputs             : ${#iso_tags[@]} final cfg ROOT file(s)"
      iso_view_internal_enabled && say "  final ROOT layout                  : ${#iso_tags[@]} cfg file(s); each contains the configured suffixed iso/cone histogram views"
    fi
  else
    say "  photon-ID modes submitted         : ${iso_submit_n} independent cfg tag(s) (fanout disabled)"
  fi
  if (( ${uepipe_in_tag:-0} )); then
    say "  clusterUEpipeline                 : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} values; tagged)"
  else
    say "  clusterUEpipeline                 : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} forced value; not tagged for ${TAG})"
  fi
  say "  upstream DST job configs          : ${BOLD}${n_cfg}${RST}"
  say "  final output cfg tags             : ${BOLD}${final_cfg_n}${RST}"
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

  say "${BOLD}Upstream DST-pass cell list (${n_cfg} entries × ${#samples[@]} samples = ${BOLD}$((n_cfg * ${#samples[@]}))${RST} submit blocks):${RST}"
  say "  These are the expensive Fun4All/DST passes. In fanout mode a representative row can write many final cfg-tag ROOT outputs."
  printf "  ${BOLD}%3s │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s${RST}\n" \
         "#" "representative_cfg_or_shard" "pt" "frac" "vz" "cone" "iso" "uepipe"
  printf "  ────┼────────────────────────────────────────────────────────────────────────┼───── ────── ─────── ─────── ────────────────── ────────────\n"
  local cfg_num=0
  for _pt in "${sim_submit_pts[@]}"; do
    for _frac in "${sim_submit_fracs[@]}"; do
      for _vz in "${sim_vzs[@]}"; do
        if iso_cone_fanout_enabled; then
          for _ue in "${uepipe_modes[@]}"; do
            local _shards _sh
            _shards="$(direct_fanout_shard_count "${#sim_cones[@]}" "${#iso_tags[@]}")"
            for (( _sh=1; _sh<=_shards; _sh++ )); do
              (( cfg_num+=1 ))
              printf "  %3d │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s\n" \
                     "$cfg_num" "isoConeFanoutShard$(printf '%03d' "$_sh")" "$_pt" "$_frac" "$_vz" "scan" "cone×iso×ID" "$_ue"
            done
          done
        else
          for _cone in "${sim_cones[@]}"; do
            for (( _ci=0; _ci<${#iso_tags[@]}; _ci++ )); do
              id_fanout_enabled && ! iso_idx_is_group_leader "$_ci" && continue
              for _ue in "${uepipe_modes[@]}"; do
                (( cfg_num+=1 ))
                local _tag
                _tag="$(matrix_cfg_tag "$_pt" "$_frac" "$_vz" "$_cone" "$_ci" "$_ue")"
                printf "  %3d │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s\n" \
                       "$cfg_num" "$_tag" "$_pt" "$_frac" "$_vz" "$_cone" "${iso_tags[$_ci]}" "$_ue"
              done
            done
          done
        fi
      done
    done
  done
  echo

  say "${BOLD}Final ROOT output cfg tags written by those upstream passes:${RST}"
  if iso_view_internal_enabled; then
    say "  Layout: one cfg ROOT file per photon-ID triplet; each file contains the configured suffixed iso/cone views:"
    say "          $(iso_view_env_value)"
  elif iso_cone_fanout_enabled; then
    say "  Layout: one cfg ROOT file per cone × iso × photon-ID output."
  elif id_fanout_enabled; then
    say "  Layout: one cfg ROOT file per photon-ID output in each upstream iso group."
  else
    say "  Layout: one scalar legacy cfg ROOT file per upstream pass."
  fi
  local _max_tags _shown_tags _total_tags_expected
  _max_tags="$(checkjobs_max_output_tags)"
  _shown_tags=0
  _total_tags_expected=0
  printf "  ${BOLD}%3s │ %-80s │ %-28s │ %-12s${RST}\n" "#" "final_cfg_tag" "photon_id_triplet" "sample dirs"
  printf "  ────┼──────────────────────────────────────────────────────────────────────────────────┼──────────────────────────────┼────────────\n"
  for _pt in "${sim_submit_pts[@]}"; do
    for _frac in "${sim_submit_fracs[@]}"; do
      for _vz in "${sim_vzs[@]}"; do
        for _cone in "${sim_cones[@]}"; do
          for (( _ci=0; _ci<${#iso_tags[@]}; _ci++ )); do
            for _ue in "${uepipe_modes[@]}"; do
              (( _total_tags_expected+=1 ))
              (( _shown_tags >= _max_tags )) && continue
              local _tag _trip
              _tag="$(matrix_cfg_tag "$_pt" "$_frac" "$_vz" "$_cone" "$_ci" "$_ue")"
              _trip="${iso_preselection[$_ci]}/${iso_tight[$_ci]}/${iso_nonTight[$_ci]}"
              (( _shown_tags+=1 ))
              printf "  %3d │ %-80s │ %-28s │ %-12s\n" "$_shown_tags" "$_tag" "$_trip" "${#samples[@]} sample(s)"
            done
          done
        done
      done
    done
  done
  if (( _total_tags_expected > _shown_tags )); then
    say "  ... showing ${_shown_tags}/${_total_tags_expected}; set RJ_CHECKJOBS_MAX_OUTPUT_TAGS=${_total_tags_expected} to print all."
  fi
  echo

  total_jobs=$(( n_cfg * per_cfg_jobs ))
  say "${BOLD}Job count summary:${RST}"
  say "  upstream configs       : ${n_cfg}"
  say "  jobs per combo (Σsamp) : ${per_cfg_jobs}"
  if id_fanout_enabled; then
    say "  final cfg outputs      : ${final_cfg_n} per sample set; written by fanout, not extra DST passes"
  else
    say "  final cfg outputs      : ${final_cfg_n}; one output cfg per upstream pass"
  fi
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
  _ex_tag="$(matrix_cfg_tag "${sim_submit_pts[0]}" "${sim_submit_fracs[0]}" "${sim_vzs[0]}" "${sim_cones[0]}" 0 "${uepipe_modes[0]}")"
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
  local -a ck_pts ck_submit_pts ck_fracs ck_submit_fracs ck_vzs ck_cones
  mapfile -t ck_pts   < <( yaml_get_values "jet_pt_min" "$data_yaml_src" )
  mapfile -t ck_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$data_yaml_src" )
  mapfile -t ck_vzs   < <( dataset_vz_values "$data_yaml_src" )
  mapfile -t ck_cones < <( yaml_get_values "coneR" "$data_yaml_src" )
  (( ${#ck_pts[@]} ))   || { err "No values found for jet_pt_min in $data_yaml_src"; exit 72; }
  (( ${#ck_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $data_yaml_src"; exit 72; }
  (( ${#ck_vzs[@]} ))   || { err "No values found for vz_cut_cm in $data_yaml_src"; exit 72; }
  (( ${#ck_cones[@]} )) || { err "No values found for coneR in $data_yaml_src"; exit 72; }
  local -a ck_view_cones=( "${ck_cones[@]}" )
  mapfile -t ck_cones < <( cone_submit_values "${ck_view_cones[@]}" )
  mapfile -t ck_submit_pts < <( jetpt_submit_values "${ck_pts[@]}" )
  mapfile -t ck_submit_fracs < <( dphi_submit_values "${ck_fracs[@]}" )
  build_iso_modes "$data_yaml_src"
  read_uepipe_modes "$data_yaml_src" "$TAG"

  local total_runs=0 total_files=0 base_jobs=0 missing=0

  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    local r8; r8="$(run8 "$rn")"
    local lf="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"

    if [[ ! -s "$lf" ]]; then
      warn "No per-run list for ${r8}; skipping"
      ((missing += 1))
      continue
    fi

    local nfiles; nfiles=$(wc -l < "$lf" | awk '{print $1}')
    local nj;     nj=$(ceil_div "$nfiles" "$gs")

    (( total_runs  += 1    ))
    (( total_files += nfiles ))
    (( base_jobs   += nj   ))
  done < "$GOLDEN"

  local iso_submit_n
  iso_submit_n="$(iso_submit_count)"
  local fanout_shards_per_axis=0
  local final_cfg_n=0
  local n_matrix=0
  if iso_cone_fanout_enabled; then
    fanout_shards_per_axis="$(direct_fanout_shard_count "${#ck_cones[@]}" "${#iso_tags[@]}")"
    n_matrix=$(( ${#ck_submit_pts[@]} * ${#ck_submit_fracs[@]} * ${#ck_vzs[@]} * fanout_shards_per_axis * ${#uepipe_modes[@]} ))
    final_cfg_n=$(( ${#ck_submit_pts[@]} * ${#ck_submit_fracs[@]} * ${#ck_vzs[@]} * ${#ck_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
  else
    n_matrix=$(( ${#ck_submit_pts[@]} * ${#ck_submit_fracs[@]} * ${#ck_vzs[@]} * ${#ck_cones[@]} * iso_submit_n * ${#uepipe_modes[@]} ))
    final_cfg_n=$(( ${#ck_submit_pts[@]} * ${#ck_submit_fracs[@]} * ${#ck_vzs[@]} * ${#ck_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
  fi
  local total_jobs=$(( n_matrix * base_jobs ))

  say "CHECKJOBS (dataset=${DATASET}, tag=${TAG})"
  say "  YAML source          : ${data_yaml_src}"
  say "  groupSize            : ${gs}"
  echo
  say "${BOLD}Matrix dimensions:${RST}"
  if jetpt_internal_enabled; then
    say "  jet_pt_min           : [${ck_pts[*]}]  (${#ck_pts[@]} values; internal hist scan)"
  else
    say "  jet_pt_min           : [${ck_pts[*]}]  (${#ck_pts[@]} cfg-tag values)"
  fi
  if dphi_internal_enabled; then
    say "  b2b_dphi_pi_fraction : [${ck_fracs[*]}]  (${#ck_fracs[@]} values; internal hist scan)"
  else
    say "  b2b_dphi_pi_fraction : [${ck_fracs[*]}]  (${#ck_fracs[@]} cfg-tag values)"
  fi
  say "  vz_cut_cm            : [${ck_vzs[*]}]  (${#ck_vzs[@]} cfg-tag value(s); $(vz_selection_reason_for_yaml "$data_yaml_src"))"
  say_vz_selection_summary "$data_yaml_src" "${ck_vzs[@]}"
  if iso_view_internal_enabled; then
    say "  coneR                : [${ck_view_cones[*]}]  (${#ck_view_cones[@]} internal iso/cone view values; submit scalar=${ck_cones[*]})"
    say "  iso/cone views       : $(iso_view_env_value)"
  else
    say "  coneR                : [${ck_cones[*]}]  (${#ck_cones[@]} values)"
  fi
  if id_fanout_enabled; then
    if iso_cone_fanout_enabled; then
      say "  iso/cone fanout      : ${fanout_shards_per_axis} upstream shard(s) per pt/dphi/vz/UE axis (cap=$(id_fanout_max_rows) outputs/pass)"
      say "  final cone×iso×ID cfgs: $(( ${#ck_cones[@]} * ${#iso_tags[@]} )) cfg tag(s) per pt/dphi/vz/UE axis"
    else
      say "  photon-ID fanout     : ${iso_submit_n} upstream shard(s), cap=${RJ_ID_FANOUT_MAX_ROWS:-15} cfg outputs/pass"
      say "  photon-ID cfg outputs: ${#iso_tags[@]} final cfg ROOT file(s)"
      iso_view_internal_enabled && say "  final ROOT layout     : ${#iso_tags[@]} cfg file(s); each contains the configured suffixed iso/cone histogram views"
    fi
  else
    say "  photon-ID cfgs       : ${iso_submit_n} independent cfg tag(s) (fanout disabled)"
  fi
  if (( ${uepipe_in_tag:-0} )); then
    say "  clusterUEpipeline    : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} values; tagged)"
  else
    say "  clusterUEpipeline    : [${uepipe_modes[*]}]  (${#uepipe_modes[@]} forced value; not tagged for ${TAG})"
  fi
  say "  upstream matrix cells: ${BOLD}${n_matrix}${RST}"
  say "  final output cfg tags: ${BOLD}${final_cfg_n}${RST}"
  echo
  say "${BOLD}Per-matrix-cell base:${RST}"
  say "  golden runs (w/lists): ${total_runs}"
  say "  total input files    : ${total_files}"
  say "  base jobs per combo  : ${base_jobs}"
  (( missing > 0 )) && warn "  runs skipped (missing lists): ${missing}"
  echo

  say "${BOLD}Upstream DST-pass cell list (${n_matrix} entries):${RST}"
  say "  These are the expensive Fun4All/DST passes. In fanout mode a representative row can write many final cfg-tag ROOT outputs."
  echo
  printf "  ${BOLD}%3s │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s${RST}\n" \
         "#" "representative_cfg_or_shard" "pt" "frac" "vz" "cone" "iso" "uepipe"
  printf "  ────┼────────────────────────────────────────────────────────────────────────┼───── ────── ─────── ─────── ────────────────── ────────────\n"
  local cfg_num=0
  for _pt in "${ck_submit_pts[@]}"; do
  for _frac in "${ck_submit_fracs[@]}"; do
  for _vz in "${ck_vzs[@]}"; do
    if iso_cone_fanout_enabled; then
      for _ue in "${uepipe_modes[@]}"; do
        local _shards _sh
        _shards="$(direct_fanout_shard_count "${#ck_cones[@]}" "${#iso_tags[@]}")"
        for (( _sh=1; _sh<=_shards; _sh++ )); do
          (( ++cfg_num ))
          printf "  %3d │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s\n" \
                 "$cfg_num" "isoConeFanoutShard$(printf '%03d' "$_sh")" "$_pt" "$_frac" "$_vz" "scan" "cone×iso×ID" "$_ue"
        done
      done
    else
      for _cone in "${ck_cones[@]}"; do
      for (( _ci=0; _ci<${#iso_tags[@]}; _ci++ )); do
      id_fanout_enabled && ! iso_idx_is_group_leader "$_ci" && continue
      for _ue in "${uepipe_modes[@]}"; do
        (( ++cfg_num ))
        local _tag
        _tag="$(matrix_cfg_tag "$_pt" "$_frac" "$_vz" "$_cone" "$_ci" "$_ue")"
        printf "  %3d │ %-70s │ %-5s %-6s %-7s %-7s %-18s %-12s\n" \
               "$cfg_num" "$_tag" "$_pt" "$_frac" "$_vz" "$_cone" "${iso_tags[$_ci]}" "$_ue"
      done
      done
      done
    fi
  done
  done
  done
  echo

  say "${BOLD}Final ROOT output cfg tags written by those upstream passes:${RST}"
  if iso_view_internal_enabled; then
    say "  Layout: one cfg ROOT file per photon-ID triplet; each file contains the configured suffixed iso/cone views:"
    say "          $(iso_view_env_value)"
  elif iso_cone_fanout_enabled; then
    say "  Layout: one cfg ROOT file per cone × iso × photon-ID output."
  elif id_fanout_enabled; then
    say "  Layout: one cfg ROOT file per photon-ID output in each upstream iso group."
  else
    say "  Layout: one scalar legacy cfg ROOT file per upstream pass."
  fi
  local _max_tags _shown_tags _total_tags_expected
  _max_tags="$(checkjobs_max_output_tags)"
  _shown_tags=0
  _total_tags_expected=0
  printf "  ${BOLD}%3s │ %-80s │ %-28s │ %-18s${RST}\n" "#" "final_cfg_tag" "photon_id_triplet" "output tree"
  printf "  ────┼──────────────────────────────────────────────────────────────────────────────────┼──────────────────────────────┼──────────────────\n"
  for _pt in "${ck_submit_pts[@]}"; do
  for _frac in "${ck_submit_fracs[@]}"; do
  for _vz in "${ck_vzs[@]}"; do
    for _cone in "${ck_cones[@]}"; do
    for (( _ci=0; _ci<${#iso_tags[@]}; _ci++ )); do
    for _ue in "${uepipe_modes[@]}"; do
      (( _total_tags_expected+=1 ))
      (( _shown_tags >= _max_tags )) && continue
      local _tag _trip
      _tag="$(matrix_cfg_tag "$_pt" "$_frac" "$_vz" "$_cone" "$_ci" "$_ue")"
      _trip="${iso_preselection[$_ci]}/${iso_tight[$_ci]}/${iso_nonTight[$_ci]}"
      (( _shown_tags+=1 ))
      printf "  %3d │ %-80s │ %-28s │ %-18s\n" "$_shown_tags" "$_tag" "$_trip" "<run8>/*.root"
    done
    done
    done
  done
  done
  done
  if (( _total_tags_expected > _shown_tags )); then
    say "  ... showing ${_shown_tags}/${_total_tags_expected}; set RJ_CHECKJOBS_MAX_OUTPUT_TAGS=${_total_tags_expected} to print all."
  fi
  echo

  say "${BOLD}Job count summary (groupSize=${gs}):${RST}"
  say "  upstream matrix cells  : ${n_matrix}"
  say "  base jobs per combo    : ${base_jobs}"
  if id_fanout_enabled; then
    say "  final cfg outputs      : ${final_cfg_n}; written by fanout, not extra DST passes"
  else
    say "  final cfg outputs      : ${final_cfg_n}; one output cfg per upstream pass"
  fi
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
  _ex_tag="$(matrix_cfg_tag "${ck_submit_pts[0]}" "${ck_submit_fracs[0]}" "${ck_vzs[0]}" "${ck_cones[0]}" 0 "${uepipe_modes[0]}")"
  say "  e.g. ${DIM}${DEST_BASE}/${_ex_tag}/00067599/*.root${RST}"
}

workflow_check() {
  local yaml_src
  yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  [[ -s "$yaml_src" ]] || { err "Workflow check YAML not found or empty: $yaml_src"; exit 72; }

  local -a pts submit_pts fracs submit_fracs vzs cones
  mapfile -t pts   < <( yaml_get_values "jet_pt_min" "$yaml_src" )
  mapfile -t fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$yaml_src" )
  mapfile -t vzs   < <( dataset_vz_values "$yaml_src" )
  mapfile -t cones < <( yaml_get_values "coneR" "$yaml_src" )
  (( ${#pts[@]} ))   || { err "No values found for jet_pt_min in $yaml_src"; exit 72; }
  (( ${#fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $yaml_src"; exit 72; }
  (( ${#vzs[@]} ))   || { err "No values found for vz_cut_cm in $yaml_src"; exit 72; }
  (( ${#cones[@]} )) || { err "No values found for coneR in $yaml_src"; exit 72; }
  mapfile -t cones < <( cone_submit_values "${cones[@]}" )
  mapfile -t submit_pts < <( jetpt_submit_values "${pts[@]}" )
  mapfile -t submit_fracs < <( dphi_submit_values "${fracs[@]}" )

  build_iso_modes "$yaml_src"
  read_uepipe_modes "$yaml_src" "$TAG"

  local iso_submit_n
  iso_submit_n="$(iso_submit_count)"
  local n_cfg=$(( ${#submit_pts[@]} * ${#submit_fracs[@]} * ${#vzs[@]} * ${#cones[@]} * iso_submit_n * ${#uepipe_modes[@]} ))
  local ex_tag
  ex_tag="$(matrix_cfg_tag "${submit_pts[0]}" "${submit_fracs[0]}" "${vzs[0]}" "${cones[0]}" 0 "${uepipe_modes[0]}")"
  local override
  override="$(sim_make_yaml_override "$yaml_src" "${submit_pts[0]}" "${submit_fracs[0]}" "${vzs[0]}" "${cones[0]}" "${iso_sliding[0]}" "${iso_fixed[0]}" "${uepipe_modes[0]}" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}" "$ex_tag" "WORKFLOWCHECK" "false")"

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
  if id_fanout_enabled; then
    say "  production mode  : direct RecoilJets histogram machinery (${#iso_tags[@]} final cfg ROOT files; ${iso_submit_n} photon-ID fanout shard(s), cap=${RJ_ID_FANOUT_MAX_ROWS:-15})"
  else
    say "  production mode  : direct legacy one-cfg-per-job validation (ID fanout disabled)"
  fi
}

# Submit a set of runs (from a round file or the whole golden list) to Condor
#   submit_condor <runs_source> [firstChunk]
submit_condor() {
  local source="$1"
  local first_chunk="${2:-}"
  local direct_max_jobs="${RJ_DIRECT_MAX_JOBS:-0}"

  [[ -s "$source" ]] || { err "Run source not found or empty: $source"; exit 5; }
  [[ "$direct_max_jobs" =~ ^[0-9]+$ ]] || {
    err "RJ_DIRECT_MAX_JOBS must be a non-negative integer, got '${direct_max_jobs}'"
    exit 2
  }
  if ppg12_pp_strict_list_coverage_enabled; then
    validate_ppg12_pp_list_coverage "$source" || exit 88
  fi

  # Clean stale .sub files for this TAG only in direct-submit mode. In
  # orchestration mode previous analysis nodes are still referenced by the
  # parent DAG being built.
  if ! dag_collect_enabled; then
    rm -f "${SUB_DIR}/RecoilJets_${TAG}_"*.sub 2>/dev/null || true
  fi

  local stamp; stamp="$(date +%Y%m%d_%H%M%S)"
  local sub="${SUB_DIR}/RecoilJets_${TAG}_${stamp}.sub"
  local args_file="${SUB_DIR}/RecoilJets_${TAG}_${stamp}.args"
  local submit_stage_dir="${STAGE_DIR}/${TAG}_${stamp}"
  local exe_to_use="${BULK_FROZEN_EXE:-${EXE}}"
  local macro_env=""
  [[ -n "${BULK_FROZEN_MACRO:-}" ]] && macro_env=";RJ_MACRO_PATH=${BULK_FROZEN_MACRO}"
  local submit_extra_env
  submit_extra_env="$(build_submit_extra_env_fragment)"
  local request_memory
  request_memory="$(memory_request_from_env_or_default "2000MB")"
  local request_memory_mb
  request_memory_mb="$(memory_request_to_mb "$request_memory")"
  local direct_nevents
  direct_nevents="$(direct_worker_nevents)"
  mkdir -p "$submit_stage_dir"
  : > "$args_file"
  say "Submit chunk-list stage: ${submit_stage_dir}"

  # Snapshot YAML at submit time so idle jobs are immune to later edits
  local yaml_src="${RJ_CONFIG_YAML:-${SIM_YAML_DEFAULT}}"
  local yaml_snap="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${stamp}.yaml"
  local source_runs
  local require_non_tiny_output
  source_runs=$(grep -cE '^[0-9]+' "$source" 2>/dev/null || true)
  require_non_tiny_output="$(effective_require_non_tiny_output)"
  mkdir -p "$SIM_YAML_OVERRIDE_DIR"
  cp -f "$yaml_src" "$yaml_snap"
  force_ppg12_photon_yield_yaml_contract "$yaml_snap"
  say "YAML snapshot: ${yaml_snap}"
  say "Submit context: source=${source}  runs=${source_runs:-0}  groupSize=${GROUP_SIZE}  nEvents=${direct_nevents}  firstChunk=${first_chunk:-none}"
  say "Submit environment: RJ_DATASET=${DATASET}  RJ_VERBOSITY=0  RJ_CONFIG_YAML=${yaml_snap}${macro_env}${submit_extra_env};RJ_PROFILE_JOB=${RJ_PROFILE_JOB:-0};RJ_PROFILE_STAGE=${RJ_PROFILE_STAGE:-direct};RJ_REQUEST_MEMORY_MB=${request_memory_mb};RJ_REQUIRE_NON_TINY_OUTPUT=${require_non_tiny_output};RJ_MIN_OUTPUT_BYTES=${RJ_MIN_OUTPUT_BYTES:-50000};RJ_FAIL_ON_MISSING_CALO_INPUT=${RJ_FAIL_ON_MISSING_CALO_INPUT:-0}"

  cat > "$sub" <<SUB
universe      = vanilla
executable    = ${exe_to_use}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/job.\$(Cluster).\$(Process).err
$(condor_auto_memory_retry_block "$request_memory_mb")
$(condor_worker_failure_hold_block)
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
# Force dataset & quiet macro on Condor (YAML frozen at submit time):
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=0;RJ_CONFIG_YAML=${yaml_snap}${macro_env}${submit_extra_env};RJ_PROFILE_JOB=${RJ_PROFILE_JOB:-0};RJ_JOB_HEARTBEAT_SECONDS=${RJ_JOB_HEARTBEAT_SECONDS:-0};RJ_PROFILE_STAGE=${RJ_PROFILE_STAGE:-direct};RJ_PROFILE_LABEL=${RJ_PROFILE_LABEL:-${TAG}};RJ_REQUEST_MEMORY_MB=${request_memory_mb};RJ_REQUIRE_NON_TINY_OUTPUT=${require_non_tiny_output};RJ_MIN_OUTPUT_BYTES=${RJ_MIN_OUTPUT_BYTES:-50000};RJ_FAIL_ON_MISSING_CALO_INPUT=${RJ_FAIL_ON_MISSING_CALO_INPUT:-0}
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

    if (( direct_max_jobs > 0 )); then
      local remaining_jobs=$(( direct_max_jobs - queued ))
      if (( remaining_jobs <= 0 )); then
        break
      fi
      if (( ${#groups[@]} > remaining_jobs )); then
        say "Capping direct DATA group list: ${#groups[@]} -> ${remaining_jobs} jobs (RJ_DIRECT_MAX_JOBS=${direct_max_jobs})"
        groups=( "${groups[@]:0:remaining_jobs}" )
      fi
    fi

    local gidx=0
    for glist in "${groups[@]}"; do
      (( gidx+=1 ))
      printf '%s %s %s $(Cluster) %d %d NONE %s\n' \
             "$r8" "$glist" "$DATASET" "$direct_nevents" "$gidx" "$DEST_BASE" >> "$args_file"
      (( queued+=1 ))
    done

    if (( direct_max_jobs > 0 && queued >= direct_max_jobs )); then
      say "Reached direct DATA job cap: ${queued}/${direct_max_jobs}"
      break
    fi

  done < "$source"

  local submit_elapsed
  submit_elapsed=$(( $(date +%s) - t0 ))
  if dag_collect_enabled; then
    say "Queueing ${BOLD}${queued}${RST} jobs as DAG node → $(basename "$sub")"
  else
    say "Submitting ${BOLD}${queued}${RST} jobs → $(basename "$sub")"
  fi
  say "Submit summary: runsProcessed=${run_counter}  groupSize=${GROUP_SIZE}  firstChunk=${first_chunk:-none}  elapsed=${submit_elapsed}s"
  say "Wrapper path: ${exe_to_use}"
  submit_or_collect_condor "$sub" "analysis_${TAG}"
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

select_explicit_stat_data_run() {
  local requested="${1:?requested data run required}"
  local out_file="${2:?output run list required}"
  local stats_file="${3:-}"
  [[ "$requested" =~ ^[0-9]+$ ]] || {
    err "RJ_SMOKE_DATA_RUN must be numeric, got '${requested}'"
    return 2
  }
  local requested8
  requested8="$(run8 "$requested")"
  local found=0 rn r8 src nfiles
  while IFS= read -r rn; do
    [[ -z "$rn" || "$rn" =~ ^# ]] && continue
    r8="$(run8 "$rn")"
    [[ "$r8" == "$requested8" ]] || continue
    if [[ -n "${TRIGGER_BIT}" ]] && ! is_trigger_active "$r8" "$TRIGGER_BIT"; then
      err "Explicit smoke run ${r8} is outside the active trigger contract."
      return 2
    fi
    src="${LIST_DIR}/${LIST_PREFIX}-${r8}.list"
    [[ -s "$src" ]] || {
      err "Explicit smoke run ${r8} has no nonempty paired input list: ${src}"
      return 2
    }
    nfiles="$(grep -Evc '^[[:space:]]*($|#)' "$src" 2>/dev/null || echo 0)"
    [[ "$nfiles" =~ ^[0-9]+$ && "$nfiles" -gt 0 ]] || {
      err "Explicit smoke run ${r8} has no usable paired input rows."
      return 2
    }
    mkdir -p "$(dirname "$out_file")"
    printf '%s\n' "$r8" > "$out_file"
    if [[ -n "$stats_file" ]]; then
      mkdir -p "$(dirname "$stats_file")"
      printf 'run=%s input_files=%d explicit=1\n' "$r8" "$nfiles" > "$stats_file"
    fi
    found=1
    break
  done < "$GOLDEN"
  (( found == 1 )) || {
    err "Explicit smoke run ${requested8} is absent from the accepted golden run list ${GOLDEN}."
    return 2
  }
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
  mapfile -t vzs   < <( dataset_vz_values "$master_yaml" )
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
  mapfile -t vzs   < <( dataset_vz_values "$yaml" )
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
  mapfile -t vzs   < <( dataset_vz_values "$master_yaml" )
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
  mapfile -t vzs   < <( dataset_vz_values "$master_yaml" )
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
  [[ "$DATASET" == "isAuAu" ]] || { err "${ACTION:-scaledTriggerStudy} is valid only for isAuAu"; exit 2; }
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
  mapfile -t data_vzs   < <( dataset_vz_values "$master_yaml" )
  mapfile -t data_cones < <( yaml_get_values "coneR" "$master_yaml" )
  (( ${#data_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
  (( ${#data_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
  (( ${#data_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
  (( ${#data_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }

  local old_disable_iso_internal="${RJ_DISABLE_ISO_CONE_INTERNALIZATION-__unset__}"
  RJ_DISABLE_ISO_CONE_INTERNALIZATION=1
  build_iso_modes "$master_yaml"
  if [[ "$old_disable_iso_internal" == "__unset__" ]]; then
    unset RJ_DISABLE_ISO_CONE_INTERNALIZATION
  else
    RJ_DISABLE_ISO_CONE_INTERNALIZATION="$old_disable_iso_internal"
  fi
  read_uepipe_modes "$master_yaml" "$TAG"

  local data_pt="${SCALED_TRIGGER_JET_PT_MIN}"
  local data_frac="${SCALED_TRIGGER_DPHI_FRAC}"
  local data_vz="${SCALED_TRIGGER_VZ_CUT_CM}"
  local data_cone="${SCALED_TRIGGER_CONE_R}"
  local uepipe="${SCALED_TRIGGER_UEPIPELINE}"
  local data_sliding="${SCALED_TRIGGER_IS_SLIDING}"
  local data_fixed="${SCALED_TRIGGER_FIXED_GEV}"
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

  local study_action="${ACTION:-scaledTriggerStudy}"
  SCALED_TRIGGER_YAML="${SIM_YAML_OVERRIDE_DIR}/analysis_config_${TAG}_${SCALED_TRIGGER_CFG_TAG}_${study_action}.yaml"
  mkdir -p "$SIM_YAML_OVERRIDE_DIR"
  sed -E \
    -e "s|^([[:space:]]*jet_pt_min:).*|\\1 ${data_pt}|" \
    -e "s|^([[:space:]]*back_to_back_dphi_min_pi_fraction:).*|\\1 ${data_frac}|" \
    -e "s|^([[:space:]]*vz_cut_cm:).*|\\1 ${data_vz}|" \
    -e "s|^([[:space:]]*coneR:).*|\\1 ${data_cone}|" \
    -e "s|^([[:space:]]*isSlidingIso:).*|\\1 ${data_sliding}|" \
    -e "s|^([[:space:]]*isSlidingAndFixed:).*|\\1 false|" \
    -e "s|^([[:space:]]*fixedGeV:).*|\\1 ${data_fixed}|" \
    -e "s|^([[:space:]]*clusterUEpipeline:).*|\\1 ${uepipe}|" \
    "$master_yaml" > "$SCALED_TRIGGER_YAML"
  pin_photon_id_scalars_in_yaml "$SCALED_TRIGGER_YAML" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}"

  say "${study_action} single config:"
  say "  cfg_tag       : ${SCALED_TRIGGER_CFG_TAG}"
  say "  output tag    : ${SCALED_TRIGGER_OUTPUT_TAG}"
  say "  YAML override : ${SCALED_TRIGGER_YAML}"
  say "  run list      : ${SCALED_TRIGGER_RUNLIST} ($(grep -cE '^[0-9]+' "$SCALED_TRIGGER_RUNLIST") runs)"
  say "  pinned knobs  : jet_pt_min=${data_pt}, dphi/pi=${data_frac}, vz=${data_vz} cm, coneR=${data_cone}, sliding=${data_sliding}, UE=${uepipe}"
  if [[ "$study_action" == "scaledTriggerCentStudy" ]]; then
    say "  cent edges    : ${SCALED_TRIGGER_CENT_EDGES}"
  fi
}

# ------------------------ Parse CLI ------------------------
[[ $# -ge 1 ]] || usage
if [[ "${1:-}" == "orchestrationSelfTest" || "${2:-}" == "orchestrationSelfTest" ]]; then
  orchestration_self_test
  exit $?
fi
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
    scaledTriggerStudy|scaledTriggerCentStudy)
      ACTION="$tok"
      ;;
    local|localTest|condorDoAll|condorDoAllSmoke|condorDoAllDirect|condorDoAllFromScratch|condorHistFromPool|resume|smokeTest|condorExtract|trainFromExtraction|trainCentInput3x3FromExtraction|trainExpandedFromExtraction|trainExpandedFromExtractionCondor|applyCheck|validateOnSim|validateSim|simValidation|validateOnSimCondor|condorValidateOnSim|validateSimCondor|smokeTestFirstPass|smokeTestSecondPass|smokeTestApplyExisting)
      if [[ "$ACTION" == trainTightBDT || "$ACTION" == trainNPB || "$ACTION" == trainJetMLResidual || "$ACTION" == trainMLAll || "$ACTION" == scaledTriggerStudy || "$ACTION" == scaledTriggerCentStudy ]]; then
        TRAIN_MODE="$tok"
      elif [[ "$ACTION" == "condor" && "$tok" == "smokeTest" ]]; then
        :  # DATA smokeTest is a condor submode consumed positionally below.
      else
        ACTION="$tok"
      fi
      ;;
    checkModels|workflowCheck|isLocalIsoPing|orchestrationSelfTest|resolveGoldenRunList|condor|splitGoldenRunList|condorTest)
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
if [[ "$ACTION" == "scaledTriggerStudy" || "$ACTION" == "scaledTriggerCentStudy" ]]; then
  [[ "$DATASET" == "isAuAu" ]] || { err "${ACTION} is valid only as: $0 isAuAu ${ACTION} <local|condorDoAll>"; exit 2; }
  [[ -n "$TRAIN_MODE" ]] || TRAIN_MODE="local"
  [[ "$TRAIN_MODE" == "local" || "$TRAIN_MODE" == "condorDoAll" ]] || { err "${ACTION} mode must be local or condorDoAll, got '${TRAIN_MODE}'"; exit 2; }
  if [[ "$ACTION" == "scaledTriggerCentStudy" && -z "${RJ_SCALED_TRIGGER_OUTPUT_SUFFIX:-}" ]]; then
    SCALED_TRIGGER_OUTPUT_SUFFIX="$SCALED_TRIGGER_CENT_OUTPUT_SUFFIX"
  fi
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
  orchestrationSelfTest)
    orchestration_self_test
    exit $?
    ;;

  resolveGoldenRunList)
    [[ "$IS_SIM" -eq 0 ]] || { err "resolveGoldenRunList is only valid for data datasets"; exit 2; }
    resolve_source="${RJ_RESOLVE_RUN_SOURCE:-${GOLDEN}}"
    resolve_stamp="$(date +%Y%m%d_%H%M%S)"
    resolve_out="${RJ_RESOLVE_RUNLIST_OUT:-${SUB_DIR}/${TAG}_resolvedGolden_${resolve_stamp}_runs.list}"
    resolve_stats="${RJ_RESOLVE_STATS_OUT:-${SUB_DIR}/${TAG}_resolvedGolden_${resolve_stamp}_run_stats.txt}"
    say "${BOLD}DATA input resolver preflight requested${RST}"
    say "  dataset      : ${DATASET}"
    say "  source runs  : ${resolve_source}"
    say "  input lists  : ${LIST_DIR}"
    say "  output runs  : ${resolve_out}"
    say "  output stats : ${resolve_stats}"
    say "  requirement  : each run list must resolve required CALOFITTING inputs before submit"
    resolve_golden_run_list "$resolve_source" "$resolve_out" "$resolve_stats"
    ;;

  scaledTriggerStudy|scaledTriggerCentStudy)
    scaled_trigger_prepare_artifacts
    scaled_trigger_prepare_single_yaml
    scaled_trigger_cent_flag=0
    if [[ "$ACTION" == "scaledTriggerCentStudy" ]]; then
      scaled_trigger_cent_flag=1
    fi
    export RJ_CONFIG_YAML="$SCALED_TRIGGER_YAML"
    export RJ_SCALED_TRIGGER_RUNLIST="$SCALED_TRIGGER_RUNLIST"
    export RJ_SCALED_TRIGGER_STUDY_ONLY=1
    export RJ_SCALED_TRIGGER_CENT_STUDY="$scaled_trigger_cent_flag"
    export RJ_SCALED_TRIGGER_CENT_EDGES="$SCALED_TRIGGER_CENT_EDGES"
    export RJ_REQUEST_MEMORY="${RJ_SCALED_TRIGGER_REQUEST_MEMORY:-1000MB}"
    if [[ "$ACTION" == "scaledTriggerCentStudy" && -z "${RJ_SCALED_TRIGGER_REQUEST_MEMORY:-}" ]]; then
      export RJ_REQUEST_MEMORY="${RJ_SCALED_TRIGGER_CENT_REQUEST_MEMORY:-1500MB}"
    fi
    export RJ_SUBMIT_EXTRA_ENV="RJ_SCALED_TRIGGER_STUDY_ONLY=1;RJ_SCALED_TRIGGER_CENT_STUDY=${scaled_trigger_cent_flag};RJ_SCALED_TRIGGER_CENT_EDGES=${SCALED_TRIGGER_CENT_EDGES};RJ_SCALED_TRIGGER_RUNLIST=${SCALED_TRIGGER_RUNLIST};RJ_SCALED_TRIGGER_VZ_MAX_CM=${RJ_SCALED_TRIGGER_VZ_MAX_CM:-${SCALED_TRIGGER_VZ_CUT}}"

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

        say "${ACTION} local smoke test"
        say "  run          : ${r8}"
        say "  groupSize    : ${GROUP_SIZE}"
        say "  list chunk   : ${glist}"
        say "  events       : ${nevt}"
        say "  DEST_BASE    : ${DEST_BASE}"
        if [[ "$ACTION" == "scaledTriggerCentStudy" ]]; then
          say "  cent edges   : ${SCALED_TRIGGER_CENT_EDGES}"
        fi
        say "  memory target: local run (condor default would be ${RJ_REQUEST_MEMORY})"
        say "Invoking wrapper locally..."

        RJ_DATASET="$DATASET" RJ_VERBOSITY="$RJV" \
        RJ_CONFIG_YAML="$SCALED_TRIGGER_YAML" \
        RJ_SCALED_TRIGGER_STUDY_ONLY=1 \
        RJ_SCALED_TRIGGER_CENT_STUDY="$scaled_trigger_cent_flag" \
        RJ_SCALED_TRIGGER_CENT_EDGES="$SCALED_TRIGGER_CENT_EDGES" \
        RJ_SCALED_TRIGGER_RUNLIST="$SCALED_TRIGGER_RUNLIST" \
        RJ_SCALED_TRIGGER_VZ_MAX_CM="${RJ_SCALED_TRIGGER_VZ_MAX_CM:-${SCALED_TRIGGER_VZ_CUT}}" \
        bash "$EXE" "$r8" "$glist" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
        ;;

      condorDoAll)
        cleanup_bulk_snapshots_for_tag
        create_pipeline_snapshot "auau" "$(date +%Y%m%d_%H%M%S)"
        say "${ACTION} condorDoAll"
        say "  runs         : $(grep -cE '^[0-9]+' "$SCALED_TRIGGER_RUNLIST")"
        say "  groupSize    : ${GROUP_SIZE}"
        say "  request mem  : ${RJ_REQUEST_MEMORY}"
        say "  DEST_BASE    : ${DEST_BASE}"
        if [[ "$ACTION" == "scaledTriggerCentStudy" ]]; then
          say "  cent edges   : ${SCALED_TRIGGER_CENT_EDGES}"
        fi
        submit_condor "$SCALED_TRIGGER_RUNLIST" ""
        ;;
    esac
    ;;

  trainTightBDT|trainNPB)
    task="tight"
    [[ "$ACTION" == "trainNPB" ]] && task="npb"
    if [[ "$ACTION" == "trainTightBDT" ]]; then
      sidecar="${BASE}/scripts/auau_tight_bdt_pipeline.sh"
      [[ -f "$sidecar" ]] || { err "Missing sidecar tight-BDT workflow: $sidecar"; exit 2; }
      sidecar_args=()
      for tok in "${tokens[@]}"; do
        case "$tok" in
          trainTightBDT|local|localTest|smokeTest|smokeTestFirstPass|condorDoAll|condorExtract|smokeTestSecondPass|trainFromExtraction|trainCentInput3x3FromExtraction|trainExpandedFromExtraction|trainExpandedFromExtractionCondor|smokeTestApplyExisting|applyCheck|validateOnSim|validateSim|simValidation|validateOnSimCondor|condorValidateOnSim|validateSimCondor)
            ;;
          *) sidecar_args+=( "$tok" ) ;;
        esac
      done
      case "${TRAIN_MODE:-local}" in
        local|localTest)
          say "trainTightBDT is now sidecar-managed; forwarding to: ${sidecar} localTest"
          exec bash "$sidecar" localTest "${sidecar_args[@]}"
          ;;
        smokeTest|smokeTestFirstPass)
          say "trainTightBDT smoke extraction is now sidecar-managed; forwarding to: ${sidecar} smokeTest"
          exec bash "$sidecar" smokeTest "${sidecar_args[@]}"
          ;;
        condorDoAll|condorExtract)
          say "trainTightBDT Condor extraction is now sidecar-managed; forwarding to: ${sidecar} condorExtract"
          exec bash "$sidecar" condorExtract "${sidecar_args[@]}"
          ;;
        smokeTestSecondPass|trainFromExtraction)
          say "trainTightBDT model training is now sidecar-managed; forwarding to: ${sidecar} trainFromExtraction"
          exec bash "$sidecar" trainFromExtraction "${sidecar_args[@]}"
          ;;
        trainCentInput3x3FromExtraction)
          say "trainTightBDT 3x3-width centrality-input training is sidecar-managed; forwarding to: ${sidecar} trainCentInput3x3FromExtraction"
          exec bash "$sidecar" trainCentInput3x3FromExtraction "${sidecar_args[@]}"
          ;;
        trainExpandedFromExtraction)
          say "trainTightBDT expanded campaign planning/training is now sidecar-managed; forwarding to: ${sidecar} trainExpandedFromExtraction"
          exec bash "$sidecar" trainExpandedFromExtraction "${sidecar_args[@]}"
          ;;
        trainExpandedFromExtractionCondor)
          say "trainTightBDT expanded campaign Condor training is now sidecar-managed; forwarding to: ${sidecar} trainExpandedFromExtractionCondor"
          exec bash "$sidecar" trainExpandedFromExtractionCondor "${sidecar_args[@]}"
          ;;
        smokeTestApplyExisting|applyCheck)
          say "trainTightBDT model apply check is now sidecar-managed; forwarding to: ${sidecar} applyCheck"
          exec bash "$sidecar" applyCheck "${sidecar_args[@]}"
          ;;
        validateOnSim|validateSim|simValidation)
          say "trainTightBDT simulation validation is now sidecar-managed; forwarding to: ${sidecar} validateOnSim"
          exec bash "$sidecar" validateOnSim "${sidecar_args[@]}"
          ;;
        validateOnSimCondor|condorValidateOnSim|validateSimCondor)
          say "trainTightBDT Condor simulation validation is now sidecar-managed; forwarding to: ${sidecar} validateOnSimCondor"
          exec bash "$sidecar" validateOnSimCondor "${sidecar_args[@]}"
          ;;
        *)
          err "trainTightBDT mode must be localTest, smokeTest, condorExtract, trainFromExtraction, trainCentInput3x3FromExtraction, trainExpandedFromExtraction, trainExpandedFromExtractionCondor, applyCheck, validateOnSim, or validateOnSimCondor; got '${TRAIN_MODE:-local}'"
          exit 2
          ;;
      esac
    fi
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
    mapfile -t sim_vzs   < <( dataset_vz_values "$master_yaml" )
    mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
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

    cleanup_local_artifacts_before_local_run

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
      mapfile -t sim_vzs   < <( dataset_vz_values "$master_yaml" )
      mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
      (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
      (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
      (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
      (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
      mapfile -t sim_cones < <( cone_submit_values "${sim_cones[@]}" )
      mapfile -t sim_submit_pts < <( jetpt_submit_values "${sim_pts[@]}" )
      mapfile -t sim_submit_fracs < <( dphi_submit_values "${sim_fracs[@]}" )
      build_iso_modes "$master_yaml"
      read_uepipe_modes "$master_yaml" "$TAG"

      samples=()
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
        case "$DATASET" in
          isSimEmbedded)          samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
          isSimEmbeddedInclusive) mapfile -t samples < <(simembeddedinclusive_sample_list) ;;
          isSimInclusive|isSimJet5) samples=( "run28_jet5" "run28_jet8" "run28_jet12" "run28_jet20" "run28_jet30" "run28_jet40" ) ;;
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

      for pt in "${sim_submit_pts[@]}"; do
        for frac in "${sim_submit_fracs[@]}"; do
          for vz in "${sim_vzs[@]}"; do
          for cone in "${sim_cones[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
          for uepipe in "${uepipe_modes[@]}"; do
          SIM_CFG_TAG="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$iso_idx" "$uepipe")"
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

              RJ_VERBOSITY="$RJV" RJ_CONFIG_YAML="$yaml_override" \
                RJ_INTERNAL_ISO_VIEWS="$(iso_view_env_value)" \
                RJ_ALLOW_FIXED_RECO_ISO_VIEWS="${RJ_ALLOW_FIXED_RECO_ISO_VIEWS:-0}" \
                bash "$EXE" "$SIM_SAMPLE" "$tmp" "$DATASET" LOCAL "$nevt" 1 NONE "$DEST_BASE"
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

                RJ_VERBOSITY="$RJV" RJ_CONFIG_YAML="$yaml_override" \
                  RJ_INTERNAL_ISO_VIEWS="$(iso_view_env_value)" \
                  RJ_ALLOW_FIXED_RECO_ISO_VIEWS="${RJ_ALLOW_FIXED_RECO_ISO_VIEWS:-0}" \
                  bash "$EXE" "$SIM_SAMPLE" "$tmp" "$DATASET" LOCAL "$nevt" "$local_file_idx" NONE "$DEST_BASE"
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
      mapfile -t data_vzs   < <( dataset_vz_values "$data_yaml_src" )
      mapfile -t data_cones < <( yaml_get_values "coneR" "$data_yaml_src" )
      (( ${#data_pts[@]} ))   || { err "No values found for jet_pt_min in $data_yaml_src"; exit 72; }
      (( ${#data_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $data_yaml_src"; exit 72; }
      (( ${#data_vzs[@]} ))   || { err "No values found for vz_cut_cm in $data_yaml_src"; exit 72; }
      (( ${#data_cones[@]} )) || { err "No values found for coneR in $data_yaml_src"; exit 72; }
      mapfile -t data_cones < <( cone_submit_values "${data_cones[@]}" )
      mapfile -t data_submit_pts < <( jetpt_submit_values "${data_pts[@]}" )
      mapfile -t data_submit_fracs < <( dphi_submit_values "${data_fracs[@]}" )
      build_iso_modes "$data_yaml_src"
      read_uepipe_modes "$data_yaml_src" "$TAG"
      DATA_DEST_BASE_SAVED="${RJ_LOCAL_DATA_OUTPUT_BASE:-$DEST_BASE}"

      for data_pt in "${data_submit_pts[@]}"; do
      for data_frac in "${data_submit_fracs[@]}"; do
      for data_vz in "${data_vzs[@]}"; do
        for data_cone in "${data_cones[@]}"; do
        for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
        for uepipe in "${uepipe_modes[@]}"; do
        data_cfg_tag="$(matrix_cfg_tag "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe")"
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
        force_ppg12_photon_yield_yaml_contract "$yaml_override"
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
        RJ_INTERNAL_ISO_VIEWS="$(iso_view_env_value)" \
        RJ_ALLOW_FIXED_RECO_ISO_VIEWS="${RJ_ALLOW_FIXED_RECO_ISO_VIEWS:-0}" \
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
    mapfile -t data_vzs   < <( dataset_vz_values "$data_yaml_src" )
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
    force_ppg12_photon_yield_yaml_contract "$yaml_override"
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
    mapfile -t sim_vzs   < <( dataset_vz_values "$master_yaml" )
    mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
    mapfile -t sim_cones < <( cone_submit_values "${sim_cones[@]}" )
    mapfile -t sim_submit_pts < <( jetpt_submit_values "${sim_pts[@]}" )
    mapfile -t sim_submit_fracs < <( dphi_submit_values "${sim_fracs[@]}" )
    build_iso_modes "$master_yaml"
    read_uepipe_modes "$master_yaml" "$TAG"
    jetpt_env_for_sub="$(jetpt_env_fragment "${sim_pts[@]}")"
    dphi_env_for_sub="$(dphi_env_fragment "${sim_fracs[@]}")"
    iso_view_env_for_sub="$(iso_view_env_fragment)"
    stitch_env_for_sub="$(embedded_inclusive_stitch_env_fragment)"
    submit_extra_env_for_sub="$(build_submit_extra_env_fragment)"

    SIM_DEST_BASE_RESOLVED="$DEST_BASE"

    pt0="${sim_submit_pts[0]}"
    frac0="${sim_submit_fracs[0]}"
    vz0="${sim_vzs[0]}"
    cone0="${sim_cones[0]}"
    SIM_CFG_TAG="$(matrix_cfg_tag "$pt0" "$frac0" "$vz0" "$cone0" 0 "${uepipe_modes[0]}")"
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
$(condor_auto_memory_retry_block "2000")
$(condor_worker_failure_hold_block)
should_transfer_files = NO
stream_output = True
stream_error  = True
environment   = RJ_VERBOSITY=10;RJ_CONFIG_YAML=${yaml_override}${jetpt_env_for_sub}${dphi_env_for_sub}${iso_view_env_for_sub}${stitch_env_for_sub}${submit_extra_env_for_sub};RJ_SIM_SAMPLE=${SIM_SAMPLE};RJ_EMBEDDED_INCLUSIVE_JET_SAMPLE=${SIM_SAMPLE};RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES:-0};RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30:-0};RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES:-0};RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40:-0}
arguments     = ${SIM_SAMPLE} ${tmp} ${DATASET} \$(Cluster) 0 1 NONE ${DEST_BASE}
queue
SUB

    say "Submitting 1 ${DATASET} condorTest job (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG}) → $(basename "$sub")"
    say "Output ROOT dir: ${DEST_BASE}/${SIM_SAMPLE}"
    say "YAML override: ${yaml_override}"
    condor_submit "$sub"
    ;;

  condorDoAllSmoke)
    [[ "$IS_SIM" -eq 1 ]] || { err "condorDoAllSmoke is only valid for isSim variants"; exit 2; }
    if [[ "${GROUP_SIZE_EXPLICIT:-0}" -eq 0 ]]; then
      case "$DATASET" in
        isSimEmbedded)          GROUP_SIZE="${RJ_SMOKE_GROUPSIZE_EMBEDDED:-2}" ;;
        isSimEmbeddedInclusive) GROUP_SIZE="${RJ_SMOKE_GROUPSIZE_EMBEDDED_INCLUSIVE:-${RJ_SMOKE_GROUPSIZE_EMBEDDED:-2}}" ;;
        *)                      GROUP_SIZE="${RJ_SMOKE_GROUPSIZE_SIM:-2}" ;;
      esac
    fi
    if [[ "${MAX_JOBS_EXPLICIT:-0}" -eq 0 ]]; then
      case "$DATASET" in
        isSimEmbedded)          MAX_JOBS="${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE_EMBEDDED:-${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE:-12}}" ;;
        isSimEmbeddedInclusive) MAX_JOBS="${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE_EMBEDDED_INCLUSIVE:-${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE_EMBEDDED:-${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE:-12}}}" ;;
        *)                      MAX_JOBS="${RJ_SMOKE_SIM_MAX_JOBS_PER_SAMPLE:-12}" ;;
      esac
    fi
    if [[ -z "${RJ_REQUEST_MEMORY:-}" ]]; then
      case "$DATASET" in
        isSimEmbedded)          RJ_REQUEST_MEMORY="${RJ_SMOKE_REQUEST_MEMORY_EMBEDDED:-2500MB}" ;;
        isSimEmbeddedInclusive) RJ_REQUEST_MEMORY="${RJ_SMOKE_REQUEST_MEMORY_EMBEDDED_INCLUSIVE:-${RJ_SMOKE_REQUEST_MEMORY_EMBEDDED:-2500MB}}" ;;
        *)                      RJ_REQUEST_MEMORY="${RJ_SMOKE_REQUEST_MEMORY_SIM:-1800MB}" ;;
      esac
      export RJ_REQUEST_MEMORY
    fi
    smoke_stamp="$(date +%Y%m%d_%H%M%S)"
    sim_smoke_out="${RJ_SMOKE_OUTPUT_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaSmoke/${TAG}_smokeTest_${smoke_stamp}}"
    export RJ_DEST_BASE_OVERRIDE="$sim_smoke_out"
    export RJ_MERGE_OUT_BASE_OVERRIDE="${RJ_MERGE_OUT_BASE_OVERRIDE:-${BASE}/outputSmoke/${TAG}_smokeTest_${smoke_stamp}}"
    export RJ_PROFILE_JOB=1
    export RJ_DIRECT_NEVENTS="${RJ_SMOKE_SIM_NEVENTS:-3000}"
    export RJ_JOB_HEARTBEAT_SECONDS="${RJ_JOB_HEARTBEAT_SECONDS:-${RJ_SMOKE_JOB_HEARTBEAT_SECONDS:-120}}"
    export RJ_PROFILE_STAGE="${RJ_PROFILE_STAGE:-directSmoke}"
    export RJ_PROFILE_LABEL="${RJ_PROFILE_LABEL:-${TAG}_smokeTest}"
    if [[ -z "${RJ_VALIDATE_SIM_INPUT_MAX_LINES:-}" && "$GROUP_SIZE" =~ ^[0-9]+$ && "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 ]]; then
      export RJ_VALIDATE_SIM_INPUT_MAX_LINES="$(( GROUP_SIZE * MAX_JOBS ))"
    fi

    say "${BOLD}SIM direct-fanout smokeTest requested${RST}"
    say "  dataset      : ${DATASET}"
    say "  output base  : ${RJ_DEST_BASE_OVERRIDE}"
    say "  merge output : ${RJ_MERGE_OUT_BASE_OVERRIDE}"
    say "  groupSize    : ${GROUP_SIZE}"
    say "  maxJobs/cfg/sample: ${MAX_JOBS}"
    say "  nEvents/job  : ${RJ_DIRECT_NEVENTS} (0 means full worker input)"
    say "  request mem  : ${RJ_REQUEST_MEMORY}"
    say "  engine       : direct RecoilJets fanout; pool replay is not used"
    if [[ "${RJ_DAG_DRYRUN:-0}" == "1" || "${RJ_DAG_DRYRUN:-0}" == "true" || "${RJ_DAG_DRYRUN:-0}" == "TRUE" ]]; then
      _smoke_check_args=( "$DATASET" CHECKJOBS groupSize "$GROUP_SIZE" )
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]]; then
        _smoke_check_args+=( "SAMPLE=${SIM_SAMPLE}" )
      fi
      "$0" "${_smoke_check_args[@]}"
      echo "RECOILJETS_SMOKETEST_DRYRUN_V1"
      echo "dataset=${DATASET}"
      echo "mode=condorDoAllSmoke"
      echo "engine=direct_legacy_fanout"
      echo "legacy_output_parity=NOT_RUN_DIRECT_RECOILJETS_ENGINE"
      echo "output_base=${RJ_DEST_BASE_OVERRIDE}"
      echo "merge_output_base=${RJ_MERGE_OUT_BASE_OVERRIDE}"
      echo "continuing_to_build_full_auto_dag=1"
    fi
    _smoke_args=( "$DATASET" condorDoAll groupSize "$GROUP_SIZE" maxJobs "$MAX_JOBS" )
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]]; then
      _smoke_args+=( "SAMPLE=${SIM_SAMPLE}" )
    fi
    "$0" "${_smoke_args[@]}"
    exit $?
    ;;

  condorDoAllFromScratch)
    [[ "$IS_SIM" -eq 1 ]] || { err "condorDoAllFromScratch is only valid for isSim variants"; exit 2; }
    warn "condorDoAllFromScratch now uses direct RecoilJets fanout; pool capture/replay was removed."
    _alias_args=( "$DATASET" condorDoAll groupSize "$GROUP_SIZE" )
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]]; then
      _alias_args+=( "SAMPLE=${SIM_SAMPLE}" )
    fi
    if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 ]]; then
      _alias_args+=( maxJobs "$MAX_JOBS" )
    fi
    "$0" "${_alias_args[@]}"
    ;;

  condorHistFromPool)
    err "condorHistFromPool was removed; use direct RecoilJets fanout: $0 $DATASET condorDoAll"
    exit 2
    ;;

  condorDoAllDirect)
    [[ "$IS_SIM" -eq 1 ]] || { err "condorDoAllDirect is only valid for isSim variants"; exit 2; }
    warn "condorDoAllDirect disables ID fanout and submits one legacy cfg-tag per DST pass for validation."
    _direct_args=( "$DATASET" condorDoAll groupSize "$GROUP_SIZE" )
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]]; then
      _direct_args+=( "SAMPLE=${SIM_SAMPLE}" )
    fi
    if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 ]]; then
      _direct_args+=( maxJobs "$MAX_JOBS" )
    fi
    RJ_DISABLE_ID_FANOUT=1 RJ_DISABLE_JET_PT_INTERNALIZATION=1 RJ_DISABLE_DPHI_INTERNALIZATION=1 RJ_DIRECT_DST_DOALL=1 "$0" "${_direct_args[@]}"
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
    if [[ "$DATASET" == "isSimEmbedded" || "$DATASET" == "isSimEmbeddedInclusive" ]]; then
      if [[ "$(basename "$master_yaml")" == *"auau_bdt_validation"* || "$master_yaml" == *"auau_bdt_validation"* || "$(basename "$master_yaml")" == *"auau_bdt_etfine"* || "$master_yaml" == *"auau_bdt_etfine"* ]]; then
        if [[ -z "${RJ_ID_FANOUT_MAX_ROWS:-}" ]]; then
          # The AuAu BDT validation matrix loads many TMVA readers. The
          # validated WP0.50 campaign used one photon-ID row per worker; keep
          # that execution shape as the default for WP scans so memory remains
          # predictable. Callers may still override for stress tests.
          export RJ_ID_FANOUT_MAX_ROWS=1
          say "AuAu BDT validation config detected: defaulting RJ_ID_FANOUT_MAX_ROWS=1 (set explicitly to override)"
        fi
      fi
    fi
    direct_nevents="$(direct_worker_nevents)"
    direct_request_memory="$(memory_request_from_env_or_default "2000MB")"
    direct_memory_floor_mb=0
    if [[ -z "${RJ_REQUEST_MEMORY:-}" && -z "${RJ_REQUEST_MEMORY_MB:-}" ]]; then
      case "$DATASET" in
        isSimEmbedded|isSimEmbeddedInclusive)
          # Embedded AuAu workers are ROOT/TMVA heavy enough that the generic
          # SIM default can trigger cgroup aborts before auto-retry has a
          # chance to help. Keep this as the safe production floor; callers can
          # still override explicitly with RJ_REQUEST_MEMORY.
          direct_memory_floor_mb=4000
          if [[ "$(basename "$master_yaml")" == *"auau_bdt_validation"* || "$master_yaml" == *"auau_bdt_validation"* || "$(basename "$master_yaml")" == *"auau_bdt_etfine"* || "$master_yaml" == *"auau_bdt_etfine"* ]]; then
            # The expanded AuAu BDT validation matrix fans out many TMVA model
            # readers in one process. SDCC evidence from WP0.80 validation
            # showed cgroup holds at 4.6 GB and then again right at the 9 GB
            # allocation boundary during ROOT/model startup. Use a higher
            # campaign-only floor with enough headroom for SDCC cgroup rounding.
            direct_memory_floor_mb=12000
            export RJ_AUTO_MEMORY_RETRY_CAP_MB="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-16000}"
          fi
          direct_request_memory="${direct_memory_floor_mb}MB"
          ;;
      esac
    else
      case "$DATASET" in
        isSimEmbedded|isSimEmbeddedInclusive)
          direct_memory_floor_mb=4000
          if [[ "$(basename "$master_yaml")" == *"auau_bdt_validation"* || "$master_yaml" == *"auau_bdt_validation"* || "$(basename "$master_yaml")" == *"auau_bdt_etfine"* || "$master_yaml" == *"auau_bdt_etfine"* ]]; then
            direct_memory_floor_mb=12000
            export RJ_AUTO_MEMORY_RETRY_CAP_MB="${RJ_AUTO_MEMORY_RETRY_CAP_MB:-16000}"
          fi
          ;;
      esac
    fi
    direct_request_memory_mb="$(memory_request_to_mb "$direct_request_memory")"
    if [[ "$DATASET" == "isSimEmbedded" || "$DATASET" == "isSimEmbeddedInclusive" ]]; then
      if (( direct_memory_floor_mb > 0 && direct_request_memory_mb > 0 && direct_request_memory_mb < direct_memory_floor_mb )); then
        warn "${DATASET} condorDoAll worker request_memory=${direct_request_memory} is below the validated floor of ${direct_memory_floor_mb}MB for this config."
        warn "Set RJ_REQUEST_MEMORY=${direct_memory_floor_mb}MB or higher unless this is an intentional memory stress test."
      fi
    fi

    mapfile -t sim_pts   < <( yaml_get_values "jet_pt_min" "$master_yaml" )
    mapfile -t sim_fracs < <( yaml_get_values "back_to_back_dphi_min_pi_fraction" "$master_yaml" )
    mapfile -t sim_vzs   < <( dataset_vz_values "$master_yaml" )
    mapfile -t sim_cones < <( yaml_get_values "coneR" "$master_yaml" )
    (( ${#sim_pts[@]} ))   || { err "No values found for jet_pt_min in $master_yaml"; exit 72; }
    (( ${#sim_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $master_yaml"; exit 72; }
    (( ${#sim_vzs[@]} ))   || { err "No values found for vz_cut_cm in $master_yaml"; exit 72; }
    (( ${#sim_cones[@]} )) || { err "No values found for coneR in $master_yaml"; exit 72; }
    sim_view_cones=( "${sim_cones[@]}" )
    mapfile -t sim_cones < <( cone_submit_values "${sim_view_cones[@]}" )
    mapfile -t sim_submit_pts < <( jetpt_submit_values "${sim_pts[@]}" )
    mapfile -t sim_submit_fracs < <( dphi_submit_values "${sim_fracs[@]}" )
    jetpt_env_for_sub="$(jetpt_env_fragment "${sim_pts[@]}")"
    dphi_env_for_sub="$(dphi_env_fragment "${sim_fracs[@]}")"
    iso_view_env_for_sub="$(iso_view_env_fragment)"
    stitch_env_for_sub="$(embedded_inclusive_stitch_env_fragment)"
    build_iso_modes "$master_yaml"
    read_uepipe_modes "$master_yaml" "$TAG"

    samples=()
    if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
      case "$DATASET" in
        isSimEmbedded)          samples=( "run28_embeddedPhoton12" "run28_embeddedPhoton20" ) ;;
        isSimEmbeddedInclusive) mapfile -t samples < <(simembeddedinclusive_sample_list) ;;
        isSimInclusive|isSimJet5) samples=( "run28_jet5" "run28_jet8" "run28_jet12" "run28_jet20" "run28_jet30" "run28_jet40" ) ;;
        isSimMB)                samples=( "run28_detroit" ) ;;
        *)
          if [[ "${RJ_SIM_SIGNAL_SAMPLE_SET:-}" == "ppg12_double" || "${RJ_SIM_SIGNAL_SAMPLE_SET:-}" == "double" ]]; then
            samples=( "run28_photonjet5_double" "run28_photonjet10_double" "run28_photonjet20_double" )
          else
            samples=( "run28_photonjet5" "run28_photonjet10" "run28_photonjet20" )
          fi
          ;;
      esac
    else
      samples=( "${SIM_SAMPLE}" )
    fi

    # If CHECKJOBS was also provided, do a dry-run count for the FULL condorDoAll matrix and exit.
    if [[ "${DRYRUN:-0}" -eq 1 ]]; then
      GROUP_SIZE="$gs_doall"
      check_jobs_sim
      exit 0
    fi

    SIM_DEST_BASE_RESOLVED="$DEST_BASE"
    cleanup_dataset_outputs_if_requested

    if auto_merge_enabled; then
      need_cmd condor_submit_dag
    else
      need_cmd condor_submit
    fi
    doall_stamp="$(date +%Y%m%d_%H%M%S)"
    if [[ -n "${RJ_SUBMISSION_NAMESPACE:-}" ]]; then
      doall_stamp="${doall_stamp}_$(printf '%s' "$RJ_SUBMISSION_NAMESPACE" | tr -c 'A-Za-z0-9_.-' '_')"
    fi
    SIM_YAML_OVERRIDE_DIR="${SIM_YAML_OVERRIDE_DIR}/${TAG}_condorDoAll_${doall_stamp}"
    mkdir -p "$SIM_YAML_OVERRIDE_DIR"
    say "SIM YAML/fanout artifact dir: ${SIM_YAML_OVERRIDE_DIR}"
    say "SIM worker memory request: ${direct_request_memory} (${direct_request_memory_mb} MB)"
    if [[ "$DATASET" == "isSimEmbedded" || "$DATASET" == "isSimEmbeddedInclusive" ]]; then
      say "SIM merge memory hint   : ${RJ_SIM_FIRSTROUND_REQUEST_MEMORY:-auto/default} (firstRound; independent of worker memory)"
    fi
    say "SIM production vz selection:"
    say_vz_selection_summary "$master_yaml" "${sim_vzs[@]}"
    # Freeze pipeline for this bulk submission
    cleanup_bulk_snapshots_for_tag
    case "$DATASET" in
      isSimEmbedded|isSimEmbeddedInclusive) create_pipeline_snapshot "auau" "$doall_stamp" ;;
      *)                                    create_pipeline_snapshot "pp" "$doall_stamp" ;;
    esac
    # Clean stale .sub files only. Keep YAML overrides/snapshots because live Condor jobs may still reference them.
    if ! auto_merge_enabled; then
      rm -f "${SUB_DIR}/RecoilJets_sim_"*.sub "${SUB_DIR}/RecoilJets_${TAG}_"*.sub 2>/dev/null || true
    fi
    SIM_STAGE_NAMESPACE="${TAG}_${doall_stamp}"
    export SIM_STAGE_NAMESPACE
    say "SIM chunk-list stage namespace: ${SIM_STAGE_NAMESPACE}"

    auto_workflow=0
    auto_dag=""
    auto_dag_dir=""
    if auto_merge_enabled; then
      auto_workflow=1
      auto_dag_dir="${SUB_DIR}/auto_workflow_${TAG}_${doall_stamp}"
      mkdir -p "$auto_dag_dir"
      auto_dag="${auto_dag_dir}/RecoilJets_auto_${TAG}_${doall_stamp}.dag"
      : > "$auto_dag"
      write_cleanup_manifest "${auto_dag_dir}/cleanup_manifest.txt" "$SIM_YAML_OVERRIDE_DIR" "$SIM_DEST_BASE_RESOLVED"
      RJ_DAG_COLLECTED_NODES=()
      export RJ_COLLECT_DAG_FILE="$auto_dag"
      export RJ_COLLECT_SUB_DIR="${auto_dag_dir}/analysis_submit_files"
      export RJ_COLLECT_NODE_PREFIX="ANALYSIS_${TAG}"
      say "${BOLD}Automatic analysis-to-merge DAG enabled${RST}"
      say "  dag          : ${auto_dag}"
      say "  escape hatch : RJ_AUTO_MERGE=0 keeps analysis-only submission"
    fi

    for pt in "${sim_submit_pts[@]}"; do
      for frac in "${sim_submit_fracs[@]}"; do
        for vz in "${sim_vzs[@]}"; do
        for cone in "${sim_cones[@]}"; do
        for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
        direct_shard_idx=0
        direct_shard_total=0
        if iso_cone_fanout_enabled; then
          iso_idx_is_group_leader "$iso_idx" || continue
          direct_shard_idx="$(direct_fanout_cell_ordinal "$cone" "$iso_idx" "${sim_cones[@]}")"
          direct_shard_total="$(direct_fanout_shard_count "${#sim_cones[@]}" "${#iso_tags[@]}")"
          (( direct_shard_idx <= direct_shard_total )) || continue
        else
          id_fanout_enabled && ! iso_idx_is_group_leader "$iso_idx" && continue
        fi
        for uepipe in "${uepipe_modes[@]}"; do
        if iso_cone_fanout_enabled; then
          SIM_CFG_TAG="$(direct_fanout_shard_tag "$pt" "$frac" "$vz" "$uepipe" "$direct_shard_idx")"
        else
          SIM_CFG_TAG="$(matrix_cfg_tag "$pt" "$frac" "$vz" "$cone" "$iso_idx" "$uepipe")"
        fi
        DEST_BASE="${SIM_DEST_BASE_RESOLVED}/${SIM_CFG_TAG}"
        fanout_dirs="${SIM_YAML_OVERRIDE_DIR}/id_fanout_${SIM_CFG_TAG}_${doall_stamp}.txt"
        if iso_cone_fanout_enabled; then
          emit_direct_fanout_shard_file "$fanout_dirs" "$SIM_DEST_BASE_RESOLVED" "$pt" "$frac" "$vz" "$uepipe" "$direct_shard_idx" "${sim_cones[@]}"
        elif id_fanout_enabled; then
          emit_id_fanout_dirs_file "$fanout_dirs" "$SIM_DEST_BASE_RESOLVED" "$pt" "$frac" "$vz" "$cone" "$iso_idx" "$uepipe"
        else
          : > "$fanout_dirs"
        fi
        yaml_override="$(sim_make_yaml_override "$master_yaml" "$pt" "$frac" "$vz" "$cone" "${iso_sliding[$iso_idx]}" "${iso_fixed[$iso_idx]}" "${uepipe}" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}" "$SIM_CFG_TAG" "$doall_stamp")"

        for samp in "${samples[@]}"; do
          SIM_SAMPLE="$samp"
          GROUP_SIZE="$gs_doall"
          sim_init

          if id_fanout_enabled; then
            while IFS='|' read -r fan_dest _fan_cfg _fan_pre _fan_tight _fan_nonTight; do
              [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
              SIM_OUT_DIR="${fan_dest}/${SIM_SAMPLE}"
              if [[ "${RJ_CLEAN_OUTPUT_BASE:-0}" == "1" ]] && dag_dryrun_enabled; then
                say "DRYRUN: would clean SIM_OUT_DIR=${SIM_OUT_DIR}"
              elif [[ "${RJ_CLEAN_OUTPUT_BASE:-0}" == "1" ]]; then
                mkdir -p "$SIM_OUT_DIR"
                rm -f "${SIM_OUT_DIR}/"*.root 2>/dev/null || true
                find "${SIM_OUT_DIR}" -maxdepth 1 -name "*_LOCAL_*.root" -delete 2>/dev/null || true
              else
                mkdir -p "$SIM_OUT_DIR"
              fi
            done < "$fanout_dirs"
          else
            SIM_OUT_DIR="${DEST_BASE}/${SIM_SAMPLE}"
            if [[ "${RJ_CLEAN_OUTPUT_BASE:-0}" == "1" ]] && dag_dryrun_enabled; then
              say "DRYRUN: would clean SIM_OUT_DIR=${SIM_OUT_DIR}"
            elif [[ "${RJ_CLEAN_OUTPUT_BASE:-0}" == "1" ]]; then
              mkdir -p "$SIM_OUT_DIR"
              rm -f "${SIM_OUT_DIR}/"*.root 2>/dev/null || true
              find "${SIM_OUT_DIR}" -maxdepth 1 -name "*_LOCAL_*.root" -delete 2>/dev/null || true
            else
              mkdir -p "$SIM_OUT_DIR"
            fi
          fi
          rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_grp"*.list 2>/dev/null || true
          rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_LOCAL_"*.list 2>/dev/null || true
          rm -f "${SIM_STAGE_DIR}/${SIM_JOB_PREFIX}_condorTest_"*.list 2>/dev/null || true

          _rj_group_limit="$MAX_JOBS"
          if [[ -n "${RJ_THE134_EXTRACTION_EXECUTION_MANIFEST:-}" ]]; then
            if [[ ! "${THE134_ADMITTED_EXPECTED_JOB_COUNT:-}" =~ ^[1-9][0-9]*$ ]]; then
              err "THE-134 admitted expected job count is missing before group materialization"
              exit 30
            fi
            _rj_group_limit="$THE134_ADMITTED_EXPECTED_JOB_COUNT"
          fi
          _rj_group_output=""
          if ! _rj_group_output="$(make_sim_groups "$GROUP_SIZE" "$_rj_group_limit")"; then
            err "SIM group materialization failed before Condor submission (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"
            exit 30
          fi
          [[ -n "$_rj_group_output" ]] || {
            err "No sim groups produced (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"
            exit 30
          }
          mapfile -t groups <<< "$_rj_group_output"
          unset _rj_group_output
          (( ${#groups[@]} )) || { err "No sim groups produced (sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG})"; exit 30; }
          if [[ -n "${RJ_THE134_EXTRACTION_EXECUTION_MANIFEST:-}" &&
                "${#groups[@]}" -ne "$THE134_ADMITTED_EXPECTED_JOB_COUNT" ]]; then
            err "THE-134 group count differs before Condor submission: expected=${THE134_ADMITTED_EXPECTED_JOB_COUNT} observed=${#groups[@]}"
            exit 30
          fi
          if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 && "${#groups[@]}" -gt "$MAX_JOBS" ]]; then
            say "Capping ${DATASET} group list for sample=${SIM_SAMPLE}, tag=${SIM_CFG_TAG}: ${#groups[@]} → ${MAX_JOBS} jobs"
            groups=( "${groups[@]:0:$MAX_JOBS}" )
          fi
          unset _rj_group_limit

          stamp="$(date +%Y%m%d_%H%M%S)"
          sub="${SUB_DIR}/RecoilJets_sim_${SIM_CFG_TAG}_${SIM_SAMPLE}_${stamp}.sub"
          args_file="${SUB_DIR}/RecoilJets_sim_${SIM_CFG_TAG}_${SIM_SAMPLE}_${stamp}.args"
          : > "$args_file"
          exe_for_sub="${BULK_FROZEN_EXE:-${EXE}}"
          macro_env_for_sub=""
          [[ -n "${BULK_FROZEN_MACRO:-}" ]] && macro_env_for_sub=";RJ_MACRO_PATH=${BULK_FROZEN_MACRO}"
          fanout_env_for_sub=""
          id_fanout_enabled && fanout_env_for_sub=";RJ_ID_FANOUT_FILE=${fanout_dirs};RJ_ID_FANOUT_DIRS_FILE=${fanout_dirs}"
          submit_extra_env_for_sub="$(build_submit_extra_env_fragment)"
          direct_require_non_tiny_output="$(effective_require_non_tiny_output)"

          cat > "$sub" <<SUB
universe      = vanilla
executable    = ${exe_for_sub}
initialdir    = ${BASE}
getenv        = $(condor_getenv_directive)
log           = ${LOG_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/${SIM_JOB_PREFIX}.job.\$(Cluster).\$(Process).err
$(condor_auto_memory_retry_block "$direct_request_memory_mb")
$(condor_worker_failure_hold_block)
$(condor_late_materialization_block)
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
environment   = RJ_VERBOSITY=0;RJ_CONFIG_YAML=${yaml_override}${macro_env_for_sub}${fanout_env_for_sub}${jetpt_env_for_sub}${dphi_env_for_sub}${iso_view_env_for_sub}${stitch_env_for_sub}${submit_extra_env_for_sub};RJ_SIM_SAMPLE=${SIM_SAMPLE};RJ_EMBEDDED_INCLUSIVE_JET_SAMPLE=${SIM_SAMPLE};RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES:-0};RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30:-0};RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES:-0};RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40:-0};RJ_PROFILE_JOB=${RJ_PROFILE_JOB:-0};RJ_JOB_HEARTBEAT_SECONDS=${RJ_JOB_HEARTBEAT_SECONDS:-0};RJ_PROFILE_STAGE=${RJ_PROFILE_STAGE:-direct};RJ_PROFILE_LABEL=${RJ_PROFILE_LABEL:-${TAG}};RJ_REQUEST_MEMORY_MB=${direct_request_memory_mb};RJ_REQUIRE_NON_TINY_OUTPUT=${direct_require_non_tiny_output};RJ_MIN_OUTPUT_BYTES=${RJ_MIN_OUTPUT_BYTES:-50000};RJ_FAIL_ON_MISSING_CALO_INPUT=${RJ_FAIL_ON_MISSING_CALO_INPUT:-0}
queue arguments from ${args_file}
SUB

          gidx=0
          for glist in "${groups[@]}"; do
            (( gidx+=1 ))
            printf '%s %s %s $(Cluster) %d %d NONE %s\n' \
                   "$SIM_SAMPLE" "$glist" "$DATASET" "$direct_nevents" "$gidx" "$DEST_BASE" >> "$args_file"
          done

          if iso_cone_fanout_enabled; then
            say "Submitting ${DATASET} condorDoAll iso/cone/ID fanout (shard=${direct_shard_idx}/${direct_shard_total}, sample=${SIM_SAMPLE}, groupSize=${GROUP_SIZE}, nEvents=${direct_nevents}, uepipe=${uepipe}) → jobs=${BOLD}${#groups[@]}${RST}"
            say "Fanout outputs: $(wc -l < "$fanout_dirs" | awk '{print $1}') cfg tags from one DST pass"
            say "Fanout dirs file: ${fanout_dirs}"
          elif id_fanout_enabled; then
            say "Submitting ${DATASET} condorDoAll RecoilJets fanout (base tag=${SIM_CFG_TAG}, sample=${SIM_SAMPLE}, groupSize=${GROUP_SIZE}, nEvents=${direct_nevents}, uepipe=${uepipe}) → jobs=${BOLD}${#groups[@]}${RST}"
            say "Fanout outputs: $(wc -l < "$fanout_dirs" | awk '{print $1}') cfg tags from one DST pass"
            say "Fanout dirs file: ${fanout_dirs}"
          else
            say "Submitting ${DATASET} condorDoAllDirect legacy cfg (tag=${SIM_CFG_TAG}, sample=${SIM_SAMPLE}, groupSize=${GROUP_SIZE}, nEvents=${direct_nevents}, uepipe=${uepipe}) → jobs=${BOLD}${#groups[@]}${RST}"
          fi
          say "YAML override: ${yaml_override}"
          submit_or_collect_condor "$sub" "analysis_${TAG}_${SIM_SAMPLE}"
        done
        done
        done
        done
        done
      done
    done
    if (( auto_workflow )); then
      unset RJ_COLLECT_DAG_FILE
      unset RJ_COLLECT_SUB_DIR
      unset RJ_COLLECT_NODE_PREFIX
      runner="${auto_dag_dir}/auto_stage_runner.sh"
      final_notify="${auto_dag_dir}/auto_final_notify.sh"
      write_auto_stage_runner "$runner"
      write_auto_final_notify "$final_notify"
      sim_merge_out_base="${RJ_MERGE_OUT_BASE_OVERRIDE:-${BASE}/output}"
      sim_merge_group_size_default=300
      case "$DATASET" in
        isSimEmbedded|isSimEmbeddedInclusive)
          sim_merge_group_size_default=75
          ;;
      esac
      sim_merge_group_size="${RJ_SIM_MERGE_GROUP_SIZE:-${sim_merge_group_size_default}}"
      first_round_args=( env "MERGE_CONFIG_YAML=${master_yaml}" "MERGE_SIM_INPUT_BASE_OVERRIDE=${SIM_DEST_BASE_RESOLVED}" "MERGE_OUT_BASE_OVERRIDE=${sim_merge_out_base}" "MERGE_CFG_MATCH=${RJ_PHOTON_ID_ROW_MATCH:-}" "RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40:-0}" "RJ_STAGE_EMAIL_MODE=none" "RJ_STAGE_EMAIL_STRICT=1" )
      while IFS= read -r stitch_kv; do
        [[ -n "$stitch_kv" ]] && first_round_args+=( "$stitch_kv" )
      done < <(embedded_inclusive_stitch_env_args)
      first_round_args+=( "${BASE}/scripts/mergeRecoilJets.sh" "$DATASET" firstRound groupSize "${sim_merge_group_size}" )
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]]; then
        first_round_args+=( "SAMPLE=${SIM_SAMPLE}" )
      fi
      add_auto_stage_node "$auto_dag" "SIM_FIRSTROUND" "$runner" "sim_firstRound_${TAG}_all" "${first_round_args[@]}"
      second_round_args=( env "MERGE_CONFIG_YAML=${master_yaml}" "MERGE_SIM_INPUT_BASE_OVERRIDE=${SIM_DEST_BASE_RESOLVED}" "MERGE_OUT_BASE_OVERRIDE=${sim_merge_out_base}" "MERGE_CFG_MATCH=${RJ_PHOTON_ID_ROW_MATCH:-}" "RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40:-0}" "RJ_STAGE_EMAIL_MODE=none" "RJ_STAGE_EMAIL_STRICT=1" )
      while IFS= read -r stitch_kv; do
        [[ -n "$stitch_kv" ]] && second_round_args+=( "$stitch_kv" )
      done < <(embedded_inclusive_stitch_env_args)
      second_round_args+=( "${BASE}/scripts/mergeRecoilJets.sh" "$DATASET" secondRound condor )
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 1 ]]; then
        second_round_args+=( "SAMPLE=${SIM_SAMPLE}" )
      fi
      add_auto_stage_node "$auto_dag" "SIM_SECONDROUND" "$runner" "sim_secondRound_${TAG}_all" "${second_round_args[@]}"
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
        final_stitch_args=( env "MERGE_CONFIG_YAML=${master_yaml}" "MERGE_SIM_INPUT_BASE_OVERRIDE=${SIM_DEST_BASE_RESOLVED}" "MERGE_OUT_BASE_OVERRIDE=${sim_merge_out_base}" "MERGE_CFG_MATCH=${RJ_PHOTON_ID_ROW_MATCH:-}" "RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_THREE_SAMPLES:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET30:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES=${RJ_SIMEMBEDDEDINCLUSIVE_FOUR_SAMPLES:-0}" "RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40=${RJ_SIMEMBEDDEDINCLUSIVE_INCLUDE_JET40:-0}" "RJ_STAGE_EMAIL_MODE=none" "RJ_STAGE_EMAIL_STRICT=1" )
        while IFS= read -r stitch_kv; do
          [[ -n "$stitch_kv" ]] && final_stitch_args+=( "$stitch_kv" )
        done < <(embedded_inclusive_stitch_env_args)
        final_stitch_args+=( "${BASE}/scripts/mergeRecoilJets.sh" "$DATASET" finalStitch condor )
        add_auto_stage_node "$auto_dag" "SIM_FINALSTITCH" "$runner" "sim_finalStitch_${TAG}_all" "${final_stitch_args[@]}"
      fi
      if (( ${#RJ_DAG_COLLECTED_NODES[@]} > 0 )); then
        printf 'PARENT' >> "$auto_dag"
        printf ' %s' "${RJ_DAG_COLLECTED_NODES[@]}" >> "$auto_dag"
        printf ' CHILD SIM_FIRSTROUND\n' >> "$auto_dag"
      fi
      printf 'PARENT SIM_FIRSTROUND CHILD SIM_SECONDROUND\n' >> "$auto_dag"
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
        printf 'PARENT SIM_SECONDROUND CHILD SIM_FINALSTITCH\n' >> "$auto_dag"
      fi
      add_auto_final_node "$auto_dag" "FINAL_NOTIFY" "$final_notify" "auto_${TAG}_final_ready" "$DATASET" "./scripts/sftp_get_recoiljets_outputs.sh ${DATASET}" "${sim_merge_out_base}/${TAG}"
      validate_auto_dag_submit_files "$auto_dag"
      say "Automatic workflow DAG built:"
      say "  analysis nodes : ${#RJ_DAG_COLLECTED_NODES[@]}"
      if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
        say "  merge stages   : SIM_FIRSTROUND -> SIM_SECONDROUND -> SIM_FINALSTITCH (quiet strict validation; stage-boundary and final emails)"
      else
        say "  merge stages   : SIM_FIRSTROUND -> SIM_SECONDROUND (sample-explicit run; final stitching skipped)"
      fi
      say "  dag            : ${auto_dag}"
      if dag_dryrun_enabled; then
        echo "RECOILJETS_AUTO_DAG_DRYRUN_V1"
        echo "dataset=${DATASET}"
        echo "dag=${auto_dag}"
        echo "analysis_nodes=${#RJ_DAG_COLLECTED_NODES[@]}"
        if [[ "${SIM_SAMPLE_EXPLICIT:-0}" -eq 0 ]]; then
          echo "sim_final_stitch=1"
        else
          echo "sim_final_stitch=0"
        fi
        echo "next_action_after_ready=./scripts/sftp_get_recoiljets_outputs.sh ${DATASET}"
        sed -n '1,240p' "$auto_dag"
      else
        condor_submit_dag -notification Never "$auto_dag"
      fi
    fi
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
    mapfile -t data_vzs   < <( dataset_vz_values "$data_yaml_src" )
    mapfile -t data_cones < <( yaml_get_values "coneR" "$data_yaml_src" )
    (( ${#data_pts[@]} ))   || { err "No values found for jet_pt_min in $data_yaml_src"; exit 72; }
    (( ${#data_fracs[@]} )) || { err "No values found for back_to_back_dphi_min_pi_fraction in $data_yaml_src"; exit 72; }
    (( ${#data_vzs[@]} ))   || { err "No values found for vz_cut_cm in $data_yaml_src"; exit 72; }
    (( ${#data_cones[@]} )) || { err "No values found for coneR in $data_yaml_src"; exit 72; }
    data_view_cones=( "${data_cones[@]}" )
    mapfile -t data_cones < <( cone_submit_values "${data_view_cones[@]}" )
    mapfile -t data_submit_pts < <( jetpt_submit_values "${data_pts[@]}" )
    mapfile -t data_submit_fracs < <( dphi_submit_values "${data_fracs[@]}" )
    jetpt_env_for_sub="$(jetpt_env_fragment "${data_pts[@]}")"
    dphi_env_for_sub="$(dphi_env_fragment "${data_fracs[@]}")"
    iso_view_env_for_sub="$(iso_view_env_fragment)"
    build_iso_modes "$data_yaml_src"
    read_uepipe_modes "$data_yaml_src" "$TAG"
    DATA_DEST_BASE_SAVED="$DEST_BASE"
    if [[ "${3:-}" == "all" || -z "${3:-}" ]]; then
      cleanup_dataset_outputs_if_requested
    fi
    say "DATA production vz selection:"
    say_vz_selection_summary "$data_yaml_src" "${data_vzs[@]}"

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
        pt0="${data_submit_pts[0]}"
        frac0="${data_submit_fracs[0]}"
        vz0="${data_vzs[0]}"
        cone0="${data_cones[0]}"
        data_cfg_tag="$(matrix_cfg_tag "$pt0" "$frac0" "$vz0" "$cone0" 0 "${uepipe_modes[0]}")"

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
        force_ppg12_photon_yield_yaml_contract "$yaml_snap"
        pin_photon_id_scalars_in_yaml "$yaml_snap" "${iso_preselection[0]}" "${iso_tight[0]}" "${iso_nonTight[0]}"
        say "YAML snapshot (pt=${pt0}, frac=${frac0}, vz=${vz0}, coneR=${cone0}, iso=${iso_tags[0]}, uepipe=${uepipe_modes[0]}): ${yaml_snap}"
        DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"
        submit_extra_env_for_sub="$(build_submit_extra_env_fragment)"
        test_require_non_tiny_output="$(effective_require_non_tiny_output)"
        test_nevents="${RJ_TEST_DATA_NEVENTS:-${RJ_DIRECT_NEVENTS:-${RJ_SMOKE_DATA_NEVENTS:-3000}}}"
        [[ "$test_nevents" == "-1" ]] && test_nevents=0
        [[ "$test_nevents" =~ ^[0-9]+$ ]] || { err "RJ_TEST_DATA_NEVENTS/RJ_DIRECT_NEVENTS must be -1 or a non-negative integer, got '${test_nevents}'"; exit 2; }
        test_verbosity="${RJ_TEST_VERBOSITY:-${RJ_CONDOR_TEST_VERBOSITY:-0}}"
        [[ "$test_verbosity" =~ ^[0-9]+$ ]] || { err "RJ_TEST_VERBOSITY must be a non-negative integer, got '${test_verbosity}'"; exit 2; }
        test_request_memory="$(memory_request_from_env_or_default "2000MB")"
        test_request_memory_mb="$(memory_request_to_mb "$test_request_memory")"

        cat > "$sub" <<SUB
universe      = vanilla
executable    = ${EXE}
initialdir    = ${BASE}
getenv        = True
log           = ${LOG_DIR}/job.\$(Cluster).\$(Process).log
output        = ${OUT_DIR}/job.\$(Cluster).\$(Process).out
error         = ${ERR_DIR}/job.\$(Cluster).\$(Process).err
$(condor_auto_memory_retry_block "$test_request_memory_mb")
$(condor_worker_failure_hold_block)
should_transfer_files = NO
stream_output = True
stream_error  = True
notification  = Never
environment   = RJ_DATASET=${DATASET};RJ_VERBOSITY=${test_verbosity};RJ_CONFIG_YAML=${yaml_snap}${jetpt_env_for_sub}${dphi_env_for_sub}${iso_view_env_for_sub}${submit_extra_env_for_sub};RJ_PROFILE_JOB=${RJ_PROFILE_JOB:-0};RJ_JOB_HEARTBEAT_SECONDS=${RJ_JOB_HEARTBEAT_SECONDS:-0};RJ_PROFILE_STAGE=${RJ_PROFILE_STAGE:-testJob};RJ_PROFILE_LABEL=${RJ_PROFILE_LABEL:-${TAG}_testJob};RJ_REQUEST_MEMORY_MB=${test_request_memory_mb};RJ_REQUIRE_NON_TINY_OUTPUT=${test_require_non_tiny_output};RJ_MIN_OUTPUT_BYTES=${RJ_MIN_OUTPUT_BYTES:-50000};RJ_FAIL_ON_MISSING_CALO_INPUT=${RJ_FAIL_ON_MISSING_CALO_INPUT:-0}
arguments     = ${r8} ${glist} ${DATASET} \$(Cluster) ${test_nevents} 1 NONE ${DEST_BASE}
queue
SUB
        say "Submitting 1 test job on run ${BOLD}${r8}${RST} (first chunk, groupSize=1, nEvents=${test_nevents}, vz=${vz0}, RJ_VERBOSITY=${test_verbosity}, requestMem=${test_request_memory_mb}MB) → $(basename "$sub")"
        condor_submit "$sub"
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
        smoke_stamp="$(date +%Y%m%d_%H%M%S)"
        smoke_run_count="${RJ_SMOKE_DATA_RUNS:-10}"
        smoke_selected_runs="${SUB_DIR}/${TAG}_directSmoke_${smoke_stamp}_runs.list"
        smoke_selected_stats="${SUB_DIR}/${TAG}_directSmoke_${smoke_stamp}_run_stats.txt"
        if [[ -n "${RJ_SMOKE_DATA_RUN:-}" ]]; then
          select_explicit_stat_data_run "$RJ_SMOKE_DATA_RUN" "$smoke_selected_runs" "$smoke_selected_stats" || { err "smokeTest could not select explicit DATA run ${RJ_SMOKE_DATA_RUN} from ${GOLDEN}"; exit 99; }
        else
          select_largest_stat_data_runs "$smoke_run_count" "$smoke_selected_runs" "$smoke_selected_stats" || { err "smokeTest could not select DATA runs from ${GOLDEN}"; exit 99; }
        fi
        smoke_out_base="${RJ_SMOKE_OUTPUT_BASE:-/sphenix/tg/tg01/bulk/jbennett/thesisAnaSmoke/${TAG}_smokeTest_${smoke_stamp}}"
        export RJ_DEST_BASE_OVERRIDE="$smoke_out_base"
        export RJ_MERGE_OUT_BASE_OVERRIDE="${RJ_MERGE_OUT_BASE_OVERRIDE:-${BASE}/outputSmoke/${TAG}_smokeTest_${smoke_stamp}}"
        export RJ_GOLDEN_OVERRIDE="$smoke_selected_runs"
        export RJ_PROFILE_JOB=1
        export RJ_DIRECT_NEVENTS="${RJ_SMOKE_DATA_NEVENTS:-3000}"
        export RJ_DIRECT_MAX_JOBS="${RJ_SMOKE_DATA_MAX_JOBS:-0}"
        export RJ_JOB_HEARTBEAT_SECONDS="${RJ_JOB_HEARTBEAT_SECONDS:-${RJ_SMOKE_JOB_HEARTBEAT_SECONDS:-120}}"
        export RJ_PROFILE_STAGE="${RJ_PROFILE_STAGE:-directSmoke}"
        export RJ_PROFILE_LABEL="${RJ_PROFILE_LABEL:-${TAG}_smokeTest}"
        say "${BOLD}DATA direct-fanout smokeTest requested${RST}"
        say "  dataset       : ${DATASET}"
        if [[ -n "${RJ_SMOKE_DATA_RUN:-}" ]]; then
          say "  selected runs : exact accepted run $(run8 "$RJ_SMOKE_DATA_RUN")"
        else
          say "  selected runs : ${smoke_run_count} largest-statistics golden runs"
        fi
        say "  run list      : ${smoke_selected_runs}"
        say "  run stats     : ${smoke_selected_stats}"
        say "  output base   : ${RJ_DEST_BASE_OVERRIDE}"
        say "  merge output  : ${RJ_MERGE_OUT_BASE_OVERRIDE}"
        say "  groupSize     : ${GROUP_SIZE}"
        say "  max jobs      : ${RJ_DIRECT_MAX_JOBS} (0 means all chunks in selected runs)"
        say "  nEvents/job   : ${RJ_DIRECT_NEVENTS} (0 means full worker input)"
        say "  request mem   : ${RJ_REQUEST_MEMORY}"
        say "  engine        : direct RecoilJets fanout; pool replay is not used"
        if [[ "${RJ_DAG_DRYRUN:-0}" == "1" || "${RJ_DAG_DRYRUN:-0}" == "true" || "${RJ_DAG_DRYRUN:-0}" == "TRUE" ]]; then
          "$0" "$DATASET" CHECKJOBS groupSize "$GROUP_SIZE"
          echo "RECOILJETS_SMOKETEST_DRYRUN_V1"
          echo "dataset=${DATASET}"
          echo "mode=condor smokeTest"
          echo "engine=direct_legacy_fanout"
          echo "legacy_output_parity=NOT_RUN_DIRECT_RECOILJETS_ENGINE"
          echo "runs_or_samples=$(grep -cE '^[0-9]+' "$smoke_selected_runs" 2>/dev/null || echo 0)"
          echo "selected_run_stats=${smoke_selected_stats}"
          echo "output_base=${RJ_DEST_BASE_OVERRIDE}"
          echo "merge_output_base=${RJ_MERGE_OUT_BASE_OVERRIDE}"
          echo "continuing_to_build_full_auto_dag=1"
        fi
        "$0" "$DATASET" condor all groupSize "$GROUP_SIZE"
        ;;
      round)
        seg="${4:?round number required}"
        firstChunk="${5:-}"
        round_file="${ROUND_DIR}/goldenRuns_${TAG}_segment${seg}.txt"
        [[ -s "$round_file" ]] || { err "Round file not found: $round_file. Run 'splitGoldenRunList' first."; exit 8; }

        # Keep existing YAML overrides/snapshots; live Condor jobs may still reference them.
        mkdir -p "$SIM_YAML_OVERRIDE_DIR"

        for data_pt in "${data_submit_pts[@]}"; do
        for data_frac in "${data_submit_fracs[@]}"; do
        for data_vz in "${data_vzs[@]}"; do
          for data_cone in "${data_cones[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
          direct_shard_idx=0
          direct_shard_total=0
          if iso_cone_fanout_enabled; then
            iso_idx_is_group_leader "$iso_idx" || continue
            direct_shard_idx="$(direct_fanout_cell_ordinal "$data_cone" "$iso_idx" "${data_cones[@]}")"
            direct_shard_total="$(direct_fanout_shard_count "${#data_cones[@]}" "${#iso_tags[@]}")"
            (( direct_shard_idx <= direct_shard_total )) || continue
          else
            id_fanout_enabled && ! iso_idx_is_group_leader "$iso_idx" && continue
          fi
          for uepipe in "${uepipe_modes[@]}"; do
          if iso_cone_fanout_enabled; then
            data_cfg_tag="$(direct_fanout_shard_tag "$data_pt" "$data_frac" "$data_vz" "$uepipe" "$direct_shard_idx")"
          else
            data_cfg_tag="$(matrix_cfg_tag "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe")"
          fi
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
          force_ppg12_photon_yield_yaml_contract "$yaml_override"
          pin_photon_id_scalars_in_yaml "$yaml_override" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}"
          if iso_cone_fanout_enabled; then
            emit_direct_fanout_shard_file "$fanout_dirs" "$DATA_DEST_BASE_SAVED" "$data_pt" "$data_frac" "$data_vz" "$uepipe" "$direct_shard_idx" "${data_cones[@]}"
          elif id_fanout_enabled; then
            emit_id_fanout_dirs_file "$fanout_dirs" "$DATA_DEST_BASE_SAVED" "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe"
          else
            : > "$fanout_dirs"
          fi
          export RJ_CONFIG_YAML="$yaml_override"
          DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"

          if iso_cone_fanout_enabled; then
            say "DATA condor round iso/cone/ID fanout (pt=${data_pt}, frac=${data_frac}, vz=${data_vz}, shard=${direct_shard_idx}/${direct_shard_total}, uepipe=${uepipe}, tag=${data_cfg_tag})"
          else
            say "DATA condor round (pt=${data_pt}, frac=${data_frac}, vz=${data_vz}, coneR=${data_cone}, iso=${iso_tags[$iso_idx]}, uepipe=${uepipe}, tag=${data_cfg_tag})"
          fi
          say "  YAML override: ${yaml_override}"
          if id_fanout_enabled; then
            say "  Fanout cfgs  : $(wc -l < "$fanout_dirs" | awk '{print $1}') output cfg tags from one DST pass"
            say "  Fanout file  : ${fanout_dirs}"
            while IFS='|' read -r fan_dest _fan_cfg _fan_pre _fan_tight _fan_nonTight; do
              [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
              fanout_dest_allowed "$fan_dest" || { err "Refusing to wipe fanout DEST_BASE='$fan_dest'"; exit 62; }
              if [[ "${RJ_CLEAN_OUTPUT_BASE:-0}" == "1" ]]; then
                mkdir -p "$fan_dest"
                find "$fan_dest" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true
              else
                mkdir -p "$fan_dest"
              fi
            done < "$fanout_dirs"
          else
            say "  Fanout cfgs  : disabled; this submits exactly ${data_cfg_tag}"
          fi
          _old_submit_extra_env="${RJ_SUBMIT_EXTRA_ENV:-}"
          if id_fanout_enabled; then
            export RJ_SUBMIT_EXTRA_ENV="${_old_submit_extra_env:+${_old_submit_extra_env};}RJ_ID_FANOUT_FILE=${fanout_dirs};RJ_ID_FANOUT_DIRS_FILE=${fanout_dirs}${jetpt_env_for_sub}${dphi_env_for_sub}${iso_view_env_for_sub}"
          else
            export RJ_SUBMIT_EXTRA_ENV="${_old_submit_extra_env}${jetpt_env_for_sub}${dphi_env_for_sub}${iso_view_env_for_sub}"
          fi
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
        warn "condor allFromScratch now uses direct RecoilJets fanout; pool capture/replay was removed."
        _alias_args=( "$DATASET" condor all groupSize "$GROUP_SIZE" )
        if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 ]]; then
          _alias_args+=( maxJobs "$MAX_JOBS" )
        fi
        "$0" "${_alias_args[@]}"
        ;;
      allDirect)
        warn "condor allDirect disables ID fanout and submits one legacy cfg-tag per DST pass for validation."
        _direct_args=( "$DATASET" condor all groupSize "$GROUP_SIZE" )
        if [[ "$MAX_JOBS" =~ ^[0-9]+$ && "$MAX_JOBS" -gt 0 ]]; then
          _direct_args+=( maxJobs "$MAX_JOBS" )
        fi
        RJ_DISABLE_ID_FANOUT=1 RJ_DISABLE_JET_PT_INTERNALIZATION=1 RJ_DISABLE_DPHI_INTERNALIZATION=1 RJ_DIRECT_DST_DOALL=1 "$0" "${_direct_args[@]}"
        ;;
      all|"")
        all_stamp="$(date +%Y%m%d_%H%M%S)"
        SIM_YAML_OVERRIDE_DIR="${SIM_YAML_OVERRIDE_DIR}/${TAG}_condorAll_${all_stamp}"
        say "DATA YAML/fanout artifact dir: ${SIM_YAML_OVERRIDE_DIR}"
        # Keep existing YAML overrides/snapshots; live Condor jobs may still reference them.
        mkdir -p "$SIM_YAML_OVERRIDE_DIR"

        # Freeze pipeline for this bulk submission
        cleanup_bulk_snapshots_for_tag
        case "$DATASET" in
          isAuAu|isOO) create_pipeline_snapshot "auau" "$all_stamp" ;;
          *)           create_pipeline_snapshot "pp" "$all_stamp" ;;
        esac

        iso_submit_n="$(iso_submit_count)"
        fanout_shards_per_axis=0
        final_cfg_n=0
        if iso_cone_fanout_enabled; then
          fanout_shards_per_axis="$(direct_fanout_shard_count "${#data_cones[@]}" "${#iso_tags[@]}")"
          n_matrix=$(( ${#data_submit_pts[@]} * ${#data_submit_fracs[@]} * ${#data_vzs[@]} * fanout_shards_per_axis * ${#uepipe_modes[@]} ))
          final_cfg_n=$(( ${#data_submit_pts[@]} * ${#data_submit_fracs[@]} * ${#data_vzs[@]} * ${#data_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
        else
          n_matrix=$(( ${#data_submit_pts[@]} * ${#data_submit_fracs[@]} * ${#data_vzs[@]} * ${#data_cones[@]} * iso_submit_n * ${#uepipe_modes[@]} ))
          final_cfg_n=$(( ${#data_submit_pts[@]} * ${#data_submit_fracs[@]} * ${#data_vzs[@]} * ${#data_cones[@]} * ${#iso_tags[@]} * ${#uepipe_modes[@]} ))
        fi
        say "═══════════════════════════════════════════════════════════════"
        say "${BOLD}CONDOR ALL: ${n_matrix} upstream matrix configuration(s), groupSize=${GROUP_SIZE}${RST}"
        say "  dataset       : ${DATASET}"
        say "  YAML source   : ${data_yaml_src}"
        say "  vz selection  : [${data_vzs[*]}] ($(vz_selection_reason_for_yaml "$data_yaml_src"))"
        say "  golden list   : $(basename "$GOLDEN") ($(grep -cE '^[0-9]+' "$GOLDEN") runs)"
        if id_fanout_enabled; then
          if iso_cone_fanout_enabled; then
            say "  engine        : direct legacy iso/cone/ID fanout (${final_cfg_n} final cfg outputs; ${fanout_shards_per_axis} shard(s) per pt/dphi/vz/UE axis; cap=$(id_fanout_max_rows))"
          else
            say "  engine        : direct RecoilJets modules (${#iso_tags[@]} final cfg ROOT files; ${iso_submit_n} photon-ID fanout shard(s), cap=${RJ_ID_FANOUT_MAX_ROWS:-15})"
          fi
        else
          say "  engine        : direct legacy one-cfg-per-pass validation (fanout disabled)"
        fi
        say "═══════════════════════════════════════════════════════════════"
        echo

        auto_workflow=0
        auto_dag=""
        auto_dag_dir=""
        if auto_merge_enabled; then
          auto_workflow=1
          need_cmd condor_submit_dag
          auto_stamp="$all_stamp"
          auto_dag_dir="${SUB_DIR}/auto_workflow_${TAG}_${auto_stamp}"
          mkdir -p "$auto_dag_dir"
          auto_dag="${auto_dag_dir}/RecoilJets_auto_${TAG}_${auto_stamp}.dag"
          : > "$auto_dag"
          write_cleanup_manifest "${auto_dag_dir}/cleanup_manifest.txt" "$SIM_YAML_OVERRIDE_DIR" "$DATA_DEST_BASE_SAVED"
          RJ_DAG_COLLECTED_NODES=()
          export RJ_COLLECT_DAG_FILE="$auto_dag"
          export RJ_COLLECT_SUB_DIR="${auto_dag_dir}/analysis_submit_files"
          export RJ_COLLECT_NODE_PREFIX="ANALYSIS_${TAG}"
          say "${BOLD}Automatic analysis-to-merge DAG enabled${RST}"
          say "  dag          : ${auto_dag}"
          say "  escape hatch : RJ_AUTO_MERGE=0 keeps analysis-only submission"
        fi

        cell_num=0
        total_queued_all=0
        t0_all="$(date +%s)"

        for data_pt in "${data_submit_pts[@]}"; do
        for data_frac in "${data_submit_fracs[@]}"; do
        for data_vz in "${data_vzs[@]}"; do
          for data_cone in "${data_cones[@]}"; do
          for (( iso_idx=0; iso_idx<${#iso_tags[@]}; iso_idx++ )); do
          direct_shard_idx=0
          direct_shard_total=0
          if iso_cone_fanout_enabled; then
            iso_idx_is_group_leader "$iso_idx" || continue
            direct_shard_idx="$(direct_fanout_cell_ordinal "$data_cone" "$iso_idx" "${data_cones[@]}")"
            direct_shard_total="$(direct_fanout_shard_count "${#data_cones[@]}" "${#iso_tags[@]}")"
            (( direct_shard_idx <= direct_shard_total )) || continue
          else
            id_fanout_enabled && ! iso_idx_is_group_leader "$iso_idx" && continue
          fi
          for uepipe in "${uepipe_modes[@]}"; do
          if iso_cone_fanout_enabled; then
            data_cfg_tag="$(direct_fanout_shard_tag "$data_pt" "$data_frac" "$data_vz" "$uepipe" "$direct_shard_idx")"
          else
            data_cfg_tag="$(matrix_cfg_tag "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe")"
          fi
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
          force_ppg12_photon_yield_yaml_contract "$yaml_override"
          pin_photon_id_scalars_in_yaml "$yaml_override" "${iso_preselection[$iso_idx]}" "${iso_tight[$iso_idx]}" "${iso_nonTight[$iso_idx]}"
          if iso_cone_fanout_enabled; then
            emit_direct_fanout_shard_file "$fanout_dirs" "$DATA_DEST_BASE_SAVED" "$data_pt" "$data_frac" "$data_vz" "$uepipe" "$direct_shard_idx" "${data_cones[@]}"
          elif id_fanout_enabled; then
            emit_id_fanout_dirs_file "$fanout_dirs" "$DATA_DEST_BASE_SAVED" "$data_pt" "$data_frac" "$data_vz" "$data_cone" "$iso_idx" "$uepipe"
          else
            : > "$fanout_dirs"
          fi
          export RJ_CONFIG_YAML="$yaml_override"
          DEST_BASE="${DATA_DEST_BASE_SAVED}/${data_cfg_tag}"

          (( ++cell_num ))
          say "───────────────────────────────────────────────────────────────"
          if iso_cone_fanout_enabled; then
            say "${BOLD}[${cell_num}/${n_matrix}] ${data_cfg_tag} (iso/cone/ID shard ${direct_shard_idx}/${direct_shard_total})${RST}"
          else
            say "${BOLD}[${cell_num}/${n_matrix}] ${data_cfg_tag}${RST}"
          fi
          say "  YAML override : ${yaml_override}"
          if id_fanout_enabled; then
            say "  Fanout cfgs   : $(wc -l < "$fanout_dirs" | awk '{print $1}') output cfg tags from one DST pass"
            say "  Fanout file   : ${fanout_dirs}"
          else
            say "  Fanout cfgs   : disabled; this submits exactly ${data_cfg_tag}"
          fi
          tmp_src="${ROUND_DIR}/ALL_${TAG}_${data_cfg_tag}_$(date +%s).txt"
          grep -E '^[0-9]+' "$GOLDEN" > "$tmp_src"

          if id_fanout_enabled; then
            while IFS='|' read -r fan_dest _fan_cfg _fan_pre _fan_tight _fan_nonTight; do
              [[ -z "${fan_dest:-}" || "${fan_dest:0:1}" == "#" ]] && continue
              fanout_dest_allowed "$fan_dest" || { err "Refusing to wipe fanout DEST_BASE='$fan_dest'"; exit 62; }
              if [[ "${RJ_CLEAN_OUTPUT_BASE:-0}" == "1" ]] && dag_dryrun_enabled; then
                say "DRYRUN: would clean fanout DEST_BASE=${fan_dest}"
              elif [[ "${RJ_CLEAN_OUTPUT_BASE:-0}" == "1" ]]; then
                mkdir -p "$fan_dest"
                find "$fan_dest" -mindepth 1 -maxdepth 1 -exec rm -rf {} + 2>/dev/null || true
              else
                mkdir -p "$fan_dest"
              fi
            done < "$fanout_dirs"
          fi

          _old_submit_extra_env="${RJ_SUBMIT_EXTRA_ENV:-}"
          if id_fanout_enabled; then
            export RJ_SUBMIT_EXTRA_ENV="${_old_submit_extra_env:+${_old_submit_extra_env};}RJ_ID_FANOUT_FILE=${fanout_dirs};RJ_ID_FANOUT_DIRS_FILE=${fanout_dirs}${jetpt_env_for_sub}${dphi_env_for_sub}${iso_view_env_for_sub}"
          else
            export RJ_SUBMIT_EXTRA_ENV="${_old_submit_extra_env}${jetpt_env_for_sub}${dphi_env_for_sub}${iso_view_env_for_sub}"
          fi
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
        say "${BOLD}CONDOR ALL complete: ${cell_num} upstream configurations submitted (${elapsed_all}s)${RST}"
        if (( auto_workflow )); then
          unset RJ_COLLECT_DAG_FILE
          unset RJ_COLLECT_SUB_DIR
          unset RJ_COLLECT_NODE_PREFIX
          runner="${auto_dag_dir}/auto_stage_runner.sh"
          final_notify="${auto_dag_dir}/auto_final_notify.sh"
          write_auto_stage_runner "$runner"
          write_auto_final_notify "$final_notify"
          data_merge_out_base="${RJ_MERGE_OUT_BASE_OVERRIDE:-${BASE}/output}"
          auto_final_addchunks="${RJ_AUTO_FINAL_ADDCHUNKS:-1}"
          if [[ "$auto_final_addchunks" != "0" ]]; then
            auto_final_stage_label="DATA_PERRUN -> DATA_SLICERUNS -> DATA_FINAL_ADDCHUNKS"
            auto_final_notify_key="auto_${TAG}_final_ready"
            auto_final_next_action="./scripts/sftp_get_recoiljets_outputs.sh ${DATASET}"
          else
            auto_final_stage_label="DATA_PERRUN -> DATA_SLICERUNS"
            auto_final_notify_key="auto_${TAG}_sliceRuns_ready"
            auto_final_next_action="MERGE_RUN_BASE_OVERRIDE=${DATA_DEST_BASE_SAVED} MERGE_OUT_BASE_OVERRIDE=${data_merge_out_base} ${BASE}/scripts/mergeRecoilJets.sh addChunks ${TAG}"
          fi

          add_auto_stage_node "$auto_dag" "DATA_PERRUN" "$runner" "data_perRun_${TAG}_all" env "MERGE_RUN_BASE_OVERRIDE=${DATA_DEST_BASE_SAVED}" "MERGE_OUT_BASE_OVERRIDE=${data_merge_out_base}" "RJ_STAGE_EMAIL_MODE=none" "RJ_STAGE_EMAIL_STRICT=1" "${BASE}/scripts/mergeRecoilJets.sh" condor "$TAG"
          add_auto_stage_node "$auto_dag" "DATA_SLICERUNS" "$runner" "data_sliceRuns_${TAG}_all" env "MERGE_RUN_BASE_OVERRIDE=${DATA_DEST_BASE_SAVED}" "MERGE_OUT_BASE_OVERRIDE=${data_merge_out_base}" "RJ_STAGE_EMAIL_MODE=none" "RJ_STAGE_EMAIL_STRICT=1" "${BASE}/scripts/mergeRecoilJets.sh" addChunks "$TAG" condor sliceRuns
          if [[ "$auto_final_addchunks" != "0" ]]; then
            add_auto_stage_node "$auto_dag" "DATA_FINAL_ADDCHUNKS" "$runner" "data_final_${TAG}_all" env "MERGE_RUN_BASE_OVERRIDE=${DATA_DEST_BASE_SAVED}" "MERGE_OUT_BASE_OVERRIDE=${data_merge_out_base}" "RJ_STAGE_EMAIL_MODE=none" "RJ_STAGE_EMAIL_STRICT=1" "${BASE}/scripts/mergeRecoilJets.sh" addChunks "$TAG" condor
          fi
          if (( ${#RJ_DAG_COLLECTED_NODES[@]} > 0 )); then
            printf 'PARENT' >> "$auto_dag"
            printf ' %s' "${RJ_DAG_COLLECTED_NODES[@]}" >> "$auto_dag"
            printf ' CHILD DATA_PERRUN\n' >> "$auto_dag"
          fi
          printf 'PARENT DATA_PERRUN CHILD DATA_SLICERUNS\n' >> "$auto_dag"
          if [[ "$auto_final_addchunks" != "0" ]]; then
            printf 'PARENT DATA_SLICERUNS CHILD DATA_FINAL_ADDCHUNKS\n' >> "$auto_dag"
          fi
          add_auto_final_node "$auto_dag" "FINAL_NOTIFY" "$final_notify" "$auto_final_notify_key" "$DATASET" "$auto_final_next_action" "${data_merge_out_base}/${TAG}"
          validate_auto_dag_submit_files "$auto_dag"
          say "Automatic workflow DAG built:"
          say "  analysis nodes : ${#RJ_DAG_COLLECTED_NODES[@]}"
          say "  merge stages   : ${auto_final_stage_label} (quiet strict validation; stage-boundary and final emails)"
          say "  dag            : ${auto_dag}"
          if dag_dryrun_enabled; then
            echo "RECOILJETS_AUTO_DAG_DRYRUN_V1"
            echo "dataset=${DATASET}"
            echo "dag=${auto_dag}"
            echo "analysis_nodes=${#RJ_DAG_COLLECTED_NODES[@]}"
            echo "auto_final_addchunks=${auto_final_addchunks}"
            echo "next_action_after_ready=${auto_final_next_action}"
            sed -n '1,240p' "$auto_dag"
          else
            condor_submit_dag -notification Never "$auto_dag"
          fi
        fi
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
