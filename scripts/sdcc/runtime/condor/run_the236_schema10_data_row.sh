#!/usr/bin/env bash
set -euo pipefail

# One immutable schema-10 sparse-data row.  The analysis writes its candidate
# directly beside the final group-bulk output as .root.part.  Only this worker,
# after ROOT/content/count validation, may atomically publish the .root name.

manifest="${1:-}"
manifest_sha="${2:-}"
system="${3:-}"
row_id="${4:-}"
records_path="${5:-}"
records_sha="${6:-}"
byte_offset="${7:-}"
byte_count="${8:-}"
chunk_sha="${9:-}"
expected_processed="${10:-}"
output_path="${11:-}"
measurement_path="${12:-}"
request_memory_mb="${13:-}"
event_limit="${14:-0}"
event_offset="${15:-0}"
source_total_events="${16:-0}"

die(){ printf '[THE236-SCHEMA10-DATA][ERROR] %s\n' "$*" >&2; exit 2; }
sha_file(){ /usr/bin/sha256sum "$1" | /usr/bin/awk '{print $1}'; }

[[ -n "${_CONDOR_SCRATCH_DIR:-}" && -d "$_CONDOR_SCRATCH_DIR" ]] || die "worker sandbox is absent"
[[ "$manifest" == /* && -f "$manifest" && ! -L "$manifest" ]] || die "manifest is absent"
[[ "$manifest_sha" =~ ^[0-9a-f]{64}$ && "$(sha_file "$manifest")" == "$manifest_sha" ]] || die "manifest differs"
[[ "$system" == pp || "$system" == auau ]] || die "system differs"
[[ "$row_id" =~ ^(pp|auau)_[a-z0-9_]+$ ]] || die "row identity differs"
[[ "$records_path" == /* && -f "$records_path" && ! -L "$records_path" ]] || die "source records are absent"
[[ "$records_sha" =~ ^[0-9a-f]{64}$ && "$(sha_file "$records_path")" == "$records_sha" ]] || die "source records differ"
[[ "$byte_offset" =~ ^[0-9]+$ && "$byte_count" =~ ^[1-9][0-9]*$ ]] || die "record range differs"
[[ "$chunk_sha" =~ ^[0-9a-f]{64}$ ]] || die "chunk identity differs"
[[ "$expected_processed" =~ ^[1-9][0-9]*$ ]] || die "exact processed-event count is required"
[[ "$output_path" == /sphenix/tg/tg01/bulk/*.root && ! -e "$output_path" && ! -L "$output_path" ]] || die "final output path is unsafe or occupied"
[[ "$measurement_path" == /sphenix/tg/tg01/bulk/*.json && ! -e "$measurement_path" && ! -L "$measurement_path" ]] || die "measurement path is unsafe or occupied"
[[ "$request_memory_mb" =~ ^[0-9]+$ && "$request_memory_mb" -le 4096 ]] || die "memory request exceeds contract"
[[ "$event_limit" =~ ^[0-9]+$ ]] || die "event limit differs"
[[ "$event_offset" =~ ^[0-9]+$ ]] || die "event offset differs"
[[ "$source_total_events" =~ ^[1-9][0-9]*$ ]] || die "source total-event count differs"
if [[ "$event_limit" -gt 0 ]]; then
  [[ "$expected_processed" -eq "$event_limit" ]] || die "bounded canary processed-event contract differs"
fi
if [[ "$system" == auau ]]; then
  [[ "$event_limit" -gt 0 && "$event_limit" -le 35000 ]] || die "AuAu DST-to-TTree row exceeds the 35000-event hard cap"
  [[ "$source_total_events" -ge $((event_offset + event_limit)) ]] || die "AuAu event range exceeds its paired source"
else
  [[ "$event_offset" -eq 0 && "$source_total_events" -eq "$expected_processed" ]] || die "pp event-range fields differ"
fi
[[ "${ClusterId:-}" =~ ^[0-9]+$ && "${ProcId:-}" =~ ^[0-9]+$ ]] || die "scheduler identity is absent"

runtime="$({ python3 - "$manifest" "$system" <<'PY'
import hashlib,json,re,sys
from pathlib import Path

payload=json.loads(Path(sys.argv[1]).read_text(encoding="utf-8")); system=sys.argv[2]
if payload.get("schema")!="THE236Schema10SparseDataProductionManifestV1" or payload.get("status")!="PASS_FROZEN_UNSUBMITTED":
    raise SystemExit("manifest state differs")
if payload.get("data_retention_profile")!="sparse_photon_analysis_v1":
    raise SystemExit("retention profile differs")
runtime=payload.get("runtime",{}); binding=runtime.get(system,{})
required=("wrapper_path","wrapper_sha256","macro_path","macro_sha256","config_path","config_sha256","library_path","library_sha256","model_sha256","install_prefix","release","offline_main")
for key in required:
    if not isinstance(binding.get(key),str) or not binding[key]: raise SystemExit(f"runtime field absent: {key}")
for key in ("wrapper","macro","config"):
    path=Path(binding[f"{key}_path"])
    if not path.is_absolute() or not path.is_file() or path.is_symlink(): raise SystemExit(f"runtime file invalid: {key}")
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest!=binding[f"{key}_sha256"]: raise SystemExit(f"runtime file differs: {key}")
library=Path(binding["library_path"])
if not library.is_absolute() or not library.is_file(): raise SystemExit("runtime library invalid")
resolved_library=library.resolve(strict=True)
if not resolved_library.is_file() or resolved_library.is_symlink(): raise SystemExit("runtime library target invalid")
if library.is_symlink() and resolved_library.parent!=library.parent.resolve(strict=True): raise SystemExit("runtime library link escapes install lib")
if hashlib.sha256(resolved_library.read_bytes()).hexdigest()!=binding["library_sha256"]: raise SystemExit("runtime library differs")
install=Path(binding["install_prefix"])
if not install.is_absolute() or not install.is_dir() or install.is_symlink(): raise SystemExit("install prefix differs")
environment=payload.get("replay_environment",{})
if not isinstance(environment,dict): raise SystemExit("replay environment differs")
protected={"RJ_SCHEMA10_DATA_RETENTION_PROFILE","RJ_REPLAY_CONFIG_SHA256","RJ_REPLAY_DATASET","RJ_REPLAY_LANE","RJ_REPLAY_SAMPLE","RJ_REPLAY_SI_DI_ROLE","RJ_REPLAY_OWNERSHIP_STATE","RJ_REPLAY_INPUT_URI_SHA256","RJ_REPLAY_INPUT_FILE_SHA256","RJ_REPLAY_SOURCE_MANIFEST_SHA256","RJ_REPLAY_SOURCE_SHA256","RJ_REPLAY_PROVENANCE_MANIFEST_SHA256","RJ_OUTPUT_ROOT_CANDIDATE","RJ_REQUEST_MEMORY_MB","RJ_DATASET","RJ_IS_SIM","RJ_PP_INSTALL_PREFIX","RJ_AUAU_INSTALL_PREFIX","RJ_MACRO_PATH","RJ_CONFIG_YAML","RJ_RUNTIME_BASE","RJ_PINNED_RELEASE_NAME","RJ_PINNED_OFFLINE_MAIN"}
for key,value in environment.items():
    if not re.fullmatch(r"RJ_[A-Z0-9_]+",key) or not isinstance(value,str) or any(c in value for c in "\t\r\n"):
        raise SystemExit("unsafe environment")
    if key in protected: raise SystemExit(f"protected environment override: {key}")
if not re.fullmatch(r"[0-9a-f]{64}",binding["model_sha256"]): raise SystemExit("model hash invalid")
print("\t".join((binding["wrapper_path"],binding["macro_path"],binding["config_path"],binding["install_prefix"],binding["release"],binding["offline_main"],binding["config_sha256"],binding["model_sha256"])))
for key,value in sorted(environment.items()): print("ENV\t"+key+"\t"+value)
PY
} 2>&1)" || die "runtime binding validation failed: $runtime"

first_line="${runtime%%$'\n'*}"
IFS=$'\t' read -r wrapper macro config install_prefix release offline_main config_sha model_sha <<<"$first_line"
environment_lines=""; [[ "$runtime" == *$'\n'* ]] && environment_lines="${runtime#*$'\n'}"

while IFS=$'\t' read -r marker key value; do
  [[ -z "$marker" ]] && continue
  [[ "$marker" == ENV && "$key" == RJ_* ]] || die "replay environment row differs"
  export "$key=$value"
done <<<"$environment_lines"

chunk_path="${_CONDOR_SCRATCH_DIR}/source_chunk.tsv"
/bin/dd if="$records_path" of="$chunk_path" iflag=skip_bytes,count_bytes skip="$byte_offset" count="$byte_count" status=none
[[ -s "$chunk_path" && "$(sha_file "$chunk_path")" == "$chunk_sha" ]] || die "source record slice differs"

gate_output_path=""
gate_candidate_path=""
gate_receipt_path=""
gate_receipt_candidate=""
gate_validator=""
gate_source_pairs=0
gate_run=0
if [[ "$system" == auau ]]; then
  gate_binding="$({ python3 - "$manifest" <<'PY'
import hashlib,json,re,sys
from pathlib import Path

manifest=Path(sys.argv[1])
payload=json.loads(manifest.read_text(encoding="utf-8"))
contract=payload.get("auau_event_gate_contract",{})
if contract.get("schema")!="AuAuInlineEventGateProductionContractV1":
    raise SystemExit("AuAu event-gate contract is absent")
if contract.get("required_for_sparse_auau_data") is not True:
    raise SystemExit("AuAu event-gate contract is not mandatory")
if contract.get("photon10_scaled_bit")!=22:
    raise SystemExit("AuAu Photon10 scaled-bit contract differs")
name=contract.get("validator_name")
digest=contract.get("validator_sha256")
if not isinstance(name,str) or not re.fullmatch(r"[A-Za-z0-9_.-]+",name):
    raise SystemExit("AuAu event-gate validator name differs")
if not isinstance(digest,str) or not re.fullmatch(r"[0-9a-f]{64}",digest):
    raise SystemExit("AuAu event-gate validator hash differs")
path=manifest.parent/name
if not path.is_file() or path.is_symlink():
    raise SystemExit("AuAu event-gate validator is absent")
if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
    raise SystemExit("AuAu event-gate validator content differs")
sharding=payload.get("auau_dst_ttree_sharding_contract",{})
required={
    "schema":"AuAuDSTTTreeShardingContractV1",
    "required_for_sparse_auau_data":True,
    "workload_profile":"AUAU_DATA_DST_TO_TTREE",
    "source_pairs_per_base_unit":1,
    "max_events_per_job":35000,
    "event_range_partitioning":"DETERMINISTIC_CONTIGUOUS_OFFSET_COUNT_V1",
    "automatic_retry":False,
}
for key,value in required.items():
    if sharding.get(key)!=value:
        raise SystemExit(f"AuAu DST-to-TTree sharding contract differs: {key}")
print(f"{path}\t{digest}")
PY
  } 2>&1)" || die "AuAu event-gate binding failed: $gate_binding"
  IFS=$'\t' read -r gate_validator gate_validator_sha <<<"$gate_binding"
  [[ -n "$gate_validator" && "$gate_validator_sha" =~ ^[0-9a-f]{64}$ ]] || die "AuAu event-gate validator binding differs"

  gate_output_path="${output_path%.root}.auau_event_gate.root"
  gate_candidate_path="${gate_output_path}.part"
  gate_receipt_path="${measurement_path%.json}.auau_event_gate.json"
  gate_receipt_candidate="${gate_receipt_path}.part"
  for path in "$gate_output_path" "$gate_candidate_path" "$gate_receipt_path" "$gate_receipt_candidate"; do
    [[ ! -e "$path" && ! -L "$path" ]] || die "AuAu event-gate path is occupied: $path"
  done
  gate_source_pairs="$(/usr/bin/awk 'NF{n++} END{print n+0}' "$chunk_path")"
  [[ "$gate_source_pairs" -eq 1 ]] || die "AuAu DST-to-TTree row must contain exactly one paired source unit"
  gate_run="$({ python3 - "$chunk_path" <<'PY'
import re,sys
from pathlib import Path
runs=set()
for raw in Path(sys.argv[1]).read_text(encoding="utf-8").splitlines():
    if not raw: continue
    fields=raw.split("\t")
    if len(fields)!=2: raise SystemExit("source pair differs")
    for value in fields:
        match=re.search(r"-([0-9]{8})-[0-9]{5}[.]root$",value)
        if not match: raise SystemExit("source run identity differs")
        runs.add(int(match.group(1)))
if len(runs)!=1: raise SystemExit("source row spans runs")
print(runs.pop())
PY
  } 2>&1)" || die "AuAu event-gate run binding failed: $gate_run"
  [[ "$gate_run" =~ ^[1-9][0-9]*$ ]] || die "AuAu event-gate run differs"
fi

candidate_path="${output_path}.part"
mkdir -p "$(dirname "$output_path")" "$(dirname "$measurement_path")"
[[ ! -e "$candidate_path" && ! -L "$candidate_path" ]] || die "candidate path is occupied"

export RJ_RUNTIME_BASE="$(dirname "$(dirname "$macro")")"
export RJ_MACRO_PATH="$macro" RJ_CONFIG_YAML="$config"
export RJ_PINNED_RELEASE_NAME="$release" RJ_PINNED_OFFLINE_MAIN="$offline_main"
export RJ_REPLAY_FOUNDATION_V1=1 RJ_SCHEMA10_DATA_RETENTION_PROFILE=sparse_photon_analysis_v1
export RJ_REPLAY_CONFIG_SHA256="$config_sha" RJ_REPLAY_DATASET="${system}_data"
export RJ_REPLAY_MODEL_SHA256="$model_sha"
# The worker bound a model hash above but never declared which shower
# definition that model was trained against. RecoilJets.cc / RecoilJets_AuAu.cc
# read RJ_REPLAY_MODEL_SHOWER_DEFINITION to populate ModelEvaluationRow, and
# RJReplayFoundationV1.h rejects the row when shower_definition_id is empty or
# shower_semantic_sha256 is not valid hex64:
#   "model evaluation shower semantics or foreign-key failure"
# That aborts the whole event transaction, which is why every candidate-bearing
# row (pp_photon_positive, auau_photon_positive) failed while the zero-candidate
# "typical" rows passed. It stayed invisible for five canary generations because
# batch output was being discarded; see RJ_VERBOSITY above.
#
# Value is forced by which models exist, per shower_definition_factorial_contract.yaml:
#   H70 role=nominal, H0 role=primary_control (predecessor THE-111)
#   model_state: pp_H70 REUSE_THE116..., auau_H70 RETRAIN_REQUIRED,
#                auau_H0 REUSE_THE111_AS_CONTROL
# The AuAu H70 model does not exist yet, so AuAu data runs the THE-111 H0
# control and must declare H0. The contract forbids transferring THE-111
# working points to H70, so labelling this H70 would be the prohibited
# cross-model conflation. pp has an H70 model and declares H70.
# Matches materialize_the121_the122_harmonized_canary.py and
# resolve_the134_full_multiview_extraction.py, which both use
# "H70" if system == "pp" else "H0".
if [[ "$system" == pp ]]; then
  export RJ_REPLAY_MODEL_SHOWER_DEFINITION=H70
else
  export RJ_REPLAY_MODEL_SHOWER_DEFINITION=H0
fi
export RJ_REPLAY_LANE="${system}_data" RJ_REPLAY_SAMPLE=data RJ_REPLAY_SI_DI_ROLE=DATA
export RJ_REPLAY_OWNERSHIP_STATE="$(basename "$(dirname "$manifest")")"
export RJ_REPLAY_INPUT_URI_SHA256="$chunk_sha" RJ_REPLAY_INPUT_FILE_SHA256="$records_sha"
export RJ_REPLAY_SOURCE_MANIFEST_SHA256="$manifest_sha" RJ_REPLAY_SOURCE_SHA256="$manifest_sha"
export RJ_REPLAY_PROVENANCE_MANIFEST_SHA256="$manifest_sha"
export RJ_OUTPUT_ROOT_CANDIDATE="$candidate_path"
export RJ_REQUIRE_NON_TINY_OUTPUT=1 RJ_MIN_OUTPUT_BYTES=1000 RJ_JOB_HEARTBEAT_SECONDS=0 RJ_PROFILE_JOB=0
# Fun4All_recoilJets_unified_impl.C defaults vlevel to 0 whenever
# _CONDOR_SCRATCH_DIR is set, and vlevel==0 enables ScopedSilence, which
# redirects BOTH std::cout and std::cerr to /dev/null for the entire run.
# Every LOG(0,...) FATAL the analysis emits is therefore discarded in batch.
# That is why cluster 2391573's candidate-bearing rows aborted with
# reason=END_NONZERO and left no explanation in either log stream. Setting
# RJ_VERBOSITY=1 keeps the silence off so fatal diagnostics reach stdout,
# while LOG(10)-level chatter stays gated behind Verbosity() >= 10.
export RJ_VERBOSITY=1
export RJ_REQUEST_MEMORY_MB="$request_memory_mb"
if [[ "$system" == pp ]]; then
  export RJ_PP_INSTALL_PREFIX="$install_prefix" RJ_DATASET=isPP RJ_IS_SIM=0
  unset RJ_AUAU_INSTALL_PREFIX
  dataset=isPP
else
  export RJ_AUAU_INSTALL_PREFIX="$install_prefix" RJ_DATASET=isAuAu RJ_IS_SIM=0
  export RJ_AUAU_EVENT_GATE_OUTPUT_CANDIDATE="$gate_candidate_path"
  export RJ_AUAU_EVENT_GATE_ROW_ID="$row_id"
  export RJ_AUAU_EVENT_GATE_SOURCE_PAIRS="$gate_source_pairs"
  export RJ_EVENT_OFFSET="$event_offset"
  unset RJ_PP_INSTALL_PREFIX
  dataset=isAuAu
fi

time_raw="${_CONDOR_SCRATCH_DIR}/resource.time"
set +e
/usr/bin/time -v -o "$time_raw" bash "$wrapper" data "$chunk_path" "$dataset" "$ClusterId" "$event_limit" "$ProcId" NONE "$(dirname "$output_path")"
wrapper_status=$?
set -e
[[ "$wrapper_status" -eq 0 ]] || die "analysis wrapper exited ${wrapper_status}; candidate preserved at ${candidate_path}"
[[ -s "$candidate_path" ]] || die "candidate ROOT is absent"

if ! root_contract="$({
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n "$release" >/dev/null
  set -u
  # Use the same PyROOT object-reading path as the certified foreground
  # validator.  The earlier Cling one-liner could falsely return an empty
  # TNamed title on valid persisted metadata and turn successful science into
  # exit 2.  The sealed external count remains an independent check below.
  python3 - "$candidate_path" <<'PY'
import re
import sys

import ROOT

path = sys.argv[1]
ROOT.gROOT.SetBatch(True)
root_file = ROOT.TFile.Open(path, "READ")
if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
    raise SystemExit("RJ_DATA_ROOT_CONTRACT_V1 invalid ROOT candidate")
directory = root_file.GetDirectory("ReplayFoundationV1")
if not directory:
    raise SystemExit("RJ_DATA_ROOT_CONTRACT_V1 missing ReplayFoundationV1")

def title(key: str) -> str:
    value = directory.Get(key)
    if not value or not value.InheritsFrom("TNamed"):
        return ""
    return str(value.GetTitle())

def count(key: str) -> int:
    value = directory.Get(key)
    if not value:
        return 0
    if not value.InheritsFrom("TTree"):
        raise SystemExit(f"RJ_DATA_ROOT_CONTRACT_V1 non-tree object key={key}")
    return int(value.GetEntries())

def parse_nonnegative(key: str) -> int:
    raw = title(key)
    if not re.fullmatch(r"[0-9]+", raw):
        raise SystemExit(f"RJ_DATA_ROOT_CONTRACT_V1 invalid integer key={key} raw={raw}")
    return int(raw)

processed = parse_nonnegative("processed_events")
retained = parse_nonnegative("retained_events")
omitted = parse_nonnegative("omitted_events")
primary = parse_nonnegative("retained_primary_events")
extension = parse_nonnegative("retained_extension_events")
primary_candidates = parse_nonnegative("retained_primary_candidates")
extension_candidates = parse_nonnegative("retained_extension_candidates")

trees = (
    "RJEventV1", "RJPhotonCandidateV1", "RJModelEvaluationV1",
    "RJShowerCellV1", "RJShowerFeatureViewV1", "RJIsolationConstituentV1",
    "RJIsolationWitnessV1", "RJJetV1", "RJJetConstituentV1",
    "RJPhotonJetPairV1", "RJTruthPhotonV1",
    "RJTruthPhotonMissOccurrenceV1", "RJEmbeddedPhotonDiagnosticOccurrenceV1",
    "RJPPG12DiagnosticOccurrenceV1", "RJTruthJetV1", "RJRecoTruthLinkV1",
    "RJWeightComponentV1", "RJEventDisplaySnapshotV1",
)
entries = {}
for tree in trees:
    expected = parse_nonnegative(f"{tree}_entries")
    actual = count(tree)
    if expected != actual:
        raise SystemExit(
            f"RJ_DATA_ROOT_CONTRACT_V1 table mismatch tree={tree} "
            f"metadata={expected} actual={actual}"
        )
    entries[tree] = actual

events = entries["RJEventV1"]
truth_photons = entries["RJTruthPhotonV1"]
truth_jets = entries["RJTruthJetV1"]
links = entries["RJRecoTruthLinkV1"]
jets = entries["RJJetV1"]
pairs = entries["RJPhotonJetPairV1"]
if (
    title("rj_replay_complete") != "1"
    or title("rj_replay_schema_version") != "10"
    or title("data_retention_profile") != "sparse_photon_analysis_v1"
    or processed != retained + omitted
    or retained != events
    or retained != primary + extension
    or primary_candidates < primary
    or extension_candidates < extension
    or truth_photons != 0
    or truth_jets != 0
    or links != 0
    or (primary == 0 and (jets != 0 or pairs != 0))
):
    raise SystemExit(
        "RJ_DATA_ROOT_CONTRACT_V1 invariant failure "
        f"processed={processed} retained={retained} omitted={omitted} events={events} "
        f"primary={primary} extension={extension} "
        f"primary_candidates={primary_candidates} extension_candidates={extension_candidates} "
        f"truthPhotons={truth_photons} truthJets={truth_jets} links={links} "
        f"jets={jets} pairs={pairs}"
    )
print(
    f"RJ_DATA_ROOT_CONTRACT_V1 processed={processed} retained={retained} "
    f"primary={primary} extension={extension}"
)
root_file.Close()
PY
} 2>&1)"; then
  printf '%s\n' "$root_contract" >&2
  die "candidate ROOT validation failed; candidate preserved at ${candidate_path}"
fi

processed_actual="$({ python3 - "$root_contract" <<'PY'
import re,sys
match=re.search(r"RJ_DATA_ROOT_CONTRACT_V1 processed=(\d+)",sys.argv[1])
if not match: raise SystemExit("processed count is absent")
print(match.group(1))
PY
} 2>&1)" || die "validated processed count projection failed: $processed_actual"
[[ "$processed_actual" =~ ^[1-9][0-9]*$ ]] || die "validated processed count differs"

if [[ "$system" == auau ]]; then
  [[ -s "$gate_candidate_path" ]] || die "mandatory AuAuEventGateV1 candidate is absent"
  if ! gate_validation="$({
    set +u
    source /opt/sphenix/core/bin/sphenix_setup.sh -n "$release" >/dev/null
    set -u
    python3 "$gate_validator" \
      --gate "$gate_candidate_path" \
      --base "$candidate_path" \
      --row-id "$row_id" \
      --run "$gate_run" \
      --photon10-bit 22 \
      --expected-processed "$processed_actual" \
      --expected-source-pairs "$gate_source_pairs" \
      --receipt "$gate_receipt_candidate"
  } 2>&1)"; then
    printf '%s\n' "$gate_validation" >&2
    die "AuAu event-gate companion validation failed; candidates preserved"
  fi
  [[ -s "$gate_receipt_candidate" ]] || die "AuAu event-gate validation receipt is absent"
fi

measurement_candidate="${measurement_path}.part"
[[ ! -e "$measurement_candidate" && ! -L "$measurement_candidate" ]] || die "measurement candidate path is occupied"
python3 - "$time_raw" "$measurement_candidate" "$system" "$row_id" "$ClusterId" "$ProcId" "$request_memory_mb" "$candidate_path" "$output_path" "$manifest_sha" "$chunk_sha" "$expected_processed" "$root_contract" "$event_limit" "$gate_candidate_path" "$gate_output_path" "$gate_receipt_candidate" "$gate_receipt_path" "$event_offset" "$source_total_events" <<'PY'
import hashlib,json,math,os,re,sys
from pathlib import Path
raw=Path(sys.argv[1]); receipt=Path(sys.argv[2]); candidate=Path(sys.argv[8]); output=Path(sys.argv[9])
peaks=[int(m.group(1)) for line in raw.read_text(encoding="utf-8").splitlines() if (m:=re.fullmatch(r"\s*Maximum resident set size \(kbytes\):\s*([0-9]+)\s*",line))]
if len(peaks)!=1 or peaks[0]<=0: raise SystemExit("peak memory is absent")
peak_mb=math.ceil(peaks[0]/1024); request=int(sys.argv[7])
if peak_mb>request: raise SystemExit("peak memory exceeds request")
digest=hashlib.sha256(candidate.read_bytes()).hexdigest()
match=re.search(r"RJ_DATA_ROOT_CONTRACT_V1 processed=(\d+) retained=(\d+) primary=(\d+) extension=(\d+)",sys.argv[13])
if not match: raise SystemExit("ROOT validation summary is absent")
expected_processed=int(sys.argv[12]); processed_events=int(match.group(1)); event_limit=int(sys.argv[14])
event_offset=int(sys.argv[19]); source_total_events=int(sys.argv[20])
shortfall=expected_processed-processed_events
# An unbounded production row can reach a clean synchronized EOF one event
# before the frozen catalog total when the paired DST streams differ only at
# their terminal boundary.  The wrapper has already required a successful
# multi-input EOF and the ROOT contract independently records the events that
# were actually processed.  The same provenance-bounded boundary can omit two
# terminal records; larger source omissions still fail closed.  Bounded
# canaries remain exact.
if event_limit>0:
    terminal_range=(event_offset+event_limit)==source_total_events
    if shortfall==0:
        processed_count_contract="EXACT_BOUNDED_EVENT_RANGE_V1"
    elif sys.argv[3]=="auau" and terminal_range and shortfall==1:
        processed_count_contract="VERIFIED_PAIRED_DST_TERMINAL_RANGE_EOF_MINUS_ONE_V1"
    elif sys.argv[3]=="auau" and terminal_range and shortfall==2:
        processed_count_contract="VERIFIED_PAIRED_DST_TERMINAL_RANGE_EOF_MINUS_TWO_V1"
    else:
        raise SystemExit(f"bounded processed-event count differs expected={expected_processed} actual={processed_events}")
elif shortfall==0:
    processed_count_contract="EXACT_EOF_V1"
elif shortfall==1:
    processed_count_contract="VERIFIED_PAIRED_DST_EOF_MINUS_ONE_V1"
elif shortfall==2:
    processed_count_contract="VERIFIED_PAIRED_DST_EOF_MINUS_TWO_V1"
else:
    raise SystemExit(f"processed-event count differs expected={expected_processed} actual={processed_events}")
payload={"schema":"ResourceMeasurementV2","status":"MEASURED_TERMINAL","system":sys.argv[3],"row_id":sys.argv[4],
 "cluster_id":int(sys.argv[5]),"process_id":int(sys.argv[6]),"request_memory_mb":request,
 "raw_peak_memory_kb":peaks[0],"peak_memory_mb":peak_mb,"scientific_output_path":str(output),
 "scientific_output_sha256":digest,"scientific_output_size_bytes":candidate.stat().st_size,
 "production_manifest_sha256":sys.argv[10],"chunk_sha256":sys.argv[11],
 "expected_processed_events":expected_processed,"processed_events":processed_events,
 "event_offset":event_offset,"event_limit":event_limit,"source_total_events":source_total_events,
 "processed_event_shortfall":shortfall,"processed_count_contract":processed_count_contract,
 "retained_events":int(match.group(2)),"retained_primary_events":int(match.group(3)),
 "retained_extension_events":int(match.group(4)),"wrapper_exit_code":0,"publication":"ATOMIC_PART_TO_ROOT_V1"}
gate_candidate=Path(sys.argv[15]) if sys.argv[15] else None
gate_output=Path(sys.argv[16]) if sys.argv[16] else None
gate_receipt_candidate=Path(sys.argv[17]) if sys.argv[17] else None
gate_receipt=Path(sys.argv[18]) if sys.argv[18] else None
if sys.argv[3]=="auau":
    if not gate_candidate or not gate_candidate.is_file() or not gate_receipt_candidate or not gate_receipt_candidate.is_file():
        raise SystemExit("mandatory AuAu event-gate publication inputs are absent")
    gate_validation=json.loads(gate_receipt_candidate.read_text(encoding="utf-8"))
    if gate_validation.get("schema")!="AuAuInlineEventGateValidationReceiptV1" or gate_validation.get("status")!="PASS":
        raise SystemExit("mandatory AuAu event-gate validation receipt differs")
    payload["auau_event_gate_companion"]={
        "status":"PASS",
        "output_path":str(gate_output),
        "output_sha256":hashlib.sha256(gate_candidate.read_bytes()).hexdigest(),
        "output_size_bytes":gate_candidate.stat().st_size,
        "validation_receipt_path":str(gate_receipt),
        "validation_receipt_sha256":hashlib.sha256(gate_receipt_candidate.read_bytes()).hexdigest(),
        "publication":"COMPANION_AND_RECEIPTS_BEFORE_BASE_COMMIT_V1",
    }
elif any(sys.argv[index] for index in range(15,19)):
    raise SystemExit("pp row unexpectedly carries an AuAu event-gate artifact")
receipt.write_text(json.dumps(payload,sort_keys=True,separators=(",",":"))+"\n",encoding="utf-8")
os.chmod(receipt,0o644)
PY

# Publish companion products and receipts first; publish the base ROOT last.
# The accepted base filename is therefore the commit marker proving that its
# mandatory event-gate witness was already validated and made durable.  No move
# crosses a filesystem boundary.
if [[ "$system" == auau ]]; then
  /bin/mv -- "$gate_candidate_path" "$gate_output_path"
  /bin/mv -- "$gate_receipt_candidate" "$gate_receipt_path"
fi
/bin/mv -- "$measurement_candidate" "$measurement_path"
/bin/mv -- "$candidate_path" "$output_path"

printf '[THE236-SCHEMA10-DATA] PASS row=%s output=%s\n' "$row_id" "$output_path"
