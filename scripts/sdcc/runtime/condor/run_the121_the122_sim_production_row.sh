#!/usr/bin/env bash
set -euo pipefail

# Immutable worker for one THE-121/THE-122 simulation production chunk.
# The login-node materializer creates only data descriptors and chunk lists;
# all Fun4All execution happens here in the Condor worker sandbox.

manifest="${1:-}"
manifest_sha="${2:-}"
sample_id="${3:-}"
tuple_records_path="${4:-}"
tuple_records_sha="${5:-}"
byte_offset="${6:-}"
byte_count="${7:-}"
chunk_sha="${8:-}"
event_upper_bound="${9:-}"
output_path="${10:-}"
measurement_path="${11:-}"
request_memory_mb="${12:-}"
original_segment_id="${13:-${ProcId:-}}"

die(){ printf '[THE121-122-PRODUCTION-WORKER][ERROR] %s\n' "$*" >&2; exit 2; }
sha_file(){ /usr/bin/sha256sum "$1" | /usr/bin/awk '{print $1}'; }

[[ -n "${_CONDOR_SCRATCH_DIR:-}" && -d "${_CONDOR_SCRATCH_DIR}" ]] || die "worker sandbox is absent"
[[ "$manifest" == /* && -f "$manifest" && ! -L "$manifest" ]] || die "production manifest is absent"
[[ "$manifest_sha" =~ ^[0-9a-f]{64}$ && "$(sha_file "$manifest")" == "$manifest_sha" ]] || die "production manifest differs"
packet_root="$(dirname "$manifest")"
[[ "$tuple_records_path" == "$packet_root"/records/*.tsv && -f "$tuple_records_path" && ! -L "$tuple_records_path" ]] || die "tuple-record file is absent"
[[ "$tuple_records_sha" =~ ^[0-9a-f]{64}$ ]] || die "tuple-record identity is invalid"
[[ "$byte_offset" =~ ^[0-9]+$ && "$byte_count" =~ ^[1-9][0-9]*$ ]] || die "tuple-record range is invalid"
[[ "$chunk_sha" =~ ^[0-9a-f]{64}$ ]] || die "chunk identity is invalid"
[[ "$event_upper_bound" =~ ^[1-9][0-9]*$ ]] || die "event-count upper bound is invalid"
[[ "$output_path" == /sphenix/tg/tg01/bulk/* && "$output_path" == *.root && ! -e "$output_path" ]] || die "output path is unsafe or occupied"
[[ "$measurement_path" == /sphenix/tg/tg01/bulk/* && "$measurement_path" == *.json && ! -e "$measurement_path" ]] || die "measurement path is unsafe or occupied"
[[ "$request_memory_mb" =~ ^[0-9]+$ && "$request_memory_mb" -le 4096 ]] || die "memory request is invalid"
[[ "${ClusterId:-}" =~ ^[0-9]+$ && "${ProcId:-}" =~ ^[0-9]+$ ]] || die "scheduler identity is absent"
[[ "$original_segment_id" =~ ^[0-9]+$ ]] || die "original segment identity is invalid"

# Persist a completed scientific candidate in the final output filesystem
# before metadata validation. A validator or receipt defect must never make a
# successfully processed ROOT file disappear with the Condor scratch sandbox.
candidate_path="${output_path}.part"
measurement_candidate="${measurement_path}.part"
[[ ! -e "$candidate_path" && ! -L "$candidate_path" ]] || die "candidate path is occupied"
[[ ! -e "$measurement_candidate" && ! -L "$measurement_candidate" ]] || die "measurement candidate path is occupied"

values="$(python3 - "$manifest" "$sample_id" <<'PY'
import hashlib, json, re, sys
from pathlib import Path

payload=json.loads(Path(sys.argv[1]).read_text(encoding="utf-8")); sample_id=sys.argv[2]
if payload.get("schema")!="THE121_THE122_SIM_PRODUCTION_MANIFEST_V1" or payload.get("status")!="PASS_FROZEN":
    raise SystemExit("production manifest state differs")
samples=payload.get("samples")
if not isinstance(samples,dict) or sample_id not in samples:
    raise SystemExit("sample is not admitted")
r=samples[sample_id]; runtime=payload.get("runtime")
if not isinstance(runtime,dict): raise SystemExit("runtime binding is absent")
for key in ("system","lane","dataset","sample","si_di_role"):
    if not isinstance(r.get(key),str) or not r[key]: raise SystemExit(f"sample field absent: {key}")
for key in ("release","durable_runtime_root","dependency_runtime_root","pp_install_prefix","auau_install_prefix",
            "pp_library_path","pp_library_sha256","auau_library_path","auau_library_sha256",
            "dependency_manifest_path","dependency_manifest_sha256",
            "runtime_seal_receipt_path","runtime_seal_receipt_sha256"):
    if not isinstance(runtime.get(key),str) or not runtime[key]: raise SystemExit(f"runtime field absent: {key}")
system=r["system"]
if system not in {"pp","auau"}: raise SystemExit("system differs")
prefix="pp" if system=="pp" else "auau"
for key in (f"{prefix}_wrapper_path",f"{prefix}_wrapper_sha256",f"{prefix}_macro_path",f"{prefix}_macro_sha256",
            f"{prefix}_config_path",f"{prefix}_config_sha256","unified_macro_path","unified_macro_sha256",
            "calo_calib_path","calo_calib_sha256"):
    if not isinstance(runtime.get(key),str) or not runtime[key]: raise SystemExit(f"runtime field absent: {key}")
checks=(
 (runtime[f"{prefix}_wrapper_path"],runtime[f"{prefix}_wrapper_sha256"]),
 (runtime[f"{prefix}_macro_path"],runtime[f"{prefix}_macro_sha256"]),
 (runtime[f"{prefix}_config_path"],runtime[f"{prefix}_config_sha256"]),
 (runtime["unified_macro_path"],runtime["unified_macro_sha256"]),
 (runtime["calo_calib_path"],runtime["calo_calib_sha256"]),
 (runtime["pp_library_path"],runtime["pp_library_sha256"]),
 (runtime["auau_library_path"],runtime["auau_library_sha256"]),
 (runtime["dependency_manifest_path"],runtime["dependency_manifest_sha256"]),
 (runtime["runtime_seal_receipt_path"],runtime["runtime_seal_receipt_sha256"]),
)
for raw,digest in checks:
    p=Path(raw)
    if not p.is_absolute() or not p.is_file() or p.is_symlink(): raise SystemExit(f"bound runtime file invalid: {raw}")
    h=hashlib.sha256()
    with p.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""): h.update(block)
    if h.hexdigest()!=digest: raise SystemExit(f"bound runtime file differs: {raw}")
for key in ("durable_runtime_root","dependency_runtime_root","pp_install_prefix","auau_install_prefix"):
    p=Path(runtime[key])
    if not p.is_absolute() or not p.is_dir() or p.is_symlink() or not str(p).startswith("/sphenix/tg/tg01/bulk/"):
        raise SystemExit(f"runtime directory invalid: {key}")
base=runtime.get(f"{prefix}_replay_environment")
if not isinstance(base,dict): raise SystemExit("replay environment is absent")
env=dict(base)
env.update({
 "RJ_REPLAY_DATASET":r["dataset"], "RJ_REPLAY_LANE":r["lane"],
 "RJ_REPLAY_SAMPLE":r["sample"], "RJ_REPLAY_SI_DI_ROLE":r["si_di_role"],
 "RJ_REPLAY_FOUNDATION_CANARY":"0", "RJ_REPLAY_OWNERSHIP_STATE":payload["campaign_tag"],
})
if r["lane"] == "auau_inclusive_embedded":
    env["RJ_EMBEDDED_INCLUSIVE_JET_SAMPLE"] = r["sample"]
else:
    # Replay-environment templates are shared across lanes.  Do not let an
    # inclusive-jet selector inherited from the template leak into photon or
    # p+p rows; RJ_SIM_SAMPLE below remains the authoritative row sample.
    env.pop("RJ_EMBEDDED_INCLUSIVE_JET_SAMPLE", None)
allowed=set(runtime.get("allowed_environment_keys", []))
for key,value in env.items():
    approved = key in allowed if allowed else key.startswith("RJ_REPLAY_")
    if not approved or not re.fullmatch(r"RJ_[A-Z0-9_]+",key) or not isinstance(value,str) or any(c in value for c in "\t\r\n"):
        raise SystemExit("unsafe replay environment")
row=(system,r["dataset"],r["sample"],runtime[f"{prefix}_wrapper_path"],runtime[f"{prefix}_macro_path"],
     runtime[f"{prefix}_config_path"],runtime["pp_install_prefix"],runtime["auau_install_prefix"],
     runtime["dependency_runtime_root"],runtime["release"])
print("\t".join(row))
for key,value in sorted(env.items()): print("ENV\t"+key+"\t"+value)
PY
)" || die "production contract validation failed"

first_line="${values%%$'\n'*}"
IFS=$'\t' read -r system dataset sample wrapper macro config pp_install auau_install dependency_root release <<<"$first_line"
replay_lines=""; [[ "$values" == *$'\n'* ]] && replay_lines="${values#*$'\n'}"

export RJ_RUNTIME_BASE="$dependency_root/shared"
export RJ_SNAPSHOT_LIB_DIR="$dependency_root/shared/lib"
export RJ_SNAPSHOT_LIB_PRECEDENCE=append
if [[ "$system" == pp ]]; then
  export RJ_PP_INSTALL_PREFIX="$pp_install"; unset RJ_AUAU_INSTALL_PREFIX RJ_SHARED_HEADER_INCLUDE
else
  export RJ_AUAU_INSTALL_PREFIX="$auau_install"; unset RJ_PP_INSTALL_PREFIX
  export RJ_SHARED_HEADER_INCLUDE="$pp_install/include/caloana"
fi
export RJ_MACRO_PATH="$macro" RJ_CONFIG_YAML="$config"
export RJ_PINNED_RELEASE_NAME="$release"
export RJ_PINNED_OFFLINE_MAIN="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_ana/${release}"
export RJ_DATASET="$dataset" RJ_IS_SIM=1 RJ_SIM_SAMPLE="$sample" RJ_REPLAY_FOUNDATION_V1=1
export RJ_REQUIRE_NON_TINY_OUTPUT=1 RJ_MIN_OUTPUT_BYTES=50000 RJ_JOB_HEARTBEAT_SECONDS=0 RJ_PROFILE_JOB=0
# Fun4All_recoilJets_unified_impl.C defaults vlevel to 0 whenever
# _CONDOR_SCRATCH_DIR is set, and vlevel==0 enables ScopedSilence, which
# redirects BOTH std::cout and std::cerr to /dev/null for the whole run.
# Every LOG(0,...) FATAL is therefore discarded in batch, which is why the
# THE-236 sparse data canary failed five generations with no diagnosable
# reason. RJ_VERBOSITY=1 keeps the silence off so fatal diagnostics reach
# stdout, while LOG(10)-level chatter stays gated behind Verbosity() >= 10.
# Applied here for the T2/T3 simulation matrix; cluster 2391708 was already
# in flight from a frozen packet copy and is unaffected.
export RJ_VERBOSITY=1
# Jet constituent retention threshold. Measured 2026-08-17 on a published T1
# row (cluster 2391708 proc 673, 3,872 events, 137,191,031 bytes):
# RJJetConstituentV1 held 10,625,974 entries = 97,611,829 zip bytes = 71.5% of
# the file, ~2,744 constituent rows per event. Across the 20,006-row matrix
# that is ~200 billion rows.
#
# Nothing in THE-114's end state consumes them. capture_domain_contract.yaml
# lists 'constituent_references_where_retained' (conditional by design) and puts
# new_jet_algorithm / new_radius / new_constituent_source / new_calibration /
# new_subtraction under direct_dst_if_unretained, i.e. a fresh DST pass
# regardless. Verified directly: RJPhotonJetPairV1 stores xjgamma as its own
# branch and RJJetV1 carries corrected_pt, raw_pt, eta, phi, area, radius and
# quality_bitmask, with no constituent reference in either tree.
#
# bundle.jets.push_back(row) is unconditional in RecoilJets.cc, so every jet row
# is still written; only the constituent payload is gated, and the
# kJetConstituentPayloadRetained quality bit records the choice per jet.
# 14.9 is the ceiling of the validator's 0 <= pT < 15 range.
# RJIsolationConstituentV1 is a DIFFERENT table and is deliberately untouched;
# isolation threshold and R03/R04 replay depend on it.
export RJ_REPLAY_JET_CONSTITUENT_PT_MIN=14.9
while IFS=$'\t' read -r marker key value; do
  [[ -z "$marker" ]] && continue
  [[ "$marker" == ENV && "$key" == RJ_* ]] || die "replay environment row differs"
  export "$key=$value"
done <<<"$replay_lines"

# A production row consumes one bounded tuple bundle, not one physical DST.
# Reuse the already-verified tuple records instead of rehashing large upstream
# files.  The sealed production manifest links these two compact identities to
# the exact source lists, software, configuration, and runtime.
export RJ_REPLAY_INPUT_URI_SHA256="$chunk_sha"
export RJ_REPLAY_INPUT_FILE_SHA256="$tuple_records_sha"
export RJ_REPLAY_SOURCE_MANIFEST_SHA256="$manifest_sha"
export RJ_REPLAY_SOURCE_SHA256="$manifest_sha"
export RJ_REPLAY_PROVENANCE_MANIFEST_SHA256="$manifest_sha"
export RJ_REPLAY_SEGMENT="$original_segment_id"

scratch_output="${_CONDOR_SCRATCH_DIR}/output"
chunk_path="${_CONDOR_SCRATCH_DIR}/input_chunk.tsv"
/bin/dd if="$tuple_records_path" of="$chunk_path" iflag=skip_bytes,count_bytes \
  skip="$byte_offset" count="$byte_count" status=none
[[ -s "$chunk_path" && "$(sha_file "$chunk_path")" == "$chunk_sha" ]] || die "tuple-record slice differs"
mkdir -p "$scratch_output" "$(dirname "$output_path")" "$(dirname "$measurement_path")"
time_raw="${_CONDOR_SCRATCH_DIR}/resource.time"
set +e
# The source manifest seals a maximum of 1000 events per tuple, not an exact
# count for every physical DST.  Run to that upper bound.  The macro accepts an
# earlier stop only when all ordinary input managers independently prove EOF,
# no abort statistics are present, End() succeeds, and the summed EOF return
# code matches the exact number of ordinary managers.
/usr/bin/time -v -o "$time_raw" bash "$wrapper" "$sample" "$chunk_path" "$dataset" "$ClusterId" "$event_upper_bound" 1 NONE "$scratch_output"
wrapper_status=$?
set -e
[[ "$wrapper_status" -eq 0 ]] || die "analysis wrapper exited ${wrapper_status}"
mapfile -t roots < <(find "$scratch_output" -type f -name '*.root' -print)
[[ "${#roots[@]}" == 1 && -s "${roots[0]}" ]] || die "scientific ROOT output is absent or empty"

# Copy the completed worker output to a durable, non-published candidate name.
# A failed validator preserves this evidence under the fresh campaign namespace.
cp -- "${roots[0]}" "$candidate_path"
chmod 0644 "$candidate_path"
[[ -s "$candidate_path" ]] || die "durable scientific ROOT candidate is absent"

# Validate the authoritative completion marker and read the actual event count
# before publication.  This is a metadata-only ROOT open and does not rescan
# DST inputs or physics branches.
if ! root_contract="$({
  export RJ_VALIDATE_ROOT="$candidate_path"
  # Site setup scripts legitimately inspect optional variables such as PGHOST.
  # Keep nounset disabled only while sourcing them, then use the same PyROOT
  # metadata reader already certified for the sparse-data publication layer.
  set +e; set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n "$release" >/dev/null
  setup_status=$?
  set -u; set -e
  [[ "$setup_status" -eq 0 ]] || exit "$setup_status"
  python3 - "$candidate_path" <<'PY'
import ROOT
import sys

path = sys.argv[1]
root_file = ROOT.TFile.Open(path, "READ")
if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
    raise SystemExit(f"RJ_PRODUCTION_ROOT_CONTRACT_V1 open failure path={path}")
replay = root_file.GetDirectory("ReplayFoundationV1")
events = replay.Get("RJEventV1") if replay else None
complete = replay.Get("rj_replay_complete") if replay else None
version = replay.Get("rj_replay_schema_version") if replay else None
if (
    not replay
    or not events
    or not events.InheritsFrom("TTree")
    or not complete
    or not complete.InheritsFrom("TNamed")
    or str(complete.GetTitle()) != "1"
    or not version
    or not version.InheritsFrom("TNamed")
    or str(version.GetTitle()) != "10"
):
    raise SystemExit(
        "RJ_PRODUCTION_ROOT_CONTRACT_V1 invariant failure "
        f"replay={int(bool(replay))} events={int(bool(events))} "
        f"complete={str(complete.GetTitle()) if complete else 'MISSING'} "
        f"version={str(version.GetTitle()) if version else 'MISSING'}"
    )
actual_events = int(events.GetEntriesFast())
print(f"RJ_PRODUCTION_ROOT_CONTRACT_V1 actual_events={actual_events}")
root_file.Close()
PY
} 2>&1)"; then
  printf '%s\n' "$root_contract" >&2
  die "scientific ROOT completion contract failed; candidate preserved at ${candidate_path}"
fi
actual_events="$(printf '%s\n' "$root_contract" | /usr/bin/awk -F= '/^RJ_PRODUCTION_ROOT_CONTRACT_V1 actual_events=[0-9]+$/ {print $2}')"
[[ "$actual_events" =~ ^[1-9][0-9]*$ ]] || die "scientific ROOT event count is absent; candidate preserved at ${candidate_path}"
(( actual_events <= event_upper_bound )) || die "scientific ROOT event count exceeds bound; candidate preserved at ${candidate_path}"

python3 - "$time_raw" "$measurement_candidate" "$sample_id" "$ClusterId" "$ProcId" "$request_memory_mb" "$candidate_path" "$output_path" "$chunk_sha" "$tuple_records_sha" "$manifest_sha" "$manifest" "$event_upper_bound" "$actual_events" "$original_segment_id" <<'PY'
import hashlib,json,math,os,re,sys
from pathlib import Path
raw=Path(sys.argv[1]); receipt=Path(sys.argv[2]); candidate=Path(sys.argv[7]); output=Path(sys.argv[8])
matches=[]
for line in raw.read_text(encoding="utf-8",errors="strict").splitlines():
    m=re.fullmatch(r"\s*Maximum resident set size \(kbytes\):\s*([0-9]+)\s*",line)
    if m: matches.append(int(m.group(1)))
if len(matches)!=1 or matches[0]<=0: raise SystemExit("peak-memory field is absent")
peak=math.ceil(matches[0]/1024); request=int(sys.argv[6])
if peak>request: raise SystemExit("peak memory exceeds request")
h=hashlib.sha256()
with candidate.open("rb") as stream:
    for block in iter(lambda:stream.read(1024*1024),b""): h.update(block)
payload={"schema":"THE121_THE122_ProductionRowReceiptV2","status":"MEASURED_TERMINAL",
 "sample_id":sys.argv[3],"cluster_id":int(sys.argv[4]),"process_id":int(sys.argv[5]),
 "request_memory_mb":request,"raw_peak_memory_kb":matches[0],"peak_memory_mb":peak,
 "tuple_records_sha256":sys.argv[10],
 "chunk_sha256":sys.argv[9],"scientific_output_path":str(output),"scientific_output_sha256":h.hexdigest(),
 "scientific_output_size_bytes":candidate.stat().st_size,
 "production_manifest_sha256":sys.argv[11],"production_manifest_path":sys.argv[12],
 "expected_event_count":int(sys.argv[14]),"actual_event_count":int(sys.argv[14]),
 "event_count_upper_bound":int(sys.argv[13]),
 "event_count_contract":"POSITIVE_ACTUAL_AT_MOST_SEALED_UPPER_BOUND_V1",
 "original_segment_id":int(sys.argv[15]),
 "wrapper_exit_code":0,"publication":"ATOMIC_PART_TO_FINAL_V1","receipt_written_before_success":True}
tmp=receipt.with_name(receipt.name+f".tmp.{os.getpid()}")
tmp.write_text(json.dumps(payload,sort_keys=True,separators=(",",":"))+"\n",encoding="utf-8")
os.chmod(tmp,0o644); os.replace(tmp,receipt)
PY

# Publish only after both candidate validation and the terminal measurement
# receipt are complete. Each rename remains within its existing filesystem.
/bin/mv -- "$candidate_path" "$output_path"
/bin/mv -- "$measurement_candidate" "$measurement_path"

printf '[THE121-122-PRODUCTION-WORKER] PASS sample=%s output=%s\n' "$sample_id" "$output_path"
