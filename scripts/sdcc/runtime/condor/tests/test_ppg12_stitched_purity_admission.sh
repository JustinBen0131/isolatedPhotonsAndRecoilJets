#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
submitter="${repo_root}/scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"
tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/ppg12-closure-admission.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

err() { printf 'ERROR: %s\n' "$*" >&2; }
say() { printf '%s\n' "$*"; }
env_truthy() {
  case "${1:-}" in 1|true|TRUE|yes|YES|on|ON) return 0 ;; *) return 1 ;; esac
}
auto_merge_enabled() {
  case "${RJ_AUTO_MERGE:-1}" in 0|false|FALSE|no|NO|off|OFF) return 1 ;; esac
  return 0
}
sim_yaml_master_path() { printf '%s\n' "$RJ_CONFIG_YAML"; }

# Exercise the production functions without running the top-level dispatcher.
eval "$(sed -n '/^PPG12_STITCHED_PURITY_CONTRACT_SHA256=/,/^# Initializes paths for isSim mode/p' "$submitter" | sed '$d')"
contract_path="${repo_root}/agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
actual_contract_sha="$(python3 - "$contract_path" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

payload = json.loads(Path(sys.argv[1]).read_text())
encoded = json.dumps(
    payload, sort_keys=True, separators=(",", ":"), allow_nan=False
).encode()
print(hashlib.sha256(encoded).hexdigest())
PY
)"
[[ "$PPG12_STITCHED_PURITY_CONTRACT_SHA256" == "$actual_contract_sha" ]] || {
  printf 'FAIL submitter closure-contract hash is stale\n' >&2
  exit 1
}

source_list="${tmpdir}/source.list"
yaml="${tmpdir}/analysis_config.yaml"
admission="${tmpdir}/admission_manifest.json"
reference_manifest="${tmpdir}/reference_manifest.json"
candidate_manifest="${tmpdir}/candidate_manifest.json"
merge_audit="${tmpdir}/merge_audit.json"
gate_report="${tmpdir}/gate_report.json"
printf '/a\t/b\t/c\t/d\t/e\n' > "$source_list"
printf 'coneR: [0.3]\n' > "$yaml"
list_sha="$(ppg12_sha256_file "$source_list")"
yaml_sha="$(ppg12_sha256_file "$yaml")"
hash_a="$(printf a | shasum -a 256 | awk '{print $1}')"
runtime_dir="${tmpdir}/runtime"
runtime_manifest="${runtime_dir}/runtime_manifest.json"
mkdir -p "$runtime_dir"

python3 - "$contract_path" "$runtime_dir" "$runtime_manifest" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

contract = json.loads(Path(sys.argv[1]).read_text())
runtime_dir = Path(sys.argv[2])
manifest_path = Path(sys.argv[3])
receipt = runtime_dir / "build_receipt.json"
receipt.write_text(json.dumps({"status": "source-locked-test-runtime"}, sort_keys=True))
roles = sorted({
    role
    for group in contract["lane_runtime_contract"]["required_roles_by_lane_field"].values()
    for role in group
    if role != "lane_config"
})
files = []
for role in roles:
    path = runtime_dir / role.replace("/", "_")
    path.write_text(f"executed bytes for {role}\n")
    files.append({
        "role": role,
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    })
manifest = {
    **contract["lane_runtime_contract"]["required_manifest_values"],
    "schema_version": contract["lane_runtime_contract"]["manifest_schema_version"],
    "build_receipt": str(receipt.resolve()),
    "build_receipt_sha256": hashlib.sha256(receipt.read_bytes()).hexdigest(),
    "files": files,
}
manifest_path.write_text(json.dumps(manifest, sort_keys=True))
PY

python3 - "$admission" "$reference_manifest" "$candidate_manifest" "$merge_audit" "$gate_report" "$source_list" "$yaml" "$list_sha" "$yaml_sha" "$hash_a" "$PPG12_STITCHED_PURITY_CONTRACT_SHA256" "$contract_path" "$runtime_manifest" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

(
    out, reference_path, candidate_path, merge_path, gate_report_path,
    source_path, config_path, source_sha, config_sha, frozen_sha, contract_sha,
    contract_path, runtime_manifest_path,
) = sys.argv[1:]

def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def payload_sha256(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()

def write(path, payload):
    Path(path).write_text(json.dumps(payload, sort_keys=True))

evidence_root = Path(out).parent
frozen_evidence = evidence_root / "frozen.evidence"
frozen_evidence.write_bytes(b"a")
frozen_link = {
    "path": str(frozen_evidence.resolve()),
    "sha256": file_sha256(frozen_evidence),
}
source_link = {"path": str(Path(source_path).resolve()), "sha256": source_sha}
config_link = {"path": str(Path(config_path).resolve()), "sha256": config_sha}

contract = json.loads(Path(contract_path).read_text())
runtime_manifest_file = Path(runtime_manifest_path).resolve()
runtime_manifest = json.loads(runtime_manifest_file.read_text())
runtime_by_role = {row["role"]: row for row in runtime_manifest["files"]}
runtime_by_role["lane_config"] = {"role": "lane_config", **config_link}
runtime_source_sets = {}
runtime_hashes = {}
for field, roles in sorted(
    contract["lane_runtime_contract"]["required_roles_by_lane_field"].items()
):
    files = [runtime_by_role[role] for role in sorted(roles)]
    portable = [{"role": row["role"], "sha256": row["sha256"]} for row in files]
    set_sha = payload_sha256(portable)
    runtime_source_sets[field] = {
        "roles": sorted(roles),
        "files": files,
        "portable_role_set_sha256": set_sha,
    }
    runtime_hashes[field] = set_sha
runtime_link = {
    "path": str(runtime_manifest_file),
    "sha256": file_sha256(runtime_manifest_file),
}
receipt_path = Path(runtime_manifest["build_receipt"]).resolve()
receipt_link = {
    "path": str(receipt_path),
    "sha256": file_sha256(receipt_path),
}
runtime_evidence = {
    "manifest": runtime_link,
    "build_receipt": receipt_link,
    "manifest_payload_sha256": payload_sha256(runtime_manifest),
}

provenance_fields = (
    "implementation_sha256", "source_set_sha256", "config_contract_sha256",
    "model_set_sha256", "reconstruction_contract_sha256",
    "ownership_contract_sha256", "weights_contract_sha256",
    "estimator_contract_sha256",
)
provenance_files = [frozen_link]
provenance_set_sha = payload_sha256(provenance_files)
provenance = {field: provenance_set_sha for field in provenance_fields}
provenance_payload = {
    "schema": "ppg12-stitched-purity-provenance/v1",
    "provenance": provenance,
    "source_sets": {
        field: {"files": provenance_files, "set_sha256": provenance_set_sha}
        for field in provenance_fields
    },
}
provenance_path = evidence_root / "provenance.json"
write(provenance_path, provenance_payload)

lanes = {}
lane_rows = []
lane_links = []
for family, samples in (("inclusive", ("jet8", "jet12", "jet20", "jet30", "jet40")),
                        ("photon", ("photon5", "photon10", "photon20"))):
    for sample in samples:
        for period in ("0mrad", "1p5mrad"):
            for interaction in ("si", "di"):
                lane_id = f"{family}:{sample}:{period}:{interaction}"
                frozen_lane = {
                    "source_list_sha256": source_sha,
                    **runtime_hashes,
                    "external_scale": 1.0,
                }
                lanes[lane_id] = frozen_lane
                lane_row = {
                    "lane_id": lane_id,
                    "family": family,
                    "sample": sample,
                    "period": period,
                    "interaction": interaction,
                    "event_set_sha256": source_sha,
                    **frozen_lane,
                }
                lane_rows.append(lane_row)
                raw_lane = {
                    **lane_row,
                    "extraction": {
                        "evidence": {
                            "source_list": source_link,
                            "event_set": source_link,
                            "config": config_link,
                        },
                        "runtime_evidence": runtime_evidence,
                        "runtime_source_sets": runtime_source_sets,
                    },
                }
                raw_lane_path = evidence_root / f"lane_{len(lane_links):02d}.json"
                write(raw_lane_path, raw_lane)
                lane_links.append({
                    "lane_id": lane_id,
                    "path": str(raw_lane_path.resolve()),
                    "sha256": file_sha256(raw_lane_path),
                })
base_manifest = {
    "schema": "ppg12-stitched-purity-manifest/v1",
    "abcd_population": "unsuffixed",
    "external_scale": 1.0,
    "random_seed": 42,
    "toy_count": 20000,
    "provenance": provenance,
    "lanes": lane_rows,
}
write(reference_path, {**base_manifest, "role": "reference"})
write(candidate_path, {
    **base_manifest,
    "role": "candidate",
    "assembly": {
        "provenance": {
            "path": str(provenance_path.resolve()),
            "sha256": file_sha256(provenance_path),
        },
        "lanes": lane_links,
    },
})

merge_payload = {
    "schema": "ppg12-stitched-purity-merge-audit/v1",
    "candidate_manifest": {
        "path": str(Path(candidate_path).resolve()),
        "sha256": file_sha256(candidate_path),
    },
    "audits": [
        {
            "family": "inclusive", "status": "PASS", "failures": [],
            "max_content_delta": 0.0, "max_sumw2_delta": 0.0,
            "inputs_fixed_order": [f"inclusive-{index}" for index in range(20)],
        },
        {
            "family": "photon", "status": "PASS", "failures": [],
            "max_content_delta": 0.0, "max_sumw2_delta": 0.0,
            "inputs_fixed_order": [f"photon-{index}" for index in range(12)],
        },
    ],
}
write(merge_path, merge_payload)

gate_payload = {
    "schema": "ppg12-stitched-purity-gate-report/v1",
    "mode": "admit",
    "status": "PASS",
    "contract": {"path": "contract.json", "sha256": contract_sha},
    "reference_manifest": str(Path(reference_path).resolve()),
    "candidate_manifest": str(Path(candidate_path).resolve()),
    "failure_count": 0,
    "failures": [],
    "summary": {"expected_lane_count": 32},
}
write(gate_report_path, gate_payload)

frozen = {
    "provenance": provenance,
    "lanes": lanes,
    "abcd_population": "unsuffixed",
    "external_scale": 1.0,
    "random_seed": 42,
    "toy_count": 20000,
}
payload = {
    "schema": "ppg12-stitched-purity-admission/v1",
    "status": "PASS",
    "contract_sha256": contract_sha,
    "reference_manifest": {
        "path": str(Path(reference_path).resolve()),
        "sha256": file_sha256(reference_path),
    },
    "candidate_manifest": {
        "path": str(Path(candidate_path).resolve()),
        "sha256": file_sha256(candidate_path),
    },
    "merge_audit": {
        "path": str(Path(merge_path).resolve()),
        "sha256": file_sha256(merge_path),
    },
    "gate_report_payload_sha256": payload_sha256(gate_payload),
    "reference_frozen": frozen,
    "candidate_frozen": frozen,
}
write(out, payload)
PY

cp "$admission" "${admission}.valid"
cp "$candidate_manifest" "${candidate_manifest}.valid"
cp "$merge_audit" "${merge_audit}.valid"
cp "$gate_report" "${gate_report}.valid"

clear_gate_env() {
  unset RJ_PPG12_CLOSURE_CANARY RJ_PPG12_CLOSURE_CANARY_ID
  unset RJ_PPG12_DIRECT_SMEAR_REPAIR RJ_PPG12_DIRECT_REPAIR_LIBRECOILJETS_PATH
  unset RJ_PPG12_PHOTON_YIELD_ET_SMEAR RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4
  unset RJ_PPG12_PPSIM_G4_ONLY RJ_PPG12_PHOTON_YIELD
  unset RJ_PPG12_CLOSURE_ADMISSION_MANIFEST RJ_PPG12_CLOSURE_ADMISSION_SHA256
  unset RJ_PPG12_CLOSURE_RUNTIME_MANIFEST
  unset RJ_PPG12_CLOSURE_IMPLEMENTATION_SHA256 RJ_PPG12_CLOSURE_SOURCE_SET_SHA256
  unset RJ_PPG12_CLOSURE_CONFIG_CONTRACT_SHA256 RJ_PPG12_CLOSURE_MODEL_SET_SHA256
  unset RJ_PPG12_CLOSURE_RECONSTRUCTION_CONTRACT_SHA256
  unset RJ_PPG12_CLOSURE_OWNERSHIP_CONTRACT_SHA256 RJ_PPG12_CLOSURE_WEIGHTS_CONTRACT_SHA256
  unset RJ_PPG12_CLOSURE_ESTIMATOR_CONTRACT_SHA256
  unset RJ_PPG12_PERIOD RJ_SUBMISSION_NAMESPACE RJ_DEST_BASE_OVERRIDE
}

set_common() {
  ACTION=condorDoAll
  SIM_SAMPLE=run28_jet8
  SIM_CLEAN_LIST="$source_list"
  RJ_CONFIG_YAML="$yaml"
  GROUP_SIZE=5
  GROUP_SIZE_EXPLICIT=1
  MAX_JOBS=1
  MAX_JOBS_EXPLICIT=1
  RJ_AUTO_MERGE=0
}

expect_pass() {
  local label="$1"
  if ! validate_ppg12_stitched_purity_admission; then
    printf 'FAIL expected pass: %s\n' "$label" >&2
    exit 1
  fi
}

expect_fail() {
  local label="$1"
  if validate_ppg12_stitched_purity_admission; then
    printf 'FAIL expected rejection: %s\n' "$label" >&2
    exit 1
  fi
}

clear_gate_env
set_common
SIM_SAMPLE=run28_jet5
expect_pass 'Jet5 is outside the 32-lane closure contract'

clear_gate_env
set_common
RJ_PPG12_CLOSURE_CANARY=1
RJ_PPG12_CLOSURE_CANARY_ID=unit-canary
RJ_SUBMISSION_NAMESPACE=unit-canary
RJ_DEST_BASE_OVERRIDE="${tmpdir}/canary-output"
expect_pass 'bounded five-file canary'

GROUP_SIZE=4
expect_fail 'undersized canary group is forbidden'
GROUP_SIZE=5
RJ_AUTO_MERGE=1
expect_fail 'canary automatic merge is forbidden'

clear_gate_env
set_common
ACTION=condorDoAllDirect
RJ_PPG12_CLOSURE_CANARY=1
RJ_PPG12_CLOSURE_CANARY_ID=unit-canary-direct
RJ_SUBMISSION_NAMESPACE=unit-canary-direct
RJ_DEST_BASE_OVERRIDE="${tmpdir}/canary-output-direct"
expect_pass 'bounded direct canary is guarded'
MAX_JOBS=2
expect_fail 'direct canary maxJobs bypass is forbidden'

clear_gate_env
set_common
expect_fail 'broad production without admission is blocked'

# The one explicitly authorized 2026-07-21 direct-repair campaign may bypass
# admission only when its isolated namespace, deployed graph, and every source
# and runtime hash are exact. Override the production pins with deterministic
# test fixtures so this remains a hermetic unit test.
direct_base="${tmpdir}/direct-base"
direct_runtime="${tmpdir}/direct-runtime/libRecoilJets.so"
mkdir -p "${direct_base}/src" "${direct_base}/macros" "$(dirname "$direct_runtime")"
printf 'recoil cc\n' > "${direct_base}/src/RecoilJets.cc"
printf 'recoil h\n' > "${direct_base}/src/RecoilJets.h"
printf 'macro\n' > "${direct_base}/macros/Fun4All_recoilJets_unified_impl.C"
printf 'runtime\n' > "$direct_runtime"
BASE="$direct_base"
PPG12_DIRECT_REPAIR_EXPECTED_RECOIL_CC_SHA256="$(ppg12_sha256_file "${direct_base}/src/RecoilJets.cc")"
PPG12_DIRECT_REPAIR_EXPECTED_RECOIL_H_SHA256="$(ppg12_sha256_file "${direct_base}/src/RecoilJets.h")"
PPG12_DIRECT_REPAIR_EXPECTED_MACRO_SHA256="$(ppg12_sha256_file "${direct_base}/macros/Fun4All_recoilJets_unified_impl.C")"
PPG12_DIRECT_REPAIR_EXPECTED_LIB_SHA256="$(ppg12_sha256_file "$direct_runtime")"

clear_gate_env
set_common
MAX_JOBS=0
MAX_JOBS_EXPLICIT=0
RJ_PPG12_DIRECT_SMEAR_REPAIR=1
RJ_SUBMISSION_NAMESPACE=the97_ppg12_deployed_smear_full_20260721_1430_v1
RJ_DEST_BASE_OVERRIDE="${tmpdir}/${RJ_SUBMISSION_NAMESPACE}/outputs"
RJ_PPG12_PHOTON_YIELD_ET_SMEAR=1
RJ_PPG12_PPSIM_REBUILD_CALO_FROM_G4=1
RJ_PPG12_PPSIM_G4_ONLY=1
RJ_PPG12_PHOTON_YIELD=1
RJ_PPG12_DIRECT_REPAIR_LIBRECOILJETS_PATH="$direct_runtime"
expect_pass 'exact hash-pinned direct deployed-smear repair permits full production'

unset RJ_PPG12_PHOTON_YIELD_ET_SMEAR
expect_fail 'direct repair without explicit ET smearing is blocked'
RJ_PPG12_PHOTON_YIELD_ET_SMEAR=1

printf 'drift\n' >> "${direct_base}/src/RecoilJets.cc"
expect_fail 'direct repair source hash drift is blocked'
printf 'recoil cc\n' > "${direct_base}/src/RecoilJets.cc"

RJ_SUBMISSION_NAMESPACE=unscoped-repair
expect_fail 'direct repair outside its canonical namespace is blocked'

clear_gate_env
set_common
RJ_PPG12_CLOSURE_ADMISSION_MANIFEST="$admission"
RJ_PPG12_CLOSURE_ADMISSION_SHA256="$(ppg12_sha256_file "$admission")"
RJ_PPG12_CLOSURE_RUNTIME_MANIFEST="$runtime_manifest"
RJ_PPG12_PERIOD=0mrad
export RJ_PPG12_CLOSURE_IMPLEMENTATION_SHA256="$hash_a"
export RJ_PPG12_CLOSURE_SOURCE_SET_SHA256="$hash_a"
export RJ_PPG12_CLOSURE_CONFIG_CONTRACT_SHA256="$hash_a"
export RJ_PPG12_CLOSURE_MODEL_SET_SHA256="$hash_a"
export RJ_PPG12_CLOSURE_RECONSTRUCTION_CONTRACT_SHA256="$hash_a"
export RJ_PPG12_CLOSURE_OWNERSHIP_CONTRACT_SHA256="$hash_a"
export RJ_PPG12_CLOSURE_WEIGHTS_CONTRACT_SHA256="$hash_a"
export RJ_PPG12_CLOSURE_ESTIMATOR_CONTRACT_SHA256="$hash_a"
expect_pass 'exact hash-bound admission permits the admitted lane'

python3 - "$admission" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload.pop("candidate_manifest")
path.write_text(json.dumps(payload, sort_keys=True))
PY
RJ_PPG12_CLOSURE_ADMISSION_SHA256="$(ppg12_sha256_file "$admission")"
expect_fail 'pass-shaped admission without linked artifacts is rejected'
cp "${admission}.valid" "$admission"
RJ_PPG12_CLOSURE_ADMISSION_SHA256="$(ppg12_sha256_file "$admission")"

printf '\n' >> "$candidate_manifest"
expect_fail 'mutated linked candidate manifest is rejected'
cp "${candidate_manifest}.valid" "$candidate_manifest"

printf '\n' >> "$merge_audit"
expect_fail 'mutated linked merge audit is rejected'
cp "${merge_audit}.valid" "$merge_audit"

python3 - "$gate_report" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
payload = json.loads(path.read_text())
payload["status"] = "FAIL"
path.write_text(json.dumps(payload, sort_keys=True))
PY
expect_fail 'mutated linked gate report is rejected'
cp "${gate_report}.valid" "$gate_report"

RJ_PPG12_CLOSURE_ADMISSION_SHA256="$(printf stale | shasum -a 256 | awk '{print $1}')"
expect_fail 'changed admission bytes are rejected'
RJ_PPG12_CLOSURE_ADMISSION_SHA256="$(ppg12_sha256_file "$admission")"

printf '/changed\t/b\t/c\t/d\t/e\n' > "$source_list"
expect_fail 'changed source list is rejected'
printf '/a\t/b\t/c\t/d\t/e\n' > "$source_list"

export RJ_PPG12_CLOSURE_MODEL_SET_SHA256="$(printf changed | shasum -a 256 | awk '{print $1}')"
expect_pass 'caller-supplied model hash cannot replace computed evidence hashes'
runtime_model="$(python3 - "$runtime_manifest" <<'PY'
import json
import sys
from pathlib import Path
payload = json.loads(Path(sys.argv[1]).read_text())
print(next(row["path"] for row in payload["files"] if row["role"] == "ppg_apply_model_base_v3E"))
PY
)"
cp "$runtime_model" "${runtime_model}.valid"
printf 'mutated executed model\n' > "$runtime_model"
expect_fail 'mutated executed model is rejected even with caller-supplied replacement hash'
mv "${runtime_model}.valid" "$runtime_model"

printf changed > "${tmpdir}/frozen.evidence"
expect_fail 'changed linked model/provenance bytes are rejected'

printf 'PASS ppg12_stitched_purity_admission\n'
