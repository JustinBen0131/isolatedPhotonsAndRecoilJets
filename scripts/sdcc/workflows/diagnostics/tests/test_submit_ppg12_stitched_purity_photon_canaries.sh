#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
driver="$(cd "${script_dir}/.." && pwd -P)/submit_ppg12_stitched_purity_photon_canaries.sh"
tmp="$(mktemp -d "${TMPDIR:-/tmp}/ppg12-paired-12-test.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT
fail() { printf 'FAIL: %s\n' "$*" >&2; exit 1; }

mkdir -p "$tmp/assets" "$tmp/lane"
for name in setup apply_bdt apply_config base_e base_v3e npb mask runtime_manifest recoil_config; do
  printf 'asset %s\n' "$name" > "$tmp/assets/$name"
done
for photon in 5 10 20; do
  for interaction in si di; do
    slug="photon${photon}_${interaction}"
    if [[ "$interaction" == si ]]; then
      interaction_upper=SI
    else
      interaction_upper=DI
    fi
    {
      printf 'INPUTREADHITS::listfile[0] = "g4.list";\n'
      printf 'INPUTREADHITS::listfile[4] = "truth.list";\n'
      if [[ "$interaction" == di ]]; then
        printf 'TruthJetInput *truth = new TruthJetInput(Jet::PARTICLE);\n'
        printf 'truth->add_embedding_flag(2);\n'
      fi
    } > "$tmp/lane/${slug}.C"
    : > "$tmp/lane/${slug}_g4.list"
    : > "$tmp/lane/${slug}_truth.list"
    for index in 0 1 2 3 4; do
      identity="pythia8_PhotonJet${photon}_${interaction_upper}_${index}.root"
      printf '/source/G4Hits_%s\n' "$identity" >> "$tmp/lane/${slug}_g4.list"
      printf '/source/DST_TRUTH_JET_%s\n' "$identity" >> "$tmp/lane/${slug}_truth.list"
    done
  done
done

source_manifest="$tmp/source.json"
python3 - "$source_manifest" "$tmp" <<'PY'
import json, sys
from pathlib import Path
out, base = Path(sys.argv[1]), Path(sys.argv[2])
lanes=[]
for photon in (5,10,20):
  for period in ("0mrad","1p5mrad"):
    for interaction in ("si","di"):
      slug=f"photon{photon}_{interaction}"
      lanes.append({
        "lane_id":f"photon:photon{photon}:{period}:{interaction}",
        "sample":f"Photon{photon}","period":period,"interaction":interaction.upper(),
        "ppg_macro":str(base/"lane"/(slug+".C")),
        "g4_full_list":str(base/"lane"/(slug+"_g4.list")),
        "truthjet_full_list":str(base/"lane"/(slug+"_truth.list")),
      })
doc={"schema":"ppg12-paired-source-manifest/v1","common":{
  "setup_script":str(base/"assets/setup"),"apply_bdt":str(base/"assets/apply_bdt"),
  "apply_config":str(base/"assets/apply_config"),"base_e_model":str(base/"assets/base_e"),
  "base_v3e_model":str(base/"assets/base_v3e"),"npb_model":str(base/"assets/npb"),
  "tower_mask":str(base/"assets/mask"),"recoil_runtime_manifest":str(base/"assets/runtime_manifest"),
  "recoil_config":str(base/"assets/recoil_config")},"lanes":lanes}
out.write_text(json.dumps(doc,indent=2)+"\n")
PY

mock="$tmp/mock_paired_driver.sh"
cat > "$mock" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
mode=plan token="" lane="" sample="" period="" interaction="" output=""
while (($#)); do
  case "$1" in
    --run) mode=run; shift;; --token) token="$2"; shift 2;;
    --lane-id) lane="$2"; shift 2;; --sample) sample="$2"; shift 2;;
    --period) period="$2"; shift 2;; --interaction) interaction="$2"; shift 2;;
    --output-dir) output="$2"; shift 2;; *) shift 2;;
  esac
done
digest="$(printf '%s' "$lane|$sample|$period|$interaction|$output" | shasum -a 256 | awk '{print $1}')"
expected="ppg12-oracle:$digest"
if [[ "$mode" == plan ]]; then
  printf '  lane: %s | %s | %s\n  lane_id: %s\n  run_token: %s\n' "$sample" "$period" "$interaction" "$lane" "$expected"
  exit 0
fi
[[ "$token" == "$expected" ]]
mkdir -p "$output/comparison"
contract="$output/paired_oracle_contract.json"
printf '{"schema_version":3,"lane":{"lane_id":"%s","sample":"%s","period":"%s","interaction":"%s","rows":5}}\n' "$lane" "$sample" "$period" "$interaction" > "$contract"
contract_hash="$(shasum -a 256 "$contract" | awk '{print $1}')"
candidate="$output/comparison/paired_oracle_candidates.csv"
printf 'lane_id,runtime_contract_sha256,candidate_identity,match_status,ppg12_tag_evidence_source\n%s,%s,%s,matched,preserved_ppg12_executable\n' "$lane" "$contract_hash" "$lane:candidate" > "$candidate"
candidate_hash="$(shasum -a 256 "$candidate" | awk '{print $1}')"
printf '{"schema_version":1,"evidence_source":"preserved_ppg12_executable_aggregate","status":"PASS","mode":"full","lane_identity":{"lane_id":"%s","runtime_contract_sha256":"%s"},"provenance":{"runtime_contract":{"path":"%s","sha256":"%s"},"candidate_csv":{"path":"%s","sha256":"%s"}}}\n' "$lane" "$contract_hash" "$contract" "$contract_hash" "$candidate" "$candidate_hash" > "$output/comparison/executable_aggregate.json"
printf 'PASS\n' > "$output/RUN_STATE"
EOF
chmod +x "$mock"

campaign=paired12_unit
evidence="$tmp/evidence"
output="$tmp/output"
env_args=(
  "RJ_PPG12_PAIRED_SOURCE_MANIFEST=$source_manifest"
  "RJ_PPG12_PAIRED_ORACLE_DRIVER=$mock"
  "RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG=$campaign"
  "RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR=$evidence"
  "RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT=$output"
)

env "${env_args[@]}" "$driver" > "$tmp/plan.out"
plan="$evidence/photon_canary_plan.json"
python3 - "$plan" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
assert d["schema"]=="ppg12-stitched-purity-photon-canary-plan/v5"
assert d["bounded_contract"]["execution"]=="foreground_source_locked_paired_executable"
assert d["bounded_contract"]["python_shadow_admissible"] is False
assert d["bounded_contract"]["archived_ppg12_root_admissible"] is False
assert len(d["lanes"])==12 and len({x["lane_id"] for x in d["lanes"]})==12
assert len({x["output_base"] for x in d["lanes"]})==12
assert all(x["execution"]=="source_locked_paired_executable" for x in d["lanes"])
assert all("ppg12_oracle_root" not in x for x in d["lanes"])
assert d["source_validator"]["path"].endswith("ppg12_paired_oracle_condor_fanout.py")
PY
token="$(sed -n 's/^run_token=//p' "$tmp/plan.out")"
[[ "$token" =~ ^RUN_PPG12_PAIRED_[0-9a-f]{64}$ ]] || fail "invalid plan token"

# Raw cardinality must be checked before lane IDs can be collapsed into a
# dictionary: a 13th duplicate row is never an admissible 12-lane manifest.
duplicate_manifest="$tmp/source_13_row_duplicate.json"
python3 - "$source_manifest" "$duplicate_manifest" <<'PY'
import json, sys
source, target = sys.argv[1:]
doc = json.load(open(source))
doc["lanes"].append(dict(doc["lanes"][0]))
with open(target, "w") as stream:
    json.dump(doc, stream, indent=2, sort_keys=True)
    stream.write("\n")
PY
if env \
  "RJ_PPG12_PAIRED_SOURCE_MANIFEST=$duplicate_manifest" \
  "RJ_PPG12_PAIRED_ORACLE_DRIVER=$mock" \
  "RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG=${campaign}_duplicate" \
  "RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR=$tmp/evidence_duplicate" \
  "RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT=$tmp/output_duplicate" \
  "$driver" >"$tmp/duplicate.out" 2>"$tmp/duplicate.err"; then
  fail "13-row duplicate source manifest was accepted"
fi
grep -q 'exactly 12 raw canonical lane rows; observed=13' "$tmp/duplicate.err" || \
  fail "13-row rejection did not report raw cardinality"
[[ ! -s "$tmp/evidence_duplicate/photon_canary_plan.json" ]] || \
  fail "13-row duplicate emitted a foreground plan"

# With cardinality held at twelve, duplicate IDs must still fail before
# mapping construction rather than silently replacing a missing lane.
duplicate_id_manifest="$tmp/source_12_row_duplicate_id.json"
python3 - "$source_manifest" "$duplicate_id_manifest" <<'PY'
import json, sys
source, target = sys.argv[1:]
doc = json.load(open(source))
doc["lanes"][-1] = dict(doc["lanes"][0])
with open(target, "w") as stream:
    json.dump(doc, stream, indent=2, sort_keys=True)
    stream.write("\n")
PY
if env \
  "RJ_PPG12_PAIRED_SOURCE_MANIFEST=$duplicate_id_manifest" \
  "RJ_PPG12_PAIRED_ORACLE_DRIVER=$mock" \
  "RJ_PPG12_PHOTON_CANARY_CAMPAIGN_TAG=${campaign}_duplicate_id" \
  "RJ_PPG12_PHOTON_CANARY_EVIDENCE_DIR=$tmp/evidence_duplicate_id" \
  "RJ_PPG12_PHOTON_CANARY_OUTPUT_ROOT=$tmp/output_duplicate_id" \
  "$driver" >"$tmp/duplicate_id.out" 2>"$tmp/duplicate_id.err"; then
  fail "12-row duplicate lane ID source manifest was accepted"
fi
grep -q 'duplicates raw lane_id values' "$tmp/duplicate_id.err" || \
  fail "duplicate lane-ID rejection did not precede mapping"

if env "${env_args[@]}" "$driver" --submit --token WRONG >"$tmp/wrong" 2>&1; then
  fail "wrong token executed"
fi
[[ ! -e "$output" ]] || fail "wrong token mutated output"

env "${env_args[@]}" "$driver" --submit --token "$token" > "$tmp/run.out"
receipt="$evidence/paired_execution_receipts.tsv"
python3 - "$receipt" <<'PY'
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1]),delimiter="\t"))
assert len(rows)==12 and len({x["lane_id"] for x in rows})==12
for key in ("output_base","runtime_contract","runtime_contract_sha256","candidate_csv","candidate_csv_sha256","executable_aggregate","executable_aggregate_sha256"):
  assert len({x[key] for x in rows})==12
assert all(x["status"]=="PASS" for x in rows)
PY

rm -rf "$output"
python3 - "$tmp/lane/photon5_si_g4.list" "$tmp/lane/photon5_si_truth.list" <<'PY'
import sys
from pathlib import Path
g4, truth = map(Path, sys.argv[1:])
g4_rows = g4.read_text().splitlines()
truth_rows = truth.read_text().splitlines()
g4_rows[-1] = "/source/G4Hits_pythia8_PhotonJet5_SI_CHANGED.root"
truth_rows[-1] = "/source/DST_TRUTH_JET_pythia8_PhotonJet5_SI_CHANGED.root"
g4.write_text("\n".join(g4_rows) + "\n")
truth.write_text("\n".join(truth_rows) + "\n")
PY
env "${env_args[@]}" "$driver" > "$tmp/changed-plan.out"
changed="$(sed -n 's/^run_token=//p' "$tmp/changed-plan.out")"
[[ "$changed" != "$token" ]] || fail "source mutation did not invalidate token"

printf 'PPG12_PAIRED_12_TEST_PASS\n'
