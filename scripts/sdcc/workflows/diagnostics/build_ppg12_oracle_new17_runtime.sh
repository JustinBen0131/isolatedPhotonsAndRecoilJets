#!/usr/bin/env bash
# Build a sealed new.17 RecoilJets runtime for the paired PPG12 oracle.
#
# Default mode is a read-only plan.  The build is foreground-only, writes to a
# new contained output root, and requires the exact token printed by plan mode.
# It never submits Condor, installs into a user prefix, or mutates the checkout.

set -euo pipefail

die() {
  printf 'PPG12_ORACLE_NEW17_BUILD_FAIL: %s\n' "$*" >&2
  exit 2
}

usage() {
  cat <<'EOF'
Usage:
  build_ppg12_oracle_new17_runtime.sh --output-dir ABS [--jobs N] \
    [--photon-source-dir ABS] [--ppg-repo ABS] [--ppg-revision SHA] \
    [--ppg-source-runtime-manifest ABS] \
    [--estimator-revision SHA] [--truth-vertex-reweight-0mrad ABS] \
    [--truth-vertex-reweight-1p5mrad ABS] [--yaml-cpp-include-dir ABS] \
    [--roounfold-library ABS] [--roounfold-include-dir ABS] \
    [--roounfold-pcm ABS]
  build_ppg12_oracle_new17_runtime.sh --build --token TOKEN \
    --output-dir ABS [--jobs N] [--photon-source-dir ABS] \
    [--ppg-repo ABS] [--ppg-revision SHA] [--estimator-revision SHA] \
    [--ppg-source-runtime-manifest ABS] \
    [--truth-vertex-reweight-0mrad ABS] \
    [--truth-vertex-reweight-1p5mrad ABS] [--yaml-cpp-include-dir ABS] \
    [--roounfold-library ABS] [--roounfold-include-dir ABS] \
    [--roounfold-pcm ABS]

Default mode prints the immutable build contract and authorization token.
--build performs the foreground build only when TOKEN exactly matches that
contract.  ABS must be a new path below REPO/.recoiljets_tmp or /tmp.

--ppg-source-runtime-manifest imports the exact source-locked libCaloAna24.so
from a previously sealed new.17 runtime.  It does not permit the archived June
binary: the origin manifest and receipt must prove the same source revision and
common-runtime rebuild contract accepted by the default build mode.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "${script_dir}/../../../.." && pwd -P)"
expected_offline="/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.17"

mode=plan
provided_token=""
output_dir=""
setup_script="/opt/sphenix/core/bin/sphenix_setup.sh"
jobs=4
photon_source_dir="${repo_root}/src"
ppg_repo="${repo_root}/ppg12codeGit"
# This is the last CaloAna24 source revision before the archived June PPG12
# production.  The working tree is deliberately ignored: git-show exports the
# committed source into the sealed build root.
ppg_revision="1c0ff86bf0ebabfba63a1abc4512cbe59fe48e31"
ppg_source_runtime_manifest=""
ppg_binary_mode="source_locked_rebuild"
ppg_origin_build_receipt=""
ppg_origin_library=""
ppg_origin_manifest_sha256=""
ppg_origin_receipt_sha256=""
ppg_origin_library_sha256=""
# Reconstruction and the downstream estimator have distinct historical
# contracts.  Never source the estimator from the older reconstruction
# revision or from the mutable checkout.
estimator_revision="29f8223bd9b36dffab07961b597afa94185bbdf1"
yaml_cpp_library="/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so"
yaml_cpp_include_dir="/sphenix/u/shuhang98/install/include"
roounfold_root="/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/analysis/roounfold"
roounfold_library="${roounfold_root}/libRooUnfold.so"
roounfold_include_dir="${roounfold_root}/src"
roounfold_pcm="${roounfold_root}/tmp/linuxx8664gcc/RooUnfoldDict_rdict.pcm"
expected_roounfold_library_sha256="d135771391ae250bcb64c0889571825abe9924649485890e7a9c64648ee99062"
expected_roounfold_pcm_sha256="2d91962a7b42acf246c7a80339eee71ca2f7e6df18ef76051d24a83bc61d4244"
expected_roounfold_header_tree_sha256="ea9b923a8f6bc57b28027b7183b10e87246810c326208b1e36bf2b4b7491a458"
roounfold_headers=(
  RooUnfold.h
  RooUnfoldResponse.h
  RooUnfoldBayes.h
  RooUnfoldBinByBin.h
  RooUnfoldErrors.h
  RooUnfoldInvert.h
  RooUnfoldParms.h
  RooUnfoldSvd.h
  RooUnfoldTUnfold.h
)
roounfold_library_overridden=0
roounfold_include_dir_overridden=0
roounfold_pcm_overridden=0
vertex_scan_data_file="/sphenix/user/shuhangli/ppg12/efficiencytool/results/data_histo_bdt_nom_vtxscan.root"
mbd_correction_file="/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr"
truth_vertex_reweight_0mrad="/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/0mrad/reweight.root"
truth_vertex_reweight_1p5mrad="/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/1p5mrad/reweight.root"
apply_model_dir="/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models"
apply_npb_model="/sphenix/user/shuhangli/ppg12/FunWithxgboost/npb_models/npb_score_split_tmva.root"
apply_model_names=(base base_vr base_v0 base_v1 base_v2 base_v3 base_E base_v0E base_v1E base_v2E base_v3E)
expected_apply_bdt_sha256="bd6e7c5bc9858ddad9bc835552d818c00290bb7de3f5036f44d8bf4804734366"
expected_apply_config_sha256="b8d1bc359a647cc913f213777fc42958b532b30c37a63bb318680130eb6e321b"
expected_recoeff_sha256="e9b25fdb6dd8a6bfbbad029cb90aaddc9489fdf2846c630ea63c8c41ac771eee"
expected_recoeff_config_sha256="42b7be1628843d5b7607ab988ffb58c6d019d8ade01b3c4d528611498db95732"
expected_recoeff_period_config_0mrad_sha256="3995033c8867f4b0e21d5ebc025d36395185671da474fceec128b20db7218be2"
expected_recoeff_period_config_1p5mrad_sha256="6d2e4cc691e2fdd49271486ef193055704da00bcbd6b50ced76fcdd99cd050b8"
expected_truth_vertex_reweight_0mrad_sha256="4c2a50fa2dd4fe6e3f8b823454f19753367876edabc9fe48164447a0d07be6b9"
expected_truth_vertex_reweight_1p5mrad_sha256="1429442b2cce368bcc4b4fd613f706f27b9c194f1e835d800238e059c0804d4d"
expected_apply_model_hashes=(
  "0ec432081df6bd5cbdc68c948c4c325331bc0220834b4720429c2049c86880e1"
  "8b70b9bda2430fa7147694ebb8f51de72122ee525a8b5b59dca0d50c4c86000f"
  "d71acf56911f4648baecf17c9bd567fbfa6f646e31ce20300ea1801fe77f15d4"
  "56f59bb0d80f47726cf425423169bd6556f15353c0b2204e3124d84f9d3c579f"
  "d863f545a2d9d8557243ddea38a3a3165c83c11d756eb65e1e998e90d68c75bf"
  "3a722fb2c16f0120d62e74710963187f8d2efad9b0868ea411cff33055052393"
  "7e2d5ed9d1216ab30b92a137482a5bfa111d78e86d4f5a459179c38bc7993997"
  "6a302d7ebece4f5a592a38edb8de934ebb3a2935fdedff4d7bf012d54c7870dd"
  "b75c9e3c3c4a6e6333c79567b8de0faacf8ece6591e085813be974af6f0571e9"
  "c5a14d44b3655516692b012f15f2f84419f1acba70479d23c82421b6eeb49f27"
  "7679e634260402fb3815b2733767182690eec7587f9e09bffc307a05d00d59df"
)
expected_apply_npb_sha256="d6086dadac534013cda15cdfb69c1683776d3456d9e439903589653e8ac19eab"

while (($#)); do
  case "$1" in
    --build) mode=build; shift ;;
    --token) [[ $# -ge 2 ]] || die "--token requires a value"; provided_token="$2"; shift 2 ;;
    --output-dir) [[ $# -ge 2 ]] || die "--output-dir requires a value"; output_dir="$2"; shift 2 ;;
    --setup-script) [[ $# -ge 2 ]] || die "--setup-script requires a value"; setup_script="$2"; shift 2 ;;
    --photon-source-dir) [[ $# -ge 2 ]] || die "--photon-source-dir requires a value"; photon_source_dir="$2"; shift 2 ;;
    --ppg-repo) [[ $# -ge 2 ]] || die "--ppg-repo requires a value"; ppg_repo="$2"; shift 2 ;;
    --ppg-revision) [[ $# -ge 2 ]] || die "--ppg-revision requires a value"; ppg_revision="$2"; shift 2 ;;
    --ppg-source-runtime-manifest) [[ $# -ge 2 ]] || die "--ppg-source-runtime-manifest requires a value"; ppg_source_runtime_manifest="$2"; ppg_binary_mode="source_locked_runtime_import"; shift 2 ;;
    --estimator-revision) [[ $# -ge 2 ]] || die "--estimator-revision requires a value"; estimator_revision="$2"; shift 2 ;;
    --yaml-cpp-library) [[ $# -ge 2 ]] || die "--yaml-cpp-library requires a value"; yaml_cpp_library="$2"; shift 2 ;;
    --yaml-cpp-include-dir) [[ $# -ge 2 ]] || die "--yaml-cpp-include-dir requires a value"; yaml_cpp_include_dir="$2"; shift 2 ;;
    --roounfold-library) [[ $# -ge 2 ]] || die "--roounfold-library requires a value"; roounfold_library="$2"; roounfold_library_overridden=1; shift 2 ;;
    --roounfold-include-dir) [[ $# -ge 2 ]] || die "--roounfold-include-dir requires a value"; roounfold_include_dir="$2"; roounfold_include_dir_overridden=1; shift 2 ;;
    --roounfold-pcm) [[ $# -ge 2 ]] || die "--roounfold-pcm requires a value"; roounfold_pcm="$2"; roounfold_pcm_overridden=1; shift 2 ;;
    --vertex-scan-data-file) [[ $# -ge 2 ]] || die "--vertex-scan-data-file requires a value"; vertex_scan_data_file="$2"; shift 2 ;;
    --mbd-correction-file) [[ $# -ge 2 ]] || die "--mbd-correction-file requires a value"; mbd_correction_file="$2"; shift 2 ;;
    --truth-vertex-reweight-0mrad) [[ $# -ge 2 ]] || die "--truth-vertex-reweight-0mrad requires a value"; truth_vertex_reweight_0mrad="$2"; shift 2 ;;
    --truth-vertex-reweight-1p5mrad) [[ $# -ge 2 ]] || die "--truth-vertex-reweight-1p5mrad requires a value"; truth_vertex_reweight_1p5mrad="$2"; shift 2 ;;
    --apply-model-dir) [[ $# -ge 2 ]] || die "--apply-model-dir requires a value"; apply_model_dir="$2"; shift 2 ;;
    --apply-npb-model) [[ $# -ge 2 ]] || die "--apply-npb-model requires a value"; apply_npb_model="$2"; shift 2 ;;
    --jobs) [[ $# -ge 2 ]] || die "--jobs requires a value"; jobs="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ "$jobs" =~ ^[1-8]$ ]] || die "--jobs must be an integer from 1 through 8"
[[ -n "$output_dir" && "$output_dir" == /* ]] || die "--output-dir must be absolute"
[[ "$output_dir" != *$'\n'* && "$output_dir" != *$'\r'* && "$output_dir" != *' '* ]] || \
  die "--output-dir contains whitespace or a newline"
case "$output_dir" in
  "${repo_root}/.recoiljets_tmp/"*|/tmp/*) ;;
  *) die "--output-dir must be below ${repo_root}/.recoiljets_tmp or /tmp" ;;
esac
[[ ! -e "$output_dir" ]] || die "output path already exists: $output_dir"
[[ "$setup_script" == /* && -f "$setup_script" && -s "$setup_script" ]] || \
  die "setup script is missing or empty: $setup_script"
[[ "$photon_source_dir" == /* && -d "$photon_source_dir" ]] || \
  die "--photon-source-dir must be an existing absolute directory: $photon_source_dir"
[[ "$ppg_repo" == /* && -d "$ppg_repo/.git" ]] || \
  die "--ppg-repo must be an existing absolute Git checkout: $ppg_repo"
[[ "$ppg_revision" =~ ^[0-9a-f]{40}$ ]] || die "--ppg-revision must be a full commit SHA"
[[ "$estimator_revision" =~ ^[0-9a-f]{40}$ ]] || \
  die "--estimator-revision must be a full commit SHA"
if [[ -n "$ppg_source_runtime_manifest" ]]; then
  [[ "$ppg_source_runtime_manifest" == /* && \
     "$ppg_source_runtime_manifest" != *$'\n'* && \
     "$ppg_source_runtime_manifest" != *$'\r'* && \
     "$ppg_source_runtime_manifest" != *$'\t'* ]] || \
    die "--ppg-source-runtime-manifest must be an absolute single-line path"
  [[ -f "$ppg_source_runtime_manifest" && -s "$ppg_source_runtime_manifest" ]] || \
    die "source runtime manifest is missing or empty: $ppg_source_runtime_manifest"
fi
ppg_repo_real="$(cd "$ppg_repo" && pwd -P)"
git_safe=(-c "safe.directory=${ppg_repo_real}")
git "${git_safe[@]}" -C "$ppg_repo_real" cat-file -e "${ppg_revision}^{commit}" 2>/dev/null || \
  die "--ppg-revision is not present in --ppg-repo: $ppg_revision"
git "${git_safe[@]}" -C "$ppg_repo_real" cat-file -e \
  "${estimator_revision}^{commit}" 2>/dev/null || \
  die "--estimator-revision is not present in --ppg-repo: $estimator_revision"
ppg_source_names=(configure.ac Makefile.am autogen.sh CaloAna24.cc CaloAna24.h)
for source_name in "${ppg_source_names[@]}"; do
  git "${git_safe[@]}" -C "$ppg_repo_real" cat-file -e \
    "${ppg_revision}:anatreemaker/source/${source_name}" 2>/dev/null || \
    die "PPG12 source revision lacks anatreemaker/source/${source_name}"
done

estimator_source_names=(
  RecoEffCalculator_TTreeReader.C
  CrossSectionWeights.h
  TruthVertexReweightLoader.h
  config_bdt_nom.yaml
  config_bdt_nom_0rad.yaml
  config_bdt_nom_1p5mrad.yaml
  CalculatePhotonYield.C
)
for source_name in "${estimator_source_names[@]}"; do
  git "${git_safe[@]}" -C "$ppg_repo_real" cat-file -e \
    "${estimator_revision}:efficiencytool/${source_name}" 2>/dev/null || \
    die "estimator revision lacks efficiencytool/${source_name}"
done
apply_stage_source_names=(apply_BDT.C config_nom.yaml)
for source_name in "${apply_stage_source_names[@]}"; do
  git "${git_safe[@]}" -C "$ppg_repo_real" cat-file -e \
    "${estimator_revision}:FunWithxgboost/${source_name}" 2>/dev/null || \
    die "estimator revision lacks FunWithxgboost/${source_name}"
done

canonical_photon_cc="${photon_source_dir}/PhotonClusterBuilder.cc"
canonical_photon_h="${photon_source_dir}/PhotonClusterBuilder.h"
recoil_source="${repo_root}/src"
macro_wrapper_source="${repo_root}/macros/Fun4All_recoilJets.C"
macro_impl_source="${repo_root}/macros/Fun4All_recoilJets_unified_impl.C"
trace_instrumenter="${repo_root}/scripts/diagnostics/pp_currentian/instrument_ppg12_recoeff_trace.py"

contract_inputs=(
  "$canonical_photon_cc"
  "$canonical_photon_h"
  "${recoil_source}/configure.ac"
  "${recoil_source}/Makefile.am"
  "${recoil_source}/autogen.sh"
  "${recoil_source}/RecoilJets.cc"
  "${recoil_source}/RecoilJets.h"
  "${recoil_source}/PPG12SimWeight.h"
  "$macro_wrapper_source"
  "$macro_impl_source"
  "$trace_instrumenter"
  "$setup_script"
)
for input in "${contract_inputs[@]}"; do
  [[ "$input" == /* && -f "$input" && -s "$input" ]] || die "missing contract input: $input"
done
for external in "$yaml_cpp_library" "$yaml_cpp_include_dir" "$roounfold_library" \
  "$roounfold_pcm" \
  "$roounfold_include_dir" "$vertex_scan_data_file" "$mbd_correction_file" \
  "$truth_vertex_reweight_0mrad" "$truth_vertex_reweight_1p5mrad" \
  "$apply_model_dir" "$apply_npb_model"; do
  [[ "$external" == /* && "$external" != *$'\n'* && "$external" != *$'\r'* ]] || \
    die "estimator runtime asset path must be absolute and single-line: $external"
done
sha256_file() {
  python3 - "$1" <<'PY'
from pathlib import Path
import hashlib
import sys

h = hashlib.sha256()
with Path(sys.argv[1]).open("rb") as stream:
    for block in iter(lambda: stream.read(1024 * 1024), b""):
        h.update(block)
print(h.hexdigest())
PY
}

roounfold_header_tree_sha256() {
  python3 - "$1" "${roounfold_headers[@]}" <<'PY'
from pathlib import Path
import hashlib
import sys

root = Path(sys.argv[1])
names = sys.argv[2:]
digest = hashlib.sha256()
for name in names:
    path = root / name
    if not path.is_file() or path.stat().st_size <= 0:
        raise SystemExit(f"historical RooUnfold header is missing or empty: {path}")
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    digest.update(name.encode())
    digest.update(b"\0")
    digest.update(bytes.fromhex(observed))
print(digest.hexdigest())
PY
}

if [[ "$roounfold_library_overridden" == 1 ]]; then
  [[ -f "$roounfold_library" && -s "$roounfold_library" ]] || \
    die "overridden RooUnfold library is missing or empty"
fi
if [[ "$roounfold_pcm_overridden" == 1 ]]; then
  [[ -f "$roounfold_pcm" && -s "$roounfold_pcm" ]] || \
    die "overridden RooUnfold PCM is missing or empty"
fi
if [[ "$roounfold_include_dir_overridden" == 1 ]]; then
  [[ -d "$roounfold_include_dir" ]] || \
    die "overridden RooUnfold include directory is missing"
fi

if [[ -f "$roounfold_library" && -s "$roounfold_library" ]]; then
  [[ "$(sha256_file "$roounfold_library")" == "$expected_roounfold_library_sha256" ]] || \
    die "historical RooUnfold library hash differs"
fi
if [[ -f "$roounfold_pcm" && -s "$roounfold_pcm" ]]; then
  [[ "$(sha256_file "$roounfold_pcm")" == "$expected_roounfold_pcm_sha256" ]] || \
    die "historical RooUnfold PCM hash differs"
fi
if [[ -d "$roounfold_include_dir" ]]; then
  [[ "$(roounfold_header_tree_sha256 "$roounfold_include_dir")" == \
     "$expected_roounfold_header_tree_sha256" ]] || \
    die "historical RooUnfold header-tree hash differs"
fi

sha256_optional() {
  if [[ -f "$1" && -s "$1" ]]; then
    sha256_file "$1"
  elif [[ -d "$1" ]]; then
    python3 - "$1" <<'PY'
from pathlib import Path
import hashlib
import sys

root = Path(sys.argv[1])
h = hashlib.sha256()
for path in sorted(item for item in root.rglob("*") if item.is_file()):
    h.update(str(path.relative_to(root)).encode())
    h.update(b"\0")
    h.update(hashlib.sha256(path.read_bytes()).digest())
print("tree:" + h.hexdigest())
PY
  else
    printf 'missing-at-plan\n'
  fi
}

git_blob_sha256() {
  git "${git_safe[@]}" -C "$ppg_repo_real" show "$1" \
    | python3 -c 'import hashlib,sys; print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())'
}

[[ "$(git_blob_sha256 "${estimator_revision}:efficiencytool/RecoEffCalculator_TTreeReader.C")" == \
  "$expected_recoeff_sha256" ]] || die "canonical RecoEff source hash differs from 29f contract"
[[ "$(git_blob_sha256 "${estimator_revision}:efficiencytool/config_bdt_nom.yaml")" == \
  "$expected_recoeff_config_sha256" ]] || die "canonical RecoEff config hash differs from 29f contract"
[[ "$(git_blob_sha256 "${estimator_revision}:efficiencytool/config_bdt_nom_0rad.yaml")" == \
  "$expected_recoeff_period_config_0mrad_sha256" ]] || \
  die "canonical 0mrad RecoEff config hash differs from 29f contract"
[[ "$(git_blob_sha256 "${estimator_revision}:efficiencytool/config_bdt_nom_1p5mrad.yaml")" == \
  "$expected_recoeff_period_config_1p5mrad_sha256" ]] || \
  die "canonical 1p5mrad RecoEff config hash differs from 29f contract"
[[ "$(git_blob_sha256 "${estimator_revision}:FunWithxgboost/apply_BDT.C")" == \
  "$expected_apply_bdt_sha256" ]] || die "canonical apply_BDT source hash differs from 29f contract"
[[ "$(git_blob_sha256 "${estimator_revision}:FunWithxgboost/config_nom.yaml")" == \
  "$expected_apply_config_sha256" ]] || die "canonical apply_BDT config hash differs from 29f contract"
if [[ -f "$apply_npb_model" && -s "$apply_npb_model" ]]; then
  [[ "$(sha256_file "$apply_npb_model")" == "$expected_apply_npb_sha256" ]] || \
    die "canonical split NPB model hash differs"
fi
for model_index in "${!apply_model_names[@]}"; do
  model_name="${apply_model_names[$model_index]}"
  model_path="${apply_model_dir}/model_${model_name}_split_single_tmva.root"
  if [[ -f "$model_path" && -s "$model_path" ]]; then
    [[ "$(sha256_file "$model_path")" == "${expected_apply_model_hashes[$model_index]}" ]] || \
      die "canonical split model hash differs: $model_name"
  fi
done

if [[ "$ppg_binary_mode" == source_locked_runtime_import ]]; then
  ppg_import_fields="$(python3 - "$ppg_source_runtime_manifest" \
    "$expected_offline" "$ppg_revision" "$estimator_revision" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

manifest_path = Path(sys.argv[1]).resolve()
expected_offline, expected_revision, expected_estimator = sys.argv[2:]

def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def load(path: Path, label: str) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"invalid {label}: {exc}") from exc
    if not isinstance(value, dict):
        raise SystemExit(f"{label} must be a JSON object")
    return value

manifest = load(manifest_path, "source runtime manifest")
for field, expected in (
    ("schema_version", 1),
    ("runtime_profile", "new.17"),
    ("offline_main", expected_offline),
    ("isolated_build", True),
    ("estimator_revision", expected_estimator),
):
    if manifest.get(field) != expected:
        raise SystemExit(
            f"source runtime manifest {field} differs: "
            f"expected {expected!r}, observed {manifest.get(field)!r}"
        )

receipt_path = Path(str(manifest.get("build_receipt", "")))
receipt_sha = str(manifest.get("build_receipt_sha256", ""))
if not receipt_path.is_absolute() or not receipt_path.is_file():
    raise SystemExit("source runtime build receipt is missing")
if not re.fullmatch(r"[0-9a-f]{64}", receipt_sha) or digest(receipt_path) != receipt_sha:
    raise SystemExit("source runtime build receipt hash differs")
receipt = load(receipt_path, "source runtime build receipt")
source = receipt.get("ppg12_source")
rewrites = receipt.get("staged_rewrites")
if not isinstance(source, dict) or not isinstance(rewrites, dict):
    raise SystemExit("source runtime receipt lacks PPG12 provenance")
if source.get("revision") != expected_revision:
    raise SystemExit("source runtime receipt uses the wrong PPG12 source revision")
if source.get("working_tree_ignored") is not True:
    raise SystemExit("source runtime receipt does not exclude mutable PPG12 edits")
if source.get("rebuilt_against_common_runtime") is not True:
    raise SystemExit("source runtime receipt was not rebuilt against the common runtime")
if rewrites.get("archived_ppg12_binary_reused") is not False:
    raise SystemExit("source runtime receipt uses the ABI-incompatible archived binary")
if rewrites.get("ppg12_source_locked_rebuild") is not True:
    raise SystemExit("source runtime receipt lacks a source-locked PPG12 rebuild")

rows = [
    row for row in manifest.get("files", [])
    if isinstance(row, dict) and row.get("role") == "libCaloAna24.so"
]
if len(rows) != 1:
    raise SystemExit("source runtime manifest must contain exactly one libCaloAna24.so role")
library_path = Path(str(rows[0].get("path", "")))
library_sha = str(rows[0].get("sha256", ""))
if not library_path.is_absolute() or not library_path.is_file():
    raise SystemExit("source runtime libCaloAna24.so is missing")
if not re.fullmatch(r"[0-9a-f]{64}", library_sha) or digest(library_path) != library_sha:
    raise SystemExit("source runtime libCaloAna24.so hash differs")
for path in (manifest_path, receipt_path, library_path):
    if any(character in str(path) for character in ("\n", "\r", "\t")):
        raise SystemExit("source runtime provenance path is not single-line")
print("\t".join((
    str(receipt_path.resolve()),
    str(library_path.resolve()),
    digest(manifest_path),
    receipt_sha,
    library_sha,
)))
PY
)" || die "source-locked PPG12 runtime import validation failed"
  IFS=$'\t' read -r ppg_origin_build_receipt ppg_origin_library \
    ppg_origin_manifest_sha256 ppg_origin_receipt_sha256 \
    ppg_origin_library_sha256 <<<"$ppg_import_fields"
  [[ -n "$ppg_origin_build_receipt" && -n "$ppg_origin_library" && \
     "$ppg_origin_manifest_sha256" =~ ^[0-9a-f]{64}$ && \
     "$ppg_origin_receipt_sha256" =~ ^[0-9a-f]{64}$ && \
     "$ppg_origin_library_sha256" =~ ^[0-9a-f]{64}$ ]] || \
    die "source-locked PPG12 runtime import metadata is malformed"
fi

contract_token="$({
  printf '%s\n' \
    'schema_version=1' \
    'purpose=ppg12_paired_oracle_recoil_runtime' \
    'runtime_profile=new.17' \
    "offline_main=${expected_offline}" \
    "output_dir=${output_dir}" \
    "jobs=${jobs}"
  printf 'ppg_repo=%s\nppg_revision=%s\n' "$ppg_repo_real" "$ppg_revision"
  printf 'ppg_binary_mode=%s\n' "$ppg_binary_mode"
  if [[ "$ppg_binary_mode" == source_locked_runtime_import ]]; then
    printf 'ppg_origin_manifest=%s sha256=%s\n' \
      "$ppg_source_runtime_manifest" "$ppg_origin_manifest_sha256"
    printf 'ppg_origin_receipt=%s sha256=%s\n' \
      "$ppg_origin_build_receipt" "$ppg_origin_receipt_sha256"
    printf 'ppg_origin_library=%s sha256=%s\n' \
      "$ppg_origin_library" "$ppg_origin_library_sha256"
  fi
  for source_name in "${ppg_source_names[@]}"; do
    printf 'ppg_source=%s sha256=%s\n' "$source_name" "$(
      git "${git_safe[@]}" -C "$ppg_repo_real" show \
        "${ppg_revision}:anatreemaker/source/${source_name}" \
        | python3 -c 'import hashlib,sys; print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())'
    )"
  done
  printf 'estimator_revision=%s\n' "$estimator_revision"
  printf 'expected_recoeff_sha256=%s\n' "$expected_recoeff_sha256"
  printf 'expected_recoeff_config_sha256=%s\n' "$expected_recoeff_config_sha256"
  printf 'expected_recoeff_period_config_0mrad_sha256=%s\n' \
    "$expected_recoeff_period_config_0mrad_sha256"
  printf 'expected_recoeff_period_config_1p5mrad_sha256=%s\n' \
    "$expected_recoeff_period_config_1p5mrad_sha256"
  printf 'expected_truth_vertex_reweight_0mrad_sha256=%s\n' \
    "$expected_truth_vertex_reweight_0mrad_sha256"
  printf 'expected_truth_vertex_reweight_1p5mrad_sha256=%s\n' \
    "$expected_truth_vertex_reweight_1p5mrad_sha256"
  printf 'expected_apply_bdt_sha256=%s\n' "$expected_apply_bdt_sha256"
  printf 'expected_apply_config_sha256=%s\n' "$expected_apply_config_sha256"
  for model_index in "${!apply_model_names[@]}"; do
    printf 'expected_apply_model=%s sha256=%s\n' \
      "${apply_model_names[$model_index]}" "${expected_apply_model_hashes[$model_index]}"
  done
  printf 'expected_apply_npb_sha256=%s\n' "$expected_apply_npb_sha256"
  printf 'expected_roounfold_library_sha256=%s\n' "$expected_roounfold_library_sha256"
  printf 'expected_roounfold_pcm_sha256=%s\n' "$expected_roounfold_pcm_sha256"
  printf 'expected_roounfold_header_tree_sha256=%s\n' \
    "$expected_roounfold_header_tree_sha256"
  for source_name in "${estimator_source_names[@]}"; do
    printf 'estimator_source=%s sha256=%s\n' "$source_name" "$(
      git "${git_safe[@]}" -C "$ppg_repo_real" show \
        "${estimator_revision}:efficiencytool/${source_name}" \
        | python3 -c 'import hashlib,sys; print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())'
    )"
  done
  for source_name in "${apply_stage_source_names[@]}"; do
    printf 'apply_stage_source=%s sha256=%s\n' "$source_name" "$(
      git "${git_safe[@]}" -C "$ppg_repo_real" show \
        "${estimator_revision}:FunWithxgboost/${source_name}" \
        | python3 -c 'import hashlib,sys; print(hashlib.sha256(sys.stdin.buffer.read()).hexdigest())'
    )"
  done
  for external in "$yaml_cpp_library" "${yaml_cpp_include_dir}/yaml-cpp" "$roounfold_library" \
    "$roounfold_pcm" \
    "$roounfold_include_dir" "$vertex_scan_data_file" "$mbd_correction_file" \
    "$truth_vertex_reweight_0mrad" "$truth_vertex_reweight_1p5mrad" \
    "$apply_model_dir" "$apply_npb_model"; do
    printf 'estimator_asset=%s sha256=%s\n' "$external" "$(sha256_optional "$external")"
  done
  for input in "${contract_inputs[@]}"; do
    printf 'input=%s sha256=%s\n' "$input" "$(sha256_file "$input")"
  done
} | python3 -c 'import hashlib,sys; print("ppg12-new17:" + hashlib.sha256(sys.stdin.buffer.read()).hexdigest())')"

if [[ "$mode" == plan ]]; then
  if [[ "$ppg_binary_mode" == source_locked_runtime_import ]]; then
    ppg_build_description="exact source-locked libCaloAna24.so import plus fresh libRecoilJets.so"
  else
    ppg_build_description="source-locked libCaloAna24.so plus libRecoilJets.so with renamed PPG12-oracle photon builder"
  fi
  cat <<EOF
PPG12_ORACLE_NEW17_BUILD_PLAN
  output_root: ${output_dir}
  runtime: new.17
  OFFLINE_MAIN: ${expected_offline}
  build: ${ppg_build_description}
  ppg_binary_mode: ${ppg_binary_mode}
  ppg_source_revision: ${ppg_revision}
  ppg_origin_manifest_sha256: ${ppg_origin_manifest_sha256:-not-applicable}
  ppg_origin_library_sha256: ${ppg_origin_library_sha256:-not-applicable}
  estimator_revision: ${estimator_revision}
  estimator: exact RecoEff/config/yield sources plus dual uninstrumented/instrumented execution
  RooUnfold: exact historical lib=${expected_roounfold_library_sha256} pcm=${expected_roounfold_pcm_sha256} headers=${expected_roounfold_header_tree_sha256}
  release copies: exact new.17 libcalo_reco.so, libclusteriso.so, libjetbase.so (cp -L)
  macro staging: PP wrapper/implementation with exact-count path rewrites
  validation: forbidden-route scan, readelf, ldd, and ROOT load/header smoke
  jobs: ${jobs}
  token: ${contract_token}
No files were created.  Re-run with --build --token '${contract_token}'.
EOF
  exit 0
fi

[[ "$provided_token" == "$contract_token" ]] || die "authorization token does not match the current build contract"
[[ "${RJ_PPG12_ORACLE_NEW17_CLEAN_ENV:-0}" != 1 ]] || clean_env=1

# The shell executing this file may carry an analysis release or private
# installation in multiple routing variables.  Re-exec with a blank
# environment before sourcing the one supported release.
if [[ "${clean_env:-0}" != 1 ]]; then
  self="${script_dir}/$(basename "${BASH_SOURCE[0]}")"
  reexec_args=(
    --build --token "$provided_token" --output-dir "$output_dir"
    --setup-script "$setup_script" --jobs "$jobs"
    --photon-source-dir "$photon_source_dir"
    --ppg-repo "$ppg_repo_real" --ppg-revision "$ppg_revision"
    --estimator-revision "$estimator_revision"
    --yaml-cpp-library "$yaml_cpp_library"
    --yaml-cpp-include-dir "$yaml_cpp_include_dir"
    --roounfold-library "$roounfold_library"
    --roounfold-include-dir "$roounfold_include_dir"
    --roounfold-pcm "$roounfold_pcm"
    --vertex-scan-data-file "$vertex_scan_data_file"
    --mbd-correction-file "$mbd_correction_file"
    --truth-vertex-reweight-0mrad "$truth_vertex_reweight_0mrad"
    --truth-vertex-reweight-1p5mrad "$truth_vertex_reweight_1p5mrad"
    --apply-model-dir "$apply_model_dir"
    --apply-npb-model "$apply_npb_model"
  )
  if [[ "$ppg_binary_mode" == source_locked_runtime_import ]]; then
    reexec_args+=(--ppg-source-runtime-manifest "$ppg_source_runtime_manifest")
  fi
  exec /usr/bin/env -i \
    HOME="${HOME:-/tmp}" \
    USER="${USER:-unknown}" \
    LOGNAME="${LOGNAME:-${USER:-unknown}}" \
    SHELL=/bin/bash \
    PATH=/usr/bin:/bin:/usr/sbin:/sbin \
    RJ_PPG12_ORACLE_NEW17_CLEAN_ENV=1 \
    /bin/bash --noprofile --norc "$self" "${reexec_args[@]}"
fi

umask 077

unset OFFLINE_MAIN MYINSTALL ROOT_INCLUDE_PATH LD_LIBRARY_PATH PYTHONPATH \
  CMAKE_PREFIX_PATH CPATH CPLUS_INCLUDE_PATH LIBRARY_PATH PKG_CONFIG_PATH
export PATH=/usr/bin:/bin:/usr/sbin:/sbin
set +e
set +u
# shellcheck disable=SC1090
source "$setup_script" -n new.17
setup_status=$?
set -u
set -e
[[ $setup_status -eq 0 ]] || die "sPHENIX new.17 setup failed"
[[ "${OFFLINE_MAIN:-}" == "$expected_offline" ]] || \
  die "setup resolved unexpected OFFLINE_MAIN: ${OFFLINE_MAIN:-<unset>}"

for command_name in python3 make aclocal automake autoconf libtoolize root root-config readelf ldd cmp git awk; do
  command -v "$command_name" >/dev/null 2>&1 || die "required build command is unavailable: $command_name"
done
for estimator_asset in "$yaml_cpp_library" "$roounfold_library" "$roounfold_pcm" \
  "$vertex_scan_data_file" "$mbd_correction_file" \
  "$truth_vertex_reweight_0mrad" "$truth_vertex_reweight_1p5mrad" \
  "$apply_npb_model"; do
  [[ -f "$estimator_asset" && -s "$estimator_asset" ]] || \
    die "estimator runtime asset is missing or empty: $estimator_asset"
done
[[ "$(sha256_file "$roounfold_library")" == "$expected_roounfold_library_sha256" ]] || \
  die "historical RooUnfold library hash differs"
[[ "$(sha256_file "$roounfold_pcm")" == "$expected_roounfold_pcm_sha256" ]] || \
  die "historical RooUnfold PCM hash differs"
[[ "$(basename "$roounfold_pcm")" == RooUnfoldDict_rdict.pcm ]] || \
  die "RooUnfold PCM must be named RooUnfoldDict_rdict.pcm"
[[ -d "${yaml_cpp_include_dir}/yaml-cpp" ]] || \
  die "yaml-cpp header tree is missing: ${yaml_cpp_include_dir}/yaml-cpp"
[[ -f "${yaml_cpp_include_dir}/yaml-cpp/yaml.h" && \
   -s "${yaml_cpp_include_dir}/yaml-cpp/yaml.h" ]] || \
  die "yaml-cpp umbrella header is missing: ${yaml_cpp_include_dir}/yaml-cpp/yaml.h"
[[ "$(sha256_file "$truth_vertex_reweight_0mrad")" == \
  "$expected_truth_vertex_reweight_0mrad_sha256" ]] || \
  die "canonical 0mrad truth-vertex reweight ROOT hash differs"
[[ "$(sha256_file "$truth_vertex_reweight_1p5mrad")" == \
  "$expected_truth_vertex_reweight_1p5mrad_sha256" ]] || \
  die "canonical 1p5mrad truth-vertex reweight ROOT hash differs"
for model_index in "${!apply_model_names[@]}"; do
  model_name="${apply_model_names[$model_index]}"
  model_path="${apply_model_dir}/model_${model_name}_split_single_tmva.root"
  [[ -f "$model_path" && -s "$model_path" ]] || \
    die "apply_BDT split model is missing or empty: $model_path"
  [[ "$(sha256_file "$model_path")" == "${expected_apply_model_hashes[$model_index]}" ]] || \
    die "apply_BDT split model hash differs: $model_name"
done
[[ "$(sha256_file "$apply_npb_model")" == "$expected_apply_npb_sha256" ]] || \
  die "apply_BDT split NPB model hash differs"
[[ -d "$roounfold_include_dir" ]] || \
  die "RooUnfold include directory is missing: $roounfold_include_dir"
for header in "${roounfold_headers[@]}"; do
  [[ -f "${roounfold_include_dir}/${header}" && -s "${roounfold_include_dir}/${header}" ]] || \
    die "RooUnfold header is missing: ${roounfold_include_dir}/${header}"
done
[[ "$(roounfold_header_tree_sha256 "$roounfold_include_dir")" == \
   "$expected_roounfold_header_tree_sha256" ]] || \
  die "historical RooUnfold header-tree hash differs"

forbidden_routes() {
  local value="${PATH:-}:${LD_LIBRARY_PATH:-}:${ROOT_INCLUDE_PATH:-}:${PYTHONPATH:-}:${CMAKE_PREFIX_PATH:-}:${CPATH:-}:${CPLUS_INCLUDE_PATH:-}:${LIBRARY_PATH:-}:${PKG_CONFIG_PATH:-}"
  [[ "$value" != *'/release/release_ana/'* ]] || return 1
  [[ "$value" != *'/sphenix/u/patsfan753/thesisAnalysis/install'* ]] || return 1
  [[ "$value" != *'/sphenix/user/patsfan753/install'* ]] || return 1
  [[ "$value" != *'/sphenix/u/patsfan753/thesisAnalysis_auau/install'* ]] || return 1
  return 0
}
forbidden_routes || die "new.17 setup retained an analysis-release or private-install route"
base_ld_library_path="${LD_LIBRARY_PATH:-}"
base_root_include_path="${ROOT_INCLUDE_PATH:-}"
root_libdir="$(root-config --libdir)"
[[ "$root_libdir" == /* && -d "$root_libdir" ]] || \
  die "ROOT library directory is missing or non-absolute: $root_libdir"
release_link_dirs=()
for link_dir in "${OFFLINE_MAIN}/lib" "${OFFLINE_MAIN}/lib64" "$root_libdir"; do
  [[ -d "$link_dir" ]] && release_link_dirs+=("$link_dir")
done
[[ ${#release_link_dirs[@]} -ge 2 ]] || \
  die "could not establish sealed new.17 and ROOT linker directories"
release_ldflags=""
for link_dir in "${release_link_dirs[@]}"; do
  release_ldflags+=" -L${link_dir}"
done
sealed_link_dirs="$(IFS=:; printf '%s' "${release_link_dirs[*]}")"

mkdir -p "$output_dir"
state_file="${output_dir}/BUILD_STATE"
printf 'RUNNING\n' > "$state_file"
completed=0
finish_state() {
  if [[ $completed -eq 0 ]]; then
    printf 'FAILED\n' > "$state_file"
  fi
}
trap finish_state EXIT

stage_root="${output_dir}/stage"
build_root="${output_dir}/build"
install_root="${output_dir}/install"
runtime_root="${output_dir}/runtime"
log_root="${output_dir}/logs"
recoil_stage="${stage_root}/RecoilJets"
ppg_stage="${stage_root}/PPG12CaloAna24"
estimator_stage="${stage_root}/PPG12Estimator"
mkdir -p "$stage_root" "$build_root/recoiljets" "$build_root/ppg12" \
  "$recoil_stage" "$install_root" "$runtime_root/lib" "$runtime_root/include/caloreco" \
  "$runtime_root/include/caloana" "$runtime_root/macros" "$log_root" "$ppg_stage" \
  "$estimator_stage" "$runtime_root/estimator/source" \
  "$runtime_root/estimator/macros" "$runtime_root/estimator/include" \
  "$runtime_root/estimator/config" "$runtime_root/estimator/data" \
  "$runtime_root/estimator/apply" "$runtime_root/estimator/apply/binned_models" \
  "$runtime_root/estimator/apply/npb_models" "$runtime_root/provenance"

for source_name in "${ppg_source_names[@]}"; do
  git "${git_safe[@]}" -C "$ppg_repo_real" show \
    "${ppg_revision}:anatreemaker/source/${source_name}" \
    > "${ppg_stage}/${source_name}"
done
for source_name in "${apply_stage_source_names[@]}"; do
  git "${git_safe[@]}" -C "$ppg_repo_real" show \
    "${estimator_revision}:FunWithxgboost/${source_name}" \
    > "${runtime_root}/estimator/apply/${source_name}"
done
for model_name in "${apply_model_names[@]}"; do
  source_model="${apply_model_dir}/model_${model_name}_split_single_tmva.root"
  target_model="${runtime_root}/estimator/apply/binned_models/model_${model_name}_split_single_tmva.root"
  cp -f "$source_model" "$target_model"
  cmp -s "$source_model" "$target_model" || \
    die "staged apply_BDT model differs from source: $model_name"
done
cp -f "$apply_npb_model" \
  "${runtime_root}/estimator/apply/npb_models/npb_score_split_tmva.root"
cmp -s "$apply_npb_model" \
  "${runtime_root}/estimator/apply/npb_models/npb_score_split_tmva.root" || \
  die "staged apply_BDT NPB model differs from source"
chmod 700 "${ppg_stage}/autogen.sh"

# Export estimator sources from their own canonical revision.  The mutable
# checkout and the older CaloAna reconstruction revision are never consulted.
for source_name in "${estimator_source_names[@]}"; do
  git "${git_safe[@]}" -C "$ppg_repo_real" show \
    "${estimator_revision}:efficiencytool/${source_name}" \
    > "${estimator_stage}/${source_name}"
done
cp -f "${estimator_stage}/RecoEffCalculator_TTreeReader.C" \
  "${runtime_root}/estimator/source/RecoEffCalculator_TTreeReader.C"
cp -f "${estimator_stage}/CalculatePhotonYield.C" \
  "${runtime_root}/estimator/source/CalculatePhotonYield.C"
cp -f "${estimator_stage}/CrossSectionWeights.h" \
  "${runtime_root}/estimator/include/CrossSectionWeights.h"
cp -f "${estimator_stage}/TruthVertexReweightLoader.h" \
  "${runtime_root}/estimator/include/TruthVertexReweightLoader.h"
cp -f "${estimator_stage}/config_bdt_nom.yaml" \
  "${runtime_root}/estimator/config/config_bdt_nom.yaml"
cp -f "${estimator_stage}/config_bdt_nom_0rad.yaml" \
  "${runtime_root}/estimator/config/config_bdt_nom_0rad.yaml"
cp -f "${estimator_stage}/config_bdt_nom_1p5mrad.yaml" \
  "${runtime_root}/estimator/config/config_bdt_nom_1p5mrad.yaml"
cp -L "$yaml_cpp_library" "${runtime_root}/lib/libyaml-cpp.so"
cp -R "${yaml_cpp_include_dir}/yaml-cpp" "${runtime_root}/estimator/include/"
source_yaml_tree_sha256="$(sha256_optional "${yaml_cpp_include_dir}/yaml-cpp")"
staged_yaml_tree_sha256="$(sha256_optional "${runtime_root}/estimator/include/yaml-cpp")"
[[ "$source_yaml_tree_sha256" == "$staged_yaml_tree_sha256" ]] || \
  die "staged yaml-cpp header tree differs from source"
yaml_cpp_header_receipt="${runtime_root}/estimator/yaml_cpp_header_tree_receipt.json"
python3 - "${yaml_cpp_include_dir}/yaml-cpp" \
  "${runtime_root}/estimator/include/yaml-cpp" "$yaml_cpp_header_receipt" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

source, staged, receipt = map(Path, sys.argv[1:])

def inventory(root: Path) -> tuple[str, list[dict[str, str]]]:
    rows = []
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = str(path.relative_to(root))
        file_digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append({"relative_path": relative, "sha256": file_digest})
        digest.update(relative.encode())
        digest.update(b"\0")
        digest.update(bytes.fromhex(file_digest))
    if not rows or not (root / "yaml.h").is_file():
        raise SystemExit("yaml-cpp header inventory is empty or lacks yaml.h")
    return digest.hexdigest(), rows

source_digest, source_rows = inventory(source)
staged_digest, staged_rows = inventory(staged)
if source_rows != staged_rows or source_digest != staged_digest:
    raise SystemExit("sealed yaml-cpp header inventory differs from source")
payload = {
    "schema_version": 1,
    "role": "ppg_recoeff_yaml_cpp_header_tree",
    "source_root": str(source.resolve()),
    "include_root": str(staged.parent.resolve()),
    "staged_tree": str(staged.resolve()),
    "tree_sha256": staged_digest,
    "files": staged_rows,
}
receipt.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
cp -L "$roounfold_library" "${runtime_root}/lib/libRooUnfold.so"
cp -f "$roounfold_pcm" "${runtime_root}/lib/RooUnfoldDict_rdict.pcm"
for header in "${roounfold_headers[@]}"; do
  cp -f "${roounfold_include_dir}/${header}" \
    "${runtime_root}/estimator/include/${header}"
  cmp -s "${roounfold_include_dir}/${header}" \
    "${runtime_root}/estimator/include/${header}" || \
    die "staged historical RooUnfold header differs from source: ${header}"
done
roounfold_header_receipt="${runtime_root}/estimator/roounfold_header_tree_receipt.json"
python3 - "$roounfold_include_dir" "${runtime_root}/estimator/include" \
  "$roounfold_header_receipt" "$expected_roounfold_header_tree_sha256" \
  "${roounfold_headers[@]}" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

source = Path(sys.argv[1]).resolve()
include_root = Path(sys.argv[2]).resolve()
receipt = Path(sys.argv[3])
expected_tree_digest = sys.argv[4]
names = sys.argv[5:]
expected = [
    "RooUnfold.h",
    "RooUnfoldResponse.h",
    "RooUnfoldBayes.h",
    "RooUnfoldBinByBin.h",
    "RooUnfoldErrors.h",
    "RooUnfoldInvert.h",
    "RooUnfoldParms.h",
    "RooUnfoldSvd.h",
    "RooUnfoldTUnfold.h",
]
if names != expected or len(set(names)) != len(expected):
    raise SystemExit("historical RooUnfold header contract differs from exact nine-file order")
rows = []
tree_digest = hashlib.sha256()
for name in names:
    source_path = source / name
    staged_path = include_root / name
    if not source_path.is_file() or not staged_path.is_file():
        raise SystemExit(f"historical RooUnfold header is missing: {name}")
    source_digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
    staged_digest = hashlib.sha256(staged_path.read_bytes()).hexdigest()
    if source_digest != staged_digest:
        raise SystemExit(f"historical RooUnfold header changed while staging: {name}")
    rows.append({"relative_path": name, "sha256": staged_digest})
    tree_digest.update(name.encode())
    tree_digest.update(b"\0")
    tree_digest.update(bytes.fromhex(staged_digest))
payload = {
    "schema_version": 1,
    "role": "ppg_recoeff_roounfold_header_tree",
    "source_root": str(source),
    "include_root": str(include_root),
    "tree_sha256": tree_digest.hexdigest(),
    "files": rows,
}
if payload["tree_sha256"] != expected_tree_digest:
    raise SystemExit("historical RooUnfold header tree differs from pinned digest")
receipt.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
PY
cp -f "$vertex_scan_data_file" \
  "${runtime_root}/estimator/data/data_histo_bdt_nom_vtxscan.root"
cp -f "$mbd_correction_file" \
  "${runtime_root}/estimator/data/MbdOut.corr"
cp -f "$truth_vertex_reweight_0mrad" \
  "${runtime_root}/estimator/data/truth_vertex_reweight_0mrad.root"
cmp -s "$truth_vertex_reweight_0mrad" \
  "${runtime_root}/estimator/data/truth_vertex_reweight_0mrad.root" || \
  die "staged 0mrad truth-vertex reweight ROOT differs from source"
cp -f "$truth_vertex_reweight_1p5mrad" \
  "${runtime_root}/estimator/data/truth_vertex_reweight_1p5mrad.root"
cmp -s "$truth_vertex_reweight_1p5mrad" \
  "${runtime_root}/estimator/data/truth_vertex_reweight_1p5mrad.root" || \
  die "staged 1p5mrad truth-vertex reweight ROOT differs from source"

for source_name in configure.ac Makefile.am autogen.sh RecoilJets.cc RecoilJets.h PPG12SimWeight.h; do
  cp -f "${recoil_source}/${source_name}" "${recoil_stage}/${source_name}"
done
cp -f "$canonical_photon_cc" "${recoil_stage}/PPG12OraclePhotonClusterBuilder.cc"
cp -f "$canonical_photon_h" "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h"

replace_exact() {
  local path="$1"
  local expected_count="$2"
  local old="$3"
  local new="$4"
  python3 - "$path" "$expected_count" "$old" "$new" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
expected = int(sys.argv[2])
old = sys.argv[3]
new = sys.argv[4]
text = path.read_text()
observed = text.count(old)
if observed != expected:
    raise SystemExit(
        f"exact rewrite count mismatch for {path}: expected {expected}, "
        f"observed {observed}, needle={old!r}"
    )
path.write_text(text.replace(old, new))
PY
}

replace_word_exact() {
  local path="$1"
  local expected_count="$2"
  local old="$3"
  local new="$4"
  python3 - "$path" "$expected_count" "$old" "$new" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
expected = int(sys.argv[2])
old = sys.argv[3]
new = sys.argv[4]
text = path.read_text()
pattern = re.compile(rf"\b{re.escape(old)}\b")
observed = len(pattern.findall(text))
if observed != expected:
    raise SystemExit(
        f"exact word-rewrite count mismatch for {path}: expected {expected}, "
        f"observed {observed}, token={old!r}"
    )
path.write_text(pattern.sub(new, text))
PY
}

recoeff_macro="${runtime_root}/estimator/macros/RecoEffCalculator_TTreeReader.C"
recoeff_trace_macro="${runtime_root}/estimator/macros/RecoEffCalculator_TTreeReader_trace.C"
recoeff_trace_receipt="${runtime_root}/estimator/recoeff_trace_transform_receipt.json"
cp -f "${estimator_stage}/RecoEffCalculator_TTreeReader.C" "$recoeff_macro"
replace_exact "$recoeff_macro" 1 \
  '/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so' \
  "${runtime_root}/lib/libyaml-cpp.so"
replace_exact "$recoeff_macro" 1 \
  '/sphenix/user/shuhangli/ppg12/efficiencytool/MbdOut.corr' \
  "${runtime_root}/estimator/data/MbdOut.corr"
python3 "$trace_instrumenter" \
  --input "$recoeff_macro" \
  --output "$recoeff_trace_macro" \
  --receipt "$recoeff_trace_receipt" \
  --source-revision "$estimator_revision"

# Default mode rebuilds the preserved PPG12 source against the same current
# new.17 headers and libraries used by the candidate.  Import mode instead
# carries forward one exact source-locked binary from a sealed runtime whose
# manifest and receipt were verified before the authorization token was made.
if [[ "$ppg_binary_mode" == source_locked_rebuild ]]; then
  replace_exact "${ppg_stage}/configure.ac" 1 \
    'CXXFLAGS="$CXXFLAGS -Wall -Werror"' \
    'CXXFLAGS="$CXXFLAGS -Wall -Wno-error"'
  replace_exact "${ppg_stage}/Makefile.am" 1 \
    $'-lcalotrigger_io \\ ' \
    $'-lcalotrigger_io \\'

  (
    cd "${build_root}/ppg12"
    export LD_LIBRARY_PATH="${install_root}/lib:${base_ld_library_path}"
    export ROOT_INCLUDE_PATH="${install_root}/include:${base_root_include_path}"
    export CPPFLAGS="-I${install_root}/include"
    export LDFLAGS="-L${install_root}/lib${release_ldflags}"
    /bin/bash "${ppg_stage}/autogen.sh" --prefix="$install_root"
    make -j "$jobs"
    make install
  ) >"${log_root}/ppg12_build.log" 2>&1 || {
    tail -n 80 "${log_root}/ppg12_build.log" >&2 || true
    die "isolated source-locked libCaloAna24 build failed"
  }
else
  printf 'PPG12 source-locked runtime import manifest=%s library_sha256=%s\n' \
    "$ppg_source_runtime_manifest" "$ppg_origin_library_sha256" \
    >"${log_root}/ppg12_binary_import.log"
fi

# These rewrites touch staged copies only.  The oracle builder is renamed so
# it can coexist with exact release libcalo_reco without an ODR/symbol
# collision.  It is compiled into libRecoilJets; release calibration and
# RawClusterBuilderTemplate remain untouched.
replace_exact "${recoil_stage}/configure.ac" 1 \
  'CXXFLAGS="$CXXFLAGS -Wall -Wextra -Werror -Wshadow"' \
  'CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Wno-error"'
replace_exact "${recoil_stage}/RecoilJets.cc" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloBase/PhotonClusterv1.h"' \
  '#include <calobase/PhotonClusterv1.h>'
replace_exact "${recoil_stage}/RecoilJets.cc" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"' \
  '// Oracle photon builder is registered by the staged Fun4All macro.'
replace_exact "${recoil_stage}/RecoilJets.h" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/coresoftware_local/offline/packages/CaloReco/PhotonClusterBuilder.h"' \
  '// Oracle photon builder is registered by the staged Fun4All macro.'
replace_word_exact "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" 4 \
  'PhotonClusterBuilder' 'PPG12OraclePhotonClusterBuilder'
replace_word_exact "${recoil_stage}/PPG12OraclePhotonClusterBuilder.cc" 31 \
  'PhotonClusterBuilder' 'PPG12OraclePhotonClusterBuilder'
replace_exact "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" 3 \
  'CALORECO_PHOTONCLUSTERBUILDER_H' \
  'CALOANA_PPG12ORACLEPHOTONCLUSTERBUILDER_H'
replace_exact "${recoil_stage}/Makefile.am" 1 \
  $'libRecoilJets_la_SOURCES = \\\n  RecoilJets.cc' \
  $'libRecoilJets_la_SOURCES = \\\n  RecoilJets.cc \\\n  PPG12OraclePhotonClusterBuilder.cc'
replace_exact "${recoil_stage}/Makefile.am" 1 \
  $'  -lg4detectors_io \\\n  -lphg4hit' \
  $'  -lg4detectors_io \\\n  -lphg4hit \\\n  -lg4eval \\\n  -lcdbobjects \\\n  -lffamodules \\\n  -lglobalvertex_io \\\n  -lTMVA \\\n  -lTMVAUtils'
replace_exact "${recoil_stage}/Makefile.am" 1 \
  'libRecoilJets_la_LDFLAGS = -no-undefined -version-info 0:0:0' \
  'libRecoilJets_la_LDFLAGS = $(AM_LDFLAGS) -L$(ROOTSYS)/lib -no-undefined -version-info 0:0:0'
replace_exact "${recoil_stage}/Makefile.am" 1 \
  $'pkginclude_HEADERS = \\\n  RecoilJets.h' \
  $'pkginclude_HEADERS = \\\n  RecoilJets.h \\\n  PPG12OraclePhotonClusterBuilder.h'

(
  cd "${build_root}/recoiljets"
  export LD_LIBRARY_PATH="${install_root}/lib:${base_ld_library_path}"
  export ROOT_INCLUDE_PATH="${install_root}/include:${base_root_include_path}"
  export CPPFLAGS="-I${install_root}/include"
  export LDFLAGS="-L${install_root}/lib${release_ldflags}"
  /bin/bash "${recoil_stage}/autogen.sh" --prefix="$install_root"
  make -j "$jobs"
  make install
) >"${log_root}/recoiljets_build.log" 2>&1 || {
  tail -n 80 "${log_root}/recoiljets_build.log" >&2 || true
  die "isolated libRecoilJets build failed"
}

resolve_installed_lib() {
  local name="$1"
  local candidate
  for candidate in "${install_root}/lib/${name}" "${install_root}/lib64/${name}"; do
    if [[ -r "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  die "isolated build did not install ${name}"
}

resolve_release_lib() {
  local name="$1"
  local candidate
  for candidate in "${OFFLINE_MAIN}/lib/${name}" "${OFFLINE_MAIN}/lib64/${name}"; do
    if [[ -r "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  die "new.17 does not provide ${name}"
}

built_recoil="$(resolve_installed_lib libRecoilJets.so)"
if [[ "$ppg_binary_mode" == source_locked_runtime_import ]]; then
  built_ppg="$ppg_origin_library"
else
  built_ppg="$(resolve_installed_lib libCaloAna24.so)"
fi
release_calo="$(resolve_release_lib libcalo_reco.so)"
release_clusteriso="$(resolve_release_lib libclusteriso.so)"
release_jetbase="$(resolve_release_lib libjetbase.so)"

cp -L "$release_calo" "${runtime_root}/lib/libcalo_reco.so"
cp -L "$built_ppg" "${runtime_root}/lib/libCaloAna24.so"
cp -L "$built_recoil" "${runtime_root}/lib/libRecoilJets.so"
cp -L "$release_clusteriso" "${runtime_root}/lib/libclusteriso.so"
cp -L "$release_jetbase" "${runtime_root}/lib/libjetbase.so"
cmp -s "$release_calo" "${runtime_root}/lib/libcalo_reco.so" || \
  die "copied libcalo_reco differs from exact new.17 source"
cmp -s "$built_ppg" "${runtime_root}/lib/libCaloAna24.so" || \
  die "copied libCaloAna24 differs from the selected source-locked binary"
cmp -s "$release_clusteriso" "${runtime_root}/lib/libclusteriso.so" || \
  die "copied libclusteriso differs from exact new.17 source"
cmp -s "$release_jetbase" "${runtime_root}/lib/libjetbase.so" || \
  die "copied libjetbase differs from exact new.17 source"
cmp -s "$yaml_cpp_library" "${runtime_root}/lib/libyaml-cpp.so" || \
  die "copied libyaml-cpp differs from sealed estimator source"
cmp -s "$roounfold_library" "${runtime_root}/lib/libRooUnfold.so" || \
  die "copied libRooUnfold differs from sealed estimator source"
cmp -s "$roounfold_pcm" "${runtime_root}/lib/RooUnfoldDict_rdict.pcm" || \
  die "copied RooUnfold dictionary PCM differs from sealed estimator source"
ppg_provenance_manifest=""
ppg_provenance_receipt=""
if [[ "$ppg_binary_mode" == source_locked_runtime_import ]]; then
  ppg_provenance_manifest="${runtime_root}/provenance/ppg_source_runtime_manifest.json"
  ppg_provenance_receipt="${runtime_root}/provenance/ppg_source_build_receipt.json"
  cp -f "$ppg_source_runtime_manifest" "$ppg_provenance_manifest"
  cp -f "$ppg_origin_build_receipt" "$ppg_provenance_receipt"
  [[ "$(sha256_file "$ppg_provenance_manifest")" == "$ppg_origin_manifest_sha256" ]] || \
    die "copied source runtime manifest hash differs"
  [[ "$(sha256_file "$ppg_provenance_receipt")" == "$ppg_origin_receipt_sha256" ]] || \
    die "copied source runtime receipt hash differs"
  [[ "$(sha256_file "${runtime_root}/lib/libCaloAna24.so")" == \
     "$ppg_origin_library_sha256" ]] || \
    die "imported libCaloAna24 hash differs after staging"
fi
cp -f "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" \
  "${runtime_root}/include/caloana/PPG12OraclePhotonClusterBuilder.h"
# The paired-oracle manifest retains its established role/path.  This file's
# contents declare only the renamed oracle class, never the release class.
cp -f "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" \
  "${runtime_root}/include/caloreco/PhotonClusterBuilder.h"
cp -f "${recoil_stage}/PPG12OraclePhotonClusterBuilder.h" \
  "${runtime_root}/include/caloreco/PPG12OraclePhotonClusterBuilder.h"
cp -f "${recoil_stage}/RecoilJets.h" "${runtime_root}/include/caloana/RecoilJets.h"

# Preserve each copied ELF's SONAME inside the sealed library directory.
for runtime_lib in "${runtime_root}"/lib/lib*.so; do
  soname="$(readelf -d "$runtime_lib" 2>/dev/null | awk -F'[][]' '/SONAME/ {print $2; exit}' || true)"
  if [[ -n "$soname" && "$soname" != "$(basename "$runtime_lib")" ]]; then
    ln -sfn "$(basename "$runtime_lib")" "${runtime_root}/lib/${soname}"
  fi
done

calo_calib="${OFFLINE_MAIN}/rootmacros/Calo_Calib.C"
[[ -f "$calo_calib" && -s "$calo_calib" ]] || \
  die "stock new.17 Calo_Calib.C is missing: $calo_calib"

runtime_wrapper="${runtime_root}/macros/Fun4All_recoilJets.C"
runtime_impl="${runtime_root}/macros/Fun4All_recoilJets_unified_impl.C"
cp -f "$macro_wrapper_source" "$runtime_wrapper"
cp -f "$macro_impl_source" "$runtime_impl"

replace_word_exact "$runtime_impl" 16 \
  'PhotonClusterBuilder' 'PPG12OraclePhotonClusterBuilder'
replace_exact "$runtime_wrapper" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Fun4All_recoilJets_unified_impl.C"' \
  "#include \"${runtime_impl}\""
replace_exact "$runtime_impl" 3 \
  '/sphenix/u/patsfan753/thesisAnalysis/install/include' \
  "${runtime_root}/include"
replace_exact "$runtime_impl" 2 \
  '/sphenix/u/patsfan753/thesisAnalysis_auau/install/include' \
  "${runtime_root}/include"
replace_exact "$runtime_impl" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src_AuAu/RecoilJets_AuAu.h"' \
  "#include \"${runtime_root}/include/caloana/RecoilJets.h\""
replace_exact "$runtime_impl" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src/RecoilJets.h"' \
  "#include \"${runtime_root}/include/caloana/RecoilJets.h\""
replace_exact "$runtime_impl" 1 \
  '#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Calo_Calib.C"' \
  "#include \"${calo_calib}\""
replace_exact "$runtime_impl" 1 \
  'R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis_auau/install/lib/libRecoilJetsAuAu.so)' \
  'R__LOAD_LIBRARY(libRecoilJets.so)'
replace_exact "$runtime_impl" 1 \
  'R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libRecoilJets.so)' \
  'R__LOAD_LIBRARY(libRecoilJets.so)'
replace_exact "$runtime_impl" 1 \
  'R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)' \
  'R__LOAD_LIBRARY(libclusteriso.so)'
replace_exact "$runtime_impl" 1 \
  'R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)' \
  'R__LOAD_LIBRARY(libjetbase.so)'

for staged_runtime_file in "$runtime_wrapper" "$runtime_impl" \
  "${runtime_root}/include/caloana/RecoilJets.h" \
  "${runtime_root}/include/caloana/PPG12OraclePhotonClusterBuilder.h"; do
  if grep -E '/release/release_ana/|/sphenix/u/patsfan753/thesisAnalysis(_auau)?/install|/sphenix/user/patsfan753/install' \
      "$staged_runtime_file" >/dev/null; then
    die "staged runtime retains a forbidden route: $staged_runtime_file"
  fi
done

runtime_ld="${runtime_root}/lib:${base_ld_library_path}"
for runtime_lib in \
  "${runtime_root}/lib/libCaloAna24.so" \
  "${runtime_root}/lib/libcalo_reco.so" \
  "${runtime_root}/lib/libRecoilJets.so" \
  "${runtime_root}/lib/libclusteriso.so" \
  "${runtime_root}/lib/libjetbase.so" \
  "${runtime_root}/lib/libyaml-cpp.so" \
  "${runtime_root}/lib/libRooUnfold.so"; do
  basename_noext="$(basename "$runtime_lib" .so)"
  readelf -d "$runtime_lib" >"${log_root}/readelf_${basename_noext}.txt" 2>&1 || \
    die "readelf failed for $runtime_lib"
  LD_LIBRARY_PATH="$runtime_ld" ldd "$runtime_lib" \
    >"${log_root}/ldd_${basename_noext}.txt" 2>&1 || die "ldd failed for $runtime_lib"
  if grep -E 'not found|/release/release_ana/|/sphenix/u/patsfan753/thesisAnalysis(_auau)?/install|/sphenix/user/patsfan753/install' \
      "${log_root}/ldd_${basename_noext}.txt" >/dev/null; then
    die "ldd exposed an unresolved or forbidden dependency for $runtime_lib"
  fi
  if grep -E '/release/release_new/new\.[0-9]+' "${log_root}/ldd_${basename_noext}.txt" \
      | grep -vF "$expected_offline" >/dev/null; then
    die "ldd exposed a mixed new-release dependency for $runtime_lib"
  fi
done

smoke_macro="${build_root}/smoke_new17_runtime.C"
cat > "$smoke_macro" <<EOF
R__LOAD_LIBRARY(${runtime_root}/lib/libRooUnfold.so)
#include <caloana/PPG12OraclePhotonClusterBuilder.h>
#include <yaml-cpp/yaml.h>
#include <RooUnfold.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldBinByBin.h>
#include <RooUnfoldErrors.h>
#include <RooUnfoldInvert.h>
#include <RooUnfoldParms.h>
#include <RooUnfoldSvd.h>
#include <RooUnfoldTUnfold.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <iostream>
void smoke_new17_runtime()
{
  const char *libraries[] = {
    "${runtime_root}/lib/libCaloAna24.so",
    "${runtime_root}/lib/libcalo_reco.so",
    "${runtime_root}/lib/libclusteriso.so",
    "${runtime_root}/lib/libjetbase.so",
    "${runtime_root}/lib/libRecoilJets.so",
    "${runtime_root}/lib/libyaml-cpp.so",
    "${runtime_root}/lib/libRooUnfold.so"
  };
  for (const char *library : libraries)
  {
    const int status = gSystem->Load(library);
    std::cout << "PPG12_ORACLE_ROOT_LOAD path=" << library
              << " status=" << status << std::endl;
    if (status < 0) gSystem->Exit(91);
  }
  PPG12OraclePhotonClusterBuilder *builder = nullptr;
  if (builder != nullptr) gSystem->Exit(92);
  const YAML::Node yaml_header_smoke = YAML::Load("ppg12_oracle_smoke: 17");
  if (!yaml_header_smoke["ppg12_oracle_smoke"] ||
      yaml_header_smoke["ppg12_oracle_smoke"].as<int>() != 17)
  {
    gSystem->Exit(93);
  }
  TH1D measured("ppg12_oracle_measured", "", 2, 0.0, 2.0);
  TH1D truth("ppg12_oracle_truth", "", 2, 0.0, 2.0);
  TH2D migration("ppg12_oracle_migration", "", 2, 0.0, 2.0, 2, 0.0, 2.0);
  RooUnfoldResponse response(
    (const TH1 *)&measured, (const TH1 *)&truth, &migration,
    "ppg12_oracle_response", "", false);
  RooUnfoldBayes bayes(
    &response, &measured, 1, false, "ppg12_oracle_bayes", "");
  (void)bayes;
  std::cout << "PPG12_ORACLE_ROOUNFOLD_API_SMOKE_PASS" << std::endl;
}
EOF
(
  export LD_LIBRARY_PATH="$runtime_ld"
  export ROOT_INCLUDE_PATH="${runtime_root}/include:${runtime_root}/estimator/include:${base_root_include_path}"
  root -l -b -q "${smoke_macro}"
) >"${log_root}/root_smoke.log" 2>&1 || {
  tail -n 80 "${log_root}/root_smoke.log" >&2 || true
  die "ROOT load/header smoke failed"
}
[[ "$(grep -c '^PPG12_ORACLE_ROOT_LOAD ' "${log_root}/root_smoke.log")" -eq 7 ]] || \
  die "ROOT smoke did not load all seven runtime libraries"
grep -Fxq 'PPG12_ORACLE_ROOUNFOLD_API_SMOKE_PASS' "${log_root}/root_smoke.log" || \
  die "ROOT smoke did not compile and instantiate the historical RooUnfold API"
if grep -E 'TCling::(LoadPCM|RegisterModule|AutoParse)|Failed to load PCM|fatal error:|no matching constructor|redefinition of|cannot open shared object' \
    "${log_root}/root_smoke.log" >/dev/null; then
  die "ROOT smoke exposed a PCM, autoparse, header, or RooUnfold API fallback error"
fi

forbidden_routes || die "build contaminated the active shell with a forbidden route"

build_receipt="${output_dir}/build_receipt.json"
runtime_manifest="${output_dir}/runtime_manifest.json"
python3 - \
  "$build_receipt" "$output_dir" "$install_root" "$runtime_root" \
  "$expected_offline" "$calo_calib" "$jobs" "$sealed_link_dirs" \
  "$canonical_photon_cc" "$canonical_photon_h" \
  "$ppg_repo_real" "$ppg_revision" "$ppg_binary_mode" \
  "$ppg_provenance_manifest" "$ppg_origin_manifest_sha256" \
  "$ppg_provenance_receipt" "$ppg_origin_receipt_sha256" \
  "$ppg_origin_library_sha256" \
  "${ppg_stage}/CaloAna24.cc" "${ppg_stage}/CaloAna24.h" \
  "${ppg_stage}/configure.ac" "${ppg_stage}/Makefile.am" \
  "${recoil_source}/RecoilJets.cc" "${recoil_source}/RecoilJets.h" \
  "$macro_wrapper_source" "$macro_impl_source" \
  "$release_calo" "$release_clusteriso" "$release_jetbase" "$log_root" \
  "$estimator_revision" \
  "${runtime_root}/estimator/source/RecoEffCalculator_TTreeReader.C" \
  "${runtime_root}/estimator/include/CrossSectionWeights.h" \
  "${runtime_root}/estimator/include/TruthVertexReweightLoader.h" \
  "${runtime_root}/estimator/config/config_bdt_nom.yaml" \
  "${runtime_root}/estimator/config/config_bdt_nom_0rad.yaml" \
  "${runtime_root}/estimator/config/config_bdt_nom_1p5mrad.yaml" \
  "${runtime_root}/estimator/source/CalculatePhotonYield.C" \
  "${runtime_root}/estimator/apply/apply_BDT.C" \
  "${runtime_root}/estimator/apply/config_nom.yaml" \
  "$recoeff_macro" "$recoeff_trace_macro" "$recoeff_trace_receipt" \
  "$trace_instrumenter" \
  "${runtime_root}/lib/libyaml-cpp.so" \
  "$yaml_cpp_header_receipt" \
  "${runtime_root}/lib/libRooUnfold.so" \
  "${runtime_root}/lib/RooUnfoldDict_rdict.pcm" \
  "$roounfold_header_receipt" \
  "${runtime_root}/estimator/include/RooUnfoldResponse.h" \
  "${runtime_root}/estimator/include/RooUnfoldBayes.h" \
  "${runtime_root}/estimator/data/data_histo_bdt_nom_vtxscan.root" \
  "${runtime_root}/estimator/data/MbdOut.corr" \
  "${runtime_root}/estimator/data/truth_vertex_reweight_0mrad.root" \
  "${runtime_root}/estimator/data/truth_vertex_reweight_1p5mrad.root" <<'PY'
from pathlib import Path
import hashlib
import json
import os
import platform
import sys

(
    receipt, output_root, install_root, runtime_root, offline_main, calo_calib,
    jobs, sealed_link_dirs, photon_cc, photon_h, ppg_repo, ppg_revision,
    ppg_binary_mode, ppg_provenance_manifest, ppg_origin_manifest_sha256,
    ppg_provenance_receipt, ppg_origin_receipt_sha256,
    ppg_origin_library_sha256,
    ppg_cc, ppg_h, ppg_configure, ppg_makefile, recoil_cc, recoil_h, macro, impl,
    calo_source, clusteriso_source, jetbase_source, log_root,
    estimator_revision, recoeff_source, cross_section_header,
    truth_vertex_header, estimator_config, estimator_config_0mrad,
    estimator_config_1p5mrad, calculate_yield_source,
    apply_bdt_source, apply_config_source,
    recoeff_macro, recoeff_trace_macro, trace_receipt, trace_instrumenter,
    yaml_cpp, yaml_cpp_header_receipt, roounfold, roounfold_pcm,
    roounfold_header_receipt, roounfold_response_header, roounfold_bayes_header,
    vertex_scan_data, mbd_correction, truth_vertex_reweight_0mrad,
    truth_vertex_reweight_1p5mrad,
) = sys.argv[1:]

def digest(path: str | Path) -> str:
    value = Path(path)
    h = hashlib.sha256()
    with value.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

roounfold_header_data = json.loads(Path(roounfold_header_receipt).read_text())
roounfold_headers = [
    str(Path(str(roounfold_header_data["include_root"])) / str(row["relative_path"]))
    for row in roounfold_header_data.get("files", [])
    if isinstance(row, dict)
]
if len(roounfold_headers) != 9:
    raise SystemExit("sealed RooUnfold header receipt does not contain exactly nine headers")
if roounfold_response_header not in roounfold_headers or roounfold_bayes_header not in roounfold_headers:
    raise SystemExit("sealed RooUnfold Response/Bayes roles differ from header receipt")

source_paths = [
    photon_cc, photon_h, ppg_cc, ppg_h, ppg_configure, ppg_makefile,
    recoil_cc, recoil_h, macro, impl,
    calo_source, clusteriso_source, jetbase_source, calo_calib,
    recoeff_source, cross_section_header, truth_vertex_header,
    estimator_config, estimator_config_0mrad, estimator_config_1p5mrad,
    calculate_yield_source, recoeff_macro,
    apply_bdt_source, apply_config_source,
    recoeff_trace_macro, trace_receipt, trace_instrumenter, yaml_cpp,
    yaml_cpp_header_receipt, roounfold, roounfold_pcm, roounfold_header_receipt,
    vertex_scan_data, mbd_correction, truth_vertex_reweight_0mrad,
    truth_vertex_reweight_1p5mrad,
]
source_paths.extend(roounfold_headers)
apply_root = Path(apply_bdt_source).parent
apply_model_names = (
    "base", "base_vr", "base_v0", "base_v1", "base_v2", "base_v3",
    "base_E", "base_v0E", "base_v1E", "base_v2E", "base_v3E",
)
apply_models = [
    apply_root / "binned_models" / f"model_{name}_split_single_tmva.root"
    for name in apply_model_names
]
apply_models.append(apply_root / "npb_models" / "npb_score_split_tmva.root")
if len(apply_models) != 12 or any(not path.is_file() for path in apply_models):
    raise SystemExit("sealed apply_BDT stage does not contain exactly 11 split models plus NPB")
source_paths.extend(str(path) for path in apply_models)
if ppg_binary_mode not in ("source_locked_rebuild", "source_locked_runtime_import"):
    raise SystemExit(f"unsupported PPG12 binary mode: {ppg_binary_mode}")
if ppg_binary_mode == "source_locked_runtime_import":
    source_paths.extend((ppg_provenance_manifest, ppg_provenance_receipt))
log_paths = sorted(str(path) for path in Path(log_root).iterdir() if path.is_file())
source_runtime_import = None
if ppg_binary_mode == "source_locked_runtime_import":
    source_runtime_import = {
        "runtime_manifest": {
            "path": ppg_provenance_manifest,
            "sha256": ppg_origin_manifest_sha256,
        },
        "build_receipt": {
            "path": ppg_provenance_receipt,
            "sha256": ppg_origin_receipt_sha256,
        },
        "library_sha256": ppg_origin_library_sha256,
        "immutable_provenance_documents": True,
    }
data = {
    "schema_version": 1,
    "purpose": "ppg12_paired_oracle_recoil_runtime",
    "runtime_profile": "new.17",
    "offline_main": offline_main,
    "isolated_build": True,
    "output_root": output_root,
    "install_root": install_root,
    "runtime_root": runtime_root,
    "jobs": int(jobs),
    "sealed_link_directories": sealed_link_dirs.split(":"),
    "platform": platform.platform(),
    "calo_calib": {"path": calo_calib, "sha256": digest(calo_calib)},
    "ppg12_source": {
        "repository": ppg_repo,
        "revision": ppg_revision,
        "subtree": "anatreemaker/source",
        "working_tree_ignored": True,
        "rebuilt_against_common_runtime": True,
        "binary_mode": ppg_binary_mode,
        "rebuilt_in_this_runtime": ppg_binary_mode == "source_locked_rebuild",
        "source_runtime_import": source_runtime_import,
    },
    "ppg12_estimator": {
        "revision": estimator_revision,
        "subtree": "efficiencytool",
        "working_tree_ignored": True,
        "reconstruction_revision_is_distinct": True,
        "canonical_source": {
            "path": recoeff_source,
            "sha256": digest(recoeff_source),
        },
        "staged_uninstrumented_macro": {
            "path": recoeff_macro,
            "sha256": digest(recoeff_macro),
            "scientific_expression_changes": False,
            "sealed_path_rewrites_only": True,
        },
        "instrumented_macro": {
            "path": recoeff_trace_macro,
            "sha256": digest(recoeff_trace_macro),
            "transform_receipt": trace_receipt,
            "transform_receipt_sha256": digest(trace_receipt),
            "candidate_trace_side_channel_only": True,
            "requires_exact_root_equivalence": True,
        },
        "canonical_config": {
            "path": estimator_config,
            "sha256": digest(estimator_config),
        },
        "period_configs": {
            "0mrad": {
                "path": estimator_config_0mrad,
                "sha256": digest(estimator_config_0mrad),
            },
            "1p5mrad": {
                "path": estimator_config_1p5mrad,
                "sha256": digest(estimator_config_1p5mrad),
            },
        },
        "truth_vertex_reweights": {
            "0mrad": {
                "path": truth_vertex_reweight_0mrad,
                "sha256": digest(truth_vertex_reweight_0mrad),
            },
            "1p5mrad": {
                "path": truth_vertex_reweight_1p5mrad,
                "sha256": digest(truth_vertex_reweight_1p5mrad),
            },
        },
        "calculate_photon_yield": {
            "path": calculate_yield_source,
            "sha256": digest(calculate_yield_source),
        },
        "apply_bdt_stage": {
            "revision": estimator_revision,
            "working_tree_ignored": True,
            "macro": {"path": apply_bdt_source, "sha256": digest(apply_bdt_source)},
            "config": {"path": apply_config_source, "sha256": digest(apply_config_source)},
            "scored_root_shared_by_both_recoeff_runs": True,
            "model_assets": [
                {"path": str(path), "sha256": digest(path)} for path in apply_models
            ],
        },
        "yaml_cpp_header_tree_receipt": {
            "path": yaml_cpp_header_receipt,
            "sha256": digest(yaml_cpp_header_receipt),
        },
        "runtime_assets": [
            {"path": path, "sha256": digest(path)}
            for path in (
                cross_section_header, truth_vertex_header, yaml_cpp,
                yaml_cpp_header_receipt, roounfold, roounfold_pcm,
                roounfold_header_receipt, *roounfold_headers,
                vertex_scan_data, mbd_correction,
            )
        ],
    },
    "custom_photon_builder": {
        "class": "PPG12OraclePhotonClusterBuilder",
        "compiled_into": "libRecoilJets.so",
        "release_photon_builder_symbol_collision": False,
    },
    "sources": [
        {"path": path, "sha256": digest(path)} for path in source_paths
    ],
    "release_copies": {
        "libcalo_reco.so": {
            "source": calo_source,
            "source_realpath": os.path.realpath(calo_source),
            "sha256": digest(calo_source),
            "contains_release_calibration_and_cluster_template": True,
        },
        "libclusteriso.so": {
            "source": clusteriso_source,
            "source_realpath": os.path.realpath(clusteriso_source),
            "sha256": digest(clusteriso_source),
        },
        "libjetbase.so": {
            "source": jetbase_source,
            "source_realpath": os.path.realpath(jetbase_source),
            "sha256": digest(jetbase_source),
        },
    },
    "staged_rewrites": {
        "exact_count_enforced": True,
        "checkout_absolute_includes_removed": True,
        "macro_private_routes_removed": True,
        "staged_configure_werror_disabled": True,
        "target_linker_inherits_am_ldflags": True,
        "calo_truth_evaluator_linked_from_release": True,
        "calo_calib_resolved_by_sealed_release_path": True,
        "contained_recoiljets_tmp_paths_allowed": True,
        "full_calo_reco_rebuild": False,
        "custom_builder_renamed": True,
        "archived_ppg12_binary_reused": False,
        "ppg12_source_locked_rebuild": ppg_binary_mode == "source_locked_rebuild",
        "ppg12_source_locked_binary_import": ppg_binary_mode == "source_locked_runtime_import",
        "estimator_revision_separate_from_reconstruction": True,
        "estimator_trace_requires_exact_root_equivalence": True,
    },
    "validation": {
        "forbidden_route_scan": "pass",
        "readelf": "pass",
        "ldd": "pass",
        "root_load_and_header_smoke": "pass",
    },
    "logs": [{"path": path, "sha256": digest(path)} for path in log_paths],
}
Path(receipt).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY

python3 - "$runtime_manifest" "$build_receipt" "$expected_offline" \
  "$runtime_wrapper" "$runtime_impl" \
  "${runtime_root}/lib/libRecoilJets.so" \
  "${runtime_root}/lib/libCaloAna24.so" \
  "${runtime_root}/lib/libcalo_reco.so" \
  "${runtime_root}/lib/libclusteriso.so" \
  "${runtime_root}/lib/libjetbase.so" \
  "${runtime_root}/include/caloreco/PhotonClusterBuilder.h" \
  "$estimator_revision" \
  "${runtime_root}/estimator/source/RecoEffCalculator_TTreeReader.C" \
  "$recoeff_macro" "$recoeff_trace_macro" "$recoeff_trace_receipt" \
  "${runtime_root}/estimator/include/CrossSectionWeights.h" \
  "${runtime_root}/estimator/include/TruthVertexReweightLoader.h" \
  "${runtime_root}/estimator/config/config_bdt_nom.yaml" \
  "${runtime_root}/estimator/config/config_bdt_nom_0rad.yaml" \
  "${runtime_root}/estimator/config/config_bdt_nom_1p5mrad.yaml" \
  "${runtime_root}/estimator/source/CalculatePhotonYield.C" \
  "${runtime_root}/estimator/apply/apply_BDT.C" \
  "${runtime_root}/estimator/apply/config_nom.yaml" \
  "${runtime_root}/lib/libyaml-cpp.so" \
  "$yaml_cpp_header_receipt" \
  "${runtime_root}/lib/libRooUnfold.so" \
  "${runtime_root}/lib/RooUnfoldDict_rdict.pcm" \
  "$roounfold_header_receipt" \
  "${runtime_root}/estimator/include/RooUnfoldResponse.h" \
  "${runtime_root}/estimator/include/RooUnfoldBayes.h" \
  "${runtime_root}/estimator/data/data_histo_bdt_nom_vtxscan.root" \
  "${runtime_root}/estimator/data/MbdOut.corr" \
  "${runtime_root}/estimator/data/truth_vertex_reweight_0mrad.root" \
  "${runtime_root}/estimator/data/truth_vertex_reweight_1p5mrad.root" \
  "$ppg_provenance_manifest" "$ppg_provenance_receipt" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

(
    manifest, receipt, offline_main, macro, impl, recoil, ppg, calo,
    clusteriso, jetbase, photon_header,
    estimator_revision, recoeff_source, recoeff_macro, recoeff_trace_macro,
    trace_receipt, cross_section_header, truth_vertex_header, estimator_config,
    estimator_config_0mrad, estimator_config_1p5mrad, calculate_yield,
    apply_bdt, apply_config, yaml_cpp, yaml_cpp_header_receipt, roounfold,
    roounfold_pcm, roounfold_header_receipt, roounfold_response_header,
    roounfold_bayes_header, vertex_scan_data, mbd_correction,
    truth_vertex_reweight_0mrad, truth_vertex_reweight_1p5mrad,
    ppg_provenance_manifest, ppg_provenance_receipt,
) = sys.argv[1:]

def digest(path: str) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

roles = [
    ("recoil_macro", macro),
    ("recoil_impl", impl),
    ("libRecoilJets.so", recoil),
    ("libCaloAna24.so", ppg),
    ("libcalo_reco.so", calo),
    ("libclusteriso.so", clusteriso),
    ("libjetbase.so", jetbase),
    ("PhotonClusterBuilder.h", photon_header),
    ("ppg_recoeff_source_macro", recoeff_source),
    ("ppg_recoeff_macro", recoeff_macro),
    ("ppg_recoeff_trace_macro", recoeff_trace_macro),
    ("ppg_recoeff_trace_transform_receipt", trace_receipt),
    ("ppg_recoeff_cross_section_header", cross_section_header),
    ("ppg_recoeff_truth_vertex_header", truth_vertex_header),
    ("ppg_recoeff_canonical_config", estimator_config),
    ("ppg_recoeff_period_config_0mrad", estimator_config_0mrad),
    ("ppg_recoeff_period_config_1p5mrad", estimator_config_1p5mrad),
    ("ppg_calculate_photon_yield", calculate_yield),
    ("ppg_apply_bdt_macro", apply_bdt),
    ("ppg_apply_bdt_config", apply_config),
    ("ppg_recoeff_yaml_cpp", yaml_cpp),
    ("ppg_recoeff_yaml_cpp_header_tree_receipt", yaml_cpp_header_receipt),
    ("ppg_recoeff_roounfold", roounfold),
    ("ppg_recoeff_roounfold_pcm", roounfold_pcm),
    ("ppg_recoeff_roounfold_header_tree_receipt", roounfold_header_receipt),
    ("ppg_recoeff_roounfold_response_header", roounfold_response_header),
    ("ppg_recoeff_roounfold_bayes_header", roounfold_bayes_header),
    ("ppg_recoeff_vertex_scan_data", vertex_scan_data),
    ("ppg_recoeff_mbd_correction", mbd_correction),
    ("ppg_recoeff_truth_vertex_reweight_0mrad", truth_vertex_reweight_0mrad),
    ("ppg_recoeff_truth_vertex_reweight_1p5mrad", truth_vertex_reweight_1p5mrad),
]
apply_root = Path(apply_bdt).parent
model_names = (
    "base", "base_vr", "base_v0", "base_v1", "base_v2", "base_v3",
    "base_E", "base_v0E", "base_v1E", "base_v2E", "base_v3E",
)
roles.extend(
    (f"ppg_apply_model_{name}", str(apply_root / "binned_models" / f"model_{name}_split_single_tmva.root"))
    for name in model_names
)
roles.append(("ppg_apply_npb_model", str(apply_root / "npb_models" / "npb_score_split_tmva.root")))
if ppg_provenance_manifest or ppg_provenance_receipt:
    if not ppg_provenance_manifest or not ppg_provenance_receipt:
        raise SystemExit("source-runtime provenance roles must be supplied together")
    roles.extend((
        ("ppg_source_runtime_manifest", ppg_provenance_manifest),
        ("ppg_source_build_receipt", ppg_provenance_receipt),
    ))
data = {
    "schema_version": 1,
    "runtime_profile": "new.17",
    "offline_main": offline_main,
    "isolated_build": True,
    "estimator_revision": estimator_revision,
    "build_receipt": receipt,
    "build_receipt_sha256": digest(receipt),
    "files": [
        {"role": role, "path": path, "sha256": digest(path)}
        for role, path in roles
    ],
}
Path(manifest).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
PY

printf 'PASS\n' > "$state_file"
completed=1
trap - EXIT
printf 'PPG12_ORACLE_NEW17_BUILD_PASS output=%s manifest=%s receipt=%s\n' \
  "$output_dir" "$runtime_manifest" "$build_receipt"
