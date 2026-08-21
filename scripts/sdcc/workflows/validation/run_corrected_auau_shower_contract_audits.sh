#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 EXTRACTION_ROOT OUTPUT_DIR [EXPECTED_FILES_PER_SAMPLE]" >&2
  exit 64
fi

source_root="$1"
output_dir="$2"
expected_files="${3:-715}"
repo_base="${RJ_REPO_BASE:-/sphenix/u/patsfan753/scratch/thesisAnalysis}"
ml_python="${RJ_ML_PYTHON:-/sphenix/u/patsfan753/.venvs/thesis-ml/bin/python}"
label_audit="${output_dir}/label_and_source_audit.json"
shower_audit="${output_dir}/shower_contract_audit.json"
acceptance="${output_dir}/corrected_extraction_acceptance.json"

[[ -d "$source_root" ]] || { echo "Missing extraction root: $source_root" >&2; exit 66; }
[[ -x "$ml_python" ]] || { echo "Missing ML Python: $ml_python" >&2; exit 69; }
mkdir -p "$output_dir"

"$ml_python" -B "${repo_base}/scripts/ml/validation/audit_the107_auau_bdt_extraction.py" \
  --root "$source_root" \
  --output "$label_audit" \
  --expected-files-per-sample "$expected_files" \
  --workers 8

"$ml_python" -B "${repo_base}/scripts/ml/validation/audit_corrected_auau_shower_contract_roots.py" \
  --root "$source_root" \
  --output "$shower_audit" \
  --expected-files-per-sample "$expected_files" \
  --workers 8

"$ml_python" -B - "$label_audit" "$shower_audit" "$acceptance" <<'PY'
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


label_path, shower_path, output_path = map(Path, sys.argv[1:])
label = json.loads(label_path.read_text())
shower = json.loads(shower_path.read_text())
passed = label.get("status") == "PASSED" and shower.get("status") == "PASSED"
same_files = label.get("total_files") == shower.get("total_files")
same_rows = label.get("total_rows") == shower.get("total_rows")
payload = {
    "schema": "CORRECTED_AUAU_SHOWER_CONTRACT_EXTRACTION_ACCEPTANCE_V1",
    "status": "PASSED" if passed and same_files and same_rows else "FAILED",
    "promotion_status": "NOT_PROMOTED",
    "extraction_root": label.get("root"),
    "total_files": label.get("total_files"),
    "total_rows": label.get("total_rows"),
    "gates": {
        "truth_label_source_audit_passed": label.get("status") == "PASSED",
        "shower_contract_audit_passed": shower.get("status") == "PASSED",
        "audits_cover_identical_file_count": same_files,
        "audits_cover_identical_row_count": same_rows,
    },
    "inputs": {
        "label_audit": str(label_path),
        "label_audit_sha256": sha256(label_path),
        "shower_audit": str(shower_path),
        "shower_audit_sha256": sha256(shower_path),
    },
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(json.dumps({"status": payload["status"], "output": str(output_path)}, sort_keys=True))
raise SystemExit(0 if payload["status"] == "PASSED" else 2)
PY

echo "CORRECTED_AUAU_SHOWER_CONTRACT_EXTRACTION_AUDIT_COMPLETE"
echo "acceptance=${acceptance}"
