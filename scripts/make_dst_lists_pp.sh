#!/bin/bash
# make_dst_lists_pp.sh — generate per-run DST lists (no side files)
# Production: ana509_2024p022_v001
# DST type   : DST_CALOFITTING

set -uo pipefail  # (intentionally NOT using -e)

IN_BASE="/sphenix/u/patsfan753/scratch/thesisAnalysis"
LIST_FILE="$IN_BASE/Full_ppGoldenRunList_Version3.list"
OUT_DIR="$IN_BASE/dst_lists_pp"

TAG="ana521_2025p007_v001"
TYPE="DST_JETCALO"

CREATE_DST_TOOL="$(command -v CreateDstList.pl || true)"
[[ -n "$CREATE_DST_TOOL" ]] || { echo "[ERROR] CreateDstList.pl not in PATH"; exit 1; }
[[ -f "$LIST_FILE" ]]      || { echo "[ERROR] Run list not found: $LIST_FILE"; exit 1; }
[[ -s "$LIST_FILE" ]]      || { echo "[ERROR] Run list is empty: $LIST_FILE"; exit 1; }

mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit 1

echo "Tag: $TAG"
echo "Type: $TYPE"
echo "Input: $LIST_FILE"
echo "Output dir: $OUT_DIR"
echo "Tool: $CREATE_DST_TOOL"
echo "----------------------------------------"

total=0
made=0

# Always start fresh: wipe previous output list files in this directory
rm -f dst_*.list DST_*.list

while IFS= read -r raw; do
  # strip comments and whitespace
  run="${raw%%#*}"
  run="$(echo -n "$run" | tr -d '[:space:]')"
  [[ -z "$run" ]] && continue

  ((total++))
  # normalize number, keep a zero-padded version for filename check
  runnum=$(printf "%d" "$run" 2>/dev/null || echo "$run")
  pad=$(printf "%08d" "$runnum")

  # Call the tool exactly as when it worked, just changing the DST type
  perl "$CREATE_DST_TOOL" --tag "$TAG" --run "$runnum" "$TYPE"
  rc=$?

  # Accept success even if tool is noisy; verify file existence (both case variants)
  if [[ -s "dst_jetcalo-$pad.list" || -s "DST_JETCALO-$pad.list" ]]; then
    ((made++))
  else
    echo "[WARN] No list file found (or file empty) for run $runnum (exit=$rc)"
  fi
done < "$LIST_FILE"

echo "----------------------------------------"
echo "Requested runs: $total"
echo "List files made: $made"
echo "Done. Files are in: $OUT_DIR"
