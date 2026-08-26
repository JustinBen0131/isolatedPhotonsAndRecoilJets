#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_dir=${1:-"${repo_dir}/outputs/example_pp"}
mkdir -p "${output_dir}"

export PYTHONPATH="${repo_dir}/python${PYTHONPATH:+:${PYTHONPATH}}"

tree_input="${repo_dir}/tests/fixtures/photonjet_trees_pp.root"

python3 -m photonjet.cli trees validate \
  --input "${tree_input}" \
  --model-input-count 11 \
  --report "${output_dir}/tree_validation.txt"

python3 -m photonjet.cli histogram \
  --input "${tree_input}" \
  --region A \
  --dataset "${repo_dir}/config/datasets/pp_engineering_fixture.yaml" \
  --output "${output_dir}/xjgamma.json"

python3 -m photonjet.cli skim \
  --input "${tree_input}" \
  --region A \
  --output "${output_dir}/recoil_pairs.root" \
  --receipt "${output_dir}/recoil_pairs.receipt.json"

python3 -m photonjet.cli purity \
  --input "${tree_input}" \
  --non-tight-definition bounded \
  --isolation-radius 0.4 \
  --output "${output_dir}/abcd_counts.json"

python3 -m photonjet.cli plot render \
  --histogram "${output_dir}/xjgamma.json" \
  --histogram-receipt "${output_dir}/xjgamma.json.receipt.json" \
  --dataset "${repo_dir}/config/datasets/pp_engineering_fixture.yaml" \
  --output "${output_dir}/xjgamma.png"

python3 -m photonjet.cli plot emit-root-header \
  --histogram "${output_dir}/xjgamma.json" \
  --histogram-receipt "${output_dir}/xjgamma.json.receipt.json" \
  --dataset "${repo_dir}/config/datasets/pp_engineering_fixture.yaml" \
  --output "${output_dir}/PlotAnnotation.generated.h"

printf 'Example outputs: %s\n' "${output_dir}"
