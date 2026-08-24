#!/bin/bash

set -euo pipefail

cd "$HOME/corHMM"
profile=${1:-full}
manifest="tests/benchmark/manifests/jobs_ml_${profile}.csv"

if [[ ! -f "$manifest" ]]; then
  echo "Missing $manifest"
  exit 1
fi

mkdir -p logs results/benchmark/simulation/ml
n_tasks=$(($(wc -l < "$manifest") - 1))
chunk_size=1000

for ((start=1; start<=n_tasks; start+=chunk_size)); do
  end=$((start + chunk_size - 1))
  if ((end > n_tasks)); then end=$n_tasks; fi
  count=$((end - start + 1))
  offset=$((start - 1))

  sbatch \
    --account=ventz001 \
    --partition=msismall \
    --job-name=corhmm_ml \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=4G \
    --time=04:00:00 \
    --array="1-${count}%50" \
    --output="logs/ml_%A_%a.out" \
    --error="logs/ml_%A_%a.err" \
    tests/benchmark/slurm/run_array.sh "$manifest" "$offset"
done
