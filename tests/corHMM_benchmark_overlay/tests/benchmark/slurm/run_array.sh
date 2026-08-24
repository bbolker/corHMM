#!/bin/bash -l

set -euo pipefail

manifest=$1
offset=${2:-0}
row=$((offset + SLURM_ARRAY_TASK_ID))

cd "$HOME/corHMM"

module purge
module load r-rich/4.6.0_msi1.2

export BAYESTRAITS_HOME="$HOME/software/BayesTraitsV5.0.3-Linux"
export BAYESTRAITS_BIN="BayesTraitsV5"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export TBB_NUM_THREADS=1

mkdir -p logs

resource_base="logs/resource_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
cat > "${resource_base}.meta" <<EOF
array_job_id=${SLURM_ARRAY_JOB_ID}
array_task_id=${SLURM_ARRAY_TASK_ID}
slurm_job_id=${SLURM_JOB_ID}
manifest=${manifest}
row=${row}
host=$(hostname)
EOF

/usr/bin/time -v \
  -o "${resource_base}.txt" \
  Rscript tests/benchmark/run_job.R "$manifest" "$row"
