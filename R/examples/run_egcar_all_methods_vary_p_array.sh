#!/bin/bash
#SBATCH --job-name=egcar_all_p
#SBATCH --output=logs/egcar_all_p_%A_%a.out
#SBATCH --error=logs/egcar_all_p_%A_%a.err
#SBATCH --array=1-60%5
#SBATCH --time=35:00:00
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --account=pi-cdonnat

set -euo pipefail

module load R/4.2.0

export R_LIBS_USER="$HOME/Rlibs"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

ROOT="$SCRATCH/$USER/CCA-experiments"
RSCRIPT="$ROOT/r/experiments/cluster/run_egcar_all_methods_vary_p.R"

# This script is chunk-safe.  The complete experiment contains 180 tasks.
# Submit three 60-task chunks with offsets 0,60,120 and the same EGCAR_RUN_ID.
# A single full-array submission is also possible by overriding --array=1-180%5
# and leaving EGCAR_TASK_OFFSET=0.
TASK_OFFSET="${EGCAR_TASK_OFFSET:-0}"
RUN_ID="${EGCAR_RUN_ID:-job_${SLURM_ARRAY_JOB_ID}}"
OUTDIR="$ROOT/egcar_all_methods_vary_p_outputs/$RUN_ID"

mkdir -p "$OUTDIR"/{metrics,fits,cv,plots,logs}
cd "$ROOT"

local_task=${SLURM_ARRAY_TASK_ID}
global_task=$(( TASK_OFFSET + local_task ))
if (( global_task < 1 || global_task > 180 )); then
  echo "Invalid global task ${global_task}; expected 1..180" >&2
  exit 2
fi

# Ten replications per configuration.
config_id=$(( (global_task - 1) / 10 + 1 ))
rep_id=$(( (global_task - 1) % 10 + 1 ))

echo "JOB=${SLURM_JOB_ID} ARRAY_JOB=${SLURM_ARRAY_JOB_ID} LOCAL_TASK=${local_task}"
echo "offset=${TASK_OFFSET} global_task=${global_task} config_id=${config_id} rep_id=${rep_id}"
echo "cpus=${SLURM_CPUS_PER_TASK} output=${OUTDIR}"

Rscript "$RSCRIPT" check_packages
Rscript "$RSCRIPT" worker "$config_id" "$rep_id" "$OUTDIR"

# The last successful task to observe all 180 metric files performs aggregation.
expected=180
completed=$(find "$OUTDIR/metrics" -maxdepth 1 -type f -name 'metrics_*.csv' | wc -l | tr -d ' ')
echo "Completed metric files: ${completed}/${expected}"

if [[ "$completed" -eq "$expected" ]]; then
  if mkdir "$OUTDIR/.aggregation_lock" 2>/dev/null; then
    echo "All worker files found; aggregating."
    Rscript "$RSCRIPT" aggregate 0 0 "$OUTDIR"
  else
    echo "Another task already acquired the aggregation lock."
  fi
fi
