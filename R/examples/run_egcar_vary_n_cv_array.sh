#!/bin/bash
#SBATCH --job-name=egcar_n_cv
#SBATCH --output=logs/egcar_n_cv_%A_%a.out
#SBATCH --error=logs/egcar_n_cv_%A_%a.err
#SBATCH --array=1-360%10
#SBATCH --time=35:00:00
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --account=pi-cdonnat

set -euo pipefail

module load R/4.2.0
export R_LIBS_USER="$HOME/Rlibs"

# Each R worker uses forked processes. Keep every BLAS process single-threaded
# so 5 CV folds do not each start an additional multithreaded BLAS pool.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Fixed simulation and algorithm settings. The signal parameter lambda is not
# the sparsity penalty; the latter is called rho in the R code.
export EGCAR_SUPPORT_SIZE=5
export EGCAR_SIGNAL_LAMBDA=0.80
export EGCAR_TOEPLITZ_ALPHA=0.30
export EGCAR_THEORY_RHO_CONSTANT=1.0
export EGCAR_CV_FOLDS=5
export EGCAR_ADMM_MU=1.0
export EGCAR_ADMM_MAX_ITER=2000
export EGCAR_ADMM_ABS_TOL=1e-5
export EGCAR_ADMM_REL_TOL=1e-4
export EGCAR_COVARIANCE_RIDGE=1e-4
export EGCAR_ORACLE_RIDGE=1e-8
export EGCAR_SCREEN_CAP=200
export EGCAR_SCREEN_SQRT_MULTIPLIER=2.5

ROOT="$SCRATCH/$USER/CCA-experiments"
SCRIPT="$ROOT/r/experiments/cluster/run_egcar_vary_n_cv.R"
OUTDIR="$ROOT/egcar_vary_n_cv_outputs/job_${SLURM_ARRAY_JOB_ID}"
EXPECTED_TASKS=360

mkdir -p "$ROOT/logs" "$OUTDIR"/{metrics,cv,snapshots,plots,config}
cd "$ROOT"

# 36 configurations = 2 p settings x 3 ranks x 6 n values.
# Each configuration has 10 independent Monte Carlo repetitions.
task=${SLURM_ARRAY_TASK_ID}
config_id=$(( (task - 1) / 10 + 1 ))
rep_id=$(( (task - 1) % 10 + 1 ))

echo "JOB=${SLURM_JOB_ID} ARRAY_JOB=${SLURM_ARRAY_JOB_ID} TASK=${task}"
echo "config_id=${config_id} rep_id=${rep_id} cpus=${SLURM_CPUS_PER_TASK}"
echo "outdir=${OUTDIR}"

Rscript "$SCRIPT" "$config_id" "$rep_id" "$OUTDIR"

# The metrics CSV is written last by each worker. The last successful worker
# acquires this directory lock and creates all mean-over-10 summary plots.
completed=$(find "$OUTDIR/metrics" -maxdepth 1 -type f \
  -name 'metrics_c*_rep*.csv' | wc -l | tr -d ' ')
echo "completed=${completed}/${EXPECTED_TASKS}"

if [[ "$completed" -eq "$EXPECTED_TASKS" ]]; then
  if mkdir "$OUTDIR/.aggregate_lock" 2>/dev/null; then
    echo "All workers completed; aggregating."
    Rscript "$SCRIPT" aggregate 0 "$OUTDIR"
    touch "$OUTDIR/AGGREGATION_COMPLETE"
  else
    echo "Another task already acquired the aggregation lock."
  fi
fi
