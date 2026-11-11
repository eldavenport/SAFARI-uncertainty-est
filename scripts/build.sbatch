#!/bin/bash
# Build U1-U6 and the checks for one winter per array task (task ID = the winter's first year), after its fetch.
#
#   sbatch --array=2014 --dependency=aftercorr:<fetch job id> scripts/build.sbatch
#
# aftercorr starts each task once the fetch task with the same ID has succeeded.
#
# Time on the laptop (2014-15, 14 cores): U1-U5 2 min, checks report 14 s, hold-out of the six error models 7 min,
# U6 2 min per window x 21 windows; about 55 min in all. The U6 and hold-out linear algebra use threaded BLAS,
# hence 8 CPUs. Memory: not measured; the largest pieces are one winter of OSTIA (~2 GB in double precision, U4)
# and U6 blocks of 4000 state points x ~5000 observations (~0.2 GB each, several at once); 32G leaves room. Record
# sacct MaxRSS and Elapsed from the pilot and tighten both.
#SBATCH --job-name=ou-build
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x.%A_%a.out

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
module purge
module load miniforge
conda activate ocean-uncertainty
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

scripts/run_winter.sh build "${SLURM_ARRAY_TASK_ID}"
