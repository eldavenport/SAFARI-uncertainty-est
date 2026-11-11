#!/bin/bash
# Download the sources for one winter per array task; the task ID is the winter's first year.
#
#   cd /projects/$USER/SAFARI-uncertainty-est && mkdir -p logs
#   sbatch --array=2014 scripts/fetch.sbatch          # pilot winter
#   sbatch --array=2011-2021 scripts/fetch.sbatch     # all winters (finished ones are skipped)
#
# Needs outbound internet on the compute node (SAFARI-DA's prefetch jobs download from compute nodes) and
# Copernicus Marine credentials: `copernicusmarine login` once, or ~/.copernicusmarine_env (chmod 600) exporting
# COPERNICUSMARINE_SERVICE_USERNAME and COPERNICUSMARINE_SERVICE_PASSWORD. Never in this repository.
#
# Size: 2.8 GB per winter, 83 MB for WOA23 (measured on the laptop for 2014-15). Memory: the largest array held is
# one month of OSTIA (0.05 degree, ~0.4 GB in double precision) plus one global WOA23 file being cut (~0.6 GB);
# 8G leaves room. Time: not measured on Alpine; record sacct Elapsed from the pilot and tighten.
#SBATCH --job-name=ou-fetch
#SBATCH --account=ucb-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x.%A_%a.out

set -eo pipefail
cd "${SLURM_SUBMIT_DIR}"
module purge
module load miniforge
conda activate ocean-uncertainty
[[ -f ~/.copernicusmarine_env ]] && source ~/.copernicusmarine_env

scripts/run_winter.sh fetch "${SLURM_ARRAY_TASK_ID}"
