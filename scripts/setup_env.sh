#!/bin/bash
# Build the ocean-uncertainty conda env on CURC Alpine. Run it by hand inside an interactive compile-node session,
# not on a login node and not through sbatch (the same practice as SAFARI-DA's setup_env.sh):
#
#   acompile
#   cd /projects/$USER/SAFARI-uncertainty-est
#   bash scripts/setup_env.sh          # ends with "ENV OK"
#
# Loading miniforge creates ~/.condarc if missing, which puts envs under /projects/$USER. pip's cache is redirected
# because its default (~/.cache/pip) is on the 2 GB /home quota.
set -eo pipefail   # no -u: conda's shell functions reference unset variables
cd "$(dirname "$0")/.."

export PIP_CACHE_DIR=/projects/$USER/.cache/pip
module purge
module load miniforge

ENV_NAME=ocean-uncertainty
if ! conda env list | grep -q "^$ENV_NAME "; then
    # Python 3.12: the version the code was developed and tested with (README).
    conda create -y -n "$ENV_NAME" -c conda-forge python=3.12 pip
fi
conda activate "$ENV_NAME"
python -m pip install -r requirements.txt -e .
python -m pytest -q tests
python -c "import copernicusmarine, argopy, gsw; print('ENV OK')"
