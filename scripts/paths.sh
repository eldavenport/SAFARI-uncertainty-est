# Storage locations on CURC Alpine. Sourced, not executed.
# Raw downloads go on scratch, which CURC purges periodically (curc.readthedocs.io): rerun the fetch job if files go
# missing. Products are kept on /projects, which is not purged.
DATA_ROOT=/scratch/alpine/$USER/ocean_uncertainty/data
PRODUCT_ROOT=/projects/$USER/ocean_uncertainty/products
