#!/bin/bash
# Runs Snakemake as a Slurm job on Triton. Snakemake then submits the
# workflow jobs. Submit from the pipeline folder:
#   sbatch snakemake_slurm.sh                            (full grid)
#   sbatch snakemake_slurm.sh --configfile config/small.yaml   (smoke test)
#SBATCH --job-name=loocv_snakemake
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G                   # the local rule combine holds all trials
#SBATCH --output=logs/snakemake_slurm-%j.out
#SBATCH --error=logs/snakemake_slurm-%j.err

# One BLAS thread per R process: the matrices have 3 columns, and more
# threads than Slurm CPUs slow the jobs down.
export OPENBLAS_NUM_THREADS=1

module load mamba
# Absolute path: mamba resolves a relative prefix inside envs_dirs.
source activate "$SLURM_SUBMIT_DIR/env"
snakemake --profile profiles/triton --cores 2 "$@"
