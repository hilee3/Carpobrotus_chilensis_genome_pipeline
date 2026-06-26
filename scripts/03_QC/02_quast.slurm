#!/usr/bin/env bash
#SBATCH --job-name=quast
#SBATCH --output=quast.%j.out
#SBATCH --error=quast.%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
ASSEMBLY="/path/to/final_nuclear_assembly.fa"
OUTDIR="quast_results"

# -----------------------------
# Software requirements
# -----------------------------
# QUAST

quast.py "${ASSEMBLY}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --output-dir "${OUTDIR}" \
    --large \
    --k-mer-stats
