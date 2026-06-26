#!/usr/bin/env bash
#SBATCH --job-name=busco
#SBATCH --output=busco.%j.out
#SBATCH --error=busco.%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
ASSEMBLY="/path/to/final_nuclear_assembly.fa"

# BUSCO lineage
LINEAGE="eudicotyledons_odb12"

# Output prefix
OUT_PREFIX="busco"

# -----------------------------
# Software requirements
# -----------------------------
# BUSCO

busco \
    -i "${ASSEMBLY}" \
    -l "${LINEAGE}" \
    -m genome \
    -o "${OUT_PREFIX}" \
    --cpu "${SLURM_CPUS_PER_TASK}"
