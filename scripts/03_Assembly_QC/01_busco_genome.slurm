#!/usr/bin/env bash
#SBATCH --job-name=busco
#SBATCH --output=busco.%j.out
#SBATCH --error=busco.%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# -----------------------------
# Inputs
# -----------------------------
# Final nuclear assembly
ASSEMBLY="/path/to/assembly_split.fa"

# BUSCO lineage
LINEAGE="eudicotyledons_odb12"

# Output prefix
OUT_PREFIX="busco"

# -----------------------------
# Run BUSCO
# -----------------------------
busco \
    -i "${ASSEMBLY}" \
    -l "${LINEAGE}" \
    -m genome \
    -o "${OUT_PREFIX}" \
    --cpu "${SLURM_CPUS_PER_TASK}"