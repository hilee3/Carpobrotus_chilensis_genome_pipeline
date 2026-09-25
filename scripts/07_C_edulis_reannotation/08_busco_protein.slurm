#!/bin/bash
#SBATCH --job-name=busco_edulis_protein
#SBATCH --output=busco_edulis_protein_%j.out
#SBATCH --error=busco_edulis_protein_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# Required software:
#   BUSCO

# -----------------------------
# Inputs
# -----------------------------

PROTEINS="/path/to/edulis_braker3/braker.longest_isoform.gffread.aa"

LINEAGE="eudicotyledons_odb12"
OUT_PREFIX="busco_edulis_annotation_longest"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# -----------------------------
# Run BUSCO in protein mode
# -----------------------------

busco \
    -i "${PROTEINS}" \
    -l "${LINEAGE}" \
    -m protein \
    -o "${OUT_PREFIX}" \
    --cpu "${THREADS}"
