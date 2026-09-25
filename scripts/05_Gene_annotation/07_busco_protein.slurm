#!/bin/bash
#SBATCH --job-name=busco_protein
#SBATCH --output=busco_protein_%j.out
#SBATCH --error=busco_protein_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail

# Required software:
#   BUSCO

# -----------------------------
# Inputs
# -----------------------------

PROTEINS="/path/to/braker.g103t2_removed.longest_isoform.gffread.aa"
LINEAGE="eudicotyledons_odb12"
OUT_PREFIX="busco_braker_annotation_longest"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

# -----------------------------
# Run BUSCO in protein mode
# -----------------------------

busco \
    -i "${PROTEINS}" \
    -l "${LINEAGE}" \
    -m protein \
    -o "${OUT_PREFIX}" \
    --cpu "${THREADS}"
