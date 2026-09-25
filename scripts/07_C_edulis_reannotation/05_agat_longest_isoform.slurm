#!/bin/bash
#SBATCH --job-name=agat_longest_edulis
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=agat_longest_edulis_%j.out
#SBATCH --error=agat_longest_edulis_%j.err

set -euo pipefail

# Required software:
#   AGAT

# -----------------------------
# Inputs
# -----------------------------

BRAKER_DIR="/path/to/edulis_braker3"

GFF3="${BRAKER_DIR}/braker.gff3"
OUT="${BRAKER_DIR}/braker.longest_isoform.agat.gff3"

agat_sp_keep_longest_isoform.pl \
    --gff "${GFF3}" \
    -o "${OUT}"

echo "[INFO] Output: ${OUT}"
