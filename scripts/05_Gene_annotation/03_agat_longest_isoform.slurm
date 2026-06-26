#!/bin/bash
#SBATCH --job-name=agat_longest
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=agat_longest_%j.out
#SBATCH --error=agat_longest_%j.err

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
BRAKER_DIR="/path/to/braker3_output"

# Curated BRAKER3 annotation.
# One transcript model (g103.t2) containing an internal stop codon
# was removed prior to longest isoform selection.
GFF3="${BRAKER_DIR}/braker.g103t2_removed.gff3"
OUT="${BRAKER_DIR}/braker.g103t2_removed.longest_isoform.agat.gff3"

# -----------------------------
# Software requirements
# -----------------------------
# AGAT must be available in PATH.

command -v agat_sp_keep_longest_isoform.pl >/dev/null || {
    echo "[ERR] AGAT not found in PATH"
    exit 1
}

[[ -r "$GFF3" ]] || {
    echo "[ERR] Input GFF3 not found: $GFF3"
    exit 1
}

agat_sp_keep_longest_isoform.pl \
    --gff "$GFF3" \
    -o "$OUT"

echo "[INFO] Output: $OUT"
