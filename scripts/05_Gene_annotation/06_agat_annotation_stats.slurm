#!/bin/bash
#SBATCH --job-name=agat_annotation_stats
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=agat_annotation_stats_%j.out
#SBATCH --error=agat_annotation_stats_%j.err

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Complete BRAKER3 annotation.
GFF="/path/to/braker.gff3"
OUTDIR="/path/to/output/agat_stats"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# -----------------------------
# Software requirements
# -----------------------------
# AGAT must be available in PATH.

[[ -r "$GFF" ]] || {
    echo "[ERR] GFF not found: $GFF"
    exit 1
}

command -v agat_sp_statistics.pl >/dev/null || {
    echo "[ERR] agat_sp_statistics.pl not in PATH"
    exit 1
}

echo "[INFO] Running AGAT statistics on: $GFF"
date

agat_sp_statistics.pl \
    --gff "$GFF" \
    -o braker_annotation_stats.txt

echo "[INFO] Done"
date
echo "[INFO] Output: $OUTDIR/braker_annotation_stats.txt"