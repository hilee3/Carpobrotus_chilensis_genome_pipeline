#!/bin/bash
#SBATCH --job-name=agat_braker_stats
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=agat_braker_stats_%j.out
#SBATCH --error=agat_braker_stats_%j.err

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.

# Curated longest-isoform BRAKER3 annotation.
# One transcript model (g103.t2) containing an internal stop codon
# was removed prior to longest isoform selection.
GFF="/path/to/braker.g103t2_removed.longest_isoform.agat.gff3"

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
    -o braker_agat_stats.txt

echo "[INFO] Done"
date
echo "[INFO] Output: $OUTDIR/braker_agat_stats.txt"