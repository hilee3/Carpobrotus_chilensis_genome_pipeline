#!/usr/bin/env bash
#SBATCH --job-name=tiara_screening
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=tiara.%j.out
#SBATCH --error=tiara.%j.err

set -euo pipefail

# -----------------------------
# Purpose
# -----------------------------
#
# Screen the final nuclear assembly using Tiara.
#
# Tiara classifies contigs into eukaryotic, organellar,
# prokaryotic, archaeal, or unknown categories.
# Contigs classified as residual organellar sequences
# were manually inspected and removed prior to the final
# assembly release.

# -----------------------------
# Inputs
# -----------------------------
# Edit these paths before running.
ASSEMBLY="/path/to/final_nuclear_assembly.fa.gz"
OUTDIR="tiara_output"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "${OUTDIR}"

# -----------------------------
# Software requirements
# -----------------------------
# Tiara

# -----------------------------
# Check dependency
# -----------------------------
command -v tiara >/dev/null || {
    echo "[ERROR] Tiara not found."
    exit 127
}

# -----------------------------
# Run Tiara
# -----------------------------

tiara \
    -i "${ASSEMBLY}" \
    -o "${OUTDIR}" \
    -t "${THREADS}" \
    --pr

echo "[DONE] Tiara screening completed."
