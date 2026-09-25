#!/usr/bin/env bash
#SBATCH --job-name=LAI
#SBATCH --output=LAI_%j.out
#SBATCH --error=LAI_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

# Required software:
#   LTR_retriever (LAI)
#   BLAST+

# -----------------------------
# Inputs
# -----------------------------

GENOME="/path/to/genome.fa"  # Copy of the final assembly (assembly_split.fa) used in the EDTA workflow.
INTACT_LTR="/path/to/genome.fa.mod.pass.list"
REPEATMASKER_OUT="/path/to/genome.fa.out"

OUTDIR="/path/to/output/lai"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# -----------------------------
# Calculate LTR Assembly Index
# -----------------------------

LAI \
    -genome "${GENOME}" \
    -intact "${INTACT_LTR}" \
    -all "${REPEATMASKER_OUT}" \
    -q \
    -t "${THREADS}"

echo "LAI analysis completed."
