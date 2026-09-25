#!/usr/bin/env bash
#SBATCH --job-name=btk_busco
#SBATCH --output=btk_busco_%j.out
#SBATCH --error=btk_busco_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail

# Required software:
#   BlobToolKit (blobtools)

# -----------------------------
# Inputs
# -----------------------------

WORKDIR="/path/to/blobtools"

BLOBDIR="${WORKDIR}/chilensis_blobdir"
BUSCO_TABLE="/path/to/run_eudicotyledons_odb12/full_table.tsv"

cd "${WORKDIR}"


# -----------------------------
# Add BUSCO results
# -----------------------------

blobtools add \
    --busco "${BUSCO_TABLE}" \
    "${BLOBDIR}"

blobtools validate \
    --basic \
    "${BLOBDIR}"
